#!/usr/bin/python3

import sys
from itertools import zip_longest
import re
import os.path
from operator import itemgetter
from itertools import groupby
from itertools import islice

miRBasestr = "reference/miRBase_miRNA_v21.str"
miRBaseHairpin = "reference/miRBase_hairpin_v21.fa"
miRBaseMature = "reference/miRBase_mature_v21.fa"
cut_site_variation="data/cut_site_variation_distribution_in_specific_regions_with_sh-lib_and_extended_out_hairpins_as_reference.txt"

# the mapping profile with extended hairpin out sequences, at least 118 nt. 
mapping_profile="data/start-end-seq-by-hairpin.all.sumbysample_with_sh-lib_and_extended_out_hairpins_as_reference"

class RNAStructure:

    def __init__(self, vienna):        

        lines = []
        for line in open(vienna):
            line = line.strip()
            lines.append(line)
        
        # resolves haripin secondary structure
        self.name = lines[0].strip() 
        self.name = re.sub('/.*RNAfold_','',self.name)
        self.sequence = lines[1].strip() 
        self.sequence_with_U = self.sequence.upper()
        self.sequence = self.sequence.upper().replace('U','T')
        self.dotbracket = lines[2].strip()
        self.energe = self.dotbracket.split(' ')[1]
        self.dotbracket = self.dotbracket.split(' ')[0]
        self.basepairs = sorted(self.parse_basepairs(self.dotbracket))
        self.annotationDict = self.parse_annotation(self.dotbracket)
        (self.index, self.annotation) = zip(*self.annotationDict.items())
        self.index = list(map(str, self.index))
        self.longesthairpin = sorted(self.longest_hairpin(self.dotbracket))
        self.longesthairpinwithgap = self.parse_longesthairpin_withgap(self.longesthairpin,self.sequence)
        self.structure = self.parse_structurefolding(self.longesthairpin,self.sequence)
        self.dotbracket_corrected = '.'*len(self.dotbracket[:self.longesthairpin[0][0]]) + self.dotbracket[self.longesthairpin[0][0]:self.longesthairpin[-1][0]+1] + '.'*len(self.dotbracket[self.longesthairpin[-1][0]+1:self.longesthairpin[-1][1]]) + self.dotbracket[self.longesthairpin[-1][1]:self.longesthairpin[0][1]+1] + '.'*len(self.dotbracket[self.longesthairpin[0][1]+1:])
        
        ## need to update the reference information
        self.ref = self.parse_ref(self.name, self.sequence, miRBaseHairpin, miRBaseMature, miRBasestr)
        self.refstructure = self.ref[3]
        self.refhairpin = self.ref[0]
        self.refmirna = self.ref[1]

        
        ## old cutting profile (rewrite)
        self.cutposition = self.parse_cut_positions(self.name, self.sequence,cut_site_variation)
        self.sequence_with_cut_position = self.cutposition[-1]
        if (type(self.ref[-1]) == str):
            self.cutposition = [self.ref[-1]] + self.cutposition[0:len(self.cutposition)-1]
        elif (type(self.cutposition[0:len(self.cutposition)-1]) == str):
            self.cutposition = self.ref[-1] + [self.cutposition[0:len(self.cutposition)-1]]
        else:
            self.cutposition = self.ref[-1] + self.cutposition[0:len(self.cutposition)-1]

        self.droshaORdicer_cutposition = self.parse_dominant_cut_position(self.refmirna, self.cutposition)
        
        
        ## new version of mapping profile
        self.mappingprofile = self.parse_mapping_profile(self.name, self.sequence)
        self.mappingprofile_statistics = self.parse_beginning_of_mapping_profile(self.mappingprofile, self.longesthairpinwithgap, self.ref, self.sequence, self.name)
        
        self.mappingprofile_print = []
        for item in self.mappingprofile:
            tmp = list(item)
            tmp[4] = str(tmp[4])
            tmp[5] = str(tmp[5])
            tmp[6] = str(tmp[6])
             
            self.mappingprofile_print.append(tmp[0] + ' ' + ';'.join((tmp[1:7]))) 

        ## bulge profile (need to write, have a look the python scirpt)
        ## motif profile (need to write, have a look the python script)

        
        
        
        
    #############
    # functions #
    #############
        
        
    def parse_basepairs(self, dotbracket): 
        stack = []
        for i, char in enumerate(dotbracket): 
            if char == '(':
                stack.append(i) 
            elif char == ')':
                j = stack.pop() 
                yield j,i


    def longest_hairpin(self, dotbracket):
        stack = []
        pairsindex = []
        maxium = 0
        state = False
        tmp = 0
        for i, char in enumerate(dotbracket):
            if char == '(':
                if state == True:
                    tmp = len(stack)
                    if maxium < len(pairsindex):
                        maxium = len(pairsindex)
                        output = pairsindex                        
                    pairsindex = []

                stack.append(i)
                state = False

            elif char == ')':
                state = True
                j = stack.pop()
                pairsindex.append((j,i))
                if len(stack) == tmp:
                    if maxium < len(pairsindex):
                        maxium = len(pairsindex)
                        output = pairsindex
                    pairsindex = []
                elif len(stack) < tmp:
                    tmp = len(stack)
                    if tmp == 0:
                        if maxium < len(pairsindex):
                            maxium = len(pairsindex)
                            output = pairsindex
                        pairsindex = []

        for i in output:
            yield i


    def parse_annotation(self, dotbracket):
        output = dict()
        tmp = []
        left = False
        for i, char in enumerate(dotbracket):
            if char == '(': 
                left = True
                output[i] = 's'
                tmp = []
            elif char == ')':
                left = False
                output[i] = 's'
                if len(tmp) != 0:
                    for n in tmp:
                        output[n] = 'h' 
                tmp = []
            elif char == '.' and left == True:
                tmp.append(i)
                output[i] = 'm'
            elif char == '.' and left == False:
                output[i] = 'm'

            elif char == ' ':
                break

        return output


    def parse_longesthairpin_withgap(self,longesthairpin,sequence):
        pairswithgap = []
        count = 0
        for i,j in longesthairpin:
            if count == 0:
                p1 = i
                p2 = j
                pairswithgap.append((p1,sequence[p1],'-',sequence[p2],p2))
                count = count + 1

            else:
                if i - p1 > p2 - j:
                    n = i - p1 - (p2 - j) + 1
                    m = p2 - j
                    for nt in range(1,n):
                        pairswithgap.append((p1 + nt, sequence[p1 + nt],'>','|', p2))
                    for nt in range(1,m):
                        pairswithgap.append((p1 + nt + (n - 1), sequence[p1 + nt + (n - 1)],'>',sequence[p2 - nt], p2 - nt))

                elif i - p1 < p2 - j:
                    n = p2 - j - (i - p1) + 1
                    m = i - p1
                    for nt in range(1,n):
                        pairswithgap.append((p1,'|','<', sequence[p2 - nt], p2 - nt))
                    for nt in range(1,m):
                        pairswithgap.append((p1 + nt, sequence[p1 + nt],'<', sequence[p2 - nt - (n -1)], p2 - nt - (n - 1)))

                elif i - p1 == p2 - j:
                    n = i - p1
                    for nt in range(1,n):
                        pairswithgap.append((p1 + nt, sequence[p1 + nt], '*', sequence[p2 - nt], p2 - nt))

                pairswithgap.append((i,sequence[i],'-',sequence[j],j))

                p1 = i
                p2 = j

        if len(sequence[:longesthairpin[0][0]]) > len(sequence[longesthairpin[0][1]+1:]):
            diff = len(sequence[:longesthairpin[0][0]])-len(sequence[longesthairpin[0][1]+1:])
            f5 = list(range(0,longesthairpin[0][0]))
            f3 = list('.' * diff) + list(reversed(list(range(longesthairpin[0][1]+1,len(sequence)))))
        else:
            diff = len(sequence[longesthairpin[0][1]+1:]) - len(sequence[:longesthairpin[0][0]])
            f5 = list('.' * diff) + list(range(0,longesthairpin[0][0]))
            f3 = list(reversed(list(range(longesthairpin[0][1]+1,len(sequence)))))

        flanking = []
        for k in range(len(f5)):
            if f5[k] != '.' and f3[k] != '.':
                flanking.append((f5[k],sequence[f5[k]],'*',sequence[f3[k]],f3[k]))
            elif f5[k] == '.' and  f3[k] != '.':
                flanking.append((f5[k],'|','*',sequence[f3[k]],f3[k]))
            elif f5[k] != '.' and f3[k] == '.':
                flanking.append((f5[k],sequence[f5[k]],'*','|',f3[k]))
            else:
                flanking.append((f5[k],'|','*','|',f3[k]))
            
        pairswithgap = flanking + pairswithgap

        loop = []
        loopindex = longesthairpin[len(longesthairpin)-1][1]-longesthairpin[len(longesthairpin)-1][0] - 1
       
        sep = int(loopindex/2)
        for nt in range(1,sep+1):
            loop.append((longesthairpin[len(longesthairpin)-1][0] + nt, sequence[longesthairpin[len(longesthairpin)-1][0] + nt],'*', sequence[longesthairpin[len(longesthairpin)-1][1] - nt],longesthairpin[len(longesthairpin)-1][1] - nt))

        if loopindex % 2 != 0:
            loop.append((longesthairpin[len(longesthairpin)-1][0] + sep + 1, sequence[longesthairpin[len(longesthairpin)-1][0] + sep + 1],'*', sequence[longesthairpin[len(longesthairpin)-1][0] + sep + 1], longesthairpin[len(longesthairpin)-1][0] + sep + 1))

        pairswithgap = pairswithgap + loop

        return pairswithgap
    


    def parse_structurefolding(self,longesthairpin,sequence):
        
        line1 = []
        line2 = []
        line3 = []
        line4 = []
        line5 = []

        count = 0
        for i,j in longesthairpin:
            if count == 0:
                p1 = i
                p2 = j
                line1.extend(' ')
                line2.extend(sequence[p1])
                line3.extend('|')
                line4.extend(sequence[p2])
                line5.extend(' ')
                count = count + 1
            else:
                if i - p1 > p2 - j:
                    n = i - p1 - (p2 - j) + 1
                    m = p2 - j
                    for nt in range(1,n):
                        line1.extend(sequence[p1 + nt])
                        line2.extend(' ')
                        line3.extend(' ')
                        line4.extend(' ')
                        line5.extend('-')
                    for nt in range(1,m):
                        line1.extend(sequence[p1 + nt + (n - 1)])
                        line2.extend(' ')
                        line3.extend(' ')
                        line4.extend(' ')
                        line5.extend(sequence[p2 - nt])

                elif i - p1 < p2 - j:
                    n = p2 - j - (i - p1) + 1
                    m = i - p1
                    for nt in range(1,n):
                        line1.extend('-')
                        line2.extend(' ')
                        line3.extend(' ')
                        line4.extend(' ')
                        line5.extend(sequence[p2 - nt])
                    for nt in range(1,m):
                        line1.extend(sequence[p1 + nt])
                        line2.extend(' ')
                        line3.extend(' ')
                        line4.extend(' ')
                        line5.extend(sequence[p2 - nt - (n - 1)])


                elif i - p1 == p2 - j:
                    n = i - p1
                    for nt in range(1,n):
                        line1.extend(sequence[p1 + nt])
                        line2.extend(' ')
                        line3.extend(' ')
                        line4.extend(' ')
                        line5.extend(sequence[p2 - nt])

                line1.extend(' ')
                line2.extend(sequence[i])
                line3.extend('|')
                line4.extend(sequence[j])
                line5.extend(' ')
        
                p1 = i
                p2 = j


        loopindex = longesthairpin[len(longesthairpin)-1][1]-longesthairpin[len(longesthairpin)-1][0] - 1

        sep = int(loopindex/2)
        for nt in range(1,sep+1):
            line1.extend(sequence[longesthairpin[len(longesthairpin)-1][0] + nt])
            line2.extend(' ')
            line3.extend(' ')
            line4.extend(' ')
            line5.extend(sequence[longesthairpin[len(longesthairpin)-1][1] - nt])

        if loopindex % 2 != 0:
            line1.extend(' ')
            line2.extend(' ')
            line3.extend(sequence[longesthairpin[len(longesthairpin)-1][0] + sep + 1])
            line4.extend(' ')
            line5.extend(' ')


        if len(sequence[:longesthairpin[0][0]]) > len(sequence[longesthairpin[0][1]+1:]):
            diff = len(sequence[:longesthairpin[0][0]])-len(sequence[longesthairpin[0][1]+1:])
            f5 = sequence[:longesthairpin[0][0]]
            f3 = ''.join(('.' * diff, ''.join(reversed(sequence[longesthairpin[0][1]+1:]))))
        else:
            diff = len(sequence[longesthairpin[0][1]+1:]) - len(sequence[:longesthairpin[0][0]])
            f5 = ''.join(('.' * diff,sequence[:longesthairpin[0][0]]))
            f3 = ''.join(reversed(sequence[longesthairpin[0][1]+1:]))

        out = []
        out.append(f5 + ''.join(line1))
        out.append(' '*len(f5) + ''.join(line2))
        out.append(' '*len(f5) + ''.join(line3))
        out.append(' '*len(f5) + ''.join(line4))
        out.append(f3 + ''.join(line5))

        return out



    def parse_ref(self, name, sequence, miRBaseHairpin, miRBaseMature, miRBasestr):
        
    
        def fasta_iter(fasta_name):
            """
            given a fasta file. yield tuples of header, sequence
            """
            fh = open(fasta_name)
            faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
            for header in faiter:
                header = next(header)[1:].split(' ')[0]
                seq = ''.join(s.strip() for s in next(faiter))
                yield header, seq
            
        out = []

        fastahairpin = dict()
        for header,seq in fasta_iter(miRBaseHairpin):
            if header in fastahairpin:
                pass
            else:
                fastahairpin[header]=seq

        fastamature = dict()
        for header,seq in fasta_iter(miRBaseMature):
            if header in fastamature:
                pass
            else:
                fastamature[header]=seq


        with open(miRBasestr) as f:
            for line in f:
                splits=line.strip().split(' ')
                f1 = splits[0]

                if (name.split(':')[0] == f1):

                    pre = '.'*len(sequence)
                    mat = '.'*len(sequence)
                    seq = '.'*len(sequence)
                    MIID = name.split(':')[1]
                    prename=splits[0].strip('>')
                
                    n = len(splits)

                    refstructure = ''.join(islice(f, 7)).strip("\n").split('\n')
                
                    refstructureout=[]
                    for i in refstructure:
                        refstructureout.append(''.join(("#refstructure:",i)))

                    refstructure = refstructureout

                    if (fastahairpin.get(prename)):
                        precursor = fastahairpin.get(prename).replace('U','T')

                        if sequence.find(precursor) >= 0:
                            precursorbegin = sequence.find(precursor)
                            precursorend = sequence.find(precursor) + len(precursor)
                            pre = ''.join((pre[:precursorbegin],precursor,pre[precursorend:]))
                
                        if precursor.find(sequence) > 0:
                            pre = sequence
                    

                    refcutposition = []
                    if (n > 4):
                        matureNameList=''
                        for number in range(4,n):
                            matureName=splits[number].strip('[').split(':')[0]
                            mature = fastamature.get(matureName).replace('U','T')                            
                            maturebegin = sequence.find(mature)
                            if maturebegin < 0:
                                mat = '.'*len(sequence)
                                refcutposition.append('.'*len(sequence) + ' (Ref;none;none)')
                            else:
                                matureNameList = ' '.join((matureNameList,matureName))
                                matureend = sequence.find(mature) + len(mature)
                                refcutposition.append(''.join((seq[:maturebegin],'R',seq[maturebegin+1:])) + ' (Ref;' + str(maturebegin) + ';none)')
                                mat = ''.join((mat[:maturebegin],mature,mat[matureend:]))

                    pre = ' '.join((pre,prename))
                    mat = ''.join((mat,matureNameList))
                    refprecursor = ' '.join((precursor,prename))
                    break
        
                else:
                
                    pre = ' '.join(('.'*len(sequence),'none'))
                    mat = ' '.join(('.'*len(sequence),'none'))
                    MIID = name.split(':')[1]
                    refcutposition = '.'*len(sequence) + ' (Ref;none;none)'
                    refstructureout = []
                    refstructureout.append("#refstructure:")
                    refstructureout.append("#refstructure:")
                    refstructureout.append("#refstructure:")
                    refstructureout.append("#refstructure:")
                    refstructureout.append("#refstructure:")
                    refstructure = refstructureout
                    refprecursor = ' '.join((pre,'none'))

        return pre, mat, MIID, refstructure, refprecursor, refcutposition



    def parse_cut_positions(self, name, sequence, cut_site_variation):
        with open(cut_site_variation, 'r') as f:
            first_line = f.readline()

        out = []

        cell_types=['HEK','NIH3T3','MEF']

        for cell_type in cell_types:
            for i,line in enumerate(first_line.split("\t")):
                if line == cell_type + '_' + name.strip('>'):
                    line_number=i

                    dominant_cut=dict()
                    fp = open(cut_site_variation)
                    for i, line in enumerate(fp):
                        if i != 0:
                            dominant_cut[int(line.split("\t")[0])] = int(line.split("\t")[line_number])

                    dominant_cut = dict((k-1, v) for k, v in dominant_cut.items() if v > 0)
                    dominant_cut = sorted(dominant_cut.items(), key=itemgetter(1), reverse=True)

                    for i,item in enumerate(dominant_cut):
                        seq = '.'*len(sequence)
                        if i == 0:
                            seq = ''.join((seq[:item[0]],'D',seq[item[0]+1:]))
                            #if cell_type == cell_types[0]:
                            sequence = ''.join((sequence[:item[0]],sequence[item[0]].lower(),sequence[item[0]+1:]))
                        else:
                            seq = ''.join((seq[:item[0]],'A',seq[item[0]+1:]))
                        seq = seq + ' ' + '(' + cell_type + ';' + str(item[0]) + ';' + str(item[1]) + ')'
                        out.append(seq)

        out.append(sequence)
        return out


    def parse_dominant_cut_position(self, refmirna, cutposition):

        cutlocation_R = [m.start() for m in re.finditer('[A,T,C,G]+', refmirna)]
        cutlocation_D = []
        cutlocation_A = []

        for i in cutposition:
            m = re.search('D',i)
            if m:
                cutlocation_D.append(int(i.split(' ')[1].split(';')[1]))

        for i in cutposition:
            m = re.search('A',i)
            if m:
                cutlocation_A.append(int(i.split(' ')[1].split(';')[1]))

        #cutlocation_D = list(set(cutlocation_D))
        x = 0
        for i,var in enumerate(cutlocation_D):
            if var == x:
                cutlocation_D.pop(i)
            else:
                pass
            x = var

        #cutlocation_D = list(set(cutlocation_D))
        x = 0
        for i,var in enumerate(cutlocation_A):
            if var == x:
                cutlocation_A.pop(i)
            else:
                pass
            x = var

        cutlocation = set(cutlocation_D + cutlocation_R + cutlocation_A)
        
        annolocation = dict()
        
        for location in cutlocation:
            if location in cutlocation_D:
                annolocation[location] = 'D'
            elif location in cutlocation_R and location not in annolocation:
                annolocation[location] = 'R'
            elif location in cutlocation_A and location not in annolocation:
                annolocation[location] = 'A'
            else:
                pass
        
                    
        return annolocation
 

    #########################
    # parse_mapping_profile #
    #########################
    
    # 1. get the detailed seqeunce mapping profile in each sample. 
    def parse_mapping_profile(self, name, sequence):

        out = ()
        with open(mapping_profile, 'r') as f:
            for line in f:
                line_split = line.split('\t')
                id = line_split[0]
                cell_type = line_split[1]
                sample = line_split[2]
                count = int(line_split[3])
                begin = int(line_split[5])-1
                end = int(line_split[6])
                seq_insert = line_split[7]
                if id == name.strip('>') and 'mock' in sample:
                    seq = '.'*len(sequence)
                    seq = ''.join((seq[:begin],seq_insert,seq[end:]))
                    seq = (seq,cell_type,'mock',sample,begin,end,count)
                    out = out + (seq,)
                elif id == name.strip('>') and 'mock' not in sample:
                    seq = '.'*len(sequence)
                    seq = ''.join((seq[:begin],seq_insert,seq[end:]))
                    seq = (seq,cell_type,'transfection',sample,begin,end,count)
                    out = out + (seq,)
                else:
                    pass
            
        out = sorted(out, key=lambda key: (key[1], key[2], key[0], key[6]),reverse=True)
    
        return out

    
    # 2. identity the 5p, 3p beginning and end positions and get the statistics of mapping distribution.
    def parse_beginning_of_mapping_profile(self, mappingprofile, longesthairpinwithgap, ref, sequence, name):    

        incorrect_structure = False
        # sum seq counts based on cell type, condition and beginning sites
        refstart = [m.start() for m in re.finditer('(?<=.)[A,T,C,G]+', ref[1])]
        refend =  [m.end() for m in re.finditer('(?<=.)[A,T,C,G]+', ref[1])]
        
        sumbysample_beginning = dict()
        for line in mappingprofile:
            replicate = line[3].split('_')[1]
            concentration = line[3].split('_')[3]
            amount = line[3].split('_')[4]
            NorOrFacs = line[3].split('_')[5]
            key = (line[1],line[2],concentration,replicate,amount,NorOrFacs,line[4])
            if key in sumbysample_beginning:
                sumbysample_beginning[key] = sumbysample_beginning[key] + line[6]
            else:
                sumbysample_beginning[key] = line[6]
                
        sumbysample_end = dict()
        for line in mappingprofile:
            replicate = line[3].split('_')[1]
            concentration = line[3].split('_')[3]
            amount = line[3].split('_')[4]
            NorOrFacs = line[3].split('_')[5]
            key = (line[1],line[2],concentration,replicate,amount,NorOrFacs,line[5])
            if key in sumbysample_end:
                sumbysample_end[key] = sumbysample_end[key] + line[6]
            else:
                sumbysample_end[key] = line[6]
                
                
        sumbysample_begin_ref = dict()
        sumbysample_end_ref = dict()
        if len(refstart) == 0:
            sumbysample_begin_ref = sumbysample_beginning
            sumbysample_end_ref = sumbysample_end
        elif len(refstart) == 1:
            sumbysample_begin_ref = sumbysample_beginning
            sumbysample_end_ref = sumbysample_end
                
        elif len(refstart) == 2:
            for line in mappingprofile:
                for start in refstart:
                    if int(start) - 10 < line[4] < int(start) + 10:
                        replicate = line[3].split('_')[1]
                        concentration = line[3].split('_')[3]
                        amount = line[3].split('_')[4]
                        NorOrFacs = line[3].split('_')[5]
                        key = (line[1],line[2],concentration,replicate,amount,NorOrFacs,line[4])
                        if key in sumbysample_begin_ref:
                            sumbysample_begin_ref[key] = sumbysample_begin_ref[key] + line[6]
                        else:
                            sumbysample_begin_ref[key] = line[6] 
                            
                for end in refend:   
                    if int(end) - 10 < line[5] < int(end) + 10:
                        replicate = line[3].split('_')[1]
                        concentration = line[3].split('_')[3]
                        amount = line[3].split('_')[4]
                        NorOrFacs = line[3].split('_')[5]
                        key = (line[1],line[2],concentration,replicate,amount,NorOrFacs,line[5])
                        if key in sumbysample_end_ref:
                            sumbysample_end_ref[key] = sumbysample_end_ref[key] + line[6]
                        else:
                            sumbysample_end_ref[key] = line[6]
                    
        else:
            sumbysample_begin_ref = sumbysample_beginning
            sumbysample_end_ref = sumbysample_end
    
        beginning_end_ref = [sumbysample_begin_ref, sumbysample_end_ref]
        beginning_end = [sumbysample_beginning, sumbysample_end]
        label = ['beginning','end']
        anno_out = []
        
        for i,var_dict in enumerate(beginning_end):
            
            # get sum counts of each hairpin
            sumbysample_hairpin = dict()
            for key, var in var_dict.items():
                key_new = (key[0],key[1],key[2],key[3],key[4],key[5])
                if key_new not in sumbysample_hairpin:
                    sumbysample_hairpin[key_new] = var
                else:
                    sumbysample_hairpin[key_new] = sumbysample_hairpin[key_new] + var
        
            if i == 0:
                for key,var in sumbysample_hairpin.items():
                    name_total = '_'.join((key[0],key[1],key[2],key[3],key[4],key[5],'all','all','all','count'))
                    anno_out.append([name_total,var]) 
                    
                    
        for i,var_dict in enumerate(beginning_end_ref):
              
            # identify dominant cutting site and counts that give rise to mature sequences
            sumbysample_mature_D = dict()
            for key,var in var_dict.items():
                if int(key[6]) < 20:
                    pass
                else:
                    key_new = (key[0],key[1],key[2],key[3],key[4],key[5])
                    if key_new not in sumbysample_mature_D:
                        sumbysample_mature_D[key_new] = (key[6], var)
                    if var > sumbysample_mature_D[key_new][1]:
                        sumbysample_mature_D[key_new] = (key[6], var)
    
            # identify dominant cutting site and counts that give rise to star sequences
            sumbysample_star_D = dict()
            for key, var in var_dict.items():
                if int(key[6]) < 20:
                    pass
                else:
                    key_in_mature_D = (key[0],key[1],key[2],key[3],key[4],key[5])
                    position = sumbysample_mature_D[key_in_mature_D][0]
                    if position <= longesthairpinwithgap[-1][0] and int(key[6]) > longesthairpinwithgap[-1][0]:
                        for pair in longesthairpinwithgap:
                            if pair[0] == position:
                                position2 = pair[4]
                                break
                        if 'position2' not in locals():
                            incorrect_structure = True
                        elif position2 == '.': 
                            incorrect_structure = True
                        else:
                            if i == 0: 
                                if int(key[6]) < position2 - 10 and key_in_mature_D not in sumbysample_star_D:
                                    sumbysample_star_D[key_in_mature_D]=(key[6],var)
                                if int(key[6]) < position2 - 10 and sumbysample_star_D[key_in_mature_D][1] < var:
                                    sumbysample_star_D[key_in_mature_D]=(key[6],var)
                            if i == 1:
                                if int(key[6]) < position2 + 30 and key_in_mature_D not in sumbysample_star_D:
                                    sumbysample_star_D[key_in_mature_D]=(key[6],var)
                                if int(key[6]) < position2 + 30 and sumbysample_star_D[key_in_mature_D][1] < var:
                                    sumbysample_star_D[key_in_mature_D]=(key[6],var)    

                    elif position > longesthairpinwithgap[-1][0] and int(key[6]) <= longesthairpinwithgap[-1][0]:                       
                        for pair in longesthairpinwithgap:
                            if pair[4] == position:
                                position2 = pair[0]
                                break
                        if 'position2' not in locals():
                            incorrect_structure = True
                        elif position2 == '.':
                            incorrect_structure = True
                        else:
                            if i == 0:  
                                if int(key[6]) < position2 - 10 and key_in_mature_D not in sumbysample_star_D:
                                    sumbysample_star_D[key_in_mature_D]=(key[6],var)
                                if int(key[6]) < position2 - 10 and sumbysample_star_D[key_in_mature_D][1] < var:
                                    sumbysample_star_D[key_in_mature_D]=(key[6],var)
                            if i == 1:
                                if int(key[6]) > position2+10 and key_in_mature_D not in sumbysample_star_D:
                                    sumbysample_star_D[key_in_mature_D]=(key[6],var)
                                if int(key[6]) > position2+10 and sumbysample_star_D[key_in_mature_D][1] < var:
                                    sumbysample_star_D[key_in_mature_D]=(key[6],var)                           
                    else:
                        pass    


            # sum seq counts that located at 5nt upstream or downstream of 5p and 3p beginning site
            sumbysample_updown5_count = dict()
            sumbysample_updown5_position = dict()

            for key, var in var_dict.items():
                key_in_maturestar_D = (key[0],key[1],key[2],key[3],key[4],key[5])
                if key_in_maturestar_D in sumbysample_mature_D:
                    position_mature = sumbysample_mature_D[key_in_maturestar_D][0]
                    if (key_in_maturestar_D+('D', position_mature,)) not in sumbysample_updown5_count:
                        sumbysample_updown5_count[key_in_maturestar_D+('D', position_mature,)] = ()
                        sumbysample_updown5_position[key_in_maturestar_D+('D', position_mature,)] = ()
                    if key[6] <= position_mature + 5 and key[6] >= position_mature - 5:
                        sumbysample_updown5_count[key_in_maturestar_D+('D', position_mature,)] = sumbysample_updown5_count[key_in_maturestar_D+('D', position_mature,)] + (var,)
                        sumbysample_updown5_position[key_in_maturestar_D+('D', position_mature,)] = sumbysample_updown5_position[key_in_maturestar_D+('D', position_mature,)] + (key[6],)
            
                if key_in_maturestar_D in sumbysample_star_D:
                    position_star = sumbysample_star_D[key_in_maturestar_D][0]
                    if (key_in_maturestar_D+('M', position_star,)) not in sumbysample_updown5_count:
                        sumbysample_updown5_count[key_in_maturestar_D+('M', position_star,)] = ()
                        sumbysample_updown5_position[key_in_maturestar_D+('M', position_star,)] = ()
                    if key[6] <= position_star + 5 and key[6] >= position_star - 5:
                        sumbysample_updown5_count[key_in_maturestar_D+('M', position_star,)] = sumbysample_updown5_count[key_in_maturestar_D+('M',position_star,)] + (var,)
                        sumbysample_updown5_position[key_in_maturestar_D+('M', position_star,)] = sumbysample_updown5_position[key_in_maturestar_D+('M',position_star,)] + (key[6],)
   
            # sum seq counts that have dominant beginning and end positions. 
            if i == 0:
                sumbysample_beginning_end = dict()
                for line in mappingprofile:
                    replicate = line[3].split('_')[1]
                    concentration = line[3].split('_')[3]
                    amount = line[3].split('_')[4]
                    NorOrFacs = line[3].split('_')[5]
                    key = (line[1],line[2],concentration,replicate,amount,NorOrFacs,line[4],line[5])
                    if key in sumbysample_beginning_end:
                        sumbysample_beginning_end[key] = sumbysample_beginning_end[key] + line[6]
                    else:
                        sumbysample_beginning_end[key] = line[6]
          
                sumbysample_beginningend_D = dict()
                for key,var in sumbysample_mature_D.items():
                    value1 = key[0]
                    value2 = key[1]
                    value3 = key[2]
                    value4 = key[3]
                    value5 = key[4]
                    value6 = key[5]                    
                    value7 = var[0]
                    var_new = 0
                    for key1,var1 in sumbysample_beginning_end.items():
                        if value1 == key1[0] and value2 == key1[1] and value3 == key1[2] and value4 == key1[3] and value5 == key1[4] and value6 == key1[5] and value7 == key1[6]:
                            if var1 > var_new:
                                key_new = key1
                                var_new = var1        
                    sumbysample_beginningend_D[key_new] = var_new
                
                sumbysample_beginningend_M = dict()
                for key,var in sumbysample_star_D.items():
                    value1 = key[0]
                    value2 = key[1]
                    value3 = key[2]
                    value4 = key[3]
                    value5 = key[4]
                    value6 = key[5]                    
                    value7 = var[0]
                    var_new = 0
                    for key1,var1 in sumbysample_beginning_end.items():
                        if value1 == key1[0] and value2 == key1[1] and value3 == key1[2] and value4 == key1[3] and value5 == key1[4] and value6 == key1[5] and value7 == key1[6]:
                            if var1 > var_new:
                                key_new = key1
                                var_new = var1        
                    sumbysample_beginningend_M[key_new] = var_new

                    
            # get the all annotataion
            
            if i == 0:                
                for key,var in sumbysample_beginningend_D.items():
                    if key[6] < longesthairpinwithgap[-1][0]:
                        strand = '5p'
                    else:
                        strand = '3p'
                    
                    name_count = '_'.join((key[0],key[1],key[2],key[3],key[4],key[5],strand,'beginningend','D','count'))
                    anno_out.append([name_count,var])
                    
                for key,var in sumbysample_beginningend_M.items():
                    if key[6] < longesthairpinwithgap[-1][0]:
                        strand = '5p'
                    else:
                        strand = '3p'
                    
                    name_count = '_'.join((key[0],key[1],key[2],key[3],key[4],key[5],strand,'beginningend','M','count'))
                    anno_out.append([name_count,var])

                refstart = [m.start() for m in re.finditer('(?<=.)[A,T,C,G]+', ref[1])]
                refend =  [m.end() for m in re.finditer('(?<=.)[A,T,C,G]+', ref[1])]
                for ref in range(len(refstart)):
                    if refstart[ref] < longesthairpinwithgap[-1][0]:
                        strand = '5p'
                    else:
                        strand = '3p'
                    name_position = '_'.join(('AREF','ref','ref','ref','ref','ref',strand,'beginning','ref','position'))
                    anno_out.append([name_position,refstart[ref]])
                    name_position = '_'.join(('AREF','ref','ref','ref','ref','ref',strand,'end','ref','position'))
                    anno_out.append([name_position,refend[ref]])    

                    
            for key,var in sumbysample_mature_D.items():
                if var[0] < longesthairpinwithgap[-1][0]:
                    strand = '5p'
                else:
                    strand = '3p'
        
                name_position = '_'.join((key[0],key[1],key[2],key[3],key[4],key[5],strand,label[i],'D','position'))
                name_count = '_'.join((key[0],key[1],key[2],key[3],key[4],key[5],strand,label[i],'D','count'))
                anno_out.append([name_position,var[0]])
                anno_out.append([name_count,var[1]])

            for key,var in sumbysample_star_D.items():
                if var[0] < longesthairpinwithgap[-1][0]:
                    strand = '5p'
                else:
                    strand = '3p'
        
                name_position = '_'.join((key[0],key[1],key[2],key[3],key[4],key[5],strand,label[i],'M','position'))
                name_count = '_'.join((key[0],key[1],key[2],key[3],key[4],key[5],strand,label[i],'M','count'))
                anno_out.append([name_position,var[0]])
                anno_out.append([name_count,var[1]])  
    
            for key,var in sumbysample_updown5_count.items():
                if key[7] < longesthairpinwithgap[-1][0]:
                    strand = '5p'
                else:
                    strand = '3p'
     
                name_count = '_'.join((key[0],key[1],key[2],key[3],key[4],key[5],strand,label[i],key[6],'updown5ntcount'))
                anno_out.append([name_count,var])  
    
    
            for key,var in sumbysample_updown5_position.items():
                if key[7] < longesthairpinwithgap[-1][0]:
                    strand = '5p'
                else:
                    strand = '3p'
     
                name_count = '_'.join((key[0],key[1],key[2],key[3],key[4],key[5],strand,label[i],key[6],'updown5ntposition'))
                anno_out.append([name_count,var]) 
        
        return sorted(anno_out)


