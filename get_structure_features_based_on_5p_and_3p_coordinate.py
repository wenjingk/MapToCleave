#!/usr/bin/python3

import sys
from itertools import zip_longest
import re
import os.path
from operator import itemgetter
from itertools import groupby
import subprocess
from scripts.RNAStructure import RNAStructure


vienna=sys.argv[1]
sample=vienna.split('/')[-1]
sample=re.sub('RNAfold_', '', sample)
sample=re.sub('.fold', '', sample)
species="HEK"
rna = RNAStructure(vienna)

# identify drosha and dicer cutting site
out = {}

for i in rna.mappingprofile_statistics:
    #print(i)
    if re.search('transfection_.*_[3,5]p_beginning_[D,M]_position',i[0]):
        if species in i[0]:
            new_pattern = re.sub("_position","_count",i[0])
            for k in rna.mappingprofile_statistics:
                if re.search(new_pattern,k[0]) and i[1] not in out:
                    out[i[1]] = k[1]
                elif re.search(new_pattern,k[0]) and i[1] in out:
                    out[i[1]] = out[i[1]] + k[1]
                else:
                    pass

for line in rna.cutposition:
    var1=line.split(' ')[1].strip('(').split(';')[0]
    var2=line.split(' ')[1].strip('(').split(';')[1]
    if var2 == 'none':
        pass
    else:
        if var1 == 'Ref':
            if int(var2) not in out:
                out[int(var2)] = 1

               
sorted_out = sorted(out.items(), key=itemgetter(1), reverse=True)


for i,var in enumerate(sorted_out):
    if i == 0:
        x = var[0]
    else:
        if x - 10 < var[0] and x + 10 > var[0]:
            sorted_out[i] = (sorted_out[i][0],0)
        elif x > rna.longesthairpinwithgap[-1][0] and 'y' not in locals() and x - 25 > var[0]:
            y = var[0]
        elif x < rna.longesthairpinwithgap[-1][0] and 'y' not in locals() and var[0] > rna.longesthairpinwithgap[-1][4]:
            y = var[0]
        elif 'y' in locals():
            if y - 10 < var[0] and y + 10 > var[0]:
                sorted_out[i] = (sorted_out[i][0],0)
        else:
            sorted_out[i] = (sorted_out[i][0],0)  
            
if 'y' in locals():
    del y

out = {}
for var in sorted_out:
    if var[1] != 0:
        out[var[0]] = var[1]
        
#######################
# Interested features #
#######################
bulge_unopened_5p={}
L1=list(range(-25,31))
L2 = ['NA'] * len(range(-25,31))
bulge_unopened_5p=dict(zip(L1,L2))

bulge_unopened_3p={}
L1=list(range(-25,31))
L2 = ['NA'] * len(range(-25,31))
bulge_unopened_3p=dict(zip(L1,L2))

bulge_unopenedExtended={}
L1=list(range(-25,31))
L2 = ['NA'] * len(range(-25,31))
bulge_unopenedExtended=dict(zip(L1,L2))

bulge_unopenedCollapsed={}
L1=list(range(-25,31))
L2 = ['NA'] * len(range(-25,31))
bulge_unopenedCollapsed=dict(zip(L1,L2))


stem_features={}
stem_features['stem_length_with_unopened_5pAsRef'] = 'NA'
stem_features['basal_stem_length_with_unopened_5pAsRef'] = 'NA'
stem_features['miRNA_duplex_length_with_unopened_5pAsRef'] = 'NA'

basic_info={}
basic_info['pre-miRNA'] = rna.sequence.replace('T', 'U')
basic_info['seq_5p'] = 'NA'
basic_info['seq_3p'] = 'NA'
basic_info['structure'] = 'NA'
basic_info['start_5p'] = 'NA'
basic_info['end_5p'] = 'NA'
basic_info['start_3p'] = 'NA'
basic_info['end_3p'] = 'NA'

start_5p = None
end_5p = None
start_3p = None
end_3p = None

motif = {}
motif['UG_atBasalJunction'] = 'NA'
motif['UGUorGUG_atApicalJunction'] = 'NA'
motif['GHG_at_3pBasalStem'] = 'NA'
motif['CNNC_at3pFlanking'] = 'NA'

GHG = {}
GHG['GHG_fang'] = 'NA'
GHG['GHG_kwon'] = 'NA'


###################
# Data processing #
###################

#############################################
# get basic information about 3p and 5p arm #
#############################################           
        
arm = "NA"
if len(out) != 0:

    if max(out.values()) > 1:
        cut_in_seq = max(out, key=out.get)
    else:
        cut_in_seq = min(out.keys())

    if rna.longesthairpin[-1][1] > 80 or rna.longesthairpin[-1][1] < 30:
        basic_info['structure'] = 'weird'
    else:
        basic_info['structure'] = 'fold'
        if cut_in_seq < rna.longesthairpin[-1][1]:
            # get start and end position of 5p arm
            p5_cut_in_seq = cut_in_seq
            start_5p = p5_cut_in_seq
            end_5p = dict()
            for i in rna.mappingprofile:
                if i[1] == species and i[2] == 'transfection' and i[4] == start_5p:
                    if i[5] in end_5p:
                        end_5p[i[5]] = end_5p[i[5]]+i[6]
                    else:
                        end_5p[i[5]] = i[6]
                        
            if not end_5p:
                end_5p = start_5p + rna.refmirna[start_5p:].find('.')
            else:
                end_5p = sorted(end_5p.items(), key=itemgetter(1), reverse=True)[0][0]
            
            
            # get start and end position of 3p arm
            p3_cut_in_seq = None
            out = {key: value for key, value in out.items() if key != p5_cut_in_seq}
            out_sorted_keys = sorted(out, key=out.get, reverse=True)
            for r in out_sorted_keys:
                if r >= rna.longesthairpin[-1][1]:
                    p3_cut_in_seq = r
                    break
                    
                    
            arm='5p'
                    
            if p3_cut_in_seq == None:
                for i in rna.longesthairpinwithgap:
                    if i[0]==end_5p:
                        start_3p = i[4] + 3
                    if i[0]==start_5p:
                        end_3p = i[4] + 3
            else:
                start_3p = p3_cut_in_seq
                end_3p = dict()
                for i in rna.mappingprofile:
                    if i[1] == species and i[2] == 'transfection' and i[4] == start_3p:
                        if i[5] in end_3p:
                            end_3p[i[5]] = end_3p[i[5]]+i[6]
                        else:
                            end_3p[i[5]] = i[6]
                if not end_3p and rna.refmirna[start_3p:].find('.') >= 0:
                    end_3p = start_3p + rna.refmirna[start_3p:].find('.')
                elif not end_3p and rna.refmirna[start_3p:].find('.') == -1:
                    end_3p = len(rna.sequence)
                else:
                    end_3p = sorted(end_3p.items(), key=itemgetter(1), reverse=True)[0][0]                        
                        
        else:
            p3_cut_in_seq = cut_in_seq
            start_3p = p3_cut_in_seq
            end_3p = dict()
            for i in rna.mappingprofile:
                if i[1] == species and i[2] == 'transfection' and i[4] == start_3p:
                    if i[5] in end_3p:
                        end_3p[i[5]] = end_3p[i[5]]+i[6]
                    else:
                        end_3p[i[5]] = i[6]
            if not end_3p and rna.refmirna[start_3p:].find('.') >= 0:
                end_3p = start_3p + rna.refmirna[start_3p:].find('.')
            elif not end_3p and rna.refmirna[start_3p:].find('.') == -1:
                end_3p = len(rna.sequence)
            else:
                end_3p = sorted(end_3p.items(), key=itemgetter(1), reverse=True)[0][0]
            
            
            p5_cut_in_seq = None
            out = {key: value for key, value in out.items() if key != p3_cut_in_seq}
            out_sorted_keys = sorted(out, key=out.get, reverse=True)
            for r in out_sorted_keys:
                if r < rna.longesthairpin[-1][1]:
                    p5_cut_in_seq = r
                    break

            arm='3p'
            
            
            if p5_cut_in_seq == None:
                for i in rna.longesthairpinwithgap:
                    if i[4]==end_3p:
                        start_5p = i[0] + 3
                    if i[4]==start_3p:
                        end_5p = i[0] + 3   
            else:
                start_5p = p5_cut_in_seq
                end_5p = dict()
                for i in rna.mappingprofile:
                    if i[1] == species and i[2] == 'transfection' and i[4] == start_5p:
                        if i[5] in end_5p:
                            end_5p[i[5]] = end_5p[i[5]]+i[6]
                        else:
                            end_5p[i[5]] = i[6]
                        
                if not end_5p:
                    end_5p = start_5p + rna.refmirna[start_5p:].find('.')
                else:
                    end_5p = sorted(end_5p.items(), key=itemgetter(1), reverse=True)[0][0]                
            
            
        if start_5p != None and end_5p != None:
            seq_5p = rna.sequence_with_U[start_5p:end_5p]  
            basic_info['seq_5p'] = seq_5p
            basic_info['start_5p'] = start_5p
            basic_info['end_5p'] = end_5p
        
        if start_3p != None and end_3p != None:
            seq_3p = rna.sequence_with_U[start_3p:end_3p]   
            basic_info['seq_3p'] = seq_3p
            basic_info['start_3p'] = start_3p
            basic_info['end_3p'] = end_3p
            
                        
        ##################################################
        # get bulge features based on unopened structure #
        ##################################################

        hairpin_with_gap = []

        for var in rna.longesthairpinwithgap:
            if var[0] != '.':
                hairpin_with_gap.append((var[0] - start_5p,var))

        ######        
        #### 1. get the bulge profile in unopened_5p hairpin structure
        ######
        unopened_5p = []
        
        for i in hairpin_with_gap:
            #if i[1][2] != '-':
            #    unopened_5p.append((i[0], i[1][2]))
            
            if i[1][2] == '*':
                unopened_5p.append((i[0], i[1][2]))
            if i[1][2] == '<' and  i[1][1] != '|':
                unopened_5p.append((i[0], i[1][2]))
            if i[1][2] == '>':
                unopened_5p.append((i[0], i[1][2]))
            
        unopened_5p = sorted(set(unopened_5p), key=lambda x: x[0])
                
        for key,var in bulge_unopened_5p.items():
            if key >= unopened_5p[0][0] and key <= unopened_5p[-1][0]:
                bulge_unopened_5p[key] = 0
             
        # profile bulge located at -25 ~ 30 with begining nucleotide of 5p arm as coordinate reference (0).
        for item in unopened_5p:
            if item[0] <= 30 and item[0] >= -25:
                bulge_unopened_5p[item[0]] = item[1]   

        # get the pair nucleotides        
        for key,var in bulge_unopened_5p.items():
            if var == 0:
                for line in hairpin_with_gap:
                    if key == line[0] and line[1][1] != '|':
                    #if key == line[0] and bulge_unopened_5p[key] == 0:
                        bulge_unopened_5p[key] = line[1][1] + line[1][2] + line[1][3]
        
        # write out
        bulge_unopened_5p_outFile="./output/"+ sample +"."+species+".unopened_5p_ref"
        with open(bulge_unopened_5p_outFile, 'w') as output:
           for key,var in bulge_unopened_5p.items():
               output.write(str(key) + '\t' + str(var) + '\n')        
                
        ######
        #### 2. get the bulge profile in unopened_3p hairpin structure
        ######
        hairpin_with_gap_4_3p = []
        
        for i,x in enumerate(rna.longesthairpinwithgap):
            if x[0] == start_5p:
                start_5p_at_3parm = x[4]
                break
            
        for var in rna.longesthairpinwithgap:
            if var[4] != '.':
                hairpin_with_gap_4_3p.append((start_5p_at_3parm - var[4],var))
        
        
        unopened_3p = []
        
        for i in hairpin_with_gap_4_3p:
            #if i[1][2] != '-':
            #    unopened_3p.append((i[0], i[1][2]))            

            if i[1][2] == '*':
                unopened_3p.append((i[0], i[1][2]))
            if i[1][2] == '<':
                unopened_3p.append((i[0], i[1][2]))
            if i[1][2] == '>' and i[1][3] != '|':
                unopened_3p.append((i[0], i[1][2]))
            
        unopened_3p = sorted(set(unopened_3p), key=lambda x: x[0])
                
        for key,var in bulge_unopened_3p.items():
            if key >= unopened_3p[0][0] and key <= unopened_3p[-1][0]:
                bulge_unopened_3p[key] = 0
             
        # profile bulge located at -25 ~ 30 with begining nucleotide of 5p arm as coordinate reference (0).
        for item in unopened_3p:
            if item[0] <= 30 and item[0] >= -25:
                bulge_unopened_3p[item[0]] = item[1]   

        # get the pair nucleotides        
        for key,var in bulge_unopened_3p.items():
            if var == 0:
                for line in hairpin_with_gap_4_3p:
                    if key == line[0] and line[1][3] != '|':
                        bulge_unopened_3p[key] = line[1][1] + line[1][2] + line[1][3]       
        
        #write out
        bulge_unopened_3p_outFile="./output/"+ sample +"."+species+".unopened_3p_ref"
        with open(bulge_unopened_3p_outFile, 'w') as output:
           for key,var in bulge_unopened_3p.items():
               output.write(str(key) + '\t' + str(var) + '\n')                     
                
                              
        #####################
        # get stem features #
        #####################
            

        # Measure the length of stem
        count = 0
        state = False
        stem_begin = None
        stem_end = None
        for i,var in enumerate(reversed(hairpin_with_gap)):
            if var[0] > -22 and var[0] < -10:
                if var[1][2] != '-':
                    if state == True:
                        count = count + 1
                    else:
                        state = True
                else:
                    state = False

            if count == 1:
                stem_begin = list(reversed(hairpin_with_gap))[i-2][0]
                break

        count = 0
        state = False
        for i,var in enumerate(hairpin_with_gap):
            if var[0] > 18 and var[0] < 33:
                if var[1][2] != '-':
                    if state == True:
                        count = count + 1
                    else:
                        state = True
                else:
                    state = False

            if count == 1:
                stem_end = hairpin_with_gap[i-2][0]
                break

        if stem_begin != None and stem_end != None:
            stem_len = abs(stem_begin) + stem_end + 1
            stem_features['stem_length_with_unopened_5pAsRef'] = stem_len
            stem_features['basal_stem_length_with_unopened_5pAsRef'] = abs(stem_begin)
            stem_features['miRNA_duplex_length_with_unopened_5pAsRef'] = stem_end + 1

        #write out
        statistics_outFile="./output/"+ sample +"."+species+".statistics"
        with open(statistics_outFile, 'w') as output:
           for key,var in basic_info.items():
               output.write(str(key) + '\t' + str(var) + '\n')   
           for key,var in stem_features.items():
               output.write(str(key) + '\t' + str(var) + '\n')   

        
        ##################################
        # get bulge detailed information #
        ##################################
        
        bulge_features_unopened_5p = []
        coordinate = []
        coordinate_filtered = []
        exist_bulge = True
        count = 1
        nt_set = set(['A', 'T', 'C', 'G', 'N'])
        for i,var in enumerate(hairpin_with_gap):
            if var[1][0] == '.' or var[1][4] == '.':
                pass
            else:
                if var[1][2] != '-' and var[0] != hairpin_with_gap[-1][0]:
                    if exist_bulge == False:
                        coordinate.append(hairpin_with_gap[i-1])
                        count = count + 1
                    exist_bulge = True
                    coordinate.append(var)
                else:
                    coordinate.append(var)
                    if len(coordinate) >= 3:
                        for j in coordinate:
                            if j[1][3] == '|' and j[1][1] == '|':
                                coordinate_filtered.append(j)
                            elif j[1][1] in nt_set:
                                coordinate_filtered.append(j)
                            else:
                                pass
                        if len(coordinate_filtered) >= 3:
                            name = 'bulge' + str(count)
                            size = len(coordinate)-2
                            if coordinate[1][1][2] == '>':
                                type = '>'
                            elif coordinate[1][1][2] == '<':
                                type = '<'
                            else:
                                type = '*'        
                
                            if var[0] == hairpin_with_gap[-1][0]:
                                position_end = coordinate_filtered[-1][0]
                            else:
                                position_end = coordinate_filtered[-2][0]
                            position_start = coordinate_filtered[1][0]                
                            p5_seq = []
                            p3_seq = []
                            for i in coordinate:
                                p5_seq.append(i[1][1])
                                p3_seq.append(i[1][3])
                            p5_seq = ''.join((p5_seq))
                            p5_seq = re.sub('[|]','-',p5_seq)
                            p3_seq = ''.join((p3_seq))
                            p3_seq = re.sub('[|]','-',p3_seq)
                            bulge_features_unopened_5p.append(rna.name.strip('>') + '\t' + name + '\t' + str(size) + '\t' + 
                                            str(position_start) + '\t' + str(position_end) + '\t' + type + '\t' + 
                                            p5_seq + '\t' + p3_seq + '\t' + 
                                            p5_seq[0] + '-' + p3_seq[0] + '\t' +
                                            p5_seq[-1] + '-' + p3_seq[-1] + '\n')
                    exist_bulge = False
                    coordinate = []  
                    coordinate_filtered = []
                    
        feature_outFile="./output/"+ sample +"."+species+".bulgefeatures_unopened_5p_ref"
        with open(feature_outFile, 'w') as output:
           for i in bulge_features_unopened_5p:
               output.write(i)                    

        
        bulge_features_unopened_3p = []
        coordinate = []
        coordinate_filtered = []
        exist_bulge = True
        count = 1
        
        for i,var in enumerate(hairpin_with_gap_4_3p):
            if var[1][0] == '.' or var[1][4] == '.':
                pass
            else:
                if var[1][2] != '-' and var[0] != hairpin_with_gap_4_3p[-1][0]:
                    if exist_bulge == False:
                        coordinate.append(hairpin_with_gap_4_3p[i-1])
                        count = count + 1
                    exist_bulge = True
                    coordinate.append(var)
                else:
                    coordinate.append(var)
                    if len(coordinate) >= 3:
                        for j in coordinate:
                            if j[1][3] == '|' and j[1][1] == '|':
                                coordinate_filtered.append(j)
                            elif j[1][3] in nt_set:
                                coordinate_filtered.append(j)
                            else:
                                pass
                                
                        if len(coordinate_filtered) >= 3:
                            name = 'bulge' + str(count)
                            size = len(coordinate)-2
                            if coordinate[1][1][2] == '>':
                                type = '>'
                            elif coordinate[1][1][2] == '<':
                                type = '<'
                            else:
                                type = '*'        
                
                            if var[0] == hairpin_with_gap_4_3p[-1][0]:
                                position_end = coordinate_filtered[-1][0]
                            else:
                                position_end = coordinate_filtered[-2][0]
                            position_start = coordinate_filtered[1][0]                
                            p5_seq = []
                            p3_seq = []
                            for i in coordinate:
                                p5_seq.append(i[1][1])
                                p3_seq.append(i[1][3])
                            p5_seq = ''.join((p5_seq))
                            p5_seq = re.sub('[|]','-',p5_seq)
                            p3_seq = ''.join((p3_seq))
                            p3_seq = re.sub('[|]','-',p3_seq)
                            bulge_features_unopened_3p.append(rna.name.strip('>') + '\t' + name + '\t' + str(size) + '\t' + 
                                            str(position_start) + '\t' + str(position_end) + '\t' + type + '\t' + 
                                            p5_seq + '\t' + p3_seq + '\t' + 
                                            p5_seq[0] + '-' + p3_seq[0] + '\t' +
                                            p5_seq[-1] + '-' + p3_seq[-1] + '\n')
                    exist_bulge = False
                    coordinate = []  
                    coordinate_filtered = []
                    
        feature_outFile="./output/"+ sample +"."+species+".bulgefeatures_unopened_3p_ref"
        with open(feature_outFile, 'w') as output:
           for i in bulge_features_unopened_3p:
               output.write(i)                    
                    
        ##############
        # get motifs #
        ##############

        #1. get the CNNC motif
        if len(rna.sequence) >= end_3p+24:
            if len(list(re.finditer('C..C', rna.sequence[end_3p+15:end_3p+24]))) == 0:
                motif['CNNC_at3pFlanking'] = 0
            else:
                motif['CNNC_at3pFlanking'] = 1


        #2. get the UG at Basal
        if start_5p-15 >= 0:
            if len(list(re.finditer('TG', rna.sequence[start_5p-15:start_5p-11]))) == 0:
                motif['UG_atBasalJunction'] = 0
            else:
                motif['UG_atBasalJunction'] = 1

        #3. get the UGU or GUG at Apical
        rx = re.compile('TGT|GTG')
        if len(list(rx.finditer(rna.sequence[start_5p+19:start_5p+27]))) == 0:
            motif['UGUorGUG_atApicalJunction'] = 0
        else:
            motif['UGUorGUG_atApicalJunction'] = 1


        #4. get GHG motif at 3p basal stem

        hairpin_with_gap = []
        for var in rna.longesthairpinwithgap:
            if var[0] != '.':
                hairpin_with_gap.append((var[0] - start_5p,var))

        ## definition based on Fang et al. 2015
        position_7 = ['C-G','T-G']
        position_6 = ['C*T','T*C','G*A','G-C','A-T','A*C','T-A']
        position_5 = ['A-T','T-A','G-C','C-G']

        ## definition based on Kwon et al. 2018
        GHGscore = dict()
        with open("./data/Kwon_GHGscore_from_suppl_table1.txt",'r') as kwon:
            for line in kwon:
                line=line.strip()
                line=line.split("\t")
                GHGscore[line[0]] = round(float(line[1]), 2)



        GHG_region=[]
        for i,var in enumerate(hairpin_with_gap):
            if var[0] >= -7 and var[0] <= -5 and var[1][1] != '|' and var[1][3] != '|':
                GHG_region.append(var)


        len(GHG_region)
        base_7 = None
        base_6 = None
        base_5 = None
        kwon_seq = None

        if len(GHG_region)==3:
            base_7 = ''.join([GHG_region[0][1][1],GHG_region[0][1][2],GHG_region[0][1][3]])
            base_6 = ''.join([GHG_region[1][1][1],GHG_region[1][1][2],GHG_region[1][1][3]])
            base_5 = ''.join([GHG_region[2][1][1],GHG_region[2][1][2],GHG_region[2][1][3]])
            kwon_seq = ''.join([GHG_region[0][1][1],GHG_region[1][1][1],GHG_region[2][1][1],
                                GHG_region[0][1][3],GHG_region[1][1][3],GHG_region[2][1][3]])

            kwon_seq = kwon_seq.replace('T','U')
            GHG['GHG_kwon'] = GHGscore[kwon_seq]

        else:
            GHG['GHG_kwon'] = 'NA'


        if base_7 != None and base_6 != None and base_5 != None:
            if base_7 in position_7 and base_6 in position_6 and base_5 in position_5:
                motif['GHG_at_3pBasalStem'] = 1
                GHG['GHG_fang'] = 1
            else:
                motif['GHG_at_3pBasalStem'] = 0
                GHG['GHG_fang'] = 0

        motif_outFile="./output/"+ sample +"."+species+".motif"
        with open(motif_outFile, 'w') as output:
          for key,var in motif.items():
              output.write(str(key) + '\t' + str(var) + '\n')


        GHG_outFile="./output/"+ sample +"."+species+".GHG"
        with open(GHG_outFile, 'w') as output:
           for key,var in GHG.items():
               output.write(str(key) + '\t' + str(var) + '\n')


                
