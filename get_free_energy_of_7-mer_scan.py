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
#vienna="../hairpin-input-data/RNAfold_data_with_extended_out_hairpins_and_constrain/mmu-mir-883a:MI0005476:in:none.fold"
sample=vienna.split('/')[-1]
sample=re.sub('RNAfold_', '', sample)
sample=re.sub('.fold', '', sample)
#species=sys.argv[2]
species="HEK"
rna = RNAStructure(vienna)
miRNA = rna.name.strip(">")

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
            
#check if the 3p and 5p are encoded in the longest hairpin
index_longesthairpin_up = list(range(rna.longesthairpin[0][0],rna.longesthairpin[-1][0]))
index_longesthairpin_down = list(range(rna.longesthairpin[-1][1],rna.longesthairpin[0][1]))
    
non_longest_hairpin = False
if start_5p in index_longesthairpin_up:
    non_longest_hairpin = False
else:
    non_longest_hairpin = True 

# check if the 3p and 5p form miRNA duplex rather than internal loop or other structures.        
if non_longest_hairpin == False:
    
    test=[]

    ###
    ## scan stem in 5-mer
    ###
    
    n_mer_size = 7
    for n in range(-16,25-n_mer_size):
        
        Gibb_scan = 'NA'
        GC_content_scan = 'NA'
        GC_pairing_percent_scan = 'NA'
        size_of_bulge_scan = 'NA'
        
        up_scan_seq = rna.sequence_with_U[start_5p+n:start_5p+n+n_mer_size]
        up_scan_structure = rna.dotbracket_corrected[start_5p+n:start_5p+n+n_mer_size]
        down_scan_coordinate = []
        scanstem=[]

        for i,var in enumerate(rna.longesthairpinwithgap):
            if var[0] in list(range(start_5p+n,start_5p+n+n_mer_size)) and var[3] != '|':
                down_scan_coordinate.append(var[4])
            if var[0] in list(range(start_5p+n,start_5p+n+n_mer_size)):
                scanstem.append((i,var))

        if len(down_scan_coordinate) != 0:
            down_scan_seq = rna.sequence_with_U[min(down_scan_coordinate):max(down_scan_coordinate)+1]  
            down_scan_structure = rna.dotbracket_corrected[min(down_scan_coordinate):max(down_scan_coordinate)+1]
        else:
            down_scan_seq = ''
            down_scan_structure = ''

        duplex_scan_seq = up_scan_seq+'&'+down_scan_seq
        duplex_scan_structure = up_scan_structure+'&'+down_scan_structure
        
        

        if len([m.start() for m in re.finditer('\(', up_scan_structure)]) == len([m.start() for m in re.finditer('\)', down_scan_structure)]) and len([m.start() for m in re.finditer('\(', down_scan_structure)]) == 0 and len([m.start() for m in re.finditer('\)', up_scan_structure)]) == 0:

            outfile_scan="./output/" + 'tmp.tmp'
            with open(outfile_scan, 'w') as out:
                out.write('>'+ miRNA + '\n' + duplex_scan_seq + '\n' + duplex_scan_structure + '\n')

            p1 = subprocess.Popen(["cat", outfile_scan], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(['RNAeval'], stdin=p1.stdout, stdout=subprocess.PIPE)
            p1.stdout.close()

            output,err = p2.communicate()
            Gibb_scan = ''.join(output.decode("utf-8").split('\n')[2].split()[1:]).strip('\)').strip('\(')

        # GC content/13 nt, GC pairing percentage among all pairs, n_bulges, total_size_of_bulges
        GC_content_scan = len([m.start() for m in re.finditer('[G,C]', duplex_scan_seq)])/(len(duplex_scan_seq)-1)*100

        count_pair_scan = 0
        count_GCpair_scan = 0
        for line in scanstem:
            if line[1][2] == '-':
                count_pair_scan = count_pair_scan + 1
                if ''.join(line[1][1:4]) == "C-G" or ''.join(line[1][1:4]) == "G-C":
                    count_GCpair_scan = count_GCpair_scan + 1
        if count_pair_scan != 0:
            GC_pairing_percent_scan = count_GCpair_scan/count_pair_scan*100

        size_of_bulge_scan = len(scanstem) - count_pair_scan

        test.append((n,duplex_scan_seq,duplex_scan_structure,Gibb_scan,GC_content_scan,GC_pairing_percent_scan,size_of_bulge_scan))
    
    outfile_test="./output/" + miRNA + '_' + str(n_mer_size) + '-mer.scan.gibbs'
    with open(outfile_test, 'w') as out:
        for row in test:
            out.write("\t".join([str(i) for i in row]) + '\n')

