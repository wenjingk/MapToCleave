#!/usr/bin/python3

import sys
from scripts.RNAStructure import RNAStructure

# input
vienna=sys.argv[1]
rna = RNAStructure(vienna)


# ### MapToCleave hairpin information

# In[2]:

print('# MapToCleave hairpin identity')
print(rna.name)


# In[3]:

print('# MapToCleave hairpin sequence (118nt)')
print(rna.sequence_with_U)


# In[4]:

print('# MapToCleave hairpin secondary structure in dotbracket format')
print(rna.dotbracket)


# In[5]:

print('# MapToCleave hairpin secondary structure in text format')
for i in rna.structure:
    print(i)


# In[6]:

print('# reference hairpin, mature and star information based on miRBase v21')
print(rna.refhairpin)
print(rna.refmirna)


# In[7]:

print('# start site of 5p and/or 3p arm identified based on the reads pooled from HEK-293T, NIH-3T3 and MEF transfection samples')
print('# the R means reference start site; D means dominant start site; A means alternative start site. The D and A sites are generated based on the mapping profile')
print('# the text in bracket means sample;position;read abundance')
print(rna.refhairpin)
print(rna.refmirna)

for i in rna.cutposition:
    print(i)
    


# In[8]:

for i in rna.mappingprofile_print:
    print(i)


# In[9]:

print('# information obtained from the mapping profile')
print('# D means mature arm and M means start arm identified based on the mapping profile')
for i in rna.mappingprofile_statistics:
    print(i)

