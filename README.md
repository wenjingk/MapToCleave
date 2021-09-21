# MapToCleave code

<!-- MarkdownTOC -->

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Input](#input)
- [Pipeline](#pipeline)
- [Output](#output)

<!-- /MarkdownTOC -->

<br><br>
## Introduction

MapToCleave is a high throughput screening method to profile hairpin features that affect miRNA processing in living cells (doi: https://doi.org/10.1101/2021.08.03.454879). The scripts used for the study are listed as following: 

1. The scripts used to profile structural features and sequence motifs of a given hairpin can be found in the [Pipeline](#pipeline) section. 
2. The R code used to generate main Figures can be found in jupyter notebook `R_code_related_to_the_main_figures.ipynb`. 


<br><br>
## Requirements
- `RNAfold` from The ViennaRNA Package ([here](https://www.tbi.univie.ac.at/RNA/))
- `partition` and `ProbabilityPlot` from ([RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html))
- Tested in python version 3.6 and 3.8. 


<br><br>
## Input

An example of input file at `data/hsa-mir-371a:MI0000779:in:none.fold` ([here](https://github.com/wenjingk/MapToCleave/blob/master/data/hsa-mir-371a:MI0000779:in:none.fold))

```
>hsa-mir-371a:MI0000779:in:none
CCGCCUUGCCGCAUCCCCUCAGCCUGUGGCACUCAAACUGUGGGGGCACUUUCUGCUCUCUGGUGAAAGUGCCGCCAUCUUUUGAGUGUUACCGCUUGAGAAGACUCAACCUGCGGAG
.................((((((..((((((((((((..((((.(((((((((.(((....)))))))))))).))))..)))))))))))).)).)))).................. (-46.70)
```

Each hairpin can be divided into four parts, comprising the flanking region, basal stem, miRNA duplex and apical loop, depending on the single and double stranded structure changes and Drosha cleavage site (illustrated by the example hairpin below). Since Drosha cleavage is the critical entry point for canonical miRNA biogenesis, we defined the coordinates along the hairpin stem loop relative to the Drosha cleavage site at the 5’ strand. For example, in the hairpin stem loop of hsa-mir-371a:

![Image of hairpin](https://github.com/wenjingk/MapToCleave/blob/master/images/hairpin.png)

Here, the nucleotides away from the Drosha cleavage site to the apical loop are counted from 0 to 28 and 0 to 27 for the 5’ and 3’ strand. The nucleotides away from the Drosha cleavage site to the flanking sequences are counted from -1 to -30 and -1 to -31 for the 5’ and 3’ strand respectively. The 5p and 3p arm sequences at the 5’ and 3’ strand are colored in red and blue respectively. In this coordinate system, all the features were counted according to the relative distance from the Drosha cleavage site.

For each hairpin, we parsed the secondary structure to profile base pairing and mismatch information of each nucleotide at the 5’ strand and at the 3’ strand of the hairpin. If the nucleotide is base paired, we record the type of base pair, including A-U, U-A, C-G, G-C, G-U or U-G; if the nucleotide is unpaired and present as part of a bulge, we record the features of the bulge, including its size, symmetry and nucleotide content. We also checked if the known sequence motifs that influence miRNA processing are present at the expected position. For example, the GU motif at the basal junction is expected to be located at position -15 to -12 of the 5’strand; the UGU or GUG motif at the apical junction is expected to be located at position + 19 to +27 of the 5’ strand; the CNNC motif is expected to be located at position -24 to -17 of the 3’ strand.


<br><br>
## Pipeline

1. Show secondary structure and mapping profile of the selected hairpin sequence.

```bash
id="hsa-mir-371a:MI0000779:in:none"
python3.6 RNAStructure_processing_example.py data/${id}.fold
```

2. Profile structural features and known sequence motifs that affect miRNA processing based on the secondary RNA structure generated in step 1. 

```bash
id="hsa-mir-371a:MI0000779:in:none"
python3.6 get_structure_features_based_on_5p_and_3p_coordinate.py data/${id}.fold
```

3. calculate local free energy in 7-mers window (The n-mer can be changed in the python script).

```bash
id="hsa-mir-371a:MI0000779:in:none"
python3.6 get_free_energy_of_7-mer_scan.py data/${id}.fold
```

4. Calculate positional shannon entropy using the same approach in `Rice et al. 2020` ([here](https://www.cell.com/molecular-cell/fulltext/S1097-2765(20)30735-8)).

```bash
id="hsa-mir-371a:MI0000779:in:none"
cat data/${id}.fold | head -2 > data/${id}_tmp.fa
partition data/${id}_tmp.fa output/${id}.pfs
ProbabilityPlot output/${id}.pfs output/${id}_base_pair_probability_from_RNAstructure.txt -t
# the positional shannon entropy is recoodinate based on the Drosha cleavage site (at position 0) and Dicer cleavage site (at position 100)
Rscript scripts/get_shannon_entropy_rice_method.R ${id}
```

<br><br>
## Output

1. outputs from pipeline step 1. 
see example in `RNAStructure_processing_example.ipynb`.

2. outputs from pipeline step 2. using id="hsa-mir-371a:MI0000779:in:none" as example. 

- output/${id}.HEK.statistics
The file provides basic information of the hairpin.

- output/${id}.HEK.unopened_5p_ref and output/${id}.HEK.unopened_5p_ref
The files provide positional information of base pairing indicated by “-” and symmetric mismatch bulge indicated by “*”, “>” asymmetric bulge with more nucleotides at the 5p strand or “<”
asymmetric bulge with more nucleotides at the 3p strand.

- output/${id}.HEK.bulgefeatures_unopened_5p_ref and output/${id}.HEK.bulgefeatures_unopened_3p_ref
The files provide detailed information on bulges with columns for identity, bulge name, maximum size of the bulge, start position of the bulge, end position of the bulge, type of bulge, nucleotides at the 5p strand, nucleotides at the 3p strand, the adjacent base pair at the 5’ side of the bulge and the adjacent base pair at the 3’ side of the bulge. The 5p and 3p in the file names mean that the features are counted respectively at the 5p and 3p strand of the hairpin. 

- output/${id}.HEK.motif and output/${id}.HEK.GHG
The files show the presence of the known sequence motifs and GHG features that are identified using different definitions from Fang et al. 2015, Kwon et al. 2019 and the MapToCleave study. 

**Note**: All the features were counted based on the 5’ side of Drosha cleavage site.

3. outputs from pipeline step 3. 

- output/${id}_7-mer.scan.gibbs
The file provide detailed information on the 7-mer windows slided across the given hairpin structure. The columns represent 5' strand position, nucleotides, structures and minimial free eenergy, GC content, percentage of GC base pair, and number of mismatches of the extracted 7-mer windows. 

4. outputs from pipeline step 4. 

- output/${i}_base_pair_probability_from_RNAstructure.txt
- output/${i}_positional_entropy_drosha_site_as_rice_paper.txt

