# AmpliSeqDesigner

This tool designs primers for multiplex PCR of environmental or similarly contaminated samples to detect and genotype target organism (ex. *Salmonella* Typhi or *Salmonella* Typhimurium). While tool can design primer of any length, their are inteded for Oxford Nanopore devices.

<br/>
AmpliSeqDesigner requires three core inputs:

A) User generated VCF files where each file is a sample of target organism mapped to user chosen reference sequence

B) User defined genotypes for each VCF file (ex. GenoTyphi for *Salmonella* Typhi or MLST for invasive non-typhoidal *Salmonella*). 

C) User chosen list of assemblies of organism that primers **should not target** - this can be easily obtained via "datasets" tool from NCBI and for *Salmonella* Typhi can be set of non-Typhi *Salmonella* and other Enterobacterales. 

<br/>
AmpliSeqDesigner works in XXX stage:

First, it uses VCFs from (A) and genotypes data from (B) to identify which SNPs are uniquely associated with genotypes

Second, it compares target organism reference sequence around SNPs identfied in previous step to the assemblies in (C) to identify homologues in non-target organism.

Third, AmpliSeqDesigner looks for sites in homologues identified in previous step where target and non-target organisms differ by one or more nucletodies. 

Fourth, it generates primers pairs where at least one 3' primer end is placed at sites identified at previous step. Mismatches at primer 3' end have high impact on primer efficincy and such placement improves detection of target organism in contaminated samples.

## Installation

At present, the tool can be downloaded via GitHub, but it will soon be available as a bioconda package. 
Dependencies:
--mafft
--blast
--minimap2

## Running
The tool requires two input arguments and can be run as follows:
```
python run.py -m Amplicon -c config.json
or
python run.py -m SNP -c config.json
```
"SNP" option only identified and reports the number of SNPs that differentiate genotypes, whereas Amplicon also attempts to design amplicons for these SNPs.

Due to large number of options and to improve reproducability most inputs are specified via a JSON file (config.json above). 
