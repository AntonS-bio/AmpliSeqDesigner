from inputs_validation import ValidateFiles
from fasta_utils import FastaUtilities
from generate_primers import PrimerGenerator
import name_converters
import time 

#### START Input validation tests ####
# file_validator=ValidateFiles()
# #fasta_file_name="/home/ubuntu/HandyAmpliconTool/test_data/test_data.fasta"
# fasta_file_name="/home/ubuntu/HandyAmpliconTool/test_data/GCF_000195995.1.fna"
# bed_file_name="/home/ubuntu/HandyAmpliconTool/test_data/test_data.bed"
# vcf_file="/home/ubuntu/HandyAmpliconTool/test_data/vcfs/8490_5#12.vcf"
# file_validator.validate_bed(bed_file_name)
# file_validator.validate_fasta(fasta_file_name)
# file_validator.check_contigs_in_fasta(bed_file_name, fasta_file_name)
# file_validator.validate_vcf(vcf_file)
#### END Input validation tests ####

#### START Identification of genotype defining SNPS #### 
from os import listdir
import metadata_utils as metadata_utils
#from metadata_utils import MetadataUtilities
from load_vcfs import VCFutilities
from hierarchy_utils import HierarchyUtilities
import pandas as pd
import numpy as np
from collections import Counter
from typing import Dict, List
import sys


vcf_dir: str="/home/ubuntu/HandyAmpliconTool/test_data/vcfs/"
vcf_files: List[str]=[f'{vcf_dir}{f}' for f in listdir(vcf_dir)]
repeat_regions_file: str="/home/ubuntu/HandyAmpliconTool/test_data/ref_repeats.bed"
file_validator=ValidateFiles()
file_validator.validate_many(vcf_files, "vcf")
file_validator.validate_bed(repeat_regions_file)
file_validator.contigs_in_vcf(repeat_regions_file,vcf_files[0])
meta_data_file="/home/ubuntu/HandyAmpliconTool/test_data/TGC_data.csv"
metadata_utils.load_metadata(meta_data_file,",")
metadata_utils.samples_in_metadata(vcf_files)
metadata_utils.genotype_column="Final_genotype" #this will be an input

vcf_utils=VCFutilities()
master_vcf=pd.DataFrame()
for vcf in vcf_files:
    vcf_to_add=vcf_utils.load_file(vcf, )
    if master_vcf.shape[1]==0:
        master_vcf=vcf_to_add.copy()
    else:
        master_vcf=master_vcf.join(vcf_to_add, how="outer")
master_vcf.fillna("REF", inplace=True)

vcf_utils.remove_repeat_regions(master_vcf,repeat_regions_file)

hierarchy_file="/home/ubuntu/HandyAmpliconTool/test_data/genotype_hierarcy.tsv"
file_validator.validate_hierarchy(hierarchy_file)
hierarchy_utils=HierarchyUtilities()
hierarchy_utils.load_hierarchy(hierarchy_file)
target_snps=hierarchy_utils.find_defining_snps(master_vcf)

#### !!!! For testing only
print('dense : {:0.0f} bytes'.format(master_vcf.memory_usage().sum() / 1e3) )
sdf = master_vcf.astype(pd.SparseDtype("str", "REF"))
print('sparse: {:0.0f} bytes'.format(sdf.memory_usage().sum() / 1e3) )
#master_vcf=master_vcf.astype(pd.SparseDtype("str", "REF"))
#print(Counter([meta_data.get_metavalue(f,"Final_genotype") for f in samples]))

gt_snp_df: pd.DataFrame=pd.DataFrame(index=[snp for gt_snps in target_snps.values() for snp in gt_snps], columns=list(target_snps.keys())).fillna(0)
for genotype in target_snps.keys():
    gt_snp_df.loc[target_snps[genotype],genotype]=1
with open("/home/ubuntu/HandyAmpliconTool/test_data/test_gt_snps.tsv", "w") as output_file:
    output_file.write("CHR\tPOS\t"+"\t".join(gt_snp_df.columns)+"\n")
    for index in gt_snp_df.index:
        output_file.write("\t".join(list(index)+[str(f) for f in gt_snp_df.loc[[index],gt_snp_df.columns] ])+"\n")
#### END Identification of genotype defining SNPS #### 


#### Primer generation section ####
#/home/ubuntu/HandyAmpliconTool/test_data/amplified_regions.fasta
# template_sequences=[]
# with open(bed_file_name) as bed_file:
#     fasta_parser=FastaUtilities()
#     for line in bed_file:
#         chr, start, end = line.split("\t")[0:3]
#         template_sequences.append(fasta_parser.get_fasta_subseq(fasta_file_name,  chr, int(start), int(end)))
# generator=PrimerGenerator()
# for template in template_sequences:
#     ### ADD test to check that template is long enough for design of primers
#     primers=generator.generate_primers(template,(100, len(template)-100))
#     break
#### Primer generation section ####