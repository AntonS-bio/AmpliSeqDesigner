from inputs_validation import ValidateFiles
from fasta_utils import FastaUtilities
from generate_primers import PrimerGenerator
import name_converters

#### Input validation tests ####
# file_validator=ValidateFiles()
# #fasta_file_name="/home/ubuntu/HandyAmpliconTool/test_data/test_data.fasta"
# fasta_file_name="/home/ubuntu/HandyAmpliconTool/test_data/GCF_000195995.1.fna"
# bed_file_name="/home/ubuntu/HandyAmpliconTool/test_data/test_data.bed"
# vcf_file="/home/ubuntu/HandyAmpliconTool/test_data/vcfs/8490_5#12.vcf"
# file_validator.validate_bed(bed_file_name)
# file_validator.validate_fasta(fasta_file_name)
# file_validator.check_contigs_in_fasta(bed_file_name, fasta_file_name)
# file_validator.validate_vcf(vcf_file)
#### Input validation tests ####

#### Identification of genotype defining SNPS #### 
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


vcf_dir="/home/ubuntu/HandyAmpliconTool/test_data/vcfs/"
vcf_files=[f'{vcf_dir}{f}' for f in listdir(vcf_dir)]
file_validator=ValidateFiles()
file_validator.validate_many(vcf_files, "vcf")

meta_data_file="/home/ubuntu/HandyAmpliconTool/test_data/TGC_data.csv"
#meta_data=MetadataUtilities(meta_data_file,",")
metadata_utils.load_metadata(meta_data_file,",")
metadata_utils.samples_in_metadata(vcf_files)
metadata_utils.genotype_column="Final_genotype" #this will be an input

#sys.exit()

vcf_utils=VCFutilities()
master_vcf=pd.DataFrame()
for vcf in vcf_files:
    vcf_to_add=vcf_utils.load_file(vcf, )
    if master_vcf.shape[1]==0:
        master_vcf=vcf_to_add.copy()
    else:
        master_vcf=master_vcf.join(vcf_to_add, how="outer")
    #print(master_vcf.shape)
master_vcf.fillna("REF", inplace=True)
print("end")

print('dense : {:0.2f} bytes'.format(master_vcf.memory_usage().sum() / 1e3) )
sdf = master_vcf.astype(pd.SparseDtype("str", "REF"))
print('sparse: {:0.2f} bytes'.format(sdf.memory_usage().sum() / 1e3) )
samples=[name_converters.get_sample(f)  for f in sdf.columns]
#print(Counter([meta_data.get_metavalue(f,"Final_genotype") for f in samples]))

hierarchy_file="/home/ubuntu/HandyAmpliconTool/test_data/genotype_hierarcy.tsv"
file_validator.validate_hierarchy(hierarchy_file)
hierarchy_utils=HierarchyUtilities()
hierarchy_utils.load_hierarchy(hierarchy_file)
target_snps=hierarchy_utils.find_defining_snps(master_vcf)
a=1

# #identify alleles that occur in all subgenotypes of given genotype, but nowhere else
# gt_specific_alleles={}
# for gt, sub_gt in gt_subgenotypes.items():
#     gt_specific_alleles[gt]=list(alleles_df.loc[ (np.sum(alleles_df[sub_gt], axis=1)==len(sub_gt)) & (np.sum(alleles_df[gt_complements[gt]], axis=1)==0), "Pos" ])

# #select the required genotypes to type
# target_gt=["1", "2","2.5","2.2.2","2.3.2","3","3.3","3.1.1","3.3.1","4",
#            "4.3.1","4.3.1.1","4.3.1.2","4.3.1.2.1","4.3.1.2.1.1","4.3.1.1.P1"]
# #target_gt=["2","3","4"]


# gt_allele_df=pd.DataFrame(index= list(set([index for sublist in  gt_specific_alleles.values() for index in sublist ])), columns=target_gt , dtype=int  ).fillna(0)
# for gt in target_gt:
#     gt_allele_df.loc[gt_specific_alleles[gt],gt]=1
# gt_allele_df=gt_allele_df.loc[np.sum(gt_allele_df,axis=1)>0]
# gt_allele_df.to_csv("gt_specific_alleles_v2.tsv",sep="\t")
# np.sum(gt_allele_df, axis=0)







#for each file in vcf_dir, get genotype from metadat




#### Identification of genotype defining SNPS ####

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


