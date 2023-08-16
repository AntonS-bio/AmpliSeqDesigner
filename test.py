from inputs_validation import ValidateFiles
from fasta_utils import FastaUtilities
#from generate_primers import PrimerGenerator
import name_converters
import time 
import sys
from multiprocessing import cpu_count, pool
from os import listdir

vcf_dir="/home/ubuntu/converted_vcfs"
vcf_files=[f'{vcf_dir}/{f}' for f in listdir(vcf_dir) ]
print(len(vcf_files))
valid=0
for file in vcf_files:
    try:
        with open(file) as temp_file:
            if temp_file.readline()[0]=="#":
                valid+=1
    except:
        pass
print(f'Valid: {valid}')
print(f'Invalid: {len(vcf_files)-valid}')

sys.exit()

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



#### START optimise the selection of SNPs ####




#### END optimise the selection of SNPs ####

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