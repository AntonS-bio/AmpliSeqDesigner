# from inputs_validation import ValidateFiles
# from fasta_utils import FastaUtilities
# #from generate_primers import PrimerGenerator
# import name_converters
# import time 
import sys
# from multiprocessing import cpu_count, pool
# from os import listdir

#### START Input validation tests ####
# from inputs_validation import ValidateFiles
# file_validator=ValidateFiles()
# #fasta_file_name="/home/ubuntu/HandyAmpliconTool/test_data/test_data.fasta"
# fasta_file_name="/home/ubuntu/HandyAmpliconTool/test_data/inputs/GCA_000195995.1_for_VCFs.fasta"
# bed_file_name="/home/ubuntu/HandyAmpliconTool/test_data/inputs/test_data.bed"
# vcf_file="/home/ubuntu/HandyAmpliconTool/test_data/vcfs/8490_5#12.vcf"
# file_validator.validate_bed(bed_file_name)
# file_validator.validate_fasta(fasta_file_name)
# file_validator.contigs_in_fasta(bed_file_name, fasta_file_name)
# file_validator.validate_vcf(vcf_file)
#### END Input validation tests ####

#### START Identify genotype defining SNPs ####
# from identify_genotype_snps import GenotypeSnpIdentifier
# import name_converters
# import pandas as pd
# from collections import Counter
# name_converters.name_stubs.add(".sorted")

# snp_identifier=GenotypeSnpIdentifier(vcf_dir="/home/ubuntu/converted_vcfs/",
#                                      hierarchy_file="/home/ubuntu/HandyAmpliconTool/test_data/inputs/genotype_hierarcy.tsv",
#                                      meta_data_file="/home/ubuntu/HandyAmpliconTool/test_data/inputs/TGC_data.csv",
#                                      genotype_column="Final_genotype",
#                                      repeat_regions_file="/home/ubuntu/HandyAmpliconTool/test_data/inputs/ref_repeats.bed",
#                                      meta_deliminter=",",
#                                      specificity=90.5,
#                                      senstivity=90.5)
# gt_snp_df: pd.DataFrame = snp_identifier.identify_snps()

# snp_to_drop=[index for index in gt_snp_df.index if True not in gt_snp_df.loc[ [index] ].values]
# gt_snp_df.drop(index=snp_to_drop, inplace=True)


# ### make output tsv pretty, as pandas can't really output or load multi index dataframes via tsv
# with open("/home/ubuntu/HandyAmpliconTool/test_data/outputs/test_gt_snps.tsv","w") as ouput_file:
#     ouput_file.write("\t".join( [f for f in ["Contig","Position"]+list(gt_snp_df.columns) ] )+"\n")
#     for i, index in enumerate(gt_snp_df.index):
#         ouput_file.write("\t".join([str(f) for f in [index[0],index[1]]+list( gt_snp_df.loc[ [index] ].values[0]) ])+"\n")

# gt_snp_df.to_pickle("/home/ubuntu/HandyAmpliconTool/test_data/outputs/test_gt_snps.pkl")

#### END Identify genotype defining SNPs ####

#sys.exit()

#### START optimise the selection of SNPs ####
# import pandas as pd
# from short_list_snps import SnpOptimiser

# gt_snp_df=pd.read_pickle("/home/ubuntu/HandyAmpliconTool/test_data/outputs/test_gt_snps.pkl")

# snp_opimiser=SnpOptimiser(gt_snp_df)
# max_iterval_len=1000
# amplicon_intervals=snp_opimiser.optimise(max_iterval_len,gt_snp_df)
# with open("/home/ubuntu/HandyAmpliconTool/test_data/outputs/multi_gt_intervals.bed", "w") as output_bed:
#     for interval in amplicon_intervals:
#         interval_start=min([f[1] for f in interval["snps"]])
#         interval_end=max([f[1] for f in interval["snps"]])
#         interval_len=interval_end-interval_start
#         chr=interval["snps"][0][0]
#         gts="_".join(interval["genotypes"])
#         output_bed.write('\t'.join([chr, str(interval_start),str(interval_end),gts])+"\n")

#### END optimise the selection of SNPs ####


#### START MSA generation section ####
from identify_species_snps import IdentifySpeciesSnps
snp_identifier=IdentifySpeciesSnps(ref_fasta="/home/ubuntu/HandyAmpliconTool/test_data/inputs/GCA_000195995.1_for_VCFs.fasta",
                                   msa_dir="/home/ubuntu/HandyAmpliconTool/test_data/msa/",
                                   negative_genomes_dir="/home/ubuntu/HandyAmpliconTool/test_data/inputs/test_negative_genomes/genomes",
                                   temp_blast_db_dir="/home/ubuntu/HandyAmpliconTool/test_data/tempBlastDB/",
                                   amplicons_bed='/home/ubuntu/HandyAmpliconTool/test_data/outputs/multi_gt_intervals.bed')
snp_identifier.identify_snps()
#### START MSA generation section ####



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



