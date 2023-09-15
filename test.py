# from inputs_validation import ValidateFiles
# from fasta_utils import FastaUtilities
# #from generate_primers import PrimerGenerator
# import name_converters
# import time 
import sys
from typing import Tuple, List, Dict
import pickle
from data_classes import Genotype, Genotypes, SNP, FlankingAmplicon, Amplicon

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
# from typing import List
# import name_converters
# import pandas as pd
# from collections import Counter
# name_converters.name_stubs.add(".sorted")

# snp_identifier=GenotypeSnpIdentifier(vcf_dir="/home/ubuntu/converted_vcfs/",
#                                      hierarchy_file="/home/ubuntu/HandyAmpliconTool/test_data/inputs/extra_genotypes_hierarchy.tsv",
#                                      meta_data_file="/home/ubuntu/HandyAmpliconTool/test_data/inputs/TGC_data.csv",
#                                      genotype_column="Final_genotype",
#                                      repeat_regions_file="/home/ubuntu/HandyAmpliconTool/test_data/inputs/ref_repeats.bed",
#                                      meta_deliminter=",",
#                                      specificity=99,
#                                      senstivity=99)
# genotypes: Genotypes = snp_identifier.identify_snps()

# genotypes.genotypes_to_snp_matrix().to_csv("/home/ubuntu/HandyAmpliconTool/test_data/outputs/test_v2.tsv", sep="\t", index=False)

# with open("/home/ubuntu/HandyAmpliconTool/test_data/outputs/test_gt_snps.pkl", "wb") as output:
#     pickle.dump(genotypes, output)

#### END Identify genotype defining SNPs ####


#### START optimise the selection of SNPs ####

import pandas as pd
from short_list_snps import SnpOptimiser

with open("/home/ubuntu/HandyAmpliconTool/test_data/outputs/test_gt_snps.pkl", "rb") as pickled_file:
    gt_snps: Genotypes = pickle.load(pickled_file)

snp_opimiser=SnpOptimiser()
max_iterval_len=1000
amplicon_intervals=snp_opimiser.optimise(max_iterval_len,gt_snps)
with open("/home/ubuntu/HandyAmpliconTool/test_data/outputs/multi_gt_intervals.bed", "w") as output_bed:
    for interval in amplicon_intervals:
        interval_start=min([f[1] for f in interval["snps"]])
        interval_end=max([f[1] for f in interval["snps"]])
        interval_len=interval_end-interval_start
        chr=interval["snps"][0][0]
        gts="_".join(interval["genotypes"])
        output_bed.write('\t'.join([chr, str(interval_start),str(interval_end),gts])+"\n")

#### END optimise the selection of SNPs ####

#### START tviD/parC species SNPs section ####
# from identify_species_snps import IdentifySpeciesSnps
# snp_identifier=IdentifySpeciesSnps(ref_fasta="/home/ubuntu/HandyAmpliconTool/test_data/inputs/GCA_000195995.1_for_VCFs.fasta",
#                                    msa_dir="/home/ubuntu/HandyAmpliconTool/test_data/msa/",
#                                    negative_genomes_dir="/home/ubuntu/HandyAmpliconTool/test_data/inputs/test_negative_genomes/genomes",
#                                    temp_blast_db_dir="/home/ubuntu/HandyAmpliconTool/test_data/tempBlastDB/",
#                                    amplicons_bed='/home/ubuntu/HandyAmpliconTool/test_data/outputs/tviD.bed')

# middle_amplicon: List[ Amplicon ]=snp_identifier.identify_insequence_snps(max_blast_length_diff=10, min_blast_identity=70)

# snps_file=open('/home/ubuntu/HandyAmpliconTool/test_data/outputs/species_snps.bed',"w")
# for amplicon in middle_amplicon:
#     for snp in amplicon.snps:
#         snp.to_file(snps_file, sep="\t")
# snps_file.close()

# print("a")
#### End tviD species SNPs section ####

#### START MSA generation section ####
# from identify_species_snps import IdentifySpeciesSnps
# snp_identifier=IdentifySpeciesSnps(ref_fasta="/home/ubuntu/HandyAmpliconTool/test_data/inputs/GCA_000195995.1_for_VCFs.fasta",
#                                    msa_dir="/home/ubuntu/HandyAmpliconTool/test_data/msa/",
#                                    negative_genomes_dir="/home/ubuntu/HandyAmpliconTool/test_data/inputs/test_negative_genomes/genomes",
#                                    temp_blast_db_dir="/home/ubuntu/HandyAmpliconTool/test_data/tempBlastDB/",
#                                    amplicons_bed='/home/ubuntu/HandyAmpliconTool/test_data/outputs/multi_gt_intervals.bed')

# flanking_amplicons: List[ FlankingAmplicon ]=snp_identifier.identify_flanking_snps(max_seq_len=1000, max_blast_length_diff=10, min_blast_identity=70)

# with open("/home/ubuntu/HandyAmpliconTool/test_data/outputs/test_species_amplicons.pkl", "wb") as output:
#     pickle.dump(flanking_amplicons, output)

# snps_file=open('/home/ubuntu/HandyAmpliconTool/test_data/outputs/species_snps.bed',"w")
# for snp in species_snps:
#     snp.to_file(snps_file, sep="\t")
# snps_file.close()
#### END MSA generation section ####

#### START list amplicons with species defining SNPs on flanks section ####

# with open("/home/ubuntu/HandyAmpliconTool/test_data/outputs/test_gt_snps.pkl", "rb") as pickled_file:
#     gts: Genotypes = pickle.load(pickled_file)

# with open("/home/ubuntu/HandyAmpliconTool/test_data/outputs/test_species_amplicons.pkl", "rb") as pickled_file:
#     flanking_amplicons: List[ FlankingAmplicon ] = pickle.load(pickled_file)

# with open("/home/ubuntu/HandyAmpliconTool/test_data/outputs/snps.bed", "w") as output_file:
#     for parent in sorted(set([f.parent for f in flanking_amplicons]), key=lambda x: x.ref_start):
#         for genotype_snp in [snp for genotype in gts.genotypes for snp in genotype.defining_snps]:
#             if parent.snp_in_amplicon(genotype_snp):
#                 parent.snps.append(genotype_snp)

#         #check the snp closest to start and end of parent amplicon
#         left_amplicon=[f for f in flanking_amplicons if f.parent.id==parent.id and f.is_left][0]
#         right_amplicon=[f for f in flanking_amplicons if f.parent.id==parent.id and not f.is_left][0]

#         if not left_amplicon.has_homologues and not right_amplicon.has_homologues:
#             #neither flank has homologue
#             line_name=f'{parent.name}_no_hom'
#             output_file.write("\t".join( [str(f) for f in  [parent.ref_contig, parent.ref_start, parent.ref_end, line_name]  ] )+"\n")

#         elif left_amplicon.has_homologues and right_amplicon.has_homologues:
#             if len(left_amplicon.snps)>0 and len(right_amplicon.snps)>0:
#                 #both flanks might SNPs
#                 left_last_snp=sorted(left_amplicon.snps, key=lambda x: x.position)[-1]
#                 right_first_snp=sorted(right_amplicon.snps, key=lambda x: x.position)[0]
#                 if right_first_snp.position-left_last_snp.position<=right_amplicon.max_len:
#                     #both flanks do have SNPs within a max_len size of amplicon
#                     line_name=f'{parent.name}_left_{left_last_snp.ref_base}/{left_last_snp.alt_base}_right_{right_first_snp.ref_base}/{right_first_snp.alt_base}'
#                     output_file.write("\t".join( [str(f) for f in  [parent.ref_contig, left_last_snp.position, right_first_snp.position, line_name]  ] )+"\n")
#                 else:
#                     #both have SNPs, but they are too far apart
#                     line_name=f'{parent.name}_left_{left_last_snp.ref_base}/{left_last_snp.alt_base}'
#                     output_file.write("\t".join( [str(f) for f in  [parent.ref_contig, left_last_snp.position, parent.ref_end, line_name]  ] )+"\n")
#                     line_name=f'{parent.name}_right_{right_first_snp.ref_base}/{right_first_snp.alt_base}'
#                     output_file.write("\t".join( [str(f) for f in  [parent.ref_contig, parent.ref_start, right_first_snp.position, line_name]  ] )+"\n")

#         elif left_amplicon.has_homologues and len(left_amplicon.snps)>0:
#             left_last_snp=sorted(left_amplicon.snps, key=lambda x: x.position)[-1]
#             line_name=f'{parent.name}_left_{left_last_snp.ref_base}/{left_last_snp.alt_base}'
#             output_file.write("\t".join( [str(f) for f in  [parent.ref_contig, left_last_snp.position, parent.ref_end, line_name]  ] )+"\n")
            
#         elif right_amplicon.has_homologues and len(right_amplicon.snps)>0:
#             right_first_snp=sorted(right_amplicon.snps, key=lambda x: x.position)[0]
#             line_name=f'{parent.name}_right_{right_first_snp.ref_base}/{right_first_snp.alt_base}'
#             output_file.write("\t".join( [str(f) for f in  [parent.ref_contig, parent.ref_start, right_first_snp.position, line_name]  ] )+"\n")


# print("a")
#### END list amplicons with species defining SNPs on flanks section ####




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



