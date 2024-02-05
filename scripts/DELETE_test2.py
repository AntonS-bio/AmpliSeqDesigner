# from inputs_validation import ValidateFiles
# from fasta_utils import FastaUtilities
# #from generate_primers import PrimerGenerator
# import name_converters
# import time 
import sys
from typing import Tuple, List, Dict
import pickle
from data_classes import Genotype, Genotypes, SNP, FlankingAmplicon, Amplicon

root_dir="/home/lshas17/"

#### START check that prerequisites exist ####
from shutil import which
if which("mafft") is None:
    raise OSError("Missing MAFFT program. It's available via Conda.")
if which("makeblastdb") is None or which("blastn") is None:
    raise OSError("Missing NCBI BLAST program. It's available via Conda.")

#### END check that prerequisites exist ####

#### START Input validation tests ####
# print("Validating input files")
# from inputs_validation import ValidateFiles
# file_validator=ValidateFiles()
# #fasta_file_name=f'{root_dir}/HandyAmpliconTool/test_data/test_data.fasta'
# fasta_file_name=f'{root_dir}/HandyAmpliconTool/test_data/inputs/GCA_000195995.1_for_VCFs.fasta'
# bed_file_name=f'{root_dir}/HandyAmpliconTool/test_data/inputs/ref_repeats.bed'
# vcf_file=f'{root_dir}/HandyAmpliconTool/test_data/vcfs/8490_5#12.vcf'
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

# snp_identifier=GenotypeSnpIdentifier(vcf_dir=f'{root_dir}/converted_vcfs/',
#                                      hierarchy_file=f'{root_dir}/HandyAmpliconTool/test_data/inputs/extra_genotypes_hierarchy.tsv',
#                                      meta_data_file=f'{root_dir}/HandyAmpliconTool/test_data/inputs/TGC_data.csv',
#                                      genotype_column="Final_genotype",
#                                      repeat_regions_file=f'{root_dir}/HandyAmpliconTool/test_data/inputs/ref_repeats.bed',
#                                      meta_deliminter=",",
#                                      specificity=99,
#                                      senstivity=99)
# genotypes: Genotypes = snp_identifier.identify_snps()

# genotypes.genotypes_to_snp_matrix().to_csv(f'{root_dir}/HandyAmpliconTool/test_data/outputs/test_v2.tsv', sep="\t", index=False)

# with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/genotype_snps.tsv', "w") as snps_file:
#     for genotype in genotypes.genotypes:
#         for snp in genotype.defining_snps:
#             if snp.passes_filters:
#                 fourth_col=f'SNP:{snp.ref_base}/{snp.alt_base}/GT:{genotype.name}/SP:{snp.specificity:.1f}/SE:{snp.sensitivity:.1f}'
#                 snps_file.write("\t".join( [str(f) for f in [snp.ref_contig_id, snp.position, snp.position+1, fourth_col] ])+"\n")

# with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/test_gt_snps.pkl', "wb") as output:
#     pickle.dump(genotypes, output)

#### END Identify genotype defining SNPs ####


#### START optimise the selection of SNPs ####

# import pandas as pd
# from short_list_snps import SnpOptimiser

# with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/test_gt_snps.pkl', "rb") as pickled_file:
#     gt_snps: Genotypes = pickle.load(pickled_file)

# snp_opimiser=SnpOptimiser()
# max_iterval_len=1000
# amplicon_intervals=snp_opimiser.optimise(max_iterval_len,gt_snps, rare_gts=["4.3.1.1.P1"])
# with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/multi_gt_intervals.bed', "w") as output_bed:
#     for interval in amplicon_intervals:
#         interval_start=min([f[1] for f in interval["snps"]])
#         interval_end=max([f[1] for f in interval["snps"]])
#         interval_len=interval_end-interval_start
#         chr=interval["snps"][0][0]
#         gts="_".join(interval["genotypes"])
#         output_bed.write('\t'.join([chr, str(interval_start),str(interval_end),gts])+"\n")

### END optimise the selection of SNPs ####

#### START MSA generation section ####
# from identify_species_snps import IdentifySpeciesSnps
# snp_identifier=IdentifySpeciesSnps(ref_fasta=f'{root_dir}/HandyAmpliconTool/test_data/inputs/GCA_000195995.1_for_VCFs.fasta',
#                                    msa_dir=f'{root_dir}/HandyAmpliconTool/test_data/msa/',
#                                    negative_genomes_dir=f'{root_dir}/HandyAmpliconTool/test_data/inputs/test_negative_genomes/genomes',
#                                    temp_blast_db_dir=f'{root_dir}/HandyAmpliconTool/test_data/tempBlastDB/',
#                                    amplicons_bed=f'{root_dir}/HandyAmpliconTool/test_data/inputs/current_amplicons.bed')

# flanking_amplicons: Genotype=snp_identifier.identify_flanking_snps(max_seq_len=0, max_blast_length_diff=90, min_blast_identity=60)

# with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/test_species_amplicons.pkl', "wb") as output:
#     pickle.dump(flanking_amplicons, output)

# sys.exit()
#### END MSA generation section ####


#### START write VCF ####

with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/test_species_amplicons.pkl', "rb") as pickled_file:
    species: Genotype = pickle.load(pickled_file)

genotypes=Genotypes()
genotypes.genotypes.append(species)

with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/snps.vcf', "w") as vcf_output_file:
    #write the vcf header
    vcf_output_file.write('##fileformat=VCFv4.2'+"\n")
    vcf_output_file.write('##FILTER=<ID=PASS,Description="All filters passed">'+"\n")
    vcf_output_file.write('##ALT=<ID=*,Description="Represents allele(s) other than observed.">'+"\n")

    for ref_contig in set([snp.ref_contig_id for genotype in genotypes.genotypes for snp in genotype.defining_snps]):
        contig_max_position=max([snp.position for genotype in genotypes.genotypes for snp in genotype.defining_snps if snp.ref_contig_id == ref_contig])
        vcf_output_file.write(f'##contig=<ID={ref_contig},length={str(contig_max_position)}>'+"\n")
    header_line="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
    gt_columns={}
    for i, gt in enumerate(genotypes.genotypes):
        header_line=header_line+gt.name+"\t"
        gt_columns[gt.name]=i
    header_line=header_line+"NonTargetSerovar\n"
    gt_columns["NonTargetSerovar"]=len(gt_columns)
    vcf_output_file.write(header_line)
   
    coordinates=sorted(set([coordinate for genotype in genotypes.genotypes for coordinate in genotype.defining_snp_coordinates]))
    species_amplicons=[genotype for genotype in genotypes.genotypes if genotype.name=="species"][0].amplicons
    #species_amplicons=[amplicon for amplicon in species_amplicons if amplicon.has_flanking]
    for contig_id, position in coordinates:
        if  len([amplicon for amplicon in species_amplicons if amplicon.coord_in_amplicon([contig_id, position])])==0:
            continue
        snps_at_coordinates=[(genotype, snp) for genotype in genotypes.genotypes for snp in genotype.defining_snps if snp.position==position and snp.ref_contig_id==contig_id]
        alt_alleles=[base for base in [snp[1].alt_base for snp in snps_at_coordinates] ][0]
        if len(alt_alleles)>1:
            print(f'Excess alleles at pos: {str(position)} contig {contig_id}')
            continue
        alt_str=f'{alt_alleles}'
        for genotype, snp in snps_at_coordinates:
            #check that snp is in multi genotype region
            if True not in set([f.snp_in_amplicon(snp) for f in species.amplicons]):
                continue

            if snp.passes_filters:
                if genotype.name!="species":
                    if genotype.get_genotype_allele(snp)==snp.alt_base:
                        suffix=["1:."]*len(gt_columns)
                        suffix[gt_columns[genotype.name]]="1:"+str(genotype.get_genotype_allele_depth(snp)) #0 is REF allele
                    else:
                        suffix=["1:."]*len(gt_columns) #set
                        suffix[gt_columns[genotype.name]]="0:"+str(genotype.get_genotype_allele_depth(snp))
                    vcf_snp_id="_".join( ["GT",genotype.name,snp.ref_contig_id,str(snp.position+1)] )
                    suffix[gt_columns["species"]]=0
                    suffix[gt_columns["NonTargetSerovar"]]=".:."
                else:
                    vcf_snp_id=[amplicon for amplicon in species_amplicons if amplicon.coord_in_amplicon([contig_id, position])][0].name
                    #vcf_snp_id="_".join( ["Serovar", snp.ref_contig_id ,str(snp.position+1), snp.alt_base] )
                    suffix=["1:."]*len(gt_columns)
                    suffix[gt_columns["NonTargetSerovar"]]="1:"+str(genotypes.genotypes[-1].get_genotype_allele_depth(snp))
                    alt_str=snp.alt_base
                vcf_output_file.write("\t".join([str(f) for f in [snp.ref_contig_id,
                                                                    snp.position+1,
                                                                    vcf_snp_id,
                                                                    snp.ref_base,
                                                                    alt_str,
                                                                    ".",
                                                                    "PASS",
                                                                    ".",
                                                                    "GT:DP",
                                                                    ]+suffix ]   ) +"\n" )

#### END write VCF ####


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



