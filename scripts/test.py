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

#### START Input validation tests ####
print("Validating input files")
from inputs_validation import ValidateFiles
file_validator=ValidateFiles()
#fasta_file_name=f'{root_dir}/HandyAmpliconTool/test_data/test_data.fasta'
fasta_file_name=f'{root_dir}/HandyAmpliconTool/test_data/inputs/GCA_000195995.1_for_VCFs.fasta'
bed_file_name=f'{root_dir}/HandyAmpliconTool/test_data/inputs/ref_repeats.bed'
#vcf_file=f'{root_dir}/HandyAmpliconTool/test_data/vcfs/8490_5#12.vcf'
vcf_dir=f'{root_dir}/converted_vcfs/'
file_validator.validate_bed(bed_file_name)
file_validator.validate_fasta(fasta_file_name)
file_validator.contigs_in_fasta(bed_file_name, fasta_file_name)
file_validator.validate_vcf(vcf_dir)
#### END Input validation tests ####

#### START Identify genotype defining SNPs ####
from identify_genotype_snps import GenotypeSnpIdentifier
from typing import List
import name_converters
import pandas as pd
from collections import Counter
name_converters.name_stubs.add(".sorted")

snp_identifier=GenotypeSnpIdentifier(vcf_dir=vcf_dir,
                                     hierarchy_file=f'{root_dir}/HandyAmpliconTool/test_data/inputs/genotype_hierarcy.tsv',
                                     meta_data_file=f'{root_dir}/HandyAmpliconTool/test_data/inputs/TGC_data.csv',
                                     genotype_column="Final_genotype",
                                     repeat_regions_file=f'{root_dir}/HandyAmpliconTool/test_data/inputs/ref_repeats.bed',
                                     meta_deliminter=",",
                                     specificity=98,
                                     senstivity=98)
genotypes: Genotypes = snp_identifier.identify_snps()

genotypes.genotypes_to_snp_matrix().to_csv(f'{root_dir}/HandyAmpliconTool/test_data/outputs/test_v2.tsv', sep="\t", index=False)

with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/genotype_snps.tsv', "w") as snps_file:
    for genotype in genotypes.genotypes:
        for snp in genotype.defining_snps:
            if snp.passes_filters:
                fourth_col=f'SNP:{snp.ref_base}/{snp.alt_base}/GT:{genotype.name}/SP:{snp.specificity:.1f}/SE:{snp.sensitivity:.1f}'
                snps_file.write("\t".join( [str(f) for f in [snp.ref_contig_id, snp.position, snp.position+1, fourth_col] ])+"\n")
sys.exit()
with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/test_gt_snps.pkl', "wb") as output:
    pickle.dump(genotypes, output)

#### END Identify genotype defining SNPs ####

#### START add a bespoke genotype (AMR) ####

# extra_genotype=Genotype("parC")
# chr="AL513382_143N_pHCM1_120N_pHCM2"
# for position, ref_base in zip([3196469, 3196470, 3196471, 3196457, 3196458, 3196459], ["G","C","T","T","T","C"]):
#     new_snp=SNP(ref_contig_id=chr, ref_base=ref_base, alt_base="A", position=position, passes_filters=True)
#     new_snp.is_genotype_snp=True
#     extra_genotype.add_genotype_allele(new_snp, "A", 1)

#### END add a bespoke genotype (AMR) ####

#### START optimise the selection of SNPs ####

import pandas as pd
from short_list_snps import SnpOptimiser

with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/test_gt_snps.pkl', "rb") as pickled_file:
    gt_snps: Genotypes = pickle.load(pickled_file)

#gt_snps.genotypes.append(extra_genotype)

snp_opimiser=SnpOptimiser()
max_iterval_len=1000
amplicon_intervals=snp_opimiser.optimise(max_iterval_len,gt_snps, rare_gts=["4.3.1.1.P1"])
with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/multi_gt_intervals.bed', "w") as output_bed:
    for interval in amplicon_intervals:
        interval_start=min([f[1] for f in interval["snps"]])
        interval_end=max([f[1] for f in interval["snps"]])
        interval_len=interval_end-interval_start
        chr=interval["snps"][0][0]
        gts="_".join(interval["genotypes"])
        output_bed.write('\t'.join([chr, str(interval_start),str(interval_end),gts])+"\n")

# # ### END optimise the selection of SNPs ####

# # #### START MSA generation section ####
from identify_species_snps import IdentifySpeciesSnps
snp_identifier=IdentifySpeciesSnps(ref_fasta=f'{root_dir}/HandyAmpliconTool/test_data/inputs/GCA_000195995.1_for_VCFs.fasta',
                                   msa_dir=f'{root_dir}/HandyAmpliconTool/test_data/msa/',
                                   negative_genomes_dir=f'{root_dir}/HandyAmpliconTool/test_data/inputs/test_negative_genomes/genomes',
                                   temp_blast_db_dir=f'{root_dir}/HandyAmpliconTool/test_data/tempBlastDB/',
                                   amplicons_bed=f'{root_dir}/HandyAmpliconTool/test_data/outputs/multi_gt_intervals.bed')

flanking_amplicons: Genotype=snp_identifier.identify_flanking_snps(max_seq_len=1000, max_blast_length_diff=90, min_blast_identity=60)

with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/test_species_amplicons.pkl', "wb") as output:
    pickle.dump(flanking_amplicons, output)


#### END MSA generation section ####


#### START write VCF ####

with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/test_species_amplicons.pkl', "rb") as pickled_file:
    species: Genotype = pickle.load(pickled_file)

with open(f'{root_dir}/HandyAmpliconTool/test_data/outputs/test_gt_snps.pkl', "rb") as pickled_file:
    genotypes: Genotypes = pickle.load(pickled_file)

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
    for contig_id, position in coordinates:
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
                    vcf_snp_id="_".join( ["Serovar", snp.ref_contig_id ,str(snp.position+1), snp.alt_base] )
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



