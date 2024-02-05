#from typing import Tuple, List, Dict
import pickle
from os.path import exists
from os import mkdir
from shutil import which
from data_classes import Genotype, Genotypes, InputConfiguration
from inputs_validation import ValidateFiles
from identify_genotype_snps import GenotypeSnpIdentifier
import name_converters
from snp_optimiser import SnpOptimiser
from identify_species_snps import IdentifySpeciesSnps
from primers_generator import PrimersGenerator

#Check dependencies

if which("mafft") is None:
    raise OSError("Missing MAFFT program. It's available via Conda.")
if which("makeblastdb") is None or which("blastn") is None:
    raise OSError("Missing NCBI BLAST program. It's available via Conda.")

#switch to argparse later
config_data = InputConfiguration("/home/lshas17/HandyAmpliconTool/test_data/configs/2.3.1.json")
#config_data = InputConfiguration("/home/lshas17/HandyAmpliconTool/test_data/configs/3.1.1.json")

#### START Input validation tests ####
print("Validating input files")
file_validator=ValidateFiles()
file_validator.validate_bed(config_data.repeats_bed_file)
file_validator.validate_fasta(config_data.reference_fasta)
file_validator.contigs_in_fasta(config_data.repeats_bed_file, config_data.reference_fasta)
file_validator.validate_vcf(config_data.vcf_dir)
file_validator.validate_hierarchy(config_data.hierarchy_file, config_data)
#### END Input validation tests ####

#create directory for output is does not exit
for value in [config_data.output_dir, config_data.msa_dir, config_data.temp_blast_db]:
    if not exists(value):
        mkdir(value)

#### START Identify genotype defining SNPs ####
for stub in config_data.name_stubs:
    name_converters.name_stubs.add(stub)

snp_identifier=GenotypeSnpIdentifier(config_data)

genotypes: Genotypes = snp_identifier.identify_snps()

genotypes.genotypes_to_snp_matrix().to_csv(config_data.genotype_snps, sep="\t", index=False)

with open(config_data.genotype_snps, "w") as snps_file:
    for genotype in genotypes.genotypes:
        for snp in genotype.defining_snps:
            if snp.passes_filters:
                fourth_col=f'SNP:{snp.ref_base}/{snp.alt_base}/GT:{genotype.name}/SP:{snp.specificity:.1f}/SE:{snp.sensitivity:.1f}'
                snps_file.write("\t".join( [str(f) for f in [snp.ref_contig_id, snp.position, snp.position+1, fourth_col] ])+"\n")

with open(config_data.genotypes_data, "wb") as output:
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


with open(config_data.genotypes_data, "rb") as pickled_file:
    gt_snps: Genotypes = pickle.load(pickled_file)

#gt_snps.genotypes.append(extra_genotype)

snp_opimiser=SnpOptimiser()
max_iterval_len=1000
amplicon_intervals=snp_opimiser.optimise(max_iterval_len,gt_snps, rare_gts=config_data.gts_with_few_snps)
with open(config_data.multi_gt_intervals, "w") as output_bed:
    for interval in amplicon_intervals:
        interval_start=min([f[1] for f in interval["snps"]])
        interval_end=max([f[1] for f in interval["snps"]])
        #interval_len=interval_end-interval_start
        contig_id=interval["snps"][0][0]
        gts="_".join(interval["genotypes"])
        output_bed.write('\t'.join([contig_id, str(interval_start),str(interval_end),gts])+"\n")

# # ### END optimise the selection of SNPs ####

# # #### START MSA generation section ####

snp_identifier=IdentifySpeciesSnps(ref_fasta=config_data.reference_fasta,
                                   msa_dir=config_data.msa_dir,
                                   negative_genomes_dir=config_data.negative_genomes,
                                   temp_blast_db_dir=config_data.temp_blast_db,
                                   amplicons_bed=config_data.multi_gt_intervals)


species_genotype: Genotype=snp_identifier.generate_flanking_amplicons()
flanking_amplicons=snp_identifier.get_bifurcating_snps(species_genotype)


with open(config_data.species_data, "wb") as output:
    pickle.dump(flanking_amplicons, output)


#### END MSA generation section ####


#### START write VCF ####

with open(config_data.species_data, "rb") as pickled_file:
    species: Genotype = pickle.load(pickled_file)

with open(config_data.genotypes_data, "rb") as pickled_file:
    genotypes: Genotypes = pickle.load(pickled_file)

genotypes.genotypes.append(species)

#### START Generate primers ####
generator=PrimersGenerator(config_data)
generator._for_testing_load_gts(config_data.species_data, config_data.genotypes_data)
generator.find_candidate_primers()
with open(config_data.output_dir+"primers.tsv","w") as output_file:
    header="\t".join(["Name", "Penalty", "Contig", "Start","End","Length",
                    "Forward","Forward Tm", "Forward GC",
                    "Reverse","Reverse Tm", "Reverse GC"])+"\n"
    output_file.write(header)
    for pair in generator.new_primer_pairs:
        output_file.write(pair.to_string()+"\n")


# with open(config_data.snps_vcf, "w") as vcf_output_file:
#     #write the vcf header
#     vcf_output_file.write('##fileformat=VCFv4.2'+"\n")
#     vcf_output_file.write('##FILTER=<ID=PASS,Description="All filters passed">'+"\n")
#     vcf_output_file.write('##ALT=<ID=*,Description="Represents allele(s) other than observed.">'+"\n")

#     for ref_contig in set([snp.ref_contig_id for genotype in genotypes.genotypes for snp in genotype.defining_snps]):
#         contig_max_position=max([snp.position for genotype in genotypes.genotypes for snp in genotype.defining_snps if snp.ref_contig_id == ref_contig])
#         vcf_output_file.write(f'##contig=<ID={ref_contig},length={str(contig_max_position)}>'+"\n")
#     header_line="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
#     gt_columns={}
#     for i, gt in enumerate(genotypes.genotypes):
#         header_line=header_line+gt.name+"\t"
#         gt_columns[gt.name]=i
#     header_line=header_line+"NonTargetSerovar\n"
#     gt_columns["NonTargetSerovar"]=len(gt_columns)
#     vcf_output_file.write(header_line)

#     coordinates=sorted(set([coordinate for genotype in genotypes.genotypes for coordinate in genotype.defining_snp_coordinates]))
#     for contig_id, position in coordinates:
#         snps_at_coordinates=[(genotype, snp) for genotype in genotypes.genotypes for snp in genotype.defining_snps if snp.position==position and snp.ref_contig_id==contig_id]
#         alt_alleles=[base for base in [snp[1].alt_base for snp in snps_at_coordinates] ][0]
#         if len(alt_alleles)>1:
#             print(f'Excess alleles at pos: {str(position)} contig {contig_id}')
#             continue
#         alt_str=f'{alt_alleles}'
#         for genotype, snp in snps_at_coordinates:
#             #check that snp is in multi genotype region
#             if True not in set([f.snp_in_amplicon(snp) for f in species.amplicons]):
#                 continue

#             if snp.passes_filters:
#                 if genotype.name!="species":
#                     if genotype.get_genotype_allele(snp)==snp.alt_base:
#                         suffix=["1:."]*len(gt_columns)
#                         suffix[gt_columns[genotype.name]]="1:"+str(genotype.get_genotype_allele_depth(snp)) #0 is REF allele
#                     else:
#                         suffix=["1:."]*len(gt_columns) #set
#                         suffix[gt_columns[genotype.name]]="0:"+str(genotype.get_genotype_allele_depth(snp))
#                     vcf_snp_id="_".join( ["GT",genotype.name,snp.ref_contig_id,str(snp.position+1)] )
#                     suffix[gt_columns["species"]]=0
#                     suffix[gt_columns["NonTargetSerovar"]]=".:."
#                 else:
#                     vcf_snp_id="_".join( ["Serovar", snp.ref_contig_id ,str(snp.position+1), snp.alt_base] )
#                     suffix=["1:."]*len(gt_columns)
#                     suffix[gt_columns["NonTargetSerovar"]]="1:"+str(genotypes.genotypes[-1].get_genotype_allele_depth(snp))
#                     alt_str=snp.alt_base
#                 vcf_output_file.write("\t".join([str(f) for f in [snp.ref_contig_id,
#                                                                     snp.position+1,
#                                                                     vcf_snp_id,
#                                                                     snp.ref_base,
#                                                                     alt_str,
#                                                                     ".",
#                                                                     "PASS",
#                                                                     ".",
#                                                                     "GT:DP",
#                                                                     ]+suffix ]   ) +"\n" )

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



