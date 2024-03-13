#from typing import Tuple, List, Dict
import pickle
from os.path import exists, expanduser
from os import mkdir
from shutil import which
from sys import exit
from data_classes import Genotype, Genotypes, InputConfiguration
from inputs_validation import ValidateFiles
from load_vcfs import VCFutilities
from identify_genotype_snps import GenotypeSnpIdentifier
import name_converters
from snp_optimiser import SnpOptimiser
from identify_species_snps import IdentifySpeciesSnps
from primers_generator import PrimersGenerator
import argparse

parser = argparse.ArgumentParser(description='Generate list of SNPs that uniquely identify one or more genotypes')
parser.add_argument('-c','--config_file', metavar='', type=str,
                    help='Config file, see sample.json for example.', required=True)
parser.add_argument('-m','--mode', metavar='', type=str, choices=['SNP', 'Amplicon'],
                    help='"SNP" to only get genotype defining SNPs, "Amplicon" to also generate amplicons', required=True)

try:
    args = parser.parse_args()
except:
    parser.print_help()
    exit(0)

run_mode=args.mode
config_file=expanduser(args.config_file)
def main():

    #Check dependencies

    if which("mafft") is None:
        raise OSError("Missing MAFFT program. It's available via Conda.")
    if which("makeblastdb") is None or which("blastn") is None:
        raise OSError("Missing NCBI BLAST program. It's available via Conda.")

    #switch to argparse later
    config_data = InputConfiguration( config_file )

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
                    fourth_col=f'SNP:{snp.ref_base}/{snp.alt_base}/GT:{genotype.name}/SP:{snp.specificity:.2f}/SE:{snp.sensitivity:.2f}'
                    snps_file.write("\t".join( [str(f) for f in [snp.ref_contig_id, snp.position, snp.position+1, fourth_col] ])+"\n")

    with open(config_data.genotypes_data, "wb") as output:
        pickle.dump(genotypes, output)


    with open(config_data.genotypes_data, "rb") as pickled_file:
        gt_snps: Genotypes = pickle.load(pickled_file)

    with open(config_data.genotypes_data, "rb") as pickled_file:
            gt_snps: Genotypes = pickle.load(pickled_file)


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

    if run_mode!="Amplicon":
        vcf_utils=VCFutilities()
        vcf_utils.output_genotypes_vcf(genotypes, config_data.snps_vcf)
        exit(0)
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

    with open(config_data.species_data, "rb") as pickled_file:
        species: Genotype = pickle.load(pickled_file)

    with open(config_data.genotypes_data, "rb") as pickled_file:
        genotypes: Genotypes = pickle.load(pickled_file)

    genotypes.genotypes.append(species)

    vcf_utils=VCFutilities()
    vcf_utils.output_species_vcf(genotypes, config_data.snps_vcf)

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


main()

