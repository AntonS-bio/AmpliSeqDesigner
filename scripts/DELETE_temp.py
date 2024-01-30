import unittest
import inputs_validation
from identify_genotype_snps import GenotypeSnpIdentifier
import name_converters
from data_classes import Genotype, Genotypes, InputConfiguration

config_data = InputConfiguration("~/HandyAmpliconTool/unit_test_data/unittest.json")

for stub in config_data.name_stubs:
    name_converters.name_stubs.add(stub)

snp_identifier=GenotypeSnpIdentifier(vcf_dir=config_data.vcf_dir,
                                    hierarchy_file=config_data.hierarchy_file,
                                    meta_data_file=config_data.meta_data_file,
                                    genotype_column=config_data.genotype_column,
                                    repeat_regions_file=config_data.repeats_bed_file,
                                    meta_deliminter=config_data.metadata_delim,
                                    specificity=config_data.snp_sensitivity,
                                    senstivity=config_data.snp_specificity)
genotypes: Genotypes = snp_identifier.identify_snps()