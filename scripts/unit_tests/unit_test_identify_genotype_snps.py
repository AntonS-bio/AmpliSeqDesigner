from os.path import expanduser, realpath, dirname
unit_test_dir = dirname(realpath(__file__))
from sys import path
path.append( f'{unit_test_dir}/..')
print(unit_test_dir)

import unittest
from identify_genotype_snps import GenotypeSnpIdentifier
import name_converters
from data_classes import InputConfiguration, Genotypes
import metadata_utils
from snp_optimiser import SnpOptimiser
from os import listdir

import pickle

class TestGenotypSnps(unittest.TestCase):
    valid_data=expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/")
    invalid_data=expanduser("~/HandyAmpliconTool/unit_test_data/invalid_data/")
    config_file=expanduser("~/HandyAmpliconTool/unit_test_data/unittest.json")

    def setUp(self) -> None:
        self.config_data = InputConfiguration(self.config_file)
        return super().setUp()

    def test_metadata_parsing(self):
        correct_gt_column=self.config_data.genotype_column
        self.config_data._config_data["metadata_parameters"]["genotype_column"]="NoneSuch"
        self.assertRaises(ValueError,metadata_utils.load_metadata, 
                        config=self.config_data) #metadata is tab delimited
        self.config_data._config_data["metadata_parameters"]["genotype_column"]=correct_gt_column
        self.assertTrue(metadata_utils.load_metadata(self.config_data) )


    def test_get_metavalue(self):
        metadata_utils.load_metadata(self.config_data)
        vcf_files=[f'{self.config_data.vcf_dir}{f}' for f in listdir(self.config_data.vcf_dir) ]
        metadata_utils.samples_in_metadata(vcf_files)
        self.assertEqual("3.1.2", metadata_utils.get_metavalue("10071_3#75",self.config_data.genotype_column))
        with self.assertWarns(Warning):
            self.assertEqual(["DummyID"], metadata_utils.samples_in_metadata(["DummyID"]))
        self.assertWarns(UserWarning, metadata_utils.samples_in_metadata, samples=["DummyID", "DummyID2"])
        self.assertRaises(ValueError, metadata_utils.get_metavalue, 
                          sample="Dummy_ID",
                          value_column=self.config_data.genotype_column)

    def test_snp_identification_init(self):
        config_data = InputConfiguration(self.config_file)
        
        for stub in config_data.name_stubs:
            name_converters.name_stubs.add(stub)

        GenotypeSnpIdentifier(config_data)


    def test_snp_identification(self):
       
        for stub in self.config_data.name_stubs:
            name_converters.name_stubs.add(stub)

        snp_identifier=GenotypeSnpIdentifier(self.config_data)
        computed_genotypes: Genotypes = snp_identifier.identify_snps()
        # This file should only be changed with a
        # major change in the way SNPs are identified. 
        # At the moment, no such change is planned.
        # with open(expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/genotype.pkl"), "wb") as output:
        #     pickle.dump(computed_genotypes, output)
        with open(f'{self.valid_data}/genotype.pkl', "rb") as pickled_file:
            valid_genotypes: Genotypes = pickle.load(pickled_file)
        for computed_genotype in computed_genotypes.genotypes:
            for valid_genotype in valid_genotypes.genotypes:
                if computed_genotype.name==valid_genotype.name:
                    self.assertTrue(sorted(computed_genotype._subgenotypes)==sorted(valid_genotype._subgenotypes) )
                    self.assertTrue(sorted(computed_genotype._alleles)==sorted(valid_genotype._alleles) )

    def test_snp_optimiser(self):
        with open(expanduser(f'{self.valid_data}/amplicon_intervals.pkl'), "rb") as pickled_file:
            valid_amplicon_intervals = pickle.load(pickled_file)

        with open(f'{self.valid_data}/genotype.pkl', "rb") as pickled_file:
            valid_genotypes: Genotypes = pickle.load(pickled_file)

        snp_opimiser=SnpOptimiser()
        max_iterval_len=1000
        amplicon_intervals=snp_opimiser.optimise(max_iterval_len,valid_genotypes, rare_gts=["4.3.1.1.P1"])
        # This file should only be changed with a major change in the way SNPs are identified. 
        # At the moment, no such change is planned.        
        # with open(expanduser(f'{self.valid_data}/amplicon_intervals.pkl'), "wb") as output:
        #     pickle.dump(amplicon_intervals, output)
        # This file is required by some downstream analysis, but should be kept fixed
        # with open(self.config_data.genotypes_data, "wb") as output:
        #     pickle.dump(valid_genotypes, output)
        self.assertTrue(amplicon_intervals==valid_amplicon_intervals)


if __name__ == '__main__':
    unittest.main(verbosity=2)