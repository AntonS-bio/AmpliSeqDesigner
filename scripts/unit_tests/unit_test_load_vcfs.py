from os.path import expanduser, realpath, dirname
from typing import List, Dict
import unittest
from sys import path
unit_test_dir = dirname(realpath(__file__))
path.append( f'{unit_test_dir}/..')
print(unit_test_dir)
from data_classes import InputConfiguration, Sample, SNP
from load_vcfs import VCFutilities
import metadata_utils as metadata_utils

class TestLoadVCF(unittest.TestCase):
    valid_data=expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/")
    invalid_data=expanduser("~/HandyAmpliconTool/unit_test_data/invalid_data/")
    config_file=expanduser("~/HandyAmpliconTool/unit_test_data/unittest.json")
    valid_vcf=f'{valid_data}/vcfs/8490_5#12.vcf'
    multisample_vcf=f'{valid_data}/multi_sample.vcf'
    invalid_vcf=f'{valid_data}/gts.tsv'

    def setUp(self) -> None:
        self.config_data = InputConfiguration(self.config_file)
        metadata_utils.load_metadata(self.config_data)
        return super().setUp()

    def test_determine_vcf_type(self):
        vcf_loader=VCFutilities()
        self.assertEqual(vcf_loader.determine_vcf_type(self.valid_vcf),"single_sample")
        self.assertEqual(vcf_loader.determine_vcf_type(self.multisample_vcf),"multi_sample")
        self.assertRaises(ValueError,vcf_loader.determine_vcf_type,filename=self.invalid_vcf)

    def test_vcf_to_snps(self):
        vcf_loader=VCFutilities()
        existing_snps: Dict[SNP, SNP]={}
        sample=Sample(self.valid_vcf.replace(".vcf",""), self.valid_vcf)
        metadata_utils.samples_in_metadata([self.valid_vcf])
        #vcf_obj.genotype=metadata_utils.get_metavalue(self.valid_vcf, metadata_utils.genotype_column)
        self.assertRaises(ValueError, vcf_loader.vcf_to_snps, 
                          filename=self.multisample_vcf,
                          existing_snps=existing_snps,
                          sample=sample)
        vcf_loader.vcf_to_snps(self.valid_vcf, existing_snps, sample)
        self.assertEqual(len(sample.snps),681)
        self.assertEqual(len(existing_snps),681)
        vcf_loader.vcf_to_snps(self.valid_vcf, existing_snps, sample) #this should make no changes to either sample or existing_snps 
        self.assertEqual(len(existing_snps),681)
        self.assertEqual(len(sample.snps),681)
        
    def test_load_repeat_regions(self):
        vcf_loader=VCFutilities()
        vcf_loader.load_repeat_regions(self.config_data.repeats_bed_file)
        self.assertEqual(len(vcf_loader.repeat_coordinates), 346834)


if __name__ == '__main__':
    unittest.main(verbosity=2)