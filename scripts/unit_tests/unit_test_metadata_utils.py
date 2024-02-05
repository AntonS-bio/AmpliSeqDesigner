from os.path import expanduser, realpath, dirname
import unittest
from sys import path
unit_test_dir = dirname(realpath(__file__))
path.append( f'{unit_test_dir}/..')
print(unit_test_dir)
from data_classes import InputConfiguration
import metadata_utils as metadata_utils

class TestMetadataUtils(unittest.TestCase):
    valid_data=expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/")
    invalid_data=expanduser("~/HandyAmpliconTool/unit_test_data/invalid_data/")
    config_file=expanduser("~/HandyAmpliconTool/unit_test_data/unittest.json")
    valid_vcf=f'{valid_data}/vcfs/8490_5#12.vcf'
    multisample_vcf=f'{valid_data}/multi_sample.vcf'
    invalid_vcf=f'{valid_data}/gts.tsv'

    def setUp(self) -> None:
        return super().setUp()

    def _reset_data(self):
        self.config_data = InputConfiguration(self.config_file)
        metadata_utils.load_metadata(self.config_data)
        metadata_utils.samples_in_metadata(list(metadata_utils.meta_data.index))


    def test_load_metadata(self):
        self._reset_data()
        self.assertEqual(metadata_utils.meta_data.shape,(109,1))
        self.config_data._config_data["metadata_parameters"]["genotype_column"]="NoneSuch"
        self.assertRaises(ValueError,metadata_utils.load_metadata, config=self.config_data)
        self._reset_data()
        self.config_data._config_data["input_files"]["meta_data_file"]=f'{self.invalid_data}/gts.tsv'
        self.assertRaises(ValueError,metadata_utils.load_metadata, config=self.config_data)

    def test_get_metavalue(self):
        self._reset_data()
        self.assertEqual(metadata_utils.get_metavalue(sample="10060_5#14", value_column=self.config_data.genotype_column), "2.1.9"    )
        self.assertRaises(ValueError, metadata_utils.get_metavalue, 
                          sample="NoneSuch",value_column="Genotype")
        self.assertRaises(ValueError, metadata_utils.get_metavalue, 
                          sample="10060_5#14",value_column="NoneSuch")
        
    def test_samples_in_metadata(self):
        self._reset_data()
        self.assertWarns(Warning,metadata_utils.samples_in_metadata, samples=list(["NoneSuch"]))
        self.assertTrue( "NoneSuch" in  metadata_utils.samples_in_metadata(["NoneSuch"])   )

    def samples_in_metadata(self):
        pass


if __name__ == '__main__':
    unittest.main(verbosity=2)