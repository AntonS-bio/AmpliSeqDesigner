from os.path import expanduser, realpath, dirname
from typing import List, Dict
import unittest
from sys import path
unit_test_dir = dirname(realpath(__file__))
path.append( f'{unit_test_dir}/..')
print(unit_test_dir)
from data_classes import InputConfiguration
import metadata_utils as metadata_utils
import name_converters

class TestNameConverters(unittest.TestCase):
    valid_data=expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/")
    invalid_data=expanduser("~/HandyAmpliconTool/unit_test_data/invalid_data/")
    config_file=expanduser("~/HandyAmpliconTool/unit_test_data/unittest.json")
    valid_vcf_full_name=f'{valid_data}/vcfs/8490_5#12.vcf'
    valid_vcf_filename='8490_5#12.vcf'
    valid_vcf_sample='8490_5#12'

    def setUp(self) -> None:
        return super().setUp()
    
    def _reset_data(self):
        self.config_data = InputConfiguration(self.config_file)
        name_converters.clear_all_names()
        metadata_utils.load_metadata(self.config_data)
        metadata_utils.samples_in_metadata(list(metadata_utils.meta_data.index))

    def test_add_value(self):
        self._reset_data()
        self.assertTrue( name_converters.add_value(self.valid_vcf_filename, set(metadata_utils.meta_data.index)) )

    def test_get_sample(self):
        self._reset_data()
        name_converters.add_value(self.valid_vcf_filename, set(metadata_utils.meta_data.index))
        self.assertEqual( name_converters.get_sample(self.valid_vcf_full_name), self.valid_vcf_sample )
        self.assertEqual( name_converters.get_sample(self.valid_vcf_filename), self.valid_vcf_sample )
        self.assertEqual( name_converters.get_sample(self.valid_vcf_filename), self.valid_vcf_sample )

    def test_remove_name_stubs(self):
        self._reset_data()
        name_converters.name_stubs.add(f'{self.valid_data}/vcfs/')
        self.assertEqual( name_converters.remove_name_stubs(self.valid_vcf_full_name), self.valid_vcf_filename)
        self.assertEqual( name_converters.remove_name_stubs(self.valid_vcf_filename), self.valid_vcf_filename)

    def test_value_exits(self):
        self._reset_data()
        self.assertFalse(name_converters.value_exists(self.valid_vcf_full_name))
        self.assertFalse(name_converters.value_exists(self.valid_vcf_filename))
        name_converters.add_value(self.valid_vcf_full_name, set(metadata_utils.meta_data.index))
        self.assertTrue(name_converters.value_exists(self.valid_vcf_full_name))
        self.assertTrue(name_converters.value_exists(self.valid_vcf_filename))

    def test_address_to_filename(self):
        self._reset_data()
        self.assertEqual( name_converters.address_to_filename(self.valid_vcf_full_name), self.valid_vcf_filename)

    def test_filename_to_prefix(self):
        self._reset_data()
        self.assertEqual( name_converters.filename_to_prefix(self.valid_vcf_full_name), self.valid_vcf_sample)
        self.assertEqual( name_converters.filename_to_prefix(self.valid_vcf_filename), self.valid_vcf_sample)


if __name__ == '__main__':
    unittest.main(verbosity=2)