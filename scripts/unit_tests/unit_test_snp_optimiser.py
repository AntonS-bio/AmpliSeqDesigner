from os.path import expanduser, realpath, dirname
import unittest
from typing import List
from sys import path
unit_test_dir = dirname(realpath(__file__))
path.append( f'{unit_test_dir}/..')
print(unit_test_dir)
from data_classes import InputConfiguration
import metadata_utils as metadata_utils

class TestSnpOptimiser(unittest.TestCase):
    valid_data=expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/")
    invalid_data=expanduser("~/HandyAmpliconTool/unit_test_data/invalid_data/")
    temp_dir=expanduser("~/HandyAmpliconTool/unit_test_data/temp_data/")
    config_file=expanduser("~/HandyAmpliconTool/unit_test_data/unittest.json")
    ref_fasta=f'{valid_data}/GCF_000195995.1_short.fna'
    long_fasta=f'{valid_data}/GCF_000195995.1_for_vcf.fna'

    def setUp(self) -> None:
        self.config_data = InputConfiguration(self.config_file)
        return super().setUp()



if __name__ == '__main__':
    unittest.main(verbosity=2)