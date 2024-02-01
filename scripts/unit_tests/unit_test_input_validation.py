from os.path import realpath, dirname, expanduser
unit_test_dir = dirname(realpath(__file__))
from sys import path
path.append( f'{unit_test_dir}/..')
import unittest
import inputs_validation
from data_classes import InputConfiguration


class TestInputValidation(unittest.TestCase):
    valid_data=expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/")
    invalid_data=expanduser("~/HandyAmpliconTool/unit_test_data/invalid_data/")
    config_file=expanduser("~/HandyAmpliconTool/unit_test_data/unittest.json")

    def test_validate_hierarchy(self):
        config=InputConfiguration(self.config_file)
        validator=inputs_validation.ValidateFiles()
        self.assertRaises(FileExistsError,validator.validate_hierarchy, hierarchy_file_name="empty", config_data=config)
        self.assertTrue(validator.validate_hierarchy(f'{self.valid_data}/genotype_hierarcy.tsv', config) )
        self.assertRaises(ValueError,validator.validate_hierarchy, hierarchy_file_name=f'{self.invalid_data}/genotype_hierarcy.tsv', config_data=config)

    def test_validate_bed(self):
        validator=inputs_validation.ValidateFiles()
        self.assertRaises(FileExistsError,validator.validate_bed, bed_file_name="empty")
        self.assertTrue(validator.validate_bed(f'{self.valid_data}/test_data.bed'))
        self.assertRaises(ValueError,validator.validate_bed, bed_file_name=f'{self.invalid_data}/test_data.bed')

    def test_contigs_in_fasta(self):
        validator=inputs_validation.ValidateFiles()
        self.assertTrue(validator.contigs_in_fasta(f'{self.valid_data}/test_data.bed', f'{self.valid_data}/GCF_000195995.1_short.fna' ))
        self.assertRaises(ValueError,validator.contigs_in_fasta,
                           bed_file_name=f'{self.invalid_data}/test_data.bed',
                           fasta_file_name=f'{self.valid_data}/GCF_000195995.1_short.fna')
        
    def test_contigs_in_vcf(self):
        validator=inputs_validation.ValidateFiles()
        self.assertTrue(validator.contigs_in_vcf(f'{self.valid_data}/test_data.bed', f'{self.valid_data}/ERR4451473.vcf' ))
        with self.assertWarns(Warning):
            validator.contigs_in_vcf(bed_file_name=f'{self.invalid_data}/test_vcf_data.bed',
                           vcf_file_name=f'{self.valid_data}/ERR4451473.vcf') 

    def test_fasta_has_dashes(self):
        validator=inputs_validation.ValidateFiles()
        self.assertFalse(validator.fasta_has_dashes(f'{self.valid_data}/GCF_000195995.1_short.fna', ))
        with self.assertWarns(Warning):
            validator.fasta_has_dashes(f'{self.invalid_data}/GCF_000195995.1_short.fna')


if __name__ == '__main__':
    unittest.main(verbosity=2)