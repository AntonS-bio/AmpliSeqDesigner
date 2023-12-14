import unittest
import inputs_validation


class TestFileValidators(unittest.TestCase):
    valid_data="/home/lshas17/HandyAmpliconTool/unit_test_data/valid_data/"
    invalid_data="/home/lshas17/HandyAmpliconTool/unit_test_data/invalid_data/"

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

    def test_validate_hierarchy(self):
        validator=inputs_validation.ValidateFiles()
        self.assertRaises(FileExistsError,validator.validate_hierarchy, hierarchy_file_name="empty")
        self.assertTrue(validator.validate_hierarchy(f'{self.valid_data}/genotype_hierarcy.tsv'))
        self.assertRaises(ValueError,validator.validate_hierarchy, hierarchy_file_name=f'{self.invalid_data}/genotype_hierarcy.tsv')

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
        self.assertFalse(validator.contigs_in_vcf(bed_file_name=f'{self.invalid_data}/test_vcf_data.bed',
                           vcf_file_name=f'{self.valid_data}/ERR4451473.vcf') )

    def test_fasta_has_dashes(self):
        validator=inputs_validation.ValidateFiles()
        self.assertFalse(validator.fasta_has_dashes(f'{self.valid_data}/GCF_000195995.1_short.fna', ))
        self.assertTrue(validator.fasta_has_dashes(f'{self.invalid_data}/GCF_000195995.1_short.fna') )


if __name__ == '__main__':
    unittest.main(verbosity=2)