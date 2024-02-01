from os.path import expanduser, realpath, dirname
from os import listdir
import unittest
from sys import path
from typing import List, Dict
import pickle
unit_test_dir = dirname(realpath(__file__))
path.append( f'{unit_test_dir}/..')
print(unit_test_dir)
from data_classes import InputConfiguration, Genotypes, Sample, SNP
from hierarchy_utils import HierarchyUtilities
import metadata_utils as metadata_utils
from load_vcfs import VCFutilities

class TestHierarchyUtils(unittest.TestCase):
    valid_data=expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/")
    invalid_data=expanduser("~/HandyAmpliconTool/unit_test_data/invalid_data/")
    config_file=expanduser("~/HandyAmpliconTool/unit_test_data/unittest.json")
    ref_fasta=f'{valid_data}/GCF_000195995.1_short.fna'
    long_fasta=f'{valid_data}/GCF_000195995.1_for_vcf.fna'

    def setUp(self) -> None:
        self.config_data = InputConfiguration(self.config_file)
        return super().setUp()
    
    def test_load_hierarchy(self):
        hierarchy=HierarchyUtilities()
        result=hierarchy.load_hierarchy(self.config_data.hierarchy_file)
        with open(self.config_data.hierarchy_file) as input_file:
            for line in input_file:
                values=line.strip().split("\t")
                self.assertTrue(values[0] in result)
                values=[f for f in values if len(f)>0]
                self.assertTrue(len(values) == len(result[values[0]].subgenotypes) )

    def test_bifurcation(self):
        pass
        #Used for debugging bifurcating SNP identification
        # hierarchy=HierarchyUtilities()
        # hierarchy.load_hierarchy("/home/lshas17/HandyAmpliconTool/test_data/inputs/2_3_1_hierarchy.tsv")
        # vcf_dir="/home/lshas17/converted_vcfs/"
        # vcf_files: List[str]=[f'{vcf_dir}{f}' for f in listdir(vcf_dir) ]
        # metadata_utils.load_metadata(self.config_data)
        # metadata_utils.samples_in_metadata(vcf_files)
        # with open(self.valid_data+"vcfs.pkl", "rb") as pickled_file:
        #     vcfs: List[Sample] = pickle.load(pickled_file)
        # genotype_bifurcating_snps: Genotypes=hierarchy.find_defining_snps(vcfs)

    def test_find_defining_snps(self):
        hierarchy=HierarchyUtilities()
        hierarchy.load_hierarchy(self.config_data.hierarchy_file)
        vcf_utils=VCFutilities()
        vcf_files: List[str]=[f'{self.config_data.vcf_dir}{f}' for f in listdir(self.config_data.vcf_dir) ]
        metadata_utils.load_metadata(self.config_data)
        metadata_utils.samples_in_metadata(vcf_files)
        vcfs: List[Sample]=[]
        print("Loading VCFs")
        all_snps: Dict[SNP, SNP]={} #this is needed for speed. List lookup is slow and sets by nature don't support indexing
        for vcf in vcf_files:
            vcf_obj=Sample(vcf.replace(".vcf",""), vcf)
            vcf_obj.genotype=metadata_utils.get_metavalue(vcf.split("/")[-1].replace(".vcf",""), self.config_data.genotype_column)
            vcf_utils.vcf_to_snps(vcf, all_snps, vcf_obj)
            vcfs.append( vcf_obj )
        genotype_bifurcating_snps: Genotypes=hierarchy.find_defining_snps(vcfs)
        
        expected_snp_counts={"0":2, "4":25, "2.3.2": 3, "2.3.1": 104, "3.1.1": 270, "4.1":26, "4.1.1": 114, "4.3.1.1.P1": 34, "4.3.1.3": 3, "4.3.1.3.Bdq": 33}
        for genotype in genotype_bifurcating_snps.genotypes:
            self.assertTrue(genotype.name in expected_snp_counts)
            self.assertTrue( len(genotype.defining_snp_coordinates)==expected_snp_counts[genotype.name] )
            expected_snp_counts.pop(genotype.name)
        self.assertTrue(len(expected_snp_counts)==0)

if __name__ == '__main__':
    unittest.main(verbosity=2)