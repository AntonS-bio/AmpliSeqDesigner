from inputs_validation import ValidateFiles

from os import listdir
import metadata_utils as metadata_utils
from load_vcfs import VCFutilities
from hierarchy_utils import HierarchyUtilities
import pandas as pd
from typing import Dict, List
from data_classes import Genotype, Genotypes, Sample, SNP, InputConfiguration
from tqdm import tqdm

class GenotypeSnpIdentifier:

    debug=True
    # def __init__(self, vcf_dir: str, hierarchy_file: str, meta_data_file: str,
    #              genotype_column: str, senstivity: float, specificity: float, **kwargs) -> None:
    def __init__(self, config: InputConfiguration) -> None:    

        self.vcf_utils=VCFutilities()
        self.file_validator=ValidateFiles()
        self.master_vcf=pd.DataFrame()
        if self.debug:
            self.vcf_files: List[str]=[f'{config.vcf_dir}{f}' for f in listdir(config.vcf_dir) ]
        else:
            self.vcf_files: List[str]=[f'{config.vcf_dir}{f}' for f in listdir(config.vcf_dir) ] #use this to limit the number of samples loaded in debugging mode

        if config.repeats_bed_file!="":
            self.file_validator.validate_bed(config.repeats_bed_file)
            self.file_validator.contigs_in_vcf(config.repeats_bed_file,self.vcf_files[0]) ###Don't check every VCF, only the first one. Unlikely to be problematic

        self.vcf_utils.load_repeat_regions(config.repeats_bed_file)
        metadata_utils.load_metadata(config)
        samples_without_metadata=metadata_utils.samples_in_metadata(self.vcf_files)
        for sample in samples_without_metadata:
            self.vcf_files.remove(sample)
        metadata_utils.genotype_column=config.genotype_column #this will be an input

        if float(config.snp_sensitivity)<1 or float(config.snp_specificity)<1:
            raise ValueError("Either sensitivity or specifity is <1, did you enter deciman instead of integer? Ex: 0.1 instead of 10.")
        self.file_validator.validate_hierarchy(config.hierarchy_file)
        self.hierarchy_utils=HierarchyUtilities( sensitivity_limit=config.snp_sensitivity/100, specificity_limit=config.snp_specificity/100)
        self.hierarchy_utils.load_hierarchy(config.hierarchy_file)


    def identify_snps(self) -> Genotypes:
        """Scans VCF files for SNPs that segregate genotypes of interest from the rest.

        """
        vcfs: List[Sample]=[]
        print("Loading VCFs")
        all_snps: Dict[SNP, SNP]={} #this is needed for speed. List lookup is slow and sets by nature don't support indexing
        with tqdm(total=len(self.vcf_files)) as progress_meter:
            for i, vcf in enumerate(self.vcf_files):
                vcf_obj=Sample(vcf.replace(".vcf",""), vcf)
                vcf_obj.genotype=metadata_utils.get_metavalue(vcf, metadata_utils.genotype_column)
                self.vcf_utils.vcf_to_snps(vcf, all_snps, vcf_obj)
                vcfs.append( vcf_obj )
                progress_meter.update(1)

        genotype_bifurcating_snps: Genotypes=self.hierarchy_utils.find_defining_snps(vcfs)

        return genotype_bifurcating_snps
    