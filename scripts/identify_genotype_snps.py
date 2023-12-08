from inputs_validation import ValidateFiles

from os import listdir
import metadata_utils as metadata_utils
from load_vcfs import VCFutilities
from hierarchy_utils import HierarchyUtilities
import pandas as pd
from typing import Dict, List
from data_classes import Genotype, Genotypes, Sample, SNP
from tqdm import tqdm

class GenotypeSnpIdentifier:

    debug=True
    def __init__(self, vcf_dir: str, hierarchy_file: str, meta_data_file: str,
                 genotype_column: str, senstivity: float, specificity: float, **kwargs) -> None:
        """Initialises class instance

        :param vcf_dir: The directory containing VCF files that will be used to determine genotype defining SNPs, defaults to ""
        :type vcf_dir: str, optional

        :param hierarchy_file: Tab delimited file specifying heirarchy of genotypes. Ex. "4 4.1 4.2 4.1.1" specifies 
        that genotype 4 has three subgenotypes:4.1, 4.2 and 4.1.1
        :type hierarchy_file: str
                
        :param meta_data_file: A file containing metadata for VCF samples. The genotypes of each sample in VCF 
        will be determined based on a column from this file
        :type meta_data_file: str

        :param genotype_column: Name of column in meta_data_file that contains the genotype labels
        :type meta_data_file: str

        :param repeat_regions_file: Bed file containing the regions (usually repeats) which will be excluded from analysis
        :type repeat_regions_file: str, optional
        
        :param metadata_delimiter: Delimiter in the metadata file, defaults to ","
        :type metadata_delimiter: str, optional

        :param senstivity: Cutoff SNP sensitivity of SNP vs target genotype ex. 99, 90. Intended to accomodate noisy data.
        :type metadata_delimiter: float

        :param specificity: Cutoff SNP specificity of SNP vs target genotype ex. 99, 90. Intended to accomodate noisy data.
        :type metadata_delimiter: float
        """
        self.vcf_utils=VCFutilities()
        self.file_validator=ValidateFiles()
        self.master_vcf=pd.DataFrame()
        if self.debug:
            self.vcf_files: List[str]=[f'{vcf_dir}{f}' for f in listdir(vcf_dir) ]
        else:
            self.vcf_files: List[str]=[f'{vcf_dir}{f}' for f in listdir(vcf_dir) ]

        if kwargs.get("repeat_regions_file","")!="":
            self.repeat_regions_file=kwargs.get("repeat_regions_file","")
            self.file_validator.validate_bed(self.repeat_regions_file)
            self.file_validator.contigs_in_vcf(self.repeat_regions_file,self.vcf_files[0]) ###Don't check every VCF, only the first one. Unlikely to be problematic

        self.vcf_utils.load_repeat_regions(self.repeat_regions_file)
        metadata_utils.load_metadata(meta_data_file, kwargs.get("metadata_delimiter",",") )
        samples_without_metadata=metadata_utils.samples_in_metadata(self.vcf_files)
        for sample in samples_without_metadata:
            self.vcf_files.remove(sample)
        metadata_utils.genotype_column=genotype_column #this will be an input

        if float(senstivity)<1 or float(specificity)<1:
            raise ValueError("Either sensitivity or specifity is <1, did you enter deciman instead of integer? Ex: 0.1 instead of 10.")
        self.file_validator.validate_hierarchy(hierarchy_file)
        self.hierarchy_utils=HierarchyUtilities( sensitivity_limit=senstivity/100, specificity_limit=specificity/100)
        self.hierarchy_utils.load_hierarchy(hierarchy_file)


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
    