from inputs_validation import ValidateFiles
import time 
from multiprocessing import cpu_count, Pool

#### START Identification of genotype defining SNPS #### 
from os import listdir
import metadata_utils as metadata_utils
from load_vcfs import VCFutilities
from hierarchy_utils import HierarchyUtilities
import pandas as pd
from typing import Dict, List
from data_classes import Genotype, Genotypes
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
            self.vcf_files: List[str]=[f'{vcf_dir}{f}' for f in listdir(vcf_dir) ][0:2000]
        else:
            self.vcf_files: List[str]=[f'{vcf_dir}{f}' for f in listdir(vcf_dir) ]

        if kwargs.get("repeat_regions_file","")!="":
            self.repeat_regions_file=kwargs.get("repeat_regions_file","")
            self.file_validator.validate_bed(self.repeat_regions_file)
            self.file_validator.contigs_in_vcf(self.repeat_regions_file,self.vcf_files[0]) ###!!!! Why 0?

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
        vcfs: List[pd.DataFrame]=[]
        print("Loading VCFs")
        with tqdm(total=len(self.vcf_files)) as progress_meter:
            for i, vcf in enumerate(self.vcf_files):
                vcfs.append( self.vcf_utils.load_file(vcf ) )
                progress_meter.update(1)

        ##get total indices from all vcfs to preallocate dataframe, this is much faster than merge and the loading of the vcfs can be parallelised
        vcf_columns=[""]*len(vcfs)
        with tqdm(total=len(vcfs)) as progress_meter:
            for i, vcf_data in enumerate(vcfs):
                vcf_columns[i]=vcf_data.columns[-1]
                progress_meter.update(1)
        all_vcf_indices=sorted(set([k for f in vcfs for k in f.index]))
        master_vcf=pd.DataFrame(index=all_vcf_indices, columns=vcf_columns)

        with tqdm(total=len(vcfs)) as progress_meter:
            for vcf_data in vcfs:
                master_vcf.loc[vcf_data.index, vcf_data.columns[-1]]=vcf_data[vcf_data.columns[-1]]
                progress_meter.update(1)

        del vcfs

        master_vcf.fillna("REF", inplace=True)
        #sys.exit()

        if self.repeat_regions_file!="":
            self.vcf_utils.remove_repeat_regions(master_vcf,self.repeat_regions_file)

        genotype_bifurcating_snps: Genotypes=self.hierarchy_utils.find_defining_snps(master_vcf)

        #### !!!! For testing only
        print('dense : {:0.0f} bytes'.format(master_vcf.memory_usage().sum() / 1e3) )
        sdf = master_vcf.astype(pd.SparseDtype("str", "REF"))
        print('sparse: {:0.0f} bytes'.format(sdf.memory_usage().sum() / 1e3) )
        #master_vcf=master_vcf.astype(pd.SparseDtype("str", "REF"))
        #print(Counter([meta_data.get_metavalue(f,"Final_genotype") for f in samples]))

        return genotype_bifurcating_snps