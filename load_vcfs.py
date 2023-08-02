import pandas as pd
from collections import Counter
import warnings
from typing import Dict
import name_converters

## Consider replacing some of this with GATKs VariantsToTable

class VCFutilities():
    def __init__(self) -> None:
        pass
        

    def load_file(self, filename: str) -> pd.DataFrame:
        #### read the VCF file until #CHROM is reached
        #### the preceeding lines are header
        vcf_file_type="Unknown"
        header_rows=0
        with open(filename) as vcf_file:
            for i, line in enumerate(vcf_file):
                if line[0:6]=="#CHROM":
                    header_rows=i
                    header_columns=line.strip().split("\t")
                    if len(header_columns)==10:
                        vcf_file_type="single_sample"
                    elif len(header_columns)>10:
                        vcf_file_type="multi_sample"
                    else:
                        raise ValueError(f'VCF file {filename} has fewer than 10 columns which is minimum required')

                    break

        if vcf_file_type=="single_sample":
            vcf_datatypes={"CHROM":"string","POS":int, "REF": "string", "ALT": "string", "FORMAT": "string"}
            vcf_data=pd.read_csv(filename, sep="\t", index_col=[0,1], header=0, usecols=[0,1,4,5,8,9],  dtype=vcf_datatypes, skiprows=header_rows)
            #check for duplicate combinations of contig+position i.e. indices
            #the tool is meant for bacteria, so non-haploid calls are a problem
            self._has_duplicate_positions(vcf_data, filename)
            #drop invariant sites, non-haploid and delecitons are not supported.
            self._drop_invariant_sites(vcf_data, filename)
            vcf_data.drop(vcf_data.columns[0:-1], axis=1, inplace=True)
            name_converters.add_address(filename)
            name_converters.add_value(vcf_data.columns[-1],  name_converters.get_sample(filename))
        else:
            raise ValueError(f'The vcf type {vcf_file_type} is not currently supported')
        return vcf_data

    def _drop_invariant_sites(self, data: pd.DataFrame, filename:str) -> None:
        #shrinks data to only show non-reference alleles
        if len(data["FORMAT"].unique())==1:
            value_columns=data["FORMAT"].unique()[0].split(":")
            if "GT" in value_columns:
                genotype_col=value_columns.index("GT")
            else:
                raise ValueError(f'VCF file {filename} does not have genotype code [GT] in column FORMAT')
        else:
            raise ValueError(f'VCF file {filename} has multiple different combinations of values in FORMAT column. Not currently supported')
        
        ##### !!!!! This needs further work. The same position may have different bases in different files i.e. three or four alleles at same position
        sample_column=data.columns[-1]
        data[sample_column]=[f.split(":")[genotype_col] for f in data[sample_column]]
        #Genotype is either a base or "." for no call
        data.drop(  data.index[data[sample_column]=="0"]  , axis="index", inplace=True)
        no_calls_mask=data[sample_column]=="."
        data.loc[~no_calls_mask,sample_column]=data.loc[~no_calls_mask,"ALT"] #this only works for single sample VCF of haploids
        ##### !!!!! This needs further work. 

    
    def _has_duplicate_positions(self, data: pd.DataFrame, filename: str) -> bool:
        if Counter(data.index.duplicated())[True]!=0:
            duplicate_indices=data.index[data.index.duplicated()]
            values_to_print="\n"+"\n".join(duplicate_indices[0:min(5, len(duplicate_indices))])
            warnings.warn(f'VCF file {filename} has {len(duplicate_indices)} duplicated positions. These will be dropped as tools is meant for bacteria. Few examples are: {values_to_print}')
            data.drop(duplicate_indices, axis="index", inplace=True)
            return True
        else:
            return False
    
    def merge_vcfs(self, master_vcf: pd.DataFrame, vcf_to_add: pd.DataFrame) -> None:
        master_vcf.join(vcf_to_add, how="outer")


