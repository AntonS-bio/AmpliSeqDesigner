import pandas as pd
from collections import Counter
import warnings
from typing import Dict, List, Tuple, Set
import name_converters
from data_classes import SNP, Sample
import metadata_utils
## Consider replacing some of this with GATKs VariantsToTable

   
class VCFutilities():

    def __init__(self) -> None:
        self._repeat_coordinates=set()
        pass

    def vcf_to_snps(self, filename: str, existing_snps:Dict[SNP, SNP], sample: Sample):
        #### read the VCF file until #CHROM is reached
        #### the preceeding lines are header
        vcf_file_type="Unknown"
        with open(filename) as vcf_file:
            for i, line in enumerate(vcf_file):
                if line[0:6]=="#CHROM":
                    header_columns=line.strip().split("\t")
                    if len(header_columns)==10:
                        vcf_file_type="single_sample"
                    elif len(header_columns)>10:
                        vcf_file_type="multi_sample"
                    else:
                        raise ValueError(f'VCF file {filename} has fewer than 10 columns which is minimum required')

                    break

        duplicate_positions=[]
        if vcf_file_type=="single_sample":
            vcf_datatypes={"CHROM":"string","POS":int, "REF": "string", "ALT": "string", "FORMAT": "string"}
            with open(filename) as vcf_file_handle:
                for line in vcf_file_handle:
                    if line[0]=="#" or line=="\n":
                        continue
                    chrom, pos, id, ref, alt, qual, filter, info, format, values=line.strip().split("\t")
                    pos=int(pos)-1  #the rest of the code is 0 indexed like BED and BAM, but VCF coordinates are 1-indexed
                    if (chrom, pos) in self.repeat_coordinates:
                        continue
                    if alt.find(",")>-1:
                        #multiploid line, skip with a warning
                        duplicate_positions.append((chrom, pos))
                        continue
                    gt_index=format.split(":").index("GT")
                    if gt_index==-1:
                        raise ValueError(f'VCF file {filename} does not have genotype code [GT] in SAMPLE colum at {chrom} {str(pos-1)}')

                    allele = ref if values.split(",")[gt_index].split(":")[0]=="0" else alt
                    snp=SNP(ref_contig_id=chrom, ref_base=ref, alt_base=allele, position=pos)
                    if snp not in existing_snps:
                        existing_snps[snp]=snp
                        sample.snps.append(snp)
                    else:
                        sample.snps.append( existing_snps[snp] )


            if len(duplicate_positions)>0:
                warnings.warn(f'VCF file {filename} has {len(duplicate_positions)} duplicated positions')

        else:
            raise ValueError(f'The vcf type {vcf_file_type} is not currently supported')

    # def load_file(self, filename: str) -> pd.DataFrame:
    #     #### read the VCF file until #CHROM is reached
    #     #### the preceeding lines are header
    #     vcf_file_type="Unknown"
    #     header_rows=0
    #     with open(filename) as vcf_file:
    #         for i, line in enumerate(vcf_file):
    #             if line[0:6]=="#CHROM":
    #                 header_rows=i
    #                 header_columns=line.strip().split("\t")
    #                 if len(header_columns)==10:
    #                     vcf_file_type="single_sample"
    #                 elif len(header_columns)>10:
    #                     vcf_file_type="multi_sample"
    #                 else:
    #                     raise ValueError(f'VCF file {filename} has fewer than 10 columns which is minimum required')

    #                 break

    #     if vcf_file_type=="single_sample":
    #         vcf_datatypes={"CHROM":"string","POS":int, "REF": "string", "ALT": "string", "FORMAT": "string"}
    #         vcf_data=pd.read_csv(filename, sep="\t", header=0, usecols=[0,1,4,5,8,9],  dtype=vcf_datatypes, skiprows=header_rows)
    #         vcf_data["POS"]=vcf_data["POS"]-1 #the rest of the code is 0 indexed like BED and BAM, but VCF coordinates are 1-indexed
    #         new_index=pd.MultiIndex.from_frame(vcf_data[ ["#CHROM","POS"] ])
    #         vcf_data.set_index(new_index, inplace=True)
    #         #check for duplicate combinations of contig+position i.e. indices
    #         #the tool is meant for bacteria, so non-haploid calls are a problem
    #         self._has_duplicate_positions(vcf_data, filename)
    #         #drop invariant sites, non-haploid and delecitons are not supported.
    #         self._drop_invariant_sites(vcf_data, filename)
    #         vcf_data.drop(vcf_data.columns[0:-1], axis=1, inplace=True)
    #         name_converters.add_value(vcf_data.columns[-1],metadata_utils.meta_data.index)
    #     else:
    #         raise ValueError(f'The vcf type {vcf_file_type} is not currently supported')
    #     return vcf_data

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

    # def remove_repeat_regions(self, vcf_data, bed_file) -> None:
    #     regions_to_exclude: List[ Tuple[str,int] ]=[]
    #     with open(bed_file) as bed_data:
    #         for line in bed_data:
    #             values=line.strip().split("\t")
    #             regions_to_exclude=regions_to_exclude+[(values[0],f) for f in range(int(values[1]), int(values[2]))]
    #     vcf_data.drop( [f for f in regions_to_exclude if f in vcf_data.index] , axis="index", inplace=True )

    def load_repeat_regions(self, bed_file: str):
        self._repeat_coordinates: Set[ Tuple[str, int] ] =set()
        if bed_file=="":
            return None
        with open(bed_file) as bed_data:
            for line in bed_data:
                values=line.strip().split("\t")
                self._repeat_coordinates.update( [(values[0],f) for f in range(int(values[1]), int(values[2]))] )

    
    @property
    def repeat_coordinates(self) -> Set[ Tuple[str, int] ]:
        return self._repeat_coordinates

    def _has_duplicate_positions(self, data: pd.DataFrame, filename: str) -> bool:
        if Counter(data.index.duplicated())[True]!=0:
            duplicate_indices=data.index[data.index.duplicated()]
            warnings.warn(f'VCF file {filename} has {len(duplicate_indices)} duplicated positions')
            values_to_print=""
            for index in duplicate_indices[0:min(5, len(duplicate_indices))]:
                values_to_print=f'{values_to_print}\n{index}' #the index is tuple (chr, pos), for clarity generating string to print is done in loop
            warnings.warn(f'VCF file {filename} has {len(duplicate_indices)} duplicated positions. These will be dropped as tools is meant for bacteria. Few examples are: {values_to_print}')
            data.drop(duplicate_indices, axis="index", inplace=True)
            return True
        else:
            return False

    master_vcf_temp: pd.DataFrame = pd.DataFrame()
    def merge_vcfs(self, additional_vcf: pd.DataFrame) -> None:
        self.master_vcf_temp.loc[additional_vcf.index, additional_vcf.columns[-1]]=additional_vcf[additional_vcf.columns[-1]]
