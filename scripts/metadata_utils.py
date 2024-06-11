import pandas as pd
from typing import List
from collections import Counter
import name_converters
from data_classes import InputConfiguration
from os.path import exists, join
import warnings

meta_data: pd.DataFrame = pd.DataFrame()
metadata_filename: str = ""
genotype_column: str=""
output_dir: str=""
warnings.formatwarning = lambda msg, *args, **kwargs: f'{msg}\n'


def load_metadata(config: InputConfiguration) -> bool:
    """Processes metadata file with genotypes information 
    :param config: config object with address of the metadata file
    :type config: InputConfiguration
    """
    global meta_data, output_dir, metadata_filename, genotype_column
    metadata_filename=config.meta_data_file
    genotype_column=config.genotype_column
    output_dir = config.output_dir

    if not exists(metadata_filename):
        warnings.warn(f'Metadata file {metadata_filename} does not exist. Please check the spelling.')
        return False

    meta_data=pd.read_csv(config.meta_data_file, sep=config.metadata_delim, index_col=0)
    if meta_data.size==0 or len(meta_data.columns)<1:
        warnings.warn(f'Metadata file has single column. Perhaps delimiter {config.metadata_delim} is incorrect? Use \\t for tab.')
        return False
    if config.genotype_column not in meta_data.columns:
        warnings.warn(f'Could not find genotype column {config.genotype_column} in metadata file {config.meta_data_file}. Did you specify correct delimiter?')
        return False
    if Counter(meta_data.index.duplicated())[True]!=0:
        duplicate_indices=meta_data.index[meta_data.index.duplicated()]
        values_to_print="\n"+"\n".join(duplicate_indices)
        warnings.warn(f'Metadata file {config.meta_data_file} has {len(duplicate_indices)} duplicated indices in column {duplicate_indices.name}: {values_to_print}')
        return False
    
    return True
    
def get_metavalue(sample: str, value_column: str) -> str:
    """Get the metadata value for provided sample 
    :param sample: sample name, name_converter is used to identify homogenous name
    :type sample: str
    :param value_column: metadata file column name from which to return value
    :type value_column: str
    :return: metadata value 
    :rtype: str
    """    
    global meta_data
    if value_column not in meta_data.columns:
        raise ValueError(f'File {metadata_filename} does not have column {value_column}')
    if name_converters.value_exists(sample):
        return str(meta_data.loc[name_converters.get_sample(sample),value_column])
    else:
        raise ValueError(f'Sample {sample} was not found among supplied VCFs')

def samples_in_metadata(samples: List[str]) -> List[str]:
    """Checks if samples are present in metadata file, is used to identify homogenous name
    :param samples:list of samples names to check against loaded metadata file
    :type samples: List[str]
    :return: List of names that were not found in metadata
    :rtype: List[str]
    """   
    global meta_data, output_dir
    samples_without_metadata:List[str]=[]
    for sample in samples:
        if not name_converters.add_value(sample, set(meta_data.index)):
            samples_without_metadata.append(sample)
    if len(samples_without_metadata)!=0:
        with open(join(output_dir,"samples_wo_metadata.txt"), "w") as output:
            for sample in samples_without_metadata:
                output.write(sample+"\n")
        #missing_to_print="\n"+"\n".join( samples_without_metadata[ 0: min(5,len(samples_without_metadata)) ] )
        warnings.warn(f'Warning! {len(samples_without_metadata)} samples are missing in metadata. Full list in {output_dir}/samples_wo_metadata.txt')
    return samples_without_metadata
