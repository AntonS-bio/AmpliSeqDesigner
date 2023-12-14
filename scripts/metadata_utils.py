import pandas as pd
from typing import List
from collections import Counter
import name_converters
from data_classes import InputConfiguration
import warnings

meta_data: pd.DataFrame = pd.DataFrame()
metadata_filename: str = ""
genotype_column: str=""
warnings.formatwarning = lambda msg, *args, **kwargs: f'{msg}\n'

#class MetadataUtilities:
def load_metadata(config: InputConfiguration) -> bool:
        #filename: str, separator: str) -> True:
    global meta_data, metadata_filename
    meta_data=pd.read_csv(config.meta_data_file, sep=config.metadata_delim, index_col=0)
    if meta_data.size==0 or len(meta_data.columns)<1:
        raise ValueError(f'Metadata file has single column. Perhaps delimiter {config.metadata_delim} is incorrect? Use \\t for tab.')
    if config.genotype_column not in meta_data.columns:
        raise ValueError(f'Could not find genotype column {config.genotype_column} in metadata file {config.meta_data_file}.')
    if Counter(meta_data.index.duplicated())[True]!=0:
        duplicate_indices=meta_data.index[meta_data.index.duplicated()]
        values_to_print="\n"+"\n".join(duplicate_indices)
        raise ValueError(f'Metadata file {config.meta_data_file} has {len(duplicate_indices)} duplicated indices in column {duplicate_indices.name}: {values_to_print}')
    metadata_filename=config.meta_data_file
    return True
    
def get_metavalue(sample: str, value_column: str) -> str:
    global meta_data
    if value_column not in meta_data.columns:
        raise ValueError(f'File {metadata_filename} does not have column {value_column}')
    if name_converters.value_exists(sample):
        return str(meta_data.loc[name_converters.get_sample(sample),value_column])
    else:
        raise ValueError(f'Sample {sample} was not found among supplied VCFs')

def samples_in_metadata(samples: List[str]) -> List[str]:
    global meta_data
    samples_without_metadata:List[str]=[]
    for sample in samples:
        if not name_converters.add_value(sample, set(meta_data.index)):
            samples_without_metadata.append(sample)
    if len(samples_without_metadata)!=0:
        missing_to_print="\n"+"\n".join( samples_without_metadata[ 0: min(5,len(samples_without_metadata)) ] )
        warnings.warn(f'Warning! {len(samples_without_metadata)} samples are missing in metadata. Few examples are: {missing_to_print}\nThese samples will be excluded.')
    return samples_without_metadata
