from typing import Dict, Set
from os.path import split

_value_to_sample: Dict[str,str]={}
converters=[_value_to_sample]
master_sample_set: Set[str]=set()


def get_sample(input_value) -> str:
    for converter in converters:
        if input_value in converter:
            return converter[input_value]
    raise ValueError(f'Value {input_value} is unknown among IDs')
    
def value_exists(input_value) -> bool:
    for converter in converters:
        if input_value in converter:
            return True
    return False

def add_address(address) -> None:
    probable_sample_name=filename_to_prefix(address_to_filename((address)))
    _value_to_sample[address]=probable_sample_name
    _value_to_sample[address_to_filename(address)]=probable_sample_name

def add_value(key, value) -> None:
    _value_to_sample[key]=value
    _value_to_sample[value]=value #sample to sample itself.

def address_to_filename(address) -> str:
    directory, filename=split(address)
    return filename
    
def filename_to_prefix(filename) -> str:
    sample=".".join(filename.split(".")[0:-1])
    return sample
