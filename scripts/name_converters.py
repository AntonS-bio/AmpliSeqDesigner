from typing import Dict, Set
from os.path import split

_value_to_sample: Dict[str,str]={}
converters=[_value_to_sample]
master_sample_set: Set[str]=set()
name_stubs: Set[str]=set()

def get_sample(input_value) -> str:
    for converter in converters:
        if input_value in converter:
            return converter[input_value]
    raise ValueError(f'Value {input_value} is unknown among IDs')

def remove_name_stubs(raw_name: str) -> str:
    cleaned_name=raw_name
    for stub in name_stubs:
        cleaned_name=cleaned_name.replace(stub,"")
    return cleaned_name

def value_exists(input_value) -> bool:
    for converter in converters:
        if input_value in converter:
            return True
    return False

def add_address(address) -> None:
    probable_sample_name=filename_to_prefix(address_to_filename((address)))
    _value_to_sample[address]=probable_sample_name
    _value_to_sample[address_to_filename(address)]=probable_sample_name

def add_value(value_to_add: str, valid_values: Set[str]) -> bool:
    values_to_add=sample_name(value_to_add)
    if values_to_add.full_address in valid_values:
        _value_to_sample[value_to_add]=values_to_add.full_address
    if values_to_add.file_name in valid_values:
        _value_to_sample[value_to_add]=values_to_add.file_name
        _value_to_sample[values_to_add.file_name]=values_to_add.file_name
    if values_to_add.file_prefix in valid_values:
        _value_to_sample[value_to_add]=values_to_add.file_prefix
        _value_to_sample[values_to_add.file_name]=values_to_add.file_prefix
        _value_to_sample[values_to_add.file_prefix]=values_to_add.file_prefix
    if values_to_add.file_prefix_wo_stub in valid_values:
        _value_to_sample[value_to_add]=values_to_add.file_prefix_wo_stub
        _value_to_sample[values_to_add.file_name]=values_to_add.file_prefix_wo_stub
        _value_to_sample[values_to_add.file_prefix]=values_to_add.file_prefix_wo_stub
        _value_to_sample[values_to_add.file_prefix_wo_stub]=values_to_add.file_prefix_wo_stub

    return value_to_add in _value_to_sample

def address_to_filename(address) -> str:
    directory, filename=split(address)
    return filename
    
def filename_to_prefix(filename) -> str:
    sample=".".join(filename.split(".")[0:-1])
    return sample

class sample_name:
    def __init__(self, value: str) -> None:
        self._full_address=value
        self._file_name=address_to_filename(value)
        self._file_prefix=filename_to_prefix(self.file_name)
        self._file_prefix_wo_stub=remove_name_stubs(self.file_prefix)

    @property
    def full_address(self):
        return self._full_address

    @property
    def file_name(self):
        return self._file_name

    @property
    def file_prefix(self):
        return self._file_prefix

    @property
    def file_prefix_wo_stub(self):
        return self._file_prefix_wo_stub