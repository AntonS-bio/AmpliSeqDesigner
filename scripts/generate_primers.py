import primer3
from typing import (
    Any,
    Dict
)

class PrimerGenerator:

    def __init__(self):
        self.global_args_generic={
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_NUM_RETURN' : 500,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'P3_FILE_FLAG' : 1
        }

    #PRIMER_TASK is a global_arg for the primer3
    # the options are:
    # generic - generic, designs primer pairs, no ordering of pairs
    # check_primers - checks primers, whatever that means. Sequence template is not essential
    # pick_primers_list - picks all primers in sequence template, but limited by SEQUENCE_INCLUDED_REGION, SEQUENCE_EXCLUDED_REGION, SEQUENCE_PRIMER_PAIR_OK_REGION_LIST
    # returns primers sorted by quality starting with best primers. Can be run on left and right separately yielding a list of unpaired primers
    # pick_sequencing_primers - selects the primers optimised for sequencing (as opposed to cloning)



    def generate_primers(self, fasta_seq: str, target_region: tuple[int, int]) -> Dict[str, Any]:
        #copy over global arguments, but some will be replaced
        global_args_to_use=dict( [ (f[0], f[1]) for f in self.global_args_generic.items()  ] )

        # a list of lists where [x,y] are min and max desired size of product. Primer3 favours shorter products
        # set to be no less than target region and up to the lenght of template
        global_args_to_use['PRIMER_PRODUCT_SIZE_RANGE']=[[target_region[1]-target_region[0], len(fasta_seq)]] 

        seq_args={"SEQUENCE_ID":"tempid", "SEQUENCE_TEMPLATE":fasta_seq, "SEQUENCE_TARGET":[target_region[0], target_region[1]-target_region[0]]}
        # SEQUENCE_TARGET has format <start>, <length>
        primers=primer3.bindings.design_primers(seq_args, global_args_to_use)
        return primers
    

