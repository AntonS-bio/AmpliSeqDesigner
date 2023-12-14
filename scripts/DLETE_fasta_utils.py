from inputs_validation import ValidateFiles
import warnings

class FastaUtilities:

    def __init__(self):
        self.input_validator=ValidateFiles()

    def get_fasta_subseq(self, fasta_file_name, seq_contig, seq_start, seq_end): ## Add unit test for cross check with bedtools getfasta output
        self.input_validator.validate_fasta(fasta_file_name)
        if seq_start>=seq_end:
            raise ValueError(f'Sequence start {seq_start} in greater or equal than end {seq_end} for file {fasta_file_name} contig {seq_contig}')
        if seq_start<0 or seq_end<0:
            raise ValueError(f'Neither start {seq_start} nor end {seq_end} can be negative for file {fasta_file_name} contig {seq_contig}')
        self.input_validator.fasta_has_dashes(fasta_file_name)
        with open(fasta_file_name) as fasta_file:
            in_target_seq=False
            seq_position=0
            target_seq=""
            for line in fasta_file:
                if line[0]==">":
                    contig_id=line.strip().split(" ")[0][1:] #the first character is ">"
                    if contig_id==seq_contig:
                        in_target_seq=True
                    else:
                        in_target_seq=False
                if in_target_seq:
                    if seq_position+len(line.strip())>=seq_start:
                        #the target sequences either starts or continues on this line
                        if target_seq=="":
                            target_seq+=line.strip()[ seq_start-seq_position : min(seq_end-seq_start,len(line.strip())) ]
                        else:
                            target_seq+=line.strip()[ 0 : min(seq_end-seq_position,len(line.strip())) ]
                    seq_position+=len(line.strip())
                if seq_position>seq_end:
                    break # exit iterating over the file lines
            return target_seq
     
        
