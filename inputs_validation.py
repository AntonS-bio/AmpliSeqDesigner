### Checks all the inputs are valid before executing the program
### this avoids program crashing mid-way due to bad inputs.

from os.path import exists
import warnings

class ValidateFiles:
    def validate_bed(self, bed_file_name):
        if not exists(bed_file_name):
            raise FileExistsError(f'File or directory {bed_file_name} does not exist')
        
        with open(bed_file_name) as bed_file:
            for line_counter, line in enumerate(bed_file):
                if line.find("\t")==-1:
                    raise ValueError(f'No tab-delimited found in file {bed_file_name} on line {line_counter}:\n {line}')
                line_values=line.split("\t")
                if len(line_values)<3:
                    raise ValueError(f'Min numpber of tab-delimited columns in 3, but only {len(line_values)} were found on line {line_counter}:\n {line} ')
                if not line_values[1].isdigit() or not line_values[2].isdigit():
                    raise ValueError(f'Expecting numeric values in column 1 and 2, but none numeric values (posibly decimal) values were found in {bed_file_name} on line {line_counter}:\n {line}')
                if int(line_values[2])<=int(line_values[1]):
                    raise ValueError(f'Value in column 2 must be greater than value in column 1. Check {bed_file_name} on line {line_counter}:\n {line}')
            if 'line_counter' not in vars():
                raise ValueError(f'{bed_file_name} is empty')

    def validate_fasta(self, fasta_file_name):
        if not exists(fasta_file_name):
            raise FileExistsError(f'File or directory {fasta_file_name} does not exist')
        with open(fasta_file_name) as fasta_file:
            first_fifty_char=fasta_file.readline()[0:50]
            if len(first_fifty_char)==0:
                raise ValueError(f'Fasta file {fasta_file_name} is empty')
            if first_fifty_char[0]!=">":
                raise ValueError(f'Fasta file must have ">" on first line in {fasta_file_name}\n {first_fifty_char}')
            
    def check_fasta_for_dashes(self, fasta_file_name):
        with open(fasta_file_name) as fasta_file:
            for line in fasta_file:
                if line[0]!=">":
                    if line.find("-")>-1:
                        warnings.warn(f'Fasta file {fasta_file_name} has "-". This is likely to cause problems')

    def check_contigs_in_fasta(self, bed_file_name, fasta_file_name):
        ### Assume that bed and fasta files have already been validated
        bed_contigs=set()
        with open(bed_file_name) as bed_file:
            for line in (bed_file):
                bed_contigs.add(line.split("\t")[0])

        with open(fasta_file_name) as fasta_file:        
            for line in fasta_file:
                if line[0]==">":
                    contig_id=line.split(" ")[0][1:]
                    if contig_id in bed_contigs:
                        bed_contigs.remove(contig_id)
                if len(bed_contigs)==0:
                    break
            if len(bed_contigs)!=0:
                missing_contigs="\n".join( list(bed_contigs)[0:min(len(bed_contigs),10)] )
                raise ValueError(f'Some bed file {bed_file_name} contig IDs not found in fasta file {fasta_file_name} (first 10):\n {missing_contigs} ')
            








        