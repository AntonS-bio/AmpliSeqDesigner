import subprocess
from os.path import exists, isfile
from os import mkdir
from data_classes import BlastResult, InputConfiguration
from typing import List


class BlastRunner:
    def __init__(self, **kwargs) -> None:
        """Class for running BLAST searches
        :key word_size: word size to use in BLAST search, int
        :key e_value: e-value threashold for BLAST search, float
        """
        self.word_size:int=kwargs.get("word_size",InputConfiguration.blast_word_size)
        self.e_value: float= kwargs.get("e_value",InputConfiguration.blast_evalue)
    
    def db_from_file(self, file_name:str, db_dir: str) -> bool:
        self.db_file_name = file_name
        self.db_dir=db_dir
        #create blast DB
        if exists(f'{self.db_dir}') and isfile(f'{self.db_dir}'):
            raise OSError(f'Cannot create directory {self.db_dir} because there is a file with such name')
        if not exists(self.db_file_name):
            raise ValueError(f'Error generating blast database. Source file {self.db_file_name} does not exist.')
        if not exists(f'{self.db_dir}'):
            mkdir(f'{self.db_dir}')
        outcome=subprocess.run(f'makeblastdb -in {self.db_file_name} -title temp -out {self.db_dir}/temp -dbtype nucl \
            -blastdb_version 4 1>/dev/null', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if outcome.returncode==1:
            raise OSError(f'Error generating blast database: {outcome.stderr}')
        return True


    def _seq_to_file(self, seq_header: str, sequence: str, temp_dir: str) -> str:
        if seq_header[0]==">":
            input_str=seq_header+"\n"+sequence
        else:
            input_str=">"+seq_header+"\n"+sequence
        return self._write_to_fasta_file(input_str, temp_dir)

    def _write_to_fasta_file(self, string_to_write, temp_dir):
        #create blast DB
        if exists(f'{temp_dir}') and isfile(f'{temp_dir}'):
            raise OSError(f'Cannot create directory {temp_dir} because there is a file with such name')
        if not exists(temp_dir):
            mkdir(f'{temp_dir}')
        temp_fasta_file=temp_dir+"/temp.fasta"
        with open(temp_fasta_file, "w") as temp_fasta:
            temp_fasta.write(string_to_write)
        return temp_fasta_file


    def db_from_string(self, seq_header:str, sequence:str, db_dir: str) -> bool:
        fasta_file=self._seq_to_file(seq_header, sequence, db_dir)
        self.db_from_file(fasta_file, db_dir)
        return True

    def db_from_multi_sequence_string(self, sequence:str, db_dir: str) -> bool:
        fasta_file=self._write_to_fasta_file(sequence, db_dir)
        self.db_from_file(fasta_file, db_dir)
        return True

    def run_from_file(self, query_file:str) -> List[BlastResult]:
        #quickblast: reference, query, prints query
        blast_results=subprocess.run(f'blastn -query {query_file} -task \'megablast\' \
                -max_target_seqs 1000000000 -db {self.db_dir}/temp \
                -num_threads 1 -evalue {self.e_value} -word_size {self.word_size} \
                -outfmt \"6 delim=  qseqid qstart qend sseqid sstart send pident evalue qseq\"'
                , shell=True, executable="/bin/bash", stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
        if blast_results.returncode==1:
            raise OSError(f'Error running BLAST: {blast_results.stderr}')
        raw_blast_results=blast_results.stdout.decode().strip().split("\n")
        blast_hits: List[BlastResult]=[]

        for blast_hit in raw_blast_results:
            if blast_hit=="":
                continue #catches the empty line at the end of results
            blast_hit=blast_hit.split("\t")
            new_hit=BlastResult()
            new_hit.qseqid=blast_hit[0]
            new_hit.qstart=blast_hit[1]
            new_hit.qend=blast_hit[2]
            new_hit.sseqid=blast_hit[3]
            new_hit.sstart=blast_hit[4]
            new_hit.send=blast_hit[5]
            new_hit.pident=blast_hit[6]
            new_hit.evalue=blast_hit[7]
            new_hit.qseq=blast_hit[8]
            new_hit.query_file_name=query_file
            blast_hits.append(new_hit)

        return blast_hits
    
    def run_from_string(self,  seq_header: str, sequence: str, temp_dir: str) -> List[BlastResult]:
        fasta_file=self._seq_to_file(seq_header, sequence, temp_dir)
        blast_hits: List[BlastResult]=self.run_from_file(fasta_file)
        return blast_hits
    
    def run_from_multi_sequence_string(self,  sequence: str, temp_dir: str) -> List[BlastResult]:
        fasta_file=self._write_to_fasta_file(sequence, temp_dir)
        blast_hits: List[BlastResult]=self.run_from_file(fasta_file)
        return blast_hits