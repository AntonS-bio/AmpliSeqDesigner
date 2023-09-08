import subprocess
from os.path import exists, isfile
from os import mkdir
from data_classes import BlastResult
from typing import List

class BlastRunner:
    def __init__(self) -> None:
        pass
    
    def db_from_file(self, file_name:str, db_dir: str):
        self.db_file_name = file_name
        self.db_dir=db_dir
        #create blast DB
        if exists(f'{self.db_dir}') and isfile(f'{self.db_dir}'):
            raise OSError(f'Cannot create directory {self.db_dir} because there is a file with such name')
        if not exists(f'{self.db_dir}'):
            mkdir(f'{self.db_dir}')
        outcome=subprocess.run(f'makeblastdb -in {self.db_file_name} -title temp -out {self.db_dir}/temp -dbtype nucl \
            -blastdb_version 4 1>/dev/null', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if outcome.returncode==1:
            raise OSError(f'Error generating blast database: {outcome.stderr}')

    def db_from_string(self, seq_header:str, sequence:str, db_dir: str):
        #Doesn't work due to <() having problems due to permissions, needs to write down the fasta file.
        self.db_dir=db_dir
        if seq_header[0]==">":
            input_str=seq_header+"\n"+sequence
        else:
            input_str=">"+seq_header+"\n"+sequence
        #create blast DB
        if exists(f'{self.db_dir}') and isfile(f'{self.db_dir}'):
            raise OSError(f'Cannot create directory {self.db_dir} because there is a file with such name')
        if not exists(f'{self.db_dir}'):
            mkdir(f'{self.db_dir}')
        outcome=subprocess.run(f'makeblastdb -in <({input_str}) -title temp -out {self.db_dir}/temp -dbtype nucl \
            -blastdb_version 4 1>/dev/null', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if outcome.returncode==1:
            raise OSError(f'Error generating BLAST database: {outcome.stderr}')

    def run_from_file(self, query_file:str) -> List[BlastResult]:
        #quickblast: reference, query, prints query
        blast_results=subprocess.run(f'blastn -query {query_file} -task \'megablast\' \
                -max_target_seqs 1000000000 -db {self.db_dir}/temp \
                -num_threads 1 -evalue 1.0E-5 -word_size 20 \
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