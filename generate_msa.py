from Bio import SeqIO
import subprocess
from os.path import exists
from os import remove

amlicon_sequences_files="/home/ubuntu/HandyAmpliconTool/test_data/multi_gt_intervals.fasta"
temp_fasta_file_dir="/home/ubuntu/HandyAmpliconTool/test_genomes/"
fasta_files_dir="/home/ubuntu/HandyAmpliconTool/test_genomes/"
msa_dir="/home/ubuntu/HandyAmpliconTool/test_data/msa/"

for record in SeqIO.parse(amlicon_sequences_files,"fasta"):
    print(record.id)
    temp_fasta_filename=f'{temp_fasta_file_dir}/{record.id}.fna'
    with open(temp_fasta_filename, "w") as temp_fasta:
        temp_fasta.write(f'>{record.id}'+"\n")
        temp_fasta.write(f'{str(record.seq)}'+"\n")
    subprocess.call(f'python ~/utilities/buildGeneTree/buildGeneTree.py -r {temp_fasta_filename} -d {fasta_files_dir} -o {msa_dir} -l 10 -i 60 -c 4', shell=True)
    if exists(temp_fasta_filename):
        remove(temp_fasta_filename)
    break