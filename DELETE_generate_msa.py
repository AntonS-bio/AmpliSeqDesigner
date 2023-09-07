from Bio import SeqIO
import subprocess
from os.path import exists
from os import remove
import pandas as pd
from collections import Counter
import sys

amlicon_sequences_files="/home/ubuntu/HandyAmpliconTool/test_data/multi_gt_intervals.fasta"
temp_fasta_file_dir="/home/ubuntu/HandyAmpliconTool/test_genomes/"
fasta_files_dir="/home/ubuntu/HandyAmpliconTool/test_genomes/"
msa_dir="/home/ubuntu/HandyAmpliconTool/test_data/msa/"
data_dir="/home/ubuntu/HandyAmpliconTool/test_data/"
ref_fasta=f'{data_dir}/GCA_000195995.1_for_VCFs.fasta'

amplicon_ids=[]
with open(f'{data_dir}/multi_gt_intervals.bed') as input_bed_file:
    for line in input_bed_file:
        print(line)
        bed_line_values = line.strip().split("\t")
        ampl_chr, ampl_start, ampl_end, ampl_target=bed_line_values[0:4]
        ampl_start=int(ampl_start)
        ampl_end=int(ampl_end)
        amplicon_id='_'.join( [str(f) for f in bed_line_values[0:4] ] )
        amplicon_ids.append(amplicon_id)
        for record in SeqIO.parse(ref_fasta,"fasta"):
            if record.id==ampl_chr:
                temp_fasta_filename=f'{temp_fasta_file_dir}/{amplicon_id}.fna'
                with open(temp_fasta_filename, "w") as temp_fasta:
                    temp_fasta.write(">"+amplicon_id+"\n")
                    temp_fasta.write(f'{str(record.seq[ampl_start:ampl_end])}'+"\n")
                # subprocess.call(f'python ~/utilities/buildGeneTree/buildGeneTree.py -r {temp_fasta_filename} -d {fasta_files_dir} -o {msa_dir} -l 10 -i 60 -c 4', shell=True)
                # if exists(temp_fasta_filename):
                #     remove(temp_fasta_filename)
                continue
#sys.exit()
#load the MSA as pandas dataframe
#preallocate df for speed
for amplicon_id in amplicon_ids:
    indices=[]
    columns=0
    for record in SeqIO.parse(f'{msa_dir}/{amplicon_id}.fna_aligned.fasta',"fasta"):
        indices.append(record.id)
        columns=len(str(record.seq))
    msa_df=pd.DataFrame(index=indices, columns=['pos_'+str(f) for f in range(0,columns)], dtype=str)
    for record in SeqIO.parse(f'{msa_dir}/{amplicon_id}.fna_aligned.fasta',"fasta"):
        msa_df.loc[record.id]=list(str(record.seq).upper())

    ##find at which positions left and right of the target SNPs the target sequence differs from the rest of genomes
    #target_index="2.2.2_2.5::AL513382_143N_pHCM1_120N_pHCM2:862-1818"
    segregating_columns=[]
    print(msa_df.shape)
    if msa_df.shape[0]==1:
        print(f'{amplicon_id} has no homologues')        
    else:
        for col in msa_df.columns:
            target_nucleotide=msa_df.loc[amplicon_id,col]
            if Counter(msa_df[col])[target_nucleotide]==1:
                segregating_columns.append(col)
        print(f'{amplicon_id} has {str(len(segregating_columns))}')