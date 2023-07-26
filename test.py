from inputs_validation import ValidateFiles
from fasta_utils import FastaUtilities
from generate_primers import PrimerGenerator

tests=ValidateFiles()
#fasta_file_name="/home/ubuntu/HandyAmpliconTool/test_data/test_data.fasta"
fasta_file_name="/home/ubuntu/HandyAmpliconTool/test_data/GCF_000195995.1.fna"
bed_file_name="/home/ubuntu/HandyAmpliconTool/test_data/test_data.bed"
tests.validate_bed(bed_file_name)
tests.validate_fasta(fasta_file_name)
tests.check_contigs_in_fasta(bed_file_name, fasta_file_name)
#/home/ubuntu/HandyAmpliconTool/test_data/amplified_regions.fasta

template_sequences=[]
with open(bed_file_name) as bed_file:
    fasta_parser=FastaUtilities()
    for line in bed_file:
        chr, start, end = line.split("\t")[0:3]
        template_sequences.append(fasta_parser.get_fasta_subseq(fasta_file_name,  chr, int(start), int(end)))


generator=PrimerGenerator()
for template in template_sequences:
    ### ADD test to check that template is long enough for design of primers
    primers=generator.generate_primers(template,(100, len(template)-100))
    break
