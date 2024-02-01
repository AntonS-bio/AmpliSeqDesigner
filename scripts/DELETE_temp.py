from data_classes import InputConfiguration
from primers_generator import PrimersGenerator


config_data = InputConfiguration("/home/lshas17/HandyAmpliconTool/test_data/configs/2.3.1.json")
generator=PrimersGenerator(config_data)
generator._for_testing_load_gts(config_data.species_data, config_data.genotypes_data)
generator.find_candidate_primers()
with open(config_data.output_dir+"primers.tsv","w") as output_file:
    header="\t".join(["Name", "Penalty", "Contig", "Target","Ref","Alt", "Start","End","Length",
                    "Forward","Forward Tm", "Forward GC",
                    "Reverse","Reverse Tm", "Reverse GC"])+"\n"
    output_file.write(header)
    for pair in generator.new_primer_pairs:
        output_file.write(pair.to_string()+"\n")