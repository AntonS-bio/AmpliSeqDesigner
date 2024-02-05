import warnings
from typing import Dict, Tuple, Set
from data_classes import SNP, Sample
## Consider replacing some of this with GATKs VariantsToTable

   
class VCFutilities():

    def __init__(self) -> None:
        self._repeat_coordinates=set()
        pass


    def determine_vcf_type(self, filename: str) -> str:
        """Checks whether VCF is single- or multisample.

        :param filename: full name to VCF file
        :type filename: str
        :return: VCF type as either single_sample or multi_sample
        """
        vcf_file_type="Unknown"
        with open(filename) as vcf_file:
            for line in vcf_file:
                if line[0:6]=="#CHROM":
                    header_columns=line.strip().split("\t")
                    if len(header_columns)==10:
                        vcf_file_type="single_sample"
                    elif len(header_columns)>10:
                        vcf_file_type="multi_sample"
                    else:
                        raise ValueError(f'VCF file {filename} has fewer than 10 columns which is minimum required')
                    break
        if vcf_file_type=="Unknown":
            raise ValueError(f'Header line not found in VCF file {filename}. Were looking for line starting with #CHROM')
        return vcf_file_type

    def vcf_to_snps(self, filename: str, existing_snps:Dict[SNP, SNP], sample: Sample):
        """Converts VCF file into SNPs. The SNPs will be added to the inputted sample object 
        
        :param filename: full name to VCF file
        :type filename: str
        :param existing_snps: Dictionary of SNPs loaded from previous VCF files
        :type existing_snps:  Dict[SNP, SNP]
        :param sample: Dictionary of SNPs loaded from previous VCF files
        :type sample: Sample that contains VCF data 
        """
        vcf_file_type=self.determine_vcf_type(filename)

        multiploid_positions=[]
        if vcf_file_type=="single_sample":
            #vcf_datatypes={"CHROM":"string","POS":int, "REF": "string", "ALT": "string", "FORMAT": "string"}
            with open(filename) as vcf_file_handle:
                for line in vcf_file_handle:
                    if line[0]=="#" or line=="\n":
                        continue
                    chrom, pos, id, ref, alt, qual, filter, info, format, values=line.strip().split("\t")
                    pos=int(pos)-1  #the rest of the code is 0 indexed like BED and BAM, but VCF coordinates are 1-indexed
                    if (chrom, pos) in self.repeat_coordinates:
                        continue
                    if alt.find(",")>-1:
                        #multiploid line, skip with a warning
                        multiploid_positions.append((chrom, pos))
                        continue
                    gt_index=format.split(":").index("GT")
                    if gt_index==-1:
                        raise ValueError(f'VCF file {filename} does not have genotype code [GT] in SAMPLE colum at {chrom} {str(pos-1)}')

                    allele = ref if values.split(",")[gt_index].split(":")[0]=="0" else alt
                    snp=SNP(ref_contig_id=chrom, ref_base=ref, alt_base=allele, position=pos)
                    if snp not in existing_snps:
                        existing_snps[snp]=snp
                        
                    if snp not in sample.snps:
                        sample.snps.append(existing_snps[snp])



            if len(multiploid_positions)>0:
                warnings.warn(f'VCF file {filename} has {len(multiploid_positions)} duplicated positions')

        else:
            raise ValueError(f'The vcf type {vcf_file_type} is not currently supported')


    def load_repeat_regions(self, bed_file: str) -> bool:
        self._repeat_coordinates: Set[ Tuple[str, int] ] =set()
        if bed_file=="":
            return True
        with open(bed_file) as bed_data:
            for line in bed_data:
                values=line.strip().split("\t")
                self._repeat_coordinates.update( [(values[0],f) for f in range(int(values[1]), int(values[2]))] )
        return True
        
    
    @property
    def repeat_coordinates(self) -> Set[ Tuple[str, int] ]:
        return self._repeat_coordinates

