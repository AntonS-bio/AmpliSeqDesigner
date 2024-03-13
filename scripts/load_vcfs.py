import warnings
from typing import Dict, Tuple, Set, List
from data_classes import SNP, Sample, Genotypes, Genotype
from io import TextIOWrapper
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
                snps_to_add: List[SNP]=[]
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
                        snps_to_add.append(existing_snps[snp])
            sample.snps.extend(snps_to_add)



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
        
    def output_genotypes_vcf(self, genotypes: Genotypes, output_file: str) -> None:
        with open(output_file, "w") as vcf_output_file:
            #write the vcf header
            gt_columns=self._write_vcf_header(genotypes, vcf_output_file)

            coordinates=sorted(set([coordinate for genotype in genotypes.genotypes for coordinate in genotype.defining_snp_coordinates]))
            for contig_id, position in coordinates:
                snps_at_coordinates=[(genotype, snp) for genotype in genotypes.genotypes for snp in genotype.defining_snps if snp.position==position and snp.ref_contig_id==contig_id]
                alt_alleles=[base for base in [snp[1].alt_base for snp in snps_at_coordinates] ][0]
                if len(alt_alleles)>1:
                    print(f'Excess alleles at pos: {str(position)} contig {contig_id}')
                    continue
                alt_str=f'{alt_alleles}'
                for genotype, snp in snps_at_coordinates:
                    #check that snp is in multi genotype region
                    if snp.passes_filters:
                        if genotype.get_genotype_allele(snp)==snp.alt_base:
                            suffix=["0:."]*len(gt_columns)
                            for gt in genotype.subgenotypes:
                                if gt in gt_columns: #if genotype X and X.1 are separate targets, X.1 needs to show all SNPs of X,
                                    #but if only X is target, X.1 will not be an output column in VCF
                                    suffix[gt_columns[gt]]="1:"+str(genotype.get_genotype_allele_depth(snp)) #0 is REF allele
                            vcf_snp_id="_".join( ["GT",genotype.name,snp.ref_contig_id,str(snp.position+1)] )
                        else:
                            suffix=["1:."]*len(gt_columns) #set
                            for gt in genotype.subgenotypes:
                                if gt in gt_columns: #if genotype X and X.1 are separate targets, X.1 needs to show all SNPs of X,
                                    #but if only X is target, X.1 will not be an output column in VCF
                                    suffix[gt_columns[gt]]="0:"+str(genotype.get_genotype_allele_depth(snp)) #0 is REF allele
                            vcf_snp_id="_".join( ["GT",genotype.name,snp.ref_contig_id,str(snp.position+1)] )

                        vcf_output_file.write("\t".join([str(f) for f in [snp.ref_contig_id,
                                                                            snp.position+1,
                                                                            vcf_snp_id,
                                                                            snp.ref_base,
                                                                            alt_str,
                                                                            ".",
                                                                            "PASS",
                                                                            ".",
                                                                            "GT:DP",
                                                                            ]+suffix ]   ) +"\n" )


    def output_species_vcf(self, genotypes: Genotypes, output_file: str) -> None:
        with open(output_file, "w") as vcf_output_file:
            #write the vcf header
            gt_columns=self._write_vcf_header(genotypes, vcf_output_file, extra_columns=["NonTargetSerovar"])

            coordinates=sorted(set([coordinate for genotype in genotypes.genotypes for coordinate in genotype.defining_snp_coordinates]))
            for contig_id, position in coordinates:
                snps_at_coordinates=[(genotype, snp) for genotype in genotypes.genotypes for snp in genotype.defining_snps if snp.position==position and snp.ref_contig_id==contig_id]
                alt_alleles=[base for base in [snp[1].alt_base for snp in snps_at_coordinates] ][0]
                if len(alt_alleles)>1:
                    print(f'Excess alleles at pos: {str(position)} contig {contig_id}')
                    continue
                alt_str=f'{alt_alleles}'
                for genotype, snp in snps_at_coordinates:
                    #check that snp is in multi genotype region

                    if snp.passes_filters:
                        if genotype.name!="species":
                            if genotype.get_genotype_allele(snp)==snp.alt_base:
                                suffix=["1:."]*len(gt_columns)
                                suffix[gt_columns[genotype.name]]="1:"+str(genotype.get_genotype_allele_depth(snp)) #0 is REF allele
                            else:
                                suffix=["1:."]*len(gt_columns) #set
                                suffix[gt_columns[genotype.name]]="0:"+str(genotype.get_genotype_allele_depth(snp))
                            vcf_snp_id="_".join( ["GT",genotype.name,snp.ref_contig_id,str(snp.position+1)] )
                            suffix[gt_columns["species"]]=0
                            suffix[gt_columns["NonTargetSerovar"]]=".:."
                        else:
                            vcf_snp_id="_".join( ["Serovar", snp.ref_contig_id ,str(snp.position+1), snp.alt_base] )
                            suffix=["1:."]*len(gt_columns)
                            suffix[gt_columns["NonTargetSerovar"]]="1:"+str(genotypes.genotypes[-1].get_genotype_allele_depth(snp))
                            alt_str=snp.alt_base
                        vcf_output_file.write("\t".join([str(f) for f in [snp.ref_contig_id,
                                                                            snp.position+1,
                                                                            vcf_snp_id,
                                                                            snp.ref_base,
                                                                            alt_str,
                                                                            ".",
                                                                            "PASS",
                                                                            ".",
                                                                            "GT:DP",
                                                                            ]+suffix ]   ) +"\n" )

    def _write_vcf_header(self, genotypes: Genotypes, file_handle: TextIOWrapper, **kwargs) -> Dict[str, int]:
        extra_columns=kwargs.get("extra_columns", [])
        #write the vcf header
        file_handle.write('##fileformat=VCFv4.2'+"\n")
        file_handle.write('##FILTER=<ID=PASS,Description="All filters passed">'+"\n")
        file_handle.write('##ALT=<ID=*,Description="Represents allele(s) other than observed.">'+"\n")

        for ref_contig in set([snp.ref_contig_id for genotype in genotypes.genotypes for snp in genotype.defining_snps]):
            contig_max_position=max([snp.position for genotype in genotypes.genotypes for snp in genotype.defining_snps if snp.ref_contig_id == ref_contig])
            file_handle.write(f'##contig=<ID={ref_contig},length={str(contig_max_position)}>'+"\n")
        header_line="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        gt_columns={}
        for i, gt in enumerate(genotypes.genotypes):
            header_line=header_line+gt.name+"\t"
            gt_columns[gt.name]=i
        header_line=header_line+"\t"+"\t".join([f for f in extra_columns])+"\n"
        for column in extra_columns:
            gt_columns[column]=len(gt_columns)
        file_handle.write(header_line)
        return gt_columns
    
    @property
    def repeat_coordinates(self) -> Set[ Tuple[str, int] ]:
        return self._repeat_coordinates

