from inputs_validation import ValidateFiles
import name_converters
import time 
from multiprocessing import cpu_count, Pool

#### START Identification of genotype defining SNPS #### 
from os import listdir
import metadata_utils as metadata_utils
from load_vcfs import VCFutilities
from hierarchy_utils import HierarchyUtilities
import pandas as pd
from typing import Dict, List
import sys
from tqdm import tqdm

##### !!Testing inputs
name_converters.name_stubs.add(".sorted")
#vcf_dir: str="/home/ubuntu/HandyAmpliconTool/test_data/vcfs/"
vcf_dir: str="/home/ubuntu/converted_vcfs/"
meta_data_file="/home/ubuntu/HandyAmpliconTool/test_data/TGC_data.csv"
meta_deliminter=","
genotype_column="Final_genotype"
hierarchy_file="/home/ubuntu/HandyAmpliconTool/test_data/genotype_hierarcy.tsv"
repeat_regions_file: str="/home/ubuntu/HandyAmpliconTool/test_data/ref_repeats.bed"
##### !!Testing inputs

vcf_files: List[str]=[f'{vcf_dir}{f}' for f in listdir(vcf_dir) ][0:3000]
#vcf_files=["/home/ubuntu/converted_vcfs/32708_1#84.vcf"]
file_validator=ValidateFiles()
#file_validator.validate_many(vcf_files, "vcf")
if repeat_regions_file!="":
    file_validator.validate_bed(repeat_regions_file)
    file_validator.contigs_in_vcf(repeat_regions_file,vcf_files[0])
metadata_utils.load_metadata(meta_data_file,meta_deliminter)
samples_without_metadata=metadata_utils.samples_in_metadata(vcf_files)
for sample in samples_without_metadata:
    vcf_files.remove(sample)
metadata_utils.genotype_column=genotype_column #this will be an input

file_validator.validate_hierarchy(hierarchy_file)
hierarchy_utils=HierarchyUtilities()
hierarchy_utils.load_hierarchy(hierarchy_file)

vcf_utils=VCFutilities()
master_vcf=pd.DataFrame()

start_time=time.time()
vcfs: List[pd.DataFrame]=[]
with tqdm(total=len(vcf_files)) as progress_meter:
    for i, vcf in enumerate(vcf_files):
        vcfs.append( vcf_utils.load_file(vcf ) )
        progress_meter.update(1)

##get total indices from all vcfs to preallocate dataframe, this is much faster than merge and the loading of the vcfs can be parallelised
vcf_columns=[""]*len(vcfs)
with tqdm(total=len(vcfs)) as progress_meter:
    for i, vcf_data in enumerate(vcfs):
        vcf_columns[i]=vcf_data.columns[-1]
        progress_meter.update(1)
all_vcf_indices=sorted(set([k for f in vcfs for k in f.index]))
master_vcf=pd.DataFrame(index=all_vcf_indices, columns=vcf_columns)
#master_vcf.sort_index(inplace=True)

if __name__ == '__main__':
    # Multiprocessing doesn't speed up this step in short (few positions) datasets. Try later on longer dataset.
    # vcf_utils.master_vcf_temp=master_vcf
    # pool = Pool(processes=max(cpu_count()-1,1))
    # pool.map( vcf_utils.merge_vcfs, vcfs)
    # pool.close()
    # pool.join()
    
    with tqdm(total=len(vcfs)) as progress_meter:
        for vcf_data in vcfs:
            master_vcf.loc[vcf_data.index, vcf_data.columns[-1]]=vcf_data[vcf_data.columns[-1]]
            progress_meter.update(1)

    del vcf

    master_vcf.fillna("REF", inplace=True)
    #sys.exit()

    if repeat_regions_file!="":
        vcf_utils.remove_repeat_regions(master_vcf,repeat_regions_file)

    # pool = Pool(processes=max(cpu_count()-1,1))
    # pool.map( vcf_utils.merge_vcfs, vcfs)
    # pool.close()
    # pool.join()    
    genotype_bifurcating_snps=hierarchy_utils.find_defining_snps(master_vcf)

    #### !!!! For testing only
    print('dense : {:0.0f} bytes'.format(master_vcf.memory_usage().sum() / 1e3) )
    sdf = master_vcf.astype(pd.SparseDtype("str", "REF"))
    print('sparse: {:0.0f} bytes'.format(sdf.memory_usage().sum() / 1e3) )
    #master_vcf=master_vcf.astype(pd.SparseDtype("str", "REF"))
    #print(Counter([meta_data.get_metavalue(f,"Final_genotype") for f in samples]))

    gt_snp_df=pd.DataFrame(index=master_vcf.index, columns=list(genotype_bifurcating_snps.keys())).fillna(False)
    for gt in genotype_bifurcating_snps:
        gt_snp_df.loc[ genotype_bifurcating_snps[gt].index,  gt ]=genotype_bifurcating_snps[gt]["Pass"]

    #snp_to_drop=[index for index in gt_snp_df.index if True not in gt_snp_df.loc[ index ].values]

    #gt_snp_df.drop(index=snp_to_drop, inplace=True)

    #### END Identification of genotype defining SNPS #### 
    gt_snp_df.to_csv("/home/ubuntu/HandyAmpliconTool/test_data/test_gt_snps.tsv", sep="\t")