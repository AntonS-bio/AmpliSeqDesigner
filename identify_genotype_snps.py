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
vcf_dir: str="/home/ubuntu/HandyAmpliconTool/test_data/vcfs/"
meta_data_file="/home/ubuntu/HandyAmpliconTool/test_data/TGC_data.csv"
meta_deliminter=","
genotype_column="Final_genotype"
hierarchy_file="/home/ubuntu/HandyAmpliconTool/test_data/genotype_hierarcy.tsv"
##### !!Testing inputs

vcf_files: List[str]=[f'{vcf_dir}{f}' for f in listdir(vcf_dir) ]
repeat_regions_file: str=""
file_validator=ValidateFiles()
#file_validator.validate_many(vcf_files, "vcf")
if repeat_regions_file!="":
    file_validator.validate_bed(repeat_regions_file)
    file_validator.contigs_in_vcf(repeat_regions_file,vcf_files[0])
metadata_utils.load_metadata(meta_data_file,meta_deliminter)
metadata_utils.samples_in_metadata(vcf_files)
metadata_utils.genotype_column=genotype_column #this will be an input

vcf_utils=VCFutilities()
master_vcf=pd.DataFrame()

start_time=time.time()
vcfs: List[pd.DataFrame]=[]
with tqdm(total=len(vcf_files)) as progress_meter:
    for i, vcf in enumerate(vcf_files):
        vcf_to_add=vcf_utils.load_file(vcf )
        if master_vcf.shape[1]==0:
            master_vcf=vcf_to_add.copy()
        else:
            vcfs.append(vcf_to_add)
        progress_meter.update(1)

##get total indices from all vcfs to preallocate dataframe, this is much faster than merge and the loading of the vcfs can be parallelised
vcf_columns=[""]*len(vcfs)
indices=[]
for i, vcf_data in enumerate(vcfs):
    indices=indices+list(vcf_data.index)
    vcf_columns[i]=vcf_data.columns[-1]
master_vcf=pd.DataFrame(index=list(set(indices)), columns=vcf_columns)
master_vcf.sort_index(inplace=True)

if __name__ == '__main__':
    # Multiprocessing doesn't speed up this step in short (few positions) datasets. Try later on longer dataset.
    # vcf_utils.master_vcf_temp=master_vcf
    # pool = Pool(processes=max(cpu_count()-1,1))
    # pool.map( vcf_utils.merge_vcfs, vcfs)
    # pool.close()
    # pool.join()
    
    for i, vcf_data in enumerate(vcfs):
        master_vcf.loc[vcf_data.index, vcf_data.columns[-1]]=vcf_data[vcf_data.columns[-1]]
    print(time.time()-start_time)
    
    master_vcf.fillna("REF", inplace=True)
    #sys.exit()

    if repeat_regions_file!="":
        vcf_utils.remove_repeat_regions(master_vcf,repeat_regions_file)

    file_validator.validate_hierarchy(hierarchy_file)
    hierarchy_utils=HierarchyUtilities()
    hierarchy_utils.load_hierarchy(hierarchy_file)

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

    unique_bifurcating_snps=list(set([a for f in genotype_bifurcating_snps.values() for a in f]))
    gt_snp_df: pd.DataFrame=pd.DataFrame(index=unique_bifurcating_snps, columns=genotype_bifurcating_snps.keys()).fillna(0)
    for genotype in genotype_bifurcating_snps.keys():
        gt_snp_df.loc[genotype_bifurcating_snps[genotype],genotype]=1
    #### END Identification of genotype defining SNPS #### 
