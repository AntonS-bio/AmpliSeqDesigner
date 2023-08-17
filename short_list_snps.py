import pandas as pd
import numpy as np
from collections import Counter

gt_snp_df=pd.read_csv("/home/ubuntu/HandyAmpliconTool/test_data/test_gt_snps.tsv", sep="\t", index_col=[0])

snp_to_drop=[index for index in gt_snp_df.index if True not in gt_snp_df.loc[ index ].values]

gt_snp_df.drop(index=snp_to_drop, inplace=True)

gt_snp_df["Chr"]=[index.split(",")[0].replace("(","").replace("'","") for index in gt_snp_df.index]
gt_snp_df["Pos"]=[int(index.split(",")[1].replace(")","")) for index in gt_snp_df.index]
gt_snp_df.set_index(["Chr", "Pos"], inplace=True)
### Collect all snps that fall within 1000nt intervals
### This is to short list SNPs for analysis

snp_interval=1000
### Take first SNP and continue checking SNPs until range exceeds 1,000
interval_snps=[]
sorted_snps=sorted(gt_snp_df.index, key=lambda x: (x[0], x[1]), reverse=False)
i=0
max_reached_index=0
while i<len(sorted_snps):
    j=i #this means the first snp is automatically added
    while j<len(sorted_snps) and sorted_snps[j][0]==sorted_snps[i][0] and sorted_snps[j][1]-sorted_snps[i][1]<snp_interval:
        j+=1
    if j>i+1:
        if j>max_reached_index: #this account for some intervals containing >2 snps. Without this, these intervals would create multiple entries in interval_snps
            interval_snps.append( {"snps": sorted_snps[i:j], "genotypes":[] } ) # j, not j-1 because python excludes the last element of index
            max_reached_index=j-1
    i+=1
# check which GTs are captured by which lists of SNPs
# Remove those that capture same GT multiple times - this is likely due to structural variant
for interval in interval_snps:
    interval["genotypes"]=set(gt_snp_df.loc[interval["snps"]].apply(lambda row: row[row==True], axis=1)) #row==True to be clear what's being tested

interval_snps=[snp_interval for snp_interval in interval_snps if len(snp_interval["genotypes"])>1 ]

with open("/home/ubuntu/HandyAmpliconTool/test_data/multi_gt_intervals.bed", "w") as output_bed:
    for inteval in interval_snps:
        interval_start=min([f[1] for f in inteval["snps"]])
        interval_end=max([f[1] for f in inteval["snps"]])
        interval_len=min(interval_end-interval_start,snp_interval)
        interval_start=int((interval_start-(snp_interval-interval_len))/2)
        interval_end=int((interval_end+(snp_interval-interval_len))/2)
        chr=inteval["snps"][0][0]
        gts="_".join(inteval["genotypes"])
        output_bed.write('\t'.join([chr, str(interval_start),str(interval_end),gts])+"\n")

print("A")