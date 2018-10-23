import pybedtools
import re
import pandas as pd
import matplotlib.pyplot as plt

#def get_gc(bed)
# decide if you want to do everything from pandas or from text files....
#calculate this for whole genome sliding windows of size N

def simulate_links(length,gc,l1,snp):
	"""
	returns some number of genomic windows of length=length
	and statistically SIMILAR gc,l1,snp
	"""
	return


links = "/Users/mike/replication_tcga/data/links_final_grch38.bed"
a = pybedtools.BedTool(links)
result = a.nucleotide_content(fi="/Users/mike/replication_tcga/hg38.nochr.fa")\
	.to_dataframe()
result.columns = result.iloc[0]
result = result.drop(0,axis="rows")
keepers = [x for x in result.columns if "user" in x or "pct_gc" in x]
result = result.loc[:,keepers]

snp = pybedtools.BedTool("/Users/mike/replication_tcga/data/snp150_ucsc_hg38_nochr.bed")
l1 = pybedtools.BedTool("/Users/mike/replication_tcga/data/L1_ucsc_hg38_nochr.bed")
b = pybedtools.BedTool.from_dataframe(result)

result = b.intersect(snp,c=True)#.to_dataframe()

result = result.coverage(l1)#.to_dataframe()

df = result.to_dataframe(names=["chrom","start","end","name","length","pct_gc",
	"num_snps","cov1","cov2","cov3","fraction_l1"])

df["snps/kb"] = df["num_snps"] / (df["length"]/1000)

final = df.loc[:,["chrom","start","end","name","length","snps/kb","pct_gc","fraction_l1"]]
#print(final)
final.to_csv("links_annotated_grch38.bed",header=True,index=False,sep="\t")

print(final.sort_values(by="fraction_l1"))
print(final.sort_values(by="snps/kb"))
print(final.sort_values(by="pct_gc"))

df.hist(["length","snps/kb","pct_gc","fraction_l1"],bins=20)
plt.show()