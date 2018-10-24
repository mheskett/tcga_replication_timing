import pybedtools
import re
import pandas as pd
import matplotlib.pyplot as plt
import time


#def get_gc(bed)
# decide if you want to do everything from pandas or from text files....
#calculate this for whole genome sliding windows of size N


#replication
def simulate_links(length):

	#/Users/heskett/tcga_replication_timing/data
	"""
	returns some number of genomic windows of length=length
	and statistically SIMILAR gc,l1,snp

	"""
	start=time.time()

	snp = "/Users/heskett/tcga_replication_timing/data/ucsc.hg38.snp.nochr.sorted.bed"
	l1 = pybedtools.BedTool("/Users/heskett/tcga_replication_timing/data/l1.hg38.ucsc.nochr.bed")

	a=pybedtools.BedTool()
	windows=a.window_maker(g="/Users/heskett/tcga_replication_timing/data/chrom1.fa.fai",w=length,s=length/2)
	windows_nuc = windows.nucleotide_content(fi="/Users/heskett/tcga_replication_timing/data/hg38.nochr.fa")\
	.to_dataframe() ## 0th column has colnames. current colname is garbage
	gc=time.time()
	print(str(gc-start)+" time for windows")
	## now widows needs to go into below procedure
	windows_nuc.columns = windows_nuc.iloc[0]
	windows_nuc = windows_nuc.drop(0,axis="rows")
	keepers = [x for x in windows_nuc.columns if "user" in x or "pct_gc" in x]
	windows_nuc = windows_nuc.loc[:,keepers]
	nuct=time.time()
	print(str(nuct-gc)+" time for gc")
	## now get snp density
	windows_nuc = pybedtools.BedTool.from_dataframe(windows_nuc)
	windows_nuc_snp = windows_nuc.intersect(snp,c=True,sorted=True)
	snpt=time.time()
	print(str(snpt-nuct)+ " time for snp ") ## pre sorted file much faster sort -k1,1 -k2,2n
	## get L1

	windows_nuc_snp_l1 = windows_nuc_snp.coverage(l1)
	l1t=time.time()
	print(str(l1t-snpt)+" time for l1 coverage")
	## manipulate data frame
	df = windows_nuc_snp_l1.to_dataframe(names=["chrom","start","end","pct_gc",
	"num_snps","cov1","cov2","cov3","fraction_l1"])
	df["snps/kb"] = df["num_snps"] / (length/1000)


	final = df.loc[:,["chrom","start","end","snps/kb","pct_gc","fraction_l1"]]
	
	return final

#print(simulate_links(10**6,3,4,5))

def get_similar_links(df,snps=False,l1=False,gc=False): ## builds distribution OR inputs a real link
	# takes in DF of fake links of size of real link, and snp,l1,gc of REAL link
	df.describe(include=["float"],percentiles=[0.45,0.55]) # or coudl use quantile
	# get quantile OF real link, then add and subtract 0.05
	"""
	works with output of simulate links
	does boolean selection of simulate links data frame
	build distribution of snps/kb,pct_gc,fracl1, then select SIMILAR windows and return them
	"""
	quantiles = df.loc[:,["snps/kb","pct_gc","fraction_l1"]].quantile([0.45,.55]) # take nearest 10% to a real link!!!

	results = df[(df["snps/kb"].between(left=quantiles.loc[0.45,"snps/kb"],right=quantiles.loc[0.55,"snps/kb"])) &
	(df["pct_gc"].between(left=quantiles.loc[0.45,"pct_gc"],right=quantiles.loc[0.55,"pct_gc"])) &
	(df["fraction_l1"].between(left=quantiles.loc[0.45,"fraction_l1"],right=quantiles.loc[0.55,"fraction_l1"]))]

	print(results)
	return 

print(get_similar_links(simulate_links(10**5)))



# links = "/Users/mike/replication_tcga/data/links_final_grch38.bed"
# a = pybedtools.BedTool(links)
# result = a.nucleotide_content(fi="/Users/mike/replication_tcga/hg38.nochr.fa")\
# 	.to_dataframe()
# result.columns = result.iloc[0]
# result = result.drop(0,axis="rows")
# keepers = [x for x in result.columns if "user" in x or "pct_gc" in x]
# result = result.loc[:,keepers]

# snp = pybedtools.BedTool("/Users/mike/replication_tcga/data/snp150_ucsc_hg38_nochr.bed")
# l1 = pybedtools.BedTool("/Users/mike/replication_tcga/data/L1_ucsc_hg38_nochr.bed")
# b = pybedtools.BedTool.from_dataframe(result)

# result = b.intersect(snp,c=True)#.to_dataframe()

# result = result.coverage(l1)#.to_dataframe()

# df = result.to_dataframe(names=["chrom","start","end","name","length","pct_gc",
# 	"num_snps","cov1","cov2","cov3","fraction_l1"])

# df["snps/kb"] = df["num_snps"] / (df["length"]/1000)

# final = df.loc[:,["chrom","start","end","name","length","snps/kb","pct_gc","fraction_l1"]]
# #print(final)
# final.to_csv("links_annotated_grch38.bed",header=True,index=False,sep="\t")

# print(final.sort_values(by="fraction_l1"))
# print(final.sort_values(by="snps/kb"))
# print(final.sort_values(by="pct_gc"))

# df.hist(["length","snps/kb","pct_gc","fraction_l1"],bins=20)
plt.show()