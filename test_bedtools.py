import pybedtools
import re
import pandas as pd
import matplotlib.pyplot as plt
import time
import scipy.stats

def simulate_links(length,window_fraction=0.1):

	"""	
	#/Users/heskett/tcga_replication_timing/data
	returns some number of genomic windows of length=length
	and statistically SIMILAR gc ,l1,snp
	"""

	start=time.time()
	# snp doesnt need to be bedtool object yet i guess
	snp = "/Users/mike/replication_tcga/data/snp150.ucsc.hg38.nochr.sorted.bed"
	l1 = pybedtools.BedTool("/Users/mike/replication_tcga/data/L1_ucsc_hg38_nochr.bed")

	# initiate bed tool and make windows
	# all file must be sorted sort -k1,1 -k2,2n
	a=pybedtools.BedTool()
	windows=a.window_maker(g="/Users/mike/replication_tcga/data/hg38.nochr.sorted.fa.fai",
							w=length,s=length*window_fraction)
	windows_nuc = windows.nucleotide_content(fi="/Users/mike/replication_tcga/data/hg38.nochr.fa")\
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

	## get snp density
	windows_nuc = pybedtools.BedTool.from_dataframe(windows_nuc)
	windows_nuc_snp = windows_nuc.intersect(snp,c=True,sorted=True)
	snpt=time.time()
	print(str(snpt-nuct)+ " time for snp ") ## pre sorted file much faster sort -k1,1 -k2,2n
	
	## get L1
	windows_nuc_snp_l1 = windows_nuc_snp.coverage(l1)
	l1t=time.time()
	print(str(l1t-snpt)+" time for l1 coverage")
	
	## manipulate data frame
	df = windows_nuc_snp_l1.to_dataframe(names=["chrom","start","end",
										"pct_gc","num_snps","cov1",
										"cov2","cov3","fraction_l1"])
	# create snps per kb column
	df["snps/kb"] = df["num_snps"] / (length/1000)

	# return columns of interest
	final = df.loc[:,["chrom","start","end",
					"snps/kb","pct_gc","fraction_l1"]]
	
	return final

def get_similar_links(df,snps=False,l1=False,gc=False,wiggle=0.05):
	# takes in DF of fake links of size of real link, and snp,l1,gc of REAL link
	"""
	works with output of simulate links
	does boolean selection of simulate links data frame
	build distribution of snps/kb,pct_gc,fracl1, then select SIMILAR windows and return them


	"""
	# get percentile score of input parameters in fake link data frame
	snps_score,gc_score,l1_score = (scipy.stats.percentileofscore(df["snps/kb"], score=snps, kind='rank'),
									scipy.stats.percentileofscore(df["pct_gc"], score=gc, kind='rank'),
									scipy.stats.percentileofscore(df["fraction_l1"], score=l1, kind='rank'))

	results = df[(df["snps/kb"]\
			.between(left=df["snps/kb"].quantile(q=(snps_score/100)-wiggle),
					right=df["snps/kb"].quantile(q=(snps_score/100)+wiggle))) 
		& (df["pct_gc"]\
			.between(left=df["pct_gc"].quantile(q=(gc_score/100)-wiggle),
					right=df["pct_gc"].quantile(q=(gc_score/100)+wiggle)))
		& (df["fraction_l1"]\
			.between(left=df["fraction_l1"].quantile(q=(l1_score/100)-wiggle),
					right=df["fraction_l1"].quantile(q=(l1_score/100)+wiggle)))]
	
	# now check to see total number, if its too small, raise the quantile and go again
	# can use bedtools merge -d to merge overlapping features with minimum distance

	return results

print(get_similar_links(simulate_links(length=579803),snps=5.179,l1=0.175,gc=.4))