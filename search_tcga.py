import os
import csv
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import pybedtools
import re
import pybedtools
import time
import scipy.stats

## question is to see which ASAR's are most significantly altered in cancer

def simulate_links(length,window_fraction=0.25):

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

def get_similar_links(df,snps=False,l1=False,gc=False,wiggle=0.05,minimum=20):

	"""
	takes in DF of fake links of size of real link, and snp,l1,gc of REAL link
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

	if len(results) < minimum:
		wiggle+=0.05
		results = get_similar_links(df,snps=snps,gc=gc,l1=l1,wiggle=wiggle,minimum=minimum) # recursive biatch
	print(wiggle)

	
	# now check to see total number, if its too small, raise the quantile and go again
	# can use bedtools merge -d to merge overlapping features with minimum distance

	return results

def unique_name(x):
	# given a line of file, create unique name for link
	return str(x[3])+":"+str(x[0])+":"+str(x[1])+"-"+str(x[2])

def search_tcga(segments_df,link_chromosome,link_start,link_end,name):
	# segments data frame must have columns
	# chromosome, start, stop, patient, cancer_type, copy_number

	types = ["gain","loss","neutral","disruption"]
	results = {name:{i:pd.DataFrame() for i in types}}

	# should just be one data frame right?
	# with asar name, type of disruption, cancer type, patient sample name, copy number? 

	results[name]["loss"] = segments_df[(segments_df['chr'] == link_chromosome) 
							& (segments_df['start'] <= link_start)
							& (segments_df['stop'] >= link_end) 
							& (segments_df['copy_number'] <2.0 )]
	
	results[name]["gain"] = segments_df[(segments_df['chr'] == link_chromosome) 
							& (segments_df['start'] <= link_start)
							& (segments_df['stop'] >= link_end) 
							& (segments_df['copy_number'] > 2.0 )]
	
	results[name]["neutral"] = segments_df[(segments_df['chr'] == link_chromosome) 
							& (segments_df['start'] <= link_start) 
							& (segments_df['stop'] >= link_end) 
							& (segments_df['copy_number'] == 2.0 )]
	# max one per sample
	results[name]["disruption"] = segments_df[(segments_df['chr'] == link_chromosome) 
								& (segments_df['start'].between(left=link_start,right=link_end)
								| segments_df['stop'].between(left=link_start,right=link_end))]
	results[name]["disruption"] = results[name]["disruption"].drop_duplicates(subset="patient",keep="first")


	return results

#####################################################################

with open ("/Users/mike/replication_tcga/data/tcga_cancer_type_dictionary.txt") as f:
	lines = f.readlines()
	lines = (x.rstrip("\n").split("\t") for x in lines)
	cancer_atlas_dictionary = dict(lines)
	f.close()

with open ("/Users/mike/replication_tcga/data/test_link.bed") as g:
	links = g.readlines()
	links = [x.rstrip("\n").split("\t") for x in links[1:]] # cols are unique_name chrom start end name length snps pct gc l1
	links = [[str(x[3])+":"+str(x[0])+":"+str(x[1])+"-"+str(x[2]),
	str(x[0]),
	int(x[1]),
	int(x[2]),
	str(x[3]),
	int(x[4]),
	float(x[5]),
	float(x[6]),
	float(x[7])] for x in links]
	g.close()
with open ("/Users/mike/replication_tcga/data/TCGA.segtabs.final.merged.bed") as h:
	segments = h.readlines()
	segments = [x.rstrip("\n").split("\t") for x in segments]
	# segments = [[str(x[0]),
	# int(x[1]),
	# int(x[2]),
	# str(x[3]),
	# float(x[4]),
	# float(x[5]),
	# float(x[6])] for x in segments]
	h.close() # fast


print("filtering")
segments = [[str(x[0]),
	int(x[1]),
	int(x[2]),
	str(x[3]),
	str(cancer_atlas_dictionary[str(x[3][:-3])]),
	float(x[4])] for x in segments if x[3][:-3] in cancer_atlas_dictionary.keys()] # v slow...2min
segments_df= pd.DataFrame(segments)
segments_df.columns = ["chr","start","stop","patient","cancer_type","copy_number"]
types = ["gain","loss","neutral","disruption"]

links_names = [x[0] for x in links]

results = {k:{i:pd.DataFrame() for i in types} for k in links_names }
#fake_results = {k:{i:pd.DataFrame() for i in types} for k in links_names }
##########################################################
# here were going to create fake links to use for the search tcga function
print(links)

print("Starting loop")
for i in range(len(links)):

	fake_links = get_similar_links(simulate_links(length=links[i][5]),
									snps=links[i][6],
									l1=links[i][8],
									gc=links[i][7]
										)

	results.update(
			search_tcga(segments_df=segments_df,
					link_chromosome=links[i][1],
					link_start=links[i][2],
					link_end=links[i][3],
					name=links[i][0]))
print("done")
print(results)
print(fake_links)
##########################################################

## just counts tables for cancer type, may not even be useful
# cancer_types = list(np.unique([x for x in cancer_atlas_dictionary.values()]))
# samples = list(np.unique([x[3][:-3] for x in segments]))

# results[links[9][0]]["loss"]["cancer_type"].value_counts().plot(kind='bar')

# counts = {k:[] for k in links_names} ### to make counts table...
# for i in range(len(links)):
# 	counts[links[i][0]] += [len(results[links[i][0]]["gain"]),
# 	len(results[links[i][0]]["loss"]),
# 	len(results[links[i][0]]["disruption"]),
# 	len(results[links[i][0]]["neutral"])]

# pd.DataFrame.from_dict(counts,orient="index",columns=["gain","loss","disruption","neutral"]).to_csv("tcga_counts.txt",sep="\t")

# answer = pd.DataFrame()
# for i in range(len(types)):
# 	for j in range(len(links)):
# 		answer = answer.append(results[links[j][0]][types[i]])
# 	fig = answer["cancer_type"].value_counts().plot(kind='bar')
# 	plt.savefig(str(types[i])+".png")
# 	plt.close()
# 	print(answer)
# 	answer = pd.DataFrame()
