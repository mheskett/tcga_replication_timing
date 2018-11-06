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


def min_max(x):
	return min(max(0,x),1)
def faster_simulate_links(length,window_fraction=0.25,snps=False,l1=False,gc=False,wiggle=0.05,minimum=30):
	"""
	returns a list of fake genes, from a length
	predicated on the fact that bedtools.nuc is the only slow step in this analysis
	"""

	start=time.time()
	# snp doesnt need to be bedtool object yet i guess
	snp_file = "/Users/mike/replication_tcga/data/snp150.ucsc.hg38.nochr.sorted.bed"
	l1_file = pybedtools.BedTool("/Users/mike/replication_tcga/data/L1_ucsc_hg38_nochr.bed")

	# initiate bed tool and make windows
	# all file must be sorted sort -k1,1 -k2,2n
	a=pybedtools.BedTool()
	windows=a.window_maker(g="/Users/mike/replication_tcga/data/hg38.nochr.sorted.fa.fai",
							w=length,s=length*window_fraction)
	wint=time.time()
	print(wint-start," window time")

	windows_snp = windows.intersect(snp_file,c=True,sorted=True) # makes one column
	snpt=time.time()
	print(snpt-wint," snp time")
	#print(str(snpt-gc)+ " time for snp ") ## pre sorted file much faster sort -k1,1 -k2,2n
	
	## get L1
	windows_snp_l1 = windows_snp.coverage(l1_file)# produces 4 columns
	l1t=time.time()
	print(l1t-snpt," l1 time")


	## always do nuc on a reduced set of windows
	df = windows_snp_l1.to_dataframe(names=["chrom","start","end",
										"num_snps","cov1",
										"cov2","cov3","fraction_l1"])
	df["snps/kb"] = df["num_snps"] / (length/1000)

	snps_score,l1_score = (scipy.stats.percentileofscore(a=df["snps/kb"], score=snps, kind='rank'), # this giving seg fault?
							scipy.stats.percentileofscore(a=df["fraction_l1"], score=l1, kind='rank'))
	df = df[
		(df["snps/kb"]\
			.between(left=df["snps/kb"].quantile(q=min_max(snps_score/100-wiggle)),
					right=df["snps/kb"].quantile(q=min_max(snps_score/100+wiggle))))\
		& (df["fraction_l1"]\
			.between(left=df["fraction_l1"].quantile(q=min_max(l1_score/100-wiggle)), 
					right=df["fraction_l1"].quantile(q=min_max(l1_score/100+wiggle))))]
	
	b = pybedtools.BedTool()
	window_snp_l1_reduced = b.from_dataframe(df)
	window_snp_l1_nuc =  window_snp_l1_reduced.nucleotide_content(fi="/Users/mike/replication_tcga/data/hg38.nochr.fa")\
							.to_dataframe() ## 0th column has colnames. current colname is garbage
	nuct=time.time()
	print(nuct-snpt," nuc time")
	window_snp_l1_nuc.columns = window_snp_l1_nuc.iloc[0]
	window_snp_l1_nuc = window_snp_l1_nuc.drop(0,axis="rows")
	keepers = [x for x in window_snp_l1_nuc.columns if "user" in x or "pct_gc" in x]
	window_snp_l1_nuc = window_snp_l1_nuc.loc[:,keepers]
	window_snp_l1_nuc.columns = ["chrom","start","end",
										"num_snps","cov1",
										"cov2","cov3","fraction_l1","snps/kb","pct_gc"
										]
	window_snp_l1_nuc[["pct_gc","fraction_l1","snps/kb","start","end"]] = \
	window_snp_l1_nuc[["pct_gc","fraction_l1","snps/kb","start","end"]].apply(pd.to_numeric)
	gc_score = scipy.stats.percentileofscore(window_snp_l1_nuc["pct_gc"], score=gc, kind='rank')
	window_snp_l1_nuc = window_snp_l1_nuc[window_snp_l1_nuc["pct_gc"]\
			.between(left=window_snp_l1_nuc["pct_gc"].quantile(q=min_max(gc_score/100-wiggle)),
					right=window_snp_l1_nuc["pct_gc"].quantile(q=min_max(gc_score/100+wiggle)))]
	########
	if len(window_snp_l1_nuc) < minimum:
		wiggle+=0.05
		results = faster_simulate_links(length=length,snps=snps,gc=gc,l1=l1,wiggle=wiggle,minimum=minimum) # recursive biatch
	print(wiggle, " wiggle")


	
	return window_snp_l1_nuc.loc[:,["chrom","start","end",
					"snps/kb","pct_gc","fraction_l1"]]



def get_similar_links(df,snps=False,l1=False,gc=False,wiggle=0.05,minimum=30):

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

	results = df[
		(df["snps/kb"]\
			.between(left=df["snps/kb"].quantile(q=min_max(snps_score/100-wiggle)),
					right=df["snps/kb"].quantile(q=min_max(snps_score/100+wiggle))))\
		& (df["pct_gc"]\
			.between(left=df["pct_gc"].quantile(q=min_max(gc_score/100-wiggle)),
					right=df["pct_gc"].quantile(q=min_max(gc_score/100+wiggle))))\
		& (df["fraction_l1"]\
			.between(left=df["fraction_l1"].quantile(q=min_max(l1_score/100-wiggle)), # this is going below zero i think....
					right=df["fraction_l1"].quantile(q=min_max(l1_score/100+wiggle))))]

	if len(results) < minimum:
		wiggle+=0.05
		results = get_similar_links(df,snps=snps,gc=gc,l1=l1,wiggle=wiggle,minimum=minimum) # recursive biatch
	print(wiggle)

	
	# now check to see total number, if its too small, raise the quantile and go again
	# can use bedtools merge -d to merge overlapping features with minimum distance

	return results

def search_tcga(segments_df,link_chromosome,link_start,link_end,name):
	types = ["gain","loss","neutral","disruption"]
	results = {name:pd.DataFrame()}

	losses = segments_df[(segments_df['chr'] == link_chromosome) # v2 just appends these
							& (segments_df['start'] <= link_start)
							& (segments_df['stop'] >= link_end) 
							& (segments_df['copy_number'] <2.0 )]

	gains = segments_df[(segments_df['chr'] == link_chromosome) 
							& (segments_df['start'] <= link_start)
							& (segments_df['stop'] >= link_end) 
							& (segments_df['copy_number'] > 2.0 )]


	neutrals = segments_df[(segments_df['chr'] == link_chromosome) 
							& (segments_df['start'] <= link_start) 
							& (segments_df['stop'] >= link_end) 
							& (segments_df['copy_number'] == 2.0 )]
	# max one per sample
	disruptions = segments_df[(segments_df['chr'] == link_chromosome) 
								& (segments_df['start'].between(left=link_start,right=link_end)
								| segments_df['stop'].between(left=link_start,right=link_end))]
	disruptions = disruptions.drop_duplicates(subset="patient",keep="first")


	
	losses.loc[:,"type"] = pd.Series(["loss"]*len(losses.index),index=losses.index)
	gains.loc[:,"type"] = pd.Series(["gain"]*len(gains.index),index=gains.index)
	disruptions.loc[:,"type"]= pd.Series(["disruption"]*len(disruptions.index),index=disruptions.index)
	neutrals.loc[:,"type"] = pd.Series(["neutral"]*len(neutrals.index),index=neutrals.index)



	results[name] = pd.concat([losses,gains,disruptions,neutrals])
	return results


## re write main program here
with open ("/Users/mike/replication_tcga/data/tcga_cancer_type_dictionary.txt") as f:
	lines = f.readlines()
	lines = (x.rstrip("\n").split("\t") for x in lines)
	cancer_atlas_dictionary = dict(lines)
	f.close()

with open ("/Users/mike/replication_tcga/data/links_annotated_grch38.bed") as g:
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
cancer_types = list(np.unique([x for x in cancer_atlas_dictionary.values()]))
links_names = [x[0] for x in links]

## this has real links as K, and a dict of disruption types with values pandas data frame as value
#fake_results = {k:{i:pd.DataFrame() for i in types} for k in links_names }
##########################################################
# here were going to create fake links to use for the search tcga function
print("Starting loop")
output=[["link_name","cancer_type","p_value_gain","p_value_loss","p_value_disruption","p_value_neutral"]]# link  



def get_answer(links,segments_df,cancer_types):
	links_names = [x[0] for x in links]
	results = {k:pd.DataFrame() for k in links_names } # each k is a real link
	output=[["link_name","cancer_type","p_value_gain","p_value_loss","p_value_disruption","p_value_neutral"]]# link name, cancer type, p value gain, p value loss, p value disruption, p value neutral
	for i in range(len(links)):
		print(links[i][0]," ",i)
		fake_links= faster_simulate_links(length=links[i][5],
									snps=links[i][6],
									l1=links[i][8],
									gc=links[i][7])
		results.update(
				search_tcga(segments_df=segments_df,
						link_chromosome=links[i][1],
						link_start=links[i][2],
						link_end=links[i][3],
						name=links[i][0]))


		# get results for the simulated links
		fake_link_names = ["fake_link_"+str(row['chrom'])+":"+str(row["start"])+"-"+str(row["end"]) for index,row in fake_links.iterrows()]
		fake_results = {k:pd.DataFrame()for k in fake_link_names} 
		for index,row in fake_links.iterrows():
			fake_results.update(
				search_tcga(segments_df=segments_df,
							link_chromosome=row['chrom'],
							link_start=row['start'],
							link_end=row['end'],
							name="fake_link_"+str(row['chrom'])+":"+str(row["start"])+"-"+str(row["end"]) # these are named here
							))
		print("done search tcga fake gnes")

		counts = {} ### to make counts table...do cancer specific
		df = results[links[i][0]]

		counts[links[i][0]] = df.groupby(['cancer_type', 'type']).size()
		print(counts[links[i][0]])

		fake_counts = {} ### make sure this matches up with fake links names
		for j in range(len(fake_link_names)):
			fake_df = fake_results[fake_link_names[j]]
			fake_counts[fake_link_names[j]] = fake_df.groupby(['cancer_type', 'type']).size()
		# might be faster? have to learn how to index through the group by object thing

	return

get_answer(links=links,segments_df=segments_df,cancer_types=cancer_types)
