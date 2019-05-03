import os
import csv
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import pybedtools
import re
import time
import scipy.stats
import multiprocessing as mp
import sklearn
#"/Users/mike/replication_tcga" for personal use ONLY
if os.path.isdir("/Users/mike"):
	base_path = "/Users/mike/"
elif os.path.isdir("Users/heskett/"):
	base_path = "Users/heskett/"
elif os.path.isdir("/home/groups/SpellmanData/heskett/"):
	base_path = "/home/groups/SpellmanData/heskett/"
else:
	print("uknown machine")
	quit()

with open ("/Users/mike/replication_tcga/data/hg19/TCGA.segtabs.final.merged.sorted.no23.bed") as h:
	segments = h.readlines()
	segments = [x.rstrip("\n").split("\t") for x in segments]
	h.close() # fast

with open ("/Users/mike/replication_tcga/data/tcga_cancer_type_dictionary.txt") as f:
	lines = f.readlines()
	lines = (x.rstrip("\n").split("\t") for x in lines)
	cancer_atlas_dictionary = dict(lines)
	f.close()

cancer_types = list(np.unique([x for x in cancer_atlas_dictionary.values()]))


segments = [[str(x[0]),
	int(x[1]),
	int(x[2]),
	str(x[3]),
	str(cancer_atlas_dictionary[str(x[3][:-3])]),
	float(x[4])] for x in segments if x[3][:-3] in cancer_atlas_dictionary.keys()]
#print(segments)


def make_random_windows(length,number):
	# all file must be sorted sort -k1,1 -k2,2n
	# this creates N windows of length length and intersects them with TCGA segments to create a data frame
	a=pybedtools.BedTool()
	windows=a.random(g="/Users/mike/replication_tcga/data/hg19/TCGA.hg19.chromosome.file.bed",
							l=length,n=number)
	tcga = pybedtools.BedTool()
	tcga = tcga.from_dataframe(pd.DataFrame(segments))
	#df = windows.to_dataframe().loc[:,["chrom","start","end"]]
	
	# c = b.coverage(a) from bedtools doc. But coverage doesn't give a list of which features overlapped...
	# return windows.intersect(tcga,wa=True,wb=True).to_dataframe(names=["chrom","start","stop","arbitrary","length",
	# 	"strand","chrom_tcga","seg_start","seg_end","sample","cancer_type","copy_number"]).sort_values("sample")

	# ## i believe this returns a multi indexed pandas series?
	# return  windows.intersect(tcga,wa=True,wb=True).to_dataframe(names=["chrom","start","stop","arbitrary","length",
	# 	"strand","chrom_tcga","seg_start","seg_end","sample","cancer_type","copy_number"])\
	# 	.groupby(["chrom","start","stop","cancer_type"])["sample"].unique() # unique number of samples covered by each window per cancer type

	# return windows.intersect(tcga,wa=True,wb=True).to_dataframe(names=["chrom","start","stop","arbitrary","length",
	# 	"strand","chrom_tcga","seg_start","seg_end","sample","cancer_type","copy_number"])

	return windows
# #### main program
# df = pd.DataFrame()
# for i in range(10**5,3*10**5,5*10**3):
# 	df = pd.concat([df,make_random_windows(length=i,number=10)],ignore_index=True)
# #print(df)

def search_tcga(segments_df,link_chromosome,link_start,link_end):
	types = ["gain","loss","neutral","disruption"]
	segments_df=segments_df[segments_df['chr']==link_chromosome] # then remove the additional operation from below...

	losses = segments_df[(segments_df['start'] <= link_start)
							& (segments_df['stop'] >= link_end) 
							& (segments_df['copy_number'] <2.0 )]

	gains = segments_df[(segments_df['start'] <= link_start)
							& (segments_df['stop'] >= link_end) 
							& (segments_df['copy_number'] > 2.0 )]


	neutrals = segments_df[(segments_df['start'] <= link_start) 
							& (segments_df['stop'] >= link_end) 
							& (segments_df['copy_number'] == 2.0 )]
	# max one per sample
	disruptions = segments_df[(segments_df['start'].between(left=link_start,right=link_end)
								| segments_df['stop'].between(left=link_start,right=link_end))]
	disruptions = disruptions[disruptions.duplicated(subset="patient",keep=False)] # marks all dupes as true cause we want to keep dupes
	disruptions = disruptions.drop_duplicates(subset="patient",keep="first") # keeps the first dupe
	
	losses.loc[:,"type"] = pd.Series(["loss"]*len(losses.index),index=losses.index).astype(str) # empty list was being assigned float64 type...
	gains.loc[:,"type"] = pd.Series(["gain"]*len(gains.index),index=gains.index).astype(str)
	disruptions.loc[:,"type"]= pd.Series(["disruption"]*len(disruptions.index),index=disruptions.index).astype(str)
	neutrals.loc[:,"type"] = pd.Series(["neutral"]*len(neutrals.index),index=neutrals.index).astype(str)

	return pd.concat([losses,gains,disruptions,neutrals])


print(make_random_windows(length=2.5*10**5,number=1))
segments_df= pd.DataFrame(segments)#.astype({0: str})
segments_df.columns = ["chr","start","stop","patient","cancer_type","copy_number"]
test = search_tcga(segments_df,link_chromosome="2",link_start=126269818,link_end=126519818)

print(test)

for k in range(len(cancer_types)):	
	counts = {}
	cancer_type_coverage = test[(test["cancer_type"]==cancer_types[k])]["patient"].nunique()
	print(cancer_types[k],cancer_type_coverage)

# for k in range(len(cancer_types)):	
# 	counts = {}
# 	cancer_type_coverage = test[(test["cancer_type"]==cancer_types[k])]["type"].value_counts()
# 	print(cancer_types[k],cancer_type_coverage)
#print(segments_df.groupby(["cancer_type"])["patient"].nunique()) # total samples
