import os
import csv
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import pybedtools
import re

### check how many samples are not in sample manifest....


def unique_name(x):
	# given a line of file, create unique name for link
	return str(x[3])+":"+str(x[0])+":"+str(x[1])+"-"+str(x[2])


########## Load necessary files
platform_name=sys.platform

# "/Users/heskett/tcga_replication_timing/data/"

if platform_name=="linux":
	working_directory =  "/home/exacloud/lustre1/SpellmanLab/heskett/tcga_replication_timing/data/"
if platform_name=="darwin":
    working_directory = "/Users/mike/replication_tcga/data/"

cancer_atlas_path = working_directory + "TCGA_mastercalls.abs_segtabs.fixed.txt"
asar_path = working_directory + "all_links.txt"
cancer_atlas_dictionary_path = working_directory + "tcga_cancer_type_dictionary.txt"

# cancer_atlas = pd.read_table(cancer_atlas_path,header=0)
# asar = pd.read_table(asar_path,names=["chr","start","end","annotation","length"],
# 		dtype={"chr":str,"start":int,"end":int,"annotation":str,"length":float})
# #asar.columns=["chr","start","end","annotation","length"]
	#dtype={"chr":str,"start":int,"end":int,"annotation":str,"length":float})
with open ("/Users/mike/replication_tcga/data/tcga_cancer_type_dictionary.txt") as f:
	lines = f.readlines()
	lines = (x.rstrip("\n").split("\t") for x in lines)
	cancer_atlas_dictionary = dict(lines)
	f.close()

## create list of cancer types

cancer_types = list(np.unique([x for x in cancer_atlas_dictionary.values()]))
print(cancer_types)
with open ("/Users/mike/replication_tcga/data/links_segments_overlap_grch38.bed") as f:
	lines = f.readlines()
	intersections = [x.rstrip("\n").split("\t") for x in lines]
	intersections = [[str(x[0]),
	int(x[1]),
	int(x[2]),
	str(x[3]),
	int(x[4]),
	str(x[5]),
	int(x[6]),
	int(x[7]),
	str(x[8]),
	float(x[9]),
	float(x[10]),
	float(x[11]),
	int(x[12])] for x in intersections]
	#print(intersections)
	f.close()

#print(intersections)

all_links = pd.read_table("/Users/mike/replication_tcga/data/links_final.bed",
	index_col=None,
	names=["chr","start","stop","name","length"])
#print(all_links)

print(all_links["length"].describe())
all_links.hist("length")
#plt.show()


links = list(np.unique([str(x[3])+":"+str(x[0])+":"+str(x[1])+"-"+str(x[2]) for x in intersections]))
samples = list(np.unique([x[8][:-3] for x in intersections]))
types = ["gain","loss","neutral"]
#print(cancer_atlas_dictionary.keys())
remove_samples = [x for x in samples if x not in cancer_atlas_dictionary.keys()]
#normals = [x[8] for x in intersections if x[-3:]=="-10"]
#print(normals)
#print(remove_samples)
#results = {k:{l : [] for l in types} for k in links } #[[sample_name,start,stop,cancertype]]
#print(results)
results = {k:[] for k in links }

for i in range(len(intersections)):
	if intersections[i][8][:-3] in remove_samples:
		continue
	if (intersections[i][11]<2.0):
		results[unique_name(intersections[i])]+=[[intersections[i][8],"loss",cancer_atlas_dictionary[intersections[i][8][:-3]],intersections[i][5],intersections[i][6],intersections[i][7]]]
	if (intersections[i][11]>2.0):
		results[unique_name(intersections[i])]+=[[intersections[i][8],"gain",cancer_atlas_dictionary[intersections[i][8][:-3]],intersections[i][5],intersections[i][6],intersections[i][7]]]
	if (intersections[i][12] == intersections[i][4]) and (intersections[i][11]==2.0):
		results[unique_name(intersections[i])]+=[[intersections[i][8],"neutral",cancer_atlas_dictionary[intersections[i][8][:-3]],intersections[i][5],intersections[i][6],intersections[i][7]]]
#
#	if (intersections[i][12] < intersections[i][4]) and (intersections[i][11]!=2.0):
#		results[unique_name(intersections[i])]+=[[intersections[i][8],"disruption",cancer_atlas_dictionary[intersections[i][8][:-3]],intersections[i][5],intersections[i][6],intersections[i][7]]]

results_df = {} # dictionary of data frames for each link
for i in range(len(links)):
	tmp=pd.DataFrame.from_dict(results[links[i]])
	tmp.columns = ["sample","type","cancer_type","chr","start","end"]
	results_df[links[i]]=tmp
#print(results_df["VLINC273:6:141034658-141219628"])
#results_df["VLINC273:6:141034658-141219628"].to_csv("vlink273.txt")

#test_df = results_df["VLINC273:6:141034658-141219628"]
#print(test_df["type"].value_counts())

counts_df = pd.DataFrame() ## this counts number of gains and losses across all
for i in range(len(links)):
	counts_df = counts_df.append(results_df[links[i]]["type"].value_counts().rename(links[i]))
counts_df.to_csv("disruption_counts.txt",sep="\t")
#print(counts_df.sum(axis="columns"))

############## create simulated links
for i in range(len(links)):

	x = pybedtools.BedTool()
#l is length, n is number
# woohoo
	length = re.split("[-:]",links[i])
	length = int(length[3])-int(length[2])

	y = x.random(l=length,n=100,g="/Users/mike/replication_tcga/data/hg38.cleaned.bed")# could change this to the file i created "hg38.cleaned.bed"
	print(y.intersect("/Users/mike/replication_tcga/data/TCGA.segtabs.bed",wao=True))







### analyze results
### make data frame with just counts of loss,gain,disruption,null for all asars, with cancer type included???? could probably do this in the loop...


# lengths = pd.DataFrame(counts,columns=["asar","loss","gain","disruption","neutral"]).set_index("asar")
# lengths.hist()
# plt.show()

### results keys are asars


#print(cancer_atlas)
#print(asar.dtypes)
#print(asar.loc[2,"start"]+3)

########
# okay do boolean indexing of cancer atlas data
# easiest way is to just make a new column with all of the
# sample id where links are disrupted
#for index,row in asar.iterrows():
	#cancer_atlas[(cancer_atlas["start"] > row["start"])]
