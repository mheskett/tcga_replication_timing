import os
import csv
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import pybedtools
import re
import pybedtools


def unique_name(x):
	# given a line of file, create unique name for link
	return str(x[3])+":"+str(x[0])+":"+str(x[1])+"-"+str(x[2])



with open ("/Users/mike/replication_tcga/data/tcga_cancer_type_dictionary.txt") as f:
	lines = f.readlines()
	lines = (x.rstrip("\n").split("\t") for x in lines)
	cancer_atlas_dictionary = dict(lines)
	f.close()

with open ("/Users/mike/replication_tcga/data/links_final_grch38.bed") as g:
	links = g.readlines()
	links = [x.rstrip("\n").split("\t") for x in links]
	links = [[str(x[3])+":"+str(x[0])+":"+str(x[1])+"-"+str(x[2]),
	str(x[0]),
	int(x[1]),
	int(x[2]),
	str(x[3]),
	int(x[4])] for x in links]
	g.close()

with open ("/Users/mike/replication_tcga/data/TCGA.segtabs.bed") as h:
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


#print(segments)
print("filtering")
segments = [[str(x[0]),
	int(x[1]),
	int(x[2]),
	str(x[3]),
	str(cancer_atlas_dictionary[str(x[3][:-3])]),
	float(x[6])] for x in segments if x[3][:-3] in cancer_atlas_dictionary.keys()] # v slow...2min
print(len(segments))

segments_df= pd.DataFrame(segments)
segments_df.columns = ["chr","start","stop","patient","cancer_type","copy_number"]
print(segments_df.dtypes)

## do this instead of looping below...?
# create 4 DF's per link
types = ["gain","loss","neutral","disruption"]
links_names = [x[0] for x in links]

results = {k:{i:pd.DataFrame() for i in types} for k in links_names }
# print(segments_df[segments_df["copy_number"]==2.0])
# print(segments_df[segments_df["chr"]=="3"])
# print(segments_df[segments_df["copy_number"]==2.0])
print("Starting loop")
for i in range(len(links)):
	link_start = links[i][2]
	link_end = links[i][3]
	link_chromosome = links[i][1]
	name = links[i][0]

	results[name]["loss"] = segments_df[(segments_df['chr'] == link_chromosome) 
	& (segments_df['start'] <= link_start) & (segments_df['stop'] >= link_end) 
	& (segments_df['copy_number'] <2.0 )]
	
	results[name]["gain"] = segments_df[(segments_df['chr'] == link_chromosome) 
	& (segments_df['start'] <= link_start) & (segments_df['stop'] >= link_end) 
	& (segments_df['copy_number'] > 2.0 )]
	
	results[name]["neutral"] = segments_df[(segments_df['chr'] == link_chromosome) 
	& (segments_df['start'] <= link_start) & (segments_df['stop'] >= link_end) 
	& (segments_df['copy_number'] == 2.0 )]

	# if either seg end points are within an asar.....
	# use series.between....... may have to do and and or and and 
	results[name]["disruption"] = segments_df[(segments_df['chr'] == link_chromosome) 
	& (segments_df['start'].between(left=link_start,right=link_end) 
	| segments_df['stop'].between(left=link_start,right=link_end))
	& (segments_df['copy_number'] != 2.0 )] ## should be max one per link per sample--so need to process this further

print("done")
print(results)
# segments_df.loc
# (segments_df[2]-segments_df[1]).hist()
# print((segments_df[2]-segments_df[1]).describe())
# plt.show()

links_names = [x[0] for x in links]
cancer_types = list(np.unique([x for x in cancer_atlas_dictionary.values()]))
samples = list(np.unique([x[3][:-3] for x in segments]))

#print(links)

#results = {k:[] for k in links_names}

# print("looping")
# for i in range(len(links)):
# 	for j in range(len(segments)):
# 		if (links[i][1]==segments[j][0]) and (links[i][2]>=segments[j][1]) and (links[i][3]<=segments[j][2]):
# 			if segments[j][6] == 2.0:
# 				print("neutral")
# 				results[links[i][0]]+=[ [segments[j][0], segments[j][1], segments[j][2], segments[j][3], cancer_atlas_dictionary[segments[j][3][:-3]], segments[j][6]] ]
# 				continue
# 			if segments[j][6] < 2.0:
# 				print("loss")
# 				continue
# 			if segments[j][6] > 2.0:
# 				print("gain")
# 				continue
# 				# 
# 		if (links[i][1]==segments[j][0]) and (links[i][2]<segments[j][1]) and (links[i][3]>segments[j][1]):
# 			print("disruption?")
# 			segments[j]
# 			#print(links[i])
# 			#print(segments[j]) ## add to a list, then check if list has more than one copy state
# 			continue
# 		if (links[i][1]==segments[j][0]) and (links[i][2]<segments[j][2]) and (links[i][3]>segments[j][2]):
# 			print("disruption?") ## add to a list, then check if list has more than one copy state, if yes its a disruption
# 			continue
# 		# if (links[i][1]==segments[j][0]) and (links[i][2]<segments[j][1]) and (links[i][3]>segments[j][2]):
# 		# 	print("disruption?") ## add to a list, then check if list has more than one copy state, if yes its a disruption
# 		# 	continue
# #print(results)
results[links[9][0]]["loss"]["cancer_type"].value_counts().plot(kind='bar')
#print(len(results[links[9][0]]["disruption"]))
#print(links[0])
#plt.show(\]

counts = {k:[] for k in links_names} ### to make counts table...
for i in range(len(links)):
	counts[links[i][0]] += [len(results[links[i][0]]["gain"]),
	len(results[links[i][0]]["loss"]),
	len(results[links[i][0]]["disruption"]),
	len(results[links[i][0]]["neutral"])]
#print(results[links[0][0]]["loss"])
pd.DataFrame.from_dict(counts,orient="index",columns=["gain","loss","disruption","neutral"]).to_csv("tcga_counts.txt",sep="\t")
# of gains losses and disruptions per each

## join all loss tables, gain tables, then do histogram by cancer type...
#print(results[links[0][0]]["loss"])
answer = pd.DataFrame()
for i in range(len(types)):
	for j in range(len(links)):
		answer = answer.append(results[links[j][0]][types[i]])
	fig = answer["cancer_type"].value_counts().plot(kind='bar')
	plt.savefig(str(types[i])+".png")
	plt.close()
	print(answer)
	answer = pd.DataFrame()

## simulated links

for i in range(len(links)):

	x = pybedtools.BedTool()
#l is length, n is number
	length = re.split("[-:]",links[i])
	length = int(length[3])-int(length[2])

	y = x.random(l=length,n=100,g="/Users/mike/replication_tcga/data/hg38.cleaned.bed")
# results_df = {} # dictionary of data frames for each link
# for i in range(len(links_names)):
# 	tmp=pd.DataFrame.from_dict(results[links_names[i]])
# 	tmp.columns = ["chr","start","stop","sample","cancer_type","copy_number"]
# 	results_df[links_names[i]]=tmp
