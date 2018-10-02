import os
import csv
import numpy as np
import pandas as pd
import sys

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
with open ("/Users/heskett/tcga_replication_timing/data/tcga_cancer_type_dictionary.txt") as f:
	lines = f.readlines()
	lines = (x.rstrip("\n").split("\t") for x in lines)
	cancer_atlas_dictionary = dict(lines)
	f.close()
with open ("/Users/heskett/tcga_replication_timing/data/links_segments_overlap.bed") as f:
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


links = list(np.unique([str(x[3])+":"+str(x[0])+":"+str(x[1])+"-"+str(x[2]) for x in intersections]))
samples = list(np.unique([x[8][:-3] for x in intersections]))
types = ["gain","loss","disruption","null"]
#print(cancer_atlas_dictionary.keys())
remove_samples = [x for x in samples if x not in cancer_atlas_dictionary.keys()]
normals = [x[8] for x in intersections if x[-3:]=="-10"]
#print(normals)
#print(remove_samples)
results = {k:{l : [] for l in types} for k in links } #[[sample_name,start,stop,cancertype]]
#print(results)

for i in range(len(intersections)):
	if intersections[i][8][:-3] in remove_samples:
		continue
	if (intersections[i][12] == intersections[i][4]) and (intersections[i][11]==1.0):
		results[unique_name(intersections[i])]["loss"]+=[[intersections[i][8],cancer_atlas_dictionary[intersections[i][8][:-3]],intersections[i][5],intersections[i][6],intersections[i][7]]]
	if (intersections[i][12] == intersections[i][4]) and (intersections[i][11]>2.0):
		results[unique_name(intersections[i])]["gain"]+=[[intersections[i][8],cancer_atlas_dictionary[intersections[i][8][:-3]],intersections[i][5],intersections[i][6],intersections[i][7]]]
	if (intersections[i][12] == intersections[i][4]) and (intersections[i][11]==2.0):
		results[unique_name(intersections[i])]["null"]+=[[intersections[i][8],cancer_atlas_dictionary[intersections[i][8][:-3]],intersections[i][5],intersections[i][6],intersections[i][7]]]
	if (intersections[i][12] < intersections[i][4]):
		results[unique_name(intersections[i])]["disruption"]+=[[intersections[i][8],cancer_atlas_dictionary[intersections[i][8][:-3]],intersections[i][5],intersections[i][6],intersections[i][7]]]
#print(results["VLINC246:6:85094762-85318808"]["loss"])

### analyze results
### make data frame with just counts of loss,gain,disruption,null for all asars, with cancer type included???? could probably do this in the loop...


print(len(results["VLINC246:6:85094762-85318808"]["loss"]))
print(len(results["VLINC246:6:85094762-85318808"]["gain"]))
print(len(results["VLINC246:6:85094762-85318808"]["disruption"]))
print(len(results["VLINC246:6:85094762-85318808"]["null"]))

test  = pd.DataFrame(columns=["asar","loss","gain","disruption","null","cancer_type"])

counts = []
for i in range(len(links)):
	counts += [[links[i],len(results[links[i]]["loss"]),
	len(results[links[i]]["gain"]),
	len(results[links[i]]["disruption"]),
	len(results[links[i]]["null"])]]

print(pd.DataFrame(counts,columns=["asar","loss","gain","disruption","null"]).set_index("asar"))

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
