import os
import csv
import numpy as np
import pandas as pd
import sys





########## Load necessary files
platform_name=sys.platform

if platform_name=="linux":
	working_directory =  "/home/exacloud/lustre1/SpellmanLab/heskett/tcga_replication_timing/data/"
if platform_name=="darwin":
    working_directory = "/Users/mike/replication_tcga/data/"

cancer_atlas_path = working_directory + "TCGA_mastercalls.abs_segtabs.fixed.txt"
asar_path = working_directory + "all_links.txt"
cancer_atlas_dictionary_path = working_directory + "tcga_cancer_type_dictionary.txt"

cancer_atlas = pd.read_table(cancer_atlas_path,header=0)
asar = pd.read_table(asar_path,names=["chr","start","end","annotation","length"],
		dtype={"chr":str,"start":int,"end":int,"annotation":str,"length":float})
#asar.columns=["chr","start","end","annotation","length"]
	#dtype={"chr":str,"start":int,"end":int,"annotation":str,"length":float})
with open (cancer_atlas_dictionary_path) as f:
	lines = f.readlines()
	lines = (x.rstrip("\n").split("\t") for x in lines)
	cancer_atlas_dictionary = dict(lines)

#print(cancer_atlas)
#print(asar.dtypes)
#print(asar.loc[2,"start"]+3)

########
# okay do boolean indexing of cancer atlas data
# easiest way is to just make a new column with all of the
# sample id where links are disrupted
#for index,row in asar.iterrows():
	#cancer_atlas[(cancer_atlas["start"] > row["start"])]
