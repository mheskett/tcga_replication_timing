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

with open ("/Users/heskett/tcga_replication_timing/data/hg19/TCGA.segtabs.final.merged.sorted.no23.bed") as h:
	segments = h.readlines()
	segments = [x.rstrip("\n").split("\t") for x in segments]
	h.close() # fast

with open ("/Users/heskett/tcga_replication_timing/data/tcga_cancer_type_dictionary.txt") as f:
	lines = f.readlines()
	lines = (x.rstrip("\n").split("\t") for x in lines)
	cancer_atlas_dictionary = dict(lines)
	f.close()

segments = [[str(x[0]),
	int(x[1]),
	int(x[2]),
	str(x[3]),
	str(cancer_atlas_dictionary[str(x[3][:-3])]),
	float(x[4])] for x in segments if x[3][:-3] in cancer_atlas_dictionary.keys()]
#print(segments)


def make_random_windows(length,number):
	# all file must be sorted sort -k1,1 -k2,2n
	# need to pick 
	a=pybedtools.BedTool()
	windows=a.random(g="/Users/heskett/tcga_replication_timing/data/hg19/TCGA.hg19.chromosome.file.bed",
							l=length,n=number)
	tcga = pybedtools.BedTool()
	tcga = tcga.from_dataframe(pd.DataFrame(segments))
	#df = windows.to_dataframe().loc[:,["chrom","start","end"]]
	
	# c = b.coverage(a) from bedtools doc. But coverage doesn't give a list of which features overlapped...
	return windows.intersect(tcga,wa=True,wb=True)

# #### main program
# df = pd.DataFrame()
# for i in range(10**5,3*10**5,5*10**3):
# 	df = pd.concat([df,make_random_windows(length=i,number=10)],ignore_index=True)
# #print(df)

print(make_random_windows(length=10**6,number=10))