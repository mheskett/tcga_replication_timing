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
def min_max(x):
	return min(max(0,x),1)

def is_real_link(x):

	b = pybedtools.BedTool("/Users/heskett/tcga_replication_timing/data/hg19/links.annotated.hg19.bed")
	a = x # returns only non-overlappers with real links
	# x is a bed tool object...

	return a.intersect(b,v=True)

def faster_simulate_links(length,window_fraction=0.25,snps=False,l1=False,gc=False,wiggle=0.1,minimum=40,maximum=100):
	"""
	returns a list of fake genes, from a length
	predicated on the fact that bedtools.nuc is the only slow step in this analysis
	need to: trim down reference genome to include regions that are not in TCGA
	then: do a check to see how many samples contain each fake gene....
	Need a file that has all non coding regions of genome.
	"""
	if length < 20000:
		window_fraction = 1

	start=time.time()
	# snp doesnt need to be bedtool object yet i guess
	snp_file = "/Users/heskett/tcga_replication_timing/data/hg19/snp150.ucsc.hg19.nochr.bed"
	l1_file = pybedtools.BedTool("/Users/heskett/tcga_replication_timing/data/hg19/l1.ucsc.hg19.nochr.bed")

	# initiate bed tool and make windows
	# all file must be sorted sort -k1,1 -k2,2n
	a=pybedtools.BedTool()
	# "/Users/heskett/tcga_replication_timing/data/hg19/TCGA.nochr.bed"
	windows=a.window_maker(b="/Users/heskett/tcga_replication_timing/data/hg19/ucsc.not.i.e.tcga.covered.bed",
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
										"cov2","cov3","fraction_l1"]) # check that crom is str or int?
	df["snps/kb"] = df["num_snps"] / (length/1000)

	snps_score,l1_score = (scipy.stats.percentileofscore(a=df["snps/kb"], score=snps, kind='rank'), # this giving seg fault?
							scipy.stats.percentileofscore(a=df["fraction_l1"], score=l1, kind='rank'))
	df = df[
		(df["snps/kb"]\
			.between(left=df["snps/kb"].quantile(q=min_max(snps_score/100-wiggle)), # could separate these for optimization if necessary
					right=df["snps/kb"].quantile(q=min_max(snps_score/100+wiggle))))\
		& (df["fraction_l1"]\
			.between(left=df["fraction_l1"].quantile(q=min_max(l1_score/100-wiggle)), 
					right=df["fraction_l1"].quantile(q=min_max(l1_score/100+wiggle))))]
		## should add another recursive call right here...
	b = pybedtools.BedTool()
	window_snp_l1_reduced = b.from_dataframe(df)
	window_snp_l1_nuc =  window_snp_l1_reduced.nucleotide_content(fi="/Users/heskett/tcga_replication_timing/data/hg19/human_g1k_v37_nochr.fasta")#\
							#.to_dataframe() ## 0th column has colnames. current colname is garbage

	# columns at this point are...
	# chrom start stop num snps, cov1, cov2, cov3, fraction l1, snps/kb, %AT, %GC, whatever...
	window_snp_l1_nuc = is_real_link(window_snp_l1_nuc)
	window_snp_l1_nuc = window_snp_l1_nuc.to_dataframe()
	window_snp_l1_nuc = window_snp_l1_nuc.loc[:,0:10]
	nuct=time.time()
	print(nuct-snpt," nuc time")
	window_snp_l1_nuc.columns = ["chrom","start","end","num_snps","cov1",
								"cov2","cov3","fraction_l1","snps/kb","AT","pct_gc"]

	window_snp_l1_nuc[["pct_gc","fraction_l1","snps/kb","start","end"]] = \
	window_snp_l1_nuc[["pct_gc","fraction_l1","snps/kb","start","end"]].apply(pd.to_numeric)
	gc_score = scipy.stats.percentileofscore(window_snp_l1_nuc["pct_gc"], score=gc, kind='rank')
	window_snp_l1_nuc = window_snp_l1_nuc[window_snp_l1_nuc["pct_gc"]\
			.between(left=window_snp_l1_nuc["pct_gc"].quantile(q=min_max(gc_score/100-wiggle)),
					right=window_snp_l1_nuc["pct_gc"].quantile(q=min_max(gc_score/100+wiggle)))]
	print(len(window_snp_l1_nuc))
	########
	## find out how many samples contain each fake link
    #######
	if len(window_snp_l1_nuc) < minimum:
		wiggle+=0.05
		window_snp_l1_nuc = faster_simulate_links(length=length,snps=snps,gc=gc,l1=l1,wiggle=wiggle,minimum=minimum) # recursive biatch
		print(wiggle, " wiggle") # not returning this!!!!
	if len(window_snp_l1_nuc) > maximum:
		window_snp_l1_nuc = window_snp_l1_nuc.sample(n=maximum,axis="rows",random_state=1) # shuld randomly sample rows!!
	final = window_snp_l1_nuc.loc[:,["chrom","start","end",
					"snps/kb","pct_gc","fraction_l1"]].astype({"chrom":str,"start":int,"end":int,
						"snps/kb":float,"pct_gc":float,"fraction_l1":float})

	final.loc[:,["chrom","start","end"]].to_csv("simulate.asars."+ str(l1)+ "."+ str(length)+"."+str(snps)+"."+str(gc)+".bed", 
		index=None, header=None, sep="\t")
	print(final)
	return final

files=[]
for i in range(5):
	l1=.25 + i*0.05 
	length=200000
	gc=0.35
	snps=3.9
	faster_simulate_links(length,window_fraction=0.25,snps=snps,l1=l1,gc=gc,wiggle=0.1,minimum=40,maximum=200)
	files+=["simulate.asars."+ str(l1)+ "."+ str(length)+"."+str(snps)+"."+str(gc)+".bed"]
	file_string = "igv " + " ".join(files)

os.system(file_string)
