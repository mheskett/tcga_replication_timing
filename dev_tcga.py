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

def plot_features(length,window_fraction=1):

	start=time.time()
	# snp doesnt need to be bedtool object yet i guess
	snp_file = "/Users/mike/replication_tcga/data/snp150.ucsc.hg38.nochr.sorted.bed"
	l1_file = pybedtools.BedTool("/Users/mike/replication_tcga/data/L1_ucsc_hg38_nochr.bed")

	# initiate bed tool and make windows
	# all file must be sorted sort -k1,1 -k2,2n
	a=pybedtools.BedTool()
	windows=a.window_maker(g="/Users/mike/replication_tcga/data/TCGA.nochr.bed",
							w=length,s=length*window_fraction)

	windows_snp = windows.intersect(snp_file,c=True,sorted=True) # makes one column
	#print(str(snpt-gc)+ " time for snp ") ## pre sorted file much faster sort -k1,1 -k2,2n
	
	## get L1
	windows_snp_l1 = windows_snp.coverage(l1_file)# produces 4 columns

	df = windows_snp_l1.to_dataframe(names=["chrom","start","end",
										"num_snps","cov1",
										"cov2","cov3","fraction_l1"])
	df["snps/kb"] = df["num_snps"] / (length/1000)
	
	b = pybedtools.BedTool()
	window_snp_l1_reduced = b.from_dataframe(df)
	window_snp_l1_nuc =  window_snp_l1_reduced.nucleotide_content(fi="/Users/mike/replication_tcga/data/hg38.nochr.fa")#\
							#.to_dataframe() ## 0th column has colnames. current colname is garbage

	# columns at this point are...
	# chrom start stop num snps, cov1, cov2, cov3, fraction l1, snps/kb, %AT, %GC, whatever...
	window_snp_l1_nuc = window_snp_l1_nuc.to_dataframe()
	window_snp_l1_nuc = window_snp_l1_nuc.loc[:,0:10]
	window_snp_l1_nuc.columns = ["chrom","start","end","num_snps","cov1","cov2","cov3","fraction_l1","snps/kb","AT","pct_gc"]
	window_snp_l1_nuc[["pct_gc","fraction_l1","snps/kb","start","end"]] = \
	window_snp_l1_nuc[["pct_gc","fraction_l1","snps/kb","start","end"]].apply(pd.to_numeric)
	window_snp_l1_nuc.plot(kind = 'bar')
	plt.show()

	return

def check_coverage(segments_df,link_chromosome,link_start,link_end):
	### check TCGA to see if a segment is covered. return a subset of patients where the gene is 100% covered?

	return 


def faster_simulate_links(length,window_fraction=0.25,snps=False,l1=False,gc=False,wiggle=0.1,minimum=40,maximum=100):
	"""
	returns a list of fake genes, from a length
	predicated on the fact that bedtools.nuc is the only slow step in this analysis
	need to: trim down reference genome to include regions that are not in TCGA
	then: do a check to see how many samples contain each fake gene....
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
	windows=a.window_maker(b="/Users/heskett/tcga_replication_timing/data/hg19/TCGA.nochr.bed",
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
	########
	## find out how many samples contain each fake link
    #######
	if len(window_snp_l1_nuc) < minimum:
		wiggle+=0.05
		window_snp_l1_nuc = faster_simulate_links(length=length,snps=snps,gc=gc,l1=l1,wiggle=wiggle,minimum=minimum) # recursive biatch
		print(wiggle, " wiggle") # not returning this!!!!
	if len(window_snp_l1_nuc) > maximum:
		window_snp_l1_nuc = window_snp_l1_nuc.sample(n=100,axis="rows",random_state=1) # shuld randomly sample rows!!
	return window_snp_l1_nuc.loc[:,["chrom","start","end",
					"snps/kb","pct_gc","fraction_l1"]].astype({"chrom":str,"start":int,"end":int,
						"snps/kb":float,"pct_gc":float,"fraction_l1":float})

def unique_name(x):
	# given a line of file, create unique name for link
	return str(x[3])+":"+str(x[0])+":"+str(x[1])+"-"+str(x[2])

def is_real_link(x):

	b = pybedtools.BedTool("/Users/heskett/tcga_replication_timing/data/hg19/links.annotated.hg19.bed")
	a = x # returns only non-overlappers with real links
	# x is a bed tool object...

	return a.intersect(b,v=True)

def search_tcga(segments_df,link_chromosome,link_start,link_end):
	types = ["gain","loss","neutral","disruption"]
	segments_df=segments_df[segments_df['chr']==link_chromosome] # then remove the additional operation from below...
	print("dododo")
	#import pdb;pdb.set_trace()
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
	disruptions = segments_df[(segments_df['start'].between(left=link_start, right=link_end)
								| segments_df['stop'].between(left=link_start, right=link_end))]
	disruptions = disruptions[disruptions.duplicated(subset="patient", keep=False)] # mcaarks all dupes as true cause we want to keep dupes
	disruptions = disruptions.drop_duplicates(subset="patient", keep="first") # keeps the first dupe
	
	losses.loc[:,"type"] = pd.Series(["loss"]*len(losses.index), index=losses.index).astype(str) # empty list was being assigned float64 type...
	gains.loc[:,"type"] = pd.Series(["gain"]*len(gains.index), index=gains.index).astype(str)
	disruptions.loc[:,"type"]= pd.Series(["disruption"]*len(disruptions.index), index=disruptions.index).astype(str)
	neutrals.loc[:,"type"] = pd.Series(["neutral"]*len(neutrals.index), index=neutrals.index).astype(str)
	print(pd.concat([losses,gains,disruptions,neutrals]))
	return pd.concat([losses,gains,disruptions,neutrals])

def cancer_specific(segments_df,links,cancer_types):

	results = { } # each k is a real link
	# output=[["link_name","cancer_type","p_value_gain","p_value_loss","p_value_disruption","p_value_neutral"]]
	output=[]

	for i in range(len(links)):
		print(links[i][0]," ",i)

		# create simulated links
		fake_links = faster_simulate_links(length=links[i][5],
										snps=links[i][6],
										l1=links[i][8],
										gc=links[i][7]
											)
		# get results for the real link
		results[links[i][0]] =\
				search_tcga(segments_df=segments_df,
						link_chromosome=links[i][1],
						link_start=links[i][2],
						link_end=links[i][3],
						)

		# get results for the simulated links
		fake_link_names = ["fake_link_"+str(row['chrom'])+":"+str(row["start"])+"-"+str(row["end"]) for index,row in fake_links.iterrows()]
		fake_results = {} 
		print("counting all disruptions in simulated links")

		for index,row in fake_links.iterrows():
			name = "fake_link_"+str(row['chrom'])+":"+str(row["start"])+"-"+str(row["end"])
			fake_results[name] = \
				search_tcga(segments_df=segments_df,
							link_chromosome=row['chrom'],
							link_start=row['start'],
							link_end=row['end'],
							) # these are named here
			# print(fake_results[name])
		print("counting cancer specific disruptions in real and simulated links")
		### Needs to have coverage = ( gain + loss + disruption + neutral ) / total samples per in a cancer type
		for k in range(len(cancer_types)):	## Get rid of low coverage cancer types here

			counts = {}
			df = results[links[i][0]]

			counts[links[i][0]] = [len(df[(df["cancer_type"]==cancer_types[k]) & (df["type"]=="gain")]) , ## just divide here by total n cancer type
									len(df[(df["cancer_type"]==cancer_types[k]) & (df["type"]=="loss")]),
									len(df[(df["cancer_type"]==cancer_types[k]) & (df["type"]=="disruption")]),
									len(df[(df["cancer_type"]==cancer_types[k]) & (df["type"]=="neutral")])]
			fake_counts = {} 
			for j in range(len(fake_link_names)):
				fake_df = fake_results[fake_link_names[j]]

				fake_counts[fake_link_names[j]] = [len(fake_df[(fake_df["cancer_type"]==cancer_types[k]) & (fake_df["type"]=="gain")]),
									len(fake_df[(fake_df["cancer_type"]==cancer_types[k]) & (fake_df["type"]=="loss")]),
									len(fake_df[(fake_df["cancer_type"]==cancer_types[k]) & (fake_df["type"]=="disruption")]),
									len(fake_df[(fake_df["cancer_type"]==cancer_types[k]) & (fake_df["type"]=="neutral")])]

			# get distribution of disruptions in fake links
			counts_df = pd.DataFrame.from_dict(counts,orient="index",columns=["gain","loss","disruption","neutral"])
			fake_counts_df = pd.DataFrame.from_dict(fake_counts,orient="index",columns=["gain","loss","disruption","neutral"])

			pvals = [links[i][0],
						cancer_types[k], # type weak means <=, type strict means only values strictly less are counted
						1-scipy.stats.percentileofscore(fake_counts_df["gain"], score=counts_df["gain"].values[0], kind='strict')/100, 
						1-scipy.stats.percentileofscore(fake_counts_df["loss"], score=counts_df["loss"].values[0], kind='strict')/100,
						1-scipy.stats.percentileofscore(fake_counts_df["disruption"], score=counts_df["disruption"].values[0], kind='strict')/100,
						1-scipy.stats.percentileofscore(fake_counts_df["neutral"], score=counts_df["neutral"].values[0], kind='strict')/100]
			pvals = [1/len(fake_links) if x==0 else x for x in pvals]
			output += [pvals]

	return output
### open and process files

with open ("/Users/heskett/tcga_replication_timing/data/hg19/tcga_cancer_type_dictionary.txt") as f:
	lines = f.readlines()
	lines = (x.rstrip("\n").split("\t") for x in lines)
	cancer_atlas_dictionary = dict(lines)
	f.close()

with open ("/Users/heskett/tcga_replication_timing/data/hg19/links.annotated.hg19.bed") as g:
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
with open ("/Users/heskett/tcga_replication_timing/data/hg19/TCGA.segtabs.final.merged.sorted.no23.bed") as h:
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

###
exclude_list = ["ACC","CHOL","DLBC","KICH","MESO","UCS","UVM"]
segments_df= pd.DataFrame(segments)#.astype({0: str})
## get rid of all segments belonging to low coverage tumor types
#### 
segments_df.columns = ["chr","start","stop","patient","cancer_type","copy_number"]
segments_df = segments_df[~segments_df.cancer_type.isin(exclude_list)]
exclude_list = ["ACC","CHOL","DLBC","KICH","MESO","UCS","UVM"]

patients_per_type = segments_df.groupby("cancer_type")["patient"].nunique()
total_patients = sum(patients_per_type)
types = ["gain","loss","neutral","disruption"]
cancer_types = list(np.unique([x for x in cancer_atlas_dictionary.values() if x not in exclude_list]))
links_names = [x[0] for x in links]


#### main program

pool = mp.Pool()
results = pool.starmap(cancer_specific,[(segments_df,[x],cancer_types) for x in links]) ##  each link needs to be list of list in parallel call

with open("links_tcga_parallel.txt", "w") as f: # should make this write each loop
    writer = csv.writer(f,delimiter="\t") #format output better
    for i in range(len(results)):
    	for j in range(len(results[i])):
    		writer.writerow(results[i][j])


# os.system("tr '    ' '\n' < links_tcga_parallel1.txt > tcga_results.txt")

# os.system(sed "s/\[//g" tcga_results.txt | sed "s/\]//g" | sed "s/'//g" | sed "s/,/      /g" > tcga_results_formatted.tsv)
