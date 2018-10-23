"""
with overlapping windows of genome of variable size
calculate L1 density<- intersect bedtools, overlap, pandas to sum the overlap
SNP density<- ucsc snp and bedtools count, GC content <- bedtools nuc

use pybedtools to generate windows, use big pandas data frame to index GTF files or bed files of genomic fatures
probably need to use seqIO from biopython, seqrecord

samtools faidx reference.fasta
samtools faidx reference.fasta lyrata:1-108
"""

import os
import subprocess
import csv
import pybedtools


def make_windows(size,slide):
	command = ["bedtools makewindows -g /Users/mike/replication_tcga/data/chromosome_info_38.txt -w "
	+str(size)+
	" -s "+
	str(slide)]

	p = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
	result = p.communicate()[0].decode('ascii')
	result = [x.split("\t") for x in result.splitlines()] ## shorten to one loop
	result = [[x[0],
	int(x[1]),
	int(x[2])] for x in result]
	return result

# for i in range(len(links)):

#         x = pybedtools.BedTool()
# #l is length, n is number
# # woohoo
#         length = re.split("[-:]",links[i])
#         length = int(length[3])-int(length[2])

#         y = x.random(l=length,n=100,g="/Users/mike/replication_tcga/data/hg38.cleaned.bed")# could change this to the file i created "hg38.cleaned.bed"
#         print(y.intersect("/Users/mike/replication_tcga/data/TCGA.segtabs.bed",wao=True))

def write_bed(x,file_name):
	## write from list of list to tab sep bed file

	with open(file_name, 'w', newline='\n') as f:
		writer = csv.writer(f,delimiter="\t",quoting=csv.QUOTE_NONE)
		writer.writerows(x)
	f.close()

	return

def read_bed(file_name):

	with open (file_name) as g:
		bed = g.readlines()
		bed = [x.rstrip("\n").split("\t") for x in links]
		bed = [[str(x[0]),int(x[1]),int(x[2])]+x[3:] for x in links]
	### read bed file of any size, return list of list. unsure how to handle multiple types
	return bed

def get_gc_L1_snp(x):

	### given a bed file, add columns with % GC, snps per kb, %L1
	a = pybedtools.BedTool(x)
	a.nucleotide_content()
	return



# a = make_windows(1000000,100000)
# print(a)
