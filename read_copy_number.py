"""
read gistic
"""
import os
import csv
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import pybedtools
import re


# table = pd.read_table("/Users/mike/replication_tcga/data/copy_number\
# gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/\
# table_amp.conf_99.txt",usecols=[1,4,5])
# table.columns=["chr","start","stop"]
name = re.split("[_-]",sys.argv[1])[3]+"_"+re.split("[/]",sys.argv[1])[-1]
name = name[:-4]+".bed"
with open (sys.argv[1]) as f:
	lines = f.readlines()[1:]
	lines = (x.rstrip("\n").split("\t") for x in lines)
	lines = [[x[1],x[4],x[5]] for x in lines]
	f.close()
pd.DataFrame(lines).to_csv("/Users/mike/replication_tcga/data/copy_number/"+name,sep="\t",index=False,header=False)