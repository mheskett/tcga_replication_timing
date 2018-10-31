import pandas as pd
import numpy as np
# with open ("/Users/mike/replication_tcga/data/TCGA.segtabs.bed") as h:
# 	segments = h.readlines()
# 	segments = [x.rstrip("\n").split("\t") for x in segments]
# 	segments = [[str(x[0]),
# 	int(x[1]),
# 	int(x[2]),
# 	str(x[3]),
# 	float(x[4]),
# 	float(x[5]),
# 	float(x[6])] for x in segments]
# 	h.close() # fast

df = pd.read_table("/Users/mike/replication_tcga/data/TCGA.segtabs.26.fixed.bed",header=None)
df = df.sort_values(by=[3,0,1])

df_m = df.as_matrix()
print(df_m[2][3])

seg_cn=0
seg_chromosome=""
seg_sample=""
seg_start=0
seg_stop=0

# for i in range(len(df_m)):
results = []
count=0
while count < len(df_m):
	seg_start=df_m[count][1]
	seg_stop=df_m[count][2]
	seg_cn=df_m[count][6]
	seg_chromosome=df_m[count][0]
	seg_sample=df_m[count][3]

	while (count < len(df_m)) and (df_m[count][0]==seg_chromosome) and \
		(df_m[count][3]==seg_sample) and \
		(df_m[count][6]==seg_cn):
			count+=1
	results +=[[int(seg_chromosome),seg_start,df_m[count-1][2],seg_sample,seg_cn]]
		 # check cn, chr, samp of next
pd.DataFrame(results).to_csv("TCGA.segtabs.26.fixed.merged.bed",sep="\t",header=None,index=None)