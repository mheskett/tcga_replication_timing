import matplotlib.pyplot as plt
import matplotlib.patches 
# set up the figure
centromere = {"1":"124535434", # hg19
"2":"95326171",
"3":"93504854",
"4":"52660117",
"5":"49405641",
"6":"61830166",
"7":"61054331",
"8":"46838887",
"9":"50367679",
"X":"61632012",
"Y":"13104553",
"10":"42254935",
"11":"54644205",
"12":"37856694",
"13":"19000000",
"14":"19000000",
"15":"20000000",
"16":"38335801",
"17":"25263006",
"18":"18460898",
"19":"27681782",
"20":"29369569",
"21":"14288129",
"22":"16000000"}

chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12",
    "13","14","15","16","17","18","19","20","21","22","X"]

stop_start = {"1":249250621,
"2":243199373,
"3":198022430,
"4":191154276,
"5":180915260,
"6":171115067,
"7":159138663,
"X":155270560,
"8":146364022,
"9":141213431,
"10":135534747,
"11":135006516,
"12":133851895,
"13":115169878,
"14":107349540,
"15":102531392,
"16":90354753,
"17":81195210,
"18":78077248,
"20":63025520,
"Y":59373566,
"19":59128983,
"22":51304566,
"21":48129895}

bedfile = "test.asar.bed"

with open(bedfile) as asars:
	lines = asars.readlines()
	lines = [x.rstrip("\n").split("\t") for x in lines] #or lines[1:] if theres a header
	asar_list = [[str(x[0]),int(x[1]),int(x[2])] for x in lines]

#asar_list = [["9",10000000,20000000],["2",3000000,7000000],["6",30000000,32000000]]
f,ax = plt.subplots(1,len(chromosomes),sharex=False,sharey=False,figsize=(14,1))
f.subplots_adjust(hspace=0)

for i in range(len(chromosomes)):
	ax[i].axhline(y=0,color="black",linestyle="-")

	#ax[i].scatter(,s=1,lw=0.05,
	          #  label=samples[j],edgecolor="black",alpha=0.6)
	#ax[i].set(xlabel=chromosomes[i]) # x axis labels or no
	ax[i].set_yticks([])
	ax[i].set_xticks([])
	ax[i].set_ylim([-10,10])
	ax[i].set_xlim([0,stop_start[chromosomes[i]]])
	ax[i].axvline(x=int(centromere[chromosomes[i]]), linestyle = "--", lw = 0.5,color="blue")
	### to plot lines over the asars
	# ax2= ax[i].twinx() 
	# x=range(0,stop_start[chromosomes[i]],100000)
	# y = [val / stop_start[chromosomes[i]] for val in x ] 
	# ax2.plot(x,y)
	# ax2.set_xticks([])
	# ax2.set_yticks([])

	
	for j in range(len(asar_list)):
		if asar_list[j][0]==chromosomes[i]:
			asar_start = asar_list[j][1]
			asar_stop = asar_list[j][2]
			ax[i].add_patch(matplotlib.patches.Rectangle((asar_start,-2),width=asar_stop - asar_start,
			height=4,color="red"))
	ax[i].set(xlabel=chromosomes[i]) # x axis labels or no



f.subplots_adjust(wspace=0, hspace=0)
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size':4}, markerscale = 3)

plt.show()