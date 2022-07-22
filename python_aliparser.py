import numpy as np
import pysam as ps
import os
import tarfile
import sys
import pickle

infile=sys.argv[1]
line_and_num=infile.split(".")[0]

file=tarfile.open(infile)
file.extractall(".")
ali_list=os.listdir(line_and_num)
perms=np.loadtxt("perms_orders_0idx.csv",delimiter=",",dtype=int)

out_arr=out_arr=np.ndarray((len(ali_list)*perms.shape[0],256))
trearr=np.repeat(int("Far_" in ali_list[0]),len(ali_list)*perms.shape[0])

charlist=['A','C','T','G']
numchars=len(charlist)

sitepatterns=np.ndarray(numchars**4,dtype=object)
idx=0
for i in range(numchars):
	for j in range(numchars):
		for k in range(numchars):
			for t in range(numchars):
				x=np.array([charlist[i],charlist[j],charlist[k],charlist[t]])
				sitepatterns[idx]=''.join(x)
				idx+=1

def site_props_calc(obs_sites,ref_sites):
    pattern_count=np.ndarray(len(ref_sites))
    for i in range(len(ref_sites)):
	    seq_length=len(obs_sites)
	    count=len(np.where(obs_sites == ref_sites[i])[0])
	    pattern_count[i]=count/seq_length
    return(pattern_count)

rowtoadd=0
for idx in np.arange(len(ali_list)):
    seq_obj=ps.FastaFile(line_and_num+"/"+ali_list[idx])
    for order in range(perms.shape[0]):
        curr_arr=np.ndarray(seq_obj.lengths[0],dtype=object)
        for site_idx in range(curr_arr.size):
            site_pattern=np.ndarray(4,dtype=object)
            for seq_num in perms[order]:
                    site_pattern[seq_num]=seq_obj.fetch(seq_obj.references[seq_num])[site_idx]
            curr_arr[site_idx]=''.join(site_pattern)
        row_props=site_props_calc(curr_arr,sitepatterns)
        print(f"Adding to row number {rowtoadd} of 120000.")
        out_arr[rowtoadd]=row_props
        rowtoadd+=1

outtuple=(out_arr,trearr)
outfilename=line_and_num+"_TT.pkl"
with open(outfilename, "wb") as fh:
        pickle.dump(outtuple, fh)

file.close()