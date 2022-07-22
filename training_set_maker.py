import numpy as np
import sys
import pickle
import os

filenames=np.loadtxt(sys.argv[1],delimiter=",",dtype=str)

for row in np.arange(filenames.shape[0]):
    num=filenames[row][0].split("_")[1].split(".")[0]
    matfile=filenames[row][0]
    mat=np.loadtxt(matfile,delimiter=",")
    trefile=filenames[row][1]
    tres=np.loadtxt(trefile,dtype=int)
    traintup=(mat,tres)
    with open(num+"_TT.pkl","wb") as fh:
        pickle.dump(traintup,fh)
    print(f"Wrote the data from line {num} to disk. Removing the source files...")
    os.remove(matfile)
    os.remove(trefile)
