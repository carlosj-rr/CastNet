import retikulos as rt
import numpy as np
founder=rt.founder_miner(0.75)
rez=rt.main(founder)
def parse_rez(rez):
     expvals=[2.49,4.98,4.98,4.98,4.98,9.96,9.96,9.96,9.96]
     for k in range(len(rez)):
         pop=rez[k]
         outvect=np.ndarray(pop.size,dtype=object)
         exp=expvals[k]
         chisqvect=np.ndarray(pop.size,dtype=object)
         df=pop.size-1
         dftab=np.ndarray(len(rez),dtype=int)
         dftab[k]=df
         chivaltab=np.ndarray(len(rez))
         for i in range(pop.size):
             obs=np.where(founder[0][1] != pop[i][1])[0].size
             outvect[i]=obs
             chisqvect[i]=((obs-exp)**2)/exp
         chisq=sum(chisqvect)
         chivaltab[k]=chisq
         meandiffs=np.mean(outvect)
         print(f"Population: {k}\nExpected avg differences: {expvals[k]}\nObserved avg differences: {meandiffs}\nChi-square p-value: {pval}")
     return
parse_rez(rez)