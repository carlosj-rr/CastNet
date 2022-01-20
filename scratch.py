import retikulos as rt
import numpy as np
from scipy.stats.distributions import chi2
founder=rt.founder_miner(0.75)
rez=rt.main(founder)
def parse_rez(rez):
     expvals=[2.49,4.98,4.98,4.98,4.98,9.96,9.96,9.96,9.96]
     print(f"Population Expected avg differences Observed avg differences Chi-square p-value")
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
         pval=1-chi2.cdf(chisq,df)
         meandiffs=np.mean(outvect)
         print(f"{k} {expvals[k]} {meandiffs} {pval}")
     return
parse_rez(rez)