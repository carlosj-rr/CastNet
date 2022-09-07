import numpy as np
import copy as cp
import retikulos as rt
import params_file as pf
import time
reps=100
founder1=rt.founder_miner()
founder2=rt.founder_miner(0.2)
allnums=[]
for i in range(reps):
    print(f"{i}/{reps}")
    start=time.time()
    pop=rt.grow_pop(founder1,pf.pop_size) # Command to time.
    stop=time.time()
    allnums.append(stop-start)

allnums=np.array(allnums)
print(f"Average time it took to run function {reps} times: {np.sum(allnums)/100}, max val: {np.max(allnums)}, min: {np.min(allnums)}.")