import retikulos as rt
np=rt.np
cp=rt.cp
founder=rt.founder_miner(0.75)
founder4x=np.ndarray((4,),dtype=object)
for i in range(4):
    founder4x[i]=cp.deepcopy(founder[0])
