import retikulos as rt
#np=rt.np
#cp=rt.cp
founder=rt.founder_miner(0.2)
randmuts=rt.randomMutations(founder[0][1].size,0.0002)
high_change=rt.codPos(randmuts,founder[0][1].shape)
a,b,muts_arr=rt.mutate_genome(founder[0][1],founder[0][2],high_change)
grn,threshs,decs,genes_on=founder[0][3],founder[0][4],founder[0][5],founder[0][8]