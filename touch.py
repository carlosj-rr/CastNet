import numpy as np
gen10=np.load("Generation_10_branch_1.npy",allow_pickle=True)
gen20=np.load("Generation_20_branch_1.npy",allow_pickle=True)

print(f"Genomes are {np.multiply(100,np.sum(gen10[0][1] != gen20[0][1])/gen[0][1].size)}% different")