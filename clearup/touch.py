import numpy as np
gen10=np.load("Generation_10_branch_1.npy",allow_pickle=True)
gen20=np.load("Generation_20_branch_1.npy",allow_pickle=True)

def calc_avg_perc_genome_diff(pop1,pop2):
    outarr=np.ndarray(len(pop1)*len(pop2))
    counter=0
    for i in range(len(pop1)):
        first_genome=pop1[i][1]
        for j in range(len(pop2)):
            second_genome=pop2[j][1]
            outarr[counter]=np.multiply(np.sum(first_genome != second_genome)/first_genome.size,100)
            counter += 1
    return np.mean(outarr)

def zeroes_in_grn_stds(pop):
    call=np.any(np.std(np.array([x[3].flatten() for x in pop[:]]),axis=0) == 0)
    if call:
        count=np.sum(np.std(np.array([x[3].flatten() for x in gen10[:]]),axis=0) == 0)
        return(call,count)
    else:
        return(call)

def calc_avg_perc_grn_diff(pop1,pop2):
    outarr=np.ndarray(len(pop1)*len(pop2))
    counter=0
    for i in range(len(pop1)):
        first_grn=pop1[i][3]
        for j in range(len(pop2)):
            second_grn=pop2[j][3]
            outarr[counter]=np.multiply(np.sum(first_grn != second_grn)/first_grn.size,100)
            counter += 1
    return np.mean(outarr)

def calc_avg_perc_theta_diff(pop1,pop2):
    outarr=np.ndarray(len(pop1)*len(pop2))
    counter=0
    for i in range(len(pop1)):
        first_theta=pop1[i][4]
        for j in range(len(pop2)):
            second_theta=pop2[j][4]
            outarr[counter]=np.multiply(np.sum(first_theta != second_theta)/first_theta.size,100)
            counter += 1
    return np.mean(outarr)

def calc_avg_perc_lambda_diff(pop1,pop2):
    outarr=np.ndarray(len(pop1)*len(pop2))
    counter=0
    for i in range(len(pop1)):
        first_lambda=pop1[i][5]
        for j in range(len(pop2)):
            second_lambda=pop2[j][5]
            outarr[counter]=np.multiply(np.sum(first_lambda != second_lambda)/first_lambda.size,100)
            counter += 1
    return np.mean(outarr)

print(f"Genomes within the gen10 pop are {calc_avg_perc_genome_diff(gen10,gen10)}% different")
print(f"Are there any links in the GRNs within the gen10 pop which do not vary?: {zeroes_in_grn_stds(gen10)} ")
print(f"GRNs within the population are on average {calc_avg_perc_grn_diff(gen10,gen10)}% different")
#print(f"thetas within the population are {np.multiply(100,np.sum(gen10[0][4] != gen10[1][4])/gen10[0][4].size)}% different")
#theta_std1=np.std(np.array([x[4] for x in gen10[:]]),axis=0)
print(f"Theta gen10 avg diffs={calc_avg_perc_theta_diff(gen10,gen10)}")
#print(f"lambdas within the population are {np.multiply(100,np.sum(gen10[0][5] != gen10[1][5])/gen10[0][5].size)}% different")
#lambda_std1=np.std(np.array([x[5] for x in gen10[:]]),axis=0).reshape(10,)
print(f"Lambdas gen10 avg diffs={calc_avg_perc_lambda_diff(gen10,gen10)}")
print(f"developments within the gen10 population are {np.multiply(100,np.sum(gen10[0][7] != gen10[1][7])/gen10[0][7].size)}% different\n\n")

print(f"Genomes 10 generations apart are {calc_avg_perc_genome_diff(gen10,gen20)}% different")
print(f"GRNs 10 generations apart are on average {calc_avg_perc_grn_diff(gen10,gen20)}% different")
print(f"thetas 10 generations apart are on average {calc_avg_perc_theta_diff(gen10,gen20)}% different")
print(f"lambdas 10 generations apart are {calc_avg_perc_lambda_diff(gen10,gen20)}% different")
print(f"developments 10 generations apart are {np.multiply(100,np.sum(gen10[0][7] != gen20[0][7])/gen10[0][7].size)}% different")