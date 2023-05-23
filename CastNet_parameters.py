# Authors: Carlos Javier Rivera Rivera and Đorđe Grbić
# License: GNU GPLv3

parallel = True # True if the run is meant to be parallelized. Note that this implies that the user has properly designed the 'main_parallel()' function

# CONSTRUCTION PARAMETERS
num_genes = 1 # How many genes in the system?
seq_length = 10000 # How long each gene's sequence is, in nucleotides? Will be automatically adjusted to the closest number that's a multiple of 3, to have complete codons.
prop_unlinked = 0.7 # GRN mask of sparseness - the proportion of GRN regulatory values that will be masked to 0 upon initialisation.
thresh_boundaries = (0.1,2) # range to initialize θs (the threshold under which each gene needs to be activated in order to be switched on).
decay_boundaries = (0,2) # range to initialize λs.
dev_steps = 15 # How many developmental steps?
pop_size = 10 # How many individuals per population?
reporting_freq=100 # Every __ generations, store a a snapshot of the population. Each snapshot will store:
                    # 1. the mean GRN values of the population,
                    # 2. the mean fitness, 
                    # 3. a plot showing the mean phenotype of the population, to be placed in an animated gif.

# SEQUENCE MUTATION PARAMETERS
seq_mutation_rate = 5.3333e-06	# Mutation prob per base, per generation.
                                    # Empirically 1.06666E-5 was a good rate for 20K gens, 10 genes, 3Kbp each.
model = "JC" # Model for base mutation, "JC", "K80", and "GTR" implemented. If using K80 or GTR, the relevant parameters below have to be uncommented and filled in.
# K80 model parameters
# k80_alpha = alpha # if using K80, uncomment and add here the probability of transitional mutations (i.e. purine-purine or pyrimidine-pyrimidine)
# k80_beta = beta # if using K80, uncomment and add here the probability of transversional mutations (i.e. purine-pyrimidine and vice versa) # IMPORTANT: k80_alpha + k80_beta MUST add up to 1

# GTR model parameters
# gtr_t2c = val # probability of T's mutating into C's, and vice versa
# gtr_t2a = val # probability of T's mutating into A's, and vice versa
# gtr_t2g = val # probability of T's mutating into G's, and vice versa
# gtr_c2a = val # probability of C's mutating into A's, and vice versa
# gtr_c2g = val # probability of C's mutating into G's, and vice versa
# gtr_a2g = val # probability of A's mutating into G's, and vice versa

# IMPORTANT: GTR PROBABILITIES OF A SINGLE BASE MUST ALL ADD UP TO 1, (e.g. gtr_t2c + gtr_t2a + gtr_t2g = gtr_t2c + gtr_c2a + gtr_c2g = gtr_t2g + gtr_c2g + gtr_a2g = 1)

new_link_bounds = (-2,2) # lower and upper limit of uniform distribution from which the initial regulatory values of the GRN will be drawn.

# SELECTION PARAMETERS
min_reproducin = 0.1 # minumum level of expression that the 'system maturity indicator' gene must have by the end of development.
prop_survivors = 0.1 # what proportion of the population passes on to the next generation.
select_strategy = "low pressure" # ""high pressure", "low pressure", and "totally relaxed"

# REPRODUCTION PARAMETERS
reproductive_strategy = "equal" # for the moment, just "equal" implemented: all surviving organisms produce the same amount of offspring, regardless of their fitness value. To include later here: "fitness bound" - organisms reproduce with a success proportional to their fitness, and "winner takes all" - the one with the highest fitness value reproduces the most, and the rest just a little, and other strategies.
# Eventually I should add recombination.
