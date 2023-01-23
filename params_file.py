# CONSTRUCTION PARAMETERS
num_genes = 5 # number of genes.
seq_length = 3000 # Number of bases per gene. Will be automatically adjusted to the closest number that's a multiple of 3, to have complete codons.
new_link_bounds = (-2,2) # values to initialize GRN.
prop_unlinked = 0.7 # Proportion of sparseness in the GRN.
thresh_boundaries = (0.1,2) # bounds of the uniform distribution used to initialize θs.
decay_boundaries = (0,2) # bounds of the uniform distribution used to initialize λs.
dev_steps = 15 #developmental steps after the t0 input vector
pop_size = 100 # size of the population in number of individuals
reporting_freq = 100 # How many generations should a snapshot be saved to disk

# MUTATION PARAMETERS
seq_mutation_rate = 5.3333e-06	# Mutation prob per base, per generation.
                                    # Empirically 1.06666E-5 was a good rate for 20K gens, 10 genes, 3Kbp each.
                                    
# SELECTION PARAMETERS
min_reproducin = 0.1 # minimum amount of the last gene that has to be found on the last developmental step in order for the system to be 'of reproductive age'
prop_survivors = 0.1 # proportion of the population that survives each generation.
select_strategy = "high pressure" # ""high pressure", "low pressure", and "totally relaxed"

# REPRODUCTION PARAMETERS
reproductive_strategy = "equal" # for the moment, just "equal": all surviving organisms produce the same amount of offspring, regardless of their fitness value. To include later: "winner takes all" - the one with the highest fitness value reproduces the most, and the rest just a little, and other strategies. Eventually I may add recombination.
#recomb_pairing = "panmictic" # recombination is still not implemented.