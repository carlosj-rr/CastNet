# CONSTRUCTION PARAMETERS
num_genes = 5
seq_length = 3000 # Number of bases per gene. Will be automatically adjusted to the closest number that's a multiple of 3, to have complete codons.
prop_unlinked = 0.7
#prop_no_threshold = 0.5
thresh_boundaries = (0.1,2) # range to initialize θs.
decay_boundaries = (0,2) # range to initialize λs.
dev_steps = 15
pop_size = 100 # How many individuals per population
reporting_freq=100 # How many generations should a snapshot be saved to disk

# SEQUENCE MUTATION PARAMETERS
seq_mutation_rate = 5.3333e-06	# Mutation prob per base, per generation.
                                    # Empirically 1.06666E-5 was a good rate for 20K gens, 10 genes, 3Kbp each.

new_link_bounds = (-2,2) # values to initialize GRN.

# SELECTION PARAMETERS
min_reproducin = 0.1
prop_survivors = 0.1
#tot_offspring = pop_size
select_strategy = "high pressure" # ""high pressure", "low pressure", and "totally relaxed"

# REPRODUCTION PARAMETERS
reproductive_strategy = "equal" # for the moment, just "equal": all surviving organisms produce the same amount of offspring, regardless of their fitness value. To include later: "winner takes all" - the one with the highest fitness value reproduces the most, and the rest just a little, and other strategies. Eventually I may add recombination.
