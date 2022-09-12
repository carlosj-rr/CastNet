# CONSTRUCTION PARAMETERS
num_genes = 10
seq_length = 3000 # Number of bases per gene. Will be automatically adjusted to the closest number that's a multiple of 3, to have complete codons.
prop_unlinked = 0.7
#prop_no_threshold = 0.5
thresh_boundaries = (0.1,2)
decay_boundaries = (0,2)
dev_steps = 15 # For the moment no more than 999 is possible
pop_size = 100 # For the moment, Multiple of 10
reporting_freq=100 # How many generations should a snapshot be saved to disk

# MUTATION PARAMETERS
#thresh_decay_mut_bounds = (-0.01,0.01)
#thresh_mutation_rate = 0 # It can also be 0.001, for example
#prob_thresh_change = 0
#decay_mutation_rate = 0
seq_mutation_rate = 1.0666666E-5	# Mutation prob per base, per generation.
                                    # Empirically 1.06666E-5 was a good rate for 20K gens, 10 genes, 3Kbp each.
#link_mutation_bounds = (-0.01,0.01)

new_link_bounds = (-2,2)

# SELECTION PARAMETERS
min_reproducin = 0.1
prop_survivors = 0.1 # For the moment, it must result in a whole number when multiplied by the pop_size
#tot_offspring = pop_size
select_strategy = "high pressure" # ""high pressure", "low pressure", and "totally relaxed""

# REPRODUCTION PARAMETERS
reproductive_strategy = "equal" # for the moment, just "none": all surviving organisms produce the same amount of offspring, regardless of their fitness value. To include later: "winner takes all" - the one with the highest fitness value reproduces the most, and the rest just a little, and other strategies. Eventually I may add recombination.
#recomb_pairing = "panmictic" # for the moment, the only option, and recombination is still not implemented.
#recomb_style = "vertical" # Options: "vertical", "horizontal", "minimal", "maximal" - how grn matrices are recombined. still non-functional
