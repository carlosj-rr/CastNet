# CONSTRUCTION PARAMETERS
num_genes = 5
seq_length = 100
prop_unlinked = 0.7
prop_no_threshold = 0.5
thresh_boundaries = (0.1,2)
decay_boundaries = (0,2)
dev_steps = 15 # For the moment no more than 999 is possible
# base_props = (0.25,0.25,0.25,0.25) # T,C,A,G
pop_size = 100 # For the moment, Multiple of 10
#**# pop_stdev = 10 # Not yet in use. Future plan to take pop size of a new generation from a random draw from a distribution

# MUTATION PARAMETERS
thresh_decay_mut_bounds = (-0.01,0.01)
thresh_mutation_rate = 0 # It can also be 0.001, for example
prob_thresh_change = 0
decay_mutation_rate = 0
seq_mutation_rate = 0.001	# Mutation likelihood per base, per generation.  Ex:
				# 1 mutation per 10,000 bases per generation: 1/10000 = 0.0001
				# For now, JC model "hardcoded" (each base is equally likely to mutate into any other base).

#grn_mutation_rate = 0		# Mutation likelihood per gene interaction, per generation.
				# Ex: 1 mutated interaction per 1,000 **interactions** 
				# --> 1/1000 = 0.001 (notice that non-interactions are not taken
				# into account, so it doesn't represent the proportion of the grn
				# that will mutate)
link_mutation_bounds = (-0.01,0.01)
#prob_grn_change = 0.000001 	# Probability that grn mutation will change grn structure
				# (i.e. it will create a new link or remove an existing one).
new_link_bounds = (-2,2)

# SELECTION PARAMETERS
min_reproducin = 0.1
prop_survivors = 0.25 # For the moment, it must result in a whole number when multiplied by the pop_size
tot_offspring = pop_size
select_strategy = "high pressure" # For the moment, just "greedy" and "random"

# REPRODUCTION PARAMETERS
reproductive_strategy = "equal" # for the moment, just "none": all surviving organisms produce the same amount of offspring, regardless of their fitness value. To include later: "winner takes all" - the one with the highest fitness value reproduces the most, and the rest just a little, and other strategies. Eventually I may add recombination.
recomb_pairing = "panmictic" # for the moment, the only option, and recombination is still not implemented.
recomb_style = "vertical" # Options: "vertical", "horizontal", "minimal", "maximal" - how grn matrices are recombined. still non-functional
