# CastNet - Authors: Carlos Javier Rivera Rivera and Đorđe Grbić
# License: GNU GPLv3

import copy as cp
import pickle
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime

import numpy as np
import scipy
from scipy import stats
import tqdm
from tqdm import tqdm
import sys
import gif

import CastNet_parameters_12 as pf
from CastNet_out_funcs import *

# A = 1, C = 2, T = 3, G = 4
#
dna_codons = np.array(
[421,
422,
423,
424,
342,
343,
412,
413,
411,
414,
332,
333,
441,
442,
443,
444,
212,
213,
131,
132,
133,
111,
114,
231,
232,
233,
234,
331,
334,
134,
112,
113,
221,
222,
223,
224,
211,
214,
141,
144,
241,
242,
243,
244,
142,
143,
321,
322,
323,
324,
121,
122,
123,
124,
431,
432,
433,
434,
344,
312,
313,
311,
314,
341],
    dtype=int,
)

trans_aas = np.array(
    [65,
65,
65,
65,
67,
67,
68,
68,
69,
69,
70,
70,
71,
71,
71,
71,
72,
72,
73,
73,
73,
75,
75,
76,
76,
76,
76,
76,
76,
77,
78,
78,
80,
80,
80,
80,
81,
81,
82,
82,
82,
82,
82,
82,
83,
83,
83,
83,
83,
83,
84,
84,
84,
84,
86,
86,
86,
86,
87,
89,
89,
95,
95,
95],
    dtype=int,
)

trans_dict = dict(zip(dna_codons,trans_aas))

args=sys.argv
if len(args) > 1:
    run_prefix=sys.argv[1]+"-"
    print(f"Run prefix provided: {run_prefix}")
    name_given=True
else:
    print("A run prefix was not provided, using default \'RUN-\'")
    run_prefix="RUN-"
    name_given=False

def founder_miner(min_fitness=0.6):
    fitness = 0
    while fitness < min_fitness:
        # Importing values for producing the genomic sequences
        n_generation = 0
        n_genes = pf.num_genes
        seq_len = pf.seq_length
        genome, proteome = make_genome_and_proteome(
            seq_len, n_genes, dna_codons, trans_aas
        )
        # Importing the values for producing all the regulatory information.
        # thresholds and decays will have the converse of this probability as 0s. See blow.
        prop_off = pf.prop_unlinked
        thresh_boundaries = pf.thresh_boundaries  # tuple of 2 values.
        decay_boundaries = pf.decay_boundaries  # tuple of 2 values.
        grn = make_grn(n_genes, prop_off)
        thresholds = random_masked_vector(
            n_genes, (1 - prop_off), min(thresh_boundaries), max(thresh_boundaries)
        )
        decays = random_masked_vector(
            n_genes, (1 - prop_off), min(decay_boundaries), max(decay_boundaries)
        )
        # Importing values for the developmental info
        dev_steps = pf.dev_steps
        start_vect = (lambda x: np.array([1] * 1 + [0] * (x - 1)))(n_genes)
        non_kos = np.array(list(map((lambda x: int(95 not in x)),proteome)))
        development = develop(start_vect, grn, decays, thresholds, dev_steps, non_kos)
        genes_on = (development.sum(axis=0) != 0).astype(int)
        fitness = calc_fitness(development)
        out_arr = np.array(
                    (
                        n_generation,
                        genome,
                        proteome,
                        grn,
                        thresholds,
                        decays,
                        start_vect,
                        development,
                        genes_on,
                        fitness,
                    ),
                    dtype=object,
                )
    return out_arr


def make_genome_and_proteome(
    seq_length, num_genes, dna_codons=dna_codons, trans_aas=trans_aas
):
    if seq_length % 3:
        seq_length = seq_length - (seq_length % 3)
        num_codons = int(seq_length / 3)
    else:
        num_codons = int(seq_length / 3)
    genome_arr = np.empty((num_genes, num_codons), dtype=int)
    proteome_arr = np.empty((num_genes, num_codons), dtype=int)
    coding_codons = dna_codons[0:61]
    amino_acids = trans_aas[0:61]
    idx_vect = np.array(range(0, len(coding_codons)))
    for i in range(0, num_genes):
        rand_codon_idx =  np.random.choice(idx_vect, (num_codons))
        
        genome_arr[i] = np.array(coding_codons[rand_codon_idx])
        proteome_arr[i] = np.array(amino_acids[rand_codon_idx])
    return genome_arr, proteome_arr


def make_grn(num_genes, prop_unlinked):
    grn = random_masked_vector(
        num_genes**2, prop_unlinked, pf.new_link_bounds[0], pf.new_link_bounds[1]
    )
    grn = grn.reshape(num_genes, num_genes)
    return grn


# Function that creates a vector of a given amount of values (within a given range), in which a certain proportion of
# the values are masked.
def random_masked_vector(num_vals, prop_zero=0, min_val=0, max_val=1):
    if min_val > max_val:
        raise ValueError(
            f"Minimum value {min_val} is larger than the maximum value {max_val}.\n"
            f"Consider revising the function call to randomMaskedVector()"
        )
    range_size = max_val - min_val
    if prop_zero == 0:
        rpv = np.array(range_size * np.random.random(num_vals) + min_val)
    else:
        mask = np.random.choice((0, 1), num_vals, p=(prop_zero, 1 - prop_zero))
        rpv = np.array(range_size * np.random.random(num_vals) + min_val)
        rpv = (rpv * mask) + 0
    return rpv


def develop(start_vect, grn, decays, thresholds, dev_steps, non_ko_vect):
    #start_vect = cp.deepcopy(start_vect)
    gene_expression_profile = np.ndarray(((pf.dev_steps + 1), pf.num_genes))
    gene_expression_profile[0] = np.array([start_vect])
    # Running the organism's development, and outputting the results
    # in an array called gene_expression_profile
    invect = start_vect
    counter = 1
    for i in range(dev_steps):
        # apply decay to all gene qties. previously: exponentialDecay(invect,decays)
        decayed_invect = (lambda x, l: x * np.exp(-l))(invect, decays)
        # calculate the regulatory effect of the decayed values.
        exp_change = np.matmul(grn, decayed_invect)
        # add the decayed amounts to the regulatory effects
        pre_thresholds = exp_change + decayed_invect
        # a vector to rectify the resulting values to their thresholds.
        thresholder = (pre_thresholds > thresholds).astype(int)
        # rectify with the thresholder vect. This step resulted in the deletion of the 'rectify()' function
        currV = pre_thresholds * thresholder * non_ko_vect
        gene_expression_profile[(i + 1)] = currV
        invect = currV
        counter = counter + 1
    return gene_expression_profile


def calc_fitness(development):
    min_reproducin = pf.min_reproducin
    is_alive = last_gene_expressed(development, min_reproducin)
    if is_alive:
        genes_on = prop_genes_on(development)
        exp_stab = 1 - expression_stability(development) # not doing what it's meant to - see issues on GitHub page.
        # added "1 -" as I realized that the R^2 value would approac 1 the more it assimilated an exponential function.
        sim_to_exp = 1 - exponential_similarity(development)
        fitness_val = np.mean([genes_on, exp_stab, sim_to_exp])
    else:
        fitness_val = 0
    return fitness_val


# is the last gene ever expressed above 'min_reproducin' level?
# AND is it expressed above 0 in the last developmental step?
def last_gene_expressed(development, min_reproducin):
    dev_steps, num_genes = development.shape
    last_col_bool = development[:, (num_genes - 1)] > min_reproducin
    last_val_last_col = development[dev_steps - 1, (num_genes - 1)]
    if last_col_bool.any() and last_val_last_col > 0:
        return_val = True
    else:
        return_val = False
    return return_val


def prop_genes_on(development):
    genes_on = development.sum(axis=0) > 0
    return genes_on.mean()


def expression_stability(development): 
    slopes = np.array(list(map((lambda x: stats.linregress(range(len(x)),x)[0]),development.T[:])))  # get the slope of a linear model for each line
    return np.mean(np.array(list(map(slope_to_angle,slopes)))/90)

def exponential_similarity(development):
    dev_steps, num_genes = development.shape
    row_means = development.mean(axis=1)
    tot_dev_steps = dev_steps
    fitted_line = stats.linregress(range(tot_dev_steps), np.log(row_means))
    r_squared = fitted_line.rvalue**2
    return r_squared

def rad_to_deg(radnum):
    return radnum*180/np.pi

def slope_to_angle(m):
    rez=np.abs(rad_to_deg(np.arctan(m)))
    return rez

class CodonError(Exception):
    pass

def translate_codon(codon):
    #if codon in trans_dict.keys():
    #    aminoac = trans_dict[codon]
    if codon in dna_codons:
        idx = np.where(dna_codons == codon)[0][0]
        aminoac = trans_aas[idx]
    else:
        raise CodonError(
            f"<{codon}> is NOT a valid codon sequence. Please make sure it is one of the following:\n{dna_codons}"
        )
    return aminoac


# Assumes input is a population (i.e. an array of organism arrays).
def grow_pop(in_orgs, out_pop_size, strategy="equal"):
    #in_orgs = cp.deepcopy(in_orgs)
    #num_in_orgs = in_orgs.shape[0]
    seed_start=datetime.now().microsecond
    if len(in_orgs.shape) == 1:
        num_in_orgs=1
    else:
        num_in_orgs=len(in_orgs)
    orgs_per_org = np.round(out_pop_size / num_in_orgs).astype(int)
    #print(f"each of the {num_in_orgs} parental organism(s) will be replicated {orgs_per_org} times")
    curr_pop_size = orgs_per_org * num_in_orgs
    seed_list=list(range(seed_start,seed_start+out_pop_size))
    #print(f"seed list is: {seed_list}")
    if strategy == "equal":
        #print("equal reproduction strategy selected...")
        orgs_per_org = np.repeat(orgs_per_org, num_in_orgs)
    elif strategy == "fitness_linked":
        #print("Reproduction is fitness bound.")
        raise NotImplementedError(
            "Fitness linked reproductive strategy is not yet implemented. Sorry!"
        )
    else:
        raise ValueError(
            f'Reproductive strategy "{strategy}" not recognized.\n'
            f"Strategy must be either 'equal' or 'fitness_linked'."
        )
    counter = 0
    out_pop = np.ndarray((curr_pop_size,), dtype=object)
    # taking each input organism and adding the requested offspring to the output population.
    for i in range(num_in_orgs):
        num_offsp = orgs_per_org[i]
        for j in range(num_offsp):
            if num_in_orgs == 1:
                out_pop[counter] = np.array([x for x in in_orgs],dtype=object)
            else:
                out_pop[counter] = np.array([x for x in in_orgs[i]],dtype=object)
            counter += 1
    #with ProcessPoolExecutor() as pool:
    #    out_pop=np.array(list(pool.map(mutation_wrapper,out_pop,np.repeat(pf.seq_mutation_rate,len(out_pop)),seed_list)))
    out_pop=np.array(list(map(mutation_wrapper,out_pop,np.repeat(pf.seq_mutation_rate,len(out_pop)),seed_list))) # ProcessPoolExecutor was used here to parallelize, but resulted in really inconsistent output.
    out_pop = cleanup_deads(out_pop)  # removing any dead organisms.
    #print(f"{out_pop.shape[0]} organisms survived")
    return out_pop


class IncorrectIndex(Exception):
    pass


# and the mutation rate of the nucleotide sequence (i.e. mutation probability per base).
def mutation_wrapper(orgarr, mut_rateseq,seed=None):
    #orgarrcp = cp.deepcopy(orgarr)
    #print(f"mutating organism using seed {seed}.") # DEBUG STATEMENT - TO REMOVE
    orgarrcp = orgarr
    in_gen_num = orgarrcp[0]
    in_genome = orgarrcp[1]
    in_proteome = orgarrcp[2]
    in_grn = orgarrcp[3]
    in_thresh = orgarrcp[4]
    in_decs = orgarrcp[5]
    start_vect = orgarrcp[6]
    in_dev = orgarrcp[7]
    in_genes_on = (in_dev.sum(axis=0) != 0).astype(int)
    in_fitness = orgarrcp[9]
    #out_genome, out_proteome, mutlocs = mutator_coordinator(in_genome,mut_rateseq)
    mutations = random_mutations(in_genome.size, mut_rateseq, seed)
    if np.any(mutations):
        #print(f"{mutations.size} mutation(s) in this reproductive event")
        mut_coords = cod_pos(mutations, in_genome.shape)
        out_genome, out_proteome, mutlocs = mutate_genome(
            in_genome, in_proteome, mut_coords
        )
        out_grn, out_thresh, out_decs = regulator_mutator(
            in_grn, in_genes_on, in_decs, in_thresh, mutlocs
        )
        non_ko_genes = np.array(list(map((lambda x: int(95 not in x)),out_proteome)))
        out_dev = develop(start_vect, out_grn, out_decs, out_thresh, pf.dev_steps, non_ko_genes)
        out_genes_on = (out_dev.sum(axis=0) != 0).astype(int)
        out_fitness = calc_fitness(out_dev)
    else:
        out_genome = in_genome
        out_proteome = in_proteome
        out_grn = in_grn
        out_thresh = in_thresh
        out_decs = in_decs
        out_dev = in_dev
        out_genes_on = (out_dev.sum(axis=0) != 0).astype(int)
        out_fitness = in_fitness
    out_gen_num = in_gen_num + 1
    out_org = np.array(
            [
                out_gen_num,
                out_genome,
                out_proteome,
                out_grn,
                out_thresh,
                out_decs,
                start_vect,
                out_dev,
                out_genes_on,
                out_fitness,
            ],
        dtype=object,
    )
    return out_org

# BELOW: Currently UNUSED function that does everything 'random_mutations()', 'cod_pos()' and 'mutate_genome()' did, in a much cleverer way, but for some reason was orders of magnitude slower.
#def mutator_coordinator(genome,mut_rateseq):
    rng=np.random.default_rng()
    new_gnome = rng.integers(111,444,genome.size).reshape(genome.shape)
    np.copyto(new_gnome,genome)
    total_bases = new_gnome.size * 3
    mutations = np.random.choice((0, 1), total_bases, p=(1 - mut_rateseq, mut_rateseq))
    m = np.array(np.where(mutations != 0)).flatten()
    if m.size:
        num_genes = new_gnome.shape[0]
        num_codons = new_gnome.shape[1]
        gene_bps = num_codons * 3
        genenum_array = np.repeat(range(num_genes), gene_bps)
        codnum_array=np.tile(np.repeat(range(num_codons),3),num_genes)
        codpos_array = np.tile([0, 1, 2], num_codons * num_genes)
        all_coords = np.array(list(zip(genenum_array,codnum_array,codpos_array)))
        mut_coords = all_coords[m]
        mut_flavor_list = []
        for row in mut_coords:
            coord = tuple(row[0:2])
            old_codon = new_gnome[coord]
            old_aa = translate_codon(old_codon)
            mut_site = row[2]
            new_codon = point_mutate_codon(old_codon,mut_site)
            new_aa = translate_codon(new_codon)
            new_gnome[coord] = new_codon
            mut_flavor_list.append(((new_aa == old_aa)*2) + (new_aa != old_aa) * (new_aa != 95))
        mut_flavor_table = np.array(list(zip(mut_coords[:,0],mut_flavor_list)))
        new_ptome = np.array(list(map(translate_codon,new_gnome.flatten()))).reshape(new_gnome.shape)
    else:
        mut_flavor_table = False
        new_ptome = np.array(list(map(translate_codon,new_gnome.flatten()))).reshape(new_gnome.shape)
    return new_gnome, new_ptome, mut_flavor_table

def random_mutations(genome_size, mut_rateseq, seed=None):  # genome_size is in CODONS
    # Each value in the genome is a codon, so the whole length (in nucleotides) is the codons times 3.
    rng = np.random.default_rng() if seed is None else np.random.default_rng(seed)
    total_bases = genome_size * 3
    mutations = rng.choice((0, 1), total_bases, p=(1 - mut_rateseq, mut_rateseq))
    m = np.array(np.where(mutations != 0)).flatten()
    if m.size:
        output = m
    else:
        output = False
    return output

def cod_pos(muts, gnome_shape):
    # base1=num+1
    num_genes = gnome_shape[0]
    num_codons = gnome_shape[1]
    if np.any(muts < 0):
        raise IncorrectIndex(
            f"There are negative index values in your mutation indices:\n{muts}.\n"
            f"This will result in untractable mutations.\nConsder double-checking the result of randomMutations()"
        )
    out_array = np.ndarray((muts.size, 3), dtype=object)
    gene_bps = num_codons * 3
    #genenum_array = np.ndarray((num_genes, gene_bps), dtype=object)
    #for i in range(num_genes):
    #    genenum_array[i, :] = i
    genenum_array = np.repeat(range(num_genes), gene_bps)
    codpos_array = np.tile([0, 1, 2], num_codons * num_genes)
    codnum_array = np.ndarray((num_genes, gene_bps), dtype=object)
    for i in range(num_genes):
        codnum_array[i, :] = np.repeat(range(num_codons), 3)
    codnum_array = codnum_array.flatten()
    for i in range(muts.size):
        basenum = muts[i]
        mut_val = np.array(
            [genenum_array[basenum], codnum_array[basenum], codpos_array[basenum]]
       )
        out_array[i, :] = mut_val
    return out_array

def mutate_genome(old_gnome, old_prome, mut_coords, seed=None):
    rng = np.random.default_rng() if seed is None else np.random.default_rng(seed)
    gnome = rng.integers(111,444,old_gnome.size).reshape(old_gnome.shape)
    np.copyto(gnome,old_gnome)
    prome = rng.integers(0,100,old_prome.size).reshape(old_prome.shape)
    np.copyto(prome,old_prome)
    #gnome,prome=old_gnome,old_prome
    # get the number of rows in the mutation coordinate array, this is the number of mutations
    mut_num = mut_coords.shape[0]
    muttype_vect = np.ndarray((mut_num, 2), dtype=object)
    if np.any(mut_coords < 0):
        raise IncorrectIndex(
            f"Some indices in the mutation coordinates are negative:\n{mut_coords}\n"
            f"This may result in untractable mutations.\nConsider examining the output of codPos()."
        )
    # DEBUG: CHANGE TO A MAP() APPLIED FUNCTION IF POSSIBLE.
    for i in range(mut_num):
        coordinates = mut_coords[i, :]
        selected_gene = coordinates[0]
        selected_codon_from_gene = coordinates[1]
        selected_codpos = coordinates[2]
        selected_codon = gnome[selected_gene, selected_codon_from_gene]
        prev_aacid = trans_dict[selected_codon] # translate_codon(selected_codon) previously
        mutated_codon = point_mutate_codon(selected_codon, selected_codpos)
        gnome[selected_gene, selected_codon_from_gene] = mutated_codon
        new_aacid = trans_dict[mutated_codon] # translate_codon(mutated_codon) previously
        if prev_aacid == new_aacid:  # Synonymous mutations are plotted as '2'
            muttype = 2
        elif new_aacid == 95:  # Nonsense mutations are plotted as '0'
            muttype = 0
        else:  # Nonsynonymous mutations are plotted as '1'
            muttype = 1
        prome[selected_gene, selected_codon_from_gene] = new_aacid # FIXED ANNOYING BUG
        muttype_vect[i] = (selected_gene, muttype)
    out_genome = gnome
    out_proteome = prome
    return out_genome, out_proteome, muttype_vect



def point_mutate_codon(codon, pos_to_mutate):
    opts=np.array([[1,2,3],[-1,1,2],[-2,-1,1],[-3,-2,-1]])
    if pos_to_mutate == 0:
        factor=100
        value=codon//100
        key=opts[value-1]
    if pos_to_mutate == 1:
        factor=10
        value=codon%100//10
        key=opts[value-1]
    if pos_to_mutate == 2:
        factor=1
        value=codon%10
        key=opts[value-1]
    #else:
    #    raise IncorrectIndex(
    #        f"ERROR: Codon position {pos_to_mutate} is not within the options {(0,1,2)}"
    #    )
    new_codon=codon+np.random.choice(key)*factor
    return new_codon

class MutationTypeError(Exception):
    pass
class NegativeIndex(Exception):
    pass

def regulator_mutator(in_grn, genes_on, in_dec, in_thresh, muttype_vect):
    rng=np.random.default_rng()
    curr_grn=rng.random(in_grn.shape[0]**2).reshape(in_grn.shape)
    np.copyto(curr_grn,in_grn)
    #curr_grn = in_grn
    #curr_thr = cp.deepcopy(in_thresh)
    curr_thr=in_thresh
    #curr_genes_on = cp.deepcopy(genes_on)
    curr_genes_on=genes_on
    #curr_dec = cp.deepcopy(in_dec)
    curr_dec=in_dec
    #curr_muttype_vect = cp.deepcopy(muttype_vect)
    curr_muttype_vect = muttype_vect
    if np.any(curr_muttype_vect < 0):
        raise NegativeIndex(
            f"There is a number in the mutations array that is negative:\n\n{curr_muttype_vect}\n\n"
            f"This may result in untractable mutations, consider inspecting the output of codPos()"
        )
    # TODO check this unused thing and remove if uneccesary
    #inactive_links = np.array(
    #    list(zip(np.where(curr_grn == 0)[0], np.where(curr_grn == 0)[1]))
    #)
    num_genes = curr_genes_on.size
    # Choosing mutations for the thresholds or decays...
    # proportion of total regulatory interactions that are thresholds
    # OR decays (simplified from "2N/(2N+N^2)", where N is the total number of genes.)
    prop = (lambda x: 2 / (2 + x))(num_genes)
    hits = np.nonzero(np.random.choice((0, 1), len(muttype_vect), p=(1 - prop, prop)))[
        0
    ]
    if hits.size > 0:
        mutsarr = curr_muttype_vect[hits]
        out_threshs, out_decs = threshs_and_decs_mutator(in_thresh, in_dec, mutsarr)
        curr_muttype_vect = np.delete(curr_muttype_vect, hits, axis=0)
    else:
        out_threshs, out_decs = curr_thr, curr_dec
    if curr_muttype_vect.size > 0:
        refmat = np.repeat(1, curr_grn.size).reshape(curr_grn.shape)
        refmat[:] = curr_genes_on # Creates a reference matrix where genes that are not on in the development have '0' in their full column, and otherwise '1'
        # For each gene that has a mutation...
        for i in curr_muttype_vect:
            gene = i[0] # get the gene index
            if gene not in range(num_genes): # raise error if index out of range - would signal a bug.
                raise IndexError(
                    f"Gene number {gene} is not within the number of genes in the system "
                    f"(integers between 0 and {num_genes-1})"
                )
            mtype = i[1] # synonymous mutation (2), nonsynonymous mutation (1), or Knockout (0)?
            exprstate = curr_genes_on[gene] # gene on (1) or off (0)?
            
            if mtype in [1, 2]: # If mutation is not KO, Do:
                if exprstate == 0 and mtype == 1: # Gene OFF, Non-Synonymous mutation (uses active sites) - if the gene is OFF, NS mutations happen in the intersection with genes that are ON.
                    active_sites = np.array(
                        list(zip(np.where(refmat == 1)[0], np.where(refmat == 1)[1]))
                    ) # returns a 2D array with the 2D coordinates of all of the non-synonymous sites - col1 = site col idx, col2 = site row idx.
                    # if in the same generation, all genes of the organism have been KO'd...impossible unless there's a bug...Do:
                    if active_sites.size == 0:
                        mutable_sites = np.array(
                            list(zip(np.repeat(gene, num_genes), range(num_genes)))
                        )
                    else:
                        mutable_sites = active_sites[
                            np.where(active_sites[:, 0] == gene)
                        ]
                    site_to_mutate = tuple(
                        mutable_sites[np.random.choice(range(mutable_sites.shape[0]))]
                    )
                    curr_grn[site_to_mutate] = weight_mut(in_grn[site_to_mutate], 0.5)
                # Gene OFF, Synonymous mutation (uses inactive sites)
                if exprstate == 0 and mtype == 2:

                    # silent sites for synonymous mutations
                    inactive_sites = np.array(
                        list(zip(np.where(refmat == 0)[0], np.where(refmat == 0)[1]))
                    )

                    gene_col = inactive_sites[np.where(inactive_sites[:, 1] == gene)]
                    gene_row = inactive_sites[np.where(inactive_sites[:, 0] == gene)]
                    mutable_sites = np.row_stack((gene_col, gene_row))
                    site_to_mutate = tuple(
                        mutable_sites[np.random.choice(range(mutable_sites.shape[0]))]
                    )
                    curr_grn[site_to_mutate] = weight_mut(in_grn[site_to_mutate], 0.5)
                # Gene ON, Non-Synonymous mutation (uses active sites)
                if exprstate == 1 and mtype == 1:

                    # non-silent sites for non-synonymous mutations
                    active_sites = np.array(
                        list(zip(np.where(refmat == 1)[0], np.where(refmat == 1)[1]))
                    )

                    gene_col = active_sites[np.where(active_sites[:, 1] == gene)]
                    gene_row = active_sites[np.where(active_sites[:, 0] == gene)]
                    mutable_sites = np.row_stack((gene_col, gene_row))
                    site_to_mutate = tuple(
                        mutable_sites[np.random.choice(range(mutable_sites.shape[0]))]
                    )

                    curr_grn[site_to_mutate] = weight_mut(in_grn[site_to_mutate], 0.5)

                # Gene ON, Synonymous mutation (uses active sites)
                if exprstate == 1 and mtype == 2:

                    # silent sites for synonymous mutations
                    inactive_sites = np.array(
                        list(zip(np.where(refmat == 0)[0], np.where(refmat == 0)[1]))
                    )

                    # If all genes are on, the 'gene_row' array will come up empty,
                    # so it switches to mutating any gene, a tiny little bit. For the future - it could simply ignore the mutation altogether - think what would reflect biological reality best.
                    if inactive_sites.size == 0:
                        mutable_sites = np.array(
                            list(zip(np.repeat(gene, num_genes), range(num_genes)))
                        )
                        site_to_mutate = tuple(
                            mutable_sites[
                                np.random.choice(range(mutable_sites.shape[0]))
                            ]
                        )
                        curr_grn[site_to_mutate] = weight_mut(
                            in_grn[site_to_mutate], 0.001
                        )
                    else:

                        gene_row = inactive_sites[
                            np.where(inactive_sites[:, 0] == gene)
                        ]
                        mutable_sites = gene_row
                        site_to_mutate = tuple(
                            mutable_sites[
                                np.random.choice(range(mutable_sites.shape[0]))
                            ]
                        )
                        curr_grn[site_to_mutate] = weight_mut(
                            in_grn[site_to_mutate], 0.5
                        )
            elif mtype == 0:
                # If mutation is KO, this will be reflected at the developmental step, so I've commented out where the column and row of the KO'd gene were multiplied by 0.
                # This means that the KO option can also be removed from threshs_and_decs_mutator(), and from the 'muttype_vect' array.
                #curr_grn[gene, :] = 0 
                #curr_grn[:, gene] = 0
                # TODO This is unused
                out_grn = curr_grn
                # Important change so the refmat is updated.
                curr_genes_on[gene] = 0
                # Important change to avoid being unable to find synonymous mutations later on.
                refmat[:] = curr_genes_on
            else:
                raise MutationTypeError(
                    f"Gene{gene}'s mutation type <{mtype}> is unclear,\nit must be one of [0,1,2]."
                )
    out_grn = curr_grn
    return out_grn, out_threshs, out_decs


def threshs_and_decs_mutator(in_thresh, in_dec, mutarr):
    rng=np.random.default_rng()
    curr_thresh = rng.random(len(in_thresh))
    np.copyto(curr_thresh,in_thresh)
    curr_dec=rng.random(len(in_dec))
    np.copyto(curr_dec,in_dec)
    # make a tuple in which the threshold array is the first value, and the decays the second.
    the_tuple = (
        curr_thresh,
        curr_dec,
    )
    # This will allow me to easily choose among them at the time of mutating, see within the for loop.
    # get the genes to be mutated from the mutarray's 1st column
    genes = mutarr[:, 0]
    # go through each gene, and decide randomly whether to make a threshold or a decay mutation in the gene.
    for i in np.arange(len(genes)):
        tuple_idx = np.random.choice((0, 1))
        # extract specific gene number that has to be mutated. This maps to the thresh and dec arrays.
        gene_num = genes[i]
        isit_ko = mutarr[i, 1] == 0 # consider removing this part - KO mutations are already implemented in the develop() function. But this would also imply removing the KO mutations from the muttypes_vect.
        if isit_ko:
            new_value = 0
        else:
            new_value = abs(weight_mut(the_tuple[tuple_idx][gene_num]))
        the_tuple[tuple_idx][gene_num] = new_value
    out_thresh, out_decs = (the_tuple[0], the_tuple[1])
    return out_thresh, out_decs


def weight_mut(value, scaler=0.01):
    # Make sure value is positive
    val = abs(value)
    if val == 0:
        """For 0, simply get 1, and then modify it by the scale
        This is for activating values that are off."""
        val = scaler / scaler
    scaled_val = val * scaler  # scale the value
    # add the scaled portion to the total value to get the final result.
    newVal = value + np.random.uniform(-scaled_val, scaled_val)
    return newVal


def store(thing, filename = None):
    now = datetime.now()
    moment = now.strftime("%m-%d-%Y_%H-%M-%S")
    if not filename:
        filename = "./EvolRun_" + moment + ".pkl"
    with open(filename, "wb") as fh:
        pickle.dump(thing, fh)
    print(f"Your session results were saved in {filename}.")


def cleanup_deads(in_pop):
    # in_pop=cp.deepcopy(in_pop)
    tot_orgs = in_pop.shape[0]
    fitnesses = np.array([x[9] for x in in_pop[:]])
    live_ones = np.nonzero(fitnesses)[0]
    if live_ones.size == tot_orgs:
        out_pop = in_pop
    elif live_ones.size != 0:
        out_pop = in_pop[live_ones]
    elif live_ones.size == 0:
        print(f"Your population went extinct. Sorry for your loss.")
        out_pop = np.array([])
    return out_pop


class UnknownSelectiveStrategy(Exception):
    pass


def select(in_pop, p=0.1, strategy="high pressure"):
    #print(in_pop)
    if (lambda x: len(x.shape) == 1)(in_pop): #if it is a single organism
        print("You have passed a single organism for selection. Confirming it is alive...")
        print(in_pop)
        if in_pop[9] > 0:
            print("It is alive, congratulations, this will be your survivor")
            return in_pop
        else:
            print("Ah well. your only survivor was dead. This lineage is now extinct.")
            return False
    pop_size = in_pop.shape[0]
    num_survivors = int(pop_size * p)
    fitnesses = np.array([x[9] for x in in_pop[:]])
    if num_survivors == 0:
        print(
            f"A proportion of {p} results in 0 survivors, you've selected your population into extinction."
        )
        return False
    else:
        out_pop=np.ndarray(num_survivors,dtype=object)
    if strategy == "high pressure":
        # returns the **indices** for the top 'num_survivors' fitnesses.
        out_idcs = np.argpartition(fitnesses, -num_survivors)[-num_survivors:]
    elif strategy == "low pressure" and p < 0.5:
        half = int(pop_size / 2)
        top_half = np.argpartition(fitnesses, -half)[-half:]
        out_idcs = np.random.choice(top_half, num_survivors, replace=False)
    elif strategy == "low pressure" and p >= 0.5:
        out_idcs = np.random.choice(range(pop_size), num_survivors, replace=False)
    elif strategy == "totally relaxed":
        out_idcs = np.random.choice(range(pop_size), num_survivors, replace=False)
    else:
        raise UnknownSelectiveStrategy(
            f"Selective strategy {strategy} is not recognized\n"
            f'This value should be any of "high pressure", "low pressure", or "totally relaxed".\n '
            f"Please double-check your input."
        )
    out_pop = in_pop[out_idcs]
    if num_survivors == 1:
        out_pop = out_pop[0]
    return out_pop

# Add a storing to disk command here
def randsplit(in_pop, out_pop_size):
    # in_pop=cp.deepcopy(in_pop)
    if len(in_pop.shape) == 1:
        inpopsize=1
    else:
        inpopsize=len(in_pop)
    if inpopsize > 2:
        idcs_lina = np.random.choice(
            range(inpopsize), int(inpopsize / 2), replace=False
        )
        idcs_linb = np.array(
            [rand for rand in np.arange(inpopsize) if rand not in idcs_lina]
        )
        lina = grow_pop(in_pop[idcs_lina], out_pop_size, "equal")
        linb = grow_pop(in_pop[idcs_linb], out_pop_size, "equal")
    elif inpopsize == 2:
        l = [1, 0]
        idx_lina = np.random.choice((0, 1), 1)[0]
        idx_linb = l[idx_lina]
        lina = grow_pop(in_pop[idx_lina], out_pop_size, "equal")
        linb = grow_pop(in_pop[idx_linb], out_pop_size, "equal")
    elif inpopsize == 1:
        print("Input population has single individual only")
        lina = grow_pop(in_pop, out_pop_size, "equal")
        linb = grow_pop(in_pop, out_pop_size, "equal")
    elif inpopsize < 1:
        raise ValueError(
            f"Input population doesn't have enough individuals: {inpopsize}."
        )
    return lina, linb

def branch_evol(parent_pop, ngens,seed=None,branch_id=0,reporting_freq=pf.reporting_freq):
    #rng = np.random.default_rng() if seed is None else np.random.default_rng(seed)
    in_pop = cp.deepcopy(parent_pop)
    #branch=np.ndarray((ngens,),dtype=object)
    #branch_key=str(np.random.randint(0,1e10))
    print(f"Size of parental population: {len(in_pop)}")
    if len(in_pop.shape) > 1:
        mean_grn=np.mean([x[3].flatten() for x in in_pop[:]],axis=0)
        flat_grns_arr=np.insert(mean_grn,0,-1,axis=0)
        mean_fitnesses=[np.mean([x[9] for x in in_pop[:]])]
    else:
        mean_grn=in_pop[3].flatten()
        flat_grns_arr=np.insert(mean_grn,0,-1,axis=0)
        mean_fitnesses=[in_pop[9]]
    #out_grns=[]
    out_devs=[]
    if parent_pop.size > 0:
        for gen in tqdm(range(ngens),desc=f"{branch_id} branch evolution progress:", ascii=False,ncols=100):
            if parent_pop.size > 0:
                #print(f"producing generation {gen+1}")
                survivors = select(parent_pop, pf.prop_survivors, pf.select_strategy)
                if type(survivors) == bool:
                    filename="Extinction1_branch"+str(branch_id)+"_generation_"+str(gen)+".pkl"
                    print(f"Branch has gone extinct, packaging and outputting a truncated branch of {gen} generation(s) into {filename}.")
                    store(parent_pop,filename) # BRANCH SELECTED INTO EXTINCTION EXIT
                    return
                #print(f"Survivor number is {len(survivors)}.")
                next_pop = grow_pop(survivors, pf.pop_size, pf.reproductive_strategy)
                if next_pop.size == 0:
                    filename="Extinction2_branch"+str(branch_id)+"_generation_"+str(gen)+".pkl"
                    print(f"The population from branch {branch_id} had no viable offspring and went extinct. Packaging and saving a truncated branch at {gen} generation(s) into {filename}.")
                    np.save(filename,parent_pop)
                    return # POPULATION UNVIABLE EXIT
                parent_pop=next_pop
                #print(f"Generation {gen+1} of {ngens} completed.")
                parent_pop = next_pop
                if (gen) % reporting_freq == 0:
                    prefix="Generation_"+str(gen)+"_branch_"+str(branch_id)
                    #grn_fig=draw_avg_grns(next_pop,gen)
                    #out_grns.append(grn_fig)
                    dev_fig=plot_avg_developments(next_pop,gen)
                    out_devs.append(dev_fig)
                    mean_grn=np.mean([ x[3].flatten() for x in next_pop[:] ],axis=0)
                    with_gennum=np.insert(mean_grn,0,gen,axis=0)
                    flat_grns_arr=np.vstack((flat_grns_arr,with_gennum))
                    mean_fitnesses.append(np.mean([x[9] for x in next_pop[:]]))
                    #np.save(filename,next_pop)
                if (gen+1) == ngens:
                    filename="TipGeneration_"+str(gen+1)+"_branch_"+str(branch_id)
                    #grn_file="Branch_"+str(branch_id)+"_AvgGRN.gif"
                    #devs_file="Branch_"+str(branch_id)+"_AvgDevs.gif"
                    #gif.save(grn_fig,grn_file,duration=1000)
                    #gif.save(dev_fig,devs_file,duration=1000)
                    #print(f"####### Saving population {gen+1} to {filename}")
                    animator(out_devs,"Branch_"+str(branch_id))
                    np.savetxt(filename+"_meanGRNs.csv",flat_grns_arr,delimiter=",")
                    np.savetxt(filename+"_fitnessProg.csv",mean_fitnesses,delimiter=",")
                    return(next_pop)
            else:
                print("Your population was extinguished.")
                return
        return(next_pop)
    else:
        print(f"Input population {parent_pop} has no individuals. Stopping simulation.")
        return

def animator(devs,file_prefix):
    gif.save(devs,file_prefix+"_avgDev.gif")
    return

def unpickle(filename):
    pickle_off = open(filename, "rb")
    output = pickle.load(pickle_off)
    return output


def gene_ali_saver(organism_array, outfile_prefix="outfile"):
    num_orgs = organism_array.size
    rand_seqs = np.random.choice(num_orgs, 10)
    num_genes = organism_array[0][1].shape[0]
    # TODO this is unused
    sequences_array = np.array([x[1] for x in organism_array])
    for i in range(num_genes):
        filename = outfile_prefix + "_gene_" + str(i) + ".fas"
        with open(filename, "w") as gene_file:  
            for j in range(rand_seqs.size):
                seq_name = ">" + outfile_prefix + "_" + str(i) + "_org" + str(j)
                sequence=''.join(organism_array[j][1][i])
                print(seq_name, file=gene_file)
                print(sequence, file=gene_file)
        print("Gene", str(i), "done")
        
def main_parallel():
    print("Running parallel flavor")
    founder = founder_miner()
    print("Founder created")
    founder_pop = grow_pop(founder, pf.pop_size, "equal")
    print("Founder pop created")
    anc1_stem, anc2_stem = randsplit(founder_pop, pf.pop_size)
    print("Founding split created")
    seed_start=datetime.now().microsecond
    seed_list=list(range(seed_start,seed_start+2))
    br_randnums = np.random.randint(0,1e10,2).astype(str)
    br_prefix=['ancestor1_','ancestor2_']
    br_lengths=np.repeat(10000,2)
    br_ids = [x+y for x,y in zip(br_prefix,br_randnums)]
    with ProcessPoolExecutor() as pool:
        anc1_tip,anc2_tip=list(pool.map(branch_evol,[anc1_stem,anc2_stem],br_lengths,seed_list,br_prefix))
    #anc1_tip = branch_evol(anc1_stem,10000,br_ids[0])
    #anc2_tip = branch_evol(anc2_stem,10000,br_ids[1])
    # Ancestor simulation completed here.
    br_randnums = np.random.randint(0,1e10,4).astype(str)
    br_prefix = [br_ids[0]+'-leafA_',br_ids[0]+'-leafB_',br_ids[1]+'-leafC_',br_ids[1]+'-leafD_']
    br_ids = [x+y for x,y in zip(br_prefix,br_randnums)]
    br_lengths=np.repeat(10000,4)
    seed_start=datetime.now().microsecond
    seed_list=list(range(seed_start,seed_start+4))
    leafa_stem, leafb_stem = randsplit(anc1_tip, pf.pop_size)
    leafc_stem, leafd_stem = randsplit(anc2_tip, pf.pop_size)
    with ProcessPoolExecutor() as pool:
        leafa_tip,leafb_tip,leafc_tip,leafd_tip = list(pool.map(branch_evol,[leafa_stem,leafb_stem,leafc_stem,leafd_stem],br_lengths,seed_list,br_prefix))
    #leafa_tip = branch_evol(leafa_stem,10000,br_ids[0])
    #leafb_tip = branch_evol(leafb_stem,10000,br_ids[1])
    
    #leafc_tip = branch_evol(leafc_stem,10000,br_ids[2])
    #leafd_tip = branch_evol(leafd_stem,10000,br_ids[3])
    
    return(np.array([founder_pop,anc1_tip,anc2_tip,leafa_tip,leafb_tip,leafc_tip,leafd_tip],dtype=object))

def main_experiment():
    founder = founder_miner()
    print("Founder created")
    founder_pop = grow_pop(founder, pf.pop_size, "equal")
    print("Founder pop created")
    anc1_stem, anc2_stem = randsplit(founder_pop, pf.pop_size)
    print("Founding split created")
    br_randnums = np.random.randint(0,1e10,2).astype(str)
    br_prefix=['ancestor1_','ancestor2_']
    br_ids = [x+y for x,y in zip(br_prefix,br_randnums)]
    anc1_tip = branch_evol(anc1_stem,10000,br_ids[0])
    anc2_tip = branch_evol(anc2_stem,10000,br_ids[1])
    # Ancestor simulation completed here.
    br_randnums = np.random.randint(0,1e10,4).astype(str)
    br_prefix = [br_ids[0]+'-leafA_',br_ids[0]+'-leafB_',br_ids[1]+'-leafC_',br_ids[1]+'-leafD_']
    br_ids = [x+y for x,y in zip(br_prefix,br_randnums)]
    leafa_stem, leafb_stem = randsplit(anc1_tip, pf.pop_size)
    leafa_tip = branch_evol(leafa_stem,10000,br_ids[0])
    leafb_tip = branch_evol(leafb_stem,10000,br_ids[1])
    leafc_stem, leafd_stem = randsplit(anc2_tip, pf.pop_size)
    leafc_tip = branch_evol(leafc_stem,10000,br_ids[2])
    leafd_tip = branch_evol(leafd_stem,10000,br_ids[3])
    
    return(np.array([founder_pop,anc1_tip,anc2_tip,leafa_tip,leafb_tip,leafc_tip,leafd_tip],dtype=object))

if __name__ == "__main__":
    if pf.parallel:
        result = main_parallel()
    else:
        result = main_experiment()
    store(result)
    tip_names=["founder","ancAB","ancCD","tip_A","tip_B","tip_C","tip_D"]
    ali_saver(run_prefix,result,tip_names)
    print("Analysis completed", result.shape)