import copy as cp
import pickle
from concurrent.futures import ProcessPoolExecutor  # for using multiple cores.
from datetime import datetime

import numpy as np
import scipy

import params_file as pf

dna_codons = np.array(
    [
        131,
        132,
        133,
        134,
        121,
        122,
        124,
        123,
        112,
        113,
        111,
        114,
        142,
        143,
        141,
        144,
        231,
        232,
        234,
        233,
        221,
        222,
        224,
        223,
        212,
        213,
        211,
        214,
        241,
        242,
        244,
        243,
        431,
        432,
        434,
        433,
        421,
        422,
        424,
        423,
        412,
        413,
        411,
        414,
        441,
        442,
        444,
        443,
        321,
        322,
        324,
        323,
        332,
        333,
        331,
        334,
        312,
        313,
        342,
        343,
        344,
        311,
        314,
        341,
    ],
    dtype=int,
)

trans_aas = np.array(
    [73,
    73,
    73,
    77,
    84,
    84,
    84,
    84,
    78,
    78,
    75,
    75,
    83,
    83,
    82,
    82,
    76,
    76,
    76,
    76,
    80,
    80,
    80,
    80,
    72,
    72,
    81,
    81,
    82,
    82,
    82,
    82,
    86,
    86,
    86,
    86,
    65,
    65,
    65,
    65,
    68,
    68,
    69,
    69,
    71,
    71,
    71,
    71,
    83,
    83,
    83,
    83,
    70,
    70,
    76,
    76,
    89,
    89,
    67,
    67,
    87,
    95,
    95,
    95
    ],
    dtype=int,
)


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
        development = develop(start_vect, grn, decays, thresholds, dev_steps)
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
    idx_vect = np.array(range(0, len(dna_codons) - 3))
    genome_arr = np.empty((num_genes, num_codons), dtype=int)
    proteome_arr = np.empty((num_genes, num_codons), dtype=int)
    for i in range(0, num_genes):
        rand_codon_idx = np.hstack(
            (
                np.random.choice(idx_vect, (num_codons - 1)),
                np.random.choice((61, 62, 63), 1),
            )
        )
        genome_arr[i] = np.array(dna_codons[rand_codon_idx])
        proteome_arr[i] = np.array(trans_aas[rand_codon_idx])
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


def develop(start_vect, grn, decays, thresholds, dev_steps):
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
        currV = pre_thresholds * thresholder
        gene_expression_profile[(i + 1)] = currV
        invect = currV
        counter = counter + 1
    return gene_expression_profile


def calc_fitness(development):
    min_reproducin = pf.min_reproducin
    is_alive = last_gene_expressed(development, min_reproducin)
    if is_alive:
        genes_on = prop_genes_on(development)
        exp_stab = expression_stability(development) # not doing what it's meant to - see issues on GitHub page.
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


def expression_stability(development):  # TODO I haven't thought deeply about this.
    row_sums = development.sum(axis=1)  # What proportion of the data range is
    # TODO the stdev? Less = better
    stab_val = row_sums.std() / (row_sums.max() - row_sums.min())
    return stab_val


def exponential_similarity(development):
    dev_steps, num_genes = development.shape
    row_means = development.mean(axis=1)
    tot_dev_steps = dev_steps
    fitted_line = scipy.stats.linregress(range(tot_dev_steps), np.log(row_means))
    r_squared = fitted_line.rvalue**2
    return r_squared


class CodonError(Exception):
    pass


def translate_codon(codon):
    if codon in dna_codons:
        idx = np.where(dna_codons == codon)[0][0]
        aminoac = trans_aas[idx]
    else:
        raise CodonError(
            f"<{codon}> is NOT a valid codon sequence. Please make sure it is one of the following:\n{dna_codons}"
        )
    return aminoac


# Assumes input is a population (i.e. an array of organism arrays), it should crash if it doesn't find 2 dimensions.
def grow_pop(in_orgs, out_pop_size, strategy="equal"):
    #in_orgs = cp.deepcopy(in_orgs)
    #num_in_orgs = in_orgs.shape[0]
    if len(in_orgs.shape) == 1:
        num_in_orgs=1
    else:
        num_in_orgs=len(in_orgs)
    orgs_per_org = np.array([np.round(out_pop_size / num_in_orgs).astype(int)])
    curr_pop_size = orgs_per_org * num_in_orgs
    if strategy == "equal":
        print("equal reproduction strategy selected...")
        orgs_per_org = np.repeat(orgs_per_org, num_in_orgs)
    elif strategy == "fitness_linked":
        print("Reproduction is fitness bound.")
        raise NotImplementedError(
            "Fitness linked reproductive strategy is not yet implemented. Sorry!"
        )
    else:
        raise ValueError(
            f'Reproductive strategy "{strategy}" not recognized.\n'
            f"Strategy must be either 'equal' or 'fitness_linked'."
        )
    counter = 0
    out_pop = np.ndarray((curr_pop_size[0],), dtype=object)
    # taking each input organism and adding the requested offspring to the output population.
    for i in range(num_in_orgs):
        num_offsp = orgs_per_org[i]
        for j in range(num_offsp):
            out_pop[counter] = cp.deepcopy(in_orgs[i])
            counter = counter + 1
    out_pop=np.array(list(map(mutation_wrapper,out_pop,np.repeat(pf.seq_mutation_rate,len(out_pop)))))
    out_pop = cleanup_deads(out_pop)  # removing any dead organisms.
    print(f"{out_pop.shape[0]} organisms survived")
    return out_pop


class NegativeIndex(Exception):
    pass


# Input is an organism array, as produced by the founder_miner() function,
# and the mutation rate of the nucleotide sequence (i.e. mutation probability per base).
def mutation_wrapper(orgarr, mut_rateseq):
    #orgarrcp = cp.deepcopy(orgarr[0])
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
    mutations = random_mutations(in_genome.size, mut_rateseq)
    if np.any(mutations):
        print(f"{mutations.size} mutation(s) in this reproductive event")
        mut_coords = cod_pos(mutations, in_genome.shape)
        out_genome, out_proteome, mutlocs = mutate_genome(
            in_genome, in_proteome, mut_coords
        )
        out_grn, out_thresh, out_decs = regulator_mutator(
            in_grn, in_genes_on, in_decs, in_thresh, mutlocs
        )
        out_dev = develop(start_vect, out_grn, out_decs, out_thresh, pf.dev_steps)
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


def random_mutations(genome_size, mut_rateseq):  # genome_size is in CODONS
    # Each value in the genome is a codon, so the whole length (in nucleotides) is the codons times 3.
    total_bases = genome_size * 3
    mutations = np.random.choice((0, 1), total_bases, p=(1 - mut_rateseq, mut_rateseq))
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
        raise NegativeIndex(
            f"There are negative index values in your mutation indices:\n{muts}.\n"
            f"This will result in untractable mutations.\nConsder double-checking the result of randomMutations()"
        )
    out_array = np.ndarray((muts.size, 3), dtype=object)
    gene_bps = num_codons * 3
    genenum_array = np.ndarray((num_genes, gene_bps), dtype=object)
    for i in range(num_genes):
        genenum_array[i, :] = i
    genenum_array = genenum_array.flatten()
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


def mutate_genome(old_gnome, old_prome, mut_coords):
    gnome = cp.deepcopy(old_gnome)
    prome = cp.deepcopy(old_prome)
    #gnome,prome=old_gnome,old_prome
    # get the number of rows in the mutation coordinate array, this is the number of mutations
    mut_num = mut_coords.shape[0]
    muttype_vect = np.ndarray((mut_num, 2), dtype=object)
    if np.any(mut_coords < 0):
        raise NegativeIndex(
            f"Some indices in the mutation coordinates are negative:\n{mut_coords}\n"
            f"This may result in untractable mutations.\nConsider examining the output of codPos()."
        )
    for i in range(mut_num):
        coordinates = mut_coords[i, :]
        selected_gene = coordinates[0]
        selected_codon_from_gene = coordinates[1]
        selected_codpos = coordinates[2]
        selected_codon = gnome[selected_gene, selected_codon_from_gene]
        prev_aacid = translate_codon(selected_codon)
        mutated_codon = point_mutate_codon(selected_codon, selected_codpos)
        gnome[selected_gene, selected_codon_from_gene] = mutated_codon
        new_aacid = translate_codon(mutated_codon)
        if prev_aacid == new_aacid:  # Synonymous mutations are plotted as '2'
            muttype = 2
        elif new_aacid == "_":  # Nonsense mutations are plotted as '0'
            muttype = 0
        else:  # Nonsynonymous mutations are plotted as '1'
            muttype = 1
        prome[selected_gene, selected_codpos] = new_aacid
        muttype_vect[i] = (selected_gene, muttype)
    out_genome = gnome
    out_proteome = prome
    return out_genome, out_proteome, muttype_vect


def point_mutate_codon(codon, pos_to_mutate): # Has to be made consisten with the new integer codification of bases.
    if pos_to_mutate < 0:
        raise NegativeIndex(
            f"Codon position {pos_to_mutate} is negative\nThis can cause intractable mutations\n"
            f"Consider verifying the output of codPos()"
        )
    bases = (1, 2, 3, 4)
    txt_codon=np.array(list(str(codon)))
    subs=np.setxor1d(bases,txt_codon[pos_to_mutate])
    def reform_txt_codon(txt_codon,pos_to_mutate,sub):
        txt_codon[pos_to_mutate]=sub
        return(int(''.join(txt_codon)))
    options=np.ndarray(3,dtype=int)
    for i in range(len(subs)):
        options[i]=reform_txt_codon(txt_codon,pos_to_mutate,subs[i])
    new_codon = np.random.choice(options)
    return new_codon

class MutationTypeError(Exception):
    pass


def regulator_mutator(in_grn, genes_on, in_dec, in_thresh, muttype_vect):
    curr_grn = cp.deepcopy(in_grn)
    curr_thr = cp.deepcopy(in_thresh)
    curr_genes_on = cp.deepcopy(genes_on)
    curr_dec = cp.deepcopy(in_dec)
    curr_muttype_vect = cp.deepcopy(muttype_vect)
    if np.any(curr_muttype_vect < 0):
        raise NegativeIndex(
            f"There is a number in the mutations array that is negative:\n\n{curr_muttype_vect}\n\n"
            f"This may result in untractable mutations, consider inspecting the output of codPos()"
        )
    # TODO check this unused thing and remove if uneccesary
    inactive_links = np.array(
        list(zip(np.where(curr_grn == 0)[0], np.where(curr_grn == 0)[1]))
    )
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
        refmat[:] = curr_genes_on

        for i in curr_muttype_vect:
            gene = i[0]
            if gene not in range(num_genes):
                raise IndexError(
                    f"Gene number {gene} is not within the number of genes in the system "
                    f"(integers between 0 and {num_genes-1})"
                )
            mtype = i[1]
            exprstate = curr_genes_on[gene]

            if mtype in [1, 2]:
                # Gene OFF, Non-Synonymous mutation (uses active sites)
                if exprstate == 0 and mtype == 1:
                    # non-silent sites for non-synonymous mutations
                    active_sites = np.array(
                        list(zip(np.where(refmat == 1)[0], np.where(refmat == 1)[1]))
                    )

                    # if in the same generation, all genes of the organism have been KO'd...
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

                # Gene OFF, Synonymous mutation (uses inactive sites)
                if exprstate == 1 and mtype == 2:

                    # silent sites for synonymous mutations
                    inactive_sites = np.array(
                        list(zip(np.where(refmat == 0)[0], np.where(refmat == 0)[1]))
                    )

                    # If all genes are on, the 'gene_row' array will come up empty,
                    # so it switches to mutating any gene, a tiny little bit.
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
                # If mutation is KO
                curr_grn[gene, :] = 0
                curr_grn[:, gene] = 0
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
    in_thresh = cp.deepcopy(in_thresh)
    in_dec = cp.deepcopy(in_dec)
    # make a tuple in which the threshold array is the first value, and the decays the second.
    the_tuple = (
        in_thresh,
        in_dec,
    )
    # This will allow me to easily choose among them at the time of mutating, see within the for loop.
    # get the genes to be mutated from the mutarray's 1st column
    genes = mutarr[:, 0]
    # go through each gene, and decide randomly whether to make a threshold or a decay mutation in the gene.
    for i in np.arange(len(genes)):
        tuple_idx = np.random.choice((0, 1))
        # extract specific gene number that has to be mutated. This maps to the thresh and dec arrays.
        gene_num = genes[i]
        isit_ko = mutarr[i, 1] == 0
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


def store(thing):
    now = datetime.now()
    moment = now.strftime("%m-%d-%Y_%H-%M-%S")
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
        out_pop = cp.deepcopy(in_pop)
    elif live_ones.size != 0:
        out_pop = cp.deepcopy(in_pop[live_ones])
    elif live_ones.size == 0:
        print(f"Your population went extinct. Sorry for your loss.")
        out_pop = np.array([])
    return out_pop


class UnknownSelectiveStrategy(Exception):
    pass


def select(in_pop, p=0.1, strategy="high pressure"):
    pop_size = in_pop.shape[0]
    num_survivors = int(pop_size * p)
    fitnesses = np.array([x[9] for x in in_pop[:]])
    if num_survivors == 0:
        raise ValueError(
            f"A proportion of {p} results in 0 survivors, you've selected your population into extinction."
        )

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
    out_pop = cp.deepcopy(in_pop[out_idcs])
    return out_pop


def randsplit(in_pop, out_pop_size):
    # in_pop=cp.deepcopy(in_pop)
    inpopsize = in_pop.shape[0]
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


def old_randsplit(in_pop, out_pop_size):
    # in_pop=cp.deepcopy(in_pop)
    inpopsize = in_pop.shape[0]
    idcs_lina = np.random.choice(range(inpopsize), int(inpopsize / 2), replace=False)
    idcs_linb = np.array(
        [rand for rand in np.arange(inpopsize) if rand not in idcs_lina]
    )
    lina = grow_pop(in_pop[idcs_lina], out_pop_size, "equal")
    linb = grow_pop(in_pop[idcs_linb], out_pop_size, "equal")
    return lina, linb


def branch_evol(in_pop, ngens):
    in_pop = cp.deepcopy(in_pop)
    if in_pop.size:
        for gen in range(ngens):
            print(f"producing generation {gen}")
            survivors = select(in_pop, pf.prop_survivors, pf.select_strategy)
            next_pop = grow_pop(survivors, pf.pop_size, pf.reproductive_strategy)
            in_pop = next_pop
    return in_pop


def unpickle(filename):
    pickle_off = open(filename, "rb")
    output = pickle.load(pickle_off)
    return output


def bloop(organism_array, outfile_prefix="outfile"):
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


#run2=unpickle(filename)

def main_serial():
    founder = founder_miner(0.3)
    results_array = np.ndarray(13, dtype=object)
    founder_pop = grow_pop(founder, pf.pop_size, "equal")
    results_array[0] = cp.deepcopy(founder_pop)
    anc1_stem, anc2_stem = randsplit(founder_pop, pf.pop_size)
    # stem_lin3,stem_lin4=randsplit(founder_pop,pf.pop_size)
    results_array[1] = cp.deepcopy(anc1_stem)
    results_array[2] = cp.deepcopy(anc2_stem)
    # results_array[3]=cp.deepcopy(stem_lin3)
    # results_array[4]=cp.deepcopy(stem_lin4)
    anc_branches = np.array([anc1_stem, anc2_stem], dtype=object)
    genslist1 = np.array([10, 10])

    for i in range(len(anc_branches)):
        results_array[i + 3] = branch_evol(anc_branches[i], genslist1[i])

    # anc1_tip,anc2_tip=list(result)
    # results_array[3]=cp.deepcopy(anc1_tip)
    # results_array[4]=cp.deepcopy(anc2_tip)

    # results_array[7]=cp.deepcopy(tip_lin3)
    # results_array[8]=cp.deepcopy(tip_lin4)
    leafa_stem, leafb_stem = randsplit(results_array[3], pf.pop_size)
    results_array[5], results_array[6] = cp.deepcopy(leafa_stem), cp.deepcopy(
        leafb_stem
    )
    leafc_stem, leafd_stem = randsplit(results_array[4], pf.pop_size)
    results_array[7], results_array[8] = cp.deepcopy(leafc_stem), cp.deepcopy(
        leafd_stem
    )

    four_leaves = np.array(
        [leafa_stem, leafb_stem, leafc_stem, leafd_stem], dtype=object
    )
    genslist2 = np.array([10, 10, 10, 10])

    for i in range(len(four_leaves)):
        results_array[i + 9] = branch_evol(four_leaves[i], genslist2[i])

    # leafa_tip,leafb_tip,leafc_tip,leafd_tip=list(result)
    # results_array[9],results_array[10],results_array[11],results_array[12]=cp.deepcopy(leafa_tip),
    # cp.deepcopy(leafb_tip),cp.deepcopy(leafc_tip),cp.deepcopy(leafd_tip)
    return results_array


def main_parallel():
    founder = founder_miner(0.3)
    results_array = np.ndarray(13, dtype=object)
    founder_pop = grow_pop(founder, pf.pop_size, "equal")
    results_array[0] = cp.deepcopy(founder_pop)
    anc1_stem, anc2_stem = randsplit(founder_pop, pf.pop_size)
    # stem_lin3,stem_lin4=randsplit(founder_pop,pf.pop_size)
    results_array[1] = cp.deepcopy(anc1_stem)
    results_array[2] = cp.deepcopy(anc2_stem)
    # results_array[3]=cp.deepcopy(stem_lin3)
    # results_array[4]=cp.deepcopy(stem_lin4)
    anc_branches = np.array([anc1_stem, anc2_stem], dtype=object)
    genslist1 = np.array([10, 10])

    with ProcessPoolExecutor() as pool:
        result = pool.map(branch_evol, anc_branches, genslist1)

    anc1_tip, anc2_tip = list(result)
    results_array[3] = cp.deepcopy(anc1_tip)
    results_array[4] = cp.deepcopy(anc2_tip)

    # results_array[7]=cp.deepcopy(tip_lin3)
    # results_array[8]=cp.deepcopy(tip_lin4)
    leafa_stem, leafb_stem = randsplit(anc1_tip, pf.pop_size)
    results_array[5], results_array[6] = cp.deepcopy(leafa_stem), cp.deepcopy(
        leafb_stem
    )
    leafc_stem, leafd_stem = randsplit(anc2_tip, pf.pop_size)
    results_array[7], results_array[8] = cp.deepcopy(leafc_stem), cp.deepcopy(
        leafd_stem
    )

    four_leaves = np.array(
        [leafa_stem, leafb_stem, leafc_stem, leafd_stem], dtype=object
    )
    genslist2 = np.array([10, 10, 10, 10])

    with ProcessPoolExecutor() as pool:
        result = pool.map(branch_evol, four_leaves, genslist2)

    leafa_tip, leafb_tip, leafc_tip, leafd_tip = list(result)
    results_array[9], results_array[10], results_array[11], results_array[12] = (
        cp.deepcopy(leafa_tip),
        cp.deepcopy(leafb_tip),
        cp.deepcopy(leafc_tip),
        cp.deepcopy(leafd_tip),
    )
    return results_array

if __name__ == "__main__":
    result = main_serial()
    print("Analysis completed", result.shape)
    store(result)
