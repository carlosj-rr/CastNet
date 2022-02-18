from multiprocessing.sharedctypes import Value
import numpy as np
import copy as cp
import scipy, os, time
import random
import params_file as pf
from scipy import stats
import math
import sys
import matplotlib.pyplot as plt
import pickle
from datetime import datetime,date
import cProfile #This is to benchmark the code. Recommended by Djole.
from concurrent.futures import ProcessPoolExecutor # for using multiple cores.

dna_codons=np.array(['ATA', 'ATC', 'ATT', 'ATG', 'ACA', 'ACC', 'ACG', 'ACT', 'AAC',
       'AAT', 'AAA', 'AAG', 'AGC', 'AGT', 'AGA', 'AGG', 'CTA', 'CTC',
       'CTG', 'CTT', 'CCA', 'CCC', 'CCG', 'CCT', 'CAC', 'CAT', 'CAA',
       'CAG', 'CGA', 'CGC', 'CGG', 'CGT', 'GTA', 'GTC', 'GTG', 'GTT',
       'GCA', 'GCC', 'GCG', 'GCT', 'GAC', 'GAT', 'GAA', 'GAG', 'GGA',
       'GGC', 'GGG', 'GGT', 'TCA', 'TCC', 'TCG', 'TCT', 'TTC', 'TTT',
       'TTA', 'TTG', 'TAC', 'TAT', 'TGC', 'TGT', 'TGG','TAA','TAG','TGA'], dtype=object)

trans_aas=np.array(['I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S',
       'S', 'R', 'R', 'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H',
       'Q', 'Q', 'R', 'R', 'R', 'R', 'V', 'V', 'V', 'V', 'A', 'A', 'A',
       'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G', 'S', 'S', 'S', 'S',
       'F', 'F', 'L', 'L', 'Y', 'Y', 'C', 'C', 'W','_','_','_'], dtype=object)

def founder_miner(min_fitness=0.6):
    fitness=0
    while fitness < min_fitness:
        # Importing values for producing the genomic sequences
        n_generation=0
        n_genes=pf.num_genes
        seq_len=pf.seq_length
        genome,proteome=makeGenomeandProteome(seq_len,n_genes,dna_codons,trans_aas)
        # Importing the values for producing all the regulatory information.
        prop_off=pf.prop_unlinked # thresholds and decays will have the converse of this probability as 0s. See blow.
        thresh_boundaries=pf.thresh_boundaries # tuple of 2 values.
        decay_boundaries=pf.decay_boundaries # tuple of 2 values.
        grn=makeGRN(n_genes,prop_off)
        thresholds=randomMaskedVector(n_genes,(1-prop_off),min(thresh_boundaries),max(thresh_boundaries))
        decays=randomMaskedVector(n_genes,(1-prop_off),min(decay_boundaries),max(decay_boundaries))
        # Importing values for the developmental info
        dev_steps=pf.dev_steps
        start_vect=(lambda x: np.array([1]*1+[0]*(x-1)))(n_genes)
        development=develop(start_vect,grn,decays,thresholds,dev_steps)
        genes_on=(development.sum(axis=0) != 0).astype(int)
        fitness=calcFitness(development)
        out_arr=np.array([np.array((n_generation,genome,proteome,grn,thresholds,decays,start_vect,development,genes_on,fitness),dtype=object)])
    return(out_arr)

def makeGenomeandProteome(seq_length,num_genes,dna_codons=dna_codons,trans_aas=trans_aas):
    if seq_length % 3:
        seq_length = seq_length - (seq_length % 3)
        num_codons = int(seq_length/3)
    else:
        num_codons=int(seq_length/3)
    idx_vect=np.array(range(0,len(dna_codons)-3))
    genome_arr=np.empty((num_genes,num_codons),dtype=object)
    proteome_arr=np.empty((num_genes,num_codons),dtype=object)
    for i in range(0,num_genes):
        rand_codon_idx=np.hstack((np.random.choice(idx_vect,(num_codons-1)),np.random.choice((61,62,63),1)))
        #len(rand_codons)
        genome_arr[i]=np.array(dna_codons[rand_codon_idx])
        proteome_arr[i]=np.array(trans_aas[rand_codon_idx])
    return(genome_arr,proteome_arr)

def makeGRN(numGenes,prop_unlinked):
    grn = randomMaskedVector(numGenes ** 2,prop_unlinked,pf.new_link_bounds[0],pf.new_link_bounds[1])
    grn = grn.reshape(numGenes,numGenes)
    return(grn)

# Function that creates a vector of a given amount of values (within a given range), in which a certain proportion of the values are masked.
def randomMaskedVector(num_vals,prop_zero=0,min_val=0,max_val=1):
    if min_val > max_val:
        raise ValueError(f"Minimum value {min_val} is larger than the maximum value {max_val}.\nConsider revising the function call to randomMaskedVector()")
    range_size = max_val - min_val
    if prop_zero == 0:
        rpv = np.array(range_size * np.random.random(num_vals) + min_val)
    else:
        mask = np.random.choice((0,1),num_vals,p=(prop_zero,1-prop_zero))
        rpv = np.array(range_size * np.random.random(num_vals) + min_val)
        rpv = (rpv * mask) + 0
    return(rpv)

#mut_kinds=np.array(["nonsense","non-synonymous","synonymous"])

def develop(start_vect,grn,decays,thresholds,dev_steps):
    start_vect = cp.deepcopy(start_vect)
    geneExpressionProfile = np.ndarray(((pf.dev_steps+1),pf.num_genes))
    geneExpressionProfile[0] = np.array([start_vect])
    #Running the organism's development, and outputting the results
    #in an array called geneExpressionProfile
    invect = start_vect
    counter=1
    for i in range(dev_steps):
        decayed_invect = (lambda x, l: x*np.exp(-l))(invect,decays) # apply decay to all gene qties. previously: exponentialDecay(invect,decays)
        exp_change = np.matmul(grn,decayed_invect) #calculate the regulatory effect of the decayed values.
        pre_thresholds = exp_change + decayed_invect # add the decayed amounts to the regulatory effects
        thresholder = (pre_thresholds > thresholds).astype(int) # a vector to rectify the resulting values to their thresholds.
        currV = pre_thresholds * thresholder # rectify with the thresholder vect. This step resulted in the deletion of the 'rectify()' function
        geneExpressionProfile[(i+1)] = currV
        invect = currV
        counter=counter+1
    return(geneExpressionProfile)

def calcFitness(development):
    min_reproducin = pf.min_reproducin
    is_alive = lastGeneExpressed(development,min_reproducin)
    if is_alive:
        genes_on = propGenesOn(development)
        exp_stab = expressionStability(development)
        sim_to_exp = 1-exponentialSimilarity(development) #added "1 -" as I realized that the R^2 value would approac 1 the more it assimilated an exponential function.
        fitness_val = np.mean([genes_on,exp_stab,sim_to_exp])
    else:
        fitness_val = 0
    return(fitness_val)
def lastGeneExpressed(development,min_reproducin): # is the last gene ever expressed above 'min_reproducin' level? AND is it expressed above 0 in the last developmental step?
    dev_steps,num_genes = development.shape
    last_col_bool = development[:,(num_genes - 1)] > min_reproducin
    last_val_last_col = development[dev_steps - 1, (num_genes - 1)]
    if last_col_bool.any() and last_val_last_col > 0:
        return_val = True
    else:
        return_val = False
    return(return_val)
def propGenesOn(development):
    genes_on = development.sum(axis=0) > 0
    return(genes_on.mean())
def expressionStability(development):  # I haven't thought deeply about this.
    row_sums = development.sum(axis=1)# What proportion of the data range is
    stab_val = row_sums.std() / (row_sums.max() - row_sums.min()) # the stdev? Less = better
    return(stab_val)
def exponentialSimilarity(development):
    dev_steps,num_genes = development.shape
    row_means = development.mean(axis=1)
    tot_dev_steps = dev_steps
    fitted_line = scipy.stats.linregress(range(tot_dev_steps),np.log(row_means))
    r_squared = fitted_line.rvalue ** 2
    return(r_squared)

class CodonError(Exception):
    pass

def translate_codon(codon):
    if np.where(dna_codons == codon)[0].size != 0:
        idx=np.where(dna_codons == codon)[0][0]
        aminoac=trans_aas[idx]
    else:
        raise CodonError(f"<{codon}> is NOT a valid codon sequence. Please make sure it is one of the following:\n{dna_codons}")
    return(aminoac)

# Assumes input is a population (i.e. an array of organism arrays), it should crash if it doesn't find 2 dimensions.
def grow_pop(in_orgs,out_pop_size,strategy='equal'):
    in_orgs=cp.deepcopy(in_orgs)
    num_in_orgs=in_orgs.shape[0]
    orgs_per_org=np.array([np.round(out_pop_size/num_in_orgs).astype(int)])
    curr_pop_size=orgs_per_org*num_in_orgs
    if strategy == 'equal':
        orgs_per_org=np.repeat(orgs_per_org,num_in_orgs)
    elif strategy == 'fitness_linked':
        print("Reproduction is fitness bound.")
        pass
    else:
        raise ValueError(f"Reproductive strategy: <{strategy}> not recognized")
    counter=0
    out_pop=np.ndarray((curr_pop_size[0],),dtype=object)
    #print(f"Shape of output population is {out_pop.shape}\nOrganisms per organisms are {orgs_per_org }")
    for k in range(num_in_orgs): # taking each input organism and adding the requested offspring to the output population.
        num_offsp=orgs_per_org[k]
        #print(f"The current organism will produce {num_offsp} offspring")
        for i in range(num_offsp):
            indiv=mutation_wrapper(in_orgs,pf.seq_mutation_rate)[0]
            out_pop[counter]=indiv
            print(f"Produced organism #{counter}")
            counter=counter+1
            #print(np.all(out_pop[counter == out_pop[(counter-1)]]))
    out_pop=cleanup_deads(out_pop) # removing any dead organisms.
    return(out_pop)

class NegativeIndex(Exception):
    pass

# Input is an organism array, as produced by the founder_miner() function, and the mutation rate of the nucleotide sequence (i.e. mutation probability per base).
def mutation_wrapper(orgarr,mut_rateseq):
    orgarrcp=cp.deepcopy(orgarr[0])
    in_gen_num=orgarrcp[0]
    in_genome=orgarrcp[1]
    in_proteome=orgarrcp[2]
    in_grn=orgarrcp[3]
    in_thresh=orgarrcp[4]
    in_decs=orgarrcp[5]
    start_vect=orgarrcp[6]
    in_dev=orgarrcp[7]
    in_genes_on=(in_dev.sum(axis=0) != 0).astype(int)
    in_fitness=orgarrcp[9]
    mutations=randomMutations(in_genome.size,mut_rateseq)
    if np.any(mutations):
        mut_coords=codPos(mutations,in_genome.shape)
        out_genome,out_proteome,mutlocs=mutate_genome(in_genome,in_proteome,mut_coords)
        out_grn,out_thresh,out_decs=regulator_mutator(in_grn,in_genes_on,in_decs,in_thresh,mutlocs)
        out_dev=develop(start_vect,out_grn,out_decs,out_thresh,pf.dev_steps)
        out_genes_on=(out_dev.sum(axis=0) != 0).astype(int)
        out_fitness=calcFitness(out_dev)
    else:
        print("No mutations this time.")
        out_genome=in_genome
        out_proteome=in_proteome
        out_grn=in_grn
        out_thresh=in_thresh
        out_decs=in_decs
        out_dev=in_dev
        out_genes_on=(out_dev.sum(axis=0) != 0).astype(int)
        out_fitness=in_fitness
    out_gen_num=in_gen_num+1
    out_org=np.array([[out_gen_num,out_genome,out_proteome,out_grn,out_thresh,out_decs,start_vect,out_dev,out_genes_on,out_fitness]],dtype=object)
    return(out_org)

def randomMutations(genome_size,mut_rateseq): # genome_size is in CODONS
    total_bases=genome_size*3 #Each value in the genome is a codon, so the whole length (in nucleotides) is the codons times 3.
    mutations=np.random.choice((0,1),total_bases,p=(1-mut_rateseq,mut_rateseq))
    m=np.array(np.where(mutations != 0)).flatten()
    if m.size:
        output=m
    else:
        output=False
    return(output)

def codPos(muts,gnome_shape):
    #base1=num+1
    num_genes=gnome_shape[0]
    num_codons=gnome_shape[1]
    if (np.any(muts < 0)):
        raise NegativeIndex(f"There are negative index values in your mutation indices:\n{muts}.\n This will result in untractable mutations.\nConsder double-checking the result of randomMutations()")
    out_array=np.ndarray((muts.size,3),dtype=object)
    gene_bps=num_codons*3
    genenum_array=np.ndarray((num_genes,gene_bps),dtype=object)
    for i in range(num_genes):
        genenum_array[i,:]=i
    genenum_array=genenum_array.flatten()
    codpos_array=np.tile([0,1,2],num_codons*num_genes)
    codnum_array=np.ndarray((num_genes,gene_bps),dtype=object)
    for i in range(num_genes):
        codnum_array[i,:]=np.repeat(range(num_codons),3)
    codnum_array=codnum_array.flatten()
    for i in range(muts.size):
        basenum=muts[i]
        mut_val=np.array([genenum_array[basenum],codnum_array[basenum],codpos_array[basenum]])
        out_array[i,:]=mut_val
    return(out_array)

def mutate_genome(old_gnome,old_prome,mut_coords):
    gnome=cp.deepcopy(old_gnome)
    prome=cp.deepcopy(old_prome)
    mut_num=mut_coords.shape[0] #get the number of rows in the mutation coordinate array, this is the number of mutations
    muttype_vect=np.ndarray((mut_num,2),dtype=object)
    if (np.any(mut_coords < 0)):
        raise NegativeIndex(f"Some indices in the mutation coordinates are negative:\n{mut_coords}\nThis may result in untractable mutations.\nConsider examining the output of codPos().")
    for i in range(mut_num):
        coordinates=mut_coords[i,:]
        selected_gene=coordinates[0]
        selected_codon_from_gene=coordinates[1]
        selected_codpos=coordinates[2]
        selected_codon=gnome[selected_gene,selected_codon_from_gene]
        prev_aacid=translate_codon(selected_codon)
        mutated_codon=pointMutateCodon(selected_codon,selected_codpos)
        gnome[selected_gene,selected_codon_from_gene]=mutated_codon
        new_aacid=translate_codon(mutated_codon)
        if prev_aacid == new_aacid: #Synonymous mutations are plotted as '2'
            muttype=2
        elif new_aacid == "_": # Nonsense mutations are plotted as '0'
            muttype=0
        else: # Nonsynonymous mutations are plotted as '1'
            muttype=1
        prome[selected_gene,selected_codpos]=new_aacid
        muttype_vect[i]=(selected_gene,muttype)
    out_genome=gnome
    out_proteome=prome
    return(out_genome,out_proteome,muttype_vect)

def pointMutateCodon(codon,pos_to_mutate):
    if pos_to_mutate < 0:
        raise NegativeIndex(f"Codon position {pos_to_mutate} is negative\nThis can cause intractable mutations\nConsider verifying the output of codPos()")
    bases=("T","C","A","G")
    base=codon[pos_to_mutate]
    change = [x for x in bases if x != base]
    new_base = np.random.choice(change)
    split_codon=np.array(list(codon))
    split_codon[pos_to_mutate]=new_base
    new_codon="".join(split_codon)
    return(new_codon)

class MutationTypeError(Exception):
    pass

def regulator_mutator(in_grn,genes_on,in_dec,in_thresh,muttype_vect):
    curr_grn=cp.deepcopy(in_grn)
    curr_thr=cp.deepcopy(in_thresh)
    curr_genes_on=cp.deepcopy(genes_on)
    curr_dec=cp.deepcopy(in_dec)
    curr_muttype_vect=cp.deepcopy(muttype_vect)
    if np.any(curr_muttype_vect < 0):
        raise NegativeIndex(f"There is a number in the mutations array that is negative:\n\n{curr_muttype_vect}\n\nThis may result in untractable mutations, consider inspecting the output of codPos()")
    inactive_links=np.array(list(zip(np.where(curr_grn == 0)[0],np.where(curr_grn == 0)[1])))
    num_genes=curr_genes_on.size
    # Choosing mutations for the thresholds or decays...
    prop=(lambda x: 2/(2+x))(num_genes) #proportion of total regulatory interactions that are thresholds OR decays (simplified from "2N/(2N+N^2)", where N is the total number of genes.)
    hits=np.nonzero(np.random.choice((0,1),len(muttype_vect),p=(1-prop,prop)))[0]
    if hits.size > 0:
        mutsarr=curr_muttype_vect[hits]
        out_threshs,out_decs=threshs_and_decs_mutator(in_thresh,in_dec,mutsarr)
        curr_muttype_vect=np.delete(curr_muttype_vect,hits,axis=0)
    else:
        out_threshs,out_decs=curr_thr,curr_dec
    if curr_muttype_vect.size > 0:
        refmat=np.repeat(1,curr_grn.size).reshape(curr_grn.shape)
        refmat[:]=curr_genes_on
        print(refmat)
        for i in curr_muttype_vect:
            gene=i[0]
            if gene not in range(num_genes):
                raise IndexError(f"Gene number {gene} is not within the number of genes in the system (integers between 0 and {num_genes-1})")
            mtype=i[1]
            exprstate=curr_genes_on[gene]
            #print(f"Expected: Gene {gene}\tExpression {exprstate}\tMutation type {mtype}.\n")
            if mtype in [1,2]:
                if exprstate == 0 and mtype == 1: # Gene OFF, Non-Synonymous mutation (uses active sites)
                    active_sites=np.array(list(zip(np.where(refmat == 1)[0],np.where(refmat == 1)[1])))  # non-silent sites for non-synonymous mutations
                    #print("Observed: Gene X\tExpression 0\tMutation type 1~~~~\n<<<<>>>>\n")
                    if active_sites.size == 0: #if in the same generation, all genes of the organism have been KO'd...
                        #print("All genes have been turned OFF, NS mutation will be on any link from the gene, and it will turn it ON.")
                        mutable_sites=np.array(list(zip(np.repeat(gene,num_genes),range(num_genes))))
                    else:
                        #print("There are enough NS sites to mutate...As you were!")
                        mutable_sites=active_sites[np.where(active_sites[:,0] == gene)]
                    site_to_mutate=tuple(mutable_sites[np.random.choice(range(mutable_sites.shape[0]))])
                    curr_grn[site_to_mutate]=weight_mut(in_grn[site_to_mutate],0.5)
                    #print(f"value {in_grn[site_to_mutate]} mutated into {curr_grn[site_to_mutate]}.")
                if exprstate == 0 and mtype == 2: #Gene OFF, Synonymous mutation (uses inactive sites)
                    #print("Observed: Gene X\tExpression 0\tMutation type 2~~~~\n<<<<>>>>\n")
                    #print("Using the inactives site matrix...")
                    inactive_sites=np.array(list(zip(np.where(refmat == 0)[0],np.where(refmat == 0)[1]))) #silent sites for synonymous mutations
                    #print(inactive_sites[0:10])
                    gene_col=inactive_sites[np.where(inactive_sites[:,1] == gene)]
                    gene_row=inactive_sites[np.where(inactive_sites[:,0] == gene)]
                    mutable_sites=np.row_stack((gene_col,gene_row))
                    site_to_mutate=tuple(mutable_sites[np.random.choice(range(mutable_sites.shape[0]))])
                    curr_grn[site_to_mutate]=weight_mut(in_grn[site_to_mutate],0.5)
                    #print(f"value {in_grn[site_to_mutate]} mutated into {curr_grn[site_to_mutate]}.")
                if exprstate == 1 and mtype == 1: #Gene ON, Non-Synonymous mutation (uses active sites)
                    #print("Observed: Gene X\tExpression 1\tMutation type 1~~~~\n<<<<>>>>\n")
                    active_sites=np.array(list(zip(np.where(refmat == 1)[0],np.where(refmat == 1)[1])))  # non-silent sites for non-synonymous mutations
                    #print(f"Gene {gene} is ON (exprstate={exprstate}, AND the mutation typs is NS (Mtype={mtype})")
                    gene_col=active_sites[np.where(active_sites[:,1] == gene)]
                    gene_row=active_sites[np.where(active_sites[:,0] == gene)]
                    mutable_sites=np.row_stack((gene_col,gene_row))
                    site_to_mutate=tuple(mutable_sites[np.random.choice(range(mutable_sites.shape[0]))])
                    #print(site_to_mutate)
                    #print(in_grn[site_to_mutate])
                    curr_grn[site_to_mutate]=weight_mut(in_grn[site_to_mutate],0.5)
                    #print(f"value {in_grn[site_to_mutate]} mutated into {curr_grn[site_to_mutate]}.")
                if exprstate == 1 and mtype == 2:#Gene OFF, Synonymous mutation (uses inactive sites)
                    #print("Observed: Gene X\tExpression 1\tMutation type 2~~~~\n<<<<>>>>\n")
                    #print(f"Gene {gene} is ON (exprstate={exprstate}, AND the mutation typs is S (Mtype={mtype}")
                    #print("Using the inactives site matrix...")
                    inactive_sites=np.array(list(zip(np.where(refmat == 0)[0],np.where(refmat == 0)[1]))) #silent sites for synonymous mutations
                    #print(inactive_sites[0:10])
                    if inactive_sites.size == 0: # If all genes are on, the 'gene_row' array will come up empty, so it switches to mutating any gene, a tiny little bit.
                        #print("All genes are on, synonymous mutation will be in an active link, but it will be very tiny...")
                        mutable_sites=np.array(list(zip(np.repeat(gene,num_genes),range(num_genes))))
                        site_to_mutate=tuple(mutable_sites[np.random.choice(range(mutable_sites.shape[0]))])
                        curr_grn[site_to_mutate]=weight_mut(in_grn[site_to_mutate],0.001)
                    else:
                        #print("there are more than 0 genes off, phew!!")
                        gene_row=inactive_sites[np.where(inactive_sites[:,0] == gene)]
                        mutable_sites=gene_row
                        site_to_mutate=tuple(mutable_sites[np.random.choice(range(mutable_sites.shape[0]))])
                        curr_grn[site_to_mutate]=weight_mut(in_grn[site_to_mutate],0.5)
            elif mtype == 0: # If mutation is KO
                #print(f"Observed: Gene X\tExpression X\tMutation 0~~~~\n<<<<>>>>\n")
                curr_grn[gene,:]=0
                curr_grn[:,gene]=0
                out_grn=curr_grn
                curr_genes_on[gene]=0 # Important change so the refmat is updated.
                #print(curr_genes_on)
                refmat[:]=curr_genes_on # Important change to avoid being unable to find synonymous mutations later on.
            else:
                raise MutationTypeError(f"Gene{gene}'s mutation type <{mtype}> is unclear,\nit must be one of [0,1,2].")
    out_grn=curr_grn
    return(out_grn,out_threshs,out_decs)


def threshs_and_decs_mutator(in_thresh,in_dec,mutarr):
    in_thresh=cp.deepcopy(in_thresh)
    in_dec=cp.deepcopy(in_dec)
    the_tuple=(in_thresh,in_dec) # make a tuple in which the threshold array is the first value, and the decays the second.
    # This will allow me to easily choose among them at the time of mutating, see within the for loop.
    genes=mutarr[:,0] # get the genes to be mutated from the mutarray's 1st column
    for i in np.arange(len(genes)): #go through each gene, and decide randomly whether to make a threshold or a decay mutation in the gene.
        tuple_idx=np.random.choice((0,1))
        gene_num=genes[i] # extract specific gene number that has to be mutated. This maps to the thresh and dec arrays.
        isit_ko=mutarr[i,1] == 0
        if isit_ko:
            new_value=0
        else:
            new_value=abs(weight_mut(the_tuple[tuple_idx][gene_num]))
        the_tuple[tuple_idx][gene_num]=new_value
    out_thresh,out_decs=(the_tuple[0],the_tuple[1])
    return(out_thresh,out_decs)

def weight_mut(value,scaler=0.01):
    val=abs(value) #Make sure value is positive
    if val == 0:
        '''For 0, simply get 1, and then modify it by the scale
        This is for activating values that are off.'''
        val=scaler/scaler
    scaled_val=val*scaler #scale the value
    newVal=value+np.random.uniform(-scaled_val,scaled_val) #add the scaled portion to the total value to get the final result.
    return(newVal)

def store(thing):
    now=datetime.now()
    moment=now.strftime("%m-%d-%Y_%H-%M-%S")
    filename="./EvolRun_"+moment+".pkl"
    with open(filename,'wb') as fh:
        pickle.dump(thing,fh)
    print(f"Your session results were saved in {filename}.")

def cleanup_deads(in_pop):
    #in_pop=cp.deepcopy(in_pop)
    tot_orgs=in_pop.shape[0]
    fitnesses=np.array([ x[9] for x in in_pop[:] ])
    live_ones=np.nonzero(fitnesses)[0]
    if live_ones.size == tot_orgs:
        out_pop=cp.deepcopy(in_pop)
    elif live_ones.size != 0:
        #print(f"{tot_orgs - live_ones.size} organisms are dead. Sorry for your loss...")
        out_pop=cp.deepcopy(in_pop[live_ones])
    elif live_ones.size == 0:
        print(f"Your population went extinct. Sorry for your loss.")
        out_pop=np.array([])
    return(out_pop)
class UnknownSelectiveStrategy(Exception):
    pass

def select(in_pop,p=0.1,strategy='high pressure'):
    pop_size=in_pop.shape[0]
    num_survivors=int(pop_size*p)
    fitnesses=np.array([ x[9] for x in in_pop[:] ])
    if num_survivors == 0:
        raise ValueError(f"A proportion of {p} results in 0 survivors, you've selected your population into extinction.")
    print(f"Proportion of survivors will be {p}, which will be {num_survivors}, out of a total of {pop_size}") # DEBUG
    if strategy == "high pressure":
        out_idcs=np.argpartition(fitnesses,-num_survivors)[-num_survivors:] # returns the **indices** for the top 'num_survivors' fitnesses.
    elif strategy == "low pressure" and p < 0.5:
        half=int(pop_size/2)
        top_half=np.argpartition(fitnesses,-half)[-half:]
        out_idcs=np.random.choice(top_half,num_survivors,replace=False)
    elif strategy == "low pressure" and p >= 0.5:
        out_idcs=np.random.choice(range(pop_size),num_survivors,replace=False)
    elif strategy == "totally relaxed":
        out_idcs=np.random.choice(range(pop_size),num_survivors,replace=False)
    else:
        raise UnknownSelectiveStrategy(f"Selective strategy {strategy} is not recognized\nThis value should be any of \"high pressure\", \"low pressure\", or \"totally relaxed\".\n Please double-check your input.")
    out_pop=cp.deepcopy(in_pop[out_idcs])
    return(out_pop)
    
def randsplit(in_pop,out_pop_size):
    #in_pop=cp.deepcopy(in_pop)
    inpopsize=in_pop.shape[0]
    idcs_lina=np.random.choice(range(inpopsize),int(inpopsize/2),replace=False)
    idcs_linb=np.array([ rand for rand in np.arange(inpopsize) if rand not in idcs_lina])
    print(f"The first random subselection of indices is of size {idcs_lina.size}, and the second of {idcs_linb.size}.")
    print(f"Do they share any number whatsoever?:\n{np.any(idcs_lina == idcs_linb)}")
    print(f"Output populations should be of {out_pop_size} individuals.")
    #lina=grow_pop(in_pop[idcs_lina],out_pop_size,'equal')
    #linb=grow_pop(in_pop[idcs_linb],out_pop_size,'equal')
    return(lina,linb)

def main(founder):
    founder=cp.deepcopy(founder)
    results_array=np.ndarray(9,dtype=object)
    founder_pop=grow_pop(founder,pf.pop_size,'equal')
    results_array[0]=cp.deepcopy(founder_pop)
    stem_lin1,stem_lin2=randsplit(founder_pop,pf.pop_size)
    stem_lin3,stem_lin4=randsplit(founder_pop,pf.pop_size)
    results_array[1]=cp.deepcopy(stem_lin1)
    results_array[2]=cp.deepcopy(stem_lin2)
    results_array[3]=cp.deepcopy(stem_lin3)
    results_array[4]=cp.deepcopy(stem_lin4)
    four_branches=np.array([stem_lin1,stem_lin2,stem_lin3,stem_lin4],dtype=object)
    n_genslist1=np.array([100,100,100,100])

    with ProcessPoolExecutor() as pool:
        result = pool.map(branch_evol,four_branches,n_genslist1)
        
    tip_lin1,tip_lin2,tip_lin3,tip_lin4=np.array(list(result),dtype=object)
    results_array[5]=cp.deepcopy(tip_lin1)
    results_array[6]=cp.deepcopy(tip_lin2)
    results_array[7]=cp.deepcopy(tip_lin3)
    results_array[8]=cp.deepcopy(tip_lin4)
    if False:
        stem_lin3,stem_lin4=randsplit(tip_lin1,pf.pop_size)
        results_array[5],results_array[6]=cp.deepcopy(stem_lin3),cp.deepcopy(stem_lin4)
        stem_lin5,stem_lin6=randsplit(tip_lin2,pf.pop_size)
        results_array[7],results_array[8]=cp.deepcopy(stem_lin5),cp.deepcopy(stem_lin6)
        
        four_branches=np.array([stem_lin3,stem_lin4, stem_lin5, stem_lin6],dtype=object)
        n_genslist2=np.array([10,10,10,10])
        
        with ProcessPoolExecutor() as pool:
            result = pool.map(branch_evol,four_branches,n_genslist2)
            
        tip_lin3,tip_lin4,tip_lin5,tip_lin6=np.array(list(result),dtype=object)
        results_array[9],results_array[10],results_array[11],results_array[12]=cp.deepcopy(tip_lin3),cp.deepcopy(tip_lin4),cp.deepcopy(tip_lin5),cp.deepcopy(tip_lin6)
    return(results_array)

def branch_evol(in_pop,ngens):
    in_pop=cp.deepcopy(in_pop)
    if in_pop.size:
        for gen in np.arange(ngens):
            print(f"producing generation {gen}")
            survivors=select(in_pop,pf.prop_survivors,pf.select_strategy)
            next_pop=grow_pop(survivors,pf.pop_size,pf.reproductive_strategy)
            in_pop=next_pop
    else:
        pass
    return(in_pop)

def unpickle(filename):
    pickle_off=open(filename,'rb')
    output=pickle.load(pickle_off)
    return(output)

if __name__ == "__main__":
    result=main()
#print(result.shape)
#store(result)

