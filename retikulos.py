import numpy as np
import copy as cp
import scipy, os, time
import random
import params_file as pf
from scipy import stats
import math
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
        #print(genome)
        # Importing the values for producing all the regulatory information.
        prop_off=pf.prop_unlinked # thresholds and decays will have the converse of this probability as 0s. See blow.
        thresh_boundaries=pf.thresh_boundaries # tuple of 2 values.
        decay_boundaries=pf.decay_boundaries # tuple of 2 values.
        grn=makeGRN(n_genes,prop_off)
        thresholds=randomMaskedVector(n_genes,(1-prop_off),thresh_boundaries[0],thresh_boundaries[1])
        decays=randomMaskedVector(n_genes,(1-prop_off),decay_boundaries[0],decay_boundaries[1])
        # Importing values for the developmental info
        dev_steps=pf.dev_steps
        start_vect=(lambda x: np.array([1]*1+[0]*(x-1)))(n_genes)
        development=develop(start_vect,grn,decays,thresholds,dev_steps)
        genes_on=(development.sum(axis=0) != 0).astype(int)
        #print(f"Current fitness {fitness} is lower than minimum {min_fitness}")
        fitness=calcFitness(development)
        out_arr=np.array([np.array((n_generation,genome,proteome,grn,thresholds,decays,start_vect,development,genes_on,fitness),dtype=object)])
    return(out_arr)

def translate_codon(codon):
    idx=np.where(dna_codons == codon)[0][0]
    aminoac=trans_aas[idx]
    return(aminoac)

def makeGRN(numGenes,prop_unlinked):
    grn = randomMaskedVector(numGenes ** 2,prop_unlinked,pf.new_link_bounds[0],pf.new_link_bounds[1])
    grn = grn.reshape(numGenes,numGenes)
    return(grn)

def makeGenomeandProteome(seq_length,num_genes,dna_codons=dna_codons,trans_aas=trans_aas):
    if seq_length % 3:
#        print("Sequence length",seq_length,"is not a multiple of 3.")
        seq_length = seq_length - (seq_length % 3)
        num_codons = int(seq_length/3)
#        print("Rounding to", seq_length,"for",num_codons,"codons")
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

# Function that creates a vector of a given amount of values (within a given range), in which a certain proportion of the values are masked.
def randomMaskedVector(num_vals,prop_zero=0,min_val=0,max_val=1):
    if min_val > max_val:
        print("Error: minimum value greater than maximum value")
        return
    range_size = max_val - min_val
    if prop_zero == 0:
        rpv = np.array(range_size * np.random.random(num_vals) + min_val)
    else:
        mask = np.random.choice((0,1),num_vals,p=(prop_zero,1-prop_zero))
        rpv = np.array(range_size * np.random.random(num_vals) + min_val)
        rpv = (rpv * mask) + 0
    return(rpv)

def mutate_genome(old_gnome,old_prome,mut_coords):
    gnome=cp.deepcopy(old_gnome)
    prome=cp.deepcopy(old_prome)
    mut_num=mut_coords.shape[0] #get the number of rows in the mutation coordinate array, this is the number of mutations
    muttype_vect=np.ndarray((mut_num,2),dtype=object)
    for i in range(mut_num):
        coordinates=mut_coords[i,:]
        #print(coordinates)
        selected_gene=coordinates[0]
        selected_codon_from_gene=coordinates[1]
        selected_codpos=coordinates[2]
        #print((selected_gene,selected_codon_from_gene),selected_codpos)
        selected_codon=gnome[selected_gene,selected_codon_from_gene]
        prev_aacid=translate_codon(selected_codon)
        mutated_codon=pointMutateCodon(selected_codon,selected_codpos)
        print(f">>>>>>>>>>{selected_codon} 8===D** {mutated_codon}<<<<<<<<<<8===D") # BUG: THIS MAY BE SMTH ELSE - Mutated codons seem to end up being the same always
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
#    out_genome=gnome
#    out_proteome=prome
    return(gnome,prome,muttype_vect)

def codPos(muts,num_genes,num_codons):
    #base1=num+1
    out_array=np.ndarray((muts.size,3),dtype=object)
    gene_bps=num_codons*3
    genome_bps=gene_bps*num_genes
    genenum_array=np.ndarray((num_genes,gene_bps),dtype=object)
    for i in range(num_genes):
        genenum_array[i,:]=i
    genenum_array=genenum_array.flatten()
    #print("genenum_array:",genenum_array)
    codpos_array=np.tile([0,1,2],num_codons*num_genes)
    #print("codpos_array:",codpos_array)
    codnum_array=np.ndarray((num_genes,gene_bps),dtype=object)
    for i in range(num_genes):
        codnum_array[i,:]=np.repeat(range(num_codons),3)
    codnum_array=codnum_array.flatten()
    #print("codnum_array:",codnum_array)
    for i in range(muts.size):
        basenum=muts[i]
        mut_val=np.array([genenum_array[basenum],codnum_array[basenum],codpos_array[basenum]])
        out_array[i,:]=mut_val
    return(out_array)
    

def randomMutations(in_genome,mut_rateseq):
    total_bases=in_genome.size*3 #Each value in the genome is a codon, so the whole length (in nucleotides) is the codons times 3
    mutations=np.random.choice((0,1),total_bases,p=(1-mut_rateseq,mut_rateseq))
    m=np.array(np.where(mutations != 0)).flatten()
    if m.size:
        output=m
    else:
        output=False
    return(output)

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
    mutations=randomMutations(in_genome,mut_rateseq)
    #print(mutations)
    if np.any(mutations):
        mut_coords=codPos(mutations,in_genome.shape[0],in_genome.shape[1])
        #print(mut_coords)
        out_genome,out_proteome,mutlocs=mutate_genome(in_genome,in_proteome,mut_coords)
        out_grn,out_thresh,out_decs=regulator_mutator(in_grn,in_genes_on,in_decs,in_thresh,mutlocs)
        out_dev=develop(start_vect,out_grn,out_decs,out_thresh,pf.dev_steps)
        out_genes_on=(out_dev.sum(axis=0) != 0).astype(int)
        out_fitness=calcFitness(out_dev)
    else:
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

def pointMutateCodon(codon,pos_to_mutate):
    bases=("T","C","A","G")
    base=codon[pos_to_mutate]
    change = [x for x in bases if x != base]
    new_base = np.random.choice(change)
    split_codon=np.array(list(codon))
    split_codon[pos_to_mutate]=new_base
    new_codon="".join(split_codon)
    return(new_codon)

def weight_mut(value,scaler=0.01):
    val=abs(value) #Make sure value is positive
    if val == 0:
        '''For values at zero, simply get 1, and then modify it by the scale
        This is for activating thresholds that are 0.'''
        val=scaler/scaler
    scaled_val=val*scaler #scale the value
    newVal=value+np.random.uniform(-scaled_val,scaled_val) #add the scaled portion to the total value to get the final result.
    return(newVal)

def threshs_and_decs_mutator(in_thresh,in_dec,mutarr):
    #print(f"Input thresholds were: {in_thresh}")
    #print(f"Input decays were: {in_dec}")
    #print(f"Input mutarr was:\n{mutarr}")
    the_tuple=(in_thresh,in_dec) # make a tuple in which the threshold array is the first value, and the decays the second.
    # This will allow me to easily choose among them at the time of mutating, see within the for loop.
    num_genes=len(in_thresh) #get the number of genes from the amount of values in the thresholds array
    genes=mutarr[:,0] # get the genes to be mutated from the mutarray's 1st column
    #print(f"The array of genes to be mutated is:\n{genes}")
    for i in np.arange(len(genes)): #go through each gene, and decide randomly whether to make a threshold or a decay mutation in the gene.'''
        tuple_idx=np.random.choice((0,1))
        #print(f"Thresholds = 0, Decays = 1, Random choice was = {tuple_idx}")
        gene_num=genes[i] # extract specific gene number that has to be mutated. This maps to the thresh and dec arrays.
        #print(f"This means that gene {gene_num} will be mutated:\nValue {the_tuple[tuple_idx][gene_num]}")
        new_value=abs(weight_mut(the_tuple[tuple_idx][gene_num]))
        the_tuple[tuple_idx][gene_num]=new_value
        #print(f"...is now {new_value}")
    out_thresh,out_decs=(the_tuple[0],the_tuple[1])
    return(out_thresh,out_decs)

def regulator_mutator(in_grn,genes_on,in_dec,in_thresh,muttype_vect):
    curr_grn=cp.deepcopy(in_grn)
    curr_thr=cp.deepcopy(in_thresh)
    curr_genes_on=cp.deepcopy(genes_on)
    curr_dec=cp.deepcopy(in_dec)
    curr_muttype_vect=cp.deepcopy(muttype_vect)
    inactive_links=np.array(list(zip(np.where(curr_grn == 0)[0],np.where(curr_grn == 0)[1])))
    num_genes=pf.num_genes
    '''I'm adding here a section that decides if any of the mutations will go to the thresholds or the decays.
    If there are any changes that have to happen in the decays and/or thresholds, we can call their mutation
    function. Otherwise we can keep on going.'''
    prop=2/(2+num_genes**2) #proportion of mutable sites that are thresholds OR decays
    hits=np.nonzero(np.random.choice((0,1),len(muttype_vect),p=(1-prop,prop)))[0]
    if hits.size > 0:
        mutsarr=curr_muttype_vect[hits]
        #print(f"Sending mutations:\n{mutsarr} to decays/thresholds")
        out_threshs,out_decs=threshs_and_decs_mutator(in_thresh,in_dec,mutsarr)
        curr_muttype_vect=np.delete(curr_muttype_vect,hits,axis=0)
    else:
        out_threshs,out_decs=curr_thr,curr_dec
    if curr_muttype_vect.size > 0:
        for i in curr_muttype_vect:
            gene=i[0]
            mtype=i[1]
            #print(f"Gene {gene} has mutation type {mut_kinds[mtype]}")
            if mtype != 0: # For all non-KO mutations (i.e. synonymous, and non-synonymous)...
                # Check this block all the way down to the "<>" below
                if curr_genes_on[gene]: # If the gene is ON...
                    active_links=np.array(list(zip(np.nonzero(curr_grn)[0],np.nonzero(curr_grn)[1])))
                    #print(f"Gene {gene} is ON ({curr_genes_on[gene]}).")
                    actives_in_gene=np.concatenate((active_links[active_links[:,1] == gene,:],active_links[active_links[:,0] == gene,:]),axis=0) # get the gene's active links
                    #print(f"Gene {gene}'s active links are:\n{actives_in_gene}, and the gene's cells show:\n {curr_grn[:,gene]} \n and {curr_grn[gene,:]}")
                    #print(f"GRN is:\n{in_grn}")
                    if mtype == 1: # And the mutation is non-synonymous...
                        #print(f"Mutation {mtype} is NS")
                        #print(f"range to be used is range({len(actives_in_gene)})")
                        rand_idx=np.random.choice(np.arange(len(actives_in_gene))) # FIXED # get a random index number for mutating a link
                        coordinates=tuple(actives_in_gene[rand_idx,:]) # get the random link's specific coordinates
                        val=curr_grn[coordinates] # Extract the value that will be mutated.
                        curr_grn[coordinates]=weight_mut(val,0.5) # mutate the value.
                        #print(f"Mutating coordinate {coordinates} of the GRN, currently showing the value {val} to {in_grn[coordinates]}")
                    elif mtype == 2: # If gene is ON, and the mutation is synonymous...
                        #print(f"Mutation {mtype} is S")
                        #print(f"range to be used is range({len(actives_in_gene)})")
                        rand_idx=np.random.choice(np.arange(len(actives_in_gene))) # FIXED # Same as above
                        coordinates=tuple(actives_in_gene[rand_idx,:]) # Same as above
                        val=curr_grn[coordinates] # Same as above
                        curr_grn[coordinates]=weight_mut(val,0.001) # mutate the value by a very small amount.
                        #print(f"Mutating coordinate {coordinates} of the GRN a tiny little only, from {val} to {in_grn[coordinates]}")
                    else:
                        #print(f"Gene {gene} is neither on nor off, its state is {curr_genes_on[gene]}")
                        None
                #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<888>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                else: # If the gene is OFF...
                    #print(f"Gene {gene} is OFF ({curr_genes_on[gene]}).")
                    if mtype == 1: # And the mutation is non-synonymous
                        #print(f"And gene{gene}'s mutation is NS")
                        inactive_links=np.array(list(zip(np.where(curr_grn == 0)[0],np.where(curr_grn == 0)[1])))
                        inactives_in_gene=np.concatenate((inactive_links[inactive_links[:,1] == gene,:],inactive_links[inactive_links[:,0] == gene,:]),axis=0)
                        rand_idx=np.random.choice(np.arange(len(inactives_in_gene))) # FIXED # Same as above, but with inactives instead
                        coordinates=tuple(inactives_in_gene[rand_idx,:]) # Same as above above             
                        mean_exp_val=np.mean(np.abs(curr_grn[np.nonzero(curr_grn)])) # Mean expression amount
                        sign=np.random.choice((-1,1)) # Randomly choose between negative or positive
                        new_val=mean_exp_val*sign
                        #print(f"Flipping inactive value at coordinate {coordinates} on at level {new_val}")
                        in_grn[coordinates]=new_val
                    elif mtype == 2: # If gene is OFF, and the mutation is synonymous...
                        # check for all active links of the gene
                        active_links=np.array(list(zip(np.nonzero(curr_grn)[0],np.nonzero(curr_grn)[1])))
                        if active_links.size == 0: #If no links are active (such as in a gene that just got KO'd)...
                            all_links=np.array(list(zip(np.where(curr_grn == 0)[0],np.where(curr_grn == 0)[1]))) # Use the inactives to mutate
                        else: # Otherwise mutate any link from that gene (since it's off, all changes are synonymous)
                            inactive_links=np.array(list(zip(np.where(curr_grn == 0)[0],np.where(curr_grn == 0)[1])))
                            all_links=np.concatenate((active_links,inactive_links),axis=0)
                        #all_links=np.concatenate((inactive_links,active_links),axis=0)
                        actives_in_gene=np.concatenate((all_links[all_links[:,1] == gene,:],all_links[all_links[:,0] == gene,:]),axis=0) # get the gene's active links
                        #print(f"range to be used is range({len(actives_in_gene)})")
                        rand_idx=np.random.choice(np.arange(len(actives_in_gene))) # FIXED # get a random index number for mutating a link
                        coordinates=tuple(actives_in_gene[rand_idx,:]) # get the random link's specific coordinates
                        val=curr_grn[coordinates] # Extract the value that will be mutated.
                        #print(f"Mutating coordinate {coordinates} of the GRN, currently showing the value {val}")
                        curr_grn[coordinates]=weight_mut(val,0.5) # mutate the value.
                    else:
                        None            
            else: # If mutation is KO
                curr_grn[gene,:]=0
                curr_grn[:,gene]=0
                out_grn=curr_grn
                curr_genes_on[gene]=0 # Important change to avoid THE bug.
                #print(f"Knocking out gene{gene}: {curr_genes_on[gene]}.")
    else:
        pass
        #print("No mutations in this round")
    out_grn=cp.deepcopy(curr_grn)
    #print("Copying the input GRN deeply")
    #out_dev=develop(in_start_vect, out_grn,out_decs,out_threshs,pf.dev_steps)
    #out_genes_on=(out_dev.sum(axis=0) != 0).astype(int)
    #out_fitness=calcFitness(out_dev)
    return(out_grn,out_threshs,out_decs)

mut_kinds=np.array(["nonsense","non-synonymous","synonymous"])

def develop(start_vect,grn,decays,thresholds,dev_steps):
    start_vect = start_vect
#    print(f"Starting with vector: {start_vect}\n and thresholds {thresholds}")
    geneExpressionProfile = np.ndarray(((pf.dev_steps+1),pf.num_genes))
    geneExpressionProfile[0] = np.array([start_vect])
    #Running the organism's development, and outputting the results
    #in an array called geneExpressionProfile
    invect = start_vect
    counter=1
    for i in range(dev_steps):
#      print(f"Development step {counter}")
        decayed_invect = (lambda x, l: x*np.exp(-l))(invect,decays) # apply decay to all gene qties. previously: exponentialDecay(invect,decays)
#        print(f"Shapes of objects to be fed to matmul:\n{grn.shape}\t{decayed_invect.shape}")
        exp_change = np.matmul(grn,decayed_invect) #calculate the regulatory effect of the decayed values.
#        exp_change = myDotProd(grn,decayed_invect) #check my bootleg dot product function
#        print(f"Output of dot product:\n{exp_change}")
        pre_thresholds = exp_change + decayed_invect # add the decayed amounts to the regulatory effects
#        print(f"Result when added:\n{pre_thresholds}")
        thresholder = (pre_thresholds > thresholds).astype(int) # a vector to rectify the resulting values to their thresholds.
#        print(f"Threshold rectifier vector:\n{thresholder}")
        currV = pre_thresholds * thresholder # rectify with the thresholder vect. This step resulted in the deletion of the 'rectify()' function
 #       print(f"Rectifying with the thresholds gives:\n{currV}")
 #      currV = currV
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
        sim_to_exp = exponentialSimilarity(development)
        fitness_val = np.mean([genes_on,exp_stab,sim_to_exp])
    else:
        fitness_val = 0
    return(fitness_val)
def lastGeneExpressed(development,min_reproducin):
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

# Assumes input is a population (i.e. an array of organism arrays), it should crash if it doesn't find 2 dimensions.
def grow_pop(in_orgs,out_pop_size,strategy='equal'):
    in_orgs=cp.deepcopy(in_orgs)
    num_in_orgs=in_orgs.shape[0]
    orgs_per_org=np.array([np.floor_divide(out_pop_size,num_in_orgs)])
    #print(f"Orgs per org is {orgs_per_org}.")
    corr_pop_size=orgs_per_org*num_in_orgs
    #in_orgs=cp.deepcopy(in_orgs)
    #print(f"Making a population out of the {num_in_orgs} organisms given, reproductive strategy is {strategy}.\nEach organism will have {orgs_per_org[0]} offspring, for a total of {corr_pop_size[0]}.")
    if strategy == 'equal':
        orgs_per_org=np.repeat(orgs_per_org,num_in_orgs)
        #print(f"Offspring/organism array:\n{orgs_per_org}")
    elif strategy == 'fitness_linked':
        print("Reproduction is fitness bound.")
        pass
    else:
        print(f"Reproductive strategy {strategy} not recognized")
        raise ValueError("Invalid reproductive strategy")
    counter=0
    out_pop=np.ndarray((corr_pop_size[0],),dtype=object)
    for k in range(num_in_orgs): # taking each input organism and adding the requested offspring to the output population.
        num_offsp=orgs_per_org[k]
        for i in range(num_offsp):
            indiv=mutation_wrapper(in_orgs,pf.seq_mutation_rate)[0]
            out_pop[counter]=indiv
            #print(f"Producing organism #{counter}")
            counter=counter+1
            #print(np.all(out_pop[counter == out_pop[(counter-1)]]))
    out_pop=cleanup_deads(out_pop) # removing any dead organisms.
    return(out_pop)

def store(thing):
    now=datetime.now()
    moment=now.strftime("%m-%d-%Y_%H-%M-%S")
    filename="./EvolRun_"+moment+".pkl"
    with open(filename,'wb') as fh:
        pickle.dump(thing,fh)
    print(f"Your session results were saved in {filename}.")

def cleanup_deads(in_pop):
    in_pop=cp.deepcopy(in_pop)
    tot_orgs=in_pop.shape[0]
    fitnesses=np.array([ x[9] for x in in_pop[:] ])
    live_ones=np.nonzero(fitnesses)[0]
    #print(f"current population has {live_ones.size} organisms alive")
    if live_ones.size == tot_orgs:
        out_pop=in_pop
    elif live_ones.size != 0:
        #print(f"{live_ones.size} organisms are dead. Sorry for your loss...")
        out_pop=in_pop[live_ones]
    elif live_ones.size == 0:
        print(f"Your population went extinct. Sorry for your loss.")
        out_pop=np.array([])
    return(out_pop)

def select(in_pop,p=0.1,strategy='high pressure'):
    in_pop=cp.deepcopy(in_pop)
    pop_size=in_pop.shape[0]
    num_survivors=int(pop_size*p)
    if strategy == "high pressure":
        fitnesses=np.array([ x[9] for x in in_pop[:] ])
        out_idcs=np.argpartition(fitnesses,-num_survivors)[-num_survivors:] # returns the **indices** for the top 'num_survivors' fitnesses.
    elif strategy == "low pressure" and p < 0.5:
        fitnesses=np.array([ x[9] for x in in_pop[:] ])
        half=np.floor_divide(pop_size,2)
        top_half=np.argpartition(fitnesses,-half)[-half:]
        out_idcs=np.random.choice(top_half,num_survivors,replace=False)
    elif strategy == "low pressure" and p >= 0.5:
        print(f"Low pressure strategy is not recommended for offspring populations\nresulting from more than half the parental population\nDefaulting to total relaxation of selection...")
        out_idcs=np.random.choice(range(pop_size),num_survivors,replace=False)
    elif strategy == "totally relaxed":
        out_idcs=np.random.choice(range(pop_size),num_survivors,replace=False)
    #print(f"Out population will have indices {out_idcs}")
    out_pop=in_pop[out_idcs]
    return(out_pop)
    

def randsplit(in_pop,out_pop_size):
    in_pop=cp.deepcopy(in_pop)
    inpopsize=in_pop.shape[0]
    idcs_lina=np.random.choice(range(inpopsize),int(inpopsize/2),replace=False)
    idcs_linb=np.array([ rand for rand in np.arange(inpopsize) if rand not in idcs_lina])
    lina=grow_pop(in_pop,out_pop_size,'equal')
    linb=grow_pop(in_pop,out_pop_size,'equal')
    return(lina,linb)

def main():
    founder=founder_miner()
    results_array=np.ndarray(7,dtype=object)
    founder_pop=grow_pop(founder,pf.pop_size,'equal')
    results_array[0]=cp.deepcopy(founder_pop)
    stem_lin1,stem_lin2=randsplit(founder_pop,pf.pop_size)
    stem_lin3,stem_lin4=randsplit(founder_pop,pf.pop_size)
    results_array[1]=cp.deepcopy(stem_lin1)
    results_array[2]=cp.deepcopy(stem_lin2)
    two_branches=np.array([stem_lin1,stem_lin2,stem_lin3,stem_lin4],dtype=object)
    n_genslist1=np.array([10,10,10,10])

    with ProcessPoolExecutor() as pool:
        result = pool.map(branch_evol,two_branches,n_genslist1)
        
    tip_lin1,tip_lin2,tip_lin3,tip_lin4=np.array(list(result),dtype=object)
    results_array[3]=tip_lin1
    results_array[4]=tip_lin2
    results_array[5]=tip_lin3
    results_array[6]=tip_lin4
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

