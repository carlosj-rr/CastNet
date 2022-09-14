from wsgiref import headers
import numpy as np

# translate an arbitrary sequence of 'integrified bases' into ACTG. Notice the input is an integer, not a string of integers.
def translate_int_seq(int_seq):
    if type(int_seq) == int:
        strint_seq=''.join(list(str(int_seq)))
    elif type(int_seq) == np.ndarray:
        strint_seq=''.join(list(int_seq.astype(str)))
    trans_dict={'1':'A','2':'C','3':'T','4':'G'}
    return(''.join(list(map((lambda x: trans_dict.get(x)),list(strint_seq)))))

def ali_saver(dset_prefix,pop_arr,tip_names,asints=True):
    num_seqs=(np.mean([x.shape[0] for x in pop_arr[:]])*0.04).astype(int)
    tot_lins=len(pop_arr)
    tot_genes=pop_arr[0][0][1].shape[0]
    for gene_idx in range(tot_genes):
        headers_list=[]
        sequences_list=[]
        gene_outfilename=dset_prefix+'_gene-'+str(gene_idx)+"_ali.fas"
        print(f"Output file is {gene_outfilename}")
        for pop_idx in range(len(pop_arr)):
            curr_pop=pop_arr[pop_idx]
            curr_popsize=len(curr_pop)
            rand_indivs=np.random.choice(range(curr_popsize),num_seqs)
            pop_prefix=">Lin"+tip_names[pop_idx]
            for indiv_idx in rand_indivs:
                seqname=pop_prefix+"-"+"org_"+str(indiv_idx)+"\n"
                with open(gene_outfilename,'w') as f:
                    f.write(seqname)
                headers_list.append(seqname)
                if not asints:
                    sequence=''.join(pop_arr[pop_idx][indiv_idx][1][gene_idx])+'\n'
                    sequences_list.append(sequence)
                else:
                    sequence=translate_int_seq(pop_arr[pop_idx][indiv_idx][1][gene_idx])+"\n"
                    sequences_list.append(sequence)
                with open(gene_outfilename,'w') as f:
                    for i in range(len(headers_list)):
                        f.write(headers_list[i])
                        f.write(sequences_list[i])   
    return