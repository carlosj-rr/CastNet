from wsgiref import headers
import numpy as np
import networkx as nx
import matplotlib.pylab as plt

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

def draw_avg_grns(in_pop,prefix):
    thresholds=[0.001,0.01,0.1,1,1.001,1.01,1,1]
    mean_grn=np.mean([x[3].flatten() for x in in_pop[:]],axis=0).reshape(in_pop[0][3].shape)
    for i in thresholds:
        outfilename="AvgGRN_"+str(i)+"_thre_"+prefix+".png"
        newnet=np.sign((np.abs(mean_grn) > i) * mean_grn)
        G=nx.from_numpy_array(newnet)
        G=nx.relabel_nodes(G,dict(zip(range(len(G.nodes())),range(len(G.nodes)))))
        G=nx.drawing.nx_agraph.to_agraph(G)
        G.node_attr.update(color="green",style="empty")
        G.edge_attr.update(color="red",width="100.0")
        #plt.title(prefix+" mean GRN with threshold: "+str(i))
        G.draw(outfilename,format='png',prog='neato')
    return
 
def plot_avg_developments(in_pop,prefix):
    mean_dev=np.mean([x[7].flatten() for x in in_pop[:]],axis=0).reshape(in_pop[0][7].shape)
    stde_dev=np.std([x[7].flatten() for x in in_pop[:]],axis=0).reshape(in_pop[0][7].shape)
    yerr=stde_dev/2
    corr_arr=np.array(list(map((lambda x: 0 if x < 0 else x),(yerr-mean_dev).flatten()))).reshape(mean_dev.shape)
    yerr_abv=yerr+corr_arr
    yerr_blow=yerr-corr_arr
    outfile=prefix+"_mean_pop_dev.png"
    fig=plt.figure(figsize=(10,5))
    for i in range(mean_dev.shape[1]):
        genevals=mean_dev[:,i]
        #geneerrs=half_std[:,i]
        gene_yerrblow=yerr_blow[:,i]
        gene_yerrabv=yerr_abv[:,i]
        plt.errorbar(range(mean_dev.shape[0]),genevals,yerr=[gene_yerrblow,gene_yerrabv],capsize=5)
        
    plt.legend(list(map((lambda x,y: x+y),np.repeat("gene ",mean_dev.shape[1]),np.array(list(range(mean_dev.shape[1])),dtype=str))))
    plt.suptitle(prefix+" mean pop dev")
    fig.savefig(outfile)
    plt.close()
    return