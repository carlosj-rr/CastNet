# Authors: Carlos Javier Rivera Rivera and Đorđe Grbić
# License: GNU GPLv3

from wsgiref import headers
import numpy as np
#import networkx as nx
import matplotlib.pylab as plt
import gif
import pickle

def unpickle(filename):
    pickle_off = open(filename, "rb")
    output = pickle.load(pickle_off)
    return output

# translate an arbitrary sequence of 'integrified bases' into ACTG. Notice the input is an integer, not a string of integers.
def translate_int_seq(int_seq):
    if type(int_seq) == int:
        strint_seq=''.join(list(str(int_seq)))
    elif type(int_seq) == np.ndarray:
        strint_seq=''.join(list(int_seq.astype(str)))
    trans_dict={'1':'A','2':'C','3':'T','4':'G'}
    return(''.join(list(map((lambda x: trans_dict.get(x)),list(strint_seq)))))

def ali_saver(dset_prefix,pop_arr,tip_names,asints=True,num_seqs=None):
    if not num_seqs:
        #num_seqs=(np.mean([x.shape[0] for x in pop_arr[:]])*0.04).astype(int)
        num_seqs=1
    #tot_lins=len(pop_arr)
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

def alimaker(pickled_file,num):
    source_pops=unpickle(pickled_file)
    lowest_count=min([ x.shape[0] for x in source_pops ])
    rand_indivs=np.random.choice(range(lowest_count),10)
    for rep in range(len(rand_indivs)):
        rand_num=rand_indivs[rep]
        for gene_idx in range(tot_genes):
            outfilename="Ali_"+str(num)+"_gene"+str(gene_idx)+"-rep"+str(rep)+".fas"
            headers_list=[ ">"+x for x in name_list[:] ]
            seqs_list=[ translate_int_seq(x[rand_num][1][gene_idx]) for x in source_pops[:] ]
            for indiv_idx in range(len(headers_list)):
                with open(outfilename,"a") as fc:
                    fc.write(headers_list[indiv_idx])
                    fc.write("\n")
                    fc.write(seqs_list[indiv_idx])
                    fc.write("\n")
    return

def floating_range(start,stop,step):
        outlist=[]
        while start < stop:
            outlist.append(start)
            start+=step
        return(outlist)

#@gif.frame
#def draw_avg_grns(in_pop,gen_n):
    #thresholds=[0,2]
    #prefix=gen_n
    #mean_grn=np.mean([x[3].flatten() for x in in_pop[:]],axis=0).reshape(in_pop[0][3].shape)
    #n_genes=mean_grn.shape[0]
    #angle_range=floating_range(0,360,(360/n_genes))
    #trans_angles=[ abs(x-90) if x < 90 else 360-(x-90) for x in angle_range]
    #rad_trans=np.multiply(trans_angles,(np.pi/180))
    #coords=list(zip(np.cos(rad_trans),np.sin(rad_trans)))
    #pos_clockwise=dict(zip(range(n_genes),coords))
    #outfilename="AvgGRN_"+str(prefix)+".png"
    #G=nx.from_numpy_array(mean_grn)
    #lablist=list(map(lambda x,y: x+y,np.repeat("G",len(G.nodes)),np.array(list(range(len(G.nodes))),dtype=str)))
    #labdict=dict(zip(range(len(lablist)),lablist))
    #pos=nx.circular_layout(G)
    #edgewidth=np.array([ d['weight'] for (u,v,d) in G.edges(data=True)])
    #edge_colors=[ 'red' if x==-1 else 'green' for x in np.sign(edgewidth)]
    #fig = plt.figure(figsize=(12,12),dpi=300); plt.clf()
    #ax=fig.gca()
    #ax.set_title("Gen "+str(prefix),fontsize=16)
    #nx.draw_networkx_nodes(G,pos_clockwise, alpha=0.5, node_size=400*20, node_shape='o',node_color='blue')
    #nx.draw_networkx_edges(G,pos_clockwise,width=edgewidth*5, edge_color=edge_colors, arrows=True,arrowsize=100)
    #nx.draw_networkx_labels(G,pos_clockwise,labels=labdict, font_size=20, font_weight='bold')
    #return(fig)

@gif.frame
def plot_avg_developments(in_pop,prefix):
    mean_dev=np.mean([x[7].flatten() for x in in_pop[:]],axis=0).reshape(in_pop[0][7].shape)
    stde_dev=np.std([x[7].flatten() for x in in_pop[:]],axis=0).reshape(in_pop[0][7].shape)
    yerr=stde_dev/2
    corr_arr=np.array(list(map((lambda x: 0 if x < 0 else x),(yerr-mean_dev).flatten()))).reshape(mean_dev.shape)
    yerr_abv=yerr+corr_arr
    yerr_blow=yerr-corr_arr
    outfile=str(prefix)+"_mean_pop_dev.png"
    fig=plt.figure(figsize=(10,5))
    for i in range(mean_dev.shape[1]):
        genevals=mean_dev[:,i]
        #geneerrs=half_std[:,i]
        gene_yerrblow=yerr_blow[:,i]
        gene_yerrabv=yerr_abv[:,i]
        plt.errorbar(range(mean_dev.shape[0]),genevals,yerr=[gene_yerrblow,gene_yerrabv],capsize=5)
    plt.legend(list(map((lambda x,y: x+y),np.repeat("gene ",mean_dev.shape[1]),np.array(list(range(mean_dev.shape[1])),dtype=str))))
    plt.suptitle("Gen "+str(prefix),fontsize=16)
    #fig.savefig(outfile)
    #plt.close()
    return(fig)