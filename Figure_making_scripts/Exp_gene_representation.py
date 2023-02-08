import numpy as np
import pickle
import seaborn as sns
import matplotlib.pylab as plt
def unpickle(filename):
	pickle_off = open(filename,'rb')
	return(pickle.load(pickle_off))

tips = unpickle("EvolRun_01-23-2023_17-05-20.pkl") # Set the correct filename for the pickled file
genes_on = np.ndarray((7,5),dtype=float)
for i in range(7):
	curr_row = np.mean([ x[8] for x in tips[i] ],axis=0)
	genes_on[i] = curr_row
	
ylabs = ["Ancestral\nPopulation","(A,B) Ancestor","(C,D) Ancestor","Lineage A","Lineage B","Lineage C", "Lineage D"]
xlabs = []
for i in range(5):
	xlabs.append("Gene "+str(i+1)) #'+1' added to have gene names match the trees on the figure.

sns.set(font_scale=0.6)
ax = sns.heatmap(genes_on, linewidth = 0.5, annot=genes_on, xticklabels=xlabs, yticklabels=ylabs, cmap = "crest")
ax.xaxis.tick_top()
plt.savefig("Heat_map.svg")