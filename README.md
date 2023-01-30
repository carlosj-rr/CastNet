# CastNet
In order to use CastNet, you need to have python 3.9+ installed, as well as the following packages:
NumPy 1.20.3

Matplotlib 3.5.2

tqdm 4.64.1

scipy 1.9.0

gif 22.5.0

networkx 2.8.4

Once all dependencies are installed, a ```git clone [URL]``` should get the three main files: 1. CastNet.py, 2. CastNet_out_funcs.py, 3. CastNet_parameters.py. Use a standard text editor to modify the parameters on the ```CastNet_parameters.py``` file to suit your needs. Also, notice that the function that runs: main_experiment() is tailor-made for the topology associated with the experimental run 

## Algorithm concept and implementation

CastNet is a sequence evolution simulator that is 'coevolution aware'. The algorithm was initially designed to simulate developmental evolution following the models of Wagner [1994](https://www.pnas.org/doi/abs/10.1073/pnas.91.10.4387) and [1996](https://academic.oup.com/evolut/article/50/3/1008/6870974), as well as elements from Espinosa-Soto's work from [2018](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006172). Briefly, a collection of 'organisms' is produced, each of which is a ragger numpy array of 9 arrays. The elements of the array, per index of the organism array is as follows:
## Organism array elements:
  [0] Generation number
  
  [1] Genome in nucleotides transformed into integers (details below).
  
  [2] Proteome in amino acids transformed into integers (details below).
  
  [3] Gene regulatory matrix in which the regulatory effect of each gene on all other genes is stored.
  
  [4] A vector of thresholds above which each gene must be expressed.
  
  [5] A vector of rates of decay for each gene.
  
  [6] A vector of starting quantities for all genes, to kick-start development.
  
  [7] A matrix showing the progression of gene quantities through a set of discrete developmental steps.
  
  [8] A vector identifying the genes that were activated at any point during development (1-activated, 0-not activated).
  
  [9] The fitness value for that organism.
  
You can access these sections individually for any organism in any population, by the given index number.
