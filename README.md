# CastNet
In order to use CastNet, you need to have python 3.9+ installed, as well as the following packages:
NumPy 1.20.3

Matplotlib 3.5.2

tqdm 4.64.1

scipy 1.9.0

gif 22.5.0

networkx 2.8.4

Once all dependencies are installed, a ```git clone [URL]``` should get the three main files: 1. CastNet.py, 2. CastNet_out_funcs.py, 3. CastNet_parameters.py. Use a standard text editor to modify the parameters on the ```CastNet_parameters.py``` file to suit your needs. Also, notice that the function that runs: main_experiment() is tailor-made for the topology associated with the experimental run in the paper presenting the algorithm. Those three files and the information contained therein are the only input CastNet needs to run.

Running:

```python CastNet.py [your_run_prefix]```

## Algorithm concept and implementation

CastNet is a sequence evolution simulator that is 'coevolution aware'. The algorithm was initially designed to simulate developmental evolution following the models of Wagner [1994](https://www.pnas.org/doi/abs/10.1073/pnas.91.10.4387) and [1996](https://academic.oup.com/evolut/article/50/3/1008/6870974), as well as elements from Espinosa-Soto's work from [2018](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006172). Briefly, a collection of 'organisms' is produced, each of which is a ragged numpy array of 9 arrays. The elements of the array, per index of the organism array is as follows:
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

### [0] Generation number
This value is initialized at 0 with the founder, and will count onwards from then for every branch.

### [1]Genome
For a user-defined number of genes _n_, and a sequence length of _s_ (which CastNet corrects to the closest number that is a multiple of 3), this is a matrix of _n_ rows (each gene) by _s/3_ columns, each containing a codon. For more efficient computing, codons and, by extension, sequences are not coded as strings, but as integers. The equivalences are as follows: ```A: 1, C: 2, T: 3, G: 4```. In order to see the sequences as a string of characters, use the ```translate_int_seq()``` function in the file ```CastNet_out_funcs.py```. Once a regular run is finished, CastNet automatically produces an alignment for each gene, in FASTA character format, using the prefix given at runtime to identify each gene. The 'ali_saver()' function in the output functions file will produce these for an abitrary number of populations, and will randomly choose one member of the population for a sequence.

### [2] Proteome
This array has the exact same dimensions as the genome, but instead of codons, it has the unicode equivalent of the single-letter code for each amino acid. Stop codons are marked with the underscore, "\_", which is unicode value 95. The translation table used is based on the Standard Genetic code (i.e. vertebrate nuclear). Other codes can be added, although it is unlikely they will affect the broader results.

### [3] Gene regulatory matrix: __R__
This is an _n x n_ matrix which contains all the gene-gene interactions that can occur. Each column contains the value by which a gene's arbitrary quantity will affect another gene's arbitrary quantity. Note that this matrix is non-reversible: values above and under the diagonal are different, and mean different things. For example, the regulatory effect of Gene A (column) on Gene C (row) is __not__ the same as the regulatory effect of Gene C (column) on Gene A (row). This is akin to a gene regulatory network (GRN) as stipulated by Wagner and Espinosa-Soto (see references above). The diagonal contains the values reflecting gene self-regulation.

### [4] Thresholds vector: __*θ*__
This is a vector of size _n_ which contains, for each gene, the minimum amount of activating signal that it has to receive in order for it to 'fire', or be activated. This is implemented as a rectified linear unit (ReLU) function which has gene-independent floors given by this vector, as opposed to rectifying to 0 as it is standard when applying ReLU as an activation function. This simulates gene activations that may require higher or lower levels of signal to be switched on.

### [5] Decay rate vectors: __*λ*__
This is also a vector of size _n_ which contains, the rate of exponential decay through discrete time steps for each gene. At each time step, the residual amount  of a gene (i.e. the amount remaining since the last developmental step) is computed by transferring the value from the previous step, and multiplying it by Euler's constant _e_ raised to the power of the negative of that gene's specific decay rate. This simulates the digestion of gene products through time.

### [6] Starting quantity vector
This is also a vector of size _n_ which is used to initialize development across all organisms, populations, and branches. The default is set to 1 for the first gene, and 0 for all other genes, meaning that the first gene kick-starts the developmental process. There is no good analog to this in biology, but it was needed computationally in order for development to start in a consistent manner across all individuals. It can be modified within the CastNet.py source code to different values, if needed.

### [7] Development: __E__
This array contains the progression of gene quantities throughout development for an organism. As can be inferred from the array element above, the first row has 1 for the first gene, and 0 for all others, the following values are all extracted from the operation described in the paper presenting this software, which takes into account all the regulatory parameters (__R__, __*θ*__, and __*λ*__), and applies them to the residual amount of each gene over the user-defined _m_ number of discrete developmental steps.


