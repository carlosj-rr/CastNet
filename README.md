# CastNet

![CastNet logo: a cast net with a phylogeny embedded in the net](https://github.com/carlosj-rr/CastNet/blob/canonical/figs/CastNet_logo_smol.png)

CastNet is platform-independent, as long as you can run python.

In order to use CastNet, you need to have python 3.9+ installed, as well as the following packages:

```
NumPy 1.20.3  
Matplotlib 3.5.2  
tqdm 4.64.1  
scipy 1.9.0  
gif 22.5.0
```
Once all dependencies are installed, running `git clone https://github.com/carlosj-rr/CastNet.git` should get the three main files: 
* CastNet.py
* CastNet_out_funcs.py
* CastNet_parameters.py
 
Use a standard text editor to modify the parameters on the ```CastNet_parameters.py``` file to suit your needs (see [wiki](https://github.com/carlosj-rr/CastNet/wiki) for help and examples! _especially_ [chapter 4](https://github.com/carlosj-rr/CastNet/wiki/4.-Making-your-own-experiment)). Also, notice on the main file ```CastNet.py```, the function main_experiment() is tailor-made for the topology used in the experimental runs in the paper presenting the algorithm. Those three files and the information contained therein are the only input CastNet needs to run. The parameter values on ```CastNet_parameters.py``` are also the ones used for the paper, specifically for the run with stringent selection.

Running:

```python CastNet.py [your_run_prefix]```
