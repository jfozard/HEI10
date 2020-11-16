## Installation instructions for Ubuntu 20.04.1LTS ( Linux 5.4.0-53-generic )

(Installation has only been tested on this system; Ryzen 3600 with 32GB )

1. Install Julia, ideally using the free version of JuliaPro (v1.4.0-1) https://juliacomputing.com/products/juliapro/ 
   Alternatively download julia `wget https://julialang-s3.julialang.org/bin/linux/x64/1.4/julia-1.4.2-linux-x86_64.tar.gz && tar zxvf julia-1.4.2-linux-x86_64.tar.gz`
2. Install julia packages - either using the Pkg REPL - see https://docs.julialang.org/en/v1/stdlib/Pkg/ 
   or through running `julia add_packages.jl`
   Required Julia packages are : ArgParse DifferentialEquations DiffEqCallbacks Random SparseArrays StatsBase Distributions
   Precise details of complete julia environment used listed in `julia_installed_packages.txt`

3. Install conda (via https://docs.conda.io/en/latest/miniconda.html ; or `https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh` )
4. Create a new conda environment `conda create -n crossover python=3.7` . This needs to be python3.7 to reproduce image analysis results correctly.
5. Activate this conda environment `conda activate crossover`
5. Install packages needed `conda install matplotlib scipy statsmodels pandas lxml imageio click`
6. Install pip-installable packages `pip install svgutils dtoolbioimage ray`
7. Install inkscape (needed for figure generation) `sudo apt-get install inkscape`

Alternatively to 4-7, the conda environment can be recreated using
`conda env create -f environment.yml`
(this specifies the precise versions of packages used by the authors)

## To recreate results of manuscript

1. Run all simulations (takes several hours) `python run_simulations.py` . Modify `my_env['JULIA_NUM_THREADS']='6'` to the appropriate value
   for your machine (optional).
2. Dowload imaging data, placing these large files in a suitable location. For review purposes, these are hosted on Dropbox:
   https://www.dropbox.com/s/bugx89pwj8y598u/morgan_ox.tar.bz2?dl=0   (5GB)
   https://www.dropbox.com/s/y58jhqhxcq5as3v/morgan_ux.tar.bz2?dl=0   (3.5GB)
   https://www.dropbox.com/s/v2o7i142cy7abrg/morgan_wt.tar.bz2?dl=0   (16GB)

3. Decompress (bunzip2) and extract (tar xvf) these imaging files
4. Edit the `data_path = 'xxxx'` line in data_preprocess.py to point to this downloaded data
5. `python data_preprocess.py` (takes an hour or so).
6. Run `bash all_plots.sh` . Figures will be placed in `output/figures` folder.


## To run short test of simulation code (c 1 min)

1. Ensure that Julia and the needed packages are installed (1. and 2. of Installation instructions)
2. `julia new_sim.jl --n 2` (The julia executable may need to be added to your PATH).
3. Output is written to `test0.1,2.1,360.0,30.0,0.0,false,0.5,1.25,1,100,6.8,1.2,1.1,2000,0.5,1.0,2.2,2.0,1.dat` in the current directory
 Expected contents:
``` #Namespace(t_L=0.1,a_nodes=2.1,dt=360.0,L=30.0,Ls=0.0,filename_base=test,t_exp=false,b_nodes=0.5,n_c=1.25,start=1,n_ts=100,C0=6.8,u0=1.2,D=1.1,m=2000,density=0.5,K=1.0,C0_noise=2.2,t_C0_ratio2=2.0,n=1)
30.0
0.2372785017168222,6.284171195942124,6.329046064756079,7.549865490959156,8.435706968571894,9.381209050082024,10.395510425758815,12.741535514853943,13.113239238288754,14.658384902385036,16.67253261973717,23.19669145372131,28.55749019507202,29.59999100696399,29.
99713976695841
0.37283953307113954,0.37283953274230386,0.37283953273767434,0.37283953256574376,0.3728395324098933,0.37283953221189337,0.3728395319667844,142.744517700139,0.3728395320097407,0.3728395347390615,0.3728395381395296,0.3728395480086411,0.372839553939854,0.372839554
6504659,0.3728395547785809
```
This performs a single simulation, with the default parameters and the starting seed 1. 
Output consists of:
l1 - the parameters used
l2 - the length of the SC simulated
l3 - the positions of the initial RI
l4 - HEI10 levels at the RI at the end of the simulation

4. Parameter values and simulation results can be modified through command-line options: `julia new_sim.jl --help` for a complete list
