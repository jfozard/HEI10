
import os
import sys
import subprocess
from copy import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


## Simulation parameters
default_params = {
    'output_base': '../output/', # Output directory (adjust as appropriate)
    'K':1, 
    'D':1.1, # HEI10 diffusion constant D (\mu m^2 s^{-1})
    'a_nodes': 2.1, # RI HEI10 absorption rate \alpha (\mu m s^{-1})
    'b_nodes': 0.5, # RI HEI10 escape rate \beta (s^{-1})
    'C0_noise_base': 2.2, # Noise in HEI10 initial RI loading \sigma (2.2 a.u.)
    't_L': 0.1,  # Telomere length as fraction of bivalent x_e (no units)
    'dt': 360.0, # Length of each timestep (s) - note that integration timestep is shorter, controlled by the estimated error in the integration
    'n_c': 1.25, # Escape rate Hill coefficient \gamma (no units)
    'n_ts': 100, # Number of integration timesteps (no units)
    'C0_base': 6.8, # HEI10 initial RI loading C_0 (a.u.)
    'u0_base': 1.2, # HEI10 initial loading on bivalent (a.u.)
    'm': 2000, # Number of 
    'K': 1.0,  # Hill function threshold K_C (a.u.)
    'n': 1000, # Number of simulations per chromosome per condition
    'density': 0.5, # Ri density \rho (\mu m^{-1})
    't_C0_ratio2': 2.0, # Increased RI loading at SC telomeres f_e (a.u.) - used for simulations without centromeres
    't_C0_ratio_L': 2.0, # Increased RI loading at LH telomere end (a.u.) - used for simulations with centromeres
    't_C0_ratio_R': 2.0, # Increased RI loading at RH telomere end (a.u.) - used for simulations with centromeres
#    'Lm': [ 37.3, 40.1, 45.3, 53.3, 59.3 ], # Mean arabidopsis SC lengths
#    'Ls': [ 4.71, 4.81, 5.52, 6.86, 7.06 ], # Standard deviation of arabidopsis SC lengths
    'Lm': [ 37.5, 40.4, 45.7, 53.8, 59.8 ], # Mean arabidopsis SC lengths
    'Ls': [ 4.93, 5.06, 5.86, 7.18, 7.13 ], # Standard deviation of arabidopsis SC lengths
    'ratio': 1.0,  # Multiplicative factor modifying HEI10 loading (on both SC and at RI) for UX / OX simulations
    'length_ratio':1.0, # Multiplicative factor used to rescale SC lengths to compare female/male meiosis
}


## Simulations with alternative escape rates
escape_params = copy(default_params)

escape_params['C0_base'] = 0.93
escape_params['C0_noise_base'] = 0.30
escape_params['u0_base'] = 0.16
escape_params['t_C0_ratio2'] = 1.5
escape_params['t_C0_ratio_L'] = 1.5
escape_params['t_C0_ratio_R'] = 1.5
escape_params['n_alpha'] = 1.0

### Changes to matplotlib parameters for plotting figure panels
mpl.rcParams.update({ 
    'figure.facecolor': 'none', 
    'figure.edgecolor': 'none',
    'font.size': 20, 
    'figure.dpi': 72, 
    'figure.subplot.bottom' : .15, 
    'axes.labelsize':28,
    'savefig.edgecolor': 'none',
    'savefig.facecolor': 'none',
})





my_env = os.environ.copy()
my_env['JULIA_NUM_THREADS']='6' # Adjust this as appropriate for your CPU
julia_path='/home/foz/JuliaPro-1.4.0-1/Julia/bin/julia' # Path to Julia executable 
julia_options = []
# Possible to precompile Julia packages to speed-up startup
#julia_options = ["--sysimage", "nsim.so" ] # Cached 


import numpy as np

def plot_ri(output_path, ri_params):
    
    x = np.linspace(0,1,100)
    y = interp1d([ri_params['c_start'],ri_params['c_mid'],ri_params['c_end'],1.0],[1,ri_params['c_dsb_ratio'],1,1])(x)
    
    y0 = np.array([y[0], ri_params['c_dsb_ratio'], 1, 1])
    x0 = [0, ri_params['c_mid'], ri_params['c_end'], 1]
    integral = np.sum(0.5*(y0[1:]+y0[:-1])*np.diff(x0))

    y = (y/integral)*ri_params['density']
    
    print(100*np.sum(y))

    plt.figure()
    plt.plot(x,y)
    plt.ylabel('RI density ($\mu m^{-1}$)')
    plt.xlabel('Relative position along bivalent', size=24)
    plt.ylim(0,1.0)
    plt.xlim(0,1)
    plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    plt.yticks([0,0.5,1])
    plt.savefig(output_path+'/ri_plot.svg')






def simulate(base_params, output_dir="", extra_flags=[], extra_params={}, julia_script="new_sim.jl"):
    # Standard simulations for wild-type (WT), HEI10 over-expresser (OX) and under-expresser (UX)

    # "length_ratio"  is a multiplicative factor modifying SC lengths to permit simulation of female meiosis.
    # "ratio" is a multiplicative factor modifying HEI10 loading (on SC and at RI), used to simulate UX and OX cases
    
    # Simulate at different RI densities to explore crossover homeostasis
    for f in [0.5, 0.75, 1.0]:
        for i in range(5): # Loop over the 5 chromosomes in Arabidopsis
            params = copy(base_params)
            params.update(extra_params)

            # Modify chromosome length and standard deviation according to length_ratio
            params["L"] = params["Lm"][i]*params["length_ratio"] 
            params["Ls"] = params["Ls"][i]*params["length_ratio"]

            # Modify initial HEI10 loading according to ratio 
            
            params["C0_noise"] = params["C0_noise_base"]*params["ratio"] 
            params["C0"]=params["C0_base"]*params["ratio"]
            params["u0"]=params["u0_base"]*params["ratio"]

            # Modify RI density (for crossover homeostasis)
            params["density"]=base_params["density"]*f

            # Construct simulation command line parameters
            simu_params = [julia_path] + julia_options + [julia_script] + sum([ [f'--{x}', f'{params[x]}'] for x in [ 'n', 'm', 'u0', 'C0', 'C0_noise', 'K', 'L', 'Ls', 'D', 'a_nodes', 'b_nodes', 'dt', 'n_c', 'n_ts', 'density', 't_L', 't_C0_ratio2'] ], [])
            simu_params += extra_flags
            filename_base = output_dir + f"/at_{f}_"
            simu_params += ['--filename_base', filename_base]

            # Check output directory exists
            os.makedirs(output_dir, exist_ok=True)

            print(simu_params)

            # Run simulation
            subprocess.run(simu_params, env=my_env)


def simulate_centromere(base_params, output_dir="", extra_flags=[], extra_params={}, julia_script="new_sim_centromeres.jl"):
    # Simulations with reduced RI density within centromeric region



    # "length_ratio"  is a multiplicative factor modifying SC lengths to permit simulation of female meiosis.
    # "ratio" is a multiplicative factor modifying HEI10 loading (on SC and at RI), used to simulate UX and OX cases
    
    # Simulate at different RI densities to explore crossover homeostasis


    for i in range(2):
        params = copy(base_params)
        params.update(extra_params)

        params["L"] = params["Lm"][i]*params["length_ratio"]
        params["Ls"] = params["Ls"][i]*params["length_ratio"]
        params["C0_noise"] = params["C0_noise_base"]*params["ratio"]
        params["C0"]=params["C0_base"]*params["ratio"]
        params["u0"]=params["u0_base"]*params["ratio"]

            # filename_base

        simu_params = [julia_path] + julia_options + [julia_script] + sum([ [f'--{x}', f'{params[x]}'] for x in [ 'n', 'm', 'u0', 'C0', 'C0_noise', 'K', 'L', 'Ls', 'D', 'a_nodes', 'b_nodes', 'dt', 'n_c', 'n_ts', 'density', 't_L', 't_C0_ratio_L', 't_C0_ratio_R'] ], [])
        simu_params += extra_flags
        filename_base = output_dir + f"/at_"
        simu_params += ['--filename_base', filename_base]

        os.makedirs(output_dir, exist_ok=True)
        
        print(simu_params)
            
        subprocess.run(simu_params, env=my_env)


def simulate_kymo(base_params, fn_base, output_dir="", extra_flags=[], extra_params={}, julia_script="kymo_sol.jl", python_script="kymo.py", python_tc_script="kymo_tc.py"):

    i = 2
    for idx, seed in [(1,3), (2,2), (3,16)]:
        params = copy(base_params)
        params.update(extra_params)

        params["L"] = params["Lm"][i]*params["length_ratio"]
        params["Ls"] = params["Ls"][i]*params["length_ratio"]
        params["C0_noise"] = params["C0_noise_base"]*params["ratio"]
        params["C0"]=params["C0_base"]*params["ratio"]
        params["u0"]=params["u0_base"]*params["ratio"]
        params["start"] = seed
        
        simu_params = [julia_path] + julia_options + [julia_script] + sum([ [f'--{x}', f'{params[x]}'] for x in [ 'm', 'u0', 'C0', 'C0_noise', 'K', 'L', 'Ls', 'D', 'a_nodes', 'b_nodes', 'dt', 'n_c', 'n_ts', 'density', 't_L', 't_C0_ratio2', 'start'] ], [])
        
        simu_params += extra_flags
        
        filename = output_dir + '/'+fn_base+ f"{idx}.dat"

        simu_params += ['--filename', filename ]

        os.makedirs(output_dir, exist_ok=True)
        
        print(simu_params)
            
        subprocess.run(simu_params, env=my_env)
    

        

def run_perturb(params):

    # Run sets of simulations with each of the default parameters values
    # increased by 10% in turn in order to assess robustness to parameters
    
    output_base = params['output_base']

    for p in ['K',
              'D',
              'a_nodes',
              'b_nodes',
              'C0_noise_base',
              't_L',
              'dt',
              'n_c',
              'C0_base',
              'u0_base',
              'density',
              't_C0_ratio2',
              ]:
        new_params = copy(params)
        new_params[p] *=1.1
        new_params['output_base'] = output_base + '/perturb_'+p+'/'
        run_simulations(new_params)
    

def run_simulations(params):
    

    output_base = params['output_base']

    ## Check that output directory exists
    os.makedirs(output_base+'/julia_plots', exist_ok=True)

    # Plot the RI distribution for the simulations with centromeres
    plot_ri(output_base+"julia_plots/", { "c_dsb_ratio":0.5,  "c_start":-0.15, "c_end":0.5, "c_mid":0.25, "density":0.5 })


    ## Simulations to make kymograph plots

    # WT
    simulate_kymo(params,
                  "kymo",
                  output_dir=output_base+'/kymo/'
    )

    # Exponentially-decaying telomeric RI loading
    simulate_kymo(params,
                  "exp_kymo",
                  output_dir=output_base+'/kymo/',
                  extra_flags=["--t_exp"]
    )

    # No additional telomere end loading
    simulate_kymo(params,
                  "no_end_kymo",
                  output_dir=output_base+'/kymo/',
                  extra_params = { "t_C0_ratio2": 1.0 },
    )

    ## Simulations to generate histograms of CO density etc
    # WT
    simulate(
        params,
        output_dir=output_base+"survey_julia_new_ends/",
    )

    # OX
    simulate(
        params,        
        extra_params = {'ratio': 4.5 },
        output_dir=output_base+"survey_julia_ox"
    )

    # UX
    simulate(
        params,
        extra_params = {'ratio': 0.6 },
        output_dir=output_base+"survey_julia_ux",
    )

    
    # Reduced RI density within centromere region
    # Centromeric density is a piecewise linear function,
    # taking the value "density" for x/L < c_start and x/L > c_end
    # and the value "density"*"c_dsb_ratio" at x/L = c_mid
    
    simulate_centromere(
        params,
        output_dir=output_base+"survey_centromere8_0.75/",
        extra_flags=["--c_dsb_ratio", "0.5",  "--c_start", "-0.15", "--c_end", "0.5", "--c_mid", "0.25" ],
        julia_script="new_sim_centromeres.jl",
    )

    # Reduced SC lengths for female meiosis 
    simulate(
        params,
        output_dir=output_base+"survey_female/",
        extra_params = {'length_ratio': 0.6 },
    )

    # Exponentially-decaying telomeric RI loading
    simulate(
        params,
        output_dir=output_base+"survey_exp/",
        extra_flags=["--t_exp"],
        extra_params = { 'ratio' :1.0 },
    )

    # No additional telomere end loading
    simulate(
        params,
        output_dir=output_base+"survey_julia_no_ends/",
        extra_params = { "t_C0_ratio2":1.0 },
    )



def run_extra_simulations(escape_params):


    print(escape_params)

    output_base = escape_params['output_base']

    ## Check that output directory exists
    os.makedirs(output_base+'/julia_plots', exist_ok=True)

    # WT
    simulate_kymo(escape_params,
                  "escape_kymo",
                  output_dir=output_base+'/kymo/',
                  julia_script="kymo_escape.jl",

    )


    ## Simulations to generate histograms of CO density etc

    # WT
    simulate(
        escape_params,
        output_dir=output_base+"survey_escape/",
        julia_script="new_sim_escape.jl",

    )


run_simulations(default_params)

# Additional simulations with modified escape rates
run_extra_simulations(escape_params)
