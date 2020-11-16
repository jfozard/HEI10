
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import matplotlib as mpl


### Script to make matplotlib kymograph plot from 


mpl.rcParams.update({ #'figure.figsize': (6.0,4.0),  
    'figure.facecolor': 'none', #(1,1,1,0), # play nicely with white background in the Qt and notebook
    'axes.facecolor': 'none',
    'figure.edgecolor': 'none',       
    'font.size': 20, # 12pt labels get cutoff on 6x4 logplots, so use 10pt.
    'figure.dpi': 72, # 72 dpi matches SVG/qtconsole
    'figure.subplot.bottom' : .15, # 10pt still needs a little more room on the xlabel
    'axes.labelsize':28,
    'savefig.edgecolor': 'none',
    'savefig.facecolor': 'none',
})

LARGE_FS=32


# Kymograph plotting (plotting output from simulation output file)

def main():
    u = []

    # Read simulation output
    with open(sys.argv[1], 'r') as f:
        try:
            h = next(f)
            L, T = map(float, h.split(','))
            h = next(f)
            x = list(map(float, h.split(',')))
            while True:
                l = next(f)
                u.append(list(map(float, l.split(','))))
            
        except StopIteration:
            pass

    u = np.array(u)
    
    n = 256
    m = 100

    T = T/60/60
    
    canvas = np.zeros((n, m))


    n_dt, N = u.shape

    t_data = np.linspace(0, T, n_dt)

    t_plot = (0.5+np.arange(n))/n*T

    # Interpolate simulation data onto a fixed grid
    for i in range(N):
        v = u[:, i]
        f = interp1d(t_data, v)
        vv = f(t_plot)
        idx = int((x[i]/L)*m)
        canvas[:, idx] = np.maximum(canvas[:,idx], vv) # Maximum intensity if two traces overlap


    canvas = np.maximum(canvas, 0.05)
    plt.imshow(canvas, extent=[0,L,T,0], aspect="auto", norm=mpl.colors.LogNorm())
    plt.xlabel('position ($\mu$m)')
    plt.ylabel('time (hr)')

    cbar = plt.colorbar()

    cbar.solids.set_edgecolor("face")

    cbar.set_label("HEI10 (a.u.)")
    
    plt.savefig(sys.argv[2])
#    plt.show()


main()
        
    
    
