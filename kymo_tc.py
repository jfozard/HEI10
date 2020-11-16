
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import matplotlib as mpl

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

# Simulation time-course plotting HEI10 amount at each RI over time

def main():
    u = []
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
    
    m = 100

    T = T/60/60
    

    n_dt, N = u.shape

    t_data = np.linspace(0, T, n_dt)


    plt.figure()
    
    for i in range(N):
        v = u[:, i]
        plt.plot(t_data, v)


    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Time (hr)')

    plt.savefig(sys.argv[2])
#    plt.show()


main()
        
    
    
