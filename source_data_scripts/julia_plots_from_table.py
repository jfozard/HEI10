import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
import random

random.seed(1234)


mpl.rcParams.update({ 
    'figure.facecolor': 'none', 
    'figure.edgecolor': 'none',       
    'font.size': 20, 
    'figure.dpi': 72, 
    'figure.subplot.bottom' : .15, 
    'axes.labelsize':28,
    'savefig.edgecolor': 'none',
    'savefig.facecolor': 'none',
    'svg.fonttype' : 'none',    
})

LARGE_FS=32

import pickle

import statsmodels.api as sm

def pc(v, alpha=0.4):
    idx = (v>=alpha*np.max(v))
    return idx


import os


import csv

def lin_fit(x, y, r2=False):
    X = np.array(x)
    Y = np.array(y)
    X2 = sm.add_constant(X)
    est = sm.OLS(y, X2)
    est2 = est.fit()
    print(est2.summary())
    X = np.linspace(np.min(X), np.max(X), 100)
    X2 = sm.add_constant(X)
    yy = est2.predict(X2)
    if not r2:
        return X, yy
    else:
        return X, yy, est2.rsquared


def plot_centre(all_data, tp=-1):
    
    c = []
    v = []
    

    return data


    
def plot_diff_all(all_data, rel=False, tp=-1):
    

    d = []
    
    for L, pos, hei10 in all_data:
        
        if len(pos) and len(hei10) and tp<len(hei10) and len(pos)==len(hei10[tp]):
          if np.max(hei10[tp])>2*np.mean(hei10[tp]):
                idx = pc(hei10[tp])
                idx = np.where(idx)[0]
                if rel:
                    p = pos[idx]/L
                else:
                    p = pos[idx]
                d += list(np.diff(p))
    plt.figure()
    plt.title('N='+str(len(d)))
    plt.xlabel('Spacing ($\mu$m)', fontsize=LARGE_FS)
    plt.ylabel('Frequency', fontsize=LARGE_FS)
    plt.hist(d, density=True, bins=np.linspace(0,60,10))
    plt.autoscale(enable=True, axis='x', tight=True)

    data = data_hist(d, "Spacing")

    return data



def to_val(v):
    try:
        return int(v)
    except ValueError:
        return float(v)

def line_hist(ax, x, data, label=None):
    data = [0] + list(np.repeat(data, 2)) + [0]
    xx = np.repeat(x, 2)
    print(data, len(data))
    print(xx, len(xx))
    ax.plot(xx, data, label=label)
   # plt.ylim(0, np.max(data)*1.05)

import pandas as pd


"""
t_L                0.1
a_nodes            2.1
dt                 360
Ls                4.93
t_exp            False
b_nodes            0.5
n_c               1.25
n_ts               100
C0                 6.8
u0                 1.2
D                  1.1
m                 2000
density            0.5
K                    1
C0_noise           2.2
t_C0_ratio2          2
n                 1000
Lm                37.5
Number COs           2
L              38.9656
position_0     11.2105
intensity_0    131.154
position_1     38.4461
intensity_1    63.2641
position_2         NaN
intensity_2        NaN
position_3         NaN
intensity_3        NaN
position_4         NaN
intensity_4        NaN
"""

def make_plots(input_csv_path, output_prefix, centre_plot=True, intensity_bins=None, max_n=4, max_k=4, do_density_plot=True):
    # Load preprocessed pkl file
    new_data = pd.read_csv(input_csv_path+'_density_0.5.csv')


    all_hist_n = np.zeros((19,))
    all_hist_nn = np.zeros((4, 10))
    all_hist = np.zeros(10,)
    
    all_intensities = []
    
    nf_mean = np.zeros(3)
    
    mid = np.linspace(0,1,11)
    mid = 0.5*(mid[:-1]+mid[1:])

    nf_tot = []

    all_data_1 = []

    all_lengths = new_data['L']

    all_nco = []

    all_spacing = []

    all_foci = []
    L_array = []

    ns = 0
    for i in new_data.index:
        s = new_data.iloc[i]
        print(s)

        nco = int(s['Number COs'])
        L = s['L']
        all_nco.append(nco)

        foci = []
        foci_raw_pos = []
        foci_raw_intensities = []
        for i in range(nco):
            foci_raw_pos.append(s[f'position_{i}'])
            foci.append(s[f'position_{i}']/L)
            foci_raw_intensities.append(s[f'intensity_{i}'])

        if len(foci_raw_pos)>1:
            all_spacing += list(np.diff(foci_raw_pos))

        all_foci.append(foci)
        foci_raw_intensities = np.array(foci_raw_intensities)
        all_intensities.append(foci_raw_intensities)
        ns += 1
        print(nco)
    
    
    foci = [ p for l in all_foci for p in l ]
        
    n_foci = [ len(l) for l in all_foci ]
    
    all_hist_n, _ = np.histogram(n_foci, range(20))

    
    for n in range(1,max_k):
        f = []
        for j in range(len(all_foci)):
            if n_foci[j]==n:
                f+=list(all_foci[j])
        all_hist_nn[n-1], _ = np.histogram(f, bins=np.linspace(0,1,11))
    f = []
    for j in range(len(all_foci)):
        if n_foci[j]>=max_k:
            f+=list(all_foci[j])
    all_hist_nn[max_k-1], _ = np.histogram(f, bins=np.linspace(0,1,11))


    all_hist, _ = np.histogram(foci, bins = np.linspace(0,1,11))
        

    
    plt.figure()

    length_by_nco = [ np.array(all_lengths)[np.array(all_nco)==n] for n in range(1, np.max(all_nco)+1) ] 

    N_data = [str(len(p)) for p in length_by_nco]
    
    print(length_by_nco, N_data)

    NN = 3

    length_by_nco = [ np.array(all_lengths)[np.array(all_nco)==n] for n in range(1, NN+1) ] 

    
    print(output_prefix, len(length_by_nco), len(np.arange(1,np.max(all_nco)+1) ))

        
    
    plt.violinplot(length_by_nco, positions=np.arange(1,4), showmeans=True, vert=False)

    for i in range(3):
        plt.annotate('n='+str(len(length_by_nco[i])), (20, i+1+0.2), fontsize=18)
    
#    plt.scatter(all_lengths, all_nco)

#    X, yy, r = lin_fit(all_lengths, all_nco, r2=True)
    plt.xlabel('SC length ($\mu$m)')
    plt.ylabel('CO number')
#    plt.plot(X, yy, 'r-')
    plt.xticks([20, 40,60,80])
    plt.yticks([1,2,3])
    plt.xlim(15,85)
    
#    plt.title(r'$r^2='+f'{r:.2f}'+r'$')

    
    plt.savefig(output_prefix+"nco_vs_length.svg")

    plt.title(" ".join(N_data))
    plt.savefig(output_prefix+"nco_vs_length.png")

    ## Output data for violin plot here
    



    plt.figure()

    plt.title('N='+str(len(all_spacing)))
    plt.xlabel('Spacing ($\mu$m)', fontsize=LARGE_FS)
    plt.ylabel('Frequency', fontsize=LARGE_FS)
    plt.hist(all_spacing, density=True, bins=np.linspace(0,60,10))
    plt.autoscale(enable=True, axis='x', tight=True)


    plt.savefig(output_prefix+'diff.svg')
    
    
    with open('../input_data/H2A.csv') as f:
        reader = csv.reader(f)
        H2A = list(reader)

    with open('../input_data/MHL1_diakinesis.csv') as f:
        reader = csv.reader(f)
        MHL = list(reader)


    H2A = np.array([float(p[1]) for p in H2A[:3]])
    MHL = np.array([float(p[1]) for p in MHL[:3]])

    print(H2A, MHL)
    
    dsb_data = H2A/H2A[0]
    CO_data = MHL

    
    plt.figure()
    
    plt.bar(np.arange(max_n+2), all_hist_n[:max_n+2])
    plt.xlabel('Number of COs')
    plt.ylabel('Number of bivalents', size=24)
    plt.xticks(np.arange(max_n+2))
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.savefig(output_prefix+'julia_all_hist_n.svg')


    # USED?
    """
    plt.figure()#figsize=(10,10))
    
    plt.bar(np.arange(max_n+2), all_hist_n[:max_n+2])
    plt.title('Mean CO number {:.02f}'.format(mean_co))
    plt.xlabel('Number of COs')
    plt.ylabel('Number of bivalents', size=24)
    plt.xticks(np.arange(max_n+2))
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.savefig(output_prefix+'julia_all_hist_n_mean.svg')
    """
    
    plt.figure()

    mid = np.linspace(0,1,11)
    mid = 0.5*(mid[:-1]+mid[1:])
    plt.hist(mid, bins=np.linspace(0,1,11), weights=all_hist, histtype='stepfilled')#, width=1.0/29)
    plt.xlabel('Relative position along bivalent', size=24)
    plt.ylabel('Number of COs')

    plt.xlim(0,1)

    
    plt.savefig(output_prefix+'julia_all_hist.svg')

    print(all_hist_n)
    print('mean peaks all', np.sum(all_hist_n*np.arange(19))/np.sum(all_hist_n), np.sum(all_hist_n*np.arange(19))/np.sum(all_hist_n)*5)

    print('MHL', MHL)
    

    for i in range(4):
        plt.figure()
        plt.hist(mid, bins=np.linspace(0,1,11), weights=all_hist_nn[i], histtype='stepfilled')

        
        plt.xlim(0,1)
        plt.xlabel('Relative position along bivalent', size=24)
        plt.ylabel('Number of COs')

        
        plt.autoscale(enable=True, axis='x', tight=True)
        plt.savefig(output_prefix+'julia_hist_'+str(i+1)+'.svg')

        plt.close()

    for i in range(2,4):
        plt.figure()
        m = [ [] for j in range(i) ]
        for c in all_foci:
            if len(c) == i:
                for j in range(i):
                    m[j].append(c[j])
        plt.hist(m, bins=np.linspace(0,1,11), stacked=True)

        
        plt.xlim(0,1)
        plt.autoscale(enable=True, axis='x', tight=True)
        plt.xlabel('Relative position along bivalent', size=24)
        plt.ylabel('Number of COs')
        plt.savefig(output_prefix+f'julia_rank_{i}.svg')
        plt.close()

    data_by_n = [ [] for i in range(max_n)]

    N_cells = len(all_intensities)//5
    data_all = []
    for i in range(N_cells):
        cell_hei10 = [ all_intensities[i+j*N_cells] for j in range(5) ]
        tot_hei10 = np.sum(np.sum(x) for x in cell_hei10)
        norm_hei10 = [ np.array(x)/tot_hei10 for x in cell_hei10]
        print('norm hei10', norm_hei10)
        for n in norm_hei10:
            if len(n)>0:
                idx = min(len(n)-1, max_n-1)
                data_by_n[idx].append(n)
                data_all+= list(n)
        
        
    data_stacked = [ [ p for u in d for p in u ] for d in data_by_n]

    plt.figure()

    if intensity_bins is None:
        bins=40
    else:
        bins=intensity_bins
    
    plt.hist(data_stacked, bins=intensity_bins, stacked=True, label=['{}'.format(i) for i in range(1, max_n)] + [str(max_n)+'+'])
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.legend()
    plt.xlabel('Relative CO intensity')
    plt.ylabel('Number of COs')
    plt.savefig(output_prefix+'julia_peak_int_stacked.svg')

    plt.close()

    plt.figure()
    plt.hist(data_all, bins=intensity_bins)
    plt.autoscale(enable=True, axis='x', tight=True)

    
    plt.xlabel('Relative CO intensity')
    plt.ylabel('Number of COs')
    plt.savefig(output_prefix+'julia_peak_int_cell.svg')
    plt.close()


    if centre_plot:
        plt.figure()

        c = []
        v = []
        for foci_rel_pos, foci_raw_hei10 in zip(all_foci, all_intensities): 
            if len(foci_rel_pos)==2:
                c.append(np.sum(foci_rel_pos)/2)
                v.append(foci_raw_hei10[0]/np.sum(foci_raw_hei10))


        idx = list(range(len(c)))
        if len(idx) < 206:
            print('too few double CO')

        random.seed(1234)
        s = random.sample(idx, 206)
        c = np.array(c)
        v = np.array(v)
    
        plt.plot(c[s], v[s], 'rx')

        

        X, yy, r = lin_fit(c[s], v[s], r2=True)
        plt.plot(X, yy, 'b-')
        plt.title("N="+str(len(c[s]))+" (sample)" + f'    $R^2$ = {r:.2f}')
        plt.xlabel("Midpoint relative position")
        plt.ylabel("Left CO relative intensity", size=22)
        plt.autoscale(enable=True, axis='x', tight=True)
        plt.savefig(output_prefix+'thinned_centre.svg')
    

    print(H2A, MHL)


    if do_density_plot:
        dsb_density = []
        CO_sim = []

        for d in ['0.25', '0.375', '0.5']:
            data = pd.read_csv(input_csv_path+f'_density_{d}.csv') 
            nf = np.sum(data['Number COs'])
            nch = len(data)
            dsb_density.append(float(d))
            CO_sim.append(nf/nch)




        dsb_density = np.array(dsb_density)
        dsb_density = dsb_density/dsb_density[-1]

        dsb_data = H2A/H2A[0]
        CO_data = MHL[:3]/5

        CO_data = CO_data/CO_data[0]


        print('nf tot', nf_tot)

        plt.figure()
        CO_sim = np.array(CO_sim)
        CO_sim = CO_sim/CO_sim[-1]


        print(dsb_density, CO_sim)

        plt.plot(dsb_density, CO_sim)
        plt.plot(dsb_data, CO_data, 'rx')
        plt.xlabel('Relative RI density')
        plt.ylabel('Relative number of \nCOs per bivalent')


        plt.savefig(output_prefix+'julia_hs.svg')
    


source_data_path = '../source_data/'
test_output_path = '../source_data_test_julia/'

os.makedirs(test_output_path, exist_ok=True)

def main():
    make_plots(source_data_path+'fig2_simulations',  test_output_path+'/new_end_', intensity_bins=np.linspace(0,0.3,30))
    make_plots(source_data_path+'figS6_simulations',  test_output_path+'/escape_',  intensity_bins=np.linspace(0,0.3,30))
    make_plots(source_data_path+'fig3_simulations',  test_output_path+'/ox_new_end_',  intensity_bins=np.linspace(0,0.12,11), max_n=7, centre_plot=False)
    make_plots(source_data_path+'fig4_simulations',  test_output_path+'/ux_new_end_',  intensity_bins=np.linspace(0,0.3,30), max_n=3)
    make_plots(source_data_path+'figS3_simulations',  test_output_path+'/no_end_',  intensity_bins=np.linspace(0,0.3,30))
    make_plots(source_data_path+'figS4b_simulations',  test_output_path+'/female_',  intensity_bins=np.linspace(0,0.3,30), do_density_plot=False)

main()
                
        
