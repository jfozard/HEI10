
## Code to generate the majority of the data analysis plots

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

import statsmodels.api as sm

import matplotlib as mpl

import json

from matplotlib.colors import to_rgba_array

from statsmodels.stats.weightstats import ttest_ind
from scipy.stats import mannwhitneyu, f_oneway, kruskal

import rpy2

from rpy2.robjects.packages import importr

utils = importr('utils')

## Uncomment if necessary

#utils.chooseCRANmirror(ind=1)
#utils.install_packages('dunn.test')
#utils.install_packages('kruskal')

import rpy2.robjects as robjects

dunn = importr('dunn.test')
#kruskal_r = importr('kruskal.test')
stats = importr('stats')

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

from matplotlib import lines


import cycler

def paler(col_str):
    rgb = to_rgba_array(col_str)[0]
    rgb[3] = 0.5
    return rgb

u = [ paler(x ['color']) for x in  mpl.rcParams['axes.prop_cycle'] ]
mpl.rcParams['axes.prop_cycle'] = cycler.cycler(color=u)
    
alpha_max = 0.4


# Criterion used for identifying peak as a CO - normalized (with mean and s.d.)
# hei10 levels being above 0.4 time maximum peak level
def pc(peaks, alpha=alpha_max, criterion='max'):
    pos, v, l = peaks
    idx = (v>=alpha*np.max(v))
    return np.array(pos[idx]), v[idx], l

# Linear regression best-fit line
# Returns best fit line on range (max(min_x, min(x), max(x))
def lin_fit(x, y, min_x=None, of=None, r2=False):
    X = np.array(x)
    Y = np.array(y)
    X2 = sm.add_constant(X)
    est = sm.OLS(y, X2)
    est2 = est.fit()
    if of is None:
        print(est2.summary())
    else:
        print(est2.summary(), file=of)
    if min_x is None:
        X = np.linspace(np.min(X), np.max(X), 100)
    else:
        X = np.linspace(min_x, np.max(X), 100)
    X2 = sm.add_constant(X)
    yy = est2.predict(X2)
    if not r2:
        return X, yy
    else:
        return X, yy, est2.rsquared

## Analysis of each dataset (similar code as for data_analysis_plots)
def analyse_data(A, data):

    all_peak_data = data['all_peak_data']
    o_hei10_traces = data['o_hei10_traces']
    orig_trace_lengths = data['orig_trace_lengths']

    #print(len(A), len(all_peak_data), len(o_hei10_traces), len(orig_trace_lengths))
    
    stages = []
    cell_id = []
    all_quality = []
    all_good = []
    chromosome_id = []
    sig_peaks = []
    new_peak_hei10 = []
    orig_peak_hei10 = []


    ### For each cell in turn, process the 

    for i in range(0, len(A), 5):
        stage = A.iloc[(i//5)*5].Stage
        if(type(stage)==str):
            stage = stage.lower()
        stages += [stage]*5
        cell_id += [i//5]*5
        all_quality += [(A.iloc[i:i+5]['quality']==1).all()]*5
        all_good += ['y' if (A.iloc[i:i+5]['good trace?']=='y').all() else 'n']*5
        chromosome_id += list(range(5))
        median_all = np.median(np.concatenate([o_hei10_traces[k] for k in range(i,i+5)]))
        std_all = np.std(np.concatenate([o_hei10_traces[k] for k in range(i,i+5)]))
        for j in range(i, i+5):
            p, h, l = all_peak_data[j]
            if len(p):
                pp = pc((p, h, l), alpha_max)

                
                h=np.array(o_hei10_traces[j])

                sig_peaks.append(pp[0])
                
                orig_peak_hei10.append(h[pp[0]])

                h = h-np.median(h)
                #h=(h-median_all) #/std_all
                #h = (h-np.mean(h))/np.std(h)

                new_peak_hei10.append(h[pp[0]])

                
            else:
                sig_peaks.append([])
                new_peak_hei10.append([])
                orig_peak_hei10.append([])

    print(len(A.index), len(stages))
    
    A['all_stage']=stages
    A['cell_id']=cell_id
    A['chromosome_id']=chromosome_id
    A['orig_trace_length']=list(orig_trace_lengths.values())
    A['trace_length']=list(v[2] for v in all_peak_data.values())
    A['peaks']=list(v[0] for v in all_peak_data.values())
    A['peak_hei10']=list(v[1] for v in all_peak_data.values())
    A['all_quality']=all_quality
    A['all_good']=all_good
    A['sig_peaks']= sig_peaks
    A['num_sig_peaks']=[len(p) for p in sig_peaks]
    A['new_peak_hei10']=new_peak_hei10     # HEI10 intensities at peaks, with median cell HEI10 levels subtracted
    A['orig_peak_hei10']=orig_peak_hei10   # Raw HEI10 intensities at peaks


    # Sort the 5 chromosomes in each cell by length

    # Fo each cell, sort the chromosomes in order
    
    def analyse_SC_lengths(df):
        ch_sort_idx = []
        for j, i in enumerate(df.index[::5]):
            p = df.iloc[5*j:5*j+5]
            s = sorted(p['SC length'])
            sidx = np.argsort(p['SC length'].to_numpy())
            rank = np.zeros((5,), dtype=np.int32)
            for k in range(5):
                rank[sidx[k]] = k
            ch_sort_idx += list(rank)
        return ch_sort_idx


    print(len(A.index), len(analyse_SC_lengths(A)))
    A['ch_sort_idx'] = analyse_SC_lengths(A)

    A.set_index(['cell_id', 'chromosome_id'])


    # Lengths normalized by the total length of all SC in the cell

    norm_lengths = []
    for i in range(0, len(A), 5):
        lengths = A['SC length'][i:i+5]
        lengths = lengths/np.sum(lengths)
        norm_lengths+= list(lengths)

    A['norm_lengths']=norm_lengths


    # Significant peak intensities (with the mean cell HEi10 trace value subtracted off)
    # normalized by the sum of all peak intensities in the cell

    new_rel_peak_int_cell = []
    for i in range(0, len(A), 5):
        new_peak_int = A['new_peak_hei10'][i:i+5]
        total_cell = np.sum([p for u in new_peak_int for p in u ])
        for j in range(5):
            new_rel_peak_int = np.array(new_peak_int[i+j])/total_cell
            new_rel_peak_int_cell.append(new_rel_peak_int)

    A['new_rel_peak_int_cell']=new_rel_peak_int_cell
    
    return A





def plot_corr(A, image_output_path, data_fn):


    # Edit CSV columns + naming
    A = A.drop(A.columns[[6,7,8,9]], axis=1)
    A = A.rename(columns={'Plant ':'Plant'})


    # Read file generated by data_preprocess.py
    with open(data_fn, 'rb') as f:
        data = pickle.load(f)
        all_peak_data, orig_trace_lengths = data['all_peak_data'], data['orig_trace_lengths']
        o_hei10_traces = data['o_hei10_traces']  # Don't have dapi data for UX.


    # Initial processing of data for each cell
    A = analyse_data(A, data)


    print(A.keys())
    # Restrict data to SCs from late cells with good quality traces
    A_late = A.loc[(A.all_stage=='late') & (A.all_quality==1)]

    print('Length of A_late', len(A_late))


    A_late2 = A.loc[(A.all_stage=='late') & (A.all_quality==1) & (A.all_good=='y')]

    print('Length of A_late2', len(A_late2))



    # Find numbers of plants and SCs used
    cs = []
    for i in A_late.index:
        cs.append((A_late.Genotype[i], A_late.Date[i], A_late.Plant[i]))

    cs = list(set(cs))
    cs.sort()

    print('Plant names', cs)
    
    print(' Data from {} Plants'.format(len(cs)))


    
    print('Data from {} SC'.format(len(A_late)))


    
    # Absolute SC length per crossover number
    
    plt.figure()
    data = []
    NN = 3



    print('NN', NN)
    for i in range(1,NN+1):
        data.append(A_late['SC length'][A_late['num_sig_peaks']==i])

    total_SC_length = [ d for p in data for d in p ]
    total_CO_number = [ n+1 for n,p in enumerate(data) for d in p ]
    
        
    plt.violinplot(data, positions=range(1,NN+1), showmeans=True, vert=False)



    for i in range(3):
        plt.annotate('n='+str(len(data[i])), (17, i+1+0.22), fontsize=18)

    plt.ylabel('CO number')
    plt.xlabel('SC length ($\mu$m)')
    plt.xticks([20,40,60,80])
    plt.yticks([1,2,3])
    plt.xlim(15,85)
    plt.savefig(output_path+'/violin_number_length.svg')

    with open(output_path+'/violin_number_length_data.csv', 'w') as f:
        f.write('num_sig_peaks, SC_length\n')
        for i in range(1, NN+1):
            d = data[i-1]
            for v in np.array(d):
                f.write(str(i) + ', ' + str(v) + '\n')
                
            
    

import sys        

input_path = '../input_data/'
data_input_path = sys.argv[1]
output_path = sys.argv[2]


for image_output_base, csv_fn, data_fn in  [
        ( output_path, input_path+'/200406.csv', data_input_path+'/test.pkl')
        ]:
    A = pd.read_csv(csv_fn)
    plot_corr(A, image_output_base, data_fn)

