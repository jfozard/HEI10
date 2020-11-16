


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

import statsmodels.api as sm

import matplotlib as mpl

from matplotlib.colors import to_rgba_array


from statsmodels.nonparametric.smoothers_lowess import lowess


mpl.rcParams.update({ #'figure.figsize': (6.0,4.0),  
    'figure.facecolor': 'none', #(1,1,1,0), # play nicely with white background in the Qt and notebook
#    'axes.facecolor': 'none', #(1,1,1,0), # play nicely with white background in the Qt and notebook
    'figure.edgecolor': 'none',
 'font.size': 20, # 12pt labels get cutoff on 6x4 logplots, so use 10pt.
 'figure.dpi': 72, # 72 dpi matches SVG/qtconsole
  'figure.subplot.bottom' : .15, # 10pt still needs a little more room on the xlabel
  'axes.labelsize':28,
    'savefig.edgecolor': 'none',
    'savefig.facecolor': 'none',
    
    
})

LARGE_FS = 32

import cycler

def paler(col_str):
    rgb = to_rgba_array(col_str)[0]
    rgb[3] = 0.5
    return rgb

u = [ paler(x ['color']) for x in  mpl.rcParams['axes.prop_cycle'] ]
mpl.rcParams['axes.prop_cycle'] = cycler.cycler(color=u)
    
alpha_max = 0.4



def pc(peaks, alpha=0.2, criterion='max'):
    pos, v, l = peaks
    idx = (v>=alpha*np.max(v))
    return np.array(pos[idx]), v[idx], l


def lin_fit(x, y, min_x=None):
    X = np.array(x)
    Y = np.array(y)
    X2 = sm.add_constant(X)
    est = sm.OLS(y, X2)
    est2 = est.fit()
    print(est2.summary())
    if min_x is None:
        X = np.linspace(np.min(X), np.max(X), 100)
    else:
        X = np.linspace(min_x, np.max(X), 100)
    X2 = sm.add_constant(X)
    yy = est2.predict(X2)
    return X, yy


from scipy.interpolate import interp1d

def interp_resample(x, n):
    l = len(x)
    return interp1d(np.arange(l)/float(l-1), x, kind='linear')(np.linspace(0,1,n))


def analyse_data_ordered(A, data):

    all_peak_data = data['all_peak_data']
    o_hei10_traces = data['o_hei10_traces']
    orig_trace_lengths = data['orig_trace_lengths']
    o_dapi_traces = data['o_dapi_traces']
    
    print(len(A), len(all_peak_data), len(o_hei10_traces), len(orig_trace_lengths))
    
    stages = []
    cell_id = []
    all_quality = []
    all_good = []
    chromosome_id = []
    sig_peaks = []
    sig_peak_hei10 = []
    new_peak_hei10 = []
    orig_peak_hei10 = []

    for i in range(0, len(A), 5):
        stage = A.iloc[(i//5)*5].Stage
        if(type(stage)==str):
            stage = stage.lower()
        stages += [stage]*5
        cell_id += [i//5]*5
        all_quality += [(A.iloc[i:i+5]['quality']==1).all()]*5
        all_good += ['y' if (A.iloc[i:i+5]['good trace?']=='y').all() else 'n']*5
        chromosome_id += list(range(5))
        median_all = np.mean(np.concatenate([o_hei10_traces[k] for k in range(i,i+5)]))
        std_all = np.std(np.concatenate([o_hei10_traces[k] for k in range(i,i+5)]))
        for j in range(i, i+5):
            p, h, l = all_peak_data[j]
            if len(p):
                pp = pc((p, h, l), alpha_max)
                sig_peaks.append(pp[0])
                sig_peak_hei10.append(pp[1])

                h=np.array(o_hei10_traces[j])


                orig_peak_hei10.append(h[pp[0]])
                
                h=(h-median_all) 

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
    A['new_peak_hei10']=new_peak_hei10
    A['orig_peak_hei10']=orig_peak_hei10


    # Sort the 5 chromosomes in each cell by length
    
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

    A = A.sort_values(['cell_id', 'ch_sort_idx'])
    print(A.index)
#    A.set_index(['cell_id', 'ch_sort_idx'])

    print(A['ch_sort_idx'])


    norm_lengths = []
    for i in range(0, len(A), 5):
        # if A.Stage[i]=='late' and A.all_good[i]=='y':
        if True:
            lengths = A['SC length'][i:i+5]
            lengths = lengths/np.sum(lengths)
            norm_lengths+= list(lengths)

    A['norm_lengths']=norm_lengths

    # Think not used for any plots in figures
    """
    rel_peak_int_ch = []
    for i in range(0, len(A)):
        peak_int = A['sig_peak_hei10'][i]

        rel_peak_int = np.array(peak_int)/np.sum(peak_int)
        rel_peak_int_ch.append(list(rel_peak_int))

    A['rel_peak_int_ch']=rel_peak_int_ch  
    
    rel_peak_int_cell = []
    for i in range(0, len(A), 5):
        peak_int = A['sig_peak_hei10'][i:i+5]
        total_cell = np.sum([p for u in peak_int for p in u ])
        for j in range(5):
            rel_peak_int = np.array(peak_int[i+j])/total_cell
            rel_peak_int_cell.append(rel_peak_int)

    A['rel_peak_int_cell']=rel_peak_int_cell
    """


    new_rel_peak_int_cell = []
    for i in range(0, len(A), 5):
        new_peak_int = A['new_peak_hei10'][i:i+5]
        total_cell = np.sum([p for u in new_peak_int for p in u ])
        for j in range(5):
            new_rel_peak_int = np.array(new_peak_int[i+j])/total_cell
            new_rel_peak_int_cell.append(new_rel_peak_int)

    A['new_rel_peak_int_cell']=new_rel_peak_int_cell


    trace_order = []
    long_peaks_dapi_resample = []
    
    for i in A.index:
        dapi = o_dapi_traces[i]
        dapi = (dapi-np.mean(dapi))/np.std(dapi)
        m = len(dapi)
        d = np.mean(dapi[:m//2])>np.mean(dapi[m//2:])
        if d:
            trace_order.append(1)
        else:
            trace_order.append(0)
        if d:
            long_peaks_dapi_resample.append(interp_resample(dapi, 1000))
        else:
            long_peaks_dapi_resample.append(interp_resample(dapi[::-1], 1000))


            
    A['trace_order'] = trace_order

    A['long_peaks_dapi_resample'] = long_peaks_dapi_resample
    
    return A


def plot_stacked(peaks, max_k=3):
    data = []
    for k in range(1,max_k):
        data.append([p for u in peaks for p in u if len(u)==k])
    data.append([p for u in peaks for p in u if len(u)>=max_k])
    
    plt.figure()    
    plt.hist(data, bins=9, range=(0,1), stacked=True, label=[str(y) for y in range(1,max_k)]+[str(max_k)+'+'])
    plt.xlabel('Relative position along bivalent', size=24)
    plt.ylabel('Number of COs')
    plt.legend(fontsize=18)
    plt.xlim(0,1)
    plt.legend()


def plot_data(A, image_output_path, data_fn, max_n, stacked_bins):

    A = A.drop(A.columns[[6,7,8,9]], axis=1)
    A = A.rename(columns={'Plant ':'Plant'})

    with open(data_fn, 'rb') as f:
        data = pickle.load(f)
        all_peak_data, orig_trace_lengths = data['all_peak_data'], data['orig_trace_lengths']
        o_hei10_traces = data['o_hei10_traces']  # Don't have dapi data for UX.


    A = analyse_data_ordered(A, data)

    A = A[(A['all_quality']==1.0) & (A['all_good']=='y') & (A['all_stage']=='late')]

    
    print(len(A), len(A)/5)

    peak_counts = [[] for i in range(5)]
    long_peaks_ordered = [[] for i in range(5) ]
    for i in A.index:
        p = np.array(A['sig_peaks'][i])/float(A['trace_length'][i])
        n = A['ch_sort_idx'][i]
        d = A['trace_order'][i]
        if d:        
            #print(p)
            long_peaks_ordered[n].append(p)
        else:
            long_peaks_ordered[n].append(1-p)


    short_data = long_peaks_ordered[0] + long_peaks_ordered[1]
    short_data_n = [len(p) for p in short_data]

    print('Number of short SC', len(short_data_n), np.sum(short_data_n))

    long_data = long_peaks_ordered[2] + long_peaks_ordered[3] + long_peaks_ordered[4]
    long_data_n = [len(p) for p in long_data]

    print(long_data_n)
    
    plt.figure()
    plot_stacked(short_data)
    plt.savefig(image_output_path+'/ordered_short.svg')
    plt.close()

    plt.figure()
    hist_n(short_data_n)
    plt.savefig(image_output_path+'/ordered_short_n.svg')
    plt.close()

    short_peak_pos = []

    def sc(x):
        return (x-np.min(x))/(np.max(x) - np.min(x) + 1e-8)

    
    n_SC = 0
    for idx in A.index:
        if A['ch_sort_idx'][idx]<2:
             n_SC += 1
             d = A['trace_order'][idx]
             L = A['trace_length'][idx]

             u = sc(data['o_dapi_traces'][idx])
             u = lowess(u, np.arange(len(u)), return_sorted=False, frac=0.1)

             m = np.argmax(u)
             if d:
                 short_peak_pos.append(m/L)
             else:
                 short_peak_pos.append(1-m/L)
                 

    print(' Short bivalents for centromere pos ', n_SC)
                 
    plt.figure()
    plt.hist(short_peak_pos, bins=20, range=(0,1))
    plt.xlim(0,1)
    plt.xlabel('Relative position along bivalent', size=24)
    plt.ylabel('Number of bivalents')
    plt.savefig(image_output_path+'short_peak_pos.svg')
    plt.close()
    print('short peak len', len(short_peak_pos))

    
def hist_n(d, r=None):
    
    h, c = np.unique(d, return_counts=True)
    if 0 not in h:
        h = [0] + list(h)
        c = [0] + list(c)
    print("mean peaks", 5*np.mean(d))
        
    plt.bar(h,c)
    plt.xlabel('Number of COs')
    plt.ylabel('Number of bivalents')
    plt.xticks(range(max(d)+1))
    if r:
        plt.xlim(r[0]-0.5,r[1]-0.5)
        plt.xticks(range(*r))

    return h, c

    


 

base_path = './' 

data_output_path = 'data_output/'

image_output_path='data_output/'


for image_output_base, csv_fn, data_fn, max_n, stacked_bins in  [
        ( image_output_path, '200406.csv', data_output_path+'test.pkl', 4, np.linspace(0, 0.3, 30)),
        ]:
    A = pd.read_csv(csv_fn)

    print(A.keys())
    plot_data(A, image_output_base, data_fn, max_n, stacked_bins)
