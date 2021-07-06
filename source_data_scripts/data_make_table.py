import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

import os
import os.path
from imageio import imread
from scipy.signal import find_peaks, find_peaks_cwt
from scipy.signal import peak_prominences
import scipy.linalg as la
from functools import reduce
import pickle

alpha_max = 0.4


def pc(peaks, alpha=alpha_max, criterion='max'):
    pos, v, l = peaks
    idx = (v>=alpha*np.max(v))
    return np.array(pos[idx]), v[idx], l



## Analysis of each dataset 
def analyse_data(A, data):

    all_peak_data = data['all_peak_data']
    o_hei10_traces = data['o_hei10_traces']
    orig_trace_lengths = data['orig_trace_lengths']

    
    
    stages = []
    cell_id = []
    all_quality = []
    all_good = []
    chromosome_id = []
    sig_peaks = []
    new_peak_hei10 = []
    orig_peak_hei10 = []
    trace_median_hei10 = []
    trace_sum_hei10 = []
    
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
            hei10=np.array(o_hei10_traces[j])
            if len(p):
                pp = pc((p, h, l), alpha_max)                

                orig_peak_hei10.append(hei10[pp[0]])
                
                sig_peaks.append(pp[0])
                
                                
            else:
                sig_peaks.append([])
                orig_peak_hei10.append([])

            trace_median_hei10.append(np.median(hei10))
            trace_sum_hei10.append(np.sum(hei10))
                
                
                
    print(len(A.index), len(stages))
    
    A['Stage']=stages
    A['original_trace_pixel_length']=list(orig_trace_lengths.values())
    A['trace_pixel_length']=list(v[2] for v in all_peak_data.values())
    A['all_quality']=all_quality
    A['all_good']=all_good
    A['num_foci']=[len(p) for p in sig_peaks]
    A['trace_median_hei10']=trace_median_hei10
    A['trace_sum_hei10']=trace_sum_hei10

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


    cell_total_SC_length = []

    for i in range(0, len(A), 5):
        lengths = A['SC length'][i:i+5]
        tot_length = np.sum(lengths)
        cell_total_SC_length += [tot_length]*5
        
    A['cell_total_SC_length'] = cell_total_SC_length

    
    norm_lengths = []
    for i in range(0, len(A), 5):
        lengths = A['SC length'][i:i+5]
        lengths = lengths/np.sum(lengths)
        norm_lengths+= list(lengths)

    A['norm_lengths']=norm_lengths

    # Now get sig peak positions and intensities


    foci_data = []
    for i in range(len(A)):
        foci_pos = sig_peaks[i]
        foci_hei10 = orig_peak_hei10[i]
        data_line = {}
        for j, (p, h) in enumerate(zip(foci_pos, foci_hei10)):
            data_line[f'pos_{j}'] = p
            data_line[f'HEI10_{j}'] = h
        foci_data.append(data_line)

    print(foci_data)
    fd = pd.DataFrame(foci_data)

    A = pd.concat([A, fd], axis=1)
    
    
    return A







def process_data(csv_fn, pkl_fn, output_file, output_trace_file, specific_traces, trace_output_paths):

    # Read CSV file with data specification
    A = pd.read_csv(csv_fn)


    # Remove some unused columns and rename one field
    A = A.drop(A.columns[[6,7,8,9]], axis=1)
    A = A.drop(columns=[x for x in ('foci count', 'foci per chromosome', 'comments') if x in A])

    A = A.rename(columns={'Plant ':'Plant'})


    data = pickle.load(open(pkl_fn, 'rb'))
    

    print(data.keys())
    
    print(data['orig_trace_lengths'])





    A = analyse_data(A, data)

    for i, v in data['all_peak_data'].items():
        print(v)

    print(A.iloc[0])

    A.to_csv(output_file)


    # Write original hei10 traces to file

    with open(output_trace_file, 'w') as f:
        for i, v in data['o_hei10_traces'].items():
            f.write(', '.join(map(str, v)) + '\n')

    for i, fn in zip(specific_traces, trace_output_paths):
        with open(fn, 'w') as f:
            t = data['o_hei10_traces'][i]
            for v in t: 
                f.write(f'{v}' + '\n')
       
        
    


            
data_output_path = '../output/data_output/'

process_data('../input_data/200406.csv', data_output_path+'test.pkl', '../source_data/fig2_cytology.csv', '../source_data/fig2_cytology_raw_traces.csv', [675, 664, 492], [f'../source_data/fig1a_HEI10_trace_{s}.csv' for s in ['upper', 'mid', 'lower']])
process_data('../input_data/OX.csv', data_output_path+'test_ox.pkl', '../source_data/fig3_cytology.csv', '../source_data/fig3_cytology_raw_traces.csv', [260], ['../source_data/fig3a_HEI10_trace.csv'])
process_data('../input_data/UX.csv', data_output_path+'test_ux.pkl', '../source_data/fig4_cytology.csv', '../source_data/fig4_cytology_raw_traces.csv', [115],  ['../source_data/fig4a_HEI10_trace.csv'])
