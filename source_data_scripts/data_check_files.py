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
    pos, v = peaks
    pos = np.array(pos)
    v = np.array(v)
    idx = (v>=alpha*np.max(v))
    return np.array(pos[idx]), v[idx]









def process_data(csv_fn, pkl_fn, peak_file, trace_file, data_csv_file):

    # Read CSV file with data specification
    A = pd.read_csv(csv_fn)

    B = pd.read_csv(data_csv_file)

    print(list(B))
    
    # Remove some unused columns and rename one field
    A = A.drop(A.columns[[6,7,8,9]], axis=1)
    A = A.drop(columns=['foci count', 'foci per chromosome', 'comments'])

    A = A.rename(columns={'Plant ':'Plant'})


    data = pickle.load(open(pkl_fn, 'rb'))


    all_peak_data = data['all_peak_data']
    o_hei10_traces = data['o_hei10_traces']
    orig_trace_lengths = data['orig_trace_lengths']


    peak_index = []
    peak_height = []
    with open(peak_file, 'r') as pf:
        for i, l in enumerate(pf):
            l = l.strip().split(',')
            length, orig_length = int(l[0]), int(l[1])
            assert(orig_length == orig_trace_lengths[i])
            assert(length == len(o_hei10_traces[i]))

            l = l[2:]
            n = 0

            pi = []
            ph = []
            
            for j in range(len(l)//3):
                idx, v, nv = int(l[j*3]), float(l[j*3+1]), float(l[j*3+2])
                assert(v == o_hei10_traces[i][abs(idx)])
                if idx>0:
                    pi.append(idx)
                    ph.append(nv)
                    assert(idx == all_peak_data[i][0][n])
                    assert(abs(nv-all_peak_data[i][1][n])<1e-10)
                    n += 1
            peak_index.append(pi)
            peak_height.append(ph)
            assert(len(all_peak_data[i][0]) == len(pi))

            
    print('peak files tested')
                    
    with open(trace_file, 'r') as tf:
        for i, l in enumerate(tf):
            t = np.array(list(map(float, l.strip().split(','))))
            u = o_hei10_traces[i]
            assert(np.max(np.abs(t-u))<1e-10)
    print('trace files tested')

        

    print('HEI10 foci compared with peak file')
    
    for i in B.index:
        nf = int(B['num_foci'][i])
        pi, ph = peak_index[i], peak_height[i]
        if len(pi)>0:
            fi, fh = pc((pi, ph))
            assert(nf == len(fi))
            for j in range(len(fi)):
                assert(int(B[f'pos_{j}'][i]) == fi[j])
                assert(abs(B[f'HEI10_{j}'][i]-o_hei10_traces[i][fi[j]])<1e-6)
                
        else:
            assert(nf == 0)
    

data_output_path = 'data_output/'

process_data('200406.csv', data_output_path+'test.pkl', 'peak_data.txt', 'traces2.txt', 'test.csv')
