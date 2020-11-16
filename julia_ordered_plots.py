

import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys

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



import pickle

import statsmodels.api as sm


def pc(v, alpha=0.4):
    idx = (v>=alpha*np.max(v))
    return idx


#@numba.njit
def delete_array(x,i):
    mask = np.zeros(x.shape, dtype=np.int32) == 0
    mask[i] = False
    return x[mask]


#@numba.njit
def winner_orig(x0):
    x = x0
    j = 0
    while True:
        M = len(x)
        if M==2:
            if abs(x[0]-0.5)<abs(x[1]-0.5):
                return x[0]
            else:
                return x[1]
        
        d = 0.5*(x[2:]-x[:-2])
        d0 = 0.5*(x[0]+x[1])
        dL = 1.0 - 0.5*(x[-1]+x[-2])
        i = np.argmin(d)
        low = d[i]
        i += 1
        if d0<low:
            i=0
            low=d0
        if dL<low:
            i=M-1
            low=dL
        x = delete_array(x, i)
        #x = np.delete(x, i)

def read_file(fn):

    print('fn', fn)
    all_data = []
    head = None
    try:
        with open(fn,'r') as f:
            try:
                nbuf = []
                head = next(f)
                print('head', head)
                head = head[head.index('(')+1:head.index(')')]
                head = head.replace(' ', '')
                while True:
                    nbuf = []
                    L = next(f)
                    L = float(L)
                    pos = next(f)
                    pos = np.array(list(map(float, pos.split(','))))
                    while True:
                        hei10 = next(f).strip()
                        if not hei10:
                            break
                        hei10 = np.array(list(map(float, hei10.split(','))))
                        nbuf.append(hei10)
                    all_data.append((L, pos, nbuf))
            except StopIteration:
                if nbuf:
                    all_data.append((L, pos, nbuf))
    except IOError:
        pass


    head_dict = dict([u.split('=') for u in head.split(',')])
    del head_dict['start']
    return head_dict, all_data


import os

def summary(all_data, tp=-1, nbins=31):


    all_foci = []
    foci_intensities = []
    L_array = []
    ns = 0
    ndsb = 0
    all_pos = []
    
    for L, pos, hei10 in all_data:
        all_pos += list(pos/L)
        if len(pos) and len(hei10) and tp<len(hei10) and len(pos)==len(hei10[tp]):
#            if np.max(hei10[tp])>2*np.mean(hei10[tp]):
                ndsb += len(hei10[tp])
                idx = pc(hei10[tp])
                all_foci.append(pos[idx]/L)
                foci_intensities.append(hei10[tp][idx])
                ns += 1
                L_array.append(L)
#            else:
#                all_foci.append([])
#                ns+=1
#                foci_intensities.append([])
#                L_array.append(L)

    
    foci = [ p for l in all_foci for p in l]

#    data_by_n = [ [ u  if c==i ] for  i in range(1,max_n)] + [ [ u for u,c in zip(A_late['rel_peak_int_cell'], A_late['num_sig_peaks']) if c>=max_n ] ]

    def flatten(l):
        return [p for u in l for p in u]

    max_n = 3



    print('all data', len(all_data))
    print('n foci', len(foci), 'n ch', len(all_foci), 'ndsb', ndsb)

    
    n_foci = [len(l) for l in all_foci]

    print('n_foci_all', n_foci)
    
    hist_n, _ = np.histogram(n_foci, range(6))

    print(hist_n)
    
    hist_nn = {}
    
    for n in range(1,max_n):
        f = []
        for j in range(len(all_foci)):
            if n_foci[j]==n:
                f+=list(all_foci[j])
        hist_nn[n], _ = np.histogram(f, bins=np.linspace(0,1,nbins))

    f = []
    for j in range(len(all_foci)):
        if n_foci[j]>=max_n:
            f+=list(all_foci[j])
    hist_nn[max_n], _ = np.histogram(f, bins=np.linspace(0,1,nbins))

    
    

    hist, _ = np.histogram(foci, bins = np.linspace(0,1,nbins))

    return hist_n, hist_nn, hist, len(foci), L_array, all_pos
    

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


def plot_double(fn1, fn2, output_path):


    print(fn1)
    print(fn2)
    data1 = read_file(fn1)

    hist_n1, hist_nn1, hist1, nf1, L_array1, all_pos1 = summary(data1[1], tp=-1)


    data2 = read_file(fn2)

    hist_n2, hist_nn2, hist2, nf2, L_array2, all_pos2 = summary(data2[1], tp=-1)




    
    NB=11

#    print(hist, len(hist))


    plt.figure()
    plt.bar(range(5), hist_n1+hist_n2)
    plt.xlabel('Number of COs')
    plt.ylabel('Number of bivalents')

    plt.savefig(output_path+'/julia_ordered_short_n.svg')

    plt.figure()
    plt.hist([(np.arange(30)+0.5)/30]*3, bins=np.linspace(0,1,NB), weights=[hist_nn1[1]+hist_nn2[1], hist_nn1[2]+hist_nn2[2], hist_nn1[3]+hist_nn2[3]], label=['1','2','3+'], stacked=True)
    plt.legend(fontsize=18)
    plt.xlim(0,1)
    plt.ylim(0,550)
    plt.xlabel('Relative position along bivalent', size=24)
    plt.ylabel('Number of COs')
    plt.savefig(output_path+'/julia_ordered_short.svg')

    plt.figure()
    plt.hist(all_pos1+all_pos2, bins=100)
    plt.xlim(0,1)
    plt.xlabel('Relative position along bivalent', size=24)
    plt.ylabel('Number of RIs')
    plt.savefig(output_path+'/julia_ordered_short_dsb.svg')
    



def main(sim_data_path, output_path):
    print(sim_data_path, os.listdir(sim_data_path+'/survey_centromere8_0.75'))
    plot_double(*[sim_data_path+'/survey_centromere8_0.75/'+f for f in os.listdir(sim_data_path+'/survey_centromere8_0.75')], output_path)

main(sys.argv[1], sys.argv[2])
                
        
