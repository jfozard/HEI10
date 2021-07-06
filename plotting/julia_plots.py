import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
import random

import pickle

import statsmodels.api as sm
import pandas as pd


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

# Peak criterion
def pc(v, alpha=0.4):
    idx = (v>=alpha*np.max(v))
    return idx


# Code to read output data file
def read_file(fn):

    all_data = []
    head = None
    try:
        with open(fn,'r') as f:
            try:
                nbuf = []
                head = next(f)
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


# Make histograms, length and foci intensity arrays for the data from a single file

def summary(data, L, K, tp=-1, max_k=4):

    head, all_data = data

    all_foci = []
    foci_intensities = []
    L_array = []

    orig_foci_intensities = []
    ns = 0
    for L, pos, hei10 in all_data:
        if len(pos) and len(hei10) and tp<len(hei10) and len(pos)==len(hei10[tp]):
            idx = pc(hei10[tp])
            if len(idx) == sum(idx):
                idx = []
            all_foci.append(pos[idx]/L) # Positions normalized to length of SC
            foci_intensities.append(hei10[tp][idx]/np.sum(hei10[tp][idx]))   # Intensities normalized to sum on the same SC
            orig_foci_intensities.append(hei10[tp][idx])
            ns += 1
            L_array.append(L)
    

    
    foci = [ p for l in all_foci for p in l]  # list of all of the (relative) foci positions 


    n_foci = [len(l) for l in all_foci] # number of foci on each SC
   
    hist_n, _ = np.histogram(n_foci, range(20)) # Histogram of relative foci numbers

    hist_nn = {} # hist_nn[k] = Histogram of relative foci positions for SCs with k foci
    
    for n in range(1,max_k):
        f = []
        for j in range(len(all_foci)):
            if n_foci[j]==n:
                f+=list(all_foci[j])
        hist_nn[n], _ = np.histogram(f, bins=np.linspace(0,1,11)) 

    f = []
    for j in range(len(all_foci)):
        if n_foci[j]>=max_k:
            f+=list(all_foci[j])
    hist_nn[max_k], _ = np.histogram(f, bins=np.linspace(0,1,11))


    hist, _ = np.histogram(foci, bins = np.linspace(0,1,11)) # hist_nn[k] = Histogram of all relative foci positions

    return hist_n, hist_nn, hist, len(foci), len(all_foci), L_array, foci_intensities, orig_foci_intensities, all_foci


### Linear regression

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


# Plot for fig 2k - for SCs with two COs, rel intensity of left peak vs position of centre
def plot_centre(all_data, tp=-1):
    
    c = []
    v = []
    
    for L, pos, hei10 in all_data:
        
        if len(pos) and len(hei10) and tp<len(hei10) and len(pos)==len(hei10[tp]):
                idx = pc(hei10[tp])
                idx = np.where(idx)[0]
                
                if len(idx)==2:
                    c.append(np.sum(pos[idx])/2/L)
                    v.append(hei10[tp][idx[0]]/np.sum(hei10[tp][idx]))


        
    idx = list(range(len(c)))
    if len(idx) < 206:
        print('too few double CO')
        return
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

    

# Plot spacing between consecutive foci
def plot_diff_all(all_data, rel=False, tp=-1):

    d = []
    
    for L, pos, hei10 in all_data:
        
        if len(pos) and len(hei10) and tp<len(hei10) and len(pos)==len(hei10[tp]):
                idx = pc(hei10[tp])
                idx = np.where(idx)[0]
                if len(idx)<len(hei10[tp]):
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



 # Load all data from a directory   
def load_data(survey_dir, survey_base):
    
    files = os.listdir(survey_dir)

    all_data_ext = {}
    for f in files:
        if len(f)>=len(survey_base) and f[:len(survey_base)] == survey_base:
            params = f[len(survey_base):-4]
            all_data_ext[params] = read_file(survey_dir+'/'+f)
    

    def to_val(v):
        try:
            return int(v)
        except ValueError:
            try:
                return float(v)
            except ValueError:
                return bool(v)
    def get_class(params):
        return tuple(to_val(f) for f in params.replace('_',',').split(','))
    
    for k in all_data_ext:
        print(get_class(k))

    new_data = {}
    for k in all_data_ext:
        p = get_class(k)[0:8] + get_class(k)[9:] # Omit the starting index from file
        if p in new_data:
            new_data[p][1] += all_data_ext[k][1]
        else:
            h, d = all_data_ext[k]
            new_data[p] = [h, d]


    return new_data


def to_val(v):
    try:
        return int(v)
    except ValueError:
        return float(v)


    

## Process all output files for WT/OX/UC etc
def make_plots(data_path, output_prefix, centre_plot=True, intensity_bins=None, max_n=4):
    # Load preprocessed pkl file
    new_data = load_data(*data_path) 
    k = next(iter(new_data))


    # Organize data into groups 
    new_data_keys = [ (idx, to_val(h['L']), to_val(h['density'])) for idx, (h, v) in enumerate(new_data.values()) ]
    idx, L, density = map(list, zip(*new_data_keys))
    
    df = pd.DataFrame(data={'idx':idx, 'L':L, 'density':density})
        
    df2 = df.sort_values(['L', 'density'])
    
    df2['density_idx'] = list(range(3))*5

    grouped = df2.groupby(['L'])
    
    new_data_list = list(new_data.values())

    all_hist_n = np.zeros((19,))
    all_hist_nn = np.zeros((4, 10))
    all_hist = np.zeros(10,)

    all_lengths = []
    
    all_intensities = []
    
    nf_mean = np.zeros(3)
    
    mid = np.linspace(0,1,11)
    mid = 0.5*(mid[:-1]+mid[1:])

    nf_tot = []

    all_nco = []

    
    all_data_1 = []
    all_co = []

    for group_idx, (name, group) in enumerate(grouped):
        s = group.iloc[-1,:]

        h, v = new_data_list[int(s['idx'])]

        all_data_1 += v
        
        L = to_val(h['L'])
        density = to_val(h['density'])
        K = to_val(h['K'])

        hist_n, hist_nn, hist, nf, nch, L_array, intensities, orig_foci_intensities, co = summary((h,v), L, K, tp=-1)
        
        all_intensities.append(orig_foci_intensities)
        
        all_hist_n += np.array(hist_n)
        all_hist += np.array(hist)

        all_hist_nn += np.array(list(hist_nn.values()))

        all_co += co

        all_lengths += list(L_array)
        
        all_nco += [len(z) for z in co]
          
        dsb_density = []
        nf_all = []
        for i in range(3):
            s = group.iloc[i,:]
            j = int(s['idx'])
            dsb_density.append(s['density'])
            _, _, _, nf, nch, _, _, _, _ = summary(new_data_list[j],  L, K, tp=-1)
            nf_all.append(nf/float(nch))
            nf_mean[i] += nf/float(nch)
        nf_tot.append(np.array(nf_all))
        



    plt.figure()

    length_by_nco = [ np.array(all_lengths)[np.array(all_nco)==n] for n in range(1, np.max(all_nco)+1) ] 

    N_data = [str(len(p)) for p in length_by_nco]
    

    NN = 3

    length_by_nco = [ np.array(all_lengths)[np.array(all_nco)==n] for n in range(1, NN+1) ] 

    
    plt.violinplot(length_by_nco, positions=np.arange(1,4), showmeans=True, vert=False)

    for i in range(3):
        plt.annotate('n='+str(len(length_by_nco[i])), (20, i+1+0.2), fontsize=18)
    
    plt.xlabel('SC length ($\mu$m)')
    plt.ylabel('CO number')
    plt.xticks([20, 40,60,80])
    plt.yticks([1,2,3])
    plt.xlim(15,85)
    
    plt.savefig(output_prefix+"nco_vs_length.svg")

    with open(output_prefix+'nco_vs_length_data.csv', 'w') as f:
        f.write('CO_number,SC_length\n')
        for i in range(3):
            for v in length_by_nco[i]:
                f.write(str(i+1) + ','+str(v) + '\n')
            
    mean_co = np.mean([len(c) for c in co])

    if centre_plot:
        plt.figure()
        plot_centre(all_data_1)
        plt.savefig(output_prefix+'thinned_centre.svg')


    plt.figure()
    plot_diff_all(all_data_1)
    plt.savefig(output_prefix+'diff.svg')

    
    
    with open('../input_data/H2A.csv') as f:
        reader = csv.reader(f)
        H2A = list(reader)

    with open('../input_data/MHL1_diakinesis.csv') as f:
        reader = csv.reader(f)
        MHL = list(reader)


    H2A = np.array([float(p[1]) for p in H2A[:3]])
    MHL = np.array([float(p[1]) for p in MHL[:3]])

    
    dsb_data = H2A/H2A[0]
    CO_data = MHL

    plt.figure()
    plt.plot(dsb_density, nf_mean)
    plt.plot(dsb_data, CO_data, 'rx')
    
    plt.figure()
    
    plt.bar(np.arange(max_n+2), all_hist_n[:max_n+2])
    plt.xlabel('Number of COs')
    plt.ylabel('Number of bivalents', size=24)
    plt.xticks(np.arange(max_n+2))
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.savefig(output_prefix+'julia_all_hist_n.svg')

    
    plt.figure()#figsize=(10,10))
    
    plt.bar(np.arange(max_n+2), all_hist_n[:max_n+2])
    plt.title('Mean CO number {:.02f}'.format(mean_co))
    plt.xlabel('Number of COs')
    plt.ylabel('Number of bivalents', size=24)
    plt.xticks(np.arange(max_n+2))
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.savefig(output_prefix+'julia_all_hist_n_mean.svg')

    
    plt.figure()

    mid = np.linspace(0,1,11)
    mid = 0.5*(mid[:-1]+mid[1:])
    plt.hist(mid, bins=np.linspace(0,1,11), weights=all_hist, histtype='stepfilled')#, width=1.0/29)
    plt.xlabel('Relative position along bivalent', size=24)
    plt.ylabel('Number of COs')

    plt.xlim(0,1)
    
    plt.savefig(output_prefix+'julia_all_hist.svg')

    
    fig, ax = plt.subplots(2,2)
    for i in range(4):
        ax[i//2, i%2].bar(np.linspace(0,1,11)[:-1], all_hist_nn[i], width=1.0/29)

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
        for c in all_co:
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
        
    data_all = []
    for cell_hei10 in zip(*all_intensities):
        tot_hei10 = np.sum(np.sum(x) for x in cell_hei10)
        norm_hei10 = [ np.array(x)/tot_hei10 for x in cell_hei10]
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



    dsb_density = np.array(dsb_density)
    dsb_density = dsb_density/dsb_density[-1]
    
    dsb_data = H2A/H2A[0]
    CO_data = MHL[:3]/5

    CO_data = CO_data/CO_data[0]

    
    plt.figure()
    nf2 = np.mean(np.array(nf_tot), axis=0)
    CO_sim = nf2/1000.0
    CO_sim = CO_sim/CO_sim[-1]

    
    plt.plot(dsb_density, CO_sim)
    plt.plot(dsb_data, CO_data, 'rx')
    plt.xlabel('Relative RI density')
    plt.ylabel('Relative number of \nCOs per bivalent')
    plt.savefig(output_prefix+'julia_hs.svg')

    


def count_CO(data_path):
    # Load preprocessed pkl file
    new_data = load_data(*data_path)
    
    new_data_keys = [ (idx, to_val(h['L']), to_val(h['density'])) for idx, (h, v) in enumerate(new_data.values()) ]
    idx, L, density = map(list, zip(*new_data_keys))
    
    df = pd.DataFrame(data={'idx':idx, 'L':L, 'density':density})
        
    df2 = df.sort_values(['L', 'density'])

    
    df2['density_idx'] = list(range(3))*5

    
    grouped = df2.groupby(['L'])
    
    new_data_list = list(new_data.values())


    co_tot = 0
    nch_tot = 0

    for group_idx, (name, group) in enumerate(grouped):
        print('name', name)
        s = group.iloc[-1,:]

        h, v = new_data_list[int(s['idx'])]

        
        L = to_val(h['L'])
        density = to_val(h['density'])
        K = to_val(h['K'])

        hist_n, hist_nn, hist, nf, nch, L_array, intensities, orig_foci_intensities, co = summary((h,v), L, K, tp=-1)

        co_tot += nf
        nch_tot += nch

    return co_tot / nch_tot


def main(sim_data_path, output_path):
    make_plots((sim_data_path+'/survey_julia_new_ends/', 'at_'), output_path+'/new_end_', intensity_bins=np.linspace(0,0.3,30))
    make_plots((sim_data_path+'/survey_escape/', 'at_'),  output_path+'/escape_',  intensity_bins=np.linspace(0,0.3,30))

    make_plots((sim_data_path+'/survey_julia_ox/', 'at_'),  output_path+'/ox_new_end_', centre_plot=False, intensity_bins=np.linspace(0,0.12,11), max_n=7)
    make_plots((sim_data_path+'/survey_julia_ux/', 'at_'),  output_path+'/ux_new_end_',  intensity_bins=np.linspace(0,0.3,10), max_n=3)
    make_plots((sim_data_path+'/survey_julia_no_ends/', 'at_'), output_path+'/no_end_', intensity_bins=np.linspace(0,0.3,30))
    make_plots((sim_data_path+'/survey_exp/', 'at_'),  output_path+'/exp_',  intensity_bins=np.linspace(0,0.3,30))

    make_plots((sim_data_path+'/survey_female/', 'at_'),  output_path+'/female_', intensity_bins=np.linspace(0,0.3,30))
    wt_co = count_CO((sim_data_path+'/survey_julia_new_ends/', 'at_'))
    female_co = count_CO((sim_data_path+'/survey_female/', 'at_'))

    with open(output_path+'/summmary_female.txt', 'w') as f:
        print('wt co, female co, ratio', wt_co, female_co, female_co/wt_co, file=f)
    
main(sys.argv[1], sys.argv[2])
                
        
