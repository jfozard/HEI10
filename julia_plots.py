import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
import random


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

LARGE_FS=32

import pickle

import statsmodels.api as sm

def pc(v, alpha=0.4):
    idx = (v>=alpha*np.max(v))
    return idx


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


import os

def summary(data, L, K, tp=-1, max_k=4):

    print('tp', tp)
    head, all_data = data

    all_foci = []
    foci_intensities = []
    L_array = []

    orig_foci_intensities = []
    ns = 0
    for L, pos, hei10 in all_data:
        if len(pos) and len(hei10) and tp<len(hei10) and len(pos)==len(hei10[tp]):
            idx = pc(hei10[tp])
            all_foci.append(pos[idx]/L)
            foci_intensities.append(hei10[tp][idx]/np.sum(hei10[tp][idx]))
            orig_foci_intensities.append(hei10[tp][idx])
            ns += 1
            L_array.append(L)

    
    foci = [ p for l in all_foci for p in l]

    print('n ch', len(all_foci), 'n f', len(foci))

    n_foci = [len(l) for l in all_foci]

    print('n_foci', n_foci)
    
    print("mean peaks", np.mean(n_foci))

   
    hist_n, _ = np.histogram(n_foci, range(20))

    print(hist_n)
    
    hist_nn = {}
    
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


    hist, _ = np.histogram(foci, bins = np.linspace(0,1,11))

    return hist_n, hist_nn, hist, len(foci), len(all_foci), L_array, foci_intensities, orig_foci_intensities, all_foci



def lin_fit(x, y):
    X = np.array(x)
    Y = np.array(y)
    X2 = sm.add_constant(X)
    est = sm.OLS(y, X2)
    est2 = est.fit()
    print(est2.summary())
    X = np.linspace(np.min(X), np.max(X), 100)
    X2 = sm.add_constant(X)
    yy = est2.predict(X2)
    return X, yy



def plot_centre(all_data, tp=-1):
    
    c = []
    v = []
    
    for L, pos, hei10 in all_data:
        
        if len(pos) and len(hei10) and tp<len(hei10) and len(pos)==len(hei10[tp]):
            if np.max(hei10[tp])>2*np.mean(hei10[tp]):
                idx = pc(hei10[tp])
                idx = np.where(idx)[0]
                
                if len(idx)==2:
                    c.append(np.sum(pos[idx])/2/L)
                    v.append(hei10[tp][idx[0]]/np.sum(hei10[tp][idx]))


    idx = list(range(len(c)))
    s = random.sample(idx, 206)
    c = np.array(c)
    v = np.array(v)
    
    plt.plot(c[s], v[s], 'rx')
    X, yy = lin_fit(c[s], v[s])
    plt.plot(X, yy, 'b-')
    plt.title("N="+str(len(c[s]))+" (sample)")
    plt.xlabel("Midpoint relative position")
    plt.ylabel("Relative intensity left CO", size=24)
    plt.autoscale(enable=True, axis='x', tight=True)

    

def plot_diff(all_data, order, tp=-1):
    

    d = []
    
    for L, pos, hei10 in all_data:
        
        if len(pos) and len(hei10) and tp<len(hei10) and len(pos)==len(hei10[tp]):
            if np.max(hei10[tp])>2*np.mean(hei10[tp]):
                idx = pc(hei10[tp])
                idx = np.where(idx)[0]
                
                if len(idx)==order:
                    p = pos[idx]/L
                    d += list(np.diff(p))
    plt.figure()
    plt.hist(d)
    plt.title("{:.4g} {:.4g}".format(np.mean(d), np.std(d)))
    plt.autoscale(enable=True, axis='x', tight=True)
    
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



def load_data(survey_dir, survey_base):
    
    files = os.listdir(survey_dir)

    all_data_ext = {}
    for f in files:
        if len(f)>=len(survey_base) and f[:len(survey_base)] == survey_base:
            print(f)
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
        p = get_class(k)[0:8] + get_class(k)[9:]
        if p in new_data:
            new_data[p][1] += all_data_ext[k][1]
        else:
            h, d = all_data_ext[k]
            new_data[p] = [h, d]

    print(list(new_data))
    
    return new_data


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


def make_plots(data_path, output_prefix, centre_plot=True, intensity_bins=None, max_n=4):
    # Load preprocessed pkl file
    new_data = load_data(*data_path) 
    k = next(iter(new_data))
    print(k, new_data[k][0])
    print('loaded data')


    print(list(new_data.values())[0])
    
    new_data_keys = [ (idx, to_val(h['L']), to_val(h['density'])) for idx, (h, v) in enumerate(new_data.values()) ]
    idx, L, density = map(list, zip(*new_data_keys))
    
    df = pd.DataFrame(data={'idx':idx, 'L':L, 'density':density})
        
    df2 = df.sort_values(['L', 'density'])

    print(df2)
    
    df2['density_idx'] = list(range(3))*5

    print(df2)


    
    grouped = df2.groupby(['L'])
    for name, group in grouped:
        print(name)
        print(group)
        print('***')
    
    new_data_list = list(new_data.values())



    all_hist_n = np.zeros((19,))
    all_hist_nn = np.zeros((4, 10))
    all_hist = np.zeros(10,)
    
    all_intensities = []
    
    nf_mean = np.zeros(3)
    
    mid = np.linspace(0,1,11)
    mid = 0.5*(mid[:-1]+mid[1:])

    nf_tot = []

    all_data_1 = []
    all_co = []

    for group_idx, (name, group) in enumerate(grouped):
        print('name', name)
        s = group.iloc[-1,:]
#        print('s', s)
        print("s['idx']", s['idx'])

        h, v = new_data_list[int(s['idx'])]

        all_data_1 += v
        
        print('h, ', h)
        
        L = to_val(h['L'])
        density = to_val(h['density'])
        K = to_val(h['K'])

        print('L, density, K', L, density, K)
                
        hist_n, hist_nn, hist, nf, nch, L_array, intensities, orig_foci_intensities, co = summary((h,v), L, K, tp=-1)
        
        all_intensities.append(orig_foci_intensities)
        
        all_hist_n += np.array(hist_n)
        all_hist += np.array(hist)

        all_hist_nn += np.array(list(hist_nn.values()))

        all_co += co

          
        print(hist, len(hist))

        
        dsb_density = []
        nf_all = []
        for i in range(3):
            s = group.iloc[i,:]
            print(s)
            j = int(s['idx'])
            dsb_density.append(s['density'])
            _, _, _, nf, nch, _, _, _, _ = summary(new_data_list[j],  L, K, tp=-1)
            nf_all.append(nf/float(nch))
            nf_mean[i] += nf/float(nch)
        nf_tot.append(np.array(nf_all))
        


    mean_co = np.mean([len(c) for c in co])

    if centre_plot:
        plt.figure()
        plot_centre(all_data_1)
        plt.savefig(output_prefix+'thinned_centre.svg')


    plt.figure()
    plot_diff_all(all_data_1)
    plt.savefig(output_prefix+'diff.svg')

    
    
    with open('H2A.csv') as f:
        reader = csv.reader(f)
        H2A = list(reader)

    with open('MHL1_diakinesis.csv') as f:
        reader = csv.reader(f)
        MHL = list(reader)


    H2A = np.array([float(p[1]) for p in H2A[:3]])
    MHL = np.array([float(p[1]) for p in MHL[:3]])

    print(H2A, MHL)
    
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

    print(all_hist_n)
    print('mean peaks all', np.sum(all_hist_n*np.arange(19))/np.sum(all_hist_n), np.sum(all_hist_n*np.arange(19))/np.sum(all_hist_n)*5)

    print('MHL', MHL)
    
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


    print(H2A, MHL)

    dsb_density = np.array(dsb_density)
    dsb_density = dsb_density/dsb_density[-1]
    
    dsb_data = H2A/H2A[0]
    CO_data = MHL[:3]/5

    CO_data = CO_data/CO_data[0]


    print('nf tot', nf_tot)
    
    plt.figure()
    nf2 = np.mean(np.array(nf_tot), axis=0)
    CO_sim = nf2/1000.0
    CO_sim = CO_sim/CO_sim[-1]


    print(dsb_density, CO_sim)
    
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
    for name, group in grouped:
        print(name)
        print(group)
        print('***')
    
    new_data_list = list(new_data.values())


    co_tot = 0
    nch_tot = 0

    for group_idx, (name, group) in enumerate(grouped):
        print('name', name)
        s = group.iloc[-1,:]
#        print('s', s)
        print("s['idx']", s['idx'])

        h, v = new_data_list[int(s['idx'])]

        print('h, ', h)
        
        L = to_val(h['L'])
        density = to_val(h['density'])
        K = to_val(h['K'])

        print('L, density, K', L, density, K)
                
        hist_n, hist_nn, hist, nf, nch, L_array, intensities, orig_foci_intensities, co = summary((h,v), L, K, tp=-1)

        co_tot += nf
        nch_tot += nch

    return co_tot / nch_tot


def main(sim_data_path, output_path):
    make_plots((sim_data_path+'/survey_julia_new_ends/', 'at_'), output_path+'/new_end_', intensity_bins=np.linspace(0,0.3,30))
    make_plots((sim_data_path+'/survey_julia_ox/', 'at_'),  output_path+'/ox_new_end_', centre_plot=False, intensity_bins=np.linspace(0,0.12,11), max_n=7)
    make_plots((sim_data_path+'/survey_julia_ux/', 'at_'),  output_path+'/ux_new_end_',  intensity_bins=np.linspace(0,0.3,10), max_n=3)
    make_plots((sim_data_path+'/survey_julia_no_ends/', 'at_'), output_path+'/no_end_', intensity_bins=np.linspace(0,0.3,30))
    make_plots((sim_data_path+'/survey_exp/', 'at_'),  output_path+'/exp_',  intensity_bins=np.linspace(0,0.3,30))

    make_plots((sim_data_path+'/survey_female/', 'at_'),  output_path+'/female_', intensity_bins=np.linspace(0,0.3,30))
    wt_co = count_CO((sim_data_path+'/survey_julia_new_ends/', 'at_'))
    female_co = count_CO((sim_data_path+'/survey_female/', 'at_'))


    print('wt co, female co, ratio', wt_co, female_co, female_co/wt_co)
    
main(sys.argv[1], sys.argv[2])
                
        
