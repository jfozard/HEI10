
## Code to generate the majority of the data analysis plots

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle

import statsmodels.api as sm

import matplotlib as mpl

from matplotlib.colors import to_rgba_array


# Update matplotlib parameters
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


# Criterion used for identifying peak as a CO - normalized (with mean and s.d.)
# hei10 levels being above 0.4 time maximum peak level
def pc(peaks, alpha=alpha_max, criterion='max'):
    pos, v, l = peaks
    idx = (v>=alpha*np.max(v))
    return np.array(pos[idx]), v[idx], l

# Linear regression best-fit line
# Returns best fit line on range (max(min_x, min(x), max(x))
def lin_fit(x, y, min_x=None, of=None):
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
    
    return X, yy


## Analysis of each dataset 
def analyse_data(A, data):

    all_peak_data = data['all_peak_data']
    o_hei10_traces = data['o_hei10_traces']
    orig_trace_lengths = data['orig_trace_lengths']

    print(len(A), len(all_peak_data), len(o_hei10_traces), len(orig_trace_lengths))
    
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
        median_all = np.mean(np.concatenate([o_hei10_traces[k] for k in range(i,i+5)]))
        std_all = np.std(np.concatenate([o_hei10_traces[k] for k in range(i,i+5)]))
        for j in range(i, i+5):
            p, h, l = all_peak_data[j]
            if len(p):
                pp = pc((p, h, l), alpha_max)

                
                h=np.array(o_hei10_traces[j])

                sig_peaks.append(pp[0])
                
                orig_peak_hei10.append(h[pp[0]])
                
                h=(h-median_all) #/std_all
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
    A['new_peak_hei10']=new_peak_hei10     # HEI10 intensities at peaks, with mean cell HEI10 levels subtracted
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





def plot_data(A, image_output_path, data_fn, max_n, stacked_bins, summary_fn):



    summary_file = open(summary_fn, 'w')

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


    # Histogram of the numbers of signficant peaks on each SC
    def hist_peak_no(df):
        d = []
        for i in df.index:
            p, h, l = df['peaks'][i], df['peak_hei10'][i], df['trace_length'][i]
            #print(l)
            if(len(p)>0):
                pp = pc((p, h, l), alpha_max)
                d.append(len(pp[0]))
            else:
                d.append(0)
        fig = plt.figure()
        h, c = np.unique(d, return_counts=True)
        if 0 not in h:
            h = [0] + list(h)
            c = [0] + list(c)
        print("mean peaks", 5*np.mean(d))
        
        plt.bar(h,c)
        plt.xlabel('Number of COs')
        plt.ylabel('Number of bivalents', size=24)
        plt.xticks(range(max(d)+1))
        return h, c, fig


    # Histograms of the positions of all of the hei10 peaks
    # (positions normalized to length of SC
    # Data can be made symmetric about 0.5
    # Also make stacked histograms (showing the distributions of the 1st / 2nd / 3rd / etc crossovers
    def hist_peak_pos(df, make_symmetric=False, stacked=0, summary_file=None):
        all = []
        for i in df.index:
            p, h, l = df['peaks'][i], df['peak_hei10'][i], df['trace_length'][i]
            #print(l)
            if(len(p)>0):
                pp = pc((p, h, l), alpha_max)
                all.append(list(pp[0]/pp[2]))

        if summary_file is not None:
            print('Peak pos histogram ' + str(len(all)) + ' ' + str(np.sum([len(c) for c in all ])), file=summary_file)

                
        if stacked:
            print(list(zip(*all)))
            d = list(zip(*all))
            fig = plt.figure()
            plt.hist(d, bins=np.linspace(0,1,11), range=(0,1), stacked=True)
            plt.xlabel('Relative position along bivalent', size=24)
            plt.ylabel('Number of COs')
            plt.xlim(0,1)
        else:
            all = [u for c in all for u in c]
            if make_symmetric:
                all = all + [1-p for p in all]

            fig = plt.figure()
            plt.hist(all, bins=np.linspace(0,1,11), range=(0,1), histtype='stepfilled')
            plt.xlabel('Relative position along bivalent', size=24)
            plt.ylabel('Number of COs')
            plt.xlim(0,1)
        return fig


    # Make CO number (e.g. fig2c) , position of all CO (fig 2d) , and positions of COs on SCs with 1/2/3/... CO (fig2 e-g)
    
    def all_summary(df, summary_file):
        all_figs = []
        h, c, fig = hist_peak_no(df)
        all_figs.append(fig)
        all_figs.append(hist_peak_pos(df, make_symmetric=False, summary_file=summary_file))
        for i, x in zip(h,c):
            if x>5:
                sel = np.array(df['num_sig_peaks']) == i
                print('Pos histogram ' + str(i) + ' ' + str(c) + ' ' + str(np.sum(sel)), file=summary_file)
                all_figs.append(hist_peak_pos(df.loc[sel], make_symmetric=False, stacked = i if i>1 else 0))
        return all_figs


    print('Num late SC ', np.sum(A.all_stage=='late'), file=summary_file)

    # Number of cells this data is drawn from for reporting summary

    nc = 0

    for i in range(0, len(A.index), 5):
        if A.all_stage[i] == 'late' and any(A.quality[j]==1 for j in range(i, i+5)):
            nc += 1
            print(A.iloc[i], file=summary_file)
        else:
            print("NOT USED", file=summary_file)
            print(A.iloc[i], file=summary_file)
            print(A.quality[i:i+5], file=summary_file)

    print('SCs from ' + str(nc) + ' cells', file=summary_file)


    # Restrict data to SCs from late cells with good quality traces
    A_late = A.loc[(A.all_stage=='late') & (A.quality==1)]

    # Histogram of spacing between significant peaks (in terms of real length)
    # e.g. Fig 2h
    pos = []
    for i in A_late.index:
        s = np.array(A['sig_peaks'][i])
        LL = A['SC length'][i]
        L = A['trace_length'][i]
        pos += list(np.diff(s)/L*LL)
        print(s)

    plt.figure()
    plt.xlabel("Spacing ($\mu$m)", size=LARGE_FS)
    plt.ylabel("Frequency", size=LARGE_FS)
    data = plt.hist(pos, bins=np.linspace(0, 60, 10), density=True)
    plt.title('N='+str(len(pos)))
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.savefig(image_output_base+'all_spacing.svg')
    plt.close()


    # Find numbers of plants and SCs used
    cs = []
    for i in A_late.index:
        cs.append((A_late.Date[i], A_late.Plant[i]))

    cs = set(cs)

    print(' Data from {} Plants'.format(len(cs)), file=summary_file)

    
    print('Data from {} SC'.format(len(A_late)), file=summary_file)

    print('Mean peaks all', np.mean(A_late.num_sig_peaks), file=summary_file)



    # Generate histograms of CO number per SC,
    # positions of all COs, positions of COs 
    
    figs = all_summary(A_late, summary_file)

    figs[0].savefig(image_output_base+'all_late_hist_n.svg')
    figs[1].savefig(image_output_base+'all_late_all.svg')
    figs[2].savefig(image_output_base+'all_late_1.svg')
    figs[3].savefig(image_output_base+'all_late_2.svg')
    figs[4].savefig(image_output_base+'all_late_3.svg')

    plt.close()

    
    # Histogram of relative CO intensities
    # Peak data from traces with cell mean hei10 level subtracted
    
    plt.figure()
    plt.hist([p for u in A_late['new_rel_peak_int_cell'] for p in u], bins=stacked_bins) 
    plt.autoscale(enable=True, axis='x', tight=True)

    plt.xlabel('Relative CO intensity')
    plt.ylabel('Number of COs')
    plt.savefig(image_output_path+'new_rel_peak_int_cell.svg')
    plt.close()


    # Make histograms of relative peak intensities (cell mean subtracted), stacked according to the number of COs on each SC

    # e.g. fig2l


    data_by_n = [ [ u for u,c in zip(A_late['new_rel_peak_int_cell'], A_late['num_sig_peaks']) if c==i ] for  i in range(1,max_n)] + [ [ u for u,c in zip(A_late['new_rel_peak_int_cell'], A_late['num_sig_peaks']) if c>=max_n ] ]

    print(list(map(len, data_by_n)))


    data_stacked = [ [ p   for u in d for p in u ] for d in data_by_n]

    bin_labels = [str(i) for i in range(1,max_n)] + [str(max_n)+'+']

    
    plt.figure()
    plt.hist(data_stacked, bins=stacked_bins, stacked=True, label=bin_labels)
    plt.legend()
    plt.xlabel('Relative CO intensity')
    plt.ylabel('Number of COs')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.savefig(image_output_base+'new_rel_peak_int_cell_stacked.svg')
    plt.close()



    # Filter to late data with two COs

    
    A_late2 = A_late[A_late.num_sig_peaks==2]


    # Fig 2k 
    
    data_c = []
    data_x = []
    for i in A_late2.index:

        c = np.sum(A_late2['sig_peaks'][i])/2/A_late2['trace_length'][i]
        x = A_late2['new_peak_hei10'][i][0]/np.sum(A_late2['new_peak_hei10'][i])
        data_c.append(c)
        data_x.append(x)

    print('Linear regression new_centre', file=summary_file)
    X, yy = lin_fit(data_c, data_x, of=summary_file)


    plt.figure()
    plt.title('N='+str(len(data_x)))
    plt.plot(data_c, data_x, 'rx')
    plt.plot(X, yy, 'b')
    plt.xlabel('Midpoint relative position')
    plt.ylabel('Relative intensity left CO', size=24)
    plt.savefig(image_output_base+'new_centre.svg')
    plt.close()



    # fig 1b right panel

    # Normalize sum of hei10 peaks for each cell
    norm_hei10 = []
    norm_peaks = dict((i,[]) for i in range(1,4))
    for i in range(0, len(A), 5):
        if A.Stage[i]=='late' and A.all_good[i]=='y':
            hei10 = []
            hei10_all = []
            for j in range(i, i+5):
                hei10.append(np.sum(A['new_peak_hei10'][j]))
                hei10_all.append(A['new_peak_hei10'][j])
            hei10 = np.array(hei10)
            sc_hei10 = hei10/np.sum(hei10)
            tot_hei10 = np.sum(hei10)
            lengths = A['SC length'][i:i+5]
            lengths = lengths/np.sum(lengths)
            num_peaks = A['num_sig_peaks'][i:i+5]
            for j, h in zip(num_peaks, hei10_all):
                #print(j, h)
                if j in norm_peaks:
                    norm_peaks[j]+=list(np.array(h)/tot_hei10)

    label_map = { 1:'single', 2:'double', 3:'triple'}

    plt.figure()
    for i in range(1,4):
        plt.hist(norm_peaks[i], histtype='step', label=label_map[i], density=True, lw=2)

        print('Single double triple new ' + str(i) + ' ' + str(len(norm_peaks[i])), file=summary_file)
        
    plt.xlabel('Relative HEI10 focus intensity')
    plt.ylabel('Frequency density')
    plt.legend()
    plt.savefig(image_output_base+'single_double_triple_peak_new_intensities.svg')

    
  
    # fig 1b middle two panels
    label_map = { 1:'single', 2:'double', 3:'triple'}

    for prefix in ['new']:
        for num_co, title, out_fn in [ ((-1, 10000), 'All late SC', '_rel_length.svg'), ((2,2), 'All late SC with double CO', '_rel_length2.svg' ) ]:
        
            # Normalize sum of hei10 peaks for each cell
            plt.figure()
            norm_lengths = []
            norm_hei10 = []
            for i in range(0, len(A), 5):
                if A.Stage[i]=='late' and A.all_good[i]=='y':
                    hei10 = []
                    for j in range(i, i+5):
                        hei10.append(np.sum(A[prefix+'_peak_hei10'][j]))
                    hei10 = np.array(hei10)
                    hei10 = hei10/np.sum(hei10)
                    lengths = A['SC length'][i:i+5]
                    lengths = lengths/np.sum(lengths)
                    num_peaks = A['num_sig_peaks'][i:i+5]
                    l, u = num_co
                    norm_lengths+= list(lengths[(l<=num_peaks) & (num_peaks<=u)])
                    norm_hei10 += list(hei10[(l<=num_peaks) & (num_peaks<=u)])

                        
            plt.plot(norm_lengths, norm_hei10, 'ro')
            print(title+' Relative bivalent length vs Relative total focus HEI10 per bivalent', file=summary_file)
            yy = lin_fit(norm_lengths, norm_hei10, min_x=0, of=summary_file)
            plt.plot(*yy,'r-')
#            plt.title(title)
            plt.xlim(xmin=0)
            if num_co!=2:
                plt.ylim(ymin=0)

            plt.xlabel('Relative bivalent length')
            plt.ylabel('Relative total focus HEI10 per bivalent')
#            plt.tight_layout()
            plt.savefig(image_output_base+prefix+out_fn)
    


    # Fig 1b left panel        
    plt.figure()

    for stage, criterion, col in [ ('late', lambda s: s=='late', 'r'),
                              ('mid', lambda s: s=='mid' or s=='Mid', 'b'),
                                   ('early', lambda s: s=='early', 'g') ]:

            
        norm_lengths = []
        norm_hei10 = []
        for i in range(0, len(A), 5):
            print(i, A.Stage[i], A.all_good[i])
            if criterion(A.Stage[i]):# and A.all_good[i]=='y':
                hei10 = []
                for j in range(i, i+5):

                    h = o_hei10_traces[j]
                    #h = np.maximum(h - np.median(h), 0)
                    h = h - np.median(h)
                    hei10.append(np.sum(h))
                hei10 = np.array(hei10)
                hei10 = hei10/np.sum(hei10)
                lengths = A['SC length'][i:i+5]
                lengths = lengths/np.sum(lengths)
                num_peaks = A['num_sig_peaks'][i:i+5]
                norm_lengths+= list(lengths)
                norm_hei10 += list(hei10)
                    
        print(stage, np.mean(norm_lengths), np.std(norm_lengths), file=summary_file)
        plt.plot(norm_lengths, norm_hei10, col+'o', markersize=2)

        print(stage + '  Relative bivalent length vs Relative total focus HEI10 per bivalent', file=summary_file)

        
        yy = lin_fit(norm_lengths, norm_hei10, min_x=0, of=summary_file)
        plt.plot(*yy,col+'-', label=stage)
#    plt.title('Total HEI10 above median')
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.xlabel('Relative bivalent SC length')
    plt.ylabel('Relative total HEI10 per SC')
    plt.legend()
 #   plt.tight_layout()

    plt.savefig(image_output_base+'rel_length_tot_compare.svg')



## Plot hei10 (mean) intensity traces along a specific SC
def plot_traces(data_fn, csv_fn, idx_list, image_output_path):
    A = pd.read_csv(csv_fn)

    with open(data_fn, 'rb') as f:
        data = pickle.load(f)
        all_peak_data, orig_trace_lengths = data['all_peak_data'], data['orig_trace_lengths']
        o_hei10_traces = data['o_hei10_traces']  # Don't have dapi data for UX.
    for i in idx_list:
        plt.figure()
        v = np.array(o_hei10_traces[i])
        n = len(v)
        plt.plot(np.linspace(0,1,n), v/65535)
        plt.savefig(image_output_path+f'trace_{i}.png')
        plt.close()
        plt.figure()
        plt.plot(np.linspace(0,1,n), v/65535)
        plt.plot(np.linspace(0,1,n), np.median(v)*np.ones_like(v)/65535)
        plt.savefig(image_output_path+f'trace_median_{i}.png')
        plt.close()
        print(A.iloc[i])

        with open(image_output_path+f'trace_{i}.txt', 'w') as f:
            for u in v:
                f.write(f'{u:.02f}\n')
    

        

data_output_path = 'data_output/'
image_output_path='data_output/' 


for image_output_base, csv_fn, data_fn, max_n, stacked_bins, summary_fn in  [
        ( image_output_path, '200406.csv', data_output_path+'test.pkl', 4, np.linspace(0, 0.3, 30), data_output_path+'summary.txt'),
        ( image_output_path+'ox_', 'OX.csv', data_output_path+'test_ox.pkl', 4, np.linspace(0, 0.12, 11), data_output_path+'summary_ox.txt'),
        ( image_output_path+'ux_', 'UX.csv', data_output_path+'test_ux.pkl', 3, np.linspace(0, 0.3, 10), data_output_path+'summary_ux.txt'),
        ]:
    A = pd.read_csv(csv_fn)
    plot_data(A, image_output_base, data_fn, max_n, stacked_bins, summary_fn)


# Plots of HEI10 traces along individual chromosomes (used for panels in Fig1A, 3A and 4A)
plot_traces(data_output_path+'test.pkl', '200406.csv', [675, 664, 492 ], image_output_path)
plot_traces(data_output_path+'test_ox.pkl', 'OX.csv', [ 260], image_output_path)
plot_traces(data_output_path+'test_ux.pkl', 'UX.csv', [ 115], image_output_path)
