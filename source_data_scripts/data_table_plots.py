
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

#utils.chooseCRANmirror(ind=1)
#utils.install_packages('dunn.test')
#utils.install_packages('kruskal')

import rpy2.robjects as robjects

dunn = importr('dunn.test')
#kruskal_r = importr('kruskal.test')
stats = importr('stats')

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

from matplotlib import lines

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

## Analysis of each dataset 
def analyse_data(A):



    #print(len(A), len(all_peak_data), len(o_hei10_traces), len(orig_trace_lengths))

    peaks = []
    orig_peak_hei10 = []
    new_peak_hei10 = []

    A = A.rename(columns={'Stage':'all_stage', 'num_foci':'num_sig_peaks', 'foci_pos':'sig_peaks', 'trace_pixel_length':'trace_length'})
    

    ### For each cell in turn, process the 

    for i in range(len(A)):
        s = A.iloc[i]
        nf = int(s['num_sig_peaks'])
        pos = [ s[f'pos_{i}'] for i in range(nf) ]
        hei10 = [ s[f'HEI10_{i}'] for i in range(nf) ]

        median_hei10 = s['trace_median_hei10']
        new_hei10 = np.array(hei10) - median_hei10


        peaks.append(pos)
        orig_peak_hei10.append(hei10)
        new_peak_hei10.append(new_hei10)


    
    A['sig_peaks']=peaks
    A['peak_hei10']=orig_peak_hei10

    A['new_peak_hei10']=new_peak_hei10     # HEI10 intensities at peaks, with median cell HEI10 levels subtracted
    A['orig_peak_hei10']=orig_peak_hei10   # Raw HEI10 intensities at peaks



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





def plot_data(A, image_output_path, max_n, stacked_bins, summary_fn, summary_fn2):

    print(list(A))


    summary_file = open(summary_fn, 'w')

    # Edit CSV columns + naming


    # Initial processing of data for each cell
    A = analyse_data(A)

    print(list(A))

    

    # Histogram of the numbers of signficant peaks on each SC
    def hist_peak_no(df):
        d = df['num_sig_peaks']
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
            p, l = df['sig_peaks'][i], df['trace_length'][i]

            if len(p)>0:
                all.append(list(np.array(p)/l))

                
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
            #print(A.iloc[i], file=summary_file)
            print("USED", i, A.all_stage[i], 'All quality:', A.all_quality[i], 'All good:', A.all_good[i], file=summary_file)
        else:
            print("NOT USED", i, A.all_stage[i], 'All quality: ', A.all_quality[i], 'All good: ', A.all_good[i], file=summary_file)
#            print(A.iloc[i], file=summary_file)
#            print(A.quality[i:i+5], file=summary_file)

    print('SCs from ' + str(nc) + ' cells', file=summary_file)



    # Restrict data to SCs from late cells with good quality traces
    A_late = A.loc[(A.all_stage=='late') & (A.quality==1)]

    print('Length of A_late', len(A_late))
                  
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
        cs.append((A_late.Genotype[i], A_late.Date[i], A_late.Plant[i]))

    cs = list(set(cs))
    cs.sort()

    print('Plant names', cs, file=summary_file)
    
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
    X, yy, r2 = lin_fit(data_c, data_x, of=summary_file, r2=True)


    plt.figure()
    plt.plot(data_c, data_x, 'rx')
    plt.plot(X, yy, 'b')
    plt.title('N='+str(len(data_x)) + f'    $R^2$ = {r2:.2f}')
    #    plt.text(0.1, 0.9, , transform = plt.gca().transAxes)
    #plt.xlabel('Bivalent fraction closest to CO')
    plt.xlabel('Midpoint relative position')
    plt.ylabel('Left CO relative intensity', size=22)
    plt.savefig(image_output_base+'new_centre.svg')
    plt.close()


    data_c = []
    data_x = []
    for i in A_late2.index:


        c = A_late2['sig_peaks'][i]/A_late2['trace_length'][i]
        c = np.minimum(c, 1-c) 
        
        x = A_late2['new_peak_hei10'][i]/np.sum(A_late2['new_peak_hei10'][i])
        data_c += list(c)
        data_x += list(x)

    print('Linear regression new_centre2', file=summary_file)
    X, yy = lin_fit(data_c, data_x, of=summary_file)


    plt.figure()
    plt.title('N='+str(len(data_x)))
    plt.plot(data_c, data_x, 'rx')
    plt.plot(X, yy, 'b')
    plt.xlabel('Bivalent fraction closest to CO')
    plt.ylabel('Left CO relative intensity', size=22)
    plt.savefig(image_output_base+'new_centre2.svg')
    plt.close()


    # fig 1b right panel

    # Normalize sum of hei10 peaks for each cell
    ncells = 0
    norm_hei10 = []
    norm_peaks = dict((i,[]) for i in range(1,4))
    for i in range(0, len(A), 5):
        if A.all_stage[i]=='late' and A.all_good[i]=='y':
            ncells +=1
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

    print('single double tripple', [len(norm_peaks[i]) for i in norm_peaks], 'ncells', ncells)
        
    plt.xlabel('Relative HEI10 focus intensity')
    plt.ylabel('Frequency density')
    plt.legend()
    plt.savefig(image_output_base+'single_double_triple_peak_new_intensities.svg')


    # Write normalized peaks into text file

    with open(image_output_base+'single_double_triple.json', 'w') as of123:
        json.dump(norm_peaks, of123)
    


    

        
    h_max = max([max(u) for u in norm_peaks.values()])
    fig, ax = plt.subplots(3,1)
    fig.subplots_adjust(left=0.15, right=0.8, bottom=0.15, top=0.9)
    y_pos = []
    for i in range(3):
        ax[i].hist(norm_peaks[i+1], range=(0, h_max))#, histtype='step', density=True, lw=2, )
        #ax[i].set_ylim(0, 20)
        #ax[i].set_yticks([0,20])
        if i<2:
            ax[i].set_xticks([])
        ax[i].tick_params('y', labelsize=15)
        ax[i].text(0.8, 0.7, label_map[i+1], size=14,  transform=ax[i].transAxes)
        ax[i].text(0.8, 0.5, 'n={}'.format(len(norm_peaks[i+1])), size=14,  transform=ax[i].transAxes)
        #ax[i].set_ylabel(label_map[i+1], size=20)
        pp = ax[i].get_position()
        y_pos.append(0.5*(pp.y0+pp.y1))
        print(pp.x1, 0.5*(pp.y0+pp.y1))


    data = rpy2.robjects.ListVector([('single', robjects.FloatVector(norm_peaks[1])),
                                     ('double', robjects.FloatVector(norm_peaks[2])),
                                     ('triple', robjects.FloatVector(norm_peaks[3])) ] )
                                     

    out = dunn.dunn_test(data, method='bonferroni', list='TRUE', altp='TRUE', kw='TRUE')
    print(out, file=summary_file)

    out2 = stats.kruskal_test(data)
    print(out2[0], file=summary_file)

    
    chi2, Z, pvals, pvals_adj, compare = out
    
        
    for i0, i1, idx in [[0,1, 0], [1,2, 2]]:
        l = lines.Line2D([0.81, 0.83, 0.83, 0.81], [y_pos[i0]-0.02, y_pos[i0]-0.02, y_pos[i1]+0.02, y_pos[i1]+0.02], transform = fig.transFigure, lw=2, color='k')

        Z_val = Z[idx]
        p_val = pvals[idx]
        p_adj = pvals_adj[idx]
        c = compare[idx]
        
#        _, p_val, _ = ttest_ind(norm_peaks[i0+1], norm_peaks[i1+1])
#        u_val, p_val = mannwhitneyu(norm_peaks[i0+1], norm_peaks[i1+1], alternative='two-sided')

#        print('ttest_ind', i0+1, i1+1, ttest_ind(norm_peaks[i0+1], norm_peaks[i1+1]), 'n=', len(norm_peaks[i0+1]), len(norm_peaks[i1+1]), file=summary_file)
#        print('mann_whitney', i0+1, i1+1, mannwhitneyu(norm_peaks[i0+1], norm_peaks[i1+1]), 'n=', len(norm_peaks[i0+1]), len(norm_peaks[i1+1]), file=summary_file)

        
        fig.text(0.87, 0.5*(y_pos[i0]+y_pos[i1]), f'p={p_adj:.1e}', fontsize=15, rotation=90, ha='center', va='center')
        #fig.text(0.88, 0.5*(y_pos[i0]+y_pos[i1]), f'Z={Z_val:.1f} ', fontsize=8, rotation=90, ha='center', va='center')
        fig.lines.append(l)

    for i0, i1, idx in [[0,2, 1]]:
        l = lines.Line2D([0.81, 0.91, 0.91, 0.81], [y_pos[i0]+0.02, y_pos[i0]+0.02, y_pos[i1]-0.02, y_pos[i1]-0.02], transform = fig.transFigure, lw=2, color='k')
        #_, p_val, _ = ttest_ind(norm_peaks[i0+1], norm_peaks[i1+1])
#        u_val, p_val = mannwhitneyu(norm_peaks[i0+1], norm_peaks[i1+1], alternative='two-sided')

        Z_val = Z[idx]
        p_val = pvals[idx]
        p_adj = pvals_adj[idx]
        c = compare[idx]

        
        fig.text(0.95, 0.5*(y_pos[i0]+y_pos[i1]), f'p={p_adj:.1e}', fontsize=15, rotation=90, ha='center', va='center')
#        fig.text(0.96, 0.5*(y_pos[i0]+y_pos[i1]), f'Z={Z_val:.2f}', fontsize=8, rotation=90, ha='center', va='center')
        fig.lines.append(l)

 #       print('ttest_ind', i0+1, i1+1, ttest_ind(norm_peaks[i0+1], norm_peaks[i1+1]), 'n=', len(norm_peaks[i0+1]), len(norm_peaks[i1+1]), file=summary_file)
#        print('mann_whitney', i0+1, i1+1, mannwhitneyu(norm_peaks[i0+1], norm_peaks[i1+1]), 'n=', len(norm_peaks[i0+1]), len(norm_peaks[i1+1]), file=summary_file)
        
    ax= fig.add_axes([0.15, 0.15, 0.7, 0.8], frame_on=False)

    plt.tick_params(labelcolor="none", bottom=False, left=False)
    #ax.set_xlabel('Relative HEI10 focus intensity', size=20)
    #ax.set_ylabel('Frequency density', size=20)
    plt.savefig(image_output_base+'single_double_triple_peak_new_intensities_stacked.svg')

    print('single double triple kruskal', kruskal(norm_peaks[1], norm_peaks[2], norm_peaks[3]), file=summary_file)
    print('single double triple anova', f_oneway(norm_peaks[1], norm_peaks[2], norm_peaks[3]), file=summary_file)

        
  
    # fig 1b middle two panels
    label_map = { 1:'single', 2:'double', 3:'triple'}

    for prefix in ['new']:
        for num_co, title, out_fn in [ ((-1, 10000), 'All late SC', '_rel_length.svg'), ((2,2), 'All late SC with double CO', '_rel_length2.svg' ) ]:
        
            # Normalize sum of hei10 peaks for each cell
            plt.figure()
            norm_lengths = []
            norm_hei10 = []
            for i in range(0, len(A), 5):
                if A.all_stage[i]=='late' and A.all_good[i]=='y':
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
            xx, yy, r2 = lin_fit(norm_lengths, norm_hei10, min_x=0, of=summary_file, r2=True)
            plt.plot(xx, yy,'r-')
#            plt.title(title)
            plt.xlim(xmin=0)
            if num_co!=2:
                plt.ylim(ymin=0)
            plt.text(0.1, 0.9, f'$R^2$ = {r2:.2f}', transform = plt.gca().transAxes)

#            plt.xlabel('Relative bivalent length')
#            plt.ylabel('Relative total focus HEI10 per bivalent')
#            plt.tight_layout()
            plt.savefig(image_output_base+prefix+out_fn)
    


    # Fig 1b left panel        
    plt.figure()


    used_sc = {}
    for stage, criterion, col in [ ('late', lambda s: s=='late', 'r'),
                              ('mid', lambda s: s=='mid' or s=='Mid', 'b'),
                                   ('early', lambda s: s=='early', 'g') ]:

            
        norm_lengths = []
        norm_hei10 = []
        used_sc[stage] = 0
        for i in range(0, len(A), 5):
            print(i, A.all_stage[i], A.all_good[i])
            if criterion(A.all_stage[i]):# and A.all_good[i]=='y':


                hei10 = []
                hei10_all = []
                
                
                for j in range(i, i+5):

#                    h = o_hei10_traces[j]
                    #h = np.maximum(h - np.median(h), 0)
#                    h = h - np.median(h) #median_all

                    #h = h-median_all
                    #h = np.maximum(h - median_all, 0)
#                    hei10.append(np.sum(h))

                    hei10_background_removed_sum = A['trace_sum_hei10'][j] - A['trace_median_hei10'][j]*A['trace_length'][j]

                    hei10.append(hei10_background_removed_sum)
                    used_sc[stage] +=1

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
        
        xx, yy, r2 = lin_fit(norm_lengths, norm_hei10, min_x=0, of=summary_file, r2=True)
#        plt.plot(xx,yy,col+'-', label=stage + f'  $R^2$ = {r2:.2f}' )
        plt.plot(xx,yy,col+'-', label=stage + f' {r2:.2f}' )

    #    plt.title('Total HEI10 above median')
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
#    plt.xlabel('Relative bivalent SC length')
#    plt.ylabel('Relative total HEI10 per SC')
    plt.legend(fontsize=18)
 #   plt.tight_layout()

    plt.savefig(image_output_base+'rel_length_tot_compare.svg')


    print('rel_length_tot_compare SC numbers', used_sc, file=summary_file)

    
    ### Calulate late total foci number per cell and total SC length per cell.


    cell_idx = []
    cell_total_SC = []
#    cell_num_peaks = []

    sc_length_all = []
    
    for i in range(0, len(A), 5):
        if A.all_stage[i]=='late' and A.all_good[i]=='y':
            total_SC = 0
#            total_peaks = 0
            sc_length = []
            for j in range(i, i+5):
                total_SC += A['SC length'][j]
#                total_peaks += A['num_sig_peaks'][j]
                sc_length.append(A['SC length'][j])
            sc_length.sort()
            cell_idx.append(i//5)
            cell_total_SC.append(total_SC)
#            cell_num_peaks.append(total_peaks)

            sc_length_all.append(sc_length)



    cell_total_SC_all = []
    
    for i in range(0, len(A), 5):
        if True: #A.Stage[i]=='late':
            total_SC = 0
#            total_peaks = 0
            for j in range(i, i+5):
                total_SC += A['SC length'][j]
            cell_total_SC_all.append(total_SC)

    cell_num_peaks = []

    for i in range(0, len(A), 5):
        if A.all_stage[i]=='late':
            total_peaks = 0
            for j in range(i, i+5):
                total_peaks += A['num_sig_peaks'][j]
            cell_num_peaks.append(total_peaks)


            
    

    summary_file2 = open(summary_fn2, 'w')
    summary_file2.write('cell_idx, total_SC, num_oeaks\n')
    for d in zip(cell_idx, cell_total_SC, cell_num_peaks):
       # summary_file2.write(','.join(map(str, d))+'\n')
         summary_file2.write('{}, {}, {}\n'.format(d[0], d[1], d[2]))


    print('Cell SC length mean / std', np.mean(cell_total_SC), np.std(cell_total_SC, ddof=1), len(cell_total_SC))
    print('Cell SC length all  include bad mean / std', np.mean(cell_total_SC_all), np.std(cell_total_SC_all, ddof=1), len(cell_total_SC_all))

    print('Cell late CO number mean / std', np.mean(cell_num_peaks),  np.std(cell_num_peaks, ddof=1), len(cell_num_peaks))
    print('Mean SC lengths', np.mean(sc_length_all, axis=0), np.var(sc_length_all, axis=0, ddof=1), len(sc_length_all))

    print('Cell SC length all include bad mean / std (using n-1)', np.mean(cell_total_SC_all), np.std(cell_total_SC_all, ddof=1), len(cell_total_SC_all), file=summary_file)

    print('Cell SC length mean / std (using n-1)', np.mean(cell_total_SC), np.std(cell_total_SC, ddof=1), len(cell_total_SC), file=summary_file)
    print('Cell late CO number mean / std (using n-1)', np.mean(cell_num_peaks),  np.std(cell_num_peaks, ddof=1), len(cell_num_peaks), file=summary_file)
    print('Mean SC lengths', np.mean(sc_length_all, axis=0), np.std(sc_length_all, axis=0, ddof=1), len(sc_length_all), file=summary_file)
    

        

input_data_path = '../input_data/'
source_data_path = '../source_data/'
test_output_path = '../source_data_test/'


os.makedirs(test_output_path, exist_ok=True)


for image_output_base, csv_fn, max_n, stacked_bins, summary_fn, summary_fn2 in  [
        ( test_output_path, source_data_path+'fig2_cytology.csv', 4, np.linspace(0, 0.3, 30), test_output_path+'summary.txt', test_output_path+'cell_summary.txt'),
	( test_output_path+'ox_', source_data_path+'fig3_cytology.csv',  4, np.linspace(0, 0.12, 11), test_output_path+'summary_ox.txt', test_output_path+'cell_summary_ox.txt'),
        ( test_output_path+'ux_', source_data_path+'fig4_cytology.csv',  3, np.linspace(0, 0.3, 10), test_output_path+'summary_ux.txt', test_output_path+'cell_summary_ux.txt'),


]:
    A = pd.read_csv(csv_fn)
    plot_data(A, image_output_base,  max_n, stacked_bins, summary_fn, summary_fn2)

