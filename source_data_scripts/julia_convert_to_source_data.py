import numpy as np
import sys
import csv


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
            if np.sum(idx) == len(hei10):
                idx = []
            all_foci.append(pos[idx]/L)
            foci_intensities.append(hei10[tp][idx]/np.sum(hei10[tp][idx]))
            orig_foci_intensities.append(hei10[tp][idx])
            ns += 1
            L_array.append(L)



    return L_array, foci_intensities, orig_foci_intensities, all_foci


def to_val(v):
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return bool(v)

### What is this doing?

        
def load_data(survey_dir, survey_base):
    
    files = os.listdir(survey_dir)

    print('file list', files)
    
    all_data = []

    for f in files:
        if len(f)>=len(survey_base) and f[:len(survey_base)] == survey_base:
            all_data.append(read_file(survey_dir+'/'+f))
    
    
    return all_data


import pandas as pd


def make_data_files(data_path, output_prefix, only_wt_density=False):

    # Load preprocessed pkl file
    new_data = load_data(*data_path) 
    print('loaded data')




    
    index_df = pd.DataFrame([h for h,v in new_data])

    index_df = index_df.sort_values(['density', 'L'])

    grouped = index_df.groupby(['density'])
    
    print(list(grouped))


    


    for group_idx, (name, group) in enumerate(grouped):

        #print('density', density_values[j])

        print('name', name)


        print('group')
        print(group)

        

        dataset = []

        for i in group.index:
            s = group.loc[i]

            
            
            h, v = new_data[s.name]

        
            print('h, ', h)

        
        
            print(index_df.loc[s.name])

        
            L = to_val(h['L'])
            density = to_val(h['density'])
            if density != 0.5 and only_wt_density:
                continue
            K = to_val(h['K'])
            
            print('L, density, K', L, density, K)

            print(s)

            timepoint = -1

            for L_sc, pos, hei10 in v:
                data_row = dict(s)
                del data_row['filename_base']

                data_row['Lm'] = data_row['L']
                del data_row['L']
            

                hei10 = hei10[timepoint]
            
                idx = pc(hei10)
                if np.sum(idx)==len(pos):
                    idx = []

#            print(np.sum(idx))

                data_row['Number COs'] = np.sum(idx)
            
                data_row['L'] = L_sc
                pos_CO = pos[idx]
                intensity_CO = hei10[idx]

                for j, _ in enumerate(pos_CO):
                    data_row[f'position_{j}'] = pos_CO[j]
                    data_row[f'intensity_{j}'] = intensity_CO[j]
            
                dataset.append(data_row)


        dataset = pd.DataFrame(dataset)

        if not only_wt_density:
            print('Zero crossover SCs', np.sum(dataset['Number COs']==0))

        dataset.to_csv(output_prefix+'density_{}.csv'.format(name), index=False)

    # Now also output files for other densities.

    

    

def main(sim_data_path, output_path):

    #make_plots((sim_data_path+'/survey_escape/', 'at_'),  output_path+'/escape_',  intensity_bins=np.linspace(0,0.3,30))

    #return
    make_data_files((sim_data_path+'/survey_julia_new_ends/', 'at_'), output_path+'/fig2_simulations_')
    make_data_files((sim_data_path+'/survey_julia_ox/', 'at_'),  output_path+'/fig3_simulations_')
    make_data_files((sim_data_path+'/survey_julia_ux/', 'at_'),  output_path+'/fig4_simulations_')
    make_data_files((sim_data_path+'/survey_julia_no_ends/', 'at_'), output_path+'/figS3_simulations_')

    make_data_files((sim_data_path+'/survey_female/', 'at_'),  output_path+'/figS4b_simulations_', only_wt_density=True)

    make_data_files((sim_data_path+'/survey_escape/', 'at_'), output_path+'/figS6_simulations_')
    
main('../output', '../source_data')
                
        
