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

from measure import measure

# Script to preprocess image data, generating .pkl files containing HEI10 data and the paths for each SC
# Edit data_path to point to downloaded and decompressed image data files

##### Edit this path to point to downloaded dataset
data_path = '/media/foz/a92bd7ce-a444-4635-8659-b03fda836a5e/JIC/'

os.makedirs('data_output', exist_ok=True)


# Combine multiple dictionaries
def join_dicts(list_of_dicts):
     return reduce(lambda a, b: {**a, **b}, list_of_dicts)

def thin_points(point_list, dmin=10):
    # point_list [ (x[3], intensity, idx ) ]
    # Remove points within dmin of each other, keeping point with highest intensity (so order of removal doesn't matter)
    # return index (third element of point tuples)

    # point list - list of location, score, key
    removed_points = []
    for i in range(len(point_list)):
        if point_list[i][2] in removed_points:
            continue
        for j in range(len(point_list)):
            if i==j:
                continue
            if point_list[j][2] in removed_points:
                continue
            d = np.array(point_list[i][0]) - np.array(point_list[j][0])
            d = la.norm(d)
            if d<dmin:
 #               print(i, j, point_list[i][0], point_list[j][0], d)
                hi = point_list[i][1]
                hj = point_list[j][1]
                if hi<hj:
                    removed_points.append(point_list[i][2])
                    break
                else:
                    removed_points.append(point_list[j][2])
                
    return removed_points



def get_path_and_image(base_path, date, genotype, plant, slide, image, path, marker='HEI10', category='wt'):

    # Find the path etc within the directory structure. There are some inconsistencies in naming within some datasets that require exceptions

    if category=='wt':
        #slide_dir = date + ' ' + 'plant ' + str(plant) + ' slide ' + str(slide)
        if 'hei10' in genotype:
            date_dir = date + ' (hei10 ++)'
        else:
            date_dir = date + ' (Col-0)'
        if os.path.exists(base_path+date_dir): # First data set
            pl = os.listdir(base_path+date_dir)
            plant_dir = None
            for p in pl:
                #print('pl entry ', p)
                if int(p.split(' ')[1])==plant:
                    plant_dir = p
            if plant_dir == None:
                print('problem finding plant data', date_dir, pl)
            sl = os.listdir(base_path+date_dir + '/' + plant_dir)
            slide_dir = 'slide '+str(slide)
            if slide_dir in sl:
                d = os.listdir(base_path+date_dir + '/' + plant_dir + '/' + slide_dir) 
                if 'skeletonized images' in d:
                    image_dir = date_dir + '/' + plant_dir + '/' + slide_dir + '/skeletonized images/image'+str(image)+'/'
                elif 'skeletonized' in d:
                    image_dir = date_dir + '/' + plant_dir + '/' + slide_dir + '/skeletonized/image'+str(image)+'/'
                else:
                    image_dir = date_dir + '/' + plant_dir + '/' + slide_dir + '/skeletonised/image'+str(image)+'/'
            else:
                slide_dir = genotype + ' ' + date + ' ' + 'plant ' + str(plant) + ' slide ' + str(slide)
                image_dir = slide_dir + '/image'+ str(image)+'/'
        else:
            slide_dir = genotype + ' ' + date + ' ' + 'plant ' + str(plant) + ' slide ' + str(slide)
            image_dir = slide_dir + '/image'+ str(image)+'/'
        print(image_dir)
    else:
        slide_dir = date + ' ' + 'plant ' + str(plant) + ' slide ' + str(slide)
        if ( (category=='ox' and slide_dir=='2.8.19 plant 22 slide 2') or
           (category=='ux' and slide_dir=='2.8.19 plant 22 slide 1') ):
            image_dir = slide_dir + '/image '+ str(image)+'/'
        else:
            image_dir = slide_dir + '/image'+ str(image)+'/'


    print (base_path+image_dir+'Path'+str(path)+'.tif', base_path+image_dir+marker+'.tif')    
    return (base_path+image_dir+'Path'+str(path)+'.tif', base_path+image_dir+marker+'.tif')    


def get_traces4(base_path, date, genotype, plant, slide, image, path, marker='HEI10', traces=False, category='wt'):
    # Extract data from skeleton and image
    # Returns skeleton points, ordered list of skeleton points, mean hei10 in a spherical region about each
    #         ordered point, max hei10 in the same spherical regions
    r = measure(*get_path_and_image(base_path, date, genotype, plant, slide, image, path, marker, category=category))
    #print(ordered_points, hei10)
    return r


# Use ray package to permit multi-threaded execution
import ray
ray.init()


@ray.remote
def process_traces(base_path, A, n, category):
    # Main function preprocessing skeleton trace data
    # base_path - directory containing data for this experiment (WT / OX / UX)
    # A - pandas dataframe containing data from appropriate .csv file
    # n - cell number to process
    # category ('wt' / 'ox' / 'ux' - used to demangle filenames)
    # Returns peaks (with position and normalized intensities), length, original (non-normalized) hei10 values,
    # and the ordered point positions along each SC
     
    i = n*5
    all_peak_data = {}
    orig_trace_lengths = {}
    o_hei10_traces = {}
    o_dapi_traces = {}
    o_points = {}
    
    cell_peaks = []

    # Process each cell (groups of five chromosomes) together
    for j in range(i, i+5):

        # Extract trace data for each SC
        img = int(A.iloc[j].Image)
        img_ext = int(round((A.iloc[j].Image - img)*10)*10)
        op, points, o_hei10, o_hei10_max = get_traces4(base_path,A.iloc[j].Date,A.iloc[j].Genotype,int(A.iloc[j].Plant),int(A.iloc[j].Slide),img,img_ext+1+(j%5), traces=True, category=category)
        o_hei10_traces[j] = o_hei10
        o_points[j] = points
        
        if category =='wt':
             # For WT also extract DAPI along each SC
             op, points, o_dapi, o_dapi_max = get_traces4(base_path,A.iloc[j].Date,A.iloc[j].Genotype,int(A.iloc[j].Plant),int(A.iloc[j].Slide),img,img_ext+1+(j%5), traces=True, marker='DAPI', category=category)
        
             o_dapi_traces[j] = o_dapi


        # For peak determination normalize each trace to have mean zero and s.d. 1
        hei10 = (o_hei10 - np.mean(o_hei10))/np.std(o_hei10)
        hei10_max = (o_hei10_max - np.mean(o_hei10_max))/np.std(o_hei10_max)
        
        orig_trace_lengths[j] = len(op)


        v = hei10
        # Find peaks - these will be further refined later
        p,_ = find_peaks(v, distance=5,  prominence=0.5*np.std(v))
     
        
        peaks = np.array(p, dtype=np.int32)

        # Store peak data - using original values, not normalized ones
        peak_mean_heights = [ o_hei10[u] for u in peaks ]
        peak_max_heights = [ o_hei10_max[u] for u in peaks ]
        peak_points = [ points[u] for u in peaks ]
        
        path_length = len(points)
        
        cell_peaks.append((peaks, peak_points, peak_mean_heights, peak_max_heights, path_length))


        
    # Eliminate peaks which have another larger peak nearby (in 3D space, on any chromosome).
    # This aims to remove small peaks in the mean intensity generated when an  SC passes close
    # to a bright peak on another SC - this is nearby in space, but brighter.
    
    to_thin = []
    for k in range(5):
        for u in range(len(cell_peaks[k][0])):
            to_thin.append((cell_peaks[k][1][u], cell_peaks[k][2][u], (k, u)))
    
    # Exclude any peak with a nearby brighter peak (on any SC)
    removed_points = thin_points(to_thin)

    # Clean up and remove these peaks
    new_cell_peaks = []
    for k in range(len(cell_peaks)):
        cc = []
        pp = list(zip(*cell_peaks[k][:4]))
        for u in range(len(pp)):
            if (k,u) not in removed_points:
                cc.append(pp[u])
        ccp = list(zip(*cc))
        if not ccp:
            ccp = [[],[],[],[]]
        new_cell_peaks.append(ccp+[cell_peaks[k][4]])
        
    cell_peaks = new_cell_peaks

    # Save peak positions, normalized HEI10 intensities, and length for each SC
    for k in range(5):
        j = i + k
        peaks = np.array(cell_peaks[k][0], dtype=int)
        hei10 = o_hei10_traces[j]
        hei10 = hei10 - np.mean(hei10)
        hei10 = hei10 / np.std(hei10)
        #print(j, peaks, hei10[peaks], len(hei10))
        all_peak_data[j] = (peaks, hei10[peaks], len(hei10))

    # Also return original trace lengths, hei10 levels, and the ordered list of points along the SC skeleton
    if category=='wt':    
         return all_peak_data, orig_trace_lengths, o_hei10_traces, o_dapi_traces, o_points
    else:
         return all_peak_data, orig_trace_lengths, o_hei10_traces, o_points
    


def process_data(base_dir, csv_fn, output_fn, category):

    # Read CSV file with data specification
    A = pd.read_csv(csv_fn)
    # Remove some unused columns and rename one field
    A = A.drop(A.columns[[6,7,8,9]], axis=1)
    A = A.rename(columns={'Plant ':'Plant'})

    # Parallel execution of HEI10 trace analysis
    results = [process_traces.remote(base_dir, A, i, category) for i in range(len(A)//5)]
    results = ray.get(results)

    # DAPI data not extracted from original stacks for OX and UX data
    if category=='wt':
        all_peak_data, orig_trace_lengths, o_hei10_traces, o_dapi_traces, o_points = (map(join_dicts, zip(*results)))
        data = { 'all_peak_data': all_peak_data, 'orig_trace_lengths': orig_trace_lengths, 'o_hei10_traces':o_hei10_traces, 'o_dapi_traces':o_dapi_traces, 'o_points':o_points }
    else:
        all_peak_data, orig_trace_lengths, o_hei10_traces, o_points = (map(join_dicts, zip(*results)))
        data = { 'all_peak_data': all_peak_data, 'orig_trace_lengths': orig_trace_lengths, 'o_hei10_traces':o_hei10_traces, 'o_points': o_points  }

    with open(output_fn, 'wb') as f:
        pickle.dump(data, f, protocol=3)


data_output_path = 'data_output/'

process_data(data_path +'WT Col-0/', '200406.csv', data_output_path+'test.pkl', 'wt')
process_data(data_path+'ox/HEI10 overexpressor/', 'OX.csv', data_output_path+'test_ox.pkl', 'ox')
process_data(data_path+'underexpressor/HEI10 underexpressor/', 'UX.csv', data_output_path+'test_ux.pkl', 'ux')
