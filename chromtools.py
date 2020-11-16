
# Code to extract ordered SC paths from skeletonized images, and to measure
# fluorescence intensities (HEI10 / DAPI) in small approximately-spherical
# regions along this path
# Note that Python 3.7 (not 3.8+) is necessary to reproduce the results in the MS

import logging

import numpy as np

from imageio import volread

# JIC utilities for image loading - essentially just numpy arrays here
# Useful in jupyter notebooks
from dtoolbioimage import Image, Image3D, scale_to_uint8

# Load tif as image
def tif_to_image3d(fpath):
    vol = volread(fpath)
    return np.transpose(vol, (1, 2, 0)).view(Image3D)

# Extract list of point coordinates from skeletonized image
def skel_to_points(skel):
    coords_list = list(zip(*np.where(skel)))
    return coords_list

def visualise_ordering(points_list, dim=(568, 568, 87)):
    rdim, cdim, _ = dim
    vis = np.zeros((rdim, cdim, 3), dtype=np.uint8)

    def get_col(i):
        r = int(255 * i/len(points_list))
        g = 255 - r
        return r, g, 0

    for n, p in enumerate(points_list):
        r, c, _ = p
        wr, wc = 5, 5
        vis[r-wr:r+wr,c-wc:c+wc] = get_col(n)

    return vis.view(Image)


# Walk along path until one end is found
def find_an_end(coords_list, p_start=None):
    available_points = set(coords_list)
    if p_start is None:
        p_start = coords_list[0]
    p = p_start

    def closest_point(p):
        v = np.array(list(available_points)) - p
        sq_dists = np.sum(v * v, axis=1)
        closest_index = np.argmin(sq_dists)
        sq_d = sq_dists[closest_index]
        return list(available_points)[closest_index], sq_d

    available_points.remove(p)

    ordered_points = []
    while len(available_points):
        ordered_points.append(p)
        p_next, sq_d = closest_point(p)
        # TODO - concealing messy stuff here
        if(sq_d) > 10:
            return p
        p = p_next

        available_points.remove(p)


def closest_point(available_points, p):
    v = np.array(list(available_points)) - p
    sq_dists = np.sum(v * v, axis=1)
    closest_index = np.argmin(sq_dists)
    sq_d = sq_dists[closest_index]
    return list(available_points)[closest_index], sq_d


# Trace along a list of points, adding the nearest point at each step
def order_points(coords_list, p_start=None):
    available_points = set(coords_list)
    if p_start is None:
        p_start = coords_list[0]
    p = p_start

    available_points.remove(p)

    ordered_points = []
    while len(available_points):
        ordered_points.append(p)
        p_next, sq_d = closest_point(available_points, p)
        if sq_d > 12:
            return ordered_points
        p = p_next
        available_points.remove(p)
    ordered_points.append(p)

    return ordered_points


# Sum of measure_stack over regin where mask==1
def measure_from_mask(mask, measure_stack):
    return np.sum(mask * measure_stack)

# Max of measure_stack over region where mask==1
def max_from_mask(mask, measure_stack):
    return np.max(mask * measure_stack)


# Translate mask to point p, treating makss near stack edges correctly
def make_mask_s(p, melem, measure_stack):
    mask = melem
    
    R = melem.shape[0] // 2
    r, c, z = p

    m_data = np.zeros(melem.shape)
    s = measure_stack.shape
    o_1, o_2, o_3 = max(R-r, 0), max(R-c, 0), max(R-z,0)
    e_1, e_2, e_3 = min(R-r+s[0], 2*R), min(R-c+s[1], 2*R), min(R-z+s[2], 2*R)
    m_data[o_1:e_1,o_2:e_2,o_3:e_3] = measure_stack[max(r-R,0):min(r+R,s[0]),max(c-R,0):min(c+R,s[1]),max(z-R,0):min(z+R, s[2])]
    return mask, m_data

# Measure the (mean/max) value of measure_stack about the point p, using
# the structuring element melem. op indicates the appropriate measurement (mean/max)
def measure_at_point(p, melem, measure_stack, op='mean'):
    if op=='mean':
        mask, m_data = make_mask_s(p, melem, measure_stack)
        melem_size = np.sum(melem)
        return float(measure_from_mask(mask, m_data) / melem_size)
    else:
        mask, m_data = make_mask_s(p, melem, measure_stack)
        return float(max_from_mask(mask, m_data))

# Generate spherical region
def make_sphere(R=5):
    x, y, z = np.ogrid[-R:R, -R:R, -R:R]
    sphere = x**2 + y**2 + (2.3 * z)**2 < R**2
    return sphere

# Measure the values of measure_stack at each of the points of points_list in turn.
# Measurement is the mean / max (specified by op) on the spherical region about each point
def measure_all_with_sphere(points_list, measure_stack, op='mean'):
    melem = make_sphere()
    measure_func = lambda p: measure_at_point(p, melem, measure_stack, op)
    return list(map(measure_func, points_list))


# Measure fluorescence levels along ordered skeleton
def measure_chrom2(single_chrom):
    # single chrom - structure containing skeleton (single_chrom.skel) and
    # fluorecence levels (single_chrom.hei10) as Image3D objects (equivalent to ndarray)
    # Returns list of coordinates in skeleton, the ordered path 
    coords_list = skel_to_points(single_chrom.skel)
    vis_unordered = visualise_ordering(coords_list, dim=single_chrom.skel.shape)
    p_end = find_an_end(coords_list)
    print('start_order_points')
    ordered_points = order_points(coords_list, p_end)
    print('end_order_points', len(ordered_points))

    measurements = measure_all_with_sphere(ordered_points, single_chrom.hei10, op='mean')
    measurements_max = measure_all_with_sphere(ordered_points, single_chrom.hei10, op='max')

    return coords_list, ordered_points, measurements, measurements_max

