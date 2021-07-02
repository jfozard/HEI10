from dataclasses import dataclass
from pathlib import Path

import click

import numpy as np
import matplotlib.pyplot as plt

from dtoolbioimage import Image, Image3D

from chromtools import (
    tif_to_image3d,
    measure_chrom2,
)



# Simple structure to contain both the skeletonized image
# and the fluorescence stack to measure
@dataclass
class SingleChromImage:
    name: str
    hei10: Image3D
    skel: Image3D

    @classmethod
    def from_ds_triplet(cls, ds, triplet):
        name, skel_idn, hei10_idn = triplet
        skel = tif_to_image3d(ds.item_content_abspath(skel_idn))
        hei10 = tif_to_image3d(ds.item_content_abspath(hei10_idn))
        return cls(name, hei10, skel)

    @classmethod
    def from_fpaths(cls, path_image_fpath, hei10_image_fpath):
        skel = tif_to_image3d(path_image_fpath)
        hei10 = tif_to_image3d(hei10_image_fpath)
        name = "path"
        return cls(name, hei10, skel)


# Measure fluorescence intensity along SC skeleton
def measure(path_image_fpath, hei10_image_fpath):
    # path_image_fpath - path to skeletonized tif stack
    # hei10_image_fpath -
    # returns - points along path, points ordered along SC, mean intensity measurements about each point
    #           maximum intensity measurement in spheres about each point
    
    single_chrom = SingleChromImage.from_fpaths(path_image_fpath, hei10_image_fpath)
    orig_points, ordered_points, measurements, measurements_max = measure_chrom2(single_chrom)

    return orig_points, ordered_points, measurements, measurements_max


if __name__ == "__main__":
    main()
