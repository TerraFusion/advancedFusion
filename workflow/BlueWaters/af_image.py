"""
This module provides a suite of functions to create images 
(.png, .jpg, etc) from datasets in an Advanced Fusion HDF5
file.
"""

import h5py
import numpy as np
from PIL import Image
from enum import IntEnum
import contextlib
import time

@contextlib.contextmanager
def timeit(ident):
    #tstart = time.time()
    yield
    #elapsed = time.time() - tstart
    #print("{}: {} s".format(ident, elapsed))

path='/u/sciteam/clipp/scratch/af_out/data/2007.01/ADVNCE_FUSE_TERRA_BF_L1B_O37661_20070116120124_F000_V001.h5'
#path="/u/sciteam/clipp/scratch/af_out/data/2007.01/ADVNCE_FUSE_TERRA_BF_L1B_O37660_20070116102231_F000_V001.h5"
#dset="Source/Data_Fields/MODIS_Radiance"
dset="Target/Data_Fields/MISR_Radiance"

misr_cam = { 
    'aa' : 0,
    'af' : 1,
    'an' : 2,
    'ba' : 3,
    'bf' : 4,
    'ca' : 5,
    'cf' : 6,
    'da' : 7,
    'df' : 8,
}

misr_band = {
    'red'   : 0,
    'green' : 1,
    'blue'  : 2,
    'nir'   : 3
}

def print_misr(file, camera, band, out_path):
    """
    From an advanced fusion file, print the specified MISR camera
    and band into an output image.
    file - Path to Advanced Fusion HDF5 file
    camera - Choose from AA, AF, AN, BA, BF, CA, CF, DA, DF
    band - Choose from red, green, blue, nir.
    out_path - output image path. Must contain a valid image extension (e.g. '.png')
    """

    camera = camera.lower()
    band = band.lower()
    if camera not in misr_cam:
        raise ValueError("choose camera from {}".format( misr_cam ))
    if band not in misr_band:
        raise ValueError("Choose band from {}".format( misr_band))
    
    dset_path="/Target/Data_Fields/MISR_Radiance"
    f = h5py.File( file, "r" )
    dataset = f[dset_path]
    _FillValue = dataset.attrs['_FillValue']

    with timeit('test1'):
        # Select specific camera, band
        dataset = dataset[ misr_cam[camera], misr_band[band] ]
   
    with timeit('test2'):
        # Clip fill values to 0
        dataset[ dataset == _FillValue ] = 0.0

    with timeit('test3'):
        #  Adjust range
        dataset *= (255.0 / dataset.max())

    with timeit('test4'):
        img = Image.fromarray( np.uint8(dataset), "L" )
        img.save( out_path )

def print_modis(file, band, dset_path, out_path):
    """
    From an advanced fusion file, print the specified MODIS band
    to an output image.
    file (str)  -- Path to Advanced Fusion HDF5 file
    band (int) -- MODIS band to select
    dset_path (str) -- Absolute path of the MODIS dataset in the HDF5 file.
        Must give in format: "/Path/To/MODIS_Radiance"
    out_path (str) -- Output image path. Must contain valid image extension,
        e.g. '.png'
    """

    f = h5py.File(file, "r")
    dataset = f[dset_path]

    _FillValue = dataset.attrs['_FillValue']
    valid_min = dataset.attrs['valid_min']

    with timeit('test1'):
        dataset = dataset[ band ]

    with timeit('test2'):
        # Clip invalid values to 0
        dataset[ np.logical_or( dataset < valid_min, dataset == _FillValue) ] = 0.0

    with timeit('test3'):
        # Adjust range
        dataset *= (255.0 / dataset.max() )

    with timeit('testconvert'):
        dataset = np.uint8(dataset)

    with timeit('test4'):
        img = Image.fromarray( dataset, "L" )
        img.save(out_path)

if __name__ == "__main__":
    #print_h5_dset(dset, None)
    #print_misr( path, 'AN', 'blue', './pictures/misr_modisonmisr.png')
    print_modis(path, 0, "/Source/Data_Fields/MODIS_Radiance", "./pictures/modis_modisonmisr_0.png")
