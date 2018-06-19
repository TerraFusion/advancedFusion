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
from skimage.measure import compare_ssim as ssim

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


def print_misr(file, camera, band, dset_path, out_path):
    """
    From an advanced fusion file, print the specified MISR camera
    and band into an output image.
    file - Path to Advanced Fusion HDF5 file
    camera (int) -- Chooses the first dimension of MISR_Radiance
    band (int) -- Chooses the second dimension of MISR_Radiance
    dset_path (str) -- Dataset Path within the HDF5 file where MISR radiance
        data is stored.
    out_path - output image path. Must contain a valid image extension (e.g. '.png')
    """

    f = h5py.File( file, "r" )
    dataset = f[dset_path]
    _FillValue = dataset.attrs['_FillValue']

    with timeit('test1'):
        # Select specific camera, band
        dataset = dataset[ camera, band ]
   
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

def print_aster(file, band, dest_path, out_path):
    """
    From an advanced fusion file, print the specified ASTER band to
    an output image.
    file (str) -- Path to an AF HDF5 file
    band (int) -- ASTER band to select
    dset_path (str) -- Absolute path of the ASTER dataset in the HDF5 file.
        Must give in format "/path/to/ASTER_Radiance"
    out_path (str) -- Output image path. Must contain a valid image
        extension.
    """

    f = h5py.File(file, "r")
    dataset = f[dset_path]
    
    _FillValue = dataset.attrs('_FillValue')

    dataset = dataset[band]
    # Clip fill values to 0
    dataset[ dataset == _FillValue ] = 0.0
    # Adjust range
    dataset *= (255.0 / dataset.max())
    dataset = np.uint8(dataset)
    img = Image.fromarray( dataset, "L")
    img.save(out_path)

def structural_similarity( img_path1, img_path2 ):

    image1 = np.asarray( Image.open( img_path1 ) )
    image2 = np.asarray( Image.open( img_path2 ) )

    return ssim( image1, image2 )


if __name__ == "__main__":
    #print_h5_dset(dset, None)
    #print_misr( path, 'AN', 'blue', './pictures/misr_modisonmisr.png')
    #print_modis(path, 0, "/Source/Data_Fields/MODIS_Radiance", "./pictures/modis_modisonmisr_0.png")
    #image1 = np.asarray( Image.open( "./pictures/modis_modisonmisr_0.png" ) )
    image1 = np.asarray( Image.open('./pictures/cosmogram136.jpg'))
    image2 = np.asarray( Image.open("./pictures/misr_modisonmisr.png") )
    print(ssim( image1, image2 ))
