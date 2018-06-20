import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import argparse
import logging
import os
import errno
import datetime
import bfutils
import sys
from mpi4py import MPI
import re
import time
from schwimmbad import MPIPool
import yaml
import subprocess
from enum import IntEnum
import af_image
from PIL import Image
import json

#============================================================================
# The following disables yaml's conversion of certain values into
# Python bools. I was having issues with it converting the string "ON"
# to True. The string "True" would then be passed to the AF tool, which 
# is invalid. This prevents that behavior.
def add_bool( self, node ):
     return self.construct_scalar(node)
yaml.constructor.Constructor.add_constructor(u'tag:yaml.org,2002:bool', add_bool)
#
#===========================================================================



# Define communication tags
Tags = IntEnum('Tags', 'START SUCCEED_AF SUCCEED_H5DUMP SUCCEED_IMG SUCCEED_ALL' )

LOG_FMT='[%(asctime)s] [%(name)12s] [%(levelname)8s] [%(filename)s:%(lineno)d] %(message)s'
LOG_DATE_FMT='%d-%m-%Y:%H:%M:%S8'

class Granule:
    """
    Class that holds all information necessary for a specific orbit of BF->AF
    data.
    """
    
    wf_params = {"COMPARE_IMAGES", "COMPARE_THRESHOLD","OUTPUT_PREFIX",
    "PRINT_ALL_IMG"}
    COMPARE_THRESHOLD_DEFAULT = 0.9
    #tags = IntEnum('tags', 'START SUCCEED FAIL_IMAGE FAIL_AF FAIL_H5DUMP FAIL')

    def __init__(self ):
        self.af_config = {}
        self.wf_config = {}
        
        # Structural similarity of the images
        self.similarity = []

        self.images = {}
        self.status = Tags.START

    def sim_above_thresh( self ):
        """
        Return True if the similarity is above the designated threshold.
        Returns False if any of the image similarities in self.similarity
        is below the COMPARE_THRESHOLD, or if the similarity list is empty.
        """

        above_thresh = True
        for sim in self.similarity:
            if sim < self.get_workflow_config()["COMPARE_THRESHOLD"]:
                above_thresh = False

        if not self.similarity and above_thresh:
            return True
        
        return False

    def set_status( self, status ):
        """
        Set the processing status. The proper tags are defined in Tags
        """
        self.status = status

    def get_status( self ):
        return self.status

    def set_input_path(self, path):
        self.af_config[in_file_key] = path
   
    def get_input_path(self):
        """
        Get the basic fusion input HDF5 file path.
        """
        return self.af_config[in_file_key]

    def set_output_path(self, path):
        """
        Set the output path of the advanced fusion file.
        """
        self.af_config[out_file_key] = path
   
    def get_output_path(self):
        """
        Get the output path of the advanced fusion file.
        """
        return self.af_config[out_file_key]

    def set_af_config_file(self, file ):
        """
        Set the path of the config file
        """
        self.af_config_file = file

    def get_af_config_file( self ):
        """
        Return path to the job config file
        """
        return self.af_config_file

    def set_config(self, config):
        """
        Stores a dict that will be placed into self.af_config_file
        (to be passed directly to the advanced fusion tool).
        For the parameters not used by the AF tool (but by this script),
        they are separated into a separate dict. 

        To retrieve parameters for the AF tool, call: get_af_config()

        To retrieve all other parameters, call: get_workflow_config() 
        """
        
        for key in config:
            if key.upper() != key:
                raise ValueError("All keys in the configuration file should \
                be upper-case!")

            if key.upper() in self.wf_params:
                val = config[key]
                if type(val) is str:
                    val = val.upper()

                self.wf_config[key.upper()] = val

            else:
                # Force all values in af_config to be strings.
                # This is done simply for consistency.
                self.af_config[key] = str(config[key])

        # Set default value for compare_threshold
        if "COMPARE_THRESHOLD" not in self.wf_params:
            self.wf_params["COMPARE_THRESHOLD"] = self.COMPARE_THRESHOLD_DEFAULT

    def get_workflow_config(self):
        """
        Return the workflow-specific parameters passed to this script.
        Returns a dict.
        """
        return self.wf_config

    def get_af_config(self):
        """
        Returns a dict that represents the configuration file passed
        to the AF tool.
        """
        return self.af_config

    def set_log_file(self, path):
        """
        Set the path of the file that will contain log output for 
        this granule.
        """
        self.log_file = path

    def get_log_file(self):
        return self.log_file

    def set_exe(self, exe):
        """
        Set path to the executable to run.
        """
        self.exe = os.path.abspath(exe)

    def get_exe(self):
        """
        Get the binary exectuable path for Advanced Fusion
        """
        return self.exe

    def set_orbit(self, orbit):
        self.orbit = orbit

    def get_orbit(self):
        """
        Return this granule's orbit (int)
        """
        return self.orbit

    def set_h5dump(self, h5dump):
        """
        Set path of h5dump file.
        """
        self.h5dump = h5dump

    def get_h5dump(self):
        """
        Return path of h5dump file
        """
        return self.h5dump

    def set_image_dir(self, dir):
        self.image_dir = dir

    def get_image_dir(self):
        return self.image_dir


    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

    
class Log(logging.Logger):
    """
    Class that contains logging utilites and metadata.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger = None
        self.run_dir = None
        self.level = None
        self.log_dirs = {}
        self.orbit = None

    def set_log_file( self, log ):
        self.log_file = log
    def get_log_file( self ):
        return self.log_file

    def set_level(self, level):
        self.level = level

    def addFileHandler(self, path, *args, **kwargs):
        """
        Wraps around logging.FileHandler and self.addHandler to add a
        file handler to the logger. Also stores the path of the log
        file to self (can be retrieved by get_log_file() at a later
        time. (path, *args, **kwargs) are the arguments passed directly
        to logging.FileHandler().
        """

        fileHandler = logging.FileHandler( path, *args, **kwargs )
        fileHandler.setFormatter(logFormatter)
        self.addHandler(fileHandler)
        self.set_log_file(path)

#=================== MPI VARIABLES =====================
mpi_comm = MPI.COMM_WORLD
mpi_rank = MPI.COMM_WORLD.Get_rank()
mpi_size = MPI.COMM_WORLD.Get_size()
rank_name = "Rank {}".format(mpi_rank)
pool = None

#================== LOGGING ======================
logging.setLoggerClass(Log)
logFormatter = logging.Formatter( LOG_FMT, LOG_DATE_FMT)
logger = logging.getLogger(name = rank_name)

#================== CONSTANTS =========================
FAILED_FILE = 'failed_granules.txt'
GRANULE_SUMMARY = 'granule_summary.txt'
out_file_key = 'OUTPUT_FILE_PATH'
in_file_key = 'INPUT_FILE_PATH'
SRC_DSET_PATH = "Source/Data_Fields/{}_Radiance"
TARGET_DSET_PATH = "Target/Data_Fields/{}_Radiance"
NUM_MODIS_BANDS = 38
IMG_EXTENSION = ".png"
LOG_SEP = "===================================================="

def make_run_dir( dir ):
    '''
    This function makes a directory called run# where # is the first integer greater than any other run# in
    the directory given by the argument dir.

    Returns the path of the newly created directory.
    '''

    if not isinstance(dir, str):
        raise TypeError("Passed argument \'dir\' is not a string!")

    greatestSeen=-1
    runRE='^run[0-9]+$'

    # Walk through every entry in dir. Find the greatest run#.
    for i in os.listdir( dir ):
        if not re.match( runRE, i ) :
            continue
        tempRunNum = int( i.replace( 'run', '') )
        if tempRunNum > greatestSeen:
            greatestSeen = tempRunNum

    newNum = greatestSeen + 1

    # We can make our new run directory
    newDir = os.path.abspath(os.path.join( dir, 'run' + str(newNum) ))
    os.mkdir( newDir )

    return newDir

def do_h5dump( granule ):
    """
    Generate h5dump of an hdf5 file
    """
    # Create the yyyy.mm dir in the h5dump directory
    try:
        os.makedirs( os.path.dirname( granule.get_h5dump() ))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    args = ['h5dump', '-H', granule.get_af_config()[out_file_key] ]
    
    with open( granule.get_h5dump(), 'w' ) as f:
        subprocess.check_call( args, stdout=f, stderr=subprocess.STDOUT)

def do_af( granule):
    """
    Generate advanced fusion output for a granule.
    """
    # Now we can call the executable
    args = [ granule.get_exe(), granule.get_af_config_file() ]
    logger.debug(' '.join(args))
    with open( granule.get_log_file(), 'a' ) as f:
        subprocess.check_call( args, stdout=f, 
            stderr=subprocess.STDOUT)
   
def class_attr_string( obj ):
    """
    Return a string that represents all of a class's attributes.

    The object must have an __iter__ method that yields attributes
    and their values.
    """
    data_str = ""
    for attr, val in obj:
        if type(val) is dict:
            # Use json library to print a pretty version of dict
            val = json.dumps(val, sort_keys=True, indent=2)

        data_str = data_str + "{}: {}\n".format(attr, val)

    return data_str

def do_img_gen_verify( data ):
    wf_config = data.get_workflow_config()
    af_config = data.get_af_config()

    if wf_config["COMPARE_IMAGES"] or wf_config["PRINT_ALL_IMG"]:
        # We need to figure out which datasets to retrieve based on
        # data's configuration.
        
        logger.info("Creating images from datasets...")

        src_dset_path = SRC_DSET_PATH.format( af_config["SOURCE_INSTRUMENT"].upper() )
        target_dset_path = TARGET_DSET_PATH.format( af_config["TARGET_INSTRUMENT"].upper() )
        
        images = []
        for instr_path in src_dset_path, target_dset_path:
            # Unfortunately, these ugly if statements are necessary because
            # each instrument function takes different parameters :(
            

            if "MODIS_Radiance" in instr_path:
                bands = af_config["MODIS_BANDS"]
                resolutions = af_config["MODIS_RESOLUTION"].upper().split()

                logger.info("Creating MODIS images...")
                
                if bands.upper() == "ALL":
                    bands = range( NUM_MODIS_BANDS )
                else:
                    # Split each band to a separate element, then convert 
                    # each element to an integer
                    bands = [ int(i) for i in bands.split() ]

                # band_idx and band are two different things. idx is 
                # index into dataset. Bands are assigned indices
                # according to which order they were given in the
                # config file.

                for resolution_idx, resolution in enumerate( resolutions ):
                    for band_idx, band in enumerate(bands):

                        img_name = instr_path.replace('/', '#') + "_{}_band_{}".format(resolution, band) + \
                            "_orbit_{}".format( data.get_orbit() ) + IMG_EXTENSION

                        img_path = os.path.join( data.get_image_dir(), img_name )
                        
                        # Save image path into data
                        dict_walker = data.images
                        modis_dset = ["MODIS", resolution ]
                        for key in modis_dset:
                            if key not in dict_walker:
                                dict_walker[str(key)] = {}

                            dict_walker = dict_walker[str(key)]

                        logger.debug("data.images['MODIS'][{}][{}] = {}".format(
                            resolution, band, img_path))

                        data.images["MODIS"][ resolution ][ str(band) ] = img_path
                        
                        images.append(img_path)
                        
                        logger.debug(LOG_SEP)
                        logger.debug("out_path: {}".format( data.get_output_path()))
                        logger.debug("instr_path: {}".format(instr_path))
                        logger.debug("img_path: {}".format(img_path))
                        logger.debug("band: {}".format(band))
                        logger.debug("band_idx: {}".format(band_idx))

                        af_image.print_modis( data.get_output_path(), band_idx, instr_path, img_path )
          
                logger.debug(LOG_SEP)
                
            elif "MISR_Radiance" in instr_path:
                # Get MISR camera angles
                angles = af_config["MISR_CAMERA_ANGLE"].upper().split()
                bands = af_config["MISR_RADIANCE"].split()
                resolutions = af_config["MISR_RESOLUTION"].upper().split()

                logger.info("Creating MISR images.")

                # We use enumerate to get index, instead of using the actual string.
                # print_misr expects integer that serves as index into dataset. The
                # order in which angles/bands are listed in the config file is the
                # order in which they are stored in the dataset.
                for res_idx, res in enumerate( resolutions ):
                    for angle_idx, angle in enumerate( angles ):
                        for band_idx, band in enumerate( bands ):
                            img_name = instr_path.replace('/', '#') + \
                                "_{}_{}_{}{}".format( angle, band.replace("_Radiance", "") , data.get_orbit(), 
                                IMG_EXTENSION )

                            img_path = os.path.join( data.get_image_dir(), img_name )

                            # Save image path into data.
                            # Recursively create a chain of dicts that
                            # describe the dataset. For instance, we would
                            # like to be able to access it like this:
                            # data.images["MODIS"]["L"]["AN"]["Red_Radiance"]
                            # 
                            # and have the path stored at that location.

                            dict_walker = data.images
                            misr_dset = ["MISR", res, angle ]
                            for key in misr_dset:
                                if key not in dict_walker:
                                    dict_walker[str(key)] = {}

                                dict_walker = dict_walker[str(key)]
                                
                            data.images["MISR"][res][angle][str(band)] = img_path

                            logger.debug(LOG_SEP)
                            logger.debug("out_path: {}".format( data.get_output_path()))
                            logger.debug("angle: {}".format(angle))
                            logger.debug("angle_idx: {}".format( angle_idx))
                            logger.debug("band: {}".format(band))
                            logger.debug("band_idx: {}".format(band_idx))
                            logger.debug("instr_path: {}".format(instr_path))
                            logger.debug("img_path: {}".format( img_path))

                            af_image.print_misr( data.get_output_path(), angle_idx, 
                                band_idx, instr_path, img_path )

                
                    logger.debug(LOG_SEP)

            elif "ASTER_Radiance" in instr_path:
                # Get ASTER parameters
                res = af_config["ASTER_RESOLUTION"].upper().split()[0]
                bands = af_config["ASTER_BANDS"].upper().split()
                
                logger.info("Creating ASTER images.")
                # Landon Clipp 6/19/2018
                # Apparently, the ASTER_Radiance dataset only stores 
                # different bands in the dimensions. So just iterate
                # over the first dimension. Note: first for loop
                for band_idx, band in enumerate( bands ):
                    img_name = "{}_{}_{}_{}{}".format( instr_path.replace('/', '#'),
                    res, band, data.get_orbit(), IMG_EXTENSION )

                    img_path = os.path.join( data.get_image_dir(), img_name )

                    dict_walker = data.images
                    aster_dset = ["ASTER", res ]
                    for key in aster_dset:
                        if key not in dict_walker:
                            dict_walker[key] = {}
                        dict_walker = dict_walker[str(key)]

                    data.images["ASTER"][res][band] = img_path

                    logger.debug(LOG_SEP)
                    logger.debug("out_path: {}".format( data.get_output_path()))
                    logger.debug("band: {}".format(band))
                    logger.debug("band_idx: {}".format(band_idx))
                    logger.debug("instr_path: {}".format( instr_path ))
                    logger.debug("img_path: {}".format(img_path))

                    af_image.print_aster( data.get_output_path(), band_idx, 
                        instr_path, img_path )

                logger.debug(LOG_SEP)

            else:
                raise ValueError("Cannot determine the instrument \
                for image generation!")

        logger.info("Images generated.")

        if wf_config["COMPARE_IMAGES"]:
            logger.info("Creating similarity metrics.")

            for idx, val in enumerate( wf_config["COMPARE_IMAGES"] ) :
                # If it's a list, retrieve the cameras from each element
                if type(val) is list:
                    dset1 = val[0].split()
                    dset2 = val[1].split()

                # Else if it's a string, the next element is our second dataset
                elif type(val) is str:
                    if len(wf_config["COMPARE_IMAGES"]) > 2:
                        raise ValueError("COMPARE_IMAGES gave more than 2 \
                        images to compare. Please refer to the README.md for \
                        this script.")

                    if idx >= 1:
                        continue

                    dset1 = val.split()
                    dset2 = wf_config["COMPARE_IMAGES"][idx + 1].split()
                else:
                    raise ValueError("Unknown value in COMPARE_IMAGES.")

                # "Recursively" iterate through the keys to retrieve the paths of the images.
                # All I'm doing is looping through the elements in dset1 and dset2 and retrieving
                # the intermediate values from data.images. For instance, 
                # say dset1 = ['MODIS', '1KM', '4']:
                # image1_path = data.images
                # image1_path = image1_path['MODIS']
                # image1_path = image1_path['1KM']
                # image1_path = image1_path['4']

                logger.debug(data.images)

                image1_path = data.images
                for key in dset1:
                    image1_path = image1_path[str(key)]

                image2_path = data.images
                for key in dset2:
                    image2_path = image2_path[str(key)]

                sim = af_image.structural_similarity( image1_path, image2_path )
                data.similarity.append( sim )

                if sim < wf_config["COMPARE_THRESHOLD"]:
                    logger.warning("Images have similarity below designated threshold!")
                    logger.warning("image 1: {}".format( image1_path ) )
                    logger.warning("image 2: {}".format( image2_path ) )
                    logger.warning("similarity: {}".format( sim ) )
                    logger.warning("similarity threshold: {}".format( wf_config["COMPARE_THRESHOLD"] ) )

                else:
                    logger.info("Images above similarity threshold.")
                    logger.debug("image 1: {}".format( image1_path ) )
                    logger.debug("image 2: {}".format( image2_path ) )
                    logger.info("similarity: {}".format( sim ) )
                    logger.info("similarity threshold: {}".format( wf_config["COMPARE_THRESHOLD"] ) )
    
def worker( data ):

    try:
        global logger

        # Remove any old file handlers from previous orbits
        logger.handlers = [ h for h in logger.handlers if (type(h) != logging.FileHandler)]

        # Set the per-rank file handler for this specific orbit
        fileHandler = logging.FileHandler( data.get_log_file(), mode='a')
        fileHandler.setFormatter( logFormatter )
        logger.addHandler( fileHandler)

        data_str = class_attr_string( data )

        logger.info("Received data: \n{}".format( data_str ))
        
        logger.info("Creating config file: {}".format(data.get_af_config_file()))
        # First need to create the config file
        with open( data.get_af_config_file(), 'w' ) as f:
            af_config = data.get_af_config()
            for key in af_config:
                f.write("{}: {}\n".format( key, af_config[key] ) )

        # Create the yyyy.mm dir in the output file path
        try:
            os.makedirs( os.path.dirname( data.get_af_config()[out_file_key] ))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        if data.get_status() != Tags.START:
            logger.info("Retrying failed steps...")

        #=======================
        # advanced fusion
        #

        if data.get_status() == Tags.START:
            logger.info("Calling the AFtool")
            try:
                do_af(data)
            except:
                logger.exception("Encountered exception when processing AF \
        granule: {}.\nSee: {}\nfor more details.".format(
                    data.get_orbit(), data.get_log_file()))

                return data

            data.set_status( Tags.SUCCEED_AF )
        
        #===============
        # h5dump
        #
        if data.get_status() == Tags.SUCCEED_AF:
            logger.info("Performing h5dump on file...") 
            try:
                do_h5dump( data )
            except:
                logger.exception("Encountered exception when processing h5dump: {}.\
        \nSee: {}\nfor more details.".format(
                    data.get_orbit(), data.get_h5dump()))

                return data
        
            logger.info("Done.")
            data.set_status( Tags.SUCCEED_H5DUMP )

        #===========================
        # image generation/verification
        #
        if data.get_status() == Tags.SUCCEED_H5DUMP:
            try:
                do_img_gen_verify( data )

#                wf_config = data.get_workflow_config()
#          
#                if wf_config["COMPARE_IMAGES"] or wf_config["PRINT_ALL_IMG"]:
#                    # We need to figure out which datasets to retrieve based on
#                    # data's configuration.
#                    
#                    logger.info("Creating images from datasets...")
#
#                    src_dset_path = SRC_DSET_PATH.format( af_config["SOURCE_INSTRUMENT"].upper() )
#                    target_dset_path = TARGET_DSET_PATH.format( af_config["TARGET_INSTRUMENT"].upper() )
#                    
#                    images = []
#                    for instr_path in src_dset_path, target_dset_path:
#                        # Unfortunately, these ugly if statements are necessary because
#                        # each instrument function takes different parameters :(
#                        
#
#                        if "MODIS_Radiance" in instr_path:
#                            bands = af_config["MODIS_BANDS"]
#                            resolutions = af_config["MODIS_RESOLUTION"].upper().split()
#
#                            logger.info("Creating MODIS images...")
#                            
#                            if bands.upper() == "ALL":
#                                bands = range( NUM_MODIS_BANDS )
#                            else:
#                                # Split each band to a separate element, then convert 
#                                # each element to an integer
#                                bands = [ int(i) for i in bands.split() ]
#
#                            # band_idx and band are two different things. idx is 
#                            # index into dataset. Bands are assigned indices
#                            # according to which order they were given in the
#                            # config file.
#
#                            for resolution_idx, resolution in enumerate( resolutions ):
#                                for band_idx, band in enumerate(bands):
#
#                                    img_name = instr_path.replace('/', '#') + "_{}_band_{}".format(resolution, band) + \
#                                        "_orbit_{}".format( data.get_orbit() ) + IMG_EXTENSION
#
#                                    img_path = os.path.join( data.get_image_dir(), img_name )
#                                    
#                                    # Save image path into data
#                                    dict_walker = data.images
#                                    modis_dset = ["MODIS", resolution ]
#                                    for key in modis_dset:
#                                        if key not in dict_walker:
#                                            dict_walker[str(key)] = {}
#
#                                        dict_walker = dict_walker[str(key)]
#
#                                    logger.debug("data.images['MODIS'][{}][{}] = {}".format(
#                                        resolution, band, img_path))
#
#                                    data.images["MODIS"][ resolution ][ str(band) ] = img_path
#                                    
#                                    images.append(img_path)
#                                    
#                                    logger.debug(LOG_SEP)
#                                    logger.debug("out_path: {}".format( data.get_output_path()))
#                                    logger.debug("instr_path: {}".format(instr_path))
#                                    logger.debug("img_path: {}".format(img_path))
#                                    logger.debug("band: {}".format(band))
#                                    logger.debug("band_idx: {}".format(band_idx))
#
#                                    af_image.print_modis( data.get_output_path(), band_idx, instr_path, img_path )
#                      
#                            logger.debug(LOG_SEP)
#                            
#                        elif "MISR_Radiance" in instr_path:
#                            # Get MISR camera angles
#                            angles = af_config["MISR_CAMERA_ANGLE"].upper().split()
#                            bands = af_config["MISR_RADIANCE"].split()
#                            resolutions = af_config["MISR_RESOLUTION"].upper().split()
#
#                            logger.info("Creating MISR images.")
#
#                            # We use enumerate to get index, instead of using the actual string.
#                            # print_misr expects integer that serves as index into dataset. The
#                            # order in which angles/bands are listed in the config file is the
#                            # order in which they are stored in the dataset.
#                            for res_idx, res in enumerate( resolutions ):
#                                for angle_idx, angle in enumerate( angles ):
#                                    for band_idx, band in enumerate( bands ):
#                                        img_name = instr_path.replace('/', '#') + \
#                                            "_{}_{}_{}{}".format( angle, band.replace("_Radiance", "") , data.get_orbit(), 
#                                            IMG_EXTENSION )
#
#                                        img_path = os.path.join( data.get_image_dir(), img_name )
#
#                                        # Save image path into data.
#                                        # Recursively create a chain of dicts that
#                                        # describe the dataset. For instance, we would
#                                        # like to be able to access it like this:
#                                        # data.images["MODIS"]["L"]["AN"]["Red_Radiance"]
#                                        # 
#                                        # and have the path stored at that location.
#
#                                        dict_walker = data.images
#                                        misr_dset = ["MISR", res, angle ]
#                                        for key in misr_dset:
#                                            if key not in dict_walker:
#                                                dict_walker[str(key)] = {}
#
#                                            dict_walker = dict_walker[str(key)]
#                                            
#                                        data.images["MISR"][res][angle][str(band)] = img_path
#
#                                        logger.debug(LOG_SEP)
#                                        logger.debug("out_path: {}".format( data.get_output_path()))
#                                        logger.debug("angle: {}".format(angle))
#                                        logger.debug("angle_idx: {}".format( angle_idx))
#                                        logger.debug("band: {}".format(band))
#                                        logger.debug("band_idx: {}".format(band_idx))
#                                        logger.debug("instr_path: {}".format(instr_path))
#                                        logger.debug("img_path: {}".format( img_path))
#
#                                        af_image.print_misr( data.get_output_path(), angle_idx, 
#                                            band_idx, instr_path, img_path )
#
#                            
#                                logger.debug(LOG_SEP)
#            
#                        elif "ASTER_Radiance" in instr_path:
#                            # Get ASTER parameters
#                            res = af_config["ASTER_RESOLUTION"].upper().split()[0]
#                            bands = af_config["ASTER_BANDS"].upper().split()
#                            
#                            logger.info("Creating ASTER images.")
#                            # Landon Clipp 6/19/2018
#                            # Apparently, the ASTER_Radiance dataset only stores 
#                            # different bands in the dimensions. So just iterate
#                            # over the first dimension. Note: first for loop
#                            for band_idx, band in enumerate( bands ):
#                                img_name = "{}_{}_{}_{}{}".format( instr_path.replace('/', '#'),
#                                res, band, data.get_orbit(), IMG_EXTENSION )
#
#                                img_path = os.path.join( data.get_image_dir(), img_name )
#
#                                dict_walker = data.images
#                                aster_dset = ["ASTER", res ]
#                                for key in aster_dset:
#                                    if key not in dict_walker:
#                                        dict_walker[key] = {}
#
#                                data.images["ASTER"][res][band] = img_path
#
#                                logger.debug(LOG_SEP)
#                                logger.debug("out_path: {}".format( data.get_output_path()))
#                                logger.debug("band: {}".format(band))
#                                logger.debug("band_idx: {}".format(band_idx))
#                                logger.debug("instr_path: {}".format( instr_path ))
#                                logger.debug("img_path: {}".format(img_path))
#
#                                af_image.print_aster( data.get_output_path(), band_idx, 
#                                    instr_path, img_path )
#
#                            logger.debug(LOG_SEP)
#
#                        else:
#                            raise ValueError("Cannot determine the instrument \
#                            for image generation!")
#
#                    logger.info("Images generated.")
#
#                    if wf_config["COMPARE_IMAGES"]:
#                        logger.info("Creating similarity metrics.")
#
#                        for idx, val in enumerate( wf_config["COMPARE_IMAGES"] ) :
#                            # If it's a list, retrieve the cameras from each element
#                            if type(val) is list:
#                                dset1 = val[0].split()
#                                dset2 = val[1].split()
#
#                            # Else if it's a string, the next element is our second dataset
#                            elif type(val) is str:
#                                dest1 = val.split()
#                                dset2 = wf_config["COMPARE_IMAGES"][idx + 1].split()
#                            else:
#                                raise ValueError("Unknown value in COMPARE_IMAGES.")
#
#                            # "Recursively" iterate through the keys to retrieve the paths of the images.
#                            # All I'm doing is looping through the elements in dset1 and dset2 and retrieving
#                            # the intermediate values from data.images. For instance, 
#                            # say dset1 = ['MODIS', '1KM', '4']:
#                            # image1_path = data.images
#                            # image1_path = image1_path['MODIS']
#                            # image1_path = image1_path['1KM']
#                            # image1_path = image1_path['4']
#
#                            logger.debug(data.images)
#
#                            image1_path = data.images
#                            for key in dset1:
#                                image1_path = image1_path[str(key)]
#
#                            image2_path = data.images
#                            for key in dset2:
#                                image2_path = image2_path[str(key)]
#
#                            sim = af_image.structural_similarity( image1_path, image2_path )
#                            data.similarity.append(sim  )
#
#                            if sim < wf_config["COMPARE_THRESHOLD"]:
#                                logger.warning("Images have similarity below designated threshold!")
#                                logger.warning("image 1: {}".format( image1_path ) )
#                                logger.warning("image 2: {}".format( image2_path ) )
#                                logger.warning("similarity: {}".format( sim ) )
#                                logger.warning("similarity threshold: {}".format( wf_config["COMPARE_THRESHOLD"] ) )
#
#                                data.set_status( Tags.FAIL_IMG_CMP )
#                                return data
#
#                            else:
#                                logger.info("Images above similarity threshold.")
#                                logger.debug("image 1: {}".format( image1_path ) )
#                                logger.debug("image 2: {}".format( image2_path ) )
#                                logger.info("similarity: {}".format( sim ) )
#                                logger.info("similarity threshold: {}".format( wf_config["COMPARE_THRESHOLD"] ) )
        
            except:
                logger.exception("Encountered exception")
                return data

            data.set_status( Tags.SUCCEED_IMG )

        logger.info("Done.")

        data.set_status( Tags.SUCCEED_ALL )
        return data

    except:
        logger.exception("Encountered exception.")
        return data

def check_h5dump():
    """
    Check if h5dump is visible.
    """

    args = ['which', 'h5dump']
    proc = subprocess.Popen( args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
    proc.wait()

    with open( logger.get_log_file(), 'a') as f:
        for line in proc.stdout:
            print( line.decode('UTF-8' ) )
            f.write( line.decode('UTF-8') )
            
    if proc.returncode != 0:
        raise RuntimeError("h5dump not visible!")


def main(pool):
    parser = argparse.ArgumentParser(description="This is an MPI \
    application that distributes the Advanced Fusion (AF) tool across \
    multiple input files using a master-slave paradigm. You may only \
    run this script compute nodes.", 
    formatter_class=argparse.ArgumentDefaultsHelpFormatter )

    parser.add_argument("input_dir", help="Directory of the input \
    Basic Fusion files. Only give the top-level directory (will \
    perform a recursive search).", type=str)

    parser.add_argument("output_dir", help="Output AF directory. \
    Will be created if it doesn't exist", type=str)

    parser.add_argument("af_tool", help="Path to the AF binary.",
    type=str)

    req_grp = parser.add_argument_group(title='Required flags')
    
    req_grp.add_argument("--range", "-r", help="Specify the inclusive \
    range of orbits to process. May specify multiple ranges. At least \
    one range is required.", nargs=2, action='append', type=int, 
    required=True)
    
    req_grp.add_argument("--config", "-c", help="""Specify one or more
    yaml-compatible files that describe the parameters of the AF run.
    The parameters are either passed directly to the AF tool, or are 
    used only by this script. The set of recognized parameters is 
    described in the README of this directory. If multiple config files 
    are given, multiple calls will be made to the AF tool, each with 
    one output file.""", type=str, required=True, action='append' )
    
    parser.add_argument("--ll", help="Define the log output level.",
    type=str, choices=["CRITICAL", "ERROR", "WARNING", "INFO",
        "DEBUG"], default="DEBUG")
   
    parser.add_argument("--retry-failed", "-f", help="""If a previous run 
    of this script resulted in some failed granules, you may pass the
    failed_granlues.txt file back to this script to have it retry ONLY 
    those granules that failed. This serves as a "checkpoint" file. The 
    script will restart processing at the last successful step. 
    For instance, if the granule succeeded in generating an AF file but
    failed in generating an image, the script will attempt to generate 
    the image again (does not attempt to regenerate AF file). Note that 
    all of the files listed in failed_granules.txt must still exist.
    NOTE: This is a BETA feature and has some undocumented behavior!""",
    type=str )

    args = parser.parse_args()

    #Define the log level. Logger has already been defined globally,
    # but we need to add a few more parameters to it.
    ll = getattr( logging, args.ll )
    rootLogger = logging.getLogger()
    global logger
    rootLogger.setLevel( ll )

    #=================================
    # SLAVE ENTRY POINT              =
    #=================================
    if not pool.is_master():        #=
        pool.wait()                 #=
        return                      #=
    #=================================

    # Console handler should only be for master rank
    consoleHandler = logging.StreamHandler(sys.stdout)
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)
   
    logger.info("Creating output directory")
    logger.debug("Output dir: {}".format(args.output_dir))
    
    def create_dirs( dirs ):
        for i in out_dirs:
            # Create directories
            try:
                os.makedirs( out_dirs[i] )
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
    
    out_dirs = {}
    out_dirs['data'] = os.path.join(args.output_dir, 'data')
    out_dirs['logs'] = os.path.join(args.output_dir, 'logs')
    create_dirs(out_dirs)
    out_dirs['h5dump'] = os.path.join( args.output_dir, 'data', 'h5dump')
    out_dirs['run'] = make_run_dir( out_dirs['logs'] )
    out_dirs['images'] = os.path.join( out_dirs['data'], 'images' )
    out_dirs['af_log'] = os.path.join( out_dirs['run'], 'af_log' )
    out_dirs['configs'] = os.path.join( out_dirs['run'], 'configs' )
    create_dirs(out_dirs)

    
    # Create empty file in our run dir whose name is the current date
    now = datetime.datetime.now()
    cur_date = "{}_{}_{}.{}hr_{}min_{}sec".format( now.year, now.month, now.day,
        now.hour, now.minute, now.second)

    with open( os.path.join( out_dirs['run'], cur_date ), 'w' ) as f:
        pass
    
    
    logger.info("run dir: {}".format(out_dirs['run']))

    # Set file handler for master
    log_file = os.path.join( out_dirs['af_log'], 'master.log')
    logger.addFileHandler( log_file )

    # Print out command-line args
    logger.info( " ".join(sys.argv) )

    logger.info("Checking for h5dump visibility...")
    check_h5dump()
    logger.info("h5dump visible.")
    
    logger.info("Parsing input directory for BF files")
  
    configs = []
    for config_file in args.config:
        # Read the config file
        with open( config_file, 'r') as f:
            conf = yaml.load(f)
        configs.append(conf)

    # Discover all of the Basic Fusion files, making sure we have
    # all the files requested by orbit_min and orbit_max.
    # We create the orbits set to create a fast way of
    # determining which orbits we are not able to find.
    orbits = set()
    for r in args.range:
        orbits = orbits | set( range(r[0], r[1] + 1 ) )

    count = 0
    jobs = []

    seen_year_month_config = {}
    seen_year_month_log = {}
    seen_year_month_image = {}

    if args.retry_failed:
        logger.info("Retrying failed granules: {}".format( args.retry_failed))

        with open(args.retry_failed, "r") as f:
            granules = f.read().split('\n\n')

            for granule_idx, granule in enumerate(granules):
                if len(granule) == 0:
                    continue

                # Load granule string into dict using yaml
                granule_dict = yaml.load( granule )
                new_granule = Granule()

                for key in granule_dict:
                    setattr( new_granule, key, granule_dict[key] )

                jobs.append(new_granule)

    else:
        for root, dirs, files in os.walk( args.input_dir ):
            for file in files:
                
                try:
                    orbit = bfutils.file.bf_file_orbit( file )
                except ValueError:
                    # if bf_file_orbit gives ValueError, then file isn't 
                    # a basic fusion file. Discard it and move on.
                    continue
                
                        
                # If this granule is in our requested orbit range
                if orbit in orbits:
                    count = count + 1
                   
                    
                    for config_idx, config in enumerate(configs):
                        

                        # Print every so often
                        if count % 500 == 0:
                            logger.info("Orbit: {}".format(orbit))

                        granule = Granule()

                        in_path = os.path.join(root, file)
                        o_start = bfutils.file.orbit_start( 
                            bfutils.file.bf_file_orbit( in_path) )

                        year_month_dir = o_start[0:4] + '.' + o_start[4:6]

                        if "OUTPUT_PREFIX" in config:
                            prefix = config["OUTPUT_PREFIX"] + "_"
                        else:
                            prefix = "ADVNCE_FUSE_{}_".format(config_idx)

                        out_path = os.path.join( out_dirs['data'], 
                            year_month_dir, prefix + file )

                        # Set the configuration for this AF run
                        af_config = config.copy()
                        af_config[in_file_key] = in_path
                        af_config[out_file_key] = out_path
                        granule.set_config( af_config )

                        # Create the config year-month directory
                        config_ym_dir = os.path.join( out_dirs['configs'], 
                        year_month_dir )
                        
                        if config_ym_dir not in seen_year_month_config:
                            try:
                                os.makedirs( config_ym_dir )
                            except OSError as e:
                                if e.errno != errno.EEXIST:
                                    raise        
                            # Keep track of which ym dirs we've seen so we don't make
                            # more unnecessary (slow) system calls
                            seen_year_month_config[ config_ym_dir ] = None

                        config_file = os.path.join( config_ym_dir, str(orbit) 
                        + '_config_{}.txt'.format(config_idx) ) 

                        granule.set_af_config_file( config_file ) 

                        # Create the log year-month directory
                        log_ym_dir = os.path.join( out_dirs['af_log'], year_month_dir )
                        if log_ym_dir not in seen_year_month_log:
                            try:
                                os.makedirs( log_ym_dir )
                            except OSError as e:
                                if e.errno != errno.EEXIST:
                                    raise        
                            seen_year_month_log[log_ym_dir] = None

                        log_file = os.path.join( log_ym_dir, '{}_config_{}.log'.format(orbit, config_idx) )
                        granule.set_log_file( log_file )

                        granule.set_orbit(orbit)

                        granule.set_exe( args.af_tool )
                
                        # Create the image directory
                        image_orbit_dir = os.path.join( out_dirs['images'], year_month_dir, str(orbit) )
                        try:
                            os.makedirs( image_orbit_dir )
                        except OSError as e:
                            if e.errno != errno.EEXIST:
                                raise        
                        granule.set_image_dir( image_orbit_dir )

                        h5_path = os.path.join( out_dirs['h5dump'], year_month_dir, os.path.basename( out_path + '.h5dump'))
                        granule.set_h5dump( h5_path )

                        jobs.append( granule )
                    
                    orbits -= {orbit}

   
        logger.info("Done.")

        if orbits:
            logger.error("Could not find Basic Fusion files for \
            specified orbit range. Printing missing orbits to debug \
            output.")
        
            logger.debug(orbits)

            raise RuntimeError("BF files not found.")

    logger.info("Sending jobs to workers...")
    logger.info("Waiting for jobs to complete...")
    results = pool.map( worker, jobs )
    pool.close()
    logger.info("Done.")

    logger.info("Checking for failed granules...")
    
    
    logger.failed_file = os.path.join( out_dirs['run'], FAILED_FILE )
    logger.granule_summary = os.path.join( out_dirs['run'], GRANULE_SUMMARY )

    fail_count = 0

    # Sort the results by orbit
    results = sorted( results, key = lambda k: k.get_orbit() )
   
    with open( logger.granule_summary, 'w' ) as g_summary:
        with open( logger.failed_file, 'w' ) as f:
            for granule in results:

                attr_string = class_attr_string(granule) + '\n'
                g_summary.write( attr_string )

                # "If granule didn't succeed or (we are doing similarity
                # comparison and any similarity is below the threshold)"
                if granule.get_status() != Tags.SUCCEED_ALL or  \
                (granule.similarity and not granule.sim_above_thresh):
                    fail_count = fail_count + 1
                    f.write( attr_string )

    logger.info(
        "{}/{} ({} %) granules failed. See: {} for more details.".format(
        fail_count, len(results), round(100 * fail_count / len(results), 2), 
        logger.failed_file))
            
if __name__ == "__main__":
    try:
        pool = MPIPool()
        main(pool)
    except SystemExit as e:
        if e.code == 0:
            pool.close()
            raise
        else:
            logger.exception("Encountered exception")
            mpi_comm.Abort()
    except:
        logger.exception("Encountered exception")
        mpi_comm.Abort()
        
