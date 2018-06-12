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

# Define communication tags
Tags = IntEnum('Tags', 'FAIL SUCCEED' )

LOG_FMT='[%(asctime)s] [%(name)12s] [%(levelname)8s] [%(filename)s:%(lineno)d] %(message)s'
LOG_DATE_FMT='%d-%m-%Y:%H:%M:%S8'

class Granule:
    """
    Class that holds all information necessary for a specific orbit of BF->AF
    data.
    """

    def __init__(self):
        self.input_path = None
        self.out_path = None
        self.job_config = {}

    def set_input_path(self, path):
        self.job_config[in_file_key] = path
   
    def get_input_path(self):
        return self.job_config[in_file_key]

    def set_output_path(self, path):
        self.job_config[out_file_key] = path
   
    def get_output_path(self):
        return self.job_config[out_file_key]

    def set_config_file(self, file ):
        """
        Set the path of the config file
        """
        self.job_config_file = file

    def get_config_file( self ):
        """
        Return path to the job config file
        """
        return self.job_config_file

    def set_config(self, config):
        """
        Stores a dict that will be placed into self.job_config_file
        (to be passed directly to the advanced fusion tool)
        """
        self.job_config = config

    def get_config(self):
        """
        Returns a dict that represents the configuration file passed
        to the AF tool.
        """
        return self.job_config

    def set_log_file(self, path):
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

    def set_run_dir(self, dir):
        self.log_dirs['run'] = dir

    def get_run_dir(self):
        return self.log_dirs['run']

    def create_log_dirs( self, out_path ):
        """
        Create and initialize all directories necessary to store
        log information for this run.
        """
        self.log_dirs['run'] = make_run_dir( out_path )
       
        dirs = ['configs', 'af_log']

        for dir in dirs:
            self.log_dirs[dir] = os.path.join( self.get_run_dir(), dir)
            os.makedirs( self.log_dirs[dir] )
        
        # Create empty file whose name is the current date
        now = datetime.datetime.now()
        cur_date = "{}_{}_{}.{}hr_{}min_{}sec".format( now.year, now.month, now.day,
            now.hour, now.minute, now.second)

        with open( os.path.join( self.get_run_dir(), cur_date ), 'w' ) as f:
            pass

    def get_log_dirs(self):
        """
        Return the dict that stores paths for all the log dirs
        """
        return self.log_dirs

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

#================== LOGGING ======================
logging.setLoggerClass(Log)
logFormatter = logging.Formatter( LOG_FMT, LOG_DATE_FMT)
logger = logging.getLogger(name = rank_name)

#================== CONSTANTS =========================
FAILED_FILE = 'failed_granules.txt'
out_file_key = 'OUTPUT_FILE_PATH'
in_file_key = 'INPUT_FILE_PATH'

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

    args = ['h5dump', '-H', granule.get_config()[out_file_key] ]
    
    with open( granule.get_h5dump(), 'w' ) as f:
        subprocess.check_call( args, stdout=f, stderr=subprocess.STDOUT)

def do_af( granule):
    """
    Generate advanced fusion output for a granule.
    """
    # Now we can call the executable
    args = [ granule.get_exe(), granule.get_config_file() ]
    logger.debug(' '.join(args))
    with open( granule.get_log_file(), 'a' ) as f:
        subprocess.check_call( args, stdout=f, 
            stderr=subprocess.STDOUT)
    
def worker( data ):

    global logger
    ret = Tags.SUCCEED

    # Remove any old file handlers from previous orbits
    logger.handlers = [ h for h in logger.handlers if (type(h) != logging.FileHandler)]

    # Set the per-rank file handler for this specific orbit
    fileHandler = logging.FileHandler( data.get_log_file(), mode='a')
    fileHandler.setFormatter( logFormatter )
    logger.addHandler( fileHandler)

    data_str = ""
    for attr, val in data:
        data_str = data_str + "{}: {}\n".format(attr, val)

    logger.info("Received data: \n{}".format( data_str ))
    logger.info("Creating config file: {}".format(data.get_config_file()))
    # First need to create the config file
    with open( data.get_config_file(), 'w' ) as f:
        config = data.get_config()
        for key in config:
            f.write("{}: {}\n".format( key, config[key] ) )

    # Create the yyyy.mm dir in the output file path
    try:
        os.makedirs( os.path.dirname( data.get_config()[out_file_key] ))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    #=======================
    # advanced fusion
    #

    logger.info("Calling the AFtool")
    try:
        do_af(data)
    except subprocess.CalledProcessError as e:
        logger.exception("Encountered exception when processing AF \
granule: {}.\nSee: {}\nfor more details.".format(
            data.get_orbit(), data.get_log_file()))

        ret = Tags.FAIL

    #===============
    # h5dump
    #
    logger.info("Performing h5dump on file...") 
    try:
        do_h5dump( data )
    except subprocess.CalledProcessError as e:
        logger.exception("Encountered exception when processing h5dump: {}.\
\nSee: {}\nfor more details.".format(
            data.get_orbit(), data.get_h5dump()))

        ret = Tags.FAIL

    logger.info("Done.")
    logger.info("Granule generated {}.".format( 
        'successfully' if ret == Tags.SUCCEED else 'unsuccessfully'))
    return (data, ret)

def check_h5dump():
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

    parser.add_argument("config", help="A yaml-compatible file that \
    describes the parameters to pass to the AF tool. The set of \
    recognized parameters is described on the AF GitHub. NOTE: \
    The OUTPUT_FILE_PATH and INPUT_FILE_PATH parameters will be \
    IGNORED by this script.", type=str, default="./AFconfig.txt")

    req_grp = parser.add_argument_group(title='Required flags')
    req_grp.add_argument("--range", "-r", help="Specify the inclusive \
    range of orbits to process. May specify multiple ranges. At least \
    one range is required.", nargs=2, action='append', type=int, 
    required=True)

    parser.add_argument("af_tool", help="Path to the AF binary.",
    type=str)

    parser.add_argument("--ll", help="Define the log output level.",
    type=str, choices=["CRITICAL", "ERROR", "WARNING", "INFO",
        "DEBUG"], default="INFO")

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
    
    out_dirs = { 'data': os.path.join(args.output_dir, 'data'), 
                 'logs': os.path.join(args.output_dir, 'logs'),
                 'h5dump': os.path.join( args.output_dir, 'data', 'h5dump')
               }
    for i in out_dirs:
        # Create directories
        try:
            os.makedirs( out_dirs[i] )
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    # Setup log directory for this run
    logger.create_log_dirs( out_dirs['logs'] )
    logger.info("run dir: {}".format(logger.get_run_dir()))

    # Set file handler for master
    log_file = os.path.join( logger.get_log_dirs()['af_log'], 'master.log')
    logger.addFileHandler( log_file )

    logger.info("Checking for h5dump visibility...")
    check_h5dump()
    logger.info("h5dump visible.")
    
    logger.info("Parsing input directory for BF files")
    
    # Read the config file
    with open( args.config, 'r') as f:
        config = yaml.load(f)
    
    # Discover all of the Basic Fusion files, making sure we have
    # all the files requested by orbit_min and orbit_max.
    # We create the orbits dictionary to create a fast way of
    # determining which orbits we are not able to find.
    orbits = {}
    for r in args.range:
        for o in range( r[0], r[1] + 1 ):
            orbits[o] = None

    count = 0
    jobs = []
    log_dirs = logger.get_log_dirs()
    for root, dirs, files in os.walk( args.input_dir ):
        for file in files:
            try:

                orbit = bfutils.file.bf_file_orbit( file )
                
                # If this granule is in our requested orbit range
                if orbit in orbits:
                    granule = Granule()

                    in_path = os.path.join(root, file)
                    o_start = bfutils.file.orbit_start( 
                        bfutils.file.bf_file_orbit( in_path) )

                    year_month_dir = o_start[0:4] + '.' + o_start[4:6]
                    out_path = os.path.join( out_dirs['data'], 
                        year_month_dir, "ADVNCE_FUSE_" + file )

                    job_config = config.copy()
                    job_config[in_file_key] = in_path
                    job_config[out_file_key] = out_path

                    granule.set_config( job_config )
                    config_file = os.path.join( log_dirs['configs'], str(orbit) + '_config.txt' ) 
                    granule.set_config_file( config_file ) 

                    log_file = os.path.join( log_dirs['af_log'], '{}.log'.format(orbit) )
                    granule.set_log_file( log_file )
                    granule.set_orbit(orbit)

                    granule.set_exe( args.af_tool )

                    h5_path = os.path.join( out_dirs['h5dump'], year_month_dir, os.path.basename( out_path + '.h5dump'))
                    granule.set_h5dump( h5_path )
                    jobs.append( granule )
                    orbits.pop(orbit)

            except ValueError:
                # Catches ValueError from bfutils call
                pass
   
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
    
    
    logger.failed_file = os.path.join( logger.get_run_dir(), FAILED_FILE )
    fail_count = 0
    with open( logger.failed_file, 'w' ) as f:
        for granule in results:
            if granule[1] == Tags.FAIL:
                fail_count = fail_count + 1
                f.write(
                    'Orbit: {}\nInput: {}\nOutput: {}\nLog file: {}\nConfig file: {}\nh5dump: {}\n\n'.format(
                    granule[0].get_orbit(), granule[0].get_input_path(),
                    granule[0].get_output_path(), granule[0].get_log_file(),
                    granule[0].get_config_file(), granule[0].get_h5dump()))


    if fail_count:
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
        
