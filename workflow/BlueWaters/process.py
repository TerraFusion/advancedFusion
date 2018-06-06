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

    def add_input_path(self, path):
        self.input_path = path
    
    def add_output_path(self, path):
        self.out_path = path
    
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
        return self.job_config

    def set_log_file(self, path):
        self.log_file = path

    def get_log_file(self):
        return self.log_file

    def set_exe(self, exe):
        """
        Set path to the executable to run.
        """
        self.exe = exe

    def get_exe(self):
        return self.exe

    def set_orbit(self, orbit):
        self.orbit = orbit

    def get_orbit(self):
        return self.orbit

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


#=================== MPI VARIABLES =====================
mpi_comm = MPI.COMM_WORLD
mpi_rank = MPI.COMM_WORLD.Get_rank()
mpi_size = MPI.COMM_WORLD.Get_size()
rank_name = "Rank {}".format(mpi_rank)

#================== LOGGING ======================
logging.setLoggerClass(Log)
logFormatter = logging.Formatter( LOG_FMT, LOG_DATE_FMT)
consoleHandler = logging.StreamHandler(sys.stdout)
consoleHandler.setFormatter(logFormatter)
rootLogger = logging.getLogger()
rootLogger.addHandler(consoleHandler)
logger = logging.getLogger(name = rank_name)

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

def worker( data ):
    logger.info("Received data.")
    logger.info("Creating config file.")
    # First need to create the config file
    with open( data.get_config_file(), 'w' ) as f:
        config = data.get_config()
        for key in config:
            f.write("{}: {}\n".format( key, config[key] ) )

    # Create the yyyy.mm dir in the output file path
    try:
        os.makedirs( os.path.dirname( data.get_config()['OUTPUT_FILE_PATH'] ))
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    logger.info("Calling the AFtool")
    # Now we can call the executable
    args = [ data.get_exe(), data.get_config_file() ]
    logger.debug(' '.join(args))
    try:
        with open( data.get_log_file(), 'w' ) as f:
            subprocess.check_call( args, stdout=f, 
                stderr=subprocess.STDOUT)
    except:
        logger.error("Encountered exception when processing AF \
granule: {}.\nSee: {}\nfor more details.".format(
            data.get_orbit(), data.get_log_file()))
        raise

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

    parser.add_argument("orbit_min", help="Inclusive lower bound \
    for the orbit range.", type=int)

    parser.add_argument("orbit_max", help="Inclusive upper bound \
    for the orbit range", type=int)

    parser.add_argument("config", help="A yaml-compatible file that \
    describes the parameters to pass to the AF tool. The set of \
    recognized parameters is described on the AF GitHub. NOTE: \
    The OUTPUT_FILE_PATH and INPUT_FILE_PATH parameters will be \
    IGNORED by this script.", type=str, default="./AFconfig.txt")

    parser.add_argument("af_tool", help="Path to the AF binary.",
    type=str)

    parser.add_argument("--ll", help="Define the log output level.",
    type=str, choices=["CRITICAL", "ERROR", "WARNING", "INFO",
        "DEBUG"], default="INFO")

    args = parser.parse_args()


    #Define the log level. Logger has already been defined globally,
    # but we need to add a few more parameters to it.
    ll = getattr( logging, args.ll )
    global rootLogger
    global logger

    rootLogger.setLevel( ll )

    logger.info("Creating output directory")
    logger.debug("Output dir: {}".format(args.output_dir))
    
    out_dirs = { 'data': os.path.join(args.output_dir, 'data'), 
                 'logs': os.path.join(args.output_dir, 'logs')}
    for i in out_dirs:
        # Create directories
        try:
            os.makedirs( out_dirs[i] )
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    # Setup log directory for this run
    logger.create_log_dirs( out_dirs['logs'] )
    logger.debug("run dir: {}".format(logger.get_run_dir()))


    logger.info("Parsing input directory for BF files")
    
    # Read the config file
    with open( args.config, 'r') as f:
        config = yaml.load(f)
    
    # Discover all of the Basic Fusion files, making sure we have
    # all the files requested by orbit_min and orbit_max.
    # We create the orbits dictionary to create a fast way of
    # determining which orbits we are not able to find.
    orbits = {}
    for i in range( args.orbit_min, args.orbit_max + 1 ):
        orbits[i] = None

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
                    job_config['INPUT_FILE_PATH'] = in_path
                    job_config['OUTPUT_FILE_PATH'] = out_path

                    granule.set_config( job_config )
                    config_file = os.path.join( log_dirs['configs'], str(orbit) + '_config.txt' ) 
                    granule.set_config_file( config_file ) 

                    log_file = os.path.join( log_dirs['af_log'], '{}.log'.format(orbit) )
                    granule.set_log_file( log_file )
                    granule.set_orbit(orbit)

                    granule.set_exe( args.af_tool )
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

    logger.info( logger.get_run_dir() )
    logger.info("Sending jobs to workers...")
    logger.info("Waiting for jobs to complete...")
    results = pool.map( worker, jobs )
    logger.info("Done.")

if __name__ == "__main__":
    try:
        pool = MPIPool()
        
        if not pool.is_master():
            pool.wait()
        else:
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
        
