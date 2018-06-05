import argparse
import logging
import os
import errno
import datetime
import bfutils
import sys
from mpi4py import MPI
import re
import mpi_master_slave as mpi_ms
import time

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

class Log(logging.Logger):
    """
    Class that contains logging utilites and metadata.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger = None
        self.run_dir = None
        self.level = None

    def set_level(self, level):
        self.level = level

    def set_run_dir(self, dir):
        self.run_dir = dir

    def get_run_dir(self):
        return self.run_dir

    def set_log_dir(self, dir):
        self.log_dir = dir

    def get_log_dir(self):
        return self.log_dir

    def get_log_struct(self):
        return self.log_struct

    def set_log_struct(self, struct):
        self.log_struct = struct

class MyApp:
    """
    This is the main application that runs the master-slave workflow. 
    This is where all of the heavy duty computation happens.
    """
    def __init__(self, slaves):
        # When creating the master we tell it whats slaves it can handle
        self.master = mpi_ms.Master(slaves)
        self.work_queue = mpi_ms.WorkQueue(self.master)

    def terminate_slaves(self):
        self.master.terminate_slaves()

    def run(self, tasks=10):
        """
        This is the core of the application. Keep starting slaves 
        as long as there is work do to.
        """

        for i in range(tasks):
            self.work_queue.add_work( data = ('Do task', i) )

        # Keep starting slaves as long as there is work to do
        while not self.work_queue.done():
            self.work_queue.do_work()
            
            # Reclaim returned data
            for return_data in self.work_queue.get_completed_work():
                done, message = return_data

                if done:
                    logger.info("Slaved finished task and says {}".format(message))

            time.sleep(0.1)

class MySlave(mpi_ms.Slave):
    """
    Slave process that extends the Slave class. Overrides the 'do_work"
    method and calls 'Slave.run'. The master will do the rest.
    """
    def __init__(self):
        super().__init__()

    def do_work( self, data ):
        logger.info("Data: {}".format(data))
        return True, "Everything is good!"

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


def main():
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

    #----------- MPI SLAVE -----------
    if mpi_rank != 0:
        try:
            MySlave().run()              
        except:
            logger.exception("Encountered exception.")
            mpi_comm.Abort()
        finally:
            return
            
    #---------------------------------

    logger.info("Creating output directory")
    logger.debug("Output dir: {}".format(args.output_dir))
    # Create directories
    try:
        os.makedirs( args.output_dir )
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    # Create run directory
    logger.set_run_dir( make_run_dir( args.output_dir ) )
    logger.debug("run dir: {}".format(logger.get_run_dir()))

    # Create subdirectories in run directory
    run_structure = { 'log': None }

    for dir in run_structure:
        path = os.path.join( logger.get_run_dir(), dir )
        os.makedirs( path )

    # Create empty file whose name is the current date
    now = datetime.datetime.now()
    cur_date = "{}_{}_{}.{}hr_{}min_{}sec".format( now.year, now.month, now.day,
        now.hour, now.minute, now.second)

    with open( os.path.join( logger.get_run_dir(), cur_date ), 'w' ) as f:
        pass

    logger.info("Parsing input directory for BF files")
    
    # Discover all of the Basic Fusion files, making sure we have
    # all the files requested by orbit_min and orbit_max.
    # We create the orbits dictionary to create a fast way of
    # determining which orbits we are not able to find.
    orbits = {}
    for i in range( args.orbit_min, args.orbit_max + 1 ):
        orbits[i] = None

    for root, dirs, files in os.walk( args.input_dir ):
        for file in files:
            try:
                orbit = bfutils.file.bf_file_orbit( file )
                
                try:
                    # Remove this key from the dict
                    orbits.pop(orbit)
                except KeyError:
                    pass

            except ValueError:
                # Catches ValueError from bfutils call
                pass
    
    if orbits:
        logger.error("Could not find Basic Fusion files for \
        specified orbit range. Printing missing orbits to debug \
        output (see log file).")

        raise RuntimeError("BF files not found.")

    app = MyApp( slaves=range(1, mpi_size) )
    
    try:
        app.run()
    except:
        logger.exception("Encountered exception")
        mpi_comm.Abort()
    
    app.terminate_slaves()

if __name__ == "__main__":
    main()
