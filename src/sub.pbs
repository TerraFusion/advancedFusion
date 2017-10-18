#!/bin/bash
### set the number of nodes
### set the number of PEs per node
#PBS -l nodes=512:ppn=32:xe
### set the wallclock time
#PBS -l walltime=36:00:00
### set the job name
#PBS -N af_test
### set the job stdout and stderr
#PBS -e $PBS_JOBID.err
#PBS -o $PBS_JOBID.out
### set email notification
##PBS -m bea
##PBS -M yllo2@illinois.edu

# NOTE: lines that begin with "#PBS" are not interpreted by the shell but ARE
# used by the batch system, wheras lines that begin with multiple # signs,
# like "##PBS" are considered "commented out" by the batch system
# and have no effect.

# If you launched the job in a directory prepared for the job to run within,
# you'll want to cd to that directory
# [uncomment the following line to enable this]
# cd $PBS_O_WORKDIR

# Alternatively, the job script can create its own job-ID-unique directory
# to run within.  In that case you'll need to create and populate that
# directory with executables and perhaps inputs
# [uncomment and customize the following lines to enable this behavior]
# mkdir -p /scratch/sciteam/$USER/$PBS_JOBID
# cd /scratch/sciteam/$USER/$PBS_JOBID
# cp /scratch/job/setup/directory/* .

# To add certain modules that you do not have added via ~/.modules
. /opt/modules/default/init/bash # NEEDED to add module commands to shell
#module load craype-hugepages2M  perftools
module load cray-hdf5
# export APRUN_XFER_LIMITS=1  # to transfer shell limits to the executable

### launch the application
### redirecting stdin and stdout if needed
### NOTE: (the "in" file must exist for input)

aprun ./testReproHDF5

### For more information see the man page for aprun