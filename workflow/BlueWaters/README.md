# Advanced Fusion parallel workflow

process.py wraps the AF tool in an MPI pool to allow for massively-parallel processing. In order to run it, you must first solve the dependencies in your Python environment by performing the following commands:

```
module load bwpy
module load bwpy-mpi
bwpy-environ

pip install --user bfutils
pip install --user schwimmbad
pip install --user globuslite
pip install --user batch4py
```

You must also resolve the dependencies of the AFtool by loading the proper environment modules. This step may differ depending on how you compiled the AFtool!

```
module swap PrgEnv-cray PrgEnv-intel
module load cray-hdf5/1.8.16
export LD_LIBRARY_PATH=/projects/sciteam/jq0/TerraFusion/gdal/lib/:$LD_LIBRARY_PATH
```

## Running in interactive mode

You may run the parallel workflow using interactive batch mode. To enter into this:

```
qsub -I -l walltime=2:00:00 -l nodes=1:ppn=32
```

After your interactive job has started, you can then run the process.py script's help message:


```
aprun -n 1 -b -- bwpy-environ -- python process -h
```

Note that the -b option is necessary to run Python scripts in Blue Waters.

An example of an actual run may look like this:

```
aprun -n 2 -d 16 -b -- bwpy-environ -- python process -r 37436 37437 ~/scratch/basicfusion ~/scratch/af_out ./config_file.txt ~/projects/advancedFusion/AFtool
```

This will run the application on 2 ranks, with 1 rank being a master and 1 being worker rank. Note that the AFtool uses a significant amount of RAM, so it's recommended not to use more than 4 ranks per node. 2 or 3 ranks per node is recommended.

## Running in batch mode

An example PBS script has been provided, "AFtest.pbs". This script requests 2 nodes and launches 4 ranks: 1 rank as the master and 3 as workers. This script may be launched on the login node as such:

```
qsub AFtest.pbs
```

Some of the parameters in the script may need to be changed to suit your environment.
