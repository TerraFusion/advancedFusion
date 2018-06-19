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

## Configuration file

process.py requires a configuration file that describes the parameters of the run. The parameters given are either consumed directly by the script, or are passed to the AF tool itself. The parameters that the AF tool accepts are described in its documentation. process.py will ignore the keys INPUT_FILE_PATH and OUTPUT_FILE_PATH, instead relying on internal logic for determining the input and output file paths. This is done because it is unreasonable for the user to specify every single input and output path.

The keys that this script accepts are:

### COMPARE_IMAGES  

Provide a list of resampled datasets to compare using a structural similarity index. The template for each instrument is:

**MISR**

MISR [MISR_RESOLUTION] [MISR_CAMERA_ANGLE] [MISR_RADIANCE] 

where each [] variable is a value you have already specified elsewhere in the config file.

**MODIS**

MODIS [MODIS_RESOLUTION] [MODIS_BANDS] 

again, where each [] variable is a value already specified in the config file.

**ASTER**

ASTER [ASTER_RESOLUTION] [ASTER BAND]

An example of a line in the config file that directs the script to compare MISR AN Red_Radiance L vs MODIS 4 1KM, and MISR AN Blue_Radiance L vs MODIS 3 1KM is:

```
COMPARE_IMAGES: [ ['MISR L AN Red_Radiance', 'MODIS 1KM 4'], ['MISR L AN Blue_Radiance', 'MODIS 1KM 3'] ]
```

If you are requesting only a single comparison, you may do something like:

```
COMPARE_IMAGES: ['MISR L AN Red_Radiance', 'MODIS 1KM 4']
```

### COMPARE_THRESHOLD

This parameter describes the threshold for the structural similarity percentage, below which an image comparison will be marked as bad. If this is not provided explicitly, the parameter will default to 0.90. Must be a value between 0 and 1.

### OUTPUT_PREFIX

This parameter sets the prefix for the output Advanced Fusion files. This is useful (read: EXTERMELY useful) if you are using multiple configuration files during the same process.py run and want to differentiate your files from each other.
