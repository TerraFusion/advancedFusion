## Advanced Fusion parallel workflow

process.py wraps the AF tool in an MPI pool to allow for massively-parallel processing. In order to run it, you must first solve the dependencies in your Python environment by performing the following commands:

```
module load bwpy
module load bwpy-mpi
bwpy-environ

pip install --user bfutils
pip install --user schwimmbad
```

After this is done, you may begin to write your PBS script. An example has been provided under AFtest.pbs. The PBS script must perform two critical steps for meeting the runtime dependencies of the AF tool itself:

1. Perform a `module swap PrgEnv-cray PrgEnv-intel` (or equivalent) in order to use the Intel libraries
2. Export the paths of your gdal and HDF5 libraries to LD_LIBRARY_PATH. 


