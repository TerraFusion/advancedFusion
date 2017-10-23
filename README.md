# Advanced Terra Fusion Resampling Tool
Terra Data Fusion Project - University of Illinois

# Contributers
Yizhao Gao (ygao29@illinois.edu)  
Ya Long Lo (yllo2@illinois.edu)

# Introduction
This is the second phase of [NASA TERRA DATA FUSION](https://earthdata.nasa.gov/community/community-data-system-programs/access-projects/terra-data-fusion-products) project. It is designed to resample the radiance values from any arbitrary Terra instrument to any other Terra instument. It takes in the output file generated by [BasicFusion](https://github.com/TerraFusion/basicFusion) and generate a new HDF5 file with all radiance values resampled to the target instrument.  
For example, this tool can be used to resample:  
* MISR 1.1km resolution data to MODIS 1km resolution pixels
* ASTER 15m resolution data to MISR 250m resolution pixels 

# Key Features
1. All calculatations are based on spherical (great circle) distances using the latitude and longtitude of pixels
2. The input pixels are assumned to be irregularily spaced, which is especially true for MODIS
3. The granule for analysis is an entire orbit (one circle of the Terra satallite around the earth)
