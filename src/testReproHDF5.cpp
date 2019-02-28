/**
 * testRepro.c
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {07/17/2017}
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>
#include <sys/time.h>
#include "reproject.h"
#include "io.h"

int main(int argc, char ** argv) 
{
	hid_t output_file = H5Fcreate("misr_modis_test_repro.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	char* file_path = "/projects/TDataFus/kent/temp/40-orbit-file/Jun15.2/TERRA_BF_L1B_O69365_F000_V000.h5";
	hid_t file;
	if(0 > (file = af_open(file_path))) {
		printf("File not found\n");
		exit(1);
	}

	int nCellMODIS;
	int nCellMISR;

	double * MISR_Lat = get_misr_lat(file, "L", &nCellMISR);
	double * MISR_Lon = get_misr_long(file, "L", &nCellMISR);

	double * MISR_Rad;

	double * MODIS_Lat = get_modis_lat(file, "_1KM", "EV_1KM_RefSB", &nCellMODIS); 
	double * MODIS_Lon = get_modis_long(file, "_1KM", "EV_1KM_RefSB", &nCellMODIS);

	printf("MISR CELLS: %d\n, MODIS CELLS: %d\n", nCellMISR, nCellMODIS);

	double * MODIS_Rad_Out;

	int * tarNNSouID;

	//MISR TO MODIS NN
	
	double ** p_MISR_Lat = &MISR_Lat;
	double ** p_MISR_Lon = &MISR_Lon;

	tarNNSouID = (int *)malloc(sizeof(int) * nCellMODIS);

	//Finding nearest points
	printf("nearest_neighbor\n");
	nearestNeighbor(p_MISR_Lat, p_MISR_Lon, nCellMISR, MODIS_Lat, MODIS_Lon, tarNNSouID, nCellMODIS, 1000);

	MISR_Lat = * p_MISR_Lat;
	MISR_Lon = * p_MISR_Lon;

	free(MISR_Lat);
	free(MISR_Lon);

	printf("writing modis geo\n");
	// modis 1KM output width is '1354'
	int lat_status =  af_write_mm_geo(output_file, 0, MODIS_Lat, nCellMODIS, 1354,-1,-1);
	int long_status = af_write_mm_geo(output_file, 1, MODIS_Lon, nCellMODIS, 1354,-1,-1);

	free(MODIS_Lat);
	free(MODIS_Lon);

	printf("getting misr\n");
	MISR_Rad = get_misr_rad(file, "AN", "L", "Blue_Radiance", &nCellMISR);
	int nCellMODIS_rad;
	double* MODIS_Rad = get_modis_rad(file, "_1KM", "EV_1KM_RefSB", &nCellMODIS_rad);
	
	////here
	MODIS_Rad_Out = (double *)malloc(sizeof(double) * nCellMODIS);
	
	printf("interpolating\n");
	nnInterpolate(MISR_Rad, MODIS_Rad_Out, tarNNSouID, nCellMODIS);
	printf("writing data fields\n");
	int data_write_status = af_write_misr_on_modis(output_file, MODIS_Rad_Out, MODIS_Rad, nCellMODIS_rad, nCellMODIS);

	free(MISR_Rad);
	free(MODIS_Rad_Out);

	free(tarNNSouID);

	herr_t ret = af_close(file);

	return 0;
}
