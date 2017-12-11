/*
 * PROGRMMER: 
 *	- Jonathan Kim (jkm@illinois.edu)
 * 
 * DESCRIPTION:
 *	This is a advance fusion tool that generates resampled data from Terra satellite.
 *	Uses I/O functions and resampling functions.
 *
 * INPUT: 
 *	A input text file, which points to a Basic Fusion orbit data file (HDF5 format).
 *
 * OUTPUT:
 *	A resampled orbit data file as HDF5 format which contains source and target 
 *	instrument data. 
 */

#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include <hdf5.h>
#include "io.h"
#include "reproject.h"
#include "AF_InputParmeterFile.h"
#include "AF_debug.h"

const hsize_t RESAMPLE_DATA_WIDTH = 1354;

void Usage(int &argc, char *argv[])
{
	std::cout	<< "Usage: \n"
				<< "   " << argv[0] << "  <parameter-input-file>\n";
}


int main(int argc, char *argv[])
{

	if (argc < 2) {
		Usage(argc, argv);
		return 1;
	}

	//----------------------------------
	// parse input parameter from file 
	AF_InputParmeterFile inputArgs;
	inputArgs.headerFileName = argv[1];
	inputArgs.ParseByLine();

	//-------------------------
	// get instrument names
	std::string srcInstrument = inputArgs.GetSourceInstrument();
	std::string targetInstrument = inputArgs.GetTargetInstrument();
	
	//----------------------
	// create output file
	std::string outputFile = inputArgs.GetOuputFilePath();
	#if DEBUG_TOOL
	std::cout << "DBG main> outputFile: " << outputFile << std::endl;
	#endif

	hid_t output_file = H5Fcreate(outputFile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


	//-----------------------------
	// handle input BF data file
	std::string inputDataPath = inputArgs.GetInputBFdataPath();
	#if DEBUG_TOOL
	std::cout << "DBG main> inputDataPath: " << inputDataPath << std::endl;
	#endif
	hid_t src_file;
	if(0 > (src_file = af_open((char*)inputDataPath.c_str()))) {
		std::cerr << "Error: File not found - " << inputDataPath << std::endl;
		exit(1);
	}
	
	//-------------------------------------------------
	// handle Source instrument latitude and longitude
	std::cout << "Getting source instrument latitude longitude  data...\n";
	std::string misr_resolution = inputArgs.GetMISR_Resolution();
	#if DEBUG_TOOL
	std::cout << "DBG main> misr_resolution: " << misr_resolution << std::endl;
	#endif
	int nCellsrc;
	double* srcLatitude = NULL;
	double* srcLongitude = NULL;
	srcLatitude = get_misr_lat(src_file, (char*)misr_resolution.c_str(), &nCellsrc);
	srcLongitude = get_misr_long(src_file, (char*)misr_resolution.c_str(), &nCellsrc);

	//-------------------------------------------------
	// handle Target instrument latitude and longitude
	std::string modis_resolution = inputArgs.GetMODIS_Resolution();
	std::cout << "Getting target instrument latitude longitude  data...\n";
	#if DEBUG_TOOL
	std::cout << "DBG main> modis_resolution: " << modis_resolution << std::endl;
	#endif
	int nCelldest;
	double* targetLatitude = NULL;
	double* targetLongitude = NULL;
	targetLatitude = get_modis_lat(src_file, (char*)modis_resolution.c_str(), &nCelldest);
	targetLongitude = get_modis_long(src_file, (char*)modis_resolution.c_str(), &nCelldest);
	
	std::cout << "Writing target geolocation data...\n";
	int lat_status =  af_write_mm_geo(output_file, 0, targetLatitude, nCelldest);
	if(lat_status < 0) {
		std::cerr << "Error: writing latitude geolocation.\n";
		return -1;
	}
	int long_status = af_write_mm_geo(output_file, 1, targetLongitude, nCelldest);
	if(long_status < 0) {
		std::cerr << "Error: writing longitude geolocation.\n";
		return -1;
	}

	int * targetNNsrcID = NULL;
	targetNNsrcID = new int [nCelldest];
	
	std::cout <<  "Running nearest neighbor block index method... \n";
	nearestNeighborBlockIndex(&srcLatitude, &srcLongitude, nCellsrc, targetLatitude, targetLongitude, targetNNsrcID, NULL, nCelldest, 1000);
	//nearestNeighbor(&srcLatitude, &srcLongitude, nCellsrc, targetLatitude, targetLongitude, targetNNsrcID, nCelldest, 1000);

	if(srcLatitude)
		free(srcLatitude);
	if(srcLongitude)
		free(srcLongitude);
	if(targetLatitude);
		free(targetLatitude);
	if(targetLongitude)
		free(targetLongitude);
	
	//-------------------------------------------------
	// Get Source instrument radiance
	std::cout << "Getting source instrument radiance data...\n";
	double* src_rad=NULL;
	std::string misr_radiance = inputArgs.GetMISR_Radiance();
	std::vector<std::string> misr_cameraAngles = inputArgs.GetMISR_CameraAngles();
	#if DEBUG_TOOL
	std::cout << "DBG main> misr_radiance: " << misr_radiance << std::endl;
	for(int i = 0; i < misr_cameraAngles.size(); i++) {
		std::cout << "DBG main> misr_cameraAngles[" << i << "]:" << misr_cameraAngles[i] << std::endl;
	}
	#endif

	src_rad = get_misr_rad(src_file, (char*)misr_cameraAngles[0].c_str(), (char*)misr_resolution.c_str(), (char*)misr_radiance.c_str(), &nCellsrc);
	if (src_rad == NULL) {
		std::cerr	<< "Please verify these MISR input values. \n" 
					<< "	   - Camera: " << misr_cameraAngles[0] << "\n" 
					<< "	   - Resolution: " << misr_resolution << "\n"
					<< "	   - Radiance: " << misr_radiance << "\n" 
					<< std::endl;
		return -1;
	}
	
	//-------------------------------------------------
	// Get Target instrument radiance
	std::cout << "Getting target instrument radiance data...\n";
	int nCelldest_rad;
	double* dest_rad = NULL;
	char bands[1][50];
	//std::vector<std::string> bands; // TODO: use this instead of above

	// handle target instrument bands
	std::vector<int>  modis_bands = inputArgs.GetMODIS_Bands();
	sprintf(bands[0], "%d", modis_bands[0]);
	#if DEBUG_TOOL
	for(int i = 0; i < modis_bands.size(); i++) {
		std::cout << "DBG main> modis_bands[" << i << "]:" << modis_bands[i] << std::endl;
	}
	std::cout << "DBG main> bands[0]:" << bands[0] << std::endl;
	#endif

	dest_rad = get_modis_rad(src_file, (char*)modis_resolution.c_str(), bands, 1, &nCelldest_rad);
	
	//-------------------------------------------------
	// handle resample method
	double* src_rad_out = new double [nCelldest];
	int new_src_size = nCelldest;
	int * nsrcPixels = NULL;
	//Interpolating
	std::string resampleMethod =  inputArgs.GetResampleMethod();
	std::cout << "Interpolating using '" << resampleMethod << "' method...\n";
	if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "nnInterpolate")) {
		nnInterpolate(src_rad, src_rad_out, targetNNsrcID, nCelldest);
	}
	else if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "summaryInterpolate")) {
		nsrcPixels = new int [nCelldest];
		summaryInterpolate(src_rad, targetNNsrcID, nCellsrc, src_rad_out, nsrcPixels, nCelldest);
		#if DEBUG_TOOL
		std::cout << "DBG> No nodata values: \n";
		for(int i = 0; i < nCelldest; i++) {
			if(nsrcPixels[i] > 0) {
				printf("%d,\t%lf\n", nsrcPixels[i], src_rad_out[i]);
			}
		}
		#endif
	}

	if (targetNNsrcID)
		delete [] targetNNsrcID;


	//-------------------------------------------------
	// Writing target and source data to output file
	hid_t group_id = H5Gcreate2(output_file, "/Data_Fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if(group_id < 0) {
		std::cerr <<  "Error: H5Gcreate2 in output file.\n";
		return -1;
	}

	//------------------------------------
	//Write target instrument data first
	std::cout << "Writing target instrument '" << targetInstrument << "' data...\n";
	hsize_t target_dim[3];
	target_dim[0] = 1;
	target_dim[2] = RESAMPLE_DATA_WIDTH;
	target_dim[1] = (nCelldest_rad)/1/RESAMPLE_DATA_WIDTH;
	hid_t target_dataspace = H5Screate_simple(3, target_dim, NULL);
	hid_t	target_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t	target_status = H5Tset_order(target_datatype, H5T_ORDER_LE);  
	hid_t target_dataset = H5Dcreate2(output_file, "/Data_Fields/modis_rad", target_datatype, target_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if(target_dataset < 0) {
		std::cerr <<  "Error: H5Dcreate2 target data in output file.\n";
		return -1;
	}
	target_status = H5Dwrite(target_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dest_rad);
	if(target_status < 0) {
		std::cerr  <<  "Error: H5Dwrite target data in output file.\n";
		return -1;
	}
	H5Sclose(target_dataspace);
	H5Tclose(target_datatype);
	H5Dclose(target_dataset);

	if (dest_rad)
		free(dest_rad);

	//--------------------------------------------
	// Write reampled source instrument data next
	std::cout << "Writing source instrument '" << srcInstrument << "' data...\n";
	hsize_t src_dim[2];
	src_dim[0] = (new_src_size) / RESAMPLE_DATA_WIDTH;
	src_dim[1] = RESAMPLE_DATA_WIDTH;
	hid_t src_dataspace = H5Screate_simple(2, src_dim, NULL);
	hid_t src_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t src_status = H5Tset_order(src_datatype, H5T_ORDER_LE);  
	hid_t src_dataset = H5Dcreate2(output_file, "/Data_Fields/misr_out", src_datatype, src_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if(src_dataset < 0){
		std::cerr <<  "Error: H5Dcreate2 source data in output file.\n";
		return -1;
	}
	src_status = H5Dwrite(src_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_rad_out);
	if(src_status < 0){
		std::cerr  <<  "Error: H5Dwrite source data in output file.\n";
		return -1;
	}
	H5Sclose(src_dataspace);
	H5Tclose(src_datatype);
	H5Dclose(src_dataset);

	herr_t grp_status = H5Gclose(group_id);
	if(grp_status < 0) {
		std::cerr <<  "Error: H5Gclose in output file.\n";
		return -1;
	}
	
	if(src_rad)
		free(src_rad);
	if(src_rad_out)
		delete [] src_rad_out;
	if (nsrcPixels)
		delete [] nsrcPixels;


	//----------------------------------------
	// Closing data files
	std::cout  <<  "Writing done. Closing files...\n";
	herr_t close_status = af_close(src_file);
	if(close_status < 0){
		std::cerr  <<  "Error: closing input data file.\n";
		return -1;
	}

	close_status = af_close(output_file);
	if(close_status < 0){
		std::cerr  <<  "Error: closing output data file.\n";
		return -1;
	}

	return 0;
}
