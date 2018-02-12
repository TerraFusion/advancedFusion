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
#include <vector>

#include <hdf5.h>
#include "io.h"
#include "reproject.h"
#include "AF_InputParmeterFile.h"
#include "AF_debug.h"

#define SUCCEED 0
#define FAILED -1

const hsize_t RESAMPLE_MODIS_DATA_WIDTH = 1354;

void DisplayTimeval (struct timeval *tv)
{
	long milliseconds;

	/* Compute milliseconds from microseconds.  */
	milliseconds = (*tv).tv_usec / 1000;
	/* Print the formatted time, in seconds, followed by a decimal point and the milliseconds.  */
	printf ("%ld.%03ld\n", (*tv).tv_sec, milliseconds);
}

void Usage(int &argc, char *argv[])
{
	std::cout	<< "Usage: \n"
				<< "   " << argv[0] << "  <parameter-input-file>\n";
}

/*###############################################
  Read & Write radiance data as target MODIS
#################################################*/
static int WriteModisRadianceAsTrg(hid_t outputFile, hid_t modisDatatype, hid_t modisFilespace, double* modisData, int modisDataSize, int bandIdx)
{
	int ret = SUCCEED;

	herr_t status;
	hid_t modis_dataset;
	if(bandIdx==0) { // means new
		modis_dataset = H5Dcreate2(outputFile, "/Data_Fields/modis_rad", modisDatatype, modisFilespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(modis_dataset < 0) {
			std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dcreate2 target data in output file.\n";
			return FAILED;
		}
	}
	else {
		modis_dataset = H5Dopen2(outputFile, "/Data_Fields/modis_rad", H5P_DEFAULT);
		if(modis_dataset < 0) {
			std::cerr <<  __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dopen2 target data in output file.\n";
			return FAILED;
		}
	}


	//------------------------------
	// select memspace
	int rankMem=2;
	hsize_t dim2dMem[2];
	hsize_t start2dMem[2];
	hsize_t count2dMem[2];
	dim2dMem[0] = modisDataSize/RESAMPLE_MODIS_DATA_WIDTH; // y
	dim2dMem[1] = RESAMPLE_MODIS_DATA_WIDTH; // x
	hid_t memspace = H5Screate_simple(rankMem, dim2dMem, NULL);

	start2dMem[0] = 0; // y
	start2dMem[1] = 0; // x
	count2dMem[0] = modisDataSize/RESAMPLE_MODIS_DATA_WIDTH; // y
	count2dMem[1] = RESAMPLE_MODIS_DATA_WIDTH; // x

	status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start2dMem, NULL, count2dMem, NULL);

	//------------------------------
	// select filespace
	hsize_t star3dFile[3];
	hsize_t count3dFile[3];
	star3dFile[0] = bandIdx;
	star3dFile[1] = 0; // y
	star3dFile[2] = 0; // x
	count3dFile[0] = 1;
	count3dFile[1] = modisDataSize/RESAMPLE_MODIS_DATA_WIDTH; // y
	count3dFile[2] = RESAMPLE_MODIS_DATA_WIDTH;  // x

	status = H5Sselect_hyperslab(modisFilespace, H5S_SELECT_SET, star3dFile, NULL, count3dFile, NULL);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Sselect_hyperslab for Modis target .\n";
		ret = -1;
		goto done;
	}

	status = H5Dwrite(modis_dataset, H5T_NATIVE_DOUBLE, memspace, modisFilespace, H5P_DEFAULT, modisData);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dwrite for Modis target .\n";
		ret = -1;
		goto done;
	}

done:
	H5Dclose(modis_dataset);
	return ret;
}


int af_ModisAsTrg_GenerateBandsOutputCumulative(hid_t outputFile,hid_t srcFile, int singleDataSize, std::string modisResolution,  std::vector<std::string> &bands)
{

	printf("af_ModisAsTrg_GenerateBandsOutputCumulative\n");

	//----------------------------------
	// prepare group of dataset
	printf("writing data fields\n");
	hid_t group_id = H5Gcreate2(outputFile, "/Data_Fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if(group_id < 0) {
		std::cerr <<  "Error: H5Gcreate2 in output file.\n";
		return FAILED;
	}
	herr_t grp_status = H5Gclose(group_id);
	if(grp_status < 0) {
		std::cerr <<  "Error: H5Gclose in output file.\n";
		return FAILED;
	}

	hid_t modisDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t	status = H5Tset_order(modisDatatype, H5T_ORDER_LE);
	if(status < 0) {
		printf("Error: MODIS write error in H5Tset_order\n");
		return FAILED;
	}


	hsize_t modis_dim[3];
	modis_dim[0] = bands.size();
	modis_dim[1] = singleDataSize/RESAMPLE_MODIS_DATA_WIDTH; // NY;
	modis_dim[2] = RESAMPLE_MODIS_DATA_WIDTH; // NX;
	hid_t modisDataspace = H5Screate_simple(3, modis_dim, NULL);

	int numCells;
	double *modisSingleData;
	std::vector<std::string> singleBandVec;
	for (int i=0; i< bands.size(); i++) {
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> bands[" << i << "]" << bands[i] << "\n";
		#endif
		// insert to only begin
		singleBandVec.insert(singleBandVec.begin(), bands[i]);
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> singleBandVec[0]: " << singleBandVec[0] << "\n";
		#endif
		//Read multi bands (instead Generate single band data for test)
		//GenerateDataDouble2D(singleData, singleDataDims2D, atoi(bands[i].c_str()));
		//DisplayDataDouble2D(singleData, singleDataDims2D);
		//---------------------------------
		// read trg radiance from BF file
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		modisSingleData = get_modis_rad(srcFile, (char*)modisResolution.c_str(), singleBandVec, 1, &numCells);
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> numCells: " << numCells << "\n";
		#endif
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG_TIME> Read target Modis single band data  DONE.");
		#endif

		//---------------------------------
		// write trg radiance to AF file
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		WriteModisRadianceAsTrg(outputFile, modisDatatype, modisDataspace,  modisSingleData, numCells, i);
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG_TIME> Write target Modis single band data  DONE.");
		#endif
		//
		// TODO: buffer is not reused. make memory allocation out of get_modis_rad() to improve performance
		free(modisSingleData);
	}

	H5Tclose(modisDatatype);
	H5Sclose(modisDataspace);
}

int   AF_GenerateTargetRadiancesOutput(AF_InputParmeterFile &inputArgs, hid_t outputFile, int singleDataSize, hid_t srcFile, std::map<std::string, strVec_t> & trgInputMultiVarsMap)
{

	strVec_t multiVarNames;
	std::string instrument = inputArgs.GetTargetInstrument();
	std::string mixType = "COMBINATION"; // Default
	//--------------------------------------------------------
	// Use this for Single multi-value Variable case. MODIS

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> trgInputMultiVarsMap.size(): " << trgInputMultiVarsMap.size() << "\n";
	#endif
	// display strBuf with array index
	if (instrument == "MODIS") {
		multiVarNames = inputArgs.GetMultiVariableNames("MODIS"); // modis_MultiVars;

		if (trgInputMultiVarsMap.size() != 1) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building target input list with MODIS. There must be only one multi-value variable.\n";
			return FAILED;
		}

		std::cout << "Display trgInputMultiVarsMap with array index\n";
		// TODO - decide for loop inside or loop outside
		#if 1 // case for looping inside
		af_ModisAsTrg_GenerateBandsOutputCumulative(outputFile, srcFile, singleDataSize, inputArgs.GetMODIS_Resolution(), trgInputMultiVarsMap[multiVarNames[0]]);
		#else // case for looping from outside
		for(int j = 0; j < trgInputMultiVarsMap[multiVarNames[0]].size(); ++j) {
			std::cout << trgInputMultiVarsMap[multiVarNames[0]][j]	<< ", ";
			AF_Method2_MODIStrg_ReadWriteSingleBandData(outputFile, srcFile, singleDataSize, trgInputMultiVarsMap[multiVarNames[0]][j] /* single Band*/);
		}
		std::cout << std::endl;
		#endif
	}
	//--------------------------------------------------------
	// Use these for two multi-value Variable case. MISR and ASTER
	else if (instrument == "MISR") {
		multiVarNames = inputArgs.GetMultiVariableNames("MISR"); // misr_MultiVars;

		if (trgInputMultiVarsMap.size() != 2) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building target input list with MISR. There must be only two multi-value variables.\n";
			return FAILED;
		}

		std::string misrCamera;
		std::string misrRadiance;
		if(mixType == "COMBINATION") {
			// JK_TODO - Try inside loop first like MODIS
			for(int i=0; i<trgInputMultiVarsMap.size()-1;i++) {
				for(int j=i; j <trgInputMultiVarsMap[multiVarNames[i]].size(); j++) {
					misrCamera = trgInputMultiVarsMap[multiVarNames[i]][j];
					for (int k=0; k<trgInputMultiVarsMap[multiVarNames[i+1]].size(); k++) {
						misrRadiance = trgInputMultiVarsMap[multiVarNames[i+1]][k];
						std::cout << misrCamera << ":" << misrRadiance << "\n";
					}
				}
			}
		}
		// Not used for now. Place holder for later.
		#if 0
		else if(mixType == "PAIR") {
			std::cout << "\nJKDBG> mixType == PAIR\n";
			int numCameras = trgInputMultiVarsMap[multiVarNames[0]].size();
			int numRadiances = trgInputMultiVarsMap[multiVarNames[1]].size();
			std::cout << "JKDBG> var0 num of cameras: " << numCameras << "\n";
			std::cout << "JKDBG> var1 num of radiances: " << numRadiances << "\n";
			int minNumVals = (numCameras < numRadiances) ? numCameras : numRadiances;
			std::cout << "JKDBG> minNumVals: " << minNumVals << "\n";
			for(int i=0; i<minNumVals ;i++) {
				//std::cout << "JKDBG trgInputMultiVarsMap[multiVarNames[0]][" << i << "]" <<  trgInputMultiVarsMap[multiVarNames[0]][i] << "\n";
				//std::cout << "JKDBG trgInputMultiVarsMap[multiVarNames[1]][" << i << "]" <<  trgInputMultiVarsMap[multiVarNames[1]][i] << "\n";
				// misr camera
				misrCamera = trgInputMultiVarsMap[multiVarNames[0]][i];
				// misr radiance
				misrRadiance = trgInputMultiVarsMap[multiVarNames[1]][i];
				std::cout << misrCamera << " : " << misrRadiance <<  "\n";
			}
		}
		#endif
	}
	return SUCCEED;
}




int main(int argc, char *argv[])
{
	int ret = SUCCEED;

	if (argc < 2) {
		Usage(argc, argv);
		return FAILED;
	}

	//----------------------------------
	// parse input parameter from file 
	AF_InputParmeterFile inputArgs;
	inputArgs.headerFileName = argv[1];
	inputArgs.ParseByLine();

	//-------------------------
	// get instrument names
	std::string srcInstrument = inputArgs.GetSourceInstrument();
	std::string trgInstrument = inputArgs.GetTargetInstrument();
	
#if 0 // TEST : multi-value variable map
	//---------------------------------------------------
	// build target instrument multi-value variable map
	std::map<std::string, strVec_t> trgInputMultiVarsMap;
	inputArgs.BuildMultiValueVariableMap(trgInstrument, trgInputMultiVarsMap);
	inputArgs.DBG_displayinputListMap(trgInstrument, trgInputMultiVarsMap, "COMBINATION");

	//---------------------------------------------------
	// build source instrument multi-value variable map 
	std::map<std::string, strVec_t> srcInputMultiVarsMap;
	inputArgs.BuildMultiValueVariableMap(srcInstrument, srcInputMultiVarsMap);
	inputArgs.DBG_displayinputListMap(srcInstrument, srcInputMultiVarsMap, "COMBINATION");
	exit(1);
#endif

	/* ---------------------------------------------------
	 * Create output file
	 */
	std::string outputFile = inputArgs.GetOuputFilePath();
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL main> outputFile: " << outputFile << std::endl;
	#endif

	hid_t output_file = H5Fcreate(outputFile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


	/* ---------------------------------------------------
	 * Handle input BF data file
	 */
	std::string inputDataPath = inputArgs.GetInputBFdataPath();
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL main> inputDataPath: " << inputDataPath << std::endl;
	#endif
	hid_t src_file;
	if(0 > (src_file = af_open((char*)inputDataPath.c_str()))) {
		std::cerr << "Error: File not found - " << inputDataPath << std::endl;
		exit(1);
	}
	
	/* ---------------------------------------------------
	 * Handle Source instrument latitude and longitude
	 */
	std::cout << "\nGetting source instrument latitude & longitude data...\n";
	std::string misr_resolution = inputArgs.GetMISR_Resolution();
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL main> misr_resolution: " << misr_resolution << std::endl;
	#endif
	int srcCellNum;
	double* srcLatitude = NULL;
	double* srcLongitude = NULL;
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	srcLatitude = get_misr_lat(src_file, (char*)misr_resolution.c_str(), &srcCellNum);
	srcLongitude = get_misr_long(src_file, (char*)misr_resolution.c_str(), &srcCellNum);
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> get misr lat/long DONE.");
	#endif

	/* ---------------------------------------------------
	 * Handle Target instrument latitude and longitude
	 */
	std::string modis_resolution = inputArgs.GetMODIS_Resolution();
	std::cout << "\nGetting target instrument latitude & longitude data...\n";
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL main> modis_resolution: " << modis_resolution << std::endl;
	#endif
	int trgCellNum;
	double* targetLatitude = NULL;
	double* targetLongitude = NULL;
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	targetLatitude = get_modis_lat(src_file, (char*)modis_resolution.c_str(), &trgCellNum);
	targetLongitude = get_modis_long(src_file, (char*)modis_resolution.c_str(), &trgCellNum);
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> get modis lat/long DONE.");
	#endif
	
	std::cout << "\nWriting target geolocation data...\n";
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	int lat_status =  af_write_mm_geo(output_file, 0, targetLatitude, trgCellNum);
	if(lat_status < 0) {
		std::cerr << "Error: writing latitude geolocation.\n";
		return FAILED;
	}
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> write geo lattitude data DONE.");
	#endif

	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	int long_status = af_write_mm_geo(output_file, 1, targetLongitude, trgCellNum);
	if(long_status < 0) {
		std::cerr << "Error: writing longitude geolocation.\n";
		return FAILED;
	}
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> write geo longitude data DONE.");
	#endif

	int * targetNNsrcID = NULL;
	targetNNsrcID = new int [trgCellNum];
	
	std::cout <<  "\nRunning nearest neighbor block index method... \n";
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	nearestNeighborBlockIndex(&srcLatitude, &srcLongitude, srcCellNum, targetLatitude, targetLongitude, targetNNsrcID, NULL, trgCellNum, 1000);
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> nearestNeighborBlockIndex DONE.");
	#endif

	if(srcLatitude)
		free(srcLatitude);
	if(srcLongitude)
		free(srcLongitude);
	if(targetLatitude)
		free(targetLatitude);
	if(targetLongitude)
		free(targetLongitude);
	
	/*------------------------------------------------------
	 * Generate target instrument radiance to output file
	 */
	// build target instrument multi-value variable map
	std::map<std::string, strVec_t> trgInputMultiVarsMap;
	inputArgs.BuildMultiValueVariableMap(trgInstrument, trgInputMultiVarsMap);
	// write target instrument radiances to output file
	ret = AF_GenerateTargetRadiancesOutput(inputArgs, output_file, trgCellNum, src_file, trgInputMultiVarsMap);
	if (ret < 0) {
		return FAILED;
	}


	//=========================================
	// START: Source radiance data write

	//-------------------------------------------------
	// Read Source instrument radiance from BF file
	std::cout << "\nGetting source instrument radiance data...\n";
	double* src_rad=NULL;
	std::vector<std::string> misr_radiance = inputArgs.GetMISR_Radiance();
	std::vector<std::string> misr_cameraAngles = inputArgs.GetMISR_CameraAngles();
	#if DEBUG_TOOL
	for(int i = 0; i < misr_radiance.size(); i++)
		std::cout << "DBG_TOOL main> misr_radiance[" << i << "]: " << misr_radiance[i] << std::endl;
	for(int i = 0; i < misr_cameraAngles.size(); i++)
		std::cout << "DBG_TOOL main> misr_cameraAngles[" << i << "]:" << misr_cameraAngles[i] << std::endl;
	#endif

	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	// TODO: [0] for only single value, change to multiple values
	src_rad = get_misr_rad(src_file, (char*)misr_cameraAngles[0].c_str(), (char*)misr_resolution.c_str(), (char*)misr_radiance[0].c_str(), &srcCellNum);
	if (src_rad == NULL) {
		std::cerr	<< "Please verify these MISR input values. \n" 
					<< "	   - Camera: " << misr_cameraAngles[0] << "\n" 
					<< "	   - Resolution: " << misr_resolution << "\n"
					<< "	   - Radiance: " << misr_radiance[0] << "\n" 
					<< std::endl;
		return FAILED;
	}
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> get_misr_rad DONE.");
	#endif
	
	
	//-------------------------------------------------
	// handle resample method
	double* src_rad_out = new double [trgCellNum];
	int new_src_size = trgCellNum;
	int * nsrcPixels = NULL;
	//Interpolating
	std::string resampleMethod =  inputArgs.GetResampleMethod();
	std::cout << "\nInterpolating using '" << resampleMethod << "' method...\n";
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "nnInterpolate")) {
		nnInterpolate(src_rad, src_rad_out, targetNNsrcID, trgCellNum);
	}
	else if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "summaryInterpolate")) {
		nsrcPixels = new int [trgCellNum];
		summaryInterpolate(src_rad, targetNNsrcID, srcCellNum, src_rad_out, nsrcPixels, trgCellNum);
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL> No nodata values: \n";
		for(int i = 0; i < trgCellNum; i++) {
			if(nsrcPixels[i] > 0) {
				printf("%d,\t%lf\n", nsrcPixels[i], src_rad_out[i]);
			}
		}
		#endif
	}
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG> nnInterpolate  DONE.");
	#endif

	if (targetNNsrcID)
		delete [] targetNNsrcID;


	//--------------------------------------------
	// Write reampled source instrument data next
	std::cout << "\nWriting source instrument '" << srcInstrument << "' data...\n";
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	hsize_t src_dim[2];
	src_dim[0] = (new_src_size) / RESAMPLE_MODIS_DATA_WIDTH;
	src_dim[1] = RESAMPLE_MODIS_DATA_WIDTH;
	hid_t src_dataspace = H5Screate_simple(2, src_dim, NULL);
	hid_t src_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t src_status = H5Tset_order(src_datatype, H5T_ORDER_LE);  
	hid_t src_dataset = H5Dcreate2(output_file, "/Data_Fields/misr_out", src_datatype, src_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if(src_dataset < 0){
		std::cerr <<  "Error: H5Dcreate2 source data in output file.\n";
		return FAILED;
	}
	src_status = H5Dwrite(src_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_rad_out);
	if(src_status < 0){
		std::cerr  <<  "Error: H5Dwrite source data in output file.\n";
		return FAILED;
	}
	H5Sclose(src_dataspace);
	H5Tclose(src_datatype);
	H5Dclose(src_dataset);

	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> Write Source data DONE.");
	#endif
	
	if(src_rad)
		free(src_rad);
	if(src_rad_out)
		delete [] src_rad_out;
	if (nsrcPixels)
		delete [] nsrcPixels;

	// END: Source
	//=========================================

	//----------------------------------------
	// Closing data files
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	std::cout  <<  "\nWriting done. Closing files...\n";
	herr_t close_status = af_close(src_file);
	if(close_status < 0){
		std::cerr  <<  "Error: closing input data file.\n";
		return FAILED;
	}

	close_status = af_close(output_file);
	if(close_status < 0){
		std::cerr  <<  "Error: closing output data file.\n";
		return FAILED;
	}
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> Closing files DONE.");
	#endif

	return 0;
}
