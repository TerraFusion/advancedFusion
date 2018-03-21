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
 *	instrument data paired with geo location data. 
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
#include "misrutil.h"

#define SUCCEED 0
#define FAILED -1


/*------------------------
 * hdf5 related names
 */
const std::string SRC_DATA_GROUP = "/Source/Data_Fields";
const std::string TRG_DATA_GROUP = "/Target/Data_Fields";
const std::string MODIS_RADIANCE_DSET = "MODIS_Radiance";
const std::string MISR_RADIANCE_DSET = "MISR_Radiance";



void DisplayTimeval (struct timeval *tv)
{
	long milliseconds;

	/* Compute milliseconds from microseconds.	*/
	milliseconds = (*tv).tv_usec / 1000;
	/* Print the formatted time, in seconds, followed by a decimal point and the milliseconds.	*/
	printf ("%ld.%03ld\n", (*tv).tv_sec, milliseconds);
}

void Usage(int &argc, char *argv[])
{
	std::cout	<< "Usage: \n"
				<< "   " << argv[0] << "  <parameter-input-file>\n";
}

/*##############################################################
 *
 * Util functions
 *
 *#############################################################*/

int AF_GetGeolocationDataFromInstrument(std::string instrument, AF_InputParmeterFile &inputArgs, hid_t inputFile, double **latitude /*OUT*/, double **longitude /*OUT*/, int &cellNum /*OUT*/)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	if (instrument == "MODIS" ) {
		std::string resolution = inputArgs.GetMODIS_Resolution();
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> Modis resolution: " << resolution << "\n";
		#endif
		*latitude = get_modis_lat(inputFile, (char*) resolution.c_str(), &cellNum);
		if (latitude == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get MODIS latitude.\n";
			return FAILED;
		}
		*longitude = get_modis_long(inputFile, (char*) resolution.c_str(), &cellNum);
		if (longitude == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get MODIS longitude.\n";
			return FAILED;
		}
	}
	else if (instrument == "MISR") {
		std::string resolution = inputArgs.GetMISR_Resolution();
		*latitude = get_misr_lat(inputFile, (char*) resolution.c_str(), &cellNum);
		if (latitude == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get MISR latitude.\n";
			return FAILED;
		}
		*longitude = get_misr_long(inputFile, (char*) resolution.c_str(), &cellNum);
		if (longitude == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get MISR longitude.\n";
			return FAILED;
		}

		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> Misr resolution: " << resolution << ", cellNum: " << cellNum << "\n";
		#endif
	}
	#if 0 // TODO LATER - place holder
	else if (instrument == "ASTER" ) {
		//Get ASTER input parameters EX: "TIR", "ImageData10"
		//*latitude = get_ast_lat(inputFile, "TIR", "ImageData10", &cellNum);
		//*longitude = get_ast_long(inputFile, "TIR", "ImageData10", &cellNum);
	}
	#endif
	else {
		std::cerr << __FUNCTION__ << "> Error: invalid instrument - " << instrument << "\n";
		return FAILED;
	}

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
	
	return SUCCEED;
}



/*=====================================
 * Get output width of an instrument
 *
 * RETURN:
 *  0 : SUCCEED
 * -1 : FAILED
 *
 * OUT parameters:
 * - crossTrackWidth : return width value.
 * - alongTrackHeight : return height value.
 *						This is only valid for Misr target shift case.
 *						If 0, caller should not use this.
 */
int af_GetWidthAndHeightForOutputDataSize(std::string instrument, AF_InputParmeterFile &inputArgs, int &crossTrackWidth /*OUT*/, int &alongTrackHeight /*OUT*/)
{
	int ret = 0;
	crossTrackWidth = 0;
	alongTrackHeight = 0;
	std::string misrShift = inputArgs.GetMISR_Shift();
	std::string trgInstrument = inputArgs.GetTargetInstrument();

	// MISR shift case. Shift always if misr is target and Shift is On for all source instrument and target geolocation data
	if(misrShift == "ON" && trgInstrument == "MISR") {
		std::string resolution = inputArgs.GetMISR_Resolution();
		getMISRFinalImageSize(&alongTrackHeight, &crossTrackWidth, (resolution=="L") ? 0 : 1);
	}
	else if (instrument == "MODIS") {
		std::string resolution = inputArgs.GetMODIS_Resolution();
		if (resolution == "_1KM") {
			crossTrackWidth = 1354;
		}
		else if (resolution == "_500m") {
			crossTrackWidth = 2708;  // 1354 * 2
		}
		else if (resolution == "_250m") {
			crossTrackWidth = 5416;  // 1354 * 4
		}
	}
	else if (instrument == "MISR") {
		std::string resolution = inputArgs.GetMISR_Resolution();
		if (resolution == "L") {  // 1.1KM
			crossTrackWidth = 512;
		}
		else if (resolution == "H") {  // 275M
			crossTrackWidth = 2048;
		}
	}
	else {
		return -1;  // fail
	}

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> misrShift: " << misrShift << ", instrument: " << instrument << ", crossTrackWidth: " << crossTrackWidth << ", alongTrackHeight: " << alongTrackHeight << "\n";
	#endif

	return 0;
}

/*##############################################################
 *
 * Generate Target instrument radiance to output file
 *
 *#############################################################*/

/* 
 * function prototype for AF_GenerateTargetRadiancesOutput()
 */
int af_GenerateOutputCumulative_ModisAsTrg(AF_InputParmeterFile &inputArgs, hid_t outputFile,hid_t srcFile, int trgCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap);
int af_GenerateOutputCumulative_MisrAsTrg(AF_InputParmeterFile &inputArgs, hid_t outputFile,hid_t srcFile, int trgCellNumOri, std::map<std::string, strVec_t> &inputMultiVarsMap);


/*=================================================
 * Generate target instrument randiance output.
 * Initial function.
 */
int   AF_GenerateTargetRadiancesOutput(AF_InputParmeterFile &inputArgs, hid_t outputFile, int trgCellNum, hid_t srcFile, std::map<std::string, strVec_t> & trgInputMultiVarsMap)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret;

	//----------------------------------
	// prepare group of dataset
	printf("Creating target group...\n");
	hid_t gcpl_id = H5Pcreate (H5P_LINK_CREATE);
	if(gcpl_id < 0) {
		std::cerr << __FUNCTION__ <<  "> Error: H5Pcreate.\n";
		return FAILED;
	}
	herr_t status = H5Pset_create_intermediate_group (gcpl_id, 1);
	if(status < 0) {
		std::cerr << __FUNCTION__ <<  "> Error: H5Pset_create_intermediate_group.\n";
		return FAILED;
	}
	hid_t group_id = H5Gcreate2(outputFile, TRG_DATA_GROUP.c_str(), gcpl_id, H5P_DEFAULT, H5P_DEFAULT);
	if(group_id < 0) {
		std::cerr << __FUNCTION__ <<  "> Error: H5Gcreate2 in output file.\n";
		return FAILED;
	}
	herr_t grp_status = H5Gclose(group_id);
	if(grp_status < 0) {
		std::cerr << __FUNCTION__ <<  "> Error: H5Gclose in output file.\n";
		return FAILED;
	}

	strVec_t multiVarNames;
	std::string instrument = inputArgs.GetTargetInstrument();
	// std::string mixType = "COMBINATION"; // Default
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

		// TODO_LATER - decide for loop inside or loop outside
		#if 1 // case for looping inside
		ret = af_GenerateOutputCumulative_ModisAsTrg(inputArgs, outputFile, srcFile, trgCellNum, trgInputMultiVarsMap);
		if(ret == FAILED) {
			std::cerr << __FUNCTION__ <<  "> Error: Generating MODIS target output.\n";
			return FAILED;
		}
		#else // case for looping from outside
		for(int j = 0; j < trgInputMultiVarsMap[multiVarNames[0]].size(); ++j) {
			std::cout << trgInputMultiVarsMap[multiVarNames[0]][j]	<< ", ";
			AF_Method2_MODIStrg_ReadWriteSingleBandData(outputFile, srcFile, trgCellNum, trgInputMultiVarsMap[multiVarNames[0]][j] /* single Band*/);
		}
		std::cout << std::endl;
		#endif
	}
	//--------------------------------------------------------
	// Use these for two multi-value Variable case. MISR and ASTER
	else if (instrument == "MISR") {

		if (trgInputMultiVarsMap.size() != 2) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building target input list with MISR. There must be only two multi-value variables.\n";
			return FAILED;
		}

		ret = af_GenerateOutputCumulative_MisrAsTrg(inputArgs, outputFile, srcFile, trgCellNum, trgInputMultiVarsMap);
		if(ret == FAILED) {
			std::cerr << __FUNCTION__ <<  "> Error: Generating MISR target output.\n";
			return FAILED;
		}
	}

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
	return SUCCEED;
}



/*============================================================
 * Target as MODIS functions for the output
 */
static int af_WriteSingleRadiance_ModisAsTrg(hid_t outputFile, hid_t modisDatatype, hid_t modisFilespace, double* modisData, int modisDataSize, int outputWidth, int bandIdx)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;

	herr_t status;
	hid_t modis_dataset;
	std::string dsetPath = TRG_DATA_GROUP + "/" + MODIS_RADIANCE_DSET;

	if(bandIdx==0) { // means new
		modis_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), modisDatatype, modisFilespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(modis_dataset < 0) {
			std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dcreate2 target data in output file.\n";
			return FAILED;
		}
	}
	else {
		modis_dataset = H5Dopen2(outputFile, dsetPath.c_str(), H5P_DEFAULT);
		if(modis_dataset < 0) {
			std::cerr <<  __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dopen2 target data in output file.\n";
			return FAILED;
		}
	}


	//------------------------------
	// select memspace
	const int ranksMem=2;
	hsize_t dim2dMem[ranksMem];
	hsize_t start2dMem[ranksMem];
	hsize_t count2dMem[ranksMem];
	dim2dMem[0] = modisDataSize/outputWidth; // y
	dim2dMem[1] = outputWidth; // x
	hid_t memspace = H5Screate_simple(ranksMem, dim2dMem, NULL);

	start2dMem[0] = 0; // y
	start2dMem[1] = 0; // x
	count2dMem[0] = modisDataSize/outputWidth; // y
	count2dMem[1] = outputWidth; // x

	status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start2dMem, NULL, count2dMem, NULL);

	//------------------------------
	// select filespace
	const int ranksFile=3; // [bands][y][x]
	hsize_t star3dFile[ranksFile];
	hsize_t count3dFile[ranksFile];
	star3dFile[0] = bandIdx;
	star3dFile[1] = 0; // y
	star3dFile[2] = 0; // x
	count3dFile[0] = 1;
	count3dFile[1] = modisDataSize/outputWidth; // y
	count3dFile[2] = outputWidth;  // x

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

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
	return ret;
}


int af_GenerateOutputCumulative_ModisAsTrg(AF_InputParmeterFile &inputArgs, hid_t outputFile,hid_t srcFile, int trgCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	std::string modisResolution = inputArgs.GetMODIS_Resolution();
	// get multi-value variable for modis
	strVec_t bands = inputMultiVarsMap[MODIS_BANDS];

	// data type
	hid_t modisDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t	status = H5Tset_order(modisDatatype, H5T_ORDER_LE);
	if(status < 0) {
		std::cerr << __FUNCTION__ <<  "> Error: MODIS write error in H5Tset_order.\n";
		return FAILED;
	}


	/*----------------------------
	 * create data space
	 */
	/* handle different data width and height. no need to care for misr-trg shift case
	 */
	int ret;
	int widthShifted = 0;
	int heightShifted = 0;
	ret = af_GetWidthAndHeightForOutputDataSize(inputArgs.GetTargetInstrument() /*target base output*/, inputArgs, widthShifted, heightShifted);
	if( (ret < 0) || (heightShifted > 0) ) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error in af_GetWidthAndHeightForOutputDataSize() \n";
		return FAILED;
	}
    int targetOutputWidth = widthShifted;
	const int rankSpace=3;  // [bands][y][x]
	hsize_t modis_dim[rankSpace];
	modis_dim[0] = bands.size();
	modis_dim[1] = trgCellNum/targetOutputWidth; // NY;
	modis_dim[2] = targetOutputWidth; // NX;
	hid_t modisDataspace = H5Screate_simple(rankSpace, modis_dim, NULL);

	int numCells;
	double *modisSingleData;
	std::vector<std::string> singleBandVec;
	for (int i=0; i< bands.size(); i++) {
		std::cout << "Processing MODIS band: " << bands[i] << "\n";
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> bands[" << i << "]: " << bands[i] << "\n";
		#endif
		// insert to only begin
		singleBandVec.insert(singleBandVec.begin(), bands[i]);
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> singleBandVec[0]: " << singleBandVec[0] << "\n";
		#endif
		#if 0 // TOOL_TEST remove later
		//Read multi bands (instead Generate single band data for test)
		//GenerateDataDouble2D(singleData, singleDataDims2D, atoi(bands[i].c_str()));
		//DisplayDataDouble2D(singleData, singleDataDims2D);
		#endif
		//---------------------------------
		// read trg radiance from BF file
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		modisSingleData = get_modis_rad(srcFile, (char*)modisResolution.c_str(), singleBandVec, 1, &numCells);
		if (modisSingleData == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get MODIS radiance.\n";
			return FAILED;
		}
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
		af_WriteSingleRadiance_ModisAsTrg(outputFile, modisDatatype, modisDataspace,  modisSingleData, numCells, targetOutputWidth, i);
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG_TIME> Write target Modis single band data  DONE.");
		#endif
		//
		// TODO: buffer is not reused. make memory allocation out of get_modis_rad() to improve performance
		free(modisSingleData);
	}

	H5Tclose(modisDatatype);
	H5Sclose(modisDataspace);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif

	return SUCCEED;
}



/*==================================================================
 * Target as MISR functions for the output
 */
static int af_WriteSingleRadiance_MisrAsTrg(hid_t outputFile, hid_t misrDatatype, hid_t misrFilespace, double* misrData, int misrDataSize, int outputWidth, int cameraIdx, int radianceIdx)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;

	herr_t status;
	hid_t misr_dataset;
	std::string dsetPath = TRG_DATA_GROUP + "/" + MISR_RADIANCE_DSET;

	if( (cameraIdx + radianceIdx) ==0 ) { // means new
		misr_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), misrDatatype, misrFilespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(misr_dataset < 0) {
			std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dcreate2 target data in output file.\n";
			return FAILED;
		}
	}
	else {
		misr_dataset = H5Dopen2(outputFile, dsetPath.c_str(), H5P_DEFAULT);
		if(misr_dataset < 0) {
			std::cerr <<  __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dopen2 target data in output file.\n";
			return FAILED;
		}
	}


	//------------------------------
	// select memspace
	const int ranksMem=2;
	hsize_t dim2dMem[ranksMem];
	hsize_t start2dMem[ranksMem];
	hsize_t count2dMem[ranksMem];
	dim2dMem[0] = misrDataSize/outputWidth; // y
	dim2dMem[1] = outputWidth; // x
	hid_t memspace = H5Screate_simple(ranksMem, dim2dMem, NULL);

	start2dMem[0] = 0; // y
	start2dMem[1] = 0; // x
	count2dMem[0] = misrDataSize/outputWidth; // y
	count2dMem[1] = outputWidth; // x

	status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start2dMem, NULL, count2dMem, NULL);

	//------------------------------
	// select filespace
	const int ranksFile=4; // [cameras][radiances][y][x]
	hsize_t star3dFile[ranksFile];
	hsize_t count3dFile[ranksFile];
	star3dFile[0] = cameraIdx;
	star3dFile[1] = radianceIdx;
	star3dFile[2] = 0; // y
	star3dFile[3] = 0; // x
	count3dFile[0] = 1;
	count3dFile[1] = 1;
	count3dFile[2] = misrDataSize/outputWidth; // y
	count3dFile[3] = outputWidth;  // x

	status = H5Sselect_hyperslab(misrFilespace, H5S_SELECT_SET, star3dFile, NULL, count3dFile, NULL);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Sselect_hyperslab for Misr target .\n";
		ret = -1;
		goto done;
	}

	status = H5Dwrite(misr_dataset, H5T_NATIVE_DOUBLE, memspace, misrFilespace, H5P_DEFAULT, misrData);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dwrite for Misr target .\n";
		ret = -1;
		goto done;
	}

done:
	H5Dclose(misr_dataset);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
	return ret;
}


int af_GenerateOutputCumulative_MisrAsTrg(AF_InputParmeterFile &inputArgs, hid_t outputFile,hid_t srcFile, int trgCellNumOri, std::map<std::string, strVec_t> &inputMultiVarsMap)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif
	int ret;

	// data type
	hid_t misrDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t	status = H5Tset_order(misrDatatype, H5T_ORDER_LE);
	if(status < 0) {
		std::cerr << __FUNCTION__ <<  "> Error: MISR write error in H5Tset_order.\n";
		return FAILED;
	}


	// get misr input variables
	std::string misrResolution = inputArgs.GetMISR_Resolution();
	// get multi-value variable for misr
	strVec_t cameras = inputMultiVarsMap[MISR_CAMERA_ANGLE];
	strVec_t radiances = inputMultiVarsMap[MISR_RADIANCE];


	/*----------------------------
	 * create data space
	 */
	int widthShifted = 0;
	int heightShifted = 0;
	int targetOutputWidth = 0;
	int trgCellNum = 0;

	/*
	 * get targetOutputWidth and trgCellNum
	 */
	ret = af_GetWidthAndHeightForOutputDataSize("MISR", inputArgs, widthShifted, heightShifted);
	if(ret < 0) {
		return FAILED;
	}
    targetOutputWidth = widthShifted;
	if(heightShifted != 0) {  // means misr-trg shif
		trgCellNum = widthShifted * heightShifted;
	}
	else { // no shift
		trgCellNum = trgCellNumOri;
	}

	const int rankSpace=4;  // [cameras][radiances][y][x]
	hsize_t misr_dim[rankSpace];
	misr_dim[0] = cameras.size();
	misr_dim[1] = radiances.size();
	misr_dim[2] = trgCellNum/targetOutputWidth; // NY;
	misr_dim[3] = targetOutputWidth; // NX;
	hid_t misrDataspace = H5Screate_simple(rankSpace, misr_dim, NULL);

	std::string misrShift = inputArgs.GetMISR_Shift();
	int numCells;
	double *misrSingleDataPtr = NULL;
	double *misrSingleData = NULL;
	std::string misrCamera;
	std::string misrRadiance;

	/*
	 * loop multi-value variables in combination to generate radiance output
	 */
	for(int i=0; i <cameras.size(); i++) {
		misrCamera = cameras[i];
		for (int j=0; j< radiances.size(); j++) {
			misrRadiance = radiances[j];
			std::cout << "Processing MISR camera: " << misrCamera << ", radiance: " << misrRadiance << "\n";
			//---------------------------------
			// read trg radiance from BF file
			#if DEBUG_ELAPSE_TIME
			StartElapseTime();
			#endif
			misrSingleData = get_misr_rad(srcFile, (char*)misrCamera.c_str(), (char*)misrResolution.c_str(),(char*) misrRadiance.c_str(), &numCells);
			if (misrSingleData == NULL) {
				std::cerr << __FUNCTION__ <<  "> Error: failed to get MISR radiance.\n";
				return FAILED;
			}
			#if DEBUG_ELAPSE_TIME
			StopElapseTimeAndShow("DBG_TIME> Read target Misr single band data  DONE.");
			#endif

			#if DEBUG_TOOL
			std::cout << "DBG_TOOL " << __FUNCTION__ << "> numCells: " << numCells << "\n";
			#endif
			//-----------------------------------------------------------------------
			// check if need to shift by MISR_TARGET_BLOCKUNSTACK==ON & MISR target case for writing
			double * misrSingleDataShifted = NULL;
			if(misrShift == "ON") {
				std::cout << "\nTarget MISR radiance block unstacking...\n";
				#if DEBUG_ELAPSE_TIME
				StartElapseTime();
				#endif
				misrSingleDataShifted = (double *) malloc(sizeof(double) * widthShifted * heightShifted);
				MISRBlockOffset(misrSingleData, misrSingleDataShifted, (misrResolution == "L") ? 0 : 1);
				#if DEBUG_ELAPSE_TIME
				StopElapseTimeAndShow("DBG_TIME> target MISR radiance block unstack DONE.");
				#endif

				misrSingleDataPtr = misrSingleDataShifted;
				numCells = widthShifted * heightShifted;
			}
			else { // OFF
				misrSingleDataPtr = misrSingleData;
			}

			//---------------------------------
			// write trg radiance to Output file
			#if DEBUG_ELAPSE_TIME
			StartElapseTime();
			#endif
			af_WriteSingleRadiance_MisrAsTrg(outputFile, misrDatatype, misrDataspace,  misrSingleDataPtr, numCells, targetOutputWidth, i, j);
			#if DEBUG_ELAPSE_TIME
			StopElapseTimeAndShow("DBG_TIME> Write target Misr single band data  DONE.");
			#endif

			// TODO: buffer is not reused. make memory allocation out of get_misr_rad() to improve performance
			if (misrSingleData)
				free(misrSingleData);

			if (misrSingleDataShifted)
				free(misrSingleDataShifted);
		}
	}

	H5Tclose(misrDatatype);
	H5Sclose(misrDataspace);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif

	return SUCCEED;
}



/*##############################################################
 *
 * Generate Source instrument radiance to output file
 *
 *#############################################################*/

/* 
 * function prototype for AF_GenerateTargetRadiancesOutput()
 */
int af_GenerateOutputCumulative_ModisAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNumNoShift, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap);
int af_GenerateOutputCumulative_MisrAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNum, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap);


/*===================================================================
 * Generate source instrument randiance output.
 * Initial function.
 */
int   AF_GenerateSourceRadiancesOutput(AF_InputParmeterFile &inputArgs, hid_t outputFile, int * targetNNsrcID, int trgCellNum, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> & srcInputMultiVarsMap)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif
	int ret = SUCCEED;

	//----------------------------------
	// prepare group of dataset
	printf("writing data fields\n");
	hid_t gcpl_id = H5Pcreate (H5P_LINK_CREATE);
	if(gcpl_id < 0) {
		std::cerr << __FUNCTION__ <<  "> Error: H5Pcreate.\n";
		return FAILED;
	}
	herr_t status = H5Pset_create_intermediate_group (gcpl_id, 1);
	if(status < 0) {
		std::cerr << __FUNCTION__ <<  "> Error: H5Pset_create_intermediate_group.\n";
		return FAILED;
	}
	hid_t group_id = H5Gcreate2(outputFile, SRC_DATA_GROUP.c_str(), gcpl_id, H5P_DEFAULT, H5P_DEFAULT);
	if(group_id < 0) {
		std::cerr << __FUNCTION__ <<  "> Error: H5Gcreate2 in output file.\n";
		return FAILED;
	}
	herr_t grp_status = H5Gclose(group_id);
	if(grp_status < 0) {
		std::cerr << __FUNCTION__ <<  "> Error: H5Gclose in output file.\n";
		return FAILED;
	}

	strVec_t multiVarNames;
	std::string instrument = inputArgs.GetSourceInstrument();
	//std::string mixType = "COMBINATION"; // Default
	//--------------------------------------------------------
	// Use this for Single multi-value Variable case. MODIS

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> srcInputMultiVarsMap.size(): " << srcInputMultiVarsMap.size() << "\n";
	#endif

	if (instrument == "MODIS") {
		if (srcInputMultiVarsMap.size() != 1) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building target input list with MODIS. There must be only one multi-value variable.\n";
			return FAILED;
		}

		ret = af_GenerateOutputCumulative_ModisAsSrc(inputArgs, outputFile, targetNNsrcID,  trgCellNum, srcFile, srcCellNum, srcInputMultiVarsMap);
		if (ret == FAILED) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> failed generating output for MISR.\n";
			ret = FAILED;
			goto done;
		}
	}
	//--------------------------------------------------------
	// Use these for two multi-value Variable case. MISR and ASTER
	else if (instrument == "MISR") {

		if (srcInputMultiVarsMap.size() != 2) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building target input list with MISR. There must be only two multi-value variables.\n";
			ret = FAILED;
			goto done;
		}

		ret = af_GenerateOutputCumulative_MisrAsSrc(inputArgs, outputFile, targetNNsrcID,  trgCellNum, srcFile, srcCellNum, srcInputMultiVarsMap);
		if (ret == FAILED) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> failed generating output for MISR.\n";
			ret = FAILED;
			goto done;
		}
	}

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
done:
	return ret;
}


/*======================================================
 * Source as MODIS functions for the output
 */
static int af_WriteSingleRadiance_ModisAsSrc(hid_t outputFile, hid_t modisDatatype, hid_t modisFilespace, double* processedData, int trgCellNum, int outputWidth, int bandIdx)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;
	herr_t status;
	hid_t modis_dataset;
	std::string dsetPath = SRC_DATA_GROUP + "/" + MODIS_RADIANCE_DSET;

	if(bandIdx==0) { // means new
		modis_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), modisDatatype, modisFilespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(modis_dataset < 0) {
			std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dcreate2 target data in output file.\n";
			return FAILED;
		}
	}
	else {
		modis_dataset = H5Dopen2(outputFile, dsetPath.c_str(), H5P_DEFAULT);
		if(modis_dataset < 0) {
			std::cerr <<  __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dopen2 target data in output file.\n";
			return FAILED;
		}
	}


	//------------------------------
	// select memspace
	int ranksMem=2;
	hsize_t dim2dMem[ranksMem];
	hsize_t start2dMem[ranksMem];
	hsize_t count2dMem[ranksMem];
	dim2dMem[0] = trgCellNum/outputWidth; // y
	dim2dMem[1] = outputWidth; // x
	hid_t memspace = H5Screate_simple(ranksMem, dim2dMem, NULL);

	start2dMem[0] = 0; // y
	start2dMem[1] = 0; // x
	count2dMem[0] = trgCellNum/outputWidth; // y
	count2dMem[1] = outputWidth; // x

	status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start2dMem, NULL, count2dMem, NULL);

	//------------------------------
	// select filespace for modis
	const int ranksFile=3; // [bands][y][x]
	hsize_t startFile[ranksFile];
	hsize_t countFile[ranksFile];
	startFile[0] = bandIdx;
	startFile[1] = 0; // y
	startFile[2] = 0; // x
	countFile[0] = 1;
	countFile[1] = trgCellNum/outputWidth; // y
	countFile[2] = outputWidth;  // x

	status = H5Sselect_hyperslab(modisFilespace, H5S_SELECT_SET, startFile, NULL, countFile, NULL);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Sselect_hyperslab for Modis target .\n";
		ret = -1;
		goto done;
	}

	status = H5Dwrite(modis_dataset, H5T_NATIVE_DOUBLE, memspace, modisFilespace, H5P_DEFAULT, processedData);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dwrite for Modis target .\n";
		ret = -1;
		goto done;
	}

done:
	H5Dclose(modis_dataset);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
	return ret;
}



int af_GenerateOutputCumulative_ModisAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNumNoShift, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;

	// strVec_t multiVarNames = inputArgs.GetMultiVariableNames("MODIS"); // modis_MultiVars;
	std::string modisResolution = inputArgs.GetMODIS_Resolution();

	// two multi-value variables are expected as this point
	strVec_t bands = inputMultiVarsMap[MODIS_BANDS];

	// data type
	hid_t modisDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t	status = H5Tset_order(modisDatatype, H5T_ORDER_LE);
	if(status < 0) {
		printf("Error: MODIS write error in H5Tset_order\n");
		return FAILED;
	}


	/*------------------------------------------
	 * create space for modis
	 */
	const int rankSpace=3;  // [bands][y][x]
	hsize_t modisDims[rankSpace];
	/* handle different data width and height, also misr-trg shift case
	 */
	int widthShifted;
	int heightShifted;
	int trgCellNum;
	ret = af_GetWidthAndHeightForOutputDataSize(inputArgs.GetTargetInstrument() /*target base output */, inputArgs, widthShifted, heightShifted);
	if(ret < 0) {
		return FAILED;
	}
    int srcOutputWidth = widthShifted;
	if(inputArgs.GetMISR_Shift() == "ON" && inputArgs.GetTargetInstrument() == "MISR") {
		trgCellNum = widthShifted * heightShifted;
	}
	else {
		trgCellNum = trgCellNumNoShift;
	}
	modisDims[0] = bands.size();
	modisDims[1] = trgCellNum/srcOutputWidth; // NY;
	modisDims[2] = srcOutputWidth; // NX;
	hid_t modisDataspace = H5Screate_simple(rankSpace, modisDims, NULL);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> trgCellNum: " << trgCellNum << ", srcCellNum: " << srcCellNum << "\n";
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> srcOutputWidth: " << srcOutputWidth <<  "\n";
	#endif
	int numCells;
	double *modisSingleData=NULL;

	std::vector<std::string> singleBandVec;
	//-----------------------------------------------------------------
	// TODO: improve by preparing these memory allocation out of loop
	// srcProcessedData , nsrcPixels
	// modisSingleData
	double * srcProcessedData = NULL;
	int * nsrcPixels = NULL;
	// Note: This is Combination case only
	for (int i=0; i< bands.size(); i++) {
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> bands[" << i << "]" << bands[i] << "\n";
		#endif
		// insert to only begin
		singleBandVec.insert(singleBandVec.begin(), bands[i]);

		//---------------------------------
		// read src band from BF file
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		modisSingleData = get_modis_rad(srcFile, (char*)modisResolution.c_str(), singleBandVec, 1, &numCells);
		if (modisSingleData == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get MODIS band.\n";
			return FAILED;
		}
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> numCells: " << numCells << "\n";
		#endif
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG_TIME> Read source Misr single band data	DONE.");
		#endif

		//-------------------------------------------------
		// handle resample method
		srcProcessedData = new double [trgCellNumNoShift];
		int new_src_size = trgCellNumNoShift;
		//Interpolating
		std::string resampleMethod =  inputArgs.GetResampleMethod();
		std::cout << "Interpolating with '" << resampleMethod << "' method on " << inputArgs.GetSourceInstrument() << " by " << bands[i] << ".\n";
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "nnInterpolate")) {
			nnInterpolate(modisSingleData, srcProcessedData, targetNNsrcID, trgCellNumNoShift);
		}
		else if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "summaryInterpolate")) {
			nsrcPixels = new int [trgCellNumNoShift];
			summaryInterpolate(modisSingleData, targetNNsrcID, srcCellNum, srcProcessedData, nsrcPixels, trgCellNumNoShift);
			#if 0 // DEBUG_TOOL
			std::cout << "DBG_TOOL> No nodata values: \n";
			for(int i = 0; i < trgCellNumNoShift; i++) {
				if(nsrcPixels[i] > 0) {
					printf("%d,\t%lf\n", nsrcPixels[i], srcProcessedData[i]);
				}
			}
			#endif
		}
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG> nnInterpolate  DONE.");
		#endif


		//-----------------------------------------------------------------------
		// check if need to shift by MISR (shift==ON & target) case before writing
		double * srcProcessedDataShifted = NULL;
		double * srcProcessedDataPtr = NULL;
		if(inputArgs.GetMISR_Shift() == "ON" && inputArgs.GetTargetInstrument() == "MISR") {
			std::cout << "\nSource MODIS radiance MISR-base shifting...\n";
			#if DEBUG_ELAPSE_TIME
			StartElapseTime();
			#endif
			srcProcessedDataShifted = new double [widthShifted * heightShifted];
			MISRBlockOffset(srcProcessedData, srcProcessedDataShifted, (inputArgs.GetMISR_Resolution() == "L") ? 0 : 1);
			#if DEBUG_ELAPSE_TIME
			StopElapseTimeAndShow("DBG_TIME> source MODIS radiance MISR-base shift DONE.");
			#endif

			srcProcessedDataPtr = srcProcessedDataShifted;
			numCells = widthShifted * heightShifted;
		}
		else { // no misr-trg shift
			srcProcessedDataPtr = srcProcessedData;
			numCells = trgCellNum;
		}

		//---------------------------------
		// write src band to AF file
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		ret = af_WriteSingleRadiance_ModisAsSrc(outputFile, modisDatatype, modisDataspace,  srcProcessedDataPtr, numCells /*processed size*/, srcOutputWidth, i /*bandIdx*/);
		if (ret == FAILED) {
			std::cerr << __FUNCTION__ << "> Error: returned fail.\n";
		}
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG_TIME> Write source Misr single band data  DONE.");
		#endif
		//
		// TODO: buffer is not reused. make memory allocation out of get_modis_rad() to improve performance

		// free memory
		if (modisSingleData)
			free(modisSingleData);
		if(srcProcessedData)
			delete [] srcProcessedData;
		if (nsrcPixels)
			delete [] nsrcPixels;
		if(srcProcessedDataShifted)
			delete [] srcProcessedDataShifted;
	} // i loop

	H5Tclose(modisDatatype);
	H5Sclose(modisDataspace);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif

	return ret;
}


/*=================================================
 * Source as MISR functions for the output
 */
static int af_WriteSingleRadiance_MisrAsSrc(hid_t outputFile, hid_t misrDatatype, hid_t misrFilespace, double* processedData, int trgCellNum, int outputWidth, int cameraIdx, int radIdx)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;
	herr_t status;
	hid_t misr_dataset;
	std::string dsetPath = SRC_DATA_GROUP + "/" + MISR_RADIANCE_DSET;

	if(cameraIdx==0 && radIdx==0) { // means new
		misr_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), misrDatatype, misrFilespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(misr_dataset < 0) {
			std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dcreate2 target data in output file.\n";
			return FAILED;
		}
	}
	else {
		misr_dataset = H5Dopen2(outputFile, dsetPath.c_str(), H5P_DEFAULT);
		if(misr_dataset < 0) {
			std::cerr <<  __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dopen2 target data in output file.\n";
			return FAILED;
		}
	}


	//------------------------------
	// select memspace
	int ranksMem=2;
	hsize_t dim2dMem[ranksMem];
	hsize_t start2dMem[ranksMem];
	hsize_t count2dMem[ranksMem];
	dim2dMem[0] = trgCellNum/outputWidth; // y
	dim2dMem[1] = outputWidth; // x
	hid_t memspace = H5Screate_simple(ranksMem, dim2dMem, NULL);

	start2dMem[0] = 0; // y
	start2dMem[1] = 0; // x
	count2dMem[0] = trgCellNum/outputWidth; // y
	count2dMem[1] = outputWidth; // x

	status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start2dMem, NULL, count2dMem, NULL);

	//------------------------------
	// select filespace for misr
	const int ranksFile=4; // [cameras][radiances][y][x]
	hsize_t startFile[ranksFile];
	hsize_t countFile[ranksFile];
	startFile[0] = cameraIdx;
	startFile[1] = radIdx;
	startFile[2] = 0; // y
	startFile[3] = 0; // x
	countFile[0] = 1;
	countFile[1] = 1;
	countFile[2] = trgCellNum/outputWidth; // y
	countFile[3] = outputWidth;  // x

	status = H5Sselect_hyperslab(misrFilespace, H5S_SELECT_SET, startFile, NULL, countFile, NULL);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Sselect_hyperslab for Misr target .\n";
		ret = -1;
		goto done;
	}

	status = H5Dwrite(misr_dataset, H5T_NATIVE_DOUBLE, memspace, misrFilespace, H5P_DEFAULT, processedData);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dwrite for Misr target .\n";
		ret = -1;
		goto done;
	}

done:
	H5Dclose(misr_dataset);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
	return ret;
}


int af_GenerateOutputCumulative_MisrAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNum, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;

	// strVec_t multiVarNames = inputArgs.GetMultiVariableNames("MISR"); // misr_MultiVars;
	std::string misrResolution = inputArgs.GetMISR_Resolution();

	// two multi-value variables are expected as this point
	strVec_t cameras = inputMultiVarsMap[MISR_CAMERA_ANGLE];
	strVec_t radiances = inputMultiVarsMap[MISR_RADIANCE];

	// data type
	hid_t misrDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t	status = H5Tset_order(misrDatatype, H5T_ORDER_LE);
	if(status < 0) {
		printf("Error: MISR write error in H5Tset_order\n");
		return FAILED;
	}


	/*------------------------
	 * create space for misr
	 */
	const int rankSpace=4;  // [cameras][radiances][y][x]
	hsize_t misrDims[rankSpace];
	/* handle different data width and height for output. no need to care for misr-trg shift in this
	 */
	int widthShifted = 0;
	int heightShifted = 0;
	ret = af_GetWidthAndHeightForOutputDataSize(inputArgs.GetTargetInstrument() /*target base output*/, inputArgs, widthShifted, heightShifted);
	if( (ret < 0) || (heightShifted > 0) ) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error in af_GetWidthAndHeightForOutputDataSize() \n";
		return FAILED;
	}
    int srcOutputWidth = widthShifted;
	misrDims[0] = cameras.size();
	misrDims[1] = radiances.size();
	misrDims[2] = trgCellNum/srcOutputWidth; // NY;
	misrDims[3] = srcOutputWidth; // NX;
	hid_t misrDataspace = H5Screate_simple(rankSpace, misrDims, NULL);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> trgCellNum: " << trgCellNum << ", srcCellNum: " << srcCellNum << "\n";
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> srcOutputWidth: " << srcOutputWidth <<  "\n";
	#endif
	int numCells;
	double *misrSingleData=NULL;

	std::string singleRad;
	std::string singleCamera;
	//-----------------------------------------------------------------
	// TODO: improve by preparing these memory allocation out of loop
	// srcProcessedData , nsrcPixels
	// misrSingleData 
	double * srcProcessedData = NULL;
	int * nsrcPixels = NULL;
	// Note: This is Combination case only
	for(int j=0; j < cameras.size(); j++) {
		singleCamera = cameras[j];
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> cameras[" << j << "]" << cameras[j] << "\n";
		#endif
		for (int i=0; i< radiances.size(); i++) {
			#if DEBUG_TOOL
			std::cout << "DBG_TOOL " << __FUNCTION__ << "> radiances[" << i << "]" << radiances[i] << "\n";
			#endif
			// insert to only begin
			singleRad =  radiances[i];
	
			//---------------------------------
			// read src radiance from BF file
			#if DEBUG_ELAPSE_TIME
			StartElapseTime();
			#endif
			misrSingleData = get_misr_rad(srcFile, (char*) singleCamera.c_str(), (char*)misrResolution.c_str(), (char*)singleRad.c_str(), &numCells);
			if (misrSingleData == NULL) {
				std::cerr << __FUNCTION__ <<  "> Error: failed to get MISR radiance.\n";
				return FAILED;
			}
			#if DEBUG_TOOL
			std::cout << "DBG_TOOL " << __FUNCTION__ << "> numCells: " << numCells << "\n";
			#endif
			#if DEBUG_ELAPSE_TIME
			StopElapseTimeAndShow("DBG_TIME> Read source Misr single band data	DONE.");
			#endif
	
			//-------------------------------------------------
			// handle resample method
			srcProcessedData = new double [trgCellNum];
			int new_src_size = trgCellNum;
			//Interpolating
			std::string resampleMethod =  inputArgs.GetResampleMethod();
			std::cout << "Interpolating with '" << resampleMethod << "' method on " << inputArgs.GetSourceInstrument() << " by " << cameras[j] << " : " << radiances[i] << ".\n";
			#if DEBUG_ELAPSE_TIME
			StartElapseTime();
			#endif
			if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "nnInterpolate")) {
				nnInterpolate(misrSingleData, srcProcessedData, targetNNsrcID, trgCellNum);
			}
			else if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "summaryInterpolate")) {
				nsrcPixels = new int [trgCellNum];
				summaryInterpolate(misrSingleData, targetNNsrcID, srcCellNum, srcProcessedData, nsrcPixels, trgCellNum);
				#if 0 // DEBUG_TOOL
				std::cout << "DBG_TOOL> No nodata values: \n";
				for(int i = 0; i < trgCellNum; i++) {
					if(nsrcPixels[i] > 0) {
						printf("%d,\t%lf\n", nsrcPixels[i], srcProcessedData[i]);
					}
				}
				#endif
			}
			#if DEBUG_ELAPSE_TIME
			StopElapseTimeAndShow("DBG> nnInterpolate  DONE.");
			#endif
	
			//---------------------------------
			// write src radiance to AF file
			#if DEBUG_ELAPSE_TIME
			StartElapseTime();
			#endif
			ret = af_WriteSingleRadiance_MisrAsSrc(outputFile, misrDatatype, misrDataspace,  srcProcessedData, trgCellNum /*processed size*/, srcOutputWidth, j /*cameraIdx*/, i /*radIdx*/);
			if (ret == FAILED) {
				std::cerr << __FUNCTION__ << "> Error: returned fail.\n";
			}
			#if DEBUG_ELAPSE_TIME
			StopElapseTimeAndShow("DBG_TIME> Write source Misr single band data  DONE.");
			#endif
			//
			// TODO: buffer is not reused. make memory allocation out of get_misr_rad() to improve performance
	
			// free memory
			if (misrSingleData)
				free(misrSingleData);
			if(srcProcessedData)
				delete [] srcProcessedData;
			if (nsrcPixels)
				delete [] nsrcPixels;
		} // i loop
	} // j loop

	H5Tclose(misrDatatype);
	H5Sclose(misrDataspace);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif

	return ret;
}





/*===========================================
 * Main 
 */

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
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL main> src instrument: " << srcInstrument << std::endl;
	std::cout << "DBG_TOOL main> target instrument: " << trgInstrument << std::endl;
	#endif
	
#if 0 // TEST : multi-value variable map , remove later
	//---------------------------------------------------
	// build target instrument multi-value variable map
	std::map<std::string, strVec_t> trgInputMultiVarsMapTEST;
	inputArgs.BuildMultiValueVariableMap(trgInstrument, trgInputMultiVarsMapTEST);
	inputArgs.DBG_displayinputListMap(trgInstrument, trgInputMultiVarsMapTEST, "COMBINATION");

	//---------------------------------------------------
	// build source instrument multi-value variable map 
	std::map<std::string, strVec_t> srcInputMultiVarsMapTEST;
	inputArgs.BuildMultiValueVariableMap(srcInstrument, srcInputMultiVarsMapTEST);
	inputArgs.DBG_displayinputListMap(srcInstrument, srcInputMultiVarsMapTEST, "COMBINATION");
	exit(1);
#endif


	/* ===================================================
	 * Create output file
	 */
	std::string outputFile = inputArgs.GetOuputFilePath();
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL main> outputFile: " << outputFile << std::endl;
	#endif

	hid_t output_file = H5Fcreate(outputFile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


	/* ===================================================
	 * Handle input BF data file
	 */
	std::string inputDataPath = inputArgs.GetInputBFdataPath();
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL main> inputDataPath: " << inputDataPath << std::endl;
	#endif
	hid_t inputFile;
	if(0 > (inputFile = af_open((char*)inputDataPath.c_str()))) {
		std::cerr << "Error: File not found - " << inputDataPath << std::endl;
		exit(1);
	}


	
	/* ===================================================
	 * Get Source instrument latitude and longitude
	 */
	std::cout << "\nGetting source instrument latitude & longitude data...\n";
	int srcCellNum;
	double* srcLatitude = NULL;
	double* srcLongitude = NULL;
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	ret = AF_GetGeolocationDataFromInstrument(srcInstrument, inputArgs, inputFile, &srcLatitude /*OUT*/, &srcLongitude /*OUT*/, srcCellNum /*OUT*/);
	if (ret == FAILED) {
		std::cerr << __FUNCTION__ << "> Error getting geolocation data from source instrument - " << srcInstrument << ".\n";
		return FAILED;
	}
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> get source lat/long DONE.");
	#endif
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL main> srcCellNum: " <<  srcCellNum << "\n";
	#endif



	/* ===================================================
	 * Get Target instrument latitude and longitude
	 */
	std::cout << "\nGetting target instrument latitude & longitude data...\n";
	int trgCellNum;
	double* targetLatitude = NULL;
	double* targetLongitude = NULL;
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	ret = AF_GetGeolocationDataFromInstrument(trgInstrument, inputArgs, inputFile, &targetLatitude /*OUT*/, &targetLongitude /*OUT*/, trgCellNum /*OUT*/);
	if (ret == FAILED) {
		std::cerr << __FUNCTION__ << "> Error: getting geolocation data from target instrument - " << trgInstrument << ".\n";
		return FAILED;
	}
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> get target lat/long DONE.");
	#endif
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL main> trgCellNum: " <<  trgCellNum << "\n";
	#endif
	


	/* ===================================================
	 * Output Target instrument latitude and longitude
	 */
	int trgCellNumNew;
	int trgOutputWidth;
	double * targetLatitudePtr = NULL;
	int widthShifted;
	int heightShifted;
	double * targetLatitudeShifted = NULL;

	ret = af_GetWidthAndHeightForOutputDataSize(inputArgs.GetTargetInstrument() /* target base output */, inputArgs, widthShifted, heightShifted);
	trgOutputWidth = widthShifted;
	// misr-target base shift
	if(inputArgs.GetMISR_Shift() == "ON" && inputArgs.GetTargetInstrument() == "MISR") {
		std::cout << "Target latitude MISR-base block unstacking...\n";
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		targetLatitudeShifted = (double *) malloc(sizeof(double) * widthShifted * heightShifted);
		std::string misrResolution = inputArgs.GetMISR_Resolution();
		MISRBlockOffset(targetLatitude, targetLatitudeShifted, (misrResolution == "L") ? 0 : 1);
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG_TIME> Target latitude MISR-base block unstacking DONE.");
		#endif

		targetLatitudePtr = targetLatitudeShifted;
		trgCellNumNew = widthShifted * heightShifted;
	}
	else { // no misr shift
		targetLatitudePtr = targetLatitude;
		trgCellNumNew = trgCellNum;
	}
	

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL main> trgOutputWidth: " <<  trgOutputWidth << "\n";
	#endif

	std::cout << "\nWriting target geolocation data...\n";
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	int lat_status =  af_write_mm_geo(output_file, 0, targetLatitudePtr, trgCellNumNew, trgOutputWidth);
	if(lat_status < 0) {
		std::cerr << "Error: writing latitude geolocation.\n";
		return FAILED;
	}
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> write geo lattitude data DONE.");
	#endif

	if (targetLatitudeShifted)
		free(targetLatitudeShifted);


	double * targetLongitudePtr = NULL;
	double * targetLongitudeShifted = NULL;

	if(inputArgs.GetMISR_Shift() == "ON" && inputArgs.GetTargetInstrument() == "MISR") {
		std::cout << "Target longitude MISR-base block unstacking...\n";
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		targetLongitudeShifted = (double *) malloc(sizeof(double) * widthShifted * heightShifted);
		std::string misrResolution = inputArgs.GetMISR_Resolution();
		MISRBlockOffset(targetLongitude, targetLongitudeShifted, (misrResolution == "L") ? 0 : 1);
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG_TIME> Target longitude MISR-base block unstacking DONE.");
		#endif

		targetLongitudePtr = targetLongitudeShifted;
	}
	else { // no misr shift
		targetLongitudePtr = targetLongitude;
	}

	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	int long_status = af_write_mm_geo(output_file, 1, targetLongitudePtr, trgCellNumNew, trgOutputWidth);
	if(long_status < 0) {
		std::cerr << "Error: writing longitude geolocation.\n";
		return FAILED;
	}
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> write geo longitude data DONE.");
	#endif

	if (targetLongitudeShifted)
		free(targetLongitudeShifted);



	/* ===========================================================
	 * Calculate nearest neighbot source over target geolocation
	 */
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
	


	/*======================================================
	 * Target instrument: Generate radiance to output file
	 */
	std::cout  <<  "\nGenerating target instrument " << trgInstrument << " radiance output...\n";
	// build target instrument multi-value variable map
	std::map<std::string, strVec_t> trgInputMultiVarsMap;
	ret = inputArgs.BuildMultiValueVariableMap(trgInstrument, trgInputMultiVarsMap);
	if (ret < 0) {
		std::cerr << __FUNCTION__ << "> Error: build multi-value variable map for " << trgInstrument << ".\n";
		return FAILED;
	}
	// write target instrument radiances to output file
	ret = AF_GenerateTargetRadiancesOutput(inputArgs, output_file, trgCellNum, inputFile, trgInputMultiVarsMap);
	if (ret < 0) {
		std::cerr << "Error: generate target radiance output.\n";
		return FAILED;
	}
	std::cout << "Writing target radiance output done.\n";



	/*======================================================
	 * Source instrument: Generate radiance to output file
	 */
	std::cout  <<  "\nGenerating source instrument " << srcInstrument << " radiance output...\n";
	// build source instrument multi-value variable map
	std::map<std::string, strVec_t> srcInputMultiVarsMap;
	ret = inputArgs.BuildMultiValueVariableMap(srcInstrument, srcInputMultiVarsMap);
	if (ret < 0) {
		std::cerr << __FUNCTION__ << "> Error: build multi-value variable map for " << srcInstrument << ".\n";
		return FAILED;
	}
	// write source instrument radiances to output file
	ret = AF_GenerateSourceRadiancesOutput(inputArgs, output_file, targetNNsrcID, trgCellNum, inputFile, srcCellNum, srcInputMultiVarsMap);
	if (ret < 0) {
		std::cerr << "Error: generate source radiance output.\n";
		return FAILED;
	}
	std::cout << "Writing source radiance output done.\n";

	if (targetNNsrcID)
		delete [] targetNNsrcID;



	//==========================================
	// Closing data files
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	std::cout  <<  "\nClosing file...\n";
	herr_t close_status = af_close(inputFile);
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
	StopElapseTimeAndShow("DBG_TIME> Closing file DONE.");
	#endif

	return 0;
}
