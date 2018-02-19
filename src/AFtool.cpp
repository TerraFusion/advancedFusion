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

#define SUCCEED 0
#define FAILED -1

const hsize_t RESAMPLE_MODIS_DATA_WIDTH = 1354; // 1KM
//const hsize_t RESAMPLE_MODIS_DATA_WIDTH = 2708; // 500m
//const hsize_t RESAMPLE_MODIS_DATA_WIDTH = 5416; // 250m
const hsize_t RESAMPLE_MISR_DATA_WIDTH = 1354;

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
		*longitude = get_modis_long(inputFile, (char*) resolution.c_str(), &cellNum);
	}
	else if (instrument == "MISR") {
		std::string resolution = inputArgs.GetMISR_Resolution();
		*latitude = get_misr_lat(inputFile, (char*) resolution.c_str(), &cellNum);
		*longitude = get_misr_long(inputFile, (char*) resolution.c_str(), &cellNum);
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

/*##############################################################
 *
 * Generate Target instrument radiance to output file
 *
 *#############################################################*/

/*==================================
 * Target output as MODIS 
 */
static int af_WriteSingleRadiance_ModisAsTrg(hid_t outputFile, hid_t modisDatatype, hid_t modisFilespace, double* modisData, int modisDataSize, int bandIdx)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

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
	const int ranksMem=2;
	hsize_t dim2dMem[2];
	hsize_t start2dMem[2];
	hsize_t count2dMem[2];
	dim2dMem[0] = modisDataSize/RESAMPLE_MODIS_DATA_WIDTH; // y
	dim2dMem[1] = RESAMPLE_MODIS_DATA_WIDTH; // x
	hid_t memspace = H5Screate_simple(ranksMem, dim2dMem, NULL);

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

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
	return ret;
}


int af_GenerateOutputCumulative_ModisAsTrg(hid_t outputFile,hid_t srcFile, int trgCellNum, std::string modisResolution,  std::vector<std::string> &bands /* multi */)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	// data type
	hid_t modisDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t	status = H5Tset_order(modisDatatype, H5T_ORDER_LE);
	if(status < 0) {
		printf("Error: MODIS write error in H5Tset_order\n");
		return FAILED;
	}


	hsize_t modis_dim[3];
	modis_dim[0] = bands.size();
	modis_dim[1] = trgCellNum/RESAMPLE_MODIS_DATA_WIDTH; // NY;
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
		af_WriteSingleRadiance_ModisAsTrg(outputFile, modisDatatype, modisDataspace,  modisSingleData, numCells, i);
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
}

int   AF_GenerateTargetRadiancesOutput(AF_InputParmeterFile &inputArgs, hid_t outputFile, int trgCellNum, hid_t srcFile, std::map<std::string, strVec_t> & trgInputMultiVarsMap)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	// JK_TODO - change to "/Target/Data_Fields"
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

		// TODO_LATER - decide for loop inside or loop outside
		#if 1 // case for looping inside
		// [0] is cause we know modis has only one multi-value valriable at this point
		af_GenerateOutputCumulative_ModisAsTrg(outputFile, srcFile, trgCellNum, inputArgs.GetMODIS_Resolution(), trgInputMultiVarsMap[multiVarNames[0]]);
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
		multiVarNames = inputArgs.GetMultiVariableNames("MISR"); // misr_MultiVars;

		if (trgInputMultiVarsMap.size() != 2) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building target input list with MISR. There must be only two multi-value variables.\n";
			return FAILED;
		}

		std::string misrCamera;
		std::string misrRadiance;
		if(mixType == "COMBINATION") {
			// JK_TODO_LATER - Try inside loop first like MODIS
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

/*==================================
 * Source output as MODIS 
 */
#if 0 // JK_TODO_LATER // place holder for Modis as source output
static int af_WriteSingleRadiance_ModisAsSrc(hid_t outputFile, hid_t modisDatatype, hid_t modisFilespace, double* modisData, int modisDataSize, int bandIdx)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

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
	int ranksMem=2;
	hsize_t dim2dMem[2];
	hsize_t start2dMem[2];
	hsize_t count2dMem[2];
	dim2dMem[0] = modisDataSize/RESAMPLE_MODIS_DATA_WIDTH; // y
	dim2dMem[1] = RESAMPLE_MODIS_DATA_WIDTH; // x
	hid_t memspace = H5Screate_simple(ranksMem, dim2dMem, NULL);

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

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
	return ret;
}


int af_GenerateOutputCumulative_ModisAsSrc(hid_t outputFile,hid_t srcFile, int singleDataSize, std::string modisResolution,  std::vector<std::string> &bands)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

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
		#if 0 // TOOL_TEST remove later
		//Read multi bands (instead Generate single band data for test)
		//GenerateDataDouble2D(singleData, singleDataDims2D, atoi(bands[i].c_str()));
		//DisplayDataDouble2D(singleData, singleDataDims2D);
		#endif
		//---------------------------------
		// read src radiance from BF file
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
		// write src radiance to AF file
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		af_WriteSingleRadiance_ModisAsSrc(outputFile, modisDatatype, modisDataspace,  modisSingleData, numCells, i);
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
}
#endif // JK_TODO_LATER


/*==================================
 * Source output as MISR 
 */
static int af_WriteSingleRadiance_MisrAsSrc(hid_t outputFile, hid_t misrDatatype, hid_t misrFilespace, double* misrData, int misrDataSize, int cameraIdx, int radIdx)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;
	herr_t status;
	hid_t misr_dataset;

	if(cameraIdx==0 && radIdx==0) { // means new
		misr_dataset = H5Dcreate2(outputFile, "/Data_Fields/misr_out", misrDatatype, misrFilespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(misr_dataset < 0) {
			std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dcreate2 target data in output file.\n";
			return FAILED;
		}
	}
	else {
		misr_dataset = H5Dopen2(outputFile, "/Data_Fields/misr_out", H5P_DEFAULT);
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
	dim2dMem[0] = misrDataSize/RESAMPLE_MISR_DATA_WIDTH; // y
	dim2dMem[1] = RESAMPLE_MISR_DATA_WIDTH; // x
	hid_t memspace = H5Screate_simple(ranksMem, dim2dMem, NULL);

	start2dMem[0] = 0; // y
	start2dMem[1] = 0; // x
	count2dMem[0] = misrDataSize/RESAMPLE_MISR_DATA_WIDTH; // y
	count2dMem[1] = RESAMPLE_MISR_DATA_WIDTH; // x

	status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start2dMem, NULL, count2dMem, NULL);

	//------------------------------
	// select filespace
	const int ranksFile=4;
	hsize_t startFile[ranksFile];
	hsize_t countFile[ranksFile];
	startFile[0] = cameraIdx;
	startFile[1] = radIdx;
	startFile[2] = 0; // y
	startFile[3] = 0; // x
	countFile[0] = 1;
	countFile[1] = 1;
	countFile[2] = misrDataSize/RESAMPLE_MISR_DATA_WIDTH; // y
	countFile[3] = RESAMPLE_MISR_DATA_WIDTH;  // x

	status = H5Sselect_hyperslab(misrFilespace, H5S_SELECT_SET, startFile, NULL, countFile, NULL);
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


int af_GenerateOutputCumulative_MisrAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNum, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &srcInputMultiVarsMap)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;

	strVec_t multiVarNames = inputArgs.GetMultiVariableNames("MISR"); // misr_MultiVars;
	std::string misrResolution = inputArgs.GetMISR_Resolution();

	// two multi-value variables are expected as this point
	strVec_t cameras = srcInputMultiVarsMap[MISR_CAMERA_ANGLE];
	strVec_t radiances = srcInputMultiVarsMap[MISR_RADIANCE];

	// data type
	hid_t misrDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t	status = H5Tset_order(misrDatatype, H5T_ORDER_LE);
	if(status < 0) {
		printf("Error: MISR write error in H5Tset_order\n");
		return FAILED;
	}


	// create space
	const int rankSpace=4;
	hsize_t misrDims[rankSpace];
	misrDims[0] = cameras.size();
	misrDims[1] = radiances.size();
	misrDims[2] = trgCellNum/RESAMPLE_MISR_DATA_WIDTH; // NY;
	misrDims[3] = RESAMPLE_MISR_DATA_WIDTH; // NX;
	hid_t misrDataspace = H5Screate_simple(rankSpace, misrDims, NULL);

	int numCells;
	double *misrSingleData=NULL;

	std::string singleRad;
	std::string singleCamera;
	//-----------------------------------------------------------------
	// TODO: improve by preparing these memory allocation out of loop
	// src_rad_out , nsrcPixels
	// misrSingleData 
	double * src_rad_out = NULL;
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
			#if DEBUG_TOOL
			std::cout << "DBG_TOOL " << __FUNCTION__ << "> numCells: " << numCells << "\n";
			#endif
			#if DEBUG_ELAPSE_TIME
			StopElapseTimeAndShow("DBG_TIME> Read target Misr single band data  DONE.");
			#endif
	
			//-------------------------------------------------
			// handle resample method
			src_rad_out = new double [trgCellNum];
			int new_src_size = trgCellNum;
			//Interpolating
			std::string resampleMethod =  inputArgs.GetResampleMethod();
			std::cout << "\nInterpolating using '" << resampleMethod << "' method...\n";
			#if DEBUG_ELAPSE_TIME
			StartElapseTime();
			#endif
			if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "nnInterpolate")) {
				nnInterpolate(misrSingleData, src_rad_out, targetNNsrcID, trgCellNum);
			}
			else if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "summaryInterpolate")) {
				nsrcPixels = new int [trgCellNum];
				summaryInterpolate(misrSingleData, targetNNsrcID, srcCellNum, src_rad_out, nsrcPixels, trgCellNum);
				#if 0 // DEBUG_TOOL
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
	
			//---------------------------------
			// write src radiance to AF file
			#if DEBUG_ELAPSE_TIME
			StartElapseTime();
			#endif
			ret = af_WriteSingleRadiance_MisrAsSrc(outputFile, misrDatatype, misrDataspace,  src_rad_out, numCells, j /*cameraIdx*/, i /*radIdx*/);
			if (ret == FAILED) {
				std::cerr << __FUNCTION__ << "> Error: returned fail.\n";
			}
			#if DEBUG_ELAPSE_TIME
			StopElapseTimeAndShow("DBG_TIME> Write target Misr single band data  DONE.");
			#endif
			//
			// TODO: buffer is not reused. make memory allocation out of get_misr_rad() to improve performance
	
			// free memory
			if (misrSingleData)
				free(misrSingleData);
	    	if(src_rad_out)
	        	delete [] src_rad_out;
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

/*===============================================
 * Generate source instrument randiance output.
 * Initial function.
 */
int   AF_GenerateSourceRadiancesOutput(AF_InputParmeterFile &inputArgs, hid_t outputFile, int * targetNNsrcID, int trgCellNum, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> & srcInputMultiVarsMap)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif
	int ret = SUCCEED;

	#if 0 // JK_TODO_LATER - change to  "/Source/Data_Fields"
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
	#endif

	strVec_t multiVarNames;
	std::string instrument = inputArgs.GetSourceInstrument();
	std::string mixType = "COMBINATION"; // Default
	//--------------------------------------------------------
	// Use this for Single multi-value Variable case. MODIS

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> srcInputMultiVarsMap.size(): " << srcInputMultiVarsMap.size() << "\n";
	#endif

	if (instrument == "MODIS") {
		multiVarNames = inputArgs.GetMultiVariableNames("MODIS"); // modis_MultiVars;

		#if 0 // TODO LATER
		if (srcInputMultiVarsMap.size() != 1) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building target input list with MODIS. There must be only one multi-value variable.\n";
			return FAILED;
		}

		// TODO - decide for loop inside or loop outside
		#if 1 // case for looping inside
		// [0] is cause we know modis has only one multi-value valriable at this point
		af_GenerateOutputCumulative_ModisAsSrc(outputFile, srcFile, trgCellNum, inputArgs.GetMODIS_Resolution(), srcInputMultiVarsMap[multiVarNames[0]]);
		#else // case for looping from outside
		for(int j = 0; j < srcInputMultiVarsMap[multiVarNames[0]].size(); ++j) {
			std::cout << srcInputMultiVarsMap[multiVarNames[0]][j]	<< ", ";
			af_GenerateOutputCumulative_ModisAsSrc(outputFile, srcFile, trgCellNum, srcInputMultiVarsMap[multiVarNames[0]][j] /* single Band*/);
		}
		std::cout << std::endl;
		#endif

		#endif // TODO LATER
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
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> returned failure.\n";
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
	
#if 0 // TEST : multi-value variable map , remove later
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
	 * Get Source instrument latitude and longitude
	 */
	std::cout << "\nGetting source instrument latitude & longitude data...\n";
	int srcCellNum;
	double* srcLatitude = NULL;
	double* srcLongitude = NULL;
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	ret = AF_GetGeolocationDataFromInstrument(srcInstrument, inputArgs, src_file, &srcLatitude /*OUT*/, &srcLongitude /*OUT*/, srcCellNum /*OUT*/);
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

	/* ---------------------------------------------------
	 * Get Target instrument latitude and longitude
	 */
	int trgCellNum;
	double* targetLatitude = NULL;
	double* targetLongitude = NULL;
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	ret = AF_GetGeolocationDataFromInstrument(trgInstrument, inputArgs, src_file, &targetLatitude /*OUT*/, &targetLongitude /*OUT*/, trgCellNum /*OUT*/);
	if (ret == FAILED) {
		std::cerr << __FUNCTION__ << "> Error getting geolocation data from target instrument - " << trgInstrument << ".\n";
		return FAILED;
	}
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> get source lat/long DONE.");
	#endif
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL main> trgCellNum: " <<  trgCellNum << "\n";
	#endif
	
	/* ---------------------------------------------------
	 * Output Target instrument latitude and longitude
	 */
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


	/* -----------------------------------------------------------
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
	

	/*------------------------------------------------------
	 * Target instrument: Generate radiance to output file
	 */
	// build target instrument multi-value variable map
	std::map<std::string, strVec_t> trgInputMultiVarsMap;
	inputArgs.BuildMultiValueVariableMap(trgInstrument, trgInputMultiVarsMap);
	// write target instrument radiances to output file
	ret = AF_GenerateTargetRadiancesOutput(inputArgs, output_file, trgCellNum, src_file, trgInputMultiVarsMap);
	if (ret < 0) {
		return FAILED;
	}


	/*------------------------------------------------------
	 * Source instrument: Generate radiance to output file
	 */
	// build source instrument multi-value variable map
	std::map<std::string, strVec_t> srcInputMultiVarsMap;
	inputArgs.BuildMultiValueVariableMap(srcInstrument, srcInputMultiVarsMap);
	// write source instrument radiances to output file
	ret = AF_GenerateSourceRadiancesOutput(inputArgs, output_file, targetNNsrcID, trgCellNum, src_file, srcCellNum, srcInputMultiVarsMap);
	if (ret < 0) {
		return FAILED;
	}

	if (targetNNsrcID)
		delete [] targetNNsrcID;


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
