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
#include <sys/time.h>
#include <vector>
#include <sstream>

#include <hdf5.h> // JK_RE may rm
#include "io.h"
#include "reproject.h" // JK_RE may rm
#include "AF_InputParmeterFile.h"
#include "AF_debug.h"
#include "misrutil.h"
#include "gdalio.h"
#include "AF_common.h" // JK_RE may rm
#include "AF_output_util.h" // JK_RE may rm
#include "AF_output_MODIS.h"
#include "AF_output_MISR.h"


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

	/*======================================================
 	 * MODIS section
 	 */
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
	/*======================================================
 	 * MISR section
 	 */
	else if (instrument == "MISR") {
		std::string resolution = inputArgs.GetMISR_Resolution();
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> Misr resolution: " << resolution << "\n";
		#endif
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

	}
	/*======================================================
 	 * ASTER section
 	 */
	else if (instrument == "ASTER" ) {
		//Get ASTER input parameters EX: "TIR", "ImageData10"
		std::string resolution = inputArgs.GetASTER_Resolution();
		strVec_t bands = inputArgs.GetASTER_Bands();
		// pass the first band (bands[0]) as it always exists and lat&lon is the same for a resolution.
		*latitude = get_ast_lat(inputFile, (char*) resolution.c_str(), (char*)bands[0].c_str(), &cellNum);
		if (latitude == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get ASTER latitude.\n";
			return FAILED;
		}
		*longitude = get_ast_long(inputFile, (char*) resolution.c_str(), (char*)bands[0].c_str(), &cellNum);
		if (longitude == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get ASTER longitude.\n";
			return FAILED;
		}
	}
	/*======================================================
	 * USER_DEFINE section
	 */
	else if (instrument == "USER_DEFINE") {
		// targetX == lon , targetY == lat, nPoints == cellNum according to test_userdefinedgrids.cpp
	    int userOuputEPSG = inputArgs.GetUSER_EPSG();
	    double userXmin = inputArgs.GetUSER_xMin();
	    double userXmax = inputArgs.GetUSER_xMax();
	    double userYmin = inputArgs.GetUSER_yMin();
	    double userYmax = inputArgs.GetUSER_yMax();
	    double userRsolution = inputArgs.GetUSER_Resolution();
		cellNum = getCellCenterLatLon(userOuputEPSG, userXmin, userYmin, userXmax, userYmax, userRsolution, longitude, latitude);

		#if DEBUG_TOOL
		double *lonPtr = (double*)*longitude;
		double *latPtr = (double*)*latitude;
		for(int i = 0; i < 10; i++) {
			printf("JKDBG> i:%d, X:%lf,\t Y:%lf\n", i, lonPtr[i], latPtr[i]);
		}
		#endif
	}
	else {
		std::cerr << __FUNCTION__ << "> Error: invalid instrument - " << instrument << "\n";
		return FAILED;
	}


	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> Instrument: " << instrument << " cellNum: " << cellNum << "\n";
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
	
	return SUCCEED;
}



/*##############################################################
 *
 * Generate Target instrument radiance to output file
 *
 *#############################################################*/

/* 
 * function prototype for AF_GenerateTargetRadiancesOutput()
 */

/*=================================================
 * Generate Target instrument randiance outputs.
 * Initial function.
 */
int   AF_GenerateTargetRadiancesOutput(AF_InputParmeterFile &inputArgs, hid_t outputFile, int trgCellNum, hid_t srcFile, std::map<std::string, strVec_t> & trgInputMultiVarsMap)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret;

	std::string instrument = inputArgs.GetTargetInstrument();

	// there is no radiance for this case. just skip.
	if (instrument == "USER_DEFINE") {
		return SUCCEED;
	}

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
	// std::string mixType = "COMBINATION"; // Default
	//--------------------------------------------------------
	// Use this for Single multi-value Variable case. MODIS

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> trgInputMultiVarsMap.size(): " << trgInputMultiVarsMap.size() << "\n";
	#endif
	/*======================================================
 	 * MODIS section
 	 */
	if (instrument == "MODIS") {
		multiVarNames = inputArgs.GetMultiVariableNames("MODIS"); // modis_MultiVars;

		// for one multi-value Variable case.
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
	/*======================================================
 	 * MISR section
 	 */
	else if (instrument == "MISR") {

		// for two multi-value Variable case.
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




/*##############################################################
 *
 * Generate Source instrument radiance to output file
 *
 *#############################################################*/

/* 
 * function prototype for AF_GenerateTargetRadiancesOutput()
 */ 
int af_GenerateOutputCumulative_AsterAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNumNoShift, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap);


/*===================================================================
 * Generate Source instrument randiance outputs.
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

	/*======================================================
 	 * MODIS section
 	 */
	if (instrument == "MODIS") {
		// for one multi-value Variable case.
		if (srcInputMultiVarsMap.size() != 1) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building source input list with MODIS. There must be only one multi-value variable.\n";
			return FAILED;
		}

		ret = af_GenerateOutputCumulative_ModisAsSrc(inputArgs, outputFile, targetNNsrcID,  trgCellNum, srcFile, srcCellNum, srcInputMultiVarsMap);
		if (ret == FAILED) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> failed generating output for MODIS.\n";
			ret = FAILED;
			goto done;
		}
	}
	/*======================================================
 	 * MISR section
 	 */
	else if (instrument == "MISR") {
		// for two multi-value Variable case.
		if (srcInputMultiVarsMap.size() != 2) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building source input list with MISR. There must be only two multi-value variables.\n";
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
	/*======================================================
 	 * ASTER section
 	 */
	else if (instrument == "ASTER") {
		// for one multi-value Variable case.
		if (srcInputMultiVarsMap.size() != 1) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building source input list with ASTER. There must be only one multi-value variable.\n";
			return FAILED;
		}

		ret = af_GenerateOutputCumulative_AsterAsSrc(inputArgs, outputFile, targetNNsrcID,  trgCellNum, srcFile, srcCellNum, srcInputMultiVarsMap);
		if (ret == FAILED) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> failed generating output for ASTER.\n";
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
 * ASTER section as Source. functions for the output.
 */
static int af_WriteSingleRadiance_AsterAsSrc(hid_t outputFile, hid_t asterDatatype, hid_t asterFilespace, double* processedData, int trgCellNum, int outputWidth, int bandIdx)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;
	herr_t status;
	hid_t aster_dataset;
	std::string dsetPath = SRC_DATA_GROUP + "/" + ASTER_RADIANCE_DSET;

	if(bandIdx==0) { // means new
		aster_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), asterDatatype, asterFilespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(aster_dataset < 0) {
			std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dcreate2 target data in output file.\n";
			return FAILED;
		}
	}
	else {
		aster_dataset = H5Dopen2(outputFile, dsetPath.c_str(), H5P_DEFAULT);
		if(aster_dataset < 0) {
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
	// select filespace for aster
	const int ranksFile=3; // [bands][y][x]
	hsize_t startFile[ranksFile];
	hsize_t countFile[ranksFile];
	startFile[0] = bandIdx;
	startFile[1] = 0; // y
	startFile[2] = 0; // x
	countFile[0] = 1;
	countFile[1] = trgCellNum/outputWidth; // y
	countFile[2] = outputWidth;  // x

	status = H5Sselect_hyperslab(asterFilespace, H5S_SELECT_SET, startFile, NULL, countFile, NULL);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Sselect_hyperslab for Aster target .\n";
		ret = -1;
		goto done;
	}

	status = H5Dwrite(aster_dataset, H5T_NATIVE_DOUBLE, memspace, asterFilespace, H5P_DEFAULT, processedData);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dwrite for Aster target .\n";
		ret = -1;
		goto done;
	}

done:
	H5Dclose(aster_dataset);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
	return ret;
}


int af_GenerateOutputCumulative_AsterAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNumNoShift, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;

	// strVec_t multiVarNames = inputArgs.GetMultiVariableNames("ASTER"); // aster_MultiVars;
	std::string asterResolution = inputArgs.GetASTER_Resolution();

	// two multi-value variables are expected as this point
	strVec_t bands = inputMultiVarsMap[ASTER_BANDS];

	// data type
	hid_t asterDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t	status = H5Tset_order(asterDatatype, H5T_ORDER_LE);
	if(status < 0) {
		printf("Error: ASTER write error in H5Tset_order\n");
		return FAILED;
	}


	/*------------------------------------------
	 * create space for aster
	 */
	const int rankSpace=3;  // [bands][y][x]
	hsize_t asterDims[rankSpace];
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
	asterDims[0] = bands.size();
	asterDims[1] = trgCellNum/srcOutputWidth; // NY;
	asterDims[2] = srcOutputWidth; // NX;
	hid_t asterDataspace = H5Screate_simple(rankSpace, asterDims, NULL);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> trgCellNum: " << trgCellNum << ", srcCellNum: " << srcCellNum << "\n";
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> srcOutputWidth: " << srcOutputWidth <<  "\n";
	#endif
	int numCells;
	double *asterSingleData=NULL;

	//-----------------------------------------------------------------
	// TODO: improve by preparing these memory allocation out of loop
	// srcProcessedData , nsrcPixels
	// asterSingleData
	double * srcProcessedData = NULL;
	int * nsrcPixels = NULL;
	// Note: This is Combination case only
	for (int i=0; i< bands.size(); i++) {
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> bands[" << i << "]" << bands[i] << "\n";
		#endif

		//---------------------------------
		// read src band from BF file
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		asterSingleData = get_ast_rad(srcFile, (char*)asterResolution.c_str(), (char*)bands[i].c_str(), &numCells);
		if (asterSingleData == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get ASTER band.\n";
			return FAILED;
		}
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> numCells: " << numCells << "\n";
		#endif
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG_TIME> Read source ASTER single band data	DONE.");
		#endif

		//-------------------------------------------------
		// handle resample method
		// Note: resample should be done with trgCellNumNoShift
		srcProcessedData = new double [trgCellNumNoShift];
		//Interpolating
		std::string resampleMethod =  inputArgs.GetResampleMethod();
		std::cout << "Interpolating with '" << resampleMethod << "' method on " << inputArgs.GetSourceInstrument() << " by " << bands[i] << ".\n";
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "nnInterpolate")) {
			nnInterpolate(asterSingleData, srcProcessedData, targetNNsrcID, trgCellNumNoShift);
		}
		else if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "summaryInterpolate")) {
			nsrcPixels = new int [trgCellNumNoShift];
			summaryInterpolate(asterSingleData, targetNNsrcID, srcCellNum, srcProcessedData, nsrcPixels, trgCellNumNoShift);
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
			std::cout << "\nSource ASTER radiance MISR-base shifting...\n";
			#if DEBUG_ELAPSE_TIME
			StartElapseTime();
			#endif
			srcProcessedDataShifted = new double [widthShifted * heightShifted];
			MISRBlockOffset(srcProcessedData, srcProcessedDataShifted, (inputArgs.GetMISR_Resolution() == "L") ? 0 : 1);
			#if DEBUG_ELAPSE_TIME
			StopElapseTimeAndShow("DBG_TIME> source ASTER radiance MISR-base shift DONE.");
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
		ret = af_WriteSingleRadiance_AsterAsSrc(outputFile, asterDatatype, asterDataspace,  srcProcessedDataPtr, numCells /*processed size*/, srcOutputWidth, i /*bandIdx*/);
		if (ret == FAILED) {
			std::cerr << __FUNCTION__ << "> Error: returned fail.\n";
		}
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG_TIME> Write source ASTER single band data  DONE.");
		#endif
		//
		// TODO: buffer is not reused. make memory allocation out of get_aster_rad() to improve performance

		// free memory
		if (asterSingleData)
			free(asterSingleData);
		if(srcProcessedData)
			delete [] srcProcessedData;
		if (nsrcPixels)
			delete [] nsrcPixels;
		if(srcProcessedDataShifted)
			delete [] srcProcessedDataShifted;
	} // i loop

	H5Tclose(asterDatatype);
	H5Sclose(asterDataspace);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif

	return ret;
}

/*##############################################################
 * Test functions
 *#############################################################*/
void Test_Parser(std::string headerFile)
{
	int ret;
    AF_InputParmeterFile inputArgs;
    inputArgs.headerFileName = headerFile;
    inputArgs.ParseByLine();
    ret = inputArgs.CheckParsedValues();
    if (ret < 0) {
		std::cout << __FUNCTION__ << ":" << __LINE__ << " > Failed inputArgs.CheckParsedValues()\n";
        return;
    }

	/*-------------------------------------------------------
	 * COMMON section
	 */
	std::string inputDataPath = inputArgs.GetInputBFdataPath();
	std::cout << "TEST Parser> INPUT_FILE_PATH: " << inputDataPath << std::endl;

	std::string outputFile = inputArgs.GetOuputFilePath();
	std::cout << "TEST Parser> OUTPUT_FILE_PATH: " << outputFile << std::endl;

	std::string resampleMethod =  inputArgs.GetResampleMethod();
    std::cout << "TEST Parser> Resample Method: " << resampleMethod << std::endl;

    // Instruments
    std::string srcInstrument = inputArgs.GetSourceInstrument();
    std::string trgInstrument = inputArgs.GetTargetInstrument();
    std::cout << "TEST Parser> SOURCE instrument: " << srcInstrument << std::endl;
    std::cout << "TEST Parser> TARGET instrument: " << trgInstrument << std::endl;
	std::cout << "\n";


	/*-------------------------------------------------------
	 * MODIS section
	 */
	if(srcInstrument == "MODIS" || trgInstrument == "MODIS") {
		std::string modisResolution = inputArgs.GetMODIS_Resolution();
		std::cout << "TEST Parser>  MODIS resolution: " << modisResolution << "\n";

		strVec_t modisBands = inputArgs.GetMODIS_Bands();
		std::cout << "TEST Parser>  MODIS bands: ";
		for(int i=0; i < modisBands.size(); i++) {
			std::cout << modisBands[i] << " ";
		}
		std::cout << "\n";
		std::cout << "\n";
	}

	/*-------------------------------------------------------
	 * MISR section
	 */
	if(srcInstrument == "MISR" || trgInstrument == "MISR") {
		std::string misrResolution = inputArgs.GetMISR_Resolution();
		std::cout << "TEST Parser>  MISR resolution: " << misrResolution << "\n";

		strVec_t misrCameras = inputArgs.GetMISR_CameraAngles();
		std::cout << "TEST Parser>  MISR cameras: ";
		for(int i=0; i < misrCameras.size(); i++) {
			std::cout << misrCameras[i] << " ";
		}
		std::cout << "\n";

		strVec_t misrRadiances = inputArgs.GetMISR_Radiance();
		std::cout << "TEST Parser>  MISR radiances: ";
		for(int i=0; i < misrRadiances.size(); i++) {
			std::cout << misrRadiances[i] << " ";
		}
		std::cout << "\n";
		std::cout << "\n";
	}

	/*-------------------------------------------------------
	 * ASTER section
	 */
	if(srcInstrument == "ASTER" || trgInstrument == "ASTER") {
		std::string asterResolution = inputArgs.GetASTER_Resolution();
		std::cout << "TEST Parser>  ASTER resolution: " << asterResolution << "\n";

		strVec_t asterBands = inputArgs.GetASTER_Bands();
		std::cout << "TEST Parser>  ASTER bands: ";
		for(int i=0; i < asterBands.size(); i++) {
			std::cout << asterBands[i] << " ";
		}
		std::cout << "\n";
		std::cout << "\n";
	}

	/*-------------------------------------------------------
	 * USER_DEFINE section
	 */
	if(srcInstrument == "USER_DEFINE" || trgInstrument == "USER_DEFINE") {
	    int userEPSG = inputArgs.GetUSER_EPSG();
	    std::cout << "TEST Parser> USER EPSG: " << userEPSG << std::endl;
	    double userXmin = inputArgs.GetUSER_xMin();
	    std::cout << "TEST Parser> USER X min: " << userXmin << std::endl;
	    double userXmax = inputArgs.GetUSER_xMax();
	    std::cout << "TEST Parser> USER X max: " << userXmax << std::endl;
	    double userYmin = inputArgs.GetUSER_yMin();
	    std::cout << "TEST Parser> USER Y min: " << userYmin << std::endl;
	    double userYmax = inputArgs.GetUSER_yMax();
	    std::cout << "TEST Parser> USER Y max: " << userYmax << std::endl;
	    int userRsolution = inputArgs.GetUSER_Resolution();
	    std::cout << "TEST Parser> USER Rsolution: " << userRsolution << std::endl;
		std::cout << "\n";
	}
}


/*===========================================
 * Main 
 */

int main(int argc, char *argv[])
{
	int ret;

	if (argc < 2) {
		Usage(argc, argv);
		return FAILED;
	}

	#if 0 // TEST : parser
	Test_Parser(argv[1]);
	exit(1);
	#endif

	//----------------------------------
	// parse input parameter from file 
	std::cout << "\nUser input handling ...\n";
	AF_InputParmeterFile inputArgs;
	inputArgs.headerFileName = argv[1];
	inputArgs.ParseByLine();
	ret = inputArgs.CheckParsedValues();
	if (ret < 0) {
		return FAILED;
	}

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
	int trgCellNumNoShift;
	double* targetLatitude = NULL;
	double* targetLongitude = NULL;
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	ret = AF_GetGeolocationDataFromInstrument(trgInstrument, inputArgs, inputFile, &targetLatitude /*OUT*/, &targetLongitude /*OUT*/, trgCellNumNoShift /*OUT*/);
	if (ret == FAILED) {
		std::cerr << __FUNCTION__ << "> Error: getting geolocation data from target instrument - " << trgInstrument << ".\n";
		return FAILED;
	}
	#if DEBUG_ELAPSE_TIME
	StopElapseTimeAndShow("DBG_TIME> get target lat/long DONE.");
	#endif
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL main> trgCellNumNoShift: " <<  trgCellNumNoShift << "\n";
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

	/*
	 * Figure out trgOutputWidth and new trgCellNum by MISR target shift case
	 */
	ret = af_GetWidthAndHeightForOutputDataSize(inputArgs.GetTargetInstrument() /* target base output */, inputArgs, widthShifted, heightShifted);
	trgOutputWidth = widthShifted;
	// MISR-target & shift case update
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
		trgCellNumNew = trgCellNumNoShift;
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

	// MISR-target & shift case update
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
	 * Calculate nearest neighbor source over target geolocation
	 * Note: use not shifted trgCellNum for this
	 */
	int * targetNNsrcID = NULL;
	
	std::cout <<  "\nRunning nearest neighbor block index method... \n";
	#if DEBUG_ELAPSE_TIME
	StartElapseTime();
	#endif
	std::string resampleMethod =  inputArgs.GetResampleMethod();
	// source is low and target is similar or high resolution case (ex: MISRtoMODIS and vice versa)
	if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "nnInterpolate")) {
		targetNNsrcID = new int [trgCellNumNoShift];
		nearestNeighborBlockIndex(&srcLatitude, &srcLongitude, srcCellNum, targetLatitude, targetLongitude, targetNNsrcID, NULL, trgCellNumNoShift, 1000);
	} 
	// source is high and target is low resolution case (ex: ASTERtoMODIS)
	else if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "summaryInterpolate")) {
		targetNNsrcID = new int [srcCellNum];
		// this need to swap source and target for high resolution to low resolution case like ASTER to MODIS
		nearestNeighborBlockIndex(&targetLatitude, &targetLongitude, trgCellNumNoShift, srcLatitude, srcLongitude, targetNNsrcID, NULL, srcCellNum, 1000);
	}
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
	// Note: pass not-shifted-trgCellNum as it will internally replace if condition met
	ret = AF_GenerateTargetRadiancesOutput(inputArgs, output_file, trgCellNumNoShift, inputFile, trgInputMultiVarsMap);
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
	// Note: pass not-shifted-trgCellNum as it will internally replace if condition met
	ret = AF_GenerateSourceRadiancesOutput(inputArgs, output_file, targetNNsrcID, trgCellNumNoShift, inputFile, srcCellNum, srcInputMultiVarsMap);
	if (ret < 0) {
		std::cerr << "Error: generate source radiance output.\n";
		return FAILED;
	}
	std::cout << "Writing source radiance output done.\n";

	if (targetNNsrcID)
		delete [] targetNNsrcID;



	//==========================================
	// Closing data file
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
