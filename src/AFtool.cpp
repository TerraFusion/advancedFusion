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

#include "reproject.h"
#include "gdalio.h"
#include "misrutil.h"
#include "io.h"
#include "AF_InputParmeterFile.h"
#include "AF_debug.h"
#include "AF_common.h"
#include "AF_output_util.h"
#include "AF_output_MODIS.h"
#include "AF_output_MISR.h"
#include "AF_output_ASTER.h"


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
	if (instrument == MODIS_STR ) {
		std::string resolution = inputArgs.GetMODIS_Resolution();
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> Modis resolution: " << resolution << "\n";
		#endif
		*latitude = get_modis_lat(inputFile, (char*) resolution.c_str(), &cellNum);
		if (*latitude == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get MODIS latitude.\n";
			return FAILED;
		}
		*longitude = get_modis_long(inputFile, (char*) resolution.c_str(), &cellNum);
		if (*longitude == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get MODIS longitude.\n";
			return FAILED;
		}
	}
	/*======================================================
 	 * MISR section
 	 */
	else if (instrument == MISR_STR) {
		std::string resolution = inputArgs.GetMISR_Resolution();
		#if DEBUG_TOOL
		std::cout << "DBG_TOOL " << __FUNCTION__ << "> Misr resolution: " << resolution << "\n";
		#endif
		*latitude = get_misr_lat(inputFile, (char*) resolution.c_str(), &cellNum);
		if (*latitude == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get MISR latitude.\n";
			return FAILED;
		}
		*longitude = get_misr_long(inputFile, (char*) resolution.c_str(), &cellNum);
		if (*longitude == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get MISR longitude.\n";
			return FAILED;
		}

	}
	/*======================================================
 	 * ASTER section
 	 */
	else if (instrument == ASTER_STR ) {
		//Get ASTER input parameters EX: "TIR", "ImageData10"
		std::string resolution = inputArgs.GetASTER_Resolution();
		strVec_t bands = inputArgs.GetASTER_Bands();
		// pass the first band (bands[0]) as it always exists and lat&lon is the same for a resolution.
		*latitude = get_ast_lat(inputFile, (char*) resolution.c_str(), (char*)bands[0].c_str(), &cellNum);
		if (*latitude == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get ASTER latitude.\n";
			return FAILED;
		}
		*longitude = get_ast_long(inputFile, (char*) resolution.c_str(), (char*)bands[0].c_str(), &cellNum);
		if (*longitude == NULL) {
			std::cerr << __FUNCTION__ <<  "> Error: failed to get ASTER longitude.\n";
			return FAILED;
		}
	}
	/*======================================================
	 * USER_DEFINE section
	 */
	else if (instrument == USERGRID_STR) {
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
 * Generate Target instrument radiance data to output file
 *
 */
int   AF_GenerateTargetRadiancesOutput(AF_InputParmeterFile &inputArgs, hid_t outputFile, int trgCellNum, hid_t srcFile, std::map<std::string, strVec_t> & trgInputMultiVarsMap)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret;

	std::string instrument = inputArgs.GetTargetInstrument();

	// there is no radiance for this case. just skip.
	if (instrument == USERGRID_STR) {
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
	if (instrument == MODIS_STR) {
		multiVarNames = inputArgs.GetMultiVariableNames(MODIS_STR); // modis_MultiVars;

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
	else if (instrument == MISR_STR) {

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
 * Generate Source instrument radiance data to output file
 *
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
	if (instrument == MODIS_STR) {
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
	else if (instrument == MISR_STR) {
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
	else if (instrument == ASTER_STR) {
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
	if(srcInstrument == MODIS_STR || trgInstrument == MODIS_STR) {
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
	if(srcInstrument == MISR_STR || trgInstrument == MISR_STR) {
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
	if(srcInstrument == ASTER_STR || trgInstrument == ASTER_STR) {
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
	if(srcInstrument == USERGRID_STR || trgInstrument == USERGRID_STR) {
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
	if(inputArgs.GetMISR_Shift() == "ON" && inputArgs.GetTargetInstrument() == MISR_STR) {
		std::cout << "Target latitude MISR-base block unstacking...\n";
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		targetLatitudeShifted = (double *) malloc(sizeof(double) * widthShifted * heightShifted);
		std::string misrResolution = inputArgs.GetMISR_Resolution();
		MISRBlockOffset<double>(targetLatitude, targetLatitudeShifted, (misrResolution == "L") ? 0 : 1);
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
	if(inputArgs.GetMISR_Shift() == "ON" && inputArgs.GetTargetInstrument() == MISR_STR) {
		std::cout << "Target longitude MISR-base block unstacking...\n";
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		targetLongitudeShifted = (double *) malloc(sizeof(double) * widthShifted * heightShifted);
		std::string misrResolution = inputArgs.GetMISR_Resolution();
		MISRBlockOffset<double>(targetLongitude, targetLongitudeShifted, (misrResolution == "L") ? 0 : 1);
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
		double maxRadius = inputArgs.GetMaxRadiusForNNeighborFunc(srcInstrument);
		nearestNeighborBlockIndex(&srcLatitude, &srcLongitude, srcCellNum, targetLatitude, targetLongitude, targetNNsrcID, NULL, trgCellNumNoShift, maxRadius);
	} 
	// source is high and target is low resolution case (ex: ASTERtoMODIS)
	else if (inputArgs.CompareStrCaseInsensitive(resampleMethod, "summaryInterpolate")) {
		// when summaryInterpolate is used, need to swap source and target. This is cases for projecting high resolution to low resolution case like ASTER to MODIS
		targetNNsrcID = new int [srcCellNum];
		// get it from src instrument of nearestNeighbor point of view, which is switched for this case, thus use target instrument.
		double maxRadius = inputArgs.GetMaxRadiusForNNeighborFunc(trgInstrument);
		nearestNeighborBlockIndex(&targetLatitude, &targetLongitude, trgCellNumNoShift, srcLatitude, srcLongitude, targetNNsrcID, NULL, srcCellNum, maxRadius);
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
