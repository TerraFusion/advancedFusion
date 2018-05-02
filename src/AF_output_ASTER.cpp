
/*********************************************************************
 * DESCRIPTION:
 *  Generate radiance data output to HDF5 for ASTER
 *
 *
 * DEVELOPERS:
 *  - Jonathan Kim (jkm@illinois.edu)
 */

#include "AF_output_ASTER.h"

#include "AF_common.h"
#include "AF_output_util.h"
#include "io.h"
#include "reproject.h"
#include "misrutil.h"



/*===============================================================================
 *
 * ASTER as Source instrument, functions to generate radiance data
 *
 *===============================================================================*/

// T: type of data type of output data
template <typename T>
static int af_WriteSingleRadiance_AsterAsSrc(hid_t outputFile, std::string outputDsetName, hid_t dataTypeH5, hid_t fileSpaceH5, T* processedData, int trgCellNum, int outputWidth, int bandIdx)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;
	herr_t status;
	hid_t aster_dataset;
	std::string dsetPath = SRC_DATA_GROUP + "/" + outputDsetName;

	if(bandIdx==0) { // means new
		aster_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), dataTypeH5, fileSpaceH5,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
	// select memory space
	int ranksMem=2;
	hsize_t dim2dMem[ranksMem];
	hsize_t start2dMem[ranksMem];
	hsize_t count2dMem[ranksMem];
	dim2dMem[0] = trgCellNum/outputWidth; // y
	dim2dMem[1] = outputWidth; // x
	hid_t memSpaceH5 = H5Screate_simple(ranksMem, dim2dMem, NULL);

	start2dMem[0] = 0; // y
	start2dMem[1] = 0; // x
	count2dMem[0] = trgCellNum/outputWidth; // y
	count2dMem[1] = outputWidth; // x

	status = H5Sselect_hyperslab(memSpaceH5, H5S_SELECT_SET, start2dMem, NULL, count2dMem, NULL);

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

	status = H5Sselect_hyperslab(fileSpaceH5, H5S_SELECT_SET, startFile, NULL, countFile, NULL);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Sselect_hyperslab for Aster target .\n";
		ret = -1;
		goto done;
	}

	status = H5Dwrite(aster_dataset, dataTypeH5, memSpaceH5, fileSpaceH5, H5P_DEFAULT, processedData);
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

	//-----------------------------------
	// define data types for hdf5 data
	hid_t dataTypeDoubleH5 = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t	status = H5Tset_order(dataTypeDoubleH5, H5T_ORDER_LE);
	if(status < 0) {
		printf("Error: ASTER write error in H5Tset_order\n");
		return FAILED;
	}

	hid_t dataTypeIntH5 = H5Tcopy(H5T_NATIVE_INT);
	status = H5Tset_order(dataTypeIntH5, H5T_ORDER_LE);
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
	// srcProcessedData, SD, nsrcPixels
	// asterSingleData
	double * srcProcessedData = NULL; // radiance
	double * SD = NULL;  // Standard Deviation
	int * nsrcPixels = NULL; // count
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
			SD = new double [trgCellNumNoShift];
			nsrcPixels = new int [trgCellNumNoShift];
			summaryInterpolate(asterSingleData, targetNNsrcID, srcCellNum, srcProcessedData, SD, nsrcPixels, trgCellNumNoShift);
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
			// use srcProcessedDataShifted instead of srcProcessedData, and free memory
			if(srcProcessedData) {
				delete [] srcProcessedData;
				srcProcessedData = NULL;
			}
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
		// output radiance dset
		ret = af_WriteSingleRadiance_AsterAsSrc<double>(outputFile, ASTER_RADIANCE_DSET, dataTypeDoubleH5, asterDataspace,  srcProcessedDataPtr, numCells /*processed size*/, srcOutputWidth, i /*bandIdx*/);
		if (ret == FAILED) {
			std::cerr << __FUNCTION__ << "> Error: returned fail.\n";
		}

		// output standard deviation dset
		ret = af_WriteSingleRadiance_AsterAsSrc<double>(outputFile, ASTER_SD_DSET, dataTypeDoubleH5, asterDataspace,  SD, numCells /*processed size*/, srcOutputWidth, i /*bandIdx*/);
		if (ret == FAILED) {
			std::cerr << __FUNCTION__ << "> Error: returned fail.\n";
		}

		// output pixels count dset
		ret = af_WriteSingleRadiance_AsterAsSrc<int>(outputFile, ASTER_COUNT_DSET, dataTypeIntH5, asterDataspace,  nsrcPixels, numCells /*processed size*/, srcOutputWidth, i /*bandIdx*/);
		if (ret == FAILED) {
			std::cerr << __FUNCTION__ << "> Error: returned fail.\n";
		}
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG_TIME> Write source ASTER data (randiance, SD, count) of single band DONE.");
		#endif


		// TODO: buffer is not reused. make memory allocation out of get_aster_rad() to improve performance

		// free memory
		if (asterSingleData)
			free(asterSingleData);
		if(srcProcessedData)
			delete [] srcProcessedData;
		if(SD)
			delete [] SD;
		if (nsrcPixels)
			delete [] nsrcPixels;
		if(srcProcessedDataShifted)
			delete [] srcProcessedDataShifted;
	} // i loop

	H5Tclose(dataTypeDoubleH5);
	H5Sclose(asterDataspace);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif

	return ret;
}
