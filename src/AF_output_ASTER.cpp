
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
