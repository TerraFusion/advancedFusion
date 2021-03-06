
/*********************************************************************
 * DESCRIPTION:
 *	Generate radiance data output to HDF5 for ASTER
 *
 *
 * DEVELOPERS:
 *	- Jonathan Kim (jkm@illinois.edu)
 *	- Hyo-Kyung (Joe) Lee (hyoklee@hdfgroup.org)
	 \author Hyo-Kyung (Joe) Lee (hyoklee@hdfgroup.org)
	 \date May 18, 2018
	 \note removed _FillValue attribute for ASTER_Count in af_WriteSingleRadiance_AsterAsSrc().
	 \date May 17, 2018
	 \note added CF attributes and cleaned up indentation in in af_WriteSingleRadiance_AsterAsSrc().
 */

#include "AF_output_ASTER.h"

#include "AF_common.h"
#include "AF_output_util.h"
#include "io.h"
#include "reproject.h"
#include "misrutil.h"
#include <iostream>
#include <string>



/*#############################################################################
 *
 * ASTER as Source instrument, functions to generate radiance data
 *
 *############################################################################*/

/*=====================================================================
 * DESCRIPTION:
 *	 Write resampled radiance output data of a single orbit of the given
 *	 band for ASTER as the source instrument.
 *	 Only called by af_GenerateOutputCumulative_AsterAsSrc().
 *
 * PARAMETER:
 *	- outputFile : HDF5 id for output file
 *	- outputDsetName : HDF5 output dset name
 *	- dataTypeH5 : HDF5 id for output datatype 
 *	- fileSpaceH5 : HDF5 id file sapce 
 *	- processedData : resmapled data pointer
 *	- trgCellNum : number of cells (pixels) in target instrument data
 *	- outputWidth : cross-track (width) size for output image
 *	- bandIdx : ASTER band index.
 * 
 * RETURN:
 *	- Success: SUCCEED	(defined in AF_common.h)
 *	- Fail : FAILED  (defined in AF_common.h)
 *
 * NOTE:
 *	- TODO: if radiance data becomes all float internally, use float directly
 *	  without converting via HDF5 
 */
// T_IN : input data type
// T_OUT : output data type
template <typename T_IN, typename T_OUT>
static int af_WriteSingleRadiance_AsterAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, std::string outputDsetName, hid_t dataTypeH5, hid_t fileSpaceH5, T_IN* processedData, int trgCellNum, int outputWidth, int bandIdx,const strVec_t bands, hid_t ctrackDset,hid_t atrackDset,hid_t bandDset)
{
#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	std::cout << "ASTER bands: ";
	for(int i=0; i < bands.size(); i++) {
		std::cout << bands[i] << " ";
	}
	std::cout << std::endl;

#endif

	int ret = SUCCEED;
	herr_t status;
	hid_t aster_dataset;
	std::string dsetPath = SRC_DATA_GROUP + "/" + outputDsetName;

	/*-------------------------------------
	 * set output data type
	 */
	hid_t dataTypeOutH5;
	if (std::is_same<T_OUT, float>::value) {
		dataTypeOutH5 = H5T_IEEE_F32LE;
	}
	else if (std::is_same<T_OUT, double>::value) {
		dataTypeOutH5 = H5T_IEEE_F64LE;
	}
	else if (std::is_same<T_OUT, int>::value) {
		dataTypeOutH5 = H5T_NATIVE_INT;
	}
	else {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: invlid output data type T_OUT specified." << std::endl;
		return FAILED;
	}

	/*-------------------------------------
	 * if first time, create dataset
	 * otherwise, open existing one
	 */
	if(bandIdx==0) { // means new
		if(true == inputArgs.GetUseH5Chunk()) {
			hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
			if(create_chunk_comp_plist(plist_id,3,(size_t)trgCellNum,(size_t)outputWidth)<0) {
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: create_chunk_comp_plist failed.\n";
				H5Pclose(plist_id);
				return FAILED;
			}
			aster_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), dataTypeOutH5, fileSpaceH5,H5P_DEFAULT, plist_id, H5P_DEFAULT);
			H5Pclose(plist_id);

		}
		else 
			aster_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), dataTypeOutH5, fileSpaceH5,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(aster_dataset < 0) {
			std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dcreate2 target data in output file.\n";
			return FAILED;
		}
		else {
			/*
			 * make compatible with CF convention (NetCDF)
			 */

			// Attach dimension scales.
			// cross track
			if(H5DSattach_scale(aster_dataset,ctrackDset,2)<0) {
				H5Dclose(aster_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for ASTER cross-track dimension.\n";
				return FAILED;
			}

			// along track
			if(H5DSattach_scale(aster_dataset,atrackDset,1)<0) {
				H5Dclose(aster_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for ASTER along-track dimension.\n";
				return FAILED;
			}

			// band
			if(H5DSattach_scale(aster_dataset,bandDset,0)<0) {
				H5Dclose(aster_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for ASTER band dimension.\n";
				return FAILED;
			}
		
			// Change units based on dset name.
			char* units = NULL;
			float _FillValue = -999.0;
			float valid_min = 0.;
			float valid_max = 0.;
			unsigned short handle_flag = 0;
			if (outputDsetName == "ASTER_Radiance") {
				units = "Watts/m^2/micrometer/steradian";
				valid_min = 0.;
				valid_max = 569.0;
			}
			if (outputDsetName == "ASTER_Count" || outputDsetName == "ASTER_SD") {
				handle_flag = 1;
			}
			// Don't add valid_min attribute by making valid_min argument
			// same as _FillValue.
			if(af_write_cf_attributes(aster_dataset, units, _FillValue, valid_min,valid_max,handle_flag) < 0) {
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<	"> Error: af_write_cf_attributes" << std::endl;
				return FAILED;
			}

			// Add CF long name
			const char* long_name = "long_name";
			if(outputDsetName == "ASTER_Count") {
				if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),long_name,"Number of ASTER pixels in a resampled cell")<0) {
					H5Dclose(aster_dataset);
					std::cerr << __FUNCTION__ << ":" << __LINE__ <<	"> Error: cannot generate long_name attribute for ASTER_Count" << std::endl;
					return FAILED;
				}
			}
			else if(outputDsetName == "ASTER_SD") {
				if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),long_name,"Standard deviation of ASTER pixels in a resampled cell")<0) {
					H5Dclose(aster_dataset);
					std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate long_name attribute for ASTER_SD" << std::endl;
					return FAILED;
				}
			}
			else {

				// Write long_name 
				std::string aster_resolution = inputArgs.GetASTER_Resolution();
				std::string long_name_value; 
				if(aster_resolution == "TIR") {
					std::string tir = "Thermal Infrared ";
					//long_name_value = "ASTER " + "Thermal Infrared ";
					long_name_value = "ASTER Level 1T " + tir;
					long_name_value += "Radiance";
				}
				else if(aster_resolution == "SWIR") {
					std::string swir = "Short-wave Infrared ";
					long_name_value = "ASTER Level 1T " + swir;
					long_name_value += "Radiance";
				}
				else {
					std::string vnir ="Visible and Near Infrared ";
					long_name_value = "ASTER Level 1T " +vnir;
					long_name_value += "Radiance";
				}

	 			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),long_name,long_name_value.c_str())<0) {
					H5Dclose(aster_dataset);
					std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate long_name attribute for ASTER Radiance" << std::endl;
					return FAILED;
				}
 				
				std::vector<std::string> band_name_vec = inputArgs.GetASTER_Orig_Bands();
				//std::string band_name_values = std::accumulate(band_name_vec.begin(),band_name_vec.end(),std::string(",")); 
//#if 0
				std::string band_name_values;
				for(int i = 0; i<band_name_vec.size();i++) 
					band_name_values = band_name_values + band_name_vec[i] +',';
				
//#endif
				band_name_values = band_name_values.erase(band_name_values.size()-1,1);
				if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"band_names",band_name_values.c_str())<0) {
					H5Dclose(aster_dataset);
					std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate band_names" << std::endl;
					return FAILED;
				}

				// Add resolution as a number.
				float aster_resolution_value = inputArgs.GetInstrumentResolutionValue(ASTER_STR);
				if(false == af_AddSrcSpatialResolutionAttrs(outputFile,dsetPath,aster_resolution_value,true)) {
					H5Dclose(aster_dataset);
					return FAILED;
				}

				if(inputArgs.GetTargetInstrument() != "USER_DEFINE") {
					float target_resolution_value = inputArgs.GetInstrumentResolutionValue(inputArgs.GetTargetInstrument());
					if(false == af_AddSrcSpatialResolutionAttrs(outputFile,dsetPath,target_resolution_value,false)) {
						H5Dclose(aster_dataset);
						return FAILED;
					}

					if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"spatial_resolution_resampled_units","meter")<0) {
						H5Dclose(aster_dataset);
						std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
						return FAILED;
					}
				}
				if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"spatial_resolution_original_units","meter")<0) {
					H5Dclose(aster_dataset);
					std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
					return FAILED;
				}

				if("USER_DEFINE" == inputArgs.GetTargetInstrument()) {
				
					std::string comments_ud = "Check the group attributes under group /Geolocation to find the user-defined ESPG code and resolution information.";
					if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"spatial_resolution_resampled_description",comments_ud.c_str())<0) {
						H5Dclose(aster_dataset);
						std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
						return FAILED;
					}	
				}

				// Add resample method
				std::string resample_method_value = "Summary Interpolation";
				if(inputArgs.GetResampleMethod()=="nnInterpolate")
					resample_method_value = "Nearest Neighbor Interpolation";

				if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"resample_method",resample_method_value.c_str())<0) {
					H5Dclose(aster_dataset);
					std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate resample_method attr." << std::endl;
					return FAILED;
				}
			
			}

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
		ret = FAILED;
		goto done;
	}

	status = H5Dwrite(aster_dataset, dataTypeH5, memSpaceH5, fileSpaceH5, H5P_DEFAULT, processedData);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dwrite for Aster target .\n";
		ret = FAILED;
		goto done;
	}

 done:
	H5Dclose(aster_dataset);


	if (outputDsetName == "ASTER_Radiance") {
		bool output_geotiff = inputArgs.GetGeoTiffOutput();

		if(true == output_geotiff && ("USER_DEFINE" == inputArgs.GetTargetInstrument())) {
			std::string op_geotiff_fname = get_gtiff_fname(inputArgs,-1,bandIdx);
			//std::string op_geotiff_fname = "test.tif";
			int userOutputEPSG = inputArgs.GetUSER_EPSG();
			double userXmin = inputArgs.GetUSER_xMin();
			double userXmax = inputArgs.GetUSER_xMax();
			double userYmin = inputArgs.GetUSER_yMin();
			double userYmax = inputArgs.GetUSER_yMax();
			double userResolution = inputArgs.GetUSER_Resolution();
			gdalIORegister();
			
			writeGeoTiff((char*)op_geotiff_fname.c_str(),(double*)processedData, userOutputEPSG, userXmin, userYmin, userXmax, userYmax, userResolution);
		}
	}

#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
#endif
	return ret;
}


/*=====================================================================
 * DESCRIPTION:
 *	 Write resampled radiance output data of a single orbit for all the 
 *	 specified bands for ASTER as the source instrument.
 *
 * PARAMETER:
 *	- inputArgs : a class object contains all the user input parameter info
 *	- outputFile : HDF5 id for output file
 *	- targetNNsrcID : got from nearestNeighborBlockIndex()
 *	- trgCellNumNoShift : number of target instrument data cells before
 *	  applying shift (if MISR is target)
 *	- srcFile : HDF5 id for input file
 *	- srcCellNum : number of source instrument data cells
 *	- inputMultiVarsMap : To obtain multiple values from a given user input
 *	  directive which allows to have multiple values.
 * 
 * RETURN:
 *	- Success: SUCCEED	(defined in AF_common.h)
 *	- Fail : FAILED  (defined in AF_common.h)
 */
int af_GenerateOutputCumulative_AsterAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNumNoShift, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap,hid_t ctrackDset,hid_t atrackDset)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;

	// strVec_t multiVarNames = inputArgs.GetMultiVariableNames(ASTER_STR); // aster_MultiVars;
	std::string asterResolution = inputArgs.GetASTER_Resolution();

	// two multi-value variables are expected as this point
	strVec_t bands = inputMultiVarsMap[ASTER_BANDS];

	// Create band dimension 
	hid_t bandDset = create_pure_dim_dataset(outputFile,(hsize_t)(bands.size()),"Band_ASTER");
	if(bandDset < 0) {
		printf("create_pure_dim_dataset for ASTER band failed.\n");
		return FAILED;
	}
    
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
	const int rankSpace=3;	// [bands][y][x]
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
	if(inputArgs.GetMISR_Shift() == "ON" && inputArgs.GetTargetInstrument() == MISR_STR) {
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
	// srcProcessedData, SD, srcPixelCount
	// asterSingleData
	double * srcProcessedData = NULL; // radiance
	double * SD = NULL;  // Standard Deviation
	int * srcPixelCount = NULL; // count
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
			srcPixelCount = new int [trgCellNumNoShift];
			summaryInterpolate(asterSingleData, targetNNsrcID, srcCellNum, srcProcessedData, SD, srcPixelCount, trgCellNumNoShift);
			#if 0 // DEBUG_TOOL
			std::cout << "DBG_TOOL> No nodata values: \n";
			for(int i = 0; i < trgCellNumNoShift; i++) {
				if(srcPixelCount[i] > 0) {
					printf("%d,\t%lf\n", srcPixelCount[i], srcProcessedData[i]);
				}
			}
			#endif
		}
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG> nnInterpolate  DONE.");
		#endif


		//-----------------------------------------------------------------------
		// check if need to shift by MISR (shift==ON & target) case before writing
		// radiance data
		double * srcRadianceDataShifted = NULL;
		double * srcRadianceDataPtr = NULL;
		// standard deviation data
		double * srcSDDataShifted = NULL;
		double * srcSDDataPtr = NULL;
		// pixel count data
		int * srcPixelCountDataShifted = NULL;
		int * srcPixelCountDataPtr = NULL;

		// if MISR is target and Shift is On
		if(inputArgs.GetMISR_Shift() == "ON" && inputArgs.GetTargetInstrument() == MISR_STR) {
			std::cout << "\nSource ASTER radiance MISR-base shifting...\n";
			#if DEBUG_ELAPSE_TIME
			StartElapseTime();
			#endif
			/*-------------------- 
			 * shift radiance data
			 */
			srcRadianceDataShifted = new double [widthShifted * heightShifted];
			MISRBlockOffset<double>(srcProcessedData, srcRadianceDataShifted, (inputArgs.GetMISR_Resolution() == "L") ? 0 : 1);
			// use srcRadianceDataShifted instead of srcProcessedData, and free memory
			if(srcProcessedData) {
				delete [] srcProcessedData;
				srcProcessedData = NULL;
			}

			/*-------------------- 
			 * shift SD data
			 */
			srcSDDataShifted = new double [widthShifted * heightShifted];
			MISRBlockOffset<double>(SD, srcSDDataShifted, (inputArgs.GetMISR_Resolution() == "L") ? 0 : 1);
			// use srcSDDataShifted instead of SD, and free memory
			if(SD) {
				delete [] SD;
				SD = NULL;
			}

			/*-------------------- 
			 * shift PixelCount data
			 */
			srcPixelCountDataShifted = new int [widthShifted * heightShifted];
			MISRBlockOffset<int>(srcPixelCount, srcPixelCountDataShifted, (inputArgs.GetMISR_Resolution() == "L") ? 0 : 1);
			// use srcPixelCountDataShifted instead of PixelCount, and free memory
			if(srcPixelCount) {
				delete [] srcPixelCount;
				srcPixelCount = NULL;
			}
			#if DEBUG_ELAPSE_TIME
			StopElapseTimeAndShow("DBG_TIME> source ASTER radiance MISR-base shift DONE.");
			#endif

			srcRadianceDataPtr = srcRadianceDataShifted;
			srcSDDataPtr = srcSDDataShifted;
			srcPixelCountDataPtr = srcPixelCountDataShifted;
			numCells = widthShifted * heightShifted;
		}
		else { // dats with no misr-trg shift
			srcRadianceDataPtr = srcProcessedData;
			srcSDDataPtr = SD;
			srcPixelCountDataPtr = srcPixelCount;	
			numCells = trgCellNum;
		}

		//---------------------------------
		// write src band to AF file
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		// output radiance dset
		ret = af_WriteSingleRadiance_AsterAsSrc<double, float>(inputArgs,outputFile, ASTER_RADIANCE_DSET, dataTypeDoubleH5, asterDataspace,  srcRadianceDataPtr, numCells /*processed size*/, srcOutputWidth, i /*bandIdx*/,bands,ctrackDset,atrackDset,bandDset);
		if (ret == FAILED) {
			std::cerr << __FUNCTION__ << "> Error: returned fail.\n";
		}

		// output standard deviation dset
		ret = af_WriteSingleRadiance_AsterAsSrc<double, float>(inputArgs,outputFile, ASTER_SD_DSET, dataTypeDoubleH5, asterDataspace,  srcSDDataPtr, numCells /*processed size*/, srcOutputWidth, i /*bandIdx*/,bands,ctrackDset,atrackDset,bandDset);
		if (ret == FAILED) {
			std::cerr << __FUNCTION__ << "> Error: returned fail.\n";
		}

		// output pixels count dset
		ret = af_WriteSingleRadiance_AsterAsSrc<int, int>(inputArgs,outputFile, ASTER_COUNT_DSET, dataTypeIntH5, asterDataspace,	srcPixelCountDataPtr, numCells /*processed size*/, srcOutputWidth, i /*bandIdx*/,bands,ctrackDset,atrackDset,bandDset);
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
		if (srcPixelCount)
			delete [] srcPixelCount;
		if(srcRadianceDataShifted)
			delete [] srcRadianceDataShifted;
		if(srcPixelCountDataShifted)
			delete [] srcPixelCountDataShifted;
		if(srcSDDataShifted)
			delete [] srcSDDataShifted;
	} // i loop

	H5Tclose(dataTypeDoubleH5);
	H5Sclose(asterDataspace);
	H5Dclose(bandDset);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif

	return ret;
}
