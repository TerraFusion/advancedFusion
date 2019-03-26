/*********************************************************************
 * DESCRIPTION:
 *  Generate radiance data output to HDF5 for MISR
 *
 *
 * DEVELOPERS:
 *  - Jonathan Kim (jkm@illinois.edu)
 *  - Hyo-Kyung (Joe) Lee (hyoklee@hdfgroup.org)
      \author Hyo-Kyung (Joe) Lee (hyoklee@hdfgroup.org)
      \date May 18, 2018
      \note added CF attributes in af_WriteSingleRadiance_MisrAsTrg().
      \date May 15, 2018
      \note added CF attributes in af_WriteSingleRadiance_MisrAsSrc().
 */

#include "AF_output_MISR.h"

#include "AF_common.h"
#include "AF_output_util.h"
#include "io.h"
#include "reproject.h"
#include "misrutil.h"


/*#############################################################################
 *
 * MISR as Target instrument, functions to generate radiance data
 *
 *############################################################################*/

/*=========================================================================
 * DESCRIPTION:
 *   Write radiance output data of a single orbit of the given band for
 *   MISR as the target instrument.
 *   Only called by af_GenerateOutputCumulative_MisrAsTrg().                              
 *
 * PARAMETER:
 *  - outputFile : HDF5 id for output file
 *  - dataTypeH5 : HDF5 id for output datatype
 *  - fileSpaceH5 : HDF5 id file sapce
 *  - misrData : MISR data pointer (not processed)
 *  - misrDataSize : number of cells (pixels) in MISR data
 *  - outputWidth : cross-track (width) size for output image
 *  - cameraIdx : MISR camera index
 *  - radianceIdx : MISR radiance index
 *
 * RETURN:
 *  - Success: SUCCEED  (defined in AF_common.h)
 *  - Fail : FAILED  (defined in AF_common.h)
 *
 * NOTE:
 *  - TODO: if radiance data becomes all float internally, use float directly
 *    without converting via HDF5
 */
// T_IN : input data type
// T_OUT : output data type
template <typename T_IN, typename T_OUT>
static int af_WriteSingleRadiance_MisrAsTrg(AF_InputParmeterFile &inputArgs, hid_t outputFile, hid_t dataTypeH5, hid_t fileSpaceH5, T_IN* misrData, int misrDataSize, int outputWidth, int cameraIdx, int radianceIdx, hid_t ctrackDset,hid_t atrackDset,hid_t cameraDset,hid_t bandDset)
{

	int ret = SUCCEED;

	herr_t status;
	hid_t misr_dataset;
	std::string dsetPath = TRG_DATA_GROUP + "/" + MISR_RADIANCE_DSET;


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
	if( (cameraIdx + radianceIdx) ==0 ) { // means new

		if(true == inputArgs.GetUseH5Chunk()) {
			hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
			if(create_chunk_comp_plist(plist_id,4,(size_t)misrDataSize,(size_t)outputWidth)<0) {
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: create_chunk_comp_plist failed.\n";
				H5Pclose(plist_id);
				return FAILED;
			}
			misr_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), dataTypeOutH5, fileSpaceH5,H5P_DEFAULT, plist_id, H5P_DEFAULT);
			H5Pclose(plist_id);
		}
		else
			misr_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), dataTypeOutH5, fileSpaceH5,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(misr_dataset < 0) {
			std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dcreate2 target data in output file.\n";
			return FAILED;
		}
		else {
			// make compatible with CF convention (NetCDF)
			// attach dimensions
			
			// cross track
			if(H5DSattach_scale(misr_dataset,ctrackDset,3)<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MISR cross-track dimension.\n";
				return FAILED;
			}

			// along track
			if(H5DSattach_scale(misr_dataset,atrackDset,2)<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MISR along-track dimension.\n";
				return FAILED;
			}

			// band
			if(H5DSattach_scale(misr_dataset,bandDset,1)<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MISR band dimension.\n";
				return FAILED;
			}

			//cameras
			if(H5DSattach_scale(misr_dataset,cameraDset,0)<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MISR cameras dimension.\n";
				return FAILED;
			}

			char* units = "Watts/m^2/micrometer/steradian";
			float valid_min = 0;
			float valid_max = 800.0;
			float _FillValue = -999.0;
			if(af_write_cf_attributes(misr_dataset, units, _FillValue, valid_min,valid_max,0) < 0) {
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: af_write_cf_attributes" << std::endl;                            
				return FAILED;                        
			}

			// Add CF long_name 
			const char* long_name = "long_name";
 			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),long_name,"MISR Level 1B2 Radiance")<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate long_name attribute for MISR Radiance" << std::endl;
				return FAILED;
			}
//#if 0 				
			//std::string misr_resolution = inputArgs.GetMISR_Resolution();
			std::string camera_angle_values;
			std::vector<std::string> camera_angle_vec = inputArgs.GetMISR_CameraAngles();
			for(int i = 0; i<camera_angle_vec.size();i++) 
				camera_angle_values = camera_angle_values + camera_angle_vec[i] +',';
			camera_angle_values = camera_angle_values.erase(camera_angle_values.size()-1,1);

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"camera_names",camera_angle_values.c_str())<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate camera_angles" << std::endl;
				return FAILED;
			}

			std::string radiance_values;
			std::vector<std::string> radiance_vec = inputArgs.GetMISR_Radiance();
			for(int i = 0; i<radiance_vec.size();i++) 
				radiance_values = radiance_values + radiance_vec[i] +',';
			radiance_values = radiance_values.erase(radiance_values.size()-1,1);
			

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"band_names",radiance_values.c_str())<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate radiances" << std::endl;
				return FAILED;
			}

			// Add resolution as a number.
			float misr_resolution_value = inputArgs.GetInstrumentResolutionValue(MISR_STR);
			if(H5LTset_attribute_float(outputFile,dsetPath.c_str(),"spatial_resolution",&misr_resolution_value, 1 ) <0) {
                std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate resolution attribute." << std::endl;
				H5Dclose(misr_dataset);
                return FAILED;
            }

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"spatial_resolution_units","meter")<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
				return FAILED;
			}

			// Add unstacked 
			std::string misr_target_block_unstack_value = inputArgs.GetMISR_Shift();
			//if(inputArgs.GetResampleMethod()=="nnInterpolate")
			//		resample_method_value = "Nearest Neighbor Interpolation";

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"MISR_target_block_unstack",misr_target_block_unstack_value.c_str())<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate MISR target block unstack attr." << std::endl;
				return FAILED;

			}


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

	status = H5Sselect_hyperslab(fileSpaceH5, H5S_SELECT_SET, star3dFile, NULL, count3dFile, NULL);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Sselect_hyperslab for Misr target .\n";
		ret = -1;
		goto done;
	}

	status = H5Dwrite(misr_dataset, dataTypeH5, memspace, fileSpaceH5, H5P_DEFAULT, misrData);
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


/*=========================================================================
 * DESCRIPTION:
 *   Write radiance output data of a single orbit for all the
 *   specified bands for MISR as the target instrument.
 *
 * PARAMETER:
 *  - inputArgs : a class object contains all the user input parameter info
 *  - outputFile : HDF5 id for output file
 *  - srcFile : HDF5 id for input file
 *  - trgCellNumOri : number of target instrument data cells (un-modified)
 *  - inputMultiVarsMap : To obtain multiple values from a given user input
 *    directive which allows to have multiple values.
 *
 * RETURN:
 *  - Success: SUCCEED  (defined in AF_common.h)
 *  - Fail : FAILED  (defined in AF_common.h)
 */
int af_GenerateOutputCumulative_MisrAsTrg(AF_InputParmeterFile &inputArgs, hid_t outputFile,hid_t srcFile, int trgCellNumOri, std::map<std::string, strVec_t> &inputMultiVarsMap,hid_t ctrackDset, hid_t atrackDset)
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

	// Create band and camera dimensions
	hid_t bandDset = create_pure_dim_dataset(outputFile,(hsize_t)(radiances.size()),"Band_MISR");
	if(bandDset < 0) {
		printf("create_pure_dim_dataset for MISR band failed.\n");
		return FAILED;
	}
	hid_t cameraDset = create_pure_dim_dataset(outputFile,(hsize_t)(cameras.size()),"Camera");
	if(bandDset < 0) {
		printf("create_pure_dim_dataset for MISR Camera failed.\n");
		return FAILED;
	}

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
	ret = af_GetWidthAndHeightForOutputDataSize(MISR_STR, inputArgs, widthShifted, heightShifted);
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
			StopElapseTimeAndShow("DBG_TIME> Read target MISR single band data  DONE.");
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
				MISRBlockOffset<double>(misrSingleData, misrSingleDataShifted, (misrResolution == "L") ? 0 : 1);
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
			af_WriteSingleRadiance_MisrAsTrg<double, float>(inputArgs,outputFile, misrDatatype, misrDataspace,  misrSingleDataPtr, numCells, targetOutputWidth, i, j,ctrackDset,atrackDset,cameraDset,bandDset);
			#if DEBUG_ELAPSE_TIME
			StopElapseTimeAndShow("DBG_TIME> Write target MISR single band data  DONE.");
			#endif

			// TODO: buffer is not reused. make memory allocation out of get_misr_rad() to improve performance
			if (misrSingleData)
				free(misrSingleData);

			if (misrSingleDataShifted)
				free(misrSingleDataShifted);
		}
	}

	H5Dclose(cameraDset);
	H5Dclose(bandDset);
	H5Tclose(misrDatatype);
	H5Sclose(misrDataspace);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif

	return SUCCEED;
}



/*#############################################################################
 *
 * MISR as Source instrument, functions to generate radiance data
 *
 *############################################################################*/

/*=====================================================================
 * DESCRIPTION:
 *   Write resampled radiance output data of a single orbit of the given
 *   camera & radiance for MISR as the source instrument.
 *   Only called by af_GenerateOutputCumulative_MisrAsSrc().
 *
 * PARAMETER:
 *  - outputFile : HDF5 id for output file
 *  - dataTypeH5 : HDF5 id for output datatype
 *  - fileSpaceH5 : HDF5 id file sapce
 *  - processedData : resmapled data pointer
 *  - trgCellNum : number of cells (pixels) in target instrument data
 *  - outputWidth : cross-track (width) size for output image
 *  - cameraIdx : MISR camera index
 *  - radIdx : MISR radiance index
 *
 * RETURN:
 *  - Success: SUCCEED  (defined in AF_common.h)
 *  - Fail : FAILED  (defined in AF_common.h)
 *
 * NOTE:
 *  - TODO: if radiance data becomes all float internally, use float directly
 *    without converting via HDF5
 */
// T_IN : input data type
// T_OUT : output data type
template <typename T_IN, typename T_OUT>
static int af_WriteSingleRadiance_MisrAsSrc(AF_InputParmeterFile &inputArgs,hid_t outputFile, hid_t dataTypeH5, hid_t fileSpaceH5, T_IN* processedData, int trgCellNum, int outputWidth, int cameraIdx, int radIdx, hid_t ctrackDset,hid_t atrackDset,hid_t cameraDset,hid_t bandDset)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;
	herr_t status;
	hid_t misr_dataset;
	std::string dsetPath = SRC_DATA_GROUP + "/" + MISR_RADIANCE_DSET;


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
	if(cameraIdx==0 && radIdx==0) { // means new
		if(true == inputArgs.GetUseH5Chunk()) {
			hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
			if(create_chunk_comp_plist(plist_id,4,(size_t)trgCellNum,(size_t)outputWidth)<0) {
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: create_chunk_comp_plist failed.\n";
				H5Pclose(plist_id);
				return FAILED;
			}
			misr_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), dataTypeOutH5, fileSpaceH5,H5P_DEFAULT, plist_id, H5P_DEFAULT);
			H5Pclose(plist_id);
		}
		else 
			misr_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), dataTypeOutH5, fileSpaceH5,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(misr_dataset < 0) {
			std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dcreate2 target data in output file.\n";
			return FAILED;
		}
		else {
			// make compatible with CF convention (NetCDF)
			// attach dimensions
			
			// cross track
			if(H5DSattach_scale(misr_dataset,ctrackDset,3)<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MISR cross-track dimension.\n";
				return FAILED;
			}

			// along track
			if(H5DSattach_scale(misr_dataset,atrackDset,2)<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MISR along-track dimension.\n";
				return FAILED;
			}

			// band
			if(H5DSattach_scale(misr_dataset,bandDset,1)<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MISR band dimension.\n";
				return FAILED;
			}

			//cameras
			if(H5DSattach_scale(misr_dataset,cameraDset,0)<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MISR cameras dimension.\n";
				return FAILED;
			}

			
			char* units = "Watts/m^2/micrometer/steradian";
			float valid_min = 0;
			float valid_max = 800.0;
			float _FillValue = -999.0;
			if(af_write_cf_attributes(misr_dataset, units, _FillValue, valid_min,valid_max,0) < 0) {
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: af_write_cf_attributes" << std::endl;
				return FAILED;
			}
			
			// Add CF long_name 
			const char* long_name = "long_name";
 			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),long_name,"MISR Level 1B2 Radiance")<0) {
					H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate long_name attribute for MISR Radiance" << std::endl;
				return FAILED;
			}
//#if 0 				
			//std::string misr_resolution = inputArgs.GetMISR_Resolution();
			std::string camera_angle_values;
			std::vector<std::string> camera_angle_vec = inputArgs.GetMISR_CameraAngles();
			for(int i = 0; i<camera_angle_vec.size();i++) 
				camera_angle_values = camera_angle_values + camera_angle_vec[i] +',';
			camera_angle_values = camera_angle_values.erase(camera_angle_values.size()-1,1);
			

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"camera_names",camera_angle_values.c_str())<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate camera_angles" << std::endl;
				return FAILED;
			}

			std::string radiance_values;
			std::vector<std::string> radiance_vec = inputArgs.GetMISR_Radiance();
			for(int i = 0; i<radiance_vec.size();i++) 
				radiance_values = radiance_values + radiance_vec[i] +',';
			radiance_values = radiance_values.erase(radiance_values.size()-1,1);
			

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"band_names",radiance_values.c_str())<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate radiances" << std::endl;
				return FAILED;
			}

			// Add resolution as a number.
			float misr_resolution_value = inputArgs.GetInstrumentResolutionValue(MISR_STR);
			if(false == af_AddSrcSpatialResolutionAttrs(outputFile,dsetPath,misr_resolution_value,true)) {
				H5Dclose(misr_dataset);
				return FAILED;
			}

			if(inputArgs.GetTargetInstrument() != "USER_DEFINE") {
				float target_resolution_value = inputArgs.GetInstrumentResolutionValue(inputArgs.GetTargetInstrument());
				if(false == af_AddSrcSpatialResolutionAttrs(outputFile,dsetPath,target_resolution_value,false)) {
					H5Dclose(misr_dataset);
					return FAILED;
				}

				if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"spatial_resolution_resampled_units","meter")<0) {
					H5Dclose(misr_dataset);
					std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
					return FAILED;
				}
			}
			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"spatial_resolution_original_units","meter")<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
				return FAILED;
			}

			if("USER_DEFINE" == inputArgs.GetTargetInstrument()) {
				
				std::string comments_ud = "Check the group attributes under group /Geolocation to find the user-defined ESPG code and resolution information.";
				if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"spatial_resolution_resampled_description",comments_ud.c_str())<0) {
					H5Dclose(misr_dataset);
					std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
					return FAILED;
				}	
			}


			// Add resample method
			std::string resample_method_value = "Summary Interpolation";
			if(inputArgs.GetResampleMethod()=="nnInterpolate")
				resample_method_value = "Nearest Neighbor Interpolation";

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"resample_method",resample_method_value.c_str())<0) {
				H5Dclose(misr_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate resample_method attr." << std::endl;
				return FAILED;

			}
//#endif
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

	status = H5Sselect_hyperslab(fileSpaceH5, H5S_SELECT_SET, startFile, NULL, countFile, NULL);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Sselect_hyperslab for Misr target .\n";
		ret = FAILED;
		goto done;
	}

	status = H5Dwrite(misr_dataset, dataTypeH5, memspace, fileSpaceH5, H5P_DEFAULT, processedData);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dwrite for Misr target .\n";
		ret = FAILED;
		goto done;
	}
        

done:
	H5Dclose(misr_dataset);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
	return ret;
}


/*=====================================================================
 * DESCRIPTION:
 *   Write resampled radiance output data of a single orbit for all the
 *   specified cameras and radiances for ASTER as the source instrument.
 *
 * PARAMETER:
 *  - inputArgs : a class object contains all the user input parameter info
 *  - outputFile : HDF5 id for output file
 *  - targetNNsrcID : got from nearestNeighborBlockIndex()
 *  - trgCellNum : number of target instrument data cells
 *  - srcFile : HDF5 id for input file
 *  - srcCellNum : number of source instrument data cells
 *  - inputMultiVarsMap : To obtain multiple values from a given user input
 *    directive which allows to have multiple values.
 *
 * RETURN:
 *  - Success: SUCCEED  (defined in AF_common.h)
 *  - Fail : FAILED  (defined in AF_common.h)
 */
int af_GenerateOutputCumulative_MisrAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNum, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap,hid_t ctrackDset, hid_t atrackDset)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;

	// strVec_t multiVarNames = inputArgs.GetMultiVariableNames(MISR_STR); // misr_MultiVars;
	std::string misrResolution = inputArgs.GetMISR_Resolution();

	// two multi-value variables are expected as this point
	strVec_t cameras = inputMultiVarsMap[MISR_CAMERA_ANGLE];
	strVec_t radiances = inputMultiVarsMap[MISR_RADIANCE];

	// Create band and camera dimensions
	hid_t bandDset = create_pure_dim_dataset(outputFile,(hsize_t)(radiances.size()),"Band_MISR");
	if(bandDset < 0) {
		printf("create_pure_dim_dataset for MISR band failed.\n");
		return FAILED;
	}
	hid_t cameraDset = create_pure_dim_dataset(outputFile,(hsize_t)(cameras.size()),"Camera");
	if(bandDset < 0) {
		printf("create_pure_dim_dataset for MISR Camera failed.\n");
		return FAILED;
	}

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
			StopElapseTimeAndShow("DBG_TIME> Read source MISR single band data	DONE.");
			#endif
	
			//-------------------------------------------------
			// handle resample method
			srcProcessedData = new double [trgCellNum];
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
				summaryInterpolate(misrSingleData, targetNNsrcID, srcCellNum, srcProcessedData, NULL, nsrcPixels, trgCellNum);
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
			ret = af_WriteSingleRadiance_MisrAsSrc<double,float>(inputArgs,outputFile, misrDatatype, misrDataspace,  srcProcessedData, trgCellNum /*processed size*/, srcOutputWidth, j /*cameraIdx*/, i /*radIdx*/,ctrackDset,atrackDset,cameraDset,bandDset);
			if (ret == FAILED) {
				std::cerr << __FUNCTION__ << "> Error: returned fail.\n";
			}
			#if DEBUG_ELAPSE_TIME
			StopElapseTimeAndShow("DBG_TIME> Write source MISR single band data  DONE.");
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

	H5Dclose(cameraDset);
	H5Dclose(bandDset);
	H5Tclose(misrDatatype);
	H5Sclose(misrDataspace);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif

	return ret;
}
