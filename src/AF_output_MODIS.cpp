/*********************************************************************
 * DESCRIPTION:
 *	Generate radiance data output to HDF5 for MODIS
 *
 *
 * DEVELOPERS:
 *	- Jonathan Kim (jkm@illinois.edu)
 *	- Hyo-Kyung (Joe) Lee (hyoklee@hdfgroup.org)
	  \author Hyo-Kyung (Joe) Lee (hyoklee@hdfgroup.org)
	  \date may 18, 2018
	  \note added cf attributes in af_WriteSingleRadiance_ModisAsSrc().
	  \date May 15, 2018
	  \note added valid_min attribute in af_WriteSingleRadiance_ModisAsTrg().
	  \date May 14, 2018
	  \note added coordinates and _FillValue attributes in af_WriteSingleRadiance_ModisAsTrg().
 */

#include "AF_output_MODIS.h"

#include "AF_common.h"
#include "AF_output_util.h"
#include "io.h"
#include "reproject.h"
#include "misrutil.h"
#include <algorithm>
#include "gdalio.h"

//const char* ref_band_list[22] = {"1","2","3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13L", "13H", "14L", "14H", "15", "16", "17", "18", "19", "26"};
strVec_t ref_band_list = {"1","2","3", "4", "5", "6", "7","8", "9", "10", "11", "12", "13L", "13H", "14L", "14H", "15", "16", "17", "18", "19", "26"};
const std::string All_bands_1km = "1,2,3,4,5,6,7,8,9,10,11,12,13L,13H,14L,14H,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36";
const std::string All_bands_500m = "1,2,3,4,5,6,7";
const std::string All_bands_250m = "1,2";
/*#############################################################################
 *
 * MODIS as Target instrument, functions to generate radiance data
 *
 *############################################################################*/

/*=====================================================================
 * DESCRIPTION:
 *	 Write radiance output data of a single orbit of the given band for
 *	 MODIS as the target instrument.
 *	 Only called by af_GenerateOutputCumulative_ModisAsTrg().
 *
 * PARAMETER:
 *	- outputFile : HDF5 id for output file
 *	- dataTypeH5 : HDF5 id for output datatype
 *	- fileSpaceH5 : HDF5 id file sapce
 *	- modisData : MODIS data pointer (not processed)
 *	- modisDataSize : number of cells (pixels) in MODIS data
 *	- outputWidth : cross-track (width) size for output image
 *	- bandIdx : MODIS band index.
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
static int af_WriteSingleRadiance_ModisAsTrg(AF_InputParmeterFile &inputArgs,hid_t outputFile, hid_t dataTypeH5, hid_t fileSpaceH5, T_IN* modisData, int modisDataSize, int outputWidth, int bandIdx,bool has_refsb,const strVec_t bands,hid_t ctrackDset, hid_t atrackDset, hid_t bandDset)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	std::cout << "MODIS bands: ";
	for(int i=0; i < bands.size(); i++) {
		std::cout << bands[i] << " ";
	}
	std::cout << std::endl;
	#endif

	int ret = SUCCEED;

	herr_t status;
	hid_t modis_dataset;
	std::string dsetPath = TRG_DATA_GROUP + "/" + MODIS_RADIANCE_DSET;


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
			if(create_chunk_comp_plist(plist_id,3,(size_t)modisDataSize,(size_t)outputWidth)<0) {
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: create_chunk_comp_plist failed.\n";
				H5Pclose(plist_id);
				return FAILED;
			}
			modis_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), dataTypeOutH5, fileSpaceH5,H5P_DEFAULT, plist_id, H5P_DEFAULT);
			H5Pclose(plist_id);
		}
		else
			modis_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), dataTypeOutH5, fileSpaceH5,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(modis_dataset < 0) {
			std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dcreate2 target data in output file.\n";
			return FAILED;
		}
		else {
			// make compatible with CF convention (NetCDF)/
			// Attach dimension scales.
			// cross track
			if(H5DSattach_scale(modis_dataset,ctrackDset,2)<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MODIS cross-track dimension.\n";
				return FAILED;
			}

			// along track
			if(H5DSattach_scale(modis_dataset,atrackDset,1)<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MODIS along-track dimension.\n";
				return FAILED;
			}

			// band
			if(H5DSattach_scale(modis_dataset,bandDset,0)<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MODIS band dimension.\n";
				return FAILED;
			}

			char* units = "Watts/m^2/micrometer/steradian";
			float valid_min = 0;
			float valid_max = 100.0;
 			float _FillValue = -999.0;
			if(true == has_refsb)
				valid_max = 900.0;
			
			if(af_write_cf_attributes(modis_dataset, units, _FillValue,valid_min,valid_max,0) < 0) {
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: af_write_cf_attributes" << std::endl;
			}

			// Add CF long_name 
			const char* long_name = "long_name";
 			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),long_name,"MODIS Level 1B Radiance")<0) {
					H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate long_name attribute for ASTER Radiance" << std::endl;
				return FAILED;
			}
 				
			std::string modis_resolution = inputArgs.GetMODIS_Resolution();
			std::string band_name_values;
			std::vector<std::string> band_name_vec = inputArgs.GetMODIS_Bands();
			//Check if 'ALL' is in bands
			bool isAllBands = inputArgs.IsMODIS_AllBands();
			if(true == isAllBands) {
				if("_1KM" == modis_resolution) 
					band_name_values = All_bands_1km;
				else if("_500m" == modis_resolution) 
					band_name_values = All_bands_500m;
				else
					band_name_values = All_bands_250m;
			}
			else {
				for(int i = 0; i<band_name_vec.size();i++) 
					band_name_values = band_name_values + band_name_vec[i] +',';
				band_name_values = band_name_values.erase(band_name_values.size()-1,1);
			}

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"band_names",band_name_values.c_str())<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate band_names" << std::endl;
				return FAILED;
			}


			// Add resolution as a number.
			float modis_resolution_value = inputArgs.GetInstrumentResolutionValue(MODIS_STR);

			if(H5LTset_attribute_float(outputFile,dsetPath.c_str(),"spatial_resolution",&modis_resolution_value, 1 ) <0) {
            	std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate resolution attribute." << std::endl;
            	return FAILED;
        	}

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"spatial_resolution_units","meter")<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
				return FAILED;
			}

			std::vector<int> modis_radiance_type_list = inputArgs.GetMODIS_Radiance_TypeList();
			if(H5LTset_attribute_int(outputFile,dsetPath.c_str(),"MODIS_radiance_type",&modis_radiance_type_list[0],modis_radiance_type_list.size())<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate MODIS radiance type" << std::endl;
				return FAILED;
			}
			
			std::string comment_for_radiance_type ="Description for Attribute MODIS_radiance_type. \n" ;
			comment_for_radiance_type += " There are seven radiance types depending on different MODIS Bands.";
			comment_for_radiance_type +=" The corresponding band number can be found in attribute band_names.";
			comment_for_radiance_type +=" They are represented by numerical number from 1 to 7. \n";
			comment_for_radiance_type +=" 1: 1KM_Emissive, 2: 1KM_RefSB, 3: 500_Aggr1km_RefSB, 4: 250_Aggr1km_RefSB \n";
			comment_for_radiance_type +=" 5: 500_RefSB, 6: 250_Aggr500_RefSB, 7: 250_RefSB";

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"comments_for_radiance_type",comment_for_radiance_type.c_str())<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
				return FAILED;
			}

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

	status = H5Sselect_hyperslab(fileSpaceH5, H5S_SELECT_SET, star3dFile, NULL, count3dFile, NULL);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Sselect_hyperslab for Modis target .\n";
		ret = FAILED;
		goto done;
	}

	status = H5Dwrite(modis_dataset, dataTypeH5, memspace, fileSpaceH5, H5P_DEFAULT, modisData);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dwrite for Modis target .\n";
		ret = FAILED;
		goto done;
	}
done:
	H5Dclose(modis_dataset);
	
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif
	return ret;
}


/*=========================================================================
 * DESCRIPTION:
 *	 Write radiance output data of a single orbit for all the
 *	 specified bands for MODIS as the target instrument.
 *
 * PARAMETER:
 *	- inputArgs : a class object contains all the user input parameter info
 *	- outputFile : HDF5 id for output file
 *	- srcFile : HDF5 id for input file
 *	- trgCellNum : number of target instrument data cells
 *	- inputMultiVarsMap : To obtain multiple values from a given user input
 *	  directive which allows to have multiple values.
 *
 * RETURN:
 *	- Success: SUCCEED	(defined in AF_common.h)
 *	- Fail : FAILED  (defined in AF_common.h)
 */
int af_GenerateOutputCumulative_ModisAsTrg(AF_InputParmeterFile &inputArgs, hid_t outputFile,hid_t srcFile, int trgCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap,hid_t ctrackDset, hid_t atrackDset)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	std::string modisResolution = inputArgs.GetMODIS_Resolution();
	// get multi-value variable for modis
	strVec_t bands = inputMultiVarsMap[MODIS_BANDS];
	// Create band dimension 
	hid_t bandDset = create_pure_dim_dataset(outputFile,(hsize_t)(bands.size()),"Band_MODIS");
	if(bandDset < 0) {
		printf("create_pure_dim_dataset for MODIS band failed.\n");
		return FAILED;
	}
	
        bool has_refsb = false;

        for (int temp_i= 0; temp_i < bands.size(); temp_i++) {
        	has_refsb = (std::find(std::begin(ref_band_list),std::end(ref_band_list),bands[temp_i])!=std::end(ref_band_list));
		if(true == has_refsb)
			break;
	}
        
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
	const int rankSpace=3;	// [bands][y][x]
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
		StopElapseTimeAndShow("DBG_TIME> Read target MODIS single band data  DONE.");
		#endif

		//---------------------------------
		// write trg radiance to AF file
		#if DEBUG_ELAPSE_TIME
		StartElapseTime();
		#endif
		af_WriteSingleRadiance_ModisAsTrg<double,float>(inputArgs,outputFile, modisDatatype, modisDataspace,modisSingleData, numCells, targetOutputWidth, i,has_refsb,bands,ctrackDset,atrackDset,bandDset);
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG_TIME> Write target MODIS single band data  DONE.");
		#endif
		//
		// TODO: buffer is not reused. make memory allocation out of get_modis_rad() to improve performance
		free(modisSingleData);
	}

	H5Dclose(bandDset);
	H5Tclose(modisDatatype);
	H5Sclose(modisDataspace);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif

	return SUCCEED;
}



/*#############################################################################
 *
 * MODIS as Source instrument, functions to generate radiance data
 *
 *############################################################################*/

/*=====================================================================
 * DESCRIPTION:
 *	 Write resampled radiance output data of a single orbit of the given
 *	 band for MODIS as the source instrument.
 *	 Only called by af_GenerateOutputCumulative_ModisAsSrc().
 *
 * PARAMETER:
 *	- outputFile : HDF5 id for output file
 *	- dataTypeH5 : HDF5 id for output datatype
 *	- fileSpaceH5 : HDF5 id file sapce
 *	- processedData : resmapled data pointer
 *	- trgCellNum : number of cells (pixels) in target instrument data
 *	- outputWidth : cross-track (width) size for output image
 *	- bandIdx : MODIS band index.
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
static int af_WriteSingleRadiance_ModisAsSrc(AF_InputParmeterFile &inputArgs,hid_t outputFile, hid_t dataTypeH5, hid_t fileSpaceH5, T_IN* processedData, int trgCellNum, int outputWidth, int bandIdx,bool has_refsb,const strVec_t bands, hid_t ctrackDset,hid_t atrackDset,hid_t bandDset)
{
#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	std::cout << "MODIS bands: ";
	for(int i=0; i < bands.size(); i++) {
		std::cout << bands[i] << " ";
	}
	std::cout << std::endl;

#endif

	int ret = SUCCEED;
	herr_t status;
	hid_t modis_dataset;
	std::string dsetPath = SRC_DATA_GROUP + "/" + MODIS_RADIANCE_DSET;


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
			modis_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), dataTypeOutH5, fileSpaceH5,H5P_DEFAULT, plist_id, H5P_DEFAULT);
			H5Pclose(plist_id);
		}
		else
			modis_dataset = H5Dcreate2(outputFile, dsetPath.c_str(), dataTypeOutH5, fileSpaceH5,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		
		if(modis_dataset < 0) {
			std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dcreate2 target data in output file.\n";
			return FAILED;
		}
		else {
			// make compatible with CF convention (NetCDF)
			// Attach dimension scales.
			// cross track
			if(H5DSattach_scale(modis_dataset,ctrackDset,2)<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MODIS cross-track dimension.\n";
				return FAILED;
			}

			// along track
			if(H5DSattach_scale(modis_dataset,atrackDset,1)<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MODIS along-track dimension.\n";
				return FAILED;
			}

			// band
			if(H5DSattach_scale(modis_dataset,bandDset,0)<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5DSattach_scale failed for MODIS band dimension.\n";
				return FAILED;
			}

			char* units = "Watts/m^2/micrometer/steradian";
			float valid_min = 0;
			float valid_max = 100.0;
			float _FillValue = -999.0;
			if(true == has_refsb)
				valid_max = 900.0;
			
			if(af_write_cf_attributes(modis_dataset, units, _FillValue,valid_min,valid_max,0) < 0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ <<	"> Error: af_write_cf_attributes" << std::endl;
				return FAILED;
			}

			// Add CF long_name 
			const char* long_name = "long_name";
 			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),long_name,"MODIS Level 1B Radiance")<0) {
					H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate long_name attribute for ASTER Radiance" << std::endl;
				return FAILED;
			}
 				
			std::string modis_resolution = inputArgs.GetMODIS_Resolution();
			std::string band_name_values;
			std::vector<std::string> band_name_vec = inputArgs.GetMODIS_Bands();
			//Check if 'ALL' is in bands
			bool isAllBands = inputArgs.IsMODIS_AllBands();
			if(true == isAllBands) {
				if("_1KM" == modis_resolution) 
					band_name_values = All_bands_1km;
				else if("_500m" == modis_resolution) 
					band_name_values = All_bands_500m;
				else
					band_name_values = All_bands_250m;
			}
			else {
				for(int i = 0; i<band_name_vec.size();i++) 
					band_name_values = band_name_values + band_name_vec[i] +',';
				band_name_values = band_name_values.erase(band_name_values.size()-1,1);
			}

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"band_names",band_name_values.c_str())<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate band_names" << std::endl;
				return FAILED;
			}


			// Add resolution as a number.
			float modis_resolution_value = inputArgs.GetInstrumentResolutionValue(MODIS_STR);
			if(false == af_AddSrcSpatialResolutionAttrs(outputFile,dsetPath,modis_resolution_value,true)) {
				H5Dclose(modis_dataset);
				return FAILED;
			}

			if(inputArgs.GetTargetInstrument() != "USER_DEFINE") {
				float target_resolution_value = inputArgs.GetInstrumentResolutionValue(inputArgs.GetTargetInstrument());
				if(false == af_AddSrcSpatialResolutionAttrs(outputFile,dsetPath,target_resolution_value,false)) {
					H5Dclose(modis_dataset);
					return FAILED;
				}
				if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"spatial_resolution_resampled_units","meter")<0) {
                	H5Dclose(modis_dataset);
                	std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
                	return FAILED;
            	}

			}

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"spatial_resolution_original_units","meter")<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
				return FAILED;
			}

			if("USER_DEFINE" == inputArgs.GetTargetInstrument()) {
				
				std::string comments_ud = "Check the group attributes under group /Geolocation to find the user-defined ESPG code and resolution information.";
				if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"spatial_resolution_resampled_description",comments_ud.c_str())<0) {
					H5Dclose(modis_dataset);
					std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
					return FAILED;
				}	
			}

			std::vector<int> modis_radiance_type_list = inputArgs.GetMODIS_Radiance_TypeList();
			if(H5LTset_attribute_int(outputFile,dsetPath.c_str(),"MODIS_radiance_type",&modis_radiance_type_list[0],modis_radiance_type_list.size())<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate MODIS radiance type" << std::endl;
				return FAILED;
			}
		
			std::string comment_for_radiance_type ="Description for Attribute MODIS_radiance_type. \n";
			comment_for_radiance_type += " There are seven radiance types depending on different MODIS Bands.";
			comment_for_radiance_type +=" The corresponding band number can be found in attribute band_names.";
			comment_for_radiance_type +=" They are represented by numerical number from 1 to 7. \n";
			comment_for_radiance_type +=" 1: 1KM_Emissive, 2: 1KM_RefSB, 3: 500_Aggr1km_RefSB, 4: 250_Aggr1km_RefSB \n";
			comment_for_radiance_type +=" 5: 500_RefSB, 6: 250_Aggr500_RefSB, 7: 250_RefSB";

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"comments_for_radiance_type",comment_for_radiance_type.c_str())<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate spatial_resolution_units" << std::endl;
				return FAILED;
			}


			// Add resample method
			std::string resample_method_value = "Summary Interpolation";
			if(inputArgs.GetResampleMethod()=="nnInterpolate")
				resample_method_value = "Nearest Neighbor Interpolation";

			if(H5LTset_attribute_string(outputFile,dsetPath.c_str(),"resample_method",resample_method_value.c_str())<0) {
				H5Dclose(modis_dataset);
				std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: cannot generate resample_method attr." << std::endl;
				return FAILED;

			}
			
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

	status = H5Sselect_hyperslab(fileSpaceH5, H5S_SELECT_SET, startFile, NULL, countFile, NULL);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Sselect_hyperslab for Modis target .\n";
		ret = FAILED;
		goto done;
	}

	status = H5Dwrite(modis_dataset, dataTypeH5, memspace, fileSpaceH5, H5P_DEFAULT, processedData);
	if(status < 0) {
		std::cerr << __FUNCTION__ << ":" << __LINE__ <<  "> Error: H5Dwrite for Modis target .\n";
		ret = FAILED;
		goto done;
	}

 done:
	H5Dclose(modis_dataset);
//#if 0
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
		writeGeoTiff((char*)op_geotiff_fname.c_str(),processedData, userOutputEPSG, userXmin, userYmin, userXmax, userYmax, userResolution);


	}
//#endif


#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
#endif
	return ret;
}


/*=====================================================================
 * DESCRIPTION:
 *	 Write resampled radiance output data of a single orbit for all the
 *	 specified bands for MODIS as the source instrument.
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
int af_GenerateOutputCumulative_ModisAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNumNoShift, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap,hid_t ctrackDset, hid_t atrackDset)
{
	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> BEGIN \n";
	#endif

	int ret = SUCCEED;

	// strVec_t multiVarNames = inputArgs.GetMultiVariableNames(MODIS_STR); // modis_MultiVars;
	std::string modisResolution = inputArgs.GetMODIS_Resolution();

	// two multi-value variables are expected as this point
	strVec_t bands = inputMultiVarsMap[MODIS_BANDS];
	// Create band dimension 
	hid_t bandDset = create_pure_dim_dataset(outputFile,(hsize_t)(bands.size()),"Band_MODIS");
	if(bandDset < 0) {
		printf("create_pure_dim_dataset for MODIS band failed.\n");
		return FAILED;
	}

        bool has_refsb = false;

        for (int temp_i= 0; temp_i < bands.size(); temp_i++) {
        	has_refsb = (std::find(std::begin(ref_band_list),std::end(ref_band_list),bands[temp_i])!=std::end(ref_band_list));
		if(true == has_refsb)
			break;
	}
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
	const int rankSpace=3;	// [bands][y][x]
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
	if(inputArgs.GetMISR_Shift() == "ON" && inputArgs.GetTargetInstrument() == MISR_STR) {
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
		StopElapseTimeAndShow("DBG_TIME> Read source MODIS single band data	DONE.");
		#endif

		//-------------------------------------------------
		// handle resample method
		srcProcessedData = new double [trgCellNumNoShift];
		// Note: resample should be done with trgCellNumNoShift
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
			summaryInterpolate(modisSingleData, targetNNsrcID, srcCellNum, srcProcessedData, NULL, nsrcPixels, trgCellNumNoShift);
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
		if(inputArgs.GetMISR_Shift() == "ON" && inputArgs.GetTargetInstrument() == MISR_STR) {
			std::cout << "\nSource MODIS radiance MISR-base shifting...\n";
			#if DEBUG_ELAPSE_TIME
			StartElapseTime();
			#endif
			srcProcessedDataShifted = new double [widthShifted * heightShifted];
			MISRBlockOffset<double>(srcProcessedData, srcProcessedDataShifted, (inputArgs.GetMISR_Resolution() == "L") ? 0 : 1);
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
		ret = af_WriteSingleRadiance_ModisAsSrc<double, float>(inputArgs,outputFile, modisDatatype, modisDataspace,  srcProcessedDataPtr, numCells /*processed size*/, srcOutputWidth, i /*bandIdx*/,has_refsb/*radiance has refSB*/,bands,ctrackDset,atrackDset,bandDset);
		if (ret == FAILED) {
			std::cerr << __FUNCTION__ << "> Error: returned fail.\n";
		}
		#if DEBUG_ELAPSE_TIME
		StopElapseTimeAndShow("DBG_TIME> Write source MODIS single band data  DONE.");
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

	H5Dclose(bandDset);
	H5Tclose(modisDatatype);
	H5Sclose(modisDataspace);

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> END \n";
	#endif

	return ret;
}

