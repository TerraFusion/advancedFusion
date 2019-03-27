/**
 * test_modis2aster.cpp
 * Authors: Yizhao Gao <yizhaotsccsj@gmail.com>
 * Date: {04/26/2018}
 * Description: test code using summary interpolation to resample from MODIS to ASTER, using nearest neighbor approach. 
 * Things to Do: now the raster pixels are output as a long 1D array. However, it makes more sense to organize aster images again into granules based on the input.
 */
#include <string>
#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>
#include <sys/time.h>
#include "reproject.h"
#include "io.h"

int main(int argc, char ** argv)
{
	hid_t output_file = H5Fcreate("/projects/sciteam/jq0/TerraFusion/yizhao/modis_2_aster_Test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	char* file_path = "/projects/sciteam/jq0/TerraFusion/yizhao/TERRA_BF_L1B_O69400_20130104000439_F000_V000.h5";
	hid_t src_file;
	if(0 > (src_file = af_open(file_path))) {
		printf("File not found\n");
		exit(1);
	}
	
	int64_t nCellsrc;
	double* src_lat;
	double* src_long;
	src_lat = get_modis_lat(src_file, "_1KM", &nCellsrc);
	src_long = get_modis_long(src_file, "_1KM", &nCellsrc);
	
	int64_t nCelldest;
	double* dest_lat;
	double* dest_long;
	
	dest_lat = get_ast_lat(src_file, "TIR", "ImageData10", &nCelldest);
	dest_long = get_ast_long(src_file, "TIR", "ImageData10", &nCelldest);
	
	printf("writing dest geo\n");
	printf("Number of total ASTER pixels: %d\n", nCelldest);

	//The af_write_mm_geo function is not good here. It is used temporarily since it is the best we have at this stage. 
	// There is no way that ASTER pixels can be written as a 2D block since the size of each granule is different. 	
	int lat_status =  af_write_mm_geo(output_file, 0, dest_lat, nCelldest, nCelldest,-1,-1);
	int long_status = af_write_mm_geo(output_file, 1, dest_long, nCelldest, nCelldest,-1,-1);
	if(lat_status < 0 || long_status < 0){
		printf("Writing dest geolocation - error\n");
		return -1;
	}

	double* projected_Rad_Out;
	int64_t * tarNNSouID;
	double** p_src_lat = &src_lat;
	double** p_src_lon = &src_long;
	//double** p_dest_lat = &dest_lat;
	//double** p_dest_lon = &dest_long;
	
	tarNNSouID = (int64_t *)malloc(sizeof(int64_t) * nCelldest);
	
	printf("nearest neighbor\n");
	//nearestNeighbor(p_src_lat, p_src_lon, nCellsrc, dest_lat, dest_long, tarNNSouID, nCelldest, 1000);
	nearestNeighborBlockIndex(p_src_lat, p_src_lon, nCellsrc, dest_lat, dest_long, tarNNSouID, NULL, nCelldest, 1000);
	src_lat = *p_src_lat;
	src_long = *p_src_lon;	
	
	free(src_lat);
	free(src_long);
	free(dest_lat);
	free(dest_long);
	
	printf("getting source rad\n");
	double* src_rad;
	
	std::vector<std::string> bands = {"8"};
	src_rad = get_modis_rad(src_file, "_1KM", bands, bands.size(), &nCellsrc);
	
	int64_t nCelldest_rad;
	double* dest_rad;
	
	dest_rad = get_ast_rad(src_file, "TIR", "ImageData10", &nCelldest_rad);
	
	double* src_rad_out = (double *)malloc(sizeof(double) * nCelldest);
	//Interpolating
	printf("interpolating\n");
	
	//Also collects the SD
	nnInterpolate(src_rad, src_rad_out, tarNNSouID, nCelldest);

//	for(int i = 0; i < nCelldest; i++) {
//		if(nsrcPixels[i] > 0) {
//			printf("%d,\t%lf\n", nsrcPixels[i], src_rad_out[i]);
//		}
//	}

	
	//Writing
	printf("writing data fields\n");
	hid_t group_id = H5Gcreate2(output_file, "/Data_Fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//Write ASTER first (as a large 1D array)
	hsize_t ast_dim[1];
	ast_dim[0] = nCelldest;
	hid_t ast_dataspace = H5Screate_simple(1, ast_dim, NULL);
	hid_t ast_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t ast_status = H5Tset_order(ast_datatype, H5T_ORDER_LE);  
	hid_t ast_dataset = H5Dcreate2(output_file, "/Data_Fields/aster_rad", ast_datatype, ast_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	ast_status = H5Dwrite(ast_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dest_rad);
	H5Sclose(ast_dataspace);
	H5Tclose(ast_datatype);
	H5Dclose(ast_dataset);
	if(ast_status < 0){
		printf("ast write error\n");
		return -1;
	}

	//Write MODIS next (also as a large 1D array)
	hsize_t modis_dim[1];
	modis_dim[0] = nCelldest;
	hid_t modis_dataspace = H5Screate_simple(1, modis_dim, NULL);
	hid_t	modis_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t  modis_status = H5Tset_order(modis_datatype, H5T_ORDER_LE);  
	hid_t modis_dataset = H5Dcreate2(output_file, "/Data_Fields/modis_rad", modis_datatype, modis_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	modis_status = H5Dwrite(modis_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_rad_out);
	H5Sclose(modis_dataspace);
	H5Tclose(modis_datatype);
	H5Dclose(modis_dataset);
	if(modis_status < 0){
		printf("MODIS write error\n");
    		return -1;
	}


	free(src_rad);
	free(src_rad_out);
	free(tarNNSouID);

	printf("Writing done\n");
	//Closing file
	herr_t ret = af_close(src_file);
	ret = af_close(output_file);
	return 0;
}
