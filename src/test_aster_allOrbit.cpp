/*


    AUTHOR:
        Yat Long Lo; Yizhao Gao

    EMAIL:
        yllo2@illinois.edu; ygao29@illinois.edu

	PROGRAM DESCRIPTION:
		This is a program that makes use of the IO module and reprojection module to produce reprojected data of ASTER onto MODIS
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
	hid_t output_file = H5Fcreate("/projects/sciteam/jq0/TerraFusion/yizhao/aster_on_modis_3N_Test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	char* file_path = "/projects/sciteam/jq0/TerraFusion/yizhao/TERRA_BF_L1B_O69400_20130104000439_F000_V000.h5";
	hid_t src_file;
	if(0 > (src_file = af_open(file_path))) {
		printf("File not found\n");
		exit(1);
	}
	
	int nCellsrc;
	double* src_lat;
	double* src_long;
	src_lat = get_ast_lat(src_file, "VNIR", "ImageData3N", &nCellsrc);
	src_long = get_ast_long(src_file, "VNIR", "ImageData3N", &nCellsrc);
	
	int nCelldest;
	double* dest_lat;
	double* dest_long;
	
	dest_lat = get_modis_lat(src_file, "_1KM", &nCelldest);
	dest_long = get_modis_long(src_file, "_1KM", &nCelldest);
	
	printf("writing dest geo\n");
	// modis 1KM output width is '1354'
	int lat_status =  af_write_mm_geo(output_file, 0, dest_lat, nCelldest, 1354);
	int long_status = af_write_mm_geo(output_file, 1, dest_long, nCelldest, 1354);
	if(lat_status < 0 || long_status < 0){
		printf("Writing dest geolocation - error\n");
		return -1;
	}

	double* projected_Rad_Out;
	int * tarNNSouID;
	//double** p_src_lat = &src_lat;
	//double** p_src_lon = &src_long;
	double** p_dest_lat = &dest_lat;
	double** p_dest_lon = &dest_long;
	
	tarNNSouID = (int *)malloc(sizeof(int) * nCellsrc);
	
	printf("nearest neighbor\n");
	//nearestNeighbor(p_src_lat, p_src_lon, nCellsrc, dest_lat, dest_long, tarNNSouID, nCelldest, 1000);
	nearestNeighborBlockIndex(p_dest_lat, p_dest_lon, nCelldest, src_lat, src_long, tarNNSouID, NULL, nCellsrc, 1000);
	dest_lat = *p_dest_lat;
	dest_long = *p_dest_lon;
	
	
	free(src_lat);
	free(src_long);
	free(dest_lat);
	free(dest_long);
	
	printf("getting source rad\n");
	double* src_rad;
	
	src_rad = get_ast_rad(src_file, "VNIR", "ImageData3N", &nCellsrc);
	
	int nCelldest_rad;
	double* dest_rad;
	std::vector<std::string> bands = {"8"};
	dest_rad = get_modis_rad(src_file, "_1KM", bands, bands.size(), &nCelldest_rad);
	
	double* src_rad_out = (double *)malloc(sizeof(double) * nCelldest);
	int new_ast_size = nCelldest;
	//Interpolating
	int * nsrcPixels; //Number of contributing ASTER pixel to each new MODIS pixel
	double * sd; //The standard deviation of all contributing ASTER cell's value 
	nsrcPixels = (int *) malloc(sizeof(int) * nCelldest);
	sd = (double *) malloc(sizeof(int) * nCelldest);
	printf("interpolating\n");
	
	//Also collects the SD
	summaryInterpolate(src_rad, tarNNSouID, nCellsrc, src_rad_out, sd, nsrcPixels, nCelldest);

	printf("No nodata values: \n");
//	for(int i = 0; i < nCelldest; i++) {
//		if(nsrcPixels[i] > 0) {
//			printf("%d,\t%lf\n", nsrcPixels[i], src_rad_out[i]);
//		}
//	}

	
	//Writing
	printf("writing data fields\n");
	hid_t group_id = H5Gcreate2(output_file, "/Data_Fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Write MODIS first
	hsize_t modis_dim[3];
	modis_dim[0] = bands.size();
	modis_dim[1] = (nCelldest_rad)/bands.size()/1354;
	modis_dim[2] = 1354;
	hid_t modis_dataspace = H5Screate_simple(3, modis_dim, NULL);
	hid_t	modis_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t  modis_status = H5Tset_order(modis_datatype, H5T_ORDER_LE);  
	hid_t modis_dataset = H5Dcreate2(output_file, "/Data_Fields/modis_rad", modis_datatype, modis_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	modis_status = H5Dwrite(modis_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dest_rad);
	H5Sclose(modis_dataspace);
	H5Tclose(modis_datatype);
	H5Dclose(modis_dataset);
	if(modis_status < 0){
		printf("MODIS write error\n");
    		return -1;
	}
	//Write ASTER next
	hsize_t ast_dim[2];
	ast_dim[0] = (new_ast_size) / 1354;
	ast_dim[1] = 1354;
	//Write ASTER Radiance
	hid_t ast_dataspace = H5Screate_simple(2, ast_dim, NULL);
	hid_t ast_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	herr_t ast_status = H5Tset_order(ast_datatype, H5T_ORDER_LE);  
	hid_t ast_dataset = H5Dcreate2(output_file, "/Data_Fields/aster_average", ast_datatype, ast_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	ast_status = H5Dwrite(ast_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_rad_out);
	H5Sclose(ast_dataspace);
	H5Tclose(ast_datatype);
	H5Dclose(ast_dataset);
	if(ast_status < 0){
		printf("ast write error\n");
		return -1;
	}

	//Write ASTER SD
	ast_dataspace = H5Screate_simple(2, ast_dim, NULL);
	ast_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
	ast_status = H5Tset_order(ast_datatype, H5T_ORDER_LE);  
	ast_dataset = H5Dcreate2(output_file, "/Data_Fields/aster_sd", ast_datatype, ast_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	ast_status = H5Dwrite(ast_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_rad_out);
	H5Sclose(ast_dataspace);
	H5Tclose(ast_datatype);
	H5Dclose(ast_dataset);
	if(ast_status < 0){
		printf("ast write error\n");
		return -1;
	}

	//Write ASTER SD
	ast_dataspace = H5Screate_simple(2, ast_dim, NULL);
	ast_datatype = H5Tcopy(H5T_NATIVE_INT);
	ast_status = H5Tset_order(ast_datatype, H5T_ORDER_LE);  
	ast_dataset = H5Dcreate2(output_file, "/Data_Fields/aster_count", ast_datatype, ast_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	ast_status = H5Dwrite(ast_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_rad_out);
	H5Sclose(ast_dataspace);
	H5Tclose(ast_datatype);
	H5Dclose(ast_dataset);
	if(ast_status < 0){
		printf("ast write error\n");
		return -1;
	}
	
	free(nsrcPixels);
	free(sd);
	
	free(src_rad);
	free(src_rad_out);
	free(tarNNSouID);

	printf("Writing done\n");
	//Closing file
	herr_t ret = af_close(src_file);
	ret = af_close(output_file);
	return 0;
}
