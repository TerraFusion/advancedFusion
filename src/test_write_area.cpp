/*


    AUTHOR:
        Yat Long Lo

    EMAIL:
        yllo2@illinois.edu

	PROGRAM DESCRIPTION:
		This is merely an area used for testing the correctness of certain functions in the IO module
*/
#if 1 // JK_WORK
#include <string>
#include <vector>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>
#include <sys/time.h>
#include "reproject.h"
#include "io.h"

int main(int argc, char ** argv) 
{
	hid_t output_file = H5Fcreate("test_write.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	char* file_path = "/projects/sciteam/jq0/TerraFusion/testFiles/TERRA_BF_L1B_O69400_20130104000439_F000_V000.h5";
	hid_t file;
	if(0 > (file = af_open(file_path))) {
		printf("File not found\n");
		exit(1);
	}
	
	/*char* m_250_list[2] = {"1", "2"};
		char* m_500_list[5] = {"3", "4", "5", "6", "7"};
	char* km_1_ref_list[15] = {"8", "9", "10", "11", "12", "13L", "13H", "14L", "14H", "15", "16", "17", "18", "19", "26"};
	char* kme_1_list[16] = {"20", "21", "22", "23", "24", "25", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36"};
	
	int nCellMODIS;
	int band_index;
	int file_size;
	int h;
	for(h = 0; h < 15; h++){
		char* d_name = get_modis_filename("_1KM", km_1_ref_list[h], &band_index);
		double * MODIS_rad = get_modis_rad_by_band(file, "_1KM", d_name, &band_index, &file_size);
		printf("MODIS_rad: %f\n", MODIS_rad[0]);
		printf("MODIS_rad: %f\n", MODIS_rad[2748620]);
	}*/
	
	/*int geo_size;
	double* modis_lat = get_modis_lat(file, "_1KM", &geo_size);
	printf("test size: %d\n", geo_size);
	return 0;*/
	int size;
    #if 1 // JK_WORK
    std::vector<std::string> bands = {"8", "9"};
	double* modis_test = get_modis_rad(file, "_1KM", bands, bands.size(), &size);
    #else
	char bands[2][50] = {"8", "9"};
	double* modis_test = get_modis_rad(file, "_1KM", bands, 2, &size);
	#endif
	printf("test modis size: %d\n", size);
	printf("test data: %f\n", modis_test[30248359]);
	hid_t group_id = H5Gcreate2(output_file, "/Data_Fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	//Write MODIS first
	hsize_t modis_dim[3];
	#if 1 // JK_WORK
	modis_dim[0] = bands.size();
	modis_dim[2] = 1354;
	modis_dim[1] = (size)/bands.size()/1354;
	#else
	modis_dim[0] = 2;
	modis_dim[2] = 1354;
	modis_dim[1] = (size)/2/1354;
	#endif
	hid_t modis_dataspace = H5Screate_simple(3, modis_dim, NULL);
	hid_t	modis_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    herr_t  modis_status = H5Tset_order(modis_datatype, H5T_ORDER_LE);  
    hid_t modis_dataset = H5Dcreate2(output_file, "/Data_Fields/modis_rad", modis_datatype, modis_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    modis_status = H5Dwrite(modis_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, modis_test);
    H5Sclose(modis_dataspace);
	H5Tclose(modis_datatype);
	H5Dclose(modis_dataset);
    if(modis_status < 0){
    	printf("MODIS write error\n");
    	return -1;
	}
	/*int rad_size;
	double* aster_rad = get_ast_rad(file, "TIR", "ImageData10", &rad_size);
	int lat_size;
	double* aster_lat = get_ast_lat(file, "TIR", "ImageData10", &lat_size);
	int long_size;
	double* aster_long = get_ast_long(file, "TIR", "ImageData10", &long_size);
	printf("test rad size: %d\n", rad_size);
	printf("test lat size: %d\n", lat_size);
	printf("test long size: %d\n", long_size);*/
	
	herr_t ret = af_close(file);


	return 0;
}
