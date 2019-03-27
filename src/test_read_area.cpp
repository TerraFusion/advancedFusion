/*


    AUTHOR:
        Yat Long Lo

    EMAIL:
        yllo2@illinois.edu


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
	hid_t output_file = H5Fcreate("misr_modis_test_repro.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	char* file_path = "/projects/TDataFus/kent/temp/40-orbit-file/Jun15.2/TERRA_BF_L1B_O69365_F000_V000.h5";
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
	int64_t size;
    std::vector<std::string> bands = {"8", "9", "12", "14L", "20"};
	double* modis_test = get_modis_rad(file, "_1KM", bands, bands.size(), &size);
	printf("test modis size: %d\n", size);
	printf("test data: %f\n", modis_test[0]);
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
