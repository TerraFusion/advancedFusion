#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "gdalio.h"
#include "reproject.h"
#include "io.h"
#include <math.h>

int main(int argc, char ** argv) {

/* Test input for wgs 84 */
	int outputEPSG = 4326;
	double xMin = -180;
	double xMax = 180;
	double yMin = -90;
	double yMax = 90;
	double cellSize = 0.1;

/* Test input for wgs 84 / utm zone 23s 
	int outputEPSG = 32723;
	double xMin = 100000;
	double xMax = 900000;
	double yMin = 7000000;
	double yMax = 10000000;
	double cellSize = 2000;
*/

/* Test input for wgs 84 / utm zone 24s
	int outputEPSG = 32724;
	double xMin = 100000;
	double xMax = 900000;
	double yMin = 7000000;
	double yMax = 10000000;
	double cellSize = 2000;
*/

	double * targetX;
	double * targetY;

	
	gdalIORegister();

	int nPoints = getCellCenterLatLon(outputEPSG, xMin, yMin, xMax, yMax, cellSize, &targetX, &targetY);
	int crossTrack = ceil((xMax - xMin) / cellSize);

	printf("%d output cells in total.\n", nPoints);
/*
	for(int i = 0; i < nPoints; i++) {
		printf("%lf,\t%lf\n", targetX[i], targetY[i]);
	}
*/
    hid_t output_file = H5Fcreate("./test_userdefinedgrids_modis2user.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    printf("writing dest geo\n");
	// calculate width from user define input value
    int lat_status =  af_write_mm_geo(output_file, 0, targetY, nPoints, crossTrack);
    int long_status = af_write_mm_geo(output_file, 1, targetX, nPoints, crossTrack);
    if(lat_status < 0 || long_status < 0){
        printf("Writing dest geolocation - error\n");
        return -1;
    }
	

	char* file_path = "/projects/sciteam/jq0/TerraFusion/yizhao/TERRA_BF_L1B_O69400_20130104000439_F000_V000.h5";
	char* outfileName = "test4326_439_test.tif";
	
	hid_t src_file;
	if(0 > (src_file = af_open(file_path))) {
		printf("File not found\n");
		exit(1);
	}
	
	int nCellsrc;
	double* src_lat;
	double* src_long;

	src_lat = get_modis_lat(src_file, "_1KM", &nCellsrc);
	printf("lat pixel size: %d\n", nCellsrc);
	src_long = get_modis_long(src_file, "_1KM", &nCellsrc);
	printf("long pixel size: %d\n", nCellsrc);


	double** p_src_lat = &src_lat;
	double** p_src_lon = &src_long;
	
	int * tarNNSouID;
	tarNNSouID = (int *)malloc(sizeof(int) * nPoints);

	printf("nearest neighbor\n");

	nearestNeighborBlockIndex(p_src_lat, p_src_lon, nCellsrc, targetY, targetX, tarNNSouID, NULL, nPoints, 5000);
	
	src_lat = *p_src_lat;
	src_long = *p_src_lon;
	
	free(src_lat);
	free(src_long);

	free(targetX);
	free(targetY);

	printf("getting source rad\n");
	double* src_rad;
	std::vector<std::string> bands = {"25"};
	src_rad = get_modis_rad(src_file, "_1KM", bands, bands.size(), &nCellsrc);
	printf("rad pixel size: %d\n", nCellsrc);

	double* src_rad_out = (double *)malloc(sizeof(double) * nPoints);
	printf("interpolating\n");

	nnInterpolate(src_rad, src_rad_out, tarNNSouID, nPoints);
	
	printf("writing data fields\n");
	
	writeGeoTiff(outfileName, src_rad_out, outputEPSG, xMin, yMin, xMax, yMax, cellSize);

    //Writing Modis source
    printf("writing source data fields\n");
    hid_t group_id = H5Gcreate2(output_file, "/Data_Fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t modis_dim[3];
    modis_dim[0] = bands.size();
    modis_dim[1] = (nPoints)/bands.size()/crossTrack;
    modis_dim[2] = crossTrack;
    hid_t modis_dataspace = H5Screate_simple(3, modis_dim, NULL);
    hid_t   modis_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    herr_t  modis_status = H5Tset_order(modis_datatype, H5T_ORDER_LE);
    hid_t modis_dataset = H5Dcreate2(output_file, "/Data_Fields/modis_rad", modis_datatype, modis_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    modis_status = H5Dwrite(modis_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_rad_out);
    H5Sclose(modis_dataspace);
    H5Tclose(modis_datatype);
    H5Dclose(modis_dataset);
	H5Gclose(group_id);
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
