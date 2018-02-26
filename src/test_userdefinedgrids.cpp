#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "gdalio.h"
#include "reproject.h"
#include "io.h"

int main(int argc, char ** argv) {

/* Test input for wgs 84
	int outputEPSG = 4326;
	double xMin = -180;
	double xMax = 180;
	double yMin = -90;
	double yMax = 90;
	double cellSize = 1;
*/
/* Test input for wgs 84 / utm zone 23s */
	int outputEPSG = 32723;
	double xMin = 100000;
	double xMax = 900000;
	double yMin = 7000000;
	double yMax = 10000000;
	double cellSize = 2000;


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

	printf("%d output cells in total.\n", nPoints);
/*
	for(int i = 0; i < nPoints; i++) {
		printf("%lf,\t%lf\n", targetX[i], targetY[i]);
	}
*/

	char* file_path = "/projects/sciteam/jq0/TerraFusion/yizhao/TERRA_BF_L1B_O69626_20130119123228_F000_V001.h5";
	char* outfileName = "test32723.tif";
	
	hid_t src_file;
	if(0 > (src_file = af_open(file_path))) {
		printf("File not found\n");
		exit(1);
	}
	
	int nCellsrc;
	double* src_lat;
	double* src_long;

	src_lat = get_modis_lat(src_file, "_1KM", &nCellsrc);
	src_long = get_modis_long(src_file, "_1KM", &nCellsrc);
	

	double** p_src_lat = &src_lat;
	double** p_src_lon = &src_long;
	
	int * tarNNSouID;
	tarNNSouID = (int *)malloc(sizeof(int) * nPoints);

	printf("nearest neighbor\n");

	nearestNeighborBlockIndex(p_src_lat, p_src_lon, nCellsrc, targetY, targetX, tarNNSouID, NULL, nPoints, 1000);
	
	src_lat = *p_src_lat;
	src_long = *p_src_lon;
	
	free(src_lat);
	free(src_long);

	free(targetX);
	free(targetY);

	printf("getting source rad\n");
	double* src_rad;
	std::vector<std::string> bands = {"8"};
	src_rad = get_modis_rad(src_file, "_1KM", bands, bands.size(), &nCellsrc);
	herr_t ret = af_close(src_file);

	double* src_rad_out = (double *)malloc(sizeof(double) * nPoints);
	printf("interpolating\n");

	nnInterpolate(src_rad, src_rad_out, tarNNSouID, nPoints);
	
	printf("writing data fields\n");
	
	writeGeoTiff(outfileName, src_rad_out, outputEPSG, xMin, yMin, xMax, yMax, cellSize);
	

	free(src_rad);
	free(src_rad_out);
	free(tarNNSouID);

	return 0;
}
