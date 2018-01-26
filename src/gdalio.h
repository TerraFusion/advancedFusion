/**
 * gdalio.h
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {01/25/2018}
 */

#ifndef GDALIOH
#define GDALIOH

/**
 * NAME:	gdalIORegister
 * DESCRIPTION:	Register drivers
 * PARAMETERS: void
 */
void gdalIORegister();


/**
 * NAME:	getCellCenterLatLon
 * DESCRIPTION:	Get the latitude and longtitude of pixel centers given a grid
 * PARAMETERS:
 * 	int outputEPSG:		EPSG code of output spatial reference system 
 *	double xMin:		west boundary of output area
 *	double yMin:		south boundary of output area
 *	double xMax:		east boundary of output area
 *	double yMax: 		north boundary of output area
 * 	double cellSize:	output raste cell size
 *	double ** px:		longitude of output pixel centers 
 *	double ** py:		latitude of ouput pixel centers
 * Output:
 *	double ** px:		longitude of output pixel centers, memory will be allocated in this function
 *	double ** py:		latitude of ouput pixel centers, memory will be allocated in this function
 */
int getCellCenterLatLon(int outputEPSG, double xMin, double yMin, double xMax, double yMax, double cellSize, double ** px, double ** py);


/**
 * NAME:	writeGeoTiff
 * DESCRIPTION:	Write the output grid as a GeoTiff
 * PARAMETERS:
 *	char * fileName:	output GeoTiff file name
 * 	double * grid:		the grid of the output radianc values
 * 	int outputEPSG:		EPSG code of output spatial reference system 
 *	double xMin:		west boundary of output area
 *	double yMin:		south boundary of output area
 *	double xMax:		east boundary of output area
 *	double yMax: 		north boundary of output area
 * 	double cellSize:	output raste cell size
 */
void writeGeoTiff(char * fileName, double * grid, int outputEPSG, double xMin, double yMin, double xMax, double yMax, double cellSize);

#endif
