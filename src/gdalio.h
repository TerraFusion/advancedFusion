/**
 * gdalio.h
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {01/14/2018}
 */

#ifndef GDALIOH
#define GDALIOH

void gdalIORegister();
int getCellCenterLatLon(int outputEPSG, double xMin, double yMin, double xMax, double yMax, double cellSize, double ** px, double ** py);
void writeGeoTiff(char * fileName, double * grid, double xMin, double yMin, double xMax, double yMax, double cellSize);

#endif
