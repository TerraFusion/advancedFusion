/**
 * gdalio.h
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {01/04/2018}
 */

#ifndef GDALIOH
#define GDALIOH

void gdalIORegister();
int getCellCenterLatLon(int inputEPSG, double xMin, double yMin, double xMax, double yMax, double cellSize, double ** px, double ** py);

#endif
