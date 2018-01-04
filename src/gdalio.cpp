/**
 * gdalio.cpp
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {01/04/2018}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gdal.h>
#include <ogr_srs_api.h>
#include <ogr_api.h>
#include <cpl_conv.h>

void gdalIORegister() 
{	
	GDALAllRegister();
	OGRRegisterAll();
}

int getCellCenterLatLon(int inputEPSG, double xMin, double yMin, double xMax, double yMax, double cellSize, double ** px, double ** py) 
{

	int nRow = ceil((yMax - yMin) / cellSize);
	int nCol = ceil((xMax - xMin) / cellSize);

	yMax = yMin + nRow * cellSize;
	xMax = xMin + nCol * cellSize;

	int nPoints = nRow * nCol;

	if(NULL == (*px = (double *)malloc(sizeof(double) * nPoints))) 
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		printf("The number of output cells : %d\n may be too large\n", nPoints);
		exit(1);	
	}
	
	if(NULL == (*py = (double *)malloc(sizeof(double) * nPoints))) 
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		printf("The number of output cells : %d\n may be too large\n", nPoints);
		exit(1);	
	}

	x = *px;
	y = *py;

	int i, j;
	double rowY;

	for(i = 0; i < nRow; i++)
	{
		rowY = yMax - cellSize * (i + 0.5);
		for(j = 0; j < nCol; j++) 
		{
			x[i * nCol + j] = xMin + cellSize * (j + 0.5);
			y[i * nCol + j] = rowY;
		}
	}

	OGRSpatialReferenceH sourceSRS, targetSRS;
	OGRCoordinateTransformationH cTransform;

	sourceSRS = OSRNewSpatialReference(NULL);
	targetSRS = OSRNewSpatialReference(NULL);

	OSRImportFromEPSG(sourceSRS, inputEPSG);
	OSRImportFromEPSG(targetSRS, 4326);

	cTransform = OCTNewCoordinateTransformation(sourceSRS, targetSRS);

	OCTTransform(cTransform, nPoints, x, y, NULL);		

	OCTDestroyCoordinateTransformation(cTransform);	

	return nPoints;

}
