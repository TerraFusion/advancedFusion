/**
 * gdalio.cpp
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {05/25/2018}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gdal.h>
#include <ogr_srs_api.h>
#include <ogr_api.h>
#include <cpl_conv.h>
#include <omp.h>

/**
 * NAME:	gdalIORegister
 * DESCRIPTION:	Register drivers
 * PARAMETERS: void
 */
void gdalIORegister() 
{	
	GDALAllRegister();
	OGRRegisterAll();
}

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
 * Return:
 *	int:	the total number of pixels
 */
int getCellCenterLatLon(int outputEPSG, double xMin, double yMin, double xMax, double yMax, double cellSize, double ** px, double ** py) 
{

	if(yMax <=yMin || xMax <=xMin || cellSize <0)
		return -1;

	int nRow = ceil((yMax - yMin) / cellSize);
	int nCol = ceil((xMax - xMin) / cellSize);

	yMax = yMin + nRow * cellSize;
	xMax = xMin + nCol * cellSize;

	int nPoints = nRow * nCol;

	if(NULL == (*px = (double *)malloc(sizeof(double) * nPoints))) 
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		printf("The number of output cells : %d\n may be too large\n", nPoints);
		return -1;
	}
	
	if(NULL == (*py = (double *)malloc(sizeof(double) * nPoints))) 
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		printf("The number of output cells : %d\n may be too large\n", nPoints);
		return -1;	
	}

	double * x = *px;
	double * y = *py;

	int i, j;
	double rowY;

#pragma omp parallel for private(j, rowY)
	for(i = 0; i < nRow; i++)
	{
/*
		if(i == 0)
		{
			printf("%d threads\n", omp_get_num_threads());
		}
*/
		rowY = yMax - cellSize * (i + 0.5);
		for(j = 0; j < nCol; j++) 
		{
			x[i * nCol + j] = xMin + cellSize * (j + 0.5);
			y[i * nCol + j] = rowY;
		}
	}

	if(outputEPSG != 4326) 
	{
#pragma omp parallel
		{
			int nThreads = omp_get_num_threads();
			int threadID = omp_get_thread_num();


			int start = nPoints / nThreads * threadID;
			int nPointsThread;
			if(threadID != nThreads - 1)
			{
			 	nPointsThread = nPoints / nThreads;
			}
			else
			{
				nPointsThread = nPoints - start;
			}

			//TODO: Need to check error handling

			OGRSpatialReferenceH sourceSRS, targetSRS;
			OGRCoordinateTransformationH cTransform;

			sourceSRS = OSRNewSpatialReference(NULL);
			targetSRS = OSRNewSpatialReference(NULL);

			OSRImportFromEPSG(sourceSRS, outputEPSG);
			OSRImportFromEPSG(targetSRS, 4326);

			cTransform = OCTNewCoordinateTransformation(sourceSRS, targetSRS);

			OCTTransform(cTransform, nPointsThread, x + start, y + start, NULL);		

			OCTDestroyCoordinateTransformation(cTransform);

            // Forgot destroy sourceSRS and targetSRS-memory leaking. destroy here
            OSRDestroySpatialReference(sourceSRS);
            OSRDestroySpatialReference(targetSRS);

		//	printf("%d of %d: %d - %d\n", threadID, nThreads, start, nPointsThread);

		}

	}

	return nPoints;

}

/**
 * NAME:	writeGeoTiff
 * DESCRIPTION:	Write the output grid as a GeoTiff
 * PARAMETERS:
 *	char * fileName:	output GeoTiff file name
 * 	double * grid:		the grid of the output radianc values
 * 	int outputEPSG:		EPSG code of output spatial reference system (negative value if unknown) 
 *	double xMin:		west boundary of output area
 *	double yMin:		south boundary of output area
 *	double xMax:		east boundary of output area
 *	double yMax: 		north boundary of output area
 * 	double cellSize:	output raste cell size
 */
void writeGeoTiff(char * fileName, double * grid, int outputEPSG, double xMin, double yMin, double xMax, double yMax, double cellSize)
{	
	int nRow = ceil((yMax - yMin) / cellSize);
	int nCol = ceil((xMax - xMin) / cellSize);

	yMax = yMin + nRow * cellSize;
	xMax = xMin + nCol * cellSize;

	GDALDriverH hDriver;
	if(NULL == (hDriver = GDALGetDriverByName("GTiff")))
	{
		printf("ERROR: cannot get driver for GTiff\n");
		exit(1);
	}

	GDALDatasetH hDstDS;
	char *papszOptions[] = {"COMPRESS=LZW",NULL};
	hDstDS = GDALCreate(hDriver, fileName, nCol, nRow, 1, GDT_Float64, papszOptions);
	
	double adfGeoTransform[6];
	adfGeoTransform[0] = xMin;
	adfGeoTransform[1] = cellSize;
	adfGeoTransform[2] = 0;
	adfGeoTransform[3] = yMax;
	adfGeoTransform[4] = 0;
	adfGeoTransform[5] = -cellSize;

	GDALSetGeoTransform(hDstDS,adfGeoTransform);

	if(outputEPSG > 0)
	{
		char *pszSRS_WKT = NULL;
		OGRSpatialReferenceH hSRS = OSRNewSpatialReference(NULL);
		OSRImportFromEPSG(hSRS, outputEPSG);
		OSRExportToWkt(hSRS,&pszSRS_WKT);
		GDALSetProjection(hDstDS,pszSRS_WKT);
		OSRDestroySpatialReference(hSRS);
		CPLFree(pszSRS_WKT);
	}

	GDALRasterBandH hBand;
	hBand=GDALGetRasterBand(hDstDS,1);
	GDALSetRasterNoDataValue(hBand,-999.0);
	GDALRasterIO(hBand, GF_Write, 0, 0, nCol, nRow, grid, nCol, nRow, GDT_Float64, 0, 0 );

	GDALClose(hDstDS);

	return;
}

/**
 * NAME:	getMaxRadiusOfUserdefine
 * DESCRIPTION:	Get the maximum distance (in meters) for user-defined-grid to be used in "nearestNeighbor" when using summary interpolate
 * PARAMETERS:
 *	int epsgCode:		EPSG code of the spatial reference system
 * 	double cellSize:	Raste cell size
 * Return:
 *	double:		the maximum distance (in meters) to be used in "nearestNeighbor"
 */
double getMaxRadiusOfUserdefine(int epsgCode, double cellSize) {

        const double earthRadius = 6367444;

        OGRSpatialReferenceH hSRS = OSRNewSpatialReference(NULL);
        OSRImportFromEPSG(hSRS, epsgCode);

        char *SRSWKT = NULL;
        OSRExportToWkt(hSRS,&SRSWKT);
        //printf("DBG Output SRS Info:\n\tEPSG CODE: %d\n\t%s\n", epsgCode, SRSWKT);

        double maxRadius;

        int isProjected = OSRIsProjected(hSRS);
        if(isProjected) {
                //printf("\tIs a Projected Cooridinate System.\n");

                char * unitName = NULL;
                double unitConversion;
                unitConversion = OSRGetLinearUnits(hSRS, &unitName);
                //printf("\tUnit Name: %s (%lf)\n", unitName, unitConversion);

                maxRadius = cellSize * unitConversion;

        }
        else {
                //printf("\tIs a Geographic Cooridinate System.\n");

                maxRadius = earthRadius * cellSize * M_PI / 180;
        }
        //printf("DBG> Max Radius for Resampling: %lf meters\n", maxRadius);
        return maxRadius;
}



