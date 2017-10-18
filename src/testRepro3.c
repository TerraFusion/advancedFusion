/**
 * testRepro3.c
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {07/14/2017}
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "reproject.h"

int main(int argc, char ** argv) {

	double * MODISLat, * MODISLon, * MODISVal;
	double * MISRLat, * MISRLon, * MISRVal;
	double ** pMODISLat, ** pMODISLon;
	double ** pMISRLat, ** pMISRLon;
	int * souNNTarID;
	int nMODIS = 1347102;
	int nMISR = 262144;

	double maxR = 1000;

	MODISLat = (double *) malloc(sizeof(double) * nMODIS);
	MODISLon = (double *) malloc(sizeof(double) * nMODIS);
	MODISVal = (double *) malloc(sizeof(double) * nMODIS);

	MISRLat = (double *) malloc(sizeof(double) * nMISR);
	MISRLon = (double *) malloc(sizeof(double) * nMISR);
	MISRVal = (double *) malloc(sizeof(double) * nMISR);
	
	pMODISLat = &MODISLat;
	pMODISLon = &MODISLon;
	pMISRLat = &MISRLat;
	pMISRLon = &MISRLon;
	
	FILE * fLat;
	FILE * fLon;
	FILE * fVal;

	fLat = fopen("data/MODIS_Latitude-sub", "r");
	fLon = fopen("data/MODIS_Longitude-sub", "r");
	fVal = fopen("data/MODIS_EV_250_Aggr500_RefSB-sub", "r");

	for(int i = 0; i < nMODIS; i++) {
		
		fscanf(fLat, "%lf,\n", MODISLat + i);
		fscanf(fLon, "%lf,\n", MODISLon + i);
		fscanf(fVal, "%lf,\n", MODISVal + i);
	}


	fclose(fLat);
	fclose(fLon);
	fclose(fVal);
	
	fLat = fopen("data/MISR_Latitude-sub", "r");
	fLon = fopen("data/MISR_Longitude-sub", "r");
	fVal = fopen("data/MISR_Blue_Radiance-sub", "r");

	for(int i = 0; i < nMISR; i++) {
		
		fscanf(fLat, "%lf,\n", MISRLat + i);
		fscanf(fLon, "%lf,\n", MISRLon + i);
		fscanf(fVal, "%lf,\n", MISRVal + i);
	}


	fclose(fLat);
	fclose(fLon);
	fclose(fVal);

	//MISR to MODIS
/*	
	printf("Lat,Lon,Val\n");
	for(int i = 0; i < nMISR; i++) {
		printf("%lf,%lf,%lf\n", MISRLat[i], MISRLon[i], MISRVal[i]);
	}

	souNNTarID = (int *) malloc(sizeof(int) * nMISR);
	
	nearestNeighbor(pMODISLat, pMODISLon, nMODIS, MISRLat, MISRLon, souNNTarID, nMISR, maxR);
	
	MODISLat = *pMODISLat;
	MODISLon = *pMODISLon;

	int * nMISRPixels = (int *) malloc(sizeof(int) * nMODIS);

	summaryInterpolate(MISRVal, souNNTarID, nMISR, MODISVal, nMISRPixels, nMODIS);

	fLat = fopen("data/MODIS_Latitude-sub", "r");
	fLon = fopen("data/MODIS_Longitude-sub", "r");

	for(int i = 0; i < nMODIS; i++) {
		
		fscanf(fLat, "%lf,\n", MODISLat + i);
		fscanf(fLon, "%lf,\n", MODISLon + i);
	}


	fclose(fLat);
	fclose(fLon);

	printf("Lat,Lon,Val\n");
	for(int i = 0; i < nMODIS; i++) {
		printf("%lf,%lf,%d,%lf\n", MODISLat[i], MODISLon[i], nMISRPixels[i], MODISVal[i]);
	}

	free(nMISRPixels);
*/
	//MODIS to MISR

	printf("Lat,Lon,Val\n");
	for(int i = 0; i < nMODIS; i++) {
		printf("%lf,%lf,%lf\n", MODISLat[i], MODISLon[i], MODISVal[i]);
	}

	souNNTarID = (int *) malloc(sizeof(int) * nMODIS);
	
	nearestNeighbor(pMISRLat, pMISRLon, nMISR, MODISLat, MODISLon, souNNTarID, nMODIS, maxR);
	
	MISRLat = *pMISRLat;
	MISRLon = *pMISRLon;
	
	int * nMODISPixels = (int *) malloc(sizeof(int) * nMISR);

	summaryInterpolate(MODISVal, souNNTarID, nMODIS, MISRVal, nMODISPixels, nMISR);
	
	fLat = fopen("data/MISR_Latitude-sub", "r");
	fLon = fopen("data/MISR_Longitude-sub", "r");

	for(int i = 0; i < nMISR; i++) {
		
		fscanf(fLat, "%lf,\n", MISRLat + i);
		fscanf(fLon, "%lf,\n", MISRLon + i);
	}


	fclose(fLat);
	fclose(fLon);

	printf("Lat,Lon,Count,Val\n");
	for(int i = 0; i < nMISR; i++) {
		printf("%lf,%lf,%d,%lf\n", MISRLat[i], MISRLon[i],  nMODISPixels[i], MISRVal[i]);
	}

	free(nMODISPixels);

	//MODIS to MISR
	

	free(*pMODISLat);
	free(*pMODISLon);
	free(MODISVal);
	free(*pMISRLat);
	free(*pMISRLon);
	free(MISRVal);
	free(souNNTarID);

	return 0;
}
