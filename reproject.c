/**
 * reproject.c
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {07/17/2017}
 */


#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

int * pointIndexOnLat(double ** plat, double ** plon,  int * oriID, int count, int nBlockY) {

	double *lat = *plat;
	double *lon = *plon;

	double blockR = M_PI/nBlockY;

	int * index;
	int * pointsInB;
	
	double * newLon;
	double * newLat;

	if(NULL == (index = (int *)malloc(sizeof(int) * (nBlockY + 1))))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	if(NULL == (pointsInB = (int *)malloc(sizeof(int) * nBlockY)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	int i;
	for(i = 0; i < nBlockY; i++)
	{
		pointsInB[i] = 0;
	}

	int blockID;
	int j;	
	for(j = 0; j < count; j++) {
	
		blockID = (int)((lat[j] + M_PI/2) / blockR);
		
		if(blockID >= 0 && blockID < nBlockY) {
			pointsInB[blockID] ++;
		}
	}

	index[0] = 0;
	int k;
	for(k = 1; k < nBlockY + 1; k++) {
		index[k] = index[k - 1] + pointsInB[k - 1];
	}

	if(NULL == (newLon = (double *)malloc(sizeof(double) * index[nBlockY]))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	} 
	if(NULL == (newLat = (double *)malloc(sizeof(double) * index[nBlockY]))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	
	pointsInB[0] = 0;
	int m;
	for(m = 1; m < nBlockY; m++) {
		pointsInB[m] = index[m];
	}
	int n;
	for(n = 0; n < count; n++) {
		
		blockID = (int)((lat[n] + M_PI/2) / blockR);
		if(blockID >= 0 && blockID < nBlockY) {
			newLon[pointsInB[blockID]] = lon[n];
			newLat[pointsInB[blockID]] = lat[n];
			oriID[pointsInB[blockID]] = n;
			pointsInB[blockID] ++;
	
		}
	}
	


	free(pointsInB);
	free(lon);
	free(lat);

	count = index[nBlockY];
	*plon = newLon;
	*plat = newLat;

//	printf("%0x\n", lat);

	return index;
}

//Finding the nearest neiboring point's ID 
void nearestNeighbor(double ** psouLat, double ** psouLon, int nSou, double * tarLat, double * tarLon, int * tarNNSouID, int nTar, double maxR) {

	//printf("%0x\n", souLat);
	double * souLat = *psouLat;
	double * souLon = *psouLon;

	const double earthRadius = 6367444;
	double radius = maxR / earthRadius;
	int nBlockY = M_PI / radius;

	double blockR = M_PI / nBlockY;
	
	int i;
	for(i = 0; i < nSou; i++) {
		souLat[i] = souLat[i] * M_PI / 180;
		souLon[i] = souLon[i] * M_PI / 180;
	}
	int j;
	for(j = 0; j < nTar; j++) {
		tarLat[j] = tarLat[j] * M_PI / 180;
		tarLon[j] = tarLon[j] * M_PI / 180;
	}

	int * souID;
	if(NULL == (souID = (int *)malloc(sizeof(double) * nSou))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	int * souIndex = pointIndexOnLat(psouLat, psouLon, souID, nSou, nBlockY);
	souLat = *psouLat;
	souLon = *psouLon;

	double tLat, tLon;
	double sLat, sLon;
	int blockID;
	int startBlock, endBlock;

	double pDis;	
	double nnDis;
	int nnSouID;
	int m;
	for(m = 0; m < nTar; m ++) {

		tLat = tarLat[m];
		tLon = tarLon[m];

		blockID = (tLat + M_PI / 2) / blockR;
		startBlock = blockID - 1;
		endBlock = blockID + 1;

		if(startBlock < 0) {
			startBlock = 0;
		}
		if(endBlock > nBlockY - 1) {
			endBlock = nBlockY - 1;
		}

		nnDis = -1;
		int n;
		for(n = souIndex[startBlock]; n < souIndex[endBlock+1]; n++) {
			
			sLat = souLat[n];
			sLon = souLon[n];

			pDis = acos(sin(tLat) * sin(sLat) + cos(tLat) * cos(sLat) * cos(tLon - sLon));
				
			if((nnDis < 0 || nnDis > pDis) && pDis <= radius) {
				nnDis = pDis;
				nnSouID = souID[n];
			}
				
		}

		if(nnDis < 0) {
			tarNNSouID[m] = -1;
		}
		else {
			tarNNSouID[m] = nnSouID;
		}
	
		 
	}
	
	return;	
}


void nnInterpolate(double * souVal, double * tarVal, int * tarNNSouID, int nTar) {

	int nnSouID;
	int i;
	for(i = 0; i < nTar; i++) {
		nnSouID = tarNNSouID[i];
		if(nnSouID < 0) {
			tarVal[i] = -999;
		}
		else {
			tarVal[i] = souVal[nnSouID];
		}
	}
}

void summaryInterpolate(double * souVal, int * souNNTarID, int nSou, double * tarVal, int * nSouPixels, int nTar) {
	int i;
	for(i = 0; i < nTar; i++) {
	
		tarVal[i] = 0;
		nSouPixels[i] = 0;
	}

	int nnTarID;
	int j;
	for(j = 0; j < nSou; j++) {
		
		nnTarID = souNNTarID[j];
		if(nnTarID > 0 && souVal[j] >= 0) {
			tarVal[nnTarID] += souVal[j];
			nSouPixels[nnTarID] ++;	
		}
	}

	int k;
	for(k = 0; k < nTar; k++) {
	
		if(nSouPixels[k] > 0) {
			tarVal[k] = tarVal[k] / nSouPixels[k];
		}
		else {
			tarVal[k] = -999;
		}
	}

}
