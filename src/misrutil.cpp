/**
 * misrutil.cpp 
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {02/11/2018}
 */

// NOTE: dimensions of misr data
// Blocks per orbit: 180
// Pixels per block:
//	Low resolution (1.1km): 128 * 512
//	High resultion (275m): 512 * 2048

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/**
 * NAME:	getMISRFinalImageSize
 * DESCRIPTION:	Get the image size of MISR for one orbit after block offseting; results only depend on resolution (low or high)
 * PARAMETERS:
 * 	int * pNRow: the number of rows in the final image to be output
 *	int * pNCol: the number of columns in the final image to be ouput
 *	int highResolution: whether the MISR image is high or low resolution
 *		0: low resolution
 *		1: high resolution
 * OUTPUT:
 * 	int * pNRow: the number of rows in the final image to be output
 *	int * pNCol: the number of columns in the final image to be ouput
 */
void getMISRFinalImageSize(int * pNRow, int * pNCol, int highResolution) 
{
	if (highResolution == 0) 
	{
		*pNRow = 180 * 128;
		*pNCol = 512 + 1580; 
		//The number 1580 comes from the 179 block offset values
		//A prefix sum of the offset value is first calculated
		//The difference be tween the highest (64) and the lowest (-1520) prefix sum is 1580
	}
	else 
	{
		*pNRow = 180 * 512;
		*pNCol = 2048 + 1580 * 4;
	}
}

/**
 * NAME:	MISRBlockOffset
 * DESCRIPTION:	Perform MISR block offsets. This needs to be done for both geolocations and radiance values. 
 *		Please use this function to generate the array for MISR output if block offset needs to apply
 * PARAMETERS:
 *	double * originalGrid:	the original input grids (radiance values in one band, latitude or longitude)
 *	double * newGrid:	the output grid (the original input grid after block offseting)
 *	int highResolution: whether the MISR image is high or low resolution
 *		0: low resolution
 *		1: high resolution
 * Output:
 *	double * newGrid:	the output grid (the original input grid after block offseting)
 */
void MISRBlockOffset(double * originalGrid, double * newGrid, int highResolution) 
{
	int nRowPerBlock;
	int nColPerBlock;

	int nRow;
	int nCol;
	getMISRFinalImageSize(&nRow, &nCol, highResolution);
	
	int offsets[180] = {1520,1520,1536,1536,1552,1552,1552,1552,1568,1568,1568,1568,1568,1584,1584,1584,1584,1584,1584,1584,1584,1584,1584,1584,1584,1584,1584,1568,1568,1568,1568,1552,1552,1552,1536,1536,1536,1520,1520,1504,1504,1488,1488,1472,1456,1456,1440,1440,1424,1408,1408,1392,1376,1360,1360,1344,1328,1312,1296,1296,1280,1264,1248,1232,1216,1200,1184,1168,1152,1136,1120,1104,1088,1072,1056,1040,1024,1008,992,976,960,944,928,912,896,864,848,832,816,800,784,768,752,736,720,704,672,656,640,624,608,592,576,560,544,528,512,496,480,464,448,432,416,400,384,368,352,336,320,304,304,288,272,256,240,224,224,208,192,176,176,160,144,144,128,128,112,96,96,80,80,64,64,64,48,48,32,32,32,16,16,16,16,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,16,16,32,32,32,48,48};
	//This list of 180 offsets number are calculated from the original 179 offset values
	//First a zero (0) is added as the first item of the list to represent the offset of the first block
	//Then prefix sum of the offset value is calculate
	//Finally, the minium prefix sum value is substracted from each element of the prefix sum to make all values non-negative
	//The resulting list is the "int offsets[180]"

	int i, j;
	for(i = 0; i < nRow; i++)
	{
		for(j = 0; j < nCol; j ++)
		{
			newGrid[i * nCol + j] = -999;
		}
	}

	if (highResolution == 0)
	{
		nRowPerBlock = 128;
		nColPerBlock = 512;
	}
	else
	{
		nRowPerBlock = 512;
		nColPerBlock = 2048;
		for(i = 0; i < 180; i++)
		{
			offsets[i] *= 4;
		}	
	}

	int blockID;

#pragma omp parallel for private(blockID, j)
	for(i = 0; i < nRow; i++)
	{
		blockID = i / nRowPerBlock;
		for(j = 0; j < nColPerBlock; j++ )
		{
			newGrid[i * nCol + j + offsets[blockID]] = originalGrid[i * nColPerBlock + j];
		}
	}
}
