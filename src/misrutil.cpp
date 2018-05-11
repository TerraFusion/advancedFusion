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

