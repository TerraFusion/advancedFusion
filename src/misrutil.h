/**
 * misrutil.h
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {02/12/2018}
 */

#ifndef MISRUTILH
#define MISRUTILH
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
void getMISRFinalImageSize(int * pNRow, int * pNCol, int highResolution);

/**
 * NAME:	MISRBlockOffset
 * DESCRIPTION:	Perform MISR block offsets. This needs to be done for both geolocations and radiance values. 
 *		Please use this function to generate the array for MISR output if block offset needs to apply
 * PARAMETERS:
 *	double * originalGrid:	the original input grids (radiance values in one band, latitude or longitude)
 *	int highResolution: whether the MISR image is high or low resolution
 *		0: low resolution
 *		1: high resolution
 * Return:
 *	double *: the grid values after block offset
 */
double * MISRBlockOffset(double * originalGrid, int highResolution);

#endif
