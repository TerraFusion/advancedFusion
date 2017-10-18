/**
 * reproject.h
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {10/10/2017}
 */

#ifndef REPROH
#define REPROH


/**
 * NAME:	nearestNeighborBlockIndex
 * DESCRIPTION:	Find the nearest neighboring source cell's ID for each target cell
 * PARAMETERS:
 *	double ** psouLat:	the pointer to the array of latitudes of source cells (the data are changed during in the function, so please do the output before this function)
 *	double ** psouLon:	the pointer to the array of longitudes of source cells (the data are changed during in the function, so please do the output before this function)
 *	int nSou:		the number of source cells
 *	double * tarLat:	the latitudes of target cells
 *	double * tarLon:	the longitudes of target cells
 *	int * tarNNSouID:	the output IDs of nearest neighboring source cells 
 *	double * tarNNDis	the output nearest distance for each target cell (input NULL if you don't need this field)
 *	int nTar:		the number of target cells
 *	double maxR:		the maximum distance (in meters) to define neighboring cells
 * Output: 	
 *	int * tarNNSouID:	the output IDs of nearest neighboring source cells 
 *	double * tarNNDis	the output nearest distance for each target cell (input NULL if you don't need this field)
 */ 

void nearestNeighborBlockIndex(double ** psouLat, double ** psouLon, int nSou, double * tarLat, double * tarLon, int * tarNNSouID, double * tarNNDis, int nTar, double maxR);

/**
 * NAME:	nearestNeighbor
 * DESCRIPTION:	Find the nearest neighboring source cell's ID for each target cell
 * PARAMETERS:
 *	double ** psouLat:	the pointer to the array of latitudes of source cells (the data are changed during in the function, so please do the output before this function)
 *	double ** psouLon:	the pointer to the array of longitudes of source cells (the data are changed during in the function, so please do the output before this function)
 *	int nSou:		the number of source cells
 *	double * tarLat:	the latitudes of target cells
 *	double * tarLon:	the longitudes of target cells
 *	int * tarNNSouID:	the output IDs of nearest neighboring source cells 
 *	double * tarNNDis	the output nearest distance for each target cell (input NULL if you don't need this field)
 *	int nTar:		the number of target cells
 *	double maxR:		the maximum distance (in meters) to define neighboring cells
 * Output: 	
 *	int * tarNNSouID:	the output IDs of nearest neighboring source cells 
 *	double * tarNNDis	the output nearest distance for each target cell (input NULL if you don't need this field)
 */ 

void nearestNeighbor(double ** psouLat, double ** psouLon, int nSou, double * tarLat, double * tarLon, int * tarNNSouID, double * tarNNDis, int nTar, double maxR);





/**
 * NAME:	nnInterpolate
 * DESCRIPTION:	Nearest neighbor interpolation
 * PARAMETERS:
 * 	double * souVal:	the input values at source cells
 * 	double * tarVal:	the output values at target cells
 * 	int * tarNNSouID:	the IDs of nearest neighboring source cells for each target cells (generated from "nearestNeighbor") 
 *	int nTar:		the number of target cells
 * Output: 	
 * 	double * tarVal:	the output values at target cells
 */ 
void nnInterpolate(double * souVal, double * tarVal, int * tarNNSouID, int nTar);

/**
 * NAME:	summaryInterpolate
 * DESCRIPTION:	Interpolation (summary) from fine resolution to coarse resolution
 * PARAMETERS:
 * 	double * souVal:	the input values at source cells
 * 	int * souNNTarID:	the IDs of nearest neighboring target cells for each source cells (generated from "nearestNeighbor")
 * 	int nSou:		the number of source cells
 * 	double * tarVal:	the output values at target cells
 * 	int * nSouPixels:	the output numbers of contributing source cells to each target cell
 *	int nTar:		the number of target cells
 * Output:
 * 	double * tarVal:	the output values at target cells
 * 	int * nSouPixels:	the output numbers of contributing source cells to each target cell
 */
void summaryInterpolate(double * souVal, int * souNNTarID, int nSou, double * tarVal, int * nSouPixels, int nTar);

#endif


