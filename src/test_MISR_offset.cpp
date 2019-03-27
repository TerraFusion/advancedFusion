/*
    AUTHOR:
        Yizhao Gao
    EMAIL:
        ygao29@illinois.edu
	PROGRAM DESCRIPTION:
		For testing MISR block offset functionalities
*/

#include <stdio.h>
#include <stdlib.h>
#include "io.h"
#include "misrutil.h"
#include "gdalio.h"

int main(int argc, char ** argv)
{
	char* file_path = "/projects/sciteam/jq0/TerraFusion/yizhao/TERRA_BF_L1B_O69626_20130119123228_F000_V001.h5";
	hid_t src_file;
	if(0 > (src_file = af_open(file_path))) {
		printf("File not found\n");
		exit(1);
	}

	int64_t nCellMISR;
	double * MISRRad;
	MISRRad = get_misr_rad(src_file, "AN", "L", "Blue_Radiance", &nCellMISR);
	herr_t ret = af_close(src_file);

	printf("Total number of MISR pixels: %d\n", nCellMISR);

	int nRow, nCol;
	getMISRFinalImageSize(&nRow, &nCol, 0);
	printf("Size of the final output image: %d * %d\n", nRow, nCol);

	double * MISRFinal;

	MISRFinal = (double *) malloc(sizeof(double) * nRow * nCol);

	MISRBlockOffset<double>(MISRRad, MISRFinal, 0); 
	
	gdalIORegister();
	writeGeoTiff("TestMISROffset.tif", MISRFinal, -1, 0, 0, (double)(nCol), (double)(nRow), 1);

	free(MISRRad);
	free(MISRFinal);
	
	return 0;
}
