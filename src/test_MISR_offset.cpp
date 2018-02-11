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

int main(int argc, char ** argv)
{
	char* file_path = "/projects/sciteam/jq0/TerraFusion/yizhao/TERRA_BF_L1B_O69626_20130119123228_F000_V001.h5";
	hid_t src_file;
	if(0 > (src_file = af_open(file_path))) {
		printf("File not found\n");
		exit(1);
	}

	int nCellMISR;
	double * MISRRad;
	MISRRad = get_misr_rad(src_file, "AN", "L", "Blue_Radiance", &nCellMISR);
	herr_t ret = af_close(src_file);

	printf("Total number of MISR pixels: %d\n", nCellMISR);

	free(MISRRad);
	
	return 0;
}
