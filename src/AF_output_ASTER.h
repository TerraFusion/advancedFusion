#ifndef _AF_OUTPUT_ASTER_H_
#define _AF_OUTPUT_ASTER_H_

/*********************************************************************
 * DESCRIPTION:
 *  Generate radiance data output to HDF5 for ASTER
 *
 *
 * DEVELOPERS:
 *  - Jonathan Kim (jkm@illinois.edu)
 */

#include "AF_InputParmeterFile.h"
#include <hdf5.h>

//  ASTER as Target instrument, generate radiance data
//  TBD


//  ASTER as Source instrument, generate radiance data
int af_GenerateOutputCumulative_AsterAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNumNoShift, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap);


#endif // _AF_OUTPUT_ASTER_H
