#ifndef _AF_OUTPUT_MISR_H_
#define _AF_OUTPUT_MISR_H_
/*********************************************************************
 * DESCRIPTION:
 *  Generate radiance data output to HDF5 for MISR
 *
 *
 * DEVELOPERS:
 *  - Jonathan Kim (jkm@illinois.edu)
 */


#include "AF_InputParmeterFile.h"
#include <hdf5.h>

//  MODIS as Target instrument, generate radiance data
int af_GenerateOutputCumulative_MisrAsTrg(AF_InputParmeterFile &inputArgs, hid_t outputFile,hid_t srcFile, int trgCellNumOri, std::map<std::string, strVec_t> &inputMultiVarsMap);


//  MODIS as Source instrument, generate radiance data
int af_GenerateOutputCumulative_MisrAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNum, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap);

#endif // _AF_OUTPUT_MISR_H_
