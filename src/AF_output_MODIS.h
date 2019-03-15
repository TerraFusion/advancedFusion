#ifndef _AF_OUTPUT_MODIS_H_
#define _AF_OUTPUT_MODIS_H_

/*********************************************************************
 * DESCRIPTION:
 *  Generate radiance data output to HDF5 for MODIS
 *
 *
 * DEVELOPERS:
 *  - Jonathan Kim (jkm@illinois.edu)
 */

#include "AF_InputParmeterFile.h"
#include <hdf5.h>
#include <hdf5_hl.h>

//  MODIS as Target instrument, generate radiance data
//int af_GenerateOutputCumulative_ModisAsTrg(AF_InputParmeterFile &inputArgs, hid_t outputFile,hid_t srcFile, int trgCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap,hid_t ctrackDset, hid_t atrackDset);
int af_GenerateOutputCumulative_ModisAsTrg(AF_InputParmeterFile &inputArgs, hid_t outputFile,hid_t srcFile, int64_t trgCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap,hid_t ctrackDset, hid_t atrackDset);

//  MODIS as Source instrument, generate radiance data
//int af_GenerateOutputCumulative_ModisAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int *targetNNsrcID,  int trgCellNumNoShift, hid_t srcFile, int srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap,hid_t ctrackDset, hid_t atrackDset);
int af_GenerateOutputCumulative_ModisAsSrc(AF_InputParmeterFile &inputArgs, hid_t outputFile, int64_t *targetNNsrcID,  int64_t trgCellNumNoShift, hid_t srcFile, int64_t srcCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap,hid_t ctrackDset, hid_t atrackDset);


#endif // _AF_OUTPUT_MODIS_H_
