#ifndef _AF_OUTPUT_MODIS_H_
#define _AF_OUTPUT_MODIS_H_

#include "AF_InputParmeterFile.h"
#include <hdf5.h>

//  MODIS as target instrument.  Generate output to hdf5 file
int af_GenerateOutputCumulative_ModisAsTrg(AF_InputParmeterFile &inputArgs, hid_t outputFile,hid_t srcFile, int trgCellNum, std::map<std::string, strVec_t> &inputMultiVarsMap);


#endif // _AF_OUTPUT_MODIS_H_
