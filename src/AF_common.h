#ifndef _AF_COMMON_H_
#define _AF_COMMON_H_

#include <iostream>

#include "AF_debug.h"

#define SUCCEED 0
#define FAILED -1

/*------------------------
 * hdf5 related names
 */
const std::string SRC_DATA_GROUP = "/Source/Data_Fields";
const std::string TRG_DATA_GROUP = "/Target/Data_Fields";
const std::string MODIS_RADIANCE_DSET = "MODIS_Radiance";
const std::string MISR_RADIANCE_DSET = "MISR_Radiance";
const std::string ASTER_RADIANCE_DSET = "ASTER_Radiance";

#endif // _AF_COMMON_H_ 
