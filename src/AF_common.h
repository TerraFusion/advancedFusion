#ifndef _AF_COMMON_H_
#define _AF_COMMON_H_
/*********************************************************************
 * DESCRIPTION:
 *  Add here things common to share among all components
 *
 * DEVELOPER:
 *  - Jonathan Kim (jkm@illinois.edu)
 *
 */

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
const std::string ASTER_SD_DSET = "ASTER_SD";  // Standard Deviation
const std::string ASTER_COUNT_DSET = "ASTER_Count";

#endif // _AF_COMMON_H_ 
