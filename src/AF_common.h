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

/*=====================================
 * Error logging macros.
 *
 * HDF_CHECK:
 *  stmt = statement to evaluate
 *  ret_type = Return type of the statement
 *  err_cond = The error condition. For instance, if function returns 
 *             a negative value on error, you would insert "< 0"
 *             (without quotes) as err_cond.
 * example:
 *      hid_t retval;
 *      HDF_CHECK( retval = H5Fopen( "nowhere", H5F_ACC_RDONLY, H5P_DEFAULT ), 
 *          hid_t, < 0 );
 *
 * FATAL_MSG:
 *  cerr_stream = C++ style err stream statement.
 *
 * example:
 *  FATAL_MSG( "Function returned " << ret_val << std::endl );
 */
#define ERR_MSG( cerr_stream ) \
    do { \
        std::cerr << "[" << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << "()] ERROR: " << std::endl; \
        std::cerr << cerr_stream ; \
        } while(0)

#define CALL_CHECK( stmt, ret_type, err_cond ) \
do { \
    ret_type __ret_val = stmt; \
    if(__ret_val err_cond ) { \
        ERR_MSG( "\tOn call: " << #stmt << std::endl << "\treturned: " << __ret_val << std::endl ); \
    } \
    } while(0)


/*======================================
 * HDF5 output group or dset names
 */
const std::string SRC_DATA_GROUP = "/Source/Data_Fields";
const std::string TRG_DATA_GROUP = "/Target/Data_Fields";
const std::string MODIS_RADIANCE_DSET = "MODIS_Radiance";
const std::string MISR_RADIANCE_DSET = "MISR_Radiance";
const std::string ASTER_RADIANCE_DSET = "ASTER_Radiance";
const std::string ASTER_SD_DSET = "ASTER_SD";  // Standard Deviation
const std::string ASTER_COUNT_DSET = "ASTER_Count"; // Pixel Count

#endif // _AF_COMMON_H_ 
