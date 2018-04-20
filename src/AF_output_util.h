#ifndef _AF_OUTPUT_UTIL_H_
#define _AF_OUTPUT_UTIL_H_
/*********************************************************************
 * DESCRIPTION:
 *   Add functions to be used from data-output component
 *
 *
 * DEVELOPERS:
 *  - Jonathan Kim (jkm@illinois.edu)
 */

#include "AF_InputParmeterFile.h"

/*=====================================
 * Get output width of an instrument
 *
 * RETURN:
 *  0 : SUCCEED
 * -1 : FAILED
 *
 * OUT parameters:
 * - crossTrackWidth : return width value.
 * - alongTrackHeight : return height value.
 *						This is only valid for Misr target shift case.
 *						If 0, caller should not use this.
 */
int af_GetWidthAndHeightForOutputDataSize(std::string instrument, AF_InputParmeterFile &inputArgs, int &crossTrackWidth /*OUT*/, int &alongTrackHeight /*OUT*/);

#endif // _AF_OUTPUT_UTIL_H_
