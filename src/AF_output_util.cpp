
/*********************************************************************
 * DESCRIPTION:
 *   Add functions to be used from data-output component
 *
 *
 * DEVELOPERS:
 *  - Jonathan Kim (jkm@illinois.edu)
 */

#include "AF_output_util.h"

#include <math.h>

#include "AF_common.h"
#include "misrutil.h"

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
int af_GetWidthAndHeightForOutputDataSize(std::string instrument, AF_InputParmeterFile &inputArgs, int &crossTrackWidth /*OUT*/, int &alongTrackHeight /*OUT*/)
{
	int ret = 0;
	crossTrackWidth = 0;
	alongTrackHeight = 0;
	std::string misrShift = inputArgs.GetMISR_Shift();
	std::string trgInstrument = inputArgs.GetTargetInstrument();

	// MISR shift case. Shift always if misr is target and Shift is On for all source instrument and target geolocation data
	if(misrShift == "ON" && trgInstrument == MISR_STR) {
		std::string resolution = inputArgs.GetMISR_Resolution();
		getMISRFinalImageSize(&alongTrackHeight, &crossTrackWidth, (resolution=="L") ? 0 : 1);
	}
	/*-------------------------------------------------------
	 * MODIS section
	 */
	else if (instrument == MODIS_STR) {
		std::string resolution = inputArgs.GetMODIS_Resolution();
		if (resolution == "_1KM") {
			crossTrackWidth = 1354;
		}
		else if (resolution == "_500m") {
			crossTrackWidth = 2708;  // 1354 * 2
		}
		else if (resolution == "_250m") {
			crossTrackWidth = 5416;  // 1354 * 4
		}
	}
	/*-------------------------------------------------------
	 * MISR section
	 */
	else if (instrument == MISR_STR) {
		std::string resolution = inputArgs.GetMISR_Resolution();
		if (resolution == "L") {  // 1.1KM
			crossTrackWidth = 512;
		}
		else if (resolution == "H") {  // 275M
			crossTrackWidth = 2048;
		}
	}
	/*-------------------------------------------------------
	 * USER_DEFINE section
	 */
	else if (instrument == USERGRID_STR) {
		double xMin = inputArgs.GetUSER_xMin();
		double xMax = inputArgs.GetUSER_xMax();
		double cellSize = inputArgs.GetUSER_Resolution();
		crossTrackWidth = ceil((xMax - xMin) / cellSize);
	}
	else {
		return -1;  // fail
	}

	#if DEBUG_TOOL
	std::cout << "DBG_TOOL " << __FUNCTION__ << "> misrShift: " << misrShift << ", instrument: " << instrument << ", crossTrackWidth: " << crossTrackWidth << ", alongTrackHeight: " << alongTrackHeight << "\n";
	#endif

	return 0;
}
