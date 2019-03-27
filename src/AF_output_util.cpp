
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

std::string get_gtiff_fname(AF_InputParmeterFile &inputArgs,int camera_index,int band_index) {
	
	std::string err_fname ="";
	if(band_index <0){
		std::cerr<<"Insrument band index is less than 0. "<<std::endl;
		return err_fname;

	} 
	std::string gfname = inputArgs.GetOuputFilePath();
	size_t suffix_pos = gfname.rfind(".h5");
	if(suffix_pos != std::string::npos) 
		if(suffix_pos == (gfname.size()-3))
			gfname = gfname.substr(0,gfname.size()-3);
			
	std::string band_camera_info;
	std::string band_info;
	std::string camera_info;
	std::vector<std::string> band_names;
	std::vector<std::string>camera_names;

	if("MODIS" == inputArgs.GetSourceInstrument()) {
		if(true == inputArgs.IsMODIS_AllBands()){
			std::string modis_resolution = inputArgs.GetMODIS_Resolution();
			if("_1KM" == modis_resolution){
				band_names = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13L", "13H", "14L", "14H", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36"};
			}
			else if("_500m" == modis_resolution) {
				band_names = {"1", "2", "3","4", "5", "6", "7"};
			}
			else if("_250m" == modis_resolution){
				band_names = {"1", "2"};
			}

		}
		else {
			band_names = inputArgs.GetMODIS_Bands();
		}
	}
	else if("MISR" == inputArgs.GetSourceInstrument()) {
		band_names = inputArgs.GetMISR_Radiance();
		camera_names = inputArgs.GetMISR_CameraAngles();
	}
	else if("ASTER" == inputArgs.GetSourceInstrument()) {
		band_names = inputArgs.GetASTER_Bands();
        }
	if(band_index<band_names.size())
		band_info = band_names[band_index]; 
	else {
		std::cerr<<"Insrument band index is beyond the limit. band index  is "<< band_index<<std::endl;
		std::cerr<<"The limit is "<<band_names.size()-1<<std::endl;
		return err_fname;
	}
	if("MISR" == inputArgs.GetSourceInstrument()) {
		if((camera_index<camera_names.size()) && (camera_index >=0))
			camera_info = camera_names[camera_index]; 
		else {
			std::cerr<<"MISR camera index is beyond the limit. camera index  is "<< camera_index<<std::endl;
			std::cerr<<"The limit is "<<camera_names.size()-1<<std::endl;
			return err_fname;
		}
		band_camera_info = "_" + camera_info;
	}
	band_camera_info +="_";
	band_camera_info +=band_info;

	gfname += band_camera_info;
	gfname += ".tif";
		
	return gfname;

}


