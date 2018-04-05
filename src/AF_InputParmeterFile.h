
#include <vector>
#include <map>
#include <string>

#include "AF_debug.h"

/*======================================================
 * Input parameter entry strings
 */
const std::string INPUT_FILE_PATH="INPUT_FILE_PATH";
const std::string OUTPUT_FILE_PATH="OUTPUT_FILE_PATH";
const std::string RESAMPLE_METHOD="RESAMPLE_METHOD";
const std::string SOURCE_INSTRUMENT="SOURCE_INSTRUMENT";
const std::string TARGET_INSTRUMENT="TARGET_INSTRUMENT";
// MISR section ---------------------
const std::string MISR_RESOLUTION="MISR_RESOLUTION";
const std::string MISR_CAMERA_ANGLE="MISR_CAMERA_ANGLE";
const std::string MISR_RADIANCE="MISR_RADIANCE";
const std::string MISR_SHIFT="MISR_TARGET_BLOCKUNSTACK";
// MODIS section -------------------
const std::string MODIS_RESOLUTION="MODIS_RESOLUTION";
const std::string MODIS_BANDS="MODIS_BANDS";
#if 1 // JK_ASTER2MODIS
// ASTER section ----------------
const std::string ASTER_RESOLUTION="ASTER_RESOLUTION";
const std::string ASTER_BANDS="ASTER_BANDS";
#endif
// USER_DEFINE section ---------------
const std::string USER_EPSG="USER_OUTPUT_EPSG";
const std::string USER_X_MIN="USER_X_MIN";
const std::string USER_X_MAX="USER_X_MAX";
const std::string USER_Y_MIN="USER_Y_MIN";
const std::string USER_Y_MAX="USER_Y_MAX";
const std::string USER_RESOLUTION="USER_RESOLUTION";

/*-------------------------
 * New types
 */
typedef std::vector<std::string> strVec_t;

class AF_InputParmeterFile
{
	public:
	AF_InputParmeterFile();
	~AF_InputParmeterFile();
	void ParseByLine();
	int CheckParsedValues();

	std::string headerFileName;


	/*===========================================
	 * Get input parameter values
	 */
	std::string GetOuputFilePath();
	std::string GetInputBFdataPath();
	std::string GetResampleMethod();
	std::string GetSourceInstrument();
	std::string GetTargetInstrument();
	// MISR section --------------------
	std::string GetMISR_Resolution();
	std::vector<std::string>  GetMISR_CameraAngles();
	std::vector<std::string>  GetMISR_Radiance();
	std::string GetMISR_Shift();
	// MODIS section -------------------
	std::string GetMODIS_Resolution();
	std::vector<std::string>  GetMODIS_Bands();
	#if 1 // JK_ASTER2MODIS
	// ASTER  section ------------------
	std::string GetASTER_Resolution();
	std::vector<std::string>  GetASTER_Bands();
	#endif
	// USER_DEFINE section ---------
	std::string GetUSER_EPSG();
	std::string GetUSER_xMin();
	std::string GetUSER_xMax();
	std::string GetUSER_yMin();
	std::string GetUSER_yMax();
	std::string GetUSER_Resolution();


	/*===========================================
	 * Handle multi-value variables
	 */
	std::vector<std::string> &GetMultiVariableNames(std::string instrument);
	int BuildMultiValueVariableMap(std::string &instrument, std::map<std::string, strVec_t> &inputMultiVarsMap);
	void DBG_displayinputListMap(std::string &instrument, std::map<std::string, strVec_t> &trgInputMultiVarsMap, const std::string &mixType);


	/*===========================================
	 * Util functions
	 */
	bool CompareStrCaseInsensitive(const std::string& s1, const std::string& s2);

	protected:
	/*===========================================
	 * Function to Validate user input values
	 */
	// Common
	bool CheckInputBFdataPath(const std::string &str);
	bool IsSourceTargetInstrumentSame();
	// MODIS
	bool CheckRevise_MODISresolution(std::string &str);
	#if 1 // JK_ASTER2MODIS
	// ASTER
	bool CheckRevise_ASTERresolution(std::string &str);
	bool CheckRevise_ASTERbands(std::vector<std::string> &strVec);
	#endif

	private:
	bool didReadHeaderFile;


	/*================================
	 * Input parameter entries
	 */
	std::string inputBFfilePath;
	std::string outputFilePath;
	std::string resampleMethod;
	std::string sourceInstrument;
	std::string targetInstrument;
	// MISR section ----------------
	std::string misr_Resolution;
	std::vector<std::string> misr_CameraAngles;
	std::vector<std::string> misr_Radiances;
	std::string misr_Shift;
	// MODIS section  --------------
	std::string modis_Resolution;
	std::vector<std::string> modis_Bands;
	#if 1 // JK_ASTER2MODIS
	// ASTER section ------------------
	std::string aster_Resolution;
	std::vector<std::string> aster_Bands;
	#endif
	// USER_DEFINE section -------
	std::string user_EPSG;
	std::string user_xMin;
	std::string user_xMax;
	std::string user_yMin;
	std::string user_yMax;
	std::string user_Resolution;


	/*=======================================
	 * Store multi-value variables names
	 */
	std::vector<std::string> modis_MultiVars;
	std::vector<std::string> misr_MultiVars;
	#if 1 // JK_ASTER2MODIS
	std::vector<std::string> aster_MultiVars;
	#endif
};
