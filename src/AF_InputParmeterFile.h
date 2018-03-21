
#include <vector>
#include <map>
#include <string>

#include "AF_debug.h"

/*---------------------------------
 * Input parameter entry strings
 */
const std::string INPUT_FILE_PATH="INPUT_FILE_PATH";
const std::string OUTPUT_FILE_PATH="OUTPUT_FILE_PATH";
const std::string RESAMPLE_METHOD="RESAMPLE_METHOD";
const std::string SOURCE_INSTRUMENT="SOURCE_INSTRUMENT";
const std::string TARGET_INSTRUMENT="TARGET_INSTRUMENT";
const std::string MISR_RESOLUTION="MISR_RESOLUTION";
const std::string MISR_CAMERA_ANGLE="MISR_CAMERA_ANGLE";
const std::string MISR_RADIANCE="MISR_RADIANCE";
const std::string MISR_SHIFT="MISR_TARGET_BLOCKUNSTACK";
const std::string MODIS_RESOLUTION="MODIS_RESOLUTION";
const std::string MODIS_BANDS="MODIS_BANDS";

typedef std::vector<std::string> strVec_t;

class AF_InputParmeterFile
{
	public:
	AF_InputParmeterFile();
	~AF_InputParmeterFile();
	void ParseByLine();

	std::string headerFileName;

	/*------------------------------
	 * Get input parameter values
	 */
	std::string GetOuputFilePath();
	std::string GetInputBFdataPath();
	std::string GetResampleMethod();
	std::string GetSourceInstrument();
	std::string GetTargetInstrument();
	// MISR
	std::string GetMISR_Resolution();
	std::vector<std::string>  GetMISR_CameraAngles();
	std::vector<std::string>  GetMISR_Radiance();
	std::string GetMISR_Shift();
	// MODIS
	std::string GetMODIS_Resolution();
	std::vector<std::string>  GetMODIS_Bands();


	/*-------------------------------
	 * Handle multi-value variables
	 */
	std::vector<std::string> &GetMultiVariableNames(std::string instrument);
	int BuildMultiValueVariableMap(std::string &instrument, std::map<std::string, strVec_t> &inputMultiVarsMap);
	void DBG_displayinputListMap(std::string &instrument, std::map<std::string, strVec_t> &trgInputMultiVarsMap, const std::string &mixType);


	/*-------------------------------
	 * Util functions
	 */
	bool CompareStrCaseInsensitive(const std::string& s1, const std::string& s2);

	protected:
	/*------------------------------
	 * Check input parameter values
	 */
	bool CheckInputBFdataPath(const std::string &str);
	bool IsSourceTargetInstrumentSame();
	bool CheckMODIS_Resolution(std::string &str);

	private:
	bool didReadHeaderFile;

	/*-------------------------
	 * Input parameter entries
	 */
	std::string inputBFfilePath;
	std::string outputFilePath;
	std::string resampleMethod;
	std::string sourceInstrument;
	std::string targetInstrument;
	// MISR
	std::string misr_Resolution;
	std::vector<std::string> misr_CameraAngles;
	std::vector<std::string> misr_Radiances;
	std::string misr_Shift;
	// MODIS
	std::string modis_Resolution;
	std::vector<std::string> modis_Bands;

	/*-----------------------------------
	 * Store multi-value variables names
	 */
	std::vector<std::string> modis_MultiVars;
	std::vector<std::string> misr_MultiVars;
};
