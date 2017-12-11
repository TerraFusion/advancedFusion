
#include <vector>
#include <string>

#include "AF_debug.h"

class AF_InputParmeterFile
{
	public:
	AF_InputParmeterFile();
	~AF_InputParmeterFile();
	void ParseByLine();

	//------------------------------
	// Get input parameter values
	std::string GetOuputFilePath();
	std::string GetInputBFdataPath();
	std::string GetResampleMethod();
	std::string GetSourceInstrument();
	std::string GetMISR_Resolution();
	std::vector<std::string>  GetMISR_CameraAngles();
	std::string GetMISR_Radiance();
	std::string GetTargetInstrument();
	std::string GetMODIS_Resolution();
	std::vector<int>  GetMODIS_Bands();		  

	std::string headerFileName;

	bool CompareStrCaseInsensitive(const std::string& s1, const std::string& s2);

	protected:
	//------------------------------
	// Check input parameter values
	bool CheckInputBFdataPath(const std::string &str);
	bool CheckMODIS_Resolution(std::string &str);

	private:
	bool didReadHeaderFile;

	// input parameter entries
	std::string inputBFfilePath;
	std::string outputFilePath;
	std::string resampleMethod;
	std::string sourceInstrument;
	std::string misr_Resolution;
	std::vector<std::string> misr_CameraAngles;
	std::string misr_Radiance;
	std::string targetInstrument;
	std::string modis_Resolution;
	std::vector<int> modis_Bands;
};
