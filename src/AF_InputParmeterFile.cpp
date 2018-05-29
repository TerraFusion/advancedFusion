/*********************************************************************
 * DEVELOPER: 
 *	- Jonathan Kim (jkm@illinois.edu)
 * 
 * DESCRIPTION:
 *	This is input parameter handling functions for advance fusion tool.
 *
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

#include <stdlib.h>
#include <stdio.h>

#include "AF_InputParmeterFile.h"
#include "AF_debug.h"
#include "gdalio.h"


/*=================================================================
 * AF_InputParmeterFile Constructor
 */
AF_InputParmeterFile::AF_InputParmeterFile()
{
	#if DEBUG_TOOL_PARSER
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Constructor AF_InputParmeterFile()\n";
	#endif
	didReadHeaderFile = false;
	misr_Shift = "ON"; // if not specified, but only effective when MISR is target

	/*------------------------------
	 * init multi-value variables
	 */
	// MODIS
	modis_MultiVars.push_back(MODIS_BANDS);
	// MISR
	misr_MultiVars.push_back(MISR_CAMERA_ANGLE);
	misr_MultiVars.push_back(MISR_RADIANCE);
	// ASTER
	aster_MultiVars.push_back(ASTER_BANDS);
}

/*=================================================================
 * AF_InputParmeterFile Destructor
 */
AF_InputParmeterFile::~AF_InputParmeterFile()
{
	#if DEBUG_TOOL_PARSER
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << " > Destructor AF_InputParmeterFile()\n";
	#endif
}

/* ################################################################
 *
 *  Util functions
 *
 * ################################################################*/
inline bool caseInsenstiveCharCompareN(char a, char b) {
	return(toupper(a) == toupper(b));
}
bool AF_InputParmeterFile::CompareStrCaseInsensitive(const std::string& s1, const std::string& s2) {
	return((s1.size() == s2.size()) &&
			equal(s1.begin(), s1.end(), s2.begin(), caseInsenstiveCharCompareN));
}


/*=================================================================
 * Parse input parameters from a input file.
 *
 * What next:
 *  Use Get...  member functions to retrive input parameters.
 *
 */
void AF_InputParmeterFile::ParseByLine()
{
	if(didReadHeaderFile)
		return;
	didReadHeaderFile = true;

	std::ifstream headerFile(headerFileName.c_str());
	if(headerFile.fail())
	{
		std::cerr << "Error: failed to read header file!\n";
	}

	std::string line;
	size_t found;
	int pos;

	while(headerFile.good())
	{
		std::getline(headerFile, line);
		#if DEBUG_TOOL_PARSER
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> line: " <<  line << std::endl;
		#endif

		// skip comment line
		if(line[0] == '#')
			continue;
		// skip empty line
		if(line.find_first_not_of(' ') == std::string::npos)
			continue;


		/*--------------------------- 
		 * INPUT_FILE_PATH
		 * parse single exact token without '\n', '\r' or space. 
		 */
		found = line.find(INPUT_FILE_PATH.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(INPUT_FILE_PATH.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact token
				inputBFfilePath = token;
			}
			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> " << INPUT_FILE_PATH << ": "  << inputBFfilePath << std::endl;
			#endif
			continue;
		}


		/*--------------------------- 
		 * OUTPUT_FILE_PATH
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(OUTPUT_FILE_PATH.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(OUTPUT_FILE_PATH.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				outputFilePath = token;
			}
			#if DEBUG_TOOL_PARSER
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> " <<  OUTPUT_FILE_PATH << ": " << outputFilePath << std::endl;
			#endif
			continue;
		}


		/*--------------------------- 
		 * RESAMPLE_METHOD
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(RESAMPLE_METHOD.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(RESAMPLE_METHOD.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				resampleMethod = token;
			}
			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> " <<  RESAMPLE_METHOD << ": " << resampleMethod << std::endl;
			#endif
			continue;
		}


		/*--------------------------- 
		 * SOURCE_INSTRUMENT
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(SOURCE_INSTRUMENT.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(SOURCE_INSTRUMENT.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				sourceInstrument = token;
			}
			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> " <<  SOURCE_INSTRUMENT << ": " << sourceInstrument << std::endl;
			#endif
			continue;
		}


		/*--------------------------- 
		 * TARGET_INSTRUMENT
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(TARGET_INSTRUMENT.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(TARGET_INSTRUMENT.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				targetInstrument = token;
			}
			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> " <<  TARGET_INSTRUMENT << ": " << targetInstrument << std::endl;
			#endif
			continue;
		}



		/*======================================================================
		 * MISR section
		 */

		/*--------------------------- 
		 * MISR_RESOLUTION
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(MISR_RESOLUTION.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(MISR_RESOLUTION.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				misr_Resolution = token;
			}
			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> " <<  MISR_RESOLUTION << ": " << misr_Resolution << std::endl;
			#endif
			continue;
		}


		/*--------------------------- 
		 * MISR_CAMERA_ANGLE
		 * parse multiple
		 */
		found = line.find(MISR_CAMERA_ANGLE.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(MISR_CAMERA_ANGLE.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {
				misr_CameraAngles.push_back(token);
			}
			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Num of cameras:" << misr_CameraAngles.size() << std::endl;
			for(int i = 0; i < misr_CameraAngles.size(); i++) {
				std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> misr_CameraAngles[" << i << "]:" << misr_CameraAngles[i] << std::endl;
			}
			#endif
			continue;
		}


		/*--------------------------- 
		 * MISR_RADIANCE
		 * parse multiple
		 */
		found = line.find(MISR_RADIANCE.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(MISR_RADIANCE.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				misr_Radiances.push_back(token);
			}
			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Num of radiances:" << misr_Radiances.size() << std::endl;
			for(int i = 0; i < misr_Radiances.size(); i++) {
				std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> misr_Radiances[" << i << "]:" << misr_Radiances[i] << std::endl;
			}
			#endif
			continue;
		}

		/*--------------------------- 
		 * MISR_TARGET_BLOCKUNSTACK
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(MISR_SHIFT.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(MISR_SHIFT.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				misr_Shift = token;
			}
			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> " <<  MISR_SHIFT << ": " << misr_Shift << std::endl;
			#endif
			continue;
		}


		/*======================================================================
		 * MODIS section
		 */

		/*-------------------- 
		 * MODIS_RESOLUTION  
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(MODIS_RESOLUTION.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(MODIS_RESOLUTION.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				modis_Resolution = token;
			}

			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> " <<  MODIS_RESOLUTION << ": " << modis_Resolution << std::endl;
			#endif
			continue;
		}


		/*-------------------- 
		 * MODIS_BANDS  
		 * parse multiple
		 */
		found = line.find(MODIS_BANDS.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(MODIS_BANDS.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {
				modis_Bands.push_back(token.c_str());
			}
			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Num of modis bands:" << modis_Bands.size() << std::endl;
			for(int i = 0; i < modis_Bands.size(); i++) {
				std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> modis_Bands[" << i << "]:" << modis_Bands[i] << std::endl;
			}
			#endif
			continue;
		}


		/*======================================================================
		 * ASTER section
		 */

		/*-------------------- 
		 * ASTER_RESOLUTION  
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(ASTER_RESOLUTION.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(ASTER_RESOLUTION.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				aster_Resolution = token;
			}

			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> " <<  ASTER_RESOLUTION << ": " << aster_Resolution << std::endl;
			#endif
			continue;
		}

		/*-------------------- 
		 * ASTER_BANDS  
		 * parse multiple
		 */
		found = line.find(ASTER_BANDS.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(ASTER_BANDS.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {
				aster_Bands.push_back(token.c_str());
			}
			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Num of aster bands:" << aster_Bands.size() << std::endl;
			for(int i = 0; i < aster_Bands.size(); i++) {
				std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> aster_Bands[" << i << "]:" << aster_Bands[i] << std::endl;
			}
			#endif
			continue;
		}


		/*======================================================================
		 * USER_DEFINE section
		 */

		/*-------------------- 
		 * USER_OUTPUT_EPSG
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(USER_EPSG.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(USER_EPSG.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				user_EPSG = token;
			}

			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> USER_EPSG : " << user_EPSG << std::endl;
			#endif
			continue;
		}

		/*-------------- 
		 * USER_X_MIN
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(USER_X_MIN.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(USER_X_MIN.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				user_xMin = token;
			}

			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> USER_X_MIN : " << user_xMin << std::endl;
			#endif
			continue;
		}

		/*------------- 
		 * USER_X_MAX
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(USER_X_MAX.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(USER_X_MAX.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				user_xMax = token;
			}

			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> USER_X_MAX : " << user_xMax << std::endl;
			#endif
			continue;
		}

		/*-------------- 
		 * USER_Y_MIN
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(USER_Y_MIN.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(USER_Y_MIN.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				user_yMin = token;
			}

			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> USER_Y_MIN : " << user_yMin << std::endl;
			#endif
			continue;
		}

		/*------------- 
		 * USER_Y_MAX
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(USER_Y_MAX.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(USER_Y_MAX.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				user_yMax = token;
			}

			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> USER_Y_MAX : " << user_yMax << std::endl;
			#endif
			continue;
		}

		/*------------- 
		 * USER_RESOLUTION
		 * parse single exact token without '\n', '\r' or space.
		 */
		found = line.find(USER_RESOLUTION.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(USER_RESOLUTION.c_str()));
			while(line[0] == ' ' || line[0] == ':')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				user_Resolution = token;
			}

			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> USER_RESOLUTION : " << user_Resolution << std::endl;
			#endif
			continue;
		}

	} // end of while
}


/* ##########################################################
 *	Functions to check input values
 */

int AF_InputParmeterFile::CheckParsedValues()
{
	bool isValidInput;

	
	/*=================================================
	 * Common section
	 */
	isValidInput = CheckInputBFdataPath(inputBFfilePath);
	if (isValidInput == false) {
		return -1; // failed
	}

	if (IsSourceTargetInstrumentSame() == true) {
		std::cerr << "Error: Source and target instrument must be different.\n";
		return -1; // failed
	}


	/*=================================================
	 * MODIS section
	 */
	if (sourceInstrument == MODIS_STR || targetInstrument == MODIS_STR) {
		isValidInput = CheckRevise_MODISresolution(modis_Resolution);
		if (isValidInput == false) {
			return -1; // failed
		}
	}


	/*=================================================
	 * ASTER section
	 */
	if (sourceInstrument == ASTER_STR || targetInstrument == ASTER_STR) {
		isValidInput = CheckRevise_ASTERresolution(aster_Resolution);
		if (isValidInput == false) {
			return -1; // failed
		}
	
		isValidInput = CheckRevise_ASTERbands(aster_Bands);
		if (isValidInput == false) {
			return -1; // failed
		}
	}

	
	return 0; // succeed
}


/*=================================================================
 * Check if BF input file path is vaild
 *
 * Parameter:
 *  - filePath [IN] : BF input file path
 *
 * Return:
 *  - true -  valid
 * -  false - not valid
 */

bool AF_InputParmeterFile::CheckInputBFdataPath(const std::string &filePath)
{
	bool ret = true;
	#if DEBUG_TOOL_PARSER
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> BF file path: " << filePath <<	".\n"; 
	#endif
	// check if 
	std::ifstream file(filePath.c_str());
	ret = (file.good() ? true : false);
	if(ret == false) {
		std::cerr << "Error: Input data file '" << filePath << "' doesn't exist.\n"; 
	}
	return ret;
}

/*=================================================================
 * Check if source instrument is same as target instrument
 *
 * Parameter:
 *  - str [IN] : resolution string
 *
 * Return:
 *  - true -  same
 * -  false - not same
 */
bool AF_InputParmeterFile::IsSourceTargetInstrumentSame()
{
	bool ret = false;
	#if DEBUG_TOOL_PARSER
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> sourceInstrument: " << sourceInstrument << ", targetInstrument: " << targetInstrument <<   ".\n";
	#endif

	if (sourceInstrument == targetInstrument) {
		std::cerr	<< "Error: Source instrument must be different from target instrument.\n";
		ret = true;
	}

	return ret;
}


/*=================================================================
 * Check Modis resolution input
 *
 * Parameter:
 *  - str [IN/OUT] : resolution string.
 *					 Check and convert it for internal notation.
 *
 * Return:
 *  - true - success
 * -  false - fail
 */
bool AF_InputParmeterFile::CheckRevise_MODISresolution(std::string &str)
{
	bool ret = true;
	#if DEBUG_TOOL_PARSER
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Before update. \n";
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> reolution: " << str <<	".\n";
	#endif

	if (CompareStrCaseInsensitive(str, "1KM")) {
		str = "_1KM";
	}
	else if  (CompareStrCaseInsensitive(str, "500M")) {
		str = "_500m";
	}
	else if (CompareStrCaseInsensitive(str, "250M")) {
		str = "_250m";
	}
	else {
		std::cerr	<< "Error: Invalid MODIS resolution '" << str << "'.\n"
					<< "Only allowed 1KM , 500M or 250M\n"
					<< std::endl;
		ret = false;
	}

	#if DEBUG_TOOL_PARSER
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> After update. \n";
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> reolution: " << str <<	".\n";
	#endif

	return ret;
}


/*=================================================================
 * Check Aster resolution input
 *
 * Parameter:
 *  - str [IN/OUT] : resolution string.
 *					 Check and convert it for internal notation.
 *
 * Return:
 *  - true - success
 * -  false - fail
 */
bool AF_InputParmeterFile::CheckRevise_ASTERresolution(std::string &str)
{
	bool ret = true;
	#if DEBUG_TOOL_PARSER
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Before update. \n";
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> reolution: " << str <<	".\n";
	#endif

	if (CompareStrCaseInsensitive(str, "90M")) {
		str = "TIR";
	}
	else if  (CompareStrCaseInsensitive(str, "30M")) {
		str = "SWIR";
	}
	else if (CompareStrCaseInsensitive(str, "15M")) {
		str = "VNIR";
	}
	else {
		std::cerr	<< "Error: Invalid ASTER resolution '" << str << "'.\n"
					<< "Only allowed 90M , 30M or 15M\n"
					<< std::endl;
		ret = false;
	}
	#if DEBUG_TOOL_PARSER
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> After update. \n";
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> reolution: " << str <<	".\n";
	#endif

	return ret;
}

/*=================================================================
 * Check Aster bands input
 *
 * Parameter:
 *  - str [IN/OUT] : band string.
 *					 Check and convert it for internal notation.
 *
 * Return:
 *  - true - success
 * -  false - fail
 */
bool AF_InputParmeterFile::CheckRevise_ASTERbands(std::vector<std::string> &strVec)
{
	bool ret = true;
	#if DEBUG_TOOL_PARSER
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Before update" << std::endl;
	for(int i = 0; i < strVec.size(); i++) {
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> asterBands[" << i << "]:" << strVec[i] << std::endl;
    }
	#endif


	#if 0 // TODO_LATER - allow range by Resolution 
	if (CompareStrCaseInsensitive(aster_Resolution, "TIR" /*90M*/)) {
		for(int i = 0; i < strVec.size(); i++) {
			if (CompareStrCaseInsensitive(strVec[i], "10") ||
	        	CompareStrCaseInsensitive(strVec[i], "11") ||
	        	CompareStrCaseInsensitive(strVec[i], "12") ||
	        	CompareStrCaseInsensitive(strVec[i], "13") ||
	        	CompareStrCaseInsensitive(strVec[i], "14")) {
				strVec[i] = "ImageData" + strVec[i];
			}
			else {
				std::cerr	<< "Error: Invalid ASTER band '" << strVec[i] << "'.\n"
							<< "Only allowed 10 to 14 for 90M.\n";
				ret = false;
			}
		}
	}
	else if (CompareStrCaseInsensitive(aster_Resolution, "SWIR" /*30M*/)) {
		// only allow 4 ~ 9
	}
	else if (CompareStrCaseInsensitive(aster_Resolution, "VNIR" /*15M*/)) {
		// only allow 1 ~ 3
	}
	#else
	for(int i = 0; i < strVec.size(); i++) {
		if (CompareStrCaseInsensitive(strVec[i], "1") ||
			CompareStrCaseInsensitive(strVec[i], "2") ||
			CompareStrCaseInsensitive(strVec[i], "4") ||
			CompareStrCaseInsensitive(strVec[i], "5") ||
			CompareStrCaseInsensitive(strVec[i], "6") ||
			CompareStrCaseInsensitive(strVec[i], "7") ||
			CompareStrCaseInsensitive(strVec[i], "8") ||
			CompareStrCaseInsensitive(strVec[i], "9") ||
			CompareStrCaseInsensitive(strVec[i], "10") ||
			CompareStrCaseInsensitive(strVec[i], "11") ||
			CompareStrCaseInsensitive(strVec[i], "12") ||
			CompareStrCaseInsensitive(strVec[i], "13") ||
			CompareStrCaseInsensitive(strVec[i], "14") ) {
			strVec[i] = "ImageData" + strVec[i];
		}
		else if  (CompareStrCaseInsensitive(strVec[i], "3")) {
			strVec[i] = "ImageData3N";
		}
		else {
			std::cerr	<< "Error: Invalid ASTER band '" << strVec[i] << "'.\n"
						<< "Only allowed 1 to 14 integer value.\n"
						<< std::endl;
			ret = false;
		}
	}
	#endif

	#if DEBUG_TOOL_PARSER
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> After update" << std::endl;
	for(int i = 0; i < strVec.size(); i++) {
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> asterBands[" << i << "]:" << strVec[i] << std::endl;
    }
	#endif

	return ret;
}

/* ######################################################################
 * Functions to get input based parameter of internal functions
 */

/*=================================================================
 * Check Aster bands input
 *
 * Parameter:
 *  - instrument : instrumane name string.
 *
 * Return:
 *  - Success :  > 0
 *  - Fail : 0
 */
double AF_InputParmeterFile::GetMaxRadiusForNNeighborFunc(std::string instrument)
{
	double maxRadius = 0.0;

	/*---------------------
	 * MODIS section
	 */
	if(instrument == MODIS_STR) {
		if(modis_Resolution == "_1KM") {
			maxRadius = 5040.0;
		}
		else if(modis_Resolution == "_500m") {
			maxRadius = 2520.0;
		}
		else if(modis_Resolution == "_250m") {
			maxRadius = 1260.0;
		}
		#if DEBUG_TOOL_PARSER
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> modis_Resolution: " << modis_Resolution << ", maxRadius: " << maxRadius <<  std::endl;
		#endif
        return maxRadius;
    }
	/*---------------------
	 * MISR section
	 */
    else if(instrument == MISR_STR) {
		if(misr_Resolution == "H") {
			maxRadius = 302.0;
		}
		else if(misr_Resolution == "L") {
			maxRadius = 1155.0;
		}
		#if DEBUG_TOOL_PARSER
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> misr_Resolution: " << misr_Resolution << ", maxRadius: " << maxRadius <<  std::endl;
		#endif
        return maxRadius;
    }
	/*---------------------
	 * ASTER section
	 */
    else if(instrument == ASTER_STR) {
		if(aster_Resolution == "TIR") { // 90M
			maxRadius = 95.0;
		}
		else if(aster_Resolution == "SWIR") { // 30M
			maxRadius = 32.0;
		}
		else if(aster_Resolution == "VNIR") { // 15M
			maxRadius = 17.0;
		}
		#if DEBUG_TOOL_PARSER
		std::cout << "DBG_PARSER " <<  __FUNCTION__ << ":" << __LINE__ << "> aster_Resolution: " << aster_Resolution << ", maxRadius: " << maxRadius <<  std::endl;
		#endif
        return maxRadius;
	}
	/*-------------------------
	 * USER_DEFINE section
	 */
	else if(instrument == USERGRID_STR) {
		int userEPSG =  GetUSER_EPSG();
		double userResolution = GetUSER_Resolution();
		maxRadius =  getMaxRadiusOfUserdefine(userEPSG, userResolution);
		#if DEBUG_TOOL_PARSER
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> userEPSG: " << userEPSG << ", userResolution: " << userResolution << ", maxRadius: " << maxRadius <<  std::endl;
		#endif
		return maxRadius;
    }
    else {
        std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: incorrect instrument name '"<< instrument << "' is specified." << "\n";
        return 0;
    }
}


/* #########################################################################
 *  Functions to get input values from the input parameter file
 */
std::string AF_InputParmeterFile::GetOuputFilePath()
{
	return outputFilePath;
}

std::string AF_InputParmeterFile::GetInputBFdataPath()
{
	return inputBFfilePath;
}

std::string AF_InputParmeterFile::GetResampleMethod()
{
	return resampleMethod;
}

std::string AF_InputParmeterFile::GetSourceInstrument()
{
	return sourceInstrument;
}

std::string AF_InputParmeterFile::GetTargetInstrument()
{
	return targetInstrument;
}

/*---------------------
 * MISR section
 */
std::string AF_InputParmeterFile::GetMISR_Resolution()
{
	return misr_Resolution;
}

std::vector<std::string>  AF_InputParmeterFile::GetMISR_CameraAngles()
{
	return misr_CameraAngles;
}

std::vector<std::string> AF_InputParmeterFile::GetMISR_Radiance()
{
	return misr_Radiances;
}

std::string AF_InputParmeterFile::GetMISR_Shift()
{
	return misr_Shift;
}


/*---------------------
 * MODIS section
 */
std::string AF_InputParmeterFile::GetMODIS_Resolution()
{
	return modis_Resolution;
}

std::vector<std::string>  AF_InputParmeterFile::GetMODIS_Bands()
{
	return modis_Bands;
}


/*---------------------
 * ASTER section
 */
std::string AF_InputParmeterFile::GetASTER_Resolution()
{
	return aster_Resolution;
}

std::vector<std::string>  AF_InputParmeterFile::GetASTER_Bands()
{
	return aster_Bands;
}


/*---------------------
 * USER section
 */
int AF_InputParmeterFile::GetUSER_EPSG()
{
	// convert string to int
	int retValue;
	std::stringstream ss(user_EPSG);
    ss >> retValue;
	return retValue;
}

double AF_InputParmeterFile::GetUSER_xMin()
{
	// convert string to double
	double retValue;
	std::stringstream ss(user_xMin);
    ss >> retValue;
	return retValue;
}

double AF_InputParmeterFile::GetUSER_xMax()
{
	// convert string to double
	double retValue;
	std::stringstream ss(user_xMax);
    ss >> retValue;
	return retValue;
}

double AF_InputParmeterFile::GetUSER_yMin()
{
	// convert string to double
	double retValue;
	std::stringstream ss(user_yMin);
    ss >> retValue;
	return retValue;
}

double AF_InputParmeterFile::GetUSER_yMax()
{
	// convert string to double
	double retValue;
	std::stringstream ss(user_yMax);
    ss >> retValue;
	return retValue;
}

double AF_InputParmeterFile::GetUSER_Resolution()
{
	// convert string to int
	double retValue;
	std::stringstream ss(user_Resolution);
    ss >> retValue;
	return retValue;
}


/* #####################################
 *
 * Handling multi-value variables
 *
 * #####################################*/

/*=========================================================
 * Get multi-value variable names for the given instrument
 *
 * Parameters:
 *  - instrument [IN] : name string of an instrument.
 * Return:
 *  - multi-value variable names as list of strings
 */
std::vector<std::string> & AF_InputParmeterFile::GetMultiVariableNames(std::string instrument)
{
    if(instrument == MODIS_STR) {
        return modis_MultiVars;
    }
    else if(instrument == MISR_STR) {
        return misr_MultiVars;
    }
    else {
        std::cerr << __FUNCTION__ << ":" << __LINE__ << "> Error: incorrect instrument name '"<< instrument << "' is specified." << "\n";
        exit(1);
    }
}


/*=============================================================
 * Debugging purpose. Show multi-value variable map contents
 * for the given instrument.
 *
 * Parameters:
 *  - instrument [IN]: name string of an instrument.
 *  - trgInputMultiVarsMap [IN]: multi-value variable map
 *  - mixType [IN]: "COMBINATION" or "PAIR"  Only support "COMBINATION" for now
 */
void AF_InputParmeterFile::DBG_displayinputListMap(std::string &instrument, std::map<std::string, strVec_t> &trgInputMultiVarsMap, const std::string &mixType)
{
    strVec_t multiVarNames;
    //--------------------------------------------------------
    // Use this for Single multi-value Variable case. MODIS

    std::cout << "JKDBG> trgInputMultiVarsMap.size(): " << trgInputMultiVarsMap.size() << "\n";
    // display strBuf with array index
	if (instrument == MODIS_STR) {
        multiVarNames = modis_MultiVars;

        if (trgInputMultiVarsMap.size() != 1) {
            std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building target input list with MODIS. There must be only one multi-value variable.\n";
        }

        std::cout << "Display trgInputMultiVarsMap with array index\n";
        for(int i = 0; i < trgInputMultiVarsMap.size(); ++i)
        {
            for(int j = 0; j < trgInputMultiVarsMap[multiVarNames[i]].size(); ++j)
                std::cout << trgInputMultiVarsMap[multiVarNames[i]][j]  << ", ";
            std::cout << std::endl;
        }
    }
    //--------------------------------------------------------
    // Use these for two multi-value Variable case. MISR and ASTER
	else if (instrument == MISR_STR) {
        multiVarNames = misr_MultiVars;

        if (trgInputMultiVarsMap.size() != 2) {
            std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building target input list with MISR. There must be only two multi-value variables.\n";
        }

        std::string misrCamera;
        std::string misrRadiance;
        if(mixType == "PAIR") {
            std::cout << "\nJKDBG> mixType == PAIR\n";
            // TODO - Add a function to determine minNumVals among mutiple variables. Assume it's 2 for test purpose
            int numCameras = trgInputMultiVarsMap[multiVarNames[0]].size();
            int numRadiances = trgInputMultiVarsMap[multiVarNames[1]].size();
            std::cout << "JKDBG> var0 num of cameras: " << numCameras << "\n";
            std::cout << "JKDBG> var1 num of radiances: " << numRadiances << "\n";
            int minNumVals = (numCameras < numRadiances) ? numCameras : numRadiances;
            std::cout << "JKDBG> minNumVals: " << minNumVals << "\n";
            for(int i=0; i<minNumVals ;i++) {
                //std::cout << "JKDBG trgInputMultiVarsMap[multiVarNames[0]][" << i << "]" <<  trgInputMultiVarsMap[multiVarNames[0]][i] << "\n";
                //std::cout << "JKDBG trgInputMultiVarsMap[multiVarNames[1]][" << i << "]" <<  trgInputMultiVarsMap[multiVarNames[1]][i] << "\n";
                // misr camera
                misrCamera = trgInputMultiVarsMap[multiVarNames[0]][i];
                // misr rad
                misrRadiance = trgInputMultiVarsMap[multiVarNames[1]][i];
                std::cout << misrCamera << " : " << misrRadiance <<  "\n";
            }
        }
        else if(mixType == "COMBINATION") {
            std::cout << "\nJKDBG> mixType == COMBINATION\n";
            for(int i=0; i<trgInputMultiVarsMap.size()-1;i++) {
                for(int j=i; j <trgInputMultiVarsMap[multiVarNames[i]].size(); j++) {
                    misrCamera = trgInputMultiVarsMap[multiVarNames[i]][j];
                    for (int k=0; k<trgInputMultiVarsMap[multiVarNames[i+1]].size(); k++) {
                        misrRadiance = trgInputMultiVarsMap[multiVarNames[i+1]][k];
                        std::cout << misrCamera << ":" << misrRadiance << "\n";
                    }
                }
            }
        }
    }
}

/*=============================================================================
 * Build multi-value variables map which is generic container among instruments
 *
 * Parameters:
 *  - instrument [IN] : name string of an instrument.
 *  - trgInputMultiVarsMap [OUT] : multi-value variable map to build
 *
 * Return:
 *  - 0 : Succeed
 *  - -1: Failed
 */
int AF_InputParmeterFile::BuildMultiValueVariableMap(std::string &instrument, std::map<std::string, strVec_t> &inputMultiVarsMap)
{
	std::vector<std::string> strVecTmp;

	/*===========================================================
	 * MODIS section
	 */
	if (instrument == MODIS_STR) {
    	#if DEBUG_TOOL_PARSER
		if (targetInstrument == MODIS_STR)
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Target is MODIS\n";
		if (sourceInstrument == MODIS_STR)
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Source is MODIS\n";
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> modis_MultiVars.size(): " << modis_MultiVars.size() << "\n";
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> modis_Bands.size(): " << modis_Bands.size() << "\n";
		#endif

		if (modis_MultiVars.size() != 1) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building input list with MODIS. There must be only one multi-value variable.\n";
			return -1;
		}


		// Note - can this be common among instruments?
		/*----------------------------
		 * Check if 'ALL' is in bands
		 */
		bool isAllbands = false;
		for(int i=0; i < modis_Bands.size(); i++) {
			if (CompareStrCaseInsensitive(modis_Bands[i], "ALL")) {
				isAllbands = true;
				#if DEBUG_TOOL_PARSER
				std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> ALL for Modis bands.\n";
				#endif
				break;
			}
        }
		/*
		 * if ALL is in bands, build all list
		 */
		std::vector<std::string> modisBandsUpdated;
		if (isAllbands) {
			if(modis_Resolution == "_1KM") {
				// 38 bands total
				modisBandsUpdated = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13L", "13H", "14L", "14H", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36"};
			}
			else if(modis_Resolution == "_500m") {
				// 8 bands total
				modisBandsUpdated = {"1", "2", "3","4", "5", "6", "7"};
			}
			else if(modis_Resolution == "_250m") {
				// 2 bands total
				modisBandsUpdated = {"1","2"};
			}
		}
		else {
			modisBandsUpdated = modis_Bands;
		}

		#if DEBUG_TOOL_PARSER
		std::cout << __FUNCTION__ << ":" << __LINE__ << "> modisBandsUpdated: ";
		for (int i=0; i< modisBandsUpdated.size(); i++) {
			std::cout << modisBandsUpdated[i] << " ";
		}
		std::cout << "\n";
		#endif


		/*----------------------------
		 * Note: used modisBandsUpdated not modis_Bands directly to handle
		 */
		if(modis_MultiVars[0] == MODIS_BANDS) {
			strVecTmp.clear();
			for(int j=0; j < modisBandsUpdated.size(); j++) {
				strVecTmp.push_back(modisBandsUpdated[j]);
			}
			inputMultiVarsMap[modis_MultiVars[0]] = strVecTmp;
			/*---------------------
			Ex:
			inputMultiVarsMap["MODIS_BANDS"] = "{"8","9","10"}
			inputMultiVarsMap.size() == 1
			----------------------*/
		}
		else {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building input list with MODIS.\n";
		}
	}
	/*===========================================================
	 * MISR section
	 */
	else if (instrument == MISR_STR) {
    	#if DEBUG_TOOL_PARSER
		if (targetInstrument == MISR_STR)
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Target is MISR\n";
		if (sourceInstrument == MISR_STR)
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Source is MISR\n";
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> misr_MultiVars.size(): " << misr_MultiVars.size() << "\n";
		#endif
		if (misr_MultiVars.size() != 2) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building input list with MISR. There must be only two multi-value variables.\n";
		}

		for(int i=0; i<misr_MultiVars.size(); i++) {
			strVecTmp.clear();
			if(misr_MultiVars[i] == MISR_CAMERA_ANGLE) {
				for(int j=0; j < misr_CameraAngles.size(); j++) {
					strVecTmp.push_back(misr_CameraAngles[j]);
				}
				inputMultiVarsMap[MISR_CAMERA_ANGLE] = strVecTmp;
				/*---------------------
				Ex:
				inputMultiVarsMap["MISR_CAMERA_ANGLE"] = "{"AN","AF", ..}
				inputMultiVarsMap.size() == 1
				----------------------*/
			}
			else if(misr_MultiVars[i] == MISR_RADIANCE) {
				for(int j=0; j < misr_Radiances.size(); j++) {
					strVecTmp.push_back(misr_Radiances[j]);
				}
				inputMultiVarsMap[MISR_RADIANCE] = strVecTmp;
				/*---------------------
				Ex:
				inputMultiVarsMap["MISR_RADIANCE"] = "{"BLUE","GREEN", ..}
				inputMultiVarsMap.size() == 2
				----------------------*/
			}
			else {
				std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building target input list with MISR.\n";
			}
		}
	}
	/*===========================================================
	 * ASTER section
	 */
	else if (instrument == ASTER_STR) {
    	#if DEBUG_TOOL_PARSER
		if (targetInstrument == ASTER_STR)
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Target is ASTER\n";
		if (sourceInstrument == ASTER_STR)
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Source is ASTER\n";
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> aster_MultiVars.size(): " << aster_MultiVars.size() << "\n";
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> aster_Bands.size(): " << aster_Bands.size() << "\n";
		#endif

		if (aster_MultiVars.size() != 1) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building input list with ASTER. There must be only one multi-value variable.\n";
			return -1;
		}

		#if 0 // TODO_LATER - handle ALL cases. Placeholder
		/*----------------------------
		 * Check if 'ALL' is in bands
		 */
		bool isAllbands = false;
		for(int i=0; i < aster_Bands.size(); i++) {
			if (CompareStrCaseInsensitive(aster_Bands[i], "ALL")) {
				isAllbands = true;
				#if DEBUG_TOOL_PARSER
				std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> ALL for Aster bands.\n";
				#endif
				break;
			}
        }

		/* 
		 * if ALL is in bands, build all list
		 */
		std::vector<std::string> asterBandsUpdated;
		if (isAllbands) {
			if(aster_Resolution == "TIR") {
				// 38 bands total
				asterBandsUpdated = {"10", "11", "12", "13", "14"};
			}
			else if(aster_Resolution == "SWIR") {
				// 8 bands total
				asterBandsUpdated = {"4", "5", "6", "7", "8", "9"};
			}
			else if(aster_Resolution == "VNIR") {
				// 2 bands total
				asterBandsUpdated = {"ImageData1","2", "3N"};
			}
		}
		else {
			asterBandsUpdated = aster_Bands;
		}
		#else
		std::vector<std::string> asterBandsUpdated;
		asterBandsUpdated = aster_Bands;
		#endif

		#if DEBUG_TOOL_PARSER
		std::cout << __FUNCTION__ << ":" << __LINE__ << "> asterBandsUpdated: ";
		for (int i=0; i< asterBandsUpdated.size(); i++) {
			std::cout << asterBandsUpdated[i] << " ";
		}
		std::cout << "\n";
		#endif

		/*----------------------------
		 * Note: used asterBandsUpdated not aster_Bands directly to handle
		 */
		if(aster_MultiVars[0] == ASTER_BANDS) {
			strVecTmp.clear();
			for(int j=0; j < asterBandsUpdated.size(); j++) {
				strVecTmp.push_back(asterBandsUpdated[j]);
			}
			inputMultiVarsMap[aster_MultiVars[0]] = strVecTmp;
			/*---------------------
			Ex:
			inputMultiVarsMap["ASTER_BANDS"] = "{"8","9","10"}
			inputMultiVarsMap.size() == 1
			----------------------*/
		}
		else {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building input list with ASTER.\n";
		}
	}

	return 0;
}
