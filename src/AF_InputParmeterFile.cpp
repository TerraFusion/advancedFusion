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


/*=================================================================
 * AF_InputParmeterFile Constructor
 */
AF_InputParmeterFile::AF_InputParmeterFile()
{
	#if DEBUG_TOOL_PARSER
	std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Constructor AF_InputParmeterFile()\n";
	#endif
	didReadHeaderFile = false;

	/*------------------------------
	 * init multi-value variables
	 */
	// MODIS
	modis_MultiVars.push_back(MODIS_BANDS);
	// MISR
	misr_MultiVars.push_back(MISR_CAMERA_ANGLE);
	misr_MultiVars.push_back(MISR_RADIANCE);
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
	bool isValidInput = true;

	while(headerFile.good() && isValidInput)
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


		// parse single exact token without '\n', '\r' or space. 
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
			isValidInput = CheckInputBFdataPath(inputBFfilePath);
			continue;
		}


		// parse single exact token without '\n', '\r' or space.
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


		// parse single exact token without '\n', '\r' or space.
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


		// parse single exact token without '\n', '\r' or space.
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


		// parse single exact token without '\n', '\r' or space.
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


		// parse multiple
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


		// parse multiple
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


		// parse single exact token without '\n', '\r' or space.
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


		// parse single exact token without '\n', '\r' or space.
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

			isValidInput = CheckMODIS_Resolution(modis_Resolution);
			#if DEBUG_TOOL_PARSER
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> " <<  MODIS_RESOLUTION << ": " << modis_Resolution << std::endl;
			#endif
			continue;
		}


		// parse multiple
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
	} // end of while

	if (isValidInput == false) {
		std::cerr << "Error: invlid input detected.\n";
		exit(1);
	}

	if (IsSourceTargetInstrumentSame() == true) {
		exit(1);
	}
}



/*===================================================
 *	Functions to check input values
 */

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
		std::cout << "Input data file '" << filePath << "' doesn't exist.\n"; 
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
 *  - str [IN] : resolution string
 *
 * Return:
 *  - true - success
 * -  false - fail
 */
bool AF_InputParmeterFile::CheckMODIS_Resolution(std::string &str)
{
	bool ret = true;
	#if DEBUG_TOOL_PARSER
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
		std::cerr	<< "Invalid MODIS resolution '" << str << "'.\n"
					<< "Only allowed 1KM , 500M or 250M\n"
					<< std::endl;
		ret = false;
	}

	return ret;
}


//=============================================================
// Functions to get input values from the input parameter file
//
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

std::string AF_InputParmeterFile::GetTargetInstrument()
{
	return targetInstrument;
}

std::string AF_InputParmeterFile::GetMODIS_Resolution()
{
	return modis_Resolution;
}

std::vector<std::string>  AF_InputParmeterFile::GetMODIS_Bands()
{
	return modis_Bands;
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
    if(instrument == "MODIS") {
        return modis_MultiVars;
    }
    else if(instrument == "MISR") {
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
	if (instrument == "MODIS") {
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
	else if (instrument == "MISR") {
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

	if (instrument == "MODIS") {
    	#if DEBUG_TOOL_PARSER
		if (targetInstrument == "MODIS")
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Target is MODIS\n";
		if (sourceInstrument == "MODIS")
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Source is MODIS\n";
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> modis_MultiVars.size(): " << modis_MultiVars.size() << "\n";
		std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> modis_Bands.size(): " << modis_Bands.size() << "\n";
		#endif

		if (modis_MultiVars.size() != 1) {
			std::cout << __FUNCTION__ << ":" << __LINE__ <<  "> Error building input list with MODIS. There must be only one multi-value variable.\n";
			return -1;
		}


		// Note - This can be common function among instruments
		/*----------------------------
		 * Check if ALL is in bands
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

		/*------------------------------------
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
	else if (instrument == "MISR") {
    	#if DEBUG_TOOL_PARSER
		if (targetInstrument == "MISR")
			std::cout << "DBG_PARSER " << __FUNCTION__ << ":" << __LINE__ << "> Target is MISR\n";
		if (sourceInstrument == "MISR")
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

	return 0;
}
