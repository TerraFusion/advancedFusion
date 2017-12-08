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


AF_InputParmeterFile::AF_InputParmeterFile()
{
	#if DEBUG_TOOL
	std::cout << "DBG> Constructor AF_InputParmeterFile()\n";
	#endif
	didReadHeaderFile = false;
}

AF_InputParmeterFile::~AF_InputParmeterFile()
{
	#if DEBUG_TOOL
	std::cout << "DBG> Destructor AF_InputParmeterFile()\n";
	#endif
}

//---------------------------
// util functions
inline bool caseInsenstiveCharCompareN(char a, char b) {
	return(toupper(a) == toupper(b));
}

bool AF_InputParmeterFile::CompareStrCaseInsensitive(const std::string& s1, const std::string& s2) {
	return((s1.size() == s2.size()) &&
			equal(s1.begin(), s1.end(), s2.begin(), caseInsenstiveCharCompareN));
}


//------------------------------------------------
// parsing input parameters from a input file, and can retrive with Get member functions. 
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

	// input parameter entry strings
	const std::string INPUT_FILE_PATH="INPUT_FILE_PATH:";
	const std::string OUTPUT_FILE_PATH="OUTPUT_FILE_PATH:";
	const std::string RESAMPLE_METHOD="RESAMPLE_METHOD:";
	const std::string SOURCE_INSTRUMENT="SOURCE_INSTRUMENT:";
	const std::string MISR_RESOLUTION="MISR_RESOLUTION:";
	const std::string MISR_CAMERA_ANGLE="MISR_CAMERA_ANGLE:";
	const std::string MISR_RADIANCE="MISR_RADIANCE:";
	const std::string TARGET_INSTRUMENT="TARGET_INSTRUMENT:";
	const std::string MODIS_RESOLUTION="MODIS_RESOLUTION:";
	const std::string MODIS_BANDS="MODIS_BANDS:";

	std::string line;
	size_t found;
	int pos;
	bool isValidInput = true;

	while(headerFile.good() && isValidInput)
	{
		std::getline(headerFile, line);
		#if DEBUG_TOOL
		std::cout << "DBG ParseByLine()> line: " <<  line << std::endl;
		std::cout << "DBG> line size: " <<  line.size() << std::endl;
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
			while(line[0] == ' ')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact token
				inputBFfilePath = token;
			}
			#if DEBUG_TOOL
			std::cout << "DBG ParseByLine()> " << INPUT_FILE_PATH << ": "  << inputBFfilePath << std::endl;
			#endif
			isValidInput = CheckInputBFdataPath(inputBFfilePath);
			continue;
		}


		// parse single exact token without '\n', '\r' or space.
		found = line.find(OUTPUT_FILE_PATH.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(OUTPUT_FILE_PATH.c_str()));
			while(line[0] == ' ')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				outputFilePath = token;
			}
			#if DEBUG_TOOL
			std::cout << "DBG ParseByLine()> " <<  OUTPUT_FILE_PATH << ": " << outputFilePath << std::endl;
			#endif
			continue;
		}


		// parse single exact token without '\n', '\r' or space.
		found = line.find(RESAMPLE_METHOD.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(RESAMPLE_METHOD.c_str()));
			while(line[0] == ' ')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				resampleMethod = token;
			}
			#if DEBUG_TOOL
			std::cout << "DBG ParseByLine()> " <<  RESAMPLE_METHOD << ": " << resampleMethod << std::endl;
			#endif
			continue;
		}


		// parse single exact token without '\n', '\r' or space.
		found = line.find(SOURCE_INSTRUMENT.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(SOURCE_INSTRUMENT.c_str()));
			while(line[0] == ' ')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				sourceInstrument = token;
			}
			#if DEBUG_TOOL
			std::cout << "DBG ParseByLine()> " <<  SOURCE_INSTRUMENT << ": " << sourceInstrument << std::endl;
			#endif
			continue;
		}


		// parse single exact token without '\n', '\r' or space.
		found = line.find(MISR_RESOLUTION.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(MISR_RESOLUTION.c_str()));
			while(line[0] == ' ')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				misr_Resolution = token;
			}
			#if DEBUG_TOOL
			std::cout << "DBG ParseByLine()> " <<  MISR_RESOLUTION << ": " << misr_Resolution << std::endl;
			#endif
			continue;
		}


		// parse multiple
		found = line.find(MISR_CAMERA_ANGLE.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(MISR_CAMERA_ANGLE.c_str()));
			while(line[0] == ' ')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {
				misr_CameraAngles.push_back(token);
			}
			#if DEBUG_TOOL
			//std::cout << "DBG ParseByLine()> Num of cameras:" << misr_CameraAngles.size() << std::endl;
			for(int i = 0; i < misr_CameraAngles.size(); i++) {
				std::cout << "DBG ParseByLine()> misr_CameraAngles[" << i << "]:" << misr_CameraAngles[i] << std::endl;
			}
			#endif
			continue;
		}


		// parse single exact token without '\n', '\r' or space.
		found = line.find(MISR_RADIANCE.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(MISR_RADIANCE.c_str()));
			while(line[0] == ' ')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				misr_Radiance = token;
			}
			#if DEBUG_TOOL
			std::cout << "DBG ParseByLine()> " <<  MISR_RADIANCE << ": " << misr_Radiance << std::endl;
			#endif
			continue;
		}


		// parse single exact token without '\n', '\r' or space.
		found = line.find(TARGET_INSTRUMENT.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(TARGET_INSTRUMENT.c_str()));
			while(line[0] == ' ')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				targetInstrument = token;
			}
			#if DEBUG_TOOL
			std::cout << "DBG ParseByLine()> " <<  TARGET_INSTRUMENT << ": " << targetInstrument << std::endl;
			#endif
			continue;
		}


		// parse single exact token without '\n', '\r' or space.
		found = line.find(MODIS_RESOLUTION.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(MODIS_RESOLUTION.c_str()));
			while(line[0] == ' ')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {  // get exact string
				modis_Resolution = token;
			}

			isValidInput = CheckMODIS_Resolution(modis_Resolution);
			#if DEBUG_TOOL
			std::cout << "DBG ParseByLine()> " <<  MODIS_RESOLUTION << ": " << modis_Resolution << std::endl;
			#endif
			continue;
		}


		// parse multiple
		found = line.find(MODIS_BANDS.c_str());
		if(found != std::string::npos)
		{
			line = line.substr(strlen(MODIS_BANDS.c_str()));
			while(line[0] == ' ')
				line = line.substr(1);
			pos = line.find_first_of(' ', 0);
			std::stringstream ss(line); // Insert the string into a stream
			std::string token;
			while (ss >> token) {
				modis_Bands.push_back(atoi(token.c_str()));
			}
			#if DEBUG_TOOL
			//std::cout << "DBG ParseByLine()> Num of modis bands:" << modis_Bands.size() << std::endl;
			for(int i = 0; i < modis_Bands.size(); i++) {
				std::cout << "DBG ParseByLine()> modis_Bands[" << i << "]:" << modis_Bands[i] << std::endl;
			}
			#endif
			continue;
		}
	} // end of while

	if (isValidInput == false) {
		std::cerr << "Error: invlid input detected.\n";
		exit(1);
	}
}



/*===================================================
 *	Functions to check input values
 */

//--------------------
// Return: 
//	 true - OK 
//	 false - Failed
bool AF_InputParmeterFile::CheckInputBFdataPath(const std::string &filePath)
{
	bool ret = true;
	#if DEBUG_TOOL
	std::cout << "DBG CheckInputBFdataPath> BF file path: " << filePath <<	".\n"; 
	#endif
	// check if 
	std::ifstream file(filePath.c_str());
	ret = (file.good() ? true : false);
	if(ret == false) {
		std::cout << "Input data file '" << filePath << "' doesn't exist.\n"; 
	}
	return ret;
}

//--------------------
// Return: 
//	 true - OK 
//	 false - Failed
bool AF_InputParmeterFile::CheckMODIS_Resolution(std::string &str)
{
	bool ret = true;
	#if DEBUG_TOOL
	std::cout << "DBG CheckMODIS_Resolution> reolution: " << str <<  ".\n"; 
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


//===================================================
// Functions to get input values
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

std::string AF_InputParmeterFile::GetMISR_Radiance()
{
	return misr_Radiance;
}

std::string AF_InputParmeterFile::GetTargetInstrument()
{
	return targetInstrument;
}

std::string AF_InputParmeterFile::GetMODIS_Resolution()
{
	return modis_Resolution;
}

std::vector<int>  AF_InputParmeterFile::GetMODIS_Bands()
{
	return modis_Bands;
}
