/*
    AUTHOR:
        Yat Long Lo

    EMAIL:
        yllo2@illinois.edu

	PROGRAM DESCRIPTION:
		This is the entry point for the advanced fusion program. The IO module is in io.c and the repojection module is in reproject.c.
		This part of the program is still in development. At the moment, it serves as a starting point to read in an input text file and
		extract the specified parameters from it. The specified parameters are then used to call the corresponding functions to perform 
		data reading, reprojection and data writing. A sample input text file should be provided in the repository.

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>
#include <sys/time.h>
#include "reproject.h"
#include "io.h"

#define MAX_MODIS_BANDS 38

void deleteSpaces(char* src, char* dst);

int main(int argc, char ** argv) 
{
	//Input variables
	char* file_path = (char*)calloc(100, sizeof(char));
	char* output_file = (char*)calloc(100, sizeof(char));
	char* project_instrument = (char*)calloc(20, sizeof(char));
	char* base_instrument = (char*)calloc(20, sizeof(char)); 
	char* method = (char*) calloc(20, sizeof(char));
	char misr_args[3][50];
	int misr_args_size = 3;
	char modis_args[1][50];
	int modis_args_size = 2;
	char modis_bands[MAX_MODIS_BANDS][50];
	int modis_band_count = 0;
	char aster_args[2][50];
	int aster_args_size = 2;
	int project_instru_args = 0;
	int base_instru_args = 0;

	if(argc < 2){
		printf("Usage: ./af_run input_parameters.txt\n");
		return -1;
	}
	
	//Open input parameters file
	FILE* file = fopen(argv[1], "r");
	if(file == NULL){
		printf("Input parameters file does not exist\n");
		return -1;
	}
	//Read in input parameters
	char line[256];
	int line_count = 0;
	while(fgets(line, sizeof(line), file)){
		//Remove trailing newline or space characte
		size_t len = strlen(line)-1;
		if(line[len] == '\n' || line[len] == ' '){
			line[len] = '\0';
		}
		
		int token_count = 0;
		char* token = strtok(line, "=");
		while(token != NULL){
			//Expect input title
			if(token_count == 0){
				printf("title: %s\n", token);
				if(line_count == 0 && strcmp(token, "file_path") != 0){
					printf("First line of input should be file_path\n");
					return -1;
				}
				else if(line_count == 1 && strcmp(token, "output_file_path") != 0){
					printf("Second line of input should be output_file_path\n");
					return -1;
				}
				else if(line_count == 2 && strcmp(token, "project_instrument") != 0){
					printf("Third line of input should be project_instrument\n");
					return -1;
				}
				else if(line_count == 3 + project_instru_args && strcmp(token, "method") != 0){
					printf("method is missing\n");
					return -1;
				}
				else if(line_count == 4 + project_instru_args && strcmp(token,"base_instrument") != 0){
					printf("base_instrument is missing\n");
					return -1;
				}
				token_count += 1;
			}
			//Expect input value
			else if(token_count == 1){
				printf("reading_value\n");
				if(line_count == 0){
					strcpy(file_path, token);
				}
				else if(line_count == 1){
					strcpy(output_file, token);
				}
				else if(line_count == 2){
					//Start reading in project instrument arguments
					strcpy(project_instrument, token);
					//MISR
					if(strstr(project_instrument, "MISR") != NULL){
						int i;
						for(i = 0; i < misr_args_size; i++){
							fgets(line, sizeof(line), file);
							size_t len = strlen(line)-1;
							if(line[len] == '\n' || line[len] == ' '){
								line[len] = '\0';
							}
							int token_count = 0;
							char* token = strtok(line, "=");
							while(token != NULL){
								if(token_count == 0){
									if(i == 0 && strcmp(token, "resolution") != 0){
										printf("First arg of MISR should be resolution\n");
										return -1;
									}
									else if(i == 1 && strcmp(token, "camera_angle") != 0){
										printf("Second arg of MISR should be camera_angle\n");
										return -1;
									}
									else if(i == 2 && strcmp(token, "radiance") != 0){
										printf("Third arg of MISR should be radiance\n");
										return -1;
									}
									token_count += 1;
								}
								else if(token_count == 1){
									strcpy(misr_args[i], token);
									token_count = 0;
								}
								else{
									printf("Something is seriously wrong with the input parameters file\n");
									return -1;
								}
								token = strtok(NULL, "=");
							}
						}
						line_count += misr_args_size;
						project_instru_args = 3;
					}
					else if(strstr(project_instrument, "ASTER") != NULL){
						int i;
						for(i = 0; i < aster_args_size; i++){
							fgets(line, sizeof(line), file);
							size_t len = strlen(line)-1;
							if(line[len] == '\n' || line[len] == ' '){
								line[len] = '\0';
							}
							int token_count = 0;
							char* token = strtok(line, "=");
							while(token != NULL){
								if(token_count == 0){
									if(i == 0 && strcmp(token, "subsystem") != 0){
										printf("First arg of ASTER should be subsystem\n");
										return -1;
									}
									else if(i == 1 && strcmp(token, "dataset_name") != 0){
										printf("Second arg of ASTER should be dataset_name\n");
										return -1;
									}
									token_count += 1;
								}
								else if(token_count == 1){
									strcpy(aster_args[i], token);
									token_count = 0;
								}
								else{
									printf("Something is seriously wrong with the input parameters file\n");
									return -1;
								}
								token = strtok(NULL, "=");
							}
						}
					}
					//Handle other instruments as project instruments
				}
				else if(line_count == 3 + project_instru_args && project_instrument != NULL && project_instrument[0] != '\0'){
					strcpy(method, token);
				}
				else if(line_count == 4 + project_instru_args && method != NULL && method[0] != '\0'){
					//Start reading in base instrument arguments
					strcpy(base_instrument, token);
					if(strstr(base_instrument, "MODIS") != NULL){
						int i;
						for(i = 0; i < modis_args_size; i++){
							fgets(line, sizeof(line), file);
							size_t len = strlen(line)-1;
							if(line[len] == '\n' || line[len] == ' '){
								line[len] = '\0';
							}
							char* token = strtok(line, "=");
							int token_count = 0;
							while(token != NULL){
								printf("test: %s\n", token);
								if(token_count == 0){
									if(i == 0 && strcmp(token, "resolution") != 0){
										printf("First arg of MODIS should be resolution\n");
										return -1;
									}
									else if(i == 1 && strcmp(token, "bands") != 0){
										printf("Second arg of MODIS should be bands\n");
										return -1;
									}
									token_count += 1;
								}
								else if(token_count == 1){
									if(base_instru_args == 0){
										strcpy(modis_args[0], token);
									}
									else{
										//Start reading bands
										char* band_token = strtok(token, ",");
										while(band_token != NULL){
											
											strcpy(modis_bands[modis_band_count], band_token);
											modis_band_count += 1;
											band_token = strtok(NULL,",");
										}
									}
									token_count = 0;
									base_instru_args += 1;
								}
								token = strtok(NULL, "=");
							}
							line_count += 2;
						}
					}
				}
				token_count = 0;
			}
			else{
				printf("Something is seriously wrong with the input parameters file\n");
				return -1;
			}
			token = strtok(NULL,"=");
		}
		line_count += 1;
	}
	
	//Clean inputs
	deleteSpaces(project_instrument, project_instrument);
	deleteSpaces(base_instrument, base_instrument);
	deleteSpaces(method, method);
	deleteSpaces(modis_args[0], modis_args[0]);
	deleteSpaces(file_path, file_path);
	deleteSpaces(output_file, output_file);

	//Test reading input values
	printf("Input parameters\n");
	printf("file_path: %s\n", file_path);
	printf("outputfile: %s\n", output_file);
	printf("project instrument: %s\n", project_instrument);
	printf("base instrument: %s\n", base_instrument);
	printf("method: %s\n", method);
	printf("misr args: %s %s %s\n", misr_args[0], misr_args[1], misr_args[2]);
	printf("modis args: %s\n", modis_args[0]);
	int i;
	for(i = 0; i < modis_band_count; i++){
		printf("modis band: %s\n", modis_bands[i]);
	}

/////////////////////////////////////////////////////////////////////////////////////	
	//Begin AF processes
	//Open source hdf5 file (BF file)
	hid_t src_file;
	if(0 > (src_file = af_open(file_path))) {
		printf("File not found\n");
		exit(1);
	}
	
	
	
	//Create output file
	hid_t output_file_pointer = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	//Get project instrument geolocation
	int nCellsrc;
	double* src_lat;
	double* src_long;
	if(strstr(project_instrument, "MISR") != NULL){
		src_lat = get_misr_lat(src_file, misr_args[0], &nCellsrc);
		src_long = get_misr_long(src_file, misr_args[0], &nCellsrc);
	}
	else if(strstr(project_instrument, "ASTER") != NULL){
		src_lat = get_ast_lat(src_file, aster_args[0], aster_args[1], &nCellsrc);
		src_long = get_ast_long(src_file, aster_args[0], aster_args[1], &nCellsrc);
	}
	//TODO: add other cases for source geolocation


	//Get base instrument geolocation
	int nCelldest;
	double* dest_lat;
	double* dest_long;
	if(strstr(base_instrument, "MODIS") != NULL){
		if(strstr(modis_args[0], "1KM") != NULL){
			strcpy(modis_args[0], "_1KM");
		}
		else if(strstr(modis_args[0], "250m") != NULL){
			strcpy(modis_args[0], "_250m");
		}
		else if(strstr(modis_args[0], "500m") != NULL){
			strcpy(modis_args[0], "_500m");
		}
		dest_lat = get_modis_lat(src_file, modis_args[0], &nCelldest);
		dest_long = get_modis_long(src_file, modis_args[0], &nCelldest);
	}
	//TODO: add other cases for destination geolocation
	
	double* projected_Rad_Out;
	int * tarNNSouID;
	double** p_src_lat = &src_lat;
	double** p_src_lon = &src_long;
	
	tarNNSouID = (int *)malloc(sizeof(int) * nCelldest);
	
	printf("nearest_neighbor\n");
	//TODO change last argument based on resolution
	if(strstr(base_instrument, "MODIS") != NULL){
		if(strstr(modis_args[0], "_1KM") != NULL){
			nearestNeighbor(p_src_lat, p_src_lon, nCellsrc, dest_lat, dest_long, tarNNSouID, NULL, nCelldest, 1000);
		}
		else if(strstr(modis_args[0], "_250m") != NULL){
			nearestNeighbor(p_src_lat, p_src_lon, nCellsrc, dest_lat, dest_long, tarNNSouID, NULL, nCelldest, 300);
		}
		else if(strstr(modis_args[0], "_500m") != NULL){
			nearestNeighbor(p_src_lat, p_src_lon, nCellsrc, dest_lat, dest_long, tarNNSouID, NULL, nCelldest, 600);
		}
	}
	
	src_lat = *p_src_lat;
	src_long = *p_src_lon;
	
	free(src_lat);
	free(src_long);
	
	printf("writing dest geo\n");
	int lat_status =  af_write_mm_geo(output_file_pointer, 0, dest_lat, nCelldest);
	int long_status = af_write_mm_geo(output_file_pointer, 1, dest_long, nCelldest);
	if(lat_status < 0 || long_status < 0){
		printf("Writing dest geolocation - error\n");
		return -1;
	}
	
	free(dest_lat);
	free(dest_long);
	
	printf("getting source rad\n");
	double* src_rad;
	if(strstr(project_instrument, "MISR") != NULL){
		src_rad = get_misr_rad(src_file, misr_args[1], misr_args[0], misr_args[2], &nCellsrc);
	}
	else if(strstr(project_instrument, "ASTER") != NULL){
		src_rad = get_ast_rad(src_file, aster_args[0], aster_args[1], &nCellsrc);
	}
	
	printf("getting dest rad\n");
	int nCelldest_rad;
	double* dest_rad;
	if(strstr(base_instrument, "MODIS") != NULL){
		dest_rad = get_modis_rad(src_file, modis_args[0], modis_bands, modis_band_count, &nCelldest_rad);		
	}
	
	double* src_rad_out = (double *)malloc(sizeof(double) * nCelldest);
	//Interpolating
	int * nsrcPixels;
	printf("interpolating\n");
	if(strstr(method, "nnInterpolate") != NULL){
		nnInterpolate(src_rad, src_rad_out, tarNNSouID, nCelldest);
	}
	else if(strstr(method, "summaryInterpolate") != NULL){
		nsrcPixels = (int *) malloc(sizeof(int) * nCelldest);
		summaryInterpolate(src_rad, tarNNSouID, nCellsrc, src_rad_out, nsrcPixels, nCelldest);
	}
	
	printf("writing data fields\n");
	//TODO: more generalize writing function
	int data_write_status = af_write_misr_on_modis(output_file_pointer, src_rad_out, dest_rad, nCelldest_rad, modis_band_count, nCelldest);
	if(data_write_status < 0){
		printf("Writing data fields - error\n");
	}	
	
	free(src_rad);
	free(src_rad_out);
	free(tarNNSouID);

	printf("Writing done\n");
	//Closing file
	herr_t ret = af_close(src_file);

	
	//Close input parameters file
	fclose(file);
	return 0;

}


//Inputs formatting utilities
void deleteSpaces(char* src, char* dst)
{
  int s, d=0;
  for (s=0; src[s] != 0; s++)
    if (src[s] != ' ') {
       dst[d] = src[s];
       d++;
    }
  dst[d] = 0;
}
