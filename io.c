/*


    AUTHOR:
        Yat Long Lo

    EMAIL:
        yllo2@illinois.edu

	PROGRAM DESCRIPTION:
		This is the IO module for the TERRA Fusion project, used for advanced fusion. TERRA data can be retrieved for any instruments by specifying
		different paraameters, which are used for reprojection. 

*/

#include "io.h"
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <assert.h>
#define FALSE   0

//Band constants for MODIS
char* m_250_list[2] = {"1", "2"};
char* m_500_list[5] = {"3", "4", "5", "6", "7"};
char* km_1_ref_list[15] = {"8", "9", "10", "11", "12", "13L", "13H", "14L", "14H", "15", "16", "17", "18", "19", "26"};
char* kme_1_list[16] = {"20", "21", "22", "23", "24", "25", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36"};


/*
						get_misr_rad
	DESCRIPTION:
		This function retrieves a particular MISR radiance dataset based on specified parameters including camera angle,
		resolution and radiance. Checks would be done to make sure the parameters are valid. Downsampling is also supported
		if the configuration allows, in which the MISR radiance data would be averaged for lower resolution.	
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. camera_angle -- A string variable that specifies the camera angle
		2. resolution(H/L) -- A string variable that specifies the resolution 
		3. radiance -- A string variable that specifies the radiance data field
		4. size -- An integer pointer that points to the size of the radiance data after the retreival
		
	EFFECT:
		Memory would be allocated according to the size of the MISR radiance data. The variable size would also
		be set to the size of the data array
		
	RETURN:
		Returns NULL upon error
		Returns down_data (1D array) if the data requires downsampling
		Returns data (1D array) in normal situations
		
*/

double* get_misr_rad(hid_t file, char* camera_angle, char* resolution, char* radiance, int* size){
	//Path to dataset proccessing 
	int down_sampling = 0;
	char* instrument = "MISR";
	char* d_fields = "Data_Fields";
	const char* arr[] = {instrument, camera_angle, d_fields, radiance};
	
	//Dataset name parsing
	char* rad_dataset_name;
	concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(camera_angle) + strlen(d_fields) + strlen(radiance) + 4, 4);
		//Check for correct specification
	if(strcmp(camera_angle, "AN") != 0 && strcmp(radiance, "Red_Radiance") != 0 && strcmp(resolution, "H") == 0){
		printf("Error: Your specification does not support high resolution.\n");
		return NULL;
	}
	else if((strcmp(camera_angle, "AN") == 0 || strcmp(radiance, "Red_Radiance") == 0) && strcmp(resolution, "L") == 0){
		//Downsampling has to be done
		down_sampling = 1;
	}
	
	printf("Reading MISR\n");
	/*Dimensions - 180 blocks, 512 x 2048 ordered in 1D Array*/
	//Retrieve radiance dataset and dataspace
	double* data = af_read(file, rad_dataset_name);
	*size = dim_sum(af_read_size(file, rad_dataset_name), 3);
	
	if(data == NULL){
		return NULL;
	}
	printf("Reading successful\n");
	//Variable containing down sampled data
	double* down_data;
	if(down_sampling == 1){
		printf("Undergoing downsampling\n");
		hsize_t* dims = af_read_size(file, rad_dataset_name);
		*size = dims[0] * (dims[1]/4) * (dims[2]/4);
		down_data = malloc(dims[0] * (dims[1]/4) * (dims[2]/4) * sizeof(double));
		int i, j, k;
		for(i = 0; i < dims[0]; i++){
			for(j = 0; j < dims[1]; j = j + 4){
				for(k = 0; k < dims[2]; k = k + 4){
					//Retrieving 4x4 window for averaging
					//Formula for converting i, j and k to index in data array
					//int index = i*dims[1]*dims[2] + j*dims[2] + k;
					int a,b;
					int max_x = j + 4;
					int max_z = k + 4; 
					int* index_array = malloc(16*sizeof(int));
					int index_iter = 0;
					for(a = j; a < max_x; a++){
						for(b = k; b < max_z; b++){
							index_array[index_iter] = i*dims[1]*dims[2] + a*dims[2] + b;
							index_iter += 1;
						}
					}
					double* window = malloc(16*sizeof(double));
					int c;
					for(c = 0; c < 16; c++){
						window[c] = data[index_array[c]];
					}
					//Window Retrieved, get average and assign to new data grid
					double average = misr_averaging(window);
					int new_index = i*dims[1]/4*dims[2]/4 + (j/4)*dims[2]/4 + k/4;
					down_data[new_index] = average;
					free(index_array);
					free(window);
				}
			}
		}
		free(data);
		printf("Downsampling done\n");
	}
	
	if(down_sampling == 1){
		printf("rad_data: %f\n", down_data[0]);
		return down_data;
	}
	else{
		printf("rad_data: %f\n", data[0]);
		return data;	
	}
	
}

/*
						get_misr_lat
	DESCRIPTION:
		This function retrieves the corresponding geological latitude data for MISR based on the resolution specified.
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. resolution(H/L) -- A string variable that specifies the resolution 
		2. size -- An integer pointer that points to the size of the latitude data after the retreival
		
	EFFECT:
		Memory would be allocated according to the size of the MISR latitude data. The variable size would also
		be set to the size of the data array
	
	RETURN:
		Returns lat_data (1D array) if successful
		Returns NULL upon error
*/


double* get_misr_lat(hid_t file, char* resolution, int* size){
	//Path to dataset proccessing 
	char* instrument = "MISR";
	char* location;
	if(strcmp(resolution, "H") == 0){
		location = "HRGeolocation";
	}
	else{
		location = "Geolocation";
	}
	char* lat = "GeoLatitude";
	const char* arr2[] = {instrument, location, lat};
	
	//Dataset names parsing
	char* lat_dataset_name;
	concat_by_sep(&lat_dataset_name, arr2, "/", strlen(instrument) + strlen(location) + strlen(lat) + 4, 3);
	
	printf("Retrieveing latitude data for MISR\n");
	//Retrieve latitude dataset and dataspace
	double* lat_data = af_read(file, lat_dataset_name);
	*size = dim_sum(af_read_size(file, lat_dataset_name), 3);
	if(lat_data == NULL){
		return NULL;
	}
	printf("lat_data: %f\n", lat_data[0]);
	return lat_data;
}

/*
						get_misr_long
	DESCRIPTION:
		This function retrieves the corresponding geological longitude data for MISR based on the resolution specified.
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. resolution(H/L) -- A string variable that specifies the resolution 
		2. size -- An integer pointer that points to the size of the latitude data after the retreival
		
	EFFECT:
		Memory would be allocated according to the size of the MISR longitude data. The variable size would also
		be set to the size of the data array
	
	RETURN:
		Returns long_data (1D array) if successful
		Returns NULL upon error
*/


double* get_misr_long(hid_t file, char* resolution, int* size){
	//Path to dataset proccessing 
	char* instrument = "MISR";
	char* location;
	if(strcmp(resolution, "H") == 0){
		location = "HRGeolocation";
	}
	else{
		location = "Geolocation";
	}
	char* longitude = "GeoLongitude";
	const char* arr3[] = {instrument, location, longitude};
	
	//Dataset names parsing
	char* long_dataset_name;
	concat_by_sep(&long_dataset_name, arr3, "/", strlen(instrument) + strlen(location) + strlen(longitude) + 4, 3);
	
	printf("Retrieveing longitude data for MISR\n");
	//Retrieve longitude dataset and dataspace
	double* long_data = af_read(file, long_dataset_name);
	*size = dim_sum(af_read_size(file, long_dataset_name), 3);
	if(long_data == NULL){
		return NULL;
	}
	printf("long_data: %f\n", long_data[0]);
	return long_data;
}

/*
						get_misr_attr
	DESCRIPTION:	
		This function retrieves the attribute of MISR's geological dataset based on specified parameters including camera angle and resolution
	
	ARGUMENTS:	
		0. file -- A hdf file variable that points to the BasicFusion file
		1. camera_angle -- A string variable that specifies the camera angle
		2. resolution(H/L) -- A string variable that specifies the resolution 
		3. radiance -- A string variable that specifies the radiance data field
		4. attr_name -- A string variable that specifies the name of the attribute
		5. geo(0/1) -- An integer that indicates the preference over geological data (0 - latitude, 1 - longitude)
		6. attr_pt -- A void pointer that points to the retrieved attribute which can by float or char*

	EFFECT:
			Memory would be allocated to hold the attribute value
			
	RETURN:
		Returns attr_pt, a void pointer that points to the attribute value
		Returns NULL upon error
		
*/


//geo - 0:not geolocation attributes, 1:lat, 2:long
void* get_misr_attr(hid_t file, char* camera_angle, char* resolution, char* radiance, char* attr_name, int geo, void* attr_pt){
	//Path variables
	char* instrument = "MISR";
	char* d_fields = "Data_Fields";
	char* location = "Geolocation";

	//Dataset name parsing
	char* rad_dataset_name;
	if(geo == 0){
		const char* arr[4] = {instrument, camera_angle, d_fields, radiance};
		concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(camera_angle) + strlen(d_fields) + strlen(radiance) + 4, 4);
	}
	else if(geo == 1){
		char* lat = "GeoLatitude";
		const char* arr[3] = {instrument, location, lat};
		concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(location) + strlen(lat) + 4, 3);
	}
	else if(geo == 2){
		char* longitude = "GeoLongitude";
		const char* arr[3] = {instrument, location, longitude};
		concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(location) + strlen(longitude) + 4, 3);
	}
	else{
		printf("Wrong geo number");
		return NULL;
	}
	
	//Get attribute
	hid_t attr = H5Aopen_by_name(file, rad_dataset_name, attr_name, H5P_DEFAULT, H5P_DEFAULT);
	if(attr < 0){
		printf("Attribute %s does not exists\n", attr_name);
	}
	hid_t attr_type = H5Aget_type(attr);
	if(strcmp(attr_name, "Units") == 0 || strcmp(attr_name, "units") == 0){
		attr_pt = malloc(sizeof(char) * 50);
		H5Aread(attr, attr_type, attr_pt);
	}
	else if(strcmp(attr_name, "_FillValue") == 0){
		attr_pt = malloc(sizeof(float));
		H5Aread(attr, attr_type, attr_pt);
	}
	H5Aclose(attr);
	
	return attr_pt;
}

/*
						get_modis_rad
	DESCRIPTION:	
		This function retrieves MODIS randiance data fields based on the specified resolution and any number of valid bands that MODIS provides.
		Only granules that consist of all resolutions, namely 1KM, 250m and 500m, are considered valid granules. The invalid ones would not be picked
		up, while the valid ones would be stitched together band after band, in the form of a 1D array. 
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. resolution(_1KM, _250m, _500m) -- A string variable that specifies the resolution 
		2. bands -- a char array that specifies the name of bands to be retrieved. There are 38 possible bands for MODIS
		3. band_size -- an integer that specifies the number of bands to be retrieved
		4. size -- An integer pointer that points to the size of the radiance data after the retreival
		
	EFFECT:
		Memmory would be allocated according to the size of the MODIS data with the variable size set to the total size of the data
		
	RETURN:
		Returns result_data upon sucessful retrieval 
		Returns NULL upon error
*/


double* get_modis_rad(hid_t file, char* resolution, char bands[38][50], int band_size, int* size){
	printf("Reading MODIS rad\n");
	
	//Path variables
	char* instrument = "MODIS";
	char* d_fields = "Data_Fields";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	char res[3][10] = {"_1KM", "_250m", "_500m"};
	int i;
	int store_count = 0;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		
		//Check if it has all resolutions
		int include_g = 0;
		int j;
		for(j = 0; j < 3; j++){
			char* res_group_name;
			const char* d_arr[] = {name, res[j]};
			concat_by_sep(&res_group_name, d_arr, "/", strlen(name) + strlen(res[j]) + 2, 2);
			memmove(&res_group_name[0], &res_group_name[1], strlen(res_group_name));
			printf("group_name: %s\n", res_group_name);
			htri_t status = H5Lexists(group, res_group_name, H5P_DEFAULT);
			if(status <= 0){
				printf("Group does not exist\n");
				include_g = 1;
				break;
			}
		}
		if(include_g == 0){
			strcpy(names[store_count], name);
			store_count += 1;
		}
		free(name);
	}
	
	printf("num granules: %d\n", store_count);
	
	//Get dataset names from bands
	printf("Retreving dataset names\n");
	char dnames[band_size][50];
	int band_indices[band_size];
	int j;
	for(j = 0; j < band_size; j++){
		char* dname = get_modis_filename(resolution, bands[j], &band_indices[j]);
		if(dname == NULL){
			printf("Band %s is not supported for %s resolution\n", bands[j], resolution);
			return NULL;
		}
		printf("dname: %s\n", dname);
		strcpy(dnames[j], dname);
	}
	
	
	//Get total data size
	printf("Get total data size\n");
	int k;
	int m;
	int total_size = 0;
	for(m = 0; m < band_size; m++){
		for(k = 0; k < store_count; k++){
			char* name = names[k];
			const char* d_arr[] = {instrument, name, resolution, d_fields, dnames[m]};
			char* dataset_name;
			concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(d_fields) + strlen(dnames[m]), 5);
			hsize_t* curr_dim = af_read_size(file, dataset_name);
			if(curr_dim == NULL){
				continue;
			}
			total_size += curr_dim[1]*curr_dim[2];
			free(curr_dim);
		}
	}
	
	double* result_data = calloc(total_size, sizeof(double));
	int start_point = 0;
	
	//Start reading data
	int n;
	for(n = 0; n < band_size; n++){
		int file_size;
		double * MODIS_rad = get_modis_rad_by_band(file, resolution, dnames[n], &band_indices[n], &file_size);
		memcpy(&result_data[start_point], MODIS_rad, file_size*sizeof(double));
		start_point += file_size;
	}
	
	if(total_size == start_point){
		printf("Final size validated\n");
	}
	
	*size = total_size;
	
	return result_data;
}


/*
						get_modis_rad_by_band
	DESCRIPTION:	
		This function retrieves MODIS data fields of one single band by specifying band index and dataset name, which is the index of that band in the granule. Similar
		to get_modis_rad, checks are done to make sure the granule has all 3 resolutions. Then, the band data of all valid granules would be stitched 
		together.
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. resolution(_1KM, _250m, _500m) -- A string variable that specifies the resolution 
		2. d_name -- a string that specifies the dataset name that the band is located, check get_modis_filename for band and dataset name matching
		3. band_index -- a pointer to an integer that indicates the location of a band within the dataset, check global variables at the top of this
						 file for reference
		4. size -- An integer pointer that points to the size of the radiance data after the retreival 
			
	EFFECT:
		Memory would be allocated to the retrieved MODIS data with the variable size set to the size of the data
		
	RETURN:
		Returns result_data upon successful retrieval 
		Returns NULL upon error
		
*/


double* get_modis_rad_by_band(hid_t file, char* resolution, char* d_name, int* band_index, int* size){
	printf("Reading MODIS rad by band\n");
	char* instrument = "MODIS";
	char* d_fields = "Data_Fields";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	char res[3][10] = {"_1KM", "_250m", "_500m"};
	int i;
	int store_count = 0;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		
		//Check if it has all resolutions
		int include_g = 0;
		int j;
		for(j = 0; j < 3; j++){
			char* res_group_name;
			const char* d_arr[] = {name, res[j]};
			concat_by_sep(&res_group_name, d_arr, "/", strlen(name) + strlen(res[j])+2, 2);
			memmove(&res_group_name[0], &res_group_name[1], strlen(res_group_name));
			printf("group_name: %s\n", res_group_name);
			htri_t status = H5Lexists(group, res_group_name, H5P_DEFAULT);
			if(status <= 0){
				printf("Group does not exist\n");
				include_g = 1;
				break;
			}
		}
		if(include_g == 0){
			strcpy(names[store_count], name);
			store_count += 1;
		}
		free(name);
	}
	
	//Get total data size
	printf("Get total data size\n");
	int k;
	int total_size = 0;
	for(k = 0; k < store_count; k++){
		char* name = names[k];
		const char* d_arr[] = {instrument, name, resolution, d_fields, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(d_fields) + strlen(d_name), 5);
		hsize_t* curr_dim = af_read_size(file, dataset_name);
		if(curr_dim == NULL){
			continue;
		}
		total_size += curr_dim[1]*curr_dim[2];
		free(curr_dim);
	}
	
	//Allocate data size
	double* result_data = calloc(total_size, sizeof(double));
	
	//Retreving data
	int h;
	int curr_size = 0;
	int read_first = -1;
	for(h = 0; h < store_count; h++){
		double* data;
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {instrument, name, resolution, d_fields, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(d_fields) + strlen(d_name), 5);
		printf("granule_name: %s\n", name);
		data = af_read(file, dataset_name);
		if(data == NULL){
			continue;
		} 
		hsize_t* curr_dim = af_read_size(file, dataset_name);
		int band_length = curr_dim[1] * curr_dim[2];
		printf("band index: %d\n", (*band_index));
		printf("band length: %d\n", band_length);
		int read_offset = (*band_index)*band_length;
		printf("test data: %f\n", data[2748619]);
		memcpy(&(result_data[curr_size]), &(data[read_offset]), band_length*sizeof(double));
		curr_size += band_length;
		free(data);
		free(curr_dim);
	}
	*size = curr_size;
	
	assert(curr_size == total_size);
	printf("Size validated\n");
	printf("test data: %f\n", result_data[curr_size-1]);
	return result_data;
}

/*
						get_modis_lat
	DESCRIPTION:	
		This function retrieves MODIS geological latitude data based on specified resolution. Similar to other MODIS data retrieval functions, it 
		checks every granule to ensure all the resolutions are present. Then, the latitude data of valid granules would be stitched together. At the
		moment, realloc is used which is not the most efficient approach. However, this does not pose any noticable efficiency issues in terms of the flow 
		of the whole program.
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. resolution(_1KM, _250m, _500m) -- A string variable that specifies the resolution 
		2. size -- An integer pointer that points to the size of the latitude data after the retreival 	
	
	EFFECT:
		Memory would be allocated for retrieved latitude data with the variable set as the size of the latitude data
			
	RETURN:
		Returns lat_data upon successful retrieval
		Returns NULL upon error
		
*/


double* get_modis_lat(hid_t file, char* resolution, int* size){
	printf("Reading MODIS lat\n");
	//Path variables
	char* instrument = "MODIS";
	char* d_fields = "Data_Fields";
	char* location = "Geolocation";
	char* lat = "Latitude";
	
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	char res[3][10] = {"_1KM", "_250m", "_500m"};
	int i;
	int store_count = 0;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		
		//Check if it has all resolutions
		int include_g = 0;
		int j;
		for(j = 0; j < 3; j++){
			char* res_group_name;
			const char* d_arr[] = {name, res[j]};
			concat_by_sep(&res_group_name, d_arr, "/", strlen(name) + strlen(res[j])+2, 2);
			memmove(&res_group_name[0], &res_group_name[1], strlen(res_group_name));
			printf("group_name: %s\n", res_group_name);
			htri_t status = H5Lexists(group, res_group_name, H5P_DEFAULT);
			if(status <= 0){
				printf("Group does not exist\n");
				include_g = 1;
				break;
			}
		}
		if(include_g == 0){
			strcpy(names[store_count], name);
			store_count += 1;
		}
		free(name);
	}
	
	printf("num granules: %d\n", store_count);
	
	int h;
	double* lat_data;
	double curr_lat_size;
	int read_first = -1;
	for(h = 0; h < store_count; h++){
		//Path formation
		char* name = names[h];
		printf("granule name: %s\n", name);
		const char* lat_arr[] = {instrument, name, resolution, location, lat};
		char* lat_dataset_name;
		concat_by_sep(&lat_dataset_name, lat_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(location) + strlen(lat), 5);
		
		if(read_first < 0){
			curr_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 2);
			lat_data = af_read(file, lat_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and its dimension
			double* adding_lat = af_read(file, lat_dataset_name);
			double new_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 2);
			
			//Reallocating arrays of data
			lat_data = realloc(lat_data, sizeof(double)*(curr_lat_size + new_lat_size));
			memcpy(&lat_data[(int)curr_lat_size], adding_lat, sizeof(double)*new_lat_size);
			curr_lat_size += new_lat_size;
			
			free(adding_lat);
		}
	}
	*size = curr_lat_size;
	
	if(lat_data != NULL){
	printf("test_lat_data: %f\n", lat_data[0]);
	printf("test_lat_data: %f\n", lat_data[2748620]);
	printf("test_lat_data: %f\n", lat_data[5510780]);
	}
	
	return lat_data;
}

/*
						get_modis_long
	DESCRIPTION:	
		This function retrieves MODIS geological longitude data based on specified resolution. Similar to other MODIS data retrieval functions, it 
		checks every granule to ensure all the resolutions are present. Then, the longitude data of valid granules would be stitched together. At the
		moment, realloc is used which is not the most efficient approach. However, this does not pose any noticable efficiency issues in terms of the flow 
		of the whole program.
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. resolution(_1KM, _250m, _500m) -- A string variable that specifies the resolution 
		2. size -- An integer pointer that points to the size of the latitude data after the retreival 	
	
	EFFECT:
		Memory would be allocated for retrieved longitude data with the variable set as the size of the longitude data
			
	RETURN:
		Returns long_data upon successful retrieval
		Returns NULL upon error
		
*/


double* get_modis_long(hid_t file, char* resolution, int* size){
	printf("Reading MODIS long\n");
	//Path variables
	char* instrument = "MODIS";
	char* d_fields = "Data_Fields";
	char* location = "Geolocation";
	char* longitude = "Longitude";
	
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	char res[3][10] = {"_1KM", "_250m", "_500m"};
	int i;
	int store_count = 0;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		
		//Check if it has all resolutions
		int include_g = 0;
		int j;
		for(j = 0; j < 3; j++){
			char* res_group_name;
			const char* d_arr[] = {name, res[j]};
			concat_by_sep(&res_group_name, d_arr, "/", strlen(name) + strlen(res[j])+2, 2);
			memmove(&res_group_name[0], &res_group_name[1], strlen(res_group_name));
			htri_t status = H5Lexists(group, res_group_name, H5P_DEFAULT);
			if(status <= 0){
				printf("Group does not exist\n");
				include_g = 1;
				break;
			}
		}
		if(include_g == 0){
			strcpy(names[store_count], name);
			store_count += 1;
		}
		free(name);
	}
	
	int h;
	int valid_granule_count = 0;
	double* long_data;
	double curr_long_size;
	int read_first = -1;
	for(h = 0; h < store_count; h++){
		//Path formation
		char* name = names[h];
		valid_granule_count += 1;
		const char* long_arr[] = {instrument, name, resolution, location, longitude};
		char* long_dataset_name;
		concat_by_sep(&long_dataset_name, long_arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(location) + strlen(longitude), 5);
		
		if(read_first < 0){
			curr_long_size = dim_sum(af_read_size(file, long_dataset_name), 2);
			long_data = af_read(file, long_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and its dimension
			double* adding_long = af_read(file, long_dataset_name);
			double new_long_size = dim_sum(af_read_size(file, long_dataset_name), 2);
			
			//Reallocating arrays of data
			long_data = realloc(long_data, sizeof(double)*(curr_long_size + new_long_size));
			memcpy(&long_data[(int)curr_long_size], adding_long, sizeof(double)*new_long_size);
			curr_long_size += new_long_size;

			free(adding_long);
		}
	}
	*size = curr_long_size;
	
	if(long_data != NULL){
		printf("test_long_data: %f\n", long_data[0]);
		printf("test_long_data: %f\n", long_data[1]);
		printf("test_long_data: %f\n", long_data[1353]);
		printf("test_long_data: %f\n", long_data[1354]);
		printf("test_long_data: %f\n", long_data[2748620]);
		printf("test_long_data: %f\n", long_data[5510780]);
	}
	
	return long_data;
}

/*
						get_modis_attr
	DESCRIPTION:	
		This function retrieves the attribute of MODIS's  dataset based on specified parameters including resolution and dataset name
	
	ARGUMENTS:	
		0. file -- A hdf file variable that points to the BasicFusion file
		1. resolution -- A string variable that specifies the resolution 
		2. radiance -- A string variable that specifies dataset name
		3. attr_name -- A string variable that specifies the name of the attribute
		4. geo(0/1/2) -- An integer that indicates the source of attribute, which can be from radiance data or geolocation data
		5. attr_pt -- A void pointer that points to the retrieved attribute which can by float or char*

	EFFECT:
			Memory would be allocated to hold the attribute value
			
	RETURN:
		Returns attr_pt, a void pointer that points to the attribute value
		Returns NULL upon error
		
*/


//geo - 0:not geolocation attributes, 1:lat, 2:long
double* get_modis_attr(hid_t file, char* resolution, char* d_name, char* attr_name, int geo, void* attr_pt){
	//Path variables
	char* instrument = "MODIS";
	char* d_fields = "Data_Fields";
	char* location = "Geolocation";
	
	//Get one group name, assuming all attributes across granules are the same
	printf("Retrieving granule group name\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char* rad_dataset_name;
	char* name = malloc(50*sizeof(char));
	int h;
	for(h = 0; h < num_groups; h++){
		H5Gget_objname_by_idx(group, (hsize_t)h, name, 50);
		const char* arr[] = {instrument, name, resolution, d_fields, d_name};
		//Dataset name parsing
		concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(d_fields) + strlen(d_name) + 5, 5);
		memmove(&rad_dataset_name[0], &rad_dataset_name[1], strlen(rad_dataset_name));
		htri_t status = H5Lexists(group, rad_dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		else{
			break;
		}
		
	}
	
	if(geo == 1){
		char* lat = "Latitude";
		const char* arr[] = {instrument, name, resolution, location, lat};
		concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(location) + strlen(lat) + 5, 5);
	}
	else if(geo == 2){
		char* longitude = "Longitude";
		const char* arr[] = {instrument, name, resolution, location, longitude};
		concat_by_sep(&rad_dataset_name, arr, "/", strlen(instrument) + strlen(name) + strlen(resolution) + strlen(location) + strlen(longitude) + 5, 5);
	}
	
	//Get attribute 	
	hid_t attr = H5Aopen_by_name(file, rad_dataset_name, attr_name, H5P_DEFAULT, H5P_DEFAULT);
	if(attr < 0){
		printf("Attribute %s does not exists\n", attr_name);
	}
	hid_t attr_type = H5Aget_type(attr);
	if(strcmp(attr_name, "units") == 0){
		attr_pt = malloc(sizeof(char) * 50);
		H5Aread(attr, attr_type, attr_pt);
	}
	else if(strcmp(attr_name, "_FillValue") == 0){
		attr_pt = malloc(sizeof(float));
		H5Aread(attr, attr_type, attr_pt);
	}
	else if(strcmp(attr_name, "valid_min") == 0){
		attr_pt = malloc(sizeof(float));
		H5Aread(attr, attr_type, attr_pt);
	}
	H5Aclose(attr);
	
	return attr_pt;
	
}

/*
						get_modis_filename
	DESCRIPTION:	
		This functions provides important information to retrieve a particular band of MODIS data, including the dataset name the band it belongs
		to and the index that the band locates within the dataset.
		
	ARGUMENTS:
		0. resolution -- A string variable that specifies the resolution 
		1. band -- A string variable that specifies the name of the band
		2. band_index -- A variable that holds the location of the band within the dataset
		
	EFFECT:
		band_index would be set accordingly with the information of dataset name
		
	RETURN:
		Returns dataset name upon successful matching 
		Returns NULL upon error
		
*/


char* get_modis_filename(char* resolution, char* band, int* band_index){
	if(strcmp(resolution, "_1KM") == 0){
		int i;
		for(i = 0; i < 15; i++){
			if(strcmp(band, km_1_ref_list[i]) == 0){
				*band_index = i;
				return "EV_1KM_RefSB";
			}
		}
		int j;
		for(j = 0; j < 16; j++){
			if(strcmp(band, kme_1_list[j]) == 0){
				*band_index = j;
				return "EV_1KM_Emissive";
			}
		}
		int k;
		for(k = 0;k < 2; k++){
			if(strcmp(band, m_250_list[k]) == 0){
				*band_index = k;
				return "EV_250_Aggr1km_RefSB";
			}
		}
		int a;
		for(a = 0; a < 5; a++){
			if(strcmp(band, m_500_list[a]) == 0){
				*band_index = a;
				return "EV_500_Aggr1km_RefSB";
			}
		}
		return NULL; 
	}
	else if(strcmp(resolution, "_250m") == 0){
		int i;
		for(i = 0; i < 2; i++){
			if(strcmp(band, m_250_list[i]) == 0){
				*band_index = i;
				return "EV_250_RefSB";
			}
		}
		return NULL;
	}
	else if(strcmp(resolution, "_500m") == 0){
		int i;
		for(i = 0; i < 5; i++){
			if(strcmp(band, m_500_list[i]) == 0){
				*band_index = i;
				return "EV_500_RefSB";
			}
		}
		int j;
		for(j = 0; j < 2; j++){
			if(strcmp(band, m_250_list[j]) == 0){
				*band_index = j;
				return "EV_250_Aggr500_RefSB";
			}
		}
		return NULL;
	}
	return NULL;
}


/*
						get_ceres_rad
	DESCRIPTION:	
		This function retrieves CERES radiance data based on specified parameters including camera and dataset name. Granules are stitched together
		and returned as a 1D array. Modification is suggested to remove the use of realloc.
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. camera (FM1/FM2) -- a string variable that specifies the camera
		2. d_name --  a string variable that specifies the name of the dataset
		3. size -- An integer pointer that points to the size of the CERES radiance data after the retreival 
		
	EFFECT:
		Memory would be allocated for the CERES radiance data with the variable size set to the size of the data
		
	RETURN:
		Returns data upon successful retrieval
		Returns NULL upon error
*/


double* get_ceres_rad(hid_t file, char* camera, char* d_name, int* size){
	printf("Reading CERES radiance\n");
	//Path variables
	char* instrument = "CERES";
	char* rad = "Radiances";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(names[i], name);
		free(name);
	}
	int h;
	double* data;
	double curr_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {instrument, name, camera, rad, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(camera) + strlen(rad) + strlen(d_name) + 5, 5);
		printf("granule_name: %s\n", name);
		if(read_first < 0){
			data = af_read(file, dataset_name);
			if(data == NULL){
				continue;
			}
			curr_size = dim_sum(af_read_size(file, dataset_name), 1); 
			read_first = 1;
		}
		else{
			//retrieve next set of data and its dimension
			double* adding_data = af_read(file, dataset_name);
			if(adding_data == NULL){
				continue;
			}
			double new_d_size = dim_sum(af_read_size(file, dataset_name), 1);
			//Reallocating arrays of data
			data = realloc(data, sizeof(double)*(curr_size + new_d_size));
			memcpy(&data[(int)curr_size], adding_data, sizeof(double)*new_d_size);
			curr_size += new_d_size;
			
			free(adding_data);
		}
	}
	*size = curr_size;
	
	//Print statements to verify data's existence	
	if(data != NULL){
		printf("test data: %f\n", data[0]);
		printf("test_data: %f\n", data[1]);
		printf("test data: %f\n", data[2]);
	}	
	return data;
}


/*
						get_ceres_lat
	DESCRIPTION:	
		This function retrieves CERES geological latitude data based on specified parameters including camera and dataset name. Granules are stitched together
		and returned as a 1D array. Modification is suggested to remove the use of realloc.
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. camera (FM1/FM2) -- a string variable that specifies the camera
		2. d_name --  a string variable that specifies the name of the dataset
		3. size -- An integer pointer that points to the size of the CERES latitude data after the retreival 
	
	EFFECT:
		Memory would be allocated to the latitude data with variable size set as the size of the data
			
	RETURN:
		Returns lat_data upon sucessful retrieval
		Returns NULL upon error
*/


double* get_ceres_lat(hid_t file, char* camera, char* d_name, int* size){
	printf("Reading CERES lat\n");
	//Path variables
	char* instrument = "CERES";
	char* rad = "Radiances";
	char* tp = "Time_and_Position";
	char* lat = "Latitude";
	
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(names[i], name);
		free(name);
	}
	int h;
	double* lat_data;
	double curr_lat_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {name, camera, rad, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(name) + strlen(camera) + strlen(rad) + strlen(d_name), 4);
		memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
		//Check if dataset exists first
		printf("granule_name: %s\n", name);
		htri_t status = H5Lexists(group, dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		
		const char* lat_arr[] = {instrument, name, camera, tp, lat};
		char* lat_dataset_name;
		concat_by_sep(&lat_dataset_name, lat_arr, "/", strlen(instrument) + strlen(name) + strlen(camera) + strlen(tp) + strlen(lat) + 5, 5);
		
		if(read_first < 0){
			curr_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 1);
			lat_data = af_read(file, lat_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and its dimension
			double* adding_lat = af_read(file, lat_dataset_name);
			double new_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 1);
			
			//Reallocating arrays of data
			lat_data = realloc(lat_data, sizeof(double)*(curr_lat_size + new_lat_size));
			memcpy(&lat_data[(int)curr_lat_size], adding_lat, sizeof(double)*new_lat_size);
			curr_lat_size += new_lat_size;
			
			free(adding_lat);
		}
	}
	*size = curr_lat_size;
	//Print statements to verify data's existence
	if(lat_data != NULL){
		printf("test_lat_data: %f\n", lat_data[0]);
		printf("test_lat_data: %f\n", lat_data[1]);
		printf("test_lat_data: %f\n", lat_data[2]);
	}
	
	return lat_data;
}

/*
						get_ceres_long
	DESCRIPTION:	
		This function retrieves CERES geological longitude data based on specified parameters including camera and dataset name. Granules are stitched together
		and returned as a 1D array. Modification is suggested to remove the use of realloc.
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. camera (FM1/FM2) -- a string variable that specifies the camera
		2. d_name --  a string variable that specifies the name of the dataset
		3. size -- An integer pointer that points to the size of the CERES longitude data after the retreival 
	
	EFFECT:
		Memory would be allocated to the longitude data with variable size set as the size of the data
			
	RETURN:
		Returns long_data upon sucessful retrieval
		Returns NULL upon error
*/


double* get_ceres_long(hid_t file, char* camera, char* d_name, int* size){
	printf("Reading CERES long\n");
	//Path variables
	char* instrument = "CERES";
	char* rad = "Radiances";
	char* tp = "Time_and_Position";
	char* longitude = "Longitude";
	
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(names[i], name);
		free(name);
	}
	int h;
	double* long_data;
	double curr_long_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {name, camera, rad, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(name) + strlen(camera) + strlen(rad) + strlen(d_name), 4);
		memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
		//Check if dataset exists first
		printf("granule_name: %s\n", name);
		htri_t status = H5Lexists(group, dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		
		const char* long_arr[] = {instrument, name, camera, tp, longitude};
		char* long_dataset_name;
		concat_by_sep(&long_dataset_name, long_arr, "/", strlen(instrument) + strlen(name) + strlen(camera) + strlen(tp) + strlen(longitude) + 5, 5);
		
		if(read_first < 0){
			curr_long_size = dim_sum(af_read_size(file, long_dataset_name), 1);
			long_data = af_read(file, long_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and its dimension
			double* adding_long = af_read(file, long_dataset_name);
			double new_long_size = dim_sum(af_read_size(file, long_dataset_name), 1);
			//Reallocating arrays of data
			long_data = realloc(long_data, sizeof(double)*(curr_long_size + new_long_size));
			memcpy(&long_data[(int)curr_long_size], adding_long, sizeof(double)*new_long_size);
			curr_long_size += new_long_size;

			free(adding_long);
		}
	}
	*size = curr_long_size;
	
	if(long_data != NULL){
	printf("test_long_data: %f\n", long_data[0]);
	printf("test_long_data: %f\n", long_data[1]);
	printf("test_long_data: %f\n", long_data[2]);
	}
	
	return long_data;
}


/*
						get_mop_rad
	DESCRIPTION:	
		This function retrieves MOPITT radiance data. All the granules are stitched together and returned as a 1D array
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. size -- An integer pointer that points to the size of the MOPITT radiance data after the retreival 
		
	EFFECT:
		Memory would be allocated for retrieved MOPITT radiance data with the size set as the data size
		
	RETURN:
		Returns data upon successful retrieval
		Returns NULL upon  error
		
*/


double* get_mop_rad(hid_t file, int* size){
	printf("Reading MOPITT radiance\n");
	//Path variables
	char* instrument = "MOPITT";
	char* d_field = "Data_Fields";
	char* rad = "MOPITTRadiances";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(names[i], name);
		free(name);
	}
	int h;
	double* data;
	double curr_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {instrument, name, d_field, rad};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(d_field) + strlen(rad) + 4, 4);
		printf("granule_name: %s\n", name);
		if(read_first < 0){
			data = af_read(file, dataset_name);
			if(data == NULL){
				continue;
			}
			curr_size = dim_sum(af_read_size(file, dataset_name), 5); 
			read_first = 1;
		}
		else{
			//retrieve next set of data and its dimension
			double* adding_data = af_read(file, dataset_name);
			if(adding_data == NULL){
				continue;
			}
			double new_d_size = dim_sum(af_read_size(file, dataset_name), 1);
			//Reallocating arrays of data
			data = realloc(data, sizeof(double)*(curr_size + new_d_size));
			memcpy(&data[(int)curr_size], adding_data, sizeof(double)*new_d_size);
			curr_size += new_d_size;
			
			free(adding_data);
		}
	}
	*size = curr_size;
	
	//Print statements to verify data's existence	
	if(data != NULL){
		printf("test data: %f\n", data[0]);
		printf("test_data: %f\n", data[1]);
		printf("test data: %f\n", data[2]);
	}	
	return data;
}


/*
						get_mop_lat
	DESCRIPTION:	
		This function retrieves MOPITT geological latitude data. All the granules are stitched together and returned as a 1D array
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. size -- An integer pointer that points to the size of the MOPITT geological latitude data after the retreival 
		
	EFFECT:
		Memory would be allocated for retrieved MOPITT latitude data with the size set as the data size
		
	RETURN:
		Returns lat_data upon successful retrieval
		Returns NULL upon  error
		
*/



double* get_mop_lat(hid_t file, int* size){
	printf("Reading MOPITT lat\n");
	//Path variables
	char* instrument = "MOPITT";
	char* d_field = "Data_Fields";
	char* rad = "MOPITTRadiances";
	char* location = "Geolocation";
	char* lat = "Latitude";
	
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(names[i], name);
		free(name);
	}
	
	int h;
	double* lat_data;
	double curr_lat_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {name, d_field, rad};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(name) + strlen(d_field) + strlen(rad), 3);
		memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
		//Check if dataset exists first
		printf("granule_name: %s\n", name);
		htri_t status = H5Lexists(group, dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		const char* lat_arr[] = {instrument, name, location, lat};
		char* lat_dataset_name;
		concat_by_sep(&lat_dataset_name, lat_arr, "/", strlen(instrument) + strlen(name) + strlen(location) + strlen(lat) + 4, 4);
		
		if(read_first < 0){
			curr_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 3);
			lat_data = af_read(file, lat_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and its dimension
			double* adding_lat = af_read(file, lat_dataset_name);
			double new_lat_size = dim_sum(af_read_size(file, lat_dataset_name), 3);
			//Reallocating arrays of data
			lat_data = realloc(lat_data, sizeof(double)*(curr_lat_size + new_lat_size));
			memcpy(&lat_data[(int)curr_lat_size], adding_lat, sizeof(double)*new_lat_size);
			curr_lat_size += new_lat_size;

			free(adding_lat);
		}
	}
	*size = curr_lat_size;
	//Print statements to verify data's existence
	if(lat_data != NULL){
		printf("test_lat_data: %f\n", lat_data[0]);
		printf("test_lat_data: %f\n", lat_data[1]);
		printf("test_lat_data: %f\n", lat_data[2]);
	}
	
	return lat_data;
}


/*
						get_mop_long
	DESCRIPTION:	
		This function retrieves MOPITT geological longitude data. All the granules are stitched together and returned as a 1D array
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. size -- An integer pointer that points to the size of the MOPITT geological longitude data after the retreival 
		
	EFFECT:
		Memory would be allocated for retrieved MOPITT longitude data with the size set as the data size
		
	RETURN:
		Returns long_data upon successful retrieval
		Returns NULL upon  error
		
*/


double* get_mop_long(hid_t file, int* size){
	printf("Reading MOPITT longitude\n");
	//Path variables
	char* instrument = "MOPITT";
	char* d_field = "Data_Fields";
	char* rad = "MOPITTRadiances";
	char* location = "Geolocation";
	char* longitude = "Longitude";
	
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(names[i], name);
		free(name);
	}
	
	int h;
	double* long_data;
	double curr_long_size;
	int read_first = -1;
	for(h = 0; h < num_groups; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {name, d_field, rad};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(name) + strlen(d_field) +strlen(rad) + 3, 3);
		memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
		//Check if dataset exists first
		printf("granule_name: %s\n", name);
		htri_t status = H5Lexists(group, dataset_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		
		const char* long_arr[] = {instrument, name, location, longitude};
		char* long_dataset_name;
		concat_by_sep(&long_dataset_name, long_arr, "/", strlen(instrument) + strlen(name) + strlen(location) + strlen(longitude) + 4, 4);
		
		if(read_first < 0){
			curr_long_size = dim_sum(af_read_size(file, long_dataset_name), 3);
			long_data = af_read(file, long_dataset_name);
			read_first = 1;
		}
		else{
			//retrieve next set of data and its dimension
			double* adding_long = af_read(file, long_dataset_name);
			double new_long_size = dim_sum(af_read_size(file, long_dataset_name), 3);
			//Reallocating arrays of data
			long_data = realloc(long_data, sizeof(double)*(curr_long_size + new_long_size));
			memcpy(&long_data[(int)curr_long_size], adding_long, sizeof(double)*new_long_size);
			curr_long_size += new_long_size;

			free(adding_long);
		}
	}
	*size = curr_long_size;
	//Print statements to verify data's existence
	if(long_data != NULL){
		printf("test_long_data: %f\n", long_data[0]);
		printf("test_long_data: %f\n", long_data[1]);
		printf("test_long_data: %f\n", long_data[2]);
	}
	return long_data;
}


/*
						get_ast_rad
	DESCRIPTION:	
		This functions retrieves ASTER radiance data based on specified parameters including subsystem and dataset name. At the moment, the granules
		are stitched together and returned as a 1D array, which is WRONG because each granule has different dimensions. Additional work has to be done
		to account for those differenes (rotations).
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. subsystem (TIR/VNIR/SWIR) -- A string variable that specifies the subsystem to retrieve
		2. d_name -- A string variable that specifies the dataset name 
		3. size -- An integer pointer that points to the size of the ASTER radiance data after the retreival 
		
	EFFECT:
		Memory would be allocated for the retrieved ASTER radiance data, with the variable size set as the size of the data
		
	RETURN:
		Returns result_data upon successful retrieval
		Returns NULL upon error
*/


double* get_ast_rad(hid_t file, char* subsystem, char* d_name, int*size){
	printf("Reading ASTER radiance\n");
	//Path variables
	char* instrument = "ASTER";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		strcpy(names[i], name);
		free(name);
	}
	
	//Get total data size
	printf("Get total data size\n");
	int k;
	int total_size = 0;
	int store_count = 0;
	for(k = 0; k < num_groups; k++){
		char* name = names[k];
		const char* d_arr[] = {instrument, name, subsystem, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(subsystem) + strlen(d_name), 4);
		hsize_t* curr_dim = af_read_size(file, dataset_name);
		if(curr_dim == NULL){
			continue;
		}
		store_count += 1;
		total_size += curr_dim[0]*curr_dim[1];
		free(curr_dim);
	}
	
	double* result_data = calloc(total_size, sizeof(double));
	
	int h;
	int curr_size = 0;
	for(h = 0; h < store_count; h++){
		//Path formation
		char* name = names[h];
		const char* d_arr[] = {instrument, name, subsystem, d_name};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(subsystem) + strlen(d_name) + 5, 4);
		printf("dataset name: %s\n", dataset_name);
		double* data = af_read(file, dataset_name);
		if(data == NULL){
				continue;
		}
		hsize_t* curr_dim = af_read_size(file, dataset_name);
		int gran_size = curr_dim[0] * curr_dim[1];
		memcpy(&result_data[curr_size], data, sizeof(double) * gran_size);
		curr_size += gran_size;
		free(data);
		free(curr_dim);
	}
	*size = curr_size;
	
	//Print statements to verify data's existence
	if(result_data != NULL){
		printf("test data: %f\n", result_data[0]);
		printf("test_data: %f\n", result_data[1]);
		printf("test data: %f\n", result_data[2]);
	}	
	
	return result_data;
}

/*
						get_ast_lat
	DESCRIPTION:	
		This functions retrieves ASTER geological latitude data based on specified parameters including subsystem and dataset name. At the moment, the granules
		are stitched together and returned as a 1D array, which is WRONG because each granule has different dimensions. Additional work has to be done
		to account for those differenes (rotations).
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. subsystem (TIR/VNIR/SWIR) -- A string variable that specifies the subsystem to retrieve
		2. d_name -- A string variable that specifies the dataset name 
		3. size -- An integer pointer that points to the size of the ASTER latitude data after the retreival 
		
	EFFECT:
		Memory would be allocated for the retrieved ASTER geological latitude data, with the variable size set as the size of the data
		
	RETURN:
		Returns lat_data upon successful retrieval
		Returns NULL upon error
*/

double* get_ast_lat(hid_t file, char* subsystem, char* d_name, int*size){
	printf("Reading ASTER lat\n");
	//Path variables
	char* instrument = "ASTER";
	char* location = "Geolocation";
	char* lat = "Latitude";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		char* rad_group_name;
		const char* d_arr[] = {name, subsystem, d_name};
		concat_by_sep(&rad_group_name, d_arr, "/", strlen(name) + strlen(subsystem) + strlen(d_name), 3);
		memmove(&rad_group_name[0], &rad_group_name[1], strlen(rad_group_name));
		htri_t status = H5Lexists(group, rad_group_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		strcpy(names[i], name);
		free(name);
	}
	
	//Get total data size
	printf("Get total data size\n");
	int k;
	int total_size = 0;
	int store_count = 0;
	for(k = 0; k < num_groups; k++){
		char* name = names[k];
		const char* d_arr[] = {instrument, name, subsystem, location, lat};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(subsystem) + strlen(location) + strlen(lat), 5);
		hsize_t* curr_dim = af_read_size(file, dataset_name);
		if(curr_dim == NULL){
			continue;
		}
		store_count += 1;
		total_size += curr_dim[0]*curr_dim[1];
		free(curr_dim);
	}
	
	int h;
	double* lat_data = calloc(total_size, sizeof(double));
	int curr_lat_size = 0;
	int read_first = -1;
	for(h = 0; h < store_count; h++){
		//Path formation
		char* name = names[h];
		const char* lat_arr[] = {instrument, name, subsystem, location, lat};
		char* lat_dataset_name;
		concat_by_sep(&lat_dataset_name, lat_arr, "/", strlen(instrument) + strlen(name) + strlen(subsystem) + strlen(location) + strlen(lat) + 5, 5);
	    printf("dataset name: %s\n", lat_dataset_name);	
		double* data = af_read(file, lat_dataset_name);
		if(data == NULL){
				continue;
		}
		hsize_t* curr_dim = af_read_size(file, lat_dataset_name);
		int gran_size = curr_dim[0] * curr_dim[1];
		memcpy(&lat_data[curr_lat_size], data, sizeof(double) * gran_size);
		curr_lat_size += gran_size;
		free(data);
		free(curr_dim);
	}
	*size = curr_lat_size;
	//Print statements to verify data's existence
	if(lat_data != NULL){	
		printf("test_lat_data: %f\n", lat_data[0]);
		printf("test_lat_data: %f\n", lat_data[1]);
		printf("test_lat_data: %f\n", lat_data[2]);
	}
	return lat_data;
}

/*
						get_ast_long
	DESCRIPTION:	
		This functions retrieves ASTER geological longitude data based on specified parameters including subsystem and dataset name. At the moment, the granules
		are stitched together and returned as a 1D array, which is WRONG because each granule has different dimensions. Additional work has to be done
		to account for those differenes (rotations).
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. subsystem (TIR/VNIR/SWIR) -- A string variable that specifies the subsystem to retrieve
		2. d_name -- A string variable that specifies the dataset name 
		3. size -- An integer pointer that points to the size of the ASTER longitude data after the retreival 
		
	EFFECT:
		Memory would be allocated for the retrieved ASTER geological longitude data, with the variable size set as the size of the data
		
	RETURN:
		Returns long_data upon successful retrieval
		Returns NULL upon error
*/


double* get_ast_long(hid_t file, char* subsystem, char* d_name, int* size){
	printf("Reading ASTER long\n");
	//Path variables
	char* instrument = "ASTER";
	char* location = "Geolocation";
	char* longitude = "Longitude";
	//Get all granule file names
	printf("Retrieving granule group names\n");
	hid_t group = H5Gopen(file, instrument, H5P_DEFAULT);
	if(group < 0){
		printf("Group not found\n");
		return NULL;
	}
	hsize_t num_groups;
	herr_t err = H5Gget_num_objs(group, &num_groups);
	char names[(int)num_groups][50];
	int i;
	for(i = 0; i < num_groups; i++){
		char* name = malloc(50*sizeof(char));
		H5Gget_objname_by_idx(group, (hsize_t)i, name, 50);
		char* rad_group_name;
		const char* d_arr[] = {name, subsystem, d_name};
		concat_by_sep(&rad_group_name, d_arr, "/", strlen(name) + strlen(subsystem) + strlen(d_name), 3);
		memmove(&rad_group_name[0], &rad_group_name[1], strlen(rad_group_name));
		htri_t status = H5Lexists(group, rad_group_name, H5P_DEFAULT);
		if(status <= 0){
			printf("Dataset does not exist\n");
			continue;
		}
		strcpy(names[i], name);
		free(name);
	}
	
	//Get total data size
	printf("Get total data size\n");
	int k;
	int total_size = 0;
	int store_count = 0;
	for(k = 0; k < num_groups; k++){
		char* name = names[k];
		const char* d_arr[] = {instrument, name, subsystem, location, longitude};
		char* dataset_name;
		concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(name) + strlen(subsystem) + strlen(location) + strlen(longitude), 5);
		hsize_t* curr_dim = af_read_size(file, dataset_name);
		if(curr_dim == NULL){
			continue;
		}
		store_count += 1;
		total_size += curr_dim[0]*curr_dim[1];
		free(curr_dim);
	}
	
	int h;
	double* long_data = calloc(total_size, sizeof(double));
	int curr_long_size = 0;
	for(h = 0; h < store_count; h++){
		//Path formation
		char* name = names[h];
		const char* long_arr[] = {instrument, name, subsystem, location, longitude};
		char* long_dataset_name;
		concat_by_sep(&long_dataset_name, long_arr, "/", strlen(instrument) + strlen(name) + strlen(subsystem) + strlen(location) + strlen(longitude) + 5, 5);
		printf("dataset name: %s\n", long_dataset_name);	
		double* data = af_read(file, long_dataset_name);
		if(data == NULL){
				continue;
		}
		hsize_t* curr_dim = af_read_size(file, long_dataset_name);
		int gran_size = curr_dim[0] * curr_dim[1];
		memcpy(&long_data[curr_long_size], data, sizeof(double) * gran_size);
		curr_long_size += gran_size;
		free(data);
		free(curr_dim);
	}
	*size = curr_long_size;
	//Print statements to verify data's existence
	if(long_data != NULL){
		printf("test_long_data: %f\n", long_data[0]);
		printf("test_long_data: %f\n", long_data[1]);
		printf("test_long_data: %f\n", long_data[2]);
	}
	
	return long_data;
}


/*
						get_ast_rad_by_gran
	DESCRIPTION:	
		This functions retrieves ASTER radiance data based on specified parameters including subsystem, dataset name and granule name. 
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. subsystem (TIR/VNIR/SWIR) -- A string variable that specifies the subsystem to retrieve
		2. d_name -- A string variable that specifies the dataset name 
		3. gran_name -- A string variable that specifies te granule name
		4. size -- An integer pointer that points to the size of the ASTER radiance data after the retreival 
		
	EFFECT:
		Memory would be allocated for the retrieved ASTER radiance data, with the variable size set as the size of the data
		
	RETURN:
		Returns result_data upon successful retrieval
		Returns NULL upon error
*/

double* get_ast_rad_by_gran(hid_t file, char* subsystem, char* d_name, char* gran_name, int*size){
	printf("Reading ASTER radiance by gran\n");
	//Path variables
	char* instrument = "ASTER";
	char* dataset_name;
	const char* d_arr[] = {instrument, gran_name, subsystem, d_name};
	concat_by_sep(&dataset_name, d_arr, "/", strlen(instrument) + strlen(gran_name) + strlen(subsystem) + strlen(d_name), 4);
	memmove(&dataset_name[0], &dataset_name[1], strlen(dataset_name));
	htri_t status = H5Lexists(file, dataset_name, H5P_DEFAULT);
	if(status <= 0){
		printf("Dataset does not exist\n");
		return NULL;
	}
	hsize_t* curr_dim = af_read_size(file, dataset_name);
	if(curr_dim == NULL){
		printf("Retrieve current dimension error\n");
	}
	*size = curr_dim[0] * curr_dim[1];
	double* result_data = calloc(curr_dim[0]*curr_dim[1], sizeof(double));
	double* data = af_read(file, dataset_name);
	if(data == NULL){
		printf("Read error\n");
		return NULL;
	}
	memcpy(&result_data[0], data, sizeof(double)*(curr_dim[0]*curr_dim[1]));
	free(data);
	free(curr_dim);
	
	if(result_data != NULL){
		printf("print test data\n");
		printf("test: %f\n", result_data[0]);
	}
	return result_data;
}


/*
						get_ast_lat_by_gran
	DESCRIPTION:	
		This functions retrieves ASTER latitude data based on specified parameters including subsystem, dataset name and granule name. 
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. subsystem (TIR/VNIR/SWIR) -- A string variable that specifies the subsystem to retrieve
		2. d_name -- A string variable that specifies the dataset name 
		3. gran_name -- A string variable that specifies te granule name
		4. size -- An integer pointer that points to the size of the ASTER latitude data after the retreival 
		
	EFFECT:
		Memory would be allocated for the retrieved ASTER latitude data, with the variable size set as the size of the data
		
	RETURN:
		Returns result_data upon successful retrieval
		Returns NULL upon error
*/

double* get_ast_lat_by_gran(hid_t file, char* subsystem, char* d_name, char* gran_name, int*size){
	printf("Reading ASTER lat by gran\n");
	//Path variables
	char* instrument = "ASTER";
	char* location = "Geolocation";
	char* lat = "Latitude";
	
	const char* lat_arr[] = {instrument, gran_name, subsystem, location, lat};
	char* lat_dataset_name;
	concat_by_sep(&lat_dataset_name, lat_arr, "/", strlen(instrument) + strlen(gran_name) + strlen(subsystem) + strlen(location) + strlen(lat) + 5, 5);
	printf("dataset name: %s\n", lat_dataset_name);	
	hsize_t* curr_dim = af_read_size(file, lat_dataset_name);
	*size = curr_dim[0] * curr_dim[1];
	double* result_data = calloc(curr_dim[0]*curr_dim[1], sizeof(double));
	double* data = af_read(file, lat_dataset_name);
	if(data == NULL){
		return NULL;
	}
	memcpy(&result_data[0], data, sizeof(double)*(curr_dim[0]*curr_dim[1]));
	free(data);
	free(curr_dim);
	
	if(result_data != NULL){
		printf("print test data\n");
		printf("test: %f\n", result_data[0]);
	}
	return result_data;
}

/*
						get_ast_long_by_gran
	DESCRIPTION:	
		This functions retrieves ASTER latitude data based on specified parameters including subsystem, dataset name and granule name. 
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. subsystem (TIR/VNIR/SWIR) -- A string variable that specifies the subsystem to retrieve
		2. d_name -- A string variable that specifies the dataset name 
		3. gran_name -- A string variable that specifies te granule name
		4. size -- An integer pointer that points to the size of the ASTER longitude data after the retreival 
		
	EFFECT:
		Memory would be allocated for the retrieved ASTER longitude data, with the variable size set as the size of the data
		
	RETURN:
		Returns result_data upon successful retrieval
		Returns NULL upon error
*/


double* get_ast_long_by_gran(hid_t file, char* subsystem, char* d_name, char* gran_name, int*size){
	printf("Reading ASTER long by gran\n");
	//Path variables
	char* instrument = "ASTER";
	char* location = "Geolocation";
	char* longitude = "Longitude";
	
	const char* long_arr[] = {instrument, gran_name, subsystem, location, longitude};
	char* long_dataset_name;
	concat_by_sep(&long_dataset_name, long_arr, "/", strlen(instrument) + strlen(gran_name) + strlen(subsystem) + strlen(location) + strlen(longitude) + 5, 5);
	printf("dataset name: %s\n", long_dataset_name);	
	hsize_t* curr_dim = af_read_size(file, long_dataset_name);
	*size = curr_dim[0] * curr_dim[1];
	double* result_data = calloc(curr_dim[0]*curr_dim[1], sizeof(double));
	double* data = af_read(file, long_dataset_name);
	if(data == NULL){
			return NULL;
	}
	memcpy(&result_data[0], data, sizeof(double)*(curr_dim[0]*curr_dim[1]));
	free(data);
	free(curr_dim);
	
	if(result_data != NULL){
		printf("print test data\n");
		printf("test: %f\n", result_data[0]);
	}
	return result_data;
}


/*
						af_read_size
	DESCRIPTION:	
		A HDF5 API wrapper for advancedFusion to read the dimension of a dataset.
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. dataset_name -- A string variable that specifies the dataset name
		
	EFFECT:
		Memory would be allocated for the hsize_t pointer that holds the dimension of the dataset
			
	RETURN:
		Returns dims if the dataset is found and opened correctly
		Returns NULL upon error
		
*/


hsize_t* af_read_size(hid_t file, char* dataset_name){
	hid_t dataset = H5Dopen2(file, dataset_name, H5P_DEFAULT);
	if(dataset < 0){
		printf("Dataset open error\n");
		return NULL; 
	}
	hid_t dataspace = H5Dget_space(dataset);
	if(dataspace < 0){
		printf("Dataspace open error\n");
		return NULL;	
	}
	const int ndims = H5Sget_simple_extent_ndims(dataspace);
	hsize_t* dims = malloc(sizeof(hsize_t) * ndims);
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	H5Dclose(dataset);	
	H5Sclose(dataspace);
	return dims;
}

/*
						af_read
	DESCRIPTION:	
		A HDF5 API wrapper for advancedFusion to read a dataset based on specified dataset name. A special case for ASTER is handled in the function
		because the data for ASTER is of 64-bit floating point numbers, unlike the rest (32-bit floating point numbers). All data are coverted to type
		of double. 
		
	ARGUMENTS:
		0. file -- A hdf file variable that points to the BasicFusion file
		1. dataset_name -- A string variable that specifies the dataset name, which should be the full path within the BasicFusion file
		
	EFFECT:
		Memory would be allocated for data that is read in.
		
	RETURN:
		Returns data/converted_data upon successful retrieval
		Returns NULL upon error
		
*/


double* af_read(hid_t file, char* dataset_name){
	hid_t dataset = H5Dopen2(file, dataset_name, H5P_DEFAULT);
	if(dataset < 0){
		printf("Dataset open error\n");
		return NULL; 
	}
	hid_t dataspace = H5Dget_space(dataset);
	if(dataspace < 0){
		printf("Dataspace open error\n");
		return NULL;	
	}
	
	const int ndims = H5Sget_simple_extent_ndims(dataspace);
	hsize_t dims[ndims];
	H5Sget_simple_extent_dims(dataspace, dims, NULL);
	hid_t memspace = H5Screate_simple(ndims,dims,NULL);
	hid_t dtype = H5Dget_type(dataset);
	hid_t ndtype = H5Tget_native_type(dtype, H5T_DIR_DESCEND);
	if(strstr(dataset_name, "ASTER") != NULL && strstr(dataset_name, "Geolocation") != NULL){
		//Special case for ASTER geolocation because they are 64bit floating point numbers
		double* data = calloc ( dim_sum(dims, sizeof(dims)/sizeof(hsize_t)) , sizeof(double) );
		herr_t status = H5Dread(dataset, ndtype, memspace, memspace, H5P_DEFAULT, data);
		H5Dclose(dataset);	
		H5Sclose(dataspace);
		H5Tclose(dtype);
		H5Tclose(ndtype);
		if(status < 0){
			printf("read error: %d\n", status);
		}
		return data;
	}
	else{
		float* data = calloc ( dim_sum(dims, sizeof(dims)/sizeof(hsize_t)) , sizeof(ndtype) );
		double* converted_data = calloc ( dim_sum(dims, sizeof(dims)/sizeof(hsize_t)) , sizeof(double) );
		herr_t status = H5Dread(dataset, ndtype, memspace, memspace, H5P_DEFAULT, data);
		int i;
		for(i=0;i < dim_sum(dims, sizeof(dims)/sizeof(hsize_t)); i++){
			converted_data[i] = (double) data[i];
		}
		free(data);
		H5Dclose(dataset);	
		H5Sclose(dataspace);
		H5Tclose(dtype);
		H5Tclose(ndtype);
		if(status < 0){
			printf("read error: %d\n", status);
		}
		return converted_data;
	}
}

/*
						af_write_misr_on_modis
	DESCRIPTION:	
		This function writes the resultant data values to a designated output HDF5 file after reprojecting MISR data to MODIS grid. The function is
		created for demo purpose, and due to limitation of time. Ideally, a generalized write function should be written for all kinds of reprojection.
		
	ARGUMENTS:
		0. output_file -- A file pointer that points the designated HDF5 file that holds the result data
		1. misr_out -- MISR's data after reprojection
		2. modis -- MODIS radiance data
		3. modis_size -- Total size of MODIS radiance data
		4. modis_band_size -- The total number of bands for MODIS
		5. misr_size -- Total size of MISR data
		
	EFFECT:
		The reprojected result would be written to output_file
		
	RETURN:
		Returns 1 upon successful writing
		Returns -1 upon error
*/


int af_write_misr_on_modis(hid_t output_file, double* misr_out, double* modis, int modis_size, int modis_band_size, int misr_size){
	//Create datafield group
	hid_t group_id = H5Gcreate2(output_file, "/Data_Fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	//Write MODIS first
	hsize_t modis_dim[3];
	modis_dim[0] = modis_band_size;
	modis_dim[2] = 1354;
	modis_dim[1] = (modis_size)/modis_band_size/1354;
	hid_t modis_dataspace = H5Screate_simple(3, modis_dim, NULL);
	hid_t	modis_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    herr_t  modis_status = H5Tset_order(modis_datatype, H5T_ORDER_LE);  
    hid_t modis_dataset = H5Dcreate2(output_file, "/Data_Fields/modis_rad", modis_datatype, modis_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    modis_status = H5Dwrite(modis_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, modis);
    H5Sclose(modis_dataspace);
	H5Tclose(modis_datatype);
	H5Dclose(modis_dataset);
    if(modis_status < 0){
    	printf("MODIS write error\n");
    	return -1;
	}
    
    //Write 
    hsize_t misr_dim[2];
	misr_dim[0] = (misr_size) / 1354;
	misr_dim[1] = 1354;
	hid_t misr_dataspace = H5Screate_simple(2, misr_dim, NULL);
	hid_t misr_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    herr_t misr_status = H5Tset_order(misr_datatype, H5T_ORDER_LE);  
    hid_t misr_dataset = H5Dcreate2(output_file, "/Data_Fields/misr_out", misr_datatype, misr_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    misr_status = H5Dwrite(misr_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, misr_out);
    H5Sclose(misr_dataspace);
	H5Tclose(misr_datatype);
	H5Dclose(misr_dataset);
	if(misr_status < 0){
		printf("MISR write error\n");
		return -1;
	}
	
	return 1;
	
}


/*
						af_write_mm_geo
	DESCRIPTION:	
		This function writes the geological data (latitude and logitude) to the designated output HDF5 after the reprojection from MISR to MODIS.
		Again, same as af_write_misr_on_modis, a generalized version of write function is advised in the future.
		
	ARGUMENTS:
		0. output_file -- A file pointer that points the designated HDF5 file that holds the result data
		1. geo_flag (0/1) -- An integer variable that specifies whether the data is latitude or longitude
		2. geo_data -- The geolocation data to be written
		3. geo_size -- Total size of the geolocation data
		
	EFFECT:
		The geolocation data would be written to the output_file.
		
	RETURN:
		Returns 1 upon successful writing
		Returns -1 upon error
		
*/


int af_write_mm_geo(hid_t output_file, int geo_flag, double* geo_data, int geo_size){
	//Check if geolocation group exists --- TODO - change it to H5Lexists
	printf("test: %f\n", geo_data[0]);
	printf("test: %f\n", geo_data[1]);
	htri_t status = H5Lexists(output_file, "Geolocation", H5P_DEFAULT);
	if(status <= 0){
		hid_t group_id = H5Gcreate2(output_file, "/Geolocation", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	char* d_name;
	if(geo_flag == 0){
		d_name = "/Geolocation/Latitude";
	}
	else if(geo_flag == 1){
		d_name = "/Geolocation/Longitude";
	}
	hsize_t     geo_dim[2];
	geo_dim[0] = (geo_size) / 1354;
	geo_dim[1] = 1354;
	hid_t geo_dataspace = H5Screate_simple(2, geo_dim, NULL);
	hid_t geo_datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    herr_t geo_status = H5Tset_order(geo_datatype, H5T_ORDER_LE);  
    hid_t geo_dataset = H5Dcreate2(output_file, d_name, geo_datatype, geo_dataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(geo_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, geo_data);
    H5Sclose(geo_dataspace);
	H5Tclose(geo_datatype);
	H5Dclose(geo_dataset);
	
	if(geo_status < 0){
		printf("Geo Data write error\n");
		return -1;
	}
	return 1;
}

/*
						af_open
	DESCRIPTION:	
		This is a HDF5 API wrapper for advancedFusion to open a HDF5 file.
		
	ARGUMENTS:
		0. file_path -- A string variable that specifies the absolute file path to the HDF5 file
		
	EFFECT:
		A HDF5 file would be opened with its file pointer returned
		
	RETURN:
		Returns f ( if f is smaller than 0, an error has occured)
		
*/


hid_t af_open(char* file_path){
	hid_t f = H5Fopen(file_path, H5F_ACC_RDONLY, H5P_DEFAULT);
	return f;
}


/*
						af_close
	DESCRIPTION:	
		This is a HDF5 API wrapper for advancedFusion to close a HDF5 file.
		
	ARGUMENTS:
		0. file -- An identifier of the opened file
		
	EFFECT:
		A HDF5 file would be closed
		
	RETURN:
		Returns ret ( if ret is smaller than 0, an error has occured)
		
*/


herr_t af_close(hid_t file){
	herr_t ret = H5Fclose(file);
	return ret;
}


//Deprecated -- the original main function for reading data, check af_run.c for the latest version
/*int main (int argc, char *argv[]){
	//Preset filename here for easy testing 
	char* file_path = "/projects/TDataFus/kent/temp/40-orbit-file/Jun15.2/TERRA_BF_L1B_O69365_F000_V000.h5";
	hid_t file;
	argv[1] = file_path;
	if(argc < 3){
		printf("Usage: %s filename instrument_name extra_arguments\n", argv[0]);
		return 0;
	}
	else if(strcmp(argv[2], "MISR") == 0){
		if(argc < 6){
			printf("MISR Usage: %s filename MISR camera_angle resolution(H/L) radiance\n", argv[0]);
			return 0; 
		}
		else{
			//MISR input requirements fulfilled 
			//Open file
			file = af_open(file_path);
			if(file < 0){
				printf("File not found\n");
				return -1;
			}
						
			int d_size = 0;
			int* data_pt = &d_size;
			double* data = get_misr_rad(file, argv[3], argv[4], argv[5], data_pt);
			printf("Data size: %d\n", *data_pt);
			
			int lat_size = 0;
			int* lat_pt = &lat_size;
			double* lat_data = get_misr_lat(file, argv[4], lat_pt);
			printf("Lat size: %d\n", *lat_pt);
			
			int long_size = 0;
			int* long_pt = &long_size;
			double* long_data = get_misr_long(file, argv[4], long_pt);
			printf("Long size: %d\n", *long_pt);
			
			if(data != NULL && lat_data != NULL, long_data != NULL){
				printf("MISR Data retrieval successful\n");
			}
			else{
				printf("MISR Data retrieval failed\n");
			}
			
			printf("Retrieving MISR attributes\n");
			void* unit_attr;
			unit_attr = get_misr_attr(file, argv[3], argv[4], argv[5], "Units", 0, unit_attr);
			printf("Unit attr: %s\n", (char*)unit_attr);
			void* fill_attr;
			fill_attr = get_misr_attr(file, argv[3], argv[4], argv[5], "_FillValue", 0, fill_attr);
			printf("FillValue: %f\n", *(float*) fill_attr);
			void* lat_attr;
			lat_attr = get_misr_attr(file, argv[3], argv[4], argv[5], "units", 1, lat_attr);
			printf("lat_units: %s\n", (char*)lat_attr);
			void* long_attr;
			long_attr = get_misr_attr(file, argv[3], argv[4], argv[5], "units", 2, long_attr);
			printf("long_units: %s\n", (char*)long_attr);
			
			herr_t ret = af_close(file);
		}
	}
	else if(strcmp(argv[2], "MODIS") == 0){
		if(argc < 4){
			printf("MODIS Usage: %s filename MODIS resolution(1KM/500m/250m)\n", argv[0]);
		}
		else{
			file = af_open(file_path);
			if(file < 0){
				printf("File not found\n");
				return -1;
			}
			char* resolution = argv[3];
			char* d_name = "";
			if(strcmp(resolution, "1KM") == 0){
				resolution = "_1KM";
				d_name = "EV_1KM_RefSB";
			}
			else if(strcmp(resolution, "250M") == 0){
				resolution = "_250m";
				d_name = "EV_250_RefSB";
			}
			else if(strcmp(resolution, "500M") == 0){
				resolution = "_500m";
				d_name = "EV_500_RefSB";
			}
			else{
				printf("Wrong resolution, choose from 1KM, 500M or 250M\n");
			}
			
			int data_size = 0;
			int* data_pt = &data_size;
			double* data = get_modis_rad(file, resolution, d_name, data_pt);
			printf("Data size: %d\n", *data_pt);
			
			int lat_size = 0;
			int* lat_pt = &lat_size;
			double* lat_data = get_modis_lat(file, resolution, d_name, lat_pt);
			printf("Lat size: %d\n", *lat_pt);
			
			int long_size = 0;
			int* long_pt = &long_size;
			double* long_data = get_modis_long(file, resolution, d_name, long_pt);
			printf("Long size: %d\n", *long_pt);

			
			if(data != NULL && lat_data != NULL, long_data != NULL){
				printf("MODIS retrieval successful\n");
			}
			else{
				printf("MODIS retrieval failed\n");
			}
			
		
			printf("Retrieving MODIS attributes\n");
			void* unit_attr;
			unit_attr = get_modis_attr(file, resolution, d_name, "units", 0, unit_attr);
			void* fill_attr;
			fill_attr = get_modis_attr(file, resolution, d_name, "_FillValue", 0, fill_attr);
			void* min_attr;
			min_attr = get_modis_attr(file, resolution, d_name, "valid_min", 0, min_attr);
			void* lat_attr;
			lat_attr = get_modis_attr(file, resolution, d_name, "units", 1, lat_attr);
			void* long_attr;
			long_attr = get_modis_attr(file, resolution, d_name, "units", 2, long_attr);
			printf("Unit attr: %s\n", (char*)unit_attr);
			printf("FillValue: %f\n", *(float*) fill_attr);
			printf("valid_min: %f\n", *(float*) min_attr);
			printf("lat unit: %s\n", (char*) lat_attr);
			printf("long unit: %s\n", (char*) long_attr); 
			
			herr_t ret = af_close(file);
		}
	}
	else if(strcmp(argv[2], "CERES") == 0){
		if(argc < 6){
			printf("CERES Usage: %s filename CERES camera radiance(LW/SW/WN/TOT) filtered/unfiltered(F/U)");
		}
		else{
				file = af_open(file_path);
				if(file < 0){
					printf("File not found\n");
					return -1;
				}
				char* d_name = calloc(30, sizeof(char));
				if(strcmp(argv[5], "F") == 0){
					char* f = "_Filtered";
					char* r = "_Radiance";
					strcpy(d_name, argv[4]);
					strncat(d_name, f, strlen(f));
					strncat(d_name, r, strlen(r));
				}
				else{
					char* r = "_Radiance";
					strcpy(d_name, argv[4]);
					strncat(d_name, r, strlen(r));
				}
				
				int data_size = 0;
				int* data_pt = &data_size;
				double* data = get_ceres_rad(file, argv[3], d_name, data_pt);
				printf("Data size: %d\n", *data_pt);
				
				int lat_size = 0;
				int* lat_pt = &lat_size;
				double* lat_data = get_ceres_lat(file, argv[3], d_name, lat_pt);
				printf("Lat size: %d\n", *lat_pt);

				int long_size = 0;
				int* long_pt = &long_size;
				double* long_data = get_ceres_long(file, argv[3], d_name, long_pt);
				printf("Long size: %d\n", *long_pt);
				
				herr_t ret = af_close(file);
		}
	}
	else if(strcmp(argv[2], "MOPITT") == 0){
		file = af_open(file_path);
		if(file < 0){
			printf("File not found\n");
			return -1;
		}
		int data_size = 0;
		int* data_pt = &data_size;
		double* data = get_mop_rad(file, data_pt);
		printf("Data size: %d\n", *data_pt);
		
		int lat_size = 0;
		int* lat_pt = &lat_size;
		double* lat_data = get_mop_lat(file, lat_pt);
		printf("Lat size: %d\n", *lat_pt);

		int long_size = 0;
		int* long_pt = &long_size;
		double* long_data = get_mop_long(file, long_pt);
		printf("Long size: %d\n", *long_pt);
		
		herr_t ret = af_close(file);
	}
	else if(strcmp(argv[2], "ASTER") == 0){
		if(argc < 5){
			printf("ASTER Usage: %s filename ASTER subsystem(TIR/VNIR/SWIR) dataset_name\n");
		}
		else{
				file = af_open(file_path);
				if(file < 0){
					printf("File not found\n");
					return -1;
				}
				
				int data_size = 0;
				int* data_pt = &data_size;
				double* data = get_ast_rad(file, argv[3], argv[4], data_pt);
				if(data != NULL){
					printf("Data size: %d\n", *data_pt);	
				}
				int lat_size = 0;
				int* lat_pt = &lat_size;
				double* lat_data = get_ast_lat(file, argv[3], argv[4], lat_pt);
				printf("Lat size: %d\n", *lat_pt);
				
				int long_size = 0;
				int* long_pt = &long_size;
				printf("going into ast_long\n");
				double* long_data = get_ast_long(file, argv[3], argv[4], long_pt);
				printf("Long size: %d\n", *long_pt);	
				
				herr_t ret = af_close(file);	
		}
	}
	else{
		printf("Invalid instrument\n");
	}
	
	return 0;
}*/


//Helper Functions
//String helper

/*
						concat_by_sep
	DESCRIPTION:	
		This is a helper function for forming paths when reading in data in a HDF5 file. Words are concatenated with a specified separator 
		
	ARGUMENTS:
		0. source -- This is where the resultant path would be stored
		1. w -- The strings to be concatenated 
		2. sep -- the separator to be used 
		3. length -- length of the resultant path
		4. arr_size -- number of strings in w
		
	EFFECT:
		Memory would be allocated to source to hold the path
		
	RETURN:
		Nothing would be returned
		
*/


void concat_by_sep(char** source, const char** w, char* sep, size_t length, int arr_size){
	int i;
	*source = calloc(length+20, sizeof(char));
	for(i = 0; i < arr_size; i++){
		if(i == 0){
			strncpy(*source, sep, strlen(sep));
			strncat(*source, w[i], strlen(w[i]));
		}
		else{
			strncat(*source, sep, strlen(sep));
			strncat(*source, w[i], strlen(w[i]));
		}
	}
}

/*
						dim_sum
	DESCRIPTION:	
		This function sums up the dimensions of a dataset based on the given dimension pointer, to be the size of the dataset
		
	ARGUMENTS:
		0. dims -- A variable that holds the size of each dimension of the dataset, usually obtained from af_read_size
		1. arr_len -- An integer that indicates the number of dimensions of that dataset
		
	EFFECT:
		A variable would be returned, containing the sum of the sizes of the dimensions
		
	RETURN:
		Returns sum 
*/


//Summing up dimensions
double dim_sum(hsize_t* dims, int arr_len){
	double sum = 0.0;
	int i;
	for(i = 0; i < arr_len; i++){
		if(i == 0){
			sum = (double)dims[i];
		}
		else{
			sum *= (double)dims[i];
		}
	}
	return sum;
}


/*
						misr_averaging
	DESCRIPTION:	
		This function is used when downsampling is needed during the retrieval of MISR data. A window of MISR data, 16 numbers, would be downsampled
		to an average. 
		
	ARGUMENTS:
		0. window -- An array of double that contains the data of a particular window in the MISR data
		
	EFFECT:
		The resultant value would be assigned to a different data array of the downsampled size
		
	RETURN:
		Returns sum/count, the average of window
		
*/


double misr_averaging(double window[16]){
	double sum = 0.0;
	double count = 0.0;
	int i;
	for(i = 0; i < 16; i++){
		if(window[i] < 0){
			return -999.0;
		}
		else{
			sum += window[i];
			count += 1;
		}
	}
	return sum/count;
}
