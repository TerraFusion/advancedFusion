/****************************************************************************
 * DESCRIPTION:
 *  IO functions for developing the Advance Fusion tool for the ACCESS TERRA
 *  Fusion project. Basic Fusion TERRA data is used as an input and retrives
 *  data from instruments by specifying desired paraameters. The data is
 *  used for resampling and reprojection.
 *
 * DEVELOPERS:
 *  - Jonathan Kim (jkm@illinois.edu)
 *  - Yat Long Lo (yllo2@illinois.edu) - Author
 *
 */
#include <vector>
#include <string>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#define FALSE   0

//HDF5 API operations wrapper
hid_t af_open(char* file_path);
herr_t af_close(hid_t file);
double* af_read(hid_t file, char* dataset_name);
double* af_read_hyperslab(hid_t file, char*dataset_name, int x_offset, int y_offset, int z_offset);
hsize_t* af_read_size(hid_t file, char* dataset_name);
//int af_write_misr_on_modis(hid_t output_file, double* misr_out, double* modis, int modis_size, int modis_band_size, int misr_size);
int af_write_misr_on_modis(hid_t output_file, double* misr_out, double* modis, uint64_t modis_size, int modis_band_size, int misr_size);
//int af_write_mm_geo(hid_t output_file, int geo_flag, double* geo_data, int geo_size, int outputWidth,hid_t ctrackDset,hid_t atrackDset);
int af_write_mm_geo(hid_t output_file, int geo_flag, double* geo_data, uint64_t geo_size, int outputWidth,hid_t ctrackDset,hid_t atrackDset);
int af_write_attr_float(hid_t dset, char* name, float val);
int af_write_attr_str(hid_t dset, char* name, char* val);
int af_write_cf_attributes(hid_t dset, char* units, float _FillValue,
                           float valid_min,float valid_max, unsigned short handle_flag);
    
//Instrument data retrieval functions
//double* get_misr_rad(hid_t file, char* camera_angle, char* resolution, char* radiance, int* size);
double* get_misr_rad(hid_t file, char* camera_angle, char* resolution, char* radiance, uint64_t* size);
//double* get_misr_lat(hid_t file, char* resolution, int* size);
double* get_misr_lat(hid_t file, char* resolution, uint64_t* size);
//double* get_misr_long(hid_t file, char* resolution, int* size);
double* get_misr_long(hid_t file, char* resolution, uint64_t* size);
void* get_misr_attr(hid_t file, char* camera_angle, char* resolution, char* radiance, char* attr_name, int geo, void* attr_pt);
//double* get_modis_rad(hid_t file, char* resolution, std::vector<std::string> &bands, int band_size, int* size);
double* get_modis_rad(hid_t file, char* resolution, std::vector<std::string> &bands, int band_size, uint64_t* size);
//double* get_modis_rad_by_band(hid_t file, char* resolution, char* d_name, int* band_index, int* size);
double* get_modis_rad_by_band(hid_t file, char* resolution, char* d_name, int* band_index, uint64_t* size);
//double* get_modis_lat(hid_t file, char* resolution, int* size);
double* get_modis_lat(hid_t file, char* resolution, uint64_t * size);
//double* get_modis_long(hid_t file, char* resolution, int* size);
double* get_modis_long(hid_t file, char* resolution, uint64_t* size);
void* get_modis_attr(hid_t file, char* resolution, char* d_name, char* attr_name, int geo, void* attr_pt);
char* get_modis_filename(char* resolution, char* band, int* band_index);
double* get_ceres_rad(hid_t file, char* camera, char* d_name, int* size);
double* get_ceres_lat(hid_t file, char* camera, char* d_name, int* size);
double* get_ceres_long(hid_t file, char* camera, char* d_name, int* size);
double* get_mop_rad(hid_t file, int* size);
double* get_mop_lat(hid_t file, int*size);
double* get_mop_long(hid_t file, int* size);
//double* get_ast_rad(hid_t file, char* subsystem, char* d_name, int*size);
double* get_ast_rad(hid_t file, char* subsystem, char* d_name, uint64_t *size);
//double* get_ast_lat(hid_t file, char* subsystem, char* d_name, int*size);
double* get_ast_lat(hid_t file, char* subsystem, char* d_name, uint64_t *size);
//double* get_ast_long(hid_t file, char* subsystem, char* d_name, int*size);
double* get_ast_long(hid_t file, char* subsystem, char* d_name, uint64_t *size);
//double* get_ast_rad_by_gran(hid_t file, char* subsystem, char* d_name, char* gran_name, int*size);
double* get_ast_rad_by_gran(hid_t file, char* subsystem, char* d_name, char* gran_name, uint64_t *size);
//double* get_ast_lat_by_gran(hid_t file, char* subsystem, char* d_name, char* gran_name, int*size);
double* get_ast_lat_by_gran(hid_t file, char* subsystem, char* d_name, char* gran_name, uint64_t*size);
//double* get_ast_long_by_gran(hid_t file, char* subsystem, char* d_name, char* gran_name, int*size);
double* get_ast_long_by_gran(hid_t file, char* subsystem, char* d_name, char* gran_name, uint64_t *size);

//Helper functions
void concat_by_sep(char** source, const char** w, char* sep, size_t length, int arr_size);
//double dim_sum(hsize_t* dims, int arr_len);
double dim_sum(hsize_t* dims, uint64_t arr_len);
//double dim_sum_free(hsize_t* dims, int arr_len);
double dim_sum_free(hsize_t* dims, uint64_t arr_len);
double float_to_double(float f);
double misr_averaging(double window[16]);
hid_t  create_pure_dim_dataset(hid_t loc_id, hsize_t dim_size,char* dim_name);
