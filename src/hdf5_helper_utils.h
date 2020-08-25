/*
* The header file for the hdf5 utilities used in the
* *util.c files
*/

#ifndef HDF5_UTILS_H_   
#define HDF5_UTILS_H_

int checkfordataset(hid_t file_id, char *datasetname);
int get_rank(hid_t dataspace_id);
int get_rank_by_name(hid_t group_id, char *name);

#endif // HDF5_UTILS_H_
