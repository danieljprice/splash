#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <hdf5.h>

/*
 *  utility function which checks whether a particular dataset
 *  is present in the HDF5 file
 */

int checkfordataset(hid_t file_id, char *datasetname)
{
  /* default is zero */
  int ispresent=0;

  /* get number of datasets in file */
  hsize_t ndatasets[1];
  H5Gget_num_objs(file_id, ndatasets);

  /* loop over all datasets looking for dataset matching datasetname
     set function value to true (1) if it is present  */
  int i;
  char name[256];
  for (i=0;i<(int)ndatasets[0];i++) {
      H5Gget_objname_by_idx(file_id, i, name, 256);
      /* printf(" dataset %s in file \n",name); */
      if (strcmp(name,datasetname)==0) {
         /*printf(" %s in file \n",datasetname);*/
         ispresent = 1;
      }
  }
  return ispresent;
}

/*
 *  utility function to get dimensionality of a dataset
 */
int get_rank(hid_t dataspace_id)
{
  int nrank    = H5Sget_simple_extent_ndims(dataspace_id);
  hsize_t dims[nrank], maxdims[nrank];
  nrank        = H5Sget_simple_extent_dims(dataspace_id,dims,maxdims);
  int rank;
  if(nrank>1)
    {
      rank = dims[1];
    } else {
      rank = 1;
    }
  return rank;
}
/*
 *  utility function to get dimensionality of a dataset
 */
int get_rank_by_name(hid_t group_id, char *name)
{

  if (!checkfordataset(group_id,name)) { return 0; }

  herr_t HDF5_error = -1;
#if H5_VERSION_GE(1,8,0)
  hid_t dataset_id   = H5Dopen2(group_id,name,H5P_DEFAULT);
#else
  hid_t dataset_id   = H5Dopen(group_id,name);
#endif
  if (dataset_id == HDF5_error) 
     { printf("ERROR opening %s data set \n",name); return 0; }

  hid_t dataspace_id = H5Dget_space(dataset_id);
  int rank           = get_rank(dataspace_id);
  H5Dclose(dataset_id);
  return rank;
}

