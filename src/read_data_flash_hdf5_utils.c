/*
 * This subroutine performs the calls to the HDF5 library for the FLASH
 * tracer particles data read
 *
 * Easier to do it this way and link with c than to try to link against
 * the Fortran interface (in the latter case the modules must
 * have been compiled with the *exact* compiler used to compile splash
 * which is a real pain).
 *
 */
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <hdf5.h>
void read_flash_hdf5_header(char *filename, float *time, int *npart, int *ncol, int *ierr)
   {
   hid_t     file_id;
   hid_t     dataset_id;
   hid_t     dataspace_id;
   herr_t    status;
   herr_t    HDF5_error = -1;

   //printf(" opening %s \n",filename);
   file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR opening %s \n",filename); *ierr = 1; return; }
   
   // READ NPART AND NCOL from file
   
   dataset_id = H5Dopen(file_id,"tracer particles");
   if (dataset_id == HDF5_error) 
      { printf("ERROR opening tracer particle data set \n"); *ierr = 2; return; }
   
   dataspace_id = H5Dget_space(dataset_id);
      
   // get dimensional information from dataspace
   hsize_t HDFxdims[4], HDFmaxdims[4];
   int rank = H5Sget_simple_extent_dims(dataspace_id, HDFxdims, HDFmaxdims);
   if (rank > 4) { printf("RANK of dataset exceeds array bounds \n"); *ierr = 3; return; }

   // from the dimensional info, calculate the size of the buffer.
   *npart = HDFxdims[0];
   *ncol = HDFxdims[1];
   //printf(" number of particles %i \n",HDFxdims[0]);
   //printf(" number of columns = %i \n",HDFxdims[1]);

   status = H5Sclose(dataspace_id);
   if (status == HDF5_error) { printf("ERROR closing dataspace \n"); *ierr = 4; }
   
   status = H5Dclose(dataset_id);
   if (status == HDF5_error) { printf("ERROR closing dataset \n"); *ierr = 4; }
   
   /* 
    * read the time from the file - this is 
    * contained in a pointlessly complicated
    * compound structure that we have to replicate here
    *
    */
   
   dataset_id = H5Dopen(file_id,"real scalars");
   if (dataset_id == HDF5_error) 
      { printf("ERROR opening real scalars data set for time \n"); *ierr = 5; return; }
   
   dataspace_id = H5Dget_space(dataset_id);
   rank = H5Sget_simple_extent_dims(dataspace_id, HDFxdims, HDFmaxdims);
   if (rank > 4) { printf("RANK of dataset exceeds array bounds \n"); *ierr = 3; return; }
   
   int lenheader = HDFxdims[0];
   //printf(" header contains %i items \n",lenheader);
   
   typedef struct h_t {
       char     name[80];
       double   value;
   } h_t;
   h_t    header[lenheader];

   hid_t strtype = H5Tcopy(H5T_C_S1);
   status = H5Tset_size (strtype, 80);
   
   hid_t compound_id = H5Tcreate(H5T_COMPOUND, sizeof(h_t));
   H5Tinsert(compound_id,"name",HOFFSET(h_t, name), strtype);
   H5Tinsert(compound_id,"value",HOFFSET(h_t, value), H5T_NATIVE_DOUBLE);
   
   status = H5Dread(dataset_id,compound_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,header);
   if (status == HDF5_error) { printf("ERROR reading header \n"); *ierr = 6; }
   *time = header[0].value;
   
   status = H5Sclose(dataspace_id);
   if (status == HDF5_error) { printf("ERROR closing dataspace \n"); *ierr = 4; }
   
   status = H5Dclose(dataset_id);
   if (status == HDF5_error) { printf("ERROR closing dataset \n"); *ierr = 4; }

   // Additionally check if an SPH_density data set is present
   if (checkfordataset(file_id,"SPH_density")==1) {
      printf(" file contains SPH-calculated densities \n");
      *ncol += 1;
   } else {
      printf(" file does not contain SPH_density data set \n");
   }
   
   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }
   
   }

void read_flash_hdf5_data(char *filename, int *npart, int *ncol, int *isrequired, int *ierr)
   {
   hid_t     file_id;
   hid_t     dataset_id, SPHdataset_id;
   hid_t     dataspace_id, SPHdataspace_id;
   hid_t     memspace_id;
   herr_t    status;
   herr_t    HDF5_error = -1;

   //printf(" re-opening %s \n",filename);
   file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR re-opening %s \n",filename); *ierr = 1; return; }
   
   // re-open tracer particle dataspace
   
   dataset_id = H5Dopen(file_id,"tracer particles");
   if (dataset_id == HDF5_error) 
      { printf("ERROR opening tracer particle data set \n"); *ierr = 2; return; }

   dataspace_id = H5Dget_space(dataset_id);

   // Additionally check if an SPH_density data set is present
    
   int ncolloop;
   int gotSPHdata = checkfordataset(file_id,"SPH_density");
   if (gotSPHdata==0) {
      //printf(" file does not contain SPH_density data set \n");
      ncolloop = *ncol;
   } else {
      //printf(" file contains SPH-calculated densities \n");
      ncolloop = *ncol - 1;
   }

   // make a temporary space to put each column as we read it
   hsize_t nparth[1];
   nparth[0] = *npart;
   memspace_id = H5Screate_simple(1,nparth,NULL);
   
   // dynamically allocate a temporary double array to store one column
   double* temp = 0;
   temp = malloc(*npart*sizeof(double));

   // read particle information into one column
   hsize_t offset[2], count[2];
   
   int i;
   count[0] = *npart;
   count[1] = 1;
   /* 
    * read particle IDs first so we can sort into ID order
    */
   offset[0] = 0;
   offset[1] = 5;  // ID in column 5: should check this but this is for the future
   status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
   if (status == HDF5_error) { printf("ERROR creating hyperslab \n"); *ierr = 4; }
   if (!H5Sselect_valid(dataspace_id)) { printf("ERROR selecting hyperslab \n"); *ierr = 5; }
   // read ID
   H5Dread(dataset_id,H5T_NATIVE_DOUBLE,memspace_id,dataspace_id,H5P_DEFAULT,temp);
   int* tempid = 0;
   tempid = malloc(*npart*sizeof(int));
   
   // convert temp (double) into integer array
   for (i=0;i<*npart;i++) {
       tempid[i] = (int)temp[i];
   }

   /* 
    * start loop from 1 because first array is a useless "tag" array that
    * we don't need to read.
    *
    */
   for (i=1;i<ncolloop;i++) {
       if (i != 5 && isrequired[i-1]==1) {  // skip column 5 which is the particle ID, above
          offset[0] = 0;
          offset[1] = i;
          //printf(" getting column %i offset %i \n",i,offset[0]);
          //printf(" getting column %i count %i \n",i,count[0]);

          status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
          if (status == HDF5_error) { printf("ERROR creating hyperslab \n"); *ierr = 4; }

          if (!H5Sselect_valid(dataspace_id)) { printf("ERROR selecting hyperslab \n"); *ierr = 5; }
          //printf(" getting column %i count %i \n",i,count[0]);

          H5Dread(dataset_id,H5T_NATIVE_DOUBLE,memspace_id,dataspace_id,H5P_DEFAULT,temp);

          // call Fortran back, sending values in temp array to fill into the main splash dat array

          receive_data_fromc(&i,&*npart,temp,tempid);
       }
   }

   status = H5Sclose(dataspace_id);
   if (status == HDF5_error) { printf("ERROR closing dataspace \n"); *ierr = 4; }
   
   status = H5Dclose(dataset_id);
   if (status == HDF5_error) { printf("ERROR closing dataset \n"); *ierr = 4; }
   
   // Additionally read SPH_density data set if it is present
   if (gotSPHdata==1) {
      SPHdataset_id = H5Dopen(file_id,"SPH_density");

      printf(" reading SPH_density data set \n");

      status = H5Dread(SPHdataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,temp);
      if (status == HDF5_error) { printf("ERROR reading SPH_density dataspace \n"); *ierr = 3; }

      // call Fortran back, sending values in temp array to fill into the main splash dat array
      receive_data_fromc(&ncolloop,&*npart,temp,tempid);

      status = H5Dclose(SPHdataset_id);
      if (status == HDF5_error) { printf("ERROR closing SPH_density dataset \n"); *ierr = 4; }
   }
   
   // deallocate memory
   free(tempid);
   free(temp);

   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }

   }

void write_tracer_particle_density(char *filename, int *npart, double *dens, int *ierr)
   {
   hid_t     file_id;
   hid_t     dataset_id;
   hid_t     memspace_id, dataspace_id;
   herr_t    status;
   herr_t    HDF5_error = -1;

   printf(" opening %s for write \n",filename);
   file_id = H5Fopen(filename,H5F_ACC_RDWR,H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR opening %s for write \n",filename); *ierr = 1; return; }
   
   if (checkfordataset(file_id,"SPH_density")==1) {
   // SPH data set already exists, we are going to overwrite it
      dataset_id = H5Dopen(file_id,"SPH_density");
      if (dataset_id == HDF5_error)
         { printf("ERROR opening SPH_density data set for overwrite \n"); *ierr = 2; return; }
   
   } else {
   // create new dataspace for tracer particle density
      hsize_t nparth[1];
      nparth[0] = *npart;
      memspace_id = H5Screate_simple(1,nparth,NULL);
 
      dataset_id = H5Dcreate(file_id,"SPH_density",H5T_NATIVE_DOUBLE,memspace_id,H5P_DEFAULT);
      if (dataset_id == HDF5_error)
         { printf("ERROR creating SPH_density data set \n"); *ierr = 2; return; }
   }

   // write dataset to file
   status = H5Dwrite(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,dens);
   if (status == HDF5_error) { printf("ERROR writing to file \n"); *ierr = 4; }

   // clean up

   status = H5Dclose(dataset_id);
   if (status == HDF5_error) { printf("ERROR closing dataset \n"); *ierr = 6; }

   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }

   }

/*
 *  utility function which checks whether a particular dataset
 *  is present in the HDF5 file
 */
int checkfordataset(hid_t file_id, char *datasetname)
{
    // default is zero

    int ispresent=0;
    
    // get number of datasets in file
    
    hsize_t ndatasets[1];
    H5Gget_num_objs(file_id, ndatasets);
    
    // loop over all datasets looking for the "SPH_density" dataset
    // set function value to true (1) if it is present
    
    int i;
    char name[256];
    for (i=0;i<ndatasets[0];i++) {
        H5Gget_objname_by_idx(file_id, i, name, 256);
        if (strcmp(name,datasetname)==0) {
           //printf(" %s in file \n",datasetname);
           ispresent = 1;
        }
    }
    return ispresent;
}

