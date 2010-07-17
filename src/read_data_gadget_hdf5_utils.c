/*
 * This subroutine performs the calls to the HDF5 library for the
 * GADGET data read
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
void read_gadget_hdf5_header(char   *filename,
                             int    maxtypes,
                             int    *npartoftype[maxtypes],
                             double *massoftype[maxtypes],
                             double *time,
                             double *redshift,
                             int *iFlagSfr,
                             int *iFlagFeedback,
                             int *Nall[maxtypes],
                             int *iFlagCool,
                             int *ndim,
                             int *ndimV,
                             int *nfiles,
                             int *ncol,
                             int nlabels,
                             char *label[nlabels],
                             int *ierr)
   {
   hid_t     file_id;
   hid_t     group_id, dataset_id;
   hid_t     attrib_id, dataspace_id;
   herr_t    status;
   herr_t    HDF5_error = -1;

   *ierr = 0;

   printf(" opening %s \n",filename);
   file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR opening %s \n",filename); *ierr = 1; return; }
   
   /*
    * Open the "Header" dataset and read the header information
    *
    */
   
   if (!checkfordataset(file_id,"Header"))
      {
         printf(" ERROR: \"Header\" dataset not found in GADGET HDF5 file\n");
         *ierr = 2;
         return;
      }

   group_id = H5Gopen(file_id,"Header");
   if (group_id == HDF5_error) 
      { printf("ERROR opening Header data set \n"); *ierr = 2; return; }

   int nattrib;
   int i;
   char name[256];
   nattrib = H5Aget_num_attrs(group_id);

   /*
    * Read through all of the attributes in the header, so we
    * can still spit out the values even if they are not used by SPLASH
    */
   for(i=0; i < nattrib; i++) {
      attrib_id = H5Aopen_idx(group_id,i);
      ssize_t  attr_status;
      attr_status = H5Aget_name(attrib_id, 256, name);
      
      hid_t  type_id;
      type_id = H5Aget_type(attrib_id);
      //type_class = H5Tget_native_type(type_id,H5T_DIR_ASCEND);
      if (strcmp(name,"Time")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,time);
      } else if (strcmp(name,"MassTable")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,massoftype);
         //printf(" Masses = %i %f \n",maxtypes,*massoftype); // %f %f %f\n",massoftype[0]);
      } else if (strcmp(name,"MassTable")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,massoftype);
      } else if (strcmp(name,"NumPart_ThisFile")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_INT,npartoftype);
      } else if (strcmp(name,"NumPart_Total")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_INT,Nall);
      } else if (strcmp(name,"Redshift")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,redshift);
      } else if (strcmp(name,"NumFilesPerSnapshot")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_INT,nfiles);
      } else if (strcmp(name,"Flag_Sfr")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_INT,iFlagSfr);
      } else if (strcmp(name,"Flag_Cooling")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_INT,iFlagCool);
      } else if (strcmp(name,"Flag_Feedback")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_INT,iFlagFeedback);
      } else {
         printf(" unknown attribute %s \n",name);
      }

      if (status==HDF5_error) {
         printf(" ERROR reading attribute %s \n",name);
      }

      status = H5Aclose(attrib_id);
   }

   status = H5Gclose(group_id);
   if (status == HDF5_error) 
      { printf("ERROR closing Header data set \n"); *ierr = 3; return; }
      
   /*
    * Now we need to get the number of data columns in the file
    * (from the number of datasets in the "PartType0" group)
    */
    
   group_id = H5Gopen(file_id,"PartType0");
   if (group_id == HDF5_error) 
      { printf("ERROR opening PartType0 data set \n"); *ierr = 2; return; }
   
   hsize_t ndatasets;
   status = H5Gget_num_objs(group_id, &ndatasets);
   printf(" number of datasets = %i \n",(int)ndatasets);

   *ncol  = 0;
   *ndim  = 0;
   *ndimV = 0;
   for(i=0; i < ndatasets; i++) {
       status       = H5Gget_objname_by_idx(group_id, i, name, 256);
       dataset_id   = H5Dopen(group_id,name);
       dataspace_id = H5Dget_space(dataset_id);
       int nrank    = H5Sget_simple_extent_ndims(dataspace_id);
       
       hsize_t dims[nrank], maxdims[nrank];
       nrank        = H5Sget_simple_extent_dims(dataspace_id,dims,maxdims);
       int ngas     = dims[0];
       int rank;
       if(nrank>1)
         {
           rank = dims[1];
         } else {
           rank = 1;
         }
       if (strcmp(name,"Coordinates")==0) { *ndim = rank; }
       if (strcmp(name,"Velocities" )==0) { *ndimV = rank; }
       if (strcmp(name,"ParticleIDs")==0) 
          { 
            printf(" ignoring %s \n",name);
          } else {
            printf(" %s x %i \n",name,rank);
            *ncol = *ncol + rank;       
          }
       status       = H5Dclose(dataset_id);
       //for (int j=0; j < rank; j++) {
          set_blocklabel(&i,&name);
       //}
       //label[i] = *name;
   }
   //printf(" number of columns = %i \n",*ncol);
   
   status = H5Gclose(group_id);

   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }
   
   }

void read_gadget_hdf5_data(char *filename, int *npart, int *ncol, int *isrequired, int *ierr)
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

          //receive_data_fromc(&i,&*npart,temp,tempid);
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
      //receive_data_fromc(&ncolloop,&*npart,temp,tempid);

      status = H5Dclose(SPHdataset_id);
      if (status == HDF5_error) { printf("ERROR closing SPH_density dataset \n"); *ierr = 4; }
   }
   
   // deallocate memory
   free(tempid);
   free(temp);

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
        printf(" dataset %s in file \n",name);
        if (strcmp(name,datasetname)==0) {
           //printf(" %s in file \n",datasetname);
           ispresent = 1;
        }
    }
    return ispresent;
}

