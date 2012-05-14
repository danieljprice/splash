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
static int debug = 0;

int checkfordataset(hid_t file_id, char *datasetname);
int read_gadgethdf5_dataset(hid_t group_id, char *datasetname, int itype, int maxtypes, int npartoftype[maxtypes],
                            int i0[maxtypes], int ncol, int isrequired[ncol], int   *j);
int get_rank(hid_t dataspace_id);
int get_rank_by_name(hid_t group_id, char *name);
void set_blocklabel(int *icol, int *irank, char *name);
void read_gadgethdf5_data_fromc(int *icol, int *npartoftypei, double temparr[*npartoftypei],
                                int *itype, int *i0);

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
                             int *ierr)
   {
   hid_t     file_id;
   hid_t     group_id, dataset_id;
   hid_t     attrib_id, dataspace_id;
   herr_t    status;
   herr_t    HDF5_error = -1;

   *ierr = 0;

   if (debug) printf("DEBUG: opening %s \n",filename);
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
   double BoxSize,HubbleParam,Omega0,OmegaLambda;
   int iFlagStellarAge,iFlagMetals;
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
      } else if (strcmp(name,"BoxSize")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,&BoxSize);
      } else if (strcmp(name,"HubbleParam")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,&HubbleParam);
      } else if (strcmp(name,"Omega0")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,&Omega0);
      } else if (strcmp(name,"OmegaLambda")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,&OmegaLambda);
      } else if (strcmp(name,"Flag_StellarAge")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_INT,&iFlagStellarAge);
      } else if (strcmp(name,"Flag_Metals")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_INT,&iFlagMetals);
      } else {
         if (debug) printf("DEBUG: unknown attribute %s \n",name);
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
   if (debug) printf("DEBUG: number of datasets = %i \n",(int)ndatasets);

   *ncol  = 0;
   *ndim  = 0;
   *ndimV = 0;
   int rank = 0;
   int j = 0;

   *ndim = get_rank_by_name(group_id,"Coordinates");
   set_blocklabel(&j,ndim,"Coordinates");
   *ncol = *ncol + *ndim;
   if (*ndim > 0) j++;

   *ndimV = get_rank_by_name(group_id,"Velocities");
   set_blocklabel(&j,ndimV,"Velocities");
   *ncol = *ncol + *ndimV;
   if (*ndimV > 0) j++;

   rank = get_rank_by_name(group_id,"Masses");
   set_blocklabel(&j,&rank,"Masses");
   *ncol = *ncol + rank;
   if (rank > 0) j++;
   
   for(i=0; i < (int)ndatasets; i++) {
       status       = H5Gget_objname_by_idx(group_id, i, name, 256);
       dataset_id   = H5Dopen(group_id,name);
       dataspace_id = H5Dget_space(dataset_id);
       rank         = get_rank(dataspace_id);

       if (strcmp(name,"ParticleIDs")&&
           strcmp(name,"Coordinates")&&
           strcmp(name,"Velocities")&&
           strcmp(name,"Masses"))
          {
            if (debug) printf("DEBUG: %s x %i \n",name,rank);
            /* Send the dataset names back to Fortran 
             * one by one, so they can be filled into 
             * the array as appropriate */
            set_blocklabel(&j,&rank,name);
            *ncol = *ncol + rank;
            if (rank > 0) j++;
          } else {
            if (debug) printf("DEBUG: ignoring %s \n",name);
          }
       status       = H5Dclose(dataset_id);       
   } 
  
   status = H5Gclose(group_id);

   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }
   
   }

void read_gadget_hdf5_data(char *filename,
                           int maxtypes,
                           int    npartoftype[maxtypes],
                           double massoftype[maxtypes],
                           int ncol,
                           int isrequired[ncol],
                           int i0[maxtypes],
                           int *ierr)
   {
   hid_t     file_id;
   hid_t     group_id;
   herr_t    status;
   herr_t    HDF5_error = -1;
   char      groupname[12];
   char      datasetname[256];
   
   int i;

   if (debug) printf("DEBUG: re-opening %s \n",filename);
   file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR re-opening %s \n",filename); *ierr = 1; return; }
   
   /* read dataset for each particle type present in dump file */
   int itype;
   for (itype=0;itype<maxtypes;itype++) {
      if (npartoftype[itype] > 0) {
         /* If npartoftype[N] > 0 in header, look for dataset of the form PartTypeN */
         sprintf(groupname,"PartType%i",itype);
         if (debug) printf("DEBUG: opening group %s\n",groupname);
         group_id = H5Gopen(file_id,groupname);
         if (group_id == HDF5_error)
            { printf("ERROR opening %s group \n",groupname); *ierr = 2; }
         else {
            hsize_t ndatasets;
            status = H5Gget_num_objs(group_id, &ndatasets);
            if (debug) printf("DEBUG: number of datasets = %i \n",(int)ndatasets);

            int j = 0;
            /* read datasets common to all particle types first */
            *ierr = read_gadgethdf5_dataset(group_id,"Coordinates",itype,maxtypes,npartoftype,i0,ncol,isrequired,&j);
            *ierr = read_gadgethdf5_dataset(group_id,"Velocities",itype,maxtypes,npartoftype,i0,ncol,isrequired,&j);
            *ierr = read_gadgethdf5_dataset(group_id,"Masses",itype,maxtypes,npartoftype,i0,ncol,isrequired,&j);
/*          *ierr = read_gadgethdf5_dataset(group_id,"SmoothingLength",itype,maxtypes,npartoftype,i0,ncol,isrequired,&j);
            *ierr = read_gadgethdf5_dataset(group_id,"Density",itype,maxtypes,npartoftype,i0,ncol,isrequired,&j);
*/          
            /* read remaining datasets in the order they appear in the file */
            for(i=0; i < (int)ndatasets; i++) {
                status       = H5Gget_objname_by_idx(group_id, i, datasetname, 256);
                if (strcmp(datasetname,"ParticleIDs")&&
                    strcmp(datasetname,"Coordinates")&&
                    strcmp(datasetname,"Velocities")&&
                    strcmp(datasetname,"Masses")) {
                   *ierr = read_gadgethdf5_dataset(group_id,datasetname,itype,maxtypes,npartoftype,i0,ncol,isrequired,&j);
                }
            }
            H5Gclose(group_id);
         }
      }
   }

   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }

   }

int read_gadgethdf5_dataset(hid_t group_id,
                            char  *datasetname,
                            int   itype,
                            int   maxtypes,
                            int   npartoftype[maxtypes],
                            int   i0[maxtypes],
                            int   ncol,
                            int   isrequired[ncol],
                            int   *j) 
{
   hid_t dataset_id, dataspace_id, memspace_id;
   herr_t    status;
   herr_t    HDF5_error = -1;
   int       ierr;
   
   if (!checkfordataset(group_id,datasetname)) { ierr = 1; return ierr; }

   dataset_id   = H5Dopen(group_id,datasetname);
   dataspace_id = H5Dget_space(dataset_id);
   int rank     = get_rank(dataspace_id);
   int k, flag;

   /* do nothing if none of the columns are required */
   flag = 0;
   for (k=0;k<rank;k++) {
      if (isrequired[k]) flag = 1;
   }
   if (!flag) {
      if (debug) printf("DEBUG: skipping %s : not required\n",datasetname);
      H5Dclose(dataset_id);
      return 0;
   }

   if (debug) printf("DEBUG: got %s rank %i \n",datasetname,rank);
   /* make a temporary array to put each column as we read it */
   hsize_t nparttmp[1];
   nparttmp[0] = npartoftype[itype];
   memspace_id = H5Screate_simple(1,nparttmp,NULL);
   double *temp = 0;
   temp = malloc(npartoftype[itype]*sizeof(double));

   if (rank==1) {
      *j = *j + 1;

      /* read column from file */
      H5Dread(dataset_id,H5T_NATIVE_DOUBLE,memspace_id,dataspace_id,H5P_DEFAULT,temp);

      /* call Fortran back, sending values in temp array to fill into the main splash dat array */
      read_gadgethdf5_data_fromc(j,&npartoftype[itype],temp,&itype,&i0[itype]);

   } else {
      hsize_t offset[2], count[2];

      for (k=1;k<=rank;k++) {
         *j = *j + 1;
         count[0] = npartoftype[itype];
         count[1] = 1;
         offset[0] = 0;
         offset[1] = k-1; /* rank */
         status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
         if (status == HDF5_error) { printf("ERROR creating hyperslab \n"); ierr = 4; }
         if (!H5Sselect_valid(dataspace_id)) { printf("ERROR selecting hyperslab \n"); ierr = 5; }

         /* read column from file */
         if (isrequired[*j]) {
            H5Dread(dataset_id,H5T_NATIVE_DOUBLE,memspace_id,dataspace_id,H5P_DEFAULT,temp);

            /* call Fortran back, sending values in temp array to fill into the main splash dat array */
            read_gadgethdf5_data_fromc(j,&npartoftype[itype],temp,&itype,&i0[itype]);
         }
      }
   }
   free(temp);

   status = H5Sclose(dataspace_id);
   if (status == HDF5_error) { printf("ERROR closing dataspace \n"); ierr = 4; }

   H5Dclose(dataset_id);
   if (status == HDF5_error) { printf("ERROR closing dataset \n"); ierr = 7; }
   return ierr;
}

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
         //printf(" %s in file \n",datasetname);
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
  hid_t dataset_id   = H5Dopen(group_id,name);
  hid_t dataspace_id = H5Dget_space(dataset_id);
  int rank           = get_rank(dataspace_id);
  H5Dclose(dataset_id);
  return rank;
}
