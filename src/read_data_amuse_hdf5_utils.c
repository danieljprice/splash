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
int read_amuse_hdf5_dataset(hid_t group_id, char *datasetname, int itype, int maxtypes, int npartoftype[maxtypes],
                            int ncol, int isrequired[ncol], int   *j);
int get_rank(hid_t dataspace_id);
int get_rank_by_name(hid_t group_id, char *name);
void set_blocklabel(int *icol, char *name);
void read_amuse_hdf5_data_fromc(int *icol, int *npartoftypei, double temparr[*npartoftypei],int *itype);
void read_amuse_hdf5_header(char   *filename,
                             int *npart,
                             int *ncol,
                             int *ndim,
                             int *ndimV,
                             double *time,
                             int *ierr)
   {
   hid_t     file_id;
   hid_t     group_id, group_id1, group_id2;
   hid_t     attrib_id;
   herr_t    status;
   herr_t    HDF5_error = -1;
   
   *ierr = 0;
   *ndim = 0;
   *ndimV = 0;

   if (debug) printf("DEBUG: opening %s \n",filename);
   file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR opening %s \n",filename); *ierr = 1; return; }
   
   char *maingroup = "particles";
   /*
    * Open the "particles" dataset and read the number of particles attribute
    *
    */
   if (!checkfordataset(file_id,maingroup))
      {
         printf(" ERROR: \"%s\" dataset not found in AMUSE HDF5 file\n",maingroup);
         *ierr = 2;
         return;
      }

#if H5_VERSION_GE(1,8,0)
   group_id1 = H5Gopen2(file_id,maingroup,H5P_DEFAULT);
#else
   group_id1 = H5Gopen(file_id,maingroup);
#endif
   if (group_id1 == HDF5_error) 
      { printf("ERROR opening %s data set \n",maingroup); *ierr = 2; return; }

#if H5_VERSION_GE(1,8,0)
   group_id = H5Gopen2(group_id1,"0000000001",H5P_DEFAULT);
#else
   group_id = H5Gopen(group_id1,"0000000001");
#endif
   if (group_id == HDF5_error) 
      { printf("ERROR opening 00000000001 data set \n"); *ierr = 2; return; }

   int nattrib;
   int i;
   char name[256];
   nattrib = H5Aget_num_attrs(group_id);
   if (debug) printf("number of attributes found = %i\n",nattrib);

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
      /*type_class = H5Tget_native_type(type_id,H5T_DIR_ASCEND);*/
      if (strcmp(name,"time")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,time);
      } else if (strcmp(name,"number_of_particles")==0) {
         status = H5Aread(attrib_id,H5T_NATIVE_INT,npart);
      } else {
         if (debug) printf("DEBUG: unknown attribute %s \n",name);
      }

      if (status==HDF5_error) {
         printf(" ERROR reading attribute %s \n",name);
      }

      status = H5Aclose(attrib_id);
   }

   /*
    * Now we need to get the number of data columns in the file
    * (from the number of datasets in the "attributes" group)
    */
#if H5_VERSION_GE(1,8,0)
   group_id2 = H5Gopen2(group_id,"attributes",H5P_DEFAULT);
#else
   group_id2 = H5Gopen(group_id,"attributes");
#endif
   if (group_id2 == HDF5_error) 
      { printf("ERROR opening %s data set \n","attributes"); *ierr = 2; return; }
   
   hsize_t ndatasets;
   status = H5Gget_num_objs(group_id2, &ndatasets);
   if (debug) printf("DEBUG: number of datasets = %i \n",(int)ndatasets);

   *ncol  = (int)ndatasets;
   int idim;

   /* check that coordinates are present in file */
   idim = get_rank_by_name(group_id2,"x");
   if (idim <= 0) { printf("ERROR: x positions not found\n"); *ierr = 3; }

   idim = get_rank_by_name(group_id2,"y");
   if (idim <= 0) { printf("ERROR: y positions not found\n"); *ierr = 3; }
 
   idim = get_rank_by_name(group_id2,"z");
   if (idim <= 0) { 
      printf("z positions not found, assuming file is 2D \n");
      *ndim = 2;
      *ndimV = 2;
   } else {
      *ndim = 3;
      *ndimV = 3;
   }

   /* finish, close all open datasets and close file */

   status = H5Gclose(group_id2);
   if (status == HDF5_error)
      { printf("ERROR closing attributes data set \n"); *ierr = 3; return; }

   status = H5Gclose(group_id);
   if (status == HDF5_error) 
      { printf("ERROR closing %s data set \n",maingroup); *ierr = 3; return; }

   status = H5Gclose(group_id1);
   if (status == HDF5_error) 
      { printf("ERROR closing 0001 data set \n"); *ierr = 3; return; }

   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }
   if (debug) printf("DEBUG: finished header read \n");
   
   }

void read_amuse_hdf5_data(char *filename,
                           int maxtypes,
                           int npartoftype[maxtypes],
                           int ncol,
                           int isrequired[ncol],
                           int *ierr)
   {
   hid_t     file_id;
   hid_t     group_id, group_id1, group_id2;
   herr_t    status;
   herr_t    HDF5_error = -1;
   char      groupname[12];
   char      datasetname[256];
   int i;

   if (debug) printf("DEBUG: re-opening %s \n",filename);
   file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR re-opening %s \n",filename); *ierr = 1; return; }

   /* open main particles group */
   char *maingroup = "particles";
#if H5_VERSION_GE(1,8,0)
   group_id1 = H5Gopen2(file_id,maingroup,H5P_DEFAULT);
#else
   group_id1 = H5Gopen(file_id,maingroup);
#endif
   if (debug) printf("DEBUG: maxtypes = %i\n",maxtypes);
   
   /* read dataset for each particle type present in dump file */
   int itype,iobjtype;
   for (itype=0;itype<maxtypes;itype++) {
      if (debug) printf("DEBUG: type %i npartoftype = %i\n",itype,npartoftype[itype]);
      if (npartoftype[itype] > 0) {
         /* If npartoftype[N] > 0 in header, look for dataset of the form 000000000N */
         sprintf(groupname,"00000000%02i",itype+1);
         if (debug) printf("DEBUG: opening group %s\n",groupname);
#if H5_VERSION_GE(1,8,0)
         group_id = H5Gopen2(group_id1,groupname,H5P_DEFAULT);
#else
         group_id = H5Gopen(group_id1,groupname);
#endif
         if (group_id == HDF5_error)
            { printf("ERROR opening %s group \n",groupname); *ierr = 2; }
         else {

#if H5_VERSION_GE(1,8,0)
            group_id2 = H5Gopen2(group_id,"attributes",H5P_DEFAULT);
#else
            group_id2 = H5Gopen(group_id,"attributes");
#endif
            if (group_id2 == HDF5_error)
               { printf("ERROR opening attributes group \n"); *ierr = 2; }
            else {

               hsize_t ndatasets;
               status = H5Gget_num_objs(group_id2, &ndatasets);
               if (debug) printf("DEBUG: number of datasets = %i \n",(int)ndatasets);

               int j = 0;
               /* always read particle positions first */
               *ierr = read_amuse_hdf5_dataset(group_id2,"x",itype,maxtypes,npartoftype,ncol,isrequired,&j);
               j = 1;
               *ierr = read_amuse_hdf5_dataset(group_id2,"y",itype,maxtypes,npartoftype,ncol,isrequired,&j);
               j = 2;
               *ierr = read_amuse_hdf5_dataset(group_id2,"z",itype,maxtypes,npartoftype,ncol,isrequired,&j);

               /* read remaining datasets in the order they appear in the file */
               for(i=0; i < (int)ndatasets; i++) {
                   status       = H5Gget_objname_by_idx(group_id2, i, datasetname, 256);
                   iobjtype     = H5Gget_objtype_by_idx(group_id2, i);
                   if (strcmp(datasetname,"x")&&
                       strcmp(datasetname,"y")&&
                       strcmp(datasetname,"z")&& (iobjtype == H5G_DATASET)) {
                      *ierr = read_amuse_hdf5_dataset(group_id2,datasetname,itype,maxtypes,npartoftype,ncol,isrequired,&j);
                   }
               }

               /* close "attributes" group */
               H5Gclose(group_id2);

            }
            
            H5Gclose(group_id);
         }
      }
   }
   H5Gclose(group_id1);

   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }

   }

int read_amuse_hdf5_dataset(hid_t group_id,
                            char  *datasetname,
                            int   itype,
                            int   maxtypes,
                            int   npartoftype[maxtypes],
                            int   ncol,
                            int   isrequired[ncol],
                            int   *j) 
{
   hid_t dataset_id, dataspace_id, memspace_id;
   herr_t    status;
   herr_t    HDF5_error = -1;
   int       ierr = 0;
   char      name[256];
   
   if (!checkfordataset(group_id,datasetname)) { ierr = 1; return ierr; }

#if H5_VERSION_GE(1,8,0)
   dataset_id   = H5Dopen2(group_id,datasetname,H5P_DEFAULT);
#else
   dataset_id   = H5Dopen(group_id,datasetname);
#endif
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
      read_amuse_hdf5_data_fromc(j,&npartoftype[itype],temp,&itype);

      strcpy(name,datasetname);
      set_blocklabel(j,name);

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
            read_amuse_hdf5_data_fromc(j,&npartoftype[itype],temp,&itype);
         }
         strcpy(name,datasetname);
         set_blocklabel(j,name);
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
