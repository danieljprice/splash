/*
 * This subroutine performs the calls to the HDF5 library for the
 * CACTUS data read
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
static int debug = 1;
int read_cactus_dataset(hid_t file_id,char *name,int *nattrib,int *ncells,int *ndim,int *n,double *time,int inheader);
int read_cactus_iteration(hid_t file_id,int iter,int *next,int *ncells,int *ndim,int *nattrib,double *time,int inheader);
int read_cactus_grid(hid_t dataset_id,hid_t dataspace_id,int ndim,int *n,int nx,int ny,int nz,int nghost[3],int ioffset[3],double orig[3],double delta[3]);
void get_ndim_ncells(hid_t dataspace_id, int *ndim, int *nx,int *ny,int *nz);
void set_blocklabel(int *icol, char *name);
void read_cactus_hdf5_data_fromc(int *icol,int *ntot,int *np,double temparr[*np]);
void read_cactus_itype_fromc(int *ntot,int *np,int itype[*np]);

void read_cactus_hdf5_header(char   *filename,
                             int *npart,
                             int *ncol,
                             int *ndim,
                             int *ndimV,
                             double *time,
                             int *ierr)
   {
   hid_t     file_id;
   herr_t    status;
   herr_t    HDF5_error = -1;
   
   *ierr = 0;
   *ndim = 0;
   *ndimV = 0;
   *ncol = 6; /* x,y,z,h,m,density */

   if (debug) printf("DEBUG: opening %s \n",filename);
   file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR opening %s \n",filename); *ierr = 1; return; }
   
   /*
    * Open the first iteration and read the number of dimensions / cells
    *
    */
   int nattrib;
   int iter = 0;
   int next;
   *ierr = read_cactus_iteration(file_id,iter,&next,npart,ndim,&nattrib,time,1);
   *ierr = 0;

   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }
   if (debug) printf("DEBUG: finished header read \n");

   *ndimV = *ndim;
   return;
   
   }

void read_cactus_hdf5_data(char *filename,
                           int maxtypes,
                           int npartoftype[maxtypes],
                           int ncol,
                           int isrequired[ncol],
                           int *ierr)
   {
   hid_t     file_id;
   herr_t    status;
   herr_t    HDF5_error = -1;
   int next,nattrib,npart,ndim,nsteps,iter;
   double time;

   if (debug) printf("DEBUG: re-opening %s \n",filename);
   file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
   if (file_id == HDF5_error)
      { printf("ERROR re-opening %s \n",filename); *ierr = 1; return; }

   iter = 0;
   nsteps = read_cactus_iteration(file_id,iter,&next,&npart,&ndim,&nattrib,&time,0);
   if (nsteps <= 0) {
      *ierr = -1;
   }

   status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 7; }

   }


/*
 *  utility function which checks whether a particular dataset
 *  is present in the HDF5 file
 */
int read_cactus_iteration(hid_t file_id,int iter,int *next,int *ncells,int *ndim,int *nattrib,double *time,int inheader)
{
  hsize_t ndatasets[1];
  /* get number of datasets in file */
  H5Gget_num_objs(file_id, ndatasets);
  int it,tl,level,cnum;

  /* loop over all datasets looking for dataset matching the desired iteration number
     set function value to true (1) if it is present  */
  int i,nsub,ierr;
  char name[256];
  char field[12];
  int nsteps = 0;
  *ndim = 0;
  *ncells = 0;
  if (inheader) {
     printf("%i datasets in file\n",(int)ndatasets[0]);
  } else {
     printf("READING DATA\n");
  }

  i = 0;
  *next = 0;
  nsub = 0;
  *nattrib = 0;
  int n = 0;
  while (i < (int)ndatasets[0]) {
      H5Gget_objname_by_idx(file_id, i, name, 256);
      /*printf(" dataset %s in file \n",name); */
      if (sscanf(name, "GRHYDRO::%s it=%i tl=%i rl=%i c=%i",
            field, &it, &tl, &level, &cnum)==5) {
          if (it == iter) {
             nsub++;
             printf("%s it=%i tl=%i rl=%i cnum=%i ",field,it,tl,level,cnum);
             ierr = read_cactus_dataset(file_id,name,nattrib,ncells,ndim,&n,time,inheader);
          } else if (nsub > 0) {
             *next = it;
             break;;
          }
      }
      i++;
  }
  if (nsub > 0) nsteps = 1;
  printf("Read %s from %i timestep, it=%i ncells=%i (next=%i)\n",field,nsteps,iter,*ncells,*next);
  
  /* send dataset name back to phantom */
  i = 6;
  set_blocklabel(&i,field);

  return nsteps;
}

int read_cactus_dataset(hid_t file_id,char *name,int *nattrib,int *ncells,int *ndim,int *n,double *time,int inheader)
{
  hid_t dataspace_id,dataset_id,attrib_id;
  herr_t    status;
  herr_t    HDF5_error = -1;

#if H5_VERSION_GE(1,8,0)
  dataset_id   = H5Dopen2(file_id,name,H5P_DEFAULT);
#else
  dataset_id   = H5Dopen(file_id,name);
#endif
  dataspace_id = H5Dget_space(dataset_id);
  /* get dimensionality */
  int nx,ny,nz;
  get_ndim_ncells(dataspace_id,ndim,&nx,&ny,&nz);
  *ncells = *ncells + nx*ny*nz;

  /* get number of attributes */
  *nattrib = H5Aget_num_attrs(dataset_id);

  printf("ndim=%i size=%i %i %i nattrib=%i\n",*ndim,nx,ny,nz,*nattrib);

  int i,isize;
  double origin[3],delta[3];
  int bbox[3],iorigin[3];
  int nghost[3];
  char name_attr[256];

  /*
   * Read through all of the attributes in the header, so we
   * can still spit out the values even if they are not used by SPLASH
   */
  for(i=0; i < *nattrib; i++) {
     attrib_id = H5Aopen_idx(dataset_id,i);
     ssize_t  attr_status;
     attr_status = H5Aget_name(attrib_id, 256, name_attr);

     hid_t  type_id;
     type_id = H5Aget_type(attrib_id);
     
     isize = (int) H5Aget_storage_size(attrib_id);
     /*type_class = H5Tget_native_type(type_id,H5T_DIR_ASCEND);*/
     if (strcmp(name_attr,"time")==0) {
        status = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,time);
     } else if (strcmp(name_attr,"origin")==0) { 
        status = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,&origin);
     } else if (strcmp(name_attr,"delta")==0) {
        status = H5Aread(attrib_id,H5T_NATIVE_DOUBLE,&delta);
     } else if (strcmp(name_attr,"cctk_bbox")==0) {
        status = H5Aread(attrib_id,H5T_NATIVE_INT,&bbox);
     } else if (strcmp(name_attr,"iorigin")==0) {
        status = H5Aread(attrib_id,H5T_NATIVE_INT,&iorigin);
     } else if (strcmp(name_attr,"cctk_nghostzones")==0) {
        status = H5Aread(attrib_id,H5T_NATIVE_INT,&nghost);
     } else {
        if (debug && inheader) printf("DEBUG: unknown attribute %s (%i bytes)\n",name_attr,isize);
     }

     if (status==HDF5_error) {
        printf(" ERROR reading attribute %s \n",name_attr);
     }
     status = H5Aclose(attrib_id);
  }
     printf("bbox: %i %i %i\n",bbox[0],bbox[1],bbox[2]);
     printf("iorigin: %i %i %i\n",iorigin[0],iorigin[1],iorigin[2]);
     printf("origin: %e %e %e\n",origin[0],origin[1],origin[2]);
     printf("delta: %e %e %e\n",delta[0],delta[1],delta[2]);
     printf("Time=%e Delta=%e %e %e nghost=%i %i %i\n",*time,delta[0],delta[1],delta[2],
            nghost[0],nghost[1],nghost[2]);

  if (inheader==1) {
     printf("Time=%e Delta=%e %e %e nghost=%i %i %i\n",*time,delta[0],delta[1],delta[2],
            nghost[0],nghost[1],nghost[2]);
  } else {  
    /* read dataset */ 
    read_cactus_grid(dataset_id,dataspace_id,*ndim,n,nx,ny,nz,nghost,iorigin,origin,delta);
  }

  /* close dataspace and dataset */
  status = H5Sclose(dataspace_id);
  if (status == HDF5_error) { 
     printf("ERROR closing dataspace \n");
  }

  status = H5Dclose(dataset_id);
  if (status == HDF5_error) { 
     printf("ERROR closing data set \n"); 
  }
  return (int) status;
}

/*
 *  utility function to get dimensionality of a dataset
 */
void get_ndim_ncells(hid_t dataspace_id, int *ndim, int *nx,int *ny,int *nz)
{
  *ndim    = H5Sget_simple_extent_ndims(dataspace_id);
  hsize_t dims[*ndim], maxdims[*ndim];
  H5Sget_simple_extent_dims(dataspace_id,dims,maxdims);

  /* get number of cells */
  *nx = dims[0];
  *ny = dims[1];
  *nz = dims[2];
}

int read_cactus_grid(hid_t dataset_id,hid_t dataspace_id,int ndim,int *n,int nx,int ny,int nz,
                     int nghost[3],int ioffset[3],double orig[3],double delta[3]) 
{
   hid_t     memspace_id;
   herr_t    status;
   herr_t    HDF5_error = -1;
   int       ierr = 0;

   int i,j,k;
   int ncells = nx*ny*nz;

   /* make a temporary array to put each column as we read it */
   hsize_t nsize[3];
   nsize[0] = nx;
   nsize[1] = ny;
   nsize[2] = nz;
   memspace_id = H5Screate_simple(ndim,nsize,NULL);
   double temp[nx][ny][nz];
   double dat[ncells],xx[ncells],yy[ncells],zz[ncells];
   int itype[ncells];

   status = H5Dread(dataset_id,H5T_NATIVE_DOUBLE,memspace_id,dataspace_id,H5P_DEFAULT,temp);
   if (status == HDF5_error) { printf("ERROR reading dataset \n"); ierr = 4; }

   /* map to one dimensional data */
   double xi,yi,zi;
   double dz = delta[2];
   double dy = dz;
   double dx = dz;
   
   int ip = 0;
   
   for (k=0;k<nz;k++) {
       zi = orig[2] + k*dz;
       for (j=0;j<ny;j++) { 
           yi = orig[1] + j*dy;
           for (i=0;i<nx;i++) {             
               xi = orig[0] + i*dx;
               dat[ip]=temp[i][j][k];
               xx[ip]=xi;
               yy[ip]=yi;
               zz[ip]=zi;
               /*printf(" rho %i %i %i = %e \n",i,j,k,temp[i][j][k]);*/
               *n = *n + 1;
               ip++;
           }
       }
   }
   
   /* tag ghost particles */
   ip = 0;
   int isghostx,isghosty,isghostz;
   for (k=0;k<nz;k++) {
       if (k < nghost[2] || k > nz-nghost[2]) { 
          isghostz = 1; 
       } else {
          isghostz = 0;
       }   
       for (j=0;j<ny;j++) {        
           if (j < nghost[1] || j > ny-nghost[1]) {
              isghosty = 1;
           } else {
              isghosty = 0;
           }
           for (i=0;i<nx;i++) {  
               if (i < nghost[1] || i > ny-nghost[1]) {
                  isghostx = 1; 
               } else {
                  isghostx = 0;
               }                  
               if (isghostx || isghosty || isghostz) {               
                 itype[ip] = 2;
               } else {
                 itype[ip] = 1;
               }
               ip++;
           }
       }
   }
   
   /* send data through to Fortran */
   int icol = 1;
   read_cactus_hdf5_data_fromc(&icol,n,&ncells,xx);
   icol = 2;
   read_cactus_hdf5_data_fromc(&icol,n,&ncells,yy);
   icol = 3;
   read_cactus_hdf5_data_fromc(&icol,n,&ncells,zz);
   /* smoothing length */
   for (i=0;i<ncells;i++) {
       xx[i] = dx;
   }
   icol = 4;
   read_cactus_hdf5_data_fromc(&icol,n,&ncells,xx);
   /* mass */
   for (i=0;i<ncells;i++) {
       xx[i] = dat[i]/(dx*dy*dz);
   }
   icol = 5;
   read_cactus_hdf5_data_fromc(&icol,n,&ncells,xx);
   icol = 6;
   read_cactus_hdf5_data_fromc(&icol,n,&ncells,dat);

   /* zone type (regular or ghost) */
   read_cactus_itype_fromc(n,&ncells,itype);

   return ierr;
}
