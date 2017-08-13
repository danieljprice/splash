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
#include <math.h>

#define MAX_DATASETS 100000
static int debug = 0;
static hid_t file_id;
static const herr_t HDF5_error = -1;
static int have_indexed = 0;

typedef struct
{
  char name[64];
} dataset_t;

static dataset_t Dataset[MAX_DATASETS];
static int iter[MAX_DATASETS];
static int iorder[MAX_DATASETS];

int read_cactus_dataset(hid_t file_id,char *name,int *nattrib,int *ncells,int *ndim,int *n,double *time,double *deltax,int inheader);
int read_cactus_iteration(hid_t file_id,int iter,int *next,int *nsteps,int *ncells,int *ndim,int *nattrib,double *time,double *deltax,int inheader);
int read_cactus_grid(hid_t dataset_id,hid_t dataspace_id,int ndim,int *n,int nx,int ny,int nz,int nghost[3],double orig[3],double delta[3]);
void get_ndim_ncells(hid_t dataspace_id, int *ndim, int *nx,int *ny,int *nz);
void set_blocklabel(int *icol, char *name);
void sort_cactus_data(int *n, int iarr[*n], int iorder[*n]);
void read_cactus_hdf5_data_fromc(int *icol,int *ntot,int *np,double temparr[*np]);
void read_cactus_itype_fromc(int *ntot,int *np,int itype[*np]);

void open_cactus_hdf5_file(char   *filename,
                             int istep,
                             int *npart,
                             int *ncol,
                             int *nsteps,
                             int *ndim,
                             int *ndimV,
                             double *time,
                             int *ierr)
{   
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
   int next;
   double dx;
   *ierr = read_cactus_iteration(file_id,istep,&next,nsteps,npart,ndim,&nattrib,time,&dx,1);
   *ierr = 0;

   *ndimV = *ndim;
   return;   
}
   
void close_cactus_hdf5_file(int *ierr)
{
   if (debug) printf("DEBUG: closing file \n");
   herr_t status = H5Fclose( file_id );
   if (status == HDF5_error) { printf("ERROR closing file \n"); *ierr = 1; }
   return;
}


void read_cactus_hdf5_data(char *filename,int istep,int *npart,double *time,double *dx,int *ierr)
{
   int next,nattrib,ndim,nstepsread,nsteps;

   if (file_id == HDF5_error)
      { printf("ERROR with file_id %s \n",filename); *ierr = 1; return; }

   nstepsread = read_cactus_iteration(file_id,istep,&next,&nsteps,npart,&ndim,&nattrib,time,dx,0);
   if (nstepsread <= 0) {
      *ierr = -1;
   }
}

/*
 *  read one iteration from the file, corresponding to all datasets with same iteration numbers
 */
int read_cactus_iteration(hid_t file_id,int istep,int *next,int *nsteps,int *ncells,int *ndim,int *nattrib,double *time,double *deltax,int inheader)
{
  hsize_t ndatasets[1];
  /* get number of datasets in file */
  H5Gget_num_objs(file_id, ndatasets);
  int it,tl,level,cnum;

  /* loop over all datasets looking for dataset matching the desired iteration number
     set function value to true (1) if it is present  */
  int i,j,nsub,ierr;
  char name[64];
  char thorn[24];
  int nstepsgot = 0;
  *ndim = 0;
  *ncells = 0;
  int n_datasets = (int)ndatasets[0];
  if (n_datasets > MAX_DATASETS) {
     if (inheader) printf("ERROR: EXCEEDED MAX NUMBER OF DATASETS: INCREASE MAX_DATASETS TO READ ALL\n");
     n_datasets = MAX_DATASETS;
  }
  if (inheader) {
     have_indexed = 0;
  }
  if (!have_indexed) {
     printf("Indexing %i datasets in file\n",n_datasets);
  }

  i = 0;
  *next = 0;
  nsub = 0;
  *nattrib = 0;
  *nsteps = 0;
  int n = 0;
  int mystep = 0;
  int iterprev = -1;
  while (i < n_datasets) {
      if (!have_indexed) {
         H5Gget_objname_by_idx(file_id, i, name, 64);
         strcpy(Dataset[i].name,name);
         if (sscanf(name, "%s it=%i tl=%i rl=%i c=%i",
               thorn, &it, &tl, &level, &cnum)>=5) {
            if (debug) printf("%s it=%i tl=%i rl=%i cnum=%i ",thorn,it,tl,level,cnum);
            iter[i] = it;
            /* send dataset name back to phantom */
            if (i==0) {
               j = 6;
               set_blocklabel(&j,thorn);
            }
         } else {
            iter[i] = -1;
         }
      } else {  /* have indexed file */
         strcpy(name,Dataset[iorder[i]-1].name);
         /*printf(" dataset %s in file \n",name);*/
         it = iter[iorder[i]-1];
      }
      if (it >= 0) {
          if (it != iterprev) {
              mystep++;
          }
          if (mystep == istep) {
             nsub++;
             ierr = read_cactus_dataset(file_id,name,nattrib,ncells,ndim,&n,time,deltax,inheader);
          } else if (mystep > istep && have_indexed) {
             break;;
          }
          iterprev = it;
      }
      i++;
  }
  have_indexed = 1;

  if (nsub > 0) {
     nstepsgot = 1;
     *nsteps = n_datasets/nsub;
  }
  if (inheader) {
     sort_cactus_data(&n_datasets,iter,iorder);
     printf("Read data from %i/%i timesteps, it=%i ncells=%i\n",istep,*nsteps,iterprev,*ncells);
  }
  
  return nstepsgot;
}

int read_cactus_dataset(hid_t file_id,char *name,int *nattrib,int *ncells,int *ndim,int *n,double *time,double *deltax,int inheader)
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

  if (debug) printf("ndim=%i size=%i %i %i nattrib=%i\n",*ndim,nx,ny,nz,*nattrib);

  int i,isize;
  double origin[3],delta[3];
  int bbox[6],iorigin[3];
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
  *deltax = delta[0];
  
  if (debug) {
     printf("bbox: x %i %i y %i %i z %i %i\n",bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);
     printf("iorigin: %i %i %i\n",iorigin[0],iorigin[1],iorigin[2]);
     printf("origin: %e %e %e\n",origin[0],origin[1],origin[2]);
     printf("delta: %e %e %e\n",delta[0],delta[1],delta[2]);
     printf("Time=%e Delta=%e %e %e nghost=%i %i %i \n",*time,delta[0],delta[1],delta[2],
            nghost[0],nghost[1],nghost[2]);
  }
  if (inheader==0) {
    /* read dataset */ 
    read_cactus_grid(dataset_id,dataspace_id,*ndim,n,nx,ny,nz,nghost,origin,delta);
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
                     int nghost[3],double orig[3],double delta[3]) 
{
   hid_t     memspace_id;
   herr_t    status;
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

   /* map to one dimensional arrays and determine x,y,z for each cell */
   double xi,yi,zi;
   double dx = delta[0];
   double dy = delta[1];
   double dz = delta[2];
   if (debug) {
      printf("PATCH origin = %f %f %f\n",orig[0],orig[1],orig[2]);
      printf(" nx=%i x=%f->%f\n",nx,orig[0],(orig[0]+(nx-1)*dx));
   }
   /* here we have to reconstruct the x, y and z positions of
      each cell. The following looks hacked but works. We tested it. */
   int ip = 0;
   for (k=0;k<nx;k++) {
       zi = orig[2] + k*dz;
       for (j=0;j<ny;j++) { 
           yi = orig[1] + j*dy;
           for (i=0;i<nz;i++) {
               xi = orig[0] + i*dx;
               dat[ip]=temp[k][j][i];
               xx[ip]=xi;
               yy[ip]=yi;
               zz[ip]=zi;
               /* printf("ip=%i rho %i %i %i = %f x=%f y=%f z=%f \n",ip,i,j,k,dat[ip],xi,yi,zi);*/
               *n = *n + 1;
               ip++;
           }
       }
   }
   
   /* tag ghost particles */
   ip = 0;
   int isghostx,isghosty,isghostz;
   for (k=0;k<nx;k++) {
       if ((k < nghost[2]) || (k > nx-nghost[2]-1)) { 
          isghostz = 1; 
       } else {
          isghostz = 0;
       }   
       for (j=0;j<ny;j++) {
           if ((j < nghost[1]) || (j > ny-nghost[1]-1)) {
              isghosty = 1;
           } else {
              isghosty = 0;
           }
           for (i=0;i<nz;i++) {  
               if ((i < nghost[0]) || (i > nz-nghost[0]-1)) {
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
   icol = 6;
   read_cactus_hdf5_data_fromc(&icol,n,&ncells,dat);

   /* zone type (regular or ghost) */
   read_cactus_itype_fromc(n,&ncells,itype);

   return ierr;
}
