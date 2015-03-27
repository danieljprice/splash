/*
 * This subroutine performs the calls to the SILO library for the
 * SILO data read
 *
 * We have to do it this way as the SILO read interface for Fortran
 * is incomplete
 *
 */
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <silo.h>
static int debug = 0;

void set_blocklabel(int *icol, char *name);
void read_silo_data_fromc(int *icol, int *npartoftypei, double temparr[*npartoftypei],int *itype);
void read_silo_header(const char *filename,
                             int *npart,
                             int *ncol,
                             int *ndim,
                             int *ndimV,
                             double *time,
                             int *ierr)
   {
      *npart = 0;
      *ierr = 0;
      *ncol = 0;
      *time = 0.;
      *npart = 0;
      DBfile *silofile = DBOpen(filename, DB_UNKNOWN, DB_READ);
      if (!silofile) {
         *ierr = 1;
         return;
      }

      if (!DBVersionGEFileVersion(silofile)) {
         const char *siloversion = DBFileVersion(silofile);
         printf(" WARNING! File was created with newer version (v%s) of silo library\n",siloversion);
      }
      
      if (!DBInqFileHasObjects(silofile)) {
         printf(" ERROR: silo file %s does not appear to contain any objects\n",filename);
         *ierr = 2;
         return;
      }

      DBtoc *toc = DBGetToc(silofile);
      if (debug) {
         printf(" DEBUG: File contains:\n %i curves\n %i multimesh\n %i nmultimeshadj\n %i multivar\n", \
                toc->ncurve,toc->nmultimesh,toc->nmultimeshadj,toc->nmultivar);
         printf(" %i multimat\n %i multimatspecies\n %i csgmesh\n %i csgvar\n", \
                toc->nmultimat,toc->nmultimatspecies,toc->ncsgmesh,toc->ncsgvar);
         printf(" %i defvars\n %i qmesh\n %i qvar\n %i ucdmesh\n %i ucdvar\n", \
                toc->ndefvars,toc->nqmesh,toc->nqvar,toc->nucdmesh,toc->nucdvar);
         printf(" %i ptmesh\n %i ptvar\n %i mat\n %i matspecies\n %i var\n %i obj\n", \
                toc->nptmesh,toc->nptvar,toc->nmat,toc->nmatspecies,toc->nvar,toc->nobj);
         printf(" %i dir\n %i array\n %i mrgtree\n %i groupelmap\n %i mrgvar\n", \
                toc->ndir,toc->narray,toc->nmrgtree,toc->ngroupelmap,toc->nmrgvar);
      }
      int nptmesh = toc->nptmesh;
      if (nptmesh <= 0) {
         printf(" ERROR: silo file %s does not appear to contain any point meshes\n",filename);
         *ierr = 3;
         return;
      }
      int i;
      for (i=0;i<nptmesh;i++) {
          if (i > 0) {
             printf(" WARNING: IGNORNING ptmesh #%i (%s)\n",i+1,toc->ptmesh_names[i]);
          }
      }

      /* open the first point mesh and get info */
      DBpointmesh *my_ptmesh = DBGetPointmesh(silofile,toc->ptmesh_names[0]);
      if (!my_ptmesh) {
         printf(" ERROR reading point mesh %s\n",toc->ptmesh_names[0]);
         *ierr = 4;
         return;
      } else {
         printf(" Reading point mesh %s\n",toc->ptmesh_names[0]);
         *time  = my_ptmesh->dtime;
         *npart = my_ptmesh->nels;
         *ndim  = my_ptmesh->ndims;
         *ndimV = *ndim;
         if (debug) {
            printf(" Got labels = %s %s %s \n",\
                   my_ptmesh->labels[0],my_ptmesh->labels[1],my_ptmesh->labels[2]);
            printf(" Got title = %s \n",my_ptmesh->title);
            printf(" Got units = %s %s %s \n",\
                   my_ptmesh->units[0],my_ptmesh->units[1],my_ptmesh->units[2]);
            printf(" max_extents = %f %f %f\n",\
                  my_ptmesh->max_extents[0],my_ptmesh->max_extents[1],my_ptmesh->max_extents[2]);
         }
         DBFreePointmesh(my_ptmesh);
      }

      /* Read the other point variables */
      int nptvar = toc->nptvar;
      if (nptvar <= 0) {
         printf(" WARNING: silo file %s does not appear to contain any point variables\n",filename);
      }
      *ncol = *ndim + nptvar;
      if (*ncol <= 0) {
         *ierr = 4;
         printf(" ERROR: ncol <= 0 from silo header\n");
         return;
      }
      DBClose(silofile);
   }

void read_silo_data(char *filename,
                    int maxtypes,
                    int npartoftype[maxtypes],
                    int ncol,
                    int isrequired[ncol],
                    int *ierr)
   {

      DBfile *silofile = DBOpen(filename, DB_UNKNOWN, DB_READ);
      if (!silofile) {
         *ierr = 1;
         return;
      }
      DBtoc *toc = DBGetToc(silofile);
      DBpointmesh *my_ptmesh = DBGetPointmesh(silofile,toc->ptmesh_names[0]);
      if (!my_ptmesh) {
         printf(" ERROR reading point mesh %s\n",toc->ptmesh_names[0]);
         DBClose(silofile);
         *ierr = 4;
         return;
      } else {
         if (debug) printf(" DEBUG: Reading data from point mesh %s\n",toc->ptmesh_names[0]);
         int ndim = my_ptmesh->ndims;
         /*int nels = my_ptmesh->nels;*/
         int i;
         int np = npartoftype[0];
         int idim;
         int particle_type = 1;
         
         /* read in float */
         float *x = 0;
         x = malloc(np*sizeof(float));
         
         /* must send double to splash */
         double *tmp_dbl = 0;
         tmp_dbl = malloc(np*sizeof(double));
         
         /* read each coordinate in turn */
         int icol = 0;
         for (idim=0;idim<ndim;idim++) {         
             icol = icol + 1;
             x = my_ptmesh->coords[idim];
             for (i=0;i<np;i++) {
                 tmp_dbl[i] = (double) x[i];
             }
             read_silo_data_fromc(&icol,&np,tmp_dbl,&particle_type);
             switch(icol) {
             case 1:
                set_blocklabel(&icol,"x");
                break;
             case 2:
                set_blocklabel(&icol,"y");
                break;
             case 3:
                set_blocklabel(&icol,"z");
                break;
             }
         }

         int nptvar = toc->nptvar;
         for (i=0;i<nptvar;i++) {
             printf(" Reading ptvar #%i: %s\n",i+1,toc->ptvar_names[i]);
             DBmeshvar *my_ptvar = DBGetPointvar(silofile,toc->ptvar_names[i]);
             if (debug) {
                printf(" Associated with mesh %s\n",my_ptvar->meshname);
                printf(" label = %s\n",my_ptvar->label);
                printf(" Nels,nvals,nspace,ndims = %i %i %i %i\n",\
                     my_ptvar->nels,my_ptvar->nvals,my_ptvar->nspace,my_ptvar->ndims);
             }
             icol = icol + 1;
             x = my_ptvar->vals[0];
             for (i=0;i<np;i++) {
                 tmp_dbl[i] = (double) x[i];
             }
             read_silo_data_fromc(&icol,&np,tmp_dbl,&particle_type);
             set_blocklabel(&icol,my_ptvar->name);
             /*DBFreeMeshvar(my_ptvar);*/
         }

         free(x);
         free(tmp_dbl);

         DBFreePointmesh(my_ptmesh);
      }
      DBClose(silofile);
   }

