/*
 * This subroutine performs the calls to the PBOB c routines
 *
 */
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <sph_types.h>
#include <pbob.h>
#include <particle.h>

extern PBOB ReadPBOB(char *file_name);
extern PARTICLE *ReadParticle(char *file_name,int time_slice,int start_indx,int N);

void set_blocklabel(int *icol, char *name);
void read_pbob_data_fromc(int icol, int np, double temparr[np],int itype[np], char *tag);
void read_pbob_header(char *filename,
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
      *ndim  = 2;
      *ndimV = 2;
      *ncol  = 15; /* hard wired */
      PBOB pbob;
      pbob = ReadPBOB(filename);
      /*if (pbob == NULL) {
         *ierr = 1;
         return;
      }*/
      if (pbob.endian_int==999) {
         *ierr = 2;
         return;
      }
      /*printf(" cluster_size    = %i \n",pbob.cluster_size);
      printf(" length_units    = %s \n",pbob.length_units);
      printf(" mass_units      = %s \n",pbob.mass_units);
      printf(" time_units      = %s \n",pbob.time_units);
      printf(" energy_units    = %s \n",pbob.internal_energy_units);
      */
      printf(" short title     = %s \n",pbob.short_title);
      printf(" total_particles = %llu \n",pbob.total_particles);
      printf(" nn_k            = %i \n",pbob.nn_k);
      /*
      printf(" offset          = %i \n",pbob.first_particle_byte_offset);
      printf(" length          = %i \n",pbob.particle_length_bytes);
      */
      printf(" endian          = %s \n",pbob.endian_str);
      printf(" version         = %s \n",pbob.pbob_version);
      
      *npart = pbob.total_particles;
      *time  = pbob.time;
/*
      printf(" ascii_header: \n");
      printf("%s\n",pbob.ascii_header);
*/
   }

void read_pbob_data(char *filename,
                    int npart,
                    int time_slice,
                    int *ierr)
   {

      PARTICLE *particle = NULL;
      *ierr = 0;
      int start_indx = 0;
      int N = npart;
      int i;
      /*printf(" time_slice = %i\n",time_slice);*/
      if ((particle=ReadParticle(filename,time_slice,start_indx,N))==NULL)
      {
        *ierr = 2;
        return;
      }

      int species[npart];
      for (i=0;i<N;i++) { 
          switch(particle[i].species)
          {
          case 2048:
             species[i] = 3;
             break;
          case 1024:
             species[i] = 2;
             break;
          default:
             species[i] = 1;
          }
      }

      /*
        printf("i = %7i x = %20.10f index = %i species = %i node=%i\n",i,(float)particle[i].x,(int)particle[i].index,(int)particle[i].species,(int)particle[i].start_node);
       */
      double *temp = 0;
      temp = malloc(npart*sizeof(double));
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].x; }
      read_pbob_data_fromc(1,N,temp,species,"x");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].z; }
      read_pbob_data_fromc(2,N,temp,species,"z");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].y; }
      read_pbob_data_fromc(3,N,temp,species,"y");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].vx; }
      read_pbob_data_fromc(4,N,temp,species,"vx");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].vz; }
      read_pbob_data_fromc(5,N,temp,species,"vz");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].vy; }
      read_pbob_data_fromc(6,N,temp,species,"vy");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].ax; }
      read_pbob_data_fromc(7,N,temp,species,"ax");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].az; }
      read_pbob_data_fromc(8,N,temp,species,"az");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].ay; }
      read_pbob_data_fromc(9,N,temp,species,"ay");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].U; }
      read_pbob_data_fromc(10,N,temp,species,"U");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].dUdt; }
      read_pbob_data_fromc(11,N,temp,species,"dUdt");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].h; }
      read_pbob_data_fromc(12,N,temp,species,"h");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].rho; }
      read_pbob_data_fromc(13,N,temp,species,"rho");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].p; }
      read_pbob_data_fromc(14,N,temp,species,"p");
      for (i=0;i<N;i++) { temp[i] = (double) particle[i].m; }
      read_pbob_data_fromc(15,N,temp,species,"m");
      free(temp);

      return;
   }

