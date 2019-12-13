#include <stdio.h>
#include <string.h>
#include "fitsio.h"

void read_fits_header(char *filename,
                    int  *npix,
                    int  *ncols,
                    int  *ierr)
{
   char card[FLEN_CARD];
   int status = 0,  nkeys, ii;  /* MUST initialize status */

   fitsfile *fptr;
   printf("%s\n",filename);
   fits_open_file(&fptr, filename, READONLY, &status);
   fits_get_hdrspace(fptr, &nkeys, NULL, &status);

   for (ii = 1; ii <= nkeys; ii++)  {
       fits_read_record(fptr, ii, card, &status); /* read keyword */
       printf("%s\n", card);
   }
   printf("END\n\n");  /* terminate listing with END */
     fits_close_file(fptr, &status);

     if (status)          /* print any error messages */
        fits_report_error(stderr, status);
     *ierr = status;
}

void read_fits_data(char *filename, int *ierr)
{
   fitsfile *fptr;
   printf("%s\n",filename);
   int status = 0;
   fits_open_file(&fptr, filename, READONLY, &status);
   printf("closing file\n");
}
