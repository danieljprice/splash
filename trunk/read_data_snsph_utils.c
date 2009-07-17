#include <stdio.h>
#include <stdlib.h>
#include <error.h>
#include <string.h>

#include <SDF.h>
#include <SDFread.h>
#include <bigmalloc.h>

typedef struct {
    double x, y, z;
    float vx, vy, vz;
    float mass;
    float u;
    float rho;
    float h;
} body;


void readsdf_(char *filename, int *len, float *dat, int *maxpart, 
	      int *maxcol, int *ierr)
{
    int xconf, yconf, zconf, vxconf, vyconf, vzconf, mconf, uconf,
	rhoconf, hconf;
    int gnobj, nobj, i;
    body *p;
    SDF *sdfp;
    char fname[128];

    strncpy(fname, filename, *len);
    fname[*len] = '\0';

    if ( (sdfp = SDFopen(NULL, fname)) == NULL ) {
	fprintf(stderr, "%s: %s", fname, SDFerrstring);
	exit(2);
    }

    SDFread(sdfp, (void **)(&p), &gnobj, &nobj, sizeof(body), 
	    "x", offsetof(body, x), &xconf, 
	    "y", offsetof(body, y), &yconf,
	    "z", offsetof(body, z), &zconf,
	    "vx", offsetof(body, vx), &vxconf,
	    "vy", offsetof(body, vy), &vyconf,
	    "vz", offsetof(body, vz), &vzconf,
	    "mass", offsetof(body, mass), &mconf,
	    "u", offsetof(body, u), &uconf,
	    "rho", offsetof(body, rho), &rhoconf,
	    "h", offsetof(body, h), &hconf,
	    NULL);
    SDFclose(sdfp);

    if (!xconf || !yconf || !zconf || !vxconf || !vyconf || !vzconf ||
	!mconf || !uconf || !rhoconf || !hconf) {
	fprintf(stderr, "No %s%s%s%s%s%s%s%s%s%s in %s\n", 
		(xconf==0)? "x " : "",
		(yconf==0)? "y " : "",
		(zconf==0)? "z " : "",
		(vxconf==0)? "vx " : "",
		(vyconf==0)? "vy " : "",
		(vzconf==0)? "vz " : "",
		(mconf==0)? "m " : "",
		(uconf==0)? "u " : "",
		(rhoconf==0)? "rho " : "",
		(hconf==0)? "h " : "",
		fname);
	exit(3);
    }

    printf("nobj = %d\n", nobj);

    for(i = 0; i < nobj; ++i) {
	dat[i] = p[i].x;
	dat[(*maxpart)+i] = p[i].y;
	dat[2* (*maxpart)+i] = p[i].z;
	dat[3* (*maxpart)+i] = p[i].vx;
        dat[4* (*maxpart)+i] = p[i].vy;
        dat[5* (*maxpart)+i] = p[i].vz;
	dat[6* (*maxpart)+i] = p[i].mass;
	dat[7* (*maxpart)+i] = p[i].u;
	dat[8* (*maxpart)+i] = p[i].rho;
	dat[9* (*maxpart)+i] = p[i].h;
    }

    Free(p);

    *ierr = xconf || yconf || zconf || vxconf || vyconf || vzconf || mconf ||
	uconf || rhoconf || hconf;
}


int getcol_()
{
    return 10;
}


void getdata_(char *filename, int *len, int *nobj, float *tpos, float *gamma)
{
    SDF *sdfp;
    char fname[128];

    strncpy(fname, filename, *len);
    fname[*len] = '\0';

    if ( (sdfp = SDFopen(NULL, fname)) == NULL ) {
	fprintf(stderr, "%s: %s", fname, SDFerrstring);
	exit(2);
    }

    SDFgetintOrDie(sdfp, "npart", nobj);
    SDFgetfloatOrDie(sdfp, "tpos", tpos);
    SDFgetfloatOrDie(sdfp, "gamma", gamma);

    SDFclose(sdfp);
}
