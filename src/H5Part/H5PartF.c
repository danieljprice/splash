/*
 * This is my modified H5PartF.c for use in SPLASH
 * Main changes are that the strings are now passed
 * directly from Fortran, already with the null character appended
 *
 * Modified from the original H5Part Fortran interface
 * by Daniel Price 08/04/10
 */
#include "H5Part.h"
#include <string.h>

char
_H5Part_flagsfor2c (
	char * flags
	) {

	char fbits = 0x00;

	flags = strtok ( flags, "," );
	while ( flags != NULL ) {
		if ( strcmp ( flags, "vfd_mpiposix" ) == 0 )
				fbits |= H5PART_VFD_MPIPOSIX;
		else if ( strcmp ( flags, "fs_lustre" ) == 0 )
		    		fbits |= H5PART_FS_LUSTRE;
		flags = strtok ( NULL, "," );
	}

	return fbits;
}

/* open/close interface */
h5part_int64_t
h5ptc_openr (
	const char *file_name
	) {

	H5PartFile* f = H5PartOpenFile ( file_name, H5PART_READ );

	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5ptc_openw (
	const char *file_name
	) {

	H5PartFile* f = H5PartOpenFile ( file_name, H5PART_WRITE );

	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5ptc_opena (
	const char *file_name
	) {
	
	H5PartFile* f = H5PartOpenFile ( file_name, H5PART_APPEND );

	return (h5part_int64_t)(size_t)f;
}

h5part_int64_t
h5ptc_openr_align (
	const char *file_name,
	const h5part_int64_t *align
	) {

	H5PartFile* f = H5PartOpenFileAlign ( file_name, H5PART_READ, *align );

	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5ptc_openw_align (
	const char *file_name,
	const h5part_int64_t *align
	) {

	H5PartFile* f = H5PartOpenFileAlign ( file_name, H5PART_WRITE, *align );

	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5ptc_opena_align (
	const char *file_name,
	const h5part_int64_t *align
	) {
	
	H5PartFile* f = H5PartOpenFileAlign ( file_name, H5PART_APPEND, *align );

	return (h5part_int64_t)(size_t)f;
}

#ifdef PARALLEL_IO
h5part_int64_t
h5ptc_openr_par (
	const char *file_name,
	MPI_Fint *fcomm
	) {

	MPI_Comm ccomm = MPI_Comm_f2c (*fcomm);

	H5PartFile* f = H5PartOpenFileParallel (
		file_name, H5PART_READ, ccomm );

	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5ptc_openw_par (
	const char *file_name,
	MPI_Fint *fcomm
	) {

	MPI_Comm ccomm = MPI_Comm_f2c (*fcomm);

	H5PartFile* f = H5PartOpenFileParallel (
		file_name, H5PART_WRITE, ccomm );

	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5ptc_opena_par (
	const char *file_name,
	MPI_Fint *fcomm
	) {
	
	MPI_Comm ccomm = MPI_Comm_f2c (*fcomm);
       
	H5PartFile* f = H5PartOpenFileParallel (
	       file_name, H5PART_APPEND, ccomm );
       
	return (h5part_int64_t)(size_t)f;
}

h5part_int64_t
h5ptc_openr_par_align (
	const char *file_name,
	MPI_Fint *fcomm,
	const h5part_int64_t *align
	) {

	MPI_Comm ccomm = MPI_Comm_f2c (*fcomm);

	H5PartFile* f = H5PartOpenFileParallelAlign (
		file_name, H5PART_READ, ccomm, *align );

	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5ptc_openw_par_align (
	const char *file_name,
	MPI_Fint *fcomm,
	const h5part_int64_t *align,
	const char *flags
	) {

	MPI_Comm ccomm = MPI_Comm_f2c (*fcomm);

	char fbits = H5PART_WRITE | _H5Part_flagsfor2c ( flags );

	H5PartFile* f = H5PartOpenFileParallelAlign (
		file_name, fbits, ccomm, *align );

	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5ptc_opena_par_align (
	const char *file_name,
	MPI_Fint *fcomm,
	const h5part_int64_t *align,
	const char *flags
	) {
	
	MPI_Comm ccomm = MPI_Comm_f2c (*fcomm);
       
	char fbits = H5PART_APPEND | _H5Part_flagsfor2c ( flags );

	H5PartFile* f = H5PartOpenFileParallelAlign (
		file_name, fbits, ccomm, *align );

	return (h5part_int64_t)(size_t)f;
}
#endif

h5part_int64_t
h5ptc_close (
	const h5part_int64_t *f
	) {
	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartCloseFile ( filehandle );
}

/*==============Writing and Setting Dataset info========*/

h5part_int64_t
h5ptc_readstep (
	const h5part_int64_t *f,
	const h5part_int64_t *step,
	h5part_float64_t *x,
	h5part_float64_t *y,
	h5part_float64_t *z,
	h5part_float64_t *px,
	h5part_float64_t *py,
	h5part_float64_t *pz,
	h5part_int64_t *id
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartReadParticleStep (
		filehandle,(*step)-1,x,y,z,px,py,pz,id);
}


h5part_int64_t
h5ptc_setnpoints (
	const h5part_int64_t *f,
	h5part_int64_t *np
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartSetNumParticles ( filehandle, *np );
}

h5part_int64_t
h5ptc_setnpoints_strided (
	const h5part_int64_t *f,
	h5part_int64_t *np,
        h5part_int64_t *stride
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartSetNumParticlesStrided ( filehandle, *np, *stride );
}

h5part_int64_t
h5ptc_setstep (
	const h5part_int64_t *f,
	h5part_int64_t *step ) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartSetStep ( filehandle, (*step)-1 );
}

h5part_int64_t
h5ptc_writedata_r8 (
	const h5part_int64_t *f,
	const char *name,
	const h5part_float64_t *data ) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteDataFloat64 (
		filehandle, name, data );

	return herr;
}

h5part_int64_t
h5ptc_writedata_r4 (
	const h5part_int64_t *f,
	const char *name,
	const h5part_float32_t *data ) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteDataFloat32 (
		filehandle, name, data );

	return herr;
}

h5part_int64_t
h5ptc_writedata_i8 (
	const h5part_int64_t *f,
	const char *name,
	const h5part_int64_t *data ) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteDataInt64 (
		filehandle, name, data );

	return herr;
}

h5part_int64_t
h5ptc_writedata_i4 (
	const h5part_int64_t *f,
	const char *name,
	const h5part_int32_t *data ) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteDataInt32 (
		filehandle, name, data );

	return herr;
}

/*==============Reading Data Characteristics============*/

h5part_int64_t
h5ptc_getnsteps (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartGetNumSteps ( filehandle );
}

h5part_int64_t
h5ptc_getndatasets (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartGetNumDatasets ( filehandle );
}

h5part_int64_t
h5ptc_getnpoints (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartGetNumParticles ( filehandle );
}

h5part_int64_t
h5ptc_getdatasetname ( 
	const h5part_int64_t *f,
	const h5part_int64_t *index,
	char *name,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr =  H5PartGetDatasetName (
		filehandle, *index, name, l_name );

	return herr;
}

h5part_int64_t
h5ptc_getdatasetinfo ( 
	const h5part_int64_t *f,
	const h5part_int64_t *index,
	char *name,
        h5part_int64_t *type,
        h5part_int64_t *nelem,
	const h5part_int64_t *l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;
        h5part_int64_t type_tmp;
	h5part_int64_t herr =  H5PartGetDatasetInfo (
		filehandle, *index, name, *l_name, &type_tmp, nelem );
        
        if (type_tmp == H5PART_INT64)
            {
               *type = 1;
            }
        else if (type_tmp == H5PART_INT32)
            {
               *type = 2;
            }
        else if (type_tmp == H5PART_FLOAT64)
            {
               *type = 3;
            }
        else if (type_tmp == H5PART_FLOAT32)
            {
               *type = 4;
            }
        else if (type_tmp == H5PART_CHAR)
            {
               *type = 5;
            }
        else if (type_tmp == H5PART_STRING)
            {
               *type = 6;
            }
        else
            {
               *type = 0;
            }

	return herr;
}

/*=============Setting and getting views================*/

h5part_int64_t
h5ptc_setview (
	const h5part_int64_t *f,
	const h5part_int64_t *start,
	const h5part_int64_t *end
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartSetView ( filehandle, *start, *end );
}

h5part_int64_t
h5ptc_setview_indices (
	const h5part_int64_t *f,
	const h5part_int64_t *indices,
	const h5part_int64_t *nelem
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartSetViewIndices ( filehandle, indices, *nelem );
}

h5part_int64_t
h5ptc_resetview (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartResetView ( filehandle );
}

h5part_int64_t
h5ptc_hasview (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartHasView ( filehandle );
}

h5part_int64_t
h5ptc_getview (
	const h5part_int64_t *f,
	h5part_int64_t *start,
	h5part_int64_t *end
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartGetView ( filehandle, start, end);
}
/*==================Reading data ============*/
h5part_int64_t
h5ptc_readdata_r8 (
	const h5part_int64_t *f,
	const char *name,
	h5part_float64_t *array
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadDataFloat64 (
		filehandle, name, array );

	return herr;
}

h5part_int64_t
h5ptc_readdata_r4 (
	const h5part_int64_t *f,
	const char *name,
	h5part_float32_t *array
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadDataFloat32 (
		filehandle, name, array );

	return herr;
}

h5part_int64_t
h5ptc_readdata_i8 (
	const h5part_int64_t *f,
	const char *name,
	h5part_int64_t *array
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadDataInt64 (
		filehandle, name, array );

	return herr;
}

h5part_int64_t
h5ptc_readdata_i4 (
	const h5part_int64_t *f,
	const char *name,
	h5part_int32_t *array
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadDataInt32 (
		filehandle, name, array );

	return herr;
}

/*
h5part_int64_t
h5ptc_set_verbosity_level (
	const h5part_int64_t *level
	) {
	return H5PartSetVerbosityLevel ( *level );
}
*/
