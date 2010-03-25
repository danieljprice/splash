#include "H5Part.h"
#include <string.h>
#include <hdf5.h>

char *
_H5Part_strdupfor2c (
	const char *s,
	const ssize_t len
	) {

	char *dup = (char*)malloc ( len + 1 );
	strncpy ( dup, s, len );
	char *p = dup + len;
	do {
		*p-- = '\0';
	} while ( *p == ' ' );
	return dup;
}

char *
_H5Part_strc2for (
	char * const str,
	const ssize_t l_str
	) {

	size_t len = strlen ( str );
	memset ( str+len, ' ', l_str-len );

	return str;
}

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
h5pt_openr (
	const char *file_name,
	const int l_file_name
	) {

	char *file_name2 = _H5Part_strdupfor2c ( file_name, l_file_name );

	H5PartFile* f = H5PartOpenFile ( file_name2, H5PART_READ );

	free ( file_name2 );
	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5pt_openw (
	const char *file_name,
	const int l_file_name
	) {

	char *file_name2 = _H5Part_strdupfor2c ( file_name, l_file_name );

	H5PartFile* f = H5PartOpenFile ( file_name2, H5PART_WRITE );

	free ( file_name2 );
	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5pt_opena (
	const char *file_name,
	const int l_file_name
	) {
	
	char *file_name2 = _H5Part_strdupfor2c ( file_name, l_file_name );

	H5PartFile* f = H5PartOpenFile ( file_name2, H5PART_APPEND );

	free ( file_name2 );
	return (h5part_int64_t)(size_t)f;
}

h5part_int64_t
h5pt_openr_align (
	const char *file_name,
	const h5part_int64_t *align,
	const int l_file_name
	) {

	char *file_name2 = _H5Part_strdupfor2c ( file_name, l_file_name );

	H5PartFile* f = H5PartOpenFileAlign ( file_name2, H5PART_READ, *align );

	free ( file_name2 );
	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5pt_openw_align (
	const char *file_name,
	const h5part_int64_t *align,
	const int l_file_name
	) {

	char *file_name2 = _H5Part_strdupfor2c ( file_name, l_file_name );

	H5PartFile* f = H5PartOpenFileAlign ( file_name2, H5PART_WRITE, *align );

	free ( file_name2 );
	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5pt_opena_align (
	const char *file_name,
	const h5part_int64_t *align,
	const int l_file_name
	) {
	
	char *file_name2 = _H5Part_strdupfor2c ( file_name, l_file_name );

	H5PartFile* f = H5PartOpenFileAlign ( file_name2, H5PART_APPEND, *align );

	free ( file_name2 );
	return (h5part_int64_t)(size_t)f;
}

#ifdef PARALLEL_IO
h5part_int64_t
h5pt_openr_par (
	const char *file_name,
	MPI_Fint *fcomm,
	const int l_file_name
	) {

	MPI_Comm ccomm = MPI_Comm_f2c (*fcomm);

	char *file_name2 = _H5Part_strdupfor2c ( file_name, l_file_name );

	H5PartFile* f = H5PartOpenFileParallel (
		file_name2, H5PART_READ, ccomm );

	free ( file_name2 );
	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5pt_openw_par (
	const char *file_name,
	MPI_Fint *fcomm,
	const int l_file_name
	) {

	MPI_Comm ccomm = MPI_Comm_f2c (*fcomm);
	char *file_name2 = _H5Part_strdupfor2c ( file_name, l_file_name );

	H5PartFile* f = H5PartOpenFileParallel (
		file_name2, H5PART_WRITE, ccomm );

	free ( file_name2 );
	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5pt_opena_par (
	const char *file_name,
	MPI_Fint *fcomm,
	const int l_file_name
	) {
	
	MPI_Comm ccomm = MPI_Comm_f2c (*fcomm);
	char *file_name2 = _H5Part_strdupfor2c ( file_name, l_file_name );
       
	H5PartFile* f = H5PartOpenFileParallel (
	       file_name2, H5PART_APPEND, ccomm );
       
	free ( file_name2 );
	return (h5part_int64_t)(size_t)f;
}

h5part_int64_t
h5pt_openr_par_align (
	const char *file_name,
	MPI_Fint *fcomm,
	const h5part_int64_t *align,
	const int l_file_name
	) {

	MPI_Comm ccomm = MPI_Comm_f2c (*fcomm);
	char *file_name2 = _H5Part_strdupfor2c ( file_name, l_file_name );

	H5PartFile* f = H5PartOpenFileParallelAlign (
		file_name2, H5PART_READ, ccomm, *align );

	free ( file_name2 );
	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5pt_openw_par_align (
	const char *file_name,
	MPI_Fint *fcomm,
	const h5part_int64_t *align,
	const char *flags,
	const int l_file_name,
	const int l_flags
	) {

	MPI_Comm ccomm = MPI_Comm_f2c (*fcomm);
	char *file_name2 = _H5Part_strdupfor2c ( file_name, l_file_name );
	char *flags2 = _H5Part_strdupfor2c ( flags, l_flags );

	char fbits = H5PART_WRITE | _H5Part_flagsfor2c ( flags2 );

	H5PartFile* f = H5PartOpenFileParallelAlign (
		file_name2, fbits, ccomm, *align );

	free ( file_name2 );
	return (h5part_int64_t)(size_t)f; 
}

h5part_int64_t
h5pt_opena_par_align (
	const char *file_name,
	MPI_Fint *fcomm,
	const h5part_int64_t *align,
	const char *flags,
	const int l_file_name,
	const int l_flags
	) {
	
	MPI_Comm ccomm = MPI_Comm_f2c (*fcomm);
	char *file_name2 = _H5Part_strdupfor2c ( file_name, l_file_name );
	char *flags2 = _H5Part_strdupfor2c ( flags, l_flags );
       
	char fbits = H5PART_APPEND | _H5Part_flagsfor2c ( flags2 );

	H5PartFile* f = H5PartOpenFileParallelAlign (
		file_name2, fbits, ccomm, *align );

	free ( file_name2 );
	return (h5part_int64_t)(size_t)f;
}
#endif

h5part_int64_t
h5pt_close (
	const h5part_int64_t *f
	) {
	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartCloseFile ( filehandle );
}

/*==============Writing and Setting Dataset info========*/

h5part_int64_t
h5pt_readstep (
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
h5pt_setnpoints (
	const h5part_int64_t *f,
	h5part_int64_t *np
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartSetNumParticles ( filehandle, *np );
}

h5part_int64_t
h5pt_setnpoints_strided (
	const h5part_int64_t *f,
	h5part_int64_t *np,
        h5part_int64_t *stride
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartSetNumParticlesStrided ( filehandle, *np, *stride );
}

h5part_int64_t
h5pt_setstep (
	const h5part_int64_t *f,
	h5part_int64_t *step ) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartSetStep ( filehandle, (*step)-1 );
}

h5part_int64_t
h5pt_writedata_r8 (
	const h5part_int64_t *f,
	const char *name,
	const h5part_float64_t *data,
	const int l_name ) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 = _H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWriteDataFloat64 (
		filehandle, name2, data );

	free ( name2 );

	return herr;
}

h5part_int64_t
h5pt_writedata_r4 (
	const h5part_int64_t *f,
	const char *name,
	const h5part_float32_t *data,
	const int l_name ) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 = _H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWriteDataFloat32 (
		filehandle, name2, data );

	free ( name2 );

	return herr;
}

h5part_int64_t
h5pt_writedata_i8 (
	const h5part_int64_t *f,
	const char *name,
	const h5part_int64_t *data,
	const int l_name ) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 = _H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWriteDataInt64 (
		filehandle, name2, data );

	free ( name2 );

	return herr;
}

h5part_int64_t
h5pt_writedata_i4 (
	const h5part_int64_t *f,
	const char *name,
	const h5part_int32_t *data,
	const int l_name ) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 = _H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWriteDataInt32 (
		filehandle, name2, data );

	free ( name2 );

	return herr;
}

/*==============Reading Data Characteristics============*/

h5part_int64_t
h5pt_getnsteps (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartGetNumSteps ( filehandle );
}

h5part_int64_t
h5pt_getndatasets (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartGetNumDatasets ( filehandle );
}

h5part_int64_t
h5pt_getnpoints (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartGetNumParticles ( filehandle );
}

h5part_int64_t
h5pt_getdatasetname ( 
	const h5part_int64_t *f,
	const h5part_int64_t *index,
	char *name,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr =  H5PartGetDatasetName (
		filehandle, *index, name, l_name );

	_H5Part_strc2for ( name, l_name );
	return herr;
}

h5part_int64_t
h5pt_getdatasetinfo ( 
	const h5part_int64_t *f,
	const h5part_int64_t *index,
	char *name,
        h5part_int64_t *type,
        h5part_int64_t *nelem,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;
        h5part_int64_t type_tmp;
	h5part_int64_t herr =  H5PartGetDatasetInfo (
		filehandle, *index, name, l_name, &type_tmp, nelem );
        
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

	_H5Part_strc2for ( name, l_name );
	return herr;
}

/*=============Setting and getting views================*/

h5part_int64_t
h5pt_setview (
	const h5part_int64_t *f,
	const h5part_int64_t *start,
	const h5part_int64_t *end
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartSetView ( filehandle, *start, *end );
}

h5part_int64_t
h5pt_setview_indices (
	const h5part_int64_t *f,
	const h5part_int64_t *indices,
	const h5part_int64_t *nelem
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartSetViewIndices ( filehandle, indices, *nelem );
}

h5part_int64_t
h5pt_resetview (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartResetView ( filehandle );
}

h5part_int64_t
h5pt_hasview (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartHasView ( filehandle );
}

h5part_int64_t
h5pt_getview (
	const h5part_int64_t *f,
	h5part_int64_t *start,
	h5part_int64_t *end
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartGetView ( filehandle, start, end);
}
/*==================Reading data ============*/
h5part_int64_t
h5pt_readdata_r8 (
	const h5part_int64_t *f,
	const char *name,
	h5part_float64_t *array,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 = _H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartReadDataFloat64 (
		filehandle, name2, array );

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_readdata_r4 (
	const h5part_int64_t *f,
	const char *name,
	h5part_float32_t *array,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 = _H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartReadDataFloat32 (
		filehandle, name2, array );

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_readdata_i8 (
	const h5part_int64_t *f,
	const char *name,
	h5part_int64_t *array,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 = _H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartReadDataInt64 (
		filehandle, name2, array );

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_readdata_i4 (
	const h5part_int64_t *f,
	const char *name,
	h5part_int32_t *array,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 = _H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartReadDataInt32 (
		filehandle, name2, array );

	free ( name2 );
	return herr;
}

/*=================== Attributes ================*/

h5part_int64_t
h5pt_writefileattrib_string (
	const h5part_int64_t *f,
	const char *attrib_name,
	const char *attrib_value,
	const int l_attrib_name,
	const int l_attrib_value
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *attrib_name2 = _H5Part_strdupfor2c (attrib_name,l_attrib_name);
	char *attrib_value2= _H5Part_strdupfor2c (attrib_value,l_attrib_value);

	h5part_int64_t herr = H5PartWriteFileAttribString (
		filehandle, attrib_name2, attrib_value2 );

	free ( attrib_name2 );
	free ( attrib_value2 );
	return herr;
}

h5part_int64_t
h5pt_writestepattrib_string (
	const h5part_int64_t *f,
	const char *attrib_name,
	const char *attrib_value,
	const int l_attrib_name,
	const int l_attrib_value
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *attrib_name2 = _H5Part_strdupfor2c (attrib_name,l_attrib_name);
	char *attrib_value2= _H5Part_strdupfor2c (attrib_value,l_attrib_value);

	h5part_int64_t herr = H5PartWriteStepAttribString (
		filehandle, attrib_name2, attrib_value2 );

	free ( attrib_name2 );
	free ( attrib_value2 );
	return herr;
}

/* Reading attributes ************************* */

h5part_int64_t
h5pt_getnstepattribs (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartGetNumStepAttribs ( filehandle );
}

h5part_int64_t
h5pt_getnfileattribs (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartGetNumFileAttribs ( filehandle );
}

h5part_int64_t
h5pt_getstepattribinfo (
	const h5part_int64_t *f,
	const h5part_int64_t *idx,
	char *name,
	h5part_int64_t *nelem,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;
	h5part_int64_t type;

	h5part_int64_t herr = H5PartGetStepAttribInfo ( 
		filehandle, *idx, name, l_name, &type, nelem);

	_H5Part_strc2for( name, l_name );
	return herr;
}

h5part_int64_t
h5pt_getfileattribinfo (
	const h5part_int64_t *f,
	const h5part_int64_t *idx,
	char *name,
	h5part_int64_t *nelem,
	const int l_name ) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;
	h5part_int64_t type;

	h5part_int64_t herr = H5PartGetFileAttribInfo ( 
		filehandle, *idx, name, l_name, &type, nelem);

	_H5Part_strc2for( name, l_name );
	return herr;
}

h5part_int64_t
h5pt_readstepattrib_string (
	const h5part_int64_t *f,
	const char *attrib_name,
	char *attrib_value,
	const int l_attrib_name,
	const int l_attrib_value
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;
	
	char * attrib_name2 = _H5Part_strdupfor2c (attrib_name,l_attrib_name);

	h5part_int64_t herr = H5PartReadStepAttrib (
		filehandle, attrib_name2, attrib_value );

	_H5Part_strc2for ( attrib_value, l_attrib_value );

	free ( attrib_name2 );
	return herr;
}

h5part_int64_t
h5pt_readfileattrib_string (
	const h5part_int64_t *f,
	const char *attrib_name,
	char *attrib_value,
	const int l_attrib_name,
	const int l_attrib_value
	) {
		
	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char * attrib_name2 = _H5Part_strdupfor2c (attrib_name,l_attrib_name);

	h5part_int64_t herr = H5PartReadFileAttrib (
		filehandle, attrib_name2, attrib_value );

	_H5Part_strc2for ( attrib_value, l_attrib_value );

	free ( attrib_name2 );
	return herr;
}

h5part_int64_t
h5pt_set_verbosity_level (
	const h5part_int64_t *level
	) {
	return H5PartSetVerbosityLevel ( *level );
}
