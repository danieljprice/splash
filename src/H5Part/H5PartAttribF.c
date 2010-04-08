/*
 * This is my modified H5PartAttribF.c for use in SPLASH
 * Main changes are that the strings are now passed
 * directly from Fortran, already with the null character appended
 *
 * Modified from the original H5Part Fortran interface
 * by Daniel Price 08/04/10
 */
#include "H5Part.h"

h5part_int64_t
h5ptc_writefileattrib_r8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float64_t *data,
	const h5part_float64_t *nelem
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteFileAttrib (
		filehandle, name, H5PART_FLOAT64, data, *nelem);

	return herr;
}

h5part_int64_t
h5ptc_readfileattrib_r8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float64_t *data
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadFileAttrib (
		filehandle, name, (void*)data);

	return herr;
}

h5part_int64_t
h5ptc_writefileattrib_r4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float32_t *data,
	const h5part_float32_t *nelem
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteFileAttrib (
		filehandle, name, H5PART_FLOAT32, data, *nelem);

	return herr;
}

h5part_int64_t
h5ptc_readfileattrib_r4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float32_t *data
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadFileAttrib (
		filehandle, name, (void*)data);

	return herr;
}

h5part_int64_t
h5ptc_writefileattrib_i8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int64_t *data,
	const h5part_int64_t *nelem
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteFileAttrib (
		filehandle, name, H5PART_INT64, data, *nelem);

	return herr;
}

h5part_int64_t
h5ptc_readfileattrib_i8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int64_t *data
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadFileAttrib (
		filehandle, name, (void*)data);

	return herr;
}

h5part_int64_t
h5ptc_writefileattrib_i4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int32_t *data,
	const h5part_int32_t *nelem
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteFileAttrib (
		filehandle, name, H5PART_INT32, data, *nelem);

	return herr;
}

h5part_int64_t
h5ptc_readfileattrib_i4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int32_t *data
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadFileAttrib (
		filehandle, name, (void*)data);

	return herr;
}

h5part_int64_t
h5ptc_writestepattrib_r8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float64_t *data,
	const h5part_float64_t *nelem
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteStepAttrib (
		filehandle, name, H5PART_FLOAT64, data, *nelem);

	return herr;
}

h5part_int64_t
h5ptc_readstepattrib_r8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float64_t *data
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadStepAttrib (
		filehandle, name, (void*)data);

	return herr;
}

h5part_int64_t
h5ptc_writestepattrib_r4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float32_t *data,
	const h5part_float32_t *nelem
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteStepAttrib (
		filehandle, name, H5PART_FLOAT32, data, *nelem);

	return herr;
}

h5part_int64_t
h5ptc_readstepattrib_r4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float32_t *data
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadStepAttrib (
		filehandle, name, (void*)data);

	return herr;
}

h5part_int64_t
h5ptc_writestepattrib_i8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int64_t *data,
	const h5part_int64_t *nelem
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteStepAttrib (
		filehandle, name, H5PART_INT64, data, *nelem);

	return herr;
}

h5part_int64_t
h5ptc_readstepattrib_i8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int64_t *data
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadStepAttrib (
		filehandle, name, (void*)data);

	return herr;
}

h5part_int64_t
h5ptc_writestepattrib_i4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int32_t *data,
	const h5part_int32_t *nelem
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteStepAttrib (
		filehandle, name, H5PART_INT32, data, *nelem);

	return herr;
}

h5part_int64_t
h5ptc_readstepattrib_i4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int32_t *data
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadStepAttrib (
		filehandle, name, (void*)data);

	return herr;
}

/*=================== Attributes ================*/

h5part_int64_t
h5ptc_writefileattrib_string (
	const h5part_int64_t *f,
	const char *attrib_name,
	const char *attrib_value
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteFileAttribString (
		filehandle, attrib_name, attrib_value );

	return herr;
}

h5part_int64_t
h5ptc_writestepattrib_string (
	const h5part_int64_t *f,
	const char *attrib_name,
	const char *attrib_value
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartWriteStepAttribString (
		filehandle, attrib_name, attrib_value );

	return herr;
}

/* Reading attributes ************************* */

h5part_int64_t
h5ptc_getnstepattribs (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartGetNumStepAttribs ( filehandle );
}

h5part_int64_t
h5ptc_getnfileattribs (
	const h5part_int64_t *f
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	return H5PartGetNumFileAttribs ( filehandle );
}

h5part_int64_t
h5ptc_getstepattribinfo (
	const h5part_int64_t *f,
	const h5part_int64_t *idx,
	char *name,
	h5part_int64_t *nelem,
	const h5part_int64_t *l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;
	h5part_int64_t type;

	h5part_int64_t herr = H5PartGetStepAttribInfo ( 
		filehandle, *idx, name, *l_name, &type, nelem);

	return herr;
}

h5part_int64_t
h5ptc_getfileattribinfo (
	const h5part_int64_t *f,
	const h5part_int64_t *idx,
	char *name,
	h5part_int64_t *nelem,
	const h5part_int64_t *l_name ) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;
	h5part_int64_t type;

	h5part_int64_t herr = H5PartGetFileAttribInfo ( 
		filehandle, *idx, name, *l_name, &type, nelem);

	return herr;
}

h5part_int64_t
h5ptc_readstepattrib_string (
	const h5part_int64_t *f,
	const char *attrib_name,
	char *attrib_value
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;
	
	h5part_int64_t herr = H5PartReadStepAttrib (
		filehandle, attrib_name, attrib_value );

	return herr;
}

h5part_int64_t
h5ptc_readfileattrib_string (
	const h5part_int64_t *f,
	const char *attrib_name,
	char *attrib_value
	) {
		
	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	h5part_int64_t herr = H5PartReadFileAttrib (
		filehandle, attrib_name, attrib_value );

	return herr;
}
