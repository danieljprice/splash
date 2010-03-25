#include "H5Part.h"

h5part_int64_t
h5pt_writefileattrib_r8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float64_t *data,
	const h5part_float64_t *nelem,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWriteFileAttrib (
		filehandle, name2, H5PART_FLOAT64, data, *nelem);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_readfileattrib_r8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float64_t *data,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartReadFileAttrib (
		filehandle, name2, (void*)data);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_writefileattrib_r4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float32_t *data,
	const h5part_float32_t *nelem,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWriteFileAttrib (
		filehandle, name2, H5PART_FLOAT32, data, *nelem);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_readfileattrib_r4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float32_t *data,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartReadFileAttrib (
		filehandle, name2, (void*)data);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_writefileattrib_i8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int64_t *data,
	const h5part_int64_t *nelem,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWriteFileAttrib (
		filehandle, name2, H5PART_INT64, data, *nelem);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_readfileattrib_i8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int64_t *data,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartReadFileAttrib (
		filehandle, name2, (void*)data);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_writefileattrib_i4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int32_t *data,
	const h5part_int32_t *nelem,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWriteFileAttrib (
		filehandle, name2, H5PART_INT32, data, *nelem);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_readfileattrib_i4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int32_t *data,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartReadFileAttrib (
		filehandle, name2, (void*)data);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_writestepattrib_r8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float64_t *data,
	const h5part_float64_t *nelem,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWriteStepAttrib (
		filehandle, name2, H5PART_FLOAT64, data, *nelem);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_readstepattrib_r8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float64_t *data,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartReadStepAttrib (
		filehandle, name2, (void*)data);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_writestepattrib_r4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float32_t *data,
	const h5part_float32_t *nelem,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWriteStepAttrib (
		filehandle, name2, H5PART_FLOAT32, data, *nelem);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_readstepattrib_r4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_float32_t *data,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartReadStepAttrib (
		filehandle, name2, (void*)data);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_writestepattrib_i8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int64_t *data,
	const h5part_int64_t *nelem,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWriteStepAttrib (
		filehandle, name2, H5PART_INT64, data, *nelem);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_readstepattrib_i8 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int64_t *data,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartReadStepAttrib (
		filehandle, name2, (void*)data);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_writestepattrib_i4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int32_t *data,
	const h5part_int32_t *nelem,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartWriteStepAttrib (
		filehandle, name2, H5PART_INT32, data, *nelem);

	free ( name2 );
	return herr;
}

h5part_int64_t
h5pt_readstepattrib_i4 (
	h5part_int64_t *f,
	const char *name,
	const h5part_int32_t *data,
	const int l_name
	) {

	H5PartFile *filehandle = (H5PartFile*)(size_t)*f;

	char *name2 =_H5Part_strdupfor2c ( name, l_name );

	h5part_int64_t herr = H5PartReadStepAttrib (
		filehandle, name2, (void*)data);

	free ( name2 );
	return herr;
}
