module h5part
 implicit none
 integer*8, parameter :: H5PART_INT64   = 1
 integer*8, parameter :: H5PART_INT32   = 2
 integer*8, parameter :: H5PART_FLOAT64 = 3
 integer*8, parameter :: H5PART_FLOAT32 = 4
 integer*8, parameter :: H5PART_CHAR    = 5
 integer*8, parameter :: H5PART_STRING  = 6
 character(len=7), dimension(6), parameter :: &
    h5pt_type = (/'INT64  ',&
                  'INT32  ',&
                  'FLOAT64',&
                  'FLOAT32',&
                  'CHAR   ',&
                  'STRING '/)

interface

INTEGER*8 FUNCTION h5pt_isvalid_( ihandle ) bind(C)
    INTEGER*8, INTENT(IN) :: ihandle   !< the filename to open for reading
END FUNCTION

! Declaration of subroutines for Fortran Bindings

!> \defgroup h5part_f90_api H5Part F90 API

!> \ingroup h5part_f90_api
!! \defgroup h5partf_open       File Opening and Closing
!<

!> \ingroup h5part_f90_api
!! \defgroup h5partf_model      Setting up the Data Model
!<

!> \ingroup h5part_f90_api
!! \defgroup h5partf_data       Reading and Writing Datasets
!<

!> \ingroup h5part_f90_api
!! \defgroup h5partf_attrib     Reading and Writing Attributes
!<


!!!!!!!! File Opening and Closing !!!!!!!!

!> \ingroup h5partf_open
!! Opens a file for reading. See \ref H5PartOpenFile
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openr( filename ) bind(C)
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for reading
END FUNCTION

!> \ingroup h5partf_open
!! Opens a file for writing in truncate mode. See \ref H5PartOpenFile
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openw ( filename ) bind(C)
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for writing
END FUNCTION

!> \ingroup h5partf_open
!! Opens a file for writing in append mode. See \ref H5PartOpenFile
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_opena ( filename ) bind(C)
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for appending
END FUNCTION

!> \ingroup h5partf_open
!! Opens a parallel file for reading.
!! See \ref H5PartOpenFileParallel
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openr_par ( filename, mpi_communicator ) bind(C)
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for reading
    INTEGER, INTENT(IN) :: mpi_communicator     !< the MPI communicator used by the program
END FUNCTION

!> \ingroup h5partf_open
!! Opens a parallel file for writing in truncate mode.
!! See \ref H5PartOpenFileParallel
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openw_par ( filename, mpi_communicator ) bind(C)
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for writing
    INTEGER, INTENT(IN) :: mpi_communicator     !< the MPI_Communicator used by the program
END FUNCTION

!> \ingroup h5partf_open
!! Opens a parallel file for writing in append mode.
!! See \ref H5PartOpenFileParallel
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_opena_par ( filename, mpi_communicator ) bind(C)
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for appending
    INTEGER, INTENT(IN) :: mpi_communicator     !< the MPI_Communicator used by the program
END FUNCTION

!> \ingroup h5partf_open
!! Opens a file for reading and specifies an HDF5 alignment.
!! See \ref H5PartOpenFileAlign
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openr_align ( filename, align ) bind(C)
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for reading
    INTEGER*8, INTENT(IN) :: align              !< alignment value in bytes
END FUNCTION

!> \ingroup h5partf_open
!! Opens a file for writing in truncate mode and specifies an HDF5 alignment.
!! See \ref H5PartOpenFileAlign
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openw_align ( filename, align ) bind(C)
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for writing
    INTEGER*8, INTENT(IN) :: align              !< alignment value in bytes
END FUNCTION
 
!> \ingroup h5partf_open
!! Opens a file for writing in append mode and specifies an HDF5 alignment.
!! See \ref H5PartOpenFileAlign
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_opena_align ( filename, align ) bind(C)
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for appending
    INTEGER*8, INTENT(IN) :: align              !< alignment value in bytes
END FUNCTION

!> \ingroup h5partf_open
!! Opens a parallel file for reading and specifies an HDF5 alignment.
!! See \ref H5PartOpenFileParallelAlign
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openr_par_align ( filename, mpi_communicator, align ) bind(C)
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for reading
    INTEGER, INTENT(IN) :: mpi_communicator     !< the MPI_Communicator used by the program
    INTEGER*8, INTENT(IN) :: align              !< alignment value in bytes
END FUNCTION

!> \ingroup h5partf_open
!! Opens a parallel file for writing in truncate mode and specifies
!! an HDF5 alignment.
!! See \ref H5PartOpenFileParallelAlign
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openw_par_align ( filename, mpi_communicator, align, flags ) bind(C)
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for writing
    INTEGER, INTENT(IN) :: mpi_communicator     !< the MPI_Communicator used by the program
    INTEGER*8, INTENT(IN) :: align              !< alignment value in bytes
    CHARACTER(LEN=*), INTENT(IN) :: flags       !< additional flags
END FUNCTION

!> \ingroup h5partf_open
!! Opens a parallel file for writing in append mode and specifies
!! an HDF5 alignment.
!! See \ref H5PartOpenFileParallelAlign
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_opena_par_align ( filename, mpi_communicator, align, flags ) bind(C)
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for appending
    INTEGER, INTENT(IN) :: mpi_communicator     !< the MPI_Communicator used by the program
    INTEGER*8, INTENT(IN) :: align              !< alignment value in bytes
    CHARACTER(LEN=*), INTENT(IN) :: flags       !< additional flags
END FUNCTION

!> \ingroup h5partf_open
!! Closes a file. See \ref H5PartCloseFile
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_close ( filehandle ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_open
!! See \ref H5PartSetVerbosityLevel
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_set_verbosity_level ( level ) bind(C)
    INTEGER*8, INTENT(IN) :: level      !< the level from 0 (no output) to 5 (most detailed)
END FUNCTION


!!!!!!!! Setting up the Data Model !!!!!!!!

!> \ingroup h5partf_model
!! See \ref H5PartSetNumParticles
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_setnpoints ( filehandle, npoints ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: npoints    !< the number of particles on *this* processor
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartSetNumParticlesStrided
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_setnpoints_strided ( filehandle, npoints, stride ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: npoints    !< the number of particles on *this* processor
    INTEGER*8, INTENT(IN) :: stride     !< the stride value (e.g. the number of fields in the particle data array)
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartSetStep
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_setstep (filehandle,step) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: step       !< a timestep value >= 1
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartGetNumSteps
!! \return the number of steps or error code
!<
INTEGER*8 FUNCTION h5pt_getnsteps (filehandle) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartGetNumDatasets
!! \return the number of datasets or error code
!<
INTEGER*8 FUNCTION h5pt_getndatasets (filehandle) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartGetNumParticles
!! \return the number of particles or error code
!<
INTEGER*8 FUNCTION h5pt_getnpoints (filehandle) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartGetDatasetName
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_getdatasetname (filehandle,index,name) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: index              !< index of dataset to query (starting from 0)
    CHARACTER(LEN=*), INTENT(OUT) :: name       !< buffer to read the dataset name into 
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartGetDatasetInfo
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_getdatasetinfo (filehandle,index,name,data_type,nelem) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: index              !< index of dataset to query (starting from 0)
    CHARACTER(LEN=*), INTENT(OUT) :: name       !< buffer to read the dataset name into 
    INTEGER*8, INTENT(OUT) :: data_type         !< type of data in dataset
    INTEGER*8, INTENT(OUT) :: nelem             !< number of elements
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartSetView
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_setview (filehandle,start,end) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: start      !< offset of the first particle in the view
    INTEGER*8, INTENT(IN) :: end        !< offset of the last particle in the view (inclusive)
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartSetViewIndices
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_setview_indices (filehandle,indices,nelem) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: indices(*) !< list of indicies to select in this view
    INTEGER*8, INTENT(IN) :: nelem      !< number of particles in the list
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartResetView
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_resetview (filehandle) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartResetView
!! \return 1 if true, 0 if false, or error code
!<
INTEGER*8 FUNCTION h5pt_hasview (filehandle) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartGetView
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_getview (filehandle,start,end) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
    INTEGER*8, INTENT(OUT) :: start     !< buffer to store the offset of the first particle in the view
    INTEGER*8, INTENT(OUT) :: end       !< buffer to store the offset of the last particle in the view (inclusive)
END FUNCTION


!!!!!!!! Reading and Writing Datasets !!!!!!!!

!> \ingroup h5partf_data
!! See \ref H5PartWriteDataFloat64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_writedata_r8 ( filehandle, name, data ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    REAL*8, INTENT(IN) :: data(*)               !< the array of float64 data to write
END FUNCTION

!> \ingroup h5partf_data
!! See \ref H5PartWriteDataFloat32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_writedata_r4 ( filehandle, name, data ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    REAL, INTENT(IN) :: data(*)                 !< the array of float32 data to write
END FUNCTION

!> \ingroup h5partf_data
!! See \ref H5PartWriteDataInt64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_writedata_i8 ( filehandle, name, data ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    INTEGER*8, INTENT(IN) :: data(*)            !< the array of int64 data to write
END FUNCTION

!> \ingroup h5partf_data
!! See \ref H5PartWriteDataInt32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_writedata_i4 ( filehandle, name, data ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    INTEGER, INTENT(IN) :: data(*)              !< the array of int32 data to write
END FUNCTION


!> \ingroup h5partf_data
!! See \ref H5PartReadDataFloat64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_readdata_r8 (filehandle,name,data) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    REAL*8, INTENT(OUT) :: data(*)              !< array to read float64 data into
END FUNCTION

!> \ingroup h5partf_data
!! See \ref H5PartReadDataFloat32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_readdata_r4 (filehandle,name,data) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    REAL, INTENT(OUT) :: data(*)                !< array to read float32 data into
END FUNCTION

!> \ingroup h5partf_data
!! See \ref H5PartReadDataInt64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_readdata_i8 (filehandle,name,data) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    INTEGER*8, INTENT(OUT) :: data(*)           !< array to read int64 data into
END FUNCTION

!> \ingroup h5partf_data
!! See \ref H5PartReadDataInt32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_readdata_i4 (filehandle,name,data) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    INTEGER, INTENT(OUT) :: data(*)             !< array to read int32 data into
END FUNCTION


!!!!!!!! Reading and Writing Attributes !!!!!!!!

!> \ingroup h5partf_attrib
!! See \ref H5PartWriteFileAttribString
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_writefilestring (filehandle,name,value) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the attribute
    CHARACTER(LEN=*), INTENT(IN) :: value       !< the string value to store
END FUNCTION

!> \ingroup h5partf_attrib
!! See \ref H5PartWriteStepAttribString
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_writestepstring (filehandle,name,value) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the attribute
    CHARACTER(LEN=*), INTENT(IN) :: value       !< the string value to store
END FUNCTION

!> \ingroup h5partf_attrib
!! Reads the attribute \c name in the file root ("/")
!! into the string buffer \c value.
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_readfilestring (filehandle,name,value) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the attribute
    CHARACTER(LEN=*), INTENT(OUT) :: value      !< buffer to read the string value into
END FUNCTION

!> \ingroup h5partf_attrib
!! Reads the attribute \c name in the current step
!! into the string buffer \c value.
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_readstepstring (filehandle,name,value) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the attribute
    CHARACTER(LEN=*), INTENT(OUT) :: value      !< buffer to read the string value into
END FUNCTION

!> \ingroup h5partf_attrib
!! See \ref H5PartGetNumStepAttribs
!! \return number of attributes or error code
!<
INTEGER*8 FUNCTION h5pt_getnstepattribs (filehandle) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_attrib
!! See \ref H5PartGetNumFileAttribs
!! \return number of attributes or error code
!<
INTEGER*8 FUNCTION h5pt_getnfileattribs (filehandle) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_attrib
!! See \ref H5PartGetStepAttribInfo
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_getstepattribinfo (filehandle,idx,name,nelem) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: idx              !< index of the attribute to query (starting from 0)
    CHARACTER(LEN=*), INTENT(OUT) :: name       !< buffer to read the attribute name into
    INTEGER*8, INTENT(OUT) :: nelem             !< number of elements in the attribute's array
END FUNCTION

!> \ingroup h5partf_attrib
!! See \ref H5PartGetFileAttribInfo
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_getfileattribinfo (filehandle,idx,name,nelem) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: idx              !< index of the attribute to query (starting from 0)
    CHARACTER(LEN=*), INTENT(OUT) :: name       !< buffer to read the attribute name into
    INTEGER*8, INTENT(OUT) :: nelem             !< number of elements in the attribute's array
END FUNCTION

end interface

!
!
!**** attributes interface starts here
!
!

interface
!
!< \ingroup h5partf_attrib
!! See \ref H5PartWritefileAttribFloat64
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writefileattrib_r8 ( filehandle, name, data, nelem) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    REAL*8, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readfileattrib_r8 ( filehandle, name, data ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    REAL*8, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritefileAttribFloat32
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writefileattrib_r4 ( filehandle, name, data, nelem) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    REAL*4, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readfileattrib_r4 ( filehandle, name, data ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    REAL*4, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritefileAttribInt64
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writefileattrib_i8 ( filehandle, name, data, nelem) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    INTEGER*8, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readfileattrib_i8 ( filehandle, name, data ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    INTEGER*8, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritefileAttribInt32
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writefileattrib_i4 ( filehandle, name, data, nelem) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    INTEGER*4, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readfileattrib_i4 ( filehandle, name, data ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    INTEGER*4, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritestepAttribFloat64
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writestepattrib_r8 ( filehandle, name, data, nelem) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    REAL*8, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readstepattrib_r8 ( filehandle, name, data ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    REAL*8, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritestepAttribFloat32
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writestepattrib_r4 ( filehandle, name, data, nelem) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    REAL*4, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readstepattrib_r4 ( filehandle, name, data ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    REAL*4, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritestepAttribInt64
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writestepattrib_i8 ( filehandle, name, data, nelem) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    INTEGER*8, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readstepattrib_i8 ( filehandle, name, data ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    INTEGER*8, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritestepAttribInt32
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writestepattrib_i4 ( filehandle, name, data, nelem) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    INTEGER*4, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readstepattrib_i4 ( filehandle, name, data ) bind(C)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    INTEGER*4, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

end interface

end module h5part
