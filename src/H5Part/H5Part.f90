!-------------------------------------------------------------
!
! Fortran 90 interface to H5Part, using
! the Fortran 2003 C interoperability module (iso_c_binding)
!
! Written by Daniel Price 08/04/10
!  daniel.price@sci.monash.edu.au
!
! We first specify the interfaces to the C interface routines
! used to handle the H5partfile container object. However, all
! string conversion is done in the Fortran, not in the C.
!
!-------------------------------------------------------------
module h5part
 use, intrinsic :: iso_c_binding, only:c_char,c_int,c_int64_t,c_double,c_float
 implicit none
 integer(kind=c_int64_t), parameter, public :: H5PART_INT64   = 1
 integer(kind=c_int64_t), parameter, public :: H5PART_INT32   = 2
 integer(kind=c_int64_t), parameter, public :: H5PART_FLOAT64 = 3
 integer(kind=c_int64_t), parameter, public :: H5PART_FLOAT32 = 4
 integer(kind=c_int64_t), parameter, public :: H5PART_CHAR    = 5
 integer(kind=c_int64_t), parameter, public :: H5PART_STRING  = 6
 character(len=7), dimension(6), parameter, public :: &
    h5part_type = (/'INT64  ',&
                    'INT32  ',&
                    'FLOAT64',&
                    'FLOAT32',&
                    'CHAR   ',&
                    'STRING '/)
!
! interfaces provided by this module
!
 public :: h5pt_openr,h5pt_openw,h5pt_opena,h5pt_close
 public :: h5pt_openr_align,h5pt_openw_align,h5pt_opena_align

#ifdef PARALLEL_IO
 public :: h5pt_openr_par,h5pt_openw_par,h5pt_opena_par
 public :: h5pt_openr_par_align,h5pt_openw_par_align,h5pt_opena_par_align
#endif
 
 public :: h5pt_setnpoints,h5pt_setnpoints_strided
 public :: h5pt_getnpoints
 public :: h5pt_setstep,h5pt_getnsteps,h5pt_getndatasets
 public :: h5pt_getdatasetname,h5pt_getdatasetinfo
 
 public :: h5pt_setview,h5pt_setview_indices
 public :: h5pt_getview
 public :: h5pt_resetview,h5pt_hasview
 
 public :: h5pt_set_verbosity_level
 
 public :: h5pt_writedata
 public :: h5pt_readdata
!
! the type-specific routines are also public
! (could make these private to allow only 
!  the generic interface to be used)
!
 public :: h5pt_writedata_r8,h5pt_writedata_r4, &
           h5pt_writedata_i8,h5pt_writedata_i4
 public :: h5pt_readdata_r8,h5pt_readdata_r4, &
           h5pt_readdata_i8,h5pt_readdata_i4
 
 private

!
! generic interface for writing data of any type
! 
 interface h5pt_writedata
    module procedure h5pt_writedata_i4,h5pt_writedata_i8, &
                     h5pt_writedata_r4,h5pt_writedata_r8
 end interface h5pt_writedata

!
! generic interface for reading data of any type
! 
 interface h5pt_readdata
    module procedure h5pt_readdata_i4,h5pt_readdata_i8, &
                     h5pt_readdata_r4,h5pt_readdata_r8
 end interface h5pt_readdata

!---------------------------
!
! interfaces to c routines
!
!---------------------------
interface
!
! Opening and closing files
!
integer(kind=c_int64_t) function h5ptc_openr( filename ) bind(C)
    import
    character(kind=c_char), dimension(1), intent(in) :: filename    !< the filename to open for reading
end function

integer(kind=c_int64_t) function h5ptc_openw ( filename ) bind(C)
    import
    character(kind=c_char), dimension(1), intent(in) :: filename    !< the filename to open for writing
end function

integer(kind=c_int64_t) function h5ptc_opena ( filename ) bind(C)
    import
    character(kind=c_char), dimension(1), intent(in) :: filename    !< the filename to open for appending
end function

integer(kind=c_int64_t) function h5ptc_openr_align ( filename, align ) bind(C)
    import
    character(kind=c_char), dimension(1), intent(in) :: filename    !< the filename to open for reading
    integer(kind=c_int64_t), intent(in) :: align              !< alignment value in bytes
end function

integer(kind=c_int64_t) function h5ptc_openw_align ( filename, align ) bind(C)
    import
    character(kind=c_char), dimension(1), intent(in) :: filename    !< the filename to open for writing
    integer(kind=c_int64_t), intent(in) :: align              !< alignment value in bytes
end function
 
integer(kind=c_int64_t) function h5ptc_opena_align ( filename, align ) bind(C)
    import
    character(kind=c_char), dimension(1), intent(in) :: filename    !< the filename to open for appending
    integer(kind=c_int64_t), intent(in) :: align              !< alignment value in bytes
end function

integer(kind=c_int64_t) function h5pt_close ( filehandle ) bind(C,name="h5ptc_close")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned during file open
end function

#ifdef PARALLEL_IO
!
! Opening files (parallel I/O)
!
integer(kind=c_int64_t) function h5ptc_openr_par ( filename, mpi_communicator ) bind(C)
    import
    character(kind=c_char), dimension(1), intent(in) :: filename    !< the filename to open for reading
    integer, intent(in) :: mpi_communicator     !< the MPI communicator used by the program
end function

integer(kind=c_int64_t) function h5ptc_openw_par ( filename, mpi_communicator ) bind(C)
    import
    character(kind=c_char), dimension(1), intent(in) :: filename    !< the filename to open for writing
    integer, intent(in) :: mpi_communicator     !< the MPI_Communicator used by the program
end function

integer(kind=c_int64_t) function h5ptc_opena_par ( filename, mpi_communicator ) bind(C)
    import
    character(kind=c_char), dimension(1), intent(in) :: filename    !< the filename to open for appending
    integer, intent(in) :: mpi_communicator     !< the MPI_Communicator used by the program
end function

integer(kind=c_int64_t) function h5ptc_openr_par_align ( filename, mpi_communicator, align ) bind(C)
    import
    character(kind=c_char), dimension(1), intent(in) :: filename    !< the filename to open for reading
    integer, intent(in) :: mpi_communicator     !< the MPI_Communicator used by the program
    integer(kind=c_int64_t), intent(in) :: align              !< alignment value in bytes
end function

integer(kind=c_int64_t) function h5ptc_openw_par_align ( filename, mpi_communicator, align, flags ) bind(C)
    import
    character(kind=c_char), dimension(1), intent(in) :: filename    !< the filename to open for writing
    integer, intent(in) :: mpi_communicator     !< the MPI_Communicator used by the program
    integer(kind=c_int64_t), intent(in) :: align              !< alignment value in bytes
    character(kind=c_char), dimension(1), intent(in) :: flags       !< additional flags
end function

integer(kind=c_int64_t) function h5ptc_opena_par_align ( filename, mpi_communicator, align, flags ) bind(C)
    import
    character(kind=c_char), dimension(1), intent(in) :: filename    !< the filename to open for appending
    integer, intent(in) :: mpi_communicator     !< the MPI_Communicator used by the program
    integer(kind=c_int64_t), intent(in) :: align              !< alignment value in bytes
    character(kind=c_char), dimension(1), intent(in) :: flags       !< additional flags
end function
#endif

!
! Setting up the data model
! (where no strings are passed these are just interfaces to the c routines that convert the filehandle)
!
integer(kind=c_int64_t) function h5pt_setnpoints ( filehandle, npoints ) bind(C,name="h5ptc_setnpoints")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: npoints    !< the number of particles on *this* processor
end function

integer(kind=c_int64_t) function h5pt_setnpoints_strided ( filehandle, npoints, stride ) &
                                 bind(C,name="h5ptc_setnpoints_strided")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: npoints    !< the number of particles on *this* processor
    integer(kind=c_int64_t), intent(in) :: stride     !< the stride value (e.g. the number of fields in the particle data array)
end function

integer(kind=c_int64_t) function h5pt_setstep (filehandle,step) bind(C,name="h5ptc_setstep")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: step       !< a timestep value >= 1
end function

integer(kind=c_int64_t) function h5pt_getnsteps (filehandle) bind(C,name="h5ptc_getnsteps")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
end function

integer(kind=c_int64_t) function h5pt_getndatasets (filehandle) bind(C,name="h5ptc_getndatasets")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
end function

integer(kind=c_int64_t) function h5pt_getnpoints (filehandle) bind(C,name="h5ptc_getnpoints")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
end function

integer(kind=c_int64_t) function h5ptc_getdatasetname (filehandle,index,name,l_name) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: index              !< index of dataset to query (starting from 0)
    character(kind=c_char), dimension(1), intent(out) :: name       !< buffer to read the dataset name into 
    integer(kind=c_int64_t), intent(in) :: l_name             !< size of name
end function

integer(kind=c_int64_t) function h5ptc_getdatasetinfo (filehandle,index,name,data_type,nelem,l_name) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: index              !< index of dataset to query (starting from 0)
    character(kind=c_char), dimension(1), intent(out) :: name       !< buffer to read the dataset name into 
    integer(kind=c_int64_t), intent(out) :: data_type         !< type of data in dataset
    integer(kind=c_int64_t), intent(out) :: nelem             !< number of elements
    integer(kind=c_int64_t), intent(in) :: l_name             !< size of name
end function

!
! Setting and getting views
! (as no strings are passed these are just interfaces to the c routines that convert the filehandle)
!
integer(kind=c_int64_t) function h5pt_setview (filehandle,start,end) bind(C,name="h5ptc_setview")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: start      !< offset of the first particle in the view
    integer(kind=c_int64_t), intent(in) :: end        !< offset of the last particle in the view (inclusive)
end function

integer(kind=c_int64_t) function h5pt_setview_indices (filehandle,indices,nelem) bind(C,name="h5ptc_setview_indices")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: indices(*) !< list of indicies to select in this view
    integer(kind=c_int64_t), intent(in) :: nelem      !< number of particles in the list
end function

integer(kind=c_int64_t) function h5pt_resetview (filehandle) bind(C,name="h5ptc_resetview")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned during file open
end function

integer(kind=c_int64_t) function h5pt_hasview (filehandle) bind(C,name="h5ptc_hasview")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned during file open
end function

integer(kind=c_int64_t) function h5pt_getview (filehandle,start,end) bind(C,name="h5ptc_getview")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned during file open
    integer(kind=c_int64_t), intent(out) :: start     !< buffer to store the offset of the first particle in the view
    integer(kind=c_int64_t), intent(out) :: end       !< buffer to store the offset of the last particle in the view (inclusive)
end function

!
! Reading and writing datasets
!
integer(kind=c_int64_t) function h5ptc_writedata_r8 ( filehandle, name, data ) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(kind=c_char), dimension(1), intent(in) :: name  !< the name of the dataset
    real(kind=c_double), intent(in) :: data(*)                !< the array of float64 data to write
end function

integer(kind=c_int64_t) function h5ptc_readdata_r8 (filehandle,name,data) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(kind=c_char), dimension(1), intent(in) :: name  !< the name of the dataset
    real(kind=c_double), intent(out) :: data(*)               !< array to read float64 data into
end function

integer(kind=c_int64_t) function h5ptc_writedata_r4 ( filehandle, name, data ) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(kind=c_char), dimension(1), intent(in) :: name  !< the name of the dataset
    real(kind=c_float), intent(in) :: data(*)                 !< the array of float32 data to write
end function

integer(kind=c_int64_t) function h5ptc_readdata_r4 (filehandle,name,data) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(kind=c_char), dimension(1), intent(in) :: name  !< the name of the dataset
    real(kind=c_float), intent(out) :: data(*)                !< array to read float32 data into
end function

integer(kind=c_int64_t) function h5ptc_writedata_i8 ( filehandle, name, data ) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(kind=c_char), dimension(1), intent(in) :: name  !< the name of the dataset
    integer(kind=c_int64_t), intent(in) :: data(*)            !< the array of int64 data to write
end function

integer(kind=c_int64_t) function h5ptc_readdata_i8 (filehandle,name,data) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(kind=c_char), dimension(1), intent(in) :: name        !< the name of the dataset
    integer(kind=c_int64_t), intent(out) :: data(*)           !< array to read int64 data into
end function

integer(kind=c_int64_t) function h5ptc_writedata_i4 ( filehandle, name, data ) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(kind=c_char), dimension(1), intent(in) :: name        !< the name of the dataset
    integer(kind=c_int), intent(in) :: data(*)              !< the array of int32 data to write
end function

integer(kind=c_int64_t) function h5ptc_readdata_i4 (filehandle,name,data) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(kind=c_char), dimension(1), intent(in) :: name        !< the name of the dataset
    integer(kind=c_int), intent(out) :: data(*)             !< array to read int32 data into
end function

!
! routines where no conversion of anything is needed: i.e. just call H5Part C routines directly
!
integer(kind=c_int64_t) function h5pt_set_verbosity_level ( level ) bind(C,name="H5PartSetVerbosityLevel")
    import
    integer(kind=c_int64_t), intent(in) :: level      !< the level from 0 (no output) to 5 (most detailed)
end function

end interface

contains
!---------------------------------------------------------------------------
!
! wrappers for functions with string arguments: 
! converts strings into C strings
!
!---------------------------------------------------------------------------
integer(kind=c_int64_t) function h5pt_openr( filename )
    implicit none
    character(len=*), intent(in) :: filename    !< the filename to open for reading
    
    h5pt_openr = h5ptc_openr( cstring(filename) )

end function

integer(kind=c_int64_t) function h5pt_openw ( filename )
    implicit none
    character(len=*), intent(in) :: filename    !< the filename to open for writing
    
    h5pt_openw = h5ptc_openw ( cstring(filename) )

end function

integer(kind=c_int64_t) function h5pt_opena ( filename )
    implicit none
    character(len=*), intent(in) :: filename    !< the filename to open for appending

    h5pt_opena = h5ptc_opena ( cstring(filename) )

end function

integer(kind=c_int64_t) function h5pt_openr_align ( filename, align )
    implicit none
    character(len=*), intent(in) :: filename    !< the filename to open for reading
    integer(kind=c_int64_t), intent(in) :: align              !< alignment value in bytes
    
    h5pt_openr_align = h5ptc_openr_align ( cstring(filename), align )
    
end function

integer(kind=c_int64_t) function h5pt_openw_align ( filename, align )
    implicit none
    character(len=*), intent(in) :: filename    !< the filename to open for writing
    integer(kind=c_int64_t), intent(in) :: align              !< alignment value in bytes

    h5pt_openw_align = h5ptc_openw_align ( cstring(filename), align )

end function
 
integer(kind=c_int64_t) function h5pt_opena_align ( filename, align )
    implicit none
    character(len=*), intent(in) :: filename    !< the filename to open for appending
    integer(kind=c_int64_t), intent(in) :: align              !< alignment value in bytes
    
    h5pt_opena_align = h5ptc_opena_align ( cstring(filename), align )
    
end function

#ifdef PARALLEL_IO
!
! opening files (parallel I/O)
!
integer(kind=c_int64_t) function h5pt_openr_par ( filename, mpi_communicator )
    implicit none
    character(len=*), intent(in) :: filename    !< the filename to open for reading
    integer, intent(in) :: mpi_communicator     !< the MPI communicator used by the program

    h5pt_openr_par = h5ptc_openr_par ( cstring(filename), mpi_communicator )

end function

integer(kind=c_int64_t) function h5pt_openw_par ( filename, mpi_communicator )
    implicit none
    character(len=*), intent(in) :: filename    !< the filename to open for writing
    integer, intent(in) :: mpi_communicator     !< the MPI_Communicator used by the program

    h5pt_openw_par = h5ptc_openw_par ( cstring(filename), mpi_communicator )

end function

integer(kind=c_int64_t) function h5pt_opena_par ( filename, mpi_communicator )
    implicit none
    character(len=*), intent(in) :: filename    !< the filename to open for appending
    integer, intent(in) :: mpi_communicator     !< the MPI_Communicator used by the program

    h5pt_opena_par = h5ptc_opena_par ( cstring(filename), mpi_communicator )

end function

integer(kind=c_int64_t) function h5pt_openr_par_align ( filename, mpi_communicator, align )
    implicit none
    character(len=*), intent(in) :: filename    !< the filename to open for reading
    integer, intent(in) :: mpi_communicator     !< the MPI_Communicator used by the program
    integer(kind=c_int64_t), intent(in) :: align              !< alignment value in bytes
    
    h5pt_openr_par_align = h5ptc_openr_par_align ( cstring(filename), mpi_communicator, align )

end function

integer(kind=c_int64_t) function h5pt_openw_par_align ( filename, mpi_communicator, align, flags )
    implicit none
    character(len=*), intent(in) :: filename    !< the filename to open for writing
    integer, intent(in) :: mpi_communicator     !< the MPI_Communicator used by the program
    integer(kind=c_int64_t), intent(in) :: align              !< alignment value in bytes
    character(len=*), intent(in) :: flags       !< additional flags
    
    h5pt_openw_par_align = h5ptc_openw_par_align ( cstring(filename), mpi_communicator, align, cstring(flags) )

end function

integer(kind=c_int64_t) function h5pt_opena_par_align ( filename, mpi_communicator, align, flags )
    implicit none
    character(len=*), intent(in) :: filename    !< the filename to open for appending
    integer, intent(in) :: mpi_communicator     !< the MPI_Communicator used by the program
    integer(kind=c_int64_t), intent(in) :: align              !< alignment value in bytes
    character(len=*), intent(in) :: flags       !< additional flags
    
    h5pt_opena_par_align = h5ptc_opena_par_align ( cstring(filename), mpi_communicator, align, cstring(flags) )
    
end function
#endif

!
! Setting up the data model
!
integer(kind=c_int64_t) function h5pt_getdatasetname (filehandle,index,name)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: index              !< index of dataset to query (starting from 0)
    character(len=*), intent(out) :: name       !< buffer to read the dataset name into 
    integer(kind=c_int64_t) :: l_name
    
    l_name = len(name)
    h5pt_getdatasetname = h5ptc_getdatasetname(filehandle,index,name,l_name)
    name = fstring(name)

end function h5pt_getdatasetname

 integer(kind=c_int64_t) function h5pt_getdatasetinfo(filehandle,idx,name,data_type,nelem)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: idx              !< index of dataset to query (starting from 0)
    character(len=*), intent(out) :: name       !< buffer to read the dataset name into 
    integer(kind=c_int64_t), intent(out) :: data_type         !< type of data in dataset
    integer(kind=c_int64_t), intent(out) :: nelem             !< number of elements
    integer(kind=c_int64_t) :: l_name
    
    l_name = len(name)
    h5pt_getdatasetinfo = h5ptc_getdatasetinfo(filehandle,idx,name,data_type,nelem,l_name)
    name = fstring(name)
    
 end function h5pt_getdatasetinfo

!
! Reading and writing datasets
!
integer(kind=c_int64_t) function h5pt_writedata_r8 ( filehandle, name, data )
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(len=*), intent(in) :: name  !< the name of the dataset
    real(kind=c_double), intent(in) :: data(*)                !< the array of float64 data to write

    h5pt_writedata_r8 = h5ptc_writedata_r8 ( filehandle, cstring(name), data )

end function

integer(kind=c_int64_t) function h5pt_readdata_r8 ( filehandle, name, data)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(len=*), intent(in) :: name  !< the name of the dataset
    real(kind=c_double), intent(out) :: data(*)               !< array to read float64 data into

    h5pt_readdata_r8 = h5ptc_readdata_r8 ( filehandle, cstring(name), data)

end function

integer(kind=c_int64_t) function h5pt_writedata_r4 ( filehandle, name, data )
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(len=*), intent(in) :: name  !< the name of the dataset
    real(kind=c_float), intent(in) :: data(*)                 !< the array of float32 data to write

    h5pt_writedata_r4 = h5ptc_writedata_r4 ( filehandle, cstring(name), data )

end function

integer(kind=c_int64_t) function h5pt_readdata_r4 ( filehandle, name, data)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(len=*), intent(in) :: name  !< the name of the dataset
    real(kind=c_float), intent(out) :: data(*)                !< array to read float32 data into

    h5pt_readdata_r4 = h5ptc_readdata_r4 ( filehandle, cstring(name), data)

end function

integer(kind=c_int64_t) function h5pt_writedata_i8 ( filehandle, name, data )
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(len=*), intent(in) :: name  !< the name of the dataset
    integer(kind=c_int64_t), intent(in) :: data(*)            !< the array of int64 data to write

    h5pt_writedata_i8 = h5ptc_writedata_i8 ( filehandle, cstring(name), data )

end function

integer(kind=c_int64_t) function h5pt_readdata_i8 ( filehandle, name, data)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(len=*), intent(in) :: name        !< the name of the dataset
    integer(kind=c_int64_t), intent(out) :: data(*)           !< array to read int64 data into

    h5pt_readdata_i8 = h5ptc_readdata_i8 ( filehandle, cstring(name), data)

end function

integer(kind=c_int64_t) function h5pt_writedata_i4 ( filehandle, name, data )
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(len=*), intent(in) :: name        !< the name of the dataset
    integer(kind=c_int), intent(in) :: data(*)              !< the array of int32 data to write

    h5pt_writedata_i4 = h5ptc_writedata_i4 ( filehandle, cstring(name), data )

end function

integer(kind=c_int64_t) function h5pt_readdata_i4 (filehandle,name,data)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(len=*), intent(in) :: name        !< the name of the dataset
    integer(kind=c_int), intent(out) :: data(*)             !< array to read int32 data into

    h5pt_readdata_i4 = h5ptc_readdata_i4 ( filehandle, cstring(name), data)

end function

!---------------------------------------------------------------------------
!
! function to safely convert a string to c format (ie. with a terminating
! ascii null character)
!
!---------------------------------------------------------------------------
 function cstring(string)
  implicit none
  character(len=*), intent(in) :: string
  character(len=len(string)+1) :: cstring

  cstring = trim(string)//char(0)

 end function cstring

!---------------------------------------------------------------------------
!
! function to safely convert a string from c format (ie. with a terminating
! ascii null character) back to a normal Fortran string
!
!---------------------------------------------------------------------------
 function fstring(string)
  implicit none
  character(len=*), intent(in) :: string  !< the name of the dataset
  character(len=len(string)) :: fstring
  integer :: idx

  idx = index(string,char(0))
  if (idx.gt.1) then
     fstring = string(1:idx-1)
  else
     fstring = ''
  endif

 end function fstring

end module h5part
