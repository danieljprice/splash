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
! This module contains the attributes interface
!
!-------------------------------------------------------------
module h5partattrib
 use, intrinsic :: iso_c_binding, only:c_char,c_int,c_int64_t,c_double,c_float
 implicit none
!
! interfaces provided by this module
!
 public :: h5pt_getnfileattribs,h5pt_getfileattribinfo
 public :: h5pt_getnstepattribs,h5pt_getstepattribinfo
 public :: h5pt_writefileattrib,h5pt_readfileattrib
 public :: h5pt_writestepattrib,h5pt_readstepattrib
!
! the type-specific routines are also public
! (could make these private to allow only 
!  the generic interface to be used)
!
 public :: h5pt_writefileattrib_i4,h5pt_writefileattrib_i8, &
           h5pt_writefileattrib_r4,h5pt_writefileattrib_r8, &
           h5pt_writefileattrib_string
 public :: h5pt_readfileattrib_i4,h5pt_readfileattrib_i8, &
           h5pt_readfileattrib_r4,h5pt_readfileattrib_r8, &
           h5pt_readfileattrib_string
 public :: h5pt_writestepattrib_i4,h5pt_writestepattrib_i8, &
           h5pt_writestepattrib_r4,h5pt_writestepattrib_r8, &
           h5pt_writestepattrib_string
 public :: h5pt_readstepattrib_i4,h5pt_readstepattrib_i8, &
           h5pt_readstepattrib_r4,h5pt_readstepattrib_r8, &
           h5pt_readstepattrib_string

 private

!
! generic interface for writing file attributes of any type
!
 interface h5pt_writefileattrib
   module procedure h5pt_writefileattrib_i4,h5pt_writefileattrib_i8, &
                    h5pt_writefileattrib_r4,h5pt_writefileattrib_r8, &
                    h5pt_writefileattrib_string

 end interface h5pt_writefileattrib
!
! generic interface for reading file attributes of any type
!
 interface h5pt_readfileattrib
   module procedure h5pt_readfileattrib_i4,h5pt_readfileattrib_i8, &
                    h5pt_readfileattrib_r4,h5pt_readfileattrib_r8, &
                    h5pt_readfileattrib_string

 end interface h5pt_readfileattrib
!
! generic interface for writing step attributes of any type
!
 interface h5pt_writestepattrib
   module procedure h5pt_writestepattrib_i4,h5pt_writestepattrib_i8, &
                    h5pt_writestepattrib_r4,h5pt_writestepattrib_r8, &
                    h5pt_writestepattrib_string

 end interface h5pt_writestepattrib
!
! generic interface for reading step attributes of any type
!
 interface h5pt_readstepattrib
   module procedure h5pt_readstepattrib_i4,h5pt_readstepattrib_i8, &
                    h5pt_readstepattrib_r4,h5pt_readstepattrib_r8, &
                    h5pt_readstepattrib_string

 end interface h5pt_readstepattrib

!---------------------------
!
! interfaces to c routines
!
!---------------------------
interface
!
!  file attributes: interfaces to c routines
!
integer(kind=c_int64_t) function h5pt_getnfileattribs (filehandle) bind(C,name="h5ptc_getnfileattribs")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
end function

integer(kind=c_int64_t) function h5ptc_getfileattribinfo (filehandle,idx,name,nelem,l_name) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: idx              !< index of the attribute to query (starting from 0)
    character(kind=c_char), dimension(1), intent(out) :: name       !< buffer to read the attribute name into
    integer(kind=c_int64_t), intent(out) :: nelem             !< number of elements in the attribute's array
    integer(kind=c_int64_t), intent(in) :: l_name             !< size of name
end function

integer(kind=c_int64_t) function h5ptc_writefileattrib_r8 ( filehandle, name, data, nelem) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name !< the name of the attribute
    real(kind=c_double), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
end function

integer(kind=c_int64_t) function h5ptc_readfileattrib_r8 ( filehandle, name, data ) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name   !< the name of the attribute
    real(kind=c_double), intent(out) :: data(*) !< buffer to read value into
end function

integer(kind=c_int64_t) function h5ptc_writefileattrib_r4 ( filehandle, name, data, nelem) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name !< the name of the attribute
    real(kind=c_float), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
end function

integer(kind=c_int64_t) function h5ptc_readfileattrib_r4 ( filehandle, name, data ) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name   !< the name of the attribute
    real(kind=c_float), intent(out) :: data(*) !< buffer to read value into
end function

integer(kind=c_int64_t) function h5ptc_writefileattrib_i8 ( filehandle, name, data, nelem) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name !< the name of the attribute
    integer(kind=c_int64_t), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
end function

integer(kind=c_int64_t) function h5ptc_readfileattrib_i8 ( filehandle, name, data ) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name   !< the name of the attribute
    integer(kind=c_int64_t), intent(out) :: data(*) !< buffer to read value into
end function

integer(kind=c_int64_t) function h5ptc_writefileattrib_i4 ( filehandle, name, data, nelem) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name !< the name of the attribute
    integer(kind=c_int), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
end function

integer(kind=c_int64_t) function h5ptc_readfileattrib_i4 ( filehandle, name, data ) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name   !< the name of the attribute
    integer(kind=c_int), intent(out) :: data(*) !< buffer to read value into
end function

integer(kind=c_int64_t) function h5ptc_writefileattrib_string (filehandle,name,value) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(kind=c_char), dimension(1), intent(in) :: name        !< the name of the attribute
    character(kind=c_char), dimension(1), intent(in) :: value       !< the string value to store
end function

integer(kind=c_int64_t) function h5ptc_readfileattrib_string (filehandle,name,value) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(kind=c_char), dimension(1), intent(in) :: name        !< the name of the attribute
    character(kind=c_char), dimension(1), intent(out) :: value      !< buffer to read the string value into
end function

!
!  step attributes: interfaces to c routines
!
integer(kind=c_int64_t) function h5pt_getnstepattribs (filehandle) bind(C,name="h5ptc_getnstepattribs")
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
end function

integer(kind=c_int64_t) function h5ptc_getstepattribinfo (filehandle,idx,name,nelem,l_name) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: idx              !< index of the attribute to query (starting from 0)
    character(kind=c_char), dimension(1), intent(out) :: name       !< buffer to read the attribute name into
    integer(kind=c_int64_t), intent(out) :: nelem             !< number of elements in the attribute's array
    integer(kind=c_int64_t), intent(in) :: l_name             !< number of elements in the attribute's array
end function

integer(kind=c_int64_t) function h5ptc_writestepattrib_r8 ( filehandle, name, data, nelem) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name !< the name of the attribute
    real(kind=c_double), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
end function

integer(kind=c_int64_t) function h5ptc_readstepattrib_r8 ( filehandle, name, data ) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name   !< the name of the attribute
    real(kind=c_double), intent(out) :: data(*) !< buffer to read value into
end function

integer(kind=c_int64_t) function h5ptc_writestepattrib_r4 ( filehandle, name, data, nelem) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name !< the name of the attribute
    real(kind=c_float), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
end function

integer(kind=c_int64_t) function h5ptc_readstepattrib_r4 ( filehandle, name, data ) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name   !< the name of the attribute
    real(kind=c_float), intent(out) :: data(*) !< buffer to read value into
end function

integer(kind=c_int64_t) function h5ptc_writestepattrib_i8 ( filehandle, name, data, nelem) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name !< the name of the attribute
    integer(kind=c_int64_t), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
end function

integer(kind=c_int64_t) function h5ptc_readstepattrib_i8 ( filehandle, name, data ) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name   !< the name of the attribute
    integer(kind=c_int64_t), intent(out) :: data(*) !< buffer to read value into
end function

integer(kind=c_int64_t) function h5ptc_writestepattrib_i4 ( filehandle, name, data, nelem) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name !< the name of the attribute
    integer(kind=c_int), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
end function

integer(kind=c_int64_t) function h5ptc_readstepattrib_i4 ( filehandle, name, data ) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(kind=c_char), dimension(1), intent(in) :: name   !< the name of the attribute
    integer(kind=c_int), intent(out) :: data(*) !< buffer to read value into
end function

!
! read/write string
!
integer(kind=c_int64_t) function h5ptc_writestepattrib_string (filehandle,name,value) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(kind=c_char), dimension(1), intent(in) :: name        !< the name of the attribute
    character(kind=c_char), dimension(1), intent(in) :: value       !< the string value to store
end function

integer(kind=c_int64_t) function h5ptc_readstepattrib_string (filehandle,name,value) bind(C)
    import
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(kind=c_char), dimension(1), intent(in) :: name        !< the name of the attribute
    character(kind=c_char), dimension(1), intent(out) :: value      !< buffer to read the string value into
end function

end interface

contains

!---------------------------------------------------------------------------
!
! wrappers for functions with string arguments: 
! converts strings into C strings
!
!---------------------------------------------------------------------------

!
! file attributes
!
integer(kind=c_int64_t) function h5pt_getfileattribinfo (filehandle,idx,name,nelem)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: idx              !< index of the attribute to query (starting from 0)
    character(len=*), intent(out) :: name       !< buffer to read the attribute name into
    integer(kind=c_int64_t), intent(out) :: nelem             !< number of elements in the attribute's array
    integer(kind=c_int64_t) :: l_name             !< size of name
    
    l_name = len(name)
    h5pt_getfileattribinfo = h5ptc_getfileattribinfo (filehandle,idx,name,nelem,l_name)
    name = fstring(name)
    
end function

integer(kind=c_int64_t) function h5pt_writefileattrib_r8 ( filehandle, name, data, nelem)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in)         :: name !< the name of the attribute
    real(kind=c_double), intent(in)     :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array

    h5pt_writefileattrib_r8 = h5ptc_writefileattrib_r8 ( filehandle, cstring(name), data, nelem)

end function h5pt_writefileattrib_r8

integer(kind=c_int64_t) function h5pt_readfileattrib_r8 ( filehandle, name, data )
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name   !< the name of the attribute
    real(kind=c_double), intent(out) :: data(*) !< buffer to read value into

    h5pt_readfileattrib_r8 = h5ptc_readfileattrib_r8 ( filehandle, cstring(name), data )

end function

integer(kind=c_int64_t) function h5pt_writefileattrib_r4 ( filehandle, name, data, nelem)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name !< the name of the attribute
    real(kind=c_float), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
    
    h5pt_writefileattrib_r4 = h5ptc_writefileattrib_r4( filehandle, cstring(name), data, nelem)    
    
end function

integer(kind=c_int64_t) function h5pt_readfileattrib_r4 ( filehandle, name, data )
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name   !< the name of the attribute
    real(kind=c_float), intent(out) :: data(*) !< buffer to read value into
    
    h5pt_readfileattrib_r4 = h5ptc_readfileattrib_r4( filehandle, cstring(name), data )
    
end function

integer(kind=c_int64_t) function h5pt_writefileattrib_i8 ( filehandle, name, data, nelem)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name !< the name of the attribute
    integer(kind=c_int64_t), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
    
    h5pt_writefileattrib_i8 = h5ptc_writefileattrib_i8( filehandle, cstring(name), data, nelem)

end function

integer(kind=c_int64_t) function h5pt_readfileattrib_i8 ( filehandle, name, data )
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name   !< the name of the attribute
    integer(kind=c_int64_t), intent(out) :: data(*) !< buffer to read value into
    
    h5pt_readfileattrib_i8 = h5ptc_readfileattrib_i8( filehandle, cstring(name), data )

end function

integer(kind=c_int64_t) function h5pt_writefileattrib_i4 ( filehandle, name, data, nelem)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name !< the name of the attribute
    integer(kind=c_int), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array

    h5pt_writefileattrib_i4 = h5ptc_writefileattrib_i4( filehandle, cstring(name), data, nelem)

end function

integer(kind=c_int64_t) function h5pt_readfileattrib_i4 ( filehandle, name, data )
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name   !< the name of the attribute
    integer(kind=c_int), intent(out) :: data(*) !< buffer to read value into
    
    h5pt_readfileattrib_i4 = h5ptc_readfileattrib_i4 ( filehandle, cstring(name), data )
     
end function

integer(kind=c_int64_t) function h5pt_writefileattrib_string ( filehandle, name, value)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(len=*), intent(in) :: name        !< the name of the attribute
    character(len=*), intent(in) :: value       !< the string value to store
    
    h5pt_writefileattrib_string = h5ptc_writefileattrib_string ( filehandle, cstring(name), cstring(value))

end function

integer(kind=c_int64_t) function h5pt_readfileattrib_string ( filehandle, name, value)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(len=*), intent(in) :: name        !< the name of the attribute
    character(len=*), intent(out) :: value      !< buffer to read the string value into
    
    h5pt_readfileattrib_string = h5ptc_readfileattrib_string ( filehandle, cstring(name), value)
    value = fstring(value)

end function

!
! step attributes:
!
integer(kind=c_int64_t) function h5pt_getstepattribinfo (filehandle,idx,name,nelem)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    integer(kind=c_int64_t), intent(in) :: idx              !< index of the attribute to query (starting from 0)
    character(len=*), intent(out) :: name       !< buffer to read the attribute name into
    integer(kind=c_int64_t), intent(out) :: nelem             !< number of elements in the attribute's array
    integer(kind=c_int64_t) :: l_name             !< number of elements in the attribute's array

    l_name = len(name)
    h5pt_getstepattribinfo = h5ptc_getstepattribinfo (filehandle,idx,name,nelem,l_name)
    name = fstring(name)

end function

integer(kind=c_int64_t) function h5pt_writestepattrib_r8 ( filehandle, name, data, nelem)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name !< the name of the attribute
    real(kind=c_double), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
    
    h5pt_writestepattrib_r8 = h5ptc_writestepattrib_r8 ( filehandle, cstring(name), data, nelem)
    
end function

integer(kind=c_int64_t) function h5pt_readstepattrib_r8 ( filehandle, name, data )
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name   !< the name of the attribute
    real(kind=c_double), intent(out) :: data(*) !< buffer to read value into
    
     h5pt_readstepattrib_r8 = h5ptc_readstepattrib_r8 ( filehandle, cstring(name), data )

end function

integer(kind=c_int64_t) function h5pt_writestepattrib_r4 ( filehandle, name, data, nelem)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name !< the name of the attribute
    real(kind=c_float), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
    
    h5pt_writestepattrib_r4 = h5ptc_writestepattrib_r4 ( filehandle, cstring(name), data, nelem)
    
end function

integer(kind=c_int64_t) function h5pt_readstepattrib_r4 ( filehandle, name, data )
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name   !< the name of the attribute
    real(kind=c_float), intent(out) :: data(*) !< buffer to read value into
    
    h5pt_readstepattrib_r4 = h5ptc_readstepattrib_r4 ( filehandle, cstring(name), data )
    
end function

integer(kind=c_int64_t) function h5pt_writestepattrib_i8 ( filehandle, name, data, nelem)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name !< the name of the attribute
    integer(kind=c_int64_t), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
    
    h5pt_writestepattrib_i8 = h5ptc_writestepattrib_i8 ( filehandle, cstring(name), data, nelem)

end function

integer(kind=c_int64_t) function h5pt_readstepattrib_i8 ( filehandle, name, data )
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name   !< the name of the attribute
    integer(kind=c_int64_t), intent(out) :: data(*) !< buffer to read value into
    
    h5pt_readstepattrib_i8 = h5ptc_readstepattrib_i8 ( filehandle, cstring(name), data )
    
end function

integer(kind=c_int64_t) function h5pt_writestepattrib_i4 ( filehandle, name, data, nelem)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name !< the name of the attribute
    integer(kind=c_int), intent(in) :: data(*) !< the array of data to write into the attribute
    integer(kind=c_int64_t), intent(in) :: nelem !< the number of elements in the array
    
    h5pt_writestepattrib_i4 = h5ptc_writestepattrib_i4 ( filehandle, cstring(name), data, nelem)
    
end function

integer(kind=c_int64_t) function h5pt_readstepattrib_i4 ( filehandle, name, data )
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle !< the handle returned at file open
    character(len=*), intent(in) :: name   !< the name of the attribute
    integer(kind=c_int), intent(out) :: data(*) !< buffer to read value into
    
    h5pt_readstepattrib_i4 = h5ptc_readstepattrib_i4 ( filehandle, cstring(name), data )
    
end function
!
! read/write string
!
integer(kind=c_int64_t) function h5pt_writestepattrib_string ( filehandle, name, value)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(len=*), intent(in) :: name        !< the name of the attribute
    character(len=*), intent(in) :: value       !< the string value to store
    
    h5pt_writestepattrib_string = h5ptc_writestepattrib_string ( filehandle, cstring(name), cstring(value))
    
end function

integer(kind=c_int64_t) function h5pt_readstepattrib_string ( filehandle, name, value)
    implicit none
    integer(kind=c_int64_t), intent(in) :: filehandle         !< the handle returned during file open
    character(len=*), intent(in) :: name        !< the name of the attribute
    character(len=*), intent(out) :: value      !< buffer to read the string value into
    
    h5pt_readstepattrib_string = h5ptc_readstepattrib_string ( filehandle, cstring(name), value)
    value = fstring(value)
    
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


end module h5partattrib
