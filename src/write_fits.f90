!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2020- Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
!----------------------------------------------------------------------
!
!  Module handling read and write of FITS files
!  With thanks to Christophe Pinte
!
!----------------------------------------------------------------------
module readwrite_fits
 use iso_fortran_env, only:real32,real64
 implicit none
 public :: read_fits_image,write_fits_image,fits_error
 public :: read_fits_cube,write_fits_cube
 public :: read_fits_header
 public :: get_floats_from_fits_header,get_from_header,get_from_header_s
 public :: append_to_fits_cube
 public :: get_velocity_from_fits_header
 public :: flatten_header

 interface write_fits_image
  module procedure write_fits_image,write_fits_image64
 end interface write_fits_image

 interface write_fits_cube
  module procedure write_fits_cube,write_fits_cube64
 end interface write_fits_cube

 interface read_fits_image
  module procedure read_fits_image,read_fits_image64
 end interface read_fits_image

 interface read_fits_cube
  module procedure read_fits_cube,read_fits_cube64
 end interface read_fits_cube

 interface append_to_fits_cube
  module procedure append_to_fits_cube
 end interface append_to_fits_cube

 private

contains

!---------------------------------------------------
! subroutine to read image from FITS file
! using cfitsio library
!---------------------------------------------------
subroutine read_fits_image(filename,image,naxes,ierr,hdr)
 character(len=*), intent(in)   :: filename
 real(kind=real32), intent(out), allocatable :: image(:,:)
 character(len=80), intent(inout), allocatable, optional :: hdr(:)
 integer, intent(out) :: naxes(2),ierr
 integer :: iunit,ireadwrite,npixels,blocksize
 integer :: firstpix,nullval,group,nfound!,bitpix
 logical :: anynull
 !
 !--open file and read header information
 !
 ierr = 0
 call ftgiou(iunit,ierr)

 ireadwrite = 0
 call ftopen(iunit,filename,ireadwrite,blocksize,ierr)

 if (ierr /= 0) then
    ierr = -1
    return
 endif
 !
 !--read fits header (this is optional)
 !
 if (present(hdr)) call read_fits_head(iunit,hdr,ierr)
 !
 !--get the essential things from the header
 !
 call ftgknj(iunit,'NAXIS',1,2,naxes,nfound,ierr)
 npixels = naxes(1)*naxes(2)
 !if (present(hdr)) bitpix = abs(get_from_header('BITPIX',hdr,ierr))
 !
 !--sanity check the header read
 !
 if (npixels <= 0) then
    !print*,' ERROR: No pixels found'
    ierr = 1
    return
 endif
 !
 ! read image
 !
 firstpix = 1
 nullval = -999
 group = 1
 allocate(image(naxes(1),naxes(2)),stat=ierr)
 if (ierr /= 0) then
    ierr = 2
    return
 endif
 ierr = 0
 call ftgpve(iunit,group,firstpix,npixels,nullval,image,anynull,ierr)
 call ftclos(iunit,ierr)
 call ftfiou(iunit,ierr)

end subroutine read_fits_image

!---------------------------------------------------
! read FITS header from file
!---------------------------------------------------
subroutine read_fits_header(filename,hdr,ierr)
 character, intent(in)  :: filename
 character(len=80), allocatable, intent(out) :: hdr(:)
 integer, intent(out) :: ierr
 integer :: ireadwrite,iunit,blocksize

 ierr = 0
 call ftgiou(iunit,ierr)

 ireadwrite = 0
 call ftopen(iunit,filename,ireadwrite,blocksize,ierr)
 if (ierr /= 0) return
 !
 !--read fits header (this is optional)
 !
 call read_fits_head(iunit,hdr,ierr)
 call ftclos(iunit,ierr)
 call ftfiou(iunit,ierr)

end subroutine read_fits_header

!---------------------------------------------------
! internal subroutine to read FITS header information
!---------------------------------------------------
subroutine read_fits_head(iunit,hdr,ierr)
 integer, intent(in)  :: iunit
 integer, intent(out) :: ierr
 character(len=80), allocatable, intent(inout) :: hdr(:)
 character(len=80) :: record
 integer :: i,nkeys,nspace

! The FTGHSP subroutine returns the number of existing keywords in the
! current header data unit (CHDU), not counting the required END keyword,
 call ftghsp(iunit,nkeys,nspace,ierr)
 !
 ! allocate memory
 !
 if (allocated(hdr)) deallocate(hdr)
 allocate(hdr(nkeys))

! Read each 80-character keyword record, and print it out.
 do i = 1, nkeys
    call ftgrec(iunit,i,record,ierr)
    hdr(i) = record
    !print *,hdr(i)
 end do

end subroutine read_fits_head

!---------------------------------------------------
! internal subroutine to write FITS header information
! excluding things we have changed
!---------------------------------------------------
subroutine write_fits_head(iunit,hdr,ierr)
 integer, intent(in) :: iunit
 character(len=80), intent(in) :: hdr(:)
 integer, intent(out) :: ierr
 integer :: i,morekeys

 ierr = 0
 morekeys = 0
 ! count non-blank header lines
 do i=1,size(hdr)
    if (len_trim(hdr(i)) > 0) then
       morekeys = morekeys + 1
    endif
 enddo
 ! write header size to fits file
 call fthdef(iunit,morekeys,ierr)
 if (ierr /= 0) return

 ! write header
 do i=1,size(hdr)
    if (len_trim(hdr(i)) > 0) then
       select case(hdr(i)(1:6))
       case('SIMPLE','BITPIX','NAXIS ','NAXIS1','NAXIS2','NAXIS3','NAXIS4','EXTEND')
          ! skip the above keywords
       case default
          call ftprec(iunit,hdr(i),ierr)
       end select
    endif
 enddo

end subroutine write_fits_head

!---------------------------------------------------
! subroutine to read spectral cube from FITS file
! using cfitsio library
!---------------------------------------------------
subroutine read_fits_cube(filename,image,naxes,ierr,hdr,hdu,velocity)
 character(len=*), intent(in)   :: filename
 real(kind=real32), intent(out), allocatable :: image(:,:,:)
 character(len=80), intent(inout), allocatable, optional :: hdr(:)
 real(kind=real32), intent(out), allocatable, optional :: velocity(:)
 integer, intent(out) :: naxes(4),ierr
 integer, intent(in), optional :: hdu ! specify which hdu to read
 integer :: iunit,ireadwrite,npixels,blocksize
 integer :: firstpix,nullval,group,hdutype
 logical :: anynull
 integer :: ndim
 !
 !--open file and read header information
 !
 ierr = 0
 call ftgiou(iunit,ierr)

 ireadwrite = 0
 call ftopen(iunit,filename,ireadwrite,blocksize,ierr)
 if (ierr /= 0) then
    ierr = -1
    return
 endif
 !
 ! switch to specified hdu, if argument is given
 !
 if (present(hdu)) then
    print*,' reading hdu ',hdu
    call ftmahd(iunit,hdu,hdutype,ierr)   ! hdutype==0
    if (ierr /= 0 .or. hdutype /= 0) then
       ierr = -2
       return
    endif
 endif

 if (present(hdr)) call read_fits_head(iunit,hdr,ierr)

 call ftgidm(iunit,ndim,ierr) ! get_img_dim
 if (ndim>=3) ndim = 3
 call ftgisz(iunit,3,naxes(1:ndim),ierr)
 if (ndim==2) naxes(3) = 1

 if (present(hdr) .and. present(velocity) .and. ndim >= 3) then
    if (.not.allocated(velocity)) allocate(velocity(naxes(3)))
    call get_velocity_from_fits_header(naxes(3),velocity,hdr,ierr)
 endif

 !if (present(hdr)) bitpix = abs(get_from_header('BITPIX',hdr,ierr))

 ! call ftgknj(iunit,'NAXIS',1,2,naxes,nfound,ierr)
 npixels = product(naxes(1:ndim))
 !
 ! sanity check the header read
 !
 if (npixels <= 0) then
    ierr = 1
    return
 endif
 !
 ! read images
 !
 firstpix = 1
 nullval = -999
 group = 1
 allocate(image(naxes(1),naxes(2),naxes(3)),stat=ierr)
 if (ierr /= 0) then
    ierr = 2
    return
 endif
 ierr = 0
 call ftgpve(iunit,group,firstpix,npixels,nullval,image,anynull,ierr)
 call ftclos(iunit,ierr)
 call ftfiou(iunit,ierr)

end subroutine read_fits_cube

!---------------------------------------------------
! error code handling
!---------------------------------------------------
character(len=50) function fits_error(ierr)
 integer, intent(in) :: ierr

 select case(ierr)
 case(3)
    fits_error = 'could not match floating point type for fits image'
 case(2)
    fits_error = 'could not allocate memory'
 case(1)
    fits_error = 'no pixels found'
 case(-2)
    fits_error = 'could not open specified hdu in fits file'
 case(-1)
    fits_error = 'could not open fits file'
 case default
    fits_error = 'unknown error'
 end select

end function fits_error

!------------------------------------------------
! Writing new fits file
!------------------------------------------------
subroutine write_fits_image(filename,image,naxes,ierr,hdr)
 character(len=*), intent(in) :: filename
 integer, intent(in)  :: naxes(2)
 real(kind=real32), intent(in) :: image(naxes(1),naxes(2))
 integer, intent(out) :: ierr
 character(len=80), intent(in), optional :: hdr(:)
 integer :: iunit,blocksize,group,firstpixel,bitpix,npixels
 logical :: simple,extend

 !  Get an unused Logical Unit Number to use to open the FITS file.
 ierr = 0
 call ftgiou(iunit,ierr)

 !  Create the new empty FITS file.
 blocksize=1
 print "(a)",' writing '//trim(filename)
 call ftinit(iunit,filename,blocksize,ierr)

 !  Initialize parameters about the FITS image
 simple=.true.
 ! data size
 bitpix=-32
 extend=.true.

 !  Write the required header keywords.
 call ftphpr(iunit,simple,bitpix,2,naxes,0,1,extend,ierr)
 !  Write additional header keywords, if present
 if (present(hdr)) call write_fits_head(iunit,hdr,ierr)

 group=1
 firstpixel=1
 npixels = naxes(1)*naxes(2)
 ! write as real*4
 call ftppre(iunit,group,firstpixel,npixels,image,ierr)

 !  Close the file and free the unit number
 call ftclos(iunit, ierr)
 call ftfiou(iunit, ierr)

end subroutine write_fits_image

!-------------------------------------------------------------
! Writing new fits file (convert from double precision input)
!-------------------------------------------------------------
subroutine write_fits_image64(filename,image,naxes,ierr,hdr)
 character(len=*), intent(in) :: filename
 integer,          intent(in) :: naxes(2)
 real(kind=real64),intent(in) :: image(naxes(1),naxes(2))
 real(kind=real32), allocatable :: img32(:,:)
 integer, intent(out) :: ierr
 character(len=80), intent(in), optional :: hdr(:)

 img32 = real(image,kind=real32)  ! copy and allocate
 if (present(hdr)) then
    call write_fits_image(filename,img32,naxes,ierr,hdr)
 else
    call write_fits_image(filename,img32,naxes,ierr)
 endif
 deallocate(img32,stat=ierr)

end subroutine write_fits_image64

!------------------------------------------------
! Writing new fits file
!------------------------------------------------
subroutine write_fits_cube(filename,image,naxes,ierr,hdr)
 character(len=*), intent(in) :: filename
 integer, intent(in)  :: naxes(3)
 real(kind=real32), intent(in) :: image(naxes(1),naxes(2),naxes(3))
 integer, intent(out) :: ierr
 character(len=80), intent(in), optional :: hdr(:)
 integer :: iunit,blocksize,group,firstpixel,bitpix,npixels
 logical :: simple,extend

 !  Get an unused Logical Unit Number to use to open the FITS file.
 ierr = 0
 call ftgiou(iunit,ierr)

 !  Create the new empty FITS file.
 blocksize=1
 call ftinit(iunit,filename,blocksize,ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR: '//trim(filename)//' already exists or is not writeable'
    return
 else
    print "(a)",' writing '//trim(filename)
 endif

 !  Initialize parameters about the FITS image
 simple=.true.
 ! data size
 bitpix=-32
 extend=.true.

 !  Write the required header keywords.
 call ftphpr(iunit,simple,bitpix,3,naxes,0,1,extend,ierr)
 !  Write additional header keywords, if present
 if (present(hdr)) call write_fits_head(iunit,hdr,ierr)

 group=1
 firstpixel=1
 npixels = product(naxes)
 ! write as real*4
 call ftppre(iunit,group,firstpixel,npixels,image,ierr)

 !  Close the file and free the unit number
 call ftclos(iunit, ierr)
 call ftfiou(iunit, ierr)

end subroutine write_fits_cube


!------------------------------------------------
! Writing new fits file
!------------------------------------------------
subroutine append_to_fits_cube(filename,image,naxes,ierr,hdr)
 character(len=*), intent(in) :: filename
 integer, intent(in)  :: naxes(3)
 real(kind=real32), intent(in) :: image(naxes(1),naxes(2),naxes(3))
 integer, intent(out) :: ierr
 character(len=80), intent(in), optional :: hdr(:)
 integer :: iunit,blocksize,group,firstpixel,bitpix,npixels,ireadwrite
 logical :: simple,extend

 !  Get an unused Logical Unit Number to use to open the FITS file.
 ierr = 0
 call ftgiou(iunit,ierr)

 !  Open the FITS file for read/write
 blocksize=1
 ireadwrite=1
 call ftopen(iunit,filename,ireadwrite,blocksize,ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR: '//trim(filename)//' does not exist or is not read/writeable'
    return
 else
    print "(a)",' appending to '//trim(filename)
 endif

 !--create new hdu in the fits file
 call ftcrhd(iunit,ierr)
 if (ierr /= 0) then
    print "(a)",'ERROR creating new hdu in '//trim(filename)
    return
 endif

 !  Initialize parameters about the FITS image
 simple=.true.
 ! data size
 bitpix=-32
 extend=.true.

 !  Write the required header keywords.
 call ftphpr(iunit,simple,bitpix,3,naxes,0,1,extend,ierr)
 !  Write additional header keywords, if present
 if (present(hdr)) call write_fits_head(iunit,hdr,ierr)

 group=1
 firstpixel=1
 npixels = product(naxes)
 ! write as real*4
 call ftppre(iunit,group,firstpixel,npixels,image,ierr)

 !  Close the file and free the unit number
 call ftclos(iunit, ierr)
 call ftfiou(iunit, ierr)

end subroutine append_to_fits_cube

!-------------------------------------------------------------
! Writing new fits file (convert from double precision input)
!-------------------------------------------------------------
subroutine write_fits_cube64(filename,image,naxes,ierr,hdr)
 character(len=*), intent(in) :: filename
 integer,          intent(in) :: naxes(3)
 real(kind=real64),intent(in) :: image(naxes(1),naxes(2),naxes(3))
 character(len=80), intent(in), optional :: hdr(:)
 integer,           intent(out) :: ierr
 real(kind=real32), allocatable :: img32(:,:,:)

 img32 = real(image,kind=real32)  ! copy and allocate
 if (present(hdr)) then
    call write_fits_cube(filename,img32,naxes,ierr,hdr)
 else
    call write_fits_cube(filename,img32,naxes,ierr)
 endif
 deallocate(img32,stat=ierr)

end subroutine write_fits_cube64

!-------------------------------------------------------------
! read fits file and convert to double precision
!-------------------------------------------------------------
subroutine read_fits_image64(filename,image,naxes,ierr,hdr)
 character(len=*), intent(in)   :: filename
 real(kind=real64), intent(out), allocatable :: image(:,:)
 character(len=80), intent(inout), allocatable, optional :: hdr(:)
 integer, intent(out) :: naxes(2),ierr
 real(kind=real32), allocatable :: img32(:,:)

 if (present(hdr)) then
    call read_fits_image(filename,img32,naxes,ierr,hdr)
 else
    call read_fits_image(filename,img32,naxes,ierr)
 endif
 image = img32  ! allocate and copy, converting real type
 deallocate(img32,stat=ierr)

end subroutine read_fits_image64

!-------------------------------------------------------------
! read fits cube and convert to double precision
!-------------------------------------------------------------
subroutine read_fits_cube64(filename,image,naxes,ierr,hdr,velocity)
 character(len=*), intent(in)   :: filename
 real(kind=real64), intent(out), allocatable :: image(:,:,:)
 character(len=80), intent(inout), allocatable, optional :: hdr(:)
 real(kind=real64), intent(inout), allocatable, optional :: velocity(:)
 integer, intent(out) :: naxes(4),ierr
 real(kind=real32), allocatable :: img32(:,:,:),velocity32(:)

 if (present(hdr)) then
    if (present(velocity)) then
       call read_fits_cube(filename,img32,naxes,ierr,hdr=hdr,velocity=velocity32)
       velocity = velocity32  ! allocate and copy, converting real type
    else
       call read_fits_cube(filename,img32,naxes,ierr,hdr)
    endif
 else
    call read_fits_cube(filename,img32,naxes,ierr)
 endif
 image = img32  ! allocate and copy, converting real type
 deallocate(img32,stat=ierr)
 if (allocated(velocity32)) deallocate(velocity32,stat=ierr)

end subroutine read_fits_cube64

!--------------------------------------------------
! read all floating point variables from fits header
!--------------------------------------------------
subroutine get_floats_from_fits_header(hdr,tags,vals)
 character(len=80), intent(in) :: hdr(:)
 character(len=*),  intent(out) :: tags(:)
 real,              intent(out) :: vals(:)
 integer :: i, n, ierr

 n = 0
 do i=1,size(hdr)
    n = n + 1
    if (n <= size(tags) .and. n <= size(vals)) then
       call get_fits_header_entry(hdr(i),tags(n),vals(n),ierr)
    endif
    if (ierr /= 0) n = n - 1
 enddo

end subroutine get_floats_from_fits_header

!------------------------------------------------
! get tag:val pairs from fits header record
! will extract anything readable as a floating
! point number
!------------------------------------------------
subroutine get_fits_header_entry(record,key,rval,ierr)
 character(len=80), intent(in) :: record
 character(len=*),  intent(out) :: key
 real, intent(out) :: rval
 integer, intent(out) :: ierr
 character(len=80) :: val

 call get_fits_header_key_val(record,key,val,ierr)
 if (ierr == 0) then
    read(val,*,iostat=ierr) rval
 endif

end subroutine get_fits_header_entry

!------------------------------------------------
! get tag:val pairs from fits header record
! returns the string value
!------------------------------------------------
subroutine get_fits_header_key_val(record,key,val,ierr)
 character(len=80), intent(in) :: record
 character(len=*),  intent(out) :: key
 character(len=*),  intent(out) :: val
 integer, intent(out) :: ierr
 integer :: ieq

 key = ''
 val = ''
 ierr = -1
 ! split on equals sign
 ieq = index(record,'=')
 if (ieq > 0) then
    key = record(1:ieq-1)
    val = record(ieq+1:)
    ierr = 0
 endif

end subroutine get_fits_header_key_val

!------------------------------------------------
! search fits header to find a particular variable
! e.g. bmaj = get_from_header('BMAJ',hdr,ierr)
!------------------------------------------------
function get_from_header(key,hdr,ierr) result(val)
 character(len=*),  intent(in) :: key
 character(len=80), intent(in) :: hdr(:)
 integer,           intent(out) :: ierr
 character(len=len(key)) :: mykey
 real :: val,myval
 integer :: i

 val = 0.
 ierr = -1
 do i=1,size(hdr)
    call get_fits_header_entry(hdr(i),mykey,myval,ierr)
    if (trim(adjustl(mykey))==trim(key) .and. ierr==0) then
       val = myval
       ierr = 0
       return
    endif
 enddo

end function get_from_header

!------------------------------------------------
! search fits header to find a particular string
! e.g. bmaj = get_from_header('BMAJ',hdr,ierr)
!------------------------------------------------
function get_from_header_s(key,hdr,ierr) result(val)
 character(len=*),  intent(in) :: key
 character(len=80), intent(in) :: hdr(:)
 integer,           intent(out) :: ierr
 character(len=len(key)) :: mykey
 character(len=80) :: val,myval
 integer :: i,i1,i2,islash

 val = ''
 ierr = -1
 do i=1,size(hdr)
    call get_fits_header_key_val(hdr(i),mykey,myval,ierr)
    if (trim(adjustl(mykey))==trim(key) .and. ierr==0) then
       i1 = index(myval,"'")
       islash = index(myval,"/")
       if (islash == 0) islash = len_trim(myval)
       i2 = index(myval(1:islash),"'",back=.true.)
       if (i1 > 0) then
          val = myval(i1+1:i2-1)  ! trim quotation marks
       else
          val = myval ! not a string variable, return whole string
       endif
       ierr = 0
       return
    endif
 enddo

end function get_from_header_s

!------------------------------------------------
! delete third dimension in the fits header
!------------------------------------------------
subroutine flatten_header(hdr)
 character(len=80), intent(inout) :: hdr(:)
 character(len=80) :: mykey,myval
 integer :: i,ierr

 do i=1,size(hdr)
    call get_fits_header_key_val(hdr(i),mykey,myval,ierr)
    ! delete anything in the header that ends in '3'
    if (len_trim(mykey) > 0) then
       if (mykey(len_trim(mykey):len_trim(mykey))=='3') then
          !print*,' deleting ',hdr(i)
          hdr(i) = ''
       endif
    endif
 enddo

end subroutine flatten_header

!------------------------------------------------
! get velocity grid from the fits file
!------------------------------------------------
subroutine get_velocity_from_fits_header(nv,vel,hdr,ierr)
 integer, intent(in)  :: nv
 real(kind=real32), intent(out) :: vel(nv)
 character(len=80), intent(in)  :: hdr(:)
 integer, intent(out) :: ierr
 integer :: k
 real :: dv,restfreq,crpix,crval,nu
 character(len=80) :: velocity_type
 real, parameter :: c = 2.997924e5

 dv = get_from_header('CDELT3',hdr,ierr)
 restfreq = get_from_header('RESTFRQ',hdr,ierr)
 velocity_type = get_from_header_s('CTYPE3',hdr,ierr)
 select case (trim(velocity_type))
 case("VELO-LSR","VRAD")
    if (dv > 10.) then
       dv = dv*1e-3 ! assume m/s
    endif
    crpix = get_from_header('CRPIX3',hdr,ierr)
    crval = get_from_header('CRVAL3',hdr,ierr)
    do k=1,nv
       vel(k) = real(crval + dv * (k - crpix),kind=real32)
    enddo
 case("FREQ")
    crpix = get_from_header('CRPIX3',hdr,ierr)
    crval = get_from_header('CRVAL3',hdr,ierr)
    do k=1,nv
       nu = crval + dv*(k - crpix)
       vel(k) = real(-(nu - restfreq)/restfreq * c,kind=real32)  ! c is in km/s
    enddo
 case default
    write(*,"(1x,a)") 'warning: velocity type '//trim(velocity_type)//' not recognised in fits header'
    dv = 1.
    do k=1,nv
       vel(k) = real(dv * (k - 0.5),kind=real32)
    enddo
 end select

end subroutine get_velocity_from_fits_header

end module readwrite_fits
