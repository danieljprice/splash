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
 public :: get_floats_from_fits_header,get_from_header

 interface write_fits_image
  module procedure write_fits_image,write_fits_image64
 end interface write_fits_image

 interface write_fits_cube
  module procedure write_fits_cube,write_fits_cube64
 end interface write_fits_cube

 private

contains

!---------------------------------------------------
! subroutine to read image from FITS file
! using cfitsio library
!---------------------------------------------------
subroutine read_fits_image(filename,image,naxes,ierr,hdr)
 character(len=*), intent(in)   :: filename
 real, intent(out), allocatable :: image(:,:)
 character(len=:), intent(inout), allocatable, optional :: hdr(:)
 integer, intent(out) :: naxes(2),ierr
 integer :: iunit,ireadwrite,npixels,blocksize
 integer :: firstpix,nullval,group,nfound,bitpix
 logical :: anynull
 real(kind=real32), allocatable :: img32(:,:)
 real(kind=real64), allocatable :: img64(:,:)
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
 bitpix = abs(get_from_header('BITPIX',hdr,ierr))
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
 !
 ! handle both 32 and 64 bit images, perform type conversions if they don't match
 !
 if (kind(image)==bitpix/8) then
    call ftgpve(iunit,group,firstpix,npixels,nullval,image,anynull,ierr)
 elseif (bitpix==64) then
    ! convert float64 image -> 32 bit internal array
    allocate(img64(naxes(1),naxes(2)),stat=ierr)
    call ftgpve(iunit,group,firstpix,npixels,nullval,img64,anynull,ierr)
    image = real(img64)
    deallocate(img64)
 elseif (bitpix==32) then ! perform type conversion
    ! convert float32 image -> 64 bit internal array
    allocate(img32(naxes(1),naxes(2)),stat=ierr)
    call ftgpve(iunit,group,firstpix,npixels,nullval,img32,anynull,ierr)
    image = real(img32)
    deallocate(img32)
 else ! should never happen, but just have a go anyway, maybe bitpix is corrupt?
    call ftclos(iunit,ierr)
    call ftfiou(iunit,ierr)
    ierr = 3
    return
 endif
 call ftgpve(iunit,group,firstpix,npixels,nullval,image,anynull,ierr)
 call ftclos(iunit,ierr)
 call ftfiou(iunit,ierr)

end subroutine read_fits_image

!---------------------------------------------------
! read FITS header from file
!---------------------------------------------------
subroutine read_fits_header(filename,hdr,ierr)
 character, intent(in)  :: filename
 character(len=:), allocatable, intent(out) :: hdr(:)
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
 character(len=:), allocatable, intent(inout) :: hdr(:)
 character(len=80) :: record
 integer :: i,nkeys,nspace

! The FTGHSP subroutine returns the number of existing keywords in the
! current header data unit (CHDU), not counting the required END keyword,
 call ftghsp(iunit,nkeys,nspace,ierr)
 !
 ! allocate memory
 !
 if (allocated(hdr)) deallocate(hdr)
 allocate(character(80) :: hdr(nkeys))

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
 morekeys = size(hdr)
 call fthdef(iunit,morekeys,ierr)
 do i=1,size(hdr)
    select case(hdr(i)(1:6))
    case('SIMPLE','BITPIX','NAXIS ','NAXIS1','NAXIS2','NAXIS3','NAXIS4','EXTEND')
       ! skip the above keywords
    case default
       call ftprec(iunit,hdr(i),ierr)
    end select
 enddo

end subroutine write_fits_head

!---------------------------------------------------
! subroutine to read spectral cube from FITS file
! using cfitsio library
!---------------------------------------------------
subroutine read_fits_cube(filename,image,naxes,ierr,hdr)
 character(len=*), intent(in)   :: filename
 real, intent(out), allocatable :: image(:,:,:)
 character(len=:), intent(inout), allocatable, optional :: hdr(:)
 integer, intent(out) :: naxes(4),ierr
 integer :: iunit,ireadwrite,npixels,blocksize
 integer :: firstpix,nullval,group,bitpix
 logical :: anynull
 integer :: ndim
 real(kind=real32), allocatable :: img32(:,:,:)
 real(kind=real64), allocatable :: img64(:,:,:)
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

 if (present(hdr)) call read_fits_head(iunit,hdr,ierr)

 call ftgidm(iunit,ndim,ierr) ! get_img_dim
 call ftgisz(iunit,3,naxes(1:ndim),ierr)
 if (ndim==2) naxes(3) = 1
 if (ndim>=3) ndim = 3

 bitpix = abs(get_from_header('BITPIX',hdr,ierr))

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
 !
 ! handle both 32 and 64 bit cubes, perform type conversions if they don't match
 !
 if (kind(image)==bitpix/8) then
    call ftgpve(iunit,group,firstpix,npixels,nullval,image,anynull,ierr)
 elseif (bitpix==64) then
    ! convert float64 image -> 32 bit internal array
    allocate(img64(naxes(1),naxes(2),naxes(3)),stat=ierr)
    call ftgpve(iunit,group,firstpix,npixels,nullval,img64,anynull,ierr)
    image = real(img64)
    deallocate(img64)
 elseif (bitpix==32) then ! perform type conversion
    ! convert float32 image -> 64 bit internal array
    allocate(img32(naxes(1),naxes(2),naxes(3)),stat=ierr)
    call ftgpve(iunit,group,firstpix,npixels,nullval,img32,anynull,ierr)
    image = real(img32)
    deallocate(img32)
 else ! should never happen, but just have a go anyway, maybe bitpix is corrupt?
    call ftclos(iunit,ierr)
    call ftfiou(iunit,ierr)
    ierr = 3
    return
 endif
 call ftclos(iunit,ierr)
 call ftfiou(iunit,ierr)

end subroutine read_fits_cube

!---------------------------------------------------
! error code handling
!---------------------------------------------------
 character(len=30) function fits_error(ierr)
  integer, intent(in) :: ierr

  select case(ierr)
  case(3)
     fits_error = 'could not match floating point type for fits image'
  case(2)
     fits_error = 'could not allocate memory'
  case(1)
     fits_error = 'no pixels found'
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

  img32 = image  ! copy and allocate
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
   print "(a)",' writing '//trim(filename)
   call ftinit(iunit,filename,blocksize,ierr)

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

 img32 = image  ! copy and allocate
 if (present(hdr)) then
    call write_fits_cube(filename,img32,naxes,ierr,hdr)
 else
    call write_fits_cube(filename,img32,naxes,ierr)
 endif
 deallocate(img32,stat=ierr)

end subroutine write_fits_cube64
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
subroutine get_fits_header_entry(record,tag,val,ierr)
 character(len=80), intent(in) :: record
 character(len=*),  intent(out) :: tag
 real, intent(out) :: val
 integer, intent(out) :: ierr
 integer :: ieq

 tag = ''
 val = 0.
 ierr = -1
 ! split on equals sign
 ieq = index(record,'=')
 if (ieq > 0) then
    tag = record(1:ieq-1)
    read(record(ieq+1:),*,iostat=ierr) val
 endif

end subroutine get_fits_header_entry

!------------------------------------------------
! search fits header to find a particular variable
! e.g. bmaj = get_from_header('BMAJ',hdr,ierr)
!------------------------------------------------
function get_from_header(tag,hdr,ierr) result(val)
 character(len=*),  intent(in) :: tag
 character(len=80), intent(in) :: hdr(:)
 integer,           intent(out) :: ierr
 character(len=len(tag)) :: mytag
 real :: val,myval
 integer :: i

 val = 0.
 ierr = -1
 do i=1,size(hdr)
    call get_fits_header_entry(hdr(i),mytag,myval,ierr)
    if (trim(adjustl(mytag))==trim(tag) .and. ierr==0) then
       val = myval
       ierr = 0
       return
    endif
 enddo

end function get_from_header

end module readwrite_fits
