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
!  Copyright (C) 2019- Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING FITS FILES
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! dat(maxpart,maxplot,maxstep) : main data array
!
! npartoftype(maxstep): number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!                      (used in calc_quantities for calculating the pressure)
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!
! Columns with the 'required' flag set to false are not read
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!
!  The routine that reads the data into splash's internal arrays
!
!-------------------------------------------------------------------------
subroutine read_data(rootname,istepstart,ipos,nstepsread)
 use particle_data,    only:dat,npartoftype,masstype,time,gamma,maxcol,maxpart
 use settings_data,    only:ndim,ndimV,ncolumns,ncalc,ipartialread,iverbose
 use mem_allocation,   only:alloc
 use readwrite_fits,   only:read_fits_image,fits_error,write_fits_image
 use imageutils,       only:image_denoise
 implicit none
 integer, intent(in)                :: istepstart,ipos
 integer, intent(out)               :: nstepsread
 character(len=*), intent(in)       :: rootname
 character(len=len(rootname)+10)    :: datfile
 integer               :: i,j,k,n,ierr,nextra,naxes(2)
 integer               :: ncolstep,npixels,nsteps_to_read,its,niter
 logical               :: iexist,reallocate
 real, dimension(:,:), allocatable :: image
 real :: dx,dy

 nstepsread = 0
 nsteps_to_read = 1

 if (len_trim(rootname) > 0) then
    datfile = trim(rootname)
 else
    print*,' **** no data read **** '
    return
 endif
!
!--check if first data file exists
!
 if (iverbose==1 .and. ipos==1) print "(1x,a)",'reading FITS format'
 inquire(file=datfile,exist=iexist)
 if (.not.iexist) then
    !
    !--append .fits on the endif not already present
    !
    datfile=trim(rootname)//'.fits'
    inquire(file=datfile,exist=iexist)
 endif
 if (.not.iexist) then
    print "(a)",' *** error: '//trim(rootname)//': file not found ***'
    return
 endif
!
!--set parameters which do not vary between timesteps
!
 ndim  = 2
 ndimV = 2
 nextra = 0
 ncolstep = 5 ! x, y, h, I, m
!
!--read data from snapshots
!
 i = istepstart
 write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
 !
 !--open file and read header information
 !
 call read_fits_image(datfile,image,naxes,ierr)
 if (ierr /= 0) then
    print*,'ERROR: '//trim(fits_error(ierr))
    if (allocated(image)) deallocate(image)
    return
 endif
 npixels = product(naxes)
 !
 !--allocate or reallocate memory for main data array
 !
 ncolumns = ncolstep + nextra
 nstepsread = 1
 reallocate = (npixels > maxpart)
 if (reallocate .or. .not.(allocated(dat))) then
    call alloc(npixels,nsteps_to_read,max(ncolumns+ncalc,maxcol),mixedtypes=.false.)
 endif
!
!--now memory has been allocated, set information that splash needs
!
 gamma = 5./3.
 time = 0.
 ipartialread = .false.
 masstype(1,i) = 0.0
 npartoftype(1,i) = npixels
!
! set x,y and other things needed for splash
!
 n = 0
 dx = 1.
 dy = 1.
 do k=1,naxes(2)
    do j=1,naxes(1)
       n = n + 1
       dat(n,1,i) = j
       dat(n,2,i) = k
       dat(n,3,i) = 1.  ! smoothing length == pixel scale
       dat(n,4,i) = image(j,k)
       dat(n,5,i) = image(j,k)*dx*dy  ! flux==equivalent of "mass"
    enddo
 enddo
!
! set smoothing length
!
! call image_denoise(naxes,image,dat(1:npixels,3,i))
!
! write smoothed fits image
!
! call write_fits_image('splash-output.fits',image,naxes,ierr)
!
! clean up
!
 deallocate(image)

end subroutine read_data

!------------------------------------------------------------
! set labels for each column of data
!------------------------------------------------------------
subroutine set_labels
 use labels,         only:label,labeltype,ix,ipmass,ih,irho
 use settings_data,  only:ndim,ndimV,ntypes,UseTypeInRenderings
 use geometry,       only:labelcoord
 implicit none

 if (ndim <= 0 .or. ndim > 3) then
    print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
    return
 endif
 if (ndimV <= 0 .or. ndimV > 3) then
    print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
    return
 endif

 ix(1) = 1
 ix(2) = 2
 ih = 3
 irho = 4
 ipmass = 5

 ! set labels of the quantities read in
 if (ix(1) > 0)   label(ix(1:ndim)) = labelcoord(1:ndim,1)
 if (irho > 0)    label(irho)       = 'intensity'
 if (ipmass > 0)  label(ipmass)     = 'flux'
 if (ih > 0)      label(ih)         = 'sigma'

 ! set labels for each particle type
 ntypes = 1
 labeltype(1) = 'pixel'
 UseTypeInRenderings(:) = .true.

end subroutine set_labels
