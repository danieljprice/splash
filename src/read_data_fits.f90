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
 use particle_data,  only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol,maxstep,iamtype
 use params,         only:doub_prec
 use settings_data,  only:ndim,ndimV,ncolumns,ncalc,ipartialread,iverbose,buffer_steps_in_file
 !use settings_page,  only:legendtext
 use mem_allocation, only:alloc
 ! use labels,         only:ih,irho,ipmass
 use asciiutils,     only:cstring
 use fitsutils,      only:read_fits_header,read_fits_data
 implicit none
 integer, intent(in)                :: istepstart,ipos
 integer, intent(out)               :: nstepsread
 character(len=*), intent(in)       :: rootname
 character(len=len(rootname)+10)    :: datfile
 integer               :: i,istep,ierr,nextra
 integer               :: nunknown,ignoretl
 integer               :: ncolstep,npart_max,nstep_max,nsteps_to_read,ntoti
 logical               :: iexist,reallocate,goterrors,ignore_time_levels
 real(doub_prec) :: timetemp,dx,vol

 nstepsread = 0
 goterrors  = .false.

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
 ndim  = 3
 ndimV = 3
 nextra = 0
!
!--read data from snapshots
!
 i = istepstart
 write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
 !
 !--open file and read header information
 !
 call read_fits_header(cstring(datfile),ntoti,ncolstep,ierr)
 print*,'got ',ntoti,ncolstep,ierr
 !
 !--sanity check the header read
 !
 if (ntoti <= 0) then
    ierr = 1
    return
 endif
 read*

 ncolumns = ncolstep + nextra
 if (iverbose >= 1) print "(3(a,1x,i10))",' npart: ',ntoti,' ncolumns: ',ncolstep,' nsteps: ',nstep_max

 !
 !--reallocate memory for main data array
 !
 if (reallocate .or. .not.(allocated(dat))) then
    call alloc(npart_max,nsteps_to_read,max(ncolumns+ncalc,maxcol),mixedtypes=.true.)
 endif
 nstepsread = 1
!
!--now memory has been allocated, set arrays which are constant for all time
!
 gamma = 5./3.
!
!--set flag to indicate that only part of this file has been read
!
 ipartialread = .false.
!
!--call set labels to identify location of smoothing length
!
 !call set_labels
 call read_fits_data(ierr)

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
 use labels,         only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass,iutherm,ih,irho
 use params
 use settings_data,  only:ndim,ndimV,ntypes,UseTypeInRenderings
 use geometry,       only:labelcoord
 use system_utils,   only:envlist,ienvironment
!  use cactushdf5read, only:blocklabel
 use asciiutils,     only:lcase
 implicit none
 integer :: i,icol

 if (ndim <= 0 .or. ndim > 3) then
    print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
    return
 endif
 if (ndimV <= 0 .or. ndimV > 3) then
    print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
    return
 endif
 !blocklabel(1:5) = (/'x ','y ','z ','dx','m '/)

 ix(1) = 1
 ix(2) = 2
 ix(3) = 3
 ih = 4
 ipmass = 5
 iutherm = 0
 irho = 0
!   !do icol=1,size(blocklabel)
! !     select case(trim(blocklabel(icol)))
!      case('vel[0]')
!         ivx = icol
!      case('dens','density')
!         if (irho==0) irho = icol
!      case('rho')
!         irho = icol
!      end select
!      label(icol) = trim(blocklabel(icol))
!   enddo

 ! set labels of the quantities read in
 if (ix(1) > 0)   label(ix(1:ndim)) = labelcoord(1:ndim,1)
 !if (irho > 0)    label(irho)       = 'density'
 !if (iutherm > 0) label(iutherm)    = 'u'
 !if (ipmass > 0)  label(ipmass)     = 'particle mass'
 !if (ih > 0)      label(ih)         = 'h'

 ! set labels for vector quantities
 if (ivx > 0) then
    iamvec(ivx:ivx+ndimV-1) = ivx
    labelvec(ivx:ivx+ndimV-1) = 'v'
    do i=1,ndimV
       label(ivx+i-1) = trim(labelvec(ivx))//'_'//labelcoord(i,1)
    enddo
 endif

 ! set labels for each particle type
 ntypes = 2
 labeltype(1) = 'gas'
 labeltype(2) = 'ghost'
 UseTypeInRenderings(:) = .true.
 UseTypeInRenderings(2) = .false.

 return
end subroutine set_labels
