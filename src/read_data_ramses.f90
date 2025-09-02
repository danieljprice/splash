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
!  Copyright (C) 2025- Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
!
! THIS VERSION IS FOR RAMSES ADAPTIVE MESH REFINEMENT OUTPUT
!
!-------------------------------------------------------------------------
module readdata_ramses
 use params, only:maxplot
 use labels, only:lenlabel
 implicit none

 public :: read_data_ramses, set_labels_ramses

 private
 character(len=lenlabel), dimension(maxplot) :: blocklabel

contains

subroutine read_data_ramses(rootname,istepstart,ipos,nstepsread)
 use particle_data,  only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol,maxstep
 use params,         only:doub_prec,maxplot
 use settings_data,  only:ndim,ndimV,ncolumns,ncalc,required,ipartialread, &
                           ntypes,iverbose,debugmode
 use mem_allocation, only:alloc
 use labels,         only:ih,irho,ipmass
 use system_utils,   only:renvironment,lenvironment,ienvironment,envlist
 use asciiutils,     only:cstring
 integer, intent(in)                :: istepstart,ipos
 integer, intent(out)               :: nstepsread
 character(len=*), intent(in)       :: rootname
 character(len=len(rootname)+10)    :: datfile
 integer               :: i,j,ierr
 integer               :: nhfac
 integer               :: ncolstep,npart_max,nstep_max,ntoti
 logical               :: iexist,reallocate,debug,goterrors
 real(doub_prec)       :: timetemp
 real :: hfact,hfactmean,pmassi
 real, parameter :: pi = 3.1415926536
 integer, dimension(maxplot) :: isrequired

 ! RAMSES-specific variables
 integer :: ncpu,nlevelmax,ngridmax,nboundary,ngrid_current
 real(doub_prec) :: boxlen
 integer :: nx,ny,nz

 nstepsread = 0
 goterrors  = .false.
 debug = debugmode

 if (len_trim(rootname) > 0) then
    datfile = trim(rootname)
 else
    print*,' **** no data read **** '
    return
 endif
!
!--check if first data file exists
!
 print "(1x,a)",'reading RAMSES AMR format'
 inquire(file=datfile,exist=iexist)
 if (.not.iexist) then
    print "(a)",' *** error: '//trim(rootname)//': file not found ***'
    return
 endif
!
!--set parameters which do not vary between timesteps
!
 ndim  = 3
 ndimV = 3
!
!--read data from snapshots
!
 i = istepstart
 write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
 !
 !--open file and read header information
 !
 if (debug) print*,'DEBUG: reading header...'
 call read_ramses_header(datfile,ntoti,ncolstep,ndim,ndimV,timetemp,ncpu,nlevelmax,ngridmax,&
                         nboundary,ngrid_current,boxlen,nx,ny,nz,ierr)
 if (ierr /= 0) then
    print "(a)", '*** ERROR READING HEADER ***'
    return
 endif
 ncolumns = ncolstep

 if (iverbose >= 1) print "(2(a,1x,i10))",' npart: ',ntoti,' ncolumns: ',ncolstep
 !
 !--now read data
 !
 reallocate = .false.
 npart_max = maxpart
 nstep_max = max(maxstep,1)

 if (ntoti > maxpart) then
    reallocate = .true.
    if (maxpart > 0) then
       ! if we are reallocating, try not to do it again
       npart_max = int(1.1*ntoti)
    else
       ! if first time, save on memory
       npart_max = int(ntoti)
    endif
 endif
 !
 !--reallocate memory for main data array
 !
 if (reallocate .or. .not.(allocated(dat))) then
    call alloc(npart_max,nstep_max,max(ncolumns+ncalc,maxcol))
 endif

 !
 !--copy header data into allocated arrays
 !
 npartoftype(1,i) = ntoti
 time(i) = real(timetemp)
 masstype(:,i) = 0. ! all masses read from file
 !
 !--read data
 !
 got_particles: if (ntoti > 0) then

    !call read_ramses_data(datfile,ntypes,npartoftype(:,i),ncolumns,required,ierr)

    nstepsread = 1
 endif got_particles
!
!--now memory has been allocated, set arrays which are constant for all time
!
 gamma = 5./3.
!
!--set flag to indicate that only part of this file has been read
!
 if (.not.all(required(1:ncolstep))) ipartialread = .true.
!
!--call set labels to identify location of smoothing length
!
 call set_labels_ramses
!
!--cover the special case where no particles have been read
!
 if (ntoti <= 0) then
    npartoftype(1,i) = 1
    dat(:,:,i) = 0.
 endif

 if (nstepsread > 0) then
    print "(a,i10,a)",' >> read ',sum(npartoftype(:,istepstart+nstepsread-1)),' particles'
 endif

end subroutine read_data_ramses

!!------------------------------------------------------------
!! read RAMSES header information from AMR and info files
!!------------------------------------------------------------

subroutine read_ramses_header(rootname,ntoti,ncolstep,ndim,ndimV,time,ncpu,nlevelmax,&
                              ngridmax,nboundary,ngrid_current,boxlen,nx,ny,nz,ierr)
 use params,       only:doub_prec
 use asciiutils,   only:cstring
 character(len=*), intent(in)    :: rootname
 integer, intent(out)            :: ntoti,ncolstep,ndim,ndimV,ierr
 integer, intent(out)            :: ncpu,nlevelmax,ngridmax,nboundary,ngrid_current
 integer, intent(out)            :: nx,ny,nz
 real(doub_prec), intent(out)    :: time,boxlen
 character(len=len(rootname)+20) :: amrfile,infofile
 character(len=5)                :: nchar
 character(len=80)               :: ordering
 character(len=13)               :: GMGM
 logical                         :: iexist
 integer                         :: i,ipos,levelmin
 real(doub_prec)                 :: aexp,scale_l,scale_d,scale_t
 
 ierr = 0
 
 ! Extract output number from filename
 ipos = INDEX(rootname,'output_')
 if (ipos > 0) then
    nchar = rootname(ipos+7:ipos+11)
 else
    ! Try to extract number from end of filename
    nchar = '00001'
 endif
 
 ! Read AMR file header
 amrfile = trim(rootname)//'/amr_'//trim(nchar)//'.out00001'
 inquire(file=amrfile,exist=iexist)
 if (.not.iexist) then
    print "(a)",' *** error: '//trim(amrfile)//': file not found ***'
    ierr = 1
    return
 endif
 
 print "(a)",' reading AMR header from '//trim(amrfile)
 open(unit=10,file=amrfile,status='old',form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",' *** error opening AMR file ***'
    return
 endif
 
 read(10,iostat=ierr) ncpu
 read(10,iostat=ierr) ndim
 read(10,iostat=ierr) nx,ny,nz
 read(10,iostat=ierr) nlevelmax
 read(10,iostat=ierr) ngridmax
 read(10,iostat=ierr) nboundary
 read(10,iostat=ierr) ngrid_current
 read(10,iostat=ierr) boxlen
 close(10)
 
 if (ierr /= 0) then
    print "(a)",' *** error reading AMR header ***'
    return
 endif
 
 print "(a,i0)",' ncpu = ',ncpu
 print "(a,i0)",' ndim = ',ndim
 print "(a,3(i0,1x))",' nx,ny,nz = ',nx,ny,nz
 print "(a,i0)",' nlevelmax = ',nlevelmax
 print "(a,g0)",' boxlen = ',boxlen
 
 ! Read info file for time information
 infofile = trim(rootname)//'/info_'//trim(nchar)//'.txt'
 inquire(file=infofile,exist=iexist)
 if (.not.iexist) then
    print "(a)",' *** warning: '//trim(infofile)//': file not found, setting time=0 ***'
    time = 0.0d0
 else
    print "(a)",' reading time from '//trim(infofile)
    open(unit=10,file=infofile,form='formatted',status='old',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' *** error opening info file ***'
       time = 0.0d0
       return
    endif
    
    ! Skip initial lines and read to time information
    do i=1,8
       read(10,*,iostat=ierr)
       if (ierr /= 0) exit
    enddo
    
    if (ierr == 0) then
       read(10,'(A13,E23.15)',iostat=ierr) GMGM,time
       if (ierr /= 0) then
          print "(a)",' *** error reading time from info file ***'
          time = 0.0d0
       else
          print "(a,g0)",' time = ',time
       endif
    endif
    
    close(10)
 endif
 
 ! Set up basic column structure for AMR cells treated as SPH particles
 ncolstep = 9  ! x,y,z,vx,vy,vz,h,rho,mass
 ndimV = 3
 
 ! Set up basic column labels
 blocklabel(1) = 'x'
 blocklabel(2) = 'y'
 blocklabel(3) = 'z'
 blocklabel(4) = 'vx'
 blocklabel(5) = 'vy'
 blocklabel(6) = 'vz'
 blocklabel(7) = 'h_smooth'
 blocklabel(8) = 'density'
 blocklabel(9) = 'mass'
 
 ! For now, set a placeholder number of particles
 ! This will be replaced when we actually read the AMR data
 ntoti = ngrid_current * 8  ! rough estimate: 8 cells per grid
 
end subroutine read_ramses_header

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels_ramses
 use labels,        only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass, &
                          ih,irho,iutherm,iax,make_vector_label
 use params
 use settings_data,  only:ndim,ndimV,UseTypeInRenderings
 use geometry,       only:labelcoord
 use system_utils,   only:envlist,ienvironment
 use asciiutils,     only:lcase
 integer :: icol

 if (ndim <= 0 .or. ndim > 3) then
    print*,'*** ERROR: ndim = ',ndim,' in set_labels_ramses ***'
    return
 endif
 if (ndimV <= 0 .or. ndimV > 3) then
    print*,'*** ERROR: ndimV = ',ndimV,' in set_labels_ramses ***'
    return
 endif

 ix = 0
 iutherm = 0
 do icol=1,size(blocklabel)
    select case(trim(lcase(blocklabel(icol))))
    case('x')
       ix(1) = icol
    case('y')
       ix(2) = icol
    case('z')
       ix(3) = icol
    case('vx')
       ivx = icol
    case('ax')
       iax = icol
    case('h_smooth')
       ih = icol
    case('mass')
       ipmass = icol
    case('density')
       irho = icol
    end select
    label(icol) = trim(blocklabel(icol))
 enddo

 ! set labels of the quantities read in
 if (ix(1) > 0)   label(ix(1:ndim)) = labelcoord(1:ndim,1)

 ! set labels for vector quantities
 call make_vector_label('v',ivx,ndimV,iamvec,labelvec,label,labelcoord(:,1))
 call make_vector_label('a',iax,ndimV,iamvec,labelvec,label,labelcoord(:,1))

 ! set labels for each particle type
 labeltype(1) = 'gas'
 UseTypeInRenderings(:) = .false.
 UseTypeInRenderings(1) = .true.

end subroutine set_labels_ramses

end module readdata_ramses
