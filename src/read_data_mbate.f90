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
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM MATTHEW BATE'S CODE
! (ie. STRAIGHT FROM THE DATA DUMP)
!
! *** CONVERTS TO SINGLE PRECISION ***
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! maxplot,maxpart,maxstep      : dimensions of main data array
! dat(maxplot,maxpart,maxstep) : main data array
!
! npartoftype(1:6,maxstep) : number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!-------------------------------------------------------------------------

module readdata_mbate
 implicit none

 public :: read_data_mbate, set_labels_mbate

 private

contains

subroutine read_data_mbate(rootname,indexstart,ipos,nstepsread)
 use particle_data
 use params
 use settings_data,  only:ndim,ndimV,ncolumns,ncalc,buffer_steps_in_file
 use settings_units, only:units,get_nearest_length_unit,get_nearest_time_unit
 use labels,         only:unitslabel,headertags
 use mem_allocation
 integer, intent(in) :: indexstart,ipos
 integer, intent(out) :: nstepsread
 character(len=*), intent(in) :: rootname
 integer, parameter :: maxptmass = 1000
 real, parameter :: pi=4.*atan(1.)
 integer :: i,j,ifile,ierr,ipart
 integer :: npart_max,nstep_max,ncolstep,npart,nptmassi,nunknown
 logical :: iexist,doubleprec

 character(len=len(rootname)+10) :: dumpfile
 integer :: nprint, n1, n2, nptmass, nstepsalloc
 integer, dimension(:), allocatable :: isteps, iphase
 integer, dimension(maxptmass) :: listpm

 !--use these lines if dump is double precision
 real(doub_prec), dimension(:,:), allocatable :: dattemp
 real(doub_prec) :: udisti,umassi,utimei
 real(doub_prec) :: timei, gammai
 real(doub_prec) :: rhozero, RK2
 real(doub_prec) :: escap,tkin,tgrav,tterm
 real(doub_prec) :: dtmax,unitx
 character(len=20) :: unitlabelx

 !--use these lines for single precision
 real, dimension(:,:), allocatable :: dattemps
 real :: timesi, gammasi
 real :: rhozeros, RK2s
 real :: escaps,tkins,tgravs,tterms
 real :: dtmaxs,tcomp

 nstepsread = 0
 nstep_max = 0
 npart_max = maxpart
 ifile = 1
 buffer_steps_in_file = .true.

 dumpfile = trim(rootname)
 !
 !--check if data file exists
 !
 inquire(file=dumpfile,exist=iexist)
 if (.not.iexist) then
    print "(a)",' *** error: '//trim(dumpfile)//': file not found ***'
    return
 endif
 !
 !--fix number of spatial dimensions
 !
 ndim = 3
 ndimV = 3
 ncolstep = 11  ! number of columns in file
 ncolumns = ncolstep
 !
 !--allocate memory initially
 !
 nstep_max = max(nstep_max,indexstart,1)

 j = indexstart
 nstepsread = 0

 print "(1x,a)",'reading Matthew Bate''s/Willy Benz''s old SPH code format (on unit 15)'
 write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
 !
 !--open the (unformatted) binary file and read the number of particles
 !
 open(unit=15,iostat=ierr,file=dumpfile,status='old',form='unformatted')
 if (ierr /= 0) then
    print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
 else
    !
    !--read the number of particles in the first step,
    !  allocate memory and rewind
    !
    read(15,end=55,iostat=ierr) udisti,umassi,utimei,nprint,n1,n2,timei,gammai,rhozero,RK2
    doubleprec = .true.
    !--try single precision if non-sensible values for time, gamma etc.
    if (ierr /= 0 .or. timei < 0. .or. timei > 1e30  &
          .or. gammai < 1. .or. gammai > 10. &
          .or. rhozero < 0. .or. RK2 < 0. .or. RK2 > 1.e10) then
       doubleprec = .false.
    endif
    nstepsalloc = 1
    do while (ierr == 0)
       npart_max = max(npart_max,nprint)
       nstepsalloc = nstepsalloc + 1
       read(15,iostat=ierr) udisti,umassi,utimei,nprint
    enddo
    ierr = 0
    if (.not.allocated(dat) .or. nprint > npart_max) then
       npart_max = max(npart_max,INT(1.1*nprint))
       call alloc(npart_max,nstepsalloc,ncolstep+ncalc)
    endif

    rewind(15)
 endif
 if (ierr /= 0) then
    print "(a)",'*** ERROR READING TIMESTEP HEADER ***'
 else
!
!--loop over the timesteps in this file
!
    over_steps_in_file: do
       !npart_max = max(npart_max,nprint)
!
!--allocate/reallocate memory if j > maxstep
!
       if (j > maxstep) then
          !if (nstepsread > 2) then
          !   nstepsalloc = j + 2*nstepsread
          !else
          !   nstepsalloc = j
          !endif
          call alloc(maxpart,nstepsalloc,maxcol)
       endif
!
!--allocate integer arrays required for data read
!
       if (allocated(isteps)) deallocate(isteps)
       allocate(isteps(npart_max),stat=ierr)
       if (ierr /= 0) print*,'not enough memory in read_data_mbate'

       if (allocated(iphase)) deallocate(iphase)
       allocate(iphase(npart_max),stat=ierr)
       if (ierr /= 0) print*,'not enough memory in read_data_mbate'
!
!--now read the timestep data in the dumpfile
!
       write(*,"(a,i5,a)",advance="no") '| step ',j,': '

       if (doubleprec) then
          print "(a)",'double precision dump'
          !
          !--allocate a temporary array for (double precision) variables
          !
          if (allocated(dattemp)) deallocate(dattemp)
          allocate(dattemp(npart_max,ncolstep),stat=ierr)
          if (ierr /= 0) print*,'not enough memory in read_data_mbate'

          read(15,end=55,iostat=ierr) udisti, umassi, utimei, &
             nprint, n1, n2, timei, gammai, rhozero, RK2, &
             (dattemp(i,7), i=1, nprint),escap, tkin, tgrav, tterm, &
             (dattemp(i,1), i=1, nprint), (dattemp(i,2), i=1, nprint), &
             (dattemp(i,3), i=1, nprint), (dattemp(i,4), i=1, nprint), &
             (dattemp(i,5), i=1, nprint), (dattemp(i,6), i=1, nprint), &
             (dattemp(i,8), i=1, nprint), (dattemp(i,9), i=1, nprint), &
             (dattemp(i,10), i=1, nprint), (dattemp(i,11),i=1,nprint), &
             dtmax, (isteps(i), i=1,nprint), (iphase(i),i=1,nprint), &
             nptmass, (listpm(i), i=1,nptmass)
       else
          !
          !--allocate a temporary array for (double precision) variables
          !
          if (allocated(dattemps)) deallocate(dattemps)
          allocate(dattemps(npart_max,ncolstep),stat=ierr)
          if (ierr /= 0) print*,'not enough memory in read_data_mbate'

          print "(a)",'single precision dump'
          read(15,end=55,iostat=ierr) udisti, umassi, utimei, &
             nprint, n1, n2, timesi, gammasi, rhozeros, RK2s, &
             (dattemps(i,7), i=1, nprint),escaps, tkins, tgravs, tterms, &
             (dattemps(i,1), i=1, nprint), (dattemps(i,2), i=1, nprint), &
             (dattemps(i,3), i=1, nprint), (dattemps(i,4), i=1, nprint), &
             (dattemps(i,5), i=1, nprint), (dattemps(i,6), i=1, nprint), &
             (dattemps(i,8), i=1, nprint), (dattemps(i,9), i=1, nprint), &
             (dattemps(i,10), i=1, nprint), (dattemps(i,11),i=1,nprint), &
             dtmaxs, (isteps(i), i=1,nprint), (iphase(i),i=1,nprint), &
             nptmass, (listpm(i), i=1,nptmass)
       endif
!
!--set transformation factors between code units/real units
!
       call get_nearest_time_unit(utimei,unitx,unitslabel(0))
       units(0) = unitx ! convert to real*4
       call get_nearest_length_unit(udisti,unitx,unitlabelx)
       units(1:3) = unitx
       unitslabel(1:3) = unitlabelx
       units(4:6) = udisti/utimei
       unitslabel(4:6) = ' [cm/s]'
       units(7) = unitx
       unitslabel(7) = unitlabelx
       units(8) = (udisti/utimei)**2
       unitslabel(8) = ' [erg/g]'
       units(9) = umassi
       unitslabel(9) = ' [g]'
       units(10) = umassi/udisti**3
       unitslabel(10) = ' [g/cm^3]'
!
!--convert to single precision and separate pt masses from normal particles
!
       ipart = 0
       do i=1,nprint
          if (iphase(i)==0) then
             ipart = ipart + 1
             if (doubleprec) then
                dat(ipart,1:ncolstep,j) = real(dattemp(i,1:ncolstep))
             else
                dat(ipart,1:ncolstep,j) = dattemps(i,1:ncolstep)
             endif
          endif
       enddo
       npart = ipart
!
!--place point masses after normal particles
!
       nptmassi = 0
       do i=1,nprint
          if (iphase(i) >= 1) then
             ipart = ipart + 1
             nptmassi = nptmassi + 1
             if (doubleprec) then
                dat(ipart,1:ncolstep,j) = real(dattemp(i,1:ncolstep))
             else
                dat(ipart,1:ncolstep,j) = dattemps(i,1:ncolstep)
             endif
          endif
       enddo
       if (nptmass > 0) print*,' Number of point masses = ',nptmass
       if (nptmassi /= nptmass) print *,'WARNING: nptmass from iphase =',nptmassi,'not equal to nptmass'
!
!--put any others as unknown
!
       nunknown = 0
       do i=1,nprint
          if (iphase(i) < 0) then
             ipart = ipart + 1
             nunknown = nunknown + 1
             if (doubleprec) then
                dat(ipart,1:ncolstep,j) = real(dattemp(i,1:ncolstep))
             else
                dat(ipart,1:ncolstep,j) = dattemps(i,1:ncolstep)
             endif
          endif
       enddo
       if (nunknown > 0) print *,nunknown,' particles of unknown type (probably dead)'

       if (allocated(dattemp)) deallocate(dattemp)
       if (allocated(dattemps)) deallocate(dattemps)
       if (allocated(isteps)) deallocate(isteps)
       if (allocated(iphase)) deallocate(iphase)

       npartoftype(1,j) = npart
       npartoftype(2,j) = nptmassi
       npartoftype(3,j) = nunknown

       headertags(1:14) = (/'udist','umass','utime','npart',&
                           'n1   ','n2   ','time ','gamma',&
                           'rho0 ','RK2  ','escap','tkin ',&
                           'tgrav','tterm'/)
       if (doubleprec) then
          gamma(j) = real(gammai)
          time(j) = real(timei)
          headervals(1:14,j) = (/real(udisti),real(umassi),real(utimei),real(nprint),&
                                real(n1),real(n2),real(timei),real(gammai),&
                                real(rhozero),real(RK2),real(escap),real(tkin), &
                                real(tgrav),real(tterm)/)
       else
          gamma(j) = gammasi
          time(j) = timesi
          headervals(1:14,j) = (/real(udisti),real(umassi),real(utimei),real(nprint),&
                                real(n1),real(n2),real(timesi),real(gammasi),&
                                real(rhozeros),real(RK2s),real(escaps),real(tkins), &
                                real(tgravs),real(tterms)/)
       endif
       print*,' time = ',time(j),' gamma = ',gamma(j)

       if (ierr /= 0) then
          print "(a)",'*** INCOMPLETE DATA ***'
          nstepsread = nstepsread + 1
          exit over_steps_in_file
       else
          nstepsread = nstepsread + 1
       endif
       j = j + 1

    enddo over_steps_in_file

 endif

55 continue
!
!--reached end of file
!
 close(15)
 if (j-1  >  0) then
    print*,'>> end of dump file: nsteps =',j-1,'ntot = ', &
         sum(npartoftype(:,j-1)),'nghost=',npartoftype(2,j-1)
 endif

end subroutine read_data_mbate

!------------------------------------------------------------
! set labels for each column of data
!------------------------------------------------------------
subroutine set_labels_mbate
 use labels
 use params
 use settings_data
 use geometry, only:labelcoord
 integer :: i

 if (ndim <= 0 .or. ndim > 3) then
    print*,'*** ERROR: ndim = ',ndim,' in set_labels_mbate ***'
    return
 endif
 if (ndimV <= 0 .or. ndimV > 3) then
    print*,'*** ERROR: ndimV = ',ndimV,' in set_labels_mbate ***'
    return
 endif

 do i=1,ndim
    ix(i) = i
 enddo
 ivx = 4
 ih = 7        !  smoothing length
 iutherm = 8  !  thermal energy
 ipmass = 9   !  particle mass
 irho = 10     ! location of rho in data array
 if (ncolumns > 10) then
    label(11) = 'dgrav'
 endif

 label(ix(1:ndim)) = labelcoord(1:ndim,1)
 do i=1,ndimV
    label(ivx+i-1) = 'v\d'//labelcoord(i,1)
 enddo
 label(irho) = 'density'
 label(iutherm) = 'u'
 label(ih) = 'h'
 label(ipmass) = 'particle mass'

 !
 !--set labels for vector quantities
 !
 iamvec(ivx:ivx+ndimV-1) = ivx
 labelvec(ivx:ivx+ndimV-1) = 'v'
 do i=1,ndimV
    label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
 enddo
 !
 !--set labels for each particle type
 !
 ntypes = 3  !!maxparttypes
 labeltype(1) = 'gas'
 labeltype(2) = 'sink'
 labeltype(3) = 'unknown'
 UseTypeInRenderings(1) = .true.
 UseTypeInRenderings(2) = .false.
 UseTypeInRenderings(3) = .true.

end subroutine set_labels_mbate

end module readdata_mbate
