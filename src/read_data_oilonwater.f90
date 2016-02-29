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

!--short module storing units information
module oilonwaterread
 use params
 implicit none
 real(doub_prec) :: udisti,umassi,utimei

end module oilonwaterread


subroutine read_data(rootname,indexstart,ipos,nstepsread)
  use particle_data
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation
  use labels, only:ix,ivx,ih,irho,ipmass
  use oilonwaterread, only:udisti,umassi,utimei
  implicit none
  integer, intent(in) :: indexstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer, parameter :: maxptmass = 1000
  real, parameter :: pi=3.141592653589
  integer :: i,j,k,ifile,ierr,ipart
  integer :: npart_max,nstep_max,ncolstep,npart,nunknown
  logical :: iexist,doubleprec

  character(len=len(rootname)+10) :: dumpfile
  integer :: nprint, n1, n2, nptmass, nstepsalloc
  integer :: npartoil, npartwater
  integer, dimension(:), allocatable :: isteps, iphase

  !--use these lines if dump is double precision
  real(doub_prec), dimension(:,:), allocatable :: dattemp
  real(doub_prec) :: timei, gammai
  real(doub_prec) :: rhozero, RK2
  real(doub_prec) :: escap,tkin,tgrav,tterm
  real(doub_prec) :: dtmax

  !--use these lines for single precision
  real :: timesi, gammasi
  real :: rhozeros, RK2s
  real :: escaps,tkins,tgravs,tterms
  real :: dtmaxs

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  ifile = 1

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
  ncolstep = 18  ! number of columns in file
  ncolumns = ncolstep
  !
  !--allocate memory initially
  !
  nstep_max = max(nstep_max,indexstart,1)

  j = indexstart
  nstepsread = 0
  doubleprec = .false.

  print "(1x,a)",'reading oil-on-water code format'
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
      print*,'npart = ',nprint
      if (.not.allocated(dat) .or. nprint.gt.npart_max) then
         npart_max = max(npart_max,INT(1.1*nprint))
         call alloc(npart_max,nstep_max,ncolstep+ncalc,mixedtypes=.true.)
      endif
      doubleprec = .true.
      !--try single precision if non-sensible values for time, gamma etc.
      if (ierr.ne.0 .or. timei.lt.0. .or. timei.gt.1e30  &
          .or. gammai.lt.1. .or. gammai.gt.10. &
          .or. rhozero.lt.0. .or. RK2.lt.0. .or. RK2.gt.1.e10) then
          doubleprec = .false.
      endif
      rewind(15)
   endif

   call set_labels

   if (ierr /= 0) then
      print "(a)",'*** ERROR READING TIMESTEP HEADER ***'
   else
!
!--loop over the timesteps in this file
!
   over_steps_in_file: do
     npart_max = max(npart_max,nprint)
!
!--allocate/reallocate memory if j > maxstep
!
     if (j.gt.maxstep) then
        if (nstepsread.ge.2) then
           nstepsalloc = 2*nstepsread
        else
           nstepsalloc = j
        endif
        call alloc(maxpart,nstepsalloc,maxcol,mixedtypes=.true.)
     endif
!
!--allocate integer arrays required for data read
!
     if (allocated(isteps)) deallocate(isteps)
     allocate(isteps(npart_max),stat=ierr)
     if (ierr /= 0) print*,'not enough memory in read_data'

     if (allocated(iphase)) deallocate(iphase)
     allocate(iphase(npart_max),stat=ierr)
     if (ierr /= 0) print*,'not enough memory in read_data'
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
        if (ierr /= 0) print*,'not enough memory in read_data'

        read(15,end=55,iostat=ierr) udisti, umassi, utimei, &
             nprint, n1, n2, timei, gammai, rhozero, RK2, &
             (dattemp(i,ih), i=1, nprint),escap, tkin, tgrav, tterm, &
             ((dattemp(i,k), i=1, nprint),k=1,ih-1), &
             ((dattemp(i,k), i=1, nprint),k=ih+1,11), &
             dtmax, &
              (isteps(i),    i=1,nprint), &
             ((dattemp(i,k), i=1, nprint),k=12,18), &
             nptmass, &
             ((dattemp(i,k),     i=nprint+1,nprint+nptmass),k=1,6), &
             (dattemp(i,ipmass), i=nprint+1,nprint+nptmass), &
             (dattemp(i,ih),     i=nprint+1,nprint+nptmass), &
             (isteps(i),         i=nprint+1,nprint+nptmass), &
             (iphase(i),         i=1,nprint)
     else
        print "(a)",'single precision dump'
        read(15,end=55,iostat=ierr) udisti, umassi, utimei, &
             nprint, n1, n2, timesi, gammasi, rhozeros, RK2s, &
             (dat(i,ih,j), i=1, nprint),escaps, tkins, tgravs, tterms, &
             ((dat(i,k,j), i=1, nprint),k=1,ih-1), &
             ((dat(i,k,j), i=1, nprint),k=ih+1,11), &
             dtmaxs, &
              (isteps(i),i=1, nprint), &
             ((dat(i,k,j), i=1, nprint),k=12,18), &
             nptmass, &
             ((dat(i,k,j),     i=nprint+1,nprint+nptmass),k=1,6), &
             (dat(i,ipmass,j), i=nprint+1,nprint+nptmass), &
             (dat(i,ih,j),     i=nprint+1,nprint+nptmass), &
             (isteps(i),          i=nprint+1,nprint+nptmass), &
             (iphase(i),          i=1,nprint)
     endif
!
!--extract time and gamma
!
     if (doubleprec) then
        gamma(j) = real(gammai)
        time(j) = real(timei)
     else
        gamma(j) = gammasi
        time(j) = timesi
     endif
     print*,' time = ',time(j),' gamma = ',gamma(j)

     npart = nprint+nptmass
!
!--convert to single precision if necessary
!
     if (doubleprec) then
        dat(1:npart,1:ncolstep,j) = real(dattemp(1:npart,1:ncolstep))
     endif
!
!--set particle types using iphase
!
     ipart = 0
     npartoil = 0
     npartwater = 0
     nunknown = 0
     do i=1,nprint
        select case(iphase(i))
        case(0)
           npartoil = npartoil + 1
           iamtype(i,j) = 1
        case(1)
           npartwater = npartwater + 1
           iamtype(i,j) = 2
        case default
           nunknown = nunknown + 1
           iamtype(i,j) = 4
        end select
     enddo
!
!--set particle type for point masses
!
     if (nptmass.gt.0) then
        iamtype(nprint+1:nprint+nptmass,j) = 3
     endif

     if (nunknown.gt.0) print *,nunknown,' particles of unknown type (probably dead)'

     if (allocated(dattemp)) deallocate(dattemp)
     if (allocated(isteps)) deallocate(isteps)
     if (allocated(iphase)) deallocate(iphase)

     npartoftype(1,j) = npartoil
     npartoftype(2,j) = npartwater
     npartoftype(3,j) = nptmass
     npartoftype(4,j) = nunknown
     print "(1x,3(a,i10))",' n(oil) = ',npartoil,' n(water) = ',npartwater,' nptmass = ',nptmass

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
if (j-1 .gt. 0) then
   print*,'>> end of dump file: nsteps =',j-1
endif

return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels
  use params
  use physcon, only:solarrcgs,solarmcgs
  use settings_data
  use geometry, only:labelcoord
  use oilonwaterread, only:udisti,umassi,utimei
  use settings_units, only:units,unitzintegration
  implicit none
  integer :: i

  if (ndim.le.0 .or. ndim.gt.3) then
     print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
     return
  endif
  if (ndimV.le.0 .or. ndimV.gt.3) then
     print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
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
  if (ncolumns.gt.10) then
     label(11) = 'dgrav'
     label(12) = 'torque t'
     label(13) = 'torque g'
     label(14) = 'torque p'
     label(15) = 'torque v'
     label(16) = 'torque c'
     label(17) = 'hecomp'
     label(18) = 'potential energy'
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
  !--set transformation factors between code units/real units
  !
  units(1:3) = udisti/solarrcgs
  unitslabel(1:3) = ' [R\d\(2281)\u]'
  units(4:6) = udisti/utimei
  unitslabel(4:6) = ' [cm/s]'
  units(ih) = units(1)
  unitslabel(ih) = unitslabel(1)
  units(iutherm) = (udisti/utimei)**2
  unitslabel(iutherm) = ' [erg/g]'
  units(ipmass) = umassi/solarmcgs
  unitslabel(ipmass) = ' [M\d\(2281)\u]'
  units(irho) = umassi/udisti**3
  unitslabel(irho) = ' [g/cm\u3\d]'

  !--unit for z integration - leave this as cm to get g/cm^2 in column density
  unitzintegration = udisti
  labelzintegration = ' [cm]'

  !--time units for legend
  units(0) = utimei/3.1536e7
  unitslabel(0) = ' yrs'

  !
  !--set labels for each particle type
  !
  ntypes = 4
  labeltype(1) = 'gas (oil)'
  labeltype(2) = 'gas (water)'
  labeltype(3) = 'point mass'
  labeltype(4) = 'unknown'
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2) = .true.
  UseTypeInRenderings(3) = .false.
  UseTypeInRenderings(4) = .true.

!-----------------------------------------------------------

  return
end subroutine set_labels
