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
!  Copyright (C) 2005-2016 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR OUTPUT FROM THE GADGET CODE
! AS MODIFIED BY JAMIE BOLTON
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
!-------------------------------------------------------------------------

subroutine read_data(rootname,istart,ipos,nstepsread)
  use particle_data
  use params
  use labels
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation
  implicit none
  integer, intent(IN) :: istart,ipos
  integer, intent(OUT) :: nstepsread
  character(LEN=*), intent(IN) :: rootname
  character(LEN=LEN(rootname)+10) :: datfile
  integer, dimension(maxparttypes) :: npartoftypei
  integer, dimension(:), allocatable :: iamtemp
  integer :: i,itype,icol,ifile,idashpos,ierr
  integer :: index1,index2,indexstart,indexend,Nmassesdumped
  integer :: ncol_max,npart_max,nstep_max,ntoti
  logical :: iexist,reallocate
  real(doub_prec) :: timetemp
  real(doub_prec), dimension(6) :: Massoftype
  real, dimension(:), allocatable :: dattemp1
  real, dimension(:,:), allocatable :: dattemp

  nstepsread = 0

  if (len_trim(rootname).gt.0) then
     datfile = trim(rootname)
  else
     print*,' **** no data read **** '
     return
  endif
!
!--check if first data file exists
!
  inquire(file=datfile,exist=iexist)
  if (.not.iexist) then
     print "(a)",' *** error: ',trim(datfile),' file not found ***'
     return
  endif
!
!--set parameters which do not vary between timesteps
!
  ndim = 3
  ndimV = 3
  ncol_max = 15 ! 3 x pos, 3 x vel, utherm, rho, Ne, h, pmass
!
!--read data from snapshots
!
  i = istart

  print "(1x,a)",'reading Jamie Bolton''s modified GADGET format'
  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  !--open data file and read data
  !
  open(11,ERR=81,file=datfile,status='old',form='unformatted')
  !
  !--read header for this timestep
  !
  read(11,ERR=70,end=80) npartoftypei,Massoftype,timetemp
  ntoti = int(sum(npartoftypei))
  print*,'time             : ',timetemp
  print*,'Npart (by type)  : ',npartoftypei
  print*,'Mass  (by type)  : ',Massoftype
  print*,'N_gas            : ',npartoftypei(1)
  print*,'N_total          : ',ntoti

  !
  !--if successfully read header, increment the nstepsread counter
  !
  nstepsread = nstepsread + 1
  !
  !--now read data
  !
  reallocate = .false.
  npart_max = maxpart
  nstep_max = max(maxstep,1)

  if (ntoti.gt.maxpart) then
     reallocate = .true.
     npart_max = int(1.1*ntoti)
  endif
  if (i.ge.maxstep .and. i.ne.1) then
     nstep_max = i + max(10,INT(0.1*nstep_max))
     reallocate = .true.
  endif
  !
  !--reallocate memory for main data array
  !
  if (reallocate .or. .not.(allocated(dat))) then
     call alloc(npart_max,nstep_max,max(ncol_max+ncalc,maxcol))
  endif
  !
  !--copy header into header arrays
  !
  npartoftype(:,i) = npartoftypei
  time(i) = real(timetemp)


  if (ntoti.gt.0) then
     if (allocated(dattemp)) deallocate(dattemp)
     allocate(dattemp(3,ntoti))
     !
     !--read positions of all particles
     !
     print*,'positions ',ntoti
     read (11, iostat=ierr) dattemp(1:3,1:ntoti)
     if (ierr /= 0) then
        print "(a)",'error encountered whilst reading positions '
        return
     else
        do icol=1,3
           dat(1:ntoti,icol,i) = dattemp(icol,1:ntoti)
        enddo
     endif
     !
     !--same for velocities
     !
     print*,'velocities ',ntoti
     read (11, iostat=ierr) dattemp(1:3,1:ntoti)
     if (ierr /= 0) then
        print "(a)",'error encountered whilst reading velocities'
     else
        do icol=4,6
           dat(1:ntoti,icol,i) = dattemp(icol-3,1:ntoti)
        enddo
     endif
     !
     !--read particle ID
     !
     print*,'particle ID ',ntoti
     if (allocated(iamtemp)) deallocate(iamtemp)
     allocate(iamtemp(npart_max))
     read (11, end=66,ERR=73) iamtemp(1:ntoti)
     deallocate(iamtemp)
     !
     !--read particle masses
     !
     !--work out total number of masses dumped
     Nmassesdumped = 0
     do itype = 1,6
        if (abs(Massoftype(itype)).lt.1.e-8) then
           Nmassesdumped = Nmassesdumped + Npartoftype(itype,i)
        endif
     enddo
     print*,'particle masses ',Nmassesdumped

     !--read this number of entries
     if (allocated(dattemp1)) deallocate(dattemp1)
     allocate(dattemp1(Nmassesdumped))
     if (Nmassesdumped.gt.0) then
        read(11,end=66,err=74) dattemp1(1:Nmassesdumped)
     endif
     !--now copy to the appropriate sections of the .dat array
     indexstart = 1
     index1 = 1

     do itype=1,6
        if (Npartoftype(itype,i).ne.0) then
           index2 = index1 + Npartoftype(itype,i) -1
           if (abs(Massoftype(itype)).lt.1.e-8) then ! masses dumped
              indexend = indexstart + Npartoftype(itype,i) - 1
              print*,'read ',Npartoftype(itype,i),' masses for type ', &
                     itype,index1,'->',index2,indexstart,'->',indexend
              dat(index1:index2,7,i) = dattemp1(indexstart:indexend)
           else  ! masses not dumped
              print*,'setting masses for type ',itype,' = ', &
                     real(Massoftype(itype)),index1,'->',index2
              dat(index1:index2,7,i) = real(Massoftype(itype))
           endif
           index1 = index2 + 1
           indexstart = indexend + 1
        endif
     enddo
     deallocate(dattemp1)
     !
     !--read other quantities for rest of particles
     !
     print*,'gas properties ',npartoftype(1,i)
     do icol=8,15
        !!print*,icol
        read (11, end=66,ERR=78) dat(1:npartoftype(1,i),icol,i)
        !
        !--for some reason the smoothing length output by GADGET is
        !  twice the usual SPH smoothing length
        !
        if (icol.eq.15) then
           dat(1:npartoftype(1,i),icol,i) = 0.5*dat(1:npartoftype(1,i),icol,i)
        endif
     enddo


  else
     ntoti = 1
     npartoftype(1,i) = 1
     dat(:,:,i) = 0.
  endif

  !!ntot(i-1) = j-1
!
!--now memory has been allocated, set arrays which are constant for all time
!
  gamma = 5./3.
  goto 68

66 continue
  print*,'*** end of file reached in ',trim(datfile),' ***'
  ! timestep there but data incomplete
  goto 68

68 continue
  !
  !--close data file and return
  !
  close(unit=11)

  ncolumns = ncol_max
  print*,'ncolumns = ',ncolumns

  print*,'>> Finished reading: steps =',nstepsread-istart+1, &
         'last step ntot =',sum(npartoftype(:,istart+nstepsread-1))
  return
!
!--errors
!
70 continue
  print*,' *** Error encountered while reading timestep header ***'
  print*,' Npartoftype = ',Npartoftype(:,i)
  print*,' Massoftype = ',Massoftype
  return

73 continue
  print*,' *** Error encountered while reading particle ID ***'
  return

74 continue
  print*,' *** Error encountered while reading particle masses ***'
  return

78 continue
  print*,' *** Error encountered while reading gas particle properties ***'
  return

80 continue
  print*,' *** data file empty, no steps read ***'
  return

81 continue
  print*,' *** Error: can''t open data file ***'
  return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels
  use params
  use settings_data
  use geometry, only:labelcoord
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
  ipmass = 7
  irho = 9        ! location of rho in data array
  ipr = 0
  iutherm = 8     !  thermal energy
  ih = 15         !  smoothing length
  !
  !--set labels of the quantities read in
  !
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(irho) = '\gr'
  label(iutherm) = 'u'
  label(10) = 'NHp'
  label(11) = 'NHep'
  label(12) = 'NHepp'
  label(13) = 'NH0'
  label(14) = 'NHe0'
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

  !--set labels for each particle type
  !
  ntypes = 6
  labeltype(1) = 'gas'
  labeltype(2) = 'dark matter'
  labeltype(3) = 'boundary 1'
  labeltype(4) = 'boundary 2'
  labeltype(5) = 'star'
  labeltype(6) = 'sink / black hole'
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2:6) = .false.

!-----------------------------------------------------------
  return
end subroutine set_labels
