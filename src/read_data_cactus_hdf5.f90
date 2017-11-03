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
!  Copyright (C) 2005-2017 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR HDF5 OUTPUT FROM THE CACTUS CODE
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
!
!  The module below contains interface routines to c functions
!  that perform the actual calls to the HDF5 libs
!
!-------------------------------------------------------------------------
module cactushdf5read
 use params, only:maxplot,doub_prec
 use labels, only:lenlabel
 use, intrinsic :: iso_c_binding, only:c_int,c_double,c_char
 implicit none
 character(len=lenlabel), dimension(maxplot) :: blocklabel
 character(len=130) :: datfileprev = ' '
 logical :: file_is_open = .false.
 integer :: ntoti_prev,ncol_prev,nstep_prev

 interface
   subroutine open_cactus_hdf5_file(filename,istep,npart,ncol,nstep_max,ndim,ndimV,time,ignoretl,ierr) bind(c)
    import
    character(kind=c_char), dimension(*), intent(in) :: filename
    integer(kind=c_int), intent(in), value :: istep,ignoretl
    integer(kind=c_int), intent(out) :: npart,ncol,nstep_max,ndim,ndimV,ierr
    real(kind=c_double), intent(out) :: time
   end subroutine open_cactus_hdf5_file

   subroutine read_cactus_hdf5_data(filename,istep,npart,time,dx,ignoretl,ierr) bind(c)
    import
    character(kind=c_char), dimension(*), intent(in)  :: filename
    integer(kind=c_int), intent(in), value :: istep,ignoretl
    integer(kind=c_int), intent(out) :: npart,ierr
    real(kind=c_double), intent(out) :: time,dx
   end subroutine read_cactus_hdf5_data

   subroutine close_cactus_hdf5_file(ierr) bind(c)
    import
    integer(kind=c_int), intent(out) :: ierr
   end subroutine close_cactus_hdf5_file

 end interface

contains

 subroutine calc_trK(gxxd,gxyd,gxzd,gyyd,gyzd,gzzd,kxxd,kxyd,kxzd,kyyd,kyzd,kzzd,trk)
   !
   ! Subroutine to calculate trace K ( trK = g^{ij} K_{ij} )
   ! Takes g_{ij} and K_{ij} at one position and time, returns trK
   !
   real, intent(in) :: gxxd,gxyd,gxzd,gyyd,gyzd,gzzd ! spatial down metric components
   real, intent(in) :: kxxd,kxyd,kxzd,kyyd,kyzd,kzzd ! down extrinsic curvature components
   real :: gxxu,gxyu,gxzu,gyyu,gyzu,gzzu             ! spatial up metric components
   real, intent(out) :: trk
   real, dimension(3,3) :: gijd, giju  ! 4x4 metric down and up respectively
   real :: det

   ! down components
   gijd(1,1) = gxxd
   gijd(1,2) = gxyd
   gijd(1,3) = gxzd
   gijd(2,1) = gxyd
   gijd(2,2) = gyyd
   gijd(2,3) = gyzd
   gijd(3,1) = gxzd
   gijd(3,2) = gyzd
   gijd(3,3) = gzzd

   call inv3x3(gijd,giju,det)

   ! up (inverse) components
   gxxu = giju(1,1)
   gxyu = giju(1,2)
   gxzu = giju(1,3)
   gyyu = giju(2,2)
   gyzu = giju(2,3)
   gzzu = giju(3,3)


   trk = (gxxu * kxxd) + (2. * gxyu * kxyd) + (2. * gxzu * kxzd) + &
        & (gyyu * kyyd) + (2. * gyzu * kyzd) + (gzzu * kzzd)

 end subroutine calc_trK

 pure subroutine inv3x3(A,B,det)
   real, intent(in), dimension(3,3) :: A
   real, intent(out), dimension(3,3) :: B ! inverse matrix
   real, intent(out) :: det

   det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - &
        & A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + &
        & A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

   B(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
   B(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
   B(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
   B(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
   B(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
   B(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
   B(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
   B(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
   B(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

   B(:,:) = B(:,:)/det

 end subroutine inv3x3

 subroutine compute_extra_columns(ncols,nextra,dat)
  integer, intent(in)  :: ncols
  integer, intent(out) :: nextra
  real, intent(inout), optional  :: dat(:,:)
  integer :: i,n,itrk
  integer :: igxx,igxy,igxz,igyy,igyz,igzz
  integer :: ikxx,ikxy,ikxz,ikyy,ikyz,ikzz

  nextra = 0
  !
  ! find gxx,gxy,gxz and kxx,kxy,kxz etc in columns
  !
  n = size(dat(:,1))
  do i=1,ncols
     select case(blocklabel(i))
     case('gxx')
        igxx = i
     case('gxy')
        igxy = i
     case('gxz')
        igxz = i
     case('gyy')
        igyy = i
     case('gyz')
        igyz = i
     case('gzz')
        igzz = i
     case('kxx')
        ikxx = i
     case('kxy')
        ikxy = i
     case('kxz')
        ikxz = i
     case('kyy')
        ikyy = i
     case('kyz')
        ikyz = i
     case('kzz')
        ikzz = i
     end select
  enddo
  if (igxx > 0 .and. igxy > 0 .and. igxz > 0 .and. igyy > 0 .and. igyz > 0 .and. igzz > 0 .and. &
      ikxx > 0 .and. ikxy > 0 .and. ikxz > 0 .and. ikyy > 0 .and. ikyz > 0 .and. ikzz > 0) then
     itrk = ncols + 1
     blocklabel(itrk) = 'tr K'
     nextra = 1
     if (present(dat)) then
        !print*,' getting trk ',igxx,igxy,igxz,igyy,igyz,igzz,ikxx,ikxy,ikxz,ikyy,ikyz,ikzz,itrk
        do i=1,n
           call calc_trK(dat(i,igxx),dat(i,igxy),dat(i,igxz),dat(i,igyy),dat(i,igyz),dat(i,igzz),&
                         dat(i,ikxx),dat(i,ikxy),dat(i,ikxy),dat(i,ikyy),dat(i,ikyz),dat(i,ikzz),dat(i,itrk))
        enddo
     endif
  endif

 end subroutine compute_extra_columns

end module cactushdf5read

!-------------------------------------------------------------------------
!
!  The routine that reads the data into splash's internal arrays
!
!-------------------------------------------------------------------------
subroutine read_data(rootname,istepstart,ipos,nstepsread)
  use particle_data,  only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol,maxstep,iamtype
  use params,         only:doub_prec
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc,ipartialread,iverbose,buffer_steps_in_file
  use settings_page,  only:legendtext
  use mem_allocation, only:alloc
  use labels,         only:ih,irho,ipmass
  use system_utils,   only:renvironment,lenvironment,ienvironment,envlist
  use asciiutils,     only:cstring
  use cactushdf5read, only:open_cactus_hdf5_file,read_cactus_hdf5_data,close_cactus_hdf5_file,&
                           datfileprev,file_is_open,ntoti_prev,ncol_prev,nstep_prev,compute_extra_columns
  use dataread_utils, only:count_types
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

  if (len_trim(rootname).gt.0) then
     datfile = trim(rootname)
  else
     print*,' **** no data read **** '
     return
  endif
!
!--check if first data file exists
!
  if (iverbose==1 .and. ipos==1) print "(1x,a)",'reading CACTUS HDF5 format'
  inquire(file=datfile,exist=iexist)
  if (.not.iexist) then
     !
     !--append .h5 on the end if not already present
     !
     datfile=trim(rootname)//'.h5'
     inquire(file=datfile,exist=iexist)
  endif
!
!--close previous file if filenames do not match
!
  if (trim(datfile)/=trim(datfileprev) .and. file_is_open) then
     call close_cactus_hdf5_file(ierr)
     file_is_open = .false.
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
  ignore_time_levels = lenvironment('CSPLASH_IGNORE_TIME_LEVELS')
  ignoretl = 0
  if (ignore_time_levels) ignoretl = 1
  nextra = 0
! 
!--read data from snapshots
!
  i = istepstart
  if (.not.file_is_open) write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  !--open file and read header information
  !
  if (file_is_open) then
     ntoti     = ntoti_prev
     ncolstep  = ncol_prev
     nstep_max = nstep_prev
  else
     call open_cactus_hdf5_file(cstring(datfile),ipos,ntoti,ncolstep,nstep_max,ndim,ndimV,timetemp,ignoretl,ierr)
     if (ierr /= 0) then
        print "(a)", '*** ERROR READING HEADER ***'
        call close_cactus_hdf5_file(ierr)
        return
     endif
     file_is_open = .true.
     datfileprev = datfile
     ntoti_prev = ntoti
     ncol_prev  = ncolstep
     nstep_prev = nstep_max
  endif
  call compute_extra_columns(ncolstep,nextra)
  ncolumns = ncolstep + nextra
  if (iverbose >= 1) print "(3(a,1x,i10))",' npart: ',ntoti,' ncolumns: ',ncolstep,' nsteps: ',nstep_max

  istep = 1 
  over_snapshots: do istep=1,nstep_max
  !
  !--now read data
  !
  reallocate = .false.
  npart_max = maxpart
  if (buffer_steps_in_file) then
     nsteps_to_read = nstep_max
  else
     nsteps_to_read = max(maxstep,1)
  endif
  if (nsteps_to_read > maxstep) reallocate = .true.

  if (ntoti.gt.maxpart) then
     reallocate = .true.
     if (maxpart.gt.0) then
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
     call alloc(npart_max,nsteps_to_read,max(ncolumns+ncalc,maxcol),mixedtypes=.true.)
  endif

  !
  !--read particle data
  !
  got_particles: if (ntoti > 0) then

     if (buffer_steps_in_file .or. ipos.eq.istep) then
        call read_cactus_hdf5_data(cstring(datfile),istep,ntoti,timetemp,dx,ignoretl,ierr)
        call set_labels
        ! set smoothing length and particle mass
        !print*,' Setting h = ',dx, 'ndim = ',ndim,' in column ',ih,' step ',i
        if (ih > 0) dat(:,ih,i) = real(dx)
        vol = dx**ndim
        if (ipmass > 0 .and. irho > 0) dat(:,ipmass,i) = dat(:,irho,i)*real(vol)
        !
        ! compute extra quantities (tr K, 3^R, etc)
        !
        call compute_extra_columns(ncolstep,nextra,dat(:,:,i))
        !
        ! get number of cells of each type (normal, ghost)
        !
        call count_types(ntoti,iamtype(:,i),npartoftype(:,i),nunknown)
        masstype(:,i) = 0. ! all masses read from file
        time(i) = real(timetemp)
        i = i + 1
     endif

     nstepsread = nstepsread + 1

  endif got_particles
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
  call set_labels
  
  enddo over_snapshots

  if (nstepsread.gt.0) then
     print "(a,i10,a)",' >> read ',sum(npartoftype(:,istepstart)),' cells'
  endif

  !if (ipos==nstep_max) then
  !   call close_cactus_hdf5_file(ierr)
  !   file_is_open = .false.
  !endif

end subroutine read_data

subroutine read_cactus_hdf5_data_fromc(icol,ntot,np,temparr) bind(c)
  use, intrinsic :: iso_c_binding, only:c_int,c_double
  use particle_data,  only:dat
  use settings_data,  only:debugmode
  use labels,         only:label
  implicit none
  integer(kind=c_int), intent(in) :: icol,ntot,np
  real(kind=c_double), intent(in) :: temparr(np)
  integer(kind=c_int) :: icolput
  integer :: nmax,i1,i2

  icolput = icol
  i1 = ntot-np+1
  i2 = ntot
  
  if (debugmode) print "(a,i2,a,i8,a,i8)",&
  'DEBUG: reading column ',icol,' -> '//trim(label(icolput))//' parts ',i1,' to ',i2
  
  ! check column is within array limits
  if (icolput.gt.size(dat(1,:,1)) .or. icolput.eq.0) then
     print "(a,i2,a)",' ERROR: column = ',icolput,' out of range in receive_data_fromc'
     return
  endif
  if (i2 > size(dat(:,1,1))) then
     print*,' ERROR with index range: ',i1,':',i2,' exceeds size ',size(dat(:,1,1)),' for column ',icol
     read*
     return
  endif
  ! ensure no array overflows
  nmax = min(i2,size(dat(:,1,1)))

  ! copy data into main splash array
  dat(i1:i2,icolput,1) = real(temparr(1:np))

  return
end subroutine read_cactus_hdf5_data_fromc

subroutine read_cactus_itype_fromc(ntot,np,itype) bind(c)
  use, intrinsic :: iso_c_binding, only:c_int
  use particle_data,  only:iamtype
  use params,         only:int1
  implicit none
  integer(kind=c_int), intent(in) :: ntot,np
  integer(kind=c_int), intent(in) :: itype(np)
  integer :: i1,i2,len_type

  i1 = ntot-np+1
  i2 = ntot
  ! set particle type
  len_type = size(iamtype(:,1))
  if (len_type.gt.1) then
     if (i2 > len_type) then
        print*,'error with itype length',i2,len_type
        return
     endif
     iamtype(i1:i2,1) = int(itype(1:np),kind=int1)
  endif

  return
end subroutine read_cactus_itype_fromc

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels,        only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass,iutherm,ih,irho
  use params
  use settings_data,  only:ndim,ndimV,ntypes,UseTypeInRenderings
  use geometry,       only:labelcoord
  use system_utils,   only:envlist,ienvironment
  use cactushdf5read, only:blocklabel
  use asciiutils,     only:lcase
  implicit none
  integer :: i,icol

  if (ndim.le.0 .or. ndim.gt.3) then
     print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
     return
  endif
  if (ndimV.le.0 .or. ndimV.gt.3) then
     print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
     return
  endif
  blocklabel(1:5) = (/'x ','y ','z ','dx','m '/)

  ix(1) = 1
  ix(2) = 2
  ix(3) = 3
  ih = 4
  ipmass = 5
  iutherm = 0
  irho = 0
  do icol=1,size(blocklabel)
     select case(trim(blocklabel(icol)))
     case('vel[0]')
        ivx = icol
     case('dens','density')
        if (irho==0) irho = icol
     case('rho')
        irho = icol
     end select
     label(icol) = trim(blocklabel(icol))
  enddo

  ! set labels of the quantities read in
  if (ix(1).gt.0)   label(ix(1:ndim)) = labelcoord(1:ndim,1)
  !if (irho.gt.0)    label(irho)       = 'density'
  !if (iutherm.gt.0) label(iutherm)    = 'u'
  !if (ipmass.gt.0)  label(ipmass)     = 'particle mass'
  !if (ih.gt.0)      label(ih)         = 'h'

  ! set labels for vector quantities
  if (ivx.gt.0) then
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
  

!-----------------------------------------------------------
  return
end subroutine set_labels

subroutine set_blocklabel(icol,name,lenname) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int, c_char
 use cactushdf5read, only:blocklabel
 use asciiutils,    only:fstring
 implicit none
 integer(kind=c_int),    intent(in) :: icol,lenname
 character(kind=c_char), intent(in) :: name(lenname)
 character(len=24) :: temp
 integer :: ivar
 
 temp = fstring(name)
 ivar = index(temp,'::')
 if (ivar > 0) temp = temp(ivar+2:)
 if (icol <= size(blocklabel)) then
    blocklabel(icol) = trim(temp)
 else
    print*,'ERROR - too many columns in file'
 endif
 !print*,icol,' name = ',trim(blocklabel(icol))

end subroutine set_blocklabel

subroutine sort_cactus_data(n,iter,iorder) bind(c)
 use, intrinsic :: iso_c_binding, only:c_int
 use sort, only:indexxi
 implicit none
 integer(kind=c_int), intent(in)  :: n
 integer(kind=c_int), intent(in)  :: iter(n)
 integer(kind=c_int), intent(out) :: iorder(n)

 call indexxi(n,iter,iorder)

end subroutine sort_cactus_data
