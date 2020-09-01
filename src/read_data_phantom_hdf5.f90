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
module readphantomhdf5

  implicit none

  integer, parameter :: itypemap_sink_phantom    = 3
  integer, parameter :: itypemap_unknown_phantom = 9

contains

!---------------------------------------------------------------------
! function mapping iphase setting in Phantom to splash particle types
!---------------------------------------------------------------------
elemental integer function itypemap_phantom(iphase)
  integer(kind=1), intent(in) :: iphase

  select case(int(iphase))
  case(1:2)
    itypemap_phantom = iphase
  case(3:7) ! put sinks as type 3, everything else shifted by one
    itypemap_phantom = iphase + 1
  case(-3) ! sink particles, either from external_binary or read from dump
    itypemap_phantom = itypemap_sink_phantom
  case default
    itypemap_phantom = itypemap_unknown_phantom
  end select

end function itypemap_phantom

subroutine get_ncolumns(file_id,ncolumns,ndusttypes,ndustsmall,nptmass)
  use hdf5_utils, only:HID_T,exist_in_hdf5,open_hdf5group,close_hdf5group
  integer(HID_T), intent(in)  :: file_id
  integer,        intent(out) :: ncolumns
  integer,        intent(in)  :: ndusttypes,ndustsmall,nptmass
  integer(HID_T) :: group_id
  integer        :: error
  logical        :: exist

  ncolumns = 0

  call open_hdf5group(file_id,'particles',group_id,error)

  !-- Count rank 1 arrays
  call exist_in_hdf5(group_id,'h',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'u',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'alpha',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'dt',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'poten',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'psi',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'divB',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'divBsymm',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'eta_OR',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'eta_HE',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'eta_AD',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'ne_on_n',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'grainsize',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'graindens',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'vrel_on_vfrag',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'st',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'temperature',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'divv',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'luminosity',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  call exist_in_hdf5(group_id,'beta',exist,error)
  if (exist) then
     ncolumns = ncolumns + 1
  endif

  !-- Count rank 2 arrays (3 x npart)
  call exist_in_hdf5(group_id,'xyz',exist,error)
  if (exist) then
     ncolumns = ncolumns + 3
  endif

  call exist_in_hdf5(group_id,'vxyz',exist,error)
  if (exist) then
     ncolumns = ncolumns + 3
  endif

  call exist_in_hdf5(group_id,'Bxyz',exist,error)
  if (exist) then
     ncolumns = ncolumns + 3
  endif

  call exist_in_hdf5(group_id,'curlBxyz',exist,error)
  if (exist) then
     ncolumns = ncolumns + 3
  endif

  call exist_in_hdf5(group_id,'curlvxyz',exist,error)
  if (exist) then
     ncolumns = ncolumns + 3
  endif

  !-- Count rank 2 arrays (ndusttypes x npart)
  call exist_in_hdf5(group_id,'dustfrac',exist,error)
  if (exist) then
    ncolumns = ncolumns + ndusttypes
  endif
  call exist_in_hdf5(group_id,'tstop',exist,error)
  if (exist) then
    ncolumns = ncolumns + ndusttypes
  endif

  !-- Count abundance array (5 x npart)
  call exist_in_hdf5(group_id,'abundance',exist,error)
  if (exist) then
     ncolumns = ncolumns + 5
  endif

  !--- Count dust deltavxyz array (3 x ndustsmall x npart)
  call exist_in_hdf5(group_id,'deltavxyz',exist,error)
  if (exist) then
     ncolumns = ncolumns + (3*ndustsmall)
  endif

  call close_hdf5group(group_id,error)

  !-- Count sink arrays --------------------
  call open_hdf5group(file_id,'sinks',group_id,error)

  ! Rank 1 arrays (nptmass)
  call exist_in_hdf5(group_id,'m',exist,error)
  call exist_in_hdf5(group_id,'h',exist,error)
  call exist_in_hdf5(group_id,'hsoft',exist,error)
  call exist_in_hdf5(group_id,'maccreted',exist,error)
  call exist_in_hdf5(group_id,'tlast',exist,error)

  ! Rank 2 (3 x nptmass arrays)
  call exist_in_hdf5(group_id,'xyz',exist,error)
  call exist_in_hdf5(group_id,'spinxyz',exist,error)
  call exist_in_hdf5(group_id,'vxyz',exist,error)

  call close_hdf5group(group_id,error)

end subroutine get_ncolumns

end module readphantomhdf5

!-------------------------------------------------------------------------
!
!  The routine that reads the data into splash's internal arrays
!
!-------------------------------------------------------------------------

module readdata_phantom_hdf5
 implicit none

 public :: read_data_phantom_hdf5, set_labels_phantom_hdf5

 private
contains

subroutine read_data_phantom_hdf5(dumpfile,ifile,ipos,nstepsread)
  use particle_data,   only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol,iamtype,maxstep
  use settings_data,   only:ndim,ndimV,ncolumns,ncalc,ntypes,ndusttypes,maxparttypes
  use labels,          only:label,labeltype,ipmass,irho,ih,ix,ivx,idivB,&
                            iBfirst,idustfrac,igraindens,igrainsize,iutherm
  use hdf5_utils,      only:open_hdf5file,open_hdf5group,close_hdf5file,&
                            close_hdf5group,read_from_hdf5,HID_T
  use readphantomhdf5, only:get_ncolumns,itypemap_phantom
  use mem_allocation,  only:alloc
  implicit none
  integer,          intent(in)  :: ifile,ipos
  integer,          intent(out) :: nstepsread
  character(len=*), intent(in)  :: dumpfile
  real,             allocatable  :: array_1d(:),array_2d(:,:),array_3d(:,:,:)
  integer(kind=1),  allocatable  :: itype(:)
  character(len=7), dimension(5) :: abundance_labels
  integer(HID_T) :: file_id,group_id
  integer        :: i,icol,nptmass,ndustsmall,npart,ndustlarge,error,ncolstep,npartoftypei(100)
  logical        :: got
  real           :: hfact,massoftypei(100),timei,gammai

  !-- Dimensions
  ndim  = 3
  ndimV = 3

  !-- Label particle type with strings
  labeltype(1) = 'gas'
  labeltype(2) = 'dust (old)'
  labeltype(3) = 'sink'
  labeltype(4) = 'ghost'
  labeltype(5) = 'star'
  labeltype(6) = 'dark matter'
  labeltype(7) = 'bulge'
  labeltype(8) = 'dust'
  labeltype(9) = 'unknown/dead'

  print "(a)",' Reading PHANTOM HDF5 data format '
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)

  call open_hdf5file(dumpfile,file_id,error,read_only=.true.)
  nstepsread = 1

  !-- Read values from header
  call open_hdf5group(file_id,'header',group_id,error)
  call read_from_hdf5(nptmass,'nptmass',group_id,got,error)
  call read_from_hdf5(ndustsmall,'ndustsmall',group_id,got,error)
  call read_from_hdf5(ndustlarge,'ndustlarge',group_id,got,error)
  call read_from_hdf5(npart,'nparttot',group_id,got,error)
  call read_from_hdf5(ntypes,'ntypes',group_id,got,error)
  call read_from_hdf5(hfact,'hfact',group_id,got,error)
  call read_from_hdf5(gammai,'gamma',group_id,got,error)
  call read_from_hdf5(timei,'time',group_id,got,error)
  call read_from_hdf5(npartoftypei(1:ntypes),'npartoftype',group_id,got,error)
  call read_from_hdf5(massoftypei(1:ntypes),'massoftype',group_id,got,error)
  call close_hdf5group(group_id,error)

  ndusttypes = ndustsmall + ndustlarge

  !-- Get number of columns in file
  call get_ncolumns(file_id,ncolstep,ndusttypes,ndustsmall,nptmass)

  !-- Add calculated quantities to data columns as well as 1 for density and 1 for mass
  ncolumns   = ncolstep + 2 + ncalc

  !--(re)allocate memory
  if (.not.allocated(dat) .or. npart+nptmass>maxpart .or. ncolumns>maxcol) then
    call alloc(npart+nptmass,1,ncolumns,mixedtypes=.true.)
  endif

  !-- Assign things now that memory has been allocated
  gamma(ifile)                = gammai
  time(ifile)                 = timei
  npartoftype(:,ifile)        = 0      ! initialise to zero
  masstype(:,ifile)           = 0.     ! I don't think this array is actually used anywhere in splash,
                                       ! since sinks all have different masses anyway

  npartoftype(1,ifile) = npartoftypei(1)                        ! gas
  npartoftype(2,ifile) = npartoftypei(2)                        ! dust (old)
  npartoftype(3,ifile) = nptmass                                ! sink
  npartoftype(4,ifile) = npartoftypei(3)                        ! boundary
  npartoftype(5,ifile) = npartoftypei(4)                        ! star
  npartoftype(6,ifile) = npartoftypei(5)                        ! dark matter
  npartoftype(7,ifile) = npartoftypei(6)                        ! bulge
  npartoftype(8,ifile) = sum(npartoftypei(7:7+ndusttypes-1))    ! dust
  npartoftype(9,ifile) = sum(npartoftypei(7+ndusttypes:ntypes)) ! unknown/dead?

  ! npart = sum(npartoftype(1:9,ifile))

  masstype(1,ifile) = massoftypei(1)                        ! gas
  masstype(2,ifile) = massoftypei(2)                        ! dust (old)
  ! masstype(3,ifile) = 0.                                    ! sink ???
  masstype(4,ifile) = massoftypei(3)                        ! boundary
  masstype(5,ifile) = massoftypei(4)                        ! star
  masstype(6,ifile) = massoftypei(5)                        ! dark matter
  masstype(7,ifile) = massoftypei(6)                        ! bulge
  ! masstype(8,ifile) = 0.                                    ! dust ???
  ! masstype(9,ifile) = 0.                                    ! unknown/dead ???

  print "(a,i10,a,es10.3,a,i2)",' npart = ',npart,' time = ',time(ifile),' ncolread = ',ncolstep

  !-- Try to read every possible particle array
  call open_hdf5group(file_id,'particles',group_id,error)

  allocate(itype(npart))
  call read_from_hdf5(itype,'itype',group_id,got,error)
  if (got) then
    iamtype(:,ifile) = 0
    iamtype(1:npart,ifile) = itypemap_phantom(itype(:))  ! Particles (remap phantom integer labels to splash labels)
    iamtype(npart+1:npart+nptmass,ifile) = 3             ! Sinks
  endif
  deallocate(itype)

  icol = 0

  allocate(array_1d(npart),array_2d(3,npart))

  call read_from_hdf5(array_2d,'xyz',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(1,:)
    ix(1) = icol
    label(icol) = 'x'
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(2,:)
    ix(2) = icol
    label(icol) = 'y'
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(3,:)
    ix(3) = icol
    label(icol) = 'z'
  endif

  !-- Set mass array for gas particles
  icol = icol + 1
  dat(1:npart,icol,ifile) = masstype(1,ifile)
  ipmass = icol
  label(icol) = 'm'

  call read_from_hdf5(array_1d,'h',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    ih = icol
    label(icol) = 'h'
  endif

  !-- Construct density array for each gas particle
  icol = icol + 1
  dat(1:npart,icol,ifile) = dat(1:npart,ipmass,ifile)*(hfact/dat(1:npart,ih,ifile))**3
  irho = icol
  label(icol) = 'density'

  call read_from_hdf5(array_2d,'vxyz',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(1,:)
    ivx = icol
    label(icol) = 'vx'
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(2,:)
    label(icol) = 'vy'
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(3,:)
    label(icol) = 'vz'
  endif

  call read_from_hdf5(array_1d,'u',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    iutherm = icol
    label(icol) = 'u'
  endif

  call read_from_hdf5(array_1d,'alpha',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'alpha'
  endif

  call read_from_hdf5(array_1d,'dt',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'dt'
  endif

  call read_from_hdf5(array_1d,'poten',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'poten'
  endif

  call read_from_hdf5(array_1d,'psi',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'psi'
  endif

  call read_from_hdf5(array_1d,'divB',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    idivB = icol
    label(icol) = 'divB'
  endif

  call read_from_hdf5(array_1d,'divBsymm',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'divBsymm'
  endif

  call read_from_hdf5(array_1d,'eta_OR',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'eta_{OR}'
  endif

  call read_from_hdf5(array_1d,'eta_HE',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'eta_{HE}'
  endif

  call read_from_hdf5(array_1d,'eta_AD',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'eta_{AD}'
  endif

  call read_from_hdf5(array_1d,'ne_on_n',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'ne/n'
  endif

  call read_from_hdf5(array_1d,'grainsize',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    igrainsize = icol
    label(icol) = 'grainsize'
  endif

  call read_from_hdf5(array_1d,'graindens',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    igraindens = icol
    label(icol) = 'graindens'
  endif

  call read_from_hdf5(array_1d,'vrel_on_vfrag',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'vrel/vfrag'
  endif

  call read_from_hdf5(array_1d,'st',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'st'
  endif

  call read_from_hdf5(array_1d,'temperature',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'temperature'
  endif

  call read_from_hdf5(array_1d,'divv',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'divv'
  endif

  call read_from_hdf5(array_1d,'luminosity',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'luminosity'
  endif

  call read_from_hdf5(array_1d,'beta',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_1d
    label(icol) = 'beta_pr'
  endif

  !-- Read other rank 2 arrays (3 x npart)
  call read_from_hdf5(array_2d,'Bxyz',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(1,:)
    iBfirst = icol
    label(icol) = 'Bx'
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(2,:)
    label(icol) = 'By'
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(3,:)
    label(icol) = 'Bz'
  endif

  call read_from_hdf5(array_2d,'curlBxyz',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(1,:)
    label(icol) = 'curlB_x'
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(2,:)
    label(icol) = 'curlB_y'
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(3,:)
    label(icol) = 'curlB_z'
  endif

  call read_from_hdf5(array_2d,'curlvxyz',group_id,got,error)
  if (got) then
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(1,:)
    label(icol) = 'curlv_x'
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(2,:)
    label(icol) = 'curlv_y'
    icol = icol + 1
    dat(1:npart,icol,ifile) = array_2d(3,:)
    label(icol) = 'curlv_z'
  endif

  deallocate(array_1d,array_2d)

!-------------

  !-- Read rank 2 arrays (ndusttypes x npart)
  ! TODO: add labels
  allocate(array_2d(ndusttypes,npart))
  call read_from_hdf5(array_2d,'dustfrac',group_id,got,error)
  if (got) then
    idustfrac = icol + 1
    do i=1,ndusttypes
      icol = icol + 1
      dat(1:npart,icol,ifile) = array_2d(i,1:npart)
    enddo
  endif
  call read_from_hdf5(array_2d,'tstop',group_id,got,error)
  if (got) then
    do i=1,ndusttypes
      icol = icol + 1
      dat(1:npart,icol,ifile) = array_2d(i,:)
    enddo
  endif
  deallocate(array_2d)

  !-- Read abundance array
  allocate(array_2d(5,npart))
  call read_from_hdf5(array_2d,'abundance',group_id,got,error)
  if (got) then
  abundance_labels = (/'h2ratio','abHIq  ','abhpq  ','abeq   ','abco   '/)
    do i=1,5
      icol = icol + 1
      dat(1:npart,icol,ifile) = array_2d(i,:)
      label(icol) = trim(abundance_labels(i))
    enddo
  endif
  deallocate(array_2d)

  !--- Read dust deltavxyz array
  ! TODO: add labels
  allocate(array_3d(3,ndustsmall,npart))
  call read_from_hdf5(array_3d,'deltavxyz',group_id,got,error)
  if (got) then
    do i=1,ndustsmall
      icol = icol + 1
      dat(1:npart,icol,ifile) = array_3d(1,i,:)
      icol = icol + 1
      dat(1:npart,icol,ifile) = array_3d(2,i,:)
      icol = icol + 1
      dat(1:npart,icol,ifile) = array_3d(3,i,:)
    enddo
  endif
  deallocate(array_3d)

  call close_hdf5group(group_id,error)

  !-- Try to read every possible sink array
  call open_hdf5group(file_id,'sinks',group_id,error)

  allocate(array_1d(nptmass))
  call read_from_hdf5(array_1d,'m',group_id,got,error)
  if (got) then
    dat(npart+1:npart+nptmass,ipmass,ifile) = array_1d(:)
  endif
  call read_from_hdf5(array_1d,'h',group_id,got,error)
  if (got) then
    dat(npart+1:npart+nptmass,ih,ifile) = array_1d(:)
  endif
  call read_from_hdf5(array_1d,'hsoft',group_id,got,error)
  call read_from_hdf5(array_1d,'maccreted',group_id,got,error)
  call read_from_hdf5(array_1d,'tlast',group_id,got,error)
  deallocate(array_1d)

  allocate(array_2d(3,nptmass))
  call read_from_hdf5(array_2d,'xyz',group_id,got,error)
  if (got) then
    dat(npart+1:npart+nptmass,ix(1),ifile) = array_2d(1,:)
    dat(npart+1:npart+nptmass,ix(2),ifile) = array_2d(2,:)
    dat(npart+1:npart+nptmass,ix(3),ifile) = array_2d(3,:)
  endif
  call read_from_hdf5(array_2d,'spinxyz',group_id,got,error)
  call read_from_hdf5(array_2d,'vxyz',group_id,got,error)
  if (got) then
    dat(npart+1:npart+nptmass,ivx  ,ifile) = array_2d(1,:)
    dat(npart+1:npart+nptmass,ivx+1,ifile) = array_2d(2,:)
    dat(npart+1:npart+nptmass,ivx+2,ifile) = array_2d(3,:)
  endif
  deallocate(array_2d)

  call close_hdf5group(group_id,error)

  call close_hdf5file(file_id,error)

end subroutine read_data_phantom_hdf5

subroutine set_labels_phantom_hdf5

end subroutine set_labels_phantom_hdf5
end module readdata_phantom_hdf5
