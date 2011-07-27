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
!  Copyright (C) 2005-2011 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!----------------------------------------------------------------------------
!
!  modules containing global variables
!
!----------------------------------------------------------------------------
!
!--global parameters
!
module params
 implicit none
 integer, parameter :: doub_prec = selected_real_kind(P=10,R=30)
 integer, parameter :: sing_prec = selected_real_kind(P=5,R=15)
 integer, parameter :: int1 = selected_int_kind(1)
 integer, parameter :: int8 = selected_int_kind(12)
 integer, parameter :: maxplot=64   ! maximum number of plots (for multiplot arrays)
 integer, parameter :: maxparttypes = 12  ! max # of different particle types

 public

end module params

module physcon
 use params, only:doub_prec
 implicit none
 real(doub_prec), parameter :: solarrcgs = 6.955d10
 real(doub_prec), parameter :: solarmcgs = 1.989d33

 public
end module physcon
!
!--particle data
!
module particle_data
 use params
 implicit none
 integer :: maxpart,maxstep,maxcol ! dimensions of dat array
 integer, allocatable, dimension(:)   :: icolourme
 integer(kind=int1), allocatable, dimension(:,:) :: iamtype
 integer, allocatable, dimension(:,:) :: npartoftype
 real, allocatable, dimension(:,:)    :: masstype
 real, allocatable, dimension(:)      :: time, gamma
 real, allocatable, dimension(:,:,:)  :: dat

 public

end module particle_data

module part_utils
 implicit none
 
 public :: igettype
 private

contains
 !
 !--utility returning the type of particle i
 !  when particles are ordered by type
 !
 integer function igettype(i,noftype)
  use params, only:maxparttypes
  implicit none
  integer, intent(in) :: i
  integer, dimension(maxparttypes), intent(in) :: noftype
  integer :: ntot,jtype

  ntot     = 0
  igettype = 1 ! so even if in error, will not lead to seg fault
  over_types: do jtype=1,maxparttypes
     ntot = ntot + noftype(jtype)
     if (i <= ntot) then
        igettype = i
        exit over_types
     endif
  enddo over_types

 end function igettype

end module part_utils
!
!--filenames
!
module filenames
 implicit none
 integer, parameter :: maxfile = 10001
 integer :: nfiles,nsteps,ifileopen
 character(len=120), dimension(maxfile) :: rootname
 character(len=100) :: fileprefix
 character(len=120) :: defaultsfile,limitsfile,animfile,unitsfile
 integer, dimension(maxfile) :: nstepsinfile
 character(len=68)  :: tagline = &
  'SPLASH: A visualisation tool for SPH data (c)2004-2011 Daniel Price'
 
 public

contains

 subroutine set_filenames(prefix)
  implicit none
  character(len=*), intent(in) :: prefix

  defaultsfile = trim(adjustl(prefix))//'.defaults'
  limitsfile   = trim(adjustl(prefix))//'.limits'
  animfile     = trim(adjustl(prefix))//'.anim'
  unitsfile    = trim(adjustl(prefix))//'.units'
  fileprefix   = trim(adjustl(prefix))
  
  return
 end subroutine set_filenames

end module filenames
!
!--labels for all plots and the locations of certain useful variables
!
module labels
 use params
 implicit none
 integer, parameter :: lenlabel = 80
 integer, parameter :: lenunitslabel = 40  ! length of units label
 character(len=lenlabel), dimension(maxplot+2) :: label,labelvec
 character(len=20), dimension(maxparttypes) :: labeltype
 character(len=6), parameter :: labeldefault = 'column'
 integer, dimension(3) :: ix
 integer, dimension(maxplot) :: iamvec
 integer :: ivx,irho,iutherm,ipr,ih,irad,iBfirst,iBpol,iBtor,iax
 integer :: ipmass,ike
 integer :: idivb,iJfirst,irhostar
 integer :: iacplane,ipowerspec
 integer :: icv,iradenergy
 integer :: isurfdens,itoomre
 integer :: ipdf,icolpixmap
 integer :: irhorestframe
 
 public

contains

 subroutine reset_columnids
  implicit none
  !
  !--array positions of specific quantities
  !  Identification is used in exact solution
  !  plotting and calculation of additional quantities
  !
  ix(:) = 0
  ivx = 0      ! vx
  irho = 0     ! density
  ipr = 0      ! pressure
  iutherm = 0  ! thermal energy
  ih = 0       ! smoothing length
  irad = 0     ! radius
  ipmass = 0   ! particle mass
  ipr = 0      ! pressure
  irad = 0     ! radius
  ipowerspec = 0 ! power spectrum
  iBfirst = 0  ! Bx
  iax = 0      ! ax (acceleration)
  iBpol = 0    ! B_polx
  iBtor = 0    ! B_torx
  iacplane = 0
  ike = 0
  idivB = 0
  iJfirst = 0
  icv = 0
  iradenergy = 0
  icolpixmap = 0
  irhorestframe = 0
  
  return
 end subroutine reset_columnids

 logical function is_coord(icol,ndim)
  implicit none
  integer, intent(in) :: icol,ndim
  integer :: i
  
  is_coord = .false.
  do i=1,ndim
     if (ix(i).eq.icol) is_coord = .true.
  enddo
  
 end function is_coord

end module labels

!------------------------------------
! modules containing plot settings
!------------------------------------
!
!--data
!
module settings_data
 use params
 implicit none
 integer :: numplot,ncalc,ncolumns,nextra
 integer :: ndataplots
 integer :: ndim, ndimv 
 integer :: icoords,icoordsnew,iformat,ntypes
 integer :: istartatstep,iendatstep,nfreq
 integer :: itrackpart
 integer, dimension(10) :: isteplist
 logical :: ivegotdata, DataIsBuffered, ipartialread
 logical :: buffer_data,iUseStepList,iCalcQuantities,iRescale
 !--required array is dimensioned 0:maxplot so that required(icol) = .true.
 !  does nothing bad if icol = 0 (much safer that way)
 logical :: lowmemorymode
 logical :: debugmode
 logical, dimension(0:maxplot) :: required
 logical, dimension(maxparttypes) :: UseTypeInRenderings
 real, dimension(3) :: xorigin

 namelist /dataopts/ buffer_data,iCalcQuantities,iRescale,xorigin

 public

end module settings_data
!
!--multiplot settings
!
module multiplot
 use params
 implicit none
 integer :: nyplotmulti 
 integer, dimension(maxplot) :: multiplotx,multiploty
 integer, dimension(maxplot) :: irendermulti,ivecplotmulti
 integer, dimension(maxplot) :: itrans,icontourmulti
 logical, dimension(maxplot) :: x_secmulti
 real, dimension(maxplot) :: xsecposmulti
 logical, dimension(maxplot) :: iusealltypesmulti
 logical, dimension(maxparttypes,maxplot) :: iplotpartoftypemulti
!
!--sort these into a namelist for input/output
!
 namelist /multi/ nyplotmulti,                           &
    itrans,multiplotx,multiploty,irendermulti,           &
    ivecplotmulti,icontourmulti,x_secmulti,xsecposmulti, &
    iusealltypesmulti,iplotpartoftypemulti

 public

end module multiplot
