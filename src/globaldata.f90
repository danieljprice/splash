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
 integer, parameter :: int8 = selected_int_kind(10)
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
 real, parameter :: time_not_read_val = -0.5*huge(0.)

 public

contains

 logical function time_was_read(t)
  real, intent(in) :: t
  
  time_was_read = .true.
  if (t <= time_not_read_val) time_was_read = .false.
  
 end function time_was_read

end module particle_data
!
!--filenames
!
module filenames
 implicit none
 integer, parameter :: maxfile = 10001
 integer :: nfiles,nsteps,ifileopen
 character(len=120), dimension(maxfile) :: rootname
 character(len=100) :: fileprefix
 character(len=120) :: defaultsfile,limitsfile,unitsfile
 integer, dimension(maxfile) :: nstepsinfile
 character(len=68)  :: tagline = &
  'SPLASH: A visualisation tool for SPH data (c)2004-2015 Daniel Price'

 public

contains

 subroutine set_filenames(prefix)
  implicit none
  character(len=*), intent(in) :: prefix

  defaultsfile = trim(adjustl(prefix))//'.defaults'
  limitsfile   = trim(adjustl(prefix))//'.limits'
  unitsfile    = trim(adjustl(prefix))//'.units'
  fileprefix   = trim(adjustl(prefix))

  return
 end subroutine set_filenames

end module filenames

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
 integer :: icoords,icoordsnew,iformat,ntypes,iexact
 integer :: istartatstep,iendatstep,nfreq
 integer :: itracktype,itrackoffset,iverbose
 integer, dimension(10) :: isteplist
 logical :: ivegotdata, DataIsBuffered, ipartialread
 logical :: buffer_data,iUseStepList,iCalcQuantities,iRescale
 logical, parameter :: buffer_steps_in_file = .false.
 !--required array is dimensioned 0:maxplot so that required(icol) = .true.
 !  does nothing bad if icol = 0 (much safer that way)
 logical :: lowmemorymode
 logical :: debugmode
 logical, dimension(0:maxplot) :: required
 logical, dimension(maxparttypes) :: UseTypeInRenderings
 real, dimension(3) :: xorigin

 namelist /dataopts/ buffer_data,iCalcQuantities,iRescale,xorigin,itracktype,itrackoffset

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
