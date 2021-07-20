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
!  Copyright (C) 2005-2020 Daniel Price. All rights reserved.
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
 integer, parameter :: maxplot=512   ! maximum number of plots (for multiplot arrays)
 integer, parameter :: maxparttypes = 24  ! max # of different particle types
 integer, parameter :: ltag = 16     ! length of header tags
 integer, parameter :: maxhdr = 256  ! maximum number of header variables stored

 public

end module params

module physcon
 use params, only:doub_prec
 implicit none
 real, parameter :: pi = 4.*atan(1.)
 real(doub_prec), parameter :: solarrcgs = 6.955d10  ! cm
 real(doub_prec), parameter :: solarmcgs = 1.989d33  ! g
 real(doub_prec), parameter :: steboltz = 5.67e-5    ! erg cm^-2 K^-4 s-1
 real(doub_prec), parameter :: radconst = 7.5646d-15 ! Radiation constant erg cm^-3 K^-4
 real(doub_prec), parameter :: kboltz = 1.38066d-16
 real(doub_prec), parameter :: mh = 1.67262158d-24   ! g
 real(doub_prec), parameter :: au = 1.496d13         ! cm
 real(doub_prec), parameter :: c = 2.997924d10       ! Speed of light cm/s
 real(doub_prec), parameter :: hplanck =   6.6260755d-27  ! Planck's Constant erg/s
 real(doub_prec), parameter :: kb_on_mh = kboltz/mh
 real(doub_prec), parameter :: Lsun = 3.839d33       ! Solar luminosity, erg/s
 real(doub_prec), parameter :: cm_to_nm = 1.d7
 real(doub_prec), parameter :: keV_to_erg = 1.6022d-9 ! k_eV to erg
 real(doub_prec), parameter :: keV_to_Hz  = keV_to_erg/hplanck ! k_eV to erg

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
 real, allocatable, dimension(:,:)    :: headervals
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
 integer :: nfiles,nsteps,ifileopen,iposopen
 character(len=120), dimension(maxfile) :: rootname
 character(len=100) :: fileprefix
 character(len=120) :: defaultsfile,limitsfile,unitsfile,coloursfile
 integer, dimension(maxfile) :: nstepsinfile
 character(len=*), parameter :: tagline = &
  'SPLASH: A visualisation tool for SPH data (c)2004-2021 Daniel Price and contributors'

 public

contains

subroutine set_filenames(prefix)
 character(len=*), intent(in) :: prefix

 fileprefix   = trim(adjustl(prefix))
 if (fileprefix(len_trim(fileprefix):len_trim(fileprefix))=='.') then
    fileprefix = fileprefix(1:len_trim(fileprefix)-1)
 endif
 defaultsfile = trim(adjustl(fileprefix))//'.defaults'
 limitsfile   = trim(adjustl(fileprefix))//'.limits'
 unitsfile    = trim(adjustl(fileprefix))//'.units'
 coloursfile  = trim(adjustl(fileprefix))//'.colours'

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
 integer :: ndusttypes
 integer :: idustfrac_plot = 0
 integer :: ideltav_plot = 0
 integer :: iautorender
 integer :: icoords,icoordsnew,iformat,ntypes,iexact
 integer :: istartatstep,iendatstep,nfreq
 integer :: itracktype,itrackoffset,iverbose
 integer, dimension(10) :: isteplist
 logical :: ivegotdata, DataIsBuffered, ipartialread
 logical :: buffer_data,iUseStepList,iCalcQuantities,iRescale
 logical :: idefaults_file_read
 logical :: buffer_steps_in_file = .false.
 !--required array is dimensioned 0:maxplot so that required(icol) = .true.
 !  does nothing bad if icol = 0 (much safer that way)
 logical :: lowmemorymode
 logical :: debugmode
 logical :: UseFakeDustParticles, UseFastRender
 logical, dimension(0:maxplot) :: required
 logical, dimension(maxparttypes) :: UseTypeInRenderings
 real, dimension(3) :: xorigin
 character(len=120) :: device

 namelist /dataopts/ buffer_data,iCalcQuantities,iRescale,xorigin, &
                     itracktype,itrackoffset,idustfrac_plot,ideltav_plot,UseFakeDustParticles

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
