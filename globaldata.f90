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
 integer, parameter :: maxplot=64   ! maximum number of plots (for multiplot arrays)
 integer, parameter :: maxparttypes = 6  ! max # of different particle types

 public

end module params
!
!--particle data
!
module particle_data
 use params
 implicit none
 integer :: maxpart,maxstep,maxcol ! dimensions of dat array
 integer, allocatable, dimension(:) :: icolourme
 integer, allocatable, dimension(:,:) :: npartoftype
 real, allocatable, dimension(:,:) :: massoftype
 real, allocatable, dimension(:) :: time, gamma
 real, allocatable, dimension(:,:,:) :: dat

 public

end module particle_data
!
!--filename
!
module filenames
 implicit none
 integer, parameter :: maxfile = 10001
 integer :: nfiles,nsteps,ifileopen
 character(len=120), dimension(maxfile) :: rootname
 character(len=100) :: fileprefix
 character(len=120) :: defaultsfile,limitsfile,animfile,unitsfile
 integer, dimension(maxfile) :: nstepsinfile
 
 public

end module filenames
!
!--labels for all plots and the locations of certain useful variables
!
module labels
 use params
 implicit none
 character(len=60), dimension(maxplot+2) :: label,labelvec
 character(len=20), dimension(maxparttypes) :: labeltype
 integer, dimension(3) :: ix
 integer, dimension(maxplot) :: iamvec
 integer :: ivx,irho,iutherm,ipr,ih,irad,iBfirst
 integer :: ipmass,ike
 integer :: idivb,iJfirst
 integer :: iacplane,ipowerspec
 integer :: icv,iradenergy
 
 public

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
 integer, dimension(maxplot) :: itrans
 logical, dimension(maxplot) :: iplotcontmulti, x_secmulti
 real, dimension(maxplot) :: xsecposmulti
!
!--sort these into a namelist for input/output
!
 namelist /multi/ nyplotmulti,                                  &
    itrans,multiplotx,multiploty,irendermulti,                  &
    ivecplotmulti,iplotcontmulti,x_secmulti,xsecposmulti

 public

end module multiplot
