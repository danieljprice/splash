!----------------------------------------------------------------------------
!
!  modules containing global variables
!
!----------------------------------------------------------------------------
!
!--global parameters  (should really allocate memory appropriately)
!
module params
 implicit none
 integer, parameter :: doub_prec = selected_real_kind(P=10,R=30)
 integer, parameter :: maxplot=40   ! maximum number of plots (for multiplot arrays)
 integer, parameter :: maxparttypes = 6  ! max # of different particle types
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
 real, allocatable, dimension(:) :: time, gamma
 real, allocatable, dimension(:,:,:) :: dat
end module particle_data
!
!--filename
!
module filenames
 implicit none
 integer, parameter :: maxfile = 501
 integer :: nfiles,nstepstotal,ifileopen
 character(len=120), dimension(maxfile) :: rootname
 integer, dimension(maxfile) :: nstepsinfile
end module filenames
!
!--labels for all plots and the locations of certain useful variables
!
module labels
 use params
 implicit none
 character(len=40), dimension(maxplot+2) :: label,labelvec
 character(len=20), dimension(maxparttypes) :: labeltype
 integer, dimension(3) :: ix
 integer, dimension(maxplot) :: iamvec
 integer :: ivx,irho,iutherm,ipr,ih,irad,iBfirst
 integer :: ipmass,ike
 integer :: idivb,iJfirst
 integer :: iacplane,ipowerspec
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
 integer :: icoords, iformat, ntypes
 integer :: nstart,n_end,nfreq
 integer, dimension(10) :: isteplist
 logical :: ivegotdata, DataIsBuffered
 logical :: buffer_data,iUseStepList, iCalcQuantities
 real, dimension(maxplot) :: units
 character(len=20), dimension(maxplot) :: unitslabel

 namelist /dataopts/ buffer_data, iCalcQuantities,units,unitslabel

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
 
end module multiplot
