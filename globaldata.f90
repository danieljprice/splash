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
 integer, allocatable, dimension(:) :: ntot, icolourme
 integer, allocatable, dimension(:,:) :: iam,npartoftype
 real, allocatable, dimension(:) :: time, gamma
 real, allocatable, dimension(:,:,:) :: dat
end module particle_data
!
!--filename
!
module filenames
 implicit none
 integer, parameter :: maxfile = 501
 integer :: nfiles,nstepstotal
 character(len=120), dimension(maxfile) :: rootname
 integer, dimension(maxfile) :: nstepsinfile
end module filenames
!
!--labels for all plots and the locations of certain useful variables
!
module labels
 use params
 implicit none
 character(len=20), dimension(maxplot+2) :: label,labelvec
 character(len=7), dimension(3,3) :: labelcoord
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
 implicit none
 integer :: numplot,ncalc,ncolumns,nextra
 integer :: ndataplots
 integer :: ndim, ndimv 
 integer :: icoords, ntypes
 integer :: nstart,n_end,nfreq
 logical :: ivegotdata, buffer_data
 logical :: imulti

 namelist /dataopts/ buffer_data

end module settings_data

!
!--limits
! 
module settings_limits
 implicit none
 integer :: itrackpart
 real :: scalemax,zoom
 real, dimension(3) :: xminoffset_track, xmaxoffset_track
end module settings_limits
!
!--particle plot options
!
module settings_part
 use params
 implicit none
 integer, dimension(maxparttypes) :: imarktype
 integer :: ncircpart, icoordsnew
 integer, dimension(10) :: icircpart
 integer :: nc
 integer :: linestylein, iexact
 logical, dimension(maxparttypes) :: iplotpartoftype
 logical :: iplotline,iplotlinein,ilabelpart

 namelist /plotopts/ iplotline,iplotlinein,linestylein,  &
   imarktype,iplotpartoftype,iexact, &
   ncircpart,icircpart

end module settings_part
!
!--page options
!
module settings_page
 implicit none
 integer :: iaxis,nacross,ndown,ipapersize
 logical :: ipagechange,tile,animate,interactive,iadapt
 real :: papersizex,aspectratio
 real :: hposlegend,vposlegend,hpostitle,vpostitle,fjusttitle

 namelist /pageopts/ iaxis,nacross,ndown,interactive,iadapt, &
   ipagechange,tile,animate,ipapersize,papersizex,aspectratio, &
   hposlegend,vposlegend,hpostitle,vpostitle,fjusttitle  

end module settings_page
!
!--rendering options
!
module settings_render
 implicit none
 integer :: ncontours,npix,icolours
 logical :: iplotcont_nomulti
 logical :: iPlotColourBar

 namelist /renderopts/ npix,icolours,ncontours,iplotcont_nomulti, &
   iPlotColourBar

end module settings_render
!
!--vector plot options
!
module settings_vecplot
 implicit none
 integer :: npixvec
 logical :: UseBackgndColorVecplot, iplotpartvec

 namelist /vectoropts/ npixvec, UseBackgndColorVecplot, iplotpartvec

end module settings_vecplot
!
!--cross section/rotation options
!
module settings_xsecrot
 implicit none
 integer :: ixsec,nxsec
 logical :: xsec_nomulti, irotate, flythru
 real :: anglex, angley, anglez
 real :: xsecpos_nomulti,xseclineX1,xseclineX2,xseclineY1,xseclineY2
 real, dimension(3) :: xorigin

 namelist /xsecrotopts/ xsec_nomulti,xsecpos_nomulti,flythru, &
          xseclineX1,xseclineX2,xseclineY1,xseclineY2, &
          irotate, anglex, angley, anglez

end module settings_xsecrot
!
!--power spectrum options
!
module settings_powerspec
 implicit none
 integer :: ipowerspecy, nfreqspec
 logical :: idisordered
 real :: wavelengthmax
 
 namelist /powerspecopts/ ipowerspecy,idisordered,wavelengthmax,nfreqspec
 
end module settings_powerspec

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

!
!--tabulated column density through the kernel 
!  (used in interpolate3D_projection)
!
module column
 implicit none
 integer, parameter :: maxcoltable = 1000
 real, dimension(maxcoltable) :: coltable
end module column
