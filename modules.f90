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
 integer, allocatable, dimension(:) :: npart,ntot,nghost,ntotplot
 integer, allocatable, dimension(:,:) :: iam,npartoftype
 real, allocatable, dimension(:) :: time, gamma
 real, allocatable, dimension(:,:,:) :: dat
 real :: hfact 
end module particle_data
!
!--filename
!
module filenames
 implicit none
 integer, parameter :: maxfile = 501
 integer :: ifile,nfilesteps,nfiles
 character(len=120), dimension(maxfile) :: rootname
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
 integer :: ivx,ivlast,irho,iutherm,ipr,ih,irad,ibfirst,iblast
 integer :: ipmass
 integer :: ientrop,ipmag,ibeta,itotpr,ike,idivb,idivberr,iJfirst
 integer :: iacplane,itimestep,ipowerspec
 integer :: irad2,ivpar,ivperp,iBpar,iBperp
end module labels
!
!--plot limits
!
module limits
 use params
 implicit none
 real, dimension(maxplot,2) :: lim
end module limits
!
!--module containing plot settings
!
module settings
 use params
 implicit none
!
!--global settings
!
 integer :: numplot,ncalc,ncolumns,nextra
 integer :: ndataplots
 integer :: ndim, ndimv 
 integer :: icoords,icoordsnew
 logical :: ishowopts, imulti
!
!--limits
! 
 integer :: itrackpart
 real :: scalemax,zoom
 real, dimension(3) :: xminoffset_track, xmaxoffset_track
 logical :: iadapt
!
!--data options
!
 integer :: nstart,n_end,nfreq
 logical :: ihavereadfilename, ivegotdata
!
!--particle plot options
!
 integer :: ntypes
 integer, dimension(maxparttypes) :: imarktype
 integer :: ncircpart
 integer, dimension(10) :: icircpart
 integer :: nbins,nc
 integer :: linestylein, iexact
 logical, dimension(maxparttypes) :: iplotpartoftype
 logical :: plotcirc,plotcircall
 logical :: iplotline,iplotlinein,iplotav,ilabelpart 
!
!--page options
!
 integer :: iaxis,nacross,ndown,ipapersize
 logical :: ipagechange,tile,animate,interactive
 real :: papersizex,aspectratio
 real :: hposlegend,vposlegend,hpostitle,vpostitle,fjusttitle
!
!--rendering options
!
 integer :: ncontours,npix,icolours,ncolours
 logical :: iplotcont_nomulti
 logical :: iPlotColourBar
!
!--vector plot options
!
 integer :: npixvec
 logical :: UseBackgndColorVecplot, iplotpartvec
!
!--cross section/rotation options
!
 integer :: ixsec,nxsec
 logical :: xsec_nomulti, irotate, flythru
 real :: anglerot, angletilt
 real :: xsecpos_nomulti,xseclineX1,xseclineX2,xseclineY1,xseclineY2
 real, dimension(3) :: xorigin
!
!--power spectrum options
!
 integer :: ipowerspecy, nfreqspec
 real :: wavelengthmax
 logical :: idisordered
!
!--sort these into namelists for input/output
!
 namelist /plotopts/ &
   iadapt,xsec_nomulti,flythru, &
   plotcirc,iplotline,iplotlinein,linestylein,          &
   imarktype,iplotpartoftype,                            &
   iexact,iplotav,nbins,                                &
   icolours,                      &
   ipowerspecy,idisordered,wavelengthmax,nfreqspec,icoordsnew, &
   ncircpart,icircpart

 namelist /pageopts/ iaxis,nacross,ndown,interactive, &
   ipagechange,tile,animate,ipapersize,papersizex,aspectratio, &
   hposlegend,vposlegend,hpostitle,vpostitle,fjusttitle  

 namelist /renderopts/ npix, ncontours,iplotcont_nomulti, &
   xsec_nomulti,iPlotColourBar,xsecpos_nomulti, &
   xseclineX1,xseclineX2,xseclineY1,xseclineY2, &
   irotate, anglerot, angletilt
 
 namelist /vectoropts/ npixvec, UseBackgndColorVecplot, iplotpartvec
     
end module settings
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
!--exact solution parameters
!
module exact_params
 implicit none
!--toy star
 integer :: norder ! for toy star
 real :: htstar,atstar,ctstar,totmass,sigma,sigma0
!--sound wave
 integer :: iwaveploty,iwaveplotx ! linear wave
 real :: ampl,lambda,period
!--sedov blast wave
 real :: rhosedov,esedov
!--polytrope
 real :: polyk
!--mhd shock solutions
 integer :: ishk
!--from file
 integer, parameter :: maxexactpts = 1001
 integer :: iexactpts, iexactplotx, iexactploty
 real, dimension(maxexactpts) :: xexact,yexact
!--shock tube
 real :: rho_L, rho_R, pr_L, pr_R, v_L, v_R
!
!--sort these into a namelist for input/output
!
 namelist /exactparams/ ampl,lambda,period,iwaveploty,iwaveplotx, &
          htstar,atstar,ctstar,sigma0,norder,rhosedov,esedov, &
	  rho_L, rho_R, pr_L, pr_R, v_L, v_R, &
	  iexactplotx,iexactploty 

end module exact_params
!
!--tabulated column density through the kernel 
!  (used in interpolate3D_projection)
!
module column
 implicit none
 integer, parameter :: maxcoltable = 1000
 real, dimension(maxcoltable) :: coltable
end module column
