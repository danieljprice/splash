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
 integer, allocatable, dimension(:,:) :: iam
 real, allocatable, dimension(:) :: time, gamma
 real, allocatable, dimension(:,:,:) :: dat
 real :: hfact 
end module particle_data
!
!--filename
!
module filenames
 implicit none
 integer, parameter :: maxfile = 10
 integer :: ifile,nfilesteps,nfiles
 character(len=20), dimension(maxfile) :: rootname
end module filenames
!
!--labels for all plots and the locations of certain useful variables
!
module labels
 use params
 implicit none
 character(len=20), dimension(maxplot+2) :: label
 character(len=7), dimension(3,3) :: labelcoord
 integer, dimension(3) :: ix
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
 real :: bmin,bmax 
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
 integer :: imark, imarkg, imarksink
 integer :: ixsec,nxsec, nbins,nc
 integer :: linestylein, iexact
 integer :: ncolours,nstart,n_end,nfreq
 integer :: icoords,icoordsnew
 integer :: ncircpart, itrackpart
 integer, dimension(10) :: icircpart
 integer :: icolours
!
!--limits
! 
 real :: scalemax,zoom
 real, dimension(3) :: xminoffset_track, xmaxoffset_track

!--plot options
 logical :: interactive
 logical :: iadapt,ihavereadfilename
 logical :: plotcirc,plotcircall,flythru,imulti
 logical :: iplotline,iplotlinein,iplotav,ilabelpart
 logical :: iplotpart,iplotghost,iplotsink
 !!logical, dimension(maxparttypes) :: iplotparttype
 logical :: ishowopts, ivegotdata
!
!--page options
!
 integer :: iaxis,nacross,ndown,ipapersize
 logical :: ipagechange,tile,animate
 real :: papersizex,aspectratio
 real :: hposlegend,vposlegend,hpostitle,vpostitle,fjusttitle
!
!--rendering options
!
 integer :: ncontours,npix
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
 logical :: xsec_nomulti, irotate
 real :: xsecpos_nomulti,xseclineX1,xseclineX2,xseclineY1,xseclineY2
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
   imark, imarkg, imarksink,                            &
   iexact,iplotav,nbins,                                &
   icolours,iplotghost,iplotsink,                       &
   ipowerspecy,idisordered,wavelengthmax,nfreqspec,icoordsnew, &
   ncircpart,icircpart

 namelist /pageopts/ iaxis,nacross,ndown, &
   ipagechange,tile,animate,ipapersize,papersizex,aspectratio, &
   hposlegend,vposlegend,hpostitle,vpostitle,fjusttitle  

 namelist /renderopts/ npix, ncontours,iplotcont_nomulti, &
   xsec_nomulti,iPlotColourBar,xsecpos_nomulti, &
   xseclineX1,xseclineX2,xseclineY1,xseclineY2
 
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
!  (used in interpolate3d_projection)
!
module column
 implicit none
 integer, parameter :: maxcoltable = 1000
 real, dimension(maxcoltable) :: coltable
end module column
