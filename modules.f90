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
 integer, parameter :: max=55000
 integer, parameter :: maxstep=100
 integer, parameter :: ndimmax = 3
 integer, parameter :: maxplot=24+2*ndimmax + 1	! maximum number of plots 
end module params
!
!--particle data
!
module particle_data
 use params
 implicit none
 integer, dimension(maxstep) :: npart,ntot,nghost,ntotplot
 integer, dimension(max) :: iam
 real, dimension(maxstep) :: time, gamma
 real, dimension(maxplot,max,maxstep) :: dat
 real, dimension(maxplot,2) :: lim
 real :: hfact 
 real :: bmin,bmax
 
end module particle_data
!
!--filename
!
module filenames
 implicit none
 integer :: ifile,nfilesteps
 character(len=20) :: rootname
end module filenames
!
!--labels for all plots and the locations of certain useful variables
!
module labels
 use params
 implicit none
 character(len=20), dimension(maxplot+2) :: label
 character(len=1), dimension(3) :: labelcoord
 integer, dimension(3) :: ix
 integer :: ivx,ivlast,irho,iutherm,ipr,ih,irad,ibfirst,iblast
 integer :: ipmass
 integer :: ientrop,irad2,ipmag,ibeta,itotpr,ike,idivb,idivberr
 integer :: iacplane,itimestep,ipowerspec
end module labels
!
!--module containing plot settings
!
module settings
 use params
 implicit none
 integer :: numplot,ncalc,ncolumns,nextra,nacross,ndown
 integer :: ndataplots
 integer :: ndim, ndimv 
 integer :: imark, imarkg, imarksink
 integer :: ixsec,nxsec, nbins,nc,icircpart
 integer :: linestylein, iexact
 integer :: irenderplot
 integer :: ncolours,nstart,n_end,nfreq
 
 integer :: ncontours_nomulti,npix_nomulti,npixvec_nomulti
 integer :: ivecplot_nomulti,irender,icolours
 integer :: ipapersize,menuitems
 
 real :: scalemax,zoom
 real :: xsecpos_nomulti
 real :: papersizex,aspectratio
!--plot options
 logical :: animate,iadapt,magfield,iexist,ihavereadfilename
 logical :: plotcirc,plotcircall,flythru,imulti,ipagechange
 logical :: iplotline,iplotlinein,iplotav,ilabelpart
 logical :: iplotpart,iplotghost,iplotsink
 logical :: ishowopts, ivegotdata,isamexaxis
 
 logical :: backgnd_vec_nomulti
 logical :: iplotcont_nomulti,xsec_nomulti,iplotpartvec_nomulti

!
!--sort these into a namelist for input/output
!
 namelist /plotopts/ animate,magfield,iadapt,xsec_nomulti,flythru, &
   plotcirc,iplotline,iplotlinein,linestylein, 		&
   imark, imarkg, imarksink,				&
   nacross,ndown,					&
   iexact,iplotav,nbins, 				&
   irender,ivecplot_nomulti,iplotpartvec_nomulti, 	&
   npix_nomulti,npixvec_nomulti,			&
   iplotcont_nomulti,ncontours_nomulti,			&
   icolours,iplotghost,iplotsink,			&
   ipapersize,papersizex,aspectratio
     
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
 integer, dimension(maxplot) :: npixmulti,npixvecmulti 
 integer, dimension(maxplot) :: ncontoursmulti,itrans
 logical, dimension(maxplot) :: iplotcontmulti, iplotpartvecmulti     
 logical, dimension(maxplot) :: x_secmulti,backgnd_vec_multi
 real, dimension(maxplot) :: xsecposmulti
!
!--sort these into a namelist for input/output
!
 namelist /multi/ nyplotmulti,					&
    itrans,multiplotx,multiploty,irendermulti, 			&
    iplotcontmulti,ncontoursmulti,ivecplotmulti,npixmulti,	&
    npixvecmulti,iplotpartvecmulti,x_secmulti,xsecposmulti
 
end module multiplot
!
!--exact solution parameters
!
module exact_params
 implicit none
 integer, parameter :: ipolycmax=1000
!--toy star
 integer :: norder
 real :: htstar,atstar,ctstar,totmass,sigma,sigma0 ! toy star parameters
!--sound wave
 real delta,lambda
!--polytrope
 integer :: ipolyc
 real mtot,maxrho,akfac,den(ipolycmax),rad(ipolycmax) 
!
!--sort these into a namelist for input/output
!
 namelist /exactparams/ delta,lambda,htstar,atstar,ctstar,sigma0,norder
     

end module exact_params
!
!--tabulated column density through the kernel 
!  (used in interpolate3d_projection)
!
module column
 implicit none
 integer, parameter :: maxcoltable = 1000
 real, parameter :: dmaxcoltable = 1./float(maxcoltable)
 real, dimension(maxcoltable) :: coltable
end module column
