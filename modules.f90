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
end module params
!
!--particle data
!
module particle_data
 use params
 implicit none
 integer :: icrap
 integer :: maxpart,maxstep,maxcol  ! dimensions of dat array
 integer, allocatable, dimension(:) :: npart,ntot,nghost,ntotplot
 integer, allocatable, dimension(:,:) :: iam
 real, allocatable, dimension(:) :: time, gamma
 real, allocatable, dimension(:,:,:) :: dat
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
 integer :: ientrop,ipmag,ibeta,itotpr,ike,idivb,idivberr
 integer :: iacplane,itimestep,ipowerspec
 integer :: irad2,ivpar,ivperp,iBpar,iBperp
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
 integer :: ixsec,nxsec, nbins,nc
 integer :: linestylein, iexact
 integer :: irenderplot
 integer :: ncolours,nstart,n_end,nfreq
 integer :: icoords
 integer :: ncircpart
 integer, dimension(10) :: icircpart
 
 integer :: ncontours_nomulti,npix_nomulti,npixvec_nomulti
 integer :: ivecplot_nomulti,icolours
 integer :: ipapersize,menuitems
 
 real :: scalemax,zoom
 real :: xsecpos_nomulti
 real :: papersizex,aspectratio
!--plot options
 logical :: axes, animate, interactive
 logical :: iadapt,ihavereadfilename
 logical :: plotcirc,plotcircall,flythru,imulti,ipagechange
 logical :: iplotline,iplotlinein,iplotav,ilabelpart
 logical :: iplotpart,iplotghost,iplotsink
 logical :: ishowopts, ivegotdata,isamexaxis
 
 logical :: backgnd_vec_nomulti
 logical :: iplotcont_nomulti,xsec_nomulti,iplotpartvec_nomulti
!
!--power spectrum options
!
 integer :: ipowerspecy, nfreqspec
 real :: wavelengthmax
 logical :: idisordered

!
!--sort these into a namelist for input/output
!
 namelist /plotopts/ axes, &
   animate,iadapt,xsec_nomulti,flythru, &
   plotcirc,iplotline,iplotlinein,linestylein,          &
   imark, imarkg, imarksink,                            &
   nacross,ndown,                                       &
   iexact,iplotav,nbins,                                &
   ivecplot_nomulti,iplotpartvec_nomulti,       &
   npix_nomulti,npixvec_nomulti,                        &
   iplotcont_nomulti,ncontours_nomulti,                 &
   icolours,iplotghost,iplotsink,                       &
   ipapersize,papersizex,aspectratio,                   &
   ipowerspecy,idisordered,wavelengthmax,nfreqspec,icoords, &
   ncircpart,icircpart
     
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
 namelist /multi/ nyplotmulti,                                  &
    itrans,multiplotx,multiploty,irendermulti,                  &
    iplotcontmulti,ncontoursmulti,ivecplotmulti,npixmulti,      &
    npixvecmulti,iplotpartvecmulti,x_secmulti,xsecposmulti
 
end module multiplot
!
!--exact solution parameters
!
module exact_params
 implicit none
 integer :: norder ! for toy star
 integer, parameter :: ipolycmax=1000
 real :: htstar,atstar,ctstar,totmass,sigma,sigma0 ! toy star parameters
!--sound wave
 real :: delta,lambda
!--sedov blast wave
 real :: rhosedov,esedov
!--polytrope
 integer :: ipolyc
 real :: mtot,maxrho,akfac
 real, dimension(ipolycmax) :: den,rad 
!
!--sort these into a namelist for input/output
!
 namelist /exactparams/ delta,lambda,htstar,atstar,ctstar,sigma0,norder,rhosedov,esedov
     

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
