!----------------------------------------------------------------------------
!
!  modules containing global variables
!
!----------------------------------------------------------------------------
!
!--global parameters  (should really allocate memory appropriately)
!
      MODULE params
       IMPLICIT NONE
       INTEGER, PARAMETER :: max=15000
       INTEGER, PARAMETER :: maxstep=350
       INTEGER, PARAMETER :: ndimmax = 3
       INTEGER, PARAMETER :: maxplot=24+2*ndimmax + 1	! maximum number of plots       
      END MODULE params
!
!--particle data
!      
      MODULE particle_data
       USE params
       IMPLICIT NONE
       INTEGER, DIMENSION(maxstep) :: npart,ntot,nghost,ntotplot
       INTEGER, DIMENSION(max) :: iam
       REAL, DIMENSION(maxstep) :: time, gamma
       REAL, DIMENSION(maxplot,max,maxstep) :: dat
       REAL, DIMENSION(maxplot,2) :: lim
       REAL :: hfact       
       REAL :: Bmin,Bmax
       
      END MODULE particle_data
!
!--filename
!      
      MODULE filenames
       IMPLICIT NONE
       INTEGER :: ifile,nfilesteps
       CHARACTER(LEN=20) :: rootname
      END MODULE filenames
!
!--labels for all plots and the locations of certain useful variables
!
      MODULE labels
       USE params
       IMPLICIT NONE
       CHARACTER(LEN=20), DIMENSION(maxplot+2) :: label
       CHARACTER(LEN=1), DIMENSION(3) :: labelcoord
       INTEGER, DIMENSION(3) :: ix
       INTEGER :: ivx,ivlast,irho,iutherm,ipr,ih,irad,iBfirst,iBlast
       INTEGER :: ipmass
       INTEGER :: ientrop,irad2,ipmag,ibeta,itotpr,ike,idivB,idivBerr
       INTEGER :: iextra,itimestep
      END MODULE labels
!
!--module containing plot settings
!
      MODULE settings
       USE params
       IMPLICIT NONE
       INTEGER :: numplot,ncalc,ncolumns,nextra,nacross,ndown
       INTEGER :: ndataplots
       INTEGER :: ndim, ndimV       
       INTEGER :: imark, imarkg, imarksink
       INTEGER :: ixsec,nxsec, nbins,nc,icircpart
       INTEGER :: linestylein, iexact
       INTEGER :: irenderplot
       INTEGER :: ncolours,nstart,n_end,nfreq
       
       INTEGER :: ncontours_nomulti,npix_nomulti,npixvec_nomulti
       INTEGER :: ivecplot_nomulti,irender,icolours
       INTEGER :: ipapersize,menuitems
       
       REAL :: scalemax,zoom
       REAL :: xsecpos_nomulti
       REAL :: papersizex,aspectratio
!--plot options
       LOGICAL :: animate,iadapt,magfield,iexist,ihavereadfilename
       LOGICAL :: plotcirc,plotcircall,flythru,imulti,ipagechange
       LOGICAL :: iplotline,iplotlinein,iplotav,ilabelpart
       LOGICAL :: iplotpart,iplotghost,iplotsink
       LOGICAL :: ishowopts, ivegotdata,isamexaxis
       
       LOGICAL :: backgnd_vec_nomulti
       LOGICAL :: iplotcont_nomulti,xsec_nomulti,iplotpartvec_nomulti

      END MODULE settings
!
!--multiplot settings
!      
      MODULE multiplot
       USE params
       IMPLICIT NONE
       INTEGER, DIMENSION(maxplot) :: multiplotx,multiploty
       INTEGER, DIMENSION(maxplot) :: irendermulti,ivecplotmulti      
       INTEGER, DIMENSION(maxplot) :: npixmulti,npixvecmulti 
       INTEGER, DIMENSION(maxplot) :: ncontoursmulti,itrans
       INTEGER :: nyplotmulti
       LOGICAL, DIMENSION(maxplot) :: iplotcontmulti, iplotpartvecmulti     
       LOGICAL, DIMENSION(maxplot) :: x_secmulti,backgnd_vec_multi
       REAL, DIMENSION(maxplot) :: xsecposmulti
      END MODULE multiplot
!      
!--exact solution parameters      
!
      MODULE exact_params
       IMPLICIT NONE
       INTEGER, PARAMETER :: ipolycmax=1000
!--toy star
       INTEGER :: norder
       REAL :: Htstar,Atstar,Ctstar,totmass,sigma,sigma0 ! toy star parameters
!--sound wave      
       REAL delta,lambda
!--polytrope
       INTEGER :: ipolyc
       REAL mtot,maxrho,akfac,den(ipolycmax),rad(ipolycmax)       
      END MODULE exact_params
!
!--tabulated column density through the kernel 
!  (used in interpolate3D_projection)
!
      MODULE column
       IMPLICIT NONE
       INTEGER, PARAMETER :: maxcoltable = 1000
       REAL, PARAMETER :: dmaxcoltable = 1./FLOAT(maxcoltable)
       REAL, DIMENSION(maxcoltable) :: coltable
      END MODULE column
