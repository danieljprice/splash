!!
!!     set initial default options
!!      these are used if no defaults file is found
!!
subroutine defaults_set
  use exact_params
  use labels
  use multiplot
  use settings
  implicit none
  integer :: i
  !
  !--set default options
  !
  numplot=maxplot 	! reset if read from file
  ncalc = 0		! number of columns to calculate(e.g. radius)
  nextra = 0	! extra plots aside from particle data
  iACplane = 0	! label position of toy star AC plane plot
  ipowerspec = 0	! label position of power spectrum plot
  ncolumns=maxplot-ncalc	! number of columns in data file
  ndim = 3		! number of coordinate dimensions
  ndimV = ndim	! default velocity same dim as coords
  nstart = 1	! timestep to start from
  n_end = 1000	! timestep to finish on
  nfreq = 1		! frequency of timesteps to read
  icoords = 1	! co-ordinate system of simulation
  icoordsnew = icoords ! co-ordinate system to plot in

  iaxis = 0	! turns axes off/on
  iadapt = .true.	! adaptive plot limits
  plotcirc = .false.	! plot circle of radius 2h around particles
  plotcircall = .false.	!  " " around all particle
  icircpart = 1		!  " " around a specific particle
  ncircpart = 1
  interactive = .true.
  animate = .false.
  tile = .false.
  itrackpart = 0
  xminoffset_track = 0.5
  xmaxoffset_track = 0.5

  xsec_nomulti = .false.		! take cross section of data / particles
  flythru = .false.		! take series of cross sections through data
  ipagechange = .true.	! if false plots graphs on top of each other
  scalemax = 1.0	! for rescaling adaptive limits
  zoom = 1.0	! for rescaling fixed limits
  imark = 1	! PGPLOT marker for particles
  imarkg = 4	! PGPLOT marker for ghost particles
  imarksink = 17	! PGPLOT marker for sink particles 
  nacross = 1	! number of plots across page
  ndown = 1		! number of plots down page
  ipapersize = 0	! paper size option
  papersizex = 0.0	! size of x paper (no call to PGPAP if zero)
  aspectratio = 0.0	! aspect ratio of paper (no call to PGPAP if zero)
  iplotline = .false.	! plot line joining the particles
  iplotlinein = .false.	! " " but on first step only
  linestylein = 4		! PGPLOT line style for above
  iexact = 0		! exact solution to plot
  iplotav = .false.		! plot average line through particles
  nbins = 24		! number of bins for this
  ilabelpart = .false.	! plot particle numbers
  iplotpart = .true.	! flag whether or not to plot actual SPH particles
  iplotpartvec = .true.	! whether to plot particles on vector plot
  iplotghost = .true.	! plot ghost particles
  iplotsink = .true.	! plot sink particles
  npix = 100		! pixels in x direction for rendering
  npixvec = 40	! pixels in x direction on vector plots
  iplotcont_nomulti = .true.	! plot contours
  xsecpos_nomulti = 0.   ! position of cross section
  ncontours = 30		! number of contours to plot
  icolours = 0		! colour scheme to use
  ncolours=10		! number of colours in colour table
  itrans(:) = 0		! no transformations (log10 etc)
  UseBackgndColorVecplot = .false. ! plot vector plot using black/white
  iPlotColourBar = .true.

  hposlegend = 0.75     ! horizontal legend position as fraction of viewport
  vposlegend = 2.0      ! vertical legend position in character heights
  hpostitle = 0.5     ! horizontal title position as fraction of viewport
  vpostitle = 1.0      ! vertical title position in character heights
  fjusttitle = 0.5      ! justification factor for title
!
!--set coordinate labels for all coordinate systems
!
  labelcoord(1,1) = 'x'
  labelcoord(2,1) = 'y'
  labelcoord(3,1) = 'z'
  labelcoord(1,2) = 'r'
  labelcoord(2,2) = '\gphi'
  labelcoord(3,2) = 'z'
  labelcoord(1,3) = 'r'
  labelcoord(2,3) = '\gtheta'
  labelcoord(3,3) = '\gphi'
  
!
!--exact solution parameters
!
  lambda = 1.0	! sound wave exact solution : wavelength
  ampl = 0.005	! sound wave exact solution : amplitude
  period = 1.0
  iwaveplot = 5
  htstar = 1.   ! toy star crap
  atstar = 1.
  ctstar = 1.
  norder = 0
  sigma0 = 0.
  rhosedov = 1.0  ! sedov blast wave
  esedov = 1.0    ! blast wave energy
  polyk = 1.0     ! polytropic k
  rho_L = 1.0     ! shock tube (default is sod problem)
  rho_R = 0.125
  pr_L = 1.0
  pr_R = 0.1
  v_L = 0.0
  v_R = 0.0
  iexactpts = maxexactpts
  iexactplotx = 0
  iexactploty = 0
  ishk = 0
!
!--multiplot
!
  nyplotmulti = 4		! number of plots in multiplot
  multiploty(:) = 0
  do i=1,4
     multiploty(i) = ndim+i 	! first plot : y axis
  enddo
  multiplotx(:) = 1		! first plot : x axis
  irendermulti(:) = 0	! rendering
  ivecplotmulti(:) = 0	! vector plot
  x_secmulti(:) = .false.	! take cross section?
  xsecposmulti(:) = 0.0	! position of cross section
  iplotcontmulti(:) = .false.
  !
  !--array positions of specific quantities
  !
  ix = 0       ! coords
  ivx = 0      ! vx
  irho = 0     ! density
  ipr = 0      ! pressure
  iutherm = 0  ! thermal energy
  ih = 0       ! smoothing length
  ipmass = 0   ! particle mass
  ipr = 0      ! pressure
  irad = 0     ! radius
  ipmag = 0    ! magnetic pressure
  ibeta = 0    ! plasma beta
  itotpr = 0   ! total pressure
  ike = 0      ! specific kinetic energy
  idivBerr = 0 ! div B error
  iBfirst = 0  ! Bx
  irad2 = 0    ! r_parallel
  ivpar = 0    ! v_parallel
  ivperp = 0   ! v_perp
  iBpar = 0    ! B_parallel
  iBperp = 0   ! B_perp
  !
  !--power spectrum options
  !      
  idisordered = .false.
  ipowerspecy = ndim+1
  wavelengthmax = 1.0
  nfreqspec = 32

  return    
end subroutine defaults_set
