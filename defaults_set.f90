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
  !
  !--set default options
  !
  numplot=maxplot 	! reset if read from file
  ncalc = 0		! number of columns to calculate(e.g. radius)
  nextra = 0	! extra plots aside from particle data
  iACplane = 0	! label position of toy star AC plane plot
  ipowerspec = 0	! label position of power spectrum plot
  ncolumns=maxplot-ncalc	! number of columns in data file
  ndim = ndimmax		! number of coordinate dimensions
  ndimV = ndim	! default velocity same dim as coords
  nstart = 1	! timestep to start from
  n_end = maxstep	! timestep to finish on
  nfreq = 1		! frequency of timesteps to read
  icoords = 0	! co-ordinate system of simulation

  axes = .true.	! turns axes off/on
  animate = .true.	! turns off/on prompt between page changes
  magfield = .true.	! historical - can be used to set whether MHD or not
  iadapt = .true.	! adaptive plot limits
  plotcirc = .false.	! plot circle of radius 2h around particles
  plotcircall = .false.	!  " " around all particle
  icircpart = 1		!  " " around a specific particle
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
  iplotpartvec_nomulti = .true.	! whether to plot particles on vector plot
  iplotghost = .true.	! plot ghost particles
  iplotsink = .true.	! plot sink particles
  npix_nomulti = 100		! pixels in x direction for rendering
  npixvec_nomulti = 20	! pixels in x direction on vector plots
  ivecplot_nomulti = 0	! choice of vector plot
  iplotcont_nomulti = .true.	! plot contours
  xsecpos_nomulti = 0.   ! position of cross section
  ncontours_nomulti = 30		! number of contours to plot
  icolours = 0		! colour scheme to use
  ncolours=10		! number of colours in colour table
  itrans(:) = 0		! no transformations (log10 etc)
  backgnd_vec_nomulti = .false. ! plot vector plot using black/white

  isamexaxis = .true.	! these are overwritten later in program
  irenderplot = 0		! this is just so it is set to something
  labelcoord(1) = 'x'
  labelcoord(2) = 'y'
  labelcoord(3) = 'z'

  lambda = 1.0	! sound wave exact solution : wavelength
  delta = 0.005	! sound wave exact solution : amplitude

  nyplotmulti = 1		! number of plots in multiplot
  multiploty(1) = ndim+1 	! first plot : y axis
  multiplotx(1) = 1		! first plot : x axis
  irendermulti(:) = 0	! rendering
  ivecplotmulti(:) = 0	! vector plot
  x_secmulti(:) = .false.	! take cross section?
  xsecposmulti(:) = 0.0	! position of cross section
  npixmulti(:) = 400	! number of pixels in render plots
  npixvecmulti(:) = 40	! number of pixels in vector plots
  ncontoursmulti(:) = 30	! number of contours to use
  iplotcontmulti(:) = .false.
  iplotpartvecmulti(:) = .false.
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

  return    
end subroutine defaults_set
