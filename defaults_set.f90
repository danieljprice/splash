!!
!!     set initial default options
!!      these are used if no defaults file is found
!!
subroutine defaults_set
  use exact, only:defaults_set_exact
  use filenames
  use labels
  use limits
  use multiplot
  use settings_limits
  use settings_data
  use settings_part
  use settings_page
  use settings_render
  use settings_vecplot
  use settings_xsecrot
  use settings_powerspec
  use particle_data
  implicit none
  integer :: i
!
!--data options (most should be set upon call to read_data)
!
  numplot=maxplot   ! reset if read from file
  ncalc = 0         ! number of columns to calculate(e.g. radius)
  nextra = 0        ! extra plots aside from particle data
  ipowerspec = 0    ! label position of power spectrum plot
  ncolumns=maxplot-ncalc        ! number of columns in data file
  ndim = 3          ! number of coordinate dimensions
  ndimV = ndim      ! default velocity same dim as coords
  nstart = 1        ! timestep to start from
  n_end = 1000      ! timestep to finish on
  nfreq = 1         ! frequency of timesteps to read
  icoords = 1       ! co-ordinate system of simulation
  icoordsnew = icoords ! co-ordinate system to plot in
  buffer_data = .false.
!
!--default for interactive mode
!
  interactive = .true.
!
!--limits options
!
  iadapt = .true.      ! adaptive plot limits
  scalemax = 1.0       ! for rescaling adaptive limits
  zoom = 1.0           ! for rescaling fixed limits
  itrans(:) = 0        ! no transformations (log10 etc)
  itrackpart = 0       ! particle to track (none)
  xminoffset_track = 0.5 ! offset of limits from tracked particle
  xmaxoffset_track = 0.5 !
!
!--page options
!
  iaxis = 0                ! turns axes off/on
  ipagechange = .true.     ! if false plots graphs on top of each other
  animate = .false.
  tile = .false.
  nacross = 1           ! number of plots across page
  ndown = 1             ! number of plots down page
  ipapersize = 0        ! paper size option
  papersizex = 0.0      ! size of x paper (no call to PGPAP if zero)
  aspectratio = 0.0     ! aspect ratio of paper (no call to PGPAP if zero)
  hposlegend = 0.75     ! horizontal legend position as fraction of viewport
  vposlegend = 2.0      ! vertical legend position in character heights
  hpostitle = 0.5       ! horizontal title position as fraction of viewport
  vpostitle = 1.0       ! vertical title position in character heights
  fjusttitle = 0.5      ! justification factor for title
!
!--particle plot options
!
  ncircpart = 0
  iplotline = .false.     ! plot line joining the particles
  iplotlinein = .false.   ! " " but on first step only
  linestylein = 4         ! PGPLOT line style for above
  iexact = 0              ! exact solution to plot
  ilabelpart = .false.    ! plot particle numbers
  iplotpartvec = .true.   ! whether to plot particles on vector plot
  
  iplotpartoftype(1) = .true. ! whether or not to plot particles of certain types
  iplotpartoftype(2:maxparttypes) = .false.
  imarktype = 1              ! PGPLOT marker for all particles
  imarktype(2) = 4           ! PGPLOT marker for ghost/dark matter particles
  imarktype(3) = 17          ! PGPLOT marker for sink particles 
  labeltype(1) = 'gas'
  labeltype(2) = 'type 2'
  labeltype(3) = 'type 3'
  labeltype(4) = 'type 4'
  labeltype(5) = 'type 5'
  labeltype(6) = 'type 6'
  
!
!--render options
!
  icolours = 1               ! colour scheme to use
  npix = 100                 ! pixels in x direction for rendering
  iPlotColourBar = .true.! whether or not to plot the colour bar
  iplotcont_nomulti = .false. ! plot contours
  ncontours = 30             ! number of contours to plot
!
!--cross section/rotation options
!  
  xsec_nomulti = .false.    ! take cross section of data / particles
  xsecpos_nomulti = 0.      ! position of cross section
  flythru = .false.         ! take series of cross sections through data
  xseclineX1 = 0.0
  xseclineX2 = 0.0
  xseclineY1 = 0.0
  xseclineY2 = 0.0
  irotate = .false.
  anglex = 0.
  angley = 0.
  anglez = 0.
  xorigin = 0.
!
!--vector plot options
!
  npixvec = 40        ! pixels in x direction on vector plots
  UseBackgndColorVecplot = .false. ! plot vector plot using black/white
  iamvec(:) = 0
  labelvec = ' '
  iVecplotLegend = .true.
  hposlegendvec = 0.1
  vposlegendvec = -1.0
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
  call defaults_set_exact
!
!--multiplot
!
  nyplotmulti = 4           ! number of plots in multiplot
  multiploty(:) = 0
  do i=1,4
     multiploty(i) = ndim+i  ! first plot : y axis
  enddo
  multiplotx(:) = 1          ! first plot : x axis
  irendermulti(:) = 0        ! rendering
  ivecplotmulti(:) = 0       ! vector plot
  x_secmulti(:) = .false.    ! take cross section?
  xsecposmulti(:) = 0.0      ! position of cross section
  iplotcontmulti(:) = .false.
  !
  !--array positions of specific quantities
  !
  ix = 0
  !do i=1,ndim
  !   ix(i) = i       ! coords
  !enddo
  ivx = 0      ! vx
  irho = 0     ! density
  ipr = 0      ! pressure
  iutherm = 0  ! thermal energy
  ih = 0       ! smoothing length
  ipmass = 0   ! particle mass
  ipr = 0      ! pressure
  irad = 0     ! radius
  !
  !--power spectrum options
  !      
  idisordered = .false.
  ipowerspecy = ndim+1
  wavelengthmax = 1.0
  nfreqspec = 32
  !
  !--filenames
  !
  rootname = ' '
  !
  !--limits
  !
  lim(:,:) = 0.
  itrans(:) = 0
  !
  !--data array sizes
  !
  maxpart = 0
  maxcol = 0
  maxstep = 0
  return    
end subroutine defaults_set
