program supersphplot
!---------------------------------------------------------------------------------
!     plotting utility for SPH data in 1, 2 and 3 dimensions.
!
!     uses PGPLOT routines to plot graphs and utilises the rendering
!     tools to plot density renderings and vector plots of 2D and 3D data.
!
!     subroutines as follows (in alphabetical order):
!
!     allocate           : allocates memory for main arrays
!     calc_quantities    : calculates additional quantities from particle data
!     colour_demo        : demonstration of colour schemes for rendering
!     colour_set	 : sets up pgplot colour table for rendering
!     coord_transform    : transforms between various coord systems
!     danpgsch           : sets character height independent of page size
!     danpgtile          : my utility for tiling plots on the pgplot page
!     danpgwedg          : my very minor modification of pgwedg
!     defaults_read	 : read default plot options from file
!     defaults_set	 : sets default plot options if not read from file
!     defaults_write	 : write default plot options to file
!     exact_fromfile     : reads an exact solution tabulated in a file
!     exact_mhdshock     : some tabulated solutions for mhd shocks 
!     exact_polytrope    : exact solution for a polytrope
!     exact_rhoh	 : plots exact relation between density and smoothing length
!     exact_sedov        : exact solution for sedov blast wave
!     exact_shock        : exact solution for hydrodynamic shocks
!     exact_swave        : exact solution for a linear sound wave
!     exact_toystar      : exact solution for the toy star problem
!     get_data           : wrapper for main data read
!     interactive_part   : interactive utilities for particle plots
!     interpolate1D	 : interpolation of 1D sph data to 1D grid using sph kernel
!     interpolate_vec    : interpolation of vector data to 2D grid by simple averaging
!     interpolate2D	 : interpolation of 2D sph data to 2D grid using sph kernel     
!     interpolate2D_xsec : oblique 1D cross section through 2D sph data using kernel
!     interpolate3D	 : interpolation of 3D sph data to 3D grid using sph kernel
!     interpolate3D_fastxsec   : fast cross section through 3D data using sph kernel
!     interpolate3D_projection : fast projection of 3D data to 2D grid using integrated sph kernel
!     interpolate3D_xsec_vec   : fast cross section of vector quantity in 3D data using kernel
!     legend		       : plots legend on plot (time)
!     limits_read              : reads plot limits from file
!     limits_save              : saves plot limits to file
!     limits_set               : calculates plot limits
!     main               : main plotting loop
!     menu               : main menu
!     modules		 : contains all shared (global) variables
!     options_data       : sets options relating to current data
!     options_exact	 : sets options and params for exact solution calculation/plotting
!     options_limits     : sets options relating to plot limits
!     options_page       : sets options relating to page setup
!     options_particleplots : sets options relating to particle plots
!     options_powerspec  : sets options for power spectrum plotting
!     options_render	 : sets options for render plots
!     options_vector	 : sets options for vector plots
!     plot_average	 : bins particles along x-axis and plots average line
!     plot_kernel_gr     : plots the kernel shape in non-cartesian co-ordinates
!     plot_powerspectrum : calls powerspectrum and plots it
!     powerspectrum_fourier : calculates power spectrum of 1D data on ordered pts
!     powerspectrum_lomb : calculates power spectrum of 1D data on disordered pts
!     read_data_dansph   : reads data from my format of data files
!     read_data_mrbsph   : reads data from matthew bate's format of data files
!     render	 	 : takes array of pixels and plots render map/contours etc
!     setpage            : sets up the PGPLOT page (replaces call to PGENV/PGLAB)
!     supersphplot	 : main program, drives menu loop
!     titles_read        : reads a list of titles to be used to label each timestep
!     transform	 	 : applies various transformations to data (log10, 1/x, etc)
!
!     file format is specified in the subroutine read_data   
!
!     written by: Daniel Price, Institute of Astronomy, Cambridge UK
!          email: dprice@ast.cam.ac.uk
!
!     this version for both ndspmhd and matthew bate's code 2003-2004
!     summary of major changes: (for a full changelog see the CVS log - or use cvs2cl)
!
!      20/08/04 - vectorplot replaced by interpolate_vec
!      19/08/04 - azimuthal rotation works, interactive limits not permanent,
!                 dat restructured, various clean ups.
!      27/07/04 - 2D cross sections work, options_xsecrotate added
!      14/07/04 - major revamp of render/vector options + defaults save
!      19/06/04 - can transform particle coords to new coordinate systems
!      10/06/04 - exact solution for shock tubes, also from file + added read_exactparams
!      02/06/04 - interactive plotting steps forward/backwards, replots etc
!      01/06/04 - saves/reads limits to/from limits file
!                 also revamped menu - uses characters for options
!      31/05/04 - particle tracking limits, reads from Matthew's code       
!      17/05/04 - reads multiple files consecutively, error catches in interpolate 
!                 also tiling of plots using danpgtile
!      21/04/04 - page setup moved out of main -> setpage, danpgtile added
!      26/03/04 - options split into submenus
!      04/03/04 - allocatable arrays 
!                 (last non-allocatable version tagged as noalloc_04_03_04)
!      23/02/04 - lots of compiler bugs fixed. 2D->1D cross section
!      19/12/03 - separate subroutine main
!      17/12/03 - 1D interpolation and crap power spectrum
!		- some options moved to separate subroutines
!      16/12/03 - labels on particle cross sections
!      15/12/03 - namelist input/output, freeform source in modules
! 		- bug fix in read_data (nghosts) and interpolation routines
!      09/12/03 - power spectrum plotting in 1D
!      24/11/03 - calc_quantities in separate subroutine, rhoh moved
!      28/10/03 - bug fix when no data 
!      14/10/03 - new colour schemes
!      09/10/03 - gamma different betw. timesteps, paper size, menuitems
!      23/09/03 - plot error bars, fast cross-section, multi-transform
!      18/09/03 - exact mhd shocks from papers
!      16/09/03 - sink particle plotting, fixed bug in render_coarse
!      12/09/03 - fast column density plots, several bug fixes
!      11/09/03 - bug in vector plots (xminrender,xmaxrender)
!      11/08/03 - prompting, bug fix in resetting options after multiplot
!      30/07/03 - 3D cross sections/projections
!      29/07/03 - split into more subroutines (menu etc)
!      18/07/03 - transformations (log, 1/x etc)
!      15/07/03 - interpolate2D,3D - much simpler than smooth_pixels
!      25/06/03 - clever makefile - makes for dansph or m. bate sph
!               - subroutines in different files
!      20/06/03 - rendering can handle zero density, prints if min=max
!      19/06/03 - multiplot with rendering, no x array
!      16/06/03 - labels in module, specified in read_data
!      09/06/03 - ndim, ndimv changeable, reads ncolumns in data file
!      26/05/03 - calculated quantities no longer in read_data
!      16/05/03 - colour schemes + rendering using sph summation
!                 output format has changed, reads pmass array also
!      13/05/03 - read polytrope in menu option
!      01/05/03 - can manually enter plot limits 
!      29/04/03 - bug in initial plot limits fixed
!
!
!      plots can be of two types: co-ordinate plots or not
!      1) co-ordinate plots have co-ordinates as x and y axis
!         these plots can be rendered with any scalar or vector array.
!         
!         the rendering routines interpolate from the particles to either
!         a 2D or 3D grid. in 3D you can either render to a 3D grid and take
!         cross sections, or render to a 2D grid using a table of the integrated
!         sph kernel. this 2D rendering results in a map of the quantity
!         integrated through the third co-ordinate. 
!         rendering to a 3D grid can be quite slow - it is only efficient
!         if many cross sections are taken all at once from the same data.
!
!      2) other plots have a variety of options, with lines joining the particles
!         and various exact solutions. plot limits can be fixed or adaptive.
!
!      multiplot enables you to set up multiple plots per page, mixing from any type.
!
!----------------------------------------------------------------------------------
  use filenames
  use labels
  use settings
  implicit none
  logical :: iquit
  integer :: i,iprev
  !
  !--print header
  !
  call print_header
  print*,'( version 5.0 )'
  !
  !--set default options
  !
  call defaults_set
  !
  !--initialise variables
  !      
  rootname = 'blank'
  ishowopts = .false.  ! do/don't initially display menu options
  ivegotdata = .false.
  
  ! ---------------------------------------------
  ! read default options from file if it exists
  !
  call defaults_read

  ! ---------------------------------------------------
  ! get filenames from command line/file and read file
  
  iprev = 1
  i = 1
  do while (rootname(iprev)(1:1).ne.' ' .and. i.le.maxfile)
     call getarg(i,rootname(i))
     !!print*,i,rootname(i)
     iprev = i
     i = i + 1
     if (i.gt.maxfile .and. rootname(iprev)(1:1).ne.' ') then
        print*,'WARNING: number of files >= array size: setting nfiles = ',maxfile
     endif
  enddo
  if (i.gt.maxfile .and. rootname(maxfile)(1:1).ne.' ') then
     nfiles = maxfile
  else
     nfiles = iprev - 1
  endif
  print*,'number of files = ',nfiles

  if (rootname(1)(1:1).ne.' ') then
     ihavereadfilename = .true.
     call get_data
  else
     call get_data
  endif
  !------------------------------------------------------------
  ! setup kernel table for fast column density plots in 3D
  call setup_integratedkernel
  
  ! ----------------------------------------------------------------
  ! menu - loop back to here once finished plotting/setting options
  !
  menuloop: do while (.not.iquit)
     
     !----------------------------------------------------------------------
     !  print menu
     !
     call menu(iquit)
     
  enddo menuloop
                
end program supersphplot
