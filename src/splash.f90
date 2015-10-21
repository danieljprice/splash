!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2015 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!  The plotting API for SPLASH 2.0 was written by James Wetter
!  wetter.j@gmail.com
!
!-----------------------------------------------------------------

program splash
!---------------------------------------------------------------------------------
!
!     SPLASH - a plotting utility for SPH data in 1, 2 and 3 dimensions
!     Copyright (C) 2005-2015 Daniel Price
!     daniel.price@monash.edu
!
!     --------------------------------------------------------------------------
!
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!     -------------------------------------------------------------------------
!     Version history/ Changelog:
!     2.6.0  : (22/10/15)
!             SILO, falcON and .pbob data reads implemented; bug fixes in gadget-hdf5 reader;
!             can recognise particle types in ascii read; more robust sphNG read;
!             dust fraction recognised in phantom data read; Toomre Q works in physical units;
!             bug fix with disappearing units labels; bug fix in shock tube exact solution;
!             added splash calc delta; splash to ascii keeps precision; better power spectra
!     2.5.1  : (29/01/15)
!             error bar style options; support for 5K displays; can plot vectors
!             and render with colours if h not read; range restrictions apply during splash to grid;
!             improved line-style legend; now up to 6 line styles; fixes to amuse-hdf5 read; 
!             phantom read handles star/dm particles; various bugs fixed
!     2.5.0  : (22/08/14)
!             instant multiplots by giving multiple columns as y axis;
!             ability to plot multiple exact solution files on same plot;
!             compiles in parallel by default; support for tagged sphNG/Phantom format;
!             AMUSE hdf5 format reader added; various bug fixes
!     2.4.1  : (01/04/14)
!             Roche-lobe plotting vastly improved; newunit= issue fixed;
!             bug fix with reading sink velocities from Phantom; other minor bug fixes.
!     2.4.0  : (21/02/14)
!             time formatting in legend can include general functions like %(t + 1000);
!             option to include sinks in opacity rendering;
!             supports one-fluid dust visualisation;
!             C-shock exact solution; better polytrope solution
!     2.3.1  : (11/11/13)
!             SPLASH_COROTATE option to plot in frame corotating with sinks;
!             bug fixes with handling of dead/accreted/boundary particles in sphNG/phantom;
!             various other bugs fixed.
!     2.3.0  : (09/08/13)
!             can customise time formatting in legend; improvements to legends;
!             less verboseness; splash can read and plot pixel maps produced with -o ascii;
!             3D vector field plotting improved; bug fix with gfortran 4.8
!     2.2.2  : (10/05/13)
!             particle tracking by type implemented;
!             can interpolate specific columns in splash to grid;
!             SPLASH_CENTRE_ON_SINK option generic to all data reads;
!             Aly Reheam format added; option for 2nd y axis on plots;
!             bug fix with X11 linking on Ubuntu; can read gadget ICs files
!     2.2.1  : (21/02/13)
!             minor bug with axes plotting fixed;
!             Wendland kernels added; bugs with exact solution plotting fixed;
!             bug fix with tracking of dark matter particles
!     2.2.0  : (16/11/12)
!             option to use different kernels for interpolation;
!             floating/inset colour bars added;
!             splash to gadget conversion implemented;
!             splash to grid works in 2D;
!             improved interfaces to shapes and animation sequences
!             automatically turns on dark matter particle plotting if no gas
!             interactive mode help displayed automatically
!     2.1.1  : (31/08/12)
!             irregular/circular particle selection using shift-left/middle click;
!             improved h5part and GADGET HDF5 data reads;
!             splash can be compiled in double precision;
!             bug fixes with calculated quantities + change of coordinate systems;
!             improved vector plot legend; option for box+numbers but no labels added
!     2.1.0  : (16/05/12)
!             3D vector field visualisation added; 
!             GADGET HDF5 read implemented;
!             page sizes can be specified in pixels;
!             limits can auto-adapt to device aspect ratio;
!             more general exact solution from file option;
!             tiling works with one colour bar per row;
!             splash calc handles different particle types
!       2.0  : (29/08/11)
!             new giza backend - antialiased lines; real fonts; pdf, eps and svg drivers;
!             fewer build dependencies (only cairo, X11);
!             support for semi-transparent text;
!             Double rendering (with transparent background) implemented.
!     1.15.0 : (29/08/11)
!             Multiplot with different particle types implemented; calculated quantities
!             list is now pre-filled automatically; preliminary support for r-phi and r-z
!             rendering; outlined solid markers implemented; better handling of multiple types;
!             manual contour levels can be specified in splash.contours; parallel splash to grid;
!             better support for non-square pixels; clipping of numbers at edge of viewport fixed
!     1.14.1 : (17/03/11)
!             SEREN data read added; dragon read updated; build follows Gnu conventions
!             on DEST and DESTDIR (needed for macports build); can have up to 12 particle types;
!             exact solutions re-ordered; dusty wave exact solution added
!     1.14.0 : (06/12/10)
!             Can flip between rendered quantities in interactive mode using 'f/F';
!             SPLASH_DEFAULTS variable can be set for system-wide defaults;
!             can plot arbitrary functions of x,t as exact solution; asplash better
!             handles blank lines in header and can specify time, gamma location with
!             env. variables; added data read for the H5PART format; GADGET read
!             across multiple files implemented; VINE read works with particle injection;
!             error bars can be plotted for both x and y axis simultaneously;
!             default rotation angles are set if 3D perspective turned on;
!             new directory layout and more helpful error messages during build;
!             PGPLOT linking is easier to get right.
!     1.13.1 : (26/02/10)
!             bugs with new calc_quantities module fixed; generic library interface
!             implemented so backend can be changed easily; bug fix with auto pixel selection;
!             simpler foreground/background colour setting; added subgrid interpolation warning
!     1.13.0 : (25/02/10)
!             function parser incorporated; calculated quantities can now be specified
!             at runtime, arbitrary function plotting implemented as an exact
!             solution; command-line SPH->grid conversion ("splash to grid")
!             implemented; ctrl-t in interactive mode adds arbitrary text box;
!             better line style/colour changing; bug fix with tiling and y-axis labels;
!             various other bug fixes.
!     1.12.2 : (15/07/09)
!             Variable marker sizes added, can plot particles as circles with
!             size proportional to h; dark matter rendering with block-labelled
!             GADGET format fixed; VINE read handles star particles; TIPSY read
!             with ifort10.0.0 works; snsph read added; splash to phantom added;
!             does not override labels for coords, vectors by default; bug fixes
!             with contouring options; stability bug fixes with older compilers;
!             more robust memory handling; bug fix with automatic pixel selection
!             causing seg fault.
!     1.12.1 : (20/04/09)
!             Can edit/delete text shapes interactively, also the colour bar label; can customise
!             the label on projection plots; contour levels better defined; SPLASH_HMIN_CODEUNITS added;
!             option for numeric labelling of contours; contour limits can be set separately
!             to render limits for same quantity; minor bug fixes.
!     1.12.0 : (22/12/08)
!             Command-line plotting implemented; ln transform added; bug fixes in GADGET read;
!             Backspace over annotation (legends,titles,axes,colour bar) in interactive mode
!             removes it; "splash calc" command line utility calculates time sequences of
!             global quantities from a sequence of dump files; bug fix causing seg fault.
!     1.11.1 : (13/10/08)
!             automatic number of pixels and exact pixel boundaries implemented;
!             mass does not have to be read from dump file; frame changes are per-page
!             not per-dump file for animation sequences; lower stacksize footprint;
!             bug fix with circles of interaction; bug fixes with block-labelled GADGET read;
!             Steve Foulkes data read added.
!     1.11.0 : (15/08/08)
!             ability to use subset of particles in restricted parameter range(s);
!             probability density function plot option; plot-hugging colour bars added;
!             ability to annotate plot with a range of shapes; v,V,w and H implemented
!             in interactive mode for >1 panel; various bug fixes (including one with vphi).
!     1.10.2 : (08/05/08)
!             disc surface density / toomre q parameter plotting added; flash colour
!             schemes added; splash to binary convert option, can change order in
!             which particle types are plotted; splash.columns file overrides
!             column label settings; vanaverbeke format read; various bug fixes.
!     1.10.1 : (11/03/08)
!             "splash to" command line option converts binary dumps to ascii format;
!             vector plots + rotation now implemented; block labelled GADGET format read;
!             ring-spreading exact solution added.
!     1.10.0 : (28/11/07)
!             horizontal colour bars implemented; -p, -o command line options;
!             can have mixed types in data reads; TIPSY and DRAGON data reads;
!             density weighted rendering; normalisation applies to column
!             density plots; improved particle tracking; save as option; various bug fixes
!     1.9.2 : (12/09/07)
!             improvements to ascii read including asplash -e option;
!             smarter foreground/background colour changing for titles;
!             min=max problem fixed (caught by splash not pgplot);
!             fixed vector arrow length option; other minor changes and bug fixes
!     1.9.1 : (11/07/07)
!             environment variables + improvements to gadget data read; better
!             prompting; 3 new colour schemes; improved legend/title options;
!             other minor changes
!     1.9.0 : (21/05/07)
!             animation sequences implemented; origin settings now affect radius
!             calculation and are relative to tracked particle; automatic line
!             width choice for postscript devices; w key adapts vector arrows;
!             vastly improved userguide
!     1.8.1 : (28/03/07)
!             option to hide vector arrows where there are no particles added;
!             smoother 3D plotting at low pixel numbers;
!             (smoother vector plots); bug fixes with a); issues with
!             round-off error with z integration of vectors fixed.
!     1.8.0 : (14/03/07)
!             hidden particles not used in rendering; units for z integration added;
!             a) & g) implemented in interactive mode for multiple-plots-per-page;
!             improved cross section using x in interactive mode
!     1.7.2 : (19/02/07)
!             Menu shortcuts implemented; bug fix/ more sensible transformation
!             of angular vector components in different co-ordinate systems;
!             improvements to interactive zoom and origin recentreing;
!             improved colour-by-type option; restrictions on page size removed;
!             minor bug fixes
!     1.7.1 : (04/01/07)
!             command line options for defaults and limits files added;
!             minor bug fixes
!     1.7.0 : (13/12/06)
!             renamed SPLASH instead of SUPERSPHPLOT; much faster data read
!             for gadget and sphNG reads (only required columns read);
!             physical units can be saved to file; new menu formats; various
!             other bug fixes.
!     1.6.2   (24/10/06)
!           : fast particle plotting and streamline plotting implemented;
!             more bug fixes with interactive mode on multiplots; various other bug fixes.
!     1.6.1   (24/8/06)
!           : bug fixes to 1.6.0, further improvements to interactive mode on multiplots.
!     1.6.0   (10/8/06)
!           : Interactive mode on multiple plots per page; highly optimised interpolation
!             + parallel version; new Makefile; various bug fixes
!     1.5.4 (06/7/06)
!           : Handles multiple SPH/non-SPH particle types; axes redrawn after rendering;
!             minor bug fixes
!     1.5.3 (27/6/06)
!           : minor bug fixes/improvements to multiple plots per page, colour bar labelling
!             tiled plots, legend. Accelerated rendering option for projections.
!     1.5.2 (11/5/06)
!           : S) option for saving limits and defaults; MUCH faster interactive
!             replotting (no unnecessary re-rendering), a few other minor things
!     1.5.1 (26/4/06)
!           : docs updated for v1.5, other minor changes
!     1.5.0 (17/3/06)
!           : 3D perspective added, 3D opacity rendering, improved rotation,
!             colour schemes, adjustable vector arrows (+legend), improved timestepping
!             behaviour, speed enhancements, physical unit rescaling
!     1.0.5 (28/9/05)
!           : error calculation for exact solutions, legend for plot markers,
!             exact_densityprofiles added, more colour schemes,
!             unit rescaling improved, other minor changes + bug fixes
!     1.0.4 (17/8/05)
!           : better colour schemes; interactive colour scheme changing;
!             various minor changes and bug fixes
!     1.0.3 (5/7/05)
!           : rescale data option; better page setup; improved zooming;
!             interactive particle tracking + various minor changes and bug fixes
!     1.0.2 : much improved ascii data read; better line plotting; zoom on
!             powerspectrum plots + various bug fixes
!     1.0.1 : bug fixes relating to colour bars on multiplots
!     1.0   : first "official" release: version given to many people at IPAM
!             meeting and put on web.
!
!     -------------------------------------------------------------------------
!
!     Modules/subroutines as follows (in alphabetical order):
!
!     allocate           : allocates memory for main arrays
!     calc_quantities    : calculates additional quantities from particle data
!     colours            : colour schemes for rendering
!     colourparts	 : colours particles
!     defaults           : writes/reads default options to/from file
!     exact              : module handling exact solution settings
!     exact_fromfile     : reads an exact solution tabulated in a file
!     exact_mhdshock     : some tabulated solutions for mhd shocks
!     exact_polytrope    : exact solution for a polytrope
!     exact_rhoh	 : exact relation between density and smoothing length
!     exact_sedov        : exact solution for sedov blast wave
!     exact_shock        : exact solution for hydrodynamic shocks
!     exact_wave         : exact solution for a propagating sine wave
!     exact_toystar1D    : exact solution for the 1D toy star problem
!     exact_toystar2D    : exact solution for the 2D toy star problem
!     fieldlines         : module handling streamline plotting
!     get_data           : wrapper for main data read
!     geometry           : module handling different coordinate systems
!     globaldata         : various modules containing "global" variables
!     interactive        : drives interactive mode
!     interpolate1D	 : interpolation of 1D SPH data to grid using kernel
!     interpolate2D	 : interpolation of 2D SPH data to grid
!     interpolate3D_xsec : 3D cross section interpolations
!     interpolate3D_projection	 : 3D interpolation integrated through domain
!     legends		       : plots (time) legend on plot
!     limits                   : sets initial plot limits and writes to/reads from limits file
!     menu               : main menu
!     options_data       : sets options relating to current data
!     options_limits     : sets options relating to plot limits
!     options_page       : sets options relating to page setup
!     options_particleplots : sets options relating to particle plots
!     options_powerspec  : sets options for power spectrum plotting
!     options_render	 : sets options for render plots
!     options_vector	 : sets options for vector plots
!     options_xsecrotate : sets options for cross sections and rotation
!     particleplot       : subroutines for particle plotting
!     plotstep           : main subroutines which drive plotting of a single timestep
!     powerspectrums     : calculates power spectrum of 1D data (2 methods)
!     read_data_dansph   : reads data from my format of data files
!     read_data_mbate    : reads data from matthew bate's format of data files
!     render	 	 : takes array of pixels and plots render map/contours etc
!     rotate             : subroutines controlling rotation of particles
!     setpage            : sets up the PGPLOT page (replaces call to PGENV/PGLAB)
!     splash             : main program, drives menu loop
!     timestepping       : controls stepping through timesteps
!     titles             : reads a list of titles to be used to label each timestep
!     transform	 	 : applies various transformations to data (log10, 1/x, etc)
!
!     File format is specified in the subroutine read_data
!
!     See the svn logs for a full ChangeLog
!
!      Plots can be of two types: co-ordinate plots or not
!
!      1) Co-ordinate plots have co-ordinates as x and y axis
!         these plots can be rendered with any scalar or vector array.
!
!         The rendering routines interpolate from the particles to either
!         a 2D or 3D grid. In 3D you can either render to a 3D grid and take
!         cross sections, or render to a 2D grid using a table of the integrated
!         SPH kernel. This 2D rendering results in a map of the quantity
!         integrated through the third co-ordinate.
!         Rendering to a full 3D grid can be quite slow - it is used only
!         if many cross sections are taken all at once from the same data.
!
!      2) other plots have a variety of options, with lines joining the particles
!         and various exact solutions. Plot limits can be fixed or adaptive.
!
!      multiplot enables you to set up multiple plots per page, mixing from any type.
!
!----------------------------------------------------------------------------------
  use filenames, only:rootname,nfiles,maxfile,defaultsfile,limitsfile, &
                      fileprefix,set_filenames
  use getdata,   only:get_data
  use geomutils, only:set_coordlabels
  use defaults,  only:defaults_set_initial,defaults_set,defaults_read
  use limits,    only:read_limits
  use kernels,   only:ikernel,select_kernel_by_name,select_kernel
  use mainmenu,  only:menu,allowrendering,set_extracols
  use mem_allocation,     only:deallocate_all
  use projections3D,      only:setup_integratedkernel
  use settings_data,      only:buffer_data,lowmemorymode,debugmode,ndim,ncolumns,ncalc,nextra,numplot,ndataplots
  use system_commands,    only:get_number_arguments,get_argument,get_environment
  use system_utils,       only:lenvironment
  use asciiutils,         only:read_asciifile,basename
  use write_pixmap,       only:isoutputformat,iwritepixmap,pixmapformat,isinputformat,ireadpixmap,readpixformat
  use convert,            only:convert_all
  use write_sphdata,      only:issphformat
  use readwrite_griddata, only:isgridformat,print_gridformats
  use analysis,           only:isanalysis
  use timestepping,       only:timestep_loop
  use settings_page,      only:interactive,device,nomenu
  implicit none
  integer :: i,ierr,nargs,ipickx,ipicky,irender,icontour,ivecplot
  logical :: ihavereadfilenames,evsplash,doconvert,useall,iexist
  character(len=120) :: string
  character(len=12)  :: convertformat
  character(len=*), parameter :: version = 'v2.6.0 [22nd Oct. 2015]'

  !
  ! initialise some basic code variables
  !
  call defaults_set_initial

  !
  !  default names for defaults file and limits file
  !
  fileprefix = 'splash'
  call set_filenames(trim(fileprefix))

  evsplash = .false.
  lowmemorymode = lenvironment('SPLASH_LOW_MEM') .or. lenvironment('SPLASH_LOWMEM')
  debugmode = lenvironment('SPLASH_DEBUG')
  !
  !  read all arguments off command line
  !
  call get_number_arguments(nargs)
  !
  !  extract command line arguments and filenames
  !
  i = 0
  nfiles = 0
  iwritepixmap = .false.
  ireadpixmap  = .false.
  doconvert = .false.
  useall = .false.
  nomenu = .false.
  ipickx = 0
  ipicky = 0
  irender = 0
  icontour = 0
  ivecplot = 0

  do while (i < nargs)
     i = i + 1
     call get_argument(i,string)

     if (string(1:1).eq.'-') then
        select case(trim(string(2:)))
        case('x')
           i = i + 1
           call get_argument(i,string)
           read(string,*,iostat=ierr) ipickx
           if (ierr /= 0 .or. ipickx <= 0) call print_usage(quit=.true.)
           nomenu = .true.
        case('y')
           i = i + 1
           call get_argument(i,string)
           read(string,*,iostat=ierr) ipicky
           if (ierr /= 0 .or. ipicky <= 0) call print_usage(quit=.true.)
           nomenu = .true.
        case('render','r','ren')
           i = i + 1
           call get_argument(i,string)
           read(string,*,iostat=ierr) irender
           if (ierr /= 0 .or. irender < 0) call print_usage(quit=.true.)
           nomenu = .true.
        case('contour','c','cont','con')
           i = i + 1
           call get_argument(i,string)
           read(string,*,iostat=ierr) icontour
           if (ierr /= 0 .or. icontour < 0) call print_usage(quit=.true.)
        case('vec','vecplot')
           i = i + 1
           call get_argument(i,string)
           read(string,*,iostat=ierr) ivecplot
           if (ierr /= 0 .or. ivecplot < 0) call print_usage(quit=.true.)
           nomenu = .true.
        case('dev','device')
           i = i + 1
           call get_argument(i,device)
        case('l')
           i = i + 1
           call get_argument(i,limitsfile)
        case('d','f')
           i = i + 1
           call get_argument(i,defaultsfile)
        case('p')
           i = i + 1
           call get_argument(i,string)
           if (len_trim(string).gt.0) then
              fileprefix = trim(string)
              call set_filenames(trim(fileprefix))
           endif
        case('o','writepix','wpix')
           i = i + 1
           call get_argument(i,string)
           if (isoutputformat(string)) then
              iwritepixmap = .true.
              pixmapformat = trim(string)
           else
              stop
           endif
        case('readpix','rpix')
           i = i + 1
           call get_argument(i,string)
           if (isinputformat(string)) then
              ireadpixmap = .true.
              readpixformat = trim(string)
           else
              stop
           endif
        case('e','ev')
           evsplash = .true.
           fileprefix = 'evsplash'
           call set_filenames(trim(fileprefix))
        case('lowmem','lm')
           lowmemorymode = .true.
        case('nolowmem','nlm')
           lowmemorymode = .false.
        case('-help')
           call print_usage
           print "(/,a)",' Basic splash usage is explained in the userguide,'
           print "(a,/)",'  located in the directory splash/docs/splash.pdf'
           stop
        case default
           call print_usage
           if (string(2:2).ne.'v') print "(a)",'unknown command line argument '''//trim(string)//''''
           stop
        end select
     elseif (trim(string).eq.'to' .or. trim(string).eq.'allto') then
     !
     !--for converting SPH formats
     !
           if (trim(string).eq.'allto') useall = .true.
           i = i + 1
           call get_argument(i,string)
           if (isgridformat(string)) then
              doconvert = .true.
              convertformat = trim(string)
           elseif (issphformat(string)) then
              doconvert = .true.
              convertformat = trim(string)
           else
              call print_gridformats()
              stop
           endif
     elseif (trim(string).eq.'calc') then
     !
     !--for performing analysis on a sequence of dump files
     !
           i = i + 1
           call get_argument(i,string)
           if (isanalysis(string)) then
              doconvert = .true.
              convertformat = trim(string)
           else
              stop
           endif
     elseif (len_trim(string).gt.0) then
        nfiles = nfiles + 1
        if (nfiles.le.maxfile) then
           rootname(nfiles) = trim(string)
        endif
     endif
  enddo

  !
  ! print header
  !
  call print_header
  !
  ! set default options (used if defaults file does not exist)
  !
  call defaults_set(evsplash)

  !
  ! read default options from file if it exists
  !
  call defaults_read(defaultsfile)
  !
  ! look for a system-wide defaults file if the environment
  ! variable SPLASH_DEFAULTS is set, no local file is present
  ! and no alternative prefix has been set.
  !
  inquire(file=defaultsfile,exist=iexist)
  if (.not.iexist .and. trim(fileprefix).eq.'splash') then
     call get_environment('SPLASH_DEFAULTS',string)
     if (len_trim(string).ne.0) then
        i = index(string,'.defaults')
        if (i.gt.0) then
           defaultsfile = trim(string)
        else
           defaultsfile = trim(string)//'.defaults'
        endif
        print "(a)",' Using SPLASH_DEFAULTS='//trim(defaultsfile)
        call defaults_read(defaultsfile)
        call set_filenames(trim(fileprefix))
     endif
  endif

  !
  ! check that we have got filenames
  !
  if (nfiles.gt.0) then
     if (nfiles.gt.maxfile) then
        print*,' WARNING: number of files >= array size: setting nfiles = ',maxfile
        nfiles = maxfile
     endif
  endif
  if (nfiles.ge.1 .and. rootname(1)(1:1).ne.' ') then
     ihavereadfilenames = .true.
     if (nfiles.gt.1) print*,nfiles,' filenames read from command line'
  else
     ihavereadfilenames = .false.
     !print "(a)",' no filenames read from command line'
     call read_asciifile(trim(fileprefix)//'.filenames',nfiles,rootname)
     !print*,nfiles,' filenames read from '//trim(fileprefix)//'.filenames file'
     if (nfiles.gt.0) then
        ihavereadfilenames = .true.
     else
        call get_argument(0,string)
        print "(/,a/,/,5x,a)",' Usage: ',trim(basename(string))//' snap_0*  (or use '&
                                //trim(fileprefix)//'.filenames to list files)'
        print "(5x,a,/)",trim(basename(string))//' --help   (for all command line options)'
        stop
     endif
  endif
  if (lowmemorymode) print "(a)",' << running in low memory mode >>'

  if (ikernel.eq.0) then
     !--if no kernel has been set
     call get_environment('SPLASH_KERNEL',string)
     if (len_trim(string).gt.0) then
        call select_kernel_by_name(string)
     else
        call select_kernel(0)
     endif
  else
     call select_kernel(ikernel)
  endif

  if (doconvert) then

     !
     ! batch convert all dump files into the output format
     !
     call convert_all(convertformat,ihavereadfilenames,useall)

  else
     !
     ! read data from file
     !
     if (buffer_data) then
        call get_data(-1,ihavereadfilenames)
     else
        call get_data(1,ihavereadfilenames,firsttime=.true.)
     endif

     !
     ! setup kernel table for fast column density plots in 3D
     !
     call setup_integratedkernel

     !
     ! read plot limits from file (overrides get_data limits settings)
     !
     call read_limits(trim(limitsfile),ierr)

     if (nomenu) then
     !
     !  initialise the things we would need if we called menu directly
     !
        call set_extracols(ncolumns,ncalc,nextra,numplot,ndataplots)
        call set_coordlabels(numplot)
        interactive = .false.
     !
     ! check command line plot invocation
     !

        if (ipicky.gt.0 .and. ipicky.le.numplot+1) then
           if (ipicky.le.numplot .and. (ipickx.eq.0 .or. ipickx.gt.numplot)) then
              print "(a)",' ERROR: x plot not set or out of bounds (use -x col)'
              stop
           endif
           if (irender.gt.0) then
              if (.not.allowrendering(ipicky,ipickx)) then
                 print "(a)",' ERROR: cannot render with x, y choice (must be coords)'
                 stop
              endif
              if (icontour.gt.numplot .or. icontour.lt.0) then
                 print "(a)",' ERROR: contour plot choice out of bounds'
                 stop
              endif
           elseif (icontour.gt.0) then
              print "(a)",' ERROR: -cont also requires -render setting'
              stop
           endif
        else
           if (irender.gt.0 .and. ndim.ge.2) then
              ipicky = 2
              ipickx = 1
              if (.not.allowrendering(ipicky,ipickx)) then
                 print "(a)",' ERROR: cannot render'
                 stop
              endif
              if (icontour.gt.numplot .or. icontour.lt.0) then
                 print "(a)",' ERROR: contour plot choice out of bounds'
                 stop
              endif
           else
              print "(a)",' ERROR: y plot not set or out of bounds (use -y col)'
              stop
           endif
        endif

        call timestep_loop(ipicky,ipickx,irender,icontour,ivecplot)
        !
        ! if we invoked an interactive device, enter the menu as usual, otherwise finish
        !
        if (interactive) call menu
     else
     !
     ! enter main menu
     !
        call menu
     endif
  endif

  !
  ! deallocate all memory (not strictly necessary)
  !
  call deallocate_all

contains

!------------------------------------------------------
! this subroutine prints the splash screen on startup
!------------------------------------------------------
subroutine print_header
 implicit none

 print 10
10 format( &
   "        _                _           _       _         _ ",/, &
   "      _(_)     ___ _ __ | | __ _ ___| |__   (_)     _ (_)",/, &
   "   _ (_)  _   / __| '_ \| |/ _` / __| '_ \      _  (_)   ",/, &
   "  (_)  _ (_)  \__ \ |_) | | (_| \__ \ | | |  _ (_)  _    ",/, &
   "      (_)  _  |___/ .__/|_|\__,_|___/_| |_| (_)  _ (_)   ",/, &
   "          (_)  (_)|_| (_) (_)  (_)(_) (_)(_) (_)(_)      ")
 print 20
20 format(/,  &
   '  ( B | y ) ( D | a | n | i | e | l ) ( P | r | i | c | e )',/)

 print "(a)",'  ( '//trim(version)//' Copyright (C) 2005-2015 )'
 print 30
30 format(/,    &
   ' * SPLASH comes with ABSOLUTELY NO WARRANTY.',/, &
   '   This is free software; and you are welcome to redistribute it ',/, &
   '   under certain conditions (see LICENCE file for details). *',/,/, &
   ' Updates/userguide: http://users.monash.edu.au/~dprice/splash ',/, &
   ' Email: daniel.price@monash.edu or splash-users@googlegroups.com',/, &
   ' Please cite Price (2007), PASA, 24, 159-173 (arXiv:0709.0832) if you ',/, &
   ' use SPLASH in print and don''t forget to send pics for the gallery.',/)

end subroutine print_header

subroutine print_usage(quit)
 use filenames, only:tagline
 implicit none
 logical, intent(in), optional :: quit
 logical :: ltemp

 print "(a)",trim(tagline)
 print "(a,/)",trim(version)
 print "(a,/)",'Usage: splash file1 file2 file3...'
 print "(a,/,a,/)",'Usage with flags: splash [-p fileprefix] [-d defaultsfile] [-l limitsfile] [-ev] ', &
               '[-lowmem] [-o format] [-x col] [-y col] [-render col] [-cont col] file1 file2 ...'

 print "(a,/)",'Command line options:'
 print "(a)",' -p fileprefix     : change prefix to ALL settings files read/written by splash '
 print "(a)",' -d defaultsfile   : change name of defaults file read/written by splash'
 print "(a)",' -l limitsfile     : change name of limits file read/written by splash'
 print "(a)",' -e, -ev           : use default options best suited to ascii evolution files (ie. energy vs time)'
 print "(a)",' -lm, -lowmem      : use low memory mode [applies only to sphNG data read at present]'
 print "(a)",' -o pixformat      : dump pixel map in specified format (use just -o for list of formats)'
 print "(/,a,/)",'Command line plotting mode:'
 print "(a)",' -x column         : specify x plot on command line (ie. do not prompt for x)'
 print "(a)",' -y column         : specify y plot on command line (ie. do not prompt for y)'
 print "(a)",' -r[ender] column  : specify rendered quantity on command line (ie. no render prompt)'
 print "(a)",'                     (will take columns 1 and 2 as x and y if -x and/or -y not specified)'
 print "(a)",' -vec[tor] column  : specify vector plot quantity on command line (ie. no vector prompt)'
 print "(a)",' -c[ontour] column : specify contoured quantity on command line (ie. no contour prompt)'
 print "(a)",' -dev device       : specify plotting device on command line (ie. do not prompt)'
 print "(a)"
 ltemp = issphformat('none')
 call print_gridformats()
 print "(a)"
 ltemp = isanalysis('none')

 if (present(quit)) then
    if (quit) stop
 endif

end subroutine print_usage

end program splash
