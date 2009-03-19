program splash
!---------------------------------------------------------------------------------
!
!     SPLASH - a plotting utility for SPH data in 1, 2 and 3 dimensions
!     Copyright (C) 2005-2009 Daniel Price 
!     daniel.price@sci.monash.edu.au
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
!     1.12.1 : (xx/03/09)
!             Can edit/delete text shapes interactively, also the colour bar label; can customise
!             the label on projection plots; contour levels better defined; SPLASH_HMIN_CODEUNITS added;
!             option for numeric labelling of contours; minor bug fixes.
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
!     Uses PGPLOT routines to plot graphs and utilises the rendering
!     tools to plot density renderings and vector plots of 2D and 3D data.
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
!     See the CVS logs for a full ChangeLog
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
  use filenames, only:rootname,nfiles,maxfile,defaultsfile,limitsfile,animfile, &
                      fileprefix,set_filenames
  use getdata, only:get_data
  use defaults, only:defaults_set_initial,defaults_set,defaults_read
  use limits, only:read_limits
  use mainmenu, only:menu,allowrendering,set_coordlabels,set_extracols
  use mem_allocation, only:deallocate_all
  use projections3D, only:setup_integratedkernel
  use settings_data, only:buffer_data,lowmemorymode,ndim,ncolumns,ncalc,nextra,numplot,ndataplots
  use settings_xsecrot, only:read_animfile
  use system_commands, only:get_number_arguments,get_argument
  use system_utils, only:lenvironment
  use asciiutils, only:read_asciifile
  use write_pixmap, only:isoutputformat,iwritepixmap,pixmapformat
  use convert, only:convert_all
  use write_sphdata, only:issphformat
  use analysis, only:isanalysis
  use timestepping, only:timestep_loop
  use settings_page, only:interactive,device,nomenu
  implicit none
  integer :: i,ierr,nargs,ipickx,ipicky,irender,icontour,ivecplot
  logical :: ihavereadfilenames,evsplash,doconvert
  character(len=120) :: string
  character(len=12) :: convertformat
  character(len=*), parameter :: version = 'v1.12.1beta [11th March ''09]'

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
  lowmemorymode = lenvironment('SPLASH_LOW_MEM')
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
  doconvert = .false.
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
           if (ierr /= 0) call print_usage(quit=.true.)
           nomenu = .true.
        case('y')
           i = i + 1
           call get_argument(i,string)
           read(string,*,iostat=ierr) ipicky
           if (ierr /= 0) call print_usage(quit=.true.)
           nomenu = .true.
        case('render','r')
           i = i + 1
           call get_argument(i,string)
           read(string,*,iostat=ierr) irender
           if (ierr /= 0) call print_usage(quit=.true.)
           nomenu = .true.
        case('contour','c','cont')
           i = i + 1
           call get_argument(i,string)
           read(string,*,iostat=ierr) icontour
           if (ierr /= 0) call print_usage(quit=.true.)
        case('vec','vecplot')
           i = i + 1
           call get_argument(i,string)
           read(string,*,iostat=ierr) ivecplot
           if (ierr /= 0) call print_usage(quit=.true.)
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
        case('o')
           i = i + 1
           call get_argument(i,string)
           if (isoutputformat(string)) then
              iwritepixmap = .true.
              pixmapformat = trim(string)
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
     elseif (trim(string).eq.'to') then
     !
     !--for converting SPH formats
     !
           i = i + 1
           call get_argument(i,string)
           if (issphformat(string)) then
              doconvert = .true.
              convertformat = trim(string)
           else
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
     print "(a)",' no filenames read from command line'
     call read_asciifile(trim(fileprefix)//'.filenames',nfiles,rootname)
     print*,nfiles,' filenames read from '//trim(fileprefix)//'.filenames file'
     print*
     if (nfiles.gt.0) then
        ihavereadfilenames = .true.
     else
        print "(a/,/,5x,a,/)",' Basic usage: ','splash dumpfile(s)'
        print "(a/,/,5x,a,/)",' e.g.: ','gsplash snap_0*'
        print "(a)",' Or write the filenames one per line in a file called ''splash.filenames'''
        print "(a)",' For a full list of command-line options, use splash --help'
        print "(a,/)",' For help on basic splash usage, consult the userguide: splash/docs/splash.pdf'
        stop
     endif
  endif
  if (lowmemorymode) print "(a)",' << running in low memory mode >>'


  if (doconvert) then

     !
     ! batch convert all dump files into the output format
     !
     call convert_all(convertformat,ihavereadfilenames)
  
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

     !
     ! read animation file if it exists
     !
     call read_animfile(animfile)

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
   "    _                                                 _  ",/, &
   "   (_)   _               _           _         _     (_)_",/, &
   "      _ (_)    ___ _ __ | | __ _ ___| |__     (_)   _  (_)",/, &
   "   _ (_)  _   / __| '_ \| |/ _` / __| '_ \       _ (_)    ",/, &
   "  (_)  _ (_)  \__ \ |_) | | (_| \__ \ | | |  _  (_) _    ",/, &
   "      (_)  _  |___/ .__/|_|\__,_|___/_| |_| (_)  _ (_)   ",/, &
   "          (_)  (_)|_| (_) (_)  (_)(_) (_)(_) (_)(_)     ")        
 print 20
20 format(/,  &
   '  ( B | y ) ( D | a | n | i | e | l ) ( P | r | i | c | e )',/)

 print "(a)",'  ( '//trim(version)//' Copyright (C) 2005-2009 )'
 print 30 
30 format(/,    &
   ' * SPLASH comes with ABSOLUTELY NO WARRANTY.',/, &
   '   This is free software; and you are welcome to redistribute it ',/, &
   '   under certain conditions (see LICENSE file for details). *',/,/, &
   ' Comments, bugs, suggestions and queries to: daniel.price@sci.monash.edu.au ',/, &
   ' Check for updates at: http://users.monash.edu.au/~dprice/splash ',/, &
   ' Please cite Price (2007), PASA, 24, 159-173 (arXiv:0709.0832) if you ',/, &
   ' use SPLASH for scientific work and if you plot something beautiful,',/, &
   ' why not send me a copy for the gallery? ',/)
      
end subroutine print_header

subroutine print_usage(quit)
 implicit none
 logical, intent(in), optional :: quit
 logical :: ltemp

 print "(a)",'SPLASH: a visualisation tool for Smoothed Particle Hydrodynamics simulations'
 print "(a)",'(c) 2005-2009 Daniel Price '
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
 print "(a)",' -dev device       : specify PGPLOT device on command line (ie. do not prompt)'
 print "(a)"
 ltemp = issphformat('none')
 print "(a)"
 ltemp = isanalysis('none')

 if (present(quit)) then
    if (quit) stop
 endif

end subroutine print_usage

end program splash
