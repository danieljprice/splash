program splash
!---------------------------------------------------------------------------------
!
!     SPLASH - a plotting utility for SPH data in 1, 2 and 3 dimensions
!     Copyright (C) 2005-2007 Daniel Price 
!     dprice@astro.ex.ac.uk
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
  use filenames, only:rootname,nfiles,maxfile,defaultsfile,limitsfile,animfile
  use getdata, only:get_data
  use defaults, only:defaults_set_initial,defaults_set,defaults_read
  use limits, only:read_limits
  use mainmenu, only:menu
  use mem_allocation, only:deallocate_all
  use projections3D, only:setup_integratedkernel
  use settings_data, only:buffer_data,lowmemorymode
  use settings_xsecrot, only:read_animfile
  use system_commands, only:get_number_arguments,get_argument
  use system_utils, only:lenvironment
  use asciiutils, only:read_asciifile
  implicit none
  integer :: i,ierr,nargs
  logical :: ihavereadfilenames,evsplash
  character(len=120) :: string
  character(len=*), parameter :: version = 'v1.9.1+ [Sep ''07]'

  !
  ! initialise some basic code variables
  !
  call defaults_set_initial
  
  !
  !  default names for defaults file and limits file
  !
  defaultsfile = 'splash.defaults'
  limitsfile = 'splash.limits'
  animfile = 'splash.anim'
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
  do while (i < nargs)
     i = i + 1
     call get_argument(i,string)

     if (string(1:1).eq.'-') then
        select case(trim(string(2:)))
        case('l')
           i = i + 1
           call get_argument(i,limitsfile)
        case('f')
           i = i + 1
           call get_argument(i,defaultsfile)
        case('e','ev')
           evsplash = .true.
           defaultsfile = 'evsplash.defaults'
           limitsfile = 'evsplash.limits'
           animfile = 'evsplash.anim'
        case('lowmem','lm')
           lowmemorymode = .true.
        case('nolowmem','nlm')
           lowmemorymode = .false.
        case default
           print "(a)",'SPLASH: a visualisation tool for Smoothed Particle Hydrodynamics simulations'
           print "(a,/)",trim(version)
           if (string(2:2).ne.'v') print "(a)",'unknown command line argument '''//trim(string)//''''
           print "(a)",'Usage: splash [-f defaultsfile] [-l limitsfile] [-ev] [-lowmem] file1 file2 ...'
           stop
        end select
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
     call read_asciifile('splash.filenames',nfiles,rootname)
     print*,nfiles,' filenames read from splash.filenames file'
     print*
     if (nfiles.gt.0) ihavereadfilenames = .true.
  endif
  if (lowmemorymode) print "(a)",' << running in low memory mode >>'

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
  
  !
  ! enter main menu
  !
  call menu
  
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

 print "(a)",'  ( '//trim(version)//' Copyright (C) 2005-2007 )'
 print 30 
30 format(/,    &
   ' * SPLASH comes with ABSOLUTELY NO WARRANTY.',/, &
   '   This is free software; and you are welcome to redistribute it ',/, &
   '   under certain conditions (see LICENSE file for details). *',/,/, &
   ' Comments, bugs, suggestions and queries to: dprice@astro.ex.ac.uk ',/, &
   ' Check for updates at: www.astro.ex.ac.uk/people/dprice/splash ',/, &
   ' Please cite Price (2007, PASA accepted, arXiv:0709.0832) if you ',/, &
   ' use SPLASH for scientific work and if you plot something beautiful,',/, &
   ' why not send me a copy for the gallery? ',/)
      
end subroutine print_header

end program splash
