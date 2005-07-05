program supersphplot
!---------------------------------------------------------------------------------
!
!     SUPERSPHPLOT - a plotting utility for SPH data in 1, 2 and 3 dimensions
!     Copyright (C) 2005 Daniel Price 
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
!     danpgsch           : sets character height independent of page size
!     danpgtile          : my utility for tiling plots on the pgplot page
!     danpgwedg          : my very minor modification of pgwedg
!     defaults           : writes/reads default options to/from file
!     exact              : module handling exact solution settings
!     exact_fromfile     : reads an exact solution tabulated in a file
!     exact_mhdshock     : some tabulated solutions for mhd shocks 
!     exact_polytrope    : exact solution for a polytrope
!     exact_rhoh	 : exact relation between density and smoothing length
!     exact_sedov        : exact solution for sedov blast wave
!     exact_shock        : exact solution for hydrodynamic shocks
!     exact_wave         : exact solution for a propagating sine wave
!     exact_toystar      : exact solution for the toy star problem
!     exact_toystar2D    : exact solution for the 2D toy star problem
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
!     supersphplot	 : main program, drives menu loop
!     timestepping       : controls stepping through timesteps
!     titles_read        : reads a list of titles to be used to label each timestep
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
  use filenames, only:rootname,nfiles,maxfile
  use getdata, only:get_data
  use defaults, only:defaults_set,defaults_read
  use mainmenu, only:menu
  use mem_allocation, only:deallocate_all
  use projections3D, only:setup_integratedkernel
  use settings_data, only:buffer_data
  use system_commands
  implicit none
  integer :: i
  logical :: ihavereadfilenames

  !
  ! print header
  !
  call print_header

  !
  ! set default options
  !
  call defaults_set
  
  !
  ! read default options from file if it exists
  !
  call defaults_read

  !
  ! get filenames from command line if possible
  !
  call get_number_arguments(nfiles)
  if (nfiles.gt.0) then
     if (nfiles.gt.maxfile) then
        print*,' WARNING: number of files >= array size: setting nfiles = ',maxfile     
        nfiles = maxfile
     endif
     do i=1,nfiles
        call get_argument(i,rootname(i))
     enddo
  endif   
  if (nfiles.ge.1 .and. rootname(1)(1:1).ne.' ') then
     ihavereadfilenames = .true.
  else
     ihavereadfilenames = .false.
     print "(a/)",' no filenames read from command line'
  endif

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
10 format(/, &
   "                                     _           _       _   ",/, &
   " ___ _   _ _ __   ___ _ __ ___ _ __ | |__  _ __ | | ___ | |_ ",/, &
   "/ __| | | | '_ \ / _ \ '__/ __| '_ \| '_ \| '_ \| |/ _ \| __|",/, &
   "\__ \ |_| | |_) |  __/ |  \__ \ |_) | | | | |_) | | (_) | |_ ",/, &
   "|___/\__,_| .__/ \___|_|  |___/ .__/|_| |_| .__/|_|\___/ \__|",/, &
   "          |_|                 |_|         |_|                ")
 print 20
20 format(  &
   '    _   _     _   _   _   _   _   _     _   _   _   _   _  ',/, &
   '   / \ / \   / \ / \ / \ / \ / \ / \   / \ / \ / \ / \ / \ ',/, &
   '  ( B | y ) ( D | a | n | i | e | l ) ( P | r | i | c | e )',/, &
   '   \_/ \_/   \_/ \_/ \_/ \_/ \_/ \_/   \_/ \_/ \_/ \_/ \_/ ',/)      

 print "(a)",' ( version 1.0.2 [01/06/05] Copyright (C) 2005)'
 print 30 
30 format(/,    &
   ' * SUPERSPHPLOT comes with ABSOLUTELY NO WARRANTY.',/, &
   '   This is free software; and you are welcome to redistribute it ',/, &
   '   under certain conditions (see LICENSE file for details). *',/,/, &
   ' Comments, bugs, suggestions and queries to: dprice@astro.ex.ac.uk ',/, &
   ' Check for updates at: www.astro.ex.ac.uk/people/dprice/supersphplot ',/, &
   ' Credits are always nice (but not essential) - However, if you plot ',/, &
   ' something *really* nice, please send me a copy. ',/)
      
end subroutine print_header
             
end program supersphplot
