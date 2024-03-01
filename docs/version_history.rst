
**3.10.2: (01/03/24)**

- reads phantom dumps with adaptive particle refinement
- improved splash to phantom conversion
- better documentation of splash calc lightcurve (thanks to Chunliang Mu)
- saving dust density limits now applies to all dust species

**3.10.1: (4/12/23)**

- bug fix with accreted particles appearing in smoothed particle plots
- automated Trad/Tgas in extra quantities from phantom dumps with radiation

**3.10.0: (30/11/23)**

- --sort flag to sort filenames for comparison plots
- --movie flag to automatically make movie from sequence of files
- giza backend now supports direct output of mp4 movies

**3.9.0: (06/11/23)**

- follow-the-label column choice, where if a label is selected for plotting from the first file, it will automagically shift to find the matching label in subsequent files, even if the column containing the quantity has changed
- implemented vtk reader capable of reading snapshots from Shamrock code
- apply transparency to smoothed particle plot only if adaptive smoothing used
- plot colour bar by default when particle colouring by quantity is used
- can read velocity array from fits header for position-position-velocity cubes
- bug fix with first page being white with smoothed particle plots on black background
- bug fix finding .comp file if underscore in the directory name

**3.8.5: (23/10/23)**

- implemented smoothed particle plots with multiple steps per page
- allow for .cols and .comp file in current directory even if the filepath is not the current dir
- sphmoments utility added
- added routine to extract velocity dimension from fits files
- bug fix with repeated string replacements giving endless backslashes in labels
- bug fix reading long filenames in denoise

**3.8.4: (18/08/23)**

- various bugs fixed in GADGET data reader
- auto-recognise GADGET block format
- improved conservation checks in splash to grid (thanks to Avi Chen)
- better handling of AREPO data
- bug fix with timestep not advancing
- no longer ask about particle types in multiplot if only one type present

**3.8.3: (05/07/23)**

- flip option (f/F in interactive mode) now persists across timesteps and works in snapshots other than the first
- bugs fixed in Tipsy data read (thanks to Alex Pettitt)
- auto-recognition of Tipsy binary formats implemented
- show units labels in calculated quantities list

**3.8.2: (12/05/23)**

- phantom data read looks for .comp file containing additional composition data
- also looks for .cols file containing any extra columns with one row per particle
- recognise opacity if extra quantity called "kappa" calculated

**3.8.1: (01/05/23)**

- seg faults in auto-magic exact solution mapping fixed
- longer line limit in determining number of columns in ascii/exact solution files
- automatically handle log in exact solution labels (e.g. logR, logT)

**3.8.0: (26/04/23)**

- plots multiple renderings with transparent background if more than one timestep per page selected
- auto-magically map exact solution columns onto splash columns
- added --exact=file1,file2 to switch on plotting of exact solution from file(s)
- added --track=maxdens and --origin=maxdens to track/recentre on maximum density
- pressing backspace over legends deletes them
- use density weighted and normalised rendering by default in projection plots of vector fields

**3.7.2: (21/02/23)**

- bug fix recognising labels like v_{phi} on command line, can now use -r vphi

**3.7.1: (09/02/23)**

- libexact build failure fixed

**3.7.0: (09/02/23)**

- splash calc extinction computes column density to all sink particles in the simulation
- bug fix with rendering vector components (e.g. vr) in non-cartesian coordinate systems
- bug fix with both quantities appearing in black and white when double rendering

**3.6.0: (31/10/22)**

- skip particles with zero weight in interpolation, large speedup in some cases (thanks to T. Bending)
- splash calc plus and splash calc minus for adding/subtracting snapshots
- added --origin=6245 flag to centre the origin on particle 6245
- added --hdu=1 flag to read from a particular hdu in a fits file
- use wcs coordinates / arcseconds for fits images if present in header
- option --dense to reset to densest clump in phantom/sphNG data read (thanks to J. Wurster)

**3.5.1: (20/06/22)**

- bug fix with autolog limits
- build failures in libexact and libread fixed and now tested
- recognise labels on command line e.g. -r density
- limits option for centred cube (thanks to J. Wurster)

**3.5.0: (17/06/22)**

- bug fix with blank lines in splash.titles
- bug fix with large line lengths in csv files
- allow blank labels in csv headers
- bug fix with display of column labels from ascii/csv files
- log colour bar by default when using -r flag if more than 3 orders of magnitude range

**3.4.0: (24/03/22)**

- density weighted interpolation now applied automatically to projection plots of quantities that are not densities
- added flags --codeunits or --code to enforce code units from command line
- successfully parse csv files where some of the fields are character strings

**3.3.5: (01/03/22)**

- bug fix with disappearing sinks in phantom MPI dumps

**3.3.4: (21/01/22)**

- improved visual appearance of normalised renderings with free boundaries
- automatically read planet-wake parameters from phantom file headers
- added --wake=1,3 flag to plot wake from sink particle 3 around star 1
- bug fix with disappearing sinks in phantom MPI dumps
- fixed seg fault in fits reader

**3.3.3: (19/11/21)**

- "splash to csv" exports to comma separated variable (.csv) format
- automatically apply -ev flag for filenames ending in .ev, .mdot or .out
- improved label recognition from ascii file headers
- additional divergent colour schemes (thanks to Sahl Rowther)
- deal with merged sink particles from phantom (thanks to James Wurster)
- bug fix with units resetting to 1
- skip blank and comment lines in splash.filenames

**3.3.2: (20/07/21)**

- bug fix with -dev flag
- silenced unnecessary dust warnings in sphNG read
- change-of-limits animation sequence works for vector plots
- automatic recognition of ndspmhd format

**3.3.1: (19/07/21)**

- f/F in interactive mode flips y axis on 2D plots to next column
- gradual transparency in double rendering rather than sharp cutoff
- removed S from main menu as now redundant
- allow longer paths with -dev flag
- added --xmin,--xmax,--ymin,--ymax flags for manual margin adjustment
- bug fix with relativistic corrections in splash calc lightcurve

**3.3.0: (20/05/21)**

- bug fix with surface density plot with physical units on
- splash calc lightcurve computes spectra from local blackbody emission if T and kappa given
- lightcurve now performs frequency-dependent ray tracing
- added "--anglex","--angley","--anglez" flags
- can add labelled arrows by typing ^ in interactive mode, also delete/edit
- capital M, 0 or ncols+1 from main menu gives multiplot
- added -multi flag for multiplot from command line

**3.2.1: (26/04/21)**

- added --xsec=1.0 and --kappa=1.0 flags to specify cross section position and opacity, respectively
- specifying --xsec automatically switches from projection to cross section
- specifying --kappa turns on opacity rendering
- bug fix in splash calc tracks
- can use --track=1,2,3 to specify list of particles

**3.2.0: (20/04/21)**

- disable ALL prompts if any command line flags set
- all environment variables can now be given as command line flags using lower case string after last underscore e.g. SPLASH_CENTRE_ON_SINK=1 becomes --sink=1 on command line
- useful options include --corotate, --sink=1, --debug and more
- splash to grid recognises flags including --periodic, --npix=100,100,100 and --convert=1,4
- added -gandalf and -f gandalf as shortcut for seren data read
- assume default xw device and disable device prompt if any command line flags set
- s/S options now do the same thing

**3.1.1: (31/03/21)**

- automatically plot y vs x given a two-column data file
- planet wake coordinate system added
- bug fix with SPLASH_COROTATE
- bug fix reading phantom dumps when number of particles of each type does not match itype array
- bug fixes in grid2pdf

**3.1.0: (16/02/21)**

- splash calc lightcurve implemented
- sink particles ON by default
- changing units rescales plot limits correctly
- further improvements to ray tracing / opacity rendering with physical opacity
- can change units temporarily without writing .units file
- auto-select closest velocity and mass unit and better default time unit in phantom/sphNG read
- error message if Inf or NaN read from .units file
- bug fix with units prompt
- floating colour bars are white not black
- automatically write copyright in Hollywood mode
- auto-render fits files
- read softening length from phantom sinks if accretion radius is zero

**3.0.2: (20/01/21)**

- opacity rendering uses physical value of kappa, can also use opacity defined on particles
- can track multiple particles with 'splash calc tracks' by specifying ids in splash.tracks file
- support for SWIFT code in gadget_hdf5
- auto-recognise format for .csv files
- improved starsmasher data read
- improved physical unit selection
- exact solution lines can be plotted in background colour
- bug fix for dead particles in phantom dumps
- seg fault in fits reader fixed
- seg fault in gadget data read fixed
- bug fix in x-menu options

**3.0.0: (26/08/20)**

- Unified splash binary with -f flag to specify format
- automated format recognition for phantom, gadget (and hdf5 variants) and fits
- cleaner d) menu
- splash is compiled in double precision by default
- rotation settings used in splash to grid to rotate particles
- bug fix in mbatesph data read
- pysplash utility for reading SPH data formats into python
- libsplash.so, libexact.so and libread.so libraries

**2.10.1: (24/06/20)**

- exact solution can appear in legend
- can also plot under data
- fits reader and denoise utility can read/write spectral cubes
- text shapes can print header variables using %(var)
- can shift cross section by precise amounts in interactive mode using number followed by u/d
- fits reader includes header quantities
- reduced verbosity for non-interactive plots
- use of fake dust particles is now via menu option, not environment variable
- max particle types = 24
- userguide in readthedocs format
- bug fix with save limits with particle tracking
- support for .pfm pixelmap format as output
- physical units are ON by default
- prompts only for particle types present in data

**2.10.0: (14/02/20)**

- much improved splash to grid - bug fixes with pixel number and roundoff error
- use Petkova (2018) method for sub-pixel rendering to 3D grid and 3D projections
- added bytestream output formats for splash to grid and splash to ascii
- can press number and -/+ to zoom out/in by that factor in interactive mode
- use SPLASH_COROTATE=1,3 to corotate with arbitrary pair of sink particles
- SPLASH_COROTATE also gives velocity field in corotating frame
- splash to ascii can write particular columns by setting SPLASH_CONVERT=1,4
- plasma beta correct in both code and physical units
- working fits reader and splash-denoise utility

**2.9.1: (08/11/19)**

- cleaner menu options for units and calculated quantities
- surface rendering allowed with 3D perspective turned off
- automatic labelling of grain sizes in density and column density plots
- adaptive limits on log colour bars show 3 dex range by default
- auto-adjust limits to device aspect ratio works with multiple panels
- bug fixes with r-z rendering
- Toomre Q prompts for mass

**2.9.0: (05/04/19)**

- general header quantities are read and available in function parser
- more robust label detection and parsing during ascii data read
- splash to grid works in non-cartesian geometries
- added flared and log-flared coordinate systems
- Doppler shift colour bar
- can customise line style and colour when plotting multiple exact solutions
- seg faults fixed
- better plot tiling decisions
- disappearing arrows bug fix
- Rafikov disc- planet exact solution added
- atan2 implemented in function parser
- various multigrain phantom read fixes (incl. seg faults)
- exact rendering implemented in 2D
- libsplash implemented for use as Python splash backend

**2.8.0: (06/04/18)**

- 360/4pi video mode added
- automatically read labels from ascii file headers
- nearest sensible unit (e.g. au or pc) used by default
- cactus hdf5 data read
- kernel-smoothed particle plots of arbitrary quantities
- Viridis, Ocean and Inferno colour schemes
- can customise line colours
- Bondi flow exact solution
- option for ticks but no labels
- correct units in surface density plots
- colour bar on top or left
- support for multi-grain dust in Phantom
- bug fix with NaNs in ascii files

**2.7.0: (03/05/17)**

- Hollywood mode added (ctrl-m in interactive mode)
- better handling of dust/gas phantom data
- added rotated cartesian geometry
- rendering implemented in r-phi coordinates
- added Fortran 2008 intrinsics to function parser
- better rectangle plotting
- better falcON data read
- Ogilvie-Lubow exact solution for planet-disc interaction
- tipsy read now works when splash compiled in double precision
- splash to gridascii2 implemented
- bugs with r-phi rendering fixed

**2.6.0: (22/10/15)**

- SILO, falcON and .pbob data reads implemented
- bug fixes in gadget-hdf5 reader
- can recognise particle types in ascii read
- more robust sphNG read
- dust fraction recognised in phantom data read
- Toomre Q works in physical units
- bug fix with disappearing units labels
- bug fix in shock tube exact solution
- added splash calc delta
- splash to ascii keeps precision
- better power spectra

**2.5.1: (29/01/15)**

- error bar style options
- support for 5K displays
- can plot vectors and render with colours if h not read
- range restrictions apply during splash to grid
- improved line-style legend
- now up to 6 line styles
- fixes to amuse-hdf5 read
- phantom read handles star/dm particles
- various bugs fixed

**2.5.0: (22/08/14)**

- instant multiplots by giving multiple columns as y axis
- ability to plot multiple exact solution files on same plot
- compiles in parallel by default
- support for tagged sphNG/Phantom format
- AMUSE hdf5 format reader added
- various bug fixes

**2.4.1: (01/04/14)**

- Roche-lobe plotting vastly improved
- newunit= issue fixed
- bug fix with reading sink velocities from Phantom
- other minor bug fixes.

**2.4.0: (21/02/14)**

- time formatting in legend can include general functions like %(t + 1000)
- option to include sinks in opacity rendering
- supports one-fluid dust visualisation
- C-shock exact solution
- better polytrope solution

**2.3.1: (11/11/13)**

- SPLASH_COROTATE option to plot in frame corotating with sinks
- bug fixes with handling of dead/accreted/boundary particles in sphNG/phantom
- various other bugs fixed.

**2.3.0: (09/08/13)**

- can customise time formatting in legend
- improvements to legends
- less verboseness
- splash can read and plot pixel maps produced with -o ascii
- 3D vector field plotting improved
- bug fix with gfortran 4.8

**2.2.2: (10/05/13)**

- particle tracking by type implemented
- can interpolate specific columns in splash to grid
- SPLASH_CENTRE_ON_SINK option generic to all data reads
- Aly Reheam format added
- option for 2nd y axis on plots
- bug fix with X11 linking on Ubuntu
- can read gadget ICs files

**2.2.1: (21/02/13)**

- minor bug with axes plotting fixed
- Wendland kernels added
- bugs with exact solution plotting fixed
- bug fix with tracking of dark matter particles

**2.2.0: (16/11/12)**

- option to use different kernels for interpolation
- floating/inset colour bars added
- splash to gadget conversion implemented
- splash to grid works in 2D
- improved interfaces to shapes and animation sequences
- automatically turns on dark matter particle plotting if no gas
- interactive mode help displayed automatically

**2.1.1: (31/08/12)**

- irregular/circular particle selection using shift-left/middle click
- improved h5part and GADGET HDF5 data reads
- splash can be compiled in double precision
- bug fixes with calculated quantities + change of coordinate systems
- improved vector plot legend
- option for box+numbers but no labels added

**2.1.0: (16/05/12)**

- 3D vector field visualisation added
- GADGET HDF5 read implemented
- page sizes can be specified in pixels
- limits can auto-adapt to device aspect ratio
- more general exact solution from file option
- tiling works with one colour bar per row
- splash calc handles different particle types

**2.0.0: (29/08/11)**

- new giza backend - antialiased lines
- real fonts
- pdf, eps and svg drivers
- fewer build dependencies (only cairo, X11)
- support for semi-transparent text
- Double rendering (with transparent background) implemented.

**1.15.0: (29/08/11)**

- Multiplot with different particle types implemented
- calculated quantities list is now pre-filled automatically
- preliminary support for r-phi and r-z rendering
- outlined solid markers implemented
- better handling of multiple types
- manual contour levels can be specified in splash.contours
- parallel splash to grid
- better support for non-square pixels
- clipping of numbers at edge of viewport fixed

**1.14.1: (17/03/11)**

- SEREN data read added
- dragon read updated
- build follows Gnu conventions on DEST and DESTDIR (needed for macports build)
- can have up to 12 particle types
- exact solutions re-ordered
- dusty wave exact solution added

**1.14.0: (06/12/10)**

- Can flip between rendered quantities in interactive mode using 'f/F'
- SPLASH_DEFAULTS variable can be set for system-wide defaults
- can plot arbitrary functions of x,t as exact solution
- asplash better handles blank lines in header and can specify time, gamma location with env. variables
- added data read for the H5PART format
- GADGET read across multiple files implemented
- VINE read works with particle injection
- error bars can be plotted for both x and y axis simultaneously
- default rotation angles are set if 3D perspective turned on
- new directory layout and more helpful error messages during build
- PGPLOT linking is easier to get right.

**1.13.1: (26/02/10)**

- bugs with new calc_quantities module fixed
- generic library interface implemented so backend can be changed easily
- bug fix with auto pixel selection
- simpler foreground/background colour setting
- added subgrid interpolation warning

**1.13.0: (25/02/10)**

- function parser incorporated
- calculated quantities can now be specified at runtime, arbitrary function plotting implemented as an exact solution
- command-line SPH->grid conversion ("splash to grid") implemented
- ctrl-t in interactive mode adds arbitrary text box
- better line style/colour changing
- bug fix with tiling and y-axis labels
- various other bug fixes.

**1.12.2: (15/07/09)**

- Variable marker sizes added, can plot particles as circles with size proportional to h
- dark matter rendering with block-labelled GADGET format fixed
- VINE read handles star particles
- TIPSY read with ifort10.0.0 works
- snsph read added
- splash to phantom added
- does not override labels for coords, vectors by default
- bug fixes with contouring options
- stability bug fixes with older compilers
- more robust memory handling
- bug fix with automatic pixel selection causing seg fault.

**1.12.1: (20/04/09)**

- Can edit/delete text shapes interactively, also the colour bar label
- can customise the label on projection plots
- contour levels better defined
- SPLASH_HMIN_CODEUNITS added
- option for numeric labelling of contours
- contour limits can be set separately to render limits for same quantity
- minor bug fixes.

**1.12.0: (22/12/08)**

- Command-line plotting implemented
- ln transform added
- bug fixes in GADGET read
- Backspace over annotation (legends,titles,axes,colour bar) in interactive mode removes it
- "splash calc" command line utility calculates time sequences of global quantities from a sequence of dump files
- bug fix causing seg fault.

**1.11.1: (13/10/08)**

- automatic number of pixels and exact pixel boundaries implemented
- mass does not have to be read from dump file
- frame changes are per-page not per-dump file for animation sequences
- lower stacksize footprint
- bug fix with circles of interaction
- bug fixes with block-labelled GADGET read
- Steve Foulkes data read added.

**1.11.0: (15/08/08)**

- ability to use subset of particles in restricted parameter range(s)
- probability density function plot option
- plot-hugging colour bars added
- ability to annotate plot with a range of shapes
- v,V,w and H implemented in interactive mode for >1 panel
- various bug fixes (including one with vphi).

**1.10.2: (08/05/08)**

- disc surface density / toomre q parameter plotting added
- flash colour schemes added
- splash to binary convert option
- can change order in which particle types are plotted
- splash.columns file overrides column label settings
- vanaverbeke format read
- various bug fixes.

**1.10.1: (11/03/08)**

- "splash to" command line option converts binary dumps to ascii format
- vector plots + rotation now implemented
- block labelled GADGET format read
- ring-spreading exact solution added.

**1.10.0: (28/11/07)**

- horizontal colour bars implemented
- -p, -o command line options
- can have mixed types in data reads
- TIPSY and DRAGON data reads
- density weighted rendering
- normalisation applies to column density plots
- improved particle tracking
- save as option
- various bug fixes

**1.9.2: (12/09/07)**

- improvements to ascii read including asplash -e option
- smarter foreground/background colour changing for titles
- min=max problem fixed (caught by splash not pgplot)
- fixed vector arrow length option
- other minor changes and bug fixes

**1.9.1: (11/07/07)**

- environment variables + improvements to gadget data read
- better prompting
- 3 new colour schemes
- improved legend/title options
- other minor changes

**1.9.0: (21/05/07)**

- animation sequences implemented
- origin settings now affect radius calculation and are relative to tracked particle
- automatic line width choice for postscript devices
- w key adapts vector arrows
- vastly improved userguide

**1.8.1: (28/03/07)**

- option to hide vector arrows where there are no particles added
- smoother 3D plotting at low pixel numbers
- smoother vector plots
- bug fixes with a)
- issues with round-off error with z integration of vectors fixed.

**1.8.0: (14/03/07)**

- hidden particles not used in rendering
- units for z integration added
- a) & g) implemented in interactive mode for multiple-plots-per-page
- improved cross section using x in interactive mode

**1.7.2: (19/02/07)**

- Menu shortcuts implemented
- bug fix/ more sensible transformation of angular vector components in different co-ordinate systems
- improvements to interactive zoom and origin recentreing
- improved colour-by-type option
- restrictions on page size removed
- minor bug fixes

**1.7.1: (04/01/07)**

- command line options for defaults and limits files added
- minor bug fixes

**1.7.0: (13/12/06)**

- renamed SPLASH instead of SUPERSPHPLOT
- much faster data read for gadget and sphNG reads (only required columns read)
- physical units can be saved to file
- new menu formats
- various other bug fixes

**1.6.2: (24/10/06)**

- fast particle plotting and streamline plotting implemented
- more bug fixes with interactive mode on multiplots
- various other bug fixes

**1.6.1: (24/8/06)**

- bug fixes to 1.6.0, further improvements to interactive mode on multiplots

**1.6.0: (10/8/06)**

- Interactive mode on multiple plots per page
- highly optimised interpolation + parallel version
- new Makefile
- various bug fixes

**1.5.4: (06/7/06)**

- Handles multiple SPH/non-SPH particle types
- axes redrawn after rendering
- minor bug fixes

**1.5.3: (27/6/06)**

- minor bug fixes/improvements to multiple plots per page
- colour bar labelling tiled plots
- legend
- Accelerated rendering option for projections.

**1.5.2: (11/5/06)**

- "S)" option for saving limits and defaults
- MUCH faster interactive replotting (no unnecessary re-rendering)
- a few other minor things

**1.5.1: (26/4/06)**

- docs updated for v1.5, other minor changes

**1.5.0: (17/3/06)**

- 3D perspective added
- 3D opacity rendering
- improved rotation, colour schemes
- adjustable vector arrows (+legend)
- improved timestepping behaviour
- speed enhancements
- physical unit rescaling

**1.0.5: (28/9/05)**

- error calculation for exact solutions
- legend for plot markers
- exact_densityprofiles added
- more colour schemes
- unit rescaling improved
- other minor changes + bug fixes

**1.0.4: (17/8/05)**

- better colour schemes
- interactive colour scheme changing
- various minor changes and bug fixes

**1.0.3: (5/7/05)**

- rescale data option
- better page setup
- improved zooming
- interactive particle tracking
- various minor changes and bug fixes

**1.0.2 :**

- much improved ascii data read
- better line plotting
- zoom on powerspectrum plots + various bug fixes

**1.0.1 :**

- bug fixes relating to colour bars on multiplots

**1.0.0 :**

- first official release
- version given to many people at IPAM meeting and put on web 
