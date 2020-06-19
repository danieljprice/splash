
Other useful information
========================

.. _sec:convert:

Converting binary dump files to ascii using splash
---------------------------------------------------

splash has a command line feature which can be used to convert binary
SPH dump files into ascii format. The syntax is

::

   splash to ascii dump001 dump002 dump???

which will convert all of the dump files listed on the command line into
ascii format (called ``dump001.ascii``, ``dump002.ascii`` etc.), with
columns as would be listed in the main menu if you opened the dump file
in splash . Note that the output *includes* calculated extra quantities
such as the radius if these have been turned on [in the d) menu] and the
settings saved to the ``splash.defaults`` file. Similarly the data will
be output in physical units if a ``splash.units`` file is present.

For other command line options, see :ref:`sec:commandline`.

.. _sec:converttogrid:

Converting SPH data files to 3D gridded data using splash
----------------------------------------------------------

splash has a command line feature which can be used to read binary SPH
dump files and output 3D gridded data in a variety of formats. The
syntax is

::

   splash to grid dump001 dump002 dump???

which will interpolate the density, velocity (if present) and magnetic
field (if present) onto a 3D grid and output the results to files (the
default output format is ascii, with one file for each quantity
interpolated). Other data columns in the SPH file can be interpolated
using the “allto” option, which interpolates *all* of the columns to the
grid:

::

   splash allto grid dump001 dump002 dump???

The grid interpolation uses the :math:`x`, :math:`y`, and :math:`z`
limits — as saved to the ``splash.limits`` file — for the box, and the
grid size is given by the “set number of pixels” option in the r)ender
menu — as saved to the ``splash.defaults`` file. Automatic pixel
determination also works (if npixels = 0) but there is a sensible upper
limit placed on the grid size determined in this manner to avoid
ridiculous memory/disk usage. Various environment variable options are
available (these are output at runtime) that can be used to change
various aspects of the grid interpolation behaviour (e.g. setting
``SPLASH_TO_GRID_PERIODIC=yes`` enforces periodic boundary conditions).

For all possible output formats, use ``splash --help`` or see the full
list of command line options in :ref:`sec:commandline`.

.. _sec:splashcalc:

Using splash to calculate global quantities as a function of time.
------------------------------------------------------------------

splash has a command line feature that can be used to calculate global
quantities on the particles as a function of time, for example kinetic,
thermal, magnetic and total energy, total linear and angular momentum.
An example to calculate the energies in a sequence of dump files is:

::

   splash calc energies dump001 dump002 dump???

Other options are given by typing ’splash calc’, which currently has the
following options:

::

     splash calc energies     : calculate KE,PE,total energy vs time
                                output to file called 'energy.out'
            calc massaboverho : mass above a series of density thresholds vs time
                                output to file called 'massaboverho.out'
            calc max          : maximum of each column vs. time
                                output to file called 'maxvals.out'
            calc min          : minimum of each column vs. time
                                output to file called 'minvals.out'
            calc diff         : (max - min) of each column vs. time
                                output to file called 'diffvals.out'
            calc amp          : 0.5*(max - min) of each column vs. time
                                output to file called 'ampvals.out'
            calc delta        : 0.5*(max - min)/mean of each column vs. time
                                output to file called 'deltavals.out'
            calc mean         : mean of each column vs. time
                                output to file called 'meanvals.out'
            calc rms          : (mass weighted) root mean square of each column vs. time
                                output to file called 'rmsvals.out'
            calc timeaverage  : time average of *all* entries for every particle
                               output to file called 'time_average.out'
            calc ratio        : ratio of *all* entries in each file compared to first
                                output to file called 'ratio.out'

For the ‘energies’ and ‘massaboverho’ options to be successful, splash
must be aware of the locations of the corresponding columns in the data
(i.e., by the column identification given in the set_labels routine
corresponding to the data read). For the ‘massaboverho’ option an input
file is required specifying the density thresholds (a default version is
written if the appropriate file is not already present).

Using splash to time average a series of files
----------------------------------------------

The ‘splash calc timeaverage’ command line option (see
:ref:`sec:splashcalc`) can be used to produce a time average of a
series of files from any splash-readable format. This computes the
time-average of every individual entry in the file as represented in
splash as a table of rows (or ‘particles’) and columns (or ‘quantities
defined on particles’). The output is an ascii file with the same rows
and columns, averaged over all the snapshots on the command line. The
number of columns is doubled in the output, giving the standard
deviation for each quantity in the corresponding column (e.g., the
standard deviation for column 1 is output in column :math:`N + 1`).

Examples of how this could be use might be to produce the time-averaged
power spectrum from a series of ascii files containing power spectra for
individual output times, or the time averaged probability density
function (PDF) from PDFs produced by splash .

The resulting ascii file, called ``time_average.out`` can be plotted
using the ascii splash binary (asplash).

For other command line options, see :ref:`sec:commandline`.

.. _sec:batchmode:

Reading/processing data into images without having to answer prompts
--------------------------------------------------------------------

Previously, the only way to run splash non-interactively was to write a
small shell script which runs splash and answers the prompts
appropriately. For example:

::

   #!/usr/bin/tcsh
   cd plot
   splash myrun* << ENDINPUT
   2
   1
   8
   0
   /png
   q
   ENDINPUT

which would plot the data in columns 2 and 1 and render the data in
column 8 with output to file ``mypostscript.ps``.

However, in more recent versions splash can be invoked with plot options
on the command line. Thus to achieve the same as in the example given
above we would simply use

::

   splash myrun* -x 1 -y 2 -render 8 -dev /png

or simply

::

   splash myrun* -r 8 -dev /png

which will assume sensible default values (2 and 1 respectively) for the
y and x axes. Similarly a vector plot can be specified with ``-vec`` and
a contour plot with ``-cont``. The full list of command-line flags is
given in :ref:`sec:commandline`.

If plotting options have been only partially specified on the command
line, then prompts will appear for only the remaining options. This can
be used for example to specify the graphics device via the ``-dev``
command line option, which means that only the device selection prompt
does not appear.

Making frames across multiple processors
----------------------------------------

Making identical plots of a series of dump files for a movie is a task
which can inherently be done in parallel. Included in the splash/scripts
directory is a perl wrapper for splash (“``splash_parallel.pl``”) which
distributes multiple instances of splash across multiple machines,
either via ssh or using Apple’s xgrid, with a common input file as
described in :ref:`sec:batchmode`. The limitation to this is that
you need to have a disk which can be mounted from all client machines
(i.e., they can read the data files) and preferably with password-less
access (e.g. using an ssh key-exchange or Kerberos authentication). The
script itself may need some slight adjustment for your particular
system.

However, with large datasets often the slowest part of the rendering
process can be reading the data file. A good way of crippling a system
is therefore to set 100 jobs going which all decide to read a large data
file from disk at the same time. To avoid this the script allows the
user to set a delay between launching jobs (preferably slightly longer
than the length of time it takes to read a single dump file), but some
care is needed to avoid disaster. You have been warned!

What about boundaries? How does the rendering work near a boundary?
-------------------------------------------------------------------

Usual practise in SPH simulations near boundaries is to introduce ghost
particles which mirror the real particles. splash does not explicitly
setup any ghost particles but will use any that are present in the data
(see next question for how to specify multiple particle types).
Additional particle types contribute to the rendering calculations but
not to the determination of the plot limits. Note, however, that splash
does *not* set up ghost particles itself, as this may depend on the type
and location of the boundary. Thus if your simulation uses ghost
particle boundaries, the ghost particles should be dumped alongside the
gas particles in the output file so that their positions, masses,
densities and smoothing lengths can be read into splash and used to
render the image appropriately.

How does splash handle multiple particle types?
-----------------------------------------------

splash can handle up to 6 different particle types. These can be turned
on and off in the particle plot o)ptions menu (:ref:`sec:opts`).
These types are be specified in the set_labels part of the read_data
routine, which contains some lines of code along the lines of:

::

   ntypes = 3
   labeltype(1) = 'gas'
   labeltype(2) = 'ghost'
   labeltype(3) = 'sink'
   UseTypeInRenderings(1) = .true.
   UseTypeInRenderings(2) = .true.
   UseTypeInRenderings(3) = .false.

which says that there are 3 particle types, with names as given, and
that types 1 and 2 are SPH particles and should be used in the rendering
where appropriate (i.e., only when plotting of this type is turned on in
the o)pts menu). Particle types which are to be used in renderings
should have masses, densities and smoothing lengths read. Non-SPH
particle types (e.g. sink particles) can be optionally plotted on top of
rendered plots.

Using special characters in the plot labels
-------------------------------------------

Several of the examples shown in this manual use special characters
(such as the :math:`\int` character) in the plot labels. In giza these
can be specified using TeX-like escape sequences, or with the escape
sequences used in pgplot. For example to plot the greek letter
:math:`\rho` we would use

::

   label = 'this would print the greek letter \rho'

or, in pgplot-style:

::

   label = 'this would print the greek letter \gr'

where ``\gr`` is the pgplot escape sequence for :math:`\rho`.

   In giza , which uses real fonts rather than the bitmapped characters
   used in pgplot, special characters are implemented with unicode
   characters. Thus, you need to select a font that has the appropriate
   characters included. The font can be changed using the ``GIZA_FONT``
   environment variable.

For other characters the procedure is similar. For example for the
integral

.. math:: \int v_x \mathrm{dx}

 we would use the TeX-like expression

::

   label = '\int v_x dx'

or equivalently, in pgplot-style

::

   label = '\(2268) v\d x \u dx'

where ``\(2268)`` is the pgplot escape sequence for the integral sign.
The ``\d`` indicates that what follows should be printed as subscript
and ``\u`` correspondingly indicates a return to normal script (or from
normal script to superscript). All of the escape sequences for special
characters are listed in the appendix to the pgplot user guide.

   WARNING: Note that the use of escape characters can be compiler
   dependent and may not therefore work on all compilers (for example
   the intel compiler needs the -nbs flag).

Making movies
-------------

See :ref:`sec:movies` and the online FAQ
(http://users.monash.edu.au/~dprice/splash/faqs.html).

.. _sec:writepixmap:

Outputting the raw pixel map to a file
--------------------------------------

The actual pixel map rendered to the graphics device (i.e., when a
quantity is rendered to pixels, not for particle plots) can be output
directly to a file, or series of files by using the ``-o`` command line
option when you invoke splash . Invoking splash with ``-o`` produces a
list of currently implemented formats (at the moment these are an ascii
dump file and ppm format). This is useful if you need to compare the
image to the output from another code (e.g. using a different
visualisation tool) or if you wish to have a “raw” rendering, that is
without annotation on the plots, but which (in the ppm case) uses more
colours. The files are given default names such as “splash_00001.dat” or
“splash_00001.ppm” where the number corresponds to the frame number as
would be rendered to the graphics device.

For other command line options, see :ref:`sec:commandline`.

