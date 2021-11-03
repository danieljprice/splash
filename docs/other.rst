
Other useful functionality
==========================

.. _sec:convert:

Conversion of binary data files to ascii
-----------------------------------------

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

Conversion of binary data files to csv
-----------------------------------------

To export data in csv format, just use:

::

   splash to csv dump001 dump002 dump???

which will convert all of the files listed on the command line into
csv format (called ``dump001.csv``, ``dump002.csv`` etc.), with
columns as would be listed in the main menu if you opened the dump file
in splash.

See also :ref:`sec:commandline`.

.. _sec:converttogrid:

Interpolation of SPH data to 2D and 3D grids
---------------------------------------------

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
ridiculous memory/disk usage. Options can be set to change
various aspects of the grid interpolation behaviour. For example you to
use periodic boundary conditions use::

  splash to grid dump001 --periodic

or alternatively you can set this as an environment variable::

  export SPLASH_TO_GRID_PERIODIC=yes
  splash to grid dump001

To change the number of pixels, you can use::

  splash to grid dump001 --periodic --npix=100,100,000

and to interpolate only columns 4 and 5 you can use::

  splash to grid dump001 --periodic --npix=100,100,000 --grid=1,5

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

For the ``energies`` and ``massaboverho`` options to be successful, splash
must be aware of the locations of the corresponding columns in the data
(i.e., by the column identification given in the set_labels routine
corresponding to the data read). For the ``massaboverho`` option an input
file is required specifying the density thresholds (a default version is
written if the appropriate file is not already present).

Using splash to time average a series of files
----------------------------------------------

The ``splash calc timeaverage`` command line option (see
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
function (PDF) from PDFs produced by splash (see :ref:`sec:pdfs:`).

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
line, no prompts will appear for the remaining options and default values will be assumed
for the remaining options. For example, the default device will be /xw giving an interactive plot.

.. _sec:pdfs:

Computing volume-weighted probability density functions from SPH data using SPLASH
-----------------------------------------------------------------------------------
The best way to compute a volume-weighted probability density function
on SPH particles is to interpolate the density field to a grid and compute
the histogram of the number of grid cells containing a given value of the desired quantity.

The grid2pdf utility included with splash can be used to compute the density PDF
from gridded data output by the ``splash to grid`` utility (see :ref:`sec:converttogrid`).

To use this feature, you will need to output grids in "binary" format, e.g::

   splash to gridbinary turb_00020

or if you want to skip the velocity interpolation (assuming density in column 6)::

   splash to gridbinary turb_00020 --grid=6

this produces a file called turb_00020.grid, then follow this with::

   cd $SPLASH_DIR; make grid2pdf
   cd -
   $SPLASH_DIR/splash/bin/grid2pdf turb_00020.grid

which produces::

   turb_00020.grid_pdf_ln_density.dat

this is just a two-column ascii file, so you can then plot this with your favourite plotting tool, e.g.::

   splash -ev turb_00020.grid_pdf_ln_density.dat


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

.. math::

   \int v_x \mathrm{dx}

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

.. _sec:writepixmap:

Outputting the raw pixel map to a file
--------------------------------------

The actual pixel map rendered to the graphics device (i.e., when a
quantity is rendered to pixels, not for particle plots) can be output
directly to a file, or series of files by using the ``-o`` command line
option when you invoke splash. This is useful if you need to compare the
image to the output from another code (e.g. using a different
visualisation tool) or if you wish to have a “raw” rendering, that is
without annotation on the plots.

Invoking splash with ``-o`` lists the currently implemented formats::

  possible formats for -o option:
  -o ppm   : dump pixel map to portable pixel map file
  -o pfm   : dump pixel map to portable float map file
  -o ascii : dump pixel map to ascii file

For example, to output the pixel map in ascii format, use::

   splash discG_00300 -o ascii -r 6 -dev /png

giving::

   > writing pixel map to file discG_00300_columndensitygcm2_proj.pix ...OK

This produces a file as follows::

  $ more discG_00300_columndensitygcm2_proj.pix
  # discG_00300_columndensitygcm2_proj.pix created by SPLASH
  # Contains 2D pixel array 512 x 512 written as
  #   do j=1,512
  #      write(*,*) dat(1:512,j)
  #   enddo
  # column density [g/cm^2]: min =   9.697428E-12 max =   7.487661E+03
  # x axis: min =  -4.000000E+03 max =   4.000000E+03
  # y axis: min =  -4.000000E+03 max =   4.000000E+03
  # 512 512
    0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+00  0.000000E+0
  ...

The number of pixels in the image can be controlled using the 'set number of pixels' option in the :ref:`sec:menu-r` (making sure you save the settings to the splash.defaults file using the :ref:`sec:menu-s`).

.. _sec:readpixmap:

Reading raw pixel maps from splash into Python
----------------------------------------------

See above for how to output the raw pixel map to a file. The resulting .pix file can be read into Python using the command::

  array = np.loadtxt('discG_00300_columndensitygcm2_proj.pix',skiprows=9)
  print (array.shape)
  plt.imshow(img)
  
A slightly more advanced script that also reads the x and y limits from the .pix file is provided in `splash/scripts/plot_pix.py <https://github.com/danieljprice/splash/blob/master/scripts/plot_pix.py>`_::

  python plot_pix.py discG_00300_columndensitygcm2_proj.pix

For other command line options, see :ref:`sec:commandline`.
