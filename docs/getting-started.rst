
Getting started
===============

.. _install:

Installing from your package manager
-------------------------------------
Pre-packaged builds of splash exist for most operating systems.

Stable version
~~~~~~~~~~~~~~
Mac OS via homebrew (recommended)::

  brew tap danieljprice/all
  brew install splash

You will also need to install `Xquartz <https://www.xquartz.org>`_ so that the X-windows server launches automatically.

Mac OS via Macports::

  sudo port install splash

Linux or Windows Linux Subsystem (Ubuntu)::

  sudo apt-get install splash

Development version
~~~~~~~~~~~~~~~~~~~

SPLASH and giza (the plotting backend) both have public repositories, so you can check out the latest and greatest code at any time. Both the splash and giza repositories are generally stable, so it is usually safe to get the very latest version. Just use:

Mac OS via homebrew::

  brew tap danieljprice/all
  brew install --HEAD splash

or compile from source following the instructions below.

Compiling from source
---------------------

Basic instructions
~~~~~~~~~~~~~~~~~~
If you have admin (super user) permissions::

   git clone https://github.com/danieljprice/giza.git; cd giza; ./configure; make; sudo make install; cd ..
   git clone https://github.com/danieljprice/splash.git
   cd splash; make SYSTEM=gfortran; sudo make install

.. _installhome:

Installing in your home space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you do not have admin permissions. That is, to install in your home space::

   git clone https://github.com/danieljprice/splash.git
   cd splash; git clone https://github.com/danieljprice/giza.git
   make SYSTEM=gfortran withgiza

.. important::
   If you have installed splash in your home space, you will need to set the following environment variables for everything to work. Put the following commands in your ~/.bashrc file or equivalent, so they are set every time you log in::

      export SPLASH_DIR=$HOME/splash
      export PATH=$PATH:$SPLASH_DIR/bin
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SPLASH_DIR/giza/lib

Advanced installation guide
---------------------------

The basic steps for installation are as follows:

#. make sure you have a recent Fortran compiler (such as gfortran)

#. make sure you have cairo on your system

#. compile splash and giza

Fortran compilers
~~~~~~~~~~~~~~~~~~~

Numerous Fortran 90/95/2003 compilers exist. The most widely
available are:

-  gfortran, the free Gnu Fortran Compiler
   http://gcc.gnu.org/wiki/GFortran

-  ifort, one of the most widely available commercial compilers (and is
   very good) with (limited) free licence for Linux.
   http://software.intel.com/en-us/articles/intel-compilers/

Both of these successfully compile splash and the giza library.

Cairo graphics library
~~~~~~~~~~~~~~~~~~~~~~~
Cairo is a low-level system library used in many applications. Thus it is highly
likely that you already have a copy on your system and already in your library path.
Look for the header file cairo.h, e.g. using::

   locate cairo.h

or have a look in the usual places (e.g. /usr/include/cairo, /usr/X11/include). If not,
you can usually use your inbuilt package manager to install cairo as follows:

   Debian/Ubuntu:
      sudo apt-get install libcairo2-dev
   Fedora/Red Hat/CentOS:
      sudo yum install cairo-devel
   OpenSUSE:
      zypper install cairo-devel
   MacPorts:
      sudo port install cairo

Alternatively, use the script provided in the root-level splash directory::

   ./install-cairo.sh

which downloads and installs both pixman and cairo into the giza/ subdirectory.
Unlike the methods above, this does not require any admin/superuser permissions.

Compiling and linking with giza
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can either install giza with your package manager, or in a subdirectory
of splash. To install in a splash subdirectory, use::

	cd splash
	git clone http://github.com/danieljprice/giza
   make withgiza

of splash.

With giza installed via your package manager (or previously compiled as below), use::

   cd splash
   make GIZA_DIR=/usr/local

where ``GIZA_DIR`` points to the directory where giza was installed.
To install giza in a splash subdirectory, use::

   cd splash
   git clone http://github.com/danieljprice/giza
   make withgiza

A successful ``make`` will produce a binary called ``splash``

Reading your data
-----------------

The most important part is getting splash to read \*your\* data format.
If you are using a publicly available code, it is reasonably likely
that I have already written a read data subroutine which will read your
dumps. If not it is best to look at some of the other examples and
change the necessary parts to suit your data files. Note that reading
directly from unformatted data files is \*much\* faster than reading
from formatted (ascii) output.

A standard ``make`` will create a binary which supports the file formats listed in
:ref:`tab:defaultreads`, plus a bunch of others (type ``splash --formats`` to see what formats your build supports).
All data formats are supported in the splash binary by default unless there 
are external library dependencies (e.g. ``HDF5``) .

The format of the file can in many cases be successfully guessed from the file extension or header, so you can simply type::

	splash disc_00000

Otherwise you can specify the data type you are reading using the ``-f`` or ``--format`` flag. For example,
the following will read a phantom dumpfile::

	splash --format phantom disc_00000

For backwards compatibility with previous version of ``splash``, one can add aliases into their `.bashrc`, or equivalent::

   alias asplash='splash ' # Alias for ascii splash
   alias ssplash='splash -f phantom '
   alias gsplash='splash -f gadget '
   alias tsplash='splash -f tipsy '

If splash is compiled with ``HDF5=yes``, the formats listed in :ref:`tab:hdf5reads` will also be available in the ``splash`` binary.
Other supported formats are listed in :ref:`tab:otherreads`, but these require additional libraries.

.. table:: Binaries and data reads compiled by default
   :name: tab:defaultreads

   +------------------------------+----------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash`` command           | Format Read                | ``read_data`` File            |  Comments                                                                                                                                                                                                                                         |
   +==============================+============================+===============================+===================================================================================================================================================================================================================================================+
   | ``splash -ascii <file>``     | ascii, csv                 | ``read_data_ascii.f90``       | Generic data read for n-column ascii formats. Automatically determines number of columns and skips header lines. Can recognise SPH particle data based on the column labels. Use ``splash -e`` to plot non-SPH data (e.g.  energy vs time files)  |
   +------------------------------+----------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash -dragon <file>``    | dragon                     | ``read_data_dragon``          | See environment variable  options.                                                                                                                                                                                                                |
   +------------------------------+----------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash -gadget <file>``    | gadget, gadget-2, gadget-3 | ``read_data_gadget.f90``      | Handles both default and block-labelled formats (see  environment variable  options).                                                                                                                                                             |
   +------------------------------+----------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash -ndspmhd <file>``   | ndspmhd                    | ``read_data_ndspmhd.f90``     | Format for the ndspmhd SPH/SPMHD code (publicly available from  my  website).                                                                                                                                                                     |
   +------------------------------+----------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash -phantom <file>``   | sphNG, Phantom             | ``read_data_sphNG.f90``       | sphNG is Matthew Bate’s SPH code. Option ``-sphng`` also  works.                                                                                                                                                                                  |
   +------------------------------+----------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash -magma <file>``     | magma                      | ``read_data_srosph.f90``      | Stephan Rosswog’s  code                                                                                                                                                                                                                           |
   +------------------------------+----------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash -seren <file>``     | seren                      | ``read_data_seren.f90``       | The SEREN SPH code (Hubber, McLeod et  al.)                                                                                                                                                                                                       |
   +------------------------------+----------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash -gasoline <file>``  | gasoline, tipsy            | ``read_data_tipsy.f90``       | Reads both binary and ascii TIPSY files (determined automatically). Option ``-tipsy`` also  works.                                                                                                                                                |
   +------------------------------+----------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash -vine <file>``      | vine                       | ``read_data_fine.f90``        | See environment variable  options.                                                                                                                                                                                                                |
   +------------------------------+----------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   |``splash -starsmasher <file>``| StarSmasher                | ``read_data_starsmasher.f90`` | The `StarSmasher <http://jalombar.github.io/starsmasher/>`_ code (Gaburov et al. 2018)                                                                                                                                                            |
   +------------------------------+----------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash file.vtk``          | vtk legacy binary          | ``read_data_vtk.f90``         | VTK legacy binary format, e.g. from Shamrock code                                                                                                                                                                                                 |
   +------------------------------+----------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Below is a list of the supported data formats that require ``HDF5``.

.. table:: Supported HDF5 data formats
   :name: tab:hdf5reads

   +--------------------------------+------------------------+-------------------------------+-----------------------------------------------------------------------------------------+
   | ``splash`` Command             | Read Format            | ``read_data`` File            |    Comments                                                                             |
   +================================+========================+===============================+=========================================================================================+
   | ``splash -gadget_hdf5 <file>`` | gadget HDF5 Files.     | ``read_data_gadget_hdf5.f90`` | Reads HDF5 format from the gadget code (automatically recognised)                       |
   +--------------------------------+------------------------+-------------------------------+-----------------------------------------------------------------------------------------+
   | ``splash -amuse <file>``       | AMUSE HDF5             | ``read_data_amuse_hdf5.f90``  | Reads HDF5 format from the AMUSE framework.                                             |
   +--------------------------------+------------------------+-------------------------------+-----------------------------------------------------------------------------------------+
   | ``splash -cactus_hdf5 <file>`` | Cactus HDF5            | ``read_data_cactus_hdf5.f90`` |                                                                                         |
   +--------------------------------+------------------------+-------------------------------+-----------------------------------------------------------------------------------------+
   | ``splash -flash_hdf5 <file>``  | FLASH tracer particles | ``read_dataflash_hdf5.f90``   | Reads tracer particle output from the FLASH code. The option ``-flash`` will also work. |
   +--------------------------------+------------------------+-------------------------------+-----------------------------------------------------------------------------------------+
   | ``splash -falcon_hdf5 <file>`` | falcON                 | ``read_data_falcON_hdf5.f90`` | Walter Dehnen’s SPH code format. The option ``-falcon`` will also work.                 |
   +--------------------------------+------------------------+-------------------------------+-----------------------------------------------------------------------------------------+

If the ``HDF5`` read files end with ``.h5``, the suffix ``_hdf5`` from the ``splash`` command can be removed.
For example::

	splash -gadget dump_000.h5

will recognise that the file ``dump_000.h5`` is in the ``HDF5`` format, and will automatically select the correct ``read_data`` routine.

Below is a list of other formats supported, but have additional library requirements.

.. table:: Other supported file formats that require external libraries
   :name: tab:otherreads

   +---------------------------+--------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash`` Command        | Read Format  | ``read_data`` File       | Comments                                                                                                                                                                   |
   +===========================+==============+==========================+============================================================================================================================================================================+
   | ``splash -pbob <file>``   | PBOB Files   | ``read_data_pbob.f90``   | Requires the PBOB Library. Compile ``splash`` with ``PBOB_DIR=/path/to/pbob/``.                                                                                            |
   +---------------------------+--------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash -h5part <file>`` | H5Part Files | ``read_data_h5part.f90`` | Requires the H5Part and HDF5 libraries. Compile ``splash`` with ``H5PART_DIR=/path/to/h5part/``.                                                                           |
   +---------------------------+--------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | ``splash file.fits``      | FITS files   | ``read_data_fits.f90``   | Requires FITS libraries. Try to compile ``splash`` with ``FITS=yes``. If this does not work, point to the location of your fits libraries with ``FITS_DIR=/path/to/fits``. |
   +---------------------------+--------------+--------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


Further details on writing your own subroutine are given in
appendix :ref:`sec:writeyourown`. The \*easiest\* way is to i)
email me a sample data file and ii) the subroutine you used to write it,
and I will happily create a data read for your file format.

.. _sec:commandline:

Command line options
--------------------

Typing ``splash --help`` gives a complete and up-to-date list of options. Currently these are:

::

   Command line options:
    
    -f format         : input file format to be read (default is ascii, --formats for full list)
    -p fileprefix     : change prefix to ALL settings files read/written by splash 
    -e, -ev           : use default options best suited for line plotting (.ev files)
    -360              : set default options suited to 360 video
    -b, --buffer      : buffer all data files into memory
    -o pixformat      : dump pixel map in specified format (use just -o for list of formats)
    
   Command line plotting mode:
    
    -x column         : x axis
    -y column         : y axis
    -r[ender] column  : column to render (will use 1 and 2 as x,y if -x,-y not specified)
    -vec[tor] column  : vector quantity to plot with arrows
    -c[ontour] column : contoured quantity
    -multi            : multiplot
    -dev device       : specify plotting device on command line (e.g. -dev /xw)
    --movie           : shortcut for -dev /mp4 to make a movie from plot sequence
    --xsec=1.0        : specify location of cross section slice
    --kappa=1.0       : specify opacity, and turn on opacity rendering
    --anglex=30       : rotate around x axis (similarly --angley, --anglez)
    --code            : enforce code units (also --codeunits)
    --sink=1          : centre on sink particle number 1
    --origin=666      : set coordinate system origin to particle number 666
    --origin=maxdens  : set coordinate system origin to particle at maximum density
    --track=666       : track particle number 666
    --track=maxdens   : track particle at maximum density
    --exact=file1,f2  : read and plot exact solution from ascii files file1 and f2
    --sort            : sort filenames for comparison (e.g. snap_000 snap1_000 snap2_000)
    
   Example data formats (type --formats for full list):
    
    -ascii,-csv          : ascii text/csv format (default)
    -phantom -sphng      : Phantom and sphNG codes (auto)
    -vtk                 : vtk legacy binary format (auto)
    -ndspmhd             : ndspmhd code (auto)
    -gandalf,-seren      : Gandalf/Seren code
    -gadget -gadget_hdf5 : Gadget code (auto)
    -falcon -falcon_hdf5 : FalcON code
    -flash  -flash_hdf5  : FLASH code
    -cactus -cactus_hdf5 : Cactus code
    -amuse  -amuse_hdf5  : AMUSE Framework
    -fits                : FITS format (auto)
    
    HDF5 files will be automatically recognised if they end with .h5, however you
    must specify a supported data format.
    add a suffix "_hdf5" to above format if your data files do not end with .h5.

    convert mode ("splash to X dumpfiles"):
    splash to ascii   : convert SPH data to ascii file dumpfile.ascii

           to binary  : convert SPH data to simple unformatted binary dumpfile.binary
                         write(1) time,npart,ncolumns
                         do i=1,npart
                            write(1) dat(1:ncolumns),itype
                         enddo
           to phantom : convert SPH data to binary dump file for PHANTOM
           to gadget  : convert SPH data to default GADGET snapshot file format

    Grid conversion mode ("splash to X dumpfiles"):
       splash to grid         : interpolate basic SPH data (density, plus velocity if present in data)
                                to 2D or 3D grid, write grid data to file (using default output=ascii)
              to gridascii    : as above, grid data written in ascii format
              to gridbinary   : as above, grid data in simple unformatted binary format:
                                   write(unit) nx,ny,nz,ncolumns,time                 [ 4 bytes each ]
                                   write(unit) (((rho(i,j,k),i=1,nx),j=1,ny),k=1,nz)  [ 4 bytes each ]
                                   write(unit) (((vx(i,j,k), i=1,nx),j=1,ny),k=1,nz)  [ 4 bytes each ]
                                   write(unit) (((vy(i,j,k), i=1,nx),j=1,ny),k=1,nz)  [ 4 bytes each ]
                                   write(unit) (((...(i,j,k),i=1,nx),j=1,ny),k=1,nz)  [ 4 bytes each ]
           allto grid         : as above, interpolating *all* columns to the grid (and output file)
           allto gridascii    : as above, with ascii output
           allto gridbinary   : as above, with binary output

    Analysis mode ("splash calc X dumpfiles") on a sequence of dump files:
     splash calc energies     : calculate KE,PE,total energy vs time
                                output to file called 'energy.out'
            calc massaboverho : mass above a series of density thresholds vs time
                                output to file called 'massaboverho.out'
            calc max          : maximum of each column vs. time
                                output to file called 'maxvals.out'
            calc min          : minimum of each column vs. time
                                output to file called 'minvals.out'
            calc diff           : (max - min) of each column vs. time
                                output to file called 'diffvals.out'
            calc amp          : 0.5*(max - min) of each column vs. time
                                output to file called 'ampvals.out'
            calc delta        : 0.5*(max - min)/mean of each column vs. time
                                output to file called 'deltavals.out'
            calc mean         : mean of each column vs. time
                                output to file called 'meanvals.out'
            calc rms          : (mass weighted) root mean square of each column vs. time
                                output to file called 'rmsvals.out'
            calc tracks       : track particle data vs time for selected particles,
               --track=1,2,3    output to tracks-1.out,tracks-2.out,tracks-3.out

     the above options all produce a small ascii file with one row per input file.
     the following option produces a file equivalent in size to one input file (in ascii format):

            calc timeaverage  : time average of *all* entries for every particle
                                output to file called 'time_average.out'

            calc ratio        : ratio of *all* entries in each file compared to first
                                output to file called 'ratio.out'

            calc plus         : add two snapshots together
                                output to file called 'plus.out'

Command-line options can be entered in any order on the command line
(even after the dump file names). For more information on the convert
utility (``splash to ascii``) see :ref:`sec:convert`. For details
of the ``-o ppm`` or ``-o ascii`` option see :ref:`sec:writepixmap`. For details of the ``-ev`` option, see :ref:`sec:evsplash`.


Options affecting all data reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Command line flags (or environment variables) that affect all data reads are:

+----------------------+-----------------------+-------------------------------------------------+
| ---defaults          | SPLASH_DEFAULTS       | gives name of system-wide ``splash.defaults``   |
|                      |                       | file (and splash.limits etc.) that will be      |
|                      |                       | used if there is none in the current dir. e.g.  |
|                      |                       | ``--defaults=$HOME/splash.defaults``            |
+----------------------+-----------------------+-------------------------------------------------+
| ---kernel='quintic'  | SPLASH_KERNEL         | changes the smoothing kernel used in the        |
|                      |                       | interpolations (e.g. ``cubic`` or ``quintic``). |
|                      |                       | Can also be changed in the :ref:`sec:menu-r`.   |
+----------------------+-----------------------+-------------------------------------------------+
| ---debug             | SPLASH_DEBUG          | if set to ``yes`` or ``true``, turns on verbose |
|                      |                       | debugging output. Useful to trace code crashes  |
|                      |                       | (but of course, this never happens…).           |
+----------------------+-----------------------+-------------------------------------------------+
| ---sink=1            | SPLASH_CENTRE_ON_SINK | if set to a number n, centres coordinates and   |
|                      |                       | velocities on the nth sink/star particle (e.g.  |
|                      |                       | ``export SPLASH_CENTRE_ON_SINK=2``)             |
+----------------------+-----------------------+-------------------------------------------------+
| ---corotate          | SPLASH_COROTATE       | plot in corotating frame based on locations of  |
|                      |                       | 2 sink particles (e.g. ``--corotate=1,3``)      |
+----------------------+-----------------------+-------------------------------------------------+
| ---origin=666        | SPLASH_ORIGIN         | recentre the coordinate origin and velocities   |
|                      |                       | to the selected particle (e.g. particle 666)    |
+----------------------+-----------------------+-------------------------------------------------+
| ---origin=maxdens    | SPLASH_ORIGIN         | reset origin to particle at maximum density     |
+----------------------+-----------------------+-------------------------------------------------+
| --dontcentrevel='y'  | SPLASH_DONTCENTREVEL  | used along with SPLASH_CENTRE_ON_SINK or        |
|                      |                       | SPLASH_ORIGIN. If true, then the velocities     |
|                      |                       | will not be made relative to the sink or        |
|                      |                       | particle.                                       |
+----------------------+-----------------------+-------------------------------------------------+
| ---track=4789        | SPLASH_TRACK          | set limits of all quantities relative to those  |
|                      |                       | of the selected particle (e.g. particle 4789)   |
+----------------------+-----------------------+-------------------------------------------------+
| ---track=maxdens     | SPLASH_TRACK          | track the particle at maximum density           |
+----------------------+-----------------------+-------------------------------------------------+
| ---vzero=1.0,1.0,1.0 | SPLASH_VZERO          | subtract reference velocity from all particles  |
|                      |                       | (velocity should be specified in code units)    |
+----------------------+-----------------------+-------------------------------------------------+
| ---exact=file1,file2 |                       | plot exact solution from files file1 and file2  |
+----------------------+-----------------------+-------------------------------------------------+
| ---beam=2.0          | SPLASH_BEAM           | if given a value :math:`>`\ 0 enforces a minimum|
|                      |                       | smoothing length, specified in code units,      |
|                      |                       | all the particles. Useful to “dumb-down” the    |
|                      |                       | resolution of SPH simulations to match          |
|                      |                       | observational resolution. If this variable is   |
|                      |                       | set the “accelerated rendering" option in the   |
|                      |                       | :ref:`sec:menu-r` is also turned on, otherwise  |
|                      |                       | slow rendering can result.                      |
+----------------------+-----------------------+-------------------------------------------------+
| ---xmin=0.1          | SPLASH_MARGIN_XMIN    | can be used to manually adjust the left page    |
| ---xmax=0.1          | SPLASH_MARGIN_XMAX    | page margin (set to fraction of viewport,       |
| ---ymin=0.1          | SPLASH_MARGIN_YMIN    | negative values are allowed).                   |
| ---ymax=0.1          | SPLASH_MARGIN_YMAX    |                                                 |
+----------------------+-----------------------+-------------------------------------------------+

.. _sec:splash:

Ascii data read
~~~~~~~~~~~~~~~~

For several data reads there are command-line flags which can be set
at runtime which are specific to the data read. For the ascii data read
(``splash -f ascii``) these are:

+-------------------------------------+-----------------------------------+
| ``--ncolumns=10``                   | if given a value :math:`>`\ 0     |
|                                     | sets the number of columns to be  |
|                                     | read from ascii data (overrides   |
|                                     | the automatic number of columns   |
|                                     | determination).                   |
+-------------------------------------+-----------------------------------+
| ``--nheaderlines=3``                | if given a value :math:`>=`\ 0    |
|                                     | sets the number of header lines   |
|                                     | to skip (overrides the automatic  |
|                                     | determination).                   |
+-------------------------------------+-----------------------------------+
| ``--columnsfile=/home/me/mylabels`` | can be used to provide the        |
|                                     | location of (path to) the default |
|                                     | ``columns`` file containing the   |
|                                     | labels for ascii data. Overridden |
|                                     | by the presence of a local        |
|                                     | ``columns`` file.                 |
+-------------------------------------+-----------------------------------+
| ``--time=3.0``                      | if given a nonzero value sets the |
|                                     | time to use in the legend (fixed  |
|                                     | for all files)                    |
+-------------------------------------+-----------------------------------+
| ``--gamma=1.667``                   | if given a nonzero value sets     |
|                                     | gamma to use in exact solution    |
|                                     | calculations (fixed for all       |
|                                     | files)                            |
+-------------------------------------+-----------------------------------+
| ``--timeheader=1``                  | sets the integer line number      |
|                                     | where the time appears in the     |
|                                     | header                            |
+-------------------------------------+-----------------------------------+
| ``--gammaheader=3``                 | sets the integer line number      |
|                                     | where gamma appears in the header |
+-------------------------------------+-----------------------------------+

The above options can be set as environment variables by prefixing them
with ASPLASH, e.g.::

  export ASPLASH_NCOLUMNS=10
  splash datafile.txt

.. _sec:splash -gadget:

GADGET data read
~~~~~~~~~~~~~~~~~

For the GADGET read (``splash -f gadget`` or just ``splash``) the options are:

+-----------------------------------+-----------------------------------+
| ``--format=2``                    | if set = 2, force read of block   |
|                                   | labelled GADGET format instead of |
|                                   | the default (non block labelled)  |
|                                   | format.                           |
+-----------------------------------+-----------------------------------+
| ``--usez``                        | if ``yes`` or ``true`` uses the   |
|                                   | redshift in the legend instead of |
|                                   | code time.                        |
+-----------------------------------+-----------------------------------+
| ``--hsoft=500.``                  | if given a value :math:`>` 0.0    |
|                                   | will assign a smoothing length to |
|                                   | dark matter particles for which   |
|                                   | rendered plots of column density  |
|                                   | can then be made.                 |
+-----------------------------------+-----------------------------------+
| ``--extracols``                   | if set to a comma separated list  |
|                                   | of column labels, will attempt to |
|                                   | read additional columns           |
|                                   | containing gas particle           |
|                                   | properties beyond the end of the  |
|                                   | file (not applicable if           |
|                                   | --format=2).                      |
+-----------------------------------+-----------------------------------+
| ``--starpartcols``                | if set to a comma separated list  |
|                                   | of column labels, will attempt to |
|                                   | read additional columns           |
|                                   | containing star particle          |
|                                   | properties beyond the end of the  |
|                                   | file (and after any extra gas     |
|                                   | particle columns) (not applicable |
|                                   | if GSPLASH_FORMAT=2).             |
+-----------------------------------+-----------------------------------+
| ``--checkids``                    | if set to ``yes`` or ``true``,    |
|                                   | reads and checks particle IDs,    |
|                                   | excluding particles with negative |
|                                   | IDs as accreted (gives them a     |
|                                   | negative smoothing length which   |
|                                   | means they are ignored in         |
|                                   | renderings).                      |
+-----------------------------------+-----------------------------------+
| ``--hcolumn``                     | if set to a positive integer,     |
|                                   | specifies the location of the     |
|                                   | smoothing length in the columns,  |
|                                   | overriding any default settings.  |
+-----------------------------------+-----------------------------------+
| ``--ignore-iflagcool``            | if set,does                       |
|                                   | not assume that extra columns are |
|                                   | present even if the cooling flag  |
|                                   | is set in the header.             |
+-----------------------------------+-----------------------------------+

For backwards compatibility, the above options can also be set as
environment variables by prefixing them with GSPLASH, e.g.::

  export GSPLASH_FORMAT=2
  splash -gadget snap_00010

which is equivalent to::

  splash -f gadget --format=2 snap_00010

For the GADGET read gsplash will also look for, and read if present,
files called ``snapshot_xxx.hsml`` and/or ``snapshot_xxx.dens`` (where
``snapshot_xxx`` is the name of the corresponding GADGET dump file)
which contain smoothing lengths and/or a density estimate for dark
matter particles (these should just be one-column ascii files).

VINE data read
~~~~~~~~~~~~~~~

For the VINE read (``splash -vine``) the options are:

+-----------------------------------+-----------------------------------+
| ``--hfac``                        | if ``yes`` or ``true`` multiplies |
|                                   | smoothing length read from the    |
|                                   | dump file by a factor of 2.8 (for |
|                                   | use with older VINE dumps where   |
|                                   | the smoothing length is defined   |
|                                   | as in a Plummer kernel rather     |
|                                   | than as the usual SPH smoothing   |
|                                   | length).                          |
+-----------------------------------+-----------------------------------+
| ``--mhd``                         | if set, reads VINE                |
|                                   | dumps containing MHD arrays       |
|                                   | (or set VINE_MHD=yes)             |
+-----------------------------------+-----------------------------------+

Phantom/sphNG data read
~~~~~~~~~~~~~~~~~~~~~~~~

For the sphNG and PHANTOM read (``splash -phantom``) the options are:

+-------------------+-------------------------------------------------------+
| ``--cm``          | resets the positions such that the centre of          |
|                   | mass is exactly at the origin.                        |
+-------------------+-------------------------------------------------------+
| ``--dense``       | resets the positions such that the centre of          |
|                   | mass of the densest clump is exactly at the origin.   |
+-------------------+-------------------------------------------------------+
| ``--omega=3.142`` | if non-zero, subtracts solid body rotation with omega |
|                   | as specified to give velocities in co-rotating frame  |
+-------------------+-------------------------------------------------------+
| ``--omegat=3.142``| same as --omega but applies to velocities also        |
+-------------------+-------------------------------------------------------+
| ``--timeunit=hrs``| sets default time units, either ’s’, ’min’, ’hrs’,    |
|                   | ’days’, ’yr’ or ’tfreefall’ (used verbatim in legend) |
+-------------------+-------------------------------------------------------+

dragon data read
~~~~~~~~~~~~~~~~~

For the dragon read (``splash -dragon``) the options are:

+-----------------------------------+-----------------------------------+
| ``--extracols``                   | specifies number of extra columns |
|                                   | present in the file which are     |
|                                   | dumped after the itype array      |
+-----------------------------------+-----------------------------------+

Stephan Rosswog data read
~~~~~~~~~~~~~~~~~~~~~~~~~~

For the srosph read (``splash``) the options are:

+-----------------------------------+-----------------------------------+
| ``--format=MHD``                  | can be ``MHD`` or ``HYDRO`` which |
|                                   | read the appropriate data format  |
|                                   | from either the MHD or            |
|                                   | hydrodynamic codes                |
+-----------------------------------+-----------------------------------+
| ``--com``                         | if set resets the                 |
|                                   | positions such that the centre of |
|                                   | mass is exactly at the origin.    |
+-----------------------------------+-----------------------------------+
| ``--corotating``                  | velocities are transformed to     |
|                                   | corotating frame                  |
+-----------------------------------+-----------------------------------+
| ``--hfact=1.2``                   | can be changed to give correct    |
|                                   | parameter in                      |
|                                   | :math:`h=h_{fact}(m/\rho)^{1/3}`  |
|                                   | used to set the particle masses   |
|                                   | when rendering minidumps (i.e.,   |
|                                   | when the mass is not dumped).     |
|                                   | Default is RSPLASH_HFACT=1.5      |
+-----------------------------------+-----------------------------------+

ndspmhd data read
~~~~~~~~~~~~~~~~~~

For the ndspmhd read (``splash -ndspmhd``) the options are:

+-----------------------------------+-----------------------------------+
| ``--barycentric``                 | plots barycentric quantities for  |
|                                   | one-fluid dust instead of         |
|                                   | creating fake second set of       |
|                                   | particles                         |
+-----------------------------------+-----------------------------------+

H5Part data read
~~~~~~~~~~~~~~~~~

For the H5PART read (``splash -f h5part``) the options are:

+-----------------------------------+------------------------------------+
| ``--ndim=3``                      | number of spatial dimensions       |
|                                   | :math:`d` (overrides value         |
|                                   | inferred from data)                |
+-----------------------------------+------------------------------------+
| ``--hfac=1.2``                    | factor to use to compute h from    |
|                                   | :math:`h = h_{fac} *(m/\rho)^{1/d}`|
|                                   | if h not present in data           |
+-----------------------------------+------------------------------------+
| ``--hsml=3.0``                    | value for global smoothing length  |
|                                   | h (if h not present in data)       |
+-----------------------------------+------------------------------------+
| ``--typeid=MatID``                | name of the dataset containing     |
|                                   | the particle type identification   |
|                                   | (default is “MatID”)               |
+-----------------------------------+------------------------------------+

.. _sec:envvariables:

Environment variables
---------------------

Several runtime options for splash can be set using environment
variables. These are variables set from your unix shell. In the bash
shell, environment variables are set from the command line using

::

   export VAR='blah'

or by putting this command in your ``.bash_profile``/``.bashrc``. In
csh, the equivalent is

::

   setenv VAR 'blah'

or by putting the above in your ``.cshrc`` file.

Changing the font
~~~~~~~~~~~~~~~~~~

Several environment variables affect the backend plotting library.
Probably the most useful is the ability to change font:

::

   export GIZA_FONT='Helvetica'

where the name is a reasonable guess as to the font you want to use (the
default is ``Times``). In particular, if you are having trouble displaying
unicode characters such as greek letters, you can just change the font
until you find one that works.

Endian changing
~~~~~~~~~~~~~~~~

On some compilers, the endian-ness (byte order) when reading unformatted
binary data files can be changed at runtime. This is useful for looking
at files on different systems to the one on which they were created
(e.g. x86 machines create little-endian files by default, whereas
IBM/powerpc machines create big-endian). Environment variables for
changing the endian-ness of the data read for some common compilers are
given below:

+-------------+----------------------------+----------------+-------------------+----------+
| Compiler    | Environment variable       | Big endian     | Little endian     | Other    |
+=============+============================+================+===================+==========+
| gfortran    | ``GFORTRAN_CONVERT_UNIT``  | ``big_endian`` | ``little_endian`` | ``swap`` |
+-------------+----------------------------+----------------+-------------------+----------+
| ifort       | ``F_UFMTENDIAN``           | ``big``        | ``little``        |          |
+-------------+----------------------------+----------------+-------------------+----------+

For compilers without this feature, almost all can change the
endian-ness at compile time, and the appropriate flags for doing so can
be set using

::

   export ENDIAN='BIG'

or LITTLE before *compiling* splash (this adds the appropriate
compile-time flags for the compiler selected using the SYSTEM
environment variable in the splash Makefile).
