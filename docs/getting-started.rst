
Getting started
===============

.. _install:

Installing from your package manager
-------------------------------------
Pre-packaged builds of splash exist for most operating systems.

Stable version
~~~~~~~~~~~~~~
Mac OS via homebrew::

  brew tap danieljprice/all
  brew install splash

Mac OS via Macports::

  sudo port install splash

Linux or Windows Linux Subsystem (Ubuntu)::

  sudo apt-get install splash

Development version
~~~~~~~~~~~~~~~~~~~

SPLASH and giza (the plotting backend) both have public repositories, so you can check out the latest and greatest code at any time. Just use:

Mac OS via homebrew::

  brew tap danieljprice/all
  brew install --HEAD splash

Compiling from source
---------------------

Basic instructions
~~~~~~~~~~~~~~~~~~
If you have admin (super user) permissions::

   git clone https://github.com/danieljprice/giza.git; cd giza; ./configure; make; sudo make install; cd ..
   git clone https://github.com/danieljprice/splash.git
   cd splash; make SYSTEM=gfortran; sudo make install

If you do not have admin permissions. That is, to install in your home space::

   git clone https://github.com/danieljprice/splash.git
   cd splash; git clone https://github.com/danieljprice/giza.git
   make SYSTEM=gfortran withgiza

Both the splash and giza repositories are generally stable, so it is usually safe to get the very latest version.

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
of splash. To install in a splash subdirectory, use
::
	cd splash
	git clone http://github.com/danieljprice/giza
   	make withgiza

of splash.

With giza installed via your package manager (or previously compiled as below), use
::
   cd splash
   make GIZA_DIR=/usr/local

where ``GIZA_DIR`` points to the directory where giza was installed.
To install giza in a splash subdirectory, use
::
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
All data formats in the splash repository that do not
have an additional dependencies (e.g. ``HDF5``) will be
supported in the splash binary as of version ``3.0.0``.
This means that the user needs to specify the data type
they are reading as a command line option. For example,
the following will read a phantom dumpfile
::
	splash --format phantom disc_00000

In some cases, the format of the file can be inferred if
the the file has a known suffix. For example, the above line can be changed if the
suffixe of the file is recognised
::
	splash disc_00000.pb
This will automatically recognise a Phantom binary dumpfile. For backwards compatibility with
previous version of ``splash``, one can add aliases into their `.bashrc`, or equivalent
::
 	alias asplash='splash ' # Alias for ascii splash
 	alias ssplash='splash -f phantom '
 	alias gsplash='splash -f gadget '
 	alias vsplash='splash -f vine '
 	alias nsplash='splash -f ndspmhd '
 	alias rsplash='splash -f srosph '
 	alias dsplash='splash -f dragon '
 	alias srsplash='splash -f seren '
 	alias tsplash='splash -f tipsy '
 	alias tsplash='splash -f tipsy '
 	alias msplash='splash -f mhutch '

If splash is compiled with ``HDF5=yes``,
the formats listed in 
:ref:`tab:hdf5reads` will also be available in the ``splash`` binary.
 Other supported formats are listed in
:ref:`tab:otherreads`, but these require additional libraries.

.. table:: Binaries and data reads compiled by default
   :name: tab:defaultreads
+------------------------------+----------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``splash`` command           | Format Read                | ``read_data`` File            | Comments                                                                                                                                                                                                                                         |
+==============================+============================+===============================+==================================================================================================================================================================================================================================================+
| ``splash -gadget <file>``    | ascii                      | ``read_data_asci.f90``        | Generic data read for n-column ascii formats. Automatically determines number of columns and skips header lines. Can recognise SPH particle data based on the column labels. Use ``splash -e`` to plot non-SPH data (e.g. energy vs time files). |
+------------------------------+----------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``splash -dragon <file>``    | dragon                     | ``read_data_dragon``          | See environment variable options.                                                                                                                                                                                                                |
+------------------------------+----------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``splash -gadget <file>``    | gadget, gadget-2, gadget-3 | ``read_data_gadget.f90``      | Handles both default and block-labelled formats (see environment variable options).                                                                                                                                                              |
+------------------------------+----------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``splash -ndspmhd <file>``   | ndspmhd                    | ``read_data_ndspmhd.f90``     | Format for the ndspmhd SPH/SPMHD code (publicly available from my website).                                                                                                                                                                      |
+------------------------------+----------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``splash -phantom <file>``   | sphNG, Phantom             | ``read_data_sphNG.f90``       | sphNG is Matthew Bate’s SPH code. Option ``-sphng``also works.                                                                                                                                                                                   |
+------------------------------+----------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``splash -magma <file>``     | magma                      | ``read_data_srosph.f90``      | Stephan Rosswog’s code                                                                                                                                                                                                                           |
+------------------------------+----------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``splash -seren <file>``     | seren                      | ``read_data_seren.f90``       | The SEREN SPH code (Hubber, McLeod et al.)                                                                                                                                                                                                       |
+------------------------------+----------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``splash -gasoline <file>``  | gasoline, tipsy            | ``read_data_tipsy.f90``       | Reads both binary and ascii TIPSY files (determined automatically). Option ``-tipsy`` also works.                                                                                                                                                |
+------------------------------+----------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``splash -vine <file>``      | vine                       | ``read_data_fine.f90``        | See environment variable options.                                                                                                                                                                                                                |
+------------------------------+----------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``splash -starsmasher <file>``| StarSmasher                | ``read_data_starsmasher.f90`` | The StarSmasher code (Gaburov et al. 2018) `<jalombar.github.io/starsmasher/>`_                                                                                                                                                                  |
+------------------------------+----------------------------+-------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. table:: Supported HDf5 data formats
   :name: tab:hdf5reads
+--------------------------------+------------------------+-------------------------------+-----------------------------------------------------------------------------------------+
| ``splash`` Command             | Read Format            | ``read_data`` File            | Comments                                                                                |
+================================+========================+===============================+=========================================================================================+
| ``splash -gadget_hdf5 <file>`` | gadget HDF5 Files.     | ``read_data_gadget_hdf5.f90`` | Reads HDF5 format from the gadget code.                                                 |
+--------------------------------+------------------------+-------------------------------+-----------------------------------------------------------------------------------------+
| ``splash -amuse <file>``       | AMUSE HDF5             | ``read_data_amuse_hdf5.f90``  | Reads HDF5 format from the AMUSE framework.                                             |
+--------------------------------+------------------------+-------------------------------+-----------------------------------------------------------------------------------------+
| ``splash -cactus_hdf5 <file>`` | Cactus HDF5            | ``read_data_cactus_hdf5.f90`` |                                                                                         |
+--------------------------------+------------------------+-------------------------------+-----------------------------------------------------------------------------------------+
| ``splash -flash_hdf5 <file>    | FLASH tracer particles | ``read_dataflash_hdf5.f90``   | Reads tracer particle output from the FLASH code. The option ``-flash`` will also work. |
+--------------------------------+------------------------+-------------------------------+-----------------------------------------------------------------------------------------+
| ``splash -falcon_hdf5 <file>   | falcON                 | ``read_data_falcON_hdf5.f90`` | Walter Dehnen’s SPH code format. The option ``-falcon`` will also work.                 |
+--------------------------------+------------------------+-------------------------------+-----------------------------------------------------------------------------------------+

If the ``HDF5`` read files end with ``.h5``, the suffix ``_hdf5`` from the ``splash`` command can be removed.
For example, 
::
	splash -gadget dump_000.h5
will recognise that the file ``dump_000.h5`` is in the ``HDF5`` format,
and will automatically select the correct ``read_data`` routine.

.. table:: Other supported file formats that require external libraries
   :name: tab:otherreads



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

    -p fileprefix     : change prefix to ALL settings files read/written by splash
    -d defaultsfile   : change name of defaults file read/written by splash
    -l limitsfile     : change name of limits file read/written by splash
    -e, -ev           : use default options best suited to ascii evolution files (ie. energy vs time)
    -lm, -lowmem      : use low memory mode [applies only to sphNG data read at present]
    -o pixformat      : dump pixel map in specified format (use just -o for list of formats)
    -f                : input file format to be read (ascii is default)

   To select data formats, use the shortcuts below, or use the -f or --format command line options
   Multiple data formats are not support in a single instance.
   Supported data formats:
    -ascii            : ascii file format (default)
    -phantom -sphng   : Phantom and sphNG codes
    -ndspmhd          : ndsphmd code
    -gadget           : Gadget code
    -seren            : Seren code
   ..plus many others. Type --formats for a full list
  
   The following formats support HDF5:
    -flash            : FLASH code
    -gadget           : Gadget code
    -cactus           : Cactus SPH code
    -falcon           : FalcON code
    -amuse            : AMUSE Framework
  
   add a suffix "_hdf5" to the above command line options if your data files do not end with .h5.

   Command line plotting mode:

    -x column         : specify x plot on command line (ie. do not prompt for x)
    -y column         : specify y plot on command line (ie. do not prompt for y)
    -r[ender] column  : specify rendered quantity on command line (ie. no render prompt)
                        (will take columns 1 and 2 as x and y if -x and/or -y not specified)
    -vec[tor] column  : specify vector plot quantity on command line (ie. no vector prompt)
    -c[ontour] column : specify contoured quantity on command line (ie. no contour prompt)
    -dev device       : specify plotting device on command line (ie. do not prompt)

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

     the above options all produce a small ascii file with one row per input file.
     the following option produces a file equivalent in size to one input file (in ascii format):

            calc timeaverage  : time average of *all* entries for every particle
                                output to file called 'time_average.out'

            calc ratio        : ratio of *all* entries in each file compared to first
                                output to file called 'ratio.out'

Command-line options can be entered in any order on the command line
(even after the dump file names). For more information on the convert
utility (``splash to ascii``) see :ref:`sec:convert`. For details
of the ``-o ppm`` or ``-o ascii`` option see :ref:`sec:writepixmap`. For details of the ``-ev`` option, see :ref:`sec:evsplash`.
