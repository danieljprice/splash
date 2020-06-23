
Getting started
===============

Compiling the code
------------------

The basic steps for installation are as follows:

#. make sure you have a recent Fortran compiler (such as gfortran)

#. compile splash and giza

#. ensure that splash can read your data format

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

Compiling and linking with giza
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can either install giza with your package manager, or in a subdirectory
of splash. To install in a splash subdirectory, use:

::
	cd splash
	git clone http://github.com/danieljprice/giza
   	make withgiza

For detailed instructions on compiling and linking with giza (or the
older pgplot library used in splash v1.x), refer to the INSTALL file in
the root directory of the splash distribution, or at:

::
   http://users.monash.edu.au/~dprice/splash/download/INSTALL.

A successful ``make`` will produce a binary called ``splash``

Reading your data
~~~~~~~~~~~~~~~~~~

The most important part is getting splash to read \*your\* data format.
If you are using a publicly available code, it is reasonably likely
that I have already written a read data subroutine which will read your
dumps. If not it is best to look at some of the other examples and
change the necessary parts to suit your data files. Note that reading
directly from unformatted data files is \*much\* faster than reading
from formatted (ascii) output.

A standard ``make`` will create a binary which supports the file formats listed in
:ref:`tab:defaultreads`. 
All data formats in the splash repository that do not
have an additional dependencies (e.g. ``HDF5``) will be
supported in the splash binary as of version ``3.0.0``.
This means that the user needs to specify the data type
they are reading as a command line option. For example,
the following will read a phantom dumpfile: 
::
	splash --format phantom disc_00000

In some cases, the format of the file can be inferred if
the the file has a known suffix. If splash is compiled with ``HDF5=yes``,
the 
:ref:`tab:otherreads` lists other data reads
implemented but not compiled by default.

.. table:: Binaries and data reads compiled by default
   :name: tab:defaultreads

   +-----------------+-----------------+-----------------+-----------------+
   | splash call     | Formats read    | read_data file  | Comments        |
   +=================+=================+=================+=================+
   | splash          | ascii           | ``read_data_asc | Generic data    |
   |                 |                 | ii.f90``        | read for        |
   |                 |                 |                 | n-column ascii  |
   |                 |                 |                 | formats.        |
   |                 |                 |                 | Automatically   |
   |                 |                 |                 | determines      |
   |                 |                 |                 | number of       |
   |                 |                 |                 | columns and     |
   |                 |                 |                 | skips header    |
   |                 |                 |                 | lines. Can      |
   |                 |                 |                 | recognise SPH   |
   |                 |                 |                 | particle data   |
   |                 |                 |                 | based on the    |
   |                 |                 |                 | column labels.  |
   |                 |                 |                 | Use ‘asplash    |
   |                 |                 |                 | -e’ to plot     |
   |                 |                 |                 | non-SPH data    |
   |                 |                 |                 | (e.g. energy vs |
   |                 |                 |                 | time files).    |
   +-----------------+-----------------+-----------------+-----------------+
   | splash -f dragon| dragon          | ``read_data_dra | see environment |
   |                 |                 | gon.f90``       | variable        |
   |                 |                 |                 | options.        |
   +-----------------+-----------------+-----------------+-----------------+
   | splash -f gadget| gadget,         | ``read_data_gad | Handles both    |
   |                 | gadget-2,       | get.f90``       | default and     |
   |                 | gadget-3        |                 | block-labelled  |
   |                 |                 |                 | formats (see    |
   |                 |                 |                 | environment     |
   |                 |                 |                 | variable        |
   |                 |                 |                 | options).       |
   +-----------------+-----------------+-----------------+-----------------+
   |splash -f ndspmhd| ndspmhd         | ``read_data_dan | Format for the  |
   |                 |                 | sph.f90``       | ndspmhd         |
   |                 |                 |                 | SPH/SPMHD code  |
   |                 |                 |                 | (publicly       |
   |                 |                 |                 | available from  |
   |                 |                 |                 | my website).    |
   +-----------------+-----------------+-----------------+-----------------+
   | splash -f magma | magma           | ``read_data_sro | Stephan         |
   |                 |                 | sph.f90``       | Rosswog’s code  |
   +-----------------+-----------------+-----------------+-----------------+
   |splash -f phantom| sphNG, phantom  | ``read_data_sph | sphNG is        |
   |                 |                 | NG.f90``        | Matthew Bate’s  |
   |                 |                 |                 | SPH code.       |
   +-----------------+-----------------+-----------------+-----------------+
   | splash -f seren | seren           | ``read_data_ser | The SEREN SPH   |
   |                 |                 | en.f90``        | code (Hubber,   |
   |                 |                 |                 | McLeod et al.)  |
   +-----------------+-----------------+-----------------+-----------------+
   | splash -f gaseoline| gasoline, tipsy | ``read_data_tip | Reads both      |
   |                 |                 | sy.f90``        | binary and      |
   |                 |                 |                 | ascii TIPSY     |
   |                 |                 |                 | files           |
   |                 |                 |                 | (determined     |
   |                 |                 |                 | automatically). |
   +-----------------+-----------------+-----------------+-----------------+
   | splash -f vine  | vine            | ``read_data_VIN | see environment |
   |                 |                 | E.f90``         | variable        |
   |                 |                 |                 | options.        |
   +-----------------+-----------------+-----------------+-----------------+

.. table:: Other data reads implemented but not compiled by default
   :name: tab:otherreads

   +-----------------+-----------------+-----------------+-----------------+
   | Format          | Binary          | read_data file  | Comments        |
   +=================+=================+=================+=================+
   | h5part          | h5splash        | ``read_data_h5p | Reads general   |
   |                 |                 | art.f90``       | files written   |
   |                 |                 |                 | with the h5part |
   |                 |                 |                 | library.        |
   |                 |                 |                 | Requires        |
   |                 |                 |                 | linking against |
   |                 |                 |                 | H5PART and HDF5 |
   |                 |                 |                 | libraries       |
   +-----------------+-----------------+-----------------+-----------------+
   | gadget HDF5     | gsplash-hdf5    | ``read_data_gad | Reads HDF5      |
   |                 |                 | get_hdf5.f90``  | format from the |
   |                 |                 |                 | gadget code.    |
   |                 |                 |                 | Requires        |
   |                 |                 |                 | linking against |
   |                 |                 |                 | HDF5 libraries  |
   +-----------------+-----------------+-----------------+-----------------+
   | amuse HDF5      | amsplash-hdf5   | ``read_data_amu | Reads HDF5      |
   |                 |                 | se_hdf5.f90``   | format from the |
   |                 |                 |                 | amuse           |
   |                 |                 |                 | framework.      |
   +-----------------+-----------------+-----------------+-----------------+
   | ``.silo``       | silosplash      | ``read_data_sil | a nice          |
   | format          |                 | o.f90``         | standardised    |
   | (particle data  |                 |                 | HDF5 particle   |
   | only)           |                 |                 | format.         |
   |                 |                 |                 | Requires silo   |
   |                 |                 |                 | libraries.      |
   +-----------------+-----------------+-----------------+-----------------+
   | SNSPH           | snsplash        | ``read_data_sns | Supernova SPH   |
   |                 |                 | ph.f90``        | (Chris Fryer et |
   |                 |                 |                 | al.). Requires  |
   |                 |                 |                 | libsw.          |
   +-----------------+-----------------+-----------------+-----------------+
   | falcON          | fsplash         | ``read_data_fal | Walter Dehnen’s |
   |                 |                 | cON.f90``       | SPH code format |
   |                 |                 |                 | (uses HDF5)     |
   +-----------------+-----------------+-----------------+-----------------+
   | Andreas         | bsplash         | ``read_data_bau |                 |
   | Bauswein’s code |                 | swein.f90``     |                 |
   +-----------------+-----------------+-----------------+-----------------+
   | Sigfried        | vsplash         | ``read_data_van |                 |
   | Vanaverbeke’s   |                 | averbeke.f90``  |                 |
   | code            |                 |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   | Regularised SPH | rsplash         | ``read_data_rsp |                 |
   | (Steinar Børve) |                 | h.f90``         |                 |
   +-----------------+-----------------+-----------------+-----------------+
   | FLASH tracer    | fsplash         | ``read_data_fla | Reads tracer    |
   | particles       |                 | sh_hdf5.f90``   | particle output |
   |                 |                 |                 | from the FLASH  |
   |                 |                 |                 | code. Requires  |
   |                 |                 |                 | linking against |
   |                 |                 |                 | HDF5 libraries  |
   +-----------------+-----------------+-----------------+-----------------+
   | Sky King/Nikos  | usplash         | ``read_data_UCL | A good example  |
   | Mastrodemos     |                 | A.f90``         | of a simple     |
   |                 |                 |                 | ascii format    |
   |                 |                 |                 | reader          |
   +-----------------+-----------------+-----------------+-----------------+
   | Jamie Bolton    | gsplash_jsb     | ``read_data_gad | Reads extra     |
   | GADGET          |                 | get_jsb.f90``   | arrays before   |
   |                 |                 |                 | the SPH         |
   |                 |                 |                 | smoothing       |
   |                 |                 |                 | length          |
   +-----------------+-----------------+-----------------+-----------------+
   | Old Matthew     | bsplash         | ``read_data_mba | similar to the  |
   | Bate code       |                 | te.f90``        | original Benz   |
   |                 |                 |                 | SPH code format |
   +-----------------+-----------------+-----------------+-----------------+
   | Foulkes/Haswell | fsplash         | ``read_data_fou | An ascii format |
   | /Murray         |                 | lkes.f90``      |                 |
   +-----------------+-----------------+-----------------+-----------------+
   | Andrea Urban    | usplash         | ``read_data_urb | An ascii format |
   | format          |                 | an.f90``        |                 |
   +-----------------+-----------------+-----------------+-----------------+
   | ``.pbob``       | psplash         | ``read_data_pbo | David Brown’s   |
   | format          |                 | b.f90``         | SPH code. To make available, you must specify |
   |                 |                 |                 | `PBOB_DIR=\path\to\PBOB\lib`
   +-----------------+-----------------+-----------------+-----------------+

Further details on writing your own subroutine are given in
appendix :ref:`sec:writeyourown`. The \*easiest\* way is to i)
email me a sample data file and ii) the subroutine you used to write it,
and I will happily create a data read for your file format.

.. _sec:commandline:

Command line options
--------------------

Typing ``splash -v`` gives a complete and up-to-date list of options. Currently these are:

::

   Command line options:

    -p fileprefix     : change prefix to ALL settings files read/written by splash
    -d defaultsfile   : change name of defaults file read/written by splash
    -l limitsfile     : change name of limits file read/written by splash
    -e, -ev           : use default options best suited to ascii evolution files (ie. energy vs time)
    -lm, -lowmem      : use low memory mode [applies only to sphNG data read at present]
    -o pixformat      : dump pixel map in specified format (use just -o for list of formats)

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
