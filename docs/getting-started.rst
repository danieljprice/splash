
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

A copy of giza is included in the splash distribution and is compiled
automatically along with splash. giza is also available as a standalone
project at:

   http://github.com/danieljprice/giza

For detailed instructions on compiling and linking with giza (or the
older pgplot library used in splash v1.x), refer to the INSTALL file in
the root directory of the splash distribution, or at:

   http://users.monash.edu.au/~dprice/splash/download/INSTALL.

A successful ‘make’ will produce a binary called `splash`

Reading your data
~~~~~~~~~~~~~~~~~~

The most important part is getting splash to read \*your\* data format.
If you are using a publically available code, it is reasonably likely
that I have already written a read data subroutine which will read your
dumps. If not it is best to look at some of the other examples and
change the necessary parts to suit your data files. Note that reading
directly from unformatted data files is \*much\* faster than reading
from formatted (ascii) output.

A standard “make” will create the binaries listed in
Table :ref:`tab:defaultreads` which read the
corresponding data formats listed in the third column.
Table :ref:`tab:otherreads` lists other data reads
implemented but not compiled by default.

.. table:: Binaries and data reads compiled by default
   :name: tab:defaultreads

   +-----------------+-----------------+-----------------+-----------------+
   | splash binary   | Formats read    | read_data file  | Comments        |
   +=================+=================+=================+=================+
   | asplash, splash | ascii           | ``read_data_asc | Generic data    |
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
   | dsplash         | dragon          | ``read_data_dra | see environment |
   |                 |                 | gon.f90``       | variable        |
   |                 |                 |                 | options.        |
   +-----------------+-----------------+-----------------+-----------------+
   | gsplash         | gadget,         | ``read_data_gad | Handles both    |
   |                 | gadget-2,       | get.f90``       | default and     |
   |                 | gadget-3        |                 | block-labelled  |
   |                 |                 |                 | formats (see    |
   |                 |                 |                 | environment     |
   |                 |                 |                 | variable        |
   |                 |                 |                 | options).       |
   +-----------------+-----------------+-----------------+-----------------+
   | nsplash         | ndspmhd         | ``read_data_dan | Format for the  |
   |                 |                 | sph.f90``       | ndspmhd         |
   |                 |                 |                 | SPH/SPMHD code  |
   |                 |                 |                 | (publicly       |
   |                 |                 |                 | available from  |
   |                 |                 |                 | my website).    |
   +-----------------+-----------------+-----------------+-----------------+
   | rsplash         | magma           | ``read_data_sro | Stephan         |
   |                 |                 | sph.f90``       | Rosswog’s code  |
   +-----------------+-----------------+-----------------+-----------------+
   | ssplash         | sphNG, phantom  | ``read_data_sph | sphNG is        |
   |                 |                 | NG.f90``        | Matthew Bate’s  |
   |                 |                 |                 | SPH code.       |
   +-----------------+-----------------+-----------------+-----------------+
   | srsplash        | seren           | ``read_data_ser | The SEREN SPH   |
   |                 |                 | en.f90``        | code (Hubber,   |
   |                 |                 |                 | McLeod et al.)  |
   +-----------------+-----------------+-----------------+-----------------+
   | tsplash         | gasoline, tipsy | ``read_data_tip | Reads both      |
   |                 |                 | sy.f90``        | binary and      |
   |                 |                 |                 | ascii TIPSY     |
   |                 |                 |                 | files           |
   |                 |                 |                 | (determined     |
   |                 |                 |                 | automatically). |
   +-----------------+-----------------+-----------------+-----------------+
   | vsplash         | vine            | ``read_data_VIN | see environment |
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
   | format          |                 | b.f90``         | SPH code        |
   +-----------------+-----------------+-----------------+-----------------+

Further details on writing your own subroutine are given in
appendix :ref:`sec:writeyourown`. The \*easiest\* way is to i)
email me a sample data file and ii) the subroutine you used to write it,
and I will happily create a data read for your file format.

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
default is ‘Times’). In particular, if you are having trouble displaying
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
| ``gfortran`` | ``GFORTRAN_CONVERT_UNIT`` | ``big_endian`` | ``little_endian`` | ``swap`` |
+-------------+----------------------------+----------------+-------------------+----------+
| ``ifort``   | ``F_UFMTENDIAN``           | ``big``        | ``little``        |          |
+-------------+----------------------------+----------------+-------------------+----------+

For compilers without this feature, almost all can change the
endian-ness at compile time, and the appropriate flags for doing so can
be set using

::

   export ENDIAN='BIG'

or LITTLE before *compiling* splash (this adds the appropriate
compile-time flags for the compiler selected using the SYSTEM
environment variable in the splash Makefile).

Variables affecting all data reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Environment variables that affect all data reads are:

+-----------------------------------+-----------------------------------+
| SPLASH_DEFAULTS                   | gives the name of a system-wide   |
|                                   | ``splash.defaults`` file (and     |
|                                   | splash.limits etc.) that will be  |
|                                   | used if there is none in the      |
|                                   | current directory. e.g.           |
|                                   | ``export SPLASH_DEFAULTS=/home/me |
|                                   | /splash.defaults``                |
+-----------------------------------+-----------------------------------+
| SPLASH_KERNEL                     | changes the smoothing kernel used |
|                                   | in the interpolations (e.g.       |
|                                   | ‘cubic’ or ‘quintic’). Can also   |
|                                   | be changed in the r)ender menu.   |
+-----------------------------------+-----------------------------------+
| SPLASH_DEBUG                      | if set to ‘yes’ or ‘true’, turns  |
|                                   | on very verbose debugging output. |
|                                   | Useful to trace code crashes (but |
|                                   | of course, this never happens…).  |
+-----------------------------------+-----------------------------------+
| SPLASH_CENTRE_ON_SINK             | if set to a number n, centres     |
|                                   | coordinates and velocities on the |
|                                   | nth sink/star particle (e.g.      |
|                                   | ``export SPLASH_CENTRE_ON_SINK=2``|
|                                   | ).                                |
+-----------------------------------+-----------------------------------+
| SPLASH_COROTATE                   | plot in corotating frame based on |
|                                   | locations of two sink particles   |
+-----------------------------------+-----------------------------------+
| SPLASH_HMIN_CODEUNITS             | if given a value :math:`>`\ 0     |
|                                   | enforces a minimum smoothing      |
|                                   | length, specified in code units   |
|                                   | as read from the dump file, on    |
|                                   | all the particles. This can be    |
|                                   | used to “dumb-down” the           |
|                                   | resolution of SPH simulations,    |
|                                   | e.g. to match observational       |
|                                   | resolution. If this variable is   |
|                                   | set it is *highly* recommended    |
|                                   | that the “use accelerated         |
|                                   | rendering” option in the r)ender  |
|                                   | menu is also turned on as quite   |
|                                   | slow rendering can otherwise      |
|                                   | result.                           |
+-----------------------------------+-----------------------------------+
| SPLASH_VZERO_CODEUNITS            | if set to a comma separated list  |
|                                   | of vector components (e.g.        |
|                                   | ``export SPLASH_VZERO_CODEUNITS=' |
|                                   | 0.0,1.0,0.0'``),                  |
|                                   | can be used to subtract a mean    |
|                                   | velocity field from all particles |
|                                   | — specified in code units as read |
|                                   | from the dump file.               |
+-----------------------------------+-----------------------------------+
| SPLASH_MARGIN_XMIN                | can be used to manually adjust    |
|                                   | the left horizontal page margin   |
|                                   | (set to fraction of viewport,     |
|                                   | negative values are allowed).     |
+-----------------------------------+-----------------------------------+
| SPLASH_MARGIN_XMAX                | right horizontal page margin (set |
|                                   | to fraction of viewport).         |
+-----------------------------------+-----------------------------------+
| SPLASH_MARGIN_YMIN                | bottom (vertical) page margin     |
|                                   | (set to fraction of viewport).    |
+-----------------------------------+-----------------------------------+
| SPLASH_MARGIN_YMAX                | top (vertical) page margin (set   |
|                                   | to fraction of viewport).         |
+-----------------------------------+-----------------------------------+

.. _sec:asplash:

Ascii data read
~~~~~~~~~~~~~~~~

For several data reads there are environment variables which can be set
at runtime which are specific to the data read. For the ascii data read
(‘asplash’) these are:

+-----------------------------------+-----------------------------------+
| ASPLASH_NCOLUMNS                  | if given a value :math:`>`\ 0     |
|                                   | sets the number of columns to be  |
|                                   | read from ascii data (overrides   |
|                                   | the automatic number of columns   |
|                                   | determination).                   |
+-----------------------------------+-----------------------------------+
| ASPLASH_NHEADERLINES              | if given a value :math:`>=`\ 0    |
|                                   | sets the number of header lines   |
|                                   | to skip (overrides the automatic  |
|                                   | determination).                   |
+-----------------------------------+-----------------------------------+
| ASPLASH_COLUMNSFILE               | can be used to provide the        |
|                                   | location of (path to) the default |
|                                   | ‘columns’ file containing the     |
|                                   | labels for ascii data (e.g.       |
|                                   | setenv ASPLASH_COLUMNSFILE        |
|                                   | ’/home/me/mylabels’). Overridden  |
|                                   | by the presence of a local        |
|                                   | ‘columns’ file.                   |
+-----------------------------------+-----------------------------------+
| ASPLASH_TIMEVAL                   | if given a nonzero value sets the |
|                                   | time to use in the legend (fixed  |
|                                   | for all files)                    |
+-----------------------------------+-----------------------------------+
| ASPLASH_GAMMAVAL                  | if given a nonzero value sets     |
|                                   | gamma to use in exact solution    |
|                                   | calculations (fixed for all       |
|                                   | files)                            |
+-----------------------------------+-----------------------------------+
| ASPLASH_HEADERLINE_TIME           | sets the integer line number      |
|                                   | where the time appears in the     |
|                                   | header                            |
+-----------------------------------+-----------------------------------+
| ASPLASH_HEADERLINE_GAMMA          | sets the integer line number      |
|                                   | where gamma appears in the header |
+-----------------------------------+-----------------------------------+

.. _sec:gsplash:

GADGET data read
~~~~~~~~~~~~~~~~~

For the GADGET read (‘gsplash’) the environment variable options are:

+-----------------------------------+-----------------------------------+
| GSPLASH_FORMAT                    | if set = 2, reads the block       |
|                                   | labelled GADGET format instead of |
|                                   | the default (non block labelled)  |
|                                   | format.                           |
+-----------------------------------+-----------------------------------+
| GSPLASH_USE_Z                     | if ‘YES’ or ‘TRUE’ uses the       |
|                                   | redshift in the legend instead of |
|                                   | code time.                        |
+-----------------------------------+-----------------------------------+
| GSPLASH_DARKMATTER_HSOFT          | if given a value :math:`>` 0.0    |
|                                   | will assign a smoothing length to |
|                                   | dark matter particles for which   |
|                                   | rendered plots of column density  |
|                                   | can then be made.                 |
+-----------------------------------+-----------------------------------+
| GSPLASH_EXTRACOLS                 | if set to a comma separated list  |
|                                   | of column labels, will attempt to |
|                                   | read additional columns           |
|                                   | containing gas particle           |
|                                   | properties beyond the end of the  |
|                                   | file (not applicable if           |
|                                   | GSPLASH_FORMAT=2).                |
+-----------------------------------+-----------------------------------+
| GSPLASH_STARPARTCOLS              | if set to a comma separated list  |
|                                   | of column labels, will attempt to |
|                                   | read additional columns           |
|                                   | containing star particle          |
|                                   | properties beyond the end of the  |
|                                   | file (and after any extra gas     |
|                                   | particle columns) (not applicable |
|                                   | if GSPLASH_FORMAT=2).             |
+-----------------------------------+-----------------------------------+
| GSPLASH_CHECKIDS                  | if set to ‘YES’ or ‘TRUE’, reads  |
|                                   | and checks particle IDs,          |
|                                   | excluding particles with negative |
|                                   | IDs as accreted (gives them a     |
|                                   | negative smoothing length which   |
|                                   | means they are ignored in         |
|                                   | renderings).                      |
+-----------------------------------+-----------------------------------+
| GSPLASH_HSML_COLUMN               | if set to a positive integer,     |
|                                   | specifies the location of the     |
|                                   | smoothing length in the columns,  |
|                                   | overriding any default settings.  |
+-----------------------------------+-----------------------------------+
| GSPLASH_IGNORE_IFLAGCOOL          | if set to ’YES’ or ‘TRUE’, does   |
|                                   | not assume that extra columns are |
|                                   | present even if the cooling flag  |
|                                   | is set in the header.             |
+-----------------------------------+-----------------------------------+

For the GADGET read gsplash will also look for, and read if present,
files called ``snapshot_xxx.hsml`` and/or ``snapshot_xxx.dens`` (where
``snapshot_xxx`` is the name of the corresponding GADGET dump file)
which contain smoothing lengths and/or a density estimate for dark
matter particles (these should just be one-column ascii files).

VINE data read
~~~~~~~~~~~~~~~

For the VINE read (‘vsplash’) the environment variable options are:

+-----------------------------------+-----------------------------------+
| VSPLASH_HFAC                      | if ‘YES’ or ‘TRUE’ multiplies the |
|                                   | smoothing length read from the    |
|                                   | dump file by a factor of 2.8 (for |
|                                   | use with older VINE dumps where   |
|                                   | the smoothing length is defined   |
|                                   | as in a Plummer kernel rather     |
|                                   | than as the usual SPH smoothing   |
|                                   | length).                          |
+-----------------------------------+-----------------------------------+
| VSPLASH_MHD                       | if ‘YES’ or ‘TRUE’ reads VINE     |
|                                   | dumps containing MHD arrays (note |
|                                   | that setting VINE_MHD also        |
|                                   | works).                           |
+-----------------------------------+-----------------------------------+

sphNG data read
~~~~~~~~~~~~~~~~

For the sphNG and PHANTOM read (‘ssplash’) the environment variable
options are:

+-----------------------------------+-----------------------------------+
| SSPLASH_RESET_CM                  | if ‘YES’ or ‘TRUE’ resets the     |
|                                   | positions such that the centre of |
|                                   | mass is exactly at the origin.    |
+-----------------------------------+-----------------------------------+
| SSPLASH_OMEGA                     | if non-zero, subtracts solid body |
|                                   | rotation with omega as specified  |
|                                   | to give velocities in co-rotating |
|                                   | frame.                            |
+-----------------------------------+-----------------------------------+
| SSPLASH_OMEGAT                    | if non-zero, subtracts solid body |
|                                   | rotation with omega as specified  |
|                                   | to give positions and velocities  |
|                                   | in co-rotating frame.             |
+-----------------------------------+-----------------------------------+
| SSPLASH_TIMEUNITS                 | sets default time units, either   |
|                                   | ’s’, ’min’, ’hrs’, ’days’, ’yrs’  |
|                                   | or ’tfreefall’ (NB: text is used  |
|                                   | verbatim in legend).              |
+-----------------------------------+-----------------------------------+

dragon data read
~~~~~~~~~~~~~~~~~

For the dragon read (‘dsplash’) the environment variable options are:

+-----------------------------------+-----------------------------------+
| DSPLASH_EXTRACOLS                 | specifies number of extra columns |
|                                   | present in the file which are     |
|                                   | dumped after the itype array      |
+-----------------------------------+-----------------------------------+

Stephan Rosswog data read
~~~~~~~~~~~~~~~~~~~~~~~~~~

For the srosph read (‘rsplash’) the environment variable options are:

+-----------------------------------+-----------------------------------+
| RSPLASH_FORMAT                    | can be ‘MHD’ or ‘HYDRO’ which     |
|                                   | read the appropriate data format  |
|                                   | from either the MHD or            |
|                                   | hydrodynamic codes                |
+-----------------------------------+-----------------------------------+
| RSPLASH_RESET_COM                 | if ‘YES’ or ‘TRUE’ resets the     |
|                                   | positions such that the centre of |
|                                   | mass is exactly at the origin.    |
+-----------------------------------+-----------------------------------+
| RSPLASH_COROTATING                | if ‘YES’ or ‘TRUE’ then           |
|                                   | velocities are transformed to     |
|                                   | corotating frame                  |
+-----------------------------------+-----------------------------------+
| RSPLASH_HFACT                     | can be changed to give correct    |
|                                   | parameter in                      |
|                                   | :math:`h=h_{fact}(m/\rho)^{1/3}`  |
|                                   | used to set the particle masses   |
|                                   | when rendering minidumps (i.e.,   |
|                                   | when the mass is not dumped).     |
|                                   | Default is RSPLASH_HFACT=1.5      |
+-----------------------------------+-----------------------------------+

ndspmhd data read
~~~~~~~~~~~~~~~~~~

For the ndspmhd read (‘nsplash’) the environment variable options are:

+-----------------------------------+-----------------------------------+
| NSPLASH_BARYCENTRIC               | plots barycentric quantities for  |
|                                   | one-fluid dust instead of         |
|                                   | creating fake second set of       |
|                                   | particles                         |
+-----------------------------------+-----------------------------------+

H5Part data read
~~~~~~~~~~~~~~~~~

For the H5PART read (‘h5splash’) the environment variable options are:

+-----------------------------------+------------------------------------+
| H5SPLASH_NDIM                     | number of spatial dimensions       |
|                                   | :math:`d` (overrides value         |
|                                   | inferred from data)                |
+-----------------------------------+------------------------------------+
| H5SPLASH_HFAC                     | factor to use to compute h from    |
|                                   | :math:`h = h_{fac} *(m/\rho)^{1/d}`|
|                                   | if h not present in data           |
+-----------------------------------+------------------------------------+
| H5SPLASH_HSML                     | value for global smoothing length  |
|                                   | h (if h not present in data)       |
+-----------------------------------+------------------------------------+
| H5SPLASH_TYPEID                   | name of the dataset containing     |
|                                   | the particle type identification   |
|                                   | (default is “MatID”)               |
+-----------------------------------+------------------------------------+

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
utility (‘splash to ascii’) see :ref:`sec:convert`. For details
of the ``-o ppm`` or ``-o ascii`` option see :ref:`sec:writepixmap`. For details of the ``-ev`` option, see :ref:`sec:evsplash`.
