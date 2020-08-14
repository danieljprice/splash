
.. _sec:writeyourown:

Data reads and environment variables
=====================================

Requesting a customised data reader
------------------------------------
The basic "splash" binary will read any ascii or csv data file where
columns correspond to different quantities and rows correspond to each particle (actually
I use splash to plot graphs for nearly all data in this form, whether SPH or not)
-- it will also sensibly skip header lines which do not have the same number of columns.

However, it is ultimately desirable to use SPLASH to directly visualise the
(binary) output of your code. If your format is not amongst those distributed,
then BEFORE you start writing your own routine, please consider whether or not a routine
to read your format would be of more general use (e.g. to other users of your code).
If so, PLEASE email me to request a new read_data routine for your format, by sending an email attaching:

a) an example dump

and

b) the source code from the routine which wrote the dump file.

Then I can write a read for your format that can be added to the SPLASH repository
and distributed in all future versions. Whilst I aim never to change the interface
to the read_data routines, it is not impossible that some changes may occur
somewhere down the line (or enhanced functionality -- for example the more advanced
data reads are able to read only the required columns for a given plot from the
file, rather than the whole file). It doesn’t matter if your code only has one user,
we are still happy to do this as it makes splash more widely useable and
saves trouble later.

Writing your own data read
---------------------------
If you *really* want to hack one yourself it is best to look at some of the
other examples and change the  necessary parts to suit your data files. Note
that reading directly from unformatted data files is *much* faster than reading
from formatted (ascii) data.

If you do end up writing your own, again, please email me the end result so I
can add it to the officially supported data reads. This also makes it much
easier for you to upgrade to newer versions as you do not require a locally
customised version.

The second best way is to attempt to modify one of the existing data
reads. Even then, there are some things to note: Most important is that,
for the rendering routines to work, the density, particle masses and
smoothing lengths for *all* of the (gas) particles *must* be read in
from the data file and their locations in the main data array labelled
using the integer parameters ``irho``, ``ipmass`` and ``ih``. Labelling
of the location of other particle quantities (e.g. ``iutherm`` for the
thermal energy) is used in order to plot the exact solutions on the
appropriate graphs and also for calculating additional quantities (e.g.
calculation of the pressure uses ``iutherm`` and ``irho``).

The positions of vector components in the data columns are indicated by
setting the variable ``iamvec`` of that column equal to the first
component of the vector of which this component is a part. So if column
4 is a vector quantity (say :math:`{\bf v}` in 3D), then
``iamvec(4) = 4``, ``iamvec(5) = 4`` and ``iamvec(6) = 4``. Similarly
the string ``labelvec`` should be set, i.e., ``labelvec = 'v'`` for
these columns.


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

Variables affecting all data reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Environment variables that affect all data reads are:

+----------------------------+-------------------------------------------------+
| SPLASH_DEFAULTS            | gives name of system-wide ``splash.defaults``   |
|                            | file (and splash.limits etc.) that will be      |
|                            | used if there is none in the current dir. e.g.  |
|                            | ``export SPLASH_DEFAULTS=$HOME/splash.defaults``|
+----------------------------+-------------------------------------------------+
| SPLASH_KERNEL              | changes the smoothing kernel used in the        |
|                            | interpolations (e.g. ``cubic`` or ``quintic``). |
|                            | Can also be changed in the :ref:`sec:menu-r`.   |
+----------------------------+-------------------------------------------------+
| SPLASH_DEBUG               | if set to ``yes`` or ``true``, turns on verbose |
|                            | debugging output. Useful to trace code crashes  |
|                            | (but of course, this never happens…).           |
+----------------------------+-------------------------------------------------+
| SPLASH_CENTRE_ON_SINK      | if set to a number n, centres coordinates and   |
|                            | velocities on the nth sink/star particle (e.g.  |
|                            | ``export SPLASH_CENTRE_ON_SINK=2``)             |
+----------------------------+-------------------------------------------------+
| SPLASH_COROTATE            | plot in corotating frame based on locations of  |
|                            | 2 sink particles (e.g. ``SPLASH_CORORATE=1,3``) |
+----------------------------+-------------------------------------------------+
| SPLASH_HMIN_CODEUNITS      | if given a value :math:`>`\ 0 enforces a minimum|
|                            | smoothing length, specified in code units,      |
|                            | all the particles. Useful to “dumb-down” the    |
|                            | resolution of SPH simulations to match          |
|                            | observational resolution. If this variable is   |
|                            | set the “accelerated rendering" option in the   |
|                            | :ref:`sec:menu-r` is also turned on as otherwise|
|                            | slow rendering can result.                      |
+----------------------------+-------------------------------------------------+
| SPLASH_VZERO_CODEUNITS     | if set to a comma separated list of vector      |
|                            | components (e.g.                                |
|                            | ``export SPLASH_VZERO_CODEUNITS='0.0,1.0,0.0'`` |
|                            | ), can be used to subtract a mean velocity from |
|                            | all particles — specified in code units as read |
|                            | from the dump file.                             |
+----------------------------+-------------------------------------------------+
| SPLASH_MARGIN_XMIN         | can be used to manually adjust the left page    |
|                            | page margin (set to fraction of viewport,       |
|                            | negative values are allowed).                   |
+----------------------------+-------------------------------------------------+
| SPLASH_MARGIN_XMAX         | right horizontal page margin (set to fraction   |
|                            | of viewport).                                   |
+----------------------------+-------------------------------------------------+
| SPLASH_MARGIN_YMIN         | bottom (vertical) page margin                   |
+----------------------------+-------------------------------------------------+
| SPLASH_MARGIN_YMAX         | top (vertical) page margin                      |
+----------------------------+-------------------------------------------------+

.. _sec:asplash:

Ascii data read
~~~~~~~~~~~~~~~~

For several data reads there are environment variables which can be set
at runtime which are specific to the data read. For the ascii data read
(``asplash``) these are:

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
|                                   | ``columns`` file containing the   |
|                                   | labels for ascii data (e.g.       |
|                                   | setenv ASPLASH_COLUMNSFILE        |
|                                   | ’/home/me/mylabels’). Overridden  |
|                                   | by the presence of a local        |
|                                   | ``columns`` file.                 |
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

For the GADGET read (``gsplash``) the environment variable options are:

+-----------------------------------+-----------------------------------+
| GSPLASH_FORMAT                    | if set = 2, reads the block       |
|                                   | labelled GADGET format instead of |
|                                   | the default (non block labelled)  |
|                                   | format.                           |
+-----------------------------------+-----------------------------------+
| GSPLASH_USE_Z                     | if ``yes`` or ``true`` uses the   |
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
| GSPLASH_CHECKIDS                  | if set to ``yes`` or ``true``,    |
|                                   | reads and checks particle IDs,    |
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
| GSPLASH_IGNORE_IFLAGCOOL          | if set to ``yes`` or ``true``,does|
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

For the VINE read (``vsplash``) the environment variable options are:

+-----------------------------------+-----------------------------------+
| VSPLASH_HFAC                      | if ``yes`` or ``true`` multiplies |
|                                   | smoothing length read from the    |
|                                   | dump file by a factor of 2.8 (for |
|                                   | use with older VINE dumps where   |
|                                   | the smoothing length is defined   |
|                                   | as in a Plummer kernel rather     |
|                                   | than as the usual SPH smoothing   |
|                                   | length).                          |
+-----------------------------------+-----------------------------------+
| VSPLASH_MHD                       | if ``yes`` or ``true`` reads VINE |
|                                   | dumps containing MHD arrays (note |
|                                   | that setting VINE_MHD also        |
|                                   | works).                           |
+-----------------------------------+-----------------------------------+

sphNG data read
~~~~~~~~~~~~~~~~

For the sphNG and PHANTOM read (``ssplash``) the environment variable
options are:

+-----------------------------------+-----------------------------------+
| SSPLASH_RESET_CM                  | if ``yes`` or ``true`` resets the |
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

For the dragon read (``dsplash``) the environment variable options are:

+-----------------------------------+-----------------------------------+
| DSPLASH_EXTRACOLS                 | specifies number of extra columns |
|                                   | present in the file which are     |
|                                   | dumped after the itype array      |
+-----------------------------------+-----------------------------------+

Stephan Rosswog data read
~~~~~~~~~~~~~~~~~~~~~~~~~~

For the srosph read (``rsplash``) the environment variable options are:

+-----------------------------------+-----------------------------------+
| RSPLASH_FORMAT                    | can be ``MHD`` or ``HYDRO`` which |
|                                   | read the appropriate data format  |
|                                   | from either the MHD or            |
|                                   | hydrodynamic codes                |
+-----------------------------------+-----------------------------------+
| RSPLASH_RESET_COM                 | if ``yes`` or ``true`` resets the |
|                                   | positions such that the centre of |
|                                   | mass is exactly at the origin.    |
+-----------------------------------+-----------------------------------+
| RSPLASH_COROTATING                | if ``yes`` or ``true`` then       |
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

For the ndspmhd read (``nsplash``) the environment variable options are:

+-----------------------------------+-----------------------------------+
| NSPLASH_BARYCENTRIC               | plots barycentric quantities for  |
|                                   | one-fluid dust instead of         |
|                                   | creating fake second set of       |
|                                   | particles                         |
+-----------------------------------+-----------------------------------+

H5Part data read
~~~~~~~~~~~~~~~~~

For the H5PART read (``h5splash``) the environment variable options are:

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
