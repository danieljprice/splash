
.. _sec:writeyourown:

Data reads and command line options
=====================================

Requesting a customised data reader
------------------------------------
The ``splash`` binary will read many formats. The default format is ascii or csv any ascii or csv data file where
columns correspond to different quantities and rows correspond to each particle (actually
I use splash to plot graphs for nearly all data in this form, whether SPH or not)
-- it will also sensibly skip header lines which do not have the same number of columns.

However, it is ultimately desirable to use SPLASH to directly visualise the
(binary) output of your code. If your format is not amongst those distributed::

   splash -f

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
