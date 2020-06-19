
.. _sec:writeyourown:

Custom data reads
=====================================

Writing your own data read is not recommended. The best way is just to email me a
sample data file and a copy of the routine that wrote it. I am very
happy to do this, will mean that your read is officially supported, will
appear in the development repository, and will be updated with new
features as necessary. It doesn’t matter if your code only has one user,
I am still happy to do this as it makes splash more widely useable and
saves trouble later.

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

