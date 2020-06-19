
Introduction
============

While many wonderful commercial software packages exist for visualising
scientific data (such as the widely used Interactive Data Language), I
found they were cumbersome for particle-based data. Much of what I
wanted to do was specific to SPH, like interpolating to an array of
pixels using the kernel. While generic routines exist for such tasks, I
could not explain how they worked, and they were slow. Also, while
interactive gizmos are handy, it was more difficult to perform the same
tasks non-interactively, as required for the production of animations.
The major work in the visualisation of SPH data is not the image
production itself but the manipulation of data prior to plotting. Much
of this manipulation makes sense within an SPH framework.

splash is designed for this specific task - to use SPH tools to analyse
SPH data. Publishable images and animations can be obtained efficiently
from the raw data with a minimum amount of user effort. The development
of powerful visualisation tools has enabled me to pick up on effects
present in my simulation results that I would not otherwise have noticed
— the difference between a raw particle plot and a rendered image can be
substantial. A key goal of splash is to eliminate the use of
crap-looking particle plots as a means of representing SPH data.

What it does
------------

splash is a utility for visualisation of output from (astrophysical)
simulations using the Smoothed Particle Hydrodynamics (SPH) method in
one, two and three dimensions. It is written in Fortran 90/95/20xx and
utilises giza , a custom-build backend graphics library to do the actual
plotting. In particular the following features are included:

-  Rendering of particle data to an array of pixels using the SPH kernel

-  Cross-sections through 2D and 3D data (as both particle plots and
   rendered images).

-  Fast projections through 3D data (i.e., column density plots, or
   integration of other quantities along the line of sight)

-  Surface renderings of 3D data.

-  Vector plots of the velocity (and other vector quantities), including
   vector plots in a cross section slice in 3D.

-  Rotation and animation sequence generation for 3D data.

-  Automatic stepping through timesteps, making animations simple to
   produce.

-  Interactive mode for detailed examination of timestep data (e.g.
   zooming, rotating, stepping forwards/backwards, log axes, adapting
   limits).

-  Remote visualisation via simple X-Windows forwarding

-  Multiple plots on page, including option to automatically tile plots
   if :math:`y-` and :math:`x-` limits are the same.

-  Plot limits can be fixed, adaptive or particle tracking.

-  Exact solutions for common SPH test problems (e.g. shock tubes, sedov
   blast wave).

-  Calculation of quantities not dumped (e.g. pressure, entropy)

-  Conversion of binary dump files to ascii format.

-  Interpolation of SPH data to 2D and 3D grids.

-  Transformation to different coordinate systems (for both coordinates
   and vector components) and rescaling of data into physical units.

-  Straightforward production of both bitmap (png) and vector (eps, pdf)
   images which can then be converted into animations or inserted into
   LaTeX documents.

What it doesn’t do
------------------

splash is geared towards gas dynamics simulations with SPH and has
basically grown out of my visualisation needs. Thus it is not
particularly useful for things like water and solids in SPH. An SPH
visualisation tool geared towards the non-gaseous side of things you may
want to have a look at is pv-meshless, by John Biddiscombe:

https://twiki.cscs.ch/twiki/bin/view/ParaViewMeshless

splash also doesn’t make coffee.

splash, the paper
------------------

The algorithms implemented in splash are not described here, but instead
described in a paper [Price07]_, available from:

http://www.publish.csiro.au/?paper=AS07022

This paper should be cited if you use splash for scientific purposes,
and please do so as it is my only form of thanks.

Version History
---------------

+-----------------------+-----------------------+-----------------------+
| 2.9.0 & 05/04/19 &    |                       |                       |
| general header        |                       |                       |
| quantities are read   |                       |                       |
| and available in      |                       |                       |
| function parser; more |                       |                       |
| robust label          |                       |                       |
| detection and parsing |                       |                       |
| during ascii data     |                       |                       |
| read; splash to grid  |                       |                       |
| works in              |                       |                       |
| non-cartesian         |                       |                       |
| geometries; added     |                       |                       |
| flared and log-flared |                       |                       |
| coordinate systems;   |                       |                       |
| Doppler shift colour  |                       |                       |
| bar; can customise    |                       |                       |
| line style and colour |                       |                       |
| when plotting         |                       |                       |
| multiple exact        |                       |                       |
| solutions; seg faults |                       |                       |
| fixed; better plot    |                       |                       |
| tiling decisions;     |                       |                       |
| disappearing arrows   |                       |                       |
| bug fix; Rafikov      |                       |                       |
| disc-planet exact     |                       |                       |
| solution added; atan2 |                       |                       |
| implemented in        |                       |                       |
| function parser;      |                       |                       |
| various multigrain    |                       |                       |
| phantom read fixes    |                       |                       |
| (incl. seg faults);   |                       |                       |
| exact rendering       |                       |                       |
| implemented in 2D;    |                       |                       |
| libsplash implemented |                       |                       |
| for use as Python     |                       |                       |
| splash backend        |                       |                       |
+-----------------------+-----------------------+-----------------------+

Licence
-------

splash - a visualisation tool for SPH data ©2004-2020 Daniel Price. This
program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version. This program is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details. You should have received a
copy of the GNU General Public License along with this program; if not,
write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
Boston, MA 02111-1307 USA.

