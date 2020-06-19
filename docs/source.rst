
Source code overview
====================

Here is a brief and outdated description of various files making up the
code:

+-----------------------------------+-----------------------------------+
| Filename                          | Description                       |
+===================================+===================================+
|                                   |                                   |
+-----------------------------------+-----------------------------------+
| allocate.f90                      | allocates memory for main arrays  |
+-----------------------------------+-----------------------------------+
| calc_quantities.f90               | calculates additional quantities  |
|                                   | from particle data                |
+-----------------------------------+-----------------------------------+
| colours.f90                       | colour schemes for rendering      |
+-----------------------------------+-----------------------------------+
| colourparts.f90                   | colours particles                 |
+-----------------------------------+-----------------------------------+
| defaults.f90                      | writes/reads default options      |
|                                   | to/from file                      |
+-----------------------------------+-----------------------------------+
| exact.f90                         | module handling exact solution    |
|                                   | settings                          |
+-----------------------------------+-----------------------------------+
| exact_densityprofiles.f90         | various :math:`N-`\ body density  |
|                                   | profiles                          |
+-----------------------------------+-----------------------------------+
| exact_fromfile.f90                | reads an exact solution tabulated |
|                                   | in a file                         |
+-----------------------------------+-----------------------------------+
| exact_mhdshock.f90                | some tabulated solutions for mhd  |
|                                   | shocks                            |
+-----------------------------------+-----------------------------------+
| exact_polytrope.f90               | exact solution for a polytrope    |
+-----------------------------------+-----------------------------------+
| exact_rhoh.f90                    | exact relation between density    |
|                                   | and smoothing length              |
+-----------------------------------+-----------------------------------+
| exact_sedov.f90                   | exact solution for sedov blast    |
|                                   | wave                              |
+-----------------------------------+-----------------------------------+
| exact_shock.f90                   | exact solution for hydrodynamic   |
|                                   | shocks                            |
+-----------------------------------+-----------------------------------+
| exact_wave.f90                    | exact solution for a propagating  |
|                                   | sine wave                         |
+-----------------------------------+-----------------------------------+
| exact_toystar.f90                 | exact solution for the toy star   |
|                                   | problem                           |
+-----------------------------------+-----------------------------------+
| exact_toystar2D.f90               | exact solution for the 2D toy     |
|                                   | star problem                      |
+-----------------------------------+-----------------------------------+
| get_data.f90                      | wrapper for main data read        |
+-----------------------------------+-----------------------------------+
| geometry.f90                      | module handling different         |
|                                   | coordinate systems                |
+-----------------------------------+-----------------------------------+
| globaldata.f90                    | various modules containing        |
|                                   | "global" variables                |
+-----------------------------------+-----------------------------------+
| interactive.f90                   | drives interactive mode           |
+-----------------------------------+-----------------------------------+
| interpolate1D.f90                 | interpolation of 1D SPH data to   |
|                                   | grid using kernel                 |
+-----------------------------------+-----------------------------------+
| interpolate2D.f90                 | interpolation of 2D SPH data to   |
|                                   | grid                              |
+-----------------------------------+-----------------------------------+
| interpolate3D_xsec.f90            | 3D cross section interpolations   |
+-----------------------------------+-----------------------------------+
| interpolate3D_projection.f90      | 3D interpolation integrated       |
|                                   | through domain                    |
+-----------------------------------+-----------------------------------+
| legends.f90                       | plots (time) legend on plot       |
+-----------------------------------+-----------------------------------+
| limits.f90                        | sets initial plot limits and      |
|                                   | writes to/reads from limits file  |
+-----------------------------------+-----------------------------------+
| menu.f90                          | main menu                         |
+-----------------------------------+-----------------------------------+
| options_data.f90                  | sets options relating to current  |
|                                   | data                              |
+-----------------------------------+-----------------------------------+
| options_limits.f90                | sets options relating to plot     |
|                                   | limits                            |
+-----------------------------------+-----------------------------------+
| options_page.f90                  | sets options relating to page     |
|                                   | setup                             |
+-----------------------------------+-----------------------------------+
| options_particleplots.f90         | sets options relating to particle |
|                                   | plots                             |
+-----------------------------------+-----------------------------------+
| options_powerspec.f90             | sets options for power spectrum   |
|                                   | plotting                          |
+-----------------------------------+-----------------------------------+
| options_render.f90                | sets options for render plots     |
+-----------------------------------+-----------------------------------+
| options_vector.f90                | sets options for vector plots     |
+-----------------------------------+-----------------------------------+
| options_xsecrotate.f90            | sets options for cross sections   |
|                                   | and rotation                      |
+-----------------------------------+-----------------------------------+
| particleplot.f90                  | subroutines for particle plotting |
+-----------------------------------+-----------------------------------+
| plotstep.f90                      | main “backbone” of the code which |
|                                   | drives plotting of a single       |
|                                   | timestep                          |
+-----------------------------------+-----------------------------------+
| read_data_dansph.f90              | reads data from my format of data |
|                                   | files                             |
+-----------------------------------+-----------------------------------+
| read_data_mbate.f90               | reads data from Matthew Bate’s    |
|                                   | format of data files              |
+-----------------------------------+-----------------------------------+
| read_data_xxx.f90                 | reads data from …                 |
+-----------------------------------+-----------------------------------+
| render.f90                        | takes array of pixels and plots   |
|                                   | render map/contours etc           |
+-----------------------------------+-----------------------------------+
| rotate.f90                        | subroutines controlling rotation  |
|                                   | of particles                      |
+-----------------------------------+-----------------------------------+
| setpage.f90                       | sets up the PGPLOT page (replaces |
|                                   | call to PGENV/PGLAB)              |
+-----------------------------------+-----------------------------------+
| splash.f90                        | main program, handles startup/    |
|                                   | command line reading              |
+-----------------------------------+-----------------------------------+
| timestepping.f90                  | controls stepping through         |
|                                   | timesteps                         |
+-----------------------------------+-----------------------------------+
| titles.f90                        | reads a list of titles to be used |
|                                   | to label each timestep            |
+-----------------------------------+-----------------------------------+
| transform.f90                     | applies various transformations   |
|                                   | to data (log10, 1/x, etc)         |
+-----------------------------------+-----------------------------------+
