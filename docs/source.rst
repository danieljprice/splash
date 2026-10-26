
Source code overview
====================

This is a brief guide to the main Fortran source files. The authoritative
list is ``src/*.f90`` in the repository; many specialised ``read_data_*``
routines exist for individual simulation codes.

+-----------------------------------+-----------------------------------+
| Filename                          | Description                       |
+===================================+===================================+
| splash.f90                        | main program; startup, version,   |
|                                   | command-line parsing              |
+-----------------------------------+-----------------------------------+
| read_data.F90                     | dispatches to format-specific     |
|                                   | read_data routines                |
+-----------------------------------+-----------------------------------+
| read_data_*.f90                   | one routine per supported data    |
|                                   | format (see ``splash --formats``) |
+-----------------------------------+-----------------------------------+
| get_data.f90                      | wrapper for main data read        |
+-----------------------------------+-----------------------------------+
| plotstep.f90                      | drives plotting of a single       |
|                                   | timestep                          |
+-----------------------------------+-----------------------------------+
| menu.f90                          | main text menu                    |
+-----------------------------------+-----------------------------------+
| interactive.f90                   | interactive mode (mouse/keyboard) |
+-----------------------------------+-----------------------------------+
| interactive_buttons.f90           | on-screen buttons in interactive  |
|                                   | plot windows (v4.0+)              |
+-----------------------------------+-----------------------------------+
| interactive_utils.f90             | coordinate transforms for         |
|                                   | interactive mouse input           |
+-----------------------------------+-----------------------------------+
| plotlib_giza.f90                  | giza plotting backend interface   |
+-----------------------------------+-----------------------------------+
| render.f90                        | render maps, contours, opacity    |
+-----------------------------------+-----------------------------------+
| interpolate*.f90                  | SPH interpolation (1D/2D/3D,      |
|                                   | projections, cross sections)      |
+-----------------------------------+-----------------------------------+
| exact*.f90, exact.f90             | exact solution test problems      |
+-----------------------------------+-----------------------------------+
| globaldata.f90                    | modules with global variables     |
+-----------------------------------+-----------------------------------+
| allocate.f90                      | memory allocation for main arrays |
+-----------------------------------+-----------------------------------+
| defaults.f90, limits.f90          | read/write ``splash.defaults``    |
|                                   | and ``splash.limits``             |
+-----------------------------------+-----------------------------------+
| options_*.f90                     | menu subsystems (render, vector,  |
|                                   | limits, page, etc.)               |
+-----------------------------------+-----------------------------------+
| write_sphdata.f90                 | ``splash to ...`` format          |
|                                   | conversion utilities                |
+-----------------------------------+-----------------------------------+
| system_utils.f90                  | strings, environment, UTF-8 text  |
+-----------------------------------+-----------------------------------+
| calc_quantities.f90               | derived quantities from particle  |
|                                   | data                              |
+-----------------------------------+-----------------------------------+
| geometry.f90                      | coordinate systems and transforms |
+-----------------------------------+-----------------------------------+
| timestepping.f90                  | stepping through dump sequences   |
+-----------------------------------+-----------------------------------+
| tests/                            | Python CLI regression tests       |
|                                   | (see ``tests/README.md``)         |
+-----------------------------------+-----------------------------------+
