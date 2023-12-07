SPLASH - an interactive visualisation tool for SPH data
=======================================================

About
-----
SPLASH is a free and open source visualisation tool for Smoothed Particle Hydrodynamics (SPH) simulations in one, two and three dimensions, developed mainly for astrophysics. It uses a command-line menu but data can be manipulated interactively in the plotting window.

Data is read *directly* from the code dump format giving rapid access to results and the visualisation is advanced forwards and backwards through timesteps by single keystrokes.

SPLASH uses the SPH smoothing kernel to render plots of density and other physical quantities, giving a smooth representation of the data. The goal is to produce beautiful plots and visualisations from SPH codes, instead of simple particle plots.

SPLASH can also be used as a standalone plotting tool for any kind of tabulated or image data from ascii, csv or .fits files.

Status
------
![build](https://github.com/danieljprice/splash/workflows/build/badge.svg)
[![docs](https://readthedocs.org/projects/splash-viz/badge/?version=latest)](https://splash-viz.readthedocs.io/en/latest/?badge=latest)

Example
-------
<img src="https://splash-viz.readthedocs.io/en/latest/_images/default-mode.png" alt="Accretion disc visualisation with SPLASH" height="300"/>

Links
-----

- Project homepage: http://users.monash.edu.au/~dprice/splash
- Code repository: https://github.com/danieljprice/splash
- Documentation: https://splash-viz.readthedocs.io
- Code paper: https://adsabs.harvard.edu/abs/2007PASA...24..159P

Install
-------
For installation instructions see the [userguide](https://splash-viz.readthedocs.io/en/latest/getting-started.html).

Usage
------------
```
   splash mydata.txt
```

Command line mode (to screen):
```
    splash -r density dump_0*
```
Command line mode (to file):
```
    splash -r density -dev myplot.pdf dump_0*
```

See the [userguide](https://splash-viz.readthedocs.io) for more.

Contributing
------------
We welcome contributions, including (but not limited to):

1. Code, via [pull request](https://github.com/danieljprice/splash/pulls). Please read developer section of user guide for guidelines.
2. Documentation, also by [pull request](https://github.com/danieljprice/splash/pulls). Docs can be edited in the docs/ directory of the main code.
3. Suggestions for features or bug reports, via the [issue tracker](https://github.com/danieljprice/splash/issues/new). Please file bugs via github rather than by email.
4. Discussion via the [mailing lists](http://users.monash.edu.au/~dprice/splash/mailinglists.html)

Citation
--------
Please cite [Price (2007)](https://adsabs.harvard.edu/abs/2007PASA...24..159P)
when using SPLASH.

License
-------
See [LICENCE](LICENCE) file for usage and distribution conditions

Copyright (c) 2004-2023 Daniel Price and contributors
