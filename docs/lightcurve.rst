
Lightcurve calculations
=======================


Since version 3.1.0, splash has implemented a **very simplistic and experimental** command line feature to get a quick and rough idea of lightcurves,
assuming blackbody radiation from each particles at their respective temperature.


.. _sec:lcrunhow:

How to run the code
-------------------

The basic syntax is

::

    splash calc lightcurve dump_00000 dump_00001 dump_????? --kappa --temperature
    

This command reads some of the settings stored in the splash.defaults, splash.limits & splash.units files.
Therefore, before running ``splash calc lightcurve``, it is necessary to first run splash with the binary files to set and save the following settings:

1. set the default limit of x, y, z (in the (l)imits submenu l2) to determine the pixel grid size (see below);

2. set the unit of length and mass to cgs units (in the (d)ata submenu d2), as the lightcurve calculation code assumes cgs units by default.

3. save the settings (enter s for the (s)ave option in the main menu)

After setting and saving the settings, you can now exit splash and run the above commandline, with the optional ``--kappa`` and ``--temperature`` flag to ask splash to compute opacity and temperature respectively (more on that later).


By default, this calculates the lightcurve seen by an observer at +z direction infinitely away from the source.
The viewing direction can be rotated by using the ``--anglex, --angley, --anglez`` :ref:`sec:commandline`.



.. _sec:lcworkhow:

How the code works
-------------------


(As implemented in source file ``interpolate3D_opacity.f90``.)

Splash first divides the image into a two dimensional grid of 'pixels' (by default 1024x1024 pixels), and then do ray-tracing with one ray per pixel.
Each SPH particles are assumed to produce a blackbody spectrum across frequencies.
The contributions from each particles are counted forwardly by default (i.e., from the furthest particles to the observer to the closest ones), using the radiative transfer equation:

.. math::

    I_{\nu, i+1} = I_{\nu, i} e^{- \tau_i} + S_{\nu, i} ( 1 - e^{-\tau_i} )

where the source function of the i-th particle :math:`S_{\nu, i}` is assumed to be blackbody,
:math:`I_{\nu, i}` is the specific intensity at just before the i-th particle,
:math:`\tau_i` is the optical depth of the i-th particle, calculated using kernel interpolation with

.. math::
    \tau_i = \int \kappa \rho dz
    = \sum_j \kappa_j m_j w_{\mathrm{col}, ij} / h_j^2

where the :math:`w_{\mathrm{col}, ij}` is the dimensionless column kernel.


The spectra were then summed across the frequencies and pixels to produce the bolometric luminosity, listed in [02  Luminosity] in the outputted .out file.


When ``--temperature`` flag is used, Temperature :math:`T` is computed from density :math:`\rho` and the specific internal energy :math:`u`,
assuming :math:`u` is comprised of gas and radiation pressure only, with them in local thermodynamic equilibrium.

When ``--kappa`` flag is used, splash calculate the kappa with the *get_opacity()* function in the *lightcurve_utils.f90* source file.
By default (i.e. without modifying the source code), it only considers Thompson electron scattering as the *only* source for opacity,
and it assumes that the gas is made of *hydrogen only*.
The latter assumption is needed for using a simplified version of Saha equation to determine free electron number density.
There are many commented code describing the opacity contribution from Kramer's law, negative hydrogen, molecules, and absorption;
However, these contents are commented and not in effect.
Opacity can also be fixed by using ``--kappa=$VALUE`` (``$VALUE`` being a float point number in cm^2/g);
Alternatively, you can neglect the --kappa option, in which case kappa is fixed at default value 0.35 cm2/g.




.. _sec:lcoutput:

Output file explanations
------------------------

Outputs are written to the file ``lightcurve.out`` by default. Note that there are some more information printed on screen as well.

Splash roughly estimates if the photosphere in the dump has enough resolution.
This information is printed on screen.

**Beware**\: when this printed surface depth is much less than the smoothing length,
it means that there might not be enough resolution at the photosphere to determine its light-related properties,
which means that the luminosity outputted could be rubbish.


1. ``lum``    (``[02  Luminosity]``) is the actual bolometric luminosity of the star, as summed from the pixels;

2. ``Tc``     (``[07         T_c]``) is the color temperature, calculated from fitting the blackbody peak (from the star's spectrum from the pixels);

3. ``lum_bb`` (``[05     L_{bol}]``) is *not* the bolometric luminosity of the star, but rather the bolometric luminosity of a blackbody with :math:`T=T_c`, scaled such that its blackbody spectrum peak matches the actual spectrum calculated from the dump file;

4. ``r_bb``   (``[06      R_{bb}]``) is the effective photospheric radius in cm implied by the fitted blackbody spectrum;

5. ``area`` (printed on screen, not in the .out file) is the 2D area of optically thick pixels in which :math:`tau>1`.

6. ``temp``   (``[04     T_{eff}]``) is the effective temperature determined from lum & area.

7. ``rphoto`` (``[03     R_{eff}]``) is the effective photospheric radius determined by lum & temp.




