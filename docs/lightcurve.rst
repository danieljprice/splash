
Lightcurve calculations
=======================


Since version 3.11.0, splash has implemented a **simple and experimental** command line feature to get a quick and rough idea of lightcurves,
assuming blackbody radiation from each particles at their respective temperature.


.. _sec:lcrunhow:

How to run the code
-------------------

The basic syntax is

::

    splash calc lightcurve dump_00000 dump_00001 dump_????? --kappa --temperature

For simulations with sink (star) particles that you want to include in the lightcurve calculation,
you will usually also need to specify their effective radius and temperature (see :ref:`sec:lcsinks` below), for example::

    splash calc lightcurve dump_* --kappa --temperature --reff=6.96e10 --teff=5772

This command reads some of the settings stored in the splash.defaults, splash.limits & splash.units files.
Therefore, before running ``splash calc lightcurve``, it is necessary to first run splash with the binary files to set and save the following settings:

1. set the default limit of x, y, z (in the (l)imits submenu l2) to determine the pixel grid size (see below);

2. set the unit of length and mass to cgs units (in the (d)ata submenu d2), as the lightcurve calculation code assumes cgs units by default.

3. if your simulation contains sink particles, enable sink rendering in the (x) menu
   (option 4, “Include sinks in opacity rendering”; see :ref:`sec:lcsinks`);

4. save the settings (enter s for the (s)ave option in the main menu)

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


When ``--temperature`` is set, temperature :math:`T` is computed at **read time** from density :math:`\rho` and
specific internal energy :math:`u` using *get_temp_from_u()* in ``lightcurve_utils.f90``, assuming gas and
radiation pressure in local thermodynamic equilibrium (see :ref:`sec:lckappa`).

Opacity :math:`\kappa` used in the ray trace is described in :ref:`sec:lckappa`.


.. _sec:lckappa:

Opacity (:math:`\kappa`)
--------------------------

The lightcurve calculation needs a mass opacity :math:`\kappa` (cm\ :sup:`2`\ /g) for every particle
(see :math:`\tau_i` in :ref:`sec:lcworkhow`). Splash obtains :math:`\kappa` in one of three ways,
evaluated in ``lightcurve.f90`` in the following order of precedence:

1. **Opacity column in the data** — if the dump (or an extra column created at read time) contains
   ``kappa`` / ``opacity``, those values are used directly.

2. **Fixed opacity** — if there is no opacity column but ``--kappa=``\ *VALUE* was given on the
   command line *and* read-time opacity generation was not triggered, splash uses a constant
   :math:`\kappa =` *VALUE* for all particles (via the ``taupartdepth`` setting).

3. **Default** — if neither of the above applies, a warning is printed and :math:`\kappa = 0.35`
   cm\ :sup:`2`\ /g is used for all particles.

If an opacity column is present, it can be scaled uniformly by setting ``taupartdepth`` in
``splash.defaults`` or via ``--kappa=``\ *factor* (the stored opacities are multiplied by this
factor when it differs from unity).

**Note:** for Phantom/sphNG, a bare ``--kappa`` or ``--kappa=``\ *number* on the command line
also sets ``SPLASH_GET_KAPPA`` and normally **creates** a per-particle opacity column at read time,
which then takes precedence over a constant ``taupartdepth``. To use a single constant :math:`\kappa`
in the lightcurve with sphNG data, omit ``--kappa`` / ``--kappatot`` so that no opacity column is
generated (the default 0.35 cm\ :sup:`2`\ /g applies), or supply opacities in the dump file.


Computing :math:`\kappa` at read time (``get_opacity``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For **Phantom / sphNG** dumps, splash can **create** an opacity column when the data are read by
calling *get_opacity()* in ``lightcurve_utils.f90``. This is triggered by command-line flags
(or the equivalent environment variables):

+---------------------------+-------------------------+-----------------------------------------------+
| ``--kappa``               | ``SPLASH_GET_KAPPA``    | compute :math:`\kappa` per particle at read   |
|                           |                         | time (electron scattering only; see below)    |
+---------------------------+-------------------------+-----------------------------------------------+
| ``--kappatot``            | ``SPLASH_GET_KAPPATOT`` | same, but use the full opacity recipe from    |
|                           |                         | Matsumoto & Metzger (2022, MM22)              |
+---------------------------+-------------------------+-----------------------------------------------+

**Preferred usage:** ``--kappa`` or ``--kappatot`` on the command line (same naming convention as
other ``SPLASH_*`` options; see :ref:`sec:commandline`).

These flags are only implemented in the sphNG data read. Other formats must supply ``kappa`` in the
file or use a fixed ``--kappa=``\ *VALUE*.

**Requirements when computing opacities at read time:**

* ``--temperature`` (``SPLASH_GET_TEMP`` / ``--temp``) is usually required first, so that a
  temperature column exists; *get_opacity* needs :math:`\rho` and :math:`T` in cgs units.
* The read must include density and specific internal energy (or an existing temperature column).
* Composition enters through hydrogen and helium mass fractions
  ``--xfrac`` / ``--yfrac`` (``SPLASH_XFRAC``, ``SPLASH_YFRAC``), defaulting to solar
  (:math:`X = 0.698`, :math:`Y = 0.287`, :math:`Z = 1 - X - Y`).

**Electron-scattering only (``--kappa``):**

* Uses *ionisation_fraction_Honly* — a hydrogen-only, analytic Saha solve for the free electron
  number density :math:`n_e`.
* Opacity is Thomson electron scattering only:

.. math::

   \kappa_{\rm es} = \frac{\sigma_{\rm T} \, n_e}{\rho}

* Absorption terms (Kramers, H\ :sup:`-`, molecules) are **not** included in this mode.

**Full opacity (``--kappatot``):**

Uses the MM22 combination (valid between :math:`T \sim 1.5\times 10^3` and :math:`10^9` K):

.. math::

   \kappa = \kappa_{\rm mol} + \kappa_{\rm es} + \frac{1}{\dfrac{1}{\kappa_{\rm H}} + \dfrac{1}{\kappa_{\rm K}}}

with

* :math:`\kappa_{\rm K} = 1.2\times 10^{26}\, Z (1+X)\, \rho\, T^{-3.5}` (free–free / Kramers),
* :math:`\kappa_{\rm H} = 1.1\times 10^{-25}\, \sqrt{Z}\, \sqrt{\rho}\, T^{7.7}` (H\ :sup:`-`),
* :math:`\kappa_{\rm mol} = 0.1\, Z`,
* :math:`\kappa_{\rm es}` as above (still using the hydrogen-only Saha routine for :math:`n_e`).


Computing temperature at read time (``get_temp_from_u``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With ``--temperature`` (``SPLASH_GET_TEMP``), sphNG/Phantom reads add a temperature column by solving

.. math::

   \rho u = \frac{3}{2}\frac{k_{\rm B}}{\mu m_{\rm H}}\rho T + a T^4

for :math:`T` (Newton–Raphson), with mean molecular weight :math:`\mu = 0.6`. This is the temperature
passed to *get_opacity* when ``--kappa`` or ``--kappatot`` is also set.


Typical command-line combinations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Gas-only, solar composition, opacities from :math:`\rho`, :math:`u` (sphNG):

  ::

      splash calc lightcurve dump_* --temperature --kappa

* Include absorption opacities (MM22):

  ::

      splash calc lightcurve dump_* --temperature --kappatot

* Fixed grey opacity (formats **without** sphNG read-time ``--kappa`` generation, or when no
  per-particle column is created):

  ::

      splash calc lightcurve dump_* --kappa=0.35

* Opacity already in dump — no ``--kappa`` / ``--kappatot`` flags; optional ``--kappa=2.0`` only
  scales the existing column if needed.


.. _sec:lcsinks:

Sink particles
--------------

Many SPH codes represent stars or compact objects as **sink particles** (particle types whose
label contains ``sink``, ``compact object`` or ``point mass``). These can be included in the
lightcurve ray trace alongside the gas.

Enabling sinks in the ray trace
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sink particles are only traced if the **rendersinks** option is enabled. This is off by default
(``RENDERSINKS=F`` in ``splash.defaults``). To turn it on, run splash interactively on a
representative dump file, open the :ref:`sec:menu-x` submenu, choose option 4
(3D surface rendering), enable 3D opacity rendering if prompted, and answer **yes** to
“Include sinks in opacity rendering”. Save your settings before running ``splash calc lightcurve``.

Without this setting, sink particles are skipped in the opacity interpolation (they receive zero
weight) and will not contribute to the synthetic lightcurve.

How sinks are treated in the calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When sinks are included, they are handled differently from SPH gas particles
(see ``lightcurve.f90`` and ``interpolate3D_opacity.f90``):

* Each sink is treated as a **grey blackbody emitter** at effective temperature :math:`T_{\rm eff}`
  (read from the temperature column, or set with ``--teff``; see below).
* The **effective radius** :math:`R_{\rm eff}` is taken from the sink particle's smoothing-length
  column (or set with ``--reff``). Internally this is converted to an SPH smoothing length
  :math:`h = R_{\rm eff}/W(0)` for the column-density kernel used in the ray trace.
* The bolometric luminosity of each sink is printed as
  :math:`L = 4\pi R_{\rm eff}^2 \sigma T_{\rm eff}^4`.
* Opacity along the sink sightline uses the same :math:`\kappa` as other particles; if :math:`\kappa`
  is missing or non-positive for a sink, a default of 0.35 cm\ :sup:`2`\ /g is used.

Gas particles continue to use kernel-interpolated :math:`\rho`, :math:`T` and :math:`\kappa` as
described in :ref:`sec:lcworkhow`.

Setting :math:`R_{\rm eff}` and :math:`T_{\rm eff}` with ``--reff`` and ``--teff``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sink particles in the dump file often do not carry meaningful stellar radii or effective
temperatures. Splash can overwrite these values for each sink **when the data are read**
(in ``adjust_data.f90``, before unit rescaling), using comma-separated lists:

+----------------------+-----------------------+-----------------------------------------------------+
| ``--reff=R1,R2,...`` | ``SPLASH_REFF``       | effective radius of sink 1, sink 2, … in code units |
+----------------------+-----------------------+-----------------------------------------------------+
| ``--teff=T1,T2,...`` | ``SPLASH_TEFF``       | effective temperature of sink 1, sink 2, … in K     |
+----------------------+-----------------------+-----------------------------------------------------+

**Preferred usage:** pass values on the command line with ``--reff`` and ``--teff``, as in::

    splash calc lightcurve dump_* --kappa --temperature --reff=6.96e10,3.5e10 --teff=5772,4000

For a binary, the first value applies to the first sink particle, the second to the second, and so on.
Only entries with values :math:`>`\ 0 are applied; omitted or zero entries leave the dump-file value
unchanged.

The same lists can be given as environment variables (``export SPLASH_REFF=...``,
``export SPLASH_TEFF=...``) for backwards compatibility. Internally, splash resolves both forms
via the same mechanism used for other ``SPLASH_*`` settings: the command-line flag name is the part
after the last underscore in lower case (see :ref:`sec:commandline`). **New scripts and workflows
should use ``--reff`` and ``--teff`` rather than environment variables.**

Some data reads (e.g. Phantom/sphNG) can map a ``Teff`` field from the dump into the temperature
column automatically; ``--teff`` still overrides per-sink values at read time when given.

Example: single star in cgs units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After configuring cgs length units in splash and saving ``splash.defaults`` with
``RENDERSINKS=T``::

    splash calc lightcurve star_00000 --kappa --temperature --reff=6.96e10 --teff=5772

where ``6.96e10`` cm is approximately one solar radius and ``5772`` K is the solar effective
temperature.


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

5. ``area`` (printed on screen, not in the .out file) is the 2D area of optically thick pixels in which :math:`\tau \ge 1`.

6. ``temp``   (``[04     T_{eff}]``) is the effective temperature determined from lum & area.

7. ``rphoto`` (``[03     R_{eff}]``) is the effective photospheric radius determined by lum & temp.

8. ``badfrac`` (``[08 bad pix %]``) is the fraction of the optically thick photosphere that is **under-resolved** by the SPH particle distribution (see below).


Under-resolved pixels (``badpix`` / ``badfrac``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

During the lightcurve ray trace, splash builds a per-pixel map ``badpix`` in
``interpolate3D_opacity.f90``. A pixel is flagged (``badpix = 1``) when, along the
backward ray from the observer, a **single particle** contributes optical depth
:math:`\tau \ge 1/3` to that pixel while the **total** optical depth at the pixel is still
:math:`\tau < 1`. In other words, the transition to optical depth unity is being set by
individual smoothing spheres that are large compared with the pixel size, rather than by a
smooth, many-particle photosphere.

Splash sums the area of flagged pixels (``badarea``) and reports:

* on screen: ``unresolved area`` in au\ :sup:`2` and the percentage ``badfrac``;
* in ``lightcurve.out``: column 8 as ``badfrac``, defined as

.. math::

   {\rm badfrac} = \frac{\rm badarea}{\rm area}

where ``area`` is the area of pixels with :math:`\tau \ge 1` used for the effective
temperature estimate. If ``area`` is zero, ``badfrac`` is set to zero.

A warning is printed when ``badfrac`` exceeds 5%, indicating that a significant fraction of the
photosphere may not be reliably resolved and that ``lum``, ``temp`` and ``rphoto`` should be
treated with caution. Increasing the number of particles is the only good way to
reduce ``badfrac``.



