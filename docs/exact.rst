
.. _sec:exact:

Exact solutions
================
splash contains a library of exact or analytic solutions to common
test problems that have been used extensively in benchmarking
codes (e.g. [PM04]_, [RP07]_, [Price12]_, [Price18]_).

Errors
------

The error norms calculated when exact solutions are plotted are as
follows: The error for each particle is given by

.. math:: e_i = f_i - f_{exact},

where the exact solution :math:`f_{exact}(x)` is the solution returned
from the exact solution subroutines (with resolution adjustable in the
exact solution options menu option) interpolated to the position of the
current particle :math:`x_i` via a simple linear interpolation. The
absolute :math:`L_1` error norm is simply the average of the errors
across the domain, calculated according to

.. math:: \Vert e \Vert_{L_1} = \frac{1}{N f_{max}} \sum_{i=1}^N \vert e_i \vert,

where :math:`f_{max}` is the maximum value of the exact solution in the
region in which the particles lie (also only particles in the current
plot are used) which is used to normalise the error estimate. A better
error norm is the :math:`L_2` or *Root Mean Square* (RMS) norm given by

.. math::

   \Vert e \Vert_{L_2} = \left[\frac{1}{N} \left( \frac{1}{f_{max}^2} \sum_{i=1}^N \vert e_i
   \vert^2 \right)\right]^{1/2}.

Finally the maximum error, or :math:`L_\infty` norm is calculated
according to

.. math:: \Vert e \Vert_{L_\infty} = \frac{1}{f_{max}} {\rm max}_i \vert e_i \vert.

which is the most stringent error norm.

The inset plot of the individual particle errors shows the fractional
deviation for each particle given by

.. math:: e_{i,frac} = (f_i - f_{exact}) / f_{exact}.

Shock tubes (Riemann problem)
-----------------------------

The subroutine ``exact_shock`` plots the exact solution for a
one-dimensional shock tube (Riemann problem). The difficult bit of the
problem is to determine the jump in pressure and velocity across the
shock front given the initial left and right states. This is performed
in a separate subroutine (riemannsolver) as there are many different
methods by which this can be done (see e.g.
[Toro92]_). The actual subroutine exact_shock
reconstructs the shock profile (consisting of a rarefaction fan, contact
discontinuity and shock, summarised in :numref:`fig:shocktube`), given the post-shock values of
pressure and velocity.

.. figure:: figs/sodshock.pdf
   :alt: exact solution for one-dimensional shock tube
   :name: fig:shocktube
   :width: 100.0%

   Example of exact solution for one-dimensional shock tube problem (red
   line) compared to the SPH solution (black line/particles), utilising
   the exact solutions incorporated in splash

The speed at which the shock travels into the ‘right’ fluid can be
computed from the post shock velocity using the relation

.. math:: v_{shock} = v_{post}\frac{(\rho_{post}/\rho_R)}{(\rho_{post}/\rho_R)- 1},

where the jump conditions imply

.. math:: \frac{\rho_{post}}{\rho_R} = \frac{(P_{post}/P_R) + \beta}{1 + \beta (P_{post}/P_R)}

with

.. math:: \beta = \frac{\gamma - 1}{\gamma + 1}.

Riemann solver
~~~~~~~~~~~~~~~

The algorithm for determining the post-shock velocity and pressure is
taken from [Toro92]_.

Polytrope
---------

The subroutine ``exact_polytrope`` computes the exact solution for a
static polytrope with arbitrary :math:`\gamma`. From Poisson’s equation

.. math:: \nabla^2 \phi = 4\pi G \rho,

assuming only radial dependence this is given by

.. math::
   :label: eq_poissonsph

   \frac{1}{r^{2}} \frac{d}{dr} \left(r^{2} \frac{d\phi}{dr} \right) = 4\pi G \rho(r).

The momentum equation assuming an equilibrium state
(:math:`{\bf v} = 0`) and a polytropic equation of state
:math:`P = K\rho^{\gamma}` gives

.. math::
   :label: eq_polyk

   \frac{d\phi}{dr} = - \frac{\gamma K}{\gamma-1}\frac{d}{dr} \left[\rho^{(\gamma -1)} \right]

Combining (:eq:`eq_poissonsph`) and
(:eq:`eq_polyk`) we obtain an equation for the density
profile

.. math::
   :label: eq:dens

   \frac{\gamma K}{4\pi G (\gamma - 1)} \frac{1}{r^{2}} \frac{d}{dr} \left[r^{2}
   \frac{d}{dr}\left( \rho^{\gamma-1} \right) \right] + \rho(r) = 0.

This equation can be rearranged to give

.. math::

   \frac{\gamma K}{4\pi G (\gamma - 1)} \frac{d^2}{dr^2}
   \left[r\rho^{\gamma-1}\right] + r\rho = 0.

The program solves this equation numerically by defining a variable

.. math:: \mathcal{E} = r \rho^{\gamma-1}

and finite differencing the equation according to

.. math::

   \frac{\mathcal{E}^{i+1} - \mathcal{E}^i + \mathcal{E}^{i-1}}{(\Delta r)^2} =
   \frac{4\pi G (\gamma - 1)}{\gamma K} r
   \left(\frac{\mathcal{E}}{r}\right)^{1/(\gamma-1)}.

Linear wave
-----------

The subroutine ``exact_wave`` simply plots a sine function on a given
graph. The function is of the form

.. math:: y = \sin{(k x - \omega t)}

where :math:`k` is the wavenumber and :math:`\omega` is the angular
frequency. These parameters are set via the input values of wavelength
:math:`\lambda = 2\pi/k` and wave period :math:`P = 2\pi/\omega`.

.. table:: Input parameters for the linear wave exact solution

   +-----------------+------------+
   | :math:`\lambda` | wavelength |
   +-----------------+------------+
   | :math:`P`       | period     |
   +-----------------+------------+

Sedov blast wave
----------------

The subroutine ``exact_sedov`` computes the self-similar Sedov solution
for a blast wave.

Toy stars
---------

The subroutine ``exact_toystar1D`` computes the exact solutions for the
‘Toy Stars’ described in [MP04]_. The system is one
dimensional with velocity :math:`v`, density :math:`\rho`, and pressure
:math:`P`. The acceleration equation is

.. math:: \frac{dv}{dt} = - \frac{1}{\rho} \frac{\partial P}{\partial x}  - \Omega^2 x,

We assume the equation of state is

.. math:: P = K \rho^\gamma,

The exact solutions provided assume the equations are scaled such that
:math:`\Omega^2 = 1`.

Static structure
~~~~~~~~~~~~~~~~~

The static structure is given by

.. math:: \bar \rho = 1- x^2,

Linear solutions
~~~~~~~~~~~~~~~~~

The linear solution for the velocity is given by

.. math:: v = 0.05 C_s G_n(x) \cos{\omega t} ).

The density is

.. math:: \rho = \bar{\rho} + \eta,

where

.. math:: \eta = 0.1 C_s \omega P_{n+1}(x) \sin{(\omega t)}.

Non-linear solution
~~~~~~~~~~~~~~~~~~~~

In this case the velocity is given by

.. math:: v = A(t) x,

while the density solution is

.. math:: \rho^{\gamma -1} = H(t) - C(t) x^2.

where the parameters A, H and C are determined by solving the ordinary
differential equations

.. math::

   \begin{aligned}
   \dot{H} & = & -AH(\gamma -1), \\
   \dot{A} & = & \frac{2K \gamma}{\gamma -1} C - 1 - A^2 \\
   \dot{C} & = & -AC(1+ \gamma),\end{aligned}

The relation

.. math::
   :label: eq:kconst

   A^2 = -1 - \frac{2 \sigma C}{\gamma -1} + kC^{\frac{2}{\gamma +1}},

is used to check the quality of the solution of the differential
equations by evaluating the constant :math:`k` (which should remain
close to its initial value).

MHD shock tubes
---------------

The subroutine ``exact_mhdshock`` computes the exact solution of the
Riemann problem in ideal magnetohydrodynamics for arbitrary left and
right states :math:`(\rho, p, v_x, v_y, v_z, B_y, B_z)` and normal field
:math:`B_x`, for an adiabatic or (with :math:`\gamma = 1`) isothermal
equation of state. The fast, rotational (Alfvén) and slow waves either
side of the contact discontinuity are each described by a single
strength parameter, and the parameters are iterated with a damped
Newton-Raphson method until the pressure, velocity and transverse
magnetic field match at the contact, following [Torrilhon02]_ and
[TY14]_. Shocks are found from the jump conditions, written as a cubic
(fast shocks) or quadratic (slow shocks) equation for the compression
ratio, and rarefactions by integrating along the right eigenvectors.

The physically relevant solution of the MHD Riemann problem is not
always unique. When the transverse fields and velocities lie in a plane
(e.g. the [BW88]_ problem) the solution returned is that of the coplanar
MHD system, in which the slow wave may be an intermediate shock or a
compound wave ([FK01]_; [TY13]_). This is the solution found by numerical
simulations. Otherwise the regular solution is returned. Problems with
zero transverse field on one side are handled using switch-on shocks and
rarefactions.

The left and right states are read automatically from the ``.setup``
file written by Phantom (``densleft``, ``prleft``, ``vxleft``, ...,
``Bxleft``, ``Byleft``, ``Bzleft`` and the corresponding right states)
if present. Otherwise they can be entered manually, or set from one of
the preset problems of [BW88]_, [RJ95]_, [Balsara98]_ and [DW94]_.

Special relativistic MHD shock tubes
------------------------------------

The subroutine ``exact_mhdshock_sr`` computes the exact solution of the
Riemann problem in special relativistic MHD, following [GR06]_. As in
that paper the solution is regular, consisting of fast, Alfvén and slow
waves either side of a contact discontinuity (or only fast waves and a
tangential discontinuity when :math:`B_x = 0`), so compound waves seen
in numerical solutions of e.g. the relativistic Brio & Wu problem are
not reproduced. Shocks and Alfvén discontinuities are found by solving
the jump conditions :math:`[[F]] = V_s [[U]]` together with one condition
fixing the wave strength, and rarefactions by integrating along the
eigenvectors of the flux Jacobian, with the magnetosonic speeds from the
quartic dispersion relation (e.g. [MB06]_). The solutions reproduce the
tables in [GR06]_ for the test problems of [Komissarov99]_ and
[Balsara01]_, which are provided as presets. Velocities and fields are
measured in the lab frame, with :math:`c = 1`, and the density may be
plotted either in the rest frame or as :math:`\rho^* = \rho W`.

h vs :math:`\rho`
-----------------

The subroutine exact_hrho simply plots the relation between smoothing
length and density, i.e.,

.. math:: h = h_{\rm fact} \left(\frac{m}{\rho}\right)^{1/\nu}

where :math:`\nu` is the number of spatial dimensions. The parameter
:math:`h_{\rm fact}` is output by the code into the header of each timestep.
For particles of different masses, a different curve is plotted for each
different mass value.
