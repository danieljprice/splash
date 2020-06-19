.. _sec:coordtransforms:

Coordinate transforms
=====================

Particle positions and vectors defined on the particles can be plotted
in non-cartesian coordinate systems. The coordinate system can be set
via the particle plot o)ptions menu, via the “change coordinate system”
option. The actual coordinate transformations are defined in a
standalone Fortran module called ``geometry.f90`` and the precise
details can be determined by looking in this file. For reference,
however the transformations are given below.

Cylindrical Polar Coordinates
------------------------------

For cylindrical coordinates the transformations are:

.. math::

   \begin{array}{lclp{1cm}lcl}
   R & = & \sqrt{x^2 + y^2}    & & x & = & R\cos\phi \\
   \phi & = & \tan^{-1}{(y/x)} &; & y & = & R\sin\phi \\
   z & = & z                             & & z & = & z\\
   \end{array}

where vectors transform according to:

.. math::

   \begin{array}{lclp{1cm}lcl}
   v_R      & = & v_x \frac{x}{R} + v_y \frac{y}{R}  & & v_x & = & v_R \cos\phi - v_\phi \sin\phi \\
   v_\phi & = & v_x \left(\frac{-y}{R}\right) + v_y \left(\frac{x}{R}\right) &; &v_y & = & v_R \sin\phi + v_\phi \cos\phi \\
   v_z      & = & v_z & & v_z & = & v_z. \\
   \end{array}

In the case where these vectors are velocities, the :math:`v_{\phi}`
component corresponds to :math:`v_{\phi} = r\dot{\phi}`. At the origin
we assume :math:`\phi = 0`, implying :math:`\cos\phi = 1` and therefore
:math:`\frac{x}{R} = 1` and :math:`\frac{y}{R} = 0`. This ensures the
transformations are reversible everywhere.

Spherical Polar Coordinates
----------------------------

For spherical coordinates the transformations are:

.. math::

   \begin{array}{lclp{1cm}lcl}
   r & = & \sqrt{x^2 + y^2 + z^{2}}    & & x & = & r\cos\phi\sin\theta\\
   \phi & = & \tan^{-1}{(y/x)}              &; & y & = & r\sin\phi\sin\theta \\
   \theta & = & \cos^{-1}(z/r)             & & z & = & r\cos\theta \\
   \end{array}

where vectors transform according to:

.. math::

   \begin{array}{lclp{1cm}lcl}
   v_r      & = & v_x \frac{x}{r} + v_y \frac{y}{r} + v_{z}\frac{z}{r}  & & v_x & = & v_r \cos\phi\sin\theta- v_\phi \sin\phi + v_\theta \cos\phi\cos\theta \\
   v_\phi & = & v_x \left(\frac{-y}{\sqrt{x^2 + y^{2}}}\right) + v_y \left(\frac{x}{\sqrt{x^2 + y^{2}}}\right) &; &v_y & = & v_r \sin\phi\sin\theta + v_\phi \cos\phi + v_{\theta} \sin\phi\cos\theta \\
   v_\theta & = & v_{x}\frac{xz}{r \sqrt{x^{2} + y^{2}}} + v_{y}\frac{yz}{r \sqrt{x^{2} + y^{2}}} - v_{z}\frac{(x^{2} + y^{2})}{r\sqrt{x^{2} + y^{2}}}  & & v_z & = & v_r \cos\theta - v_\theta \sin\theta. \\
   \end{array}

In the case where these vectors are velocities, the components
:math:`v_{\phi}` and :math:`v_{\theta}` correspond to
:math:`v_{\phi} = r\sin{\theta}\dot{\phi}` and
:math:`v_{\theta} = r\dot{\theta}` respectively.

Toroidal Coordinates
---------------------

Toroidal coordinates represent a local frame of reference inside a
torus. The coordinate transformations are given by

.. math::

   \begin{array}{lclp{1cm}lcl}
   r & = & \sqrt{[(x^2 + y^2)^{1/2} - R]^{2} + z^{2}}    & & x & = & (r\cos\theta + R) \cos\phi \\
   \theta & = & \tan^{-1} \left[\frac{z}{(\sqrt{x^{2} + y^{2}} - R)}\right]              &; & y & = & (r\cos\theta + R)\sin\phi \\
   \phi & = & \tan^{-1}(y/x)             & & z & = & r\sin\theta \\
   \end{array}

where :math:`R` is the radius of the torus. The use of the inverse
tangent in :math:`\theta` instead of :math:`\theta = \sin^{-1}(z/r)` is
necessary to get the quadrant correct (via the ``atan2`` function).
Vectors transform according to:

.. math::

   \begin{array}{lclp{2cm}lcl}
   v_r      & = & v_x \frac{x(r_{cyl} - R)}{r r_{cyl}} + v_y \frac{y(r_{cyl} - R)}{r r_{cyl}} + v_{z} \frac{z}{r}  & & v_x & = & v_r \cos\theta\cos\phi- v_\theta \sin\theta\cos\phi - v_\phi\sin\phi \\
   v_\theta & = & v_x \frac{-zx}{r r_{cyl}}  + v_y\frac{-zy}{r r_{cyl}}  + v_{z}\frac{(r_{cyl} - R)}{r} &; &v_y & = & v_r \cos\theta\sin\phi - v_\theta \sin\theta\sin\phi + v_\phi\cos\phi \\
   v_\phi & = & v_{x} \left(\frac{-y}{r_{cyl}}\right) + v_{y} \left(\frac{x}{r_{cyl}}\right) & & v_z & = & v_{r}\sin\theta + v_{\theta} \cos\theta \\
   \end{array}

where we have defined, for convenience,

.. math:: r_{\rm cyl} = \sqrt{x^{2} + y^{2}} = r\cos\theta + R. \nonumber

and the velocities :math:`v_\theta` and :math:`v_\phi` correspond to
:math:`r \dot{\theta}` and :math:`r_{\rm cyl} \dot{\phi}`, respectively.
The torus radius :math:`R` is a parameter in the ``geometry`` module and
is set to :math:`1` by default.

Flared Cylindrical Polar Coordinates
-------------------------------------

For flared cylindrical coordinates the transformations are:

.. math::

   \begin{array}{lclp{1cm}lcl}
   R & = & \sqrt{x^2 + y^2}    & & x & = & R\cos\phi \\
   \phi & = & \tan^{-1}{(y/x)} &; & y & = & R\sin\phi \\
   z' & = & z \left(R_{\rm ref}/{R}\right)^\beta  & & z & = & z' (R/R_{\rm ref})^\beta \\
   \end{array}

where :math:`R_{\rm ref}` is the reference radius and :math:`\beta` is
the flaring index, both parameters that can be set by the user. Vectors
transform according to:

.. math::

   \begin{array}{lclp{1cm}lcl}
   v_R      & = & v_x \frac{x}{R} + v_y \frac{y}{R}  & & v_x & = & v_R \cos\phi - v_\phi \sin\phi \\
   v_\phi & = & v_x \left(\frac{-y}{R}\right) + v_y \left(\frac{x}{R}\right) &; &v_y & = & v_R \sin\phi + v_\phi \cos\phi \\
   v_{z'}      & = & -\beta \frac{xz}{R^2} \left(\frac{R_{\rm ref}}{R}\right)^\beta v_x  -\beta \frac{yz}{R^2} \left(\frac{R_{\rm ref}}{R}\right)^\beta v_y + \left(\frac{R_{\rm ref}}{R}\right)^\beta v_z & & v_z & = & \beta \frac{z'}{R} \left(\frac{R}{R_{\rm ref}}\right)^\beta v_R + \left(\frac{R}{R_{\rm ref}}\right)^\beta  v_{z'}. \\
   \end{array}

In the case where these vectors are velocities, the :math:`v_{\phi}`
component corresponds to :math:`v_{\phi} = r\dot{\phi}`.

Logarithmic Flared Cylindrical Polar Coordinates
-------------------------------------------------

For logarithmic flared cylindrical coordinates the transformations are:

.. math::

   \begin{array}{lclp{1cm}lcl}
   d & = & \log_{10} ( \sqrt{x^2 + y^2} )   & & x & = & R\cos\phi \\
   \phi & = & \tan^{-1}{(y/x)} &; & y & = & R\sin\phi \\
   z' & = & z \left(R_{\rm ref}/{R}\right)^\beta  & & z & = & z' (R/R_{\rm ref})^\beta \\
   \end{array}

where :math:`R_{\rm ref}` is the reference radius and :math:`\beta` is
the flaring index, both parameters that can be set by the user. Vectors
transform according to:

.. math::

   \begin{array}{lclp{1cm}lcl}
   v_d      & = & v_x \frac{x}{R}f^{-1} + v_y \frac{y}{R}f^{-1}  & & v_x & = & f v_d \cos\phi - v_\phi \sin\phi \\
   v_\phi & = & v_x \left(\frac{-y}{R}\right) + v_y \left(\frac{x}{R}\right) &; &v_y & = & f v_d \sin\phi + v_\phi \cos\phi \\
   v_{z'}      & = & -\beta \frac{xz}{R^2} \left(\frac{R_{\rm ref}}{R}\right)^\beta v_x  -\beta \frac{yz}{R^2} \left(\frac{R_{\rm ref}}{R}\right)^\beta v_y + \left(\frac{R_{\rm ref}}{R}\right)^\beta v_z & & v_z & = & \beta \frac{z'}{R} \left(\frac{R}{R_{\rm ref}}\right)^\beta f v_d + \left(\frac{R}{R_{\rm ref}}\right)^\beta  v_{z'}. \\
   \end{array}

where :math:`R \equiv 10^d` and correspondingly :math:`d = \log_{10} R`
and :math:`f \equiv R \ln (10)`.

