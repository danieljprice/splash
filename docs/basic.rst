.. _sec:basic:

Basic splash usage
==================

Simple two column plot
----------------------

Once you have successfully compiled splash with a read data file that
will read your data format, splash is invoked with the name of the data
file(s) on the command line, e.g.

::

   splash myrun*.dat

| where splash should be replaced with ‘asplash’, ‘gsplash’ etc.
  depending on the data format.
| After a successful data read, the menu should appear as something like
  the following (the example given is for a “minidump” from Stephan
  Rosswog’s SPH code):

::

   dprice$ rsplash minidump.00001

::

       _                                                 _
      (_)   _               _           _         _     (_)_
         _ (_)    ___ _ __ | | __ _ ___| |__     (_)   _  (_)
      _ (_)  _   / __| '_ \| |/ _` / __| '_ \       _ (_)
     (_)  _ (_)  \__ \ |_) | | (_| \__ \ | | |  _  (_) _
         (_)  _  |___/ .__/|_|\__,_|___/_| |_| (_)  _ (_)
             (_)  (_)|_| (_) (_)  (_)(_) (_)(_) (_)(_)

     ( B | y ) ( D | a | n | i | e | l ) ( P | r | i | c | e )

   ...etc...

::

    You may choose from a delectable sample of plots
   -------------------------------------------------------
     1) x                     7) particle mass
     2) y                     8) B\dx
     3) z                     9) B\dy
     4) h                    10) B\dz
     5) \gr                  11) div B
     6) T
   -------------------------------------------------------
    12) multiplot [  4 ]      m) set multiplot
   -------------------------------------------------------
    d(ata) p(age) o(pts) l(imits) le(g)end h(elp)
    r(ender) v(ector) x(sec/rotate) s,S(ave) q(uit)
   -------------------------------------------------------
   Please enter your selection now (y axis or option):

The simplest plot is of two quantities which are not both coordinates.
For example, to plot density vs smoothing length, type

::

   Please enter your selection now (y axis or option): 5
   (x axis) (default=1): 4
    Graphics device/type (? to see list, default /xwin): /xw

The ``default=`` refers to the default value assigned if you just press
the return key. The last prompt asks for the device to which output
should be directed. A full list of available graphics devices is given
by typing ‘?’ at the prompt. Some of the most useful devices are given
in :ref:`tab:devices`. In the above we have selected
the X-window driver which means that the output is sent to the screen
(provided X-windows is running), as demonstrated in the screenshot shown
in :numref:`fig:rhoh`.

Many useful tasks can now be achieved by moving the mouse to the plot
window and selecting areas or pressing keystrokes – this is “interactive
mode”. Pressing ‘h’ in the plot window shows (in the terminal) the full
list of commands. Of the more useful ones are: pressing ‘l’ with the
mouse over the colour bar to use a logarithmic axis, press ’a’ on either
the colour bar or inside the plot to adapt the plot limits, select an
area with the mouse to zoom. See also :ref:`sec:interactive`.

To exit the plot, move the mouse to the plot window and press ’q’
(quit). To exit splash altogether press ’q’ again from the splash main
menu (in the terminal).

.. figure:: figs/rhoh.jpg
   :alt: Screenshot of simple two column plot to an X-window
   :name: fig:rhoh
   :width: 80.0%

   Screenshot of simple two column plot to an X-window

.. table:: Commonly used graphics devices available in giza
   :name: tab:devices

   +-----------------+-----------------+-----------------+-----------------+
   | ``/xw``,        | X-Window        | ``/png``        | Portable        |
   | ``/xwin``       | (interactive)   |                 | Network         |
   |                 |                 |                 | Graphics        |
   |                 |                 |                 | (bitmap)        |
   +-----------------+-----------------+-----------------+-----------------+
   | ``/eps``        | Encapsulated    | ``/svg``        | Scalable Vector |
   |                 | postscript (one |                 | Graphics        |
   |                 | file per page)  |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+
   | ``/pdf``        | PDF             | ``/null``       | null device (no |
   |                 |                 |                 | output)         |
   +-----------------+-----------------+-----------------+-----------------+
   | ``/ps``         | Postscript (all |                 |                 |
   |                 | pages in one    |                 |                 |
   |                 | file)           |                 |                 |
   +-----------------+-----------------+-----------------+-----------------+

.. _sec:renderplot:

Rendered plots
--------------

A more complicated plot is where both the :math:`x-` and :math:`y-` axes
refer to coordinates. For example

::

   Please enter your selection now (y axis or option):2
   (x axis) (default=1): 1
   (render) (0=none) ([0:11], default=0):5
   (vector plot) (0=none, 8=B) ([0:8], default=0):0
   Graphics device/type (? to see list, default /xwin): /xw

Notice that in this case that options appeared for rendered and vector
plots. Our choice of “5” at the (render) prompt corresponds to column 5,
which in this case is the density, producing the plot shown in the
screenshot in :numref:`fig:renderplot`.

.. figure:: figs/renderplot.jpg
   :alt: Screenshot of 3D column density plot to an X-window
   :name: fig:renderplot
   :width: 80.0%

   Screenshot of 3D column density plot to an X-window

Note that the render prompts only appear if, in the read_data
subroutine, values are set for the integer parameters irho, ipmass and
ih corresponding to the locations of density, particle mass and
smoothing length in the data arrays and provided the number of
coordinate dimensions is 2 or greater (splash can be used for SPH codes
in 1, 2 and 3 dimensions and even for plotting ascii data where there
are no “coordinates”).

Cross section slice
-------------------

To plot a cross section slice instead of a projection in 3D, type ’x’ at
the main menu to open the ’cross section/3D plotting options’ menu and
choose option 1 “switch between cross section and projection”. Then
re-plot the rendered plot again (exactly as in the previous example
:ref:`sec:renderplot`), setting the slice position at the prompt:

::

   enter z position for cross-section slice: ([-8.328:8.327], default=0.000):

which produces the plot shown in the screenshot in :numref:`fig:renderplot_xsec`

.. figure:: figs/renderplot_xsec.jpg
   :alt: Screenshot of 3D cross section slice plot to an X-window
   :name: fig:renderplot_xsec
   :width: 80.0%

   Screenshot of 3D cross section slice plot to an X-window

Vector plots
------------

A prompt to plot vector arrows on top of rendered plots (or on top of
particle plots) appears whenever vectors are present in the data (for
details of how to specify this in your data read, see
:ref:`sec:writeyourown`), taking the form:

::

   (vector plot) (0=none, 8=B) ([0:8], default=0):0

where the number refers to the column of the first component of the
vector quantity.

Vector plots in 3D show either the integral of each component along the
line of sight or, for cross sections, the vector arrows in a cross
section slice (depending on whether a projection or cross section has
been selected for 3D plots – see the rendering examples given
previously). In 2D vector plots simply show the vector arrows mapped to
a pixel array using the SPH kernel.

Settings related to vector plots can be changed via the v)ector plot
submenu (:ref:`sec:vectorplots`). The size of the arrows is set by
the maximum plot limit over all of the vector components. Alternatively
the arrow size can be changed interactively using ’v’, ’V’ (to decrease
and increase the arrow size respectively) and ’w’ (to automatically
adjust the arrow size so that the longest arrow is of order one pixel
width).

Contour plots
-------------

To plot contours of a quantity instead of a rendered plot, simply set
the colour scheme used for rendering to 0 (contours only) via the
“change colour scheme” option in the r)ender menu (type “r2” from the
main menu as a shortcut to option 2 in the render menu).

Contours of an additional quantity can also be plotted on top of a
render plot. However the prompt for an additional contour plot does not
appear by default – it can be turned on via the “plot contours” option
in the r)ender menu (type “r3” at the main menu as a shortcut). With
this option set *and a non-zero response to the render prompt*, a prompt
appears below the render prompt:

::

   (render) (0=none) ([0:11], default=0):5
   (contours) (0=none) ([0:11], default=0):6

Entering the column to use in the contour plot at this prompt (e.g.
column 6 in the above example would correspond to the temperature) gives
a rendered plot with overlaid contours.

Entering the same quantity used in the rendering at this prompt (e.g.
column 5 in the above example) triggers a subsequent prompt for the
contour limits which can then be set differently to those used in the
render plot. In this way it is possible to make a plot where the density
of one particle type is shown by the rendered plot and the density of
another particle type (with different limits) is shown by contours. This
can be achieved because once contour plotting is turned on, the
contribution of a given particle type to either the contours or rendered
plots can be turned on or off via the “turn on/off particles by type”
option in the particle plot o)ptions menu.

Moving forwards and backwards through data files
------------------------------------------------

If you have put more than one file on the command line (or alternatively
the file contains more than one dump), it is then possible to move
forwards and backwards through the data by pressing the space bar with
the cursor in the plot window (this is “interactive mode”). To see the
keystrokes for moving backwards or moving forwards/backwards by a
specified number of steps, press ’h’ in interactive mode. If you plot to
a non-interactive device, splash simply cycles through all the files on
the command line automatically.

Zooming in and out / changing plot limits
-----------------------------------------

Having plotted to an interactive device (e.g. /xw), tasks such as
zooming in and out, selecting, colouring and hiding particles, changing
the limits of both the plot and the colour bar and many other things can
be achieved using either the mouse (i.e., selecting an area on which to
zoom in) or by a combination of the mouse and a keystroke (e.g. move the
mouse over a particle and press ’c’ to see the size of the smoothing
circle for that particle). One of the most useful commands in
interactive mode is ’a’ (adapt plot limits) which can be used to restore
the plot limits to the maximum values for the data currently plotted
(similarly pressing ’a’ on the colour bar resets the colour bar limits
to the minimum and maximum values of the rendered quantity). Pressing
’h’ in interactive mode (that is, with your mouse in the plotting
window) gives the full list of interactive commands (note that the text
appears in the terminal from which splash was invoked). Press ’s’ in the
plot window to save changes between timesteps, otherwise the settings
will revert when you move to the next timestep.

These tasks can also be achieved non-interactively by a series of
drop-down submenus invoked from the main menu by typing a single
character. For example limits changing options are contained in the
l)imits submenu, so to manually set plot limits we would type “l” from
the main menu, then “2” for option 2 (set manual limits) and follow the
prompts to set the limits for a particular data column.

.. _sec:postscript:

Producing an encapsulated postscript figure for a paper
-------------------------------------------------------

Producing a postscript plot suitable for inclusion in a LaTeX file is
simple: at the device prompt, type

::

    Graphics device/type (? to see list, default /xw): /eps

that is, instead of “/xw” (for an X-window), simply type “/eps” or
“.eps” to use the encapsulated postscript driver. This produces a file
which by default is called ``splash.eps``, or if multiple files have
been read, a sequence of files called ``splash_0000.eps``,
``splash_0001.eps``, etc. To specify both the device and filename, type
the full filename (e.g. ``myfile.eps``) as the device. Files produced in
this way can be directly incorporated into LaTeX using standard packages
such as graphicx, psfig or epsfig.

Note that postscript devices do not have a ‘background’ colour, so plots
with a ‘black’ background and ‘white’ foreground will have invisible
axes labels when viewed in (e.g.) gv (actually, they are there in white
but the background is transparent - try inserting the figure into
Keynote or Powerpoint with a dark background). For plots in papers you
will therefore need to use a ‘black’ or similarly dark foreground colour
(set via the p)age submenu). When setting the foreground and background
colours an option appears such that annotation drawn over the rendered
region can be drawn in the opposite colour - thus enabling black axes
labels (off the plot) but white text in the legend (over the rendered
area).

.. _sec:movies:

Producing a sequence of plots for a movie
-----------------------------------------

To make a movie of your simulation, first specify all of the files you
want to use on the command line:

::

   > splash dump_*

and use an interactive device to adjust options until it looks right
(hint: for the nicest movies, best thing is to delete nearly all of the
annotation, e.g. using the backspace key in interactive mode). If in
interactive mode type ’s’ to save the current settings, then plot the
same thing again but to a non-interactive device. For example, to
generate a sequence of png files:

::

    Graphics device/type (? to see list, default /xw): /png

This will generate a series of images named ``splash_0000.png``,
``splash_0001.png``, ``splash_0002.png`` corresponding to each new
plotting page generated (or enter “``myfile.png``” at the device prompt
to generate ``myfile_0000.png``, ``myfile_0001.png``,
``myfile_0002.png``\ …).

Having obtained a sequence of images there are a variety of ways to make
these into an animation using both free and commercial software.
Suggestions on software packages to use for Mac, Linux and Windows can
be found in the online faq
(http://users.monash.edu.au/~dprice/splash/faqs.html). I generally use
the application “graphic converter” on Mac OS/X which makes quicktime
movies from a sequence of images.

Ten quick hints for producing good-looking plots
------------------------------------------------

In this section I have listed ten quick suggestions for simple changes
to settings which can improve the look of a visualisation substantially
compared to the default options. These are as follows:

#. Log the colour bar. To do this simply move the cursor over the colour
   bar and hit “l” (for log). Or non-interactively via the “apply log or
   inverse transformations to columns” option in the l)imits menu.

#. Adjust the colour bar limits. Position the mouse over the colour bar
   and left-click. To revert to the widest max/min possible for the data
   plotted, press ‘a’ with the cursor positioned over the colour bar.
   Limits can also be set manually in the l)imits submenu.

#. Try changing the colour scheme. Press ‘m’ or ‘M’ in interactive mode
   to cycle forwards or backwards through the available colour schemes.

#. Change the paper size. To produce high-resolution images/movies, use
   the “change paper size” option in the p)age menu to set the paper
   size in pixels.

#. Try using normalised interpolations. If your simulation does *not*
   involve free surfaces (or alternatively if the free surfaces are not
   visible in the figure), turning the “normalise interpolations” option
   on (in the r)ender submenu) may improve the smoothness of the
   rendering. This is turned off by default because it leads to
   funny-looking edges.

#. Remove annotation/axes. For movies, often axes are unnecessary and
   detract from the visual appeal. Axes, the colour bar and the various
   legends can be turned off in interactive mode by positioning the
   cursor appropriately and pressing backspace. Alternatively each can
   be turned off manually – axes via the “axes options” option in the
   p)age submenu; the colour bar by the “colour bar options” entry in
   the r)ender menu and the legends via options in the leg)end menu.

#. Change axes/page colours. The background colour (colour of the page)
   and foreground colour (used for axes etc) can be changed vie the “set
   foreground/background colours” option in the p)age submenu.

#. Move the legend or turn it off. The time legend can be moved by
   positioning the mouse and pressing ‘G’ in interactive mode. The
   legend can be turned off in the le(g)end submenu or by pressing
   backspace in interactive mode. Similarly the vector plot legend can
   be turned on/off in the v)ector submenu and moved by positioning the
   cursor and pressing ‘H’.

#. Use physical units on the axes. These can be set via the d)ata
   submenu. See :ref:`sec:changingunits` for more details.

#. Save settings to disk! Don’t waste your effort without being able to
   reproduce the plot you have been working on. Pressing ‘s’ in
   interactive mode only saves the current settings for subsequent
   timesteps. Pressing ‘s’ from the main menu saves these settings to
   disk. Pressing ‘S’ from the main menu saves both the plot options
   *and* the plot limits, so that the current plot can be reproduced
   exactly when splash is next invoked. Adding an “a”, as in “SA”, “SA”
   or “sa” to the save options gives a prompt for a different prefix to
   the filenames (e.g. ``splash.defaults`` becomes ``myplot.defaults``),
   which splash can be invoked to use via the ``-p`` command line option
   (e.g. ``splash -p myplot file1 file2...``).
