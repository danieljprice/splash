Frequently Asked Questions
==========================
Here is a collection of the most common problems with using SPLASH. You may also wish to browse the archives of the `splash-users <http://groups.google.com/group/splash-users>`_ discussion list. See also the old `pgplot faq <http://users.monash.edu.au/~dprice/splash/pgplot.html>`__.

.. _sec:moviemaking:

How do I make a movie from splash output?
-----------------------------------------

First, make sure ffmpeg is installed using your package manager (e.g. ``brew install ffmpeg`` or ``sudo apt install ffmpeg``)

Since v3.10.0 splash now has a direct-to-mp4 backend. Hence to make a movie just give a filename like ``movie.mp4``
or select ``/mp4`` at the device prompt which will produce a file called ``splash.mp4``

You can also try a completely hands-free approach from the command line using the ``--movie`` flag::

   splash --movie disc_0*

That's it!

Using ffmpeg to make a movie from png files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With ffmpeg installed, you can use the script provided in splash/scripts/movie.sh, e.g.::

   ~/splash/scripts/movie.sh

This should produce a movie called ``movie.mp4`` from your splash_*.png files. The contents of this script are just some options for ffmpeg that produce a nice movie::

   #!/bin/bash
   #
   # script to make an mpeg4 movie from splash_*.png files
   # requires ffmpeg utility
   #
   # DJP, Feb 2014
   #
   opts='-r 10 -vb 50M -bt 100M -vf setpts=4.*PTS'
   ffmpeg -i splash_%04d.png $opts movie.mp4


The main options to play with are the bitrate (``-bt 100M``) that controls the quality and the ``-vf setpts=4.*PTS`` which controls how fast the movie plays (``setpts=1.*PTS`` gives 10 frames per second, this is usually too fast so we slow it down).

The industry-standard compression codec is the ``H.264`` MPEG codec which gives a nice small movie that is still good quality. The best way to reduce the movie size further is to restrict the data rate.

.. important::
   The .mp4 files produced by ffmpeg as above will not play directly in the Slack app. To ensure this you need to mandate that the pixel format is 
   "yuv420p", by adding the flag ``-pix_fmt yuv420p``. You must also ensure that the page dimensions are given in even (not odd) numbers of pixels.

Can I make animated gifs?
~~~~~~~~~~~~~~~~~~~~~~~~~

Another option is to produce an animated gif using the
``convert`` utility that comes with the `ImageMagick <http://www.imagemagick.org>`_ package. As an example (kindly provided by John Mansour), use::

  convert -delay 10 -loop 0 splash_*.png animatedgif.gif

where delay gives a delay in ms between frames, and loop determines whether the animated gif will loop (0= infinite, 1=1
loop, etc). The problem with an animated GIF is that they are slow unless you have *lots* of memory. However it is easy to convert an animated GIF into other movie formats.

I'm getting a segmentation fault, is splash broken?
----------------------------------------------------

This could be a stacksize issue (particularly with large data reads). Try::

   ulimit -s unlimited

before invoking splash.

How do I write a postscript/eps file for LaTeX?
-----------------------------------------------

See :ref:`sec:postscript`. At the graphics device prompt, type::

   Graphics device/type (? to see list, default /xw): /eps

choose ``/eps``. This will write the output to a file called splash.eps in the current directory.
Alternatively you can specify the filename explicitly by typing "myfile.eps". The
file can be incorporated directly into LaTeX documents.

How can I get a colour scheme like the one Joe Bloggs uses?
-----------------------------------------------------------

Send me an image (any sensible format will do) containing a colour bar using the scheme you
want to rip off and I will add it into splash. Also feel free to play with `my scripts for grabbing colour maps <https://github.com/danieljprice/extractcmap>`_ yourself.


What about boundaries? How does the rendering work near a boundary?
-------------------------------------------------------------------

Usual practice in SPH simulations near boundaries is to introduce ghost
particles which mirror the real particles. splash does not explicitly
setup any ghost particles but will use any that are present in the data
(see next question for how to specify multiple particle types).
Additional particle types contribute to the rendering calculations but
not to the determination of the plot limits. Note, however, that splash
does *not* set up ghost particles itself, as this may depend on the type
and location of the boundary. Thus if your simulation uses ghost
particle boundaries, the ghost particles should be dumped alongside the
gas particles in the output file so that their positions, masses,
densities and smoothing lengths can be read into splash and used to
render the image appropriately.

How does splash handle multiple particle types?
-----------------------------------------------

splash can handle up to 6 different particle types. These can be turned
on and off in the :ref:`sec:menu-o`.
These types are be specified in the set_labels part of the read_data
routine, which contains some lines of code along the lines of:

::

   ntypes = 3
   labeltype(1) = 'gas'
   labeltype(2) = 'ghost'
   labeltype(3) = 'sink'
   UseTypeInRenderings(1) = .true.
   UseTypeInRenderings(2) = .true.
   UseTypeInRenderings(3) = .false.

which says that there are 3 particle types, with names as given, and
that types 1 and 2 are SPH particles and should be used in the rendering
where appropriate (i.e., only when plotting of this type is turned on in
the :ref:`sec:menu-o`). Particle types which are to be used in renderings
should have masses, densities and smoothing lengths read. Non-SPH
particle types (e.g. sink particles) can be optionally plotted on top of
:ref:`sec:renderplot`.

What does SPLASH stand for?
----------------------------
Urrmmm... it has SPH in it and it sounded good. I thought of:

- "Some Pretty Little Application for Smoothed (particle) Hydrodynamics"
- "Smoothed Particles Look Amazingly Stunning Here"
- "So People Love Analysing Simulations of Hydrodynamics"
- "Simulating Particles Like A Superfast Horse"

Your suggestions on a postcard please.

SPLASH is so great. Can I send you loads of money?
--------------------------------------------------
I accept donations in the form of citations to the
`SPLASH paper <https://ui.adsabs.harvard.edu/abs/2007PASA...24..159P/abstract>`_ (Price, 2007, PASA, 24, 159-173). Just like sending cash, only... not.
This may change if I am flooded with requests from people wanting to send large
sums of money.
