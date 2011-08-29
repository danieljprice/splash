To compile SPLASH with the PGPLOT backend (the default in SPLASH v1.x.x, but
still an option with SPLASH 2.x) you will need the following on your system,
both of which are freely available:

- The PGPLOT graphics subroutine library
- A Fortran 95 compiler

The basic steps for installation are:
1) make sure you have a Fortran 90/95 compiler (such as g95 or gfortran).
2) make sure you have the PGPLOT libraries installed.
3) compile SPLASH and link with PGPLOT.
4) if desired/necessary write a read_data subroutine so that SPLASH can read
   your data format.
5) make pretty pictures.

For troubleshooting of some common installation problems,
 have a look at the online FAQ.

1) ---------------- Fortran 95 compilers ---------------------------

 By now, many Fortran 90/95 compilers exist. In terms of free ones, both Intel
and Sun have non-commercial versions available for Linux and the g95 compiler,
downloadable from:

http://www.g95.org

successfully compiles SPLASH and if necessary the PGPLOT libraries.

Gfortran is also free and, as of version 4.2.0, works. The latest version
can be downloaded from:

http://gcc.gnu.org/wiki/GFortran

I strongly recommend downloading a more recent version of gfortran rather than
relying on any pre-installed version (use gfortran -v to check the version number).
In particular versions 4.1.0 and lower *do not* compile splash. Later versions also
have openMP, so you can compile and run SPLASH in parallel.

2) ------------------- PGPLOT -----------------------------------------
The PGPLOT graphics subroutine library is freely downloadable from

http://www.astro.caltech.edu/~tjp/pgplot/

or by ftp from

ftp://ftp.astro.caltech.edu/pub/pgplot/pgplot5.2.tar.gz

however check to see if it is already installed on your system (if so, the
libraries are usually located in /usr/local/pgplot).

If PGPLOT is already installed, make sure that the environment variable PGPLOT_DIR
is set to the location of the PGPLOT installation directory (e.g. /usr/local/pgplot).
Check this by typing "echo $PGPLOT_DIR".

If instead you are following the steps below, set the PGPLOT_DIR environment
variable to the directory to which you will install PGPLOT (e.g. export PGPLOT_DIR=$HOME/pgplot).

It is a good idea to add the setting of PGPLOT_DIR into your .profile/.bashrc or .tcshrc file
along with a setting for PGPLOT_DEV (e.g. to "/xw" which sets the default device 
to be the X-windows device).

 --- installing PGPLOT yourself ---

Whilst detailed installation instructions are given in the PGPLOT distribution,
the general procedure for installing your own version (if necessary) is given
below (otherwise skip to part 3).

Note: to compile PGPLOT with the png and X-windows drivers may require some packages
to be installed from your linux distribution. In Ubuntu these are "libpng-dev" and "libX11-dev"
which can be installed with "sudo apt-get libpng-dev" and "sudo apt-get libX11-dev". Otherwise
you will encounter errors regarding missing header files -- e.g. "cannot find png.h" and a whole bunch
of errors.

a) untar the pgplot5.2.tar.gz file (e.g. in your home space):
   "tar xvfz pgplot5.2.tar.gz"
b) rename the directory pgplot to something else
   (e.g. "mv pgplot pgplotsrc").
c) make a directory called pgplot and enter it:
   "mkdir pgplot; cd pgplot"
d) copy the drivers.list file from the pgplotsrc directory
   "cp ../pgplotsrc/drivers.list ."
e) edit the drivers.list file and uncomment the following drivers:
   /NULL /XW /VPS /CPS /VCPS /PS and either /PNG /TPNG or /GIF /VGIF.
   Optionally also /PPM /VPPM.
f) make a makefile using the makemake utility:
   "../pgplotsrc/makemake ../pgplotsrc linux g77_gcc"
   (I generally always use "linux" and "g77_gcc" regardless of the actual operating
   system). You should now have a file called "makefile" in your directory.
g) edit the makefile, replacing the compiler "g77" with "gfortran" or "g95" (or
   any other compiler) as appropriate. For some compilers, you may need to change
   the FFLAGC= to just -O2.
h) type "make": this should compile through all the .f files and may or may not
   compile the GIF drivers (if not, go back to step d and remove the /GIF drivers
   from the list), but will definitely die for the png driver with an error like:
   make: *** No rule to make target `png.h', needed by `pndriv.o'.  Stop.
i) edit the makefile and initially try just commenting out the line:
   pndriv.o : ./png.h ./pngconf.h ./zlib.h ./zconf.h
   by adding a preceding hash as follows:
   #pndriv.o : ./png.h ./pngconf.h ./zlib.h ./zconf.h
   If this does not work, try giving the correct paths for these files (use "locate png.h" to find where
   they are) -- usually this will be the following:
   pndriv.o : /usr/include/png.h /usr/include/pngconf.h /usr/include/zlib.h /usr/include/zconf.h
   
   Also add the above include path to the CFLAGC= variable, by amending the line
   "CFLAGC= xxx" to "CFLAGC= xxx -I/usr/include/"
   (e.g. on a mac you may need -I/sw/include/ or -I/opt/include/).
   
   Also, for the time being, comment out the line SHARED_LIB= (stuff) and
   replace it with a blank setting:
   SHARED_LIB= 
   so that only the static library (libpgplot.a) is built. If you succeed with
   this you can try going back and uncommenting this line and typing "make"
   again to build the dynamic library.

   You may further need to add the path for the png library (libpng.a) to the
   LIBS= flag (e.g. LIBS=stuff -L/sw/lib/).
j) now type "make" again, and you should obtain a successful build. Check that
   your installation is OK by running the demo programs pgdemo1, pgdemo2 etc.

3) ------------------- Compiling the code  --------------------------------

Once you have a copy of PGPLOT installed, you can proceed with the SPLASH
installation. Prior to doing so you should make sure that the PGPLOT_DIR environment
variable is set (check using "echo $PGPLOT_DIR"). If not, add the following lines
to your ~/.bashrc file (or the tcsh equivalent):

export PGPLOT_DIR=/me/whereeveriputit/pgplot
export PGPLOT_DEV=/xw         (this sets the default PGPLOT device).
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PGPLOT_DIR

Now, having untarred the SPLASH bundle: "tar xvfz splash-x.x.x.tar.gz",
you should have a directory called splash/. Enter this directory: "cd splash"
and have a look in the Makefile.

Preset options for the most common Fortran compilers are given in the Makefile
provided the variable SYSTEM is set appropriately (either on the command line or as
an environment variable). On the command line, type

"make SYSTEM=xxx BACKEND=pgplot"

Where the SYSTEM corresponds to one of those listed in the Makefile, some of the most
commonly used of which are:

gfortran -- settings for the gfortran compiler
g95 -- settings for the g95 compiler
nagf95 -- settings for the NAG f95 compiler
sunf95 -- settings for the Sun f95 compiler
ifort -- settings for the Intel Fortran Compiler
pgf90 -- settings for the Portland Group Fortran 90 compiler

Options which compile and also link PGPLOT on specific machines are:

mymac -- settings for a Mac using g95 with PGPLOT installed via fink
ukaff1a -- settings for the ukaff1a supercomputer

If you have the PGPLOT_DIR environment variable set then linking with PGPLOT and
associated libraries (libpng, libX11) *might just work* and you should find a whole
bunch of splash binaries (asplash, gsplash, ssplash etc.) have been created. If so,
then you are done with compilation and can skip directly to step 4 or 5. If not, read
on.

  ---- porting to a new SYSTEM ------

 If none of the SYSTEM variables corresponds to your local Fortran compiler, it
should be reasonably straightforward to add your own. For example you will need
to set the Fortran compiler and flags to your local version, e.g..

F90C = g95
F90FLAGS = -O

and importantly, on some compilers you will need to make sure that backslashes
(\) are treated as normal characters. For example on the following compilers you
should use:

intel fortran compiler:
F90C = ifc/ifort
F90FLAGS = -O -nbs

portland group fortran:
F90C = pgf90
F90FLAGS = -O -Mbackslash

Secondly, you will need to modify the system-dependent routines for your
compiler. These are specified via the settings:

SYSTEMFILE = system_f2003.f90 

which uses Fortran 2003 standard calls (supported by almost all recent
compilers). Alternatively the file system_unix.f90 should also work for older
(and newer) unix-based compilers. A file system_unix_NAG is included for the NAG
f95 compiler.

 The whole idea of SPLASH is that filenames should be read off the command line,
though sometimes there can be library clashes (e.g. two libraries defining the
same function)  which make these calls not work. In this case there are some
slightly convoluted workarounds given in the online FAQ.

 If you have ported to a new compiler, please send me an email with your new
SYSTEM variable and I will add it to the SPLASH Makefile, both for you and for
others (then you can just update directly). Remember it is always likely that
someone else in the same department may download SPLASH one day...

  ------------------- linking with PGPLOT libraries ---------------------

 Secondly the compiler must be able to link to the PGPLOT and X11 libraries on
your system. The settings for these are given at the top of the Makefile by the
settings:

X11LIBS= -lX11
PGPLOTLIBS= -lpgplot

If that works at a first attempt, take a moment to think several happy thoughts
about your system administrator. If these libraries are not found, you will need
to enter the library paths by hand. On most systems this is something like:

X11LIBS= -L/usr/X11R6/lib -lX11
PGPLOTLIBS= -L/usr/local/pgplot -lpgplot

(assuming the PGPLOT libraries are in the /usr/local/pgplot directory and the
X11 libraries are in /usr/X11R6/lib). If this does not work, try using the
"locate" command to find the libraries on your system:

user> locate libpgplot

user> locate libX11

 If, having found the PGPLOT and X11 libraries, the program still won't compile,
it is usually because the PGPLOT on your system has been compiled with a
different compiler to the one you are using, and the libraries from that
compiler must also be linked. For g77-compiled PGPLOT (very common) the relevant
library is g2c, so use:

PGPLOTLIBS = -L/usr/local/pgplot -lpgplot -lg2c

similarly for gfortran-compiled PGPLOT the appropriate library is libgfortran,
so use -lgfortran, and for g95-compiled PGPLOT, libg95, so use -lg95 (where
again you may need the -L flag to specify the location of the libxxx.a or
libxxx.so file).

If the PNG drivers are incorporated into the PGPLOT installation, the -lpng
libraries must also be added.

***
 A good, failsafe way of working out exactly what flags are required to link to
PGPLOT on your system is to look in the PGPLOT makefile itself, at exactly which
flags were used to compile the pgdemo programs (pgdemo1, pgdemo2 -- you should
also run these to check that PGPLOT has been installed correctly). For example
in the PGPLOT makefile on my Mac (located in /sw/lib/pgplot), the flags are 

LIBS=-L/usr/X11R6/lib -lX11 -L/sw/lib -laquaterm -Wl,-framework -Wl,Foundation

so these are the flags needed to link to PGPLOT, PLUS fink had used g77 to
compile pgplot, so I also need to add the -lg2c flag (see above).

Obviously a way round having to work out which compiler libraries to add is to
simply make sure that PGPLOT has been compiled with the same compiler you are
using to compile SPLASH.

 It is worth remembering also, that if you work in an Astronomy department, it
is almost certain that there will be someone in your department who uses PGPLOT
on your system and knows how to link to it, so it is worth asking around.
***

 Another last-resort option for linking to PGPLOT is to compile with the static
libraries explictly on the command line. To do this simply set the variable
STATICLIBS, e.g.

STATICLIBS=/usr/local/libpgplot.a

then no -lpgplot is needed and the library is treated like a normal .o file at
compile time.

Have a look at the online FAQ for some tips on common problems with linking to
PGPLOT (e.g. font problems on 64-bit machines).

4) -------------- reading your data format -------------------

The basic "splash" binary is quite general and will read any data where columns 
correspond to different quantities and rows correspond to each particle (actually
I use splash to plot graphs for nearly all data in this form, whether SPH or not)
-- it will also sensibly skip header lines which do not have the same number of columns.

However, it is ultimately desirable to use SPLASH to directly visualise the
(binary) output of your code. If you are using a widely used SPH code (e.g. GADGET,
GASOLINE, VINE, DRAGON), it is reasonably  likely that I have already written a 
read data subroutine which will read your dumps. If your format is not amongst those 
distributed, then BEFORE you start writing your own routine, please consider whether 
or not a routine to read your format would be of more general use (e.g. to other users
of your code). If so, PLEASE email me to request a new read_data routine for your 
format, by sending an email attaching:
a) an example dump 
and
b) the source code from the routine which wrote the dump file.

Then I can write a read for your format that can be added to the SPLASH repository
and distributed in all future versions. Whilst I aim never to change the interface
to the read_data routines, it is not impossible that some changes may occur
somewhere down the line (or enhanced functionality -- for example the more advanced
data reads are able to read only the required columns for a given plot from the
file, rather than the whole file).

If you *really* want to hack one yourself it is best to look at some of the
other examples and change the  necessary parts to suit your data files. Note
that reading directly from unformatted data files is *much* faster than reading
from formatted (ascii) output. Just to get started you can use the
read_data_ascii.f90 which reads from ascii files, but this will not enable the
full rendering capabilities until you specify the location of the density, h and
particle mass in the arrays (via the parameters ih, irho and ipmass in the
set_labels subroutine which is part of the read_data file). 

If you do end up writing your own, again, please email me the end result so I
can add it to the officially supported data reads. This also makes it much
easier for you to upgrade to newer versions as you do not require a locally
customised version.

5) ----- running splash/ making pretty pictures -----

For detailed help on how to use SPLASH, refer to the (quite extensive) userguide
in the /docs directory or on the splash web page.

Have fun! And remember, if you get stuck you can always email me... 
(it doesn't hurt).

Daniel Price
daniel.price@monash.edu
