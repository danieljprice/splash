##-------------------------------------------------------------------##
##     Makefile for compiling SPLASH and linking with PGPLOT         ##
##     Written by Daniel Price					     ##
##     University of Exeter, UK, 2006                    	     ##
##                                                                   ##
##     requires GNU make (on some systems this is 'gmake' instead    ##
##                        of 'make')                                 ##
##                                                                   ##
##     see the INSTALL file for detailed installation instructions   ##
##-------------------------------------------------------------------##


.KEEP_STATE:
KNOWN_SYSTEM=no
#
# some default settings for unix systems
#
F90C= 
F90FLAGS=
ifeq ($(findstring 64,$(MACHTYPE)),64)
   X11LIBS= -L/usr/X11R6/lib64 -lX11
else
   X11LIBS= -L/usr/X11R6/lib -lX11
endif
#
# change the line below depending on where/how you have installed PGPLOT
# (some settings of the SYSTEM variable for specific machines overwrite this)
#
PGPLOTLIBS = -L$(PGPLOT_DIR) -lpgplot -lpng
#
# add one of the lines below if PGPLOT was compiled with a different
# compiler to the one you are using. May also need -L/dir/ for the directory
# where the corresponding library is located (e.g. -L/usr/local/gfortran/lib -lgfortran)
#
# g77-compiled PGPLOT
#PGPLOTLIBS += -lg2c
#
# gfortran-compiled PGPLOT
#PGPLOTLIBS += -lgfortran
#
# g95-compiled PGPLOT
#PGPLOTLIBS += -lg95
#
# this file contains system-dependent routines like getarg, iargc etc.
#
SYSTEMFILE= system_unix.f90
#
# this can be used to static link the pgplot libraries if all else fails
#
STATICLIBS=
#
# set the parallel flag to 'yes' to compile with openMP directives
#
#PARALLEL=no
#
# the openMP flags should be set in the lines defining your system type
# (ie. leave line below blank)
OMPFLAGS=
#
# the endian flag can be used to compile the code to read BIG or LITTLE endian data
# some compilers also allow this to be done at runtime (e.g. g95, ifort) by setting
# an environment variable appropriately (e.g. G95_ENDIAN or F_UFMTENDIAN)
#
#ENDIAN=
#ENDIAN='BIG'
#ENDIAN='LITTLE'

#--------------------------------------------------------------
#  the following are general settings for particular compilers
#
#  set the environment variable 'SYSTEM' to one of those
#  listed to use the appropriate settings
#
#  e.g. in tcsh use
#  setenv SYSTEM 'g95'
#
#  in bash the equivalent is
#  export SYSTEM='g95'
#
#--------------------------------------------------------------

ifeq ($(SYSTEM),g95)
#  using the g95 compiler
   F90C= g95
   F90FLAGS= -O3 -ffast-math
   SYSTEMFILE= system_f2003.f90 # this is for Fortran 2003 compatible compilers
   DEBUGFLAG= -Wall -Wextra -Wno=165 -g -fbounds-check -ftrace=full
   ENDIANFLAGBIG= -fendian='BIG'
   ENDIANFLAGLITTLE= -fendian='LITTLE'
   PARALLEL= no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), gfortran)
#  gfortran compiler (part of gcc 4.x.x)
   F90C= gfortran
   F90FLAGS= -O3 -Wall -frecord-marker=4
   SYSTEMFILE= system_unix.f90
   DEBUGFLAG= -g -frange-check
   OMPFLAGS= -fopenmp
   ENDIANFLAGBIG= -fconvert=big-endian
   ENDIANFLAGLITTLE= -fconvert=little-endian
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),nagf95)
#  NAG f95 compiler
   F90C= f95
   F90FLAGS= -O3
   SYSTEMFILE= system_unix_NAG.f90
   PARALLEL= no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),sunf95)
#  sun f95 compiler on linux
   F90C= sunf95
   F90FLAGS= -fast -ftrap=%none
   OMPFLAGS= -openmp
   DEBUGFLAG= -g -C -w4 -errtags -erroff=COMMENT_1582,COMMENT_1744 -ftrap=%all
   SYSTEMFILE= system_f2003.f90
   ENDIANFLAGBIG= -xfilebyteorder=big16:%all
   ENDIANFLAGLITTLE= -xfilebyteorder=little16:%all
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),ifort)
#  this is for the intel fortran compiler (version 10)
   F90C= ifort
   F90FLAGS= -O3 -nbs -i_dynamic
   OMPFLAGS= -openmp
   DEBUGFLAG= -C -g
   SYSTEMFILE= system_f2003.f90
   ENDIANFLAGBIG= -convert big_endian
   ENDIANFLAGLITTLE= -convert little_endian
# or use setenv F_UFMTENDIAN=big or little at runtime
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),ifort8)
#  this is for the intel fortran compiler (version 8)
   F90C= ifort
   F90FLAGS= -O3 -Vaxlib -nbs
   OMPFLAGS= -openmp
   DEBUGFLAG= -C -g
   SYSTEMFILE= system_unix.f90
   ENDIANFLAGBIG= -convert big_endian
   ENDIANFLAGLITTLE= -convert little_endian
# or use setenv F_UFMTENDIAN=big or little at runtime
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),pgf90)
#  this is for the Portland Group Fortran 90 compiler (tested with version 7.2-5)
   F90C= pgf90
   F90FLAGS= -fast -mcmodel=medium -Mbackslash -Ktrap=none
   DEBUGFLAG= -C -g -gopt -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mcoff \
              -Mdwarf1 -Mdwarf2 -Melf -Mpgicoff -traceback
   OMPFLAGS= -mp
   SYSTEMFILE= system_unix.f90
   ENDIANFLAGBIG= -Mbyteswapio  # only works on a little-endian machine
   ENDIANFLAGLITTLE=
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),pathf95)
#  this is for the Pathscale f95 compiler
   F90C= pathf95
   F90FLAGS= -Ofast -mcmodel=medium
   DEBUGFLAG= -C -g
   OMPFLAGS= -openmp
   SYSTEMFILE= system_f2003.f90
   ENDIANFLAGBIG= -convert big_endian
   ENDIANFLAGLITTLE= -convert little_endian
   KNOWN_SYSTEM=yes
endif

#--------------------------------------------------------------
#
# the following presets are machine-specific
# (ie. relate to both the compiler and a specific
#      installation of pgplot)
#
#--------------------------------------------------------------

ifeq ($(SYSTEM),ukaff1a)
#  this is for ukaff1a
   F90C= xlf90_r ##-Wl,-t
   F90FLAGS= -O3 -qnoescape -qsuffix=f=f90 -qextname 
   OMPFLAGS= -qsmp=omp
   DEBUGFLAG= -C -g
   SYSTEMFILE= system_f2003.f90
   X11LIBS= -L/usr/X11R6/lib -lX11
   PGPLOTLIBS=   
   STATICLIBS= /home/dprice/pgplot/libpgplot.a /home/dprice/plot/libg2cmodified.a 
   KNOWN_SYSTEM=yes
#   EXT=_temp
endif

ifeq ($(SYSTEM),ukaff1b)
#  this is for the new nodes on ukaff 
   F90C= pathf95
   F90FLAGS= -Ofast -mcmodel=medium
   DEBUGFLAG= -C -g
   OMPFLAGS= -openmp
   SYSTEMFILE= system_f2003.f90
   ENDIANFLAGBIG= -convert big_endian
   ENDIANFLAGLITTLE= -convert little_endian
   X11LIBS= -L/home/dprice/lib -lX11 -lpng
   PGPLOTLIBS= -L/home/dprice/pgplot2 -lpgplot
   EXT=2   # call the binary something different
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),zen)
#  this is for the intel fortran compiler
   F90C= ifort
   F90FLAGS= -O3 -mcmodel=medium -axT -ipo -warn all #-assume nounderscore
   OMPFLAGS= -openmp
   DEBUGFLAG= -C -g
   SYSTEMFILE= system_f2003.f90
   ENDIANFLAGBIG= -convert big_endian
   ENDIANFLAGLITTLE= -convert little_endian
   X11LIBS=-L/usr/X11R6/lib64 -lX11
   PGPLOTLIBS= -L${PGPLOT_DIR} -lpgplot -lpng
# or use setenv F_UFMTENDIAN=big or little at runtime
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),mymac)
#  these are the settings for a Mac G4 running Panther 
#  using g95 with pgplot installed via fink
   F90C= g95
   F90FLAGS= -O3 -ffast-math
   DEBUGFLAG= -Wall -Wextra -Wno=165 -g -fbounds-check -ftrace=full
   PGPLOTLIBS= -L/sw/lib/pgplot -lpgplot -lg2c -L/sw/lib -lpng \
          -laquaterm -lcc_dynamic -Wl,-framework -Wl,Foundation
   SYSTEMFILE= system_f2003.f90
   PARALLEL= no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),intelmac)
#  these are the settings for an Intel Macbook running Tiger
#  using gfortran with pgplot installed via fink
   F90C= gfortran
   F90FLAGS= -O3 -ffast-math
   DEBUGFLAG= -Wall -Wextra -Wno=165 -g -frange-check
   PGPLOTLIBS= -L/sw/lib/pgplot -lpgplot -L/sw/lib -lpng \
          -laquaterm -Wl,-framework -Wl,Foundation -lSystemStubs
   SYSTEMFILE= system_f2003.f90
   PARALLEL= no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), gfortran_macosx)
#   gfortran with pgplot installed via fink
    F90C= gfortran
    F90FLAGS= -O3 -Wall
    PGPLOTLIBS= -L/sw/lib/pgplot -lpgplot -lg2c -L/sw/lib -lpng \
          -laquaterm # -lSystemStubs use this on OS/X Tiger
    SYSTEMFILE= system_unix.f90
    DEBUGFLAG= -g -frange-check
    OMPFLAGS= -fopenmp
    ENDIANFLAGBIG= -fconvert=big-endian
    ENDIANFLAGLITTLE= -fconvert=little-endian
    KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),myg95)
#  using the g95 compiler
   F90C= myg95
   F90FLAGS= -O3 -ffast-math
   X11LIBS= -L/usr/X11R6/lib64 -lX11
   PGPLOTLIBS = -L/usr/local64/pgplot -lpgplot -lpng -lg2c
   SYSTEMFILE= system_f2003.f90 # this is for Fortran 2003 compatible compilers
   ENDIANFLAGBIG= -fendian='BIG'
   ENDIANFLAGLITTLE= -fendian='LITTLE'
   PARALLEL= no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),maccluster)
#  using the g95 compiler on the iMac cluster in Exeter
   ifeq ($(MACHTYPE),i386)
   #  this stuff is for building universal binaries across intel/ppc macs
      F90C= myg95
      EXT= _i386
   else
      F90C= g95
      EXT= _ppc
   endif
   F90FLAGS= -O3 -ffast-math -Wall -Wextra -Wno=165 
   X11LIBS= -L/usr/X11R6/lib -lX11
   PGPLOTLIBS= -lSystemStubs
   STATICLIBS= /AstroUsers/djp212/pgplot/libpgplot.a
   SYSTEMFILE= system_f2003.f90
   ENDIANFLAGBIG= -fendian='BIG'
   ENDIANFLAGLITTLE= -fendian='LITTLE'
   DEBUGFLAG=-fbounds-check -ftrace=full
   PARALLEL= no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),astromac)
#  using the g95 compiler on the Astrophysics cluster in Exeter
   F90C= g95
   F90FLAGS= -O3 -ffast-math -Wall -Wextra -Wno=165
   X11LIBS= -L/usr/X11R6/lib -lX11
   PGPLOTLIBS= -L/opt/local/lib -lpng -lg2c
   STATICLIBS= /usr/local/pgplot/libpgplot.a
   SYSTEMFILE= system_f2003.f90
   ENDIANFLAGBIG= -fendian='BIG'
   ENDIANFLAGLITTLE= -fendian='LITTLE'
   DEBUGFLAG=-fbounds-check -ftrace=full
   PARALLEL= no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),spectrum)
#  sun f95 compiler on linux
   F90C= f95
   F90FLAGS= -fast -ftrap=%none
   OMPFLAGS= -openmp
   DEBUGFLAG= -g -C -w4 -errtags -erroff=COMMENT_1582,COMMENT_1744 -ftrap=%all
   SYSTEMFILE= system_f2003.f90
   ENDIANFLAGBIG= -xfilebyteorder=big16:%all
   ENDIANFLAGLITTLE= -xfilebyteorder=little16:%all
   X11LIBS = -L/usr/X11R6/lib64 -lX11
   KNOWN_SYSTEM=yes
endif
#
# these are the flags used for linking
#
LDFLAGS= $(X11LIBS) $(PGPLOTLIBS)

#
# this is an option to change the endian-ness at compile time
# (provided the appropriate flags are specified for the compiler)
#
ifeq ($(ENDIAN), BIG)
    F90FLAGS += ${ENDIANFLAGBIG}
endif

ifeq ($(ENDIAN), LITTLE)
    F90FLAGS += ${ENDIANFLAGLITTLE}
endif

# compile in parallel
ifeq ($(PARALLEL),yes)
    F90FLAGS += $(OMPFLAGS)
endif

# add debugging flags at compile time
ifeq ($(DEBUG),yes)
    F90FLAGS += $(DEBUGFLAG)
endif

# rename the executable for development work
ifeq ($(DEV),yes)
    EXT= -dev
endif

# Fortran flags same as F90
FC = $(F90C)
FFLAGS = $(F90FLAGS)

# define the implicit rule to make a .o file from a .f90/.f95 file
# (some Make versions don't know this)

%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@
%.o : %.f90
	$(F90C) $(F90FLAGS) -c $< -o $@
%.o : %.F90
	$(F90C) $(FPPFLAGS) $(F90FLAGS) -c $< -o $@
%.o : %.f95
	$(F90C) $(F90FLAGS) -c $< -o $@

# modules must be compiled in the correct order to check interfaces
# really should include all dependencies but I am lazy

SOURCESF90= globaldata.f90 asciiutils.f90 setpage.f90 transform.f90 \
         prompting.f90 geometry.f90 plotutils.f90 colourbar.f90 \
         colours.f90 colourparts.f90 shapes.f90 units.f90 \
         write_pixmap.f90 write_sphdata.f90 analysis.f90 discplot.f90 \
         exact_fromfile.f90 exact_mhdshock.f90 \
         exact_polytrope.f90 exact_rhoh.f90 \
         exact_sedov.f90 exact_shock.f90 exact_shock_sr.f90 exact_wave.f90 \
         exact_toystar1D.f90 exact_toystar2D.f90 \
         exact_densityprofiles.f90 exact_torus.f90 \
         exact_ringspread.f90 \
         exact.f90 limits.f90 \
         allocate.f90 titles.f90 \
         options_particleplots.f90 \
         calc_quantities.f90 get_data.f90 convert.f90 \
         options_data.f90 \
         options_limits.f90 options_page.f90 \
	 options_powerspec.f90 options_render.f90 \
	 options_vecplot.f90 options_xsecrotate.f90 pdfs.f90 \
         rotate.f90 interpolate1D.f90 interpolate2D.f90 \
         interpolate3D.f90 interpolate3D_xsec.f90 \
         interpolate3D_projection.F90 \
         interpolate3D_opacity.f90 interpolate_vec.f90 \
         interactive.f90 \
         fieldlines.f90 legends.f90 particleplot.f90 \
         powerspectrums.f90 render.f90 \
         plotstep.f90 timestepping.f90 \
         defaults.f90 menu.f90 \
         $(SYSTEMFILE) system_utils.f90 splash.f90

# these are `external' f77 subroutines
SOURCESF= 

OBJECTS1 = $(SOURCESF:.f=.o) $(SOURCESF90:.f90=.o) $(STATICLIBS)
OBJECTS= $(OBJECTS1:.F90=.o)

#
# Now compile with the appropriate data read file
# (move yours to the top so that you can simply type "make")
#
all: ascii gadget vine sphNG srosph dragon tipsy

ascii: checksystem $(OBJECTS) read_data_ascii.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o asplash$(EXT) $(OBJECTS) read_data_ascii.o
#--build universal binary on mac cluster
   ifeq ($(SYSTEM), maccluster)
	lipo -create asplash_ppc asplash_i386 -output asplash || cp asplash$(EXT) asplash
   endif

mbatesph: checksystem $(OBJECTS) read_data_mbate.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o bsplash $(OBJECTS) read_data_mbate.o

gadget: checksystem $(OBJECTS) read_data_gadget.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o gsplash $(OBJECTS) read_data_gadget.o

gadget_jsb: checksystem $(OBJECTS) read_data_gadget_jsb.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o gsplash_jsb $(OBJECTS) read_data_gadget_jsb.o

bauswein: checksystem $(OBJECTS) read_data_bauswein.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o bsplash $(OBJECTS) read_data_bauswein.o 

dragon: checksystem $(OBJECTS) read_data_dragon.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o dsplash $(OBJECTS) read_data_dragon.o

vine: checksystem $(OBJECTS) read_data_VINE.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o vsplash $(OBJECTS) read_data_VINE.o

ndspmhd: checksystem $(OBJECTS) read_data_dansph.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o splash $(OBJECTS) read_data_dansph.o

dansph: checksystem $(OBJECTS) read_data_dansph_old.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o dsplash $(OBJECTS) read_data_dansph_old.o

foulkes: checksystem $(OBJECTS) read_data_foulkes.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o fsplash $(OBJECTS) read_data_foulkes.o 

jjm: checksystem $(OBJECTS) read_data_jjm.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o jsplash $(OBJECTS) read_data_jjm.o 

jules: checksystem $(OBJECTS) read_data_jules.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o jsplash $(OBJECTS) read_data_jules.o 

scwsph: checksystem $(OBJECTS) read_data_scw.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o wsplash $(OBJECTS) read_data_scw.o

srosph: checksystem $(OBJECTS) read_data_sro.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o rsplash $(OBJECTS) read_data_sro.o 

spyros: checksystem $(OBJECTS) read_data_spyros.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ssplash $(OBJECTS) read_data_spyros.o

sphNG: checksystem $(OBJECTS) read_data_sphNG.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ssplash$(EXT) $(OBJECTS) read_data_sphNG.o
#--build universal binary on mac cluster
   ifeq ($(SYSTEM), maccluster)
	lipo -create ssplash_ppc ssplash_i386 -output ssplash || cp ssplash$(EXT) ssplash
   endif

RSPH: rsph

rsph: checksystem $(OBJECTS) read_data_rsph.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o rsplash $(OBJECTS) read_data_rsph.o

tipsy: checksystem $(OBJECTS) read_data_tipsy.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o tsplash $(OBJECTS) read_data_tipsy.o

vanaverbeke: checksystem $(OBJECTS) read_data_vanaverbeke.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o vsplash $(OBJECTS) read_data_vanaverbeke.o 

ucla: checksystem $(OBJECTS) read_data_UCLA.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o usplash $(OBJECTS) read_data_UCLA.o 

sky: ucla

steve: foulkes

sigfried: vanaverbeke

gasoline: tipsy

myall: ndspmhd dansph sphNG srosph gadget mbatesph tipsy ascii

checksystem:
   ifeq ($(KNOWN_SYSTEM), yes)
	@echo ""
	@echo "Compiling splash for $(SYSTEM) system..........."
	@echo ""
        ifeq ($(ENDIAN), BIG)
	     @echo "Flags set for conversion to BIG endian"
        endif
        ifeq ($(ENDIAN), LITTLE)
	     @echo "Flags set for conversion to LITTLE endian"
        endif
        ifeq ($(PARALLEL), yes)
	     @echo "Compiling the PARALLEL code"
        else
	     @echo "Compiling the SERIAL code"
        endif
   else
	@echo ""
	@echo " Makefile for splash by Daniel Price "
	@echo " -- see INSTALL file for detailed instructions"
	@echo ""
	@echo " make: ERROR: value of SYSTEM = $(SYSTEM) not recognised..."
	@echo " =>set the environment variable SYSTEM to one listed "
	@echo "   in the Makefile and try again"
	@echo ""
	quit
   endif

## other stuff

doc:
	cd docs; pdflatex splash; pdflatex splash

tar:
	tar cf splash.tar Makefile $(SOURCESF90) $(SOURCESF) read_data*.f90

targz:
	tar cf splash.tar Makefile $(SOURCESF90) $(SOURCESF) read_data*.f90
	gzip splash.tar

## unit tests of various modules as I write them

tests: test1 test2 test3

test1: interpolate3D_projection.o interpolate3D_xsec.o ./tests/test_interpolate3D.o 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o test_interpolation3D ./tests/test_interpolate3D.o interpolate3D_projection.o interpolate3D_xsec.o
test2: transform.o ./tests/test_transform.o 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o test_transform ./tests/test_transform.o transform.o
test3: fieldlines.o ./tests/test_fieldlines.o 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o test_fieldlines ./tests/test_fieldlines.o fieldlines.o

clean:
	rm *.o *.mod
