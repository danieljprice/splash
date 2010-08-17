##-------------------------------------------------------------------##
##     Makefile for compiling SPLASH and linking with                ##
##      required libraries                                           ##
##				                                     ##
##     Written by Daniel Price                                       ##
##     University of Exeter, UK,     2004-2008                       ##
##     Monash University, Australia, 2008-                           ##
##                                                                   ##
##     requires GNU make (on some systems this is 'gmake' instead    ##
##                        of 'make')                                 ##
##                                                                   ##
##     see the INSTALL file for detailed installation instructions   ##
##-------------------------------------------------------------------##

.KEEP_STATE:

KNOWN_SYSTEM=no
SRCDIR=../src
BINDIR=../bin
VPATH=$(SRCDIR) $(SRCDIR)/H5Part ../giza/interface ../giza/src ../utils ../tests
SHELL=/bin/bash

#
# some default settings for unix systems
#
F90C= 
F90FLAGS=
#
# change the line below if SPLASH does not find the X11 libraries
# (some settings of the SYSTEM variable for specific machines overwrite this)
#
# X11LIBS= -L/usr/X11R6/lib64 -lX11
X11LIBS= -L/usr/X11R6/lib -lX11

#--------------------------------------------------------------
# here we choose the backend plotting library for splash
# (default is PGPLOT for splash v1.x)
#--------------------------------------------------------------
ifeq ($(BACKEND), giza)
#-- C A I R O --
   PLOTLIB= ${GIZAOBJ} giza-fortran.F90 plotlib_giza.f90
   PGPLOTLIBS = -lcairo
   CFLAGS=-ggdb -Wall -Wextra -I../giza/src/ -I/usr/X11/include
   CC = gcc
   GIZASRC=${shell cd ../giza/src/; ls *.c}
   GIZAOBJ=${GIZASRC:.c=.o}
else
#-- P G P L O T --
   PLOTLIB= plotlib_pgplot.f90
#
# change the line below depending on where/how you have installed PGPLOT
# (some settings of the SYSTEM variable for specific machines overwrite this)
#
   PGPLOTLIBS = -L$(PGPLOT_DIR) -lpgplot -lpng
endif
#--------------------------------------------------------------

#
# If you need the HDF5 libraries, edit the lines below
# possibly adding a -L/libpath/ and a -I/includepath/
#
HDF5LIBS    = -L$(HDF5ROOT)/lib -lhdf5
HDF5INCLUDE = -I$(HDF5ROOT)/include
H5PART_LIBS    = -L$(H5PART_DIR)/lib -lpH5Part
H5PART_INCLUDE = -I$(H5PART_DIR)/include
#
# this file contains system-dependent routines like getarg, iargc etc.
#
SYSTEMFILE= system_f2003.f90
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
#
# destination for installed binaries
#
DEST=/usr/local/bin
#
# default C compiler
#
CC = gcc
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
   DEBUGFLAG= -Wall -Wextra -Wno=165 -g -fbounds-check -ftrace=full -freal=NaN
   ENDIANFLAGBIG= -fendian='BIG'
   ENDIANFLAGLITTLE= -fendian='LITTLE'
   PARALLEL= no
   CC = gcc
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), gfortran)
#  gfortran compiler (part of gcc 4.x.x)
   F90C= gfortran
   F90FLAGS= -O3 -Wall -frecord-marker=4
   SYSTEMFILE= system_f2003.f90
   DEBUGFLAG= -g -frange-check -fbacktrace
   OMPFLAGS= -fopenmp
   ENDIANFLAGBIG= -fconvert=big-endian
   ENDIANFLAGLITTLE= -fconvert=little-endian
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), gfortran41)
#  gfortran compiler (v4.1)
   F90C= gfortran
   F90FLAGS= -O3 -Wall -frecord-marker=4
   SYSTEMFILE= system_unix.f90
   DEBUGFLAG= -g -frange-check
   OMPFLAGS=
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
   CC=icc
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
   F90C= xlf90_r #-Wl,-t
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
   F90FLAGS= -O3 -mcmodel=medium -axT -warn all #-ipo #-assume nounderscore
   OMPFLAGS= -openmp
   DEBUGFLAG= -check all -traceback -g
   SYSTEMFILE= system_f2003.f90
   ENDIANFLAGBIG= -convert big_endian
   ENDIANFLAGLITTLE= -convert little_endian
   X11LIBS=-L/home/djp212/lib -lX11
   PGPLOTLIBS= -lpng #-L${PGPLOT_DIR} -lpgplot -lpng
   STATICLIBS= /home/djp212/pgplot/libpgplot.a
# or use setenv F_UFMTENDIAN=big or little at runtime
  # EXT=-dev
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

ifeq ($(SYSTEM),cody)
    F90C= gfortran
    F90FLAGS= -O3 -Wall
    SYSTEMFILE= system_f2003.f90
    DEBUGFLAG= -g -frange-check
    OMPFLAGS= -fopenmp
    ENDIANFLAGBIG= -fconvert=big-endian
    ENDIANFLAGLITTLE= -fconvert=little-endian
    SNFLAGS= -L$(HOME)/tree16/Objfiles/g5 -lsw
    CFLAGS = -g -O2 -Wall -I$(HOME)/tree16/include -fbounds-check
    CC = gcc
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
else
   ifeq ($(OPENMP),yes)
      F90FLAGS += $(OMPFLAGS)
   endif
endif

# add debugging flags at compile time
ifeq ($(DEBUG),yes)
    F90FLAGS += $(DEBUGFLAG)
endif

# rename the executable for development work
ifeq ($(DEV),yes)
    EXT= -dev
endif

# link with hdf5 libraries
ifeq ($(HDF5),yes)
    CFLAGS  += $(HDF5INCLUDE)
    LDFLAGS += $(HDF5LIBS)
endif

#
# MPI... no splash doesn't use it
# but sometimes you need to compile
# with the mpi compiler (e.g. to link
# correctly to MPI-HDF5 libraries)
# This just changes the compiler name
# whilst keeping the flags the same
# if MPI is set to "yes"
#
ifeq ($(MPI),yes)
    F90C= mpif90 -DPARALLEL_IO
    CC=mpicc -DPARALLEL_IO
endif

#
# If PGPLOT was compiled with a different compiler to the one used here,
# need to link to the libraries for that compiler. We attempt to do this
# automatically below by looking in the PGPLOT makefile.
#
# If the relevant library is not found, may also need -L/dir/ for the directory
# where the corresponding library is located (e.g. -L/usr/local/gfortran/lib -lgfortran)
#
# (information about what is done here is printed via the checkpgplot target, below)
#
ifneq ($(BACKEND),giza)
   ifdef PGPLOT_DIR

      PGPLOT_COMP=${shell grep 'FCOMPL=' $$PGPLOT_DIR/makefile | cut -d= -f2}

      ifneq (X$(F90C), X)  # make sure it is not just accidentally blank
      ifneq ($(PGPLOT_COMP), $(F90C))
# g77-compiled PGPLOT
        ifeq ($(PGPLOT_COMP), g77)
           PGPLOTLIBS+=-lg2c
        endif
# gfortran-compiled PGPLOT
        ifeq ($(PGPLOT_COMP), gfortran)
           PGPLOTLIBS+=-lgfortran
        endif
# g95-compiled PGPLOT
        ifeq ($(PGPLOT_COMP), g95)
           PGPLOTLIBS+=-lg95
        endif
      endif
      endif
   endif
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
%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@

# modules must be compiled in the correct order to check interfaces
# really should include all dependencies but I am lazy

SOURCESF90= $(PLOTLIB) globaldata.f90 asciiutils.f90 setpage.f90 transform.f90 \
         prompting.f90 geometry.f90 interpolation.f90 plotutils.f90 colourbar.f90 \
         colours.f90 colourparts.f90 pagecolours.f90 units.f90 \
         write_data_phantom.f90 write_pixmap.f90 \
         write_griddata.f90 write_sphdata.f90 \
         $(SYSTEMFILE) system_utils.f90 \
         discplot.f90 \
         fparser.f90 exact_function.f90 \
         exact_fromfile.f90 exact_mhdshock.f90 \
         exact_polytrope.f90 exact_rhoh.f90 \
         exact_sedov.f90 exact_shock.f90 exact_shock_sr.f90 exact_wave.f90 \
         exact_toystar1D.f90 exact_toystar2D.f90 \
         exact_densityprofiles.f90 exact_torus.f90 \
         exact_ringspread.f90 exact.f90 limits.f90 shapes.f90 \
         allocate.f90 titles.f90 \
         options_render.f90 options_particleplots.f90 \
         calc_quantities.f90 get_data.f90 \
         options_data.f90 \
         options_limits.f90 options_page.f90 \
	 options_powerspec.f90 \
	 options_vecplot.f90 options_xsecrotate.f90 pdfs.f90 \
         rotate.f90 interpolate1D.f90 interpolate2D.f90 \
         interpolate3D.f90 interpolate3D_xsec.f90 \
         interpolate3D_projection.F90 \
         interpolate3D_opacity.f90 interpolate_vec.f90 \
         interactive.f90 analysis.f90 \
         convert_grid.f90 convert.f90 \
         fieldlines.f90 legends.f90 particleplot.f90 \
         powerspectrums.f90 render.f90 \
         plotstep.f90 timestepping.f90 \
         defaults.f90 menu.f90 \
         splash.f90

# these are `external' f77 subroutines
SOURCESF= 

OBJECTS1 = $(SOURCESF:.f=.o) $(SOURCESF90:.f90=.o) $(STATICLIBS)
OBJECTS= $(OBJECTS1:.F90=.o)

#
# Now compile with the appropriate data read file
# (move yours to the top so that you can simply type "make")
#
all: ascii gadget vine sphNG srosph dragon tipsy
	@echo; echo ' SPLASH successfully compiled! ';
	@echo; echo ' Use "sudo make install" to copy the binaries to $(DEST)'; echo;

ascii: checksystem $(OBJECTS) read_data_ascii.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/asplash$(EXT) $(OBJECTS) read_data_ascii.o $(LDFLAGS) || ${MAKE} linkerror;
#--build universal binary on mac cluster
   ifeq ($(SYSTEM), maccluster)
	lipo -create $(BINDIR)/asplash_ppc $(BINDIR)/asplash_i386 -output $(BINDIR)/asplash || cp $(BINDIR)/asplash$(EXT) $(BINDIR)/asplash
   endif
	cp $(BINDIR)/asplash $(BINDIR)/splash

mbatesph: checksystem $(OBJECTS) read_data_mbate.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/bsplash $(OBJECTS) read_data_mbate.o $(LDFLAGS)


#---h5part reader---
H5PARTSRCC=H5PartF.c H5PartAttribF.c
H5PARTSRCF90=H5Part.f90 H5PartAttrib.f90
H5PARTOBJ=$(H5PARTSRCC:.c=.o) $(H5PARTSRCF90:.f90=.o)

h5part: checksystem checkh5part $(OBJECTS) $(H5PARTOBJ) read_data_h5part.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/h5splash $(OBJECTS) $(H5PARTOBJ) read_data_h5part.o $(LDFLAGS) $(H5PART_LIBS) $(HDF5LIBS)

read_data_h5part.o: read_data_h5part.f90
	$(F90C) $(F90FLAGS) $(H5PART_INCLUDE) -o $@ -c $<

H5PartF.o: H5PartF.c
	$(CC) -c $(CFLAGS) $(HDF5INCLUDE) $(H5PART_INCLUDE) $< -o $@

H5PartAttribF.o: H5PartAttribF.c
	$(CC) -c $(CFLAGS) $(HDF5INCLUDE) $(H5PART_INCLUDE) $< -o $@

#------------------

gadget: checksystem $(OBJECTS) read_data_gadget.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/gsplash $(OBJECTS) read_data_gadget.o $(LDFLAGS)

gadgetdualendian: checksystem $(OBJECTS) read_data_gadget.o read_data_gadget_otherendian.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/gsplash $(OBJECTS) read_data_gadget.o read_data_gadget_otherendian.o $(LDFLAGS)

read_data_gadget_otherendian.o: read_data_gadget.o
	cat read_data_gadget.f90 | awk "/subroutine read_data/,/end subroutine read_data/ { print }" | sed 's/subroutine read_data/subroutine read_data_otherendian/' > read_data_gadget_otherendian.f90 
	$(F90C) $(F90FLAGS) $(ENDIANFLAGBIG) -c read_data_gadget_otherendian.f90 -o read_data_gadget_otherendian.o

gadget_jsb: checksystem $(OBJECTS) read_data_gadget_jsb.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/gsplash-jsb $(OBJECTS) read_data_gadget_jsb.o $(LDFLAGS)

#---HDF5 gadget read ---
gadgethdf5: gadget_hdf5

gadget_hdf5: checksystem checkhdf5 $(OBJECTS) read_data_gadget_hdf5_utils.o read_data_gadget_hdf5.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/gsplash-hdf5 $(OBJECTS) read_data_gadget_hdf5_utils.o read_data_gadget_hdf5.o $(LDFLAGS) $(HDF5LIBS)

read_data_gadget_hdf5_utils.o: read_data_gadget_hdf5_utils.c
	$(CC) -c $(CFLAGS) $(HDF5INCLUDE) $(H5PART_INCLUDE) $< -o $@

#------------------

bauswein: checksystem $(OBJECTS) read_data_bauswein.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/bsplash $(OBJECTS) read_data_bauswein.o $(LDFLAGS)

dragon: checksystem $(OBJECTS) read_data_dragon.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/dsplash $(OBJECTS) read_data_dragon.o $(LDFLAGS)

vine: checksystem $(OBJECTS) read_data_VINE.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/vsplash $(OBJECTS) read_data_VINE.o $(LDFLAGS)

kitp: checksystem $(OBJECTS) read_data_kitp.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $(BINDIR)/ksplash $(OBJECTS) read_data_kitp.o

ndspmhd: checksystem $(OBJECTS) read_data_dansph.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/nsplash $(OBJECTS) read_data_dansph.o $(LDFLAGS)

dansph: checksystem $(OBJECTS) read_data_dansph_old.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/dsplash $(OBJECTS) read_data_dansph_old.o $(LDFLAGS)

foulkes: checksystem $(OBJECTS) read_data_foulkes.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/fsplash $(OBJECTS) read_data_foulkes.o $(LDFLAGS)

flash_hdf5: checksystem checkhdf5 $(OBJECTS) read_data_flash_hdf5.o read_data_flash_hdf5_utils.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/fsplash $(OBJECTS) read_data_flash_hdf5_utils.o read_data_flash_hdf5.o $(LDFLAGS)

flashhdf5:
	$(MAKE) flash_hdf5 HDF5=yes

read_data_flash_hdf5_utils.o: read_data_flash_hdf5_utils.c
	$(CC) -c $(CFLAGS) $(HDF5INCLUDE) $(H5PART_INCLUDE) $< -o $@

flash: flashhdf5

jjm: checksystem $(OBJECTS) read_data_jjm.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/jsplash $(OBJECTS) read_data_jjm.o $(LDFLAGS)

jules: checksystem $(OBJECTS) read_data_jules.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/jsplash $(OBJECTS) read_data_jules.o $(LDFLAGS)

oilonwater: checksystem $(OBJECTS) read_data_oilonwater.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/osplash $(OBJECTS) read_data_oilonwater.o $(LDFLAGS)

RSPH: rsph

rsph: checksystem $(OBJECTS) read_data_rsph.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/rsplash $(OBJECTS) read_data_rsph.o $(LDFLAGS)

scwsph: checksystem $(OBJECTS) read_data_scw.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/wsplash $(OBJECTS) read_data_scw.o $(LDFLAGS)

snsplash: snsph
snsph: $(OBJECTS) read_data_snsph.o read_data_snsph_utils.o
	$(FC) $(FFLAGS) -o $(BINDIR)/snsplash $(OBJECTS) read_data_snsph.o read_data_snsph_utils.o $(LDFLAGS) $(SNFLAGS) 

srosph: checksystem $(OBJECTS) read_data_sro.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/rsplash $(OBJECTS) read_data_sro.o $(LDFLAGS)

spyros: checksystem $(OBJECTS) read_data_spyros.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/ssplash $(OBJECTS) read_data_spyros.o $(LDFLAGS)

sphNG: checksystem $(OBJECTS) read_data_sphNG.o read_data_sphNG_otherendian.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/ssplash$(EXT) $(OBJECTS) read_data_sphNG.o read_data_sphNG_otherendian.o $(LDFLAGS)
#--build universal binary on mac cluster
   ifeq ($(SYSTEM), maccluster)
	lipo -create $(BINDIR)/ssplash_ppc $(BINDIR)/ssplash_i386 -output $(BINDIR)/ssplash || cp $(BINDIR)/ssplash$(EXT) $(BINDIR)/ssplash
   endif

read_data_sphNG_otherendian.o: read_data_sphNG.o
	@echo "!--This file is automatically generated during the make: do not edit" > $(SRCDIR)/${@:.o=.f90}
	cat $(SRCDIR)/read_data_sphNG.f90 | awk "/subroutine read_data/,/end subroutine read_data/ { print }" | sed 's/subroutine read_data/subroutine read_data_otherendian/' >> $(SRCDIR)/${@:.o=.f90}
	$(F90C) $(F90FLAGS) $(ENDIANFLAGBIG) -c $(SRCDIR)/${@:.o=.f90} -o $@

jjmmulti: checksystem $(OBJECTS) read_data_multiphase.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/jsplash $(OBJECTS) read_data_jjm_multiphase.o $(LDFLAGS)

tipsy: checksystem $(OBJECTS) read_data_tipsy.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/tsplash $(OBJECTS) read_data_tipsy.o $(LDFLAGS)

vanaverbeke: checksystem $(OBJECTS) read_data_vanaverbeke.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/vsplash $(OBJECTS) read_data_vanaverbeke.o $(LDFLAGS) 

ucla: checksystem $(OBJECTS) read_data_UCLA.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/usplash $(OBJECTS) read_data_UCLA.o $(LDFLAGS)

urban: checksystem $(OBJECTS) read_data_urban.o
	$(F90C) $(F90FLAGS) -o $(BINDIR)/usplash $(OBJECTS) read_data_urban.o $(LDFLAGS)

sky: ucla

steve: foulkes

sigfried: vanaverbeke

gasoline: tipsy

myall: ndspmhd dansph sphNG srosph gadget mbatesph tipsy ascii

SRCGRID2PDF= $(SYSTEMFILE) globaldata.f90 prompting.f90 transform.f90 \
             asciiutils.f90 write_griddata.f90 pdfs.f90 grid2pdf.f90
OBJGRID2PDF= ${SRCGRID2PDF:.f90=.o}
grid2pdf: checksystem $(OBJGRID2PDF)
	$(F90C) $(F90FLAGS) -o $(BINDIR)/$@ $(OBJGRID2PDF)


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
	@echo " make: ERROR: value of SYSTEM=$(SYSTEM) not recognised..."
	@echo " => set the environment variable SYSTEM to one listed "
	@echo "    in build/Makefile and try again"
	@echo ""
	@${MAKE} compilers
	@$(MAKE) err;
   endif

compilers:
	@echo "I suggest one of the following, based on detected Fortran compilers..."; echo;
	@if type -p ifort > /dev/null; then echo "make SYSTEM=ifort"; fi;
	@if type -p pathf90 > /dev/null; then echo "make SYSTEM=pathf90"; fi;
	@if type -p pgf90 > /dev/null; then echo "make SYSTEM=pgf90"; fi;
	@if type -p xlf90_r > /dev/null; then echo "make SYSTEM=ukaff1a [uses xlf90_r]"; fi;
	@if type -p gfortran > /dev/null; then echo "make SYSTEM=gfortran"; fi;
	@if type -p g95 > /dev/null; then echo "make SYSTEM=g95"; fi;
	@echo "(end of possible selections)"; echo;

checkpgplot:
   ifeq (X${PGPLOT_DIR}, X)
	@echo; echo "ERROR: PGPLOT_DIR should be set before compiling splash"; echo; ${MAKE} err;
   else
	@if [ -d $$PGPLOT_DIR ]; then echo; echo "PGPLOT_DIR=$$PGPLOT_DIR"; echo; else echo; echo "ERROR: Directory given by PGPLOT_DIR=$$PGPLOT_DIR does not exist"; echo; ${MAKE} err; fi;
   endif
   ifneq ($(PGPLOT_COMP),$(F90C))
	@echo; echo "*** WARNING: PGPLOT appears to have been compiled with a different Fortran"; echo "    compiler (${PGPLOT_COMP}) to the one you are using to compile SPLASH (${F90C}),"; echo "    so may need to link to the relevant compiler libraries ***";
# g77-compiled PGPLOT
	@if [ "${PGPLOT_COMP}" = "g77" ]; then echo "    [Adding -lg2c to the link flags for g77-compiled PGPLOT]"; fi;
# gfortran-compiled PGPLOT
	@if [ "${PGPLOT_COMP}" = "gfortran" ]; then echo "    [Adding -lgfortran to the link flags for gfortran-compiled PGPLOT]"; fi;
# g95-compiled PGPLOT
	@if [ "${PGPLOT_COMP}" = "g95" ]; then echo "    [Adding -lg95 to the link flags for g95-compiled PGPLOT]"; fi;
	@echo
#   else
#	@echo "PGPLOT was compiled with ${PGPLOT_COMP}";
   endif

checkh5part: checkhdf5
   ifeq (X${H5PART_DIR}, X)
	@echo; echo "ERROR: H5PART_DIR should be set before compiling splash with h5part read/write"; echo; ${MAKE} err;
   else
	@if [ -d $$H5PART_DIR ]; then echo; echo "H5PART_DIR=$$H5PART_DIR"; echo; else echo; echo "ERROR: Directory given by H5PART_DIR=$$H5PART_DIR does not exist"; echo; ${MAKE} err; fi;
   endif

checkhdf5:
   ifeq (X${HDF5ROOT}, X)
	@echo; echo "ERROR: HDF5ROOT should be set before compiling splash with HDF5 utilities"; echo; ${MAKE} err;
   else
	@if [ -d $$HDF5ROOT ]; then echo; echo "HDF5ROOT=$$HDF5ROOT"; echo; else echo; echo "ERROR: Directory given by HDF5ROOT=$$HDF5ROOT does not exist"; echo; ${MAKE} err; fi;
   endif

linkerror:
	@echo;
	@echo '------------- failed to link to relevant libraries required for SPLASH ------------------';
	@echo " Here are some common problems:";
	@echo;
        ifneq ($(BACKEND),giza)
	@echo " \"ld: library not found for -lpgplot\" ";
	@echo "       => Your PGPLOT_DIR is not set correctly, or PGPLOT is not installed correctly";
	@echo "          (linker needs to find one of libpgplot.a, libpgplot.so or libpgplot.dylib (OS/X) in PGPLOT_DIR)";
	@echo "          (check your PGPLOT installation by running the demo programs pgdemo1, pgdemo2 etc.)";
	@echo;
	@echo " \"ld: library not found for -lpng\" ";
	@echo "       => The png library was not found (most likely your PGPLOT does not use them either)";
	@echo "          (try removing -lpng from the PGPLOTLIBS= variable in build/Makefile)";
	@echo;
        endif
	@echo " \"ld: symbol(s) not found: "__g95_...", referenced from: blah.o \" ";
	@echo "       => Either you have changed compiler and need to \"make clean\", then \"make\" again";
	@echo "          or you need to add -lg95 to the PGPLOTLIBS= variable";
	@echo;
	@echo " \"ld: symbol(s) not found: "__gfortran_...", referenced from: blah.o \" ";
	@echo "       => Either you have changed compiler and need to \"make clean\", then \"make\" again";
	@echo "          or you need to add -lgfortran to the PGPLOTLIBS= variable";
	@echo;
	@echo " Undefined symbols:";
	@echo "   \"_main\", referenced from:";
	@echo "       start in crt1.10.5.o";
	@echo "       => You are compiling with gfortran but -lgfortran has also been added to the link flags";
	@echo "          (remove -lgfortran from the PGPLOTLIBS= variable)";
	@echo;
	@echo '-----------------------------------------------------------------------------------------';
	@${MAKE} err;

#
# install option, copies any binaries compiled to /usr/local/bin/
# run `make' first, then `make install'. Could in principle
# have install compile it as well, but environment variables
# will not be defined if "make" is run using sudo, so better
# to do the two separately
#
install: destcheck installcheck
	@for x in $(BINDIR)/?splash splash; do echo "copying $$x -> $(DEST)/$$x"; cp $$x $(DEST); done;
	@echo; echo 'installation complete';

installcheck:
	@if test -e $(BINDIR)/asplash && test -e $(BINDIR)/gsplash; then echo; \
        echo 'compiled binaries install to $(DEST)'; \
        echo '(use "sudo make install" if Permission denied)'; \
        echo; else echo;\
        echo 'run "make" first, followed by "sudo make install"'; echo;\
        $(MAKE) err; fi

destcheck: installcheck
	@if test -d $(DEST); then echo $(DEST) exists and is a directory; else \
        echo; echo "*** ERROR in make install ***"; echo "$(DEST) is not a valid directory"; echo;\
        $(MAKE) err; fi;

installclean: destcheck
	@for x in $(DEST)/?splash; do rm $$x; done;

distclean: installclean

cleanall: clean installclean

err:
	$(error aborting);

## other stuff

plotlib_pgplot.o: checkpgplot

giza-fortran.o : compilegiza

compilegiza: ${GIZAOBJ}
	@echo "Compiling giza... $<";

doc:
	cd docs; pdflatex splash; pdflatex splash

tar:
	tar cf splash.tar Makefile $(SOURCESF90) $(SOURCESF) read_data*.f90

targz:
	tar cf splash.tar Makefile $(SOURCESF90) $(SOURCESF) read_data*.f90
	gzip splash.tar

## unit tests of various modules as I write them

tests: test1 test2 test3 test_slicer test_prompt

test1: interpolate3D_projection.o interpolate3D_xsec.o test_interpolate3D.o 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $(BINDIR)/test_interpolation3D test_interpolate3D.o interpolate3D_projection.o interpolate3D_xsec.o
test2: transform.o ./tests/test_transform.o 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $(BINDIR)/test_transform ./tests/test_transform.o transform.o
test3: fieldlines.o ./tests/test_fieldlines.o 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $(BINDIR)/test_fieldlines ./tests/test_fieldlines.o fieldlines.o
test_slicer: interpolate3D_projection.o interpolate3D_xsec.o test_slicer3D.o 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $(BINDIR)/test_slicer3D test_slicer3D.o interpolate3D_projection.o interpolate3D_xsec.o
test_prompt: prompting.o test_prompting.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $(BINDIR)/test_prompt prompting.o test_prompting.o

#
#--code dependencies: MAY BE INCOMPLETE - I generate this automatically every so often
#
globaldata.o multiplot.mod settings_data.mod labels.mod filenames.mod particle_data.mod params.mod: globaldata.f90
asciiutils.o asciiutils.mod: asciiutils.f90
setpage.o pagesetup.mod: setpage.f90
transform.o transforms.mod: transform.f90
prompting.o prompting.mod: prompting.f90
geometry.o geometry.mod: geometry.f90
plotutils.o plotutils.mod: plotutils.f90
colourbar.o colourbar.mod: pagesetup.mod colourbar.f90
colours.o colours.mod: colours.f90
colourparts.o colourparts.mod: colourparts.f90
shapes.o shapes.mod: transforms.mod prompting.mod params.mod shapes.f90
units.o settings_units.mod: settings_data.mod labels.mod prompting.mod params.mod units.f90
write_pixmap.o write_pixmap.mod: colours.mod filenames.mod write_pixmap.f90
write_sphdata.o write_sphdata.mod: params.mod settings_units.mod labels.mod write_sphdata.f90
analysis.o analysis.mod: prompting.mod params.mod filenames.mod asciiutils.mod labels.mod analysis.f90
discplot.o disc.mod: transforms.mod discplot.f90
exact_fromfile.o exactfromfile.mod: exact_fromfile.f90
exact_mhdshock.o mhdshock.mod: exact_mhdshock.f90
exact_polytrope.o polytrope.mod: exact_polytrope.f90
exact_rhoh.o rhoh.mod: exact_rhoh.f90
exact_sedov.o sedov.mod: exact_sedov.f90
exact_shock.o shock.mod: exact_shock.f90
exact_shock_sr.o shock_sr.mod: exact_shock_sr.f90
exact_wave.o wave.mod: exact_wave.f90
exact_toystar1D.o toystar1d.mod: exact_toystar1D.f90
exact_toystar2D.o toystar2d.mod: exact_toystar2D.f90
exact_densityprofiles.o densityprofiles.mod: exact_densityprofiles.f90
exact_torus.o torus.mod: exact_torus.f90
exact_ringspread.o ringspread.mod: exact_ringspread.f90
exact.o exact.mod: transforms.mod ringspread.mod densityprofiles.mod wave.mod toystar2d.mod toystar1d.mod torus.mod shock_sr.mod shock.mod sedov.mod rhoh.mod polytrope.mod mhdshock.mod exactfromfile.mod labels.mod filenames.mod prompting.mod settings_data.mod exact.f90
limits.o limits.mod: asciiutils.mod settings_data.mod particle_data.mod geometry.mod labels.mod params.mod limits.f90
allocate.o mem_allocation.mod: particle_data.mod allocate.f90
titles.o titles.mod: filenames.mod asciiutils.mod titles.f90
system_f2003.o system_commands.mod: system_f2003.f90
system_utils.o system_utils.mod: system_commands.mod system_utils.f90
options_render.o settings_render.mod: params.mod prompting.mod colours.mod labels.mod colourbar.mod options_render.f90
options_particleplots.o settings_part.mod: multiplot.mod geometry.mod prompting.mod particle_data.mod settings_render.mod limits.mod labels.mod exact.mod settings_data.mod params.mod options_particleplots.f90
calc_quantities.o calcquantities.mod: settings_units.mod mem_allocation.mod settings_part.mod settings_data.mod particle_data.mod labels.mod calc_quantities.f90
get_data.o getdata.mod: system_utils.mod asciiutils.mod settings_units.mod calcquantities.mod geometry.mod labels.mod prompting.mod particle_data.mod settings_part.mod settings_data.mod limits.mod filenames.mod exact.mod get_data.f90
convert.o convert.mod: prompting.mod getdata.mod analysis.mod write_sphdata.mod filenames.mod settings_data.mod particle_data.mod convert.f90
options_data.o options_data.mod: settings_units.mod labels.mod limits.mod calcquantities.mod getdata.mod prompting.mod filenames.mod params.mod settings_data.mod options_data.f90
options_limits.o settings_limits.mod: transforms.mod labels.mod limits.mod prompting.mod calcquantities.mod settings_data.mod filenames.mod multiplot.mod options_limits.f90
options_page.o settings_page.mod: filenames.mod prompting.mod params.mod shapes.mod settings_limits.mod options_page.f90
options_powerspec.o settings_powerspec.mod: prompting.mod labels.mod limits.mod settings_data.mod options_powerspec.f90
options_vecplot.o settings_vecplot.mod: limits.mod labels.mod settings_data.mod prompting.mod options_vecplot.f90
options_xsecrotate.o settings_xsecrot.mod: transforms.mod multiplot.mod calcquantities.mod settings_data.mod prompting.mod limits.mod labels.mod filenames.mod options_xsecrotate.f90
pdfs.o pdfs.mod: asciiutils.mod filenames.mod plotutils.mod transforms.mod prompting.mod pdfs.f90
rotate.o rotation.mod: rotate.f90
interpolate1D.o interpolations1d.mod: interpolate1D.f90
interpolate2D.o interpolations2d.mod: interpolate2D.f90
interpolate3D.o interpolations3d.mod: interpolate3D.f90
interpolate3D_xsec.o xsections3d.mod: interpolate3D_xsec.f90
interpolate3D_projection.o projections3d.mod: interpolate3D_projection.F90
interpolate3D_opacity.o interpolate3d_opacity.mod: colours.mod projections3d.mod interpolate3D_opacity.f90
interpolate_vec.o interpolate_vec.mod: interpolate_vec.f90
interactive.o interactive_routines.mod: colours.mod settings_part.mod labels.mod calcquantities.mod filenames.mod transforms.mod settings_limits.mod settings_data.mod limits.mod settings_vecplot.mod settings_page.mod multiplot.mod shapes.mod settings_xsecrot.mod colourbar.mod interactive.f90
fieldlines.o fieldlines.mod: fieldlines.f90
legends.o legends.mod: legends.f90
particleplot.o particleplots.mod: geometry.mod interpolations2d.mod settings_part.mod settings_data.mod labels.mod params.mod particleplot.f90
powerspectrums.o powerspectrums.mod: interpolations3d.mod powerspectrums.f90
render.o render.mod: settings_vecplot.mod legends.mod plotutils.mod colourbar.mod render.f90
plotstep.o timestep_plotting.mod: asciiutils.mod rotation.mod fieldlines.mod interpolate_vec.mod geometry.mod shapes.mod legends.mod colourbar.mod write_pixmap.mod exactfromfile.mod disc.mod pagesetup.mod render.mod xsections3d.mod interpolate3d_opacity.mod interpolations3d.mod interpolations2d.mod interpolations1d.mod powerspectrums.mod particleplots.mod interactive_routines.mod transforms.mod colourparts.mod settings_units.mod settings_vecplot.mod settings_limits.mod toystar2d.mod toystar1d.mod exact.mod filenames.mod pdfs.mod projections3d.mod particle_data.mod settings_powerspec.mod settings_xsecrot.mod settings_render.mod settings_part.mod settings_page.mod settings_data.mod titles.mod prompting.mod multiplot.mod limits.mod labels.mod colours.mod params.mod plotstep.f90
timestepping.o timestepping.mod: params.mod settings_part.mod getdata.mod timestep_plotting.mod settings_page.mod settings_data.mod particle_data.mod filenames.mod timestepping.f90
defaults.o defaults.mod: shapes.mod titles.mod settings_units.mod settings_powerspec.mod settings_xsecrot.mod settings_vecplot.mod settings_render.mod settings_part.mod options_data.mod settings_limits.mod multiplot.mod exact.mod settings_page.mod settings_data.mod particle_data.mod limits.mod labels.mod filenames.mod defaults.f90
menu.o mainmenu.mod: settings_units.mod params.mod timestepping.mod getdata.mod geometry.mod defaults.mod transforms.mod prompting.mod multiplot.mod settings_xsecrot.mod settings_vecplot.mod settings_render.mod settings_page.mod settings_part.mod settings_limits.mod settings_data.mod options_data.mod limits.mod labels.mod filenames.mod menu.f90
splash.o: settings_page.mod timestepping.mod analysis.mod write_sphdata.mod convert.mod write_pixmap.mod asciiutils.mod system_utils.mod system_commands.mod settings_xsecrot.mod settings_data.mod projections3d.mod mem_allocation.mod mainmenu.mod limits.mod defaults.mod getdata.mod filenames.mod splash.f90
read_data_sphNG.o sphngread.mod: settings_units.mod geometry.mod calcquantities.mod labels.mod system_utils.mod mem_allocation.mod settings_data.mod particle_data.mod params.mod read_data_sphNG.f90


clean:
	rm -f *.o *.mod
