##-------------------------------------------------------------------##
##     Makefile for compiling supersphplot (uses PGPLOT)             ##
##     Written by Daniel Price					     ##
##     University of Exeter, UK, 2006                    	     ##
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
X11LIBS= -L/usr/X11R6/lib -lX11
#
# change the line below depending on where/how you have installed PGPLOT
# (some settings of the SYSTEM variable for specific machines overwrite this)
#
PGPLOTLIBS = -L$(PGPLOT_DIR) -lpgplot -lpng -lg2c
#
# this file contains system-dependent routines like getarg, iargc etc.
#
SYSTEMFILE= system_unix.f90
#
# this can be used to statically link the pgplot libraries if all else fails
#
STATICLIBS=
#
# set the parallel flag to yes to compile with openMP directives
#
PARALLEL= no
#
# the openMP flags should be set in the lines defining your system type
# (ie. leave line below blank)
OMPFLAGS=
#
# the endian flag can be used to compile the code to read BIG or LITTLE endian data
# some compilers also allow this to be done at runtime (e.g. g95, ifort) by setting
# an environment variable appropriately (e.g. G95_ENDIAN or F_UFMTENDIAN)
#
ENDIAN=
#ENDIAN='BIG'
#ENDIAN='LITTLE'

ifeq ($(SYSTEM),mymac)
#  these are the settings for a Mac with pgplot installed via fink
   F90C= g95
   F90FLAGS= -O3 -ffast-math
   PGPLOTLIBS= -L/sw/lib/pgplot -lpgplot -lg2c -L/sw/lib -lpng \
          -laquaterm -lcc_dynamic -Wl,-framework -Wl,Foundation
   SYSTEMFILE= system_f2003.f90
   PARALLEL= no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),g95)
#  using the g95 compiler
   F90C= g95
   F90FLAGS= -03 -ffast-math
   SYSTEMFILE= system_f2003.f90 # this is for Fortran 2003 compatible compilers
   ENDIANFLAGBIG= -fendian='BIG'
   ENDIANFLAGLITTLE= -fendian='LITTLE'
   PARALLEL= no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),nagf95)
#  NAG f95 compiler
   F90C= f95
   F90FLAGS= -03
   SYSTEMFILE= system_unix_NAG.f90
   PARALLEL= no
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),sunf95)
#  sun f95 compiler on linux
   F90C= sunf95
   F90FLAGS= -fast
   OMPFLAGS= -openmp
   SYSTEMFILE= system_f2003.f90
   ENDIANFLAGBIG= -xfilebyteorder=big16:%all
   ENDIANFLAGLITTLE= -xfilebyteorder=little16:%all
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),ifort)
#  this is for the intel fortran compiler
   F90C= ifort
   F90FLAGS= -O3 -Vaxlib -nbs
   OMPFLAGS= -openmp
   SYSTEMFILE= system_f2003.f90
   ENDIANFLAGBIG= -convert big_endian
   ENDIANFLAGLITTLE= -convert little_endian
# or use setenv F_UFMTENDIAN=big or little at runtime
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),pgf90)
#  this is for the Portland Group Fortran 90 compiler
   F90C= pgf90
   F90FLAGS= -O -Mbackslash
   SYSTEMFILE= system_unix.f90
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM),ukaff1a)
#  this is for ukaff1a
   F90C= xlf90_r
   F90FLAGS= -O5 -qnoescape -qsuffix=f=f90 -qextname
   OMPFLAGS= -qsmp=noauto
   SYSTEMFILE= system_f2003.f90
   PGPLOTLIBS= -L/home/dprice/pgplot -lpgplot  
   STATICLIBS= /home/dprice/pgplot/libpgplot.a /home/dprice/plot/libg2cmodified.a
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

ifeq ($(PARALLEL),yes)
    F90FLAGS += $(OMPFLAGS)
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
%.o : %.f95
	$(F90C) $(F90FLAGS) -c $< -o $@

#
# use either the parallel or serial versions of some routines
#
ifeq ($(PARALLEL),yes)
   INTERPROUTINES= interpolate3D_projection_P.f90
else
   INTERPROUTINES= interpolate3D_projection.f90
endif

# modules must be compiled in the correct order to check interfaces
# really should include all dependencies but I am lazy

SOURCESF90= globaldata.f90 transform.f90 \
         prompting.f90 geometry.f90 \
         colours.f90 colourparts.f90 \
         danpgutils.f90 \
         exact_fromfile.f90 exact_mhdshock.f90 \
         exact_polytrope.f90 exact_rhoh.f90 \
         exact_sedov.f90 exact_shock.f90 exact_wave.f90 \
         exact_toystar1D.f90 exact_toystar2D.f90 \
         exact_densityprofiles.f90 exact_torus.f90 \
         limits.f90 options_limits.f90 \
         exact.f90 options_page.f90 \
         options_particleplots.f90 \
         allocate.f90 titles.f90 \
         calc_quantities.f90 get_data.f90\
         options_data.f90 \
	 options_powerspec.f90 options_render.f90 \
	 options_vecplot.f90 options_xsecrotate.f90 \
         rotate.f90 interpolate1D.f90 \
         interpolate2D.f90 interpolate3D_xsec.f90 \
         $(INTERPROUTINES) \
         interpolate3D_opacity.f90\
         interactive.f90 \
         fieldlines.f90 legends.f90 particleplot.f90 \
         powerspectrums.f90 render.f90 setpage.f90 \
         plotstep.f90 timestepping.f90 \
         defaults.f90 menu.f90 \
         $(SYSTEMFILE) supersphplot.f90

# these are `external' f77 subroutines
SOURCESF=

OBJECTS = $(SOURCESF:.f=.o) $(SOURCESF90:.f90=.o) $(STATICLIBS)

#
# Now compile with the appropriate data read file
# (move yours to the top so that you can simply type "make")
#
ascii: checksystem $(OBJECTS) read_data_ascii.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o asupersphplot $(OBJECTS) read_data_ascii.o

mbatesph: checksystem $(OBJECTS) read_data_mbate.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o bsupersphplot $(OBJECTS) read_data_mbate.o

gadget: checksystem $(OBJECTS) read_data_gadget.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o gsupersphplot $(OBJECTS) read_data_gadget.o

vine: checksystem $(OBJECTS) read_data_VINE.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o vsupersphplot $(OBJECTS) read_data_vine.o

ndspmhd: checksystem $(OBJECTS) read_data_dansph.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o supersphplot $(OBJECTS) read_data_dansph.o

dansph: checksystem $(OBJECTS) read_data_dansph_old.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o dsupersphplot $(OBJECTS) read_data_dansph_old.o

scwsph: checksystem $(OBJECTS) read_data_scw.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o wsupersphplot $(OBJECTS) read_data_scw.o

srosph: checksystem $(OBJECTS) read_data_sro.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o rsupersphplot $(OBJECTS) read_data_sro.o 

spyros: checksystem $(OBJECTS) read_data_spyros.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ssupersphplot $(OBJECTS) read_data_spyros.o

sphNG: checksystem $(OBJECTS) read_data_sphNG.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ssupersphplot $(OBJECTS) read_data_sphNG.o

RSPH: rsph

rsph: checksystem $(OBJECTS) read_data_rsph.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o rsupersphplot $(OBJECTS) read_data_rsph.o


## some dependencies as I get around to it
defaults.o: ./titles.mod ./particle_data.mod ./settings_powerspec.mod \
            ./settings_xsecrot.mod ./settings_vecplot.mod ./settings_render.mod \
            ./settings_page.mod ./settings_part.mod ./settings_data.mod \
            ./options_data.mod ./settings_limits.mod ./multiplot.mod \
            ./limits.mod ./labels.mod ./filenames.mod ./exact.mod defaults.f90

exact.o: ./transforms.mod ./densityprofiles.mod ./wave.mod ./toystar2d.mod \
         ./toystar1d.mod ./shock.mod ./sedov.mod ./rhoh.mod ./polytrope.mod \
         ./mhdshock.mod ./exactfromfile.mod ./labels.mod ./filenames.mod \
         ./prompting.mod ./settings_data.mod exact.f90
         
options_render.o: colours.mod
interactive.o: colours.mod

checksystem:
   ifeq ($(KNOWN_SYSTEM), yes)
	@echo ""
	@echo "Compiling supersphplot for $(SYSTEM) system..........."
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
	@echo " Makefile for SUPERSPHPLOT by Daniel Price "
	@echo " -- see INSTALL file for detailed instructions"
	@echo ""
	@echo " make: WARNING: value of SYSTEM = $(SYSTEM) not recognised..."
	@echo " =>set the environment variable SYSTEM to one listed "
	@echo "   in the Makefile and try again"
	@echo ""
	quit
   endif

## other stuff

doc:
	cd docs; latex supersphplot; latex supersphplot; dvips supersphplot -o supersphplot.ps; ps2pdf13 supersphplot.ps

tar:
	tar cf supersphplot.tar Makefile $(SOURCESF90) $(SOURCESF) read_data*.f90

targz:
	tar cf supersphplot.tar Makefile $(SOURCESF90) $(SOURCESF) read_data*.f90
	gzip supersphplot.tar

## unit tests of various modules as I write them

tests: interpolate3D_projection.o interpolate3D_xsec.o ./tests/test_interpolate3D.o 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o test_interpolation3D ./tests/test_interpolate3D.o interpolate3D_projection.o interpolate3D_xsec.o
test2: transform.o ./tests/test_transform.o 
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o test_transform ./tests/test_transform.o transform.o

clean:
	rm *.o *.mod
