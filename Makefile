##-------------------------------------------------------------------##
##     Makefile for compiling supersphplot (uses PGPLOT)             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

## Compiler options
F90C =  g95
#F90C = f95
F90FLAGS =  -O -Wall  -fbounds-check
#LDFLAGS = -L/usr/X11R6/lib -lX11 -lpgplot \
#         -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2/ -lg2c \
#         -lpng

LDFLAGS =  -L/usr/X11R6/lib -lX11 -L/sw/lib -lpng -laquaterm -lcc_dynamic -Wl,-framework -Wl,Foundation -L/sw/lib/pgplot95 -lpgplot
SYSTEMFILE = system_unix.f90
#SYSTEMFILE = system_unix_NAG.f90

# Fortran flags same as F90
FC = $(F90C)
FFLAGS = $(F90FLAGS)

# define the implicit rule to make a .o file from a .f90 file
# (some Make versions don't know this)

%.o : %.f90
	$(F90C) $(F90FLAGS) -c $< -o $@

# modules must be compiled in the correct order to check interfaces

SOURCESF90= globaldata.f90 transform.f90 \
         prompting.f90 geometry.f90 \
         colours.f90 colourparts.f90 \
         exact_fromfile.f90 exact_mhdshock.f90 \
         exact_polytrope.f90 exact_rhoh.f90 \
         exact_sedov.f90 exact_shock.f90 exact_wave.f90 \
         exact_toystar1D.f90 exact_toystar2D.f90 \
         limits.f90 options_limits.f90 \
         exact.f90 options_page.f90 \
         options_particleplots.f90 \
         allocate.f90 \
         calc_quantities.f90 get_data.f90\
         options_data.f90 \
	 options_powerspec.f90 options_render.f90 \
	 options_vecplot.f90 options_xsecrotate.f90 \
         rotate.f90 interpolate1D.f90 \
         interpolate2D.f90 interpolate3D_xsec.f90 \
         interpolate3D_projection.f90\
         interactive.f90 \
         fieldlines.f90 legends.f90 particleplot.f90 \
         powerspectrums.f90 render.f90 setpage.f90 \
	 titles.f90 \
         plotstep.f90 timestepping.f90 \
         defaults.f90 menu.f90 \
         $(SYSTEMFILE) supersphplot.f90

# these are `external' f77 subroutines
SOURCESF= danpgsch.f danpgtile.f danpgwedg.f

OBJECTS = $(SOURCESF:.f=.o) $(SOURCESF90:.f90=.o) 

#
# Now compile with the appropriate data read file
# (move yours to the top so that you can simply type "make")
#

dansph: $(OBJECTS) read_data_dansph.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o supersphplot $(OBJECTS) read_data_dansph.o

jjmsph: $(OBJECTS) read_data_jjmsph.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o jsupersphplot $(OBJECTS) read_data_jjmsph.o

mbatesph: $(OBJECTS) read_data_mbate.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o hsupersphplot $(OBJECTS) read_data_mbate.o

spmhd: $(OBJECTS) read_data_mbate_mhd.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o msupersphplot $(OBJECTS) read_data_mbate_mhd.o

sinksph: $(OBJECTS) read_data_mbate_hydro.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o hdsupersphplot $(OBJECTS) read_data_mbate_hydro.o

scwsph: $(OBJECTS) read_data_scw.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o wsupersphplot $(OBJECTS) read_data_scw.o

srosph: $(OBJECTS) read_data_sro.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o rsupersphplot $(OBJECTS) read_data_sro.o

gadget: $(OBJECTS) read_data_gadget.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o gsupersphplot $(OBJECTS) read_data_gadget.o

vine: $(OBJECTS) read_data_VINE.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o vsupersphplot $(OBJECTS) read_data_vine.o

## other crap

tar:
	tar cf supersphplot.tar Makefile $(MODULES) $(SOURCES) read_data*.f90

targz:
	tar cf supersphplot.tar Makefile $(MODULES) $(SOURCES) read_data*.f90
	gzip supersphplot.tar

clean:
	rm *.o *.mod
