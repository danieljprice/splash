##-------------------------------------------------------------------##
##     Makefile for compiling supersphplot (uses PGPLOT)             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

## Compiler options
F90C =  g95
F90FLAGS =  -O2

LDFLAGS = -L/usr/X11R6/lib -lX11 -L/sw/lib/pgplot -lpgplot -lg2c -L/sw/lib -lpng \
          -laquaterm -lcc_dynamic -Wl,-framework -Wl,Foundation
#LDFLAGS = -L/usr/X11R6/lib -lX11 -lpgplot \
#         -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2/ -lg2c \
#         -lpng

# system file (top one uses Fortran 2003 system calls, as in g95)
SYSTEMFILE = system_f2003.f90 # this is for Fortran 2003 compatible compilers
#SYSTEMFILE = system_unix.f90
#SYSTEMFILE = system_unix_NAG.f90

# Fortran flags same as F90
FC = $(F90C)
FFLAGS = $(F90FLAGS)

# define the implicit rule to make a .o file from a .f90/.f95 file
# (some Make versions don't know this)

%.o : %.f90
	$(F90C) $(F90FLAGS) -c $< -o $@
%.o : %.f95
	$(F90C) $(F90FLAGS) -c $< -o $@

# modules must be compiled in the correct order to check interfaces

SOURCESF90= globaldata.f90 transform.f90 \
         prompting.f90 geometry.f90 \
         colours.f90 colourparts.f90 \
         exact_fromfile.f90 exact_mhdshock.f90 \
         exact_polytrope.f90 exact_rhoh.f90 \
         exact_sedov.f90 exact_shock.f90 exact_wave.f90 \
         exact_toystar1D.f90 exact_toystar2D.f90 \
         exact_densityprofiles.f90 \
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
         interpolate3D_projection.f90\
         interpolate3D_opacity.f90\
         interactive.f90 \
         fieldlines.f90 legends.f90 particleplot.f90 \
         powerspectrums.f90 render.f90 setpage.f90 \
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
ascii: $(OBJECTS) read_data_ascii.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o asupersphplot $(OBJECTS) read_data_ascii.o

mbatesph: $(OBJECTS) read_data_mbate.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o bsupersphplot $(OBJECTS) read_data_mbate.o

gadget: $(OBJECTS) read_data_gadget.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o gsupersphplot $(OBJECTS) read_data_gadget.o

vine: $(OBJECTS) read_data_VINE.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o vsupersphplot $(OBJECTS) read_data_vine.o

ndspmhd: $(OBJECTS) read_data_dansph.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o supersphplot $(OBJECTS) read_data_dansph.o

dansph: $(OBJECTS) read_data_dansph_old.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o dsupersphplot $(OBJECTS) read_data_dansph_old.o

scwsph: $(OBJECTS) read_data_scw.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o wsupersphplot $(OBJECTS) read_data_scw.o

srosph: $(OBJECTS) read_data_sro.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o supersphplot $(OBJECTS) read_data_sro.o 

spyros: $(OBJECTS) read_data_spyros.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o ssupersphplot $(OBJECTS) read_data_spyros.o

sphNG: $(OBJECTS) read_data_sphNG.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o ssupersphplot $(OBJECTS) read_data_sphNG.o


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

## other stuff

doc:
	cd docs; latex supersphplot; latex supersphplot; dvips supersphplot -o supersphplot.ps; ps2pdf13 supersphplot.ps

tar:
	tar cf supersphplot.tar Makefile $(MODULES) $(SOURCES) read_data*.f90

targz:
	tar cf supersphplot.tar Makefile $(MODULES) $(SOURCES) read_data*.f90
	gzip supersphplot.tar

## unit tests of various modules as I write them

tests: ./tests/test_interpolate3D.o interpolate3D_projection.o interpolate3D_xsec.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o test_interpolation3D ./tests/test_interpolate3D.o interpolate3D_projection.o interpolate3D_xsec.o

clean:
	rm *.o *.mod
