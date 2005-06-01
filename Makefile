##-------------------------------------------------------------------##
##     Makefile for compiling supersphplot (uses PGPLOT)             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

## Compiler options
F95C =  g95
F95FLAGS =  -O -Wall -fbounds-check

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
FC = $(F95C)
FFLAGS = $(F95FLAGS)

# define the implicit rule to make a .o file from a .f90/.f95 file
# (some Make versions don't know this)

%.o : %.f90
	$(F95C) $(F95FLAGS) -c $< -o $@
%.o : %.f95
	$(F95C) $(F95FLAGS) -c $< -o $@

# modules must be compiled in the correct order to check interfaces

# these need to be f95 compatible
SOURCESF95= exact_shock.f95
# most are just f90 
SOURCESF90= globaldata.f90 transform.f90 \
         prompting.f90 geometry.f90 \
         colours.f90 colourparts.f90 \
         exact_fromfile.f90 exact_mhdshock.f90 \
         exact_polytrope.f90 exact_rhoh.f90 \
         exact_sedov.f90 exact_wave.f90 \
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

OBJECTS = $(SOURCESF:.f=.o) $(SOURCESF95:.f95=.o) $(SOURCESF90:.f90=.o) 

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

dansph: $(OBJECTS) read_data_dansph.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o supersphplot $(OBJECTS) read_data_dansph.o

nina: $(OBJECTS) read_data_nina.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o nsupersphplot $(OBJECTS) read_data_nina.o

spmhd: $(OBJECTS) read_data_mbate_mhd.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o msupersphplot $(OBJECTS) read_data_mbate_mhd.o

sinksph: $(OBJECTS) read_data_mbate_hydro.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o hdsupersphplot $(OBJECTS) read_data_mbate_hydro.o

scwsph: $(OBJECTS) read_data_scw.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o wsupersphplot $(OBJECTS) read_data_scw.o

srosph: $(OBJECTS) read_data_sro.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o rsupersphplot $(OBJECTS) read_data_sro.o

spyros: $(OBJECTS) read_data_spyros.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o ssupersphplot $(OBJECTS) read_data_spyros.o


## other stuff

doc:
	cd docs; latex supersphplot; latex supersphplot; dvips supersphplot -o supersphplot.ps; ps2pdf13 supersphplot.ps

tar:
	tar cf supersphplot.tar Makefile $(MODULES) $(SOURCES) read_data*.f90

targz:
	tar cf supersphplot.tar Makefile $(MODULES) $(SOURCES) read_data*.f90
	gzip supersphplot.tar

clean:
	rm *.o *.mod
