##-------------------------------------------------------------------##
##     Makefile for compiling supersphplot (uses PGPLOT)             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

## Compiler options
F90C =  g95
F90FLAGS =  -O -Wall -fbounds-check
#LDFLAGS = -L/usr/X11R6/lib -lX11 -lpgplot \
#         -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2/ -lg2c \
#         -lpng

LDFLAGS =  -L/usr/X11R6/lib -lX11 -L/sw/lib -lpng -laquaterm -lcc_dynamic -Wl,-framework -Wl,Foundation -L/sw/lib/pgplot95 -lpgplot
#           -lcc_dynamic -Wl,-framework -Wl,Foundation \
#           -L/sw/lib/pgplot -lpgplot -lg2c
#LDFLAGS = -laquaterm -L/usr/X11R6/lib -lX11 -L/sw/lib -lpng -L/sw/lib/pgplot -lpgplot -lg2c
SYSTEMFILE = system_unix.f90

# Fortran flags same as F90
FC = $(F90C)
FFLAGS = $(F90FLAGS)

# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(F90C) $(F90FLAGS) -c $< -o $@

DANSPH = read_data_dansph.f90 
MRBSPH = read_data_mbate.f90
SCWSPH = read_data_scw.f90
SROSPH = read_data_sro.f90
JJMSPH = read_data_jjm.f90
GADGETSPH = read_data_gadget.f90

# put modules separately as these must be compiled before the others
MODULES= globaldata.f90 transform.f90 prompting.f90 \
         geometry.f90 colours.f90 colourparts.f90 limits.f90 rotate.f90 \
         interactive.f90 allocate.f90 \
         fieldlines.f90 legends.f90 particleplot.f90 \
         powerspectrums.f90 \
         toystar2D_utils.f90 exact_toystar2D.f90 \
         exact.f90 defaults.f90 plotstep.f90 timestepping.f90 \
         $(SYSTEMFILE)

# these are the normal `external' subroutines
SOURCES= supersphplot.f90 \
         calc_quantities.f90 \
	 danpgsch.f danpgtile.f danpgwedg.f \
	 exact_fromfile.f90 exact_rhoh.f90 \
	 exact_sedov.f90 exact_shock.f90 exact_wave.f90 \
	 exact_toystar.f90 \
	 exact_toystar_ACplane.f exact_mhdshock.f90 \
	 exact_polytrope.f \
	 get_data.f90 integratedkernel.f90 \
	 interpolate1D.f90 interpolate_vec.f90 \
         interpolate2D.f90 interpolate2D_xsec.f90 \
	 interpolate3D.f90 interpolate3D_fastxsec.f90 \
         interpolate3D_proj_vec.f90 \
	 interpolate3D_projection.f90 interpolate3D_xsec_vec.f90 \
	 menu.f90 options_data.f90 \
	 options_limits.f90 \
	 options_page.f90 options_particleplots.f90 \
	 options_powerspec.f90 options_render.f90 \
	 options_vecplot.f90 options_xsecrotate.f90 \
	 plot_kernel_gr.f90 \
	 print_header.f90\
         render.f90 render_vec.f90 \
	 setpage.f90 \
	 titles_read.f90 \

SOURCESALL = $(MODULES:.f90=.o) $(SOURCES:.f90=.o)

OBJJJMSPH = $(SOURCESALL:.f=.o) $(JJMSPH:.f90=.o)
OBJDANSPH = $(SOURCESALL:.f=.o) $(DANSPH:.f90=.o)
OBJMRBSPH = $(SOURCESALL:.f=.o) $(MRBSPH:.f90=.o)
OBJSCWSPH = $(SOURCESALL:.f=.o) $(SCWSPH:.f90=.o)
OBJSROSPH = $(SOURCESALL:.f=.o) $(SROSPH:.f90=.o)
OBJGADGETSPH = $(SOURCESALL:.f=.o) $(GADGETSPH:.f90=.o)

dansph: $(OBJDANSPH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o ../supersphplot $(OBJDANSPH)

jjmsph: $(OBJJJMSPH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o jsupersphplot $(OBJJJMSPH)

mrbsph: $(OBJMRBSPH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o supersphplot_mrb $(OBJMRBSPH)

scwsph: $(OBJSCWSPH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o supersphplot_scw $(OBJSCWSPH)

srosph: $(OBJSROSPH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o supersphplot_sro $(OBJSROSPH)

gadget: $(OBJGADGETSPH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o supersphplot_gadget $(OBJGADGETSPH)

## sort out dependencies on modules
defaults.o: exact.o

exact.o: toystar2D_utils.o exact_toystar2D.o
## other crap

tar:
	tar cf supersphplot.tar Makefile $(MODULES) $(SOURCES) read_data*.f90

targz:
	tar cf supersphplot.tar Makefile $(MODULES) $(SOURCES) read_data*.f90
	gzip supersphplot.tar

clean:
	rm *.o *.mod
