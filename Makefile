##-------------------------------------------------------------------##
##     Makefile for compiling supersphplot (uses PGPLOT)             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

FORTRAN = f90  
FC = lf95
F90C = lf95

FFLAGS = -O
F90FLAGS = -O
LDFLAGS = -lpgplot -L/usr/X11R6/lib -lX11

# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(F90C) -c $(F90FLAGS) $(FPPFLAGS) $< -o $@

DANSPH = read_data_dansph.f 
MRBSPH = read_data_mbate.f 
PGXTAL = read_data_dansph.f

SOURCES= modules.f90 prompting.f90 \
	 supersphplot.f90 main.f90 \
         calc_quantities.f90 colour_demo.f danpgwedg.f \
	 exact_rhoh.f90 \
	 exact_sedov.f exact_swave.f exact_toystar.f \
	 exact_toystar_ACplane.f exact_mhdshock.f90 \
	 integratedkernel.f90 \
	 interpolate1D.f \
         interpolate2D.f interpolate3D.f \
	 interpolate3D_fastxsec.f \
	 interpolate3D_projection.f \
	 int_from_string.f90 \
	 legend.f \
	 menu_actions.f \
	 options_exact.f90 options_powerspec.f90 \
	 options_render.f90 \
	 plot_average.f plot_powerspectrum.f90 \
	 powerspectrum_fourier1D.f90 \
	 powerspectrum_lomb1D.f90 \
	 print_header.f90\
	 print_menu.f read_defaults.f90 \
         render_coarse.f render.f90 \
	 setcolours.f set_defaults.f \
	 transform.f write_defaults.f90

SOURCESALL = $(SOURCES:.f90=.o)

OBJDANSPH = $(SOURCESALL:.f=.o) $(DANSPH:.f=.o)
OBJMRBSPH = $(SOURCESALL:.f=.o) $(MRBSPH:.f=.o)
OBJPGXTAL = /data/cass30c/dprice/pgxtal/src/pgxtal.o \
	    $(SOURCESALL:.f=.o) $(PGXTAL:.f=.o)

dansph: $(OBJDANSPH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o ../supersphplot $(OBJDANSPH)

mrbsph: $(OBJMRBSPH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o supersphplot_mrb $(OBJMRBSPH)

pgxtal: $(OBJPGXTAL)
	$(FC) $(FFLAGS) $(LDFLAGS) -o superpgxtal $(OBJPGXTAL)
	
tar:
	tar cf supersphplot.tar Makefile $(SOURCES) read_data*.f

targz:
	tar cf supersphplot.tar Makefile $(SOURCES) read_data*.f
	gzip supersphplot.tar

clean:
	rm *.o *.mod
