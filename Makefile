##-------------------------------------------------------------------##
##     Makefile for compiling supersphplot (uses PGPLOT)             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

FORTRAN = f90  
FC = f90

FFLAGS = -C -O4
F90FLAGS = -C -O4
LDFLAGS = -lpgplot -lX11 -lF77

DANSPH = read_data_dansph.f render_smooth.f
MRBSPH = read_data_mbate.f render_smooth.f
PGXTAL = read_data_dansph.f render_smooth_pgxtal.f

SOURCES= modules.f90 prompting.f90  supersphplot.f \
         calc_quantities.f90 colour_demo.f danpgwedg.f \
	 exact_rhoh.f90 \
	 exact_sedov.f exact_swave.f exact_toystar.f \
	 exact_toystar_ACplane.f exact_mhdshock.f90 \
	 get_render_options.f \
	 integratedkernel.f90 \
         interpolate2D.f interpolate3D.f \
	 interpolate3D_fastxsec.f \
	 interpolate3D_projection.f \
	 int_from_string.f90 \
	 legend.f \
	 lomb_powerspectrum1D.f90 \
	 menu_actions.f \
	 plot_average.f plot_powerspectrum.f90 \
	 print_header.f90\
	 print_menu.f read_defaults.f90 \
         render_coarse.f render.f \
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
