##-------------------------------------------------------------------##
##     Makefile for compiling supersphplot (uses PGPLOT)             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

## Compiler options (uncomment ONE set)
## ------------------------------------------------------------------ ##
## IoA compiler (Sun fortran)
#FC = f90
#F90C = f90
#FFLAGS = -O4
#F90FLAGS = -O4
#LDFLAGS = -lpgplot -L/usr/X11R6/lib -lX11
## ------------------------------------------------------------------ ##
## Monash compiler (Lahey-Fujitsu f95)
FC = lf95
F90C = lf95
FFLAGS = -O ##--chk aesux --chkglobal  --warn
F90FLAGS = -O ##--chk aesux --chkglobal --warn
LDFLAGS = -lpgplot -L/usr/X11R6/lib -lX11
## ------------------------------------------------------------------ ##
## Swinburne compiler (Intel fortran)
#FC = ifort
#F90C = ifort
#FFLAGS = -O -CB -std90 -check all
#F90FLAGS = -O -CB -std90 -check all
#LDFLAGS =  -lg2c -L/usr/local/pgplot -lpgplot -L/usr/X11R6/lib -lX11
## ------------------------------------------------------------------ ##

# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(F90C) -c $(F90FLAGS) $(FPPFLAGS) $< -o $@

DANSPH = read_data_dansph.f 
MRBSPH = read_data_mbate.f 
PGXTAL = read_data_dansph.f

SOURCES= modules.f90 prompting.f90 \
	 supersphplot.f90 main.f90 \
         calc_quantities.f90 \
	 colour_demo.f colour_set.f90 ../src/coord_transform.f90 \
	 danpgwedg.f \
	 defaults_read.f90 defaults_set.f90 defaults_write.f90 \
	 exact_rhoh.f90 \
	 exact_sedov.f exact_swave.f exact_toystar.f90 \
	 exact_toystar_ACplane.f exact_mhdshock.f90 \
	 integratedkernel.f90 \
	 interpolate1D.f90 \
         interpolate2D.f90 interpolate2D_xsec.f90 \
	 interpolate3D.f90 interpolate3D_fastxsec.f90 \
	 interpolate3D_projection.f90 \
	 int_from_string.f90 \
	 legend.f \
	 options.f90 options_exact.f90 \
	 options_powerspec.f90 options_render.f90 \
	 plot_average.f plot_kernel_gr.f90 \
	 plot_powerspectrum.f90 \
	 powerspectrum_fourier1D.f90 \
	 powerspectrum_lomb1D.f90 \
	 print_header.f90\
	 print_menu.f90 \
         render_coarse.f render.f90 \
	 transform.f90

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
