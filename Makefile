##-------------------------------------------------------------------##
##     Makefile for compiling supersphplot (uses PGPLOT)             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

## Compiler options (uncomment ONE set)
## ------------------------------------------------------------------ ##
## IoA compiler (Sun fortran)
FC = f90
F90C = f90
FFLAGS = -O4
F90FLAGS = -O4
LDFLAGS = -lpgplot -lX11 -lF77
## ------------------------------------------------------------------ ##
## NAGware f95 compiler
#FC = f95
#F90C = f95
#FFLAGS = -O4
#F90FLAGS = -O4
#LDFLAGS = -lpgplot -L/usr/X11R6/lib -lX11 -lg2c -lpng
## ------------------------------------------------------------------ ##
## Monash compiler (Lahey-Fujitsu f95)
#FC = lf95
#F90C = lf95
#FFLAGS = -O ##--chk aesux --chkglobal  --warn
#F90FLAGS = -O ##--chk aesux --chkglobal --warn
#LDFLAGS = -lpgplot -L/usr/X11R6/lib -lX11
## ------------------------------------------------------------------ ##
## Intel fortran compiler
#FC = ifc
#F90C = ifc
#FFLAGS = -O -C
#F90FLAGS = -O -C
#LDFLAGS = -lpgplot -L/usr/X11R6/lib -lX11 -lg2c -lpng -Vaxlib
## ------------------------------------------------------------------ ##

# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(F90C) -c $(F90FLAGS) $(FPPFLAGS) $< -o $@

DANSPH = read_data_dansph.f 
MRBSPH = read_data_mbate.f 

SOURCES= modules.f90 prompting.f90 \
	 supersphplot.f90 main.f90 \
         allocate.f90 calc_quantities.f90 \
	 colour_demo.f colour_set.f90 \
	 danpgsch.f danpgtile.f danpgwedg.f \
	 defaults_read.f90 defaults_set.f90 defaults_write.f90 \
	 exact_fromfile.f90 exact_rhoh.f90 \
	 exact_sedov.f90 exact_shock.f90 exact_wave.f90 \
	 exact_toystar.f90 exact_toystar2D.f90 \
	 exact_toystar_ACplane.f exact_mhdshock.f90 \
	 exact_polytrope.f \
	 get_data.f90 integratedkernel.f90 \
	 interactive_part.f90 \
	 interpolate1D.f90 \
         interpolate2D.f90 interpolate2D_xsec.f90 \
	 interpolate3D.f90 interpolate3D_fastxsec.f90 \
	 interpolate3D_projection.f90 \
	 int_from_string.f90 \
	 legend.f \
	 limits_read.f90 limits_save.f90 limits_set.f90 \
	 menu.f90 options_exact.f90 options_limits.f90 \
	 options_page.f90 options_particleplots.f90 \
	 options_powerspec.f90 options_render.f90 \
	 plot_average.f plot_kernel_gr.f90 \
	 plot_powerspectrum.f90 \
	 powerspectrum_fourier1D.f90 \
	 powerspectrum_lomb1D.f90 \
	 print_header.f90\
         read_exactparams.f90 riemannsolver.f90 \
         render.f90 setpage.f90 \
	 titles_read.f90 \
	 transform.f90 vectorplot.f90 \

SOURCESALL = $(SOURCES:.f90=.o)

OBJDANSPH = $(SOURCESALL:.f=.o) $(DANSPH:.f=.o) coord_transform.o
OBJMRBSPH = $(SOURCESALL:.f=.o) $(MRBSPH:.f=.o) coord_transform.o

dansph: $(OBJDANSPH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o ../supersphplot $(OBJDANSPH)

mrbsph: $(OBJMRBSPH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o supersphplot_mrb $(OBJMRBSPH)

## files from the main code

coord_transform.o:
	$(FC) -c $(FFLAGS) $(LDFLAGS) ../src/coord_transform.f90 -o $@

## other crap

tar:
	tar cf supersphplot.tar Makefile $(SOURCES) read_data*.f

targz:
	tar cf supersphplot.tar Makefile $(SOURCES) read_data*.f
	gzip supersphplot.tar

clean:
	rm *.o *.mod
