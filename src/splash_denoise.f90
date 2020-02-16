!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2020- Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
program denoise
 use readwrite_fits, only:read_fits_image,write_fits_image
 use imageutils,     only:image_denoise
 implicit none
 character(len=120) :: filein,fileout,tagline
 integer :: ierr, nargs, its, naxes(2), npixels
 real, allocatable :: image(:,:),hh(:)
 real :: fac

 nargs = command_argument_count()

 tagline = 'splash-denoise: a SPLASH imaging utility (c) 2020 Daniel Price'
 if (nargs < 2) then
    print "(a)",trim(tagline)
    print "(/,a,/)",'Usage: splash-denoise infile.fits outfile.fits [options]'
    stop
 endif

 call get_command_argument(1,filein)
 call get_command_argument(2,fileout)
 if (index(filein,'.fits')==0) stop 'input file is not .fits'
 if (index(fileout,'.fits')==0) stop 'output file is not .fits'

 ! read original file
 call read_fits_image(filein,image,naxes,ierr)
 if (ierr /= 0) stop 'error reading file'

 fac = 1.
 its = 4
 npixels = naxes(1)*naxes(2)
 allocate(hh(npixels))
 call image_denoise(naxes,image,hh,iterations=its,fac=fac)

 call write_fits_image(fileout,image,naxes,ierr)
 if (ierr /= 0) stop 'error writing output file'

 deallocate(image,hh)

end program denoise
