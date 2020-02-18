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
 use readwrite_fits,  only:read_fits_image,write_fits_image
 use imageutils,      only:image_denoise
 use iso_fortran_env, only:stderr=>error_unit, stdout=>output_unit
 implicit none
 character(len=120) :: file1,file2,fileout,tagline
 integer :: ierr, nargs, its, naxes(2), npixels
 real, allocatable :: image(:,:),hh(:),image1(:,:),image2(:,:),image_old(:,:),image_residuals(:,:)
 real :: fac,fluxold,fluxnew
 character(len=12) :: string
 logical :: iexist

 nargs = command_argument_count()

 tagline = 'splash-denoise: a SPLASH imaging utility (c) 2020 Daniel Price'
 if (nargs < 2) then
    print "(a)",trim(tagline)
    print "(/,a,/)",'Usage: splash-denoise infile.fits outfile.fits [options]'
    stop
 endif

 call get_command_argument(1,file1)
 if (index(file1,'.fits')==0) stop 'input file is not .fits'
 if (nargs >= 3) then
    call get_command_argument(2,file2)
    call get_command_argument(3,fileout)
 else
    call get_command_argument(2,fileout)
 endif
 if (index(fileout,'.fits')==0) stop 'output file is not .fits'

 ! check output file does not already exist
 inquire(file=fileout,exist=iexist)
 if (iexist) then
    write(stderr,*) 'ERROR: '//trim(fileout)//' already exists: please move or rename and try again'
    stop
 endif

 ! read original file
 call read_fits_image(file1,image1,naxes,ierr)
 if (ierr /= 0) stop 'error reading file'

 if (nargs >= 3) then
    call read_fits_image(file2,image2,naxes,ierr)
    if (ierr /= 0) stop 'error reading file'
    ! polarisation images
    write(stdout,"(a)") '>> Polarisation mode: sqrt(p**2 + q**2) using p='//trim(file1)//' q='//trim(file2)
    image = sqrt(image1**2 + image2**2)
 else
    image = image1
 endif

 image_old = image

 fac = 1.
 its = 4
 npixels = naxes(1)*naxes(2)
 allocate(hh(npixels))
 
 ! compute smoothing lengths and denoise total intensity
 call image_denoise(naxes,image,hh,iterations=its,fac=fac)

 ! for polarised images, denoise individual polarisations
 ! using h computed from the total intensity
 ! and stitch the denoised images back together
 if (nargs >= 3) then
    write(stdout,'(/,a)') '>> de-noising p...'
    call image_denoise(naxes,image1,hh,iterations=0)
    write(stdout,'(/,a)') '>> de-noising q...'
    call image_denoise(naxes,image2,hh,iterations=0)
    write(stdout,'(/,a)') '>> recombining to get total intensity...'
    image = sqrt(image1**2 + image2**2)
 endif

 fluxold = sum(image_old)
 fluxnew = sum(image)
 print "(a,g16.8,g16.8,a,1pg10.2,a)",'>> total intensity =',fluxold,fluxnew,&
                       ' err=',100.*abs(fluxold-fluxnew)/fluxold,' %'

 image_residuals = image - image_old
 call write_fits_image(fileout,image,naxes,ierr)
 call write_fits_image('old.fits',image_old,naxes,ierr)
 call write_fits_image('residuals.fits',image_residuals,naxes,ierr)

 if (ierr /= 0) write(stderr,*) 'error writing output file(s)'

 deallocate(image,hh,image_residuals)
 if (allocated(image1)) deallocate(image1)
 if (allocated(image2)) deallocate(image2)

end program denoise
