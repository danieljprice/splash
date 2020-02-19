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
 use readwrite_fits,  only:read_fits_cube,write_fits_cube,write_fits_image
 use imageutils,      only:image_denoise,image_denoise3D
 use iso_fortran_env, only:stderr=>error_unit, stdout=>output_unit
 use system_utils,    only:get_command_option,count_matching_args
 implicit none
 character(len=120) :: file1,file2,fileout,tagline
 integer :: ierr, nfiles, its, naxes(3), npixels, k, iarglist(3)
 real, allocatable :: hh(:)
 real, allocatable :: image(:,:,:),image1(:,:,:),image2(:,:,:)
 real, allocatable :: image_old(:,:,:),image_residuals(:,:,:)
 real :: fac,fluxold,fluxnew,imax,imagemax,image_max,use3D
 character(len=12) :: string
 logical :: iexist,use_3D

 nfiles = count_matching_args('.fits',iarglist)

 tagline = 'splash-denoise: a SPLASH imaging utility (c) 2020 Daniel Price'
 if (nfiles < 2) then
    print "(a)",trim(tagline)
    print "(/,a)",'Usage: splash-denoise [options] infile.fits outfile.fits'
    print "(/,a)",'Options:  --imax=3.4e-2     [intensity value above which no smoothing is applied]'
    print "(a)",  '          --fac=1.0         [scaling factor to increase/decrease smoothing scale]'
    print "(a,/)",'          --use3D=1         [denoise in 3D for spectral cubes]'
    stop
 endif

 call get_command_argument(iarglist(1),file1)
 if (nfiles >= 3) then
    call get_command_argument(iarglist(2),file2)
    call get_command_argument(iarglist(3),fileout)
 else
    call get_command_argument(iarglist(2),fileout)
 endif
 imax = get_command_option('imax',default=0.)
 fac  = get_command_option('fac',default=1.)
 use3D  = get_command_option('use3D',default=0.)

 !print*,' GOT imax=',imax,' fac =  ',fac,trim(file1),trim(file2),trim(fileout)

 ! check output file does not already exist
 inquire(file=fileout,exist=iexist)
 if (iexist) then
    write(stderr,*) 'ERROR: '//trim(fileout)//' already exists: please move or rename and try again'
    stop
 endif

 ! read original file
 call read_fits_cube(file1,image1,naxes,ierr)
 if (ierr /= 0) stop 'error reading file'

 if (nfiles >= 3) then
    call read_fits_cube(file2,image2,naxes,ierr)
    if (ierr /= 0) stop 'error reading file'
    ! polarisation images
    write(stdout,"(a)") '>> Polarisation mode: sqrt(p**2 + q**2) using p='//trim(file1)//' q='//trim(file2)
    image = sqrt(image1**2 + image2**2)
 else
    image = image1
 endif

 ! eliminate NaNs
 where (image /= image)
    image = 0.
 end where

 image_old = image

 its = 4
 use_3D = (naxes(3) > 1) .and. (use3D > 0.)
 if (use_3D) then
    npixels = naxes(1)*naxes(2)*naxes(3)
 else
    npixels = naxes(1)*naxes(2)
 endif
 allocate(hh(npixels))

 image_max = maxval(image)
 if (imax > 0.) then
    imagemax = imax
 else
    imagemax = image_max
 endif
 print "(3(a,g16.8))",'>> max intensity = ',imagemax,' of ',image_max,', using scaling factor ',fac
 imagemax = fac*imagemax

 if (use_3D) then
    !
    ! denoise all channels simultaneously in 3D
    !
    call image_denoise3D(naxes,image,hh,iterations=its,imax=imagemax)
 else
    !
    ! denoise channel by channel in 2D
    !
    do k=1,naxes(3)
       if (naxes(3) > 1) print "(a,i3)",'>> channel ',k
       call image_denoise(naxes(1:2),image(:,:,k),hh,iterations=its,imax=imagemax)
       print*,'min, max,mean h = ',minval(hh),maxval(hh),sum(hh)/real(npixels)

       ! for polarised images, denoise individual polarisations
       ! using h computed from the total intensity
       ! and stitch the denoised images back together
       if (nfiles >= 3) then
          write(stdout,'(/,a)') '>> de-noising p...'
          call image_denoise(naxes,image1,hh,iterations=0)
          write(stdout,'(/,a)') '>> de-noising q...'
          call image_denoise(naxes,image2,hh,iterations=0)
          write(stdout,'(/,a)') '>> recombining to get total intensity...'
          image = sqrt(image1**2 + image2**2)
       endif

       fluxold = sum(image_old(:,:,k))
       fluxnew = sum(image(:,:,k))
       print "(a,g16.8,g16.8,a,1pg10.2,a)",'>> total intensity =',fluxold,fluxnew,&
                       ' err=',100.*abs(fluxold-fluxnew)/fluxold,' %'
    enddo
 endif

 !call write_fits_image('slice.fits',image(:,:,70),naxes(1:2),ierr)
 !call write_fits_image('sliceold.fits',image_old(:,:,70),naxes(1:2),ierr)

 image_residuals = image - image_old
 call write_fits_cube(fileout,image,naxes,ierr)
 !call write_fits_cube('old.fits',image_old,naxes,ierr)
 call write_fits_cube('residuals.fits',image_residuals,naxes,ierr)

 if (ierr /= 0) write(stderr,*) 'error writing output file(s)'

 deallocate(image,hh)
 if (allocated(image_residuals)) deallocate(image_residuals)
 if (allocated(image1)) deallocate(image1)
 if (allocated(image2)) deallocate(image2)

end program denoise
