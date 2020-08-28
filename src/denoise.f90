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
 use readwrite_fits,  only:read_fits_cube,write_fits_cube,write_fits_image,get_from_header
 use imageutils,      only:image_denoise,image_denoise3D
 use iso_fortran_env, only:stderr=>error_unit, stdout=>output_unit
 use system_utils,    only:get_command_option,count_matching_args,get_command_flag
 implicit none
 character(len=120) :: file1,file2,fileout,tagline,filek
 integer :: ierr,jext,nfiles,its,naxes(4),npixels,k,i,iarglist(3),kstart,kend,k1,k2
 real, allocatable :: hh(:)
 real, allocatable :: image(:,:,:),image1(:,:,:),image2(:,:,:)
 real, allocatable :: image_old(:,:,:),image_residuals(:,:,:)
 character(len=:), allocatable :: fitsheader(:)
 real :: beam,fluxold,fluxnew,imax,imagemax,image_max,use3D,bmaj,bmin,cdelt1
 logical :: iexist,use_3D,skip

 nfiles = count_matching_args('.fits',iarglist)

 tagline = 'denoise: a SPLASH imaging utility (c) 2020 Daniel Price'
 if (nfiles < 1) then
    print "(a)",trim(tagline)
    print "(/,a)",'Usage: denoise [options] infile.fits [outfile.fits]'
    print "(/,a)",'Options:  --imax=3.4e-2     [intensity value above which no smoothing is applied]'
    print "(a)",  '          --beam=1.0        [beam size in pixels at max intensity]'
    print "(a)",  '          --its=4           [maximum number of smoothing length iterations]'
    print "(a)",  '          --use3D=1         [denoise in 3D for spectral cubes]'
    print "(a)",  '          --start=1         [denoise from channel 1 onwards]'
    print "(a,/)",'          --end=10          [denoise only up to channel 10]'
    print "(a)",  '          --nosubsample     [do not subsample if beam size >> pixel scale]'
    stop
 endif

 !
 ! read filenames from the command line
 !
 call get_command_argument(iarglist(1),file1)
 if (nfiles >= 3) then
    call get_command_argument(iarglist(2),file2)
    call get_command_argument(iarglist(3),fileout)
 elseif (nfiles >= 2) then
    call get_command_argument(iarglist(2),fileout)
 else
    fileout = 'denoise.fits' ! if no .fits extension found in original file
    jext = index(file1,'.fits')
    fileout = file1(1:jext-1)//'_denoise.fits' ! default if no output filename given
 endif
 !
 ! get options from the command line
 !
 imax = get_command_option('imax',default=0.)
 beam   = get_command_option('beam',default=0.)
 use3D  = get_command_option('use3D',default=0.)
 kstart = nint(get_command_option('start',default=0.))
 kend   = nint(get_command_option('end',default=0.))
 its    = nint(get_command_option('its',default=4.))
 skip   = .not.get_command_flag('nosubsample')

 !print*,' GOT imax=',imax,' fac =  ',fac,trim(file1),trim(file2),trim(fileout)

 ! check output file does not already exist
 inquire(file=fileout,exist=iexist)
 if (iexist) then
    write(stderr,*) 'ERROR: '//trim(fileout)//' already exists: please move or rename and try again'
    stop
 endif

 ! read original file
 call read_fits_cube(file1,image1,naxes,ierr,hdr=fitsheader)
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
 where (isnan(image))
    image = 0.
 end where

 image_old = image

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
 !
 ! get beamsize from fits header, in units of pixels
 !
 if (beam <= tiny(beam)) then
    bmaj = get_from_header('BMAJ',fitsheader,ierr)
    bmin = get_from_header('BMIN',fitsheader,ierr)
    cdelt1 = get_from_header('CDELT1',fitsheader,ierr)
!    beam = sqrt(bmaj*bmin)/abs(cdelt1)
    if (abs(cdelt1) > 0. .and. ierr == 0) then
       beam = bmaj/abs(cdelt1) !sqrt(bmin*bmaj)/abs(cdelt1)
       print*,' bmin = ',bmin/abs(cdelt1),' bmaj = ',bmaj/abs(cdelt1), ' mean = ',sqrt(bmin*bmaj)/abs(cdelt1)
    endif
    if (beam < 1. .or. beam > 20.) beam = 1.
 endif
 print "(3(a,1pg16.4))",'>> max intensity = ',imagemax,' of ',image_max,', min pix per beam = ',beam
 !imagemax = beam**2*imagemax

 if (use_3D) then
    !
    ! denoise all channels simultaneously in 3D
    !
    call image_denoise3D(naxes,image,hh,iterations=its,imax=imagemax,beam=beam)
 else
    !
    ! denoise channel by channel in 2D
    !
    k1 = 1
    k2 = naxes(3)
    if (kstart > 0) k1 = kstart
    if (kend > 0 .and. kend <= naxes(3)) k2 = kend
    do k=k1,k2
       if (naxes(3) > 1) print "(a,i5,a,i5)",'>> channel ',k, ' of ',naxes(3)
       call image_denoise(naxes(1:2),image(:,:,k),hh,iterations=its,imax=imagemax,beam=beam,skip=skip)
       !print*,'min, max,mean h = ',minval(hh),maxval(hh),sum(hh)/real(npixels)

       ! for polarised images, denoise individual polarisations
       ! using h computed from the total intensity
       ! and stitch the denoised images back together
       if (nfiles >= 3) then
          write(stdout,'(/,a)') '>> de-noising p...'
          call image_denoise(naxes,image1,hh,iterations=0,skip=skip)
          write(stdout,'(/,a)') '>> de-noising q...'
          call image_denoise(naxes,image2,hh,iterations=0,skip=skip)
          write(stdout,'(/,a)') '>> recombining to get total intensity...'
          image = sqrt(image1**2 + image2**2)
       endif

       fluxold = sum(image_old(:,:,k))
       fluxnew = sum(image(:,:,k))
       print "(a,g16.8,a,g16.8,a,1pg10.2,a)",'>> total intensity =',fluxnew,' was ',fluxold,&
                       ' err=',100.*(fluxnew-fluxold)/abs(fluxold),' %'
       if (naxes(3) > 1 .and. k2==k1) then
          i = index(fileout,'.fits')
          write(filek,*) k
          filek = fileout(1:i-1)//'-c'//trim(adjustl(filek))//'.fits'
          call write_fits_image(filek,image(:,:,k),naxes(1:2),ierr,hdr=fitsheader)
       endif
    enddo
 endif

 call write_fits_cube(fileout,image,naxes,ierr,hdr=fitsheader)
 image_residuals = image - image_old
 !call write_fits_cube('old.fits',image_old,naxes,ierr)
 !call write_fits_cube('residuals.fits',image_residuals,naxes,ierr)

 if (ierr /= 0) write(stderr,*) 'error writing output file(s)'

 deallocate(image,hh)
 if (allocated(image_residuals)) deallocate(image_residuals)
 if (allocated(image1)) deallocate(image1)
 if (allocated(image2)) deallocate(image2)

end program denoise
