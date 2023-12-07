!---------------------------------------------------------------------------------
!
!     sphmoments - a utility for computing kernel-smoothed moment maps 
!                  of observational data cubes
!
!     Copyright (C) 2023 Daniel Price
!     daniel.price@monash.edu
!
!     --------------------------------------------------------------------------
!
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!     -------------------------------------------------------------------------
!     Version history/ Changelog:
!     0.9.7   : (06/11/23)
!             fix velocity array read if comments in the fits header;
!             bug fix with moment calculation if velocity array in reverse order
!     0.9.1   : (09/10/23)
!             fix homebrew install
!     0.9.0   : (09/10/23)
!             velocity array read correctly from fits header
!     0.8.0   : (09/10/23)
!             write correct x and y units to output files;
!             copy and flatten fits header from original cube
!     0.7.0   : (06/10/23)
!             basic working version implemented
!
!     -------------------------------------------------------------------------
!
!     Written by DJP initially as part of the exoALMA project
!     
!---------------------------------------------------------------------------------
program sph_moments
 use readwrite_fits,  only:read_fits_cube,write_fits_image,get_from_header,flatten_header
 use iso_fortran_env, only:stderr=>error_unit, stdout=>output_unit
 use system_utils,    only:get_command_option,count_matching_args,ienvlist
 use moments,         only:get_moments
 use kernels,         only:select_kernel,kernelname
 !$ use omp_lib, only:omp_get_num_threads
 implicit none
 character(len=256) :: file1,fileout,fileout_m(5),tagline
 integer :: ierr,jext,nfiles,naxes(4),k,iarglist(2),nm,ns,ikernel
 real, allocatable :: cube(:,:,:),moment(:,:,:),zvals(:)
 character(len=80), allocatable :: fitsheader(:)
 character(len=12) :: string
 real :: dv,hfac
 logical :: iexist
 integer :: moment_type(5)

 nfiles = count_matching_args('.fits',iarglist)

 tagline = 'sph_moments v0.9.7: a SPLASH imaging utility (c) 2023 Daniel Price'
 if (nfiles < 1) then
    print "(a)",trim(tagline)
    print "(/,a)",'Usage: sph_moments [options] infile.fits [outfile.fits]'
    print "(/,a)",'Options:  --moments=0,1,9   [which moments to take, default=0,1,2,8,9]'
    print "(a)",  '          --kernel=0        [which smoothing kernel to use: 0=cubic spline 2=quartic 3=quintic]'
    print "(a)",  '          --hfac=3.         [ratio of smoothing length to channel spacing in spectral dimension]'
    print "(a)",  '          --sample=10       [factor by which to oversample line profiles, not very important]'
    stop
 endif

 !
 ! read filenames from the command line
 !
 call get_command_argument(iarglist(1),file1)
 if (nfiles >= 2) then
    call get_command_argument(iarglist(2),fileout)
 else
    fileout = 'moment.fits' ! if no .fits extension found in original file
    jext = index(file1,'.fits')
    fileout = file1(1:jext-1)//'_moment.fits' ! default if no output filename given
 endif
 !
 ! get options from the command line
 !
 ns = nint(get_command_option('sample'))
 hfac = get_command_option('hfac',default=3.)
 ikernel = nint(get_command_option('kernel',default=0.))
 call select_kernel(ikernel)

 !
 ! list of moments to compute
 !
 moment_type = ienvlist('moments',size(moment_type),errval=-1)
 nm = 0
 do k=1,size(moment_type)
    if (k > 1) then
       if (moment_type(k)==moment_type(k-1)) exit
    endif
    if (moment_type(k) >= 0 .and. moment_type(k) <= 9) then
       nm = nm + 1
    endif
 enddo
 ! default behaviour if no moments specified is to compute all
 if (nm==0) then
    nm = 5
    moment_type = (/0,1,2,8,9/)
 endif

 ! check output files do not already exist
 iexist = .false.
 do k=1,nm
    fileout_m(k) = fileout
    jext = index(fileout_m(k),'.fits')
    write(string,"(i0)") moment_type(k)
    fileout_m(k) = fileout_m(k)(1:jext-1)//trim(string)//'.fits'
    inquire(file=fileout_m(k),exist=iexist)
    if (iexist) then
       write(stderr,*) 'ERROR: '//trim(fileout_m(k))//' already exists: please move or rename and try again'
    endif
 enddo
 if (iexist) stop

 ! read original file
 write(stdout,*) '>> reading '//trim(file1)
 call read_fits_cube(file1,cube,naxes,ierr,hdr=fitsheader,velocity=zvals)
 if (ierr /= 0) stop 'error reading file'

 ! eliminate NaNs
 where (isnan(cube))
    cube = 0.
 end where

 if (naxes(3) <= 1) stop 'this is a 2D image but require 3D fits cube for moment generation'
 
 !
 ! find moments using either all channels or range of channels
 !
 write(stdout,"(1x,a,f5.2)") '>> computing moments using '//trim(kernelname(ikernel))//' kernel with h/dz = ',hfac
 ! warn about running in serial
 !$omp parallel
 !$ if (omp_get_num_threads()<=1) write(stderr,"(/,1x,a)") &
 !$    'WARNING: running in serial: type "export OMP_NUM_THREADS=8" to use 8 threads'
 !$omp end parallel
 call get_moments(cube,zvals,hfac,moment_type(1:nm),moment)
 !
 ! write results to a fits file
 !
 call flatten_header(fitsheader)
 do k=1,nm
    call write_fits_image(fileout_m(k),moment(:,:,k),naxes(1:2),ierr,hdr=fitsheader)
 enddo

 deallocate(cube,moment,zvals)

end program sph_moments
