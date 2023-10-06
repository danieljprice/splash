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
!     0.1.0   : (06/10/23)
!             basic working version implemented
!
!     -------------------------------------------------------------------------
!
!     Written by DJP initially as part of the exoALMA project
!     
!---------------------------------------------------------------------------------
program sph_moments
 use readwrite_fits,  only:read_fits_cube,write_fits_image,get_from_header
 use iso_fortran_env, only:stderr=>error_unit, stdout=>output_unit
 use system_utils,    only:get_command_option,count_matching_args
 use moments,         only:get_moments
 use kernels,         only:select_kernel
 !use timing,          only:timer_start,timer_stop
 implicit none
 character(len=240) :: file1,fileout,fileout1,tagline
 integer :: ierr,jext,nfiles,naxes(4),k,iarglist(2),k1,k2,m,ns,ikernel
 real, allocatable :: cube(:,:,:),moment(:,:,:),zvals(:)
 character(len=80), allocatable :: fitsheader(:)
 character(len=12) :: string
 real :: dv,hfac
 logical :: iexist
 integer, parameter :: moment_type(5) = (/0,1,2,8,9/)

 nfiles = count_matching_args('.fits',iarglist)

 tagline = 'sph_moments: a SPLASH imaging utility (c) 2023 Daniel Price'
 if (nfiles < 1) then
    print "(a)",trim(tagline)
    print "(/,a)",'Usage: sph_moments [options] infile.fits [outfile.fits]'
    print "(/,a)",'Options:  --moments=0,1,9   [which moments to take, default=0,1,2,8,9]'
    print "(a)",  '          --kernel=0        [which smoothing kernel to use: 0=cubic spline 2=quartic 3=quintic]'
    print "(a)",  '          --hfac=5.         [ratio of smoothing length to channel spacing in spectral dimension]'
    print "(a)",  '          --sample=10       [factor by which to oversample line profiles]'
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
 m = nint(get_command_option('moments',default=0.))
 ns = nint(get_command_option('sample'))
 hfac = get_command_option('hfac',default=5.)
 ikernel = nint(get_command_option('kernel',default=0.))

 call select_kernel(ikernel)

 ! check output file does not already exist
 inquire(file=fileout,exist=iexist)
 if (iexist) then
    write(stderr,*) 'ERROR: '//trim(fileout)//' already exists: please move or rename and try again'
    stop
 endif

 ! read original file
 call read_fits_cube(file1,cube,naxes,ierr,hdr=fitsheader)
 if (ierr /= 0) stop 'error reading file'

 !print*,' got fitsheader = ',fitsheader

 ! eliminate NaNs
 where (isnan(cube))
    cube = 0.
 end where

 if (naxes(3) <= 1) stop 'this is a 2D image but require 3D fits cube for moment generation'

 dv = get_from_header('CDELT3',fitsheader,ierr)
 print*,' got dv = ',dv
 allocate(zvals(naxes(3)))
 do k=1,naxes(3)
    zvals(k) = 0.5 + (k-1)
 enddo
 !
 ! find moments using either all channels or range of channels
 !
 k1 = 1
 k2 = naxes(3)
 write(stdout,*) '>> computing SPH kernel moments ',moment_type,'...'
 call get_moments(cube,zvals,hfac,moment_type,moment)
 !
 ! write results to a fits file
 !

 do k=1,size(moment,dim=3)
    fileout1 = fileout
    jext = index(fileout1,'.fits')
    write(string,"(i0)") moment_type(k)
    fileout1 = fileout1(1:jext-1)//trim(string)//'.fits'
    call write_fits_image(fileout1,moment(:,:,k),naxes(1:2),ierr)
 enddo

 deallocate(cube,moment)

end program sph_moments
