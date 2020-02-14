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
!----------------------------------------------------------------------
!
!  Some utilities for image manipulation using SPH algorithms
!
!----------------------------------------------------------------------
module imageutils
 implicit none

 public :: image_denoise
 private

contains

!----------------------------------------------------------------------
!
!  denoise an image using adaptive kernel convolution
!  i.e. variable smoothing lengths, where the smoothing scale
!  is inversely proportional to the image intensity
!
!----------------------------------------------------------------------
subroutine image_denoise(naxes,image,hh,iterations,fac,err)
 use interpolations2D, only:interpolate2D
 integer, intent(in)    :: naxes(2)
 real,    intent(inout) :: image(naxes(1),naxes(2))
 real,    intent(out)   :: hh(naxes(1)*naxes(2))
 integer, intent(in), optional :: iterations
 real,    intent(in), optional :: fac
 integer, intent(out), optional :: err
 real, allocatable      :: x(:),y(:),weight(:),dat(:)
 integer, allocatable   :: mask(:)
 integer :: npixels,its,niter,i,j,n,ierr
 real :: fluxold,fluxnew,imagemax,dx,dy,scalefac,xmin,ymin

 ! allocate memory for temporary arrays
 npixels = product(naxes)
 allocate(x(npixels),y(npixels),weight(npixels),dat(npixels),mask(npixels),stat=ierr)
 if (present(err)) err = ierr
 if (ierr /= 0) return

 xmin = 0.
 ymin = 0.
 dy = 1.
 dx = 1.
 mask(:) = 1 ! do not mask any pixels
 imagemax = maxval(image)
 fluxold  = sum(image)
 print*,' total intensity =',fluxold

 niter = 4
 if (present(iterations)) niter = iterations

 scalefac = 1.0
 if (present(fac)) scalefac = fac

 ! find the convolution length by iteration
 h_iterations: do its=1,niter
    print "(' Iteration: ',i2,' of ',i2)",its,niter
    n = 0
    do j=1,naxes(2)
       do i=1,naxes(1)
          n = n + 1
          x(n) = i
          y(n) = j
          if (abs(image(i,j)) > 0.) then
             hh(n) = min(max(0.4*sqrt(scalefac*imagemax/abs(image(i,j))),1.),5.*(its+1))
          else
             hh(n) = 0.
          endif
          if (its==1) dat(n) = image(i,j)
       enddo
    enddo

    weight(1:npixels) = dx*dy/hh(1:npixels)**2
    call interpolate2D(x,y,hh,weight,dat,mask,npixels, &
         xmin,ymin,image,naxes(1),naxes(2),dx,dy,&
         normalise=.false.,exact=.true.,periodicx=.false.,periodicy=.false.)

    fluxnew = sum(image)
    print "(a,g16.8,a,1pg10.2,a)",' total intensity =',fluxnew,&
                          ' err=',100.*abs(fluxold-fluxnew)/fluxold,' %'

 enddo h_iterations

 ! deallocate memory
 deallocate(x,y,mask,weight,dat)

end subroutine image_denoise

end module imageutils
