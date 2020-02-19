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

 public :: image_denoise, image_denoise3D
 private

contains

!----------------------------------------------------------------------
!
!  denoise an image using adaptive kernel convolution
!  i.e. variable smoothing lengths, where the smoothing scale
!  is inversely proportional to the image intensity
!
!----------------------------------------------------------------------
subroutine image_denoise(naxes,image,hh,iterations,imax,fac,err)
 use interpolations2D, only:interpolate2D
 use kernels,          only:select_kernel
 integer, intent(in)    :: naxes(2)
 real,    intent(inout) :: image(naxes(1),naxes(2))
 real,    intent(inout) :: hh(naxes(1)*naxes(2))
 integer, intent(in), optional :: iterations
 real,    intent(in), optional :: fac,imax
 integer, intent(out), optional :: err
 real, allocatable      :: x(:),y(:),weight(:),dat(:)
 integer, allocatable   :: mask(:)
 integer :: npixels,its,niter,i,j,n,ierr
 real :: imagemax,dx,dy,scalefac,xmin,ymin

 call select_kernel(0)

 ! allocate memory for temporary arrays
 npixels = product(naxes)
 allocate(x(npixels),y(npixels),weight(npixels),dat(npixels),mask(npixels),stat=ierr)
 if (present(err)) err = ierr
 if (ierr /= 0) return

 dy = 1.
 dx = 1.
 xmin = 0.0
 ymin = 0.0
 mask(:) = 1 ! do not mask any pixels
 if (present(imax)) then
    imagemax = imax
    print*,' imagemax = ',imagemax,' was ',maxval(image)
 else
    imagemax = maxval(image)
 endif
 !fluxold  = sum(image)
 !print*,' total intensity =',fluxold

 niter = 4
 if (present(iterations)) niter = iterations

 scalefac = 1.0
 if (present(fac)) scalefac = fac

 ! find the convolution length by iteration
 h_iterations: do its=1,max(niter,1)
    if (niter > 0) print "(' Iteration: ',i2,' of ',i2)",its,niter
    n = 0
    do j=1,naxes(2)
       do i=1,naxes(1)
          n = n + 1
          x(n) = i - 0.5*dx
          y(n) = j - 0.5*dy
          if (niter > 0) then
             if (abs(image(i,j)) > 0.) then
                hh(n) = min(max(0.4*sqrt(scalefac*imagemax/abs(image(i,j))),1.),5.*(its+1))
             else
                hh(n) = 5.*(its+1)
             endif
          endif
          if (its==1) dat(n) = image(i,j)
       enddo
    enddo

    weight(1:npixels) = dx*dy/hh(1:npixels)**2
    call interpolate2D(x,y,hh,weight,dat,mask,npixels, &
         xmin,ymin,image,naxes(1),naxes(2),dx,dy,&
         normalise=.false.,exact=.true.,periodicx=.false.,periodicy=.false.)

 enddo h_iterations

 ! deallocate memory
 deallocate(x,y,mask,weight,dat)

end subroutine image_denoise

!----------------------------------------------------------------------
!
!  denoise an image using adaptive kernel convolution in 3D
!
!----------------------------------------------------------------------
subroutine image_denoise3D(naxes,image,hh,iterations,imax,err)
 use interpolations3D, only:interpolate3D
 use kernels,          only:select_kernel
 integer, intent(in)    :: naxes(3)
 real,    intent(inout) :: image(naxes(1),naxes(2),naxes(3))
 real,    intent(inout) :: hh(naxes(1)*naxes(2)*naxes(3))
 integer, intent(in), optional :: iterations
 real,    intent(in), optional :: imax
 integer, intent(out), optional :: err
 real, allocatable      :: x(:),y(:),z(:),weight(:),dat(:)
 integer, allocatable   :: mask(:)
 integer :: npixels,its,niter,i,j,k,n,ierr
 real :: imagemax,dx,dy,dz,xmin,ymin,zmin
 real(8), allocatable :: image_dp(:,:,:)

 call select_kernel(0)

 ! allocate memory for temporary arrays
 npixels = product(naxes)
 allocate(x(npixels),y(npixels),z(npixels),weight(npixels),dat(npixels),mask(npixels),stat=ierr)
 if (present(err)) err = ierr
 if (ierr /= 0) return

 if (size(hh) /= size(x)) then
    print*,' error: smoothing length array too small'
    ierr = 3
    return
 endif

 dz = 1.
 dy = 1.
 dx = 1.
 xmin = 0.0
 ymin = 0.0
 zmin = 0.0
 mask(:) = 1 ! do not mask any pixels
 if (present(imax)) then
    imagemax = imax
 else
    imagemax = maxval(image)
 endif
 
 niter = 4
 if (present(iterations)) niter = iterations

 ! find the convolution length by iteration
 h_iterations: do its=1,max(niter,1)
    if (niter > 0) print "(' Iteration: ',i2,' of ',i2)",its,niter
    n = 0

    do k=1,naxes(3)
       do j=1,naxes(2)
          do i=1,naxes(1)
             n = n + 1
             !print*,n,size(z),size(hh),size(dat),i,j,k
             x(n) = i - 0.5*dx
             y(n) = j - 0.5*dy
             z(n) = k - 0.5*dz
             if (niter > 0) then
                if (abs(image(i,j,k)) > 0.) then
                   hh(n) = min(max(0.4*sqrt(imagemax/abs(image(i,j,k))),1.),5.*(its+1))
                else
                   hh(n) = 0.
                endif
             endif
             if (its==1) dat(n) = image(i,j,k)
          enddo
       enddo
    enddo
    allocate(image_dp(naxes(1),naxes(2),naxes(3)))
    weight(1:npixels) = dx*dy*dz/hh(1:npixels)**3
    call interpolate3D(x,y,z,hh,weight,dat,mask,npixels, &
         xmin,ymin,zmin,image_dp,naxes(1),naxes(2),naxes(3),dx,dy,dz,&
         normalise=.false.,periodicx=.false.,periodicy=.false.,periodicz=.false.)

    image = real(image_dp)
 enddo h_iterations

 ! deallocate memory
 deallocate(x,y,z,mask,weight,dat)
 if (allocated(image_dp)) deallocate(image_dp)

end subroutine image_denoise3D

end module imageutils
