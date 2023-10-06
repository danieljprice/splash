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
!  Copyright (C) 2023- Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
module moments
 implicit none

contains

subroutine get_moments(cube,zvals,hfac,moment_type,moments)
 use interpolations1D, only:interpolate1D
 real,    intent(in)  :: cube(:,:,:),zvals(:),hfac
 integer, intent(in)  :: moment_type(:)
 real,    intent(out), allocatable :: moments(:,:,:)
 real, allocatable :: profile(:),profile_smooth(:),hh(:),weight(:),zvals_smooth(:)
 integer, allocatable :: mask(:)
 real :: dz,dz_smooth,zmin
 integer :: i,j,k,m,nx,ny,nz,nm,nz_smooth,iu,jcount,ktmp(1),iprint,jprint

 ! compute moment
 nx = size(cube,1)
 ny = size(cube,2)
 nz = size(cube,3)
 nm = size(moment_type)
 nz_smooth = 10*nz
 allocate(profile(nz),profile_smooth(nz_smooth),zvals_smooth(nz_smooth),hh(nz),weight(nz),mask(nz))
 allocate(moments(nx,ny,nm))

 dz = zvals(2)-zvals(1)
 dz_smooth = (dz*nz)/nz_smooth
 zmin = zvals(1) - 0.5*dz
 do k=1,nz_smooth
    zvals_smooth(k) = zmin + (k-0.5)*dz_smooth
 enddo
 mask = 1
 hh = hfac*dz
 weight = dz/hh

 jprint = ny/2
 iprint = nx/2 + 10

 write(*,"(/,a,45x,a)") ' 0% ',' 100%'
 write(*,"(1x,a)",advance='no') '|'
 jcount = 0
 !$omp parallel do default(none) shared(cube,moments,moment_type,ny,nx,zvals,hh,weight,mask) &
 !$omp shared(nz,nz_smooth,nm,zmin,dz_smooth,zvals_smooth,jcount,iu,iprint,jprint) &
 !$omp private(j,i,k,m,ktmp,profile,profile_smooth)
 do j=1,ny
    !$omp atomic
    jcount = jcount + 1
    if (mod(jcount,ny/50)==0) write(*,"('=')",advance='no')
    do i=1,nx
       profile = cube(i,j,:)
       call interpolate1D(zvals,hh,weight,profile,mask,nz,zmin,profile_smooth,nz_smooth,dz_smooth,normalise=.false.,iverbose=0)
       do m=1,nm
          select case(moment_type(m))
          case(9) ! peak velocity
             ktmp = maxloc(profile_smooth)
             moments(i,j,m) = zvals_smooth(ktmp(1))
          case(8) ! peak intensity
             moments(i,j,m) = maxval(profile_smooth)
          case(1:2)
             moments(i,j,m) = sum(profile_smooth*zvals_smooth**moment_type(m))
          case default ! moment 0
             moments(i,j,m) = sum(profile_smooth)
          end select
       enddo
       if (j==jprint .and. i==iprint) then
          open(newunit=iu,file='profile.dat',status='replace')
          do k=1,nz
             write(iu,*) zvals(k),profile(k)
          enddo
          close(iu)

          open(newunit=iu,file='profile_smooth.dat',status='replace')
          do k=1,nz_smooth
             write(iu,*) zvals_smooth(k),profile_smooth(k)
          enddo
          close(iu)
       endif
    enddo
 enddo
 !$omp end parallel do
 write(*,"(a,/)") '|'

end subroutine get_moments

end module moments
