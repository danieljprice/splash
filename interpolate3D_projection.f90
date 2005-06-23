!----------------------------------------------------------------------
!
!  Module containing all of the routines required for 3D projections
!  (where rendered quantity is integrated along the line of sight)
!
!----------------------------------------------------------------------

module projections3D
 implicit none

 integer, parameter, private :: maxcoltable = 1000
 real, parameter, private :: pi = 3.1415926536      
 real, dimension(maxcoltable), private :: coltable

 public :: setup_integratedkernel
 public :: interpolate3D_projection
 public :: interpolate3D_proj_vec

contains

subroutine setup_integratedkernel
!-------------------------------------------------------------
!     tabulates the integral through the cubic spline kernel
!     subroutine originally by Matthew Bate
!-------------------------------------------------------------
 implicit none
 integer :: i,j
 real :: r, dist, step, ypos, v, v2, val
 real :: coldens, v2m

 print "(a)",' setting up integrated kernel table...'

 do i=1,maxcoltable
    r=(i-1)/500.
    dist=sqrt(4.0-r*r)
    step=dist/4000.0
    ypos=0.
         
    coldens=0.0
    do j=1,4000
       v=sqrt(r*r+ypos*ypos)
       if (v.lt.1.0) then
          v2=v*v
          val=1.0-1.5*v2+0.75*v2*v
          coldens=coldens+val*step
       else
          v2m=2.0-v
          val=0.25*v2m*v2m*v2m
          coldens=coldens+val*step        
       endif
       ypos=ypos+step
    end do
    coltable(i)=2.0*coldens/pi
 end do

!      check=0.0
!      do i=1,1000
!         r=(i-1)/500.
!         write(*,*) r, table(i)
!         check=check+2*pi*r*table(i)*(1.0/500.)
!      end do
!      write(*,*)check

 return
end subroutine setup_integratedkernel

!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_b m_b dat_b/rho_b W(r-r_b, h_b)
! 
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline
!
!     ** In this version 3D data is interpolated to a 2D grid by use of an
!     ** integrated form of the kernel (that is W_ab in this case is
!     ** the integral through the 3D kernel to give a 2D kernel)
!     ** This results in a column density map of the interpolated quantity
!     ** From a similar routine by Matthew Bate.
!
!     Input: particle coordinates  : x,y   (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price September 2003
!--------------------------------------------------------------------------

subroutine interpolate3D_projection(x,y,pmass,rho,hh,dat,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidth)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,pmass,rho,hh,dat
  real, intent(in) :: xmin,ymin,pixwidth
  real, intent(out), dimension(npixx,npixy) :: datsmooth

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: index, index1
  integer :: iprintinterval, iprintnext, iprogress, itmin
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: term,rho1i,dx,dy,xpix,ypix
  real :: dxx,dwdx
  real :: dmaxcoltable
  real :: t_start,t_end,t_used,tsec
  logical :: iprintprogress

  datsmooth = 0.
  term = 0.
  dmaxcoltable = 1./real(maxcoltable)
  print "(a)",'projecting from particles to pixels...'
  if (pixwidth.le.0.) then
     print "(a)",'interpolate3D_proj: error: pixel width <= 0'
     return
  endif
  !
  !--print a progress report if it is going to take a long time
  !  (a "long time" is, however, somewhat system dependent)
  !
  iprintprogress = (npart .ge. 100000) .or. (npixx*npixy .gt.100000)
  !
  !--loop over particles
  !
  iprintinterval = 25
  if (npart.ge.1e6) iprintinterval = 10
  iprintnext = iprintinterval
!
!--get starting CPU time
!
  call cpu_time(t_start)
  
  over_particles: do i=1,npart
     !
     !--report on progress
     !
     if (iprintprogress) then
        iprogress = 100*i/npart
        if (iprogress.ge.iprintnext) then
           write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
           iprintnext = iprintnext + iprintinterval
        endif
     endif
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) then
        print*,'interpolate3D_proj: error: h <= 0 ',i,hi
        return
     endif
     hi1 = 1./hi
     radkern = 2.*hi  !radius of the smoothing kernel
     !         const = 10./(7.*pi*h2)  ! normalisation constant
     const = hi1*hi1
     if (rho(i).gt.0.) then
        rho1i = 1./rho(i)
     else
        rho1i = 0.
     endif
     term = const*pmass(i)*dat(i)*rho1i
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((x(i) - radkern - xmin)/pixwidth)
     jpixmin = int((y(i) - radkern - ymin)/pixwidth)
     ipixmax = int((x(i) + radkern - xmin)/pixwidth)
     jpixmax = int((y(i) + radkern - ymin)/pixwidth)

     if (ipixmin.lt.1) ipixmin = 1  ! make sure they only contribute
     if (jpixmin.lt.1) jpixmin = 1  ! to pixels in the image
     if (ipixmax.gt.npixx) ipixmax = npixx ! (note that this optimises
     if (jpixmax.gt.npixy) jpixmax = npixy !  much better than using min/max)
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        ypix = ymin + (jpix)*pixwidth - 0.5*pixwidth
        dy = ypix - y(i)
        do ipix = ipixmin,ipixmax
           xpix = xmin + (ipix)*pixwidth - 0.5*pixwidth
           dx = xpix - x(i)
           rab = sqrt(dx**2 + dy**2)
           qq = rab*hi1
           !
           !--SPH kernel - integral through cubic spline
           !  interpolate from a pre-calculated table
           !                     
           !!--find nearest index in table
           index = int(0.5*qq*maxcoltable) + 1
           if (index.lt.maxcoltable) then
              index1 = index + 1
              !--find increment along from this index
              dxx = qq - index*maxcoltable
              !--find gradient                  
              dwdx = (coltable(index1) - coltable(index))*dmaxcoltable
              !--compute value of integrated kernel
              wab = coltable(index) + dwdx*dxx
              !
              !--calculate data value at this pixel using the summation interpolant
              !                  
              datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab          
           endif

        enddo
     enddo

  enddo over_particles
!
!--get ending CPU time
!
  call cpu_time(t_end)
  t_used = t_end - t_start
  if (t_used.gt.60) then
     itmin = int(t_used/60)
     tsec = t_used - (itmin*60)
     print*,'completed in ',itmin,' min ',tsec,'s'
  else
     print*,'completed in ',t_used,'s'
  endif
  
  return

end subroutine interpolate3D_projection

!--------------------------------------------------------------------------
!
!     Same as previous but for a vector quantity
!
!     Input: particle coordinates  : x,y   (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            vector data to smooth : vecx  (npart)
!                                    vecy  (npart)
!
!     Output: smoothed vector field   : vecsmoothx (npixx,npixy)
!                                     : vecsmoothy (npixx,npixy)
!
!     Daniel Price 23/12/04
!--------------------------------------------------------------------------

subroutine interpolate3D_proj_vec(x,y,pmass,rho,hh,vecx,vecy,npart,&
     xmin,ymin,vecsmoothx,vecsmoothy,npixx,npixy,pixwidth)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,pmass,rho,hh,vecx,vecy
  real, intent(in) :: xmin,ymin,pixwidth
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx, vecsmoothy

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: index,index1
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: rho1i,termx,termy,dx,dy,xpix,ypix
  real :: dxx,dwdx,dmaxcoltable

  vecsmoothx = 0.
  vecsmoothy = 0.
  termx = 0.
  termy = 0.
  dmaxcoltable = 1./real(maxcoltable)
  print "(a)",'projecting vector from particles to pixels...'
  if (pixwidth.le.0.) then
     print "(a)",'interpolate3D_proj_vec: error: pixel width <= 0'
     return
  endif
  !
  !--loop over particles
  !      
  over_particles: do i=1,npart
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) then
        print*,'interpolate3D_proj_vec: error: h <= 0 ',i,hi
        return
     endif
     hi1 = 1./hi
     radkern = 2.*hi    ! radius of the smoothing kernel
     const = hi1*hi1
     if (rho(i).ne.0.) then
        rho1i = 1./rho(i)
     else
        rho1i = 0.
     endif
        
     termx = const*pmass(i)*vecx(i)*rho1i
     termy = const*pmass(i)*vecy(i)*rho1i
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((x(i) - radkern - xmin)/pixwidth)
     jpixmin = int((y(i) - radkern - ymin)/pixwidth)
     ipixmax = int((x(i) + radkern - xmin)/pixwidth)
     jpixmax = int((y(i) + radkern - ymin)/pixwidth)

     ! PRINT*,'particle ',i,' x, y, z = ',x(i),y(i),z(i),dat(i),rho(i),hi
     ! PRINT*,'pixels = ',ipixmin,ipixmax,jpixmin,jpixmax

     if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
     if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
     if (ipixmax.gt.npixx) ipixmax = npixx
     if (jpixmax.gt.npixy) jpixmax = npixy
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        ypix = ymin + (jpix)*pixwidth - 0.5*pixwidth
        dy = ypix - y(i)
        do ipix = ipixmin,ipixmax
           xpix = xmin + (ipix)*pixwidth - 0.5*pixwidth
           dx = xpix - x(i)
           rab = sqrt(dx**2 + dy**2)
           qq = rab*hi1
           !
           !--SPH kernel - integral through cubic spline
           !  interpolate from a pre-calculated table
           !                     
           !!--find nearest index in table
           index = int(0.5*qq*maxcoltable) + 1
           if (index.lt.maxcoltable) then
              index1 = index + 1
              !--find increment along from this index
              dxx = qq - index*maxcoltable
              !--find gradient                  
              dwdx = (coltable(index1) - coltable(index))*dmaxcoltable
              !--compute value of integrated kernel
              wab = coltable(index) + dwdx*dxx
              !
              !--calculate data value at this pixel using the summation interpolant
              !  
              vecsmoothx(ipix,jpix) = vecsmoothx(ipix,jpix) + termx*wab          
              vecsmoothy(ipix,jpix) = vecsmoothy(ipix,jpix) + termy*wab          

           endif

        enddo
     enddo

  enddo over_particles

  return

end subroutine interpolate3D_proj_vec

end module projections3D
