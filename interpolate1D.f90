!----------------------------------------------------------------------
!
!  Module containing all of the routines required for 1D interpolation
!
!----------------------------------------------------------------------

module interpolations1D
 implicit none
 real, parameter, private :: pi = 3.1415926536      
 public :: interpolate1D
 
contains

!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even 1D grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_b m_b dat_b/rho_b W(r-r_b, h_b)
! 
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline
!
!     Input: particle coordinates  : x     (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx)
!
!     Daniel Price, Institute of Astronomy, Cambridge, Dec 2003
!--------------------------------------------------------------------------

subroutine interpolate1D(x,pmass,rho,hh,dat,npart,  &
     xmin,datsmooth,npixx,pixwidth)

  implicit none
  integer, intent(in) :: npart,npixx
  real, intent(in), dimension(npart) :: x,pmass,rho,hh,dat
  real, intent(in) :: xmin,pixwidth
  real, intent(out), dimension(npixx) :: datsmooth

  integer :: i,ipix,ipixmin,ipixmax
  real :: hi,hi1,radkern,qq,wab,rab,const
  real :: term,dx,xpix

  datsmooth = 0.
  term = 0.
  print*,'interpolating from particles to 1D grid: npix,xmin,max=',npixx,xmin,xmin+npixx*pixwidth
  if (pixwidth.le.0.) then
     print*,'interpolate2D: error: pixel width <= 0'
     return
  endif
  !
  !--loop over particles
  !      
  do i=1,npart
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) then
        print*,'interpolate1D: error: h <= 0 ',i,hi
        return
     endif
     hi1 = 1./hi
     radkern = 2.*hi   ! radius of the smoothing kernel
     const = 2./(3.*hi)   ! normalisation constant
     if (rho(i).ne.0.) term = pmass(i)*dat(i)/rho(i) 
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((x(i) - radkern - xmin)/pixwidth)
     ipixmax = int((x(i) + radkern - xmin)/pixwidth)

     if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
     if (ipixmax.gt.npixx) ipixmax = npixx ! to pixels in the image
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do ipix = ipixmin,ipixmax
        xpix = xmin + (ipix)*pixwidth - 0.5*pixwidth
        dx = xpix - x(i)
        rab = abs(dx)
        qq = rab*hi1
        !
        !--SPH kernel - standard cubic spline
        !     
        if (qq.lt.1.0) then
           wab = const*(1.-1.5*qq**2 + 0.75*qq**3)
        elseif (qq.lt.2.0) then
           wab = const*0.25*(2.-qq)**3
        else
           wab = 0.
        endif
        !
        !--calculate data value at this pixel using the summation interpolant
        !  
        datsmooth(ipix) = datsmooth(ipix) + term*wab          
     enddo

  enddo

  return

end subroutine interpolate1D

end module interpolations1D
