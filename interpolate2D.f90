!--------------------------------------------------------------------------
!     program to interpolate from particle data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_b m_b dat_b/rho_b W(r-r_b, h_b)
! 
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline
!
!     Input: particle coordinates  : x,y   (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data 	   : datsmooth (npixx,npixy)
!
!     Daniel Price, Institute of Astronomy, Cambridge, July 2003
!--------------------------------------------------------------------------

subroutine interpolate2D(x,y,pmass,rho,hh,dat,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidth)

  implicit none
  real, parameter :: pi = 3.1415926536      
  integer, intent(IN) :: npart,npixx,npixy
  real, intent(IN), dimension(npart) :: x,y,pmass,rho,hh,dat
  real, intent(IN) :: xmin,ymin,pixwidth
  real, intent(OUT), dimension(npixx,npixy) :: datsmooth

  integer :: i,j,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,h2,radkern,qq,wab,rab,const
  real :: term,dx,dy,xpix,ypix

  datsmooth = 0.
  term = 0.
  print*,'interpolating from particles to 2D grid...'
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
        print*,'interpolate2D: error: h <= 0 ',i,hi
	return
     endif
     hi1 = 1./hi
     h2 = hi*hi
     radkern = 2.*hi  ! radius of the smoothing kernel
     const = 10./(7.*pi*h2)  ! normalisation constant
     if (rho(i).ne.0.) term = pmass(i)*dat(i)/rho(i) 
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((x(i) - radkern - xmin)/pixwidth)
     jpixmin = int((y(i) - radkern - ymin)/pixwidth)
     ipixmax = int((x(i) + radkern - xmin)/pixwidth)
     jpixmax = int((y(i) + radkern - ymin)/pixwidth)

     if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
     if (ipixmax.lt.1) ipixmax = 1
     if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
     if (jpixmax.lt.1) jpixmax = 1
     if (ipixmax.gt.npixx) ipixmax = npixx
     if (ipixmin.gt.npixx) ipixmin = npixx
     if (jpixmax.gt.npixy) jpixmax = npixy
     if (jpixmin.gt.npixy) jpixmin = npixy
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        do ipix = ipixmin,ipixmax

           ypix = ymin + (jpix)*pixwidth - 0.5*pixwidth
           xpix = xmin + (ipix)*pixwidth - 0.5*pixwidth
           dy = ypix - y(i)
           dx = xpix - x(i)
           rab = sqrt(dx**2 + dy**2)
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
           datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab

        enddo
     enddo

  enddo

  return

end subroutine interpolate2D
