!--------------------------------------------------------------------------
!     program to interpolate from particle data to even 1D grid of pixels
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
!     Output: smoothed data 	   : datsmooth (npixx)
!
!     Daniel Price, Institute of Astronomy, Cambridge, Dec 2003
!--------------------------------------------------------------------------

      SUBROUTINE interpolate1D(x,pmass,rho,hh,dat,npart,
     &    xmin,datsmooth,npixx,pixwidth)
      
      IMPLICIT NONE
      REAL, PARAMETER :: pi = 3.1415926536      
      INTEGER, INTENT(IN) :: npart,npixx
      REAL, INTENT(IN), DIMENSION(npart) :: x,pmass,rho,hh,dat
      REAL, INTENT(IN) :: xmin,pixwidth
      REAL, INTENT(OUT), DIMENSION(npixx) :: datsmooth

      INTEGER :: i,j,ipix,ipixmin,ipixmax
      REAL :: hi,hi1,radkern,qq,wab,rab,const
      REAL :: term,dx,xpix

      datsmooth = 0.
      term = 0.
      PRINT*,'interpolating from particles to 1D grid...'
!
!--loop over particles
!      
      DO i=1,npart
!
!--set kernel related quantities
!
         hi = hh(i)
	 hi1 = 1./hi
	 radkern = 2.*hi	! radius of the smoothing kernel
         const = 2./(3.*hi)	! normalisation constant
	 IF (rho(i).NE.0.) term = pmass(i)*dat(i)/rho(i) 
!
!--for each particle work out which pixels it contributes to
!               
	 ipixmin = INT((x(i) - radkern - xmin)/pixwidth)
	 ipixmax = INT((x(i) + radkern - xmin)/pixwidth)
	 
         IF (ipixmin.LT.1) ipixmin = 1		! make sure they only contribute
	 IF (ipixmax.GT.npixx) ipixmax = npixx ! to pixels in the image
!
!--loop over pixels, adding the contribution from this particle
!
         DO ipix = ipixmin,ipixmax
  	    xpix = xmin + (ipix)*pixwidth - 0.5*pixwidth
 	    dx = xpix - x(i)
	    rab = abs(dx)
 	    qq = rab*hi1
!
!--SPH kernel - standard cubic spline
!		     
            IF (qq.LT.1.0) THEN
	       wab = const*(1.-1.5*qq**2 + 0.75*qq**3)
	    ELSEIF (qq.LT.2.0) THEN
	       wab = const*0.25*(2.-qq)**3
	    ELSE
               wab = 0.
            ENDIF
!
!--calculate data value at this pixel using the summation interpolant
!		  
	    datsmooth(ipix) = datsmooth(ipix) + term*wab		          
         ENDDO

      ENDDO
      
      RETURN
      
      END SUBROUTINE interpolate1D
