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
!     this program is *pretty slow* for any reasonable number 
!     of particles/pixels - only efficient if you want *lots* of 
!     cross sections through the same data
!
!     Input: particle coordinates  : x,y   (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data 	   : datsmooth (npixx,npixy,npixz)
!
!     Daniel Price, Institute of Astronomy, Cambridge 16/7/03
!--------------------------------------------------------------------------

      SUBROUTINE interpolate3D(x,y,z,pmass,rho,hh,dat,npart,
     &    xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidth)
      
      IMPLICIT NONE
      REAL, PARAMETER :: pi = 3.1415926536      
      INTEGER, INTENT(IN) :: npart,npixx,npixy,npixz
      REAL, INTENT(IN), DIMENSION(npart) :: x,y,z,pmass,rho,hh,dat
      REAL, INTENT(IN) :: xmin,ymin,zmin,pixwidth
      REAL, INTENT(OUT), DIMENSION(npixx,npixy,npixz) :: datsmooth

      INTEGER :: i,j,ipix,jpix,kpix
      INTEGER :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
      REAL :: hi,hi1,h3,radkern,qq,wab,rab,const
      REAL :: term,dx,dy,dz,xpix,ypix,zpix

      datsmooth = 0.
      term = 0.
      PRINT*,'interpolating from particles to 3D grid...'
!
!--loop over particles
!      
      DO i=1,npart
!
!--set kernel related quantities
!
         hi = hh(i)
	 hi1 = 1./hi
	 h3 = hi*hi*hi
	 radkern = 2.*hi	! radius of the smoothing kernel
         const = 1./(pi*h3)	! normalisation constant (3D)
	 term = 0.
	 IF (rho(i).NE.0.) term = pmass(i)*dat(i)/rho(i) 
!
!--for each particle work out which pixels it contributes to
!               
	 ipixmin = INT((x(i) - radkern - xmin)/pixwidth)
	 jpixmin = INT((y(i) - radkern - ymin)/pixwidth)
	 kpixmin = INT((z(i) - radkern - zmin)/pixwidth)
	 ipixmax = INT((x(i) + radkern - xmin)/pixwidth)
	 jpixmax = INT((y(i) + radkern - ymin)/pixwidth)
	 kpixmax = INT((z(i) + radkern - zmin)/pixwidth)
	 
!	 PRINT*,'particle ',i,' x, y = ',x(i),y(i),dat(i),rho(i),hi
!	 PRINT*,'pixels = ',ipixmin,ipixmax,jpixmin,jpixmax
!	 PRINT*,'xmin,ymin = ',xmin,ymin

         IF (ipixmin.LT.1) ipixmin = 1	! make sure they only contribute
	 IF (jpixmin.LT.1) jpixmin = 1  ! to pixels in the image
	 IF (kpixmin.LT.1) kpixmin = 1
	 IF (ipixmax.GT.npixx) ipixmax = npixx
	 IF (jpixmax.GT.npixy) jpixmax = npixy
	 IF (kpixmax.GT.npixz) kpixmax = npixz
!
!--loop over pixels, adding the contribution from this particle
!
         DO kpix = kpixmin,kpixmax
	    DO jpix = jpixmin,jpixmax
	       DO ipix = ipixmin,ipixmax
      
                  zpix = zmin + (kpix)*pixwidth - 0.5*pixwidth
		  ypix = ymin + (jpix)*pixwidth - 0.5*pixwidth
		  xpix = xmin + (ipix)*pixwidth - 0.5*pixwidth
 		  dz = zpix - z(i)
		  dy = ypix - y(i)
		  dx = xpix - x(i)
		  rab = SQRT(dx**2 + dy**2 + dz**2)
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
		  datsmooth(ipix,jpix,kpix) = datsmooth(ipix,jpix,kpix)
     &		                            + term*wab		          
      
               ENDDO
	    ENDDO
         ENDDO
      ENDDO
      
      RETURN
      
      END SUBROUTINE interpolate3D
