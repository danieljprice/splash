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
!     ** In this version 3D data is interpolated to a single 2D cross section
!     ** This is much faster than interpolating to a 3D grid
!     ** and is efficient if only one or two cross sections are needed.
!
!     ** Note that the cross section is always taken in the z co-ordinate
!     ** so should submit the appropriate arrays as x, y and z.
!
!     Input: particle coordinates  : x,y,z (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!            cross section location: zslice
!
!     Output: smoothed data 	   : datsmooth (npixx,npixy)
!
!     Daniel Price, Institute of Astronomy, Cambridge, 23/9/03
!--------------------------------------------------------------------------

      SUBROUTINE interpolate3D_fastxsec(
     &    x,y,z,pmass,rho,hh,dat,npart,
     &    xmin,ymin,zslice,datsmooth,npixx,npixy,pixwidth)
      
      IMPLICIT NONE
      REAL, PARAMETER :: pi = 3.1415926536      
      INTEGER, INTENT(IN) :: npart,npixx,npixy
      REAL, INTENT(IN), DIMENSION(npart) :: x,y,z,pmass,rho,hh,dat
      REAL, INTENT(IN) :: xmin,ymin,pixwidth,zslice
      REAL, INTENT(OUT), DIMENSION(npixx,npixy) :: datsmooth

      INTEGER :: i,j,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
      REAL :: hi,hi1,h3,radkern,qq,wab,rab,const
      REAL :: term,dx,dy,dz,dz2,xpix,ypix

      datsmooth = 0.
      term = 0.
      PRINT*,'taking fast cross section...',zslice
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
!
!--for each particle, work out distance from the cross section slice.
!
         dz = zslice - z(i)
	 dz2 = dz**2
!
!--if this is < 2h then add the particle's contribution to the pixels
!  otherwise skip all this and start on the next particle
!
	 IF (abs(dz) .LT. radkern) THEN

            const = 1./(pi*h3)	! normalisation constant (3D)
	    term = 0.
	    IF (rho(i).NE.0.) term = const*pmass(i)*dat(i)/rho(i) 
!
!--for each particle work out which pixels it contributes to
!               
 	    ipixmin = INT((x(i) - radkern - xmin)/pixwidth)
	    jpixmin = INT((y(i) - radkern - ymin)/pixwidth)
	    ipixmax = INT((x(i) + radkern - xmin)/pixwidth)
	    jpixmax = INT((y(i) + radkern - ymin)/pixwidth)

!	 PRINT*,'particle ',i,' x, y, z = ',x(i),y(i),z(i),dat(i),rho(i),hi
!	 PRINT*,'pixels = ',ipixmin,ipixmax,jpixmin,jpixmax
	 
            IF (ipixmin.LT.1) ipixmin = 1  ! make sure they only contribute
	    IF (jpixmin.LT.1) jpixmin = 1  ! to pixels in the image
	    IF (ipixmax.GT.npixx) ipixmax = npixx
	    IF (jpixmax.GT.npixy) jpixmax = npixy
!
!--loop over pixels, adding the contribution from this particle
!
	    DO jpix = jpixmin,jpixmax
	       ypix = ymin + (jpix)*pixwidth - 0.5*pixwidth
	       dy = ypix - y(i)
	       DO ipix = ipixmin,ipixmax
		  xpix = xmin + (ipix)*pixwidth - 0.5*pixwidth
		  dx = xpix - x(i)
		  rab = SQRT(dx**2 + dy**2 + dz2)
		  qq = rab*hi1
!
!--SPH kernel - standard cubic spline
!		     
                  IF (qq.LT.2.0) THEN
	           IF (qq.LT.1.0) THEN
		      wab = (1.-1.5*qq**2 + 0.75*qq**3)
		   ELSE
		      wab = 0.25*(2.-qq)**3
		   ENDIF
!
!--calculate data value at this pixel using the summation interpolant
!		  
 		   datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab		          
		  
		  ENDIF   
      
               ENDDO
	    ENDDO
	      
         ENDIF	! if particle within 2h of slice
      ENDDO ! over particles
      
      RETURN
      
      END SUBROUTINE interpolate3D_fastxsec
