!-------------------------------------------------------------------
!  This subroutine takes an array of particle positions in x,y
!  and a scalar or 2D vector array 'dat' defined on the particles
!  and interpolates to a grid which can then be rendered.
!
!  interpolation to coarser grid is just by averaging
!------------------------------------------------------------------

      SUBROUTINE coarse_render(x,y,xmin,xmax,dat,datmin,datmax,
     &                  ntot,npix,ivec,icolours,iplotcont)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ntot,npix,ivec,icolours
      INTEGER i,j,k
      INTEGER, PARAMETER :: nc=30	! number of contours
      INTEGER ihoc(npix,npix),numcell(npix,npix),ll(ntot),ix,iy
      LOGICAL, INTENT(IN) :: iplotcont
      REAL, INTENT(IN), DIMENSION(ntot) :: x,y
      REAL dat(ivec,ntot)
      REAL xmin,xmax,dx,scale,zero
      REAL datpix(ivec,npix,npix)
      REAL trans(6),levels(nc),dcont,datmin,datmax
      
      zero= 0.0E0
!     set up grid for rendering
      
      xmax = xmax + 0.000001
      xmin = xmin - 0.000001
      dx = (xmax-xmin)/REAL(npix)

      trans(1) = xmin-0.5*dx		! this is for the PGIMAG call
      trans(2) = dx			! see help for PGIMAG/PGGRAY/PGCONT
      trans(3) = 0.0
      trans(4) = xmin-0.5*dx
      trans(5) = 0.0
      trans(6) = dx
!
!--if doing vector plot, interpolation is to a coarser grid, so just average      
!
!     set up link list for which particles in a particular cell      
      ihoc(:,:) = -1
      numcell(:,:) = 0
      DO i=1,ntot
         ix = INT((x(i)-xmin)/dx)+1
	 iy = INT((y(i)-xmin)/dx)+1
         IF ((ix.lt.1).or.(ix.gt.npix)) THEN
!	    PRINT*,'particle ',i,' not in x domain, x =',x(i)
	 ELSEIF ((iy.lt.1).or.(iy.gt.npix)) THEN
!	    PRINT*,'particle ',i,' not in y domain, y =',y(i)
	 ELSE
  	    ll(i)=ihoc(ix,iy) 
	    ihoc(ix,iy) = i	   
	 ENDIF   
      ENDDO

!     now work out average in each cell      
      datpix(:,:,:) = 0.0
      DO i=1,npix
         DO j=1,npix
   	    k = ihoc(i,j)
 	    DO WHILE (k.ne.-1)	 
               datpix(:,i,j) = datpix(:,i,j) + dat(:,k)
	       numcell(i,j) = numcell(i,j) + 1
	       k = ll(k)
	    ENDDO
	 ENDDO	 
      ENDDO
      
      DO i=1,npix
         DO j=1,npix
	    IF (numcell(i,j).ne.0) THEN
	       datpix(:,i,j) = datpix(:,i,j)/float(numcell(i,j))
	    ENDIF
	 ENDDO
      ENDDO

      IF (ivec.eq.1) THEN
         PRINT*,'rendering...'
         dcont = (datmax-datmin)/REAL(nc+1)   ! even contour levels
         DO i=1,nc
            levels(i) = datmin + REAL(i)*dcont
         ENDDO
!
!--nb: plots use my modification of PGWEDG which plots vertical numbers on axes
!	 
         IF (icolours.EQ.1) THEN	! greyscale
            CALL DANPGWEDG('RG',1.0,4.0,datmin,datmax,' ')
            CALL PGGRAY(datpix,npix,npix,1,npix,1,npix,
     &                                   datmin,datmax,trans)
     
         ELSEIF (icolours.GT.1) THEN	! colour
            CALL DANPGWEDG('RI',1.0,4.0,datmin,datmax,' ')
!         CALL PGWEDG('RI',2.0,4.0,datmin,datmax,' ')
!         CALL PGPIXL(datpix,npix,npix,1,npix,1,npix,xmin,xmax,ymin,ymax)
            CALL PGIMAG(datpix,npix,npix,1,npix,1,npix,
     &                                   datmin,datmax,trans)
	 
	 ENDIF

         IF (iplotcont)
     &	    CALL PGCONT(datpix,npix,npix,1,npix,1,npix,levels,nc,trans)
     
      ELSEIF (ivec.eq.2) THEN
         CALL PGSAH(2,45.0,0.7)   ! arrow style
	 CALL PGSCH(0.3)	  ! size of arrow head
	 IF (datmax.eq.0.0) THEN
	    scale=0.0
	 ELSE
	    scale=0.1/datmax
	 ENDIF
	 PRINT*,'vector map...',scale
	 CALL PGVECT(datpix(1,:,:),datpix(2,:,:),npix,npix,
     &        1,npix,1,npix,scale,0,trans,-1000.0)   
         CALL PGSCH(1.0)
      ENDIF

      RETURN
      
      END SUBROUTINE coarse_render
