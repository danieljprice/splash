!--------------------------------------------------------------------------
!     calculates a smoothed array on a fine grid given a list of particle
!     co-ordinates and a scalar field on the particles
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth = sum_b m_b dat_b/rho_b W(r-r_b, h_b)
! 
!     where _b is the quantity at the neighbouring particle b.
!
!--------------------------------------------------------------------------

      SUBROUTINE smooth_pixels(x,y,pmass,rho,hh,dat,npart,ntot,
     &                         xpix,ypix,datsmooth,npixx,npixy,
     &                         xmin,xmax,ymin,ymax)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: npixx,npixy,npart,ntot
      INTEGER, DIMENSION(100,100) :: ifirstincell
      INTEGER, DIMENSION(ntot) :: ll
      INTEGER, DIMENSION(ntot) :: neighlist
      INTEGER, DIMENSION(9) :: neighcellx,neighcelly
      INTEGER :: i,j,icellx,icelly,ncellsx,ncellsy,icell
      INTEGER :: ipart,ineigh,nneigh,ipix,jpix
      INTEGER :: ipixmin,ipixmax,jpixmin,jpixmax
      INTEGER :: nneighcell

      REAL, PARAMETER :: pi = 3.1415926536
      REAL, INTENT(IN), DIMENSION(ntot) :: x,y,pmass,rho,hh,dat
      REAL, INTENT(IN), DIMENSION(npixx,npixy) :: xpix,ypix
      REAL, INTENT(OUT), DIMENSION(npixx,npixy) :: datsmooth
      REAL :: xmin,xmax,ymin,ymax
      REAL :: dxcell,dycell,hhmax
      REAL :: wab,rab,qq,dx,dy,hhi,h2,const
      REAL :: pixpercellx,pixpercelly

!------------------------------------------------------------
!     construct the link list to find particle neighbours
!------------------------------------------------------------
      PRINT*,'Constructing link list npart,ntot =',npart,ntot
      hhmax = MAXVAL(hh(1:npart))
      dxcell = 2.*hhmax		! link list cell size 2*h
      IF (dxcell.LE.0) STOP 'Error: link: max h <=0 :'
!
!--find max/min of particle distribution
!  
!      xmin = MINVAL(x(1:npart)) - 0.00001
!      xmax = MAXVAL(x(1:npart)) + 0.00001
!      ymin = MINVAL(y(1:npart)) - 0.00001
!      ymax = MAXVAL(y(1:npart)) + 0.00001
!
!--work out number of link list cells
! 
      ncellsx = INT((xmax - xmin)/dxcell)	! round off to lower
      ncellsy = INT((ymax - ymin)/dxcell)
!
!--adjust so that ncellsx,y exact division of xmax-xmin,ymax-ymin
!      
      dxcell = (xmax-xmin)/FLOAT(ncellsx)
      dycell = (ymax-ymin)/FLOAT(ncellsy)

      pixpercellx = npixx/FLOAT(ncellsx) ! pixels in each linklist cell
      pixpercelly = npixy/FLOAT(ncellsy) ! (not ghost cells)

      print*,' pixpercell = ',pixpercellx,pixpercelly,ncellsx,ncellsy     
!
!--now include ghost cells
!      
      IF (ntot.GT.npart) THEN
         ncellsx = ncellsx + 2
         ncellsy = ncellsy + 2
         xmin = xmin - dxcell
         xmax = xmax + dxcell
         ymin = ymin - dycell
         ymax = ymax + dycell
      ENDIF
      
      PRINT*,' ncells x,y hhmax = ',ncellsx,ncellsy,hhmax
      IF ((ncellsx .EQ. 0).OR.(ncellsy.EQ.0)) THEN
         PRINT*,'Error: link: number of cells=0:'
         PRINT*,'xmin,xmax,dxcell,hhmax=',xmin,xmax,dxcell,dycell,hhmax
         STOP
      ELSEIF ((ncellsx .GT. 100).OR.(ncellsy .GT. 100)) THEN	 
         PRINT*,'Error: link: number of cells exceeds dimensions '
	 PRINT*,' ncellsx, ncellsy = ',ncellsx,ncellsy
	 STOP
      ENDIF 

      ifirstincell(:,:) = -1 		! set all head of chains to -1
!
!--construct the link list
!  
      DO i=1,ntot		! including ghosts
         icellx = INT((x(i)-xmin)/dxcell) + 1
         icelly = INT((y(i)-ymin)/dycell) + 1 
         ll(i) = ifirstincell(icellx,icelly) ! link to previous start of chain
         ifirstincell(icellx,icelly) = i     ! set head of chain to current particle	 
      ENDDO

!      print*,' xmin,xmax,dx,dycell,hhmax',xmin,xmax,dxcell,dycell,hhmax
!      read*

!------------------------------------------------------------
!     now use the link list to find the smoothed array
!------------------------------------------------------------

      PRINT*,'Calculating smoothed image...'
!
!--set datsmooth to zero for all pixels
!      
      datsmooth(:,:) = 0.
!
!--loop over all cells containing particles (ie not ghosts)
!      
      overcellsy: DO icelly=2,ncellsy-1
         overcellsx: DO icellx=2,ncellsx-1
!
!--construct list of neighbouring cells for this particle
!	    
	    nneighcell = 9		! nine cells within 2h (incl. this one)
            neighcelly(1:3) = icelly	! this row of cells
	    neighcellx(1) = icellx      !  this cell
	    neighcellx(2) = icellx-1    !  cell to left
	    neighcellx(3) = icellx+1    !  cell to right
	    neighcelly(4:6) = icelly+1	! row above
	    neighcellx(4) = icellx	!  centre
	    neighcellx(5) = icellx-1	!  left
	    neighcellx(6) = icellx+1    !  right
	    neighcelly(7:9) = icelly-1	! row below
	    neighcellx(7) = icellx	!  centre
	    neighcellx(8) = icellx-1    !  left
	    neighcellx(9) = icellx+1    !  right	   
!
!--construct list of neighbouring particles for current cell
!	    
	    nneigh = 0		! start off with zero neighbours
	    DO icell = 1,nneighcell
	       ipart = ifirstincell(neighcellx(icell),neighcelly(icell))
	       DO WHILE (ipart.NE.-1)
	          nneigh = nneigh + 1	          
		  neighlist(nneigh) = ipart
!		  PRINT*,' neighlist(',nneigh,'), next = ',ipart,ll(ipart)
		  ipart = ll(ipart) 	       
	       ENDDO
	    ENDDO 
!	    PRINT*,' cell ',icellx,icelly,', neighbours = ',nneigh
!
!--compute datsmooth for all the pixels in this cell
!	    
            jpixmin = NINT(pixpercelly*(icelly-2))+1
	    jpixmax = MIN(NINT(pixpercelly*(icelly-1)),npixy)
	    ipixmin = NINT(pixpercellx*(icellx-2))+1
	    ipixmax = MIN(NINT(pixpercellx*(icellx-1)),npixx)

!	    PRINT*,' loop over pixels j=',jpixmin,jpixmax
!	    PRINT*,' loop over pixels i=',ipixmin,ipixmax

!            PRINT*,' this cell x = ',xmin + (icellx-1)*dxcell,
!     &	                            xmin + (icellx)*dxcell
!            PRINT*,' this cell y = ',ymin + (icelly-1)*dycell,
!     &	                            ymin + (icelly)*dycell

!	    ipart = ifirstincell(icellx,icelly)
!	    PRINT*,'first particle this cell = ',x(ipart),y(ipart)
!    ipart = ll(ipart)
!	    PRINT*,'next = ',x(ipart),y(ipart)
!            PRINT*,' pixels this cell x = ',xpix(ipixmin,jpixmin),
!     &	                            xpix(ipixmax,jpixmax)
!            PRINT*,' pixels this cell y = ',ypix(ipixmin,jpixmin),
!     &	                            ypix(ipixmax,jpixmax)
!            READ*
	    overpixy: DO jpix = jpixmin,jpixmax
	       overpixx: DO ipix = ipixmin,ipixmax
	       
!                  PRINT*,'current pixel = ',ipix,jpix,
!     &		          xpix(ipix,jpix),ypix(ipix,jpix)
!
!--sum over all the particles
!
	          overneigh: DO ineigh = 1,nneigh
		     ipart = neighlist(ineigh)
		     dy = ypix(ipix,jpix) - y(ipart)
		     dx = xpix(ipix,jpix) - x(ipart)
		     rab = SQRT(dx**2 + dy**2)
		     hhi = hh(ipart)
		     h2 = hhi*hhi
		     qq = rab/hhi
!		     PRINT*,'neighbour part ',ipart,x(ipart),y(ipart),qq
!
!--SPH kernel - standard cubic spline
!		     
                     const = 10./(7.*pi*h2)
		     IF (qq.LT.1.0) THEN
		        wab = const*(1.-1.5*qq**2 + 0.75*qq**3)
		     ELSEIF (qq.LT.2.0) THEN
		        wab = const*0.25*(2.-qq)**3.
		     ELSE
                        wab = 0.
                     ENDIF
!
!--calculate data value at this pixel using the summation interpolant
!		     
		     datsmooth(ipix,jpix) = datsmooth(ipix,jpix)
     &                + pmass(ipart)*dat(ipart)/rho(ipart)*wab		     
		     
!		     print*,'datsmooth = ',datsmooth(ipix,jpix),ipix,jpix
!		     READ*
		  ENDDO overneigh
!		  print*,'datsmooth = ',datsmooth(ipix,jpix),ipix,jpix		  
!		  read*
	       ENDDO overpixx
	    ENDDO overpixy
	    
	 ENDDO overcellsx
      ENDDO overcellsy
      
      END SUBROUTINE smooth_pixels
