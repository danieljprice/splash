!!
!!   calculates various additional quantities from the input data
!!
SUBROUTINE calc_quantities
 USE labels
 USE particle_data
 USE settings
 IMPLICIT NONE
 INTEGER :: i,j
 REAL :: Bmag
 REAL, PARAMETER :: pi = 3.1415926536
 REAL :: angledeg,anglexy,runit(ndimmax)	! to plot r at some angle

 PRINT*,'calculating ',ncalc,' additional quantities...'
 numplot = ncolumns + ncalc
 DO i=nstart,n_end
    DO j=1,npart(i)
!!--pressure if not in data array
       IF ((ipr.GT.ncolumns)                            &
	   .AND.(irho.NE.0).AND.(iutherm.NE.0))         &
 	   dat(ipr,j,i) = dat(irho,j,i)*dat(iutherm,j,i)*(gamma(i)-1.)
!!--entropy
       IF (ientrop.NE.0) THEN
  	  IF (dat(irho,j,i).GT.1.e-10) THEN
             dat(ientrop,j,i) = dat(ipr,j,i)/dat(irho,j,i)**gamma(i)
          ELSE  
	     dat(ientrop,j,i) = 0.
	  ENDIF		   
       ENDIF   
!!--radius	 
       IF (irad.NE.0) &    
          dat(irad,j,i) = SQRT(DOT_PRODUCT(dat(ix(1:ndim),j,i),   &
                                           dat(ix(1:ndim),j,i)))
!!--distance along the line with angle 30deg w.r.t. x axis
       IF (irad2.NE.0) THEN
          angledeg = 30.
	  anglexy = angledeg*pi/180.
	  runit(1) = COS(anglexy)
	  IF (ndim.GE.2) runit(2) = SIN(anglexy)
	  dat(irad2,j,i) = DOT_PRODUCT(dat(ix(1:ndim),j,i),runit(1:ndim))
       ENDIF
!!--specific KE
       IF ((ike.NE.0).AND.(ivx.NE.0)) &
          dat(ike,j,i) = 0.5*DOT_PRODUCT(dat(ivx:ivlast,j,i),dat(ivx:ivlast,j,i))
!!--magnetic pressure
       IF ((ipmag.NE.0).AND.(iBfirst.NE.0))   &
          dat(ipmag,j,i) = 0.5*DOT_PRODUCT(dat(iBfirst:iBlast,j,i), &
                                dat(iBfirst:iBlast,j,i))
!!--plasma beta
       IF ((ibeta.NE.0).AND.(ipmag.NE.0)) THEN
	  IF (abs(dat(ipmag,j,i)).GT.1e-10) THEN
             dat(ibeta,j,i) = dat(ipr,j,i)/dat(ipmag,j,i)
	  ELSE  
	     dat(ibeta,j,i) = 0.
	  ENDIF
       ENDIF   
!!--total pressure (gas + magnetic)     
       IF ((itotpr.NE.0).AND.(ipr.NE.0).AND.(ipmag.NE.0)) THEN
	  dat(itotpr,j,i) = dat(ipr,j,i) + dat(ipmag,j,i)
       ENDIF
!!--div B error	(h*divB / abs(B))	
       IF ((idivBerr.NE.0).AND.    &
          (idivB.NE.0).AND.(ih.NE.0).AND.(iBfirst.NE.0)) THEN
	  Bmag = SQRT(DOT_PRODUCT(dat(iBfirst:iBlast,j,i),dat(iBfirst:iBlast,j,i)))
          IF (Bmag.GT.0.) THEN
	     dat(idivBerr,j,i) = abs(dat(idivB,j,i))*dat(ih,j,i)/Bmag
	  ELSE
   	     dat(idivBerr,j,i) = 0.
	  ENDIF
       ENDIF
!!--h*SQRT(rho)/abs(B)	(timestep)
       IF ((itimestep.NE.0).AND.   &
       	  (ih.NE.0).AND.(irho.NE.0).AND.(iBfirst.NE.0)) THEN
	  Bmag = SQRT(DOT_PRODUCT(dat(iBfirst:iBlast,j,i),  &
                                  dat(iBfirst:iBlast,j,i)))
          IF (Bmag.NE.0.) THEN
	     dat(itimestep,j,i) = dat(ih,j,i)*SQRT(dat(irho,j,i))/Bmag
	  ELSE
	     dat(itimestep,j,i) = 0.
	  ENDIF
       ENDIF

    ENDDO
  ENDDO 
 RETURN
END SUBROUTINE calc_quantities
