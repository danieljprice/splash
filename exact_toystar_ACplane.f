      SUBROUTINE exact_toystar_ACplane(Astart,Cstart,sigma,gamma)
      IMPLICIT NONE
      REAL, INTENT(IN) :: Astart,Cstart,sigma,gamma
      REAL :: func,func2,constK,gam1,gamm1,gamp1,fact
      REAL xstart,xend,xcentre,C,Cnew,k,kmin,kmax
      REAL funct,fderiv,ymin,ymax,extra
      INTEGER i
      EXTERNAL func,func2
      COMMON /kconst/ k,fact,gam1,gamp1
      
      PRINT*,' plotting A-C plane...'
                   
      gamp1 = gamma + 1.
      gamm1 = gamma - 1.
      gam1 = 1./gamm1
      constK = 0.25
!      PRINT*,' K, Kdash = ',constK,constK + 0.5*sigma**2
      fact = 2.*(constK + 0.5*sigma**2)*gamma*gam1

      k = (Astart**2 + 1. + 2.*fact*Cstart*gam1)*Cstart**(-2./gamp1)
      
!      PRINT*,' k,fact = ',k,fact
!
!--find limits of plot (ie. where A = 0)
!
      C = 1.e6
      Cnew = 0.25
      
      DO WHILE (abs(C-Cnew).GT.1.e-5)
      	 
	 C = Cnew

	 funct = k*C**(2./gamp1) - 2.*fact*C*gam1 - 1.
	 fderiv = 2.*k/gamp1*C**(-gamm1/gamp1) - 2.*fact*gam1
	 
	 Cnew = C - funct/fderiv
	 
	 IF (Cnew.LT.0.) PRINT*,'eek C < 0'
	 	 
      ENDDO	 
      
      xstart = Cnew
      
      C = 1.e6
      Cnew = 6.37935
      
      DO WHILE (abs(C-Cnew).GT.1.e-5)
      
         C = Cnew
      
	 funct = k*C**(2./gamp1) - 2.*fact*C*gam1 - 1.
	 fderiv = 2.*k/gamp1*C**(-gamm1/gamp1) - 2.*fact*gam1        
	 
	 Cnew = C - funct/fderiv
	 
      ENDDO	       
      
      xend = Cnew
     
!      PRINT*,'plotting k = ',k,' Cstart = ',Cstart,' Astart = ',Astart
!      PRINT*,'min C = ',xstart,' max C = ',xend
      
      xstart = xstart + 0.000001
      xend = xend - 0.000001
      
      extra = 0.1*(xend-xstart)
      xcentre = 0.5*(xstart + xend)
      ymax = 1.5*func(xcentre)
      ymin = 1.5*func2(xcentre)

      CALL PGSWIN(xstart-extra,xend+extra,ymin,ymax,0,1)
      CALL PGBOX('BCNST',0.0,0,'1BVCNST',0.0,0)      
      CALL PGFUNX(func,10000,xstart,xend,1)
      CALL PGFUNX(func2,10000,xstart,xend,1)

      CALL PGLABEL ('C','A',' ')
      
      RETURN
      
      END SUBROUTINE exact_toystar_ACplane
      
      FUNCTION func(x)
      REAL x,func,k,term
      COMMON /kconst/ k,fact,gam1,gamp1
      
!      print*,'k = ',k
      
      term = -1 -2.*fact*x*gam1 + k*x**(2./gamp1)
      IF (term.LE.0.) THEN
!         PRINT*,' warning: func < 0 ',term
         func = 0.
      ELSE	 
         func = SQRT(term)
      ENDIF	 
      
      END 
      
      FUNCTION func2(x)
      REAL x,func2,k
      COMMON /kconst/ k,fact,gam1,gamp1
      
!      print*,' k = ',k
      
      term = -1 -2.*fact*x*gam1 + k*x**(2./gamp1)
      IF (term.LE.0.) THEN
         func2 = 0.
      ELSE
         func2 = -SQRT(term)
      ENDIF	 
      
      END 
