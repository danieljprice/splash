! ----------------------------------------------------------------------
! Plot exact solution for a magnetohydrodynamic shock
! (ie. one dimensional MHD Riemann problem)
!
! At the moment these are just taken from the tables in 
! Ryu & Jones (1995), ApJ 442, 228 or from ruler and pencil
! on the results in Balsara (1998)
!
! ----------------------------------------------------------------------

SUBROUTINE exact_mhdshock(iplot,ishk,time,gamma,xmin,xmax)
!     &     rho_L,rho_R,pr_L,pr_R,vx_L,vx_R,vy_L,vy_R,vz_L,vz_R,
!     &     By_L,By_R,Bx)
 IMPLICIT NONE
 INTEGER, PARAMETER :: maxpts=30
 INTEGER, INTENT(IN) :: iplot,ishk
 REAL, INTENT(in) :: time,gamma,xmin,xmax
 REAL :: rho_L,pr_L,vx_L,vy_L,vz_L,By_L,Bz_L,Bx
 REAL :: rho_R,pr_R,vx_R,vy_R,vz_R,By_R,Bz_R
 REAL, DIMENSION(maxpts) :: xplot,yplot
 REAL, DIMENSION(maxpts) :: rho,pr,vx,vy,vz,By,Bz
 REAL :: dx,x1,x2,x3,x4,x5,x6,x7, const
 REAL :: vsound2_L,vs_L,valf2_L,valf_L,vterm_L,vfast_L,vslow_L
 REAL :: vsound2_R,vs_R,valf2_R,valf_R,vterm_R,vfast_R,vslow_R

 INTEGER :: i,j,npts
 CHARACTER(LEN=20) :: filename
 LOGICAL :: iexist
      
 PRINT*,' Plotting exact mhd shock #',ishk,' at t = ',time
!
!--read parameters from .shk file
!
 filename = 'mshk'//ACHAR(ishk + 48)//'.shk'
 INQUIRE (file=filename, EXIST=iexist)
 IF (iexist) THEN
    PRINT*,' reading parameters from file ',filename
    OPEN (UNIT=55,FILE=filename,STATUS='old',FORM='formatted')
      READ(55,*) rho_L,rho_R
      READ(55,*) pr_L,pr_R
      READ(55,*) vx_L,vx_R
      READ(55,*) vy_L,vy_R
      READ(55,*) vz_L,vz_R
      READ(55,*) Bx
      READ(55,*) By_L,By_R
      READ(55,*) Bz_L,Bz_R
    CLOSE (UNIT=55)
 ENDIF
 WRITE(*,10) rho_L,rho_R,pr_L,pr_R,vx_L,vx_R,   &
                  vy_L,vy_R,vz_L,vz_R,		&
		  Bx,By_L,By_R,Bz_L,Bz_R
10 FORMAT( ' 1D shock: rho L: ',f8.3,' R: ',f8.3,/,   &
           '           pr  L: ',f8.3,' R: ',f8.3,/,   &
           '           vx  L: ',f8.3,' R: ',f8.3,/,   &
           '           vy  L: ',f8.3,' R: ',f8.3,/,   &
           '           vz  L: ',f8.3,' R: ',f8.3,/,   &
           '           Bx   : ',f8.3,/,   &	   	   	   
           '           By  L: ',f8.3,' R: ',f8.3,/,   &
           '           Bz  L: ',f8.3,' R: ',f8.3,/)

!
!--work out wave speeds in left and right medium
!
 IF (rho_L.LE.0. .OR. rho_R.LE.0.) THEN
    PRINT*,' Error: rho zero ',rho_L,rho_R
    RETURN
 ENDIF
 
 vsound2_L = gamma*pr_L/rho_L
 vsound2_R = gamma*pr_R/rho_R
 vs_L = SQRT(vsound2_L)
 vs_R = SQRT(vsound2_R)
 valf2_L = Bx**2/rho_L
 valf2_R = Bx**2/rho_R
 valf_L = SQRT(valf2_L)
 valf_R = SQRT(valf2_R)
 vterm_L = vsound2_L + (Bx**2 + By_L**2 + Bz_L**2)/rho_L
 vterm_R = vsound2_R + (Bx**2 + By_R**2 + Bz_R**2)/rho_R

 vfast_L = SQRT(0.5*(vterm_L + SQRT(vterm_L**2 - 4.*vsound2_L*valf2_L)))
 vfast_R = SQRT(0.5*(vterm_R + SQRT(vterm_R**2 - 4.*vsound2_R*valf2_R)))

 vslow_L = SQRT(0.5*(vterm_L - SQRT(vterm_L**2 - 4.*vsound2_L*valf2_L)))
 vslow_R = SQRT(0.5*(vterm_R - SQRT(vterm_R**2 - 4.*vsound2_R*valf2_R)))

 xplot(1) = xmin
 xplot(2) = (vx_L - vfast_L)*time
 xplot(3) = (vx_L - valf_L)*time
 xplot(4) = (vx_L - vslow_L)*time
 xplot(5) = (vx_L)*time
 xplot(6) = (vx_R)*time
 xplot(7) = (vx_R + vslow_R)*time
 xplot(8) = (vx_R + valf_R)*time
 xplot(9) = (vx_R + vfast_R)*time
 xplot(10) = (vx_L + vslow_L)*time
 xplot(11) = (vx_L + valf_L)*time
 xplot(12) = (vx_L + vfast_L)*time
 PRINT*, 'xplot 1-9 = ',xplot(1:12)
!
!--set up grid for exact solution
!          
 const = 1./SQRT(4.*3.1415926536)

 SELECT CASE(ishk)	! which solution to plot
    CASE(1)
!
!--Brio & Wu problem with gamma = 2
!
       vz = 0.
       Bz = 0.
       npts = 14
       xplot(1) = xmin
       xplot(2) = -0.18
       xplot(3) = -0.08
       xplot(4:6) = -0.03
       xplot(7) = -0.005
       xplot(8:9) = 0.06
       xplot(10:11) = 0.147
       xplot(12) = 0.33
       xplot(13) = 0.36
       xplot(14) = xmax
       PRINT*, 'xplot 1-14 = ',xplot(1:14)

       rho(1:2) = 1.0
       rho(3:4) = 0.675
       rho(5) = 0.815
       rho(6) = 0.775
       rho(7:8) = 0.7
       rho(9:10) = 0.234
       rho(11:12) = 0.115
       rho(13:14) = 0.125
       
       pr(1:2) = 1.0
       pr(3:4) = 0.45
       pr(5) = 0.7
       pr(6) = 0.64
       pr(7:10) = 0.52
       pr(11:12) = 0.08
       pr(13:14) = 0.1
       
       vx(1:2) = 0.0
       vx(3:4) = 0.635
       vx(5) = 0.48
       vx(6) = 0.52
       vx(7:10) = 0.595
       vx(11:12) = -0.24
       vx(13:14) = 0.0
       
       vy(1:2) = 0.0
       vy(3:4) = -0.23
       vy(5) = -1.3
       vy(6) = -1.4
       vy(7:10) = -1.58
       vy(11:12) = -0.165
       vy(13:14) = 0.
       
       By(1:2) = 1.0
       By(3:4) = 2.1*const
       By(5) = -1.2*const
       By(6) = -1.3*const
       By(7:10) = -1.9*const
       By(11:12) = -3.25*const
       By(13:14) = -1.0

    CASE(4)
!
!--isothermal MHD problem from Balsara (1998)
!    
       npts = 14
       xplot(1) = xmin
       xplot(2:3) = -0.15
       xplot(4:5) = 0.035
       xplot(6:7) = 0.07
       xplot(8:9) = 0.17 
       xplot(10:11) = 0.2 
       xplot(12:13) = 0.41
       xplot(14) = xmax
              
       rho(1:2) = 1.08
       rho(3:6) = 1.515
       rho(7:8) = 1.745
       rho(9:12) = 1.36
       rho(13:14) = 1.0
              
       vx(1:2) = 1.2 
       vx(3:6) = 0.65
       vx(7:8) = 0.62
       vx(9:12) = 0.54
       vx(13:14) = 0.0
       
       vy(1:2) = 0.01
       vy(3:4) = 0.13
       vy(5:6) = 0.24
       vy(7:8) = 0.071
       vy(9:10) = -0.215
       vy(11:12) = -0.125
       vy(13:14) = 0.0
       
       vz(1:2) = 0.5
       vz(3:4) = 0.57
       vz(5:6) = 0.31
       vz(7:8) = 0.255
       vz(9:10) = 0.165
       vz(11:12) = -0.06
       vz(13:14) = 0.0
       
       By(1:2) = 3.6*const
       By(3:4) = 5.2*const
       By(5:6) = 5.7*const
       By(7:8) = 5.22*const
       By(9:10) = 5.96*const
       By(11:12) = 5.58*const
       By(13:14) = 4.0*const
       
       Bz(1:2) = 2.0*const
       Bz(3:4) = 2.885*const
       Bz(5:6) = 1.76*const
       Bz(7:8) = 1.62*const
       Bz(9:10) = 1.85*const
       Bz(11:12) = 2.79*const
       Bz(13:14) = 2.0*const
       
       pr = rho
       
 END SELECT
!
!--determine which parameter to plot
!
 SELECT CASE(iplot)
    CASE(1)
       yplot = rho
    CASE(2)
       yplot = pr
    CASE(3)
       yplot = vx
    CASE(4)
       yplot = vy
    CASE(5)
       yplot = vz
    CASE(6)
       yplot = By
    CASE(7)
       yplot = Bz
    CASE(8)
       PRINT*,'gamma = ',gamma
       IF (abs(gamma-1.).GT.1.e-5) THEN
          WHERE (abs(rho) > 0.)
             yplot = pr / ((gamma-1.)*rho)
	  END WHERE
       ELSE
          PRINT*,' ***isothermal: utherm solution not valid'
          yplot = 0.	  
       ENDIF
  END SELECT     
!
!--plot exact line using PGPLOT
!
 CALL PGLINE(npts,xplot(1:npts),yplot(1:npts))
            
 RETURN
END SUBROUTINE exact_mhdshock
