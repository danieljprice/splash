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
!!     &           rho_L,rho_R,pr_L,pr_R,vx_L,vx_R,xmin,xmax)
 IMPLICIT NONE
 INTEGER, PARAMETER :: maxpts=16
 INTEGER, INTENT(IN) :: iplot,ishk
 REAL, INTENT(in) :: time,gamma,xmin,xmax
!! REAL, INTENT(in) :: rho_L,rho_R,pr_L,pr_R,vx_L,vx_R
 REAL, DIMENSION(maxpts) :: xplot,yplot
 REAL, DIMENSION(maxpts) :: rho,pr,vx,vy,vz,By,Bz
 REAL :: dx,const
 REAL :: machno,vs
 INTEGER :: i,j,npts
      
 PRINT*,' Plotting exact mhd shock #',ishk,' at t = ',time
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

    CASE(2)
!
!--fast/slow shock from RJ95
!    
       vz = 0.
       Bz = 0.
       npts = 12
       xplot(1) = xmin
       xplot(2) = -0.27
       xplot(3) = -0.09
       xplot(4) = -0.03
       xplot(5) = -0.01
       xplot(6:7) = 0.135
       xplot(8:9) = 0.25
       xplot(10:11) = 0.35
       xplot(12) = xmax

       rho(1:2) = 1.0
       rho(3:4) = 0.5955
       rho(5:6) = 0.55151
       rho(7:8) = 0.41272
       rho(9:10) = 0.2337
       rho(11:12) = 0.2
       
       pr(1:2) = 1.0
       pr(3:4) = 0.42629
       pr(5:8) = 0.37090
       pr(9:10) = 0.12402
       pr(11:12) = 0.1
       
       vx(1:2) = 0.0
       vx(3:4) = 0.81237
       vx(5:8) = 0.89416
       vx(9:10) = 0.24722
       vx(11:12) = 0.0
       
       vy(1:2) = 0.0
       vy(3:4) = -0.59961
       vy(5:8) = -0.5447
       vy(9:10) = -0.91164
       vy(11:12) = 0.
       
       By(1:2) = 1.0
       By(3:4) = 0.28431
       By(5:8) = 0.31528
       By(9:10) = 0.43086
       By(11:12) = 0.0

    CASE(3)
!
!--problem with 7 discontinuities from RJ95
!    
       npts = 16
       xplot(1) = xmin
       xplot(2:3) = -0.19
       xplot(4:5) = 0.03
       xplot(6:7) = 0.051
       xplot(8:9) = 0.12	! contact discontinuity
       xplot(10:11) = 0.18 
       xplot(12:13) = 0.205 
       xplot(14:15) = 0.45
       xplot(16) = xmax
              
       rho(1:2) = 1.08
       rho(3:4) = 1.4903
       rho(5:8) = 1.6343
       rho(9:10) = 1.4735
       rho(11:14) = 1.3090
       rho(15:16) = 1.0
              
       pr(1:2) = 0.95
       pr(3:4) = 1.6558
       pr(5:10) = 1.9317
       pr(11:14) = 1.5844
       pr(15:16) = 1.0
       
       vx(1:2) = 1.2 
       vx(3:5) = 0.60588
       vx(6:11) = 0.57538
       vx(12:14) = 0.53432
       vx(15:16) = 0.0
       
       vy(1:2) = 0.01
       vy(3:4) = 0.11235
       vy(5:6) = 0.22157
       vy(7:8) = 0.047602
       vy(9:10) = 0.047601
       vy(11:12) = -0.18411
       vy(13:14) = -0.094572
       vy(15:16) = 0.0
       
       vz(1:2) = 0.5
       vz(3:4) = 0.55686
       vz(5:6) = 0.30125
       vz(7:10) = 0.24734
       vz(11:12) = 0.17554
       vz(13:14) = -0.047286
       vz(15:16) = 0.0
       
       By(1:2) = 1.0155
       By(3:4) = 1.4383
       By(5:6) = 1.5716
       By(7:10) = 1.4126
       By(11:12) = 1.6103
       By(13:14) = 1.5078
       By(15:16) = 1.1284
       
       Bz(1:2) = 0.56419
       Bz(3:4) = 0.79907
       Bz(5:6) = 0.48702
       Bz(7:10) = 0.43772
       Bz(11:12) = 0.49899
       Bz(13:14) = 0.75392
       Bz(15:16) = 0.56419

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
 
    CASE(5)
!
!--rarefaction from RJ95
!
       npts = 6
       vy = 0.
       vz = 0.
       Bz = 0.
       xplot(1) = xmin
       xplot(2) = -0.27
       xplot(3) = -0.12
       xplot(4) = 0.12
       xplot(5) = 0.27
       xplot(6) = xmax
       
       rho(1:2) = 1.0
       rho(3:4) = 0.49653
       rho(5:6) = 1.0

       pr(1:2) = 1.0
       pr(3:4) = 0.31134
       pr(5:6) = 1.0
       
       vx(1:2) = -1.0
       vx(3:4) = 0.	! this is approximate (to 10-7)
       vx(5:6) = 1.0
       
       By(1:2) = 1.0
       By(3:4) = 0.49638
       By(5:6) = 1.0
       
    CASE(6)
!
!--mach 25 shocks from Dai and Woodward (1994)
!
       npts = 6
       rho(1:2) = 1.0
       rho(3:4) = 3.982
       rho(5:6) = 0.1
       
       pr(1:2) = 1.0
       pr(3:4) = 1806.0
       pr(5:6) = 1.0

!       machno = 0.5*25.5       
!       vs = SQRT(gamma*pr(1)/rho(1))
!
!      in this case we know the positions at all times
!      because of the Mach no's
!
       xplot(1) = xmin
       xplot(2:3) = -0.35	!-machno*vs*time
       xplot(4:5) = 0.35	!machno*vs*time
       xplot(6) = xmax
       
!       PRINT*,'speed, pos = ',vs*machno,vs,machno*vs*time,0.33/time
       
       vx(1:2) = 36.87
       vx(3:4) = 0.0
       vx(5:6) = -36.87
       
       vy(1:2) = -0.1546
       vy(3:4) = -0.07727
       vy(5:6) = 0.0
       
       vz(1:2) = -0.03864
       vz(3:4) = -0.01932
       vz(5:6) = 0.0
       
       By(1:2) = 4.0*const
       By(3:4) = 15.95*const
       By(5:6) = 4.0*const
        
       Bz(1:2) = 1.0*const
       Bz(3:4) = 3.988*const
       Bz(5:6) = 1.0*const
       
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
