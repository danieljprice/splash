C     ------------
CN    NAME: R I E M A N N 
C     ------------

CP    PURPOSE:
CP    THIS PROGRAM COMPUTES THE SOLUTION OF A 1D   
CP    RELATIVISTIC RIEMANN PROBLEM (FOR CONSTANT-GAMMA IDEAL GASES) WITH  
CP    INITIAL DATA UL IF X<0.5 AND UR IF X>0.5  
CP    IN THE WHOLE SPATIAL DOMAIN [0, 1] 
C

CC    COMMENTS:
CC    SEE MARTI AND MUELLER, JFM, 1994
CC
CC    WRITTEN BY:     Jose-Maria Marti
CC                    Departamento de Astronomia y Astrofisica 
CC                    Universidad de Valencia 
CC                    46100 Burjassot (Valencia), Spain
CC                    jose-maria.marti@uv.es
CC    AND
CC                    Ewald Mueller
CC                    Max-Planck-Institut fuer Astrophysik
CC                    Karl-Schwarzschild-Str. 1
CC                    85741 Garching, Germany
CC                    emueller@mpa-garching.mpg.de
C

      SUBROUTINE RIEMANN(MN,RAD,RHOA,PA,VELA,UA,RHOLIN,RHORIN,PLIN,PRIN,
     &                   VLIN,VRIN,GAMIN,TIN,X0)
      IMPLICIT NONE

C     ------- 
C     COMMON BLOCKS 
C     -------

      DOUBLE PRECISION RHOL, PL, UL, HL, CSL, VELL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, WR 
      COMMON /STATES/  RHOL, PL, UL, HL, CSL, VELL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, WR

      DOUBLE PRECISION RHOLS, ULS, HLS, CSLS, VELLS, VSHOCKL 
      COMMON /LS/      RHOLS, ULS, HLS, CSLS, VELLS, VSHOCKL

      DOUBLE PRECISION RHORS, URS, HRS, CSRS, VELRS, VSHOCKR 
      COMMON /RS/      RHORS, URS, HRS, CSRS, VELRS, VSHOCKR

      DOUBLE PRECISION GAMMA 
      COMMON /ADIND/   GAMMA

C     --------- 
C     INTERNAL VARIABLES 
C     ---------

      INTEGER         MN, N, I, ILOOP 
C      PARAMETER      (MN = 400)

      DOUBLE PRECISION TOL, PMIN, PMAX, DVEL1, DVEL2, CHECK
      DOUBLE PRECISION RHOLIN,RHORIN,PLIN,PRIN,VLIN,VRIN,GAMIN,TIN,XMIN

      DOUBLE PRECISION PS, VELS

      DOUBLE PRECISION RHOA(MN), PA(MN), VELA(MN), UA(MN)

      DOUBLE PRECISION XI, X0

      DOUBLE PRECISION RAD(MN), X1, X2, X3, X4, X5, T

C     ------- 
C     INITIAL STATES 
C     -------

C      WRITE(*,*) ' ADIABATIC INDEX OF THE GAS: ' 
C     READ (*,*)   GAMMA
      GAMMA = GAMIN

C      WRITE(*,*) ' TIME FOR THE SOLUTION: ' 
C      READ (*,*)   T
      T = TIN
      XMIN = X0 - 0.5D0

C     -----
C     LEFT STATE
C     -----

C      WRITE(*,*) ' -- LEFT STATE -- ' 
C      WRITE(*,*) '      PRESSURE     : ' 
C      READ (*,*) PL 
C      WRITE(*,*) '      DENSITY      : ' 
C      READ (*,*) RHOL 
C      WRITE(*,*) '      FLOW VELOCITY: ' 
C      READ (*,*) VELL
      PL = PLIN
      RHOL = RHOLIN
      VELL = VLIN
      
      PR = PRIN
      RHOR = RHORIN
      VELR = VRIN

C     ------ 
C     RIGHT STATE 
C     ------

C      WRITE(*,*) ' -- RIGHT STATE -- ' 
C      WRITE(*,*) '      PRESSURE     : ' 
C      READ (*,*) PR 
C      WRITE(*,*) '      DENSITY      : ' 
C      READ (*,*) RHOR 
C      WRITE(*,*) '      FLOW VELOCITY: ' 
C      READ (*,*) VELR

C     ------------------------------ 
C     SPECIFIC INTERNAL ENERGY, SPECIFIC ENTHALPY, SOUND SPEED AND  
C     FLOW LORENTZ FACTORS IN THE INITIAL STATES 
C     ------------------------------ 

      UL  = PL/(GAMMA-1.D0)/RHOL 
      UR  = PR/(GAMMA-1.D0)/RHOR

      HL  = 1.D0+UL+PL/RHOL 
      HR  = 1.D0+UR+PR/RHOR

      CSL = DSQRT(GAMMA*PL/RHOL/HL) 
      CSR = DSQRT(GAMMA*PR/RHOR/HR)

      WL  = 1.D0/DSQRT(1.D0-VELL**2) 
      WR  = 1.D0/DSQRT(1.D0-VELR**2)

C     -------- 
C     NUMBER OF POINTS 
C     --------

      N   = MN

C     ------------- 
C     TOLERANCE FOR THE SOLUTION 
C     -------------

      TOL = 0.D0

C

      ILOOP = 0

      PMIN  = (PL + PR)/2.D0 
      PMAX  = PMIN

 5    ILOOP = ILOOP + 1

      PMIN  = 0.5D0*MAX(PMIN,0.D0) 
      PMAX  = 2.D0*PMAX

      CALL GETDVEL(PMIN, DVEL1)

      CALL GETDVEL(PMAX, DVEL2)

      CHECK = DVEL1*DVEL2 
      IF (CHECK.GT.0.D0) GOTO 5

C     --------------------------- 
C     PRESSURE AND FLOW VELOCITY IN THE INTERMEDIATE STATES 
C     ---------------------------

      CALL GETP(PMIN, PMAX, TOL, PS)

      VELS = 0.5D0*(VELLS + VELRS)

C     --------------- 
C     SOLUTION ON THE NUMERICAL MESH 
C     ---------------

C     ----------- 
C     POSITIONS OF THE WAVES 
C     -----------

      IF (PL.GE.PS) THEN

        X1 = X0 + (VELL - CSL )/(1.D0 - VELL*CSL )*T 
        X2 = X0 + (VELS - CSLS)/(1.D0 - VELS*CSLS)*T

      ELSE

        X1 = X0 + VSHOCKL*T 
        X2 = X1

      END IF

      X3 = X0 + VELS*T

      IF (PR.GE.PS) THEN

        X4 = X0 + (VELS + CSRS)/(1.D0 + VELS*CSRS)*T 
        X5 = X0 + (VELR + CSR )/(1.D0 + VELR*CSR )*T

      ELSE

        X4 = X0 + VSHOCKR*T 
        X5 = X4

      END IF

C     ---------- 
C     SOLUTION ON THE MESH 
C     ----------

      DO 100 I=1,N

        RAD(I) = XMIN + DFLOAT(I-1)/DFLOAT(N-1)

 100  CONTINUE

      DO 120 I=1,N

        IF (RAD(I).LE.X1) THEN

          PA(I)   = PL 
          RHOA(I) = RHOL 
          VELA(I) = VELL 
          UA(I)   = UL

        ELSE IF (RAD(I).LE.X2) THEN

          XI = (RAD(I) - X0)/T

          CALL RAREF(XI, RHOL,    PL,    UL,    CSL,  VELL,  'L', 
     &                   RHOA(I), PA(I), UA(I),       VELA(I))

        ELSE IF (RAD(I).LE.X3) THEN

          PA(I)   = PS 
          RHOA(I) = RHOLS 
          VELA(I) = VELS 
          UA(I)   = ULS

        ELSE IF (RAD(I).LE.X4) THEN

          PA(I)   = PS 
          RHOA(I) = RHORS 
          VELA(I) = VELS 
          UA(I)   = URS

        ELSE IF (RAD(I).LE.X5) THEN

          XI = (RAD(I) - X0)/T

          CALL RAREF(XI, RHOR,    PR,    UR,    CSR,  VELR,  'R', 
     &                   RHOA(I), PA(I), UA(I),       VELA(I))

        ELSE

          PA(I)   = PR 
          RHOA(I) = RHOR 
          VELA(I) = VELR 
          UA(I)   = UR

        END IF

 120  CONTINUE

C      OPEN (3,FILE='solution.dat',FORM='FORMATTED',STATUS='NEW')

C      WRITE(3,150) N, T 
C 150  FORMAT(I5,1X,F10.5)

C      DO 60 I=1,N 
C        WRITE(3,200) RAD(I),PA(I),RHOA(I),VELA(I),UA(I) 
C 60   CONTINUE

C 200  FORMAT(5(E15.8,1X))

C      CLOSE(3)

      RETURN
      END

C     ---------- 
CN    NAME: G E T D V E L 
C     ----------

CP    PURPOSE: 
CP    COMPUTE THE DIFFERENCE IN FLOW SPEED BETWEEN LEFT AND RIGHT INTERMEDIATE 
CP    STATES FOR GIVEN LEFT AND RIGHT STATES AND PRESSURE 
C 

CC    COMMENTS
CC    NONE

      SUBROUTINE GETDVEL( P, DVEL )

      IMPLICIT NONE

C     ----- 
C     ARGUMENTS 
C     -----

      DOUBLEPRECISION P, DVEL

C     ------- 
C     COMMON BLOCKS 
C     -------

      DOUBLE PRECISION RHOLS,ULS,HLS,CSLS,VELLS,VSHOCKL 
      COMMON /LS/      RHOLS,ULS,HLS,CSLS,VELLS,VSHOCKL

      DOUBLE PRECISION RHORS,URS,HRS,CSRS,VELRS,VSHOCKR 
      COMMON /RS/      RHORS,URS,HRS,CSRS,VELRS,VSHOCKR

      DOUBLE PRECISION RHOL, PL, UL, HL, CSL, VELL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, WR 
      COMMON /STATES/  RHOL, PL, UL, HL, CSL, VELL, WL, 
     &                 RHOR, PR, UR, HR, CSR, VELR, WR

      DOUBLE PRECISION GAMMA 
      COMMON /ADIND/   GAMMA

C     ----- 
C     LEFT WAVE 
C     -----

      CALL GETVEL(P, RHOL, PL, UL,  HL,  CSL,  VELL,  WL, 'L', 
     &               RHOLS,    ULS, HLS, CSLS, VELLS, VSHOCKL )

C     ----- 
C     RIGHT WAVE 
C     -----

      CALL GETVEL(P, RHOR, PR, UR,  HR,  CSR,  VELR,  WR, 'R', 
     &               RHORS,    URS, HRS, CSRS, VELRS, VSHOCKR )

      DVEL = VELLS - VELRS

      RETURN 
      END

C     ------- 
CN    NAME: G E T P 
C     -------

CP    PURPOSE: 
CP    FIND THE PRESSURE IN THE INTERMEDIATE STATE OF A RIEMANN PROBLEM IN 
CP    RELATIVISTIC HYDRODYNAMICS 
C 

CC    COMMENTS: 
CC    THIS ROUTINE USES A COMBINATION OF INTERVAL BISECTION AND INVERSE 
CC    QUADRATIC INTERPOLATION TO FIND THE ROOT IN A SPECIFIED INTERVAL. 
CC    IT IS ASSUMED THAT DVEL(PMIN) AND DVEL(PMAX) HAVE OPPOSITE SIGNS WITHOUT 
CC    A CHECK. 
CC    ADAPTED FROM "COMPUTER METHODS FOR MATHEMATICAL COMPUTATION", 
CC    BY G. E. FORSYTHE, M. A. MALCOLM, AND C. B. MOLER, 
CC    PRENTICE-HALL, ENGLEWOOD CLIFFS N.J. 
C
      SUBROUTINE GETP( PMIN, PMAX, TOL, PS )

      IMPLICIT NONE

C     ----- 
C     ARGUMENTS 
C     -----

      DOUBLEPRECISION PMIN, PMAX, TOL, PS

C     ------- 
C     COMMON BLOCKS 
C     -------

      DOUBLEPRECISION GAMMA 
      COMMON /ADIND/  GAMMA

      DOUBLEPRECISION RHOL, PL, UL, HL, CSL, VELL, WL, 
     &                RHOR, PR, UR, HR, CSR, VELR, WR 
      COMMON /STATES/ RHOL, PL, UL, HL, CSL, VELL, WL, 
     &                RHOR, PR, UR, HR, CSR, VELR, WR

C     --------- 
C     INTERNAL VARIABLES 
C     ---------

      DOUBLEPRECISION A, B, C, D, E, EPS, FA, FB, FC, TOL1, 
     & XM, P, Q, R, S

C     ------------- 
C     COMPUTE MACHINE PRECISION 
C     -------------

      EPS  = 1.D0 
10    EPS  = EPS/2.D0 
      TOL1 = 1.D0 + EPS 
      IF( TOL1 .GT. 1.D0 ) GO TO 10

C     ------- 
C     INITIALIZATION 
C     -------

      A  = PMIN 
      B  = PMAX 
      CALL GETDVEL(A,FA) 
      CALL GETDVEL(B,FB)

C     ----- 
C     BEGIN STEP 
C     -----

20    C  = A 
      FC = FA 
      D  = B - A 
      E  = D 
30    IF( DABS(FC) .GE. DABS(FB) )GO TO 40 
      A  = B 
      B  = C 
      C  = A 
      FA = FB 
      FB = FC 
      FC = FA

C     -------- 
C     CONVERGENCE TEST 
C     --------

40    TOL1 = 2.D0*EPS*DABS(B) + 0.5D0*TOL 
      XM   = 0.5D0*(C - B) 
      IF( DABS(XM) .LE. TOL1 ) GO TO 90 
      IF( FB .EQ. 0.D0 ) GO TO 90

C     ------------ 
C     IS BISECTION NECESSARY? 
C     ------------

      IF( DABS(E) .LT. TOL1 ) GO TO 70 
      IF( DABS(FA) .LE. DABS(FB) ) GO TO 70

C     ------------------ 
C     IS QUADRATIC INTERPOLATION POSSIBLE? 
C     ------------------

      IF( A .NE. C ) GO TO 50

C     ---------- 
C     LINEAR INTERPOLATION 
C     ----------

      S = FB/FA 
      P = 2.D0*XM*S 
      Q = 1.D0 - S 
      GO TO 60

C     ---------------- 
C     INVERSE QUADRATIC INTERPOLATION 
C     ----------------

50    Q = FA/FC 
      R = FB/FC 
      S = FB/FA 
      P = S*(2.D0*XM*Q*(Q - R) - (B - A)*(R - 1.D0)) 
      Q = (Q - 1.D0)*(R - 1.D0)*(S - 1.D0)

C     ------ 
C     ADJUST SIGNS 
C     ------

60    IF( P .GT. 0.D0 ) Q = -Q 
      P = DABS(P)

C     -------------- 
C     IS INTERPOLATION ACCEPTABLE? 
C     --------------

      IF( (2.D0*P) .GE. (3.D0*XM*Q-DABS(TOL1*Q)) ) GO TO 70 
      IF( P .GE. DABS(0.5D0*E*Q) ) GO TO 70 
      E = D 
      D = P/Q 
      GO TO 80

C     ----- 
C     BISECTION 
C     -----

70    D = XM 
      E = D

C     ------- 
C     COMPLETE STEP 
C     -------

80    A  = B 
      FA = FB 
      IF( DABS(D) .GT. TOL1 ) B = B+D 
      IF( DABS(D) .LE. TOL1 ) B = B+DSIGN(TOL1,XM) 
      CALL GETDVEL(B,FB) 
      IF( (FB*(FC/DABS(FC))) .GT. 0.D0) GO TO 20 
      GO TO 30

C     -- 
C     DONE 
C     --

90    PS = B

      RETURN 
      END

C     --------- 
CN    NAME: G E T V E L 
C     ---------

CP    PURPOSE: 
CP    COMPUTE THE FLOW VELOCITY BEHIND A RAREFACTION OR SHOCK IN TERMS OF THE 
CP    POST-WAVE PRESSURE FOR A GIVEN STATE AHEAD THE WAVE IN A RELATIVISTIC 
CP    FLOW 
C 

CC    COMMENTS: 
CC    THIS ROUTINE CLOSELY FOLLOWS THE EXPRESSIONS IN MARTI AND MUELLER, 
CC    J. FLUID MECH., (1994)

      SUBROUTINE GETVEL( P, RHOA, PA, UA, HA, CSA, VELA, WA, S, 
     & RHO,      U,  H,  CS,  VEL,  VSHOCK )

      IMPLICIT NONE

C     ----- 
C     ARGUMENTS 
C     -----

      DOUBLE PRECISION P, RHOA, PA, UA, HA, CSA, VELA, WA  
      CHARACTER*1      S 
      DOUBLE PRECISION RHO, U, H, CS, VEL, VSHOCK

C     ------- 
C     COMMON BLOCKS 
C     -------

      DOUBLE PRECISION GAMMA 
      COMMON /ADIND/   GAMMA

C     --------- 
C     INTERNAL VARIABLES 
C     ---------

      DOUBLE PRECISION A, B, C, SIGN 
      DOUBLE PRECISION J, WSHOCK 
      DOUBLE PRECISION K, SQGL1

C     --------------- 
C     LEFT OR RIGHT PROPAGATING WAVE 
C     ---------------

      IF (S.EQ.'L') SIGN = -1.D0

      IF (S.EQ.'R') SIGN =  1.D0

C

      IF (P.GT.PA) THEN

C       --- 
C       SHOCK 
C       ---

        A  = 1.D0+(GAMMA-1.D0)*(PA-P)/GAMMA/P 
        B  = 1.D0-A 
        C  = HA*(PA-P)/RHOA-HA**2

C       ---------------- 
C       CHECK FOR UNPHYSICAL ENTHALPIES 
C       ----------------

        IF (C.GT.(B**2/4.D0/A)) STOP 
     & 'GETVEL: UNPHYSICAL SPECIFIC ENTHALPY IN INTERMEDIATE STATE'

C       ----------------------------- 
C       SPECIFIC ENTHALPY IN THE POST-WAVE STATE 
C       (FROM THE EQUATION OF STATE AND THE TAUB ADIABAT, 
C       EQ.(74), MM94) 
C       -----------------------------

        H = (-B+DSQRT(B**2-4.D0*A*C))/2.D0/A

C       --------------- 
C       DENSITY IN THE POST-WAVE STATE 
C       (FROM EQ.(73), MM94) 
C       ---------------

        RHO = GAMMA*P/(GAMMA-1.D0)/(H-1.D0)

C       ------------------------ 
C       SPECIFIC INTERNAL ENERGY IN THE POST-WAVE STATE 
C       (FROM THE EQUATION OF STATE) 
C       ------------------------

        U = P/(GAMMA-1.D0)/RHO

C       -------------------------- 
C       MASS FLUX ACROSS THE WAVE  
C       (FROM THE RANKINE-HUGONIOT RELATIONS, EQ.(71), MM94) 
C       --------------------------

        J = SIGN*DSQRT((P-PA)/(HA/RHOA-H/RHO))

C       ---------- 
C       SHOCK VELOCITY 
C       (FROM EQ.(86), MM94 
C       ----------

        A      = J**2+(RHOA*WA)**2 
        B      = -VELA*RHOA**2*WA**2 
        VSHOCK = (-B+SIGN*J**2*DSQRT(1.D0+RHOA**2/J**2))/A 
        WSHOCK = 1.D0/DSQRT(1.D0-VSHOCK**2)

C       ------------------- 
C       FLOW VELOCITY IN THE POST-SHOCK STATE 
C       (FROM EQ.(67), MM94) 
C       -------------------

        A = WSHOCK*(P-PA)/J+HA*WA*VELA 
        B = HA*WA+(P-PA)*(WSHOCK*VELA/J+1.D0/RHOA/WA)

        VEL = A/B

C       --------------------- 
C       LOCAL SOUND SPEED IN THE POST-SHOCK STATE 
C       (FROM THE EQUATION OF STATE) 
C       ---------------------

        CS = DSQRT(GAMMA*P/RHO/H)

      ELSE

C       ------ 
C       RAREFACTION 
C       ------

C       --------------------------- 
C       POLITROPIC CONSTANT OF THE GAS ACROSS THE RAREFACTION 
C       ---------------------------

        K = PA/RHOA**GAMMA

C       --------------- 
C       DENSITY BEHIND THE RAREFACTION 
C       ---------------

        RHO = (P/K)**(1.D0/GAMMA)

C       ------------------------ 
C       SPECIFIC INTERNAL ENERGY BEHIND THE RAREFACTION 
C       (FROM THE EQUATION OF STATE) 
C       ------------------------

        U = P/(GAMMA-1.D0)/RHO

C       -------------------- 
C       LOCAL SOUND SPEED BEHIND THE RAREFACTION 
C       (FROM THE EQUATION OF STATE) 
C       --------------------

        CS = DSQRT(GAMMA*P/(RHO+GAMMA*P/(GAMMA-1.D0)))

C       ------------------ 
C       FLOW VELOCITY BEHIND THE RAREFACTION 
C       ------------------

        SQGL1 = DSQRT(GAMMA-1.D0) 
        A   = (1.D0+VELA)/(1.D0-VELA)* 
     & ((SQGL1+CSA)/(SQGL1-CSA)* 
     & (SQGL1-CS )/(SQGL1+CS ))**(-SIGN*2.D0/SQGL1)

        VEL = (A-1.D0)/(A+1.D0)

      END IF

      END

C     -------- 
CN    NAME: R A R E F 
C     --------

CP    PURPOSE: 
CP    COMPUTE THE FLOW STATE IN A RAREFACTION FOR GIVEN PRE-WAVE STATE 
C
 
CC    COMMENTS: 
CC    THIS ROUTINE CLOSELY FOLLOWS THE EXPRESSIONS IN MARTI AND MUELLER, 
CC    J. FLUID MECH., (1994)

      SUBROUTINE RAREF( XI, RHOA, PA, UA, CSA, VELA, S, RHO, P, U, VEL )

      IMPLICIT NONE

C     ----- 
C     ARGUMENTS 
C     -----

      DOUBLE PRECISION XI

      DOUBLE PRECISION RHOA, PA, UA, CSA, VELA

      CHARACTER        S

      DOUBLE PRECISION RHO, P, U, VEL

C     ------- 
C     COMMON BLOCKS 
C     -------

      DOUBLE PRECISION GAMMA 
      COMMON /ADIND/   GAMMA

C     --------- 
C     INTERNAL VARIABLES 
C     ---------

      DOUBLE PRECISION B, C, D, K, L, V, OCS2, FCS2, DFDCS2, CS2, SIGN

C     --------------- 
C     LEFT OR RIGHT PROPAGATING WAVE 
C     ---------------

      IF (S.EQ.'L') SIGN =  1.D0

      IF (S.EQ.'R') SIGN = -1.D0

      B    = DSQRT(GAMMA - 1.D0) 
      C    = (B + CSA)/(B - CSA) 
      D    = -SIGN*B/2.D0 
      K    = (1.D0 + XI)/(1.D0 - XI) 
      L    = C*K**D 
      V    = ((1.D0 - VELA)/(1.D0 + VELA))**D

      OCS2 = CSA

25    FCS2   = L*V*(1.D0 + SIGN*OCS2)**D*(OCS2 - B) + 
     & (1.D0 - SIGN*OCS2)**D*(OCS2 + B)

      DFDCS2 = L*V*(1.D0 + SIGN*OCS2)**D* 
     & (1.D0 + SIGN*D*(OCS2 - B)/(1.D0 + SIGN*OCS2)) + 
     & (1.D0 - SIGN*OCS2)**D*
     & (1.D0 - SIGN*D*(OCS2 + B)/(1.D0 - SIGN*OCS2))

      CS2 = OCS2 - FCS2/DFDCS2

      IF (ABS(CS2 - OCS2)/OCS2.GT.5.E-7)THEN 
        OCS2 = CS2 
        GOTO 25 
      END IF

      VEL = (XI + SIGN*CS2)/(1.D0 + SIGN*XI*CS2)

      RHO = RHOA*((CS2**2*(GAMMA - 1.D0 - CSA**2))/ 
     & (CSA**2*(GAMMA - 1.D0 - CS2**2))) 
     & **(1.D0/(GAMMA - 1.D0))

      P   = CS2**2*(GAMMA - 1.D0)*RHO/(GAMMA - 1.D0 - CS2**2)/GAMMA

      U   = P/(GAMMA - 1.D0)/RHO

      RETURN  
      END 
