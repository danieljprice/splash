SUBROUTINE setup_integratedkernel
!-------------------------------------------------------------
!     tabulates the integral through the cubic spline kernel
!     subroutine originally by Matthew Bate
!-------------------------------------------------------------
 USE column
 IMPLICIT NONE
 INTEGER :: i,j
 REAL, PARAMETER :: pi = 3.1415926536
 REAL :: r, dist, step, ypos, v, v2, val
 REAL :: coldens, v2m

 PRINT*,'setting up integrated kernel table...'

 DO i=1,maxcoltable
    r=(i-1)/500.
    dist=SQRT(4.0-r*r)
    step=dist/4000.0
    ypos=0.
         
    coldens=0.0
    DO j=1,4000
       v=SQRT(r*r+ypos*ypos)
       IF (v.LT.1.0) THEN
          v2=v*v
          val=1.0-1.5*v2+0.75*v2*v
          coldens=coldens+val*step
       ELSE
          v2m=2.0-v
          val=0.25*v2m*v2m*v2m
          coldens=coldens+val*step        
       ENDIF
       ypos=ypos+step
    END DO
    coltable(i)=2.0*coldens/pi
 END DO

!      check=0.0
!      DO i=1,1000
!         r=(i-1)/500.
!         WRITE(*,*) r, table(i)
!         check=check+2*pi*r*table(i)*(1.0/500.)
!      END DO
!      WRITE(*,*)check

 RETURN
END SUBROUTINE setup_integratedkernel
