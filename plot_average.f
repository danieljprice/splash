!-----------------------------
! plot average line

      SUBROUTINE plot_average(xdata,ydata,npart,nbins)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: npart,nbins
      REAL, DIMENSION(npart), INTENT(IN) :: xdata,ydata
      REAL, DIMENSION(nbins) :: xplot,ymean
      REAL dxbin,xmin,xmax
      INTEGER, DIMENSION(nbins) :: ihoc,num
      INTEGER, DIMENSION(npart) :: ll
      INTEGER ibin,i,j
      
      xmax = MAXVAL(xdata) + 0.0001
      xmin = MINVAL(xdata) - 0.0001

      dxbin = (xmax-xmin)/FLOAT(nbins)
      
      DO ibin=1,nbins
         ihoc(ibin) = -1
         ymean(ibin) = 0.0
         num(ibin) = 0
      ENDDO

      DO i=1,npart
c
cc----bin particles in x to plot average density
c      
         ibin = INT((xdata(i)-xmin)/dxbin) + 1
         ll(i) = ihoc(ibin)
         ihoc(ibin) = i
         num(ibin) = num(ibin) + 1
      ENDDO

      DO ibin=1,nbins
         PRINT*,'num(',ibin,') = ',num(ibin)
         j = ihoc(ibin)
         DO i=1,num(ibin)
            ymean(ibin) = ymean(ibin) + ydata(j)
            j = ll(j)
         ENDDO
         IF (num(ibin).ne.0) THEN
            ymean(ibin) = ymean(ibin)/float(num(ibin))
         ELSE
            ymean(ibin) = 0.0
         ENDIF
         PRINT*,'ymean = ',ymean(ibin)
      ENDDO

      DO i=1,nbins
         xplot(i) = xmin+(i-1)*dxbin + 0.5*dxbin
      ENDDO

      CALL PGPT(nbins,xplot,ymean,5)
      CALL PGSLS(2)
      CALL PGLINE(nbins,xplot,ymean)
      CALL PGSLS(1)            
      
      RETURN
      END SUBROUTINE plot_average
