!------------------------------------
! plot h \propto (pmass/rho)^(1/ndim)
!------------------------------------
subroutine exact_rhoh(hfact,ndim,pmass,npart,xmin,xmax)
 implicit none
 integer, parameter :: npts = 200
 integer, intent(in) :: ndim,npart
 integer i,j
 real, dimension(0:npts) :: xplot,yplot
 real, intent(in) :: hfact,xmin,xmax
 real, intent(in), dimension(npart) :: pmass
 real dx,pmassprev
      
 ! hfact = 1.5*massp
 if (hfact.gt.0.01) then
    write(*,"(a,f5.2,a,i1,a)") 'plotting relation h = ',hfact,'*(m/rho)**(1/',ndim,')'
    dx = (xmax-xmin)/float(npts)
    xplot(0) = xmin
!
!--plot one line for each different mass value
!
    pmassprev = 0.
    do j=1,npart
       if (abs(pmass(j)-pmassprev).gt.1.e-9) then
          do i=1,npts
             xplot(i) = xplot(0)+dx*(i-1)
             yplot(i) = hfact*(pmass(j)/xplot(i))**(1./FLOAT(ndim))
             !  print*,i,' x,y = ',xplot(i),yplot(i)
          enddo
          !!print*,' plotting h = ',hfact,'*(',pmass(j),'/rho)^(1/ndim)'
          call pgline(npts+1,xplot,yplot)
       endif
       pmassprev = pmass(j)
    enddo
 endif  
 
 return
end subroutine exact_rhoh
