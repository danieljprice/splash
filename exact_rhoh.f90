!------------------------------------
! plot h \propto (pmass/rho)^(1/ndim)
!------------------------------------
subroutine exact_rhoh(hfact,ndim,pmass,npart,xmin,xmax)
 implicit none
 integer, intent(in) :: ndim,npart
 real, intent(in) :: hfact,xmin,xmax
 real, intent(in), dimension(npart) :: pmass
 integer, parameter :: npts = 200
 integer :: i
 real, dimension(0:npts) :: xplot,yplot
 real :: dx,pmassmin,pmassmax
      
 ! hfact = 1.5*massp
 if (hfact.gt.0.01) then
    xplot(0) = max(xmin,1.e-3)
    dx = (xmax-xplot(0))/float(npts)
!
!--plot a line for minimum and maximum masses
!
    pmassmin = minval(pmass)
    pmassmax = maxval(pmass)

    do i=1,npts
      xplot(i) = xplot(0)+dx*(i-1)
      yplot(i) = hfact*(pmassmax/xplot(i))**(1./FLOAT(ndim))
    enddo
    write(*,"(a,f5.2,a,1pe8.2,a,i1,a)") ' plotting h = ',hfact, &
                               '*(',pmassmax,'/rho)**(1/',ndim,')'
    call pgline(npts,xplot(1:npts),yplot(1:npts))
    yplot = 0.
    
    if (abs(pmassmin-pmassmax).gt.1.e-10 .and. pmassmin.gt.1.e-10) then
       do i=1,npts
         yplot(i) = hfact*(pmassmin/xplot(i))**(1./FLOAT(ndim))
       enddo
       write(*,"(a,f5.2,a,1pe8.2,a,i1,a)") ' plotting h = ',hfact, &
                               '*(',pmassmin,'/rho)**(1/',ndim,')'
       call pgline(npts,xplot(1:npts),yplot(1:npts))    
    endif

 endif  
 
 return
end subroutine exact_rhoh
