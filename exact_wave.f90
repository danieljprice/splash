! ----------------------------------------------------------------------
! compute exact solution for a linear wave
! plots a sine function with a given amplitude, period and wavelength
! ----------------------------------------------------------------------

subroutine exact_wave(time,ampl,period,lambda,xmin,xmax,ymean)
  implicit none
  integer, parameter :: npts = 100
  integer :: i
  real, parameter :: pi = 3.1415926536
  real, intent(in) :: time, ampl, period, lambda, xmin, xmax, ymean 
  real, dimension(npts) :: xplot,yplot
  real :: omega, dx

  print*,'plotting sine wave... mean = ',ymean
  print*,' lambda = ',lambda,' ampl = ',ampl,' period = ',period
!
! check for errors
!
  if (lambda.le.0.) then
     print*,'error: lambda <= 0'
     return
  endif
  if (period.gt.0.) then
     omega = 2.*pi/period
  else
     print*,'warning: period <= 0'
     omega = 0.
  endif

  dx = (xmax-xmin)/float(npts-1)
  do i=1,npts
     xplot(i) = xmin+(i-1)*dx
     yplot(i) = ymean + ampl*sin(2.*pi/lambda*(xplot(i)-xmin) - omega*time)
  enddo
  
  call PGLINE(npts,xplot,yplot)
  
  return
end subroutine exact_wave
