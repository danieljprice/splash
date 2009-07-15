! ----------------------------------------------------------------------
! compute exact solution for a linear wave
! plots a sine function with a given amplitude, period and wavelength
! ----------------------------------------------------------------------
module wave
  implicit none

contains

subroutine exact_wave(time,ampl,period,lambda,x0,ymean,xplot,yplot,ierr)
  implicit none
  integer :: i
  real, parameter :: pi = 3.1415926536
  real, intent(in) :: time, ampl, period, lambda, x0, ymean 
  real, intent(in), dimension(:) :: xplot
  real, intent(out), dimension(size(xplot)) :: yplot
  integer, intent(out) :: ierr
  real :: omega

  print*,'plotting sine wave... mean = ',ymean
  print*,' lambda = ',lambda,' ampl = ',ampl,' period = ',period
!
! check for errors
!
  ierr = 0
  if (lambda.le.0.) then
     print*,'error: lambda <= 0'
     ierr = 1
     return
  endif
  if (abs(period).gt.tiny(period)) then
     omega = 2.*pi/period
  else
     print*,'warning: period <= 0'
     omega = 0.
  endif

  do i=1,size(xplot)
     if (abs(ymean).le.0.) then
        yplot(i) = ymean + ampl*sin(2.*pi/lambda*(xplot(i)-x0) - omega*time)
     else
        yplot(i) = ymean*(1. + ampl*sin(2.*pi/lambda*(xplot(i)-x0) - omega*time))
     endif
  enddo
    
  return
end subroutine exact_wave

end module wave
