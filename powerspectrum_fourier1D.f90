!-------------------------------------------------------
! subroutine to compute the power spectrum 
! of evenly sampled data via a (slow) fourier transform
!--------------------------------------------------------

subroutine powerspectrum_fourier(npts,x,dat,nfreq,freq,freqmin,freqmax,power)
 implicit none
 real, parameter :: pi = 3.1415926536
 integer, intent(in) :: npts, nfreq
 integer :: i,ifreq
 real, intent(in), dimension(npts) :: x, dat
 real, intent(out), dimension(nfreq) :: freq, power
 real, intent(in) :: freqmin, freqmax
 real :: sum1,sum2
 real :: omega, domega, omegamin, omegamax, d2pi

 print*,' evaluating fourier transform'
!
!--work out range of angular frequencies (wavenumbers)
!
 omegamin = 2.*pi*freqmin
 omegamax = 2.*pi*freqmax
 domega = (omegamax-omegamin)/nfreq
 omega = omegamin
 d2pi = 1./(2.*pi)
!
!--now compute (slow!) fourier transform
! 
 do ifreq=1,nfreq
    omega = ifreq*omegamin  ! + domega
    freq(ifreq) = omega*d2pi
    power(ifreq) = 0.
    sum1 = 0.
    sum2 = 0.
    do i=1,npts
       sum1 = sum1 + dat(i)*COS(-omega*x(i)) 
       sum2 = sum2 + dat(i)*SIN(-omega*x(i)) 
    enddo
    power(ifreq) = power(ifreq) + sqrt(sum1**2 + sum2**2)/REAL(npts)
!    print*,ifreq,': freq = ',freq(ifreq),' power = ',power(ifreq)
 enddo
 
 
 
 print*,'done'
 
 return
end subroutine powerspectrum_fourier
