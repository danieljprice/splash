!
! contains subroutines for taking power spectrums on particle data
!
module powerspectrums
 implicit none

contains

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

!----------------------------------------------------------
! Subroutine to compute the power spectrum (periodogram)
! of unevenly sampled data via the Lomb (1976) method
! (algorithm described in Press et al, Numerical Recipes, sec 13.8, p569)
!
! Given the data (dat) on a set of points (x), 
! returns an array of nfreq frequencies (freq) between freqmin and freqmax
! together with the power at each frequency (power)
!----------------------------------------------------------
subroutine powerspectrum_lomb(npts,x,dat,nfreq,freq,freqmin,freqmax,power)
 implicit none
 real, parameter :: pi = 3.1415926536
 integer, intent(in) :: npts
 integer, intent(in) :: nfreq         ! number of frequencies to calculate
 integer :: i,ifreq
 real, intent(in) :: freqmin, freqmax
 real, intent(in), dimension(npts) :: x, dat 
 real, intent(out), dimension(nfreq) :: freq, power
 real :: datmean, datvar         ! mean, variance
 real :: ddat
 real :: tau, tau_numerator, tau_denominator
 real :: term1_numerator, term1_denominator
 real :: term2_numerator, term2_denominator
 real :: omega, omega_dx, cos_term, sin_term
 real :: omegamin,omegamax,domega
 real :: d2pi

 print*,'evaluating lomb periodogram...'
!
!--calculate the mean and variance of the data
!
 call mean_variance(dat,npts,datmean,datvar)
 print*,'data mean = ',datmean,' std. dev = ',sqrt(datvar)
 if (datvar.le.0.) then 
    print*,'error: variance = 0'
    power = 0.
    return
 endif 
!
!--work out range of angular frequencies (wavenumbers)
! 
 omegamin = 2.*pi*freqmin
 omegamax = 2.*pi*freqmax
 domega = (omegamax-omegamin)/nfreq
 d2pi = 1./(2.*pi)
 print*,'omega min = ',omegamin,' max = ',omegamax
!
!--loop over frequencies
!
 omega = omegamin
 
 over_frequency: do ifreq = 1,nfreq
!
!--save frequency (wavenumber) to array
!
    omega = (ifreq)*omegamin
    freq(ifreq) = omega*d2pi
!
!--calculate tau for this frequency
!
    tau_numerator = 0.
    tau_denominator = 0.
    do i=1,npts
       tau_numerator = tau_numerator + SIN(2.*omega*x(i))
       tau_denominator = tau_denominator + COS(2.*omega*x(i))
    enddo
    tau = ATAN(tau_numerator/tau_denominator)/(2.*omega)
!
!--calculate the terms in the power
! 
    term1_numerator = 0.
    term1_denominator = 0.
    term2_numerator = 0.
    term2_denominator = 0.
    do i=1,npts
       ddat = dat(i) - datmean
       omega_dx = omega*(x(i) - tau)
       cos_term = COS(omega_dx)
       sin_term = SIN(omega_dx)
       term1_numerator = term1_numerator + ddat*cos_term
       term1_denominator = term1_denominator + cos_term**2
       term2_numerator = term2_numerator + ddat*sin_term
       term2_denominator = term2_denominator + sin_term**2
    enddo
!
!--calculate the power at this frequency
!    
    power(ifreq) = 1./(2.*datvar)*(term1_numerator**2/term1_denominator + &
                                   term2_numerator**2/term2_denominator)

!    print*,ifreq,' freq = ',freq(ifreq),omega,' power = ', power(ifreq)

!
!--next frequency
!
    omega = omega + domega
 
 enddo over_frequency
 return
end subroutine powerspectrum_lomb

!-------------------------------------------------
! Subroutine to calculate the mean and variance
! of a set of data points
! Mean is trivial but variance uses a special
! formula to reduce round-off error
! see Press et al Numerical Recipes, section 14.2
! this is similar to their subroutine avevar
!-------------------------------------------------
subroutine mean_variance(x,npts,xmean,xvariance)
 implicit none
 integer, intent(in) :: npts
 real, intent(in), dimension(npts) :: x
 real, intent(out) :: xmean, xvariance
 real :: roundoff, delta
 integer :: i
!
!--calculate average
!
 xmean = 0.
 do i=1,npts
    xmean = xmean + x(i)
 enddo
 xmean = xmean/real(npts)
!
!--calculate variance using the corrected two-pass formula
!
!    var = 1/(n-1)*( sum (x-\bar{x}) - 1/n * (sum(x-\bar{x}) )^2 )
!  
!  where the last term corrects for the roundoff error 
!  in the first term
!
 xvariance = 0.
 roundoff = 0.
 
 do i=1,npts
    delta = x(i) - xmean
    roundoff = roundoff + delta
    xvariance = xvariance + delta*delta
 enddo
 xvariance = (xvariance - roundoff**2/npts)/real(npts-1)
 
 return
end subroutine mean_variance

end module powerspectrums
