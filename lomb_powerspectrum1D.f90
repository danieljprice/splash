!----------------------------------------------------------
! Subroutine to compute the power spectrum (periodogram)
! of unevenly sampled data via the Lomb (1976) method
! (see Press et al, Numerical Recipes, sec 13.8, p569)
!
! Given the data (dat) on a set of points (x), 
! returns an array of frequencies (freq) together
! with the power at each frequency (power)
!----------------------------------------------------------
subroutine lomb_powerspectrum(npts,x,dat,nfreq,freq,power,maxfreq)
 implicit none
 real, parameter :: pi = 3.1415926536
 integer, parameter :: ioversamplingfactor = 4	! oversample by 4
 integer, parameter :: iabovenyquist = 2	! up to multiple of nyquist frequency
 integer, intent(in) :: npts
 integer, intent(in) :: maxfreq	! dimensions of freq & power arrays
 integer, intent(out) :: nfreq	! actual number of freq's calculated
 integer :: i,ifreq
 real, intent(in), dimension(npts) :: x, dat 
 real, intent(out), dimension(maxfreq) :: freq, power
 real :: datmean, datvar	! mean, variance
 real :: xmin, xmax, dx, ddat
 real :: tau, tau_numerator, tau_denominator
 real :: term1_numerator, term1_denominator
 real :: term2_numerator, term2_denominator
 real :: omega, omega_dx, cos_term, sin_term
 real :: wavelengthmin,wavelengthmax,omegamin,omegamax,domega

 print*,'evaluating lomb periodogram...'
!
!--calculate the mean and variance of the data
!
 call mean_variance(dat,npts,datmean,datvar)
!
!--get max/min of data
! 
 xmin = MINVAL(x)
 xmax = MAXVAL(x)
!
!--set frequency interval between evaluations of power
! 
 nfreq = ioversamplingfactor*iabovenyquist*npts
 if (nfreq.gt.maxfreq) then
    print*,'warning: # of frequencies > array size, setting nfreq=',maxfreq
    nfreq = maxfreq
 endif
!
!--work out minimum (angular) frequency (max wavelength = xmax-xmin)
!
 wavelengthmax = 2.*(xmax-xmin)
 wavelengthmin = wavelengthmax/(REAL(nfreq))	! nyquist frequency
 
 omegamin = 2.*pi/wavelengthmax
 omegamax = iabovenyquist*2.*pi/wavelengthmin	! to some multiple of nyquist

 domega = (omegamax-omegamin)/REAL(nfreq)
!
!--loop over frequencies
!
 omega = omegamin
 
 over_frequency: do ifreq = 1,nfreq
!
!--save frequency (wavenumber) to array
!
    freq(ifreq) = omega
    print*,'freq(',ifreq,') = ',omega
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
!
!--next frequency
!
    omega = omega + domega
 
 enddo over_frequency
 return
end subroutine lomb_powerspectrum

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
 integer :: i,j
!
!--calculate average
!
 xmean = 0.
 do i=1,npts
    xmean = xmean + x(i)
 enddo
 xmean = xmean/npts
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
 xvariance = (xvariance - roundoff**2/npts)/(npts-1)
 
 return
end subroutine mean_variance
