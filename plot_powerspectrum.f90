!---------------------------------------------------
! Subroutine calls the power spectrum calculation
! and plots the results
!---------------------------------------------------
subroutine plot_powerspectrum(npts,nfreq,xlength,x,dat,idisordered,itrans)
 implicit none
 integer, parameter :: ioversamplingfactor = 4	! oversample by 4
 integer, intent(in) :: npts,nfreq,itrans
 integer :: ierr, ifreq
 real, dimension(npts), intent(in) :: x, dat
 real, dimension(nfreq) :: freq,freqplot,power
 real :: xlength
 real :: freqmin,freqmax,freqminplot,freqmaxplot,powermin,powermax
 real :: wavelengthmin,wavelengthmax,zero
 real, external :: theoretical_power
 logical, intent(in) :: idisordered

 zero = 1.e-10
!
!--get max/min of data
! 
! xmin = MINVAL(x)
! xmax = MAXVAL(x)
!
!--set frequency interval between evaluations of power
! 
! nfreq = 32*iabovenyquist	!128	!ioversamplingfactor*iabovenyquist*npts
 if (nfreq.lt.1) then
    print*,'error: nfreq = ',nfreq
    return
 endif
!
!--work out range of frequencies to compute
!
 if (xlength.le.0.) then
    print*,'error: max wavelength = ',xlength
    return 
 endif
 wavelengthmax = xlength
 wavelengthmin = wavelengthmax/(REAL(nfreq))	! nyquist frequency
 
 freqmin = 1./wavelengthmax
 freqmax = nfreq*freqmin   !1./wavelengthmin	! to some multiple of nyquist
 print*,'wavelengths from lambda = ',wavelengthmin,' to ',wavelengthmax
 print*,'frequencies range f = ',freqmin,' to ',freqmax

 if (.not.idisordered) then
!
!--if data points are evenly distributed in x
!  (e.g. after interpolation), then take the fourier series
!
!  call to power spectrum returns an array of nfreq frequencies (freq) 
!  between freqmin and freqmax, together with the power at each frequency (power) 
!
    call powerspectrum_fourier(npts,x,dat,nfreq,freq,freqmin,freqmax,power)
!
!  or evaluate the power spectrum (periodogram) on a set of disordered points
!  via the Lomb algorithm
!
 else
    call powerspectrum_lomb(npts,x,dat,nfreq,freq,freqmin,freqmax,power)
 endif
!
!--work out plot limits
! 
! freqmin = MINVAL(freq(1:nfreq))
! freqmax = MAXVAL(freq(1:nfreq))
 powermin = MINVAL(power(1:nfreq))
 powermax = MAXVAL(power(1:nfreq))
!
!--normalise power spectrum to 1
!
 if (powermax.gt.zero) then
    power = power/powermax
 else
    powermax = 0.  ! avoid numerical problems
 endif
!
!--take logarithms if appropriate
!
 if (itrans.eq.1) then
    freqminplot = LOG10(freqmin)
    freqmaxplot = LOG10(freqmax)
    where (freq > 0.)
       freqplot = LOG10(freq)
    end where
    where (power > 0.)
       power = LOG10(power)
    end where
    powermin = MINVAL(power(1:nfreq))
    powermax = MAXVAL(power(1:nfreq))   
 else	
    freqminplot = freqmin
    freqmaxplot = freqmax
    freqplot = freq    
 endif

!
!--set up plotting page
!
 print*,'plotting power spectrum...',freqminplot,freqmaxplot,powermin,powermax
 call PGSWIN(freqminplot,freqmaxplot,min(powermin,0.0),powermax,0,1)
 if (itrans.eq.1) then
    call PGBOX('BCNSTL',0.0,0,'1BVCNSTL',0.0,0)      
 else
    call PGBOX('BCNST',0.0,0,'1BVCNST',0.0,0)      
 endif
 call PGLABEL ('frequency','Power',' 1D Power Spectrum ')
!
!--plot power spectrum
!
 call PGLINE(nfreq,freqplot,power)			! as line
! call PGBIN(nfreq,freq,power,.true.)		! as histogram
!
!--plot theoretical power spectrum
! 
 print*,'plotting theoretical power '
 call PGSLS(2)	! dashed line
 power = 0.
 do ifreq = 1,nfreq
    if (freq(ifreq).gt.zero) then   ! in case there has been an error in powerspec
       if (itrans.eq.1) then
          if (theoretical_power(freq(ifreq)).gt.zero) then
	     !!print*,theoretical_power(freq(ifreq))
             power(ifreq) = log10(theoretical_power(freq(ifreq)))
          endif  
       else
          power(ifreq) = theoretical_power(freq(ifreq))
       endif
    endif
 enddo

! power = power/power(1)
 call PGLINE(nfreq,freqplot,power)
    
! if (itrans.eq.1) then
!    call PGFUNX(theoretical_power,10000,freqmin,freqmax,1) 
! else
!    call PGFUNX(theoretical_power,10000,freqmin,freqmax,1)
! endif
 call PGSLS(1)
 
 return
end subroutine plot_powerspectrum
!
!--function giving the theoretical power spectrum
!
function theoretical_power(k)
 implicit none
 integer :: nindex
 real :: k,theoretical_power
 
 nindex = -2
 theoretical_power = k**nindex
 
end function theoretical_power
