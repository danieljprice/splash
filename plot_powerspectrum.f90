!---------------------------------------------------
! Subroutine calls the power spectrum calculation
! and plots the results
!---------------------------------------------------
subroutine plot_powerspectrum(npts,x,dat,idisordered)
 implicit none
 integer, parameter :: ioversamplingfactor = 4	! oversample by 4
 integer, parameter :: iabovenyquist = 2	! up to multiple of nyquist frequency
 integer, intent(in) :: npts
 integer :: nfreq, ierr
 real, dimension(npts), intent(in) :: x, dat
 real, dimension(:), allocatable :: freq,power
 real :: xmin,xmax
 real :: freqmin,freqmax,powermin,powermax
 real :: wavelengthmin,wavelengthmax
 real, external :: theoretical_power
 logical, intent(in) :: idisordered
!
!--get max/min of data
! 
 xmin = MINVAL(x)
 xmax = MAXVAL(x)
!
!--set frequency interval between evaluations of power
! 
 nfreq = 32	!128	!ioversamplingfactor*iabovenyquist*npts
!
!--allocate arrays for frequency and power accordingly
! 
 allocate(freq(nfreq),power(nfreq),STAT=ierr)
 if (ierr.ne.0) then
    print*,'error allocating memory for frequency arrays, returning'
    return
 endif
!
!--work out range of frequencies to compute
!
 wavelengthmax = 2.*(xmax-xmin)
 wavelengthmin = wavelengthmax/(REAL(nfreq))	! nyquist frequency
 
 freqmin = 1./wavelengthmax
 freqmax = iabovenyquist/wavelengthmin	! to some multiple of nyquist
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
 power = power/powermax
!
!--set up plotting page
!
 print*,'plotting power spectrum...'
 call PGSWIN(freqmin,freqmax,0.0,1.0,0,1)
 call PGBOX('BCNST',0.0,0,'1BVCNST',0.0,0)      
 call PGLABEL ('frequency','Power',' 1D Power Spectrum ')
!
!--plot power spectrum
!
! call PGLINE(nfreq,freq,power)			! as line
 call PGBIN(nfreq,freq,power,.true.)		! as histogram
!
!--plot theoretical power spectrum
! 
 print*,'plotting theoretical power ',theoretical_power(freqmin)
 call PGSLS(2)	! dashed line
 call PGFUNX(theoretical_power,10000,freqmin,freqmax,1)
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
