!---------------------------------------------------
! Subroutine calls the power spectrum calculation
! and plots the results
!---------------------------------------------------
subroutine plot_powerspectrum(npart,x,dat)
 implicit none
 integer, parameter :: maxfreq = 100
 integer, intent(in) :: npart
 integer :: nfreq
 real, dimension(npart), intent(in) :: x, dat
 real, dimension(maxfreq) :: freq,power
 real :: freqmin,freqmax,powermin,powermax
 real, external :: theoretical_power
!
!--call subroutine to evaluate the power spectrum
!  of the data (dat) on a set of disordered points
!  (the particles), with co-ordinates x
!  returns an array of nfreq frequencies (freq) 
!  between freqmin and freqmax, together
!  with the power at each frequency (power) 
!
 call lomb_powerspectrum(npart,x,dat,nfreq,freq,power,maxfreq)
!
!--work out plot limits
! 
 freqmin = MINVAL(freq(1:nfreq))
 freqmax = MAXVAL(freq(1:nfreq))
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
 call PGLABEL ('f','Power',' 1D Lomb periodogram ')
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
