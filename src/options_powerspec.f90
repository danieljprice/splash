!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2010 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! Module containing settings and options relating to power spectrum
!  and Probability Distribution Function plots
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_powerspec
 implicit none
 integer :: ipowerspecy, ipowerspecx, nfreqspec
 integer :: nwavelengths,npdfbins
 logical :: idisordered
 real :: freqmax,freqmin

 namelist /powerspecopts/ ipowerspecy,idisordered,nwavelengths,nfreqspec,npdfbins,freqmin,freqmax

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_powerspec
  use settings_data, only:ndim
  implicit none

  idisordered = .true.
  ipowerspecy = max(ndim+1,2)
  ipowerspecx = 0 ! reset later
  nwavelengths = 128
  freqmin = 1.0
  freqmax = nwavelengths*freqmin
  nfreqspec = 1
  npdfbins  = 0

  return
end subroutine defaults_set_powerspec

!----------------------------------------------------------------------
! sets options and parameters for power spectrum calculation/plotting
!----------------------------------------------------------------------
subroutine options_powerspec
 use settings_data, only:ndim,ndataplots,numplot
 use limits,        only:lim
 use labels,        only:ipowerspec
 use prompting,     only:prompt
 implicit none
 real :: boxsize

 if (ipowerspecy.lt.ndim+1) ipowerspecy = ndim+1
 if (ipowerspecy.gt.ndataplots) ipowerspecy = ndataplots
 call prompt('enter data to take power spectrum of',ipowerspecy,ndim+1,ndataplots)
 if (ipowerspecx.ne.1) then
    if (ipowerspecx.lt.1) ipowerspecx = 1
    if (ipowerspecx.gt.ndataplots) ipowerspecx = ndataplots
    call prompt('enter column to use as "time" or "space"',ipowerspecx,1,ndataplots)
 endif
!
!--if box size has not been set then use x limits
!
 if (abs(freqmin-1.0).lt.tiny(1.)) then
    boxsize = abs(lim(1,2) - lim(1,1))
    if (boxsize.gt.tiny(boxsize)) freqmin = 1./boxsize
 endif
 call prompt('enter min frequency (default=1/box size)',freqmin,0.0)
 call prompt('enter max frequency ',freqmax,min=freqmin)

 if (ipowerspec.le.ndataplots .or. ipowerspec.gt.numplot) then
    !--this should never happen
    print*,'*** ERROR: something wrong in powerspectrum limit setting'
 else
    print*,' wavelength range ',1./freqmax,'->',1./freqmin
    lim(ipowerspec,1) = freqmin
    lim(ipowerspec,2) = freqmax
    print*,' frequency range ',lim(ipowerspec,1),'->',lim(ipowerspec,2)
    if (nfreqspec.le.1) nfreqspec = 2*nwavelengths
    call prompt('how many frequency points between these limits? ',nfreqspec,nwavelengths)
 endif

!! call prompt('use Lomb periodogram? (no=interpolate and fourier) ',idisordered)

 return
end subroutine options_powerspec

!-----------------------------------------------------------------
!
! settings for PDF calculation
!
!-----------------------------------------------------------------
subroutine options_pdf
 use prompting, only:prompt
 implicit none

 call prompt(' Enter number of bins between min and max of plot (0=auto)',npdfbins,0)

end subroutine options_pdf

end module settings_powerspec
