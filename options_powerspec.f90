!----------------------------------------------------------------------
! sets options and parameters for power spectrum calculation/plotting
!----------------------------------------------------------------------
subroutine options_powerspec
 use settings
 use prompting
 implicit none

 call prompt('enter data to take power spectrum of',ipowerspecy,ndim+1,numplot-nextra) 

 call prompt('enter box size (max wavelength)',wavelengthmax,0.0)
 
 call prompt('enter number of frequencies to sample ',nfreqspec,1)

 call prompt('use Lomb periodogram? (no=interpolate and fourier) ',idisordered)

 return
end subroutine options_powerspec
