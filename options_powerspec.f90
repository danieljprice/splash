!-------------------------------------------------------------------------
! Module containing settings and options relating to power spectrum plots
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_powerspec
 implicit none
 integer :: ipowerspecy, nfreqspec
 logical :: idisordered
 real :: wavelengthmax
 
 namelist /powerspecopts/ ipowerspecy,idisordered,wavelengthmax,nfreqspec

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_powerspec
  use settings_data, only:ndim
  implicit none

  idisordered = .false.
  ipowerspecy = ndim+1
  wavelengthmax = 1.0
  nfreqspec = 32

  return
end subroutine defaults_set_powerspec

!----------------------------------------------------------------------
! sets options and parameters for power spectrum calculation/plotting
!----------------------------------------------------------------------
subroutine options_powerspec
 use settings_data ! for ndim, numplot
 use prompting
 implicit none

 call prompt('enter data to take power spectrum of',ipowerspecy,ndim+1,numplot-nextra) 

 call prompt('enter box size (max wavelength)',wavelengthmax,0.0)
 
 call prompt('enter number of frequencies to sample ',nfreqspec,1)

 call prompt('use Lomb periodogram? (no=interpolate and fourier) ',idisordered)

 return
end subroutine options_powerspec

end module settings_powerspec
