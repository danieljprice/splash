!-------------------------------------------------------------------------
! Module containing settings and options relating to power spectrum plots
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_powerspec
 implicit none
 integer :: ipowerspecy, ipowerspecx, nfreqspec
 logical :: idisordered
 real :: wavelengthmax, oversamplefactor
 
 namelist /powerspecopts/ ipowerspecy,idisordered,wavelengthmax,nfreqspec, &
                          oversamplefactor

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_powerspec
  use settings_data, only:ndim
  implicit none

  idisordered = .true.
  ipowerspecy = ndim+1
  ipowerspecx = ndim
  wavelengthmax = 1.0 ! reset later
  nfreqspec = 32
  oversamplefactor = 2.0

  return
end subroutine defaults_set_powerspec

!----------------------------------------------------------------------
! sets options and parameters for power spectrum calculation/plotting
!----------------------------------------------------------------------
subroutine options_powerspec
 use settings_data, only:ndim,ndataplots
 use limits, only:lim
 use prompting
 implicit none
 real :: boxsize

 call prompt('enter data to take power spectrum of',ipowerspecy,ndim+1,ndataplots)
 if (ipowerspecx.ne.1) then
    ipowerspecx = 1
    call prompt('enter column to use as "time" or "space"',ipowerspecx,1,ndataplots) 
 endif
!
!--if box size has not been set then use x limits
!
 if (abs(wavelengthmax-1.0).lt.tiny(wavelengthmax)) then
    boxsize = abs(lim(1,2) - lim(1,1))
    if (boxsize.gt.tiny(boxsize)) wavelengthmax = boxsize
 endif
 call prompt('enter box size (max wavelength)',wavelengthmax,0.0)
 
 call prompt('enter number of frequencies to sample ',nfreqspec,1)
 
 call prompt('enter oversampling factor ',oversamplefactor,1.0)

!! call prompt('use Lomb periodogram? (no=interpolate and fourier) ',idisordered)

 return
end subroutine options_powerspec

end module settings_powerspec
