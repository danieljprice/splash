!----------------------------------------------------------------------
! sets options and parameters for power spectrum calculation/plotting
!----------------------------------------------------------------------
subroutine options_powerspec
 use settings
 use prompting
 implicit none

 call prompt('enter data to take power spectrum of',ipowerspecy,ndim+1,numplot-nextra) 

 call prompt('use lomb periodogram? (no=interpolate and fourier) ',idisordered)

 return
end subroutine options_powerspec
