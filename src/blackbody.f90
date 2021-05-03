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
!  Copyright (C) 2021- Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
module blackbody
 use physcon, only:c,hplanck,kboltz,cm_to_nm
 use params,  only:doub_prec
 implicit none

contains
!---------------------------------------------------------
! Planck function
! INPUT:
!    temp - temperature [K]
!    nu - frequency [Hz]
! OUTPUT:
!    B_nu - Planck function erg/s/cm^2/Hz/steradian
!---------------------------------------------------------
real elemental function B_nu(temp,nu)
 real, intent(in) :: temp,nu
 real(doub_prec) :: hnu_on_kT,hnu3_on_c2

 hnu_on_kT  = hplanck*nu/(kboltz*temp)
 hnu3_on_c2 = hplanck*nu**3/c**2

 if (hnu_on_kT < 300.) then
    B_nu = 2.*hnu3_on_c2/(exp(hnu_on_kT) - 1.d0)
 else
    B_nu = tiny(0.)
 endif

end function B_nu

!---------------------------------------------------------
! Planck function derivative
! INPUT:
!    temp - temperature [K]
!    nu - frequency [Hz]
! OUTPUT:
!    dB_nu/dT - Planck function derivative w.r.t. T
!---------------------------------------------------------
real elemental function dBnu_dT(temp,nu)
 real, intent(in) :: temp,nu
 real(doub_prec) :: hnu_on_kT,hnu3_on_c2
 real :: Bnu,expterm

 hnu_on_kT  = hplanck*nu/(kboltz*temp)
 hnu3_on_c2 = hplanck*nu**3/c**2

 Bnu = B_nu(temp,nu)
 if (hnu_on_kT < 300.) then
    expterm = exp(hnu_on_kT)
    dBnu_dT = Bnu*expterm/(expterm - 1.)*hnu_on_kT/temp
 else
    dBnu_dT = -tiny(0.)
 endif

end function dBnu_dT

!---------------------------------------------------------
! Wien's displacement law, in frequency
!---------------------------------------------------------
real function Wien_T_from_nu(nu) result(T)
 real, intent(in) :: nu

 T = nu/5.88e10

end function Wien_T_from_nu

!---------------------------------------------------------
! Wien's displacement law, frequency from temperature
!---------------------------------------------------------
real function Wien_nu_from_T(T) result(nu)
 real, intent(in) :: T

 nu = 5.88e10*T

end function Wien_nu_from_T

!---------------------------------------------------------
! Convert frequency in Hz to wavelength in nm
!---------------------------------------------------------
real function nu_to_lam(nu) result(lam)
 real, intent(in) :: nu

 lam = (c/nu)*cm_to_nm

end function nu_to_lam

!---------------------------------------------------------
! Convert frequency in Hz to wavelength in nm
!---------------------------------------------------------
real function lam_to_nu(lam) result(nu)
 real, intent(in) :: lam

 nu = (c/lam)*cm_to_nm

end function lam_to_nu

!---------------------------------------------------------
! colour temperature, from fitting peak of blackbody
!---------------------------------------------------------
subroutine get_colour_temperature(spectrum,freq,Tc,freq_max,bb_scale)
 real, intent(in) :: spectrum(:),freq(:)
 real, intent(out) :: Tc,freq_max,bb_scale
 integer :: imax(1),j

 imax = maxloc(spectrum)
 j = imax(1)
 freq_max = freq(j)
 Tc = Wien_T_from_nu(freq_max)

 ! scaling factor by which to shift Blackbody
 ! in y direction to match the peak flux
 bb_scale = spectrum(j)/B_nu(Tc,freq_max)

end subroutine get_colour_temperature

!--------------------------------------------------------
!+
!  Fit blackbody to spectrum over a frequency range
!  by least squares minimisation
!+
!--------------------------------------------------------
subroutine bb_fit(spec,freq,fmin,fmax,T,bbscale,min,max)
 real, intent(in) :: fmin,fmax
 real, intent(in) :: freq(:),spec(:)
 real, intent(inout) :: T
 real, intent(out) :: bbscale
 real, intent(in), optional :: min,max
 real :: Tmax,Tmin,func
 integer :: its

 Tmin = 1e3
 Tmax = 1e6
 if (present(min)) Tmin = min
 if (present(max)) Tmax = max
 print*,'initial min/max = ',Tmin,Tmax,0.5*(Tmin + Tmax)

 ! rootfind, using bisection method
 its = 0
 func = huge(func)
 do while(abs(func) > 1.e-20 .and. its < 100)
    T = 0.5*(Tmin + Tmax)
    func = dchisq_dT(spec,freq,fmin,fmax,T,bbscale)
    if (func > 0.) then
       Tmax = T
    else
       Tmin = T
    endif
    its = its + 1
    !print*,its,func,'T = ',T,' min/max = ',Tmin,Tmax,bbscale
 enddo
 print*,its,func,'T = ',T,' min/max = ',Tmin,Tmax,bbscale

end subroutine bb_fit

!--------------------------------------------------------
!+
!  derivative of the chisquared error between the
!  input data (spec) and the blackbody spectrum at T
!  this is the function which we find the root of
!+
!--------------------------------------------------------
real function dchisq_dT(spec,freq,fmin,fmax,T,bbscale)
 real, intent(in) :: fmin,fmax
 real, intent(in) :: freq(:),spec(:)
 real, intent(in) :: T
 real, intent(out) :: bbscale
 integer :: i,iloc(1)
 real :: err,Bnu

 ! fit bb_scale to match maximum value of the spectrum
 ! in the specified frequency range
 iloc = maxloc(spec,mask=(freq > fmin .and. freq < fmax))
 i = max(iloc(1),1)
 !print*,'matching at temperature ',i,Wien_T_from_nu(freq(i)),' spec = ', &
!        spec(i),' b_nu = ',B_nu(T,freq(i)),' scale=',bbscale

 bbscale = 1.e-5
 if (B_nu(T,freq(i)) > epsilon(0.)) bbscale = spec(i)/B_nu(T,freq(i))
 bbscale = max(bbscale,1e-5)

 dchisq_dT = 0.
 do i=1,size(freq)
    if (freq(i) > fmin .and. freq(i) < fmax) then
       Bnu = B_nu(T,freq(i))
       err = spec(i)/bbscale - Bnu
       dchisq_dT = dchisq_dT - 2.*err/dBnu_dT(T,freq(i))
    endif
 enddo

end function dchisq_dT

!--------------------------------------------------------
!+
!  Function to fill an array with equally log-spaced points
!+
!--------------------------------------------------------
function logspace(n,xmin,xmax) result(x)
 integer, intent(in) :: n
 real, intent(in)    :: xmin,xmax
 real, allocatable   :: x(:)
 integer :: i
 real    :: dx

 allocate(x(n))
 dx = log10(xmax/xmin)/real(n-1)
 do i=1,n
    x(i) = log10(xmin) + (i-1)*dx
 enddo

 x = 10.**x

end function logspace

!--------------------------------------------------------
!+
!  Integrate function on evenly spaced logarithmic grid
!  i.e. \int f(x) dx = \int x f(x) d(ln x)
!  using trapezoidal rule
!+
!--------------------------------------------------------
real function integrate_log(f,x,xmin,xmax) result(fint)
 real, intent(in) :: xmin,xmax
 real, intent(in) :: x(:),f(:)
 real :: dlogx
 integer :: n,i

 n = size(f)
 dlogx = log(xmax/xmin)/(n-1)
 fint = 0.
 do i=2,n
    if (x(i-1) > xmin .and. x(i) < xmax) then
       fint = fint + 0.5*(f(i)*x(i) + f(i-1)*x(i-1))*dlogx
    endif
 enddo

end function integrate_log

end module blackbody
