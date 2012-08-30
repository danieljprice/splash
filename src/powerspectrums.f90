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
!  Copyright (C) 2005-2012 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!
! contains subroutines for taking power spectrums on particle data
!
module powerspectrums
 implicit none
 real, parameter, private :: pi = 3.141592653589
 real, parameter, private :: twopi = 2.*pi
 public :: powerspectrum,powerspec3D_sph

 private

contains

subroutine powerspectrum(npts,x,dat,nfreqpts,freq,power,idisordered)
 implicit none
 integer, intent(in) :: npts, nfreqpts
 real, intent(in), dimension(npts) :: x
 real, intent(in), dimension(npts) :: dat
 real, intent(in), dimension(nfreqpts) :: freq
 real, intent(out), dimension(nfreqpts) :: power
 logical, intent(in) :: idisordered
 integer :: ifreq
 real :: datmean, datvar, omega

 if (.not.idisordered) then
    print*,' evaluating fourier transform'
    do ifreq=1,nfreqpts
       omega = twopi*freq(ifreq)
       !--get power at this frequency
       call power_fourier(npts,x,dat,omega,power(ifreq))
    enddo
 else
    print*,'evaluating lomb periodogram...'
!
!--calculate the mean and variance of the data
!
    call mean_variance(dat,npts,datmean,datvar)
    print*,'data mean = ',datmean,' std. dev = ',sqrt(datvar)
    if (datvar.le.0.) then
       print*,'error: variance = 0'
       power = 0.
       return
    endif
    do ifreq=1,nfreqpts
       omega = twopi*freq(ifreq)
       call power_lomb(npts,x,dat,datmean,datvar,omega,power(ifreq))
    enddo

 endif

end subroutine powerspectrum

!-------------------------------------------------------
! subroutine to compute the power spectrum
! of evenly sampled data via a (slow) fourier transform
!--------------------------------------------------------

subroutine power_fourier(npts,x,dat,omega,power)
 implicit none
 integer, intent(in) :: npts
 real, intent(in), dimension(npts) :: x, dat
 real, intent(in) :: omega
 real, intent(out) :: power
 integer :: i
 real :: sum1,sum2

 power = 0.
 sum1 = 0.
 sum2 = 0.
 do i=1,npts
    sum1 = sum1 + dat(i)*COS(-omega*x(i))
    sum2 = sum2 + dat(i)*SIN(-omega*x(i))
 enddo
 power= sqrt(sum1**2 + sum2**2)/REAL(npts)

 return
end subroutine power_fourier

!----------------------------------------------------------
! Subroutine to compute the power spectrum (periodogram)
! of unevenly sampled data via the Lomb (1976) method
! (algorithm described in Press et al, Numerical Recipes, sec 13.8, p569)
!
! Given the data (dat) on a set of points (x),
! returns an array of nfreq frequencies (freq) between freqmin and freqmax
! together with the power at each frequency (power)
!----------------------------------------------------------
subroutine power_lomb(npts,x,dat,datmean,datvar,omega,power)
 implicit none
 integer, intent(in) :: npts
 real, intent(in), dimension(npts) :: x, dat
 real, intent(in) :: datmean,datvar,omega
 real, intent(out) :: power
 integer :: i
 real :: ddat
 real :: tau, tau_numerator, tau_denominator
 real :: term1_numerator, term1_denominator
 real :: term2_numerator, term2_denominator
 real :: omega_dx, cos_term, sin_term
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
 power = 1./(2.*datvar)*(term1_numerator**2/term1_denominator + &
                         term2_numerator**2/term2_denominator)

 return
end subroutine power_lomb

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
 integer :: i
!
!--calculate average
!
 xmean = 0.
 do i=1,npts
    xmean = xmean + x(i)
 enddo
 xmean = xmean/real(npts)
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
 xvariance = (xvariance - roundoff**2/npts)/real(npts-1)

 return
end subroutine mean_variance

!
!  interface to 3D powerspectrum calculation on particles
!  assumes box size is the same in all directions
!
subroutine powerspec3D_sph(x,y,z,dat,hh,weight,icolours,npart, &
 ngrid,xmin,xmax,freq,power,normalise)
 use interpolations3D, only:interpolate3D
 implicit none
 integer, intent(in) :: npart,ngrid
 real, dimension(npart), intent(in) :: x,y,z,dat,hh,weight
 integer, dimension(npart), intent(in) :: icolours
 real, intent(in) :: xmin,xmax
 real, dimension(ngrid), intent(out) :: freq,power
 logical, intent(in) :: normalise
 real, dimension(ngrid,ngrid,ngrid) :: dat3D
 real :: dx
 integer :: logngrid,ik
 logical :: periodicx,periodicy,periodicz
!
!--make sure than ngrid is a factor of 2
!
 logngrid = int(log(real(ngrid))/log(2.))
 if (2**logngrid.ne.ngrid) then
    print*,' ERROR: ngrid not a power of 2 in powerspectrum interpolation ',2**logngrid,ngrid
 endif

 dx = (xmax - xmin)/real(ngrid)
!
!--interpolate (normalised) from particles to 3D grid suitable for FFT
!
 print*,'ngrid = ',ngrid
 periodicx = .false.
 periodicy = .false.
 periodicz = .false.
 call interpolate3D(x,y,z,hh,weight,dat,icolours,npart, &
      xmin,xmin,xmin,dat3D,ngrid,ngrid,ngrid,dx,dx,normalise,&
      periodicx,periodicy,periodicz)
!
!--setup grid of frequencies for plotting
!
 freq(1) = 0.
 do ik=2,ngrid
    freq(ik) = ik - 1.
 enddo
!
!--calculate powerspectrum using fft
!
 call power3d_fft(dat3D,ngrid,ngrid,ngrid,power,ngrid)

 return
end subroutine powerspec3D_sph

!
!--power spectrum routine using Fast Fourier Transform
!
subroutine power3d_fft(dat,nx,ny,nz,power,nk)
 implicit none
! include 'fftw3.f'
 integer, intent(in) :: nx,ny,nz,nk
 real, intent(in),  dimension(nx,ny,nz) :: dat
 real, intent(out), dimension(nk) :: power
 integer, dimension(nk) :: numk
 complex :: dati(nx,ny,nz)
 real :: ddenom,ptot
 integer :: ierr,k,j,i,kz,ky,kx,kk,ik
!--this is for ACML
! complex :: comm(nx*ny*nz+5*(nx+ny+nz))
!--this is for FFTW
! integer*8 :: plan

 ierr = 0
!
!--convert data to complex
!
 dati = cmplx(dat,0.0)

 print*,' starting 3D fft...'
!
!--do fast fourier transform via AMD Core Math Library function
!
! call cfft3d(-1,nx,ny,nz,dati,comm,ierr)
!
!--do fft via fftw
!
! call fftwf_plan_dft_3d(plan,nx,ny,nz,dati,dati,FFTW_FORWARD,FFTW_ESTIMATE)
! call fftwf_execute(plan)
! call fftwf_destroy_plan(plan)

 if (ierr /= 0) then
    write(*,*) 'error on powerspectrum output!'
 endif
 power = 0.
 numk = 0
!
!--get power from fourier coefficients
!
 do k=1,nz
    kz = min(k-1,nz-k+1)
    do j=1,ny
       ky = min(j-1,ny-j+1)
       do i=1,nx
          kx = min(i-1,nx-i+1)
          kk = sqrt(real(kx**2 + ky**2 + kz**2))
          ik = 1.5 + kk
          !--only return requested number of frequencies
          if (ik .le. nk) then
             power(ik) = power(ik) + abs(dati(i,j,k))**2
             numk(ik) = numk(ik) + 1 ! sum contributions at that frequency
          endif
       enddo
    enddo
 enddo

 ddenom = 1./(real(nx)*real(ny)*real(nz))
 power = power*ddenom**2
 ptot =  sum(power)
!
!--normalise according to power in k-space "shells"
!  and number of contributions in that shell from kx,ky and kz
!
 do ik=1,nk
    power(ik) = power(ik)/(numk(ik) + 1.e-8)*4./3.*pi*((ik-0.5)**3 - (ik-1.5)**3)
 enddo
!--rescale so that it has the same total power as before
 power = power*ptot/sum(power)

 return
end subroutine power3d_fft



end module powerspectrums
