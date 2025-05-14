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
!  Copyright (C) 2005-2023 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------
!----------------------------------------------------------------------------
!
!  modules containing physical constants
!
!----------------------------------------------------------------------------
module physcon
 use params, only:doub_prec
 implicit none
 real, parameter :: pi = 4.*atan(1.)
 real(doub_prec), parameter :: solarrcgs = 6.955d10  ! cm
 real(doub_prec), parameter :: solarmcgs = 1.989d33  ! g
 real(doub_prec), parameter :: steboltz = 5.67e-5    ! erg cm^-2 K^-4 s-1
 real(doub_prec), parameter :: radconst = 7.5646d-15 ! Radiation constant erg cm^-3 K^-4
 real(doub_prec), parameter :: kboltz = 1.38066d-16
 real(doub_prec), parameter :: mh = 1.67262158d-24   ! g
 real(doub_prec), parameter :: au = 1.496d13         ! cm
 real(doub_prec), parameter :: c = 2.997924d10       ! Speed of light cm/s
 real(doub_prec), parameter :: hplanck =   6.6260755d-27  ! Planck's Constant erg/s
 real(doub_prec), parameter :: kb_on_mh = kboltz/mh
 real(doub_prec), parameter :: Lsun = 3.839d33       ! Solar luminosity, erg/s
 real(doub_prec), parameter :: cm_to_nm = 1.d7
 real(doub_prec), parameter :: keV_to_erg = 1.6022d-9 ! k_eV to erg
 real(doub_prec), parameter :: keV_to_Hz  = keV_to_erg/hplanck ! k_eV to erg
 real(doub_prec), parameter :: mass_electron_cgs = 9.10938291d-28  !Electron mass in g
 real(doub_prec), parameter :: qe = 4.8032068d-10     !charge on electron        esu
 real(doub_prec), parameter :: sigma_e = 6.652e-25    ! Thomson cross section

 public
end module physcon
