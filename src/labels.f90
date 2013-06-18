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
!  Copyright (C) 2005-2013 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------
!
!  routines to do with storing and handling of plot labels
!
!-----------------------------------------------------------
module labels
 use params, only:maxplot, maxparttypes
 implicit none
 integer, parameter :: lenlabel = 80
 integer, parameter :: lenunitslabel = 40  ! length of units label
 character(len=lenlabel), dimension(maxplot+2) :: label,labelvec
 character(len=20), dimension(maxparttypes) :: labeltype
 character(len=6), parameter :: labeldefault = 'column'
 integer, dimension(3) :: ix
 integer, dimension(maxplot) :: iamvec
 integer :: ivx,irho,iutherm,ipr,ih,irad,iBfirst,iBpol,iBtor,iax
 integer :: ipmass,ike,ispsound
 integer :: idivb,iJfirst,irhostar
 integer :: iacplane,ipowerspec
 integer :: icv,iradenergy
 integer :: isurfdens,itoomre
 integer :: ipdf,icolpixmap
 integer :: irhorestframe,idustfrac,ideltav

 public

contains

 subroutine reset_columnids
  implicit none
  !
  !--array positions of specific quantities
  !  Identification is used in exact solution
  !  plotting and calculation of additional quantities
  !
  ix(:) = 0
  ivx = 0      ! vx
  irho = 0     ! density
  ipr = 0      ! pressure
  iutherm = 0  ! thermal energy
  ih = 0       ! smoothing length
  irad = 0     ! radius
  ipmass = 0   ! particle mass
  ipr = 0      ! pressure
  irad = 0     ! radius
  ipowerspec = 0 ! power spectrum
  iBfirst = 0  ! Bx
  iax = 0      ! ax (acceleration)
  iBpol = 0    ! B_polx
  iBtor = 0    ! B_torx
  iacplane = 0
  ike = 0
  idivB = 0
  iJfirst = 0
  icv = 0
  iradenergy = 0
  icolpixmap = 0
  irhorestframe = 0

  return
 end subroutine reset_columnids

 logical function is_coord(icol,ndim)
  implicit none
  integer, intent(in) :: icol,ndim
  integer :: i

  is_coord = .false.
  do i=1,ndim
     if (ix(i).eq.icol) is_coord = .true.
  enddo

 end function is_coord

end module labels
