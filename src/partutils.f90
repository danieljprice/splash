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
!  utility routines to do with particle identification
!
!-----------------------------------------------------------
module part_utils
 implicit none

 public :: igettype,get_tracked_particle
 private

contains
 !
 !--utility returning the type of particle i
 !  when particles are ordered by type
 !
 pure integer function igettype(i,noftype)
  use params, only:maxparttypes
  implicit none
  integer, intent(in) :: i
  integer, dimension(maxparttypes), intent(in) :: noftype
  integer :: ntot,ntot1,jtype

  ntot     = 0
  igettype = 1 ! so even if in error, will not lead to seg fault
  over_types: do jtype=1,maxparttypes
     ntot1 = ntot + noftype(jtype)
     if (i.gt.ntot .and. i.le.ntot1) then
        igettype = jtype
        exit over_types
     endif
     ntot = ntot1
  enddo over_types

 end function igettype
 
 integer function get_tracked_particle(itype,ioffset,noftype,iamtype)
  use params, only:int1,maxparttypes
  implicit none
  integer, intent(in) :: itype,ioffset
  integer, dimension(maxparttypes), intent(in) :: noftype
  integer(kind=int1), dimension(:), intent(in) :: iamtype
  integer :: i,n
  
  if (itype.le.0 .or. itype.gt.size(noftype)) then
     !--type not set, itrackpart = itrackoffset
     get_tracked_particle = ioffset
  else
     !--want to select nth particle of a particular type
     if (size(iamtype(:)).eq.1) then
        get_tracked_particle = sum(noftype(1:itype-1)) + ioffset
     else
        get_tracked_particle = 0
        i = 0
        n = 0
        do while (get_tracked_particle.eq.0 .and. i.lt.size(iamtype))
           i = i + 1
           if (iamtype(i).eq.itype) n = n + 1
           if (n.eq.ioffset) get_tracked_particle = i
        enddo
     endif
  endif

 end function get_tracked_particle

end module part_utils
