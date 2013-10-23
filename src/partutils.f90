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
 use params, only:int1
 implicit none

 public :: igettype,get_tracked_particle
 public :: locate_nth_particle_of_type
 public :: locate_first_two_of_type
 private

contains

!---------------------------------------------
!  utility returning the type of particle i
!  when particles are ordered by type
!---------------------------------------------
pure integer function igettype(i,noftype)
 use params, only:maxparttypes
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

!-------------------------------------------------------------------
! routine to find which particle is being tracked, when it is
! given in the form of type:offset
!-------------------------------------------------------------------
integer function get_tracked_particle(itype,ioffset,noftype,iamtype)
 use params, only:maxparttypes
 integer, intent(in) :: itype,ioffset
 integer, dimension(:), intent(in) :: noftype
 integer(kind=int1), dimension(:), intent(in) :: iamtype
 integer :: ntot

 if (itype.le.0 .or. itype.gt.size(noftype)) then
    !--type not set, itrackpart = itrackoffset
    get_tracked_particle = ioffset
 else
    !--want to select nth particle of a particular type
    call locate_nth_particle_of_type(ioffset,get_tracked_particle, &
         itype,iamtype,noftype,ntot)
 endif

end function get_tracked_particle

!-------------------------------------------------------------------
! routine to locate first two particles of a given type in the data
!-------------------------------------------------------------------
subroutine locate_first_two_of_type(i1,i2,itype,iamtype,noftype,ntot)
 integer, intent(out) :: i1,i2,ntot
 integer, intent(in)  :: itype
 integer(kind=int1), dimension(:), intent(in) :: iamtype
 integer, dimension(:), intent(in) :: noftype
 integer :: i,nfound
 
 !--locate first two sink particles in the data
 ntot = sum(noftype)
 if (size(iamtype(:)).eq.1) then
    i1 = sum(noftype(1:itype-1)) + 1
    i2 = i1 + 1
 else
    i1 = 0
    i2 = 0
    i = 0
    nfound = 0
    do while ((i1.eq.0 .or. i2.eq.0) .and. i.le.ntot)
       i = i + 1
       if (iamtype(i).eq.itype) nfound = nfound + 1
       if (nfound.eq.1) i1 = i
       if (nfound.eq.2) i2 = i
    enddo
 endif

end subroutine locate_first_two_of_type

!-------------------------------------------------------------
! routine to locate nth particle of a given type in the data
!-------------------------------------------------------------
pure subroutine locate_nth_particle_of_type(n,ipos,itype,iamtype,noftype,ntot)
 integer, intent(out) :: ipos,ntot
 integer, intent(in)  :: n,itype
 integer(kind=int1), dimension(:), intent(in) :: iamtype
 integer, dimension(:), intent(in) :: noftype
 integer :: i,nfound

 ntot = sum(noftype)
 if (size(iamtype(:)).eq.1) then
    ipos = sum(noftype(1:itype-1)) + n
 else
    ipos = 0
    i = 0
    nfound = 0
    do while (ipos.eq.0 .and. i.le.ntot)
       i = i + 1
       if (iamtype(i).eq.itype) nfound = nfound + 1
       if (nfound.eq.n) ipos = i
    enddo
 endif

end subroutine locate_nth_particle_of_type

end module part_utils
