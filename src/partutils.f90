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
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
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
 public :: get_binary
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

!----------------------------------------------------------
! routine to get properties of particle binary system
! INPUT:
!   i1, i2 : indexes of two particles to use
!   dat(npart,ncolumns) : particle data
!   ix(ndim) : columns containing positions
!   ivx      : column of first velocity component
!   ipmass   : column containing mass
! OUTPUT:
!   x0 : centre of mass position
!   v0 : velocity of centre of mass
!   angle : angle of binary about centre of mass (radians)
!----------------------------------------------------------
subroutine get_binary(i1,i2,dat,x0,v0,angle,ndim,ndimV,ncolumns,ix,ivx,ipmass,iverbose,ierr)
 integer,                  intent(in)  :: i1,i2,ndim,ndimV,ncolumns,ivx,ipmass,iverbose
 integer, dimension(ndim), intent(in)  :: ix
 real, dimension(:,:),     intent(in)  :: dat
 real, dimension(ndim),    intent(out) :: x0,v0
 real,                     intent(out) :: angle
 integer,                  intent(out) :: ierr
 integer               :: max
 real, dimension(ndim) :: x1,x2,v1,v2,dx
 real                  :: m1,m2,dmtot

 ierr = 0
 max = size(dat(:,1))
 if (i1 <= 0 .or. i2 <= 0 .or. i1 > max .or. i2 > max) then
    if (iverbose >= 2) print*,' star 1 = ',i1,' star 2 = ',i2
    print "(a)",' ERROR locating sink particles in the data'
    ierr = 1
    return
 endif

 x1 = 0.
 x2 = 0.
 x1(1:ndim) = dat(i1,ix(1:ndim))
 x2(1:ndim) = dat(i2,ix(1:ndim))
 !--get centre of mass
 if (ipmass > 0 .and. ipmass <= ncolumns) then
    m1 = dat(i1,ipmass)
    m2 = dat(i2,ipmass)
 else
    m1 = 1.
    m2 = 1.
 endif
 dmtot = 1./(m1 + m2)
 x0 = (m1*x1 + m2*x2)*dmtot
 if (iverbose >= 1) then
    print "(a,3(1x,es10.3),a,es10.3)",' :: star 1 pos =',x1(1:ndim),' m = ',m1
    print "(a,3(1x,es10.3),a,es10.3)",' :: star 2 pos =',x2(1:ndim),' m = ',m2
    print "(a,3(1x,es10.3))",' :: c. of mass =',x0(1:ndim)
 endif
 !--work out angle needed to rotate into corotating frame
 dx   = x0 - x1
 angle = -atan2(dx(2),dx(1))
 
 !--get velocities
 if (ivx > 0 .and. ivx + ndimV <= ncolumns) then
    v1 = dat(i1,ivx:ivx+ndimV-1)
    v2 = dat(i2,ivx:ivx+ndimV-1)
    v0 = (m1*v1 + m2*v2)*dmtot
    if (iverbose >= 1) print "(a,3(1x,es10.3))",' :: vel c of m =',v0(1:ndimV)
 else
    v0 = 0.
 endif

end subroutine get_binary

end module part_utils
