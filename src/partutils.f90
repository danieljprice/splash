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

!-----------------------------------------------------------
!
!  utility routines to do with particle identification
!
!-----------------------------------------------------------
module part_utils
 use params, only:int1
 implicit none

 public :: igettype,get_tracked_particle,get_itrackpart
 public :: locate_nth_particle_of_type
 public :: locate_first_two_of_type
 public :: locate_particle_from_string
 public :: get_binary,got_particles_of_type
 public :: get_positions_of_type,is_trackstring
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
    if (i > ntot .and. i <= ntot1) then
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
integer function get_tracked_particle(string,noftype,iamtype,dat,irho)
 character(len=*), intent(in) :: string
 integer, dimension(:), intent(in) :: noftype
 integer(kind=int1), dimension(:), intent(in) :: iamtype
 integer, intent(in) :: irho
 real, dimension(:,:), intent(in) :: dat
 integer :: ntot,itype,ioffset,ierr

 call get_itrackpart(string,itype,ioffset,ierr)
 if ((itype <= 0 .or. itype > size(noftype)) .and. ioffset > 0) then
    !--type not set, itrackpart = itrackoffset
    get_tracked_particle = ioffset
 elseif (ierr == 0 .and. ioffset > 0) then
    !--want to select nth particle of a particular type
    call locate_nth_particle_of_type(ioffset,get_tracked_particle, &
         itype,iamtype,noftype,ntot)
 elseif (ierr == 0) then
    ntot = sum(noftype)
    get_tracked_particle = locate_particle_from_string(string,ntot,dat,irho)
 else
    get_tracked_particle = 0
 endif

end function get_tracked_particle

!-------------------------------------------------------------------
! routine to find which particle is being tracked, when it is
! given in the form of type:offset
!-------------------------------------------------------------------
subroutine get_itrackpart(string,itracktype,itrackpart,ierr)
 character(len=*), intent(in)  :: string
 integer,          intent(out) :: itracktype,itrackpart,ierr
 integer :: ic

 ic = index(string,':')
 if (ic > 0) then
    read(string(1:ic-1),*,iostat=ierr) itracktype
    read(string(ic+1:),*,iostat=ierr) itrackpart
    if (itrackpart==0) itracktype = 0
 else
    itracktype = 0
    read(string,*,iostat=ierr) itrackpart
    if (itrackpart < 0 .or. ierr /= 0) itrackpart = 0
    ! return ierr = 0 if the string is not a particle id but is still valid
    ! e.g. string='maxdens'
    if (ierr /= 0 .and. is_trackstring(string)) ierr = 0
 endif

end subroutine get_itrackpart

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
 if (size(iamtype(:))==1) then
    i1 = sum(noftype(1:itype-1)) + 1
    i2 = i1 + 1
 else
    i1 = 0
    i2 = 0
    i = 0
    nfound = 0
    do while ((i1==0 .or. i2==0) .and. i <= ntot)
       i = i + 1
       if (iamtype(i)==itype) nfound = nfound + 1
       if (nfound==1) i1 = i
       if (nfound==2) i2 = i
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
 if (size(iamtype(:))==1) then
    ipos = sum(noftype(1:itype-1)) + n
 else
    ipos = 0
    i = 0
    nfound = 0
    do while (ipos==0 .and. i < ntot)
       i = i + 1
       if (iamtype(i)==itype) nfound = nfound + 1
       if (nfound==n) ipos = i
    enddo
 endif

end subroutine locate_nth_particle_of_type

!-------------------------------------------------------------
! locate the particle corresponding to various strings
!  e.g. maxdens = particle of maximum density
!-------------------------------------------------------------
pure integer function locate_particle_from_string(string,ntot,dat,irho) result(ipos)
 character(len=*), intent(in)  :: string
 integer,          intent(in)  :: ntot
 real,             intent(in)  :: dat(:,:)
 integer,          intent(in)  :: irho
 integer :: ipos_tmp(1),ntoti

 ipos = 0
 ipos_tmp = 0
 ntoti = min(ntot,size(dat(:,1))) ! prevent bounds overflow
 select case(string(1:7))
 case('maxdens')
    if (irho > 0 .and. irho <= size(dat(1,:))) ipos_tmp = maxloc(dat(1:ntoti,irho))
    ipos = ipos_tmp(1)
 end select

end function locate_particle_from_string

!-------------------------------------------------------------
! validate strings for function above
!-------------------------------------------------------------
logical function is_trackstring(string)
 character(len=*), intent(in) :: string

 select case(string(1:7))
 case('maxdens')
    is_trackstring = .true.
 case default
    is_trackstring = .false.
 end select

end function is_trackstring

!-------------------------------------------------------------
! check if any particles of type 'mytype' exist
!-------------------------------------------------------------
pure logical function got_particles_of_type(mytype,labeltype,npartoftype)
 character(len=*), intent(in) :: mytype,labeltype(:)
 integer, intent(in) :: npartoftype(:,:)
 integer :: itype

 got_particles_of_type = .false.
 do itype=1,size(npartoftype(:,1))
    if (index(labeltype(itype),mytype) > 0) then
       if (any(npartoftype(itype,:) > 0)) got_particles_of_type = .true.
    endif
 enddo

end function got_particles_of_type

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
!   omega : angular velocity of binary around centre of mass
!----------------------------------------------------------
subroutine get_binary(i1,i2,dat,x0,v0,angle,omega,ndim,ndimV,ncolumns,ix,ivx,ipmass,iverbose,ierr)
 integer,                  intent(in)  :: i1,i2,ndim,ndimV,ncolumns,ivx,ipmass,iverbose
 integer, dimension(ndim), intent(in)  :: ix
 real, dimension(:,:),     intent(in)  :: dat
 real, dimension(ndim),    intent(out) :: x0,v0
 real,                     intent(out) :: angle,omega
 integer,                  intent(out) :: ierr
 integer               :: max
 real, dimension(ndim) :: x1,x2,v1,v2,dx,dv
 real                  :: m1,m2,dmtot,dR2
 real, parameter       :: pi = 4.*atan(1.)

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
    dv = v2 - v0
    dx = x2 - x0
    dR2 = 1./dot_product(dx,dx)
    omega = dv(1)*(-dx(2)*dR2) + dv(2)*(dx(1)*dR2)
    if (iverbose >= 1) print "(a,3(1x,es10.3))",' :: vel c of m =',v0(1:ndimV)
    if (iverbose >= 1) print "(a,1x,es10.3,a,g12.4,a)",' ::      omega =',omega,&
                             ' angle =',-angle*180./pi,' degrees'
 else
    v0 = 0.
    omega = 0.
 endif

end subroutine get_binary

!----------------------------------------------------------
! routine to extract positions of sink particles
! INPUT:
!   dat(npart,ncolumns) : particle data
!   ix(3) : columns containing positions
!   itype : type of each particle
!   ntypes : number of particle types
!   npartoftype : number of particles of each type
! OUTPUT:
!   nsinks : number of sink particles
!   xpts : x positions of sink particles
!   ypts : y positions of sink particles
!   zpts : z positions of sink particles
!   ierr : error code, ierr=0 means no error was encountered
!----------------------------------------------------------
subroutine get_positions_of_type(dat,npartoftype,itype,iget_type,ix,n,xpts,ypts,zpts,ierr)
 integer, intent(in) :: npartoftype(:),ix(3),iget_type
 integer(kind=int1), intent(in) :: itype(:)
 real, intent(in)     :: dat(:,:)
 integer, intent(out) :: n,ierr
 real, intent(out), allocatable :: xpts(:),ypts(:),zpts(:)
 integer :: i,j
 integer, allocatable :: ilist(:)

 ierr = 0
 n = 0
 if (iget_type > 0) n = npartoftype(iget_type)

 allocate(xpts(n),ypts(n),zpts(n),ilist(n),stat=ierr)
 if (ierr /= 0) then
    print*,' ERROR allocating memory for sink particle positions, aborting...'
    return
 endif
 j = 0
 do i=1,sum(npartoftype)
    if (itype(i)==iget_type) then
       j = j + 1
       ilist(j) = i
    endif
 enddo
 if (j /= n) print*,' WARNING: found ',j,' particles but expecting ',n
 if (j < n) n = j

 xpts = dat(ilist(1:n),ix(1))
 ypts = dat(ilist(1:n),ix(2))
 zpts = dat(ilist(1:n),ix(3))

 deallocate(ilist)

end subroutine get_positions_of_type

end module part_utils
