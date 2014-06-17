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
 character(len=lenunitslabel), dimension(0:maxplot), public :: unitslabel
 character(len=lenunitslabel), public :: labelzintegration
 integer, dimension(3)       :: ix
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

!--------------------------------------------------------------
!
!  utility to reset default settings for column identification
!
!--------------------------------------------------------------
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

!--------------------------------------------------------------
!
!  query function for whether column is a spatial coordinate
!
!--------------------------------------------------------------
logical function is_coord(icol,ndim)
 implicit none
 integer, intent(in) :: icol,ndim
 integer :: i

 is_coord = .false.
 do i=1,ndim
    if (ix(i).eq.icol) is_coord = .true.
 enddo

end function is_coord

!-----------------------------------------------------------------
!
!  utility to strip spaces, escape sequences and
!  units labels from strings (this can be called for both
!  function strings and variable labels)
!
!-----------------------------------------------------------------
elemental function shortstring(string,unitslab)
 use asciiutils, only:string_delete
 implicit none
 character(len=lenlabel), intent(in)           :: string
 character(len=*),        intent(in), optional :: unitslab
 character(len=lenlabel)                       :: shortstring
 integer :: ipos

 shortstring = string
 !--strip off the units label
 if (present(unitslab)) then
    if (len_trim(unitslab).gt.0) then
    !--remove units label (only do this once)
       ipos = index(trim(shortstring),trim(unitslab))
       if (ipos.ne.0) then
          shortstring = shortstring(1:ipos-1)//&
                        shortstring(ipos+len_trim(unitslab)+1:len_trim(shortstring))
       endif
    endif
 endif

 !--remove spaces
 call string_delete(shortstring,' ')
 !--remove escape sequences (\d etc.)
 call string_delete(shortstring,'\d')
 call string_delete(shortstring,'\u')
 call string_delete(shortstring,'\g')
 call string_delete(shortstring,'\')
 call string_delete(shortstring,'_')

end function shortstring

!------------------------------------------------------------------
!
! Same as shortstring, but also strips any arithmetic operators
! should be applied to variable names, but not function strings
!
!-----------------------------------------------------------------
elemental function shortlabel(string,unitslab)
 use asciiutils, only:string_delete
 implicit none
 character(len=lenlabel), intent(in)           :: string
 character(len=*),        intent(in), optional :: unitslab
 character(len=lenlabel)                       :: shortlabel

 if (present(unitslab)) then
    shortlabel = shortstring(string,unitslab)
 else
    shortlabel = shortstring(string)
 endif
 !--remove arithmetic operators from labels
 call string_delete(shortlabel,'**')
 call string_delete(shortlabel,'/')
 call string_delete(shortlabel,'*')
 call string_delete(shortlabel,'+')
 call string_delete(shortlabel,'-')
 call string_delete(shortlabel,'^')
 call string_delete(shortlabel,'sqrt(')
 call string_delete(shortlabel,'(')
 call string_delete(shortlabel,')')
 call string_delete(shortlabel,'{')
 call string_delete(shortlabel,'}')
 call string_delete(shortlabel,'[')
 call string_delete(shortlabel,']')
 call string_delete(shortlabel,'<')
 call string_delete(shortlabel,'>')
 call string_delete(shortlabel,'\(2268)')

end function shortlabel

!---------------------------------------------------------------
! interface for adjusting the label for column-integrated plots
!---------------------------------------------------------------
function integrate_label(labelin,iplot,izcol,normalise,iRescale,labelzint,&
                         projlabelformat,iapplyprojformat)
  use asciiutils,      only:string_replace
  implicit none
  character(len=*), intent(in) :: labelin,labelzint,projlabelformat
  integer, intent(in) :: iplot,izcol,iapplyprojformat
  logical, intent(in) :: normalise,iRescale
  character(len=len(label)+20) :: integrate_label

  if (len_trim(projlabelformat).ne.0 .and. (iapplyprojformat.eq.0 .or. iapplyprojformat.eq.iplot)) then
     integrate_label = projlabelformat
     call string_replace(integrate_label,'%l',trim(labelin))
     if (iRescale) then
        call string_replace(integrate_label,'%z',trim(label(izcol)(1:index(label(izcol),unitslabel(izcol))-1)))
        call string_replace(integrate_label,'%uz',trim(unitslabel(izcol)))
     else
        call string_replace(integrate_label,'%z',trim(label(izcol)))
     endif
  else
     if (normalise) then
        integrate_label = '< '//trim(labelin)//' >'
     else
        if (iRescale) then
           integrate_label = '\(2268) '//trim(labelin)//' d'// &
              trim(label(izcol)(1:index(label(izcol),unitslabel(izcol))-1))//trim(labelzint)
        else
           integrate_label = '\(2268) '//trim(labelin)//' d'//trim(label(izcol))
        endif
        if (iplot.eq.irho .and. (index(labelin,'density').ne.0 .or. index(labelin,'rho').ne.0)) then
           integrate_label = 'column density'
           !--try to get units label right for column density
           !  would be nice to have a more robust way of knowing what the units mean
           if (iRescale .and. index(labelzint,'cm').gt.0  &
                        .and. trim(adjustl(unitslabel(irho))).eq.'[g/cm\u3\d]') then
              integrate_label = trim(integrate_label)//' [g/cm\u2\d]'
           endif
        endif
     endif
  endif
end function integrate_label

!-----------------------------------------------------------------
!
! utility to "guess" which particle type contains sink particles
! from the label
!
!-----------------------------------------------------------------
integer function get_sink_type(ntypes)
 implicit none
 integer, intent(in) :: ntypes
 integer :: i

 get_sink_type = 0
 do i=1,ntypes
    if (get_sink_type.eq.0 .and. index(labeltype(i),'sink').ne.0) get_sink_type = i
 enddo

end function get_sink_type

!-----------------------------------------------------------------
!
!  utility to neatly print number of particles by type
!
!-----------------------------------------------------------------
subroutine print_types(noftype,ltype)
 integer,          dimension(:), intent(in) :: noftype
 character(len=*), dimension(:), intent(in) :: ltype
 integer :: itype,n,i
 character(len=1) :: sp

 i = 0
 sp = ' '
 do itype=1,size(noftype)
    n = noftype(itype)
    if (n > 0) then
       i = i + 1
       if (i > 1) sp = ','
       if (n < 10000) then
          write(*,"(a,i4)",advance='no') trim(sp)//' n('//trim(ltype(itype))//') = ',n
       elseif (n < 1000000) then
          write(*,"(a,i6)",advance='no') trim(sp)//' n('//trim(ltype(itype))//') = ',n
       elseif (n < 100000000) then
          write(*,"(a,i8)",advance='no') trim(sp)//' n('//trim(ltype(itype))//') = ',n
       else
          write(*,"(a,i10)",advance='no') trim(sp)//' n('//trim(ltype(itype))//') = ',n
       endif
    endif
 enddo
 write(*,*)

end subroutine print_types

end module labels
