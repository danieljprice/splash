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
 use params,     only:maxplot,maxparttypes,maxhdr,ltag
 use asciiutils, only:count_non_blank
 implicit none
 integer, parameter :: lenlabel = 80
 integer, parameter :: lenunitslabel = 40  ! length of units label
 character(len=lenlabel), dimension(maxplot+2) :: label,labelvec
 character(len=ltag), dimension(maxhdr)     :: headertags
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
 integer :: idustfracsum,ideltavsum
 integer :: igrainsize,igraindens,ivrel

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
 igrainsize = 0 ! grainsize
 igraindens = 0 ! graindens
 ivrel = 0      ! relative velocity
 iacplane = 0
 ike = 0
 idivB = 0
 iJfirst = 0
 icv = 0
 iradenergy = 0
 icolpixmap = 0
 irhorestframe = 0
 idustfrac = 0
 idustfracsum = 0
 ideltav = 0
 ideltavsum = 0
 headertags = ''

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
    if (ix(i)==icol) is_coord = .true.
 enddo

end function is_coord

!--------------------------------------------------------------
!
!  query function for whether column is a spatial coordinate
!  returns the location of the third dimension in the
!  list of columns
!
!--------------------------------------------------------------
integer function get_z_dir(ndim,iplotx,iploty) result(iplotz)
 implicit none
 integer, intent(in) :: ndim,iplotx,iploty
 integer :: i

 iplotz = 0
 do i=1,ndim
    if ((iplotx /= iploty.and. &
        (ix(i) /= iplotx).and.(ix(i) /= iploty))) iplotz = ix(i)
 enddo

end function get_z_dir

!--------------------------------------------------------------
!
!  query function for which coordinate is z
!  returns an integer between 1 and ndim
!
!--------------------------------------------------------------
integer function get_z_coord(ndim,iplotx,iploty) result(iplotz)
 implicit none
 integer, intent(in) :: ndim,iplotx,iploty
 integer :: i

 i = get_z_dir(ndim,iplotx,iploty)
 iplotz = i - ix(1) + 1
 if (iplotz < 0)    iplotz = 1
 if (iplotz > ndim) iplotz = ndim

end function get_z_coord

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
    if (len_trim(unitslab) > 0) then
       !--remove units label (only do this once)
       ipos = index(trim(shortstring),trim(unitslab))
       if (ipos /= 0) then
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
 use asciiutils,      only:string_replace,string_delete
 implicit none
 character(len=*), intent(in) :: labelin,labelzint,projlabelformat
 integer, intent(in) :: iplot,izcol,iapplyprojformat
 logical, intent(in) :: normalise,iRescale
 character(len=len(label)+20) :: integrate_label

 if (len_trim(projlabelformat) /= 0 .and. (iapplyprojformat==0 .or. iapplyprojformat==iplot)) then
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
          ! use composite units label e.g. \int rho_d [g/cm^3 pc]
          integrate_label = '\int '// &
             trim(labelin(1:index(labelin,unitslabel(iplot))-1))//' d'// &
             trim(label(izcol)(1:index(label(izcol),unitslabel(izcol))-1))// &
             get_unitlabel_coldens(iRescale,labelzint,unitslabel(iplot))
          ! use mix of units labels \int rho_d [g/cm^3] dz [pc]
          !integrate_label = '\int '//trim(labelin)//' d'// &
          !  trim(label(izcol)(1:index(label(izcol),unitslabel(izcol))-1))//trim(labelzint)
       else
          integrate_label = '\int '//trim(labelin)//' d'//trim(label(izcol))
       endif
       if (index(labelin,'\rho_{d,') > 0) then
          integrate_label = labelin(1:index(labelin,unitslabel(iplot))-1)
          call string_delete(integrate_label,'\rho_{d,')
          call string_delete(integrate_label,'}')
          integrate_label = trim(integrate_label)//' dust surface density'
          integrate_label = trim(integrate_label)//get_unitlabel_coldens(iRescale,labelzint,unitslabel(iplot))
       endif
       if (iplot==irho .and. (index(labelin,'density') /= 0 .or. index(labelin,'rho') /= 0)) then
          integrate_label = 'column density'
          integrate_label = trim(integrate_label)//get_unitlabel_coldens(iRescale,labelzint,unitslabel(irho))
       endif
    endif
 endif
end function integrate_label

!-----------------------------------------------------------------
!
! utility to convert the units label to g/cm^2 where appropriate
! would be nice to have a more robust way of knowing what the units mean
!
!-----------------------------------------------------------------
function get_unitlabel_coldens(iRescale,labelzint,unitlabel)
 use asciiutils, only:string_delete,string_replace
 logical, intent(in) :: iRescale
 character(len=*), intent(in) :: labelzint,unitlabel
 character(len=lenunitslabel) :: get_unitlabel_coldens

 if (iRescale) then
    get_unitlabel_coldens = trim(unitlabel)//trim(labelzint)
    call string_delete(get_unitlabel_coldens,']')
    call string_delete(get_unitlabel_coldens,'[')
    get_unitlabel_coldens = ' ['//trim(adjustl(get_unitlabel_coldens))//']'
    call string_replace(get_unitlabel_coldens,'/cm\u3\d cm','/cm^2')
    call string_replace(get_unitlabel_coldens,'/cm^3 cm','/cm^2')
 else
    get_unitlabel_coldens = ' '
 endif

end function get_unitlabel_coldens

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
    if (get_sink_type==0 .and. index(labeltype(i),'sink') /= 0) get_sink_type = i
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

!-----------------------------------------------------------------
!
!  utility to make labels for vector quantities
!  these change depending on the coordinate system
!  e.g. v_x, v_y, v_z; B_x, B_y, B_z
!
!-----------------------------------------------------------------
subroutine make_vector_label(lvec,ivec,nvec,iamveci,labelveci,labeli,labelx)
 character(len=*), intent(in)    :: lvec
 integer,          intent(in)    :: ivec,nvec
 integer,          intent(inout) :: iamveci(:)
 character(len=*), intent(inout) :: labelveci(:),labeli(:)
 character(len=*), intent(in)    :: labelx(3)
 integer :: i

 if (ivec > 0 .and. ivec+nvec <= size(labeli)) then
    iamveci(ivec:ivec+nvec-1)   = ivec
    labelveci(ivec:ivec+nvec-1) = lvec
    do i=1,nvec
       labeli(ivec+i-1) = trim(lvec)//'_'//labelx(i)
    enddo
 endif

end subroutine make_vector_label

!-----------------------------------------------------------------
!
!  utility to neatly format a grain size for use in plot label
!
!-----------------------------------------------------------------
function get_label_grain_size(sizecm) result(string)
 use asciiutils, only:string_delete
 real, intent(in) :: sizecm
 character(len=16) :: string
 character(len=6) :: ulab

 if (sizecm > 1000.) then
    write(string,"(1pg10.3)") sizecm*0.001
    ulab = 'km'
 elseif (sizecm > 100.) then
    write(string,"(1pg10.3)") sizecm*0.01
    ulab = 'm'
 elseif (sizecm > 1.) then
    write(string,"(1pg10.3)") sizecm
    ulab = 'cm'
 elseif (sizecm > 0.1) then
    write(string,"(1pg10.3)") sizecm*10.
    ulab = 'mm'
 elseif (sizecm > 1.e-4) then
    write(string,"(1pg10.3)") sizecm*1.e4
    ulab = '\gmm'
 elseif (sizecm > 1.e-7) then
    write(string,"(1pg10.3)") sizecm*1.e7
    ulab = 'nm'
 endif
 string = adjustl(string)
 call string_delete(string,'.0 ')
 call string_delete(string,'. ')
 string = trim(string)//trim(ulab)

end function get_label_grain_size

end module labels
