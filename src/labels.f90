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
!  routines to do with storing and handling of plot labels
!
!-----------------------------------------------------------
module labels
 use params,     only:maxplot,maxparttypes,maxhdr,ltag
 use asciiutils, only:count_non_blank
 implicit none
 integer, parameter :: lenlabel = 80
 integer, parameter :: lenunitslabel = 40  ! length of units label
 character(len=lenlabel), dimension(maxplot+2) :: label,labelvec,labelreq,labelorig
 character(len=ltag), dimension(maxhdr)     :: headertags
 character(len=20), dimension(maxparttypes) :: labeltype
 character(len=6), parameter :: labeldefault = 'column'
 character(len=lenunitslabel), dimension(0:maxplot), public :: unitslabel,unitslabel_default
 character(len=lenunitslabel) :: labelzintegration,labelzintegration_default
 integer, dimension(3)       :: ix
 integer, dimension(maxplot) :: iamvec
 integer :: ivx,irho,iutherm,ipr,ih,irad,iBfirst,iBpol,iBtor,iax
 integer :: ipmass,ike,ispsound,itemp,ikappa
 integer :: idivb,iJfirst,irhostar,ipmomx
 integer :: iacplane,ipowerspec
 integer :: icv,iradenergy
 integer :: isurfdens,itoomre
 integer :: ipdf,icolpixmap
 integer :: irhorestframe,idustfrac,ideltav
 integer :: idustfracsum,ideltavsum
 integer :: igrainsize,igraindens,ivrel
 integer :: irhodust_start,irhodust_end
 integer :: nreq

 public

contains

!--------------------------------------------------------------
!
!  utility to reset default settings for column identification
!
!--------------------------------------------------------------
subroutine reset_columnids
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
 itemp = 0      ! temperature
 ikappa = 0     ! opacity
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
 ipmomx = 0
 irhodust_start = 0
 irhodust_end = 0
 headertags = ''

end subroutine reset_columnids

!--------------------------------------------------------------
!
!  set default column labels
!
!--------------------------------------------------------------
subroutine set_default_labels(mylabel)
 character(len=*), intent(out) :: mylabel(:)
 integer :: i

 do i=1,size(mylabel)
    write(mylabel(i),"(a,1x,i3)") trim(labeldefault),i
 enddo

end subroutine set_default_labels

!--------------------------------------------------------------
!
!  query function for whether column is a spatial coordinate
!
!--------------------------------------------------------------
logical function is_coord(icol,ndim)
 integer, intent(in) :: icol,ndim
 integer :: i

 is_coord = .false.
 do i=1,ndim
    if (ix(i)==icol) is_coord = .true.
 enddo

end function is_coord

!--------------------------------------------------------------
!
!  query function for whether column is a density
!  mainly used to decide whether or not to use
!  mass-weighted interpolation
!
!--------------------------------------------------------------
logical function is_density(icol)
 integer, intent(in) :: icol

 is_density = .false.

 ! if we match known density columns
 if (icol==irho .or. icol==irhorestframe) is_density = .true.

 ! if label contains rho or dens
 if (icol > 0 .and. icol < size(label)) then
    if (index(label(icol),'rho') > 0) is_density = .true.
    if (index(label(icol),'dens') > 0) is_density = .true.
 else
    is_density = .true.  ! icol = 0 should be true
 endif

end function is_density

!----------------------------------------------------------------
!
!  function returning the synonym of a label
!  so two different labels can be identified as the same
!  physical quantity,  e.g. "radius" and "r" are both the radius
!
!----------------------------------------------------------------
elemental function label_synonym(string)
 use asciiutils, only:lcase,string_delete
 character(len=*), intent(in) :: string
 character(len=max(len(string),8)) :: label_synonym
 character(len=len(string)) :: labeli
 integer :: k

 ! remove leading spaces and make lower case
 labeli = trim(adjustl(lcase(string)))

 ! split at whitespace or [ to avoid unit labels
 k = max(index(labeli,' '),index(labeli,'['))
 if (k > 1) then
    labeli = labeli(1:k-1)
 endif
 ! also remove special characters
 call string_delete(labeli,'\d')
 call string_delete(labeli,'\u')
 call string_delete(labeli,'\')
 call string_delete(labeli,'_')

 ! remove log from the front of the label
 if (labeli(1:3)=='log') labeli = labeli(4:)

 if (labeli(1:3)=='den' .or. index(labeli,'rho') /= 0 .or. labeli(1:2)=='gr' .or. &
     (index(labeli,'density') /= 0 .and. irho==0)) then
    label_synonym = 'density'
 elseif (labeli(1:5)=='pmass' .or. labeli(1:13)=='particle mass') then
    label_synonym = 'pmass'
 elseif (trim(labeli)=='h' .or. labeli(1:6)=='smooth') then
    label_synonym = 'h'
 elseif (trim(labeli)=='u' .or. labeli(1:6)=='utherm' .or. labeli(1:5)=='eint' &
     .or.(index(labeli,'internal energy') /= 0)) then
    label_synonym = 'u'
 elseif (labeli(1:2)=='pr' .or. trim(labeli)=='p' .or. &
        (index(labeli,'pressure') /= 0 .and. ipr==0)) then
    label_synonym = 'pressure'
 elseif (trim(labeli)=='r' .or. labeli(1:3)=='rad') then
    label_synonym = 'radius'
 else
    ! by default the label synonym is the same as the original label
    label_synonym = labeli
 endif

end function label_synonym

!--------------------------------------------------------------
!
!  query function for whether column is a spatial coordinate
!  returns the location of the third dimension in the
!  list of columns
!
!--------------------------------------------------------------
integer function get_z_dir(ndim,iplotx,iploty) result(iplotz)
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
 integer, intent(in) :: ndim,iplotx,iploty
 integer :: i

 i = get_z_dir(ndim,iplotx,iploty)
 iplotz = i - ix(1) + 1
 if (iplotz < 0)    iplotz = 1
 if (iplotz > ndim) iplotz = ndim

end function get_z_coord

!-----------------------------------------------------------------
!
!  utility to strip the units label from a label string
!
!-----------------------------------------------------------------
elemental function strip_units(string,unitslab)
 character(len=lenlabel), intent(in) :: string
 character(len=*),        intent(in) :: unitslab
 character(len=lenlabel)             :: strip_units
 integer :: ipos

 strip_units = string
 if (len_trim(unitslab) > 0) then
    !--remove units label (only do this once)
    ipos = index(trim(strip_units),trim(unitslab))
    if (ipos /= 0) then
       strip_units = strip_units(1:ipos-1)//&
                     strip_units(ipos+len_trim(unitslab)+1:len_trim(strip_units))
    endif
 endif

end function strip_units

!-----------------------------------------------------------------
!
!  utility to strip spaces, escape sequences and
!  units labels from strings (this can be called for both
!  function strings and variable labels)
!
!-----------------------------------------------------------------
elemental function shortstring(string,unitslab)
 use asciiutils, only:string_delete
 character(len=lenlabel), intent(in)           :: string
 character(len=*),        intent(in), optional :: unitslab
 character(len=lenlabel)                       :: shortstring

 shortstring = string
 !--strip off the units label
 if (present(unitslab)) shortstring = strip_units(shortstring,unitslab)

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
elemental function shortlabel(string,unitslab,lc)
 use asciiutils, only:string_delete,lcase
 character(len=lenlabel), intent(in)           :: string
 character(len=*),        intent(in), optional :: unitslab
 character(len=lenlabel)                       :: shortlabel
 logical, intent(in), optional :: lc

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

 if (present(lc)) then
    if (lc) shortlabel = lcase(shortlabel)
 endif

end function shortlabel

!---------------------------------------------------------------
! interface for adjusting the label for column-integrated plots
!---------------------------------------------------------------
function integrate_label(labelin,iplot,izcol,normalise,iRescale,labelzint,&
                         projlabelformat,iapplyprojformat)
 use asciiutils,      only:string_replace,string_delete
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
          integrate_label = '\int '//trim(labelin)//' d'//trim(shortlabel(label(izcol),unitslabel(izcol)))
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

 if (iRescale .and. len_trim(labelzint) > 0) then
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
 integer, intent(in) :: ntypes
 integer :: i

 get_sink_type = 0
 do i=1,ntypes
    if (get_sink_type==0 .and. index(labeltype(i),'sink') /= 0) get_sink_type = i
    if (get_sink_type==0 .and. index(labeltype(i),'compact object') /= 0) get_sink_type = i
    if (get_sink_type==0 .and. index(labeltype(i),'point mass') /= 0) get_sink_type = i
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

 if (sizecm >= 1000.) then
    write(string,"(1pg10.3)") sizecm*1.e-5
    ulab = 'km'
 elseif (sizecm >= 100.) then
    write(string,"(1pg10.3)") sizecm*0.01
    ulab = 'm'
 elseif (sizecm >= 1.) then
    write(string,"(1pg10.3)") sizecm
    ulab = 'cm'
 elseif (sizecm >= 0.09) then
    write(string,"(1pg10.3)") sizecm*10.
    ulab = 'mm'
 elseif (sizecm >= 1.e-4) then
    write(string,"(1pg10.3)") sizecm*1.e4
    ulab = '\gmm'
 elseif (sizecm >= 1.e-7) then
    write(string,"(1pg10.3)") sizecm*1.e7
    ulab = 'nm'
 else
    write(string,"(1pg10.3)") sizecm
    ulab = 'cm'
 endif
 string = adjustl(string)
 call string_delete(string,'.000 ')
 call string_delete(string,'.00 ')
 call string_delete(string,'.0 ')
 call string_delete(string,'. ')
 string = trim(string)//trim(ulab)

end function get_label_grain_size

!-----------------------------------------------------------------
!
!  save the list of labels that are actually used for plotting
!
!-----------------------------------------------------------------
subroutine set_required_labels(required)
 logical, intent(in) :: required(0:)
 integer :: icol

 ! save original list of labels
 labelorig = label

 nreq = 0
 labelreq = ''
 do icol=1,size(required)-1
    nreq = nreq + 1
    if (required(icol)) then
       labelreq(icol) = shortlabel(label(icol),unitslabel(icol))
    endif
 enddo

 !print*,'DEBUG: required labels:'
 !do icol=1,nreq
 !   if (len_trim(labelreq(icol)) > 0) print*,trim(labelreq(icol))
 !enddo

end subroutine set_required_labels

!-----------------------------------------------------------------
!
!  see if a column has shifted in the actual data
!
!-----------------------------------------------------------------
integer function check_for_shifted_column(icol,labelnew) result(inew)
 use asciiutils, only:match_tag
 integer, intent(in) :: icol
 character(len=*), intent(out), optional :: labelnew
 character(len=lenlabel) :: labelcol

 inew = icol
 labelcol = shortlabel(label(icol),unitslabel(icol))

 if (trim(labelcol) /= trim(labelreq(icol)) .and. len_trim(labelreq(icol)) > 0) then
    if (.not.present(labelnew)) write(*,"(1x,a,i3,a)",advance='no') 'column ',icol,' has shifted: was '//&
          trim(labelreq(icol))//' but got '//trim(labelcol)
    inew = match_tag(shortlabel(label(1:maxplot),unitslabel(1:maxplot)),labelreq(icol))
    if (inew > 0) then
       if (present(labelnew)) then
          labelnew = trim(shortlabel(label(inew),unitslabel(inew)))//trim(unitslabel(icol))
       else
          print "(1x,a,i3)",': found '//trim(shortlabel(label(inew),unitslabel(inew)))//' in col ',inew
       endif
    else
       inew = icol
    endif
 endif

end function check_for_shifted_column

!-----------------------------------------------------------------
!
!  compute the backwards map from inew -> icol based on where
!  a column has been shifted to. This enables lookup of original
!  units etc which can be used to scale the data
!
!-----------------------------------------------------------------
function map_shifted_columns() result(imap)
 integer :: imap(maxplot)
 integer :: i,icol

 do i=1,size(imap)
    imap(i) = i
 enddo

 do i=1,size(imap)
    icol = i
    if (len_trim(label(i)) > 0) icol = check_for_shifted_column(i)
    if (icol /= i) then
       !print*,i,' setting imap=',icol
       imap(icol) = i
    endif
 enddo

end function map_shifted_columns

!-----------------------------------------------------------------
!
!  set labels for columns which have been tagged as vectors
!  using the iamvec label
!
!-----------------------------------------------------------------
subroutine set_vector_labels(ncolumns,ndimV,iamveci,labelveci,labeli,labelcoordi)
 integer,                 intent(in)    :: ncolumns,ndimV
 integer,                 intent(inout) :: iamveci(:)
 character(len=*),        intent(in)    :: labelcoordi(:)
 character(len=lenlabel), intent(inout) :: labelveci(:),labeli(:)
 character(len=lenlabel) :: tmplabel
 integer :: i
 !
 !--set labels for vector quantities
 !
 i = 1
 do while (i <= ncolumns)
    if (iamvec(i) > 0) then
       tmplabel = labeli(i)
       call make_vector_label(tmplabel,i,ndimV,iamveci,labelveci,labeli,labelcoordi)
       i = i + ndimV
    else
       i = i + 1
    endif
 enddo

end subroutine set_vector_labels

end module labels
