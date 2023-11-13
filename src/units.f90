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
!  Copyright (C) 2005-2018 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!--------------------------------------------------------------------
! module containing subroutines to do with setting of physical units
!--------------------------------------------------------------------
module settings_units
 use params
 use labels, only:unitslabel,unitslabel_default,lenlabel,&
                  labelzintegration,labelzintegration_default
 implicit none
 real, dimension(0:maxplot), public :: units,units_default,units_old
 real, public :: unitzintegration,unitzintegration_default
 real(doub_prec), public :: unit_interp
 public :: set_units,read_unitsfile,write_unitsfile,defaults_set_units
 public :: get_nearest_length_unit,get_nearest_time_unit
 public :: get_nearest_mass_unit,get_nearest_velocity_unit

 integer, parameter :: nx = 9
 real(doub_prec), parameter :: unit_length(nx) = &
    (/1.d0,    &
      1.d2,    &
      2.01168d4, &
      1.d5,    &
      6.96d10, &
      1.496d13,&
      3.086d18,&
      3.086d21,&
      3.086d24/)

 character(len=*), parameter :: unit_labels_length(nx) = &
    (/' [cm]     ',&
      ' [m]      ',&
      ' [furlong]',&
      ' [km]     ',&
      ' [R_{Sun}]',&
      ' [au]     ',&
      ' [pc]     ',&
      ' [kpc]    ',&
      ' [Mpc]    '/)

 integer, parameter :: nt = 7
 real(doub_prec), parameter :: unit_time(nt) = &
    (/1.d-3,        &
      1.d0,         &
      3.6d3,        &
      8.64d4,       &
      3.15568926d7, &
      3.15568926d13,&
      3.15568926d16/)

 character(len=*), parameter :: unit_labels_time(nt) = &
    (/' ms  ',&
      ' s   ',&
      ' hrs ',&
      ' days',&
      ' yrs ',&
      ' Myr ',&
      ' Gyr '/)

 integer, parameter :: nv = 3
 real(doub_prec), parameter :: unit_vel(nv) = &
    (/1.d0,  &
      1.d2,  &
      1.d5/)

 character(len=*), parameter :: unit_labels_vel(nv) = &
    (/' [cm/s] ',&
      ' [m/s]  ',&
      ' [km/s] '/)

 integer, parameter :: nm = 6
 real(doub_prec), parameter :: unit_mass(nm) = &
    (/1.d0,  &
      1.d3,  &
      8.958d23,  &  ! Ceres mass
      5.979d27,  &  ! Earth mass
      1.898d30,  &  ! Jupiter mass
      1.989d33/)    ! Solar mass

 character(len=*), parameter :: unit_labels_mass(nm) = &
    (/' [g]         ',&
      ' [kg]        ',&
      ' [M_{Ceres}] ',&
      ' [M_{Earth}] ',&
      ' [M_{Jup}]   ',&
      ' [M_{Sun}]   '/)

 private

contains
!---------------------------------------------------------------
!
!  initialise the units arrays to some harmless default values
!
!---------------------------------------------------------------
subroutine defaults_set_units

 units(:) = 1.0
 units_default(:) = 1.0
 unitzintegration = 1.0
 unitzintegration_default = 1.0
 unit_interp      = 1.0d0
 unitslabel(:) = ' '
 unitslabel_default(:) = ' '
 labelzintegration = ' '
 labelzintegration_default = ' '

end subroutine defaults_set_units

!-------------------------------------------
!
!  find the nearest 'sensible' length unit
!
!-------------------------------------------
subroutine get_nearest_length_unit(udist,unit,unitlabel)
 real(doub_prec),  intent(in)  :: udist
 real(doub_prec),  intent(out) :: unit
 character(len=*), intent(out) :: unitlabel

 call get_nearest_unit(nx,unit_length,unit_labels_length,udist,unit,unitlabel)

end subroutine get_nearest_length_unit

!-------------------------------------------
!
!  find the nearest 'sensible' time unit
!
!-------------------------------------------
subroutine get_nearest_time_unit(utime,unit,unitlabel)
 real(doub_prec),  intent(in)  :: utime
 real(doub_prec),  intent(out) :: unit
 character(len=*), intent(out) :: unitlabel

 call get_nearest_unit(nt,unit_time,unit_labels_time,utime,unit,unitlabel)

end subroutine get_nearest_time_unit

!-------------------------------------------
!
!  find the nearest 'sensible' velocity unit
!
!-------------------------------------------
subroutine get_nearest_velocity_unit(uvel,unit,unitlabel)
 real(doub_prec),  intent(in)  :: uvel
 real(doub_prec),  intent(out) :: unit
 character(len=*), intent(out) :: unitlabel

 call get_nearest_unit(nv,unit_vel,unit_labels_vel,uvel,unit,unitlabel)

end subroutine get_nearest_velocity_unit

!-------------------------------------------
!
!  find the nearest 'sensible' mass unit
!
!-------------------------------------------
subroutine get_nearest_mass_unit(umass,unit,unitlabel)
 real(doub_prec),  intent(in)  :: umass
 real(doub_prec),  intent(out) :: unit
 character(len=*), intent(out) :: unitlabel

 call get_nearest_unit(nm,unit_mass,unit_labels_mass,umass,unit,unitlabel)

end subroutine get_nearest_mass_unit

!-------------------------------------------------------
!
!  find the nearest unit from a list of possible units
!
!-------------------------------------------------------
subroutine get_nearest_unit(nu,units,unit_labels,unit_in,unit_out,unitlabel_out)
 integer, intent(in) :: nu
 real(doub_prec),  intent(in) :: units(nu),unit_in
 character(len=*), intent(in) :: unit_labels(nu)
 real(doub_prec),  intent(out) :: unit_out
 character(len=*), intent(out) :: unitlabel_out
 real(doub_prec) :: err,erri
 integer :: i

 err = huge(err)
 do i = 1,nu
    ! find nearest unit in log space
    erri = abs(log10(unit_in)-log10(units(i)))
    if (erri < err) then
       unit_out      = unit_in/units(i)
       unitlabel_out = unit_labels(i)
       err = erri
    endif
 enddo

end subroutine get_nearest_unit

!-------------------------------------------------------
!
!  set units
!
!-------------------------------------------------------
subroutine set_units(ncolumns,numplot,UnitsHaveChanged)
 use prompting,     only:prompt
 use labels,        only:label,ix,ih,ivx,ipmass,idivB,iamvec,labelvec,headertags,strip_units
 use settings_data, only:ndim,ndimV,ivegotdata
 use particle_data, only:headervals,maxstep
 use asciiutils,    only:get_value
 integer, intent(in) :: ncolumns,numplot
 logical, intent(out) :: UnitsHaveChanged
 integer :: icol,ibs,ibc
 real(doub_prec) :: unitsprev
 real(doub_prec) :: udist,utime,umass
 logical :: applytoall,got_label
 character(len=lenlabel) :: mylabel

 icol = 1
 ! try to extract code units from the file headers, these only exist in some data reads
 utime = 0.d0
 udist = 0.d0
 umass = 0.d0
 if (maxstep > 0 .and. ivegotdata) then
    utime = get_value('utime',headertags,headervals(:,1))
    udist = get_value('udist',headertags,headervals(:,1))
    umass = get_value('umass',headertags,headervals(:,1))
 endif

 do while(icol >= 0)
    icol = -1
    call prompt('enter column to change units (-2=reset all,-1=quit,0=time)',icol,-2,numplot)
    if (icol >= 0) then
       unitsprev = units(icol)
       if (icol > ncolumns) then
          print "(a)",' WARNING: calculated quantities are automatically calculated in physical units '
          print "(a)",' this means that units set here will be re-scalings of these physical values'
       endif
       got_label = .true.
       if (icol==0 .and. utime > 0.) then
          ! give hints for possible time units, if utime is read from data file
          call choose_unit_from_list(units(icol),unitslabel(icol),&
                                     nt,unit_labels_time,unit_time,utime,'time')
       elseif (any(ix==icol) .or. icol==ih .and. udist > 0.) then
          ! give hints for possible length units, if udist is read from data file
          call choose_unit_from_list(units(icol),unitslabel(icol),&
                                     nx,unit_labels_length,unit_length,udist,'length')
       elseif (ivx==icol .and. udist > 0. .and. utime > 0.) then
          ! give hints for velocity units, if both udist and utime read from data file
          call choose_unit_from_list(units(icol),unitslabel(icol),&
                                     nv,unit_labels_vel,unit_vel,udist/utime,'velocity')
       elseif (ipmass==icol .and. umass > 0.) then
          ! give hints for velocity units, if both udist and utime read from data file
          call choose_unit_from_list(units(icol),unitslabel(icol),&
                                     nm,unit_labels_mass,unit_mass,umass,'mass')
       elseif (icol > 0) then
          mylabel = strip_units(label(icol),unitslabel(icol))
          call prompt('enter '//trim(mylabel)//' units (new=old*units)',units(icol))
          got_label = .false.
       else ! icol = 0
          call prompt('enter time units (new=old*units)',units(icol))
          got_label = .false.
       endif

       if (abs(units(icol)) > tiny(units)) then ! sanity check the units
          if (abs(units(icol) - unitsprev) > tiny(units)) UnitsHaveChanged = .true.

          !--prompt for the units label
          if (.not.got_label) then ! skip if label set in choose_unit_from_list
             !--suggest a label amendment if none already set
             if (len_trim(unitslabel(icol))==0) &
                call suggest_units_label(units(icol), unitslabel(icol))

             call prompt('enter label amendment ',unitslabel(icol))
          endif
       else
          UnitsHaveChanged = .true.
          units(icol) = 1.0
          unitslabel(icol) = ' '
       endif
       if (UnitsHaveChanged .and. icol > 0) then
          !
          !--prompt to apply same units to coordinates and h for consistency
          !
          if (any(ix(1:ndim)==icol) .or. (icol==ih .and. ih > 0)) then
             applytoall = .true.
             !--try to make prompts apply to whichever situation we have
             if (ndim==1) then
                if (icol==ix(1) .and. ih > 0) then
                   call prompt(' Apply these units to h?',applytoall)
                else
                   call prompt(' Apply these units to '//trim(label(ix(1)))//'?',applytoall)
                endif
             elseif (any(ix(1:ndim)==icol) .and. ih > 0) then
                call prompt(' Apply these units to all coordinates and h?',applytoall)
             else
                call prompt(' Apply these units to all coordinates?',applytoall)
             endif
             if (applytoall) then
                units(ix(1:ndim)) = units(icol)
                unitslabel(ix(1:ndim)) = unitslabel(icol)
                if (ih > 0) then
                   if (idivb > 0) then
                      ! amend units of div B and curl B
                      ibs = index(unitslabel(idivB),'/')
                      ibc = index(unitslabel(icol),'[')
                      units(idivB:idivB+3) = units(idivB:idivB+3)*units(ih)/units(icol)
                      unitslabel(idivB:idivB+3) = unitslabel(idivB)(1:ibs)//unitslabel(icol)(ibc+1:)
                   endif
                   units(ih) = units(icol)
                   unitslabel(ih) = unitslabel(icol)
                endif
             endif
             !
             !--set units for z integration in 3D
             !  so for example can have x,y,z in kpc but column density in g/cm^2
             !
             if (abs(unitzintegration-1.0) <= tiny(unitzintegration)) then
                unitzintegration = units(icol)
                labelzintegration = unitslabel(icol)
             endif
             if (ndim==3) then
                call prompt(' Enter unit for ''z'' in 3D column integrated plots ',unitzintegration)
                call prompt(' Enter label for z integration unit (e.g. [cm])',labelzintegration)
             endif
          endif
          !
          !--also sensible to apply same units to all components of a vector
          !
          if (ndimV > 1 .and. iamvec(icol) > 0) then
             applytoall = .true.
             call prompt(' Apply these units to all components of '//trim(labelvec(icol))//'?',applytoall)
             if (applytoall) then
                where (iamvec(1:ncolumns)==iamvec(icol))
                   units(1:ncolumns) = units(icol)
                   unitslabel(1:ncolumns) = unitslabel(icol)
                end where
             endif
          endif
       endif

    elseif (icol==-2) then
       UnitsHaveChanged = .true.
       print "(/a)",' resetting all units to unity...'
       units = 1.0
       unitslabel = ' '
       unitzintegration = 1.0
       labelzintegration = ' '
    endif
    print*
 enddo

end subroutine set_units

!-------------------------------------------------------
!
!  select unit from a sensible list of presets
!
!-------------------------------------------------------
subroutine choose_unit_from_list(unit,unitlabel,n,unit_labels,unit_vals,unit_code,tag)
 use prompting, only:prompt
 integer, intent(in) :: n
 real,             intent(inout) :: unit
 character(len=*), intent(inout) :: unitlabel
 real(doub_prec), intent(in)     :: unit_code,unit_vals(n)
 character(len=*), intent(in)    :: unit_labels(n),tag
 integer :: i,iselect
 real :: unitsprev

 iselect = n+1
 ! get default value of iselect if one of the units is already chosen
 do i=1,n
    if (abs(unit-unit_code/unit_vals(i)) < tiny(0.)) iselect = i
 enddo
 do i=1,n
    print "(i2,')',a,' ( x ',1pg11.4,')')",i,unit_labels(i),unit_code/unit_vals(i)
 enddo
 print "(i2,') custom')",n+1
 call prompt('Enter choice of '//trim(tag)//' unit ',iselect,1,n+1)
 if (iselect==n+1) then
    unitsprev = unit
    call prompt('Enter custom unit (new=old*unit)',unit)
    if (unit_code > 0.d0) call suggest_label(unitsprev,unit,unit_code, &
                               unit_vals,unit_labels,unitlabel)
    call prompt('enter label amendment ',unitlabel)
 else
    unit = real(unit_code/unit_vals(iselect))
    unitlabel = unit_labels(iselect)
    print "(a,1pg11.4)", ' => '//trim(tag)//' unit is now '//trim(unitlabel)//&
                         ', i.e. '//trim(tag)//' = '//trim(tag)//'_code * ',unit
 endif

end subroutine choose_unit_from_list

!-------------------------------------------------------
!
!  suggest the units label based on matching list
!  of known scalings
!
!-------------------------------------------------------
subroutine suggest_label(unitsprev,unit,ucode,unit_suggest,unit_label_suggest,unitlabel)
 real, intent(in) :: unitsprev,unit
 real(doub_prec), intent(in) :: ucode,unit_suggest(:)
 character(len=*), intent(in) :: unit_label_suggest(:)
 character(len=*), intent(inout) :: unitlabel
 integer :: i

 ! do nothing if unit has not changed
 if (abs(unit - unitsprev) < epsilon(0.)) return

 ! suggest corresponding label to the unit chosen from list
 do i=1,size(unit_suggest)
    if (abs(unit/(ucode/unit_suggest(i)) - 1.) < 0.01) unitlabel = unit_label_suggest(i)
 enddo
 !print*,' suggested unit label = ',unitlabel

end subroutine suggest_label

!-------------------------------------------------------
!
!  suggest sensible units label for an arbitrary scaling
!
!-------------------------------------------------------
subroutine suggest_units_label(unit,label)
 real, intent(in) :: unit
 character(len=*), intent(out) :: label
 real(doub_prec) :: dunits

 dunits = 1./unit
 if (dunits > 100 .or. dunits < 1.e-1) then
    write(label,"(1pe8.1)") dunits
 else
    write(label,"(f5.1)") dunits
 endif
 label = ' [ x '//trim(adjustl(label))//' ]'

end subroutine suggest_units_label

!-------------------------------------------------------
!
!  save units for all columns to a file
!
!-------------------------------------------------------
subroutine write_unitsfile(unitsfile,ncolumns)
 character(len=*), intent(in) :: unitsfile
 integer, intent(in) :: ncolumns
 integer :: i,ierr

 print "(1x,a)",'saving units to '//trim(unitsfile)

 open(unit=77,file=unitsfile,status='replace',form='formatted',iostat=ierr)
 if (ierr /=0) then
    print "(1x,a)",'ERROR: cannot write units file'
 else
    write(77,*,iostat=ierr) units(0),';',trim(unitslabel(0)),' ;',unitzintegration,';',trim(labelzintegration)
    do i=1,ncolumns
       write(77,*,iostat=ierr) units(i),';',trim(unitslabel(i))
       if (ierr /= 0) then
          print "(1x,a)",'ERROR whilst writing units file'
          close(unit=77)
          return
       endif
    enddo
 endif
 close(unit=77)

 return

end subroutine write_unitsfile

!-------------------------------------------------------
!
!  read units for all columns from a file
!
!-------------------------------------------------------
subroutine read_unitsfile(unitsfile,ncolumns,ierr,iverbose)
 use settings_data, only:ndim
 character(len=*), intent(in) :: unitsfile
 integer, intent(in) :: ncolumns
 integer, intent(out) :: ierr
 integer, intent(in), optional :: iverbose
 character(len=2*len(unitslabel)+40) :: line
 integer :: i,itemp,isemicolon,isemicolon2,isemicolon3,isverbose
 logical :: ierrzunits,iexist

 if (present(iverbose)) then
    isverbose= iverbose
 else
    isverbose = 1
 endif

 ierr = 0
 ierrzunits = .false.
 inquire(file=unitsfile,exist=iexist)
 if (.not.iexist) then
    if (isverbose > 1) print "(1x,a)",trim(unitsfile)//' not found'
    ierr = 1
    return
 endif

 open(unit=78,file=unitsfile,status='old',form='formatted',err=997)
 if (isverbose > 0) print "(a)",' read '//trim(unitsfile)
 do i=0,maxplot  ! read all units possibly present in file
!
!    read a line from the file
!
    read(78,"(a)",err=998,end=999) line
!
!    now get units from the first part of the line
!
    read(line,*,iostat=itemp) units(i)
    if (itemp /= 0) print*,'ERROR reading units for column ',i
    if (units(i) > huge(units)) print "(/,a,i2)",' ERROR: UNITS ARE INFINITE FOR COLUMN ',i
    if (isnan(units(i))) print "(/,a,i2)",' ERROR: UNITS ARE NaN FOR COLUMN ',i
!
!    units label is what comes after the semicolon
!
    isemicolon = index(line,';')
    if (i==0) then
       !--time line also contains unit of z integration
       isemicolon2 = index(line(isemicolon+1:),';')
       if (isemicolon2 > 0) then
          isemicolon2 = isemicolon + isemicolon2
          unitslabel(i) = trim(line(isemicolon+1:isemicolon2-1))
          isemicolon3 = isemicolon2 + index(line(isemicolon2+1:),';')
          read(line(isemicolon2+1:isemicolon3-1),*,iostat=itemp) unitzintegration
          if (itemp /= 0) then
             print*,'error reading unit for z integration'
             ierrzunits = .true.
          else
             labelzintegration = trim(line(isemicolon3+1:))
          endif
       else
          ierrzunits = .true.
          print*,'error: could not read z integration unit from units file'
          if (isemicolon > 0) then
             unitslabel(i) = trim(line(isemicolon+1:))
          else
             print*,'error reading units label for column ',i
          endif
       endif
    else
       if (isemicolon > 0) then
          unitslabel(i) = trim(line(isemicolon+1:))
       else
          print*,'error reading units label for column ',i
       endif
    endif
!     print*,i,'units = ',units(i),'label = ',unitslabel(i)
 enddo
 if (ierrzunits .and. ndim==3) then
    unitzintegration = units(3)
    labelzintegration = unitslabel(3)
 endif

 close(unit=78)

 return

997 continue
 if (isverbose > 0) print*,trim(unitsfile),' not found'
 ierr = 1
 return
998 continue
 print*,'*** error reading units from '//trim(unitsfile)
 ierr = 2
 if (ierrzunits .and. ndim==3) then
    unitzintegration = units(3)
    labelzintegration = unitslabel(3)
 endif
 close(unit=78)
 return
999 continue
 !--only give error if we really do not have enough columns
 !  (on first call nextra is not set)
 if (i <= ncolumns) then
    print "(1x,a,i2)",'end of file in '//trim(unitsfile)//': units read to column ',i
    ierr = -1
 endif
 if (ierrzunits .and. ndim==3) then
    unitzintegration = units(3)
    labelzintegration = unitslabel(3)
 endif

 close(unit=78)
 return

end subroutine read_unitsfile

end module settings_units
