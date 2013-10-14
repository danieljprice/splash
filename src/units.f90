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

!--------------------------------------------------------------------
! module containing subroutines to do with setting of physical units
!--------------------------------------------------------------------
module settings_units
 use params
 use labels, only:unitslabel,labelzintegration
 implicit none
 real, dimension(0:maxplot), public :: units
 real, public :: unitzintegration
 real(doub_prec), public :: unit_interp
 public :: set_units,read_unitsfile,write_unitsfile,defaults_set_units

 private

contains
!
!--initialise the units arrays to some harmless default values
!
subroutine defaults_set_units
 implicit none

 units(:) = 1.0
 unitzintegration = 1.0
 unit_interp      = 1.0d0
 unitslabel(:) = ' '
 labelzintegration = ' '

 return
end subroutine defaults_set_units
!
!--set units
!
subroutine set_units(ncolumns,numplot,UnitsHaveChanged)
  use prompting, only:prompt
  use labels, only:label,ix,ih,iamvec,labelvec
  use settings_data, only:ndim,ndimV
  implicit none
  integer, intent(in) :: ncolumns,numplot
  logical, intent(out) :: UnitsHaveChanged
  integer :: icol
  real :: unitsprev,dunits
  logical :: applytoall

  icol = 1
  do while(icol.ge.0)
     icol = -1
     call prompt('enter column to change units (-2=reset all,-1=quit,0=time)',icol,-2,numplot)
     if (icol.ge.0) then
        unitsprev = units(icol)
        if (icol.gt.ncolumns) then
           print "(a)",' WARNING: calculated quantities are automatically calculated in physical units '
           print "(a)",' this means that units set here will be re-scalings of these physical values'
        endif
        if (icol.eq.0) then
           call prompt('enter time units (new=old*units)',units(icol))
        else
           call prompt('enter '//trim(label(icol))//' units (new=old*units)',units(icol))
        endif
        if (abs(units(icol)).gt.tiny(units)) then
           if (abs(units(icol) - unitsprev).gt.tiny(units)) UnitsHaveChanged = .true.
           if (len_trim(unitslabel(icol)).eq.0) then
           !--suggest a label amendment if none already set
              dunits = 1./units(icol)
              if (dunits.gt.100 .or. dunits.lt.1.e-1) then
                 write(unitslabel(icol),"(1pe8.1)") dunits
              else
                 write(unitslabel(icol),"(f5.1)") dunits
              endif
              unitslabel(icol) = ' [ x '//trim(adjustl(unitslabel(icol)))//' ]'
           endif
           !--label amendment can be overwritten
           call prompt('enter label amendment ',unitslabel(icol))
        else
           UnitsHaveChanged = .true.
           units(icol) = 1.0
           unitslabel(icol) = ' '
        endif
        if (UnitsHaveChanged .and. icol.gt.0) then
           !
           !--prompt to apply same units to coordinates and h for consistency
           !
           if (any(ix(1:ndim).eq.icol) .or. (icol.eq.ih .and. ih.gt.0)) then
              applytoall = .true.
              !--try to make prompts apply to whichever situation we have
              if (ndim.eq.1) then
                 if (icol.eq.ix(1) .and. ih.gt.0) then
                    call prompt(' Apply these units to h?',applytoall)
                 else
                    call prompt(' Apply these units to '//trim(label(ix(1)))//'?',applytoall)
                 endif
              elseif (any(ix(1:ndim).eq.icol) .and. ih.gt.0) then
                 call prompt(' Apply these units to all coordinates and h?',applytoall)
              else
                 call prompt(' Apply these units to all coordinates?',applytoall)
              endif
              if (applytoall) then
                 units(ix(1:ndim)) = units(icol)
                 unitslabel(ix(1:ndim)) = unitslabel(icol)
                 if (ih.gt.0) then
                    units(ih) = units(icol)
                    unitslabel(ih) = unitslabel(icol)
                 endif
              endif
              !
              !--set units for z integration in 3D
              !  so for example can have x,y,z in kpc but column density in g/cm^2
              !
              if (abs(unitzintegration-1.0).le.tiny(unitzintegration)) then
                 unitzintegration = units(icol)
                 labelzintegration = unitslabel(icol)
              endif
              if (ndim.eq.3) then
                 call prompt(' Enter unit for ''z'' in 3D column integrated plots ',unitzintegration)
                 call prompt(' Enter label for z integration unit (e.g. [cm])',labelzintegration)
              endif
           endif
           !
           !--also sensible to apply same units to all components of a vector
           !
           if (ndimV.gt.1 .and. iamvec(icol).gt.0) then
              applytoall = .true.
              call prompt(' Apply these units to all components of '//trim(labelvec(icol))//'?',applytoall)
              if (applytoall) then
                 where (iamvec(1:ncolumns).eq.iamvec(icol))
                    units(1:ncolumns) = units(icol)
                    unitslabel(1:ncolumns) = unitslabel(icol)
                 end where
              endif
           endif
        endif

     elseif (icol.eq.-2) then
        UnitsHaveChanged = .true.
        print "(/a)",' resetting all units to unity...'
        units = 1.0
        unitslabel = ' '
     endif
     print*
  enddo

end subroutine set_units
!
!--save units for all columns to a file
!
subroutine write_unitsfile(unitsfile,ncolumns)
  implicit none
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
!
!--read units for all columns from a file
!
subroutine read_unitsfile(unitsfile,ncolumns,ierr,iverbose)
  use settings_data, only:ndim
  implicit none
  character(len=*), intent(in) :: unitsfile
  integer, intent(in) :: ncolumns
  integer, intent(out) :: ierr
  integer, intent(in), optional :: iverbose
  character(len=2*len(unitslabel)+40) :: line
  integer :: i,itemp,isemicolon,isemicolon2,isemicolon3
  logical :: ierrzunits,iexist,verbose
  
  if (present(iverbose)) then
     verbose= (iverbose.gt.0)
  else
     verbose = .true.
  endif

  ierr = 0
  ierrzunits = .false.
  inquire(file=unitsfile,exist=iexist)
  if (.not.iexist) then
     if (verbose) print "(1x,a)",trim(unitsfile)//' not found'
     ierr = 1
     return
  endif

  open(unit=78,file=unitsfile,status='old',form='formatted',err=997)
  if (verbose) print "(a)",' read '//trim(unitsfile)
  do i=0,maxplot  ! read all units possibly present in file
!
!    read a line from the file
!
     read(78,"(a)",err=998,end=999) line
!
!    now get units from the first part of the line
!
     read(line,*,iostat=itemp) units(i)
     if (itemp /= 0) print*,'error reading units for column ',i
!
!    units label is what comes after the semicolon
!
     isemicolon = index(line,';')
     if (i.eq.0) then
        !--time line also contains unit of z integration
        isemicolon2 = index(line(isemicolon+1:),';')
        if (isemicolon2.gt.0) then
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
           if (isemicolon.gt.0) then
              unitslabel(i) = trim(line(isemicolon+1:))
           else
              print*,'error reading units label for column ',i
           endif
        endif
     else
        if (isemicolon.gt.0) then
           unitslabel(i) = trim(line(isemicolon+1:))
        else
           print*,'error reading units label for column ',i
        endif
     endif
!     print*,i,'units = ',units(i),'label = ',unitslabel(i)
  enddo
  if (ierrzunits .and. ndim.eq.3) then
     unitzintegration = units(3)
     labelzintegration = unitslabel(3)
  endif

  close(unit=78)

  return

997 continue
  if (verbose) print*,trim(unitsfile),' not found'
  ierr = 1
  return
998 continue
  print*,'*** error reading units from '//trim(unitsfile)
  ierr = 2
  if (ierrzunits .and. ndim.eq.3) then
     unitzintegration = units(3)
     labelzintegration = unitslabel(3)
  endif
  close(unit=78)
  return
999 continue
  !--only give error if we really do not have enough columns
  !  (on first call nextra is not set)
  if (i.le.ncolumns) then
     print "(1x,a,i2)",'end of file in '//trim(unitsfile)//': units read to column ',i
     ierr = -1
  endif
  if (ierrzunits .and. ndim.eq.3) then
     unitzintegration = units(3)
     labelzintegration = unitslabel(3)
  endif

  close(unit=78)
  return

end subroutine read_unitsfile

end module settings_units
