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

!-------------------------------------------------
!
! Module containing subroutines relating to
! setting/saving default options
!
!-------------------------------------------------
module defaults
 implicit none

contains
  !!
  !! set initial default options
  !! these are used if no defaults file is found
  !!
  subroutine defaults_set(use_evdefaults)
   use exact,              only:defaults_set_exact
   use multiplot
   use settings_limits,    only:defaults_set_limits
   use options_data,       only:defaults_set_data
   use settings_part,      only:defaults_set_part,defaults_set_part_ev
   use settings_page,      only:defaults_set_page,defaults_set_page_ev
   use settings_render,    only:defaults_set_render
   use settings_vecplot,   only:defaults_set_vecplot
   use settings_xsecrot,   only:defaults_set_xsecrotate
   use settings_powerspec, only:defaults_set_powerspec
   use settings_units,     only:defaults_set_units
   use settings_data,      only:iverbose
   use titles,             only:pagetitles,steplegend
   logical, intent(in) :: use_evdefaults
   integer :: i,iverbose_save
  !
  !--set defaults for submenu options
  !
   iverbose_save = iverbose
   call defaults_set_data
   call defaults_set_limits
   call defaults_set_page
   call defaults_set_part
   call defaults_set_render
   call defaults_set_xsecrotate
   call defaults_set_vecplot
   call defaults_set_exact
   call defaults_set_powerspec
   call defaults_set_units
  !
  !--if using evsplash, override some default options
  !
   if (use_evdefaults) then
      if (iverbose_save >= 0) print "(a)",' ** ev mode: using default settings for .ev files **'
      call defaults_set_page_ev
      call defaults_set_part_ev
   endif

   itrans(:) = 0
  !
  !--multiplot
  !
   nyplotmulti = 4           ! number of plots in multiplot
   multiploty(:) = 0
   do i=1,4
      multiploty(i) = 1+i       ! first plot : y axis
   enddo
   multiplotx(:) = 1          ! first plot : x axis
   irendermulti(:) = 0        ! rendering
   ivecplotmulti(:) = 0       ! vector plot
   x_secmulti(:) = .false.    ! take cross section?
   xsecposmulti(:) = 0.0      ! position of cross section
   icontourmulti(:) = 0       ! contour plot
   iplotpartoftypemulti(:,:) = .false.
   do i=1,size(iplotpartoftypemulti(:,1))
      iplotpartoftypemulti(i,i) = .true.
   enddo
   iusealltypesmulti(:) = .true.
   !
   !--titles
   !
   pagetitles = ' '
   steplegend = ' '

   return
  end subroutine defaults_set

!
! set defaults for 360 video mode
!
subroutine defaults_set_360()
 use settings_page,   only:defaults_set_page_360
 use settings_render, only:defaults_set_render_360
 use settings_data,   only:icoordsnew

 icoordsnew = 3
 call defaults_set_render_360
 call defaults_set_page_360()

end subroutine defaults_set_360

!
!     writes default options to file (should match defaults_read)
!
subroutine defaults_write(filename)
 use filenames, only:rootname,nfiles
 character(len=*), intent(in) :: filename
 integer :: i,ierr
 integer, parameter :: iunit = 1

 open(unit=iunit,file=trim(adjustl(filename)),status='replace',form='formatted',iostat=ierr)
 if (ierr /= 0) then
    print*,'ERROR: cannot write file '//trim(filename)
    close(unit=iunit)
    return
 endif
 call write_namelists(iunit)
 do i=1,nfiles
    write(iunit,"(a)") trim(rootname(i))
 enddo
 close(unit=iunit)
 print "(a)",'default options saved to file '//trim(filename)

end subroutine defaults_write

!-----------------------------------------------
! print all settings to the terminal, in the
! same format as the defaults file
!-----------------------------------------------
subroutine defaults_print()
 use, intrinsic :: iso_fortran_env, only:output_unit

 call write_namelists(output_unit)

end subroutine defaults_print

!-----------------------------------------------
! write all namelist groups to an open unit
! (should match defaults_read)
!-----------------------------------------------
subroutine write_namelists(iunit)
 use exact,              only:exactopts,exactparams
 use settings_data,      only:dataopts
 use settings_part,      only:plotopts
 use settings_page,      only:pageopts
 use settings_render,    only:renderopts
 use settings_vecplot,   only:vectoropts
 use settings_xsecrot,   only:xsecrotopts,animopts
 use settings_powerspec, only:powerspecopts
 use multiplot,          only:multi
 use shapes,             only:shapeopts
 use calcquantities,     only:calcopts
 integer, intent(in) :: iunit

 ! without delim namelists may not be readable
 write(iunit,NML=dataopts,delim='apostrophe')
 write(iunit,NML=plotopts,delim='apostrophe')
 write(iunit,NML=pageopts,delim='apostrophe')
 write(iunit,NML=renderopts,delim='apostrophe')
 write(iunit,NML=vectoropts,delim='apostrophe')
 write(iunit,NML=xsecrotopts,delim='apostrophe')
 write(iunit,NML=powerspecopts,delim='apostrophe')
 write(iunit,NML=exactopts,delim='apostrophe')
 write(iunit,NML=exactparams,delim='apostrophe')
 write(iunit,NML=multi,delim='apostrophe')
 write(iunit,NML=shapeopts,delim='apostrophe')
 write(iunit,NML=calcopts,delim='apostrophe')
 write(iunit,NML=animopts,delim='apostrophe')

end subroutine write_namelists
!-----------------------------------------------
! reads default options from file
! uses namelist input to group the options
! these are specified in the modules
!-----------------------------------------------
subroutine defaults_read(filename)
 use filenames,          only:rootname,maxfile
 use multiplot,          only:multi
 use settings_data,      only:dataopts,idustfrac_plot,idefaults_file_read,iverbose
 use settings_part,      only:plotopts
 use settings_page,      only:pageopts
 use settings_render,    only:renderopts
 use settings_vecplot,   only:vectoropts
 use settings_xsecrot,   only:xsecrotopts,animopts
 use settings_powerspec, only:powerspecopts
 use exact,              only:exactopts,exactparams
 use shapes,             only:shapeopts
 use calcquantities,     only:calcopts
 character(len=*), intent(in) :: filename
 logical :: iexist
 integer :: ierr,i,nerr
 integer, parameter :: iunit = 1

 nerr = 0
 inquire (exist=iexist, file=filename)
 if (iexist) then
    open(unit=iunit,file=filename,status='old',form='formatted',delim='apostrophe',err=88)

    ierr = 0
    read(iunit,NML=dataopts,iostat=ierr)
    idefaults_file_read = .true.
    if (ierr /= 0) nerr = nerr + 1

    ierr = 0
    read(iunit,NML=plotopts,iostat=ierr)
    if (ierr /= 0) nerr = nerr + 1

    ierr = 0
    read(iunit,NML=pageopts,iostat=ierr)
    if (ierr /= 0) nerr = nerr + 1

    ierr = 0
    read(iunit,NML=renderopts,iostat=ierr)
    if (ierr /= 0) nerr = nerr + 1

    ierr = 0
    read(iunit,NML=vectoropts,iostat=ierr)
    if (ierr /= 0) nerr = nerr + 1

    ierr = 0
    read(iunit,NML=xsecrotopts,iostat=ierr)
    if (ierr /= 0) nerr = nerr + 1

    ierr = 0
    read(iunit,NML=powerspecopts,iostat=ierr)
    if (ierr /= 0) nerr = nerr + 1

    ierr = 0
    read(iunit,NML=exactopts,iostat=ierr)
    if (ierr /= 0) nerr = nerr + 1

    ierr = 0
    read(iunit,NML=exactparams,iostat=ierr)
    if (ierr /= 0) nerr = nerr + 1

    ierr = 0
    read(iunit,NML=multi,iostat=ierr)
    if (ierr /= 0) nerr = nerr + 1

    ierr = 0
    read(iunit,NML=shapeopts,iostat=ierr)
    if (ierr /= 0) nerr = nerr + 1

    ierr = 0
    read(iunit,NML=calcopts,iostat=ierr)
    if (ierr /= 0) nerr = nerr + 1

    ierr = 0
    read(iunit,NML=animopts,iostat=ierr)
    if (ierr /= 0) nerr = nerr + 1

    if (len_trim(rootname(1))==0) then
       do i=1,maxfile
          read(iunit,*,end=66,iostat=ierr) rootname(i)
       enddo
    endif
66  continue

    close(unit=iunit)
    if (iverbose >= 0) then
       if (nerr > 0) then
          print "(a)",' WARNING: '//trim(filename)//' incomplete (from old code version)'
       else
          print*,'read '//trim(filename)
       endif
    endif
    return
 else
    if (iverbose >= 0) print*,trim(filename)//' not found: using default settings'
    return
 endif

88 continue
 if (iverbose >= 0) print "(a)",' *** error opening defaults file '//trim(filename)//': using default settings'

 return
end subroutine defaults_read

!-----------------------------------------------
! override any setting in the defaults file
! from the command line, e.g. --npix=500 or
! --iaxis=-2, by reading each --name=value
! argument as namelist input. Character
! variables can be given with or without quotes
!-----------------------------------------------
subroutine defaults_read_command_line()
 use settings_data, only:iverbose
 character(len=*), parameter :: groups(13) = &
    (/'dataopts     ','plotopts     ','pageopts     ','renderopts   ', &
      'vectoropts   ','xsecrotopts  ','powerspecopts','exactopts    ', &
      'exactparams  ','multi        ','shapeopts    ','calcopts     ', &
      'animopts     '/)
 character(len=512) :: arg
 integer :: i,j,ieq,ierr

 do i=1,command_argument_count()
    call get_command_argument(i,arg)
    ieq = index(arg,'=')
    if (arg(1:2) /= '--' .or. ieq <= 3) cycle
    ! try each namelist in turn: only the one containing the variable reads cleanly
    do j=1,size(groups)
       call read_namelist_from_string(trim(groups(j)),trim(arg(3:)),ierr)
       ! retry with the value in quotes in case it is a character variable
       if (ierr /= 0 .and. index(arg(ieq+1:),"'")==0) &
          call read_namelist_from_string(trim(groups(j)),arg(3:ieq)//"'"//trim(arg(ieq+1:))//"'",ierr)
       if (ierr == 0) then
          if (iverbose >= 0) print "(a)",' setting '//trim(arg(3:))//' from command line'
          exit
       endif
    enddo
 enddo

end subroutine defaults_read_command_line

!-----------------------------------------------
! read a string of the form 'name=value' into
! the namelist group with the given name
!-----------------------------------------------
subroutine read_namelist_from_string(group,string,ierr)
 use multiplot,          only:multi
 use settings_data,      only:dataopts
 use settings_part,      only:plotopts
 use settings_page,      only:pageopts
 use settings_render,    only:renderopts
 use settings_vecplot,   only:vectoropts
 use settings_xsecrot,   only:xsecrotopts,animopts
 use settings_powerspec, only:powerspecopts
 use exact,              only:exactopts,exactparams
 use shapes,             only:shapeopts
 use calcquantities,     only:calcopts
 character(len=*), intent(in)  :: group,string
 integer,          intent(out) :: ierr
 character(len=len(group)+len(string)+4) :: nmlstring

 nmlstring = '&'//group//' '//string//' /'
 select case(group)
 case('dataopts');      read(nmlstring,nml=dataopts,iostat=ierr)
 case('plotopts');      read(nmlstring,nml=plotopts,iostat=ierr)
 case('pageopts');      read(nmlstring,nml=pageopts,iostat=ierr)
 case('renderopts');    read(nmlstring,nml=renderopts,iostat=ierr)
 case('vectoropts');    read(nmlstring,nml=vectoropts,iostat=ierr)
 case('xsecrotopts');   read(nmlstring,nml=xsecrotopts,iostat=ierr)
 case('powerspecopts'); read(nmlstring,nml=powerspecopts,iostat=ierr)
 case('exactopts');     read(nmlstring,nml=exactopts,iostat=ierr)
 case('exactparams');   read(nmlstring,nml=exactparams,iostat=ierr)
 case('multi');         read(nmlstring,nml=multi,iostat=ierr)
 case('shapeopts');     read(nmlstring,nml=shapeopts,iostat=ierr)
 case('calcopts');      read(nmlstring,nml=calcopts,iostat=ierr)
 case('animopts');      read(nmlstring,nml=animopts,iostat=ierr)
 case default;          ierr = 1
 end select

end subroutine read_namelist_from_string

end module defaults
