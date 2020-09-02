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
   use titles,             only:pagetitles,steplegend
   implicit none
   logical, intent(in) :: use_evdefaults
   integer :: i
  !
  !--set defaults for submenu options
  !
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
      print "(a)",' ** ev mode: using default settings for .ev files **'
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
 use exact,              only:exactopts,exactparams
 use filenames,          only:rootname,nfiles
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
 implicit none
 character(len=*), intent(in) :: filename
 integer :: i,ierr
 integer, parameter :: iunit = 1

 open(unit=iunit,file=trim(adjustl(filename)),status='replace',form='formatted', &
      delim='apostrophe',iostat=ierr) ! without delim namelists may not be readable
 if (ierr /= 0) then
    print*,'ERROR: cannot write file '//trim(filename)
    close(unit=iunit)
    return
 endif
 write(iunit,NML=dataopts)
 write(iunit,NML=plotopts)
 write(iunit,NML=pageopts)
 write(iunit,NML=renderopts)
 write(iunit,NML=vectoropts)
 write(iunit,NML=xsecrotopts)
 write(iunit,NML=powerspecopts)
 write(iunit,NML=exactopts)
 write(iunit,NML=exactparams)
 write(iunit,NML=multi)
 write(iunit,NML=shapeopts)
 write(iunit,NML=calcopts)
 write(iunit,NML=animopts)
 do i=1,nfiles
    write(iunit,"(a)") trim(rootname(i))
 enddo
 close(unit=iunit)
 print "(a)",'default options saved to file '//trim(filename)

 return
end subroutine defaults_write
!-----------------------------------------------
! reads default options from file
! uses namelist input to group the options
! these are specified in the modules
!-----------------------------------------------
subroutine defaults_read(filename)
 use filenames,          only:rootname,maxfile
 use multiplot,          only:multi
 use settings_data,      only:dataopts,idustfrac_plot,iRescale_has_been_set
 use settings_part,      only:plotopts
 use settings_page,      only:pageopts
 use settings_render,    only:renderopts
 use settings_vecplot,   only:vectoropts
 use settings_xsecrot,   only:xsecrotopts,animopts
 use settings_powerspec, only:powerspecopts
 use exact,              only:exactopts,exactparams
 use shapes,             only:shapeopts
 use calcquantities,     only:calcopts
 implicit none
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
    iRescale_has_been_set = .true.
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
    if (nerr > 0) then
       print "(a)",' WARNING: '//trim(filename)//' incomplete (from old code version)'
    else
       print*,'read '//trim(filename)
    endif
    return
 else
    print*,trim(filename)//' not found: using default settings'
    return
 endif

88 continue
 print "(a)",' *** error opening defaults file '//trim(filename)//': using default settings'

 return
end subroutine defaults_read

end module defaults
