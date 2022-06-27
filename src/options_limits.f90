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

!-------------------------------------------------------------------------
! Module containing settings and options related to the plot limits
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_limits
 implicit none
 logical :: iadapt, iadaptcoords, adjustlimitstodevice
 real    :: scalemax
 real, dimension(3) :: xminoffset_track, xmaxoffset_track

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_limits
 use multiplot, only:itrans

 iadapt           = .true.  ! adaptive plot limits
 iadaptcoords     = .false.
 adjustlimitstodevice = .false.
 scalemax         = 1.0     ! for rescaling adaptive limits
 itrans(:)        = 0       ! no transformations (log10 etc)
 xminoffset_track = 0.5     ! offset of limits from tracked particle
 xmaxoffset_track = 0.5     !

end subroutine defaults_set_limits

!----------------------------------------------------------------------
! submenu with options relating to plot limits
!----------------------------------------------------------------------
subroutine submenu_limits(ichoose)
 use filenames,      only:nsteps,nstepsinfile,ifileopen
 use settings_data,  only:ndataplots,numplot,ndim,iCalcQuantities, &
                          DataIsBuffered,itracktype,itrackoffset,ntypes
 use calcquantities, only:calc_quantities
 use multiplot,      only:itrans
 use prompting,      only:prompt,print_logical
 use limits,         only:lim,range,rangeset,anyrangeset,print_rangeinfo
 use labels,         only:label,ix,irad,is_coord,labeltype
 use transforms,     only:ntrans,transform_label
 integer, intent(in) :: ichoose
 integer             :: iaction,ipick,i,index,ierr
 integer             :: itracktypeprev,itrackoffsetprev
 character(len=120)  :: transprompt,fmtstring
 character(len=12)   :: string,string2
 character(len=20)   :: pstring,pstring2

! zoom = 1.0

 iaction = ichoose
 if (iadapt) then
    string = 'ADAPT'
 else
    string = 'FIXED'
 endif
 if (iadaptcoords) then
    string2 = 'ADAPT'
 else
    string2 = 'FIXED'
 endif
 write(pstring,"(i12)") itrackoffset
 pstring = adjustl(pstring)
 if (itracktype > 0) then
    write(pstring2,"(i12)") itracktype
    pstring=trim(adjustl(pstring2))//':'//trim(pstring)
 else
    pstring=trim(pstring)
 endif

 print "(a)",'------------------ limits options ---------------------'
10 format( &
        ' 0) exit ',/,                 &
        ' 1) use adaptive/fixed limits        ( ',a,', ',a,' )   ',/,  &
        ' 2) set limits manually ',/,     &
        ' 3) set particle tracking                      ( ',a,' )',/,   &
        ' 4) subset data by parameter range             ( ',a,')',/, &
        ' 5) apply log/other transformations to columns ')
 if (iaction <= 0 .or. iaction > 5) then
    print 10,trim(string),trim(string2),trim(pstring),print_logical(anyrangeset())
    call prompt('enter option ',iaction,0,5)
 endif
!
!--limits
!
 select case(iaction)
!------------------------------------------------------------------------
 case(1)

!+ With limits set to adaptive, plot limits are minimum
!+ and maximum of quantities at current timestep.
!+ However, the co-ordinate limits are not adapted
!+ in the case of rendered plots. With fixed limits, the
!+ plot limits retain their default values for all timesteps.

    call prompt('Use adaptive plot limits?',iadapt)
    call prompt('Use adaptive plot limits on coordinate axes?',iadaptcoords)
    call prompt('Adjust limits to aspect ratio of device?',adjustlimitstodevice)
    print "(a)",'adaptive plot limits = '//print_logical(iadapt)// &
                ' on coords = '//print_logical(iadaptcoords)

!------------------------------------------------------------------------
 case(2)

!+ Manually sets the plot limits for each column of data

    ipick = 1
    do while (ipick > 0)
       ipick = 0
       !write(*,*)
       if (ndim >= 1) then
          write(fmtstring,'(a,i3,a)') 'Enter column number to set limits (0=quit;',numplot+1,'=centred cube in xyz)'
          call prompt(trim(fmtstring),ipick,0,numplot+1)
       else
          write(fmtstring,'(a,i3,a)') 'Enter column number to set limits (0=quit)'
          call prompt(trim(fmtstring),ipick,0,numplot)
       endif
       if (ipick==numplot+1 .and. ndim >= 1) then
          call prompt(trim(label(ix(1)))//' max ',lim(ix(1),2))
          lim(ix(1),1) = -lim(ix(1),2)
          print*,'>> '//trim(label(ix(1)))//' limits set (min,max) = ',lim(ix(1),:)
          if (ndim >= 2) then
             lim(ix(2:ndim),1) = lim(ix(1),1)
             lim(ix(2:ndim),2) = lim(ix(1),2)
             print*,'>> '//trim(label(ix(2)))//' limits set (min,max) = ',lim(ix(2),:)
             if (ndim >= 3) print*,'>> '//trim(label(ix(3)))//' limits set (min,max) = ',lim(ix(3),:)
          endif
       elseif (ipick > 0) then
          call prompt(trim(label(ipick))//' min ',lim(ipick,1))
          call prompt(trim(label(ipick))//' max ',lim(ipick,2))
          print*,'>> '//trim(label(ipick))//' limits set (min,max) = ',lim(ipick,1),lim(ipick,2)
          if (is_coord(ipick,ndim)) then
             iadaptcoords = .false.
          elseif (ipick <= numplot) then
             iadapt = .false.
          endif
       endif
    enddo
!------------------------------------------------------------------------
 case(3)

!+ Co-ordinate limits are centred on the selected
!+ particle for all timesteps, with offsets as input by the user.
!+ This effectively gives the `Lagrangian' perspective.

    itrackoffsetprev = itrackoffset
    itracktypeprev   = itracktype
    print "(a,/,a,/)",'To track particle 4923, enter 4923', &
                      'To track the 43rd particle of type 3, enter 3:43'

    call prompt('Enter particle to track: ',pstring,noblank=.true.)
    call get_itrackpart(pstring,itracktype,itrackoffset,ierr)
    do while (ierr /= 0 .or. itracktype < 0 .or. itracktype > ntypes .or. itrackoffset < 0)
       if (itracktype < 0 .or. itracktype > ntypes) print "(a)",'invalid particle type'
       if (itrackoffset < 0)                         print "(a)",'invalid particle index'
       if (ierr /= 0)                                 print "(a)",'syntax error in string'
       pstring = '0'
       call prompt('Enter particle to track: ',pstring,noblank=.true.)
       call get_itrackpart(pstring,itracktype,itrackoffset,ierr)
    enddo
    if (itracktype > 0) then
       write(string,"(i12)") itrackoffset
       string = adjustl(string)
       print "(a)",'=> tracking '//trim(labeltype(itracktype))//' particle #'//trim(string)
    elseif (itrackoffset > 0) then
       write(string,"(i12)") itrackoffset
       string = adjustl(string)
       print "(a)",'=> tracking particle '//trim(string)
    else
       print "(a)",'=> particle tracking limits OFF'
    endif

    if (itrackoffset > 0) then
       do i=1,ndim
          call prompt('Enter offset for '//trim(label(ix(i)))//'min:', &
                      xminoffset_track(i))
          call prompt('Enter offset for '//trim(label(ix(i)))//'max :', &
                      xmaxoffset_track(i))
       enddo
       if ((itrackoffset /= itrackoffsetprev .or. itracktype /= itracktypeprev) &
           .and. iCalcQuantities .and. irad > 0 .and. irad <= numplot) then
          !--radius calculation is relative to tracked particle
          print "(a)",' recalculating radius relative to tracked particle '
          if (DataIsBuffered) then
             call calc_quantities(1,nsteps)
          else
             call calc_quantities(1,nstepsinfile(ifileopen))
          endif
       endif
    endif
!------------------------------------------------------------------------
 case(4)

!+ Plot subset of data by restricting parameter range

    ipick = 1
    do while (ipick > 0)
       ipick = 0
       call print_rangeinfo()

       call prompt('Enter column to use to restrict data set (-1=none/unset all,0=quit)',ipick,-1,ndataplots)
       if (ipick > 0) then
          print*,'current plot limits for '//trim(label(ipick))//': (min,max) = ',lim(ipick,1),lim(ipick,2)
          call prompt(trim(label(ipick))//' min value ',range(ipick,1))
          call prompt(trim(label(ipick))//' max value ',range(ipick,2),range(ipick,1))
          if (.not.rangeset(ipick)) then
             print*,'>> min=max: no restriction set'
          endif
       elseif (ipick==-1) then
          print "(a)",'>> removing all range restrictions on data set'
          range(:,:) = 0.
       endif
       write(*,*)
    enddo
!------------------------------------------------------------------------
 case(5)

!+ Applies log, inverse and other transformations to data columns

    index = 1
    do i=1,ntrans
       write(transprompt(index:),"(1x,i1,'=',a,',')") i,trim(transform_label('x',i))
       index = len_trim(transprompt) + 1
    enddo

    ipick = 1
    do while (ipick > 0 .and. ipick <= numplot)
       ipick = 0
       call prompt('Enter column to apply transform (0=quit,-1=all) ',ipick)
       if (ipick <= numplot .and. ipick /= 0) then
          print "(a)", trim(transprompt)
          if (ipick < 0) then
             ipick = 0
             call prompt('Which transform (or multiple e.g. 321)?',ipick,0)
             itrans(:) = ipick
             ipick = -99
          else
             call prompt('Which transform (or multiple e.g. 321)?',itrans(ipick),0)
          endif
       endif
    enddo

 end select

end subroutine submenu_limits

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
 endif

end subroutine get_itrackpart

end module settings_limits
