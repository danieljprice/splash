!-------------------------------------------------------------------------
! Module containing settings and options related to the plot limits
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_limits
 implicit none
 integer :: itrackpart
 logical :: iadapt, iadaptcoords
 real :: scalemax
 real, dimension(3) :: xminoffset_track, xmaxoffset_track
 
contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_limits
  use multiplot, only:itrans
  implicit none

  iadapt = .true.      ! adaptive plot limits
  iadaptcoords = .false.
  scalemax = 1.0       ! for rescaling adaptive limits
  itrans(:) = 0        ! no transformations (log10 etc)
  itrackpart = 0       ! particle to track (none)
  xminoffset_track = 0.5 ! offset of limits from tracked particle
  xmaxoffset_track = 0.5 !
 
  return
end subroutine defaults_set_limits

!----------------------------------------------------------------------
! submenu with options relating to plot limits
!----------------------------------------------------------------------
subroutine submenu_limits(help)
 use filenames, only:rootname,nsteps,nstepsinfile,ifileopen
 use settings_data, only:ndataplots,numplot,ndim,ivegotdata,DataIsBuffered
 !!use settings_page, only:nstepsperpage
 use multiplot, only:itrans
 use prompting
 use limits
 use labels, only:label,ix
 use transforms, only:ntrans,transform_label
 implicit none
 logical, intent(in), optional :: help
 integer :: iaction,ipick,i,ierr,index
 real :: diff, mid, zoom
 character(len=120) :: transprompt
 character(len=len(rootname)+7) :: limitsfile
 logical :: helpmode
 
 helpmode = .false.
 if (present(help)) helpmode=help
 zoom = 1.0

 iaction = 0
 if (iadapt) then
    print 10,iadapt,iadaptcoords,itrackpart,scalemax
 else
    print 10,iadapt,iadaptcoords,itrackpart,zoom
 endif
10 format(' 0) exit ',/,                 &
        ' 1) set adaptive/fixed limits  ( ',L1,L1,' )   ',/,  &
        ' 2) set manual limits ',/,     &
        ' 3) xy limits track particle      ( ',i8,' )   ',/,   &
        ' 4) zoom in/out                   ( ',f4.2,' ) ',/,   &
        ' 5) apply transformations (log10,1/x) ',/, &
        ' 6) save current limits to file ',/, &
        ' 7) re-read limits file         ',/, &
        ' 8) reset limits for all plots  ')
 if (helpmode) then
    call prompt('choose option for help on a specific item ',iaction,0,8)
 else
    call prompt('enter option ',iaction,0,8)
 endif
!
!--limits
!
 select case(iaction)
!------------------------------------------------------------------------
 case(1)
    if (helpmode) then
       print "(5(/a))",' With limits set to adaptive, plot limits are minimum', &
                   ' and maximum of quantities at current timestep. ', &
                   ' However, the co-ordinate limits are not adapted ', &
                   ' in the case of rendered plots. With fixed limits, the', &
                   ' plot limits retain their default values for all timesteps.'    
    else
       call prompt('Use adaptive plot limits?',iadapt)
       call prompt('Use adaptive plot limits on coordinate axes?',iadaptcoords)
       print*,'adaptive plot limits = ',iadapt,' on coords = ',iadaptcoords
    !if (nstepsperpage.gt.1 .and. (iadapt .or. iadaptcoords)) then
    !   print*,'WARNING: adaptive limits and multiple steps per page don''t mix'
    !endif
    endif
!------------------------------------------------------------------------
 case(2)
    if (helpmode) then
       print "(/a)",' Manually sets the plot limits for each column of data'
    else
       ipick = 1
       do while (ipick.gt.0)
          ipick = 0
          !write(*,*)
          call prompt('Enter plot number to set limits (0=quit)',ipick,0,numplot)
          if (ipick.gt.0) then
             call prompt(trim(label(ipick))//' min ',lim(ipick,1))
             call prompt(trim(label(ipick))//' max ',lim(ipick,2))
             print*,'>> '//trim(label(ipick))//' limits set (min,max) = ',lim(ipick,1),lim(ipick,2)
             if (ipick.le.ndim) then
                iadaptcoords = .false.
             elseif (ipick.le.numplot) then
                iadapt = .false.
             endif
          endif
       enddo
       return
    endif
!------------------------------------------------------------------------
 case(3)
    if (helpmode) then
       print "(4(/a))",'Co-ordinate limits are centred on the selected ', &
                       'particle for all timesteps, with offsets as input',&
                       'by the user.',&
                       'This effectively gives the `Lagrangian'' perspective.'
    else
       call prompt('Enter particle to track: ',itrackpart,0)
       print*,'tracking particle ',itrackpart
       if (itrackpart.gt.0) then
          do i=1,ndim
             call prompt('Enter offset for '//trim(label(ix(i)))//'min:', &
                         xminoffset_track(i))
             call prompt('Enter offset for '//trim(label(ix(i)))//'max :', &
                         xmaxoffset_track(i))
          enddo
       endif
    endif
!------------------------------------------------------------------------
 case(4)
    if (helpmode) then
       print "(/a)",' Zooms in/out (alternatively do this in interactive mode)'
    else
       if (.not.iadapt) then
          call prompt('Enter zoom factor for fixed limits',zoom,0.0)
          do i=1,numplot
             diff = lim(i,2)- lim(i,1)
             mid = 0.5*(lim(i,1) + lim(i,2))
             lim(i,1) = mid - 0.5*zoom*diff
             lim(i,2) = mid + 0.5*zoom*diff
          enddo
       else
          call prompt('Enter scale factor (adaptive limits)',scalemax,0.0)
       endif
    endif
!------------------------------------------------------------------------
  case(5)
     index = 1
     do i=1,ntrans
        write(transprompt(index:),"(1x,i1,'=',a,',')") i,trim(transform_label('x',i))
        index = len_trim(transprompt) + 1
     enddo

     if (helpmode) then
        print "(/a)",' Applies log, inverse and other transformations to data columns'
     else
        ipick = 1
        do while (ipick.gt.0 .and. ipick.le.numplot)
           ipick = 0
           !!print*,'Enter plot number to apply transformation '
           call prompt('Enter column to apply transform (0=quit,-1=all) ',ipick)
           if (ipick.le.numplot .and. ipick.ne.0) then
              print "(a)", trim(transprompt)
              if (ipick.lt.0) then
                 ipick = 0
                 call prompt('Which transform (or multiple e.g. 321)?',ipick,0)
                 itrans(:) = ipick
                 ipick = -99
              else
                 call prompt('Which transform (or multiple e.g. 321)?',itrans(ipick),0)
              endif
           endif
        enddo
        return
     endif
  case(6)
     if (helpmode) then
        print "(9(/a))",' Saves the current values of the fixed plot limits to', &
                     ' file. By default this file is called `filename.limits''',&
                     ' where filename is the name of the *first* data file.', &
                     ' If you use the default name, this file will be automatically ',&
                     ' read upon the next invocation of supersphplot. Using', &
                     ' another filename the limits file can be read by selecting ',&
                     ' option (7) from this menu. Note that limits read from file',&
                     ' will only apply when fixed (ie. not adaptive) limits are used.'  
     else
        limitsfile = trim(rootname(1))//'.limits'
        call prompt('Enter name of limits file to write ',limitsfile)
        !--append .limits if necessary
        !!!if (index(limitsfile,'.limits').eq.0) limitsfile = trim(limitsfile)//'.limits'
        call write_limits(limitsfile)
     endif
  case(7)
     if (helpmode) then
        print "(5(/a))",'Re-reads the plot limits from a file ',&
                    '(see help for write limits file for format)',&
                    'Note that this means that the limits contained in this file',&
                    'can be manually changed by the user whilst the program is ',&
                    'still running.'
     else
        limitsfile = trim(rootname(1))//'.limits'
        call prompt('Enter name of limits file to read ',limitsfile)
        !--append .limits if necessary
        !!!if (index(limitsfile,'.limits').eq.0) limitsfile = trim(limitsfile)//'.limits'
        call read_limits(limitsfile,ierr)
     endif
  case(8)
     if (helpmode) then
        print "(2(/a))",'Resets plot limits using all data currently in memory', &
                    'Note that these limits will only apply when fixed limits are used'
     else
        if (ivegotdata) then
           if (DataIsBuffered) then
              call set_limits(1,nsteps,1,ndataplots)
           else
              call set_limits(1,nstepsinfile(ifileopen),1,ndataplots)
           endif
        else
           print*,'no data with which to set limits!!'
        endif
     endif
  end select
 
 return
end subroutine submenu_limits

end module settings_limits
