!-------------------------------------------------------------------------
! Module containing settings and options related to the plot limits
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_limits
 implicit none
 integer :: itrackpart
 logical :: iadapt, iadaptcoords
 real :: scalemax,zoom
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
  zoom = 1.0           ! for rescaling fixed limits
  itrans(:) = 0        ! no transformations (log10 etc)
  itrackpart = 0       ! particle to track (none)
  xminoffset_track = 0.5 ! offset of limits from tracked particle
  xmaxoffset_track = 0.5 !
 
  return
end subroutine defaults_set_limits

!----------------------------------------------------------------------
! submenu with options relating to plot limits
!----------------------------------------------------------------------
subroutine submenu_limits
 use filenames, only:rootname
 use settings_data, only:nstart,n_end,ndataplots,numplot,ndim,ivegotdata
 !!use settings_page, only:nstepsperpage
 use multiplot, only:itrans
 use prompting
 use limits
 use labels
 use transforms, only:ntrans,transform_label
 implicit none
 integer :: iaction,ipick,i,ierr,index
 real :: diff, mid
 character(len=120) :: transprompt
 character(len=len(rootname)+7) :: limitsfile
 
 index = 1
 do i=1,ntrans
    write(transprompt(index:),"(1x,i1,'=',a,',')") i,trim(transform_label('x',i))
    index = len_trim(transprompt) + 1
 enddo

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
 call prompt('enter option ',iaction,0,8)
!
!--limits
!
 select case(iaction)
!------------------------------------------------------------------------
 case(1)
    call prompt('Use adaptive plot limits?',iadapt)
    call prompt('Use adaptive plot limits on coordinate axes?',iadaptcoords)
    print*,'adaptive plot limits = ',iadapt,' on coords = ',iadaptcoords
    !if (nstepsperpage.gt.1 .and. (iadapt .or. iadaptcoords)) then
    !   print*,'WARNING: adaptive limits and multiple steps per page don''t mix'
    !endif
!------------------------------------------------------------------------
 case(2)
    ipick = 1
    iadapt = .false.
    do while (ipick.gt.0)
       ipick = 0
       !write(*,*)
       call prompt('Enter plot number to set limits (0=quit)',ipick,0,numplot)
       if (ipick.gt.0) then
          call prompt(trim(label(ipick))//' min ',lim(ipick,1))
          call prompt(trim(label(ipick))//' max ',lim(ipick,2))
          print*,'>> '//trim(label(ipick))//' limits set (min,max) = ',lim(ipick,1),lim(ipick,2)
          !print "(a)", trim(transprompt)
          !call prompt('Which transform? (21 for 2 then 1)',itrans(ipick),0)
       endif
    enddo
    return
!------------------------------------------------------------------------
 case(3)
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
!------------------------------------------------------------------------
 case(4)
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
!------------------------------------------------------------------------
  case(5)
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
  case(6)
     limitsfile = trim(rootname(1))//'.limits'
     call prompt('Enter name of limits file to write ',limitsfile)
     !--append .limits if necessary
     !!!if (index(limitsfile,'.limits').eq.0) limitsfile = trim(limitsfile)//'.limits'
     call write_limits(limitsfile)
  case(7)
     limitsfile = trim(rootname(1))//'.limits'
     call prompt('Enter name of limits file to read ',limitsfile)
     !--append .limits if necessary
     !!!if (index(limitsfile,'.limits').eq.0) limitsfile = trim(limitsfile)//'.limits'
     call read_limits(limitsfile,ierr)
  case(8)
     if (ivegotdata) then
        call set_limits(nstart,n_end,1,ndataplots)
     else
        print*,'no data with which to set limits!!'
     endif
  end select
 
 return
end subroutine submenu_limits

end module settings_limits
