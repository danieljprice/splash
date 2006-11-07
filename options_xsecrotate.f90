!-------------------------------------------------------------------------
! Module containing settings and options relating to cross sections,
! rotations and 3D plotting.
! Includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_xsecrot
 implicit none
 integer :: nxsec,irotateaxes
 logical :: xsec_nomulti, irotate, flythru, use3Dperspective, use3Dopacityrendering
 real :: anglex, angley, anglez, zobserver, dzscreenfromobserver
 real :: taupartdepth, rkappa
 real :: xsecpos_nomulti,xseclineX1,xseclineX2,xseclineY1,xseclineY2
 real, dimension(3) :: xorigin,xminrotaxes,xmaxrotaxes

 namelist /xsecrotopts/ xsec_nomulti,xsecpos_nomulti,flythru, &
          xseclineX1,xseclineX2,xseclineY1,xseclineY2, &
          irotate,irotateaxes,anglex, angley, anglez, &
          xminrotaxes,xmaxrotaxes,use3Dperspective, &
          use3Dopacityrendering,zobserver,dzscreenfromobserver, &
          taupartdepth

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_xsecrotate
  implicit none

  xsec_nomulti = .false.    ! take cross section of data / particles
  xsecpos_nomulti = 0.      ! position of cross section
  flythru = .false.         ! take series of cross sections through data
  xseclineX1 = 0.0
  xseclineX2 = 0.0
  xseclineY1 = 0.0
  xseclineY2 = 0.0
  irotate = .false.
  irotateaxes = 0
  anglex = 0.
  angley = 0.
  anglez = 0.
  xorigin = 0.
  xminrotaxes = 0.
  xmaxrotaxes = 0.
  use3Dperspective = .false.
  use3Dopacityrendering = .false.
  zobserver = 0.
  dzscreenfromobserver = 0.
  rkappa = 0. ! rkappa is set from taupartdepth later
  taupartdepth = 2.

  return
end subroutine defaults_set_xsecrotate

!----------------------------------------------------------------------
! sets options relating to cross sectioning / rotation
!----------------------------------------------------------------------
subroutine submenu_xsecrotate
 use labels, only:label,ix
 use limits, only:lim
 use prompting
 use settings_data, only:ndim
 implicit none
 integer :: ians,i
 character(len=1) :: char
 character(len=4) :: text
 logical :: interact
 
 if (ndim.eq.1) print*,' WARNING: none of these options have any effect in 1D'
 ians = 0
 interact = .true.
 if (xsec_nomulti) then
    text = 'xsec'
 else
    text = 'proj'
 endif
 print 10,text,xsecpos_nomulti,print_logical(irotate), &
          print_logical(use3Dperspective),print_logical(use3Dopacityrendering), &
          irotateaxes
10  format('---------- cross section / 3D plotting options --------',/, &
           ' 0) exit ',/,       &
           ' 1) switch between cross section/projection     ( ',a4,' )',/, &
           ' 2) set cross section position                  (',f5.2,' )',/, &
           ' 3) rotation settings/options                   ( ',a,' )',/, &
           ' 4) 3D perspective on/off                       ( ',a,' )',/, &
           ' 5) 3D surface rendering on/off/options         ( ',a,' )',/, &
           ' 6) set axes for rotated/3D plots               ( ',i2,' )')
 call prompt('enter option',ians,0,6)
!
!--options
!
 select case(ians)
!------------------------------------------------------------------------
 case(1)
    xsec_nomulti = .not.xsec_nomulti 
    print *,' Cross section = ',xsec_nomulti
!------------------------------------------------------------------------
 case(2)
    flythru = .false.
    if (ndim.eq.3) then
       call prompt('Do you want a fly-through',flythru)
       if (.not.flythru) then
          call prompt('enter co-ordinate location of cross section slice', &
               xsecpos_nomulti)
       endif
    elseif (ndim.eq.2) then
       call prompt('set cross section position interactively?',interact)
       
       if (interact) then
       !
       !--set cross section position interactively
       !
          call pgbegin(0,'/xw',1,1)
          call pgenv(lim(1,1),lim(1,2),lim(2,1),lim(2,2),1,0)
          call pgcurs(xseclineX1,xseclineY1,char)
          print*,'please select cross section line'
          call pgband(1,1,xseclineX1,xseclineY1,xseclineX2,xseclineY2,char)
          print*,'cross section line: xmin = ',xseclineX1,' xmax = ',xseclineX2
          print*,'                    ymin = ',xseclineY1,' ymax = ',xseclineY2
          call pgend       
       else
       !
       !--set position manually
       !
          if (abs(xseclineX2-xseclineX1).lt.1.e-5 .and. &
              abs(xseclineY2-xseclineY1).lt.1.e-5) then
          !--if not already set (ie. if all = 0.0)
          !  then set default line to diagonal across the domain
             xseclineX1 = lim(1,1)
             xseclineX2 = lim(1,2)
             xseclineY1 = lim(2,1)
             xseclineY2 = lim(2,2)
          endif
          print*,'please set position of cross section through 2D data:'
          call prompt('enter xmin of cross section line',xseclineX1)
          call prompt('enter xmax of cross section line',xseclineX2)
          call prompt('enter ymin of cross section line',xseclineY1)
          call prompt('enter ymax of cross section line',xseclineY2)
       endif
    else
       print*,'WARNING: this option has no effect in 1D'
    endif
!------------------------------------------------------------------------
 case(3)
    call prompt('use rotation?',irotate)
    print*,'rotate = ',irotate
    if (irotate) then
       print*,'note that rotations are done in the order z-y-x '
       print*,'this means the y and x rotations are done about the *new* y and x axes'
       print*,'if in doubt, set the angles interactively in this order'
       call prompt('enter rotation angle about z axis (deg)',anglez,0.,360.)
       if (ndim.eq.3) then
          call prompt('enter rotation angle about y axis (deg)',angley,0.,360.)
          call prompt('enter rotation angle about x axis (deg)',anglex,0.,360.)
       endif
       !xorigin(1:ndim) = 0.5*(lim(1:ndim,1) + lim(1:ndim,2))
       do i=1,ndim
          call prompt('enter location of origin '//trim(label(ix(i))),xorigin(i))
       enddo
    endif
!------------------------------------------------------------------------
 case(4)
    use3Dperspective = .not.use3Dperspective
    print "(a,L1)",' 3D perspective = ',use3Dperspective
    if (.not.use3Dperspective) use3Dopacityrendering = .false.
!------------------------------------------------------------------------
 case(5)
    use3Dopacityrendering = .not.use3Dopacityrendering
    print "(a,L1)",' 3D opacity rendering = ',use3Dopacityrendering
    if (use3Dopacityrendering .and..not.use3Dperspective) then
       print "(a)",' also turning on 3D perspective (which must be set for this to work)'
       use3Dperspective = .true.
    endif
!------------------------------------------------------------------------
 case(6)
    print*,'0 : do not plot rotated axes'
    print*,'1 : plot rotated axes'
    print*,'2 : plot rotated box'
    print*,'3 : plot gridded x-y plane'
    call prompt('enter type of axes to plot',irotateaxes,0,3)
    if (irotateaxes.gt.0) then
       !--if not previously set, use current plot limits
       if (all(abs(xminrotaxes).le.tiny(xminrotaxes))) then
          xminrotaxes(:) = lim(ix(:),1)
          xmaxrotaxes(:) = lim(ix(:),2)
       endif
       do i=1,ndim
          call prompt('enter '//trim(label(ix(i)))//'min:',xminrotaxes(i))
          call prompt('enter '//trim(label(ix(i)))//'max:',xmaxrotaxes(i))
       enddo
    endif
 end select

 return
end subroutine submenu_xsecrotate

end module settings_xsecrot
