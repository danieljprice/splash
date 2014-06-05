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
!  Copyright (C) 2005-2012 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! Module containing settings and options relating to cross sections,
! rotations and 3D plotting.
! Includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_xsecrot
 implicit none
 !--public variables
 integer, public :: nframes,nseq
 integer, public :: nxsec,irotateaxes
 logical, public :: xsec_nomulti, irotate, flythru, use3Dperspective, use3Dopacityrendering
 logical, public :: writeppm, rendersinks
 real, public    :: anglex, angley, anglez, zobserver, dzscreenfromobserver
 real, public    :: taupartdepth,xsecwidth
 real, public    :: xsecpos_nomulti,xseclineX1,xseclineX2,xseclineY1,xseclineY2
 real, public, dimension(3) :: xorigin,xminrotaxes,xmaxrotaxes

 !--private variables related to animation sequences
 integer, parameter, private         :: maxseq = 6
 integer, dimension(maxseq), public :: iseqstart,iseqend,iseqtype
 integer, public :: icolchange
 real, public    :: xminseqend,xmaxseqend,yminseqend,ymaxseqend
 real, public    :: anglezend,angleyend,anglexend,zobserverend,taupartdepthend
 real, public    :: xmincolend,xmaxcolend,xsecpos_nomulti_end
 logical, private :: ihavesetsequence
 character(len=*), dimension(maxseq), parameter, private :: labelseqtype = &
    (/'steady zoom on x and y axes                       ', &
      'steady rotation                                   ', &
      'steady change of limits (e.g. for colour bar)     ', &
      'steady movement of 3D observer                    ', &
      'sequence of cross section slices through a 3D box ', &
      'steady change of opacity for 3D surface plots     '/)

 !--namelists for writing to defaults file and .anim file
 public :: xsecrotopts
 namelist /xsecrotopts/ xsec_nomulti,xsecpos_nomulti,flythru, &
          xseclineX1,xseclineX2,xseclineY1,xseclineY2, &
          irotate,irotateaxes,anglex, angley, anglez, &
          xminrotaxes,xmaxrotaxes,use3Dperspective, &
          use3Dopacityrendering,zobserver,dzscreenfromobserver, &
          taupartdepth,writeppm,xsecwidth,rendersinks

 public :: animopts
 namelist /animopts/ nseq,nframes,iseqstart,iseqend,iseqtype, &
          xminseqend,xmaxseqend,yminseqend,ymaxseqend, &
          anglezend,angleyend,anglexend,zobserverend,taupartdepthend, &
          icolchange,xmincolend,xmaxcolend,xsecpos_nomulti_end

 !--public procedure names
 public :: defaults_set_xsecrotate,submenu_xsecrotate,getsequencepos,insidesequence
 public :: setsequenceend
 
 procedure(add_sequence), pointer :: addseq => null()
 procedure(delete_sequence), pointer :: delseq => null()
 procedure(check_sequences), pointer :: checkseq => null()

 private

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
  xsecwidth = 0.0   ! width of xsec slices - zero means suggest better value to user
  irotate = .false.
  irotateaxes = 0
  anglex = 0.
  angley = 0.
  anglez = 0.
  xminrotaxes = 0.
  xmaxrotaxes = 0.
  use3Dperspective = .false.
  use3Dopacityrendering = .false.
  zobserver = 0.
  dzscreenfromobserver = 0.
  taupartdepth = 2.
  writeppm = .true.

  !--defaults for animation sequences
  nseq = 0
  nframes = 0
  iseqstart(:) = 0
  iseqend(:)   = 0
  iseqtype(:)  = 0
  xminseqend = 0.
  xmaxseqend = 0.
  yminseqend = 0.
  ymaxseqend = 0.
  anglezend = 360.
  angleyend = 0.
  anglexend = 0.
  icolchange = 0
  xmincolend = 0.
  xmaxcolend = 0.
  zobserverend = 0.
  taupartdepthend = 2000.0
  xsecpos_nomulti_end = 0.
  ihavesetsequence = .false.
  rendersinks = .false.

  return
end subroutine defaults_set_xsecrotate

!----------------------------------------------------------------------
! sets options relating to cross sectioning / rotation
!----------------------------------------------------------------------
subroutine submenu_xsecrotate(ichoose)
 use filenames,      only:nsteps,nstepsinfile,ifileopen
 use labels,         only:label,ix,irad,get_sink_type
 use limits,         only:lim
 use prompting,      only:prompt,print_logical
 use promptlist,     only:prompt_list
 use settings_data,  only:ndim,xorigin,iCalcQuantities,DataIsBuffered,ntypes
 use calcquantities, only:calc_quantities
 use plotlib,        only:plotlib_supports_alpha
 implicit none
 integer, intent(in) :: ichoose
 integer :: ians,i
 logical :: ichangedorigin
 character(len=4) :: text
 real, dimension(3) :: xorigintemp

 print "(a)",'---------- cross section / 3D plotting options --------'
 if (ndim.eq.1) print*,' WARNING: none of these options have any effect in 1D'
 ians = ichoose
 if (xsec_nomulti) then
    text = 'xsec'
 else
    text = 'proj'
 endif

 if (ians.le.0 .or. ians.gt.6) then
    print 10,text,print_logical(irotate),anglex,angley,anglez, &
             print_logical(use3Dperspective),print_logical(use3Dopacityrendering), &
             irotateaxes,nseq
10  format( &
              ' 0) exit ',/,       &
              ' 1) switch between cross section/projection      ( ',a4,' )',/, &
              ' 2) rotation on/off/settings (incl. origin pos)  ( ',a,3(1x,f5.1),' )',/, &
              ' 3) 3D perspective on/off                        ( ',a,' )',/, &
              ' 4) 3D surface rendering on/off                  ( ',a,' )',/, &
              ' 5) set axes for rotated/3D plots                ( ',i2,' )',/, &
              ' 6) set animation sequence (rotate,flythru etc.) ( ',i2,' )')
    call prompt('enter option',ians,0,6)
 endif
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
    call prompt('use rotation?',irotate)
    print "(a)",' rotation is '//trim(print_logical(irotate))
    if (irotate) then
       print*,'note that rotations are done in the order z-y-x '
       print*,'this means the y and x rotations are done about the *new* y and x axes'
       print*,'if in doubt, set the angles interactively in this order'
       call prompt('enter rotation angle about z axis (deg)',anglez,0.,360.)
       if (ndim.eq.3) then
          call prompt('enter rotation angle about y axis (deg)',angley,0.,360.)
          call prompt('enter rotation angle about x axis (deg)',anglex,0.,360.)
       endif
    endif

    !xorigin(1:ndim) = 0.5*(lim(1:ndim,1) + lim(1:ndim,2))
    xorigintemp(1:ndim) = xorigin(1:ndim)
    ichangedorigin = .false.
    print "(a)",' Note that origin settings affect both rotation and radius calculations'
    do i=1,ndim
       call prompt('enter location of origin '//trim(label(ix(i))),xorigin(i))
       if (abs(xorigin(i)-xorigintemp(i)).gt.tiny(0.)) then
          ichangedorigin = .true.
       endif
    enddo
    !--recalculate radius if origin settings have changed
    if (ichangedorigin .and. iCalcQuantities .and. irad.gt.0) then
       if (DataIsBuffered) then
          call calc_quantities(1,nsteps)
       else
          call calc_quantities(1,nstepsinfile(ifileopen))
       endif
    endif
!------------------------------------------------------------------------
 case(3)
    call prompt(' Use 3D perspective? ',use3Dperspective)
    if (use3Dperspective) then
       if (.not.irotate) irotate = .true.
       if (abs(anglez).lt.tiny(anglez) .and. &
           abs(angley).lt.tiny(angley) .and. &
           abs(anglex).lt.tiny(anglex)) then
           anglez = 30.
           anglex = 60.
           print "(a)",' setting default rotation angles to 30,0,60'
       endif
    else ! turn off opacity rendering if 3D perspective has been turned off
       use3Dopacityrendering = .false.
    endif
!------------------------------------------------------------------------
 case(4)
    call prompt(' Use 3D opacity rendering? ',use3Dopacityrendering)
    if (use3Dopacityrendering .and..not.use3Dperspective) then
       print "(a)",' also turning on 3D perspective (which must be set for this to work)'
       use3Dperspective = .true.
    endif
    if (use3Dopacityrendering) then
       if (.not.plotlib_supports_alpha) then
          print "(/,a)",' Warning: 3D opacity rendering sends only an approximate version '
          print "(a,/)",' to the PGPLOT device (not corrected for brightness) '
          call prompt(' Do you want to write a ppm file in addition to PGPLOT output?',writeppm)
       else
          writeppm = .false.
          !call prompt(' Do you want to apply the brightness correction?',writeppm)
       endif
    endif
    if (use3Dopacityrendering .and. get_sink_type(ntypes) > 0) then
       call prompt('Include sinks in opacity rendering (no=plot on top)?',rendersinks)
    endif
!------------------------------------------------------------------------
 case(5)
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
!------------------------------------------------------------------------
 case(6)
    call submenu_animseq()

 end select

 return
end subroutine submenu_xsecrotate

!----------------------------------------------------------------------
! sets up animation sequences
!----------------------------------------------------------------------
subroutine submenu_animseq()
 use promptlist, only:prompt_list
 use prompting,  only:prompt

 addseq => add_sequence
 checkseq => check_sequences
 delseq => delete_sequence
 call prompt_list(nseq,maxseq,'sequence',checkseq,addseq,delseq)

end subroutine submenu_animseq

!-----------------------------------
! print the current list of shapes
!-----------------------------------
subroutine check_sequences(n)
 implicit none
 integer, intent(in) :: n
 integer :: iseq

 print "(/,a)", ' Current list of animation sequences:'
 if (n.gt.0) then
    do iseq=1,n
       print "(i2,') ',a)",iseq,labelseqtype(iseqtype(iseq))
    enddo
 else
    print "(a)",' (none)'
 endif

end subroutine check_sequences

subroutine delete_sequence(iseq,n)
 implicit none
 integer, intent(in)    :: iseq
 integer, intent(inout) :: n
 integer :: i
 
 if (iseq.gt.0 .and. n.gt.0 .and. iseq.le.maxseq) then
    if (iseqtype(iseq).gt.0 .and. iseqtype(iseq).le.maxseq) then
       print "(a,i1,': ',a)",' deleting sequence ',iseq,trim(labelseqtype(iseqtype(iseq)))
    endif
    iseqtype(iseq) = 0
    do i=iseq+1,n
       iseqtype(i-1) = iseqtype(i)
    enddo
    n = n - 1
 endif

end subroutine delete_sequence

subroutine add_sequence(istart,iend,n)
 use prompting,     only:prompt
 use limits,        only:lim
 use labels,        only:ix,irho
 use settings_data, only:ndim,istartatstep,iendatstep,numplot
 use filenames,     only:nsteps
 implicit none
 integer, intent(in)    :: istart,iend
 integer, intent(inout) :: n
 integer :: i,j,ierr

 i = istart + 1
 over_sequences: do while (i.le.iend .and. i.le.maxseq)

    if (i.gt.n) n = i
    if (n.gt.0) then
       !--set sensible default value for number of frames
       if (nframes.eq.0) then
          if (nsteps.gt.1) then
             nframes = 1
          else
             nframes = 10
          endif
       endif
       call prompt('Enter number of frames generated between dumps (applies to all sequences)',nframes,1,500)
       !call prompt('Use same sequence position for all plots on the page?',imultiframeseq)
    endif

    print "(a,i2,a)",'----------------- sequence ',i,' ----------------------'
    if (iseqstart(i).eq.0) iseqstart(i) = max(istartatstep,1)
    if (iseqend(i).eq.0) iseqend(i) = max(1,iendatstep,istartatstep)
    if (nsteps.gt.1) then
       call prompt('Enter starting dump for sequence ',iseqstart(i),1,nsteps)
       call prompt('Enter finishing dump for sequence ',iseqend(i),1,nsteps)
    endif

    ierr = 1
    do while (ierr /= 0)
       print "(7(/,1x,i1,1x,':',1x,a))",0,'none (remove sequence) ', &
                                       (j,labelseqtype(j),j=1,maxseq)

       call prompt('Enter type of sequence ',iseqtype(i),0,maxseq)
       !--allow only one sequence of each type
       ierr = 0
       if (i.gt.0) then
          do j=1,n
             if (i.ne.j .and. (iseqtype(j).eq.iseqtype(i)) .and. (iseqtype(i).gt.0)) ierr = 2
          enddo
          if (ierr.eq.2) print "(/,a)",' Error: can only have one sequence of each type '
       endif
    end do
    select case(iseqtype(i))
    case(1)
       print "(a)",'Note: zoom sequence starts using current fixed x,y plot limits'
       if (abs(xminseqend).lt.tiny(xminseqend) .and. abs(xmaxseqend).lt.tiny(xmaxseqend)) then
          xminseqend = lim(1,1)
          xmaxseqend = lim(1,2)
       endif
       call prompt(' Enter finishing xmin ',xminseqend)
       call prompt(' Enter finishing xmax ',xmaxseqend)
       if (abs(yminseqend).lt.tiny(yminseqend) .and. abs(ymaxseqend).lt.tiny(ymaxseqend)) then
          yminseqend = lim(2,1)
          ymaxseqend = lim(2,2)
       endif
       call prompt(' Enter finishing ymin ',yminseqend)
       call prompt(' Enter finishing ymax ',ymaxseqend)
    case(2)
       if (ndim.lt.2) then
          print "(a)",' ERROR: cannot use this sequence in 1D'
          iseqtype(i) = 0
       endif
       if (.not.irotate) then
          print "(a)",' Turning rotation on...'
          irotate = .true.
       endif
       print "(a)",'Note: rotation sequence starts using current rotation settings'
       call prompt(' Enter finishing rotation angle (z axis) ',anglezend)
       call prompt(' Enter finishing rotation angle (y axis) ',angleyend)
       call prompt(' Enter finishing rotation angle (x axis) ',anglexend)
    case(3)
       if (icolchange.le.0 .or. icolchange.gt.numplot) then
          if (irho.gt.0 .and. irho.le.numplot) then
             icolchange = irho
          else
             icolchange = 1
          endif
       endif
       call prompt(' Enter column to change limits ',icolchange,1,numplot)
       print "(a)",'Note: limits start from current fixed plot limits for this column'
       if (abs(xmincolend).lt.tiny(xmincolend) .and. abs(xmaxcolend).lt.tiny(xmaxcolend)) then
          xmincolend = lim(icolchange,1)
          xmaxcolend = lim(icolchange,2)
       endif
       call prompt(' Enter finishing minimum value ',xmincolend)
       call prompt(' Enter finishing maximum value ',xmaxcolend)
    case(4)
       if (ndim.ne.3) then
          print "(a)",' ERROR: cannot use this sequence in < 3D'
          iseqtype(i) = 0
       endif
       if (.not.use3Dperspective) then
          print "(a)",'Turning 3D perspective on...'
          use3Dperspective = .true.
       endif
       print "(a)",'Note: observer starts at current observer settings '
       print "(a)",'      (screen height does not change)'
       !--try to give sensible default values
       if (abs(zobserverend).lt.tiny(zobserverend)) then
          if (abs(zobserver).gt.tiny(zobserver)) then
             zobserverend = 5.*zobserver
          elseif (ix(3).gt.0 .and. ix(3).le.numplot) then
             zobserverend = 10.*lim(ix(3),2)
          endif
       endif
       call prompt(' Enter finishing 3D observer height ',zobserverend)
    case(5)
       if (ndim.ne.3) then
          print "(a)",' ERROR: cannot use this sequence in < 3D'
          iseqtype(i) = 0
       endif
       if (.not.xsec_nomulti) then
          print "(a)",'Changing from projection to cross-section'
          xsec_nomulti = .true.
          if (use3Dperspective .and. .not.use3Dopacityrendering) then
             print "(a)",'Turning 3D perspecitve off'
             use3Dperspective = .false.
          endif
       endif
       print "(a)",'Note: slice position starts from value set at initial prompt'
       call prompt(' Enter finishing slice position ',xsecpos_nomulti_end)
    case(6)
       if (ndim.ne.3) then
          print "(a)",' ERROR: cannot use this sequence in < 3D'
          iseqtype(i) = 0
       endif
       if (.not.use3Dperspective .or. .not.use3Dopacityrendering) then
          print "(a)",'Turning 3D opacity rendering and 3D perspective on...'
          use3Dopacityrendering = .true.
          use3Dperspective = .true.
       endif
       print "(3(a,/))",'Note: opacity sequence starts from current opacity value ', &
                      '      and that logarithmic steps are used if finishing value is', &
                      '      set to more than 1000 times the starting value (or vice-versa) '
       call prompt('Enter finishing opacity in units of average smoothing length ',taupartdepthend)
    case default
       call delete_sequence(i,n)
       exit over_sequences
    end select
    i = i + 1
 enddo over_sequences

 if (all(iseqtype(1:n).eq.0)) then
    n = 0
 else
    ihavesetsequence = .true.
 endif

 return
end subroutine add_sequence

!----------------------------------------------------------------------
!
!  subroutine called from interactive mode which sets the current
!  plot settings as the end point to an animation sequence
!
!----------------------------------------------------------------------
subroutine setsequenceend(ipos,iplotx,iploty,irender,rotation, &
                          anglexi,angleyi,anglezi,zobserveri,use3Dopacity,taupartdepthi, &
                          x_sec,xsecposi,xmin,xmax,ymin,ymax,rendermin,rendermax)
 use limits,        only:lim
 use multiplot,     only:itrans
 use settings_data, only:ndim,numplot
 use transforms,    only:transform_limits,transform_limits_inverse,transform_label
 implicit none
 integer, intent(in) :: ipos,iplotx,iploty,irender
 real, intent(in)    :: anglexi,angleyi,anglezi,zobserveri,taupartdepthi,xsecposi
 real, intent(in)    :: xmin,xmax,ymin,ymax,rendermin,rendermax
 logical, intent(in) :: rotation, use3Dopacity,x_sec
 integer :: i
 real    :: xminfixed,xmaxfixed,yminfixed,ymaxfixed,renderminfixed,rendermaxfixed

 nseq = 0
 iseqtype(:) = 0
!
!--compare transformed limits
!
 xminfixed = lim(iplotx,1)
 xmaxfixed = lim(iplotx,2)
 call transform_limits(xminfixed,xmaxfixed,itrans(iplotx))

 yminfixed = lim(iploty,1)
 ymaxfixed = lim(iploty,2)
 call transform_limits(yminfixed,ymaxfixed,itrans(iploty))

 if (irender.gt.0 .and. irender.le.numplot) then
    renderminfixed = lim(irender,1)
    rendermaxfixed = lim(irender,2)
    call transform_limits(renderminfixed,rendermaxfixed,itrans(irender))
 endif

!--set however many sequences are required to capture the change in parameters
!
 !--change of x-y limits
 if (  notequal(xmin,xminfixed) .or. notequal(xmax,xmaxfixed) &
   .or.notequal(ymin,yminfixed) .or. notequal(ymax,ymaxfixed)) then
    nseq = nseq + 1
    iseqtype(nseq) = 1
    xminseqend = xmin
    xmaxseqend = xmax
    yminseqend = ymin
    ymaxseqend = ymax
    print*,trim(transform_label('xmin,max',itrans(iplotx)))//' start = ',xminfixed,xmaxfixed, &
          ' end = ',xminseqend,xmaxseqend
    print*,trim(transform_label('ymin,max',itrans(iploty)))//' start = ',yminfixed,ymaxfixed, &
          ' end = ',yminseqend,ymaxseqend
    !--always store untransformed limits
    call transform_limits_inverse(xminseqend,xmaxseqend,itrans(iplotx))
    call transform_limits_inverse(yminseqend,ymaxseqend,itrans(iploty))
 endif
 !--change of rotation angles
 if (ndim.ge.2 .and. rotation .and. &
    (notequal(anglexi,anglex).or.notequal(angleyi,angley).or.notequal(anglezi,anglez))) then
    nseq = nseq + 1
    iseqtype(nseq) = 2
    anglexend = anglexi
    angleyend = angleyi
    anglezend = anglezi
    print*,'angle x start = ',anglex,' end = ',anglexend
    print*,'angle y start = ',angley,' end = ',angleyend
    print*,'angle z start = ',anglez,' end = ',anglezend
 endif
 !--change of render limits
 if (ndim.gt.1 .and. irender.gt.0 .and. irender.le.numplot) then
    if (notequal(rendermin,renderminfixed) .or. notequal(rendermax,rendermaxfixed)) then
       nseq = nseq + 1
       iseqtype(nseq) = 3
       icolchange = irender
       xmincolend = rendermin
       xmaxcolend = rendermax
       print*,trim(transform_label('rendermin,max',itrans(irender)))//' start = ',renderminfixed,rendermaxfixed, &
              ' end = ',xmincolend,xmaxcolend
       !--always store untransformed limits
       call transform_limits_inverse(xmincolend,xmaxcolend,itrans(irender))
    endif
 endif
 !--change of observer position
 if (ndim.eq.3 .and. notequal(zobserveri,zobserver)) then
    nseq = nseq + 1
    iseqtype(nseq) = 4
    zobserverend = zobserveri
 endif
 !--change of cross section position
 if (ndim.eq.3 .and. x_sec .and. notequal(xsecpos_nomulti,xsecposi)) then
    nseq = nseq + 1
    iseqtype(nseq) = 5
    xsecpos_nomulti_end = xsecposi
 endif
 !--change of opacity
 if (use3Dopacity .and. notequal(taupartdepthi,taupartdepth)) then
    nseq = nseq + 1
    iseqtype(nseq) = 6
    taupartdepthend = taupartdepthi
 endif

 !--all sequences start from 1 and end at current dump position
 iseqstart(1:nseq) = 1
 iseqend(1:nseq) = ipos

 if (nseq.gt.0) then
    print "(1x,a,i1,a)",'total of ',nseq,' sequences set:'
    do i=1,nseq
       print "(1x,i1,': ',a)",i,trim(labelseqtype(iseqtype(i)))
    enddo
    print "(a,i5)",' sequences start at dump 1 and end at dump ',ipos
    if (nframes.le.0) then
       if (ipos.eq.1) then
          nframes = 10
       else
          nframes = 1
       endif
       print "(a,i3)",' setting number of frames = ',nframes
    endif
 else
    print "(a)",' no sequences set (no change in parameters)'
 endif

 return
end subroutine setsequenceend

!----------------------------------------------------------------------
!  utility function for comparing real numbers
!----------------------------------------------------------------------
logical function notequal(r1,r2)
 implicit none
 real, intent(in) :: r1,r2

 if (abs(r1-r2).gt.epsilon(r1)) then
    notequal = .true.
 else
    notequal = .false.
 endif

end function notequal

!----------------------------------------------------------------------
! query function determining whether or not a given timestep
! is inside an animation sequence or not
! (and thus whether or not to generate extra frames)
!----------------------------------------------------------------------
logical function insidesequence(ipos)
 implicit none
 integer, intent(in) :: ipos
 integer :: i

 insidesequence = .false.
 do i=1,nseq
    if (iseqtype(i).gt.0 .and. iseqstart(i).le.ipos .and. iseqend(i).ge.ipos) then
       insidesequence = .true.
    endif
 enddo

 return
end function insidesequence

!----------------------------------------------------------------------
! query function which returns the current plot parameters
! based on the position in each sequence
! (given the current frame & dump position)
!----------------------------------------------------------------------
subroutine getsequencepos(ipos,iframe,iplotx,iploty,irender, &
                          anglexi,angleyi,anglezi,zobserveri,dzscreen,taupartdepthi, &
                          xsecposi,xmin,xmax,ymin,ymax,rendermin,rendermax,isetrenderlimits)
 use limits,     only:lim
 use multiplot,  only:itrans
 use transforms, only:transform_limits
 implicit none
 integer, intent(in)  :: ipos,iframe,iplotx,iploty,irender
 real, intent(out)    :: anglexi,angleyi,anglezi,zobserveri,dzscreen,taupartdepthi,xsecposi
 real, intent(out)    :: xmin,xmax,ymin,ymax,rendermin,rendermax
 logical, intent(out) :: isetrenderlimits
 logical :: logtaudepth
 integer :: i,iposinseq,iposend
 real    :: xfrac,xminstart,xmaxstart,xminend,xmaxend,yminstart,ymaxstart,yminend,ymaxend

 isetrenderlimits = .false.

 do i=1,nseq
    !--set starting values based on first position
    if (ipos.ge.iseqstart(i)) then
       iposinseq = (ipos-iseqstart(i))*nframes + iframe
       iposend = (iseqend(i)-iseqstart(i))*nframes + nframes
       xfrac = (iposinseq-1)/real(iposend-1)
       xfrac = min(xfrac,1.0)

       if (iposinseq.gt.iposend) then
          print "(1x,a)",'-->  '//trim(labelseqtype(iseqtype(i)))//' finished : frac = 1.0'
       else
          print "(1x,a,i3,a,i3,a,f5.2)",'-->  frame ', &
                 iposinseq,' / ',iposend,' of '//trim(labelseqtype(iseqtype(i)))//': frac = ',xfrac
       endif
       select case(iseqtype(i))
       case(1)
          xminstart = lim(iplotx,1)
          xmaxstart = lim(iplotx,2)
          yminstart = lim(iploty,1)
          ymaxstart = lim(iploty,2)
          call transform_limits(xminstart,xmaxstart,itrans(iplotx))
          call transform_limits(yminstart,ymaxstart,itrans(iploty))
          xminend = xminseqend
          xmaxend = xmaxseqend
          yminend = yminseqend
          ymaxend = ymaxseqend
          call transform_limits(xminend,xmaxend,itrans(iplotx))
          call transform_limits(yminend,ymaxend,itrans(iploty))
          !--steps are linear in the transformed space
          !  and limits returned are *already transformed*
          xmin = xminstart + xfrac*(xminend - xminstart)
          xmax = xmaxstart + xfrac*(xmaxend - xmaxstart)
          ymin = yminstart + xfrac*(yminend - yminstart)
          ymax = ymaxstart + xfrac*(ymaxend - ymaxstart)
       case(2)
          anglexi = anglex + xfrac*(anglexend - anglex)
          angleyi = angley + xfrac*(angleyend - angley)
          anglezi = anglez + xfrac*(anglezend - anglez)
       case(3)
          !--steps are linear in the transformed space
          !  and limits returned are *already transformed*
          if (iplotx.eq.icolchange) then
             xminstart = lim(iplotx,1)
             xmaxstart = lim(iplotx,2)
             call transform_limits(xminstart,xmaxstart,itrans(iplotx))
             xminend = xmincolend
             xmaxend = xmaxcolend
             call transform_limits(xminend,xmaxend,itrans(iplotx))
             xmin = xminstart + xfrac*(xminend - xminstart)
             xmax = xmaxstart + xfrac*(xmaxend - xmaxstart)
          elseif (iploty.eq.icolchange) then
             yminstart = lim(iploty,1)
             ymaxstart = lim(iploty,2)
             call transform_limits(yminstart,ymaxstart,itrans(iploty))
             yminend = xmincolend
             ymaxend = xmaxcolend
             call transform_limits(yminend,ymaxend,itrans(iploty))
             ymin = yminstart + xfrac*(yminend - yminstart)
             ymax = ymaxstart + xfrac*(ymaxend - ymaxstart)
          elseif (irender.eq.icolchange) then
             xminstart = lim(irender,1)
             xmaxstart = lim(irender,2)
             call transform_limits(xminstart,xmaxstart,itrans(irender))
             xminend = xmincolend
             xmaxend = xmaxcolend
             call transform_limits(xminend,xmaxend,itrans(irender))
             rendermin = xminstart + xfrac*(xminend - xminstart)
             rendermax = xmaxstart + xfrac*(xmaxend - xmaxstart)
             isetrenderlimits = .true.
          endif
       case(4)
          zobserveri = zobserver + xfrac*(zobserverend - zobserver)
          dzscreen = zobserveri
       case(5)
          xsecposi = xsecpos_nomulti + xfrac*(xsecpos_nomulti_end - xsecpos_nomulti)
       case(6)
          logtaudepth = (taupartdepthend .gt. 1.001e3*taupartdepth) &
                    .or.(taupartdepthend .lt. 1.001e-3*taupartdepth)
          if (logtaudepth) then
             print "(a)",'     (incrementing optical depth logarithmically)'
             taupartdepthi = taupartdepth*(taupartdepthend/taupartdepth)**xfrac
          else
             taupartdepthi = taupartdepth + xfrac*(taupartdepthend - taupartdepth)
          endif
       end select
    endif
 enddo

 return
end subroutine getsequencepos

end module settings_xsecrot
