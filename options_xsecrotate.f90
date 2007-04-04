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
 logical, public :: writeppm
 real, public :: anglex, angley, anglez, zobserver, dzscreenfromobserver
 real, public :: taupartdepth
 real, public :: xsecpos_nomulti,xseclineX1,xseclineX2,xseclineY1,xseclineY2
 real, public, dimension(3) :: xorigin,xminrotaxes,xmaxrotaxes
 
 !--private variables related to animation sequences
 integer, parameter, private :: maxseq = 6
 integer, dimension(maxseq), private :: iseqstart,iseqend,iseqtype
 integer, private :: icolchange
 real, private :: xminseqend,xmaxseqend,yminseqend,ymaxseqend
 real, private :: anglezend,angleyend,anglexend,zobserverend,taupartdepthend
 real, private :: xmincolend,xmaxcolend,xsecpos_nomulti_end
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
          taupartdepth,writeppm

 private :: animopts
 namelist /animopts/ nseq,nframes,iseqstart,iseqend,iseqtype, &
          xminseqend,xmaxseqend,yminseqend,ymaxseqend, &
          anglezend,angleyend,anglexend,zobserverend,taupartdepthend, &
          icolchange,xmincolend,xmaxcolend,xsecpos_nomulti_end

 !--public procedure names
 public :: defaults_set_xsecrotate,submenu_xsecrotate,getsequencepos,insidesequence
 public :: write_animfile,read_animfile
 
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
  taupartdepth = 2.
  writeppm = .true.

  !--defaults for animation sequences
  nseq = 0
  nframes = 0
  iseqstart(:) = 0
  iseqend(:) = 0
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

  return
end subroutine defaults_set_xsecrotate

!----------------------------------------------------------------------
! sets options relating to cross sectioning / rotation
!----------------------------------------------------------------------
subroutine submenu_xsecrotate(ichoose)
 use filenames, only:nsteps
 use labels, only:label,ix
 use limits, only:lim
 use prompting
 use settings_data, only:ndim
 implicit none
 integer, intent(in) :: ichoose
 integer :: ians,i
 logical :: iyes
 character(len=4) :: text
 
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
              ' 2) rotation on/off/options                      ( ',a,3(1x,f5.2),' )',/, &
              ' 3) 3D perspective on/off                        ( ',a,' )',/, &
              ' 4) 3D surface rendering on/off/options          ( ',a,' )',/, &
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
 case(3)
    use3Dperspective = .not.use3Dperspective
    call prompt(' Use 3D perspective? ',use3Dperspective)
    if (.not.use3Dperspective) use3Dopacityrendering = .false.
!------------------------------------------------------------------------
 case(4)
    use3Dopacityrendering = .not.use3Dopacityrendering
    call prompt(' Use 3D opacity rendering? ',use3Dopacityrendering)
    if (use3Dopacityrendering .and..not.use3Dperspective) then
       print "(a)",' also turning on 3D perspective (which must be set for this to work)'
       use3Dperspective = .true.
    endif
    if (use3Dopacityrendering) then
       print "(/,a)",' Warning: 3D opacity rendering sends only an approximate version '
       print "(a,/)",' to the PGPLOT device (not corrected for brightness) '
       call prompt(' Do you want to write a ppm file in addition to PGPLOT output?',writeppm)
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
    print "(a,i1,a)",'Note: Up to ',maxseq,' sequences (1 of each type) can be set '
    call prompt('Enter number of sequences to use (0=none)',nseq,0,maxseq)

    if (nseq.gt.0) then
       !--set sensible default value for number of frames
       if (nframes.eq.0) then
          if (nsteps.gt.1) then
             nframes = 1
          else
             nframes = 10
          endif
       endif
       call prompt('Enter number of frames generated between dumps (applies to all sequences)',nframes,1,500)

       iyes = .true.
       if (ihavesetsequence) call prompt('Change sequence settings?',iyes)
       if (iyes) call submenu_animation()
    endif    

 end select

 return
end subroutine submenu_xsecrotate

!----------------------------------------------------------------------
! sets up animation sequences
!----------------------------------------------------------------------
subroutine submenu_animation()
 use prompting, only:prompt
 use limits, only:lim
 use labels, only:ix
 use settings_data, only:ndim,istartatstep,iendatstep,numplot
 use filenames, only:nsteps
 implicit none
 integer :: i,j,ierr

 do i = 1, nseq
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
       if (i.gt.1) then
          if (any(iseqtype(1:i-1).eq.iseqtype(i)).and.iseqtype(i).gt.0) then
             print "(a)",' Error: can only have one sequence of each type '
             ierr = 2
          endif
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
       print "(a)",'Note: opacity sequence starts from current opacity value '
       call prompt('Enter finishing opacity in units of average smoothing length ',taupartdepthend)
    end select
 enddo

 if (all(iseqtype(1:nseq).eq.0)) then
    print "(a)",' No sequences set!'
    nseq = 0
 else
    ihavesetsequence = .true.
    print "(3(/,a))",'Note: these sequences are saved to file using the S)ave option from the main menu', &
                   ' This will save the splash.defaults file, the splash.limits file', &
                   ' and a file called ''splash.anim'' which contains the animation sequence info'
    print "(/,a)",' press any key to continue...'
    read*
 endif
 
 return
end subroutine submenu_animation

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
                          anglexi,angleyi,anglezi,zobserveri,taupartdepthi, &
                          xsecposi,xmin,xmax,ymin,ymax,rendermin,rendermax)
 use limits, only:lim
 implicit none
 integer, intent(in) :: ipos,iframe,iplotx,iploty,irender
 real, intent(out) :: anglexi,angleyi,anglezi,zobserveri,taupartdepthi,xsecposi
 real, intent(out) :: xmin,xmax,ymin,ymax,rendermin,rendermax
 integer :: i,iposinseq,iposend
 real :: xfrac
 
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
          print "(1x,a,i3,a,i3,a,f4.2)",'-->  frame ', &
                 iposinseq,' / ',iposend,' of '//trim(labelseqtype(iseqtype(i)))//': frac = ',xfrac
       endif
       select case(iseqtype(i))
       case(1)
          xmin = lim(iplotx,1) + xfrac*(xminseqend - lim(iplotx,1))
          xmax = lim(iplotx,2) + xfrac*(xmaxseqend - lim(iplotx,2))
          ymin = lim(iploty,1) + xfrac*(yminseqend - lim(iploty,1))
          ymax = lim(iploty,2) + xfrac*(ymaxseqend - lim(iploty,2))
       case(2)
          anglexi = anglex + xfrac*(anglexend - anglex)
          angleyi = angley + xfrac*(angleyend - angley)
          anglezi = anglez + xfrac*(anglezend - anglez)
       case(3)
          if (iplotx.eq.icolchange) then
             xmin = lim(iplotx,1) + xfrac*(xmincolend - lim(iplotx,1))
             xmax = lim(iplotx,2) + xfrac*(xmaxcolend - lim(iplotx,2))
          elseif (iploty.eq.icolchange) then
             ymin = lim(iploty,1) + xfrac*(xmincolend - lim(iploty,1))
             ymax = lim(iploty,2) + xfrac*(xmaxcolend - lim(iploty,2))
          elseif (irender.eq.icolchange) then
             rendermin = lim(irender,1) + xfrac*(xmincolend - lim(irender,1))
             rendermax = lim(irender,2) + xfrac*(xmaxcolend - lim(irender,2))
          endif
       case(4)
          zobserveri = zobserver + xfrac*(zobserverend - zobserver)
       case(5)
          xsecposi = xsecpos_nomulti + xfrac*(xsecpos_nomulti_end - xsecpos_nomulti)
       case(6)
          taupartdepthi = taupartdepth + xfrac*(taupartdepthend - taupartdepth)
       end select
    endif
 enddo
 
 return
end subroutine getsequencepos

!-----------------------------------------------
! writes animation sequence options to file 
! (should match read_animfile)
!-----------------------------------------------
subroutine write_animfile(filename)
 use prompting, only:prompt
 implicit none
 character(len=*), intent(in) :: filename
 integer :: ierr
 logical :: iexist, idelete
 
 if (nseq.gt.0) then
    open(unit=15,file=filename,status='replace',form='formatted', &
         delim='apostrophe',iostat=ierr) ! without delim namelists may not be readable
       if (ierr /= 0) then 
          print*,'ERROR: cannot write file '//trim(filename)
          close(unit=15)
          return
       endif
       write(15,NML=animopts)
    close(unit=15)
    print*,'animation sequences saved to file '//trim(filename)
 elseif (nseq.eq.0) then
    inquire(file=trim(filename),exist=iexist)
    if (iexist) then
       idelete = .true.
       call prompt(' delete '//trim(filename)//' file? ',idelete)
       if (idelete) then
          open(unit=15,status='replace',file=filename,iostat=ierr)
          close(unit=15,status='delete',iostat=ierr)
          if (ierr /= 0) then
             print "(a)",' Error deleting '//trim(filename)
          else
             print "(a)",trim(filename)//' deleted'
          endif
       endif
    endif
 endif
    
 return              
end subroutine write_animfile

!-----------------------------------------------
! reads animation sequence options from file 
! (should match write_animfile)
!-----------------------------------------------
subroutine read_animfile(filename)
 use filenames, only:nsteps
 implicit none
 character(len=*), intent(in) :: filename
 logical :: iexist
 integer :: ierr
 
 inquire (exist=iexist, file=filename)
 if (iexist) then
    open(unit=15,file=filename,status='old',form='formatted')

    ierr = 0
    read(15,NML=animopts,end=77,iostat=ierr)
    if (ierr /= 0) print "(a)",'error reading animation sequences from '//trim(filename)

    close(unit=15)
    print "(1x,a)",'read animation sequences from '//trim(filename)
    ihavesetsequence = .true.
    return
 else
    return
 endif
 
77 continue
 print*,'**** warning: end of file in '//trim(filename)//' ****'
 close(unit=15)
 
 if (nseq.gt.0 .and. all(iseqstart(1:nseq).gt.nsteps)) then
    print "(a)",' WARNING: animation sequences have no effect!! (not enough dumps)'
 endif

 return
end subroutine read_animfile

end module settings_xsecrot
