!-------------------------------------------------------------------------
! Module containing settings and options relating to renderings
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_render
 implicit none
 integer :: ncontours,npix,icolours
 logical :: iplotcont_nomulti
 logical :: iPlotColourBar,icolour_particles
 real :: ColourBarDisp
 !--colour bar width is set here as it must be known for page setup
 !  in principle it could be user-changeable but this adds pointless options
 real, parameter :: ColourBarWidth = 5.5

 namelist /renderopts/ npix,icolours,ncontours,iplotcont_nomulti, &
   iPlotColourBar,icolour_particles,ColourBarDisp

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_render
  implicit none

  icolours = 2               ! colour scheme to use
  npix = 100                 ! pixels in x direction for rendering
  iPlotColourBar = .true.! whether or not to plot the colour bar
  iplotcont_nomulti = .false. ! plot contours
  icolour_particles = .false. ! colour particles instead of using pixels
  ncontours = 30             ! number of contours to plot
  ColourBarDisp = 3.2

  return
end subroutine defaults_set_render

!-----------------------------------------------------------------------------
! options for rendering / vector plots
!-----------------------------------------------------------------------------
subroutine submenu_render
  use colours
  use prompting
  implicit none
  integer :: ians,i
!
!--rendering options
!
  ians = 0
  print 10,npix,icolours,iplotcont_nomulti,ncontours, &
        iPlotColourBar,icolour_particles
10 format(' 0) exit ',/,                      &
           ' 1) change number of pixels           (',i5,' )',/, &
           ' 2) change colour scheme              (',i2,' )',/,    &
           ' 3) toggle plot contours              ( ',L1,' )',/, &
           ' 4) change number of contours         (',i3,' )',/, &
           ' 5) colour bar options                ( ',L1,' )',/, &
           ' 6) use particle colours not pixels   ( ',L1,' )' )
  call prompt('enter option',ians,0,6)
!
!--options
!
 select case(ians)
!------------------------------------------------------------------------
    case(1)
       call prompt('enter number of pixels along x axis',npix,1,10000)
!------------------------------------------------------------------------
    case(2)
!       promptloop: do
!          if (icolours.lt.0) icolours = 1
       write(*,"(i2,a,1x)") (i,': '//trim(schemename(i)),i=1,ncolourschemes)
       print "(a)",'(-ve = inverse, 0 = contours only)'
       call prompt('enter colour scheme for rendering ',icolours,-ncolourschemes,ncolourschemes)
!          if (icolours.lt.0) then
!             call colour_demo
!             cycle promptloop
!          else
!             exit promptloop
!          endif
!       enddo promptloop
       !
       ! by default, plot contours if no colour scheme and don't if a colour scheme chosen
       !
       if (icolours.eq.0) then
          iplotcont_nomulti = .true.
       else
          iplotcont_nomulti = .false.
       endif
!------------------------------------------------------------------------
    case(3)
       iplotcont_nomulti = .not.iplotcont_nomulti
       print*,'plot contours = ',iplotcont_nomulti
!------------------------------------------------------------------------
    case(4)
       call prompt(' enter number of contours between min,max',ncontours,1,500)
!------------------------------------------------------------------------
    case(5)
       call prompt(' plot colour bar? ',iPlotColourBar)
       print*,'plot colour bar = ',iPlotColourBar
       if (iPlotColourBar) then
          call prompt(' enter displacement of text from edge (character heights) ', &
                      ColourBarDisp)
       endif
!------------------------------------------------------------------------
    case(6)
       icolour_particles = .not.icolour_particles
       print*,'particles colouring = ',icolour_particles

  end select
    
 return
end subroutine submenu_render

end module settings_render
