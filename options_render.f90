!-----------------------------------------------------------------------------
! options for rendering / vector plots
!-----------------------------------------------------------------------------

subroutine options_render
  use colours
  use prompting
  use settings_render
  implicit none
  integer :: ians,i
!
!--rendering options
!
  ians = 0
  print 10,npix,icolours,iplotcont_nomulti,ncontours, &
        iPlotColourBar
10 format(' 0) exit ',/,                      &
           ' 1) change number of pixels          (',i5,' )',/, &
           ' 2) change colour scheme             (',i2,' )',/,    &
           ' 3) toggle plot contours             ( ',L1,' )',/, &
           ' 4) change number of contours        (',i3,' )',/, &
           ' 5) toggle colour bar                ( ',L1,' )',/, &
           ' 6) toggle pixels or particles       ( ',L1,' )' )
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
       promptloop: do
          if (icolours.lt.0) icolours = 1
          write(*,"(6(5(i2,a),/))") (i,'='//trim(schemename(i)),i=1,ncolourschemes)
          print *,'(-ve = demo, 0 = contours only)'
          call prompt('enter colour scheme for rendering ',icolours,max=ncolourschemes)
          if (icolours.lt.0) then
             call colour_demo
             cycle promptloop
          else
             exit promptloop
          endif
       enddo promptloop
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
       iPlotColourBar = .not.iPlotColourBar
       print*,'plot colour bar = ',iPlotColourBar
!------------------------------------------------------------------------
    case(6)
       icolour_particles = .not.icolour_particles
       print*,'particles colouring = ',icolour_particles

  end select
    
 return
end subroutine options_render
