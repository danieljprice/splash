!-----------------------------------------------------------------------------
! options for rendering / vector plots
!-----------------------------------------------------------------------------

subroutine options_render
  use prompting
  use settings
  implicit none
  integer :: ians
!
!--rendering options
!
 print 10,npix,icolours,iplotcont_nomulti,ncontours, &
      iPlotColourBar
10  format(' 0) exit ',/, 		     &
           ' 1) change number of pixels          (',i5,' )',/, &
           ' 2) change colour scheme             (',i2,' )',/,    &
           ' 3) toggle plot contours             ( ',L1,' )',/, &
           ' 4) change number of contours        (',i3,' )',/, &
           ' 5) toggle colour bar                ( ',L1,' )')
 call prompt('enter option',ians,0,5)
!
!--options
!
 select case(ians)
!------------------------------------------------------------------------
    case(1)
       call prompt('enter number of pixels along x axis',npix,1,10000)
!------------------------------------------------------------------------
    case(2)
100    continue
       print *,'(-ve = demo, 0 = contours only)'
       call prompt('enter colour scheme for rendering ',icolours,max=5)
       if (icolours.lt.0) then
          call colour_demo
          goto 100
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
  end select
    
 return
end subroutine options_render
