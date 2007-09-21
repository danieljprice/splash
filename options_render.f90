!-------------------------------------------------------------------------
! Module containing settings and options relating to renderings
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_render
 implicit none
 integer :: ncontours,npix,icolours
 logical :: iplotcont_nomulti
 logical :: iPlotColourBar,icolour_particles,inormalise_interpolations
 logical :: ifastrender,idensityweightedinterpolation
 real :: ColourBarDisp
 !--colour bar width is set here as it must be known for page setup
 !  in principle it could be user-changeable but this adds pointless options
 real, parameter :: ColourBarWidth = 5.5

 namelist /renderopts/ npix,icolours,ncontours,iplotcont_nomulti, &
   iPlotColourBar,icolour_particles,ColourBarDisp,inormalise_interpolations, &
   ifastrender,idensityweightedinterpolation

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_render
  implicit none

  icolours = 2               ! colour scheme to use
  npix = 200                 ! pixels in x direction for rendering
  iPlotColourBar = .true.! whether or not to plot the colour bar
  iplotcont_nomulti = .false. ! plot contours
  icolour_particles = .false. ! colour particles instead of using pixels
  ncontours = 30             ! number of contours to plot
  ColourBarDisp = 5.0 ! displacement of colour bar label from edge
  inormalise_interpolations = .false.       ! do not normalise interpolations
  ifastrender = .false. ! use accelerated rendering
  idensityweightedinterpolation = .false.

  return
end subroutine defaults_set_render

!-----------------------------------------------------------------------------
! options for rendered plots
!-----------------------------------------------------------------------------
subroutine submenu_render(ichoose)
  use colours, only:schemename,ncolourschemes,colour_demo
  use prompting, only:prompt,print_logical
  implicit none
  integer, intent(in) :: ichoose
  integer :: ians,i,ierr,icolourprev
!
!--rendering options
!
  ians = ichoose
  print "(a)",'----------------- rendering options -------------------'
  
  if (ians.le.0 .or. ians.gt.8) then
     print 10,npix,icolours,print_logical(iplotcont_nomulti),ncontours, &
           print_logical(iPlotColourBar),print_logical(icolour_particles), &
           print_logical(inormalise_interpolations),print_logical(ifastrender),&
           print_logical(idensityweightedinterpolation)
10   format( &
          ' 0) exit ',/,                      &
          ' 1) change number of pixels            (',i5,' )',/, &
          ' 2) change colour scheme               (',i2,' )',/,    &
          ' 3) plot contours                      ( ',a,' )',/, &
          ' 4) change number of contours          (',i3,' )',/, &
          ' 5) colour bar options                 ( ',a,' )',/,&
          ' 6) use particle colours not pixels    ( ',a,' )',/,& 
          ' 7) normalise interpolations           ( ',a,' )',/,&
          ' 8) use accelerated rendering          ( ',a,' )',/,&
          ' 9) use density weighted interpolation ( ',a,' )')
     call prompt('enter option',ians,0,9)
  endif
!
!--options
!
 select case(ians)
!------------------------------------------------------------------------
    case(1)
       call prompt('enter number of pixels along x axis',npix,1,10000)
!------------------------------------------------------------------------
    case(2)
       ierr = 1
       icolourprev = icolours
       write(*,"(i2,a,1x)") (i,': '//trim(schemename(i)),i=1,ncolourschemes)
       write(*,"(i2,a,1x)") ncolourschemes+1,': demo plot of all schemes'
       print "(a)",'(-ve = inverse, 0 = contours only)'
       promptloop: do while (ierr /= 0)
          call prompt('enter colour scheme for rendering ',icolours,-ncolourschemes,ncolourschemes+1)
          !
          ! demonstration plot of all colour schemes
          !      
          ierr = 0
          if (abs(icolours).eq.ncolourschemes+1) then
             call colour_demo
             icolours = icolourprev
             ierr = 1
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
       call prompt(' enter number of contours between min,max',ncontours,0,500)
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
!------------------------------------------------------------------------
    case(7)
       inormalise_interpolations = .not.inormalise_interpolations
       print*,'normalisation of interpolations = ',inormalise_interpolations
!------------------------------------------------------------------------
    case(8)
       ifastrender = .not.ifastrender
       print*,'accelerated rendering = ',ifastrender
       if (ifastrender) then
          print*,' Warning: this is slightly approximate (particle position'
          print*,'          assumed to be at centre of pixel)'
       endif
!------------------------------------------------------------------------
    case(9)
       idensityweightedinterpolation = .not.idensityweightedinterpolation
       print "(a)",' density weighted interpolation is '// &
                   print_logical(idensityweightedinterpolation)
  end select
    
 return
end subroutine submenu_render

end module settings_render
