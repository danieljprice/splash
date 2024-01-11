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
!  Copyright (C) 2005-2022 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! Module containing settings and options relating to renderings
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_render
 use colourbar, only:ColourBarDisp,iplotcolourbarlabel,ColourBarPosx,ColourBarPosy,&
                     ColourBarLen,ColourBarFmtStr,ColourBarWidth
 use labels,    only:lenlabel
 use kernels,   only:ikernel
 implicit none
 integer :: ncontours,npix,icolours,iColourBarStyle,iColourBarPos
 logical :: iplotcont_nomulti,ilabelcont
 logical :: icolour_particles,inormalise_interpolations
 logical :: ifastrender,idensityweightedinterpolation
 logical :: double_rendering,exact_rendering,iauto_densityweighted
 character(len=lenlabel+20) :: projlabelformat
 integer :: iapplyprojformat
 character(len=120) :: rgbfile

 namelist /renderopts/ npix,icolours,ncontours,iplotcont_nomulti, &
   icolour_particles,ColourBarDisp,inormalise_interpolations, &
   ifastrender,idensityweightedinterpolation,iColourBarStyle, &
   iplotcolourbarlabel,ilabelcont,projlabelformat,iapplyprojformat, &
   double_rendering,ikernel,ColourBarPosx,ColourBarPosy,ColourBarLen,&
   ColourBarFmtStr,ColourBarWidth,iColourBarPos,rgbfile,exact_rendering,&
   iauto_densityweighted

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_render

 icolours = 2               ! colour scheme to use
 npix = 0                 ! pixels in x direction for rendering
 iColourBarStyle = 1        ! whether or not to plot the colour bar and style
 iplotcont_nomulti = .false. ! plot contours
 icolour_particles = .false. ! colour particles instead of using pixels
 ncontours = 30             ! number of contours to plot
 ColourBarDisp = 5.0 ! displacement of colour bar label from edge
 inormalise_interpolations = .false.       ! do not normalise interpolations
 ifastrender = .false. ! use accelerated rendering
 idensityweightedinterpolation = .false.
 iauto_densityweighted = .true.
 iplotcolourbarlabel = .true.
 ilabelcont = .false.   ! print numeric labels on contours
 projlabelformat = ' '
 iapplyprojformat = 0
 double_rendering = .false.
 ikernel = 0 ! just take default kernel
 ColourBarPosx = 0.75 ! default values used for floating colour bars
 ColourBarPosy = 0.7
 ColourBarLen  = 0.25
 ColourBarWidth = 2.
 ColourBarFmtStr = 'BCMSTV'
 iColourBarPos = 3
 rgbfile = 'colours.ctb'
 exact_rendering = .false.

 return
end subroutine defaults_set_render

!---------------------------------------------
! set default values for 360 video
!---------------------------------------------
subroutine defaults_set_render_360
 use colourbar, only:set_floating_bar_style

 iColourBarStyle = 7
 call set_floating_bar_style(iColourBarStyle,4)

end subroutine defaults_set_render_360

!-----------------------------------------------------------------------------
! options for rendered plots
!-----------------------------------------------------------------------------
subroutine submenu_render(ichoose)
 use colourbar, only:maxcolourbarstyles,labelcolourbarstyles,barisvertical,&
                      isfloating,iscustombar,labelfloatingstyles,&
                      set_floating_bar_style,maxfloatingstyles,ilogcolourbaraxis
 use colours,   only:schemename,ncolourschemes,ncoltable,rgbtable,icustom
 use prompting, only:prompt,print_logical
 use params,    only:maxplot
 use plotlib,   only:plotlib_supports_alpha
 use filenames, only:fileprefix
 use kernels,   only:select_kernel,kernelname,nkernels
 use projections3D, only:setup_integratedkernel
 use asciiutils,    only:read_asciifile
 integer, intent(in) :: ichoose
 character(len=5)    :: string,stringd
 character(len=20)   :: kname
 integer :: ians,i,ierr,icolourprev
!
!--rendering options
!
 ians = ichoose
 print "(a)",'----------------- rendering options -------------------'

 if (ians <= 0 .or. ians > 10) then
    if (npix > 0) then
       write(string,"(i5)") npix
    else
       string = 'AUTO'
    endif
    if (iauto_densityweighted) then
       stringd = 'AUTO'
    else
       stringd = print_logical(idensityweightedinterpolation)
    endif
    kname = ''
    if (ikernel >= 0 .and. ikernel <= nkernels) kname = trim(kernelname(ikernel))
    print 10,trim(string),icolours,print_logical(iplotcont_nomulti), &
           iColourBarStyle,print_logical(icolour_particles), &
           print_logical(inormalise_interpolations),print_logical(ifastrender),&
           trim(stringd),trim(projlabelformat),&
           trim(kname)
10  format( &
          ' 0) exit ',/,                      &
          ' 1) set number of pixels               ( ',a,' )',/, &
          ' 2) change colour scheme               (',i2,' )',/,    &
          ' 3) contour plot settings              ( ',a,' )',/, &
          ' 4) change colour bar style            ( ',i2,' )',/,&
          ' 5) use particle colours not pixels    ( ',a,' )',/,&
          ' 6) normalise interpolations           ( ',a,' )',/,&
          ' 7) use accelerated rendering          ( ',a,' )',/,&
          ' 8) use density weighted interpolation ( ',a,' )',/, &
          ' 9) customize label on projection plots ( ',a,' )',/,&
          ' 10) set kernel             ( ',a,' )')
    call prompt('enter option',ians,0,10)
 endif
!
!--options
!
 select case(ians)
!------------------------------------------------------------------------
 case(1)
    print "(5(/,a),/)",' Note: setting number of pixels = 0 means that ', &
                        '       the number of pixels will be automatically ', &
                        '       chosen to match the device used for plotting.', &
                        '       The number of pixels is then determined by ', &
                        '       the page size (set in the p)age menu).'
    call prompt('enter number of pixels along x axis (0=auto)',npix,0,100000)
!------------------------------------------------------------------------
 case(2)
    ierr = 1
    icolourprev = icolours
    write(*,"(i2,a,1x)") (i,': '//trim(schemename(i)),i=1,ncolourschemes)
    write(*,"(i3,a,1x)") icustom,': custom'
    print "(a)",'(-ve = inverse, 0 = contours only)'
    promptloop: do while (ierr /= 0)
       call prompt('enter colour scheme for rendering ',icolours,&
                      -ncolourschemes+1,ncolourschemes+1,icustom,icustom)
       !
       ! custom colour map from file !demonstration plot of all colour schemes
       !
       ierr = 0
       if (icolours==ncolourschemes+1)    icolours=icustom
       if (icolours==-(ncolourschemes+1)) icolours=-icustom
       if (abs(icolours)==icustom) then
          call prompt('enter filename for rgb table',rgbfile,noblank=.true.)
          call read_asciifile(rgbfile,ncoltable,rgbtable,ierr)
          if (ierr /= 0 .or. ncoltable <= 0) then
             print "(a)",'ERROR: could not read colours from '//trim(rgbfile)
          else
             print "(a,i3,a)",'read ',ncoltable,' colours from '//trim(rgbfile)
          endif
          if (ierr /= 0) icolours = icolourprev
       endif
    enddo promptloop
!------------------------------------------------------------------------
 case(3)
    if (icolours==0) then
       print "(2(/,a),/)",' Warning: this option has no effect if colour scheme 0 is set', &
                             '          (cannot plot contours on top of contours)'
    endif
    if (plotlib_supports_alpha) then
       call prompt(' allow contour/double render prompt?',iplotcont_nomulti)
       print "(3(/,a),/)",' Double rendering renders the first quantity in black and white', &
                             ' and the second in colour with a transparent background ', &
                             ' (such that data below the colour bar minimum appears transparent)'
       call prompt('use double rendering instead of contours?',double_rendering)
    else
       call prompt('allow contour plotting prompt?',iplotcont_nomulti)
    endif

    if (double_rendering) then
       print "(a)",' Second render prompt is '//trim(print_logical(iplotcont_nomulti))
    else
       print "(a)",' Contour plotting prompt is '//trim(print_logical(iplotcont_nomulti))
    endif
    if ((iplotcont_nomulti .or. icolours==0) .and. .not.double_rendering) then
       print "(5(/,a),/)",&
             ' To set contour levels and level labels manually, create a file called', &
             '  '''//trim(fileprefix)//'.contours'' in the working directory, with the following format:',&
             '  1.0   label1 ', &
             '  2.0   label2 ', &
             '  ...'
       call prompt('enter number of contours between min,max',ncontours,0,500)
       !call prompt('plot contour labels?',ilabelcont)
    endif
!------------------------------------------------------------------------
 case(4)
    do i=0,maxcolourbarstyles
       print "(i2,')',1x,a)",i,trim(labelcolourbarstyles(i))
    enddo
    call prompt('enter colour bar style to use ',iColourBarStyle,0,maxcolourbarstyles)
    print "(a,/)",'colour bar style = '//trim(labelcolourbarstyles(iColourBarStyle))

    if (iColourBarStyle > 0) then
       if (isfloating(iColourBarStyle)) then
          print "(5(a,/),a)",' Positioning of floating colour bar: ', &
                              (trim(labelfloatingstyles(i)),i=1,maxfloatingstyles)
          call prompt('enter option',iColourBarPos,1,maxfloatingstyles)
          if (iColourBarPos >= 1 .and. iColourBarPos < maxfloatingstyles) ColourBarLen = 0.25
          if (iColourBarPos == maxfloatingstyles) then
             call prompt('enter x position of colour bar as fraction of viewport',ColourBarPosx,-1.,1.5)
             call prompt('enter y position of colour bar as fraction of viewport',ColourBarPosy,-1.,1.5)
             call prompt('enter length of colour bar as fraction of viewport',ColourBarLen,0.,1.)
          else
             call set_floating_bar_style(iColourBarStyle,iColourBarPos)
          endif
       endif
       call prompt('plot colour bar label?',iplotcolourbarlabel)
       if (barisvertical(iColourBarStyle) .and. iplotcolourbarlabel) then
          call prompt('enter displacement of text from edge (character heights) ', &
                         ColourBarDisp)
       endif
       call prompt('use log axis on colour bar?',ilogcolourbaraxis)
       if (iscustombar(iColourBarStyle)) then
          call prompt('enter width of colour bar in character heights',ColourBarWidth,0.,20.)
          if (barisvertical(iColourBarstyle)) then
             print "(a)",' A=axis,B=bottom,C=top,T=major ticks,S=minor ticks,N=labels,V=vertical,L=log,M=labels above'
          else
             print "(a)",' B=left,C=right,T=major ticks,S=minor ticks,N=labels,L=log,M=labels to left'
             if (ColourBarFmtStr=='BCMSTV') then
                ColourBarFmtStr='BCNST' ! use default string for horizontal bars instead
             endif
          endif
          call prompt('enter format code string for colour bar ticks/numbering',ColourBarFmtStr)
       endif
    endif
!------------------------------------------------------------------------
 case(5)
    icolour_particles = .not.icolour_particles
    print "(a)",' Use particle colours instead of pixels is ' &
                   //trim(print_logical(icolour_particles))
!------------------------------------------------------------------------
 case(6)
    call prompt('Use normalised interpolations?',inormalise_interpolations)
    print "(a)",' Normalisation of interpolations is ' &
                   //trim(print_logical(inormalise_interpolations))
    if (.not.inormalise_interpolations) then
       call prompt('Use exact rendering? (slower but more accurate)',exact_rendering)
       print "(a)",' Exact rendering is ' &
                   //trim(print_logical(exact_rendering))
    endif
!------------------------------------------------------------------------
 case(7)
    ifastrender = .not.ifastrender
    print "(a)",' Accelerated rendering is '//trim(print_logical(ifastrender))
    if (ifastrender) then
       print*,' Warning: this is slightly approximate (particle position'
       print*,'          assumed to be at centre of pixel)'
    endif
!------------------------------------------------------------------------
 case(8)
    idensityweightedinterpolation = .not.idensityweightedinterpolation
    print "(a)",' Density weighted interpolation is '// &
                   print_logical(idensityweightedinterpolation)
    iauto_densityweighted = .not.idensityweightedinterpolation

    if (.not.idensityweightedinterpolation) &
       call prompt('Switch density weighting on/off automatically?',iauto_densityweighted)
!------------------------------------------------------------------------
 case(9)
    print "(5(/,a),/,4(/,a),/)", &
                        ' Example format strings: ', &
                        '  \(2268) %l d%z %uz       : this is the default format "\int rho [g/cm^3] dz [cm]"', &
                        '   column %l               : would print "column density" for density', &
                        '  surface %l               : would print "surface density"', &
                        '  %l integrated through %z : would print "density integrated through z"', &
                        ' Format codes: ', &
                        ' %l  : label for rendered quantity ', &
                        ' %z  : label for ''z'' ', &
                        ' %uz : units label for z (only if physical units applied)'

    call prompt(' enter label format for projection plots: ',projlabelformat)
    call prompt(' enter which column to apply format to (0=all) ',iapplyprojformat,0,maxplot)
!------------------------------------------------------------------------
 case(10)
    do i=0,nkernels
       print "(1x,i1,')',1x,a)",i,trim(kernelname(i))
    enddo
    call prompt(' enter kernel to use for interpolations (0=default)',ikernel,0,nkernels)
    call select_kernel(ikernel)
    call setup_integratedkernel  ! need to redo the kernel table if kernel has changed
 end select

 return
end subroutine submenu_render

end module settings_render
