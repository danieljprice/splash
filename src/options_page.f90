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
!  Copyright (C) 2005-2015 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! Module containing page settings and options
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_page
 use settings_limits, only:iadapt,iadaptcoords,adjustlimitstodevice,xminoffset_track,xmaxoffset_track
 use labels,          only:lenlabel
 implicit none
 integer :: iaxis,nacross,ndown,ipapersize,nstepsperpage,linewidth,iscalepanel
 integer :: iPlotLegendOnlyOnPanel,modlinestyle,modcolour,maxlinestyle,maxcolour
 integer :: iPageColours,ipapersizeunits
 logical :: iColourEachStep,iChangeStyles,tile,interactive,nomenu
 logical :: iPlotLegend,iPlotStepLegend,iPlotTitles,usecolumnorder
 logical :: iPlotScale,iUseBackgroundColourForAxes,usesquarexy
 real    :: papersizex,aspectratio
 real    :: hposlegend,vposlegend,fjustlegend,hpostitle,vpostitle,fjusttitle
 real    :: charheight,alphalegend
 real    :: dxscale,hposscale,vposscale,yscalealt
 real    :: xminpagemargin,xmaxpagemargin,yminpagemargin,ymaxpagemargin
 character(len=lenlabel) :: legendtext, scaletext
 character(len=60)       :: device
 character(len=lenlabel) :: labelyalt

 namelist /pageopts/ iaxis,nacross,ndown,interactive,iadapt,iadaptcoords, &
   nstepsperpage,iColourEachStep,iChangeStyles,tile,ipapersize,papersizex,aspectratio, &
   iPlotLegend,iPlotStepLegend,hposlegend,vposlegend,iPlotTitles,hpostitle, &
   vpostitle,fjusttitle,legendtext,iPageColours,charheight,linewidth,&
   fjustlegend,iPlotLegendOnlyOnPanel, &
   iPlotScale,dxscale,scaletext,hposscale,vposscale,iscalepanel,iUseBackgroundColourForAxes, &
   usesquarexy,maxlinestyle,modlinestyle,maxcolour,modcolour,usecolumnorder,ipapersizeunits,&
   adjustlimitstodevice,alphalegend,yscalealt,labelyalt,xminoffset_track,xmaxoffset_track, &
   xminpagemargin,xmaxpagemargin,yminpagemargin,ymaxpagemargin

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_page
  use shapes,  only:defaults_set_shapes
  use plotlib, only:plotlib_maxlinecolour
  implicit none

  interactive = .true.     ! default for interactive mode
  iaxis = 0                ! turns axes off/on
  nstepsperpage = 1
  iColourEachStep = .true. ! change colours if nstepsperpage > 1
  iChangeStyles = .false.  ! change marker/ line styles if nstepsperpage > 1
  tile = .true.
  usecolumnorder = .true.
  nacross = 1           ! number of plots across page
  ndown = 1             ! number of plots down page
  ipapersize = 0        ! paper size option
  papersizex = 0.0      ! size of x paper (no call to PGPAP if zero)
  aspectratio = 0.0     ! aspect ratio of paper (no call to PGPAP if zero)
  ipapersizeunits = 1   ! units in which the paper size is set

  iPlotLegend = .true.  ! whether or not to plot legend
  iPlotStepLegend = .false. ! timestep legend
  hposlegend = 0.95     ! horizontal legend position as fraction of viewport
  vposlegend = 2.0      ! vertical legend position in character heights
  fjustlegend = 1.0    ! justification factor for legend
  alphalegend = 0.5    ! transparency of overlaid annotation
  legendtext = 't='
  iPlotLegendOnlyOnPanel = 0

  iPlotTitles = .false.  ! whether or not to plot titles
  hpostitle = 0.5       ! horizontal title position as fraction of viewport
  vpostitle = 1.0       ! vertical title position in character heights
  fjusttitle = 0.5      ! justification factor for title
  iPageColours = 0
  charheight = 1.0    ! character height
  linewidth = 0       ! line width

  iPlotScale = .false.
  hposscale = 0.5
  vposscale = 1.0
  dxscale = 1.0
  scaletext = '1 unit'
  iscalepanel = 0
  maxlinestyle = 5
  modlinestyle = 1
  modcolour = 1
  maxcolour = plotlib_maxlinecolour
  
  yscalealt = 1.
  labelyalt = ' '

  usesquarexy = .true. ! spatial dimensions have same scale
  call defaults_set_shapes

  xminpagemargin = 0.
  xmaxpagemargin = 0.
  yminpagemargin = 0.
  ymaxpagemargin = 0.

  return
end subroutine defaults_set_page

!---------------------------------------------
! changed default values for evsplash
!---------------------------------------------
subroutine defaults_set_page_ev
  implicit none

  nstepsperpage   = 1000
  iColourEachStep = .true.   ! change colours if nstepsperpage > 1
  iChangeStyles   = .true.   ! change marker/ line styles if nstepsperpage > 1
  iPlotLegend     = .true.   ! whether or not to plot legend
  iPlotStepLegend = .true.   ! timestep legend
  hposlegend      = 0.1      ! horizontal legend position as fraction of viewport
  vposlegend      = 2.0      ! vertical legend position in character heights
  fjustlegend     = 0.0      ! justification factor for legend

  return
end subroutine defaults_set_page_ev

!----------------------------------------------------------------------
! submenu with options relating to page setup
!----------------------------------------------------------------------
subroutine submenu_page(ichoose)
 use params,      only:maxplot
 use prompting,   only:prompt,print_logical
 use pagecolours, only:pagecolourscheme,colour_fore,colour_back,maxpagecolours
 use plotlib,     only:plotlib_supports_alpha,plotlib_maxlinecolour,plotlib_maxlinestyle,plotlib_is_pgplot
 implicit none
 integer, intent(in) :: ichoose
 integer             :: iaction,i,iunitsprev,ierr
 real                :: papersizey,mnraslength
 character(len=15)   :: paperfmtstr
 character(len=3)    :: string

 iaction = ichoose
 
 papersizey = papersizex*aspectratio
 if (ipapersizeunits.gt.0) then
    write(paperfmtstr,"(f5.2,1x,f5.2)") papersizex,papersizey
 else
    write(paperfmtstr,"(i5,' x ',i4)") nint(papersizex),nint(papersizey)
 endif

 print "(a)",'---------------- page setup options -------------------'

 if (iaction.le.0 .or. iaction.gt.9) then
    print "( "// &
          "' 0) exit ',/, "// &
          "' 1) plot n steps on top of each other   (n =',i4,')',/, "// &
          "' 2) axes options                        (',i2,')',/, "// &
          "' 3) change paper size                   ("//trim(paperfmtstr)//" )',/, "// &
          "' 4) subdivide page into panels          (',i2,1x,'x',1x,i2,', tiling is ',a,')',/, "// &
          "' 5) spatial dimensions have same scale  ( ',a,' )',/,"// &
          "' 6) set character height                (',f4.1,')',/,"// &
          "' 7) adjust line width                   (',i2, ')',/,"// &
          "' 8) adjust page margins                 ( ',f4.1,' )',/,"// &
          "' 9) set foreground/background colours   ( ',a,' )')", &
          nstepsperpage,iaxis,nacross,ndown,print_logical(tile), &
          trim(print_logical(usesquarexy)),charheight,linewidth,&
          xminpagemargin, &
          trim(pagecolourscheme(iPageColours,short=.true.))

    call prompt('enter option ',iaction,0,9)
 endif

 select case(iaction)
!------------------------------------------------------------------------
  case(1)
     call prompt('Enter number of timesteps per panel ',nstepsperpage,0)
     print*,'Plotting up to ',nstepsperpage,' timesteps per panel'
     if (nstepsperpage.gt.1) then
        if (iadapt .or. iadaptcoords) then
           print "(a)",'(note that adaptive plot limits are now off)'
           iadapt = .false.
           iadaptcoords = .false.
        endif
        if (nstepsperpage.gt.14) then
           print "(a)",'(warning: steps per panel > number of colours, ie. colours will repeat)'
        endif
        call prompt('Use different colours for each step?',iColourEachStep)
        if (iColourEachStep) then
           call prompt('How often to change colour? (1=every step, 2=every 2nd step etc.)',modcolour,1)
           call prompt('Enter max number of colours to use before repeating (16=plot lib max)',&
                       maxcolour,1,plotlib_maxlinecolour)
        endif
!!        if (.not.iColourEachStep) icolourthisstep = 1
        call prompt('Use different markers/line style for each step? ',iChangeStyles)
        if (iChangeStyles) then
           call prompt('How often to change line style (1=every step, 2=every 2nd step etc.)',modlinestyle,1)
           write(string,"(i3)",iostat=ierr) plotlib_maxlinestyle
           call prompt('Enter max number of line styles to cycle through before repeating ('// &
                       trim(adjustl(string))//'=plot lib max)',maxlinestyle,1,plotlib_maxlinestyle)
        endif

        if (iColourEachStep .or. iChangeStyles) then
           print "(/,a,/,a)",' (to change the legend text, create a file called', &
                       '  ''legend'' in the working directory, with one label per line)'
           call prompt('Plot legend of marker styles/colours?',iPlotStepLegend)
        endif
     endif
     return
!------------------------------------------------------------------------
  case(2)
     print*,'-4 : draw box and major tick marks only;'
     print*,'-3 : draw box and tick marks (major and minor) only;'
     print*,'-2 : draw no box, axes or labels;'
     print*,'-1 : draw box only;'
     print*,' 0 : draw box and label it with coordinates;'
     print*,' 1 : same as AXIS=0, but also draw the coordinate axes (X=0, Y=0);'
     print*,' 2 : same as AXIS=1, but also draw grid lines at major increments of the coordinates;'
     print*,' 3 : draw box, ticks and numbers but no axes labels;'
     print*,' 4 : same as AXIS=0, but with a second y-axis scaled and labelled differently'
     print*,'10 : draw box and label X-axis logarithmically;'
     print*,'20 : draw box and label Y-axis logarithmically;'
     print*,'30 : draw box and label both axes logarithmically.'
     call prompt('enter axis option ',iaxis,-4,30)
     if (iaxis.eq.4) then
        call prompt('enter scale factor for alternative y axis',yscalealt,0.)
        call prompt('enter label for alternative y axis',labelyalt)
     endif
     print *,' axis = ',iaxis
     return
!------------------------------------------------------------------------
  case(3)
     print*,' 0) plotting library default '
     print*,' 1) small square         :  2.92 x 2.92 inches'
     print*,' 2) medium square        :  5.85 x 5.85 inches'
     print*,' 3) large square         :  8.00 x 8.00 inches'
     print*,' 4) single small graph   :  5.85 x 4.13 inches'
     print*,' 5) duo small graph      : 11.70 x 4.13 inches'
     print*,' 6) duo graph            : 11.70 x 6.00 inches'
     if (plotlib_is_pgplot) then
        print*,'  7) Custom size '
        call prompt(' Enter option for paper size ',ipapersize,0,7)
     else
        print*,' 7) 800  x 600  pixels'
        print*,' 8) 640  x 360  pixels (360p)'
        print*,' 9) 1280 x 720  pixels (720p)'
        print*,'10) 1920 x 1080 pixels (1080p/Full HD)'
        print*,'11) 1024 x 768  pixels'
        print*,'12) 1440 x 900  pixels'
        print*,'13) 2560 x 1440 pixels'
        print*,'14) 2560 x 1600 pixels'
        print*,'15) 3840 x 2160 pixels (4KTV/Ultra HD)'
        print*,'16) 4096 x 2160 pixels (Cinema 4K)'
        print*,'17) 5120 x 2880 pixels (5K)'
        print*,'18) 27320 x 3072 pixels (CAVE-2)'
        print*,'19) 1/3 of A4 journal page  79mm x 180 mm'
        print*,'20) 1/2 of A4 journal page 118mm x 180 mm'
        print*,'21) full A4 journal page   236mm x 180 mm'
        print*,'22) Custom size '
        call prompt(' Enter option for paper size ',ipapersize,0,22)
     endif
     
     select case(ipapersize)
     case(1)
        ipapersizeunits = 1
        papersizex = 0.25*11.7
        aspectratio = 1.0
     case(2)
        ipapersizeunits = 1
        papersizex = 0.5*11.7
        aspectratio = 1.0
     case(3)
        ipapersizeunits = 1
        papersizex = 8.0
        aspectratio = 1.0
     case(4)
        ipapersizeunits = 1
        papersizex = 0.5*11.7
        aspectratio = 1./sqrt(2.)
     case(5)
        ipapersizeunits = 1
        papersizex = 11.7
        aspectratio = 0.5/sqrt(2.)
     case(6)
        ipapersizeunits = 1
        papersizex = 11.7
        papersizey = 6.0
        aspectratio = papersizey/papersizex
     case(7)
        if (plotlib_is_pgplot) then
           ipapersizeunits = 1
           call prompt(' x size (inches) ',papersizex,0.0)
           call prompt(' y size (inches) or aspect ratio (-ve)',papersizey)
           if (papersizey.lt.0.0) then
              aspectratio = abs(papersizey)
           else
              aspectratio = papersizey/papersizex
           endif
        else
           ipapersizeunits = 0
           papersizex = 800.
           papersizey = 600.
           aspectratio = papersizey/papersizex
        endif
     case(8:21)
        if (plotlib_is_pgplot) then
           ipapersizeunits = 1
           papersizex  = 0.  ! use PGPLOT default
           aspectratio = 0.
        else
           ipapersizeunits = 0
           mnraslength = 56.*0.42175 ! 56pc = 672pts converted to cm
           select case(ipapersize)
           case(8)
              papersizex = 640.
              papersizey = 360.
           case(9)
              papersizex = 1280.
              papersizey = 720.
           case(10)
              papersizex = 1920.
              papersizey = 1080.
           case(11)
              papersizex = 1024.
              papersizey = 768.
           case(12)
              papersizex = 1440.
              papersizey = 900.
           case(13)
              papersizex = 2560.
              papersizey = 1440.
           case(14)
              papersizex = 2560.
              papersizey = 1600.
           case(15)
              papersizex = 3840.
              papersizey = 2160.
           case(16)
              papersizex = 4096.
              papersizey = 2160.
           case(17)
              papersizex = 5120.
              papersizey = 2880.
           case(18)
              papersizex = 27320.
              papersizey = 3072.
           case(19)  ! 1/3 of MNRAS page
              ipapersizeunits = 2
              papersizex = 18.
              papersizey = mnraslength/3.
           case(20)  ! 1/2 of MNRAS page
              ipapersizeunits = 2
              papersizex = 18.
              papersizey = 0.5*mnraslength
           case(21)  ! full MNRAS page, allowing 2cm for caption
              ipapersizeunits = 2
              papersizex = 18.
              papersizey = mnraslength - 2.
           end select
           aspectratio = papersizey/papersizex
        endif
     case(22)
        if (plotlib_is_pgplot) then
           ipapersizeunits = 1
           papersizex  = 0.  ! use PGPLOT default
           aspectratio = 0.
        else
           iunitsprev = ipapersizeunits
           print*,' 0) pixels '
           print*,' 1) inches '
           print*,' 2) cm '
           call prompt(' choose units for paper size',ipapersizeunits,0,2)
           if (ipapersizeunits.ne.iunitsprev) then
              select case(ipapersizeunits)
              case(2)
                 papersizex = 29.7
                 papersizey = 21.0
              case(1)
                 papersizex = 11.0
                 papersizey = 8.5
              case(0)
                 papersizex = 800.
                 papersizey = -0.75
              case default
                 papersizex = 0.
                 papersizey = 0.
              end select
           endif
           call prompt(' x size in above units ',papersizex,1.)
           call prompt(' y size or aspect ratio (-ve)',papersizey)
           if (papersizey.lt.0.0) then
              aspectratio = abs(papersizey)
           else
              aspectratio = papersizey/papersizex
           endif
       endif 
     case default
        ipapersizeunits = 1
        papersizex = 0.0         ! no call to PGPAP if they are zero
        aspectratio = 0.0
     end select
     
     call prompt('Adjust plot limits to match device aspect ratio?',adjustlimitstodevice)

     return
!------------------------------------------------------------------------
  case(4)
     call prompt('Enter number of plots across (columns):',nacross,1,maxplot)
     call prompt('Enter number of plots down   (rows):',ndown,1,maxplot/nacross)
     if (nacross*ndown.gt.1) then
        call prompt('Tile plots on the page where possible?',tile)
        call prompt('Plot panels across-then-down? (no=down-then-across)',usecolumnorder)
     endif
     return
!------------------------------------------------------------------------
  case(5)
     usesquarexy = .not.usesquarexy
     print "(a)",' Same scale for spatial dimensions is '//print_logical(usesquarexy)
!------------------------------------------------------------------------
  case(6)
     call prompt('Enter character height ',charheight,0.1,10.)
     return
!------------------------------------------------------------------------
  case(7)
     print "(3(/,a))",' Setting line width to 0 means automatic line width choice:', &
                    ' This gives width = 2 for vector devices (/ps,/cps etc)', &
                    '        and width = 1 elsewhere (e.g. for pixel devices)'
     print*
     call prompt('Enter line width (0=auto)',linewidth,0)
     return
!------------------------------------------------------------------------
  case(8)
     call prompt('Enter xmin page margin as fraction of viewport ',xminpagemargin,-1.,1.)
     call prompt('Enter xmax page margin as fraction of viewport ',xmaxpagemargin,-1.,1.)
     call prompt('Enter ymin page margin as fraction of viewport ',yminpagemargin,-1.,1.)
     call prompt('Enter ymax page margin as fraction of viewport ',ymaxpagemargin,-1.,1.)
     !interactive = .not.interactive
     !print "(a)",' Interactive mode is '//print_logical(interactive)
!------------------------------------------------------------------------
  case(9)
     print "(3(/,i1,')',1x,a))",(i,pagecolourscheme(i),i=0,maxpagecolours)
     call prompt(' Choose page colour scheme ',iPageColours,0,maxpagecolours)

     write(*,"(3(/,a))",advance='no') &
       ' Overlaid (that is, drawn inside the plot borders) axis ',&
       ' ticks, legend text and titles are by default plotted ', &
       ' in the foreground colour'

     if (iPageColours.gt.0) then
        print "(a,/)",' [i.e. '//trim(colour_fore(iPageColours))//'].'
        call prompt('Do you want to plot these in background colour [i.e. '&
                    //trim(colour_back(iPageColours))//'] instead?',&
                    iUseBackgroundColourForAxes)
     else
        print "(a,/)",'.'
        call prompt('Use background colour for these? ',iUseBackgroundColourForAxes)
     endif
     if (iUseBackgroundColourForAxes .and. plotlib_supports_alpha) then
        call prompt('Enter opacity for overlaid text and annotation ',alphalegend,0.0,1.0)
     endif

     return
  end select

 return
end subroutine submenu_page

!----------------------------------------------------------------------
! submenu with options relating to legend and title settings
!----------------------------------------------------------------------
subroutine submenu_legend(ichoose)
 use filenames, only:fileprefix
 use prompting, only:prompt,print_logical
 use shapes,    only:nshapes,labelshapetype,shape,submenu_shapes
 use legends,   only:prompt_panelselect
 implicit none
 integer, intent(in) :: ichoose
 integer             :: iaction,i,ierr,i1,i2
 character(len=50)   :: string

 iaction = ichoose
 print "(a)",'---------------- legend and title options -------------------'

 if (iPlotStepLegend) then
    print "(/,a,/,a,/)",' Hint: to change the step legend text, create a file called', &
             '  '''//trim(fileprefix)//'.legend'' in the working directory, with one label per line'
 endif
 if (iPlotTitles) then
    if (.not.iPlotStepLegend) print "(a)"
    print "(a,/,a,/)",' To set the plot titles, create a file called', &
             '  '''//trim(fileprefix)//'.titles'' in the working directory, with one title per line'
 endif

 if (iaction.le.0 .or. iaction.gt.6) then
    !--format shape settings string
    if (nshapes.gt.0) then
       i2 = 2
       string = ': '
    else
       i2 = 0
       string = ' '
    endif
    do i=1,nshapes
       i1 = i2 + 1
       i2 = min(i1 + len_trim(labelshapetype(shape(i)%itype)),len(string))
       write(string(i1:i2),"(a)",iostat=ierr) trim(labelshapetype(shape(i)%itype)(1:i2-i1))
       if (i.lt.nshapes .and. i2.lt.len(string)) then
          write(string(i2:i2+1),"(', ')",iostat=ierr)
          i2 = i2 + 1
       endif
    enddo
    !--print menu
    print 20,print_logical(iPlotLegend),hposlegend,vposlegend,fjustlegend,trim(legendtext), &
             print_logical(iPlotTitles),hpostitle,vpostitle,fjusttitle, &
             print_logical(iPlotStepLegend), print_logical(iPlotScale),iPlotLegendOnlyOnPanel, &
             nshapes,trim(string)

20  format(' 0) exit ',/,                   &
        ' 1) time legend on/off/settings                   (',1x,a,1x,f5.2,1x,f5.2,1x,f5.2,1x,'"',a,'")',/, &
        ' 2) titles on/off/settings                        (',1x,a,1x,f5.2,1x,f5.2,1x,f5.2,')',/, &
        ' 3) legend for multiple steps per page on/off     (',1x,a,1x,')',/, &
        ' 4) plot scale |---| on coordinate plots          (',1x,a,1x,')',/, &
        ' 5) legend only on nth panel/first row/column     (',1x,i2,1x,')',/, &
        ' 6) annotate plot (e.g. arrow,square,circle,text) (',1x,i2,a,')')

    iaction = 0
    call prompt(' Enter option ',iaction,0,6)
 endif

 select case(iaction)
 case(1)
    call prompt('Plot time legend? ',iPlotLegend)
    print "(a)",'Time legend is '//print_logical(iPlotLegend)
    if (iPlotLegend) then
       print "(7(/,a),/)", &
       ' Example format strings: ', &
       '  t =              : this is the default format "t = 0.1 years"', &
       '  t = %t.5         : with time to 5 significant figures', &
       '  Time: %t dog-%ut : gives "Time: 0.1 dog-years"', &
       '  %(t + 2013)      : prints time offset by 2013', &
       '  %(t + 2013).5    : as above, to 5 sig. figs.', &
       '  %(t*100)         : multiplied by 100'

       call prompt('Enter legend text ',legendtext)

       print "(a)",'------ set legend position (can also be done interactively) --------'
       call prompt('Enter horizontal position as fraction of viewport', &
            hposlegend,0.0,1.0)
       call prompt('Enter vertical position in character heights from top',vposlegend)
       call prompt('Enter justification factor (0.0=left 1.0=right)',fjustlegend,0.0,1.0)
       call prompt_panelselect('legend',iPlotLegendOnlyOnPanel)
    endif
 case(2)
    print "(/,a,/,a,/)",' To set the plot titles, create a file called', &
             '  '''//trim(fileprefix)//'.titles'' in the working directory, with one title per line'
    call prompt('Use plot titles? ',iPlotTitles)
    print "(a)",'Titles are '//print_logical(iPlotTitles)
    if (iPlotTitles) then
       print "(a)",'------ set title position (can also be done interactively) --------'
       call prompt('Enter horizontal position as fraction of viewport', &
            hpostitle,0.0,1.0)
       call prompt('Enter vertical position in character heights above top',vpostitle)
       call prompt('Enter justification factor (0.0=left 1.0=right)',fjusttitle,0.0,1.0)
    endif
 case(3)
    iPlotStepLegend = .not.iPlotStepLegend
    print "(a)",'Step legend is '//print_logical(iPlotStepLegend)
    if (iPlotStepLegend) then
       print "(/,a,/,a,/)",' Hint: to change the step legend text, create a file called', &
                '  '''//trim(fileprefix)//'.legend'' in the working directory, with one label per line'
       print*,' press return to continue '
       read*
    endif
 case(4)
    call prompt('Plot scale on co-ordinate plots? ',iPlotScale)
    if (iPlotScale) then
       call prompt('Enter length of scale in the current x,y,z units ',dxscale)
       call prompt('Enter text to appear below scale (e.g. ''10 AU'')',scaletext)
       call prompt('Enter horizontal position as fraction of viewport', &
                   hposscale,0.0,1.0)
       call prompt('Enter vertical position in character heights above bottom',vposscale)
       if (nacross*ndown.gt.1) then
          call prompt('Enter which panel on the plotting page the scale should appear on '// &
               '(0=all co-ordinate plots)',iscalepanel,0,nacross*ndown)
       endif
    endif
 case(5)
    call prompt_panelselect('legend',iPlotLegendOnlyOnPanel)
 case(6)
    call submenu_shapes()
 end select

 return
end subroutine submenu_legend

end module settings_page
