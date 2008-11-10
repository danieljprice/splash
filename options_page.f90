!-------------------------------------------------------------------------
! Module containing page settings and options
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_page
 use settings_limits, only:iadapt,iadaptcoords
 implicit none
 integer :: iaxis,nacross,ndown,ipapersize,nstepsperpage,linewidth,iscalepanel
 integer :: iPlotLegendOnlyOnPanel
 logical :: iColourEachStep,iChangeStyles,tile,interactive,nomenu
 logical :: iPlotLegend,iPlotStepLegend,iPlotTitles
 logical :: iPlotScale,iUseBackgroundColourForAxes,usesquarexy
 real :: papersizex,aspectratio
 real :: hposlegend,vposlegend,fjustlegend,hpostitle,vpostitle,fjusttitle
 real :: charheight
 real :: dxscale,hposscale,vposscale
 character(len=20) :: colour_fore, colour_back, legendtext, scaletext
 character(len=60) :: device

 namelist /pageopts/ iaxis,nacross,ndown,interactive,iadapt,iadaptcoords, &
   nstepsperpage,iColourEachStep,iChangeStyles,tile,ipapersize,papersizex,aspectratio, &
   iPlotLegend,iPlotStepLegend,hposlegend,vposlegend,iPlotTitles,hpostitle, &
   vpostitle,fjusttitle,legendtext,colour_fore,colour_back,charheight,linewidth,&
   fjustlegend,iPlotLegendOnlyOnPanel, &
   iPlotScale,dxscale,scaletext,hposscale,vposscale,iscalepanel,iUseBackgroundColourForAxes, &
   usesquarexy

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_page
  use shapes, only:defaults_set_shapes
  implicit none

  interactive = .true.     ! default for interactive mode
  iaxis = 0                ! turns axes off/on
  nstepsperpage = 1
  iColourEachStep = .true. ! change colours if nstepsperpage > 1
  iChangeStyles = .false.  ! change marker/ line styles if nstepsperpage > 1
  tile = .true.
  nacross = 1           ! number of plots across page
  ndown = 1             ! number of plots down page
  ipapersize = 0        ! paper size option
  papersizex = 0.0      ! size of x paper (no call to PGPAP if zero)
  aspectratio = 0.0     ! aspect ratio of paper (no call to PGPAP if zero)
  
  iPlotLegend = .true.  ! whether or not to plot legend
  iPlotStepLegend = .false. ! timestep legend
  hposlegend = 0.95     ! horizontal legend position as fraction of viewport
  vposlegend = 2.0      ! vertical legend position in character heights
  fjustlegend = 1.0    ! justification factor for legend
  legendtext = 't='
  iPlotLegendOnlyOnPanel = 0
  
  iPlotTitles = .false.  ! whether or not to plot titles
  hpostitle = 0.5       ! horizontal title position as fraction of viewport
  vpostitle = 1.0       ! vertical title position in character heights
  fjusttitle = 0.5      ! justification factor for title
  colour_fore = ' '
  colour_back = ' '
  charheight = 1.0    ! PGPLOT character height
  linewidth = 0       ! PGPLOT line width

  iPlotScale = .false.
  hposscale = 0.5
  vposscale = 1.0
  dxscale = 1.0
  scaletext = '1 unit'
  iscalepanel = 0
  
  usesquarexy = .true. ! spatial dimensions have same scale
  call defaults_set_shapes

  return
end subroutine defaults_set_page

!---------------------------------------------
! changed default values for evsplash
!---------------------------------------------
subroutine defaults_set_page_ev
  implicit none

  nstepsperpage = 1000
  iColourEachStep = .true. ! change colours if nstepsperpage > 1
  iChangeStyles = .true.   ! change marker/ line styles if nstepsperpage > 1
  iPlotLegend = .true.  ! whether or not to plot legend
  iPlotStepLegend = .true. ! timestep legend
  hposlegend = 0.1     ! horizontal legend position as fraction of viewport
  vposlegend = 2.0      ! vertical legend position in character heights
  fjustlegend = 0.0    ! justification factor for legend
  
  return
end subroutine defaults_set_page_ev

!----------------------------------------------------------------------
! submenu with options relating to page setup
!----------------------------------------------------------------------
subroutine submenu_page(ichoose)
 use params, only:maxplot
 use prompting, only:prompt,print_logical
 implicit none
 integer, intent(in) :: ichoose
 integer :: iaction,ierr,ntries
 real :: papersizey

 iaction = ichoose
 
 papersizey = papersizex*aspectratio
 print "(a)",'---------------- page setup options -------------------'
 
 if (iaction.le.0 .or. iaction.gt.9) then
    print 10,nstepsperpage,iaxis,papersizex,papersizey,nacross,ndown,print_logical(tile), &
             trim(print_logical(usesquarexy)),charheight,linewidth,trim(print_logical(interactive))
10 format( &
          ' 0) exit ',/,                   &
          ' 1) plot n steps on top of each other   (n =',i4,')',/, &
          ' 2) axes options                        (',i2,')',/, &
          ' 3) change paper size                   (',f5.2,1x,f5.2,')',/, &
          ' 4) subdivide page into panels          (',i2,1x,'x',1x,i2,', tiling is ',a,')',/, &
          ' 5) spatial dimensions have same scale  ( ',a,' )',/,&
          ' 6) set character height                (',f4.1,')',/,&
          ' 7) adjust line width                   (',i2, ')',/,&
          ' 8) interactive mode on/off             ( ',a,' )',/,&
          ' 9) set foreground/background colours ')
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
!!        if (.not.iColourEachStep) icolourthisstep = 1
        call prompt('Use different markers/line style for each? ',iChangeStyles)
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
     print*,'10 : draw box and label X-axis logarithmically;'
     print*,'20 : draw box and label Y-axis logarithmically;'
     print*,'30 : draw box and label both axes logarithmically.'
     call prompt('enter axis option ',iaxis,-4,30)
     print *,' axis = ',iaxis
     return
!------------------------------------------------------------------------
  case(3)
     print*,' 0) PGPLOT default '
     print*,' 1) small square         :  2.92 x 2.92 inches'
     print*,' 2) medium square        :  5.85 x 5.85'
     print*,' 3) large square         :  8.00 x 8.00'
     print*,' 4) single small graph   :  5.85 x 4.13'
     print*,' 5) duo small graph      : 11.70 x 4.13'
     print*,' 6) duo graph            : 11.70 x 6.00'
     print*,' 7) Custom size '
     call prompt(' Enter option for paper size ',ipapersize,0,7)
     select case(ipapersize)
     case(1) 
        papersizex = 0.25*11.7
        aspectratio = 1.0
     case(2)
        papersizex = 0.5*11.7
        aspectratio = 1.0
     case(3)
        papersizex = 8.0
        aspectratio = 1.0     
     case(4) 
        papersizex = 0.5*11.7 
        aspectratio = 1./sqrt(2.)
     case(5)
        papersizex = 11.7
        aspectratio = 0.5/sqrt(2.)         
     case(6)
        papersizex = 11.7
         papersizey = 6.0
        aspectratio = papersizey/papersizex
     case(7)
        call prompt(' x size (inches) ',papersizex,0.0)
        call prompt(' y size (inches) or aspect ratio (-ve)',papersizey)
        if (papersizey.lt.0.0) then
           aspectratio = abs(papersizey)
        else
           aspectratio = papersizey/papersizex
        endif
     case DEFAULT
        papersizex = 0.0         ! no call to PGPAP if they are zero
        aspectratio = 0.0         
     end select
     return            
!------------------------------------------------------------------------
  case(4)
     call prompt('Enter number of plots across (columns):',nacross,1,maxplot)
     call prompt('Enter number of plots down   (rows):',ndown,1,maxplot/nacross)
     if (nacross*ndown.gt.1) then
        call prompt('Tile plots on the page where possible?',tile)
     endif
     return
!------------------------------------------------------------------------
  case(5)
     usesquarexy = .not.usesquarexy
     print "(a)",' Same scale for spatial dimensions is '//print_logical(usesquarexy)
!------------------------------------------------------------------------
  case(6)
     call prompt('Enter PGPLOT character height ',charheight,0.1,10.)
     return
!------------------------------------------------------------------------
  case(7)
     print "(3(/,a))",' Setting line width to 0 means automatic line width choice:', &
                    ' This gives width = 2 for vector devices (/ps,/cps etc)', &
                    '        and width = 1 elsewhere (e.g. for pixel devices)'
     print*
     call prompt('Enter PGPLOT line width (0=auto)',linewidth,0)
     return
!------------------------------------------------------------------------
  case(8)
     interactive = .not.interactive
     print "(a)",' Interactive mode is '//print_logical(interactive)
!------------------------------------------------------------------------
  case(9)
     ierr = 1
     ntries = 1
     !--open null device so that colours can be recognised
     call pgopen('/null')
     do while (ierr /= 0 .and. ntries.le.3)
        call prompt('Enter background colour (by name, e.g. "black") ',colour_back)
        call pgscrn(0,colour_back,ierr)
        ntries = ntries + 1
     enddo
     ierr = 1
     ntries = 1
     do while (ierr /= 0 .and. ntries.le.3)
        call prompt('Enter foreground colour (by name, e.g. "white") ',colour_fore)
        call pgscrn(1,colour_fore,ierr)
        ntries = ntries + 1
     enddo
     call pgclos

     print "(4(/,a))",' Overlaid (that is, drawn inside the plot borders) axis ', &
                      ' ticks, legend text and titles are by default plotted in ', &
                      ' the foreground colour [ie. '//trim(colour_fore)//'].'
     
     if (len_trim(colour_back).gt.0 .and. len_trim(colour_fore).gt.0) then
        call prompt('Do you want to plot these in background colour [ie. '//trim(colour_back)//'] instead ?',&
                    iUseBackgroundColourForAxes)     
     else
        call prompt('Do you want to plot these in background colour instead? ',iUseBackgroundColourForAxes)
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
 use shapes, only:nshapes,labelshapetype,shape,submenu_shapes
 implicit none
 integer, intent(in) :: ichoose
 integer :: iaction,i,ierr,i1,i2
 character(len=50) :: string

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
       if (i.lt.nshapes) then
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
       call prompt('Enter legend text ',legendtext)

       print "(a)",'------ set legend position (can also be done interactively) --------'
       call prompt('Enter horizontal position as fraction of viewport', &
            hposlegend,0.0,1.0)
       call prompt('Enter vertical position in character heights from top',vposlegend)
       call prompt('Enter justification factor (0.0=left 1.0=right)',fjustlegend,0.0,1.0)         

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
    print "(/,'  0 : plot legend on every panel ',/,"// &
           "' -1 : plot legend on first row only ',/,"// &
           "' -2 : plot legend on first column only ',/,"// &
           "'  n : plot legend on nth panel only ')"

    call prompt('Enter selection ',iPlotLegendOnlyOnPanel,-2)
 case(6)
    call submenu_shapes()
 end select

 return
end subroutine submenu_legend

end module settings_page
