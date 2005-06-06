!-------------------------------------------------------------------------
! Module containing page settings and options
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_page
 use settings_limits, only:iadapt,iadaptcoords
 implicit none
 integer :: iaxis,nacross,ndown,ipapersize,nstepsperpage
 logical :: iColourEachStep,tile,interactive
 logical :: iPlotLegend,iPlotTitles
 real :: papersizex,aspectratio
 real :: hposlegend,vposlegend,hpostitle,vpostitle,fjusttitle
 real :: charheightmm
 character(len=20) :: colour_fore, colour_back, legendtext

 namelist /pageopts/ iaxis,nacross,ndown,interactive,iadapt,iadaptcoords, &
   nstepsperpage,iColourEachStep,tile,ipapersize,papersizex,aspectratio, &
   iPlotLegend,hposlegend,vposlegend,iPlotTitles,hpostitle,vpostitle,fjusttitle,legendtext, &
   colour_fore, colour_back, charheightmm

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_page
  implicit none

  interactive = .true.     ! default for interactive mode
  iaxis = 0                ! turns axes off/on
  nstepsperpage = 1
  tile = .false.
  nacross = 1           ! number of plots across page
  ndown = 1             ! number of plots down page
  ipapersize = 0        ! paper size option
  papersizex = 0.0      ! size of x paper (no call to PGPAP if zero)
  aspectratio = 0.0     ! aspect ratio of paper (no call to PGPAP if zero)
  
  iPlotLegend = .true.  ! whether or not to plot legend
  hposlegend = 0.75     ! horizontal legend position as fraction of viewport
  vposlegend = 2.0      ! vertical legend position in character heights
  legendtext = 't='
  
  iPlotTitles = .true.  ! whether or not to plot titles
  hpostitle = 0.5       ! horizontal title position as fraction of viewport
  vpostitle = 1.0       ! vertical title position in character heights
  fjusttitle = 0.5      ! justification factor for title
  colour_fore = ' '
  colour_back = ' '
  charheightmm = 4.0    ! character height in mm

  return
end subroutine defaults_set_page

!----------------------------------------------------------------------
! submenu with options relating to page setup
!----------------------------------------------------------------------
subroutine submenu_page
 use settings_data, only:numplot
 use prompting
 implicit none
 integer :: iaction,ierr,ntries
 real :: papersizey

 iaction = 0
 papersizey = papersizex*aspectratio
 print 10,nstepsperpage,iaxis,papersizex,papersizey,nacross,ndown,tile, &
          iPlotTitles,hpostitle,vpostitle,iPlotLegend,hposlegend,vposlegend, &
          charheightmm
10 format(' 0) exit ',/,                   &
        ' 1) change steps per page  (',i2,')',/, &
        ' 2) axes options           (',i2,')',/, &
        ' 3) change paper size      (',f5.2,1x,f5.2,')',/, &
        ' 4) change plots per page  (',i2,1x,i2,')',/, &
         ' 5) toggle plot tiling     (',L1,')',/, & 
         ' 6) title options          (',L1,1x,f5.2,1x,f4.1,')',/, &
         ' 7) legend options         (',L1,1x,f5.2,1x,f4.1,')',/, &
         ' 8) set character height   (',f4.1,')',/,&
         ' 9) set foreground/background colours ')
 call prompt('enter option ',iaction,0,9)

 select case(iaction)
!------------------------------------------------------------------------
  case(1)
     call prompt('Enter number of timesteps per page',nstepsperpage,1)
     print*,'Plotting up to ',nstepsperpage,' timesteps per page'
     if (nstepsperpage.gt.1) then
        if (iadapt) then
           print*,'(note that adaptive plot limits are now off)'
           iadapt = .false.
        endif
        if (nstepsperpage.le.14) then
           call prompt('Use different colours for each step?',iColourEachStep)
        else
           iColourEachStep = .false.
           print*,'(and that there are too many steps per page to change colours)'
        endif
     endif
     return          
!------------------------------------------------------------------------
  case(2)
     print*,'-3 : same as AXIS=-1, but also draw tick marks;'
     print*,'-2 : draw no box, axes or labels;'
     print*,'-1 : draw box only;'
     print*,' 0 : draw box and label it with coordinates;'
     print*,' 1 : same as AXIS=0, but also draw the coordinate axes (X=0, Y=0);'
     print*,' 2 : same as AXIS=1, but also draw grid lines at major increments of the coordinates;'


!     print*,'10 : draw box and label X-axis logarithmically;'
!     print*,'20 : draw box and label Y-axis logarithmically;'
!     print*,'30 : draw box and label both axes logarithmically.'
     call prompt('enter axis option ',iaxis,-3,2)
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
        call prompt(' x size (inches) ',papersizex,0.0,12.0)
        call prompt(' y size (inches) or aspect ratio (-ve)', &
             papersizey,-12.0,12.0)
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
     call prompt('Enter number of plots across:',nacross,1,numplot)
     call prompt('Enter number of plots down  :',ndown,1,numplot)
     return
!------------------------------------------------------------------------
  case(5)
     tile = .not.tile
     print*,'tile plots = ',tile
!------------------------------------------------------------------------
  case(6)
     call prompt('Plot titles?',iPlotTitles)
     if (iPlotTitles) then
        call prompt('Enter horizontal position as fraction of viewport', &
             hpostitle,0.0,1.0)
        call prompt('Enter vertical position in character heights above top',vpostitle)
        call prompt('Enter justification factor (0.0=left 1.0=right)',fjusttitle,0.0,1.0)
     endif
     return
!------------------------------------------------------------------------
  case(7)
     call prompt('Plot time legend?',iPlotLegend)
     if (iPlotLegend) then
        call prompt('Enter horizontal position as fraction of viewport', &
             hposlegend,0.0,1.0)
        call prompt('Enter vertical position in character heights from top',vposlegend)
        call prompt('Enter legend text ',legendtext)
     endif
     return
!------------------------------------------------------------------------
  case(8)
     call prompt('Enter page character height in mm ',charheightmm,1.,50.)
     return
!------------------------------------------------------------------------
  case(9)
     ierr = 1
     ntries = 1
     !--open null device so that colours can be recognised
     call pgopen('/null')
     do while (ierr /= 0 .and. ntries.le.3)
        call prompt('Enter background colour (by name) ',colour_back)
        call pgscrn(0,colour_back,ierr)
        ntries = ntries + 1
     enddo
     ierr = 1
     ntries = 1
     do while (ierr /= 0 .and. ntries.le.3)
        call prompt('Enter foreground colour (by name) ',colour_fore)
        call pgscrn(1,colour_fore,ierr)
        ntries = ntries + 1
     enddo
     call pgclos
    
     return
  end select
 
 return
end subroutine submenu_page

end module settings_page
