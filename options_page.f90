!----------------------------------------------------------------------
! submenu with options relating to page setup
!----------------------------------------------------------------------
subroutine options_page
 use settings
 use prompting
 implicit none
 integer :: iaction
 real :: papersizey

 iaction = 0
 papersizey = papersizex*aspectratio
 print 10,ipagechange,iaxis,papersizex,papersizey,nacross,ndown,tile, &
          hpostitle,vpostitle,hposlegend,vposlegend,animate
10 format(' 0) exit ',/,                   &
        ' 1) toggle page change     (',L1,')',/, &
        ' 2) toggle axes            (',i2,')',/, &
        ' 3) change paper size      (',f5.2,1x,f5.2,')',/, &
        ' 4) change plots per page  (',i2,1x,i2,')',/, &
         ' 5) toggle plot tiling     (',L1,')',/, & 
         ' 6) adjust title position  (',f5.2,1x,f4.1,')',/, &
         ' 7) adjust legend position (',f5.2,1x,f4.1,')',/, &
         ' 8) toggle animate         (',L1,')')         
 call prompt('enter option ',iaction,0,8)

 select case(iaction)
!------------------------------------------------------------------------
  case(1)
     ipagechange=.not.ipagechange
     print*,' Page changing = ',ipagechange
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
     call prompt('Enter horizontal position as fraction of viewport', &
          hpostitle,0.0,1.0)
     call prompt('Enter vertical position in character heights above top',vpostitle)
     call prompt('Enter justification factor (0.0=left 1.0=right)',fjusttitle,0.0,1.0)
     return
!------------------------------------------------------------------------
  case(7)
     call prompt('Enter horizontal position as fraction of viewport', &
          hposlegend,0.0,1.0)
     call prompt('Enter vertical position in character heights from top',vposlegend)
     return
!------------------------------------------------------------------------
  case(8)
     animate = .not.animate
     print*,'animate = ',animate  
  end select
 
 return
end subroutine options_page
