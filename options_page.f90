!----------------------------------------------------------------------
! submenu with options relating to page setup
!----------------------------------------------------------------------
subroutine options_page
 use settings
 use prompting
 implicit none
 integer :: iaction,i
 real :: papersizey
 logical :: ians
 
 iaction = 0
 papersizey = papersizex*aspectratio
 print 10,ipagechange,axes,papersizex,papersizey,nacross,ndown
10 format(' 0) exit ',/, 		&
        ' 1) toggle page change     (',L1,')',/, &
        ' 2) toggle axes            (',L1,')',/, &
        ' 3) change paper size      (',f5.2,1x,f5.2,')',/, &
        ' 4) change plots per page  (',i2,1x,i2,')')
 call prompt('enter option ',iaction,0,4)

 select case(iaction)
!------------------------------------------------------------------------
  case(1)
     ipagechange=.not.ipagechange
     print*,' Page changing = ',ipagechange
     return 	
!------------------------------------------------------------------------
  case(2)
     axes=.not.axes
     print *,' Axes = ',axes
     return
!------------------------------------------------------------------------
  case(3)
     print*,' 0) PGPLOT default '
     print*,' 1) small square movie   :  2.92 x 2.92 inches'
     print*,' 2) large/multiple movie :  5.85 x 5.85'
     print*,' 3) single small graph   :  5.85 x 4.13'
     print*,' 4) duo small graph      : 11.70 x 4.13'
     print*,' 5) duo graph            : 11.70 x 6.00'
     print*,' 6) Custom size '
     call prompt(' Enter option for paper size ',ipapersize,0,5)
     select case(ipapersize)
     case(1) 
        papersizex = 0.25*11.7
        aspectratio = 1.0
     case(2)
        papersizex = 0.5*11.7
        aspectratio = 1.0
     case(3) 
        papersizex = 0.5*11.7 
        aspectratio = 1./sqrt(2.)
     case(4)
        papersizex = 11.7
        aspectratio = 0.5/sqrt(2.)	
     case(5)
        papersizex = 11.7
	papersizey = 6.0
        aspectratio = papersizey/papersizex	
     case(6)
        call prompt(' x size (inches) ',papersizex,0.0,12.0)
        call prompt(' y size (inches) or aspect ratio (-ve)', &
             papersizey,-12.0,12.0)
        if (papersizey.lt.0.0) then
           aspectratio = abs(papersizey)
        else
           aspectratio = papersizey/papersizex
        endif
     case DEFAULT
        papersizex = 0.0	! no call to PGPAP if they are zero
        aspectratio = 0.0	
     end select
     return 	  
!------------------------------------------------------------------------
  case(4)
     call prompt('Enter number of plots across:',nacross,1,numplot)
     call prompt('Enter number of plots down  :',ndown,1,numplot)
     return	 
     
  end select
 
 return
end subroutine options_page
