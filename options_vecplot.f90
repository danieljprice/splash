!----------------------------------------------------------------------
! sets options relating to vector plots
!----------------------------------------------------------------------
subroutine options_vecplot
 use prompting
 use settings_vecplot
 implicit none
 integer :: ians
  
 ians = 0
 print 10,npixvec,UseBackgndColorVecplot,iVecplotLegend
10  format(' 0) exit ',/, &
           ' 1) change number of pixels                   (',i4,' )',/, &
           ' 2) toggle background/foreground colour       (',L1,' )',/, &
           ' 3) vector plot legend settings               (',L1,' )')
 call prompt('enter option',ians,0,3)
!
!--options
!
 select case(ians)
!------------------------------------------------------------------------
 case(1)
    call prompt('enter number of pixels',npixvec,1,1000) 
!------------------------------------------------------------------------
 case(2)
    UseBackgndColorVecplot = .not.UseBackgndColorVecplot
    print*,'use background colour on vector plots = ',UseBackgndColorVecplot
!------------------------------------------------------------------------
 case(3)
    call prompt('plot vector legend?',iVecplotLegend)
    if (iVecplotLegend) then
       print*,'note that the following settings can also be changed interactively'
       call prompt('Enter horizontal position as fraction of viewport', &
                   hposlegendvec,-0.1,1.1)
       call prompt('Enter vertical position in character heights from top', &
                    vposlegendvec)
    endif
 end select

 return
end subroutine options_vecplot
