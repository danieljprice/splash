!----------------------------------------------------------------------
! sets options relating to vector plots
!----------------------------------------------------------------------
subroutine options_vecplot
 use prompting
 use settings
 implicit none
 integer :: ians
  
 print 10,npixvec,iplotpartvec,UseBackgndColorVecplot
10  format(' 0) exit ',/, 		&
           ' 1) change number of pixels                   (',i4,' )',/, &
           ' 2) toggle particle plotting if no rendering  (',L1,' )',/, &
	   ' 3) toggle background/foreground colour       (',L1,' )')
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
    iplotpartvec = .not.iplotpartvec
    print*,'particle plotting = ',iplotpartvec,' if no renderings'
!------------------------------------------------------------------------
 case(3)
    UseBackgndColorVecplot = .not.UseBackgndColorVecplot
    print*,'use background colour on vector plots = ',UseBackgndColorVecplot
 end select

 return
end subroutine options_vecplot
