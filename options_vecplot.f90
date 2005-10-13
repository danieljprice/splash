!-------------------------------------------------------------------------
! Module containing settings and options relating to vector plots
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_vecplot
 implicit none
 integer :: npixvec
 logical :: UseBackgndColorVecplot, iplotpartvec
 logical :: iVecplotLegend
 real :: hposlegendvec,vposlegendvec

 namelist /vectoropts/ npixvec, UseBackgndColorVecplot,iplotpartvec,&
          iVecplotLegend,hposlegendvec,vposlegendvec

contains

!---------------------------------------------
! set default values for these options
!---------------------------------------------
subroutine defaults_set_vecplot
  implicit none

  npixvec = 40        ! pixels in x direction on vector plots
  UseBackgndColorVecplot = .false. ! plot vector plot using black/white
  iplotpartvec = .true.   ! whether to plot particles on vector plot
  iVecplotLegend = .true.
  hposlegendvec = 0.05
  vposlegendvec = 2.0

  return
end subroutine defaults_set_vecplot

!----------------------------------------------------------------------
! sets options relating to vector plots
!----------------------------------------------------------------------
subroutine submenu_vecplot
 use prompting
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
                   hposlegendvec,0.0,1.0)
       call prompt('Enter vertical position in character heights from top', &
                    vposlegendvec,0.0)
    endif
 end select

 return
end subroutine submenu_vecplot

end module settings_vecplot
