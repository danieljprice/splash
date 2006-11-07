!-------------------------------------------------------------------------
! Module containing settings and options relating to vector plots
! includes default values of these options and submenu for changing them
!-------------------------------------------------------------------------
module settings_vecplot
 implicit none
 integer :: npixvec
 logical :: UseBackgndColorVecplot, iplotpartvec
 logical :: iVecplotLegend,iplotstreamlines
 real :: hposlegendvec,vposlegendvec

 namelist /vectoropts/ npixvec, UseBackgndColorVecplot,iplotpartvec,&
          iVecplotLegend,hposlegendvec,vposlegendvec,iplotstreamlines

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
  hposlegendvec = 0.02
  vposlegendvec = -1.5
  iplotstreamlines = .false. ! plot stream lines instead of arrows

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
 print 10,npixvec,print_logical(UseBackgndColorVecplot), &
          print_logical(iVecplotLegend),print_logical(iplotstreamlines)
10  format('--------------- vector plot options -------------------',/,&
           ' 0) exit ',/, &
           ' 1) change number of pixels                   (',i4,' )',/, &
           ' 2) use background colour for arrows          ( ',a,' )',/, &
           ' 3) vector plot legend settings               ( ',a,' )',/, &
           ' 4) plot stream/field lines instead of arrows ( ',a,' )')
 call prompt('enter option',ians,0,4)
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
    print*,'use background colour on vector plots = ', &
           print_logical(UseBackgndColorVecplot)
!------------------------------------------------------------------------
 case(3)
    call prompt('plot vector legend?',iVecplotLegend)
    if (iVecplotLegend) then
       print*,'note that the following settings can also be changed interactively'
       call prompt('Enter horizontal position as fraction of viewport', &
                   hposlegendvec,0.0,1.0)
       call prompt('Enter vertical position in character heights from top', &
                    vposlegendvec)
    endif
!------------------------------------------------------------------------
 case(4)
    iplotstreamlines = .not.iplotstreamlines
    print "(2(a,/))",' Note: the number of stream lines plotted is determined by', &
                     ' the "change number of contours" option in the r)ender menu'
    call prompt('use stream lines instead of arrows? ',iplotstreamlines)
 end select

 return
end subroutine submenu_vecplot

end module settings_vecplot
