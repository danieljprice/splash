module pagesetup
 implicit none
 public :: setpage
 
 private

contains
!
!--this subroutine determines the setup of the PGPLOT page
!  sorts out labelling of axes, positioning of windows etc
!  can be used as a replacement for PGENV and PGLABEL
!
!  inputs:
!         iplot  : position of current plot on page
!         nx     : number of plots across page
!         ny     : number of plots down page
!         xmin, xmax, ymin, ymax : plot limits
!         labelx, labely, title  : axes labels, plot title
!         just   : just=1 gives equal aspect ratios (same as in PGENV)
!         axis   : axes options (same as in PGENV)
!         pagechange : change the physical page between plots
!
subroutine setpage(iplot,nx,ny,xmin,xmax,ymin,ymax,labelx,labely,title,  &
     just,axis,isamexaxis,isameyaxis,ipagechange)
  implicit none
  integer, intent(in) :: iplot, nx, ny, just, axis
  real, intent(in) :: xmin, xmax, ymin, ymax
  character(len=*), intent(in) :: labelx, labely, title
  character(len=10) :: xopts, yopts
  logical, intent(in) :: ipagechange, isamexaxis, isameyaxis

  if (ipagechange) then

     call pgpage
     if (nx*ny.gt.1) then
        if (axis.lt.0) then
           call pgsvp(0.02,0.98,0.02,0.98) ! if no axes use full viewport
        else
           !call pgsvp(0.25,0.98,0.15,0.95)
           call pgsvp(0.25,0.98,0.17,0.93)
        endif
     else
        if (axis.lt.0) then
           call pgsvp(0.01,0.99,0.04,0.96) ! if no axes use full viewport
        else
           call pgsvp(0.1,0.9,0.1,0.9)
        endif
     endif

     if (just.eq.1) then
        call pgwnad(xmin,xmax,ymin,ymax) ! pgwnad does equal aspect ratios
     else
        call pgswin(xmin,xmax,ymin,ymax)
     endif
!
!--set plot axes (options are exactly as in PGENV, with axis=-3 added)
!
     yopts = '*'
     select case(axis)
     case(-3)
        xopts = 'BCST'
     case(-2)
        xopts = ' '
     case(-1)
        xopts = 'BC'
     case(0)
        xopts = 'BCNST'
     case(1)
        xopts = 'ABCNST'
     case(2)
        xopts = 'ABCGNST'
     case(10)
        xopts = 'BCNSTL'
        yopts = 'BCNST'
     case(20)
        xopts = 'BCNST'
        yopts = 'BCNSTL'
     case(30)
        xopts = 'BCNSTL'
        yopts = 'BCNSTL'
     case default
        CALL GRWARN('PGENV: illegal AXIS argument.')
        xopts = 'BCNST'
     end select
     if (yopts.eq.'*') yopts = xopts

     call pgbox(xopts,0.0,0,yopts,0.0,0)

  elseif (iplot.eq.1) then ! if would be changing page, instead go back to
     call pgpanl(1,1)      !                                   first panel
  elseif (nx*ny.gt.1) then ! change to next panel, regardless of ipagechange
     call pgpage
  endif
  
  !---------------------------------
  ! set plot limits and label plot
  !---------------------------------
    
  if (just.eq.1) then
     call pgwnad(xmin,xmax,ymin,ymax)  ! repeated for when not called above
  else
     call pgswin(xmin,xmax,ymin,ymax)
  endif
  
  !--label plot
  if (axis.ge.0) then
     !
     !--label x axis only if on last row
     !  or if x axis quantities are different
     !
     if (((ny*nx-iplot).lt.nx).or.(.not.isamexaxis)) then
        call pglabel(trim(labelx),' ',trim(title))
     endif
     !
     !--always label y axis
     !
     call pgmtxt('l',3.0,0.5,0.5,labely)
     !!             call pglabel(' ',labely,trim(titlex))     
  endif
  
  return
end subroutine setpage

end module pagesetup
