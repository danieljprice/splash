module pagesetup
 implicit none
 public :: setpage
 real, parameter, private :: xlabeloffset = 3.0, ylabeloffset = 3.0
 
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
!         colourbarwidth : colour bar width in character heights
!         pagechange : change the physical page between plots
!
subroutine setpage(iplot,nx,ny,xmin,xmax,ymin,ymax,labelx,labely,title,  &
     just,axis,colourbarwidth,titleoffset,isamexaxis,ipagechange)
  implicit none
  integer, intent(in) :: iplot, nx, ny, just, axis
  real, intent(in) :: xmin, xmax, ymin, ymax, colourbarwidth, titleoffset
  character(len=*), intent(in) :: labelx, labely, title
  character(len=10) :: xopts, yopts
  logical, intent(in) :: ipagechange, isamexaxis
  real :: vptxmin,vptxmax,vptymin,vptymax,xch,ych

  if (ipagechange) then

     call pgpage
     !
     !--query the character height as fraction of viewport
     !
     call pgqcs(0,xch,ych)
     !
     !--default is to use whole viewport
     !
     vptxmin = 0.001
     vptxmax = 0.999
     vptymin = 0.001
     vptymax = 0.999
     !
     !--leave room for axes labels if necessary
     !
     if (axis.GE.0) then
        !
        !--PGPLOT seems to leave a (small) buffer around all plots
        !  so if only one plot label per page, can leave less space
        !  for the axes labels since they can go into the buffer region
        !
        if (nx.gt.1) then
           vptxmin = (ylabeloffset+2.0)*xch
        else
           vptxmin = (ylabeloffset+1.0)*xch
        endif
        if (ny.gt.1 .and. .not.isamexaxis) then
           vptymin = (xlabeloffset+2.0)*ych
        elseif (ny.gt.1) then
           vptymin = (xlabeloffset+0.25)*ych
        else
           vptymin = (xlabeloffset+1.0)*ych
        endif
     endif
     !--also leave room for title if necessary
     vptymax = vptymax - titleoffset*ych
     
     !--also leave room for colour bar if necessary
     if (colourbarwidth.GT.0.) then
        vptxmax = vptxmax - (colourbarwidth + 0.25)*ych
     endif

     call pgsvp(vptxmin,vptxmax,vptymin,vptymax)

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
        call pgmtxt('B',xlabeloffset,0.5,0.5,labelx)
     endif
     !
     !--always label y axis
     !
     call pgmtxt('L',ylabeloffset,0.5,0.5,labely)
     !
     !--always plot title
     !
     call pgmtxt('T',-titleoffset,0.5,0.5,title)

  endif
  
  return
end subroutine setpage

end module pagesetup
