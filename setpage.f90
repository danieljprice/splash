module pagesetup
 implicit none
 public :: setpage, redraw_axes, setpage2, set_exactpixelboundaries
 real, parameter, public :: xlabeloffset = 3.0, ylabeloffset = 4.5
 
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
        !--leave a bit of buffer space if more than one plot on page
        !
        if (nx.gt.1) then
           vptxmin = (ylabeloffset+1.5)*xch
        else
           vptxmin = (ylabeloffset+1.0)*xch
        endif
        if (ny.gt.1 .and. .not.isamexaxis) then
           vptymin = (xlabeloffset+1.5)*ych
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
        vptxmax = vptxmax - (colourbarwidth + 1.6)*xch
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
     case(-4)
        xopts = 'BCT'
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
     case(3)
        xopts = 'BCNST'
     case(4)
        xopts = 'ABCNST'
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
  if (axis.ge.0 .and. axis.ne.3 .and. axis.ne.4) then
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

!
!--this subroutine is a cut down version of the above, which ONLY redraws the axes
!  (so that axes can be redrawn on *top* of what has been plotted).
!
!  inputs:
!         axis   : axes options (same as in PGENV, with axis=-3 added)
!

subroutine redraw_axes(iaxis)
  implicit none
  integer, intent(in) :: iaxis
  character(len=10) :: xopts, yopts
!
!--set plot axes (options are exactly as in PGENV, with axis=-3 added)
!
  yopts = '*'
  select case(iaxis)
  case(-4)
     xopts = 'BCT'
  case(-3)
     xopts = 'BCST'
  case(-2)
     xopts = ' '
  case(-1)
     xopts = 'BC'
  case(0)
     xopts = 'BCST'
  case(1)
     xopts = 'ABCST'
  case(2)
     xopts = 'ABCGST'
  case(10)
     xopts = 'BCSTL'
     yopts = 'BCST'
  case(20)
     xopts = 'BCST'
     yopts = 'BCSTL'
  case(30)
     xopts = 'BCSTL'
     yopts = 'BCSTL'
  case default
     print*,'redraw_axes: illegal AXIS argument.'
     xopts = 'BCST'
  end select
  if (yopts.eq.'*') yopts = xopts

  call pgbox(xopts,0.0,0,yopts,0.0,0)

  return
end subroutine redraw_axes

!
!--this subroutine determines the setup of the PGPLOT page
!  sorts out labelling of axes, positioning of windows etc
!  can be used as a replacement for PGENV and PGLABEL
!  
!  divides up a single page into subpanels
!  
!
!  option to tile graphs appropriately on a page in pgplot
!  divides up a single panel into subpanels, with a margin at the edge
!  should replace the call to pgenv and pglabel
!
!  for tiled plots the page setup looks like this:
!
!    |   |   |   | 
!  --+---+---+---+--
!    | 1 | 2 | 3 |
!  --+---+---+---+--
!    | 4 | 5 | 6 |
!  --+---+---+---+--
!    |   |   |   |
!
!  (ie. with margins in x and y)
!  note that we divide up a single panel, so pgbeg should be called with nx=1,ny=1
!
!  arguments:
!   iplot  : current plot number
!   nx     : number of panels in x direction
!   ny     : number of panels in y direction
!   xmin,  : xmax, ymin, ymax : plot limits (if tiled should be same for all plots)
!   labelx : x axis label (if tiled should be same for all plots)
!   labely : y axis label (if tiled should be same for all plots)
!   title  : current plot title (can differ between plots)
!   just   : just=1 gives equal aspect ratios (same as in pgenv)
!   axis   : axes options (same as in pgenv with a few extra)
!   vmarginleft,right,bottom,top : initial margin sizes (% of page (if tiled) or panel (if not))
!            (default should be zero for these)
!   tile   : assumes all plots can be tiled
!
!  This version by Daniel Price, July 2006
!
  subroutine setpage2(iplotin,nx,ny,xmin,xmax,ymin,ymax,labelx,labely,title,just,axis, &
                      vmarginleftin,vmarginrightin,vmarginbottomin,vmargintopin, &
                      colourbarwidth,titleoffset,isamexaxis,tile) 
  implicit none
  integer, intent(in) :: iplotin,nx,ny,just,axis
  real, intent(in) :: xmin, xmax, ymin, ymax, colourbarwidth, titleoffset
  real, intent(in) :: vmarginleftin,vmarginrightin,vmargintopin,vmarginbottomin
  character(len=*), intent(in) :: labelx,labely,title
  logical, intent(in) :: isamexaxis,tile
  integer iplot,ix,iy
  real vptsizeeffx,vptsizeeffy,panelsizex,panelsizey
  real vmargintop,vmarginbottom,vmarginleft,vmarginright
  real vptxmin,vptxmax,vptymin,vptymax
  real aspectratio,devaspectratio,x1,x2,y1,y2
  real xch,ych
  character(len=10)  :: xopts, yopts
  logical, parameter :: useexactpixelboundaries = .true.
!
! new page if iplot > number of plots on page
!
  if (iplotin.gt.nx*ny) then
     if (mod(iplotin,nx*ny).eq.1) call pgpage
     iplot = iplotin - (nx*ny)*((iplotin-1)/(nx*ny))
  elseif (iplotin.le.0) then
     return
  else
     iplot = iplotin
  endif
!
! check for errors in input
!      
  if (nx.le.0 .or. ny.le.0) return
!
! for tiled plots, adjust effective viewport size if just=1 and graphs are not square
!      
  if (tile .and. just.eq.1) then
     if (ymax.eq.ymin) then
        print*,'setpage: error tiling plots: ymax=ymin'
        return
     endif
!
! query the current aspect ratio of the device and set aspect ratio appropriately
!
     call pgqvsz(3,x1,x2,y1,y2)
     devaspectratio = (x2-x1)/(y2-y1)
     aspectratio = ((xmax-xmin)*nx)/((ymax-ymin)*ny)/devaspectratio
  else
     aspectratio = 1.0
  endif
!
! set positions of x and y labels in units of character height from edge
!
!  xlabeloffset = 3.0
!  ylabeloffset = 4.5
!
! query the character height as fraction of viewport
!
  call pgqcs(0,xch,ych)
!
! set margin size in units of viewport dimensions
! allow enough room for the plot labels if they are drawn
! nb: pgplot sets the character height as some fraction of the smallest
!     dimension
!
! for tiled plots, these margins apply to the whole page
! otherwise, these are applied to each panel individually
!
  vmargintop = vmargintopin
  vmarginright = vmarginrightin
  if (axis.ge.0) then
     vmarginleft = vmarginleftin + (ylabeloffset+1.0)*xch
     vmarginbottom = vmarginbottomin + (xlabeloffset+1.0)*ych

     if (.not.tile) then
        if (ny.gt.1 .and. .not.isamexaxis) then
           vmarginbottom = vmarginbottom + 0.5*ych
        elseif (ny.gt.1) then
           vmarginbottom = vmarginbottom + 0.25*ych
        endif
     endif
  else
     vmarginleft = vmarginleftin
     vmarginbottom = vmarginbottomin
  endif
!
!--set size of each panel
!
  ix = iplot - ((iplot-1)/nx)*nx
  iy = (iplot-1)/nx + 1

  if (tile) then
     !--also leave room for title if necessary
     if (titleoffset.ge.0.) then
        vmargintop = vmargintop + (titleoffset+0.75)*ych
     endif

     !
     ! effective viewport size = size - margins (only used for tiled
     !
     vptsizeeffx = 1.0 - vmarginright - vmarginleft
     vptsizeeffy = 1.0 - vmargintop - vmarginbottom
     !     reduce x or y size if just=1 to get right aspect ratio
     if (aspectratio.le.1.0) then
        if (aspectratio*vptsizeeffy.lt.vptsizeeffx) then
           vptsizeeffx = aspectratio*vptsizeeffy
        !  but this could still be bigger than the margins allow...
        else
           vptsizeeffy = vptsizeeffx/aspectratio
        endif
     elseif (aspectratio.gt.1.0) then
        if (vptsizeeffx/aspectratio.lt.vptsizeeffy) then
           vptsizeeffy = vptsizeeffx/aspectratio     
        !  but this could still be bigger than the margins allow...
        else
           vptsizeeffx = vptsizeeffy*aspectratio
        endif
     endif

     panelsizex = vptsizeeffx/nx
     panelsizey = vptsizeeffy/ny 
!         print*,ix,iy,nx,ny
!         print*,panelsizex,panelsizey,vptsizeeffx,vptsizeeffy

!     print*,'margins = ',vmarginleft,vmarginright
     vptxmin = vmarginleft + (ix-1)*panelsizex
     vptxmax = vptxmin + panelsizex
     vptymax = (1.0 - vmargintop) - (iy-1)*panelsizey
     vptymin = vptymax - panelsizey
  else
     !--use full page for non-tiled plots, then set margins inside each panel
     panelsizex = 1.0/nx
     panelsizey = 1.0/ny
     vptxmin = (ix-1)*panelsizex + vmarginleft
     vptxmax = ix*panelsizex - vmarginright
     vptymax = 1.0 - (iy-1)*panelsizey - vmargintop
     vptymin = 1.0 - iy*panelsizey + vmarginbottom


     !--also leave room for title if necessary
     if (titleoffset.ge.0.) then
        vptymax = vptymax - (titleoffset+0.75)*ych
     endif
     
     !--also leave room for colour bar if necessary
     if (colourbarwidth.GT.0.) then
        vptxmax = vptxmax - (colourbarwidth + 1.6)*xch
     endif

  endif
  !    print*,vptxmin,vptxmax,vptymin,vptymax
!
! set viewport
!
 !print*,'setting ',vptxmin,vptxmax,vptymin,vptymax
 call pgsvp(vptxmin,vptxmax,vptymin,vptymax)
!
! set axes
!
 if (just.eq.1) then
    call pgwnad(xmin,xmax,ymin,ymax)
 else
    call pgswin(xmin,xmax,ymin,ymax)
 endif
!
! adjust viewport to lie exactly on pixel boundaries
!
 if (useexactpixelboundaries) call set_exactpixelboundaries()
!
! option to return before actually doing anything
!      
  if (trim(title).eq.'NOPGBOX') return
!
! set options for call to pgbox (draws axes) and label axes where appropriate
! (options are exactly as in pgenv apart from axis=-3,-4 which i have added)
!
  yopts = '*'
  select case(axis)
  case(-4)
     xopts = 'BCT'
  case(-3)
     xopts = 'BCST'
  case(-2)
     xopts = ' '
  case(-1)
    xopts = 'BC'
  case(0)
    xopts = 'BCST'
  case(1)
    xopts = 'ABCST'
  case(2)
    xopts = 'ABCGST'
  case(10)
    xopts = 'BCSTL'
    yopts = 'BCST'
  case(20)
    xopts = 'BCST'
    yopts = 'BCSTL'
  case(30)
    xopts = 'BCSTL'
    yopts = 'BCSTL'
  case default
    print*,'setpage: illegal axis argument.'
    xopts = 'BCNST'
  end select
  if (yopts.eq.'*') yopts = xopts
!
! label plot
!
  if (tile) then
     !
     ! decide whether to number and label the y axis
     !      
     if (ix.eq.1 .and. axis.ge.0) then
        yopts = '1VN'//trim(yopts)
        call pgmtxt('L',ylabeloffset,0.5,0.5,labely)
     elseif (axis.ge.0) then
        yopts = trim(yopts)//'N'
     endif  
     !
     ! decide whether to number and label the x axis
     !      
     if (iy.eq.ny .and. axis.ge.0) then
        xopts = 'N'//trim(xopts)
        call pgmtxt('B',xlabeloffset,0.5,0.5,labelx)
     endif
     !
     ! plot the title if inside the plot boundaries
     !
     if (titleoffset.lt.0.) call pgmtxt('t',-titleoffset,0.96,1.0,title)

  elseif (axis.ge.0) then
     !
     !--label x axis only if on last row
     !  or if x axis quantities are different
     !
     if (((ny*nx-iplot).lt.nx).or.(.not.isamexaxis)) then
        call pgmtxt('B',xlabeloffset,0.5,0.5,labelx)
     endif
     !--always plot numbers
     xopts = 'N'//trim(xopts)
     !
     !--always label y axis
     !
     yopts = '1VN'//trim(yopts)
     call pgmtxt('L',ylabeloffset,0.5,0.5,labely)
     !
     !--always plot title
     !
     call pgmtxt('T',-titleoffset,0.5,0.5,title)

  endif

  call pgbox(xopts,0.0,0,yopts,0.0,0)

  return      
end subroutine

!
!--this subroutine can be called after PGSVP to
!  make sure that the viewport lies exactly on
!  pixel boundaries.
!
!  Queries PGPLOT routines directly so no need
!  for input/output
!

subroutine set_exactpixelboundaries()
 implicit none
 real :: xminpix,xmaxpix,yminpix,ymaxpix
 real :: vptxmin,vptxmax,vptymin,vptymax
 real :: dv
 real, parameter :: tol = 1.e-6
 !
 ! setting axes adjusts the viewport, so query to get adjusted settings
 !
 call pgqvp(0,vptxmin,vptxmax,vptymin,vptymax)
 !print*,'got ',vptxmin,vptxmax,vptymin,vptymax
 !
 ! adjust viewport on pixel devices so that 
 ! boundaries lie exactly on pixel boundaries
 !
 !  query viewport size in pixels
 call pgqvp(3,xminpix,xmaxpix,yminpix,ymaxpix)
 !print*,' in pixels = ',xminpix,xmaxpix,yminpix,ymaxpix

 !  work out how many viewport coords/pixel
 dv = (vptymax - vptymin)/(ymaxpix-yminpix)

 !  adjust viewport min/max to lie on pixel boundaries
 vptymin = max((nint(yminpix)-tol)*dv,0.)
 vptymax = min((nint(ymaxpix)-tol)*dv,1.0-epsilon(1.0)) ! be careful of round-off errors

 !  same for x
 dv = (vptxmax - vptxmin)/(xmaxpix-xminpix)
 vptxmin = max((nint(xminpix)-tol)*dv,0.)
 vptxmax = min((nint(xmaxpix)-tol)*dv,1.0-epsilon(1.0)) ! be careful of round-off errors

 !  adjust viewport
 !print*,'adjusting ',vptxmin,vptxmax,vptymin,vptymax
 call pgsvp(vptxmin,vptxmax,vptymin,vptymax)

 !call pgqvp(3,xminpix,xmaxpix,yminpix,ymaxpix)
 !print*,' in pixels = ',xminpix,xmaxpix,yminpix,ymaxpix

 return
end subroutine set_exactpixelboundaries

end module pagesetup
