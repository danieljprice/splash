!
! module containing some utilities for use with PGPLOT
!
module danpgutils
 implicit none
 public :: danpgsch, danpgtile
 
 private

contains
!
!  simple subroutine to set the pgplot character height
!  in a variety of units, independent of the page size
!
! arguments:
!  size   (input)  : character height in the appropriate units
!  units  (input)  : used to specify the units of the output value:
!                    units = 0 : normalized device coordinates
!                    units = 1 : inches
!                    units = 2 : millimeters
!                    units = 3 : pixels
!                    units = 4 : world coordinates
!                    other values give an error message, and are
!                    treated as 0.
!
!  Daniel Price, Institute of Astronomy, Cambridge, 2004.
!
  subroutine danpgsch(size,units)
   implicit none
   integer units
   real size
   real xch, ych, scalefac, charheight
!
! query current character height in the appropriate units
!
   call pgqch(charheight)
   call pgqcs(units,xch,ych)
!
! scale the character height appropriately to the desired value
!
   scalefac = size/ych
!
! reset character height
!
   call pgsch(scalefac*charheight)
   call pgqch(scalefac)

   return
  end subroutine danpgsch

!
!  simple subroutine to tile graphs appropriately on a page in pgplot
!  divides up a single panel into subpanels, with a margin at the edge
!  should replace the call to pgenv and pglabel
!
!  the page setup looks like this:
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
!   xmin,  : xmax, ymin, ymax : plot limits (should be same for all plots)
!   labelx : x axis label (should be same for all plots)
!   labely : y axis label (should be same for all plots)
!   title  : current plot title (can differ between plots)
!   just   : just=1 gives equal aspect ratios (same as in pgenv)
!   axis   : axes options (same as in pgenv)
!   vmarginleft,right,bottom,top : initial margin sizes (% of page)
!            (default should be zero for these)
!
!  Daniel Price, Institute of Astronomy, Cambridge, 2004.
!
  subroutine danpgtile(iplotin,nx,ny,xmin,xmax,ymin,ymax, &
                       labelx,labely,title,just,axis, &
                       vmarginleftin,vmarginrightin, &
                       vmarginbottomin,vmargintopin)
  implicit none
  integer iplotin,nx,ny,just,axis
  integer iplot,ix,iy
  real xmin,xmax,ymin,ymax,vmarginleftin,vmarginrightin
  real vmargintopin,vmarginbottomin
  real vptsizeeffx,vptsizeeffy,panelsizex,panelsizey
  real vmargintop,vmarginbottom,vmarginleft,vmarginright
  real vptxmin,vptxmax,vptymin,vptymax
  real aspectratio,devaspectratio,x1,x2,y1,y2
  real xch,ych,xlabeloffset,ylabeloffset
  character xopts*10, yopts*10
  character*(*) labelx,labely,title
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
! adjust effective viewport size if just=1 and graphs are not square
!      
  if (just.eq.1) then
     if (ymax.eq.ymin) then
        print*,'danpgtile: error: ymax=ymin'
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
  xlabeloffset = 3.0
  ylabeloffset = 4.5
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
  vmargintop = vmargintopin
  vmarginright = vmarginrightin
  if (axis.ge.0) then
     vmarginleft = vmarginleftin + (ylabeloffset+1.0)*xch
     vmarginbottom = vmarginbottomin + (xlabeloffset+1.0)*ych
  else
     vmarginleft = vmarginleftin
     vmarginbottom = vmarginbottomin
  endif
!
! effective viewport size = size - margins
!
  vptsizeeffx = 1.0 - vmarginright - vmarginleft
  vptsizeeffy = 1.0 - vmargintop - vmarginbottom
!     reduce x or y size if just=1 to get right aspect ratio
  if (aspectratio.lt.1.0) then
     vptsizeeffx = aspectratio*vptsizeeffy
  elseif (aspectratio.gt.1.0) then
     vptsizeeffy = vptsizeeffx/aspectratio
  endif
!
!--set size of each panel
!      
  panelsizex = vptsizeeffx/nx
  panelsizey = vptsizeeffy/ny 
  ix = iplot - ((iplot-1)/nx)*nx
  iy = (iplot-1)/nx + 1
!      print*,i,ix,iy
!      print*,panelsizex,panelsizey,vptsizeeffx,vptsizeeffy

  vptxmin = vmarginleft + (ix-1)*panelsizex
  vptxmax = vptxmin + panelsizex
  vptymax = (1.0 - vmargintop) - (iy-1)*panelsizey
  vptymin = vptymax - panelsizey
!      print*,vptxmin,vptxmax,vptymin,vptymax
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
! option to return before actually doing anything
!      
  if (title.eq.'nopgbox') return
!
! set options for call to pgbox (draws axes) and label axes where appropriate
! (options are exactly as in pgenv apart from axis=-3 which i have added)
!
  yopts = '*'
  if (axis.eq.-3) then
     xopts = 'bcst'
  elseif (axis.eq.-2) then
     xopts = ' '
  elseif (axis.eq.-1) then
    xopts = 'bc'
  elseif (axis.eq.0) then
    xopts = 'bcst'
  elseif (axis.eq.1) then
    xopts = 'abcst'
  elseif (axis.eq.2) then
    xopts = 'abcgst'
  elseif (axis.eq.10) then
    xopts = 'bcstl'
    yopts = 'bcst'
  elseif (axis.eq.20) then
    xopts = 'bcst'
    yopts = 'bcstl'
  elseif (axis.eq.30) then
    xopts = 'bcstl'
    yopts = 'bcstl'
  else
    call grwarn('danpgtile: illegal axis argument.')
    xopts = 'bcnst'
  endif
  if (yopts.eq.'*') yopts = xopts
!
! decide whether to number and label the y axis
!      
  if (ix.eq.1 .and. axis.ge.0) then
     yopts = '1vn'//yopts
     call pgmtxt('l',ylabeloffset,0.5,0.5,labely)
  elseif (axis.ge.0) then
     yopts = yopts//'n'
  endif  
!
! decide whether to number and label the x axis
!      
  if (iy.eq.ny .and. axis.ge.0) then
     xopts = 'n'//xopts
     call pgmtxt('b',xlabeloffset,0.5,0.5,labelx)
  endif

  call pgbox(xopts,0.0,0,yopts,0.0,0)
!
! plot the title inside the plot boundaries
!
  call pgmtxt('t',-1.5,0.96,1.0,title)

  return      
end subroutine


end module danpgutils
