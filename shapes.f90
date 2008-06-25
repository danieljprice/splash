!-----------------------------------------------------------------
! module which handles plotting of arbitrary shapes
! written by Daniel Price 2008
! as part of the SPLASH SPH visualisation package
!-----------------------------------------------------------------
module shapes
 implicit none
 integer, parameter, private :: maxshapes = 10
 integer, parameter, private :: maxshapetype = 6
 integer :: nshapes
 
 type shapedef
   integer :: itype
   integer :: icolour
   integer :: linestyle
   integer :: linewidth
   integer :: ifillstyle
   real :: xpos
   real :: ypos
   real :: xlen,ylen
   real :: angle,fjust
   character(len=20) :: text
 end type
 type(shapedef), dimension(maxshapes), public :: shape

 character(len=9), dimension(maxshapetype), &
 parameter, public :: labelshapetype = &
    (/'square   ', &
      'rectangle', &
      'arrow    ', &
      'circle   ', &
      'line     ', &
      'text     '/)
 namelist /shapeopts/ nshapes,shape
 
contains
!-----------------------------------------------------------------
! shape submenu
!-----------------------------------------------------------------
subroutine defaults_set_shapes
 implicit none
 
 nshapes = 0
 shape(:)%itype = 0
 shape(:)%icolour = 1
 shape(:)%linestyle = 1
 shape(:)%linewidth = 1
 shape(:)%ifillstyle = 2
 shape(:)%xpos = 0.5
 shape(:)%ypos = 0.5
 shape(:)%xlen = 1.
 shape(:)%ylen = 1.
 shape(:)%angle = 0.
 shape(:)%text = ' '
 shape(:)%fjust = 0.
 
 return
end subroutine defaults_set_shapes

!-----------------------------------------------------------------
! shape submenu
!-----------------------------------------------------------------
subroutine submenu_shapes()
 use prompting, only:prompt
 implicit none
 integer :: i,ishape,itype,indexi
 character(len=8) :: poslabel
 character(len=80) :: string

 itype = 1
 ishape = 0
 do while (itype.ne.0)
    indexi = 1
    do i=1,maxshapetype
       if (i.lt.10) then
          write(string(indexi:),"(1x,i1,') ',a)") i,trim(labelshapetype(i))
       else
          write(string(indexi:),"(1x,i2,') ',a)") i,trim(labelshapetype(i))
       endif
       indexi = len_trim(string) + 1
    enddo
    print "(/,a)",trim(string)
    !print "(i2,a)",(i,') '//trim(labelshapetype(i)),i=1,maxshapetype)
    ishape = ishape + 1
    call prompt('enter shape to plot (0 = finish) ',shape(ishape)%itype,0,maxshapetype)
    itype = shape(ishape)%itype
    if (itype.gt.0) then
       print "(a,i1,a)",'shape ',ishape,': type = '//trim(labelshapetype(itype))
       select case(itype)
       case(1) ! square
          call prompt('enter length of side (in x units of plot)',shape(ishape)%xlen,0.)
          shape(ishape)%ylen = shape(ishape)%xlen
          poslabel = 'centre'
       case(2) ! rectangle
          call prompt('enter length of side (in x units of plot)',shape(ishape)%xlen,0.)
          call prompt('enter length of side (in y units of plot)',shape(ishape)%ylen,0.)
          poslabel = 'centre'
       case(3) ! arrow
          call prompt('enter arrow length (in x units of plot)',shape(ishape)%xlen,0.)
          call prompt('enter angle in degrees (0 = horizontal) ',shape(ishape)%angle)
          poslabel = 'head'
       case(4) ! circle
          call prompt('enter radius ',shape(ishape)%xlen,0.)
          poslabel = 'centre'
       case(5) ! line
          call prompt('enter line length ',shape(ishape)%xlen,0.)
          call prompt('enter angle of line in degrees (0 = horizontal) ',shape(ishape)%angle)
          poslabel = 'starting'
       case(6) ! text
          call prompt('enter text string ',shape(ishape)%text)
          call prompt('enter angle for text in degrees (0 = horizontal) ',shape(ishape)%angle)
          call prompt('enter justification factor (0.0=left 1.0=right)',shape(ishape)%fjust)
          poslabel = 'starting'
       end select
       call prompt('enter '//trim(poslabel)//' x position (in x units of plot) ',shape(ishape)%xpos)
       call prompt('enter '//trim(poslabel)//' y position (in y units of plot) ',shape(ishape)%ypos)
       if (itype.eq.1 .or. itype.eq.2 .or. itype.eq.4) then
          call prompt('enter PGPLOT fill style (1=solid,2=outline,3=hatch,4=crosshatch) for '// &
                      trim(labelshapetype(itype)),shape(ishape)%ifillstyle,0,5)       
       endif
       if (itype.ne.6) then
          call prompt('enter PGPLOT line style (1=solid,2=dash,3=dotdash,4=dot,5=dashdot) for '// &
                      trim(labelshapetype(itype)),shape(ishape)%linestyle,0,5)
       endif
       if (itype.ne.6) then
          call prompt('enter PGPLOT line width for '//trim(labelshapetype(itype)),shape(ishape)%linewidth,0)
       endif
       call prompt('enter '//trim(labelshapetype(itype))//' colour (0=background, 1=foreground, 2-16=pgplot colour indices)', &
                   shape(ishape)%icolour,0,16)
    endif
 enddo
 nshapes = ishape - 1
 if (nshapes.gt.0) then
    print "(/,i2,a,/,15('-'),10(/,i2,')',1x,a10,' (x,y) = (',1pe10.2,',',1pe10.2,')'))",nshapes,' SHAPES SET: ', &
           (ishape,labelshapetype(shape(ishape)%itype),shape(ishape)%ypos,shape(ishape)%ypos,ishape=1,nshapes)
 else
    print "(a)",' NO SHAPES SET '
 endif
 
 return
end subroutine submenu_shapes

subroutine plot_shapes
 implicit none
 integer :: icolourprev,linestyleprev,linewidthprev,ifillstyle,i
 real :: xmin,xmax,ymin,ymax,dxplot,dyplot
 real :: xpos,ypos,xlen,ylen,angle,dx,dy
 real, dimension(2) :: xline,yline
 real, parameter :: pi = 3.1415926536
!
!--store current settings
!
 call pgqci(icolourprev)
 call pgqls(linestyleprev)
 call pgqlw(linewidthprev)
 call pgqfs(ifillstyle)
 !
!--convert hpos and vpos to x, y to plot arrow
!
 call pgqwin(xmin,xmax,ymin,ymax)
 dxplot = xmax - xmin
 dyplot = ymax - ymin
 
 do i=1,nshapes
    call pgsci(shape(i)%icolour)
    call pgsls(shape(i)%linestyle)
    call pgslw(shape(i)%linewidth)
    call pgsfs(shape(i)%ifillstyle)

    xpos = shape(i)%xpos
    ypos = shape(i)%ypos
    xlen = shape(i)%xlen
    ylen = shape(i)%ylen
    angle = shape(i)%angle*(pi/180.)
    
    print "(a)",'> plotting shape: '//trim(labelshapetype(shape(i)%itype))
    select case(shape(i)%itype)
    case(1,2) ! square, rectangle
       if (xlen.gt.dxplot .or. ylen.gt.dyplot) then
          print "(2x,a)",'Error: shape size exceeds plot dimensions: not plotted'
       else
          call pgrect(xpos-0.5*xlen,xpos+0.5*xlen,ypos-0.5*ylen,ypos + 0.5*ylen)   
       endif
    case(3) ! arrow
       dx = xlen*cos(angle)
       dy = xlen*sin(angle)
       !--do not plot if length > size of plot
       if (dx.gt.dxplot .or. dy.gt.dyplot) then
          print "(2x,a)",'Error: arrow length exceeds plot dimensions: arrow not plotted'
       else
          call pgarro(xpos-dx,ypos-dy,xpos,ypos)
       endif
    case(4) ! circle
       if (xlen.gt.dxplot .or. xlen.gt.dyplot) then
          print "(2x,a)",'Error: circle radius exceeds plot dimensions: circle not plotted'
       else
          call pgcirc(xpos,ypos,xlen)
       endif
    case(5) ! line
       xline(1) = xpos
       yline(1) = ypos
       xline(2) = xpos + xlen*cos(angle)
       yline(2) = ypos + xlen*sin(angle)
       call pgline(2,xline,yline)
    case(6) ! text
       call pgptext(xpos,ypos,angle,shape(i)%fjust,trim(shape(i)%text))
    end select
 enddo
 
 call pgsci(icolourprev)
 call pgsls(linestyleprev)
 call pgslw(linewidthprev)
 call pgsfs(ifillstyle)
 
end subroutine plot_shapes

end module shapes
