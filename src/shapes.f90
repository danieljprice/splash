!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2009 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! module which handles plotting of arbitrary shapes
! written by Daniel Price 2008
! as part of the SPLASH SPH visualisation package
!-----------------------------------------------------------------
module shapes
 implicit none
 integer, parameter, private :: maxshapes = 10
 integer, parameter, private :: maxshapetype = 7
 integer :: nshapes

 type shapedef
   integer :: itype
   integer :: icolour
   integer :: linestyle
   integer :: linewidth
   integer :: ifillstyle
   integer :: iunits
   integer :: iplotonpanel
   real :: xpos
   real :: ypos
   real :: xlen,ylen
   real :: angle,fjust
   character(len=120) :: text
 end type
 type(shapedef), dimension(maxshapes), public :: shape

 character(len=9), dimension(maxshapetype), &
 parameter, public :: labelshapetype = &
    (/'square   ', &
      'rectangle', &
      'arrow    ', &
      'circle   ', &
      'line     ', &
      'text     ', &
      'f(x)     '/)
 namelist /shapeopts/ nshapes,shape

 real, parameter, private :: pi = 3.1415926536

contains
!-----------------------------------------------------------------
! shape default settings
!-----------------------------------------------------------------
subroutine defaults_set_shapes
 implicit none

 nshapes = 0
 shape(:)%itype = 0
 shape(:)%icolour = 1
 shape(:)%linestyle = 1
 shape(:)%linewidth = 1
 shape(:)%ifillstyle = 2
 shape(:)%iunits = 1
 shape(:)%iplotonpanel = 0
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
 use params,        only:maxplot
 use prompting,     only:prompt
 use exactfunction, only:check_function
 implicit none
 integer            :: i,ishape,itype,indexi,iunits,ierr,itry
 character(len=8)   :: poslabel
 character(len=80)  :: string
 integer, parameter :: maxunits = 2
 character(len=20), dimension(maxunits), &
    parameter :: labelunits = &
    (/'units of plot       ', &
      'viewport coordinates'/)
!      'inches             ', &
!      'millimeters        ', &
!      'pixels             '/)

 itype = 1
 ishape = 0
 over_shapes: do while (itype.ne.0)
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
       print "(a,i1,a,/)",'shape ',ishape,': type = '//trim(labelshapetype(itype))

       if (itype.eq.7) then
          shape(ishape)%iunits = 1
       else
          !--choose units
          do i=1,maxunits
             print "(i1,')',1x,a)",i,trim(labelunits(i))
          enddo
          call prompt('enter units to use for plotting shape',shape(ishape)%iunits,0,maxunits)
       endif
       iunits = shape(ishape)%iunits

       select case(itype)
       case(1) ! square
          call prompt('enter length of side (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen,0.)
          shape(ishape)%ylen = shape(ishape)%xlen
          poslabel = 'centre'
       case(2) ! rectangle
          call prompt('enter x length of side (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen,0.)
          call prompt('enter y length of side (in '//trim(labelunits(iunits))//')',shape(ishape)%ylen,0.)
          poslabel = 'centre'
       case(3) ! arrow
          call prompt('enter arrow length (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen,0.)
          call prompt('enter angle in degrees (0 = horizontal) ',shape(ishape)%angle)
          poslabel = 'head'
       case(4) ! circle
          call prompt('enter radius (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen,0.)
          poslabel = 'centre'
       case(5) ! line
          call prompt('enter line length (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen,0.)
          call prompt('enter angle of line in degrees (0.0 = horizontal) ',shape(ishape)%angle)
          poslabel = 'starting'
       case(6) ! text
          call prompt('enter text string ',shape(ishape)%text)
          call prompt('enter angle for text in degrees (0 = horizontal) ',shape(ishape)%angle)
          call prompt('enter justification factor (0.0=left 1.0=right)',shape(ishape)%fjust)
          poslabel = 'starting'
       case(7) ! arbitrary function
          ierr = 1
          itry = 1
          do while(ierr /= 0 .and. itry.le.3)
             if (itry.gt.1) print "(a,i1,a)",'attempt ',itry,' of 3:'
             print "(a,6(/,11x,a),/)",' Examples: sin(2*pi*x)','sqrt(0.5*x)','x^2', &
             'exp(-2*x**2)','log10(x/2)','exp(p),p=sin(pi*x)','cos(z/d),z=acos(d),d=x^2'
             call prompt('enter function f(x) to plot ',shape(ishape)%text)
             call check_function(shape(ishape)%text,ierr)
             if (ierr /= 0 .and. len(shape(ishape)%text).eq.len_trim(shape(ishape)%text)) then
                print "(a,i3,a)",' (errors are probably because string is too long, max length = ',len(shape(ishape)%text),')'
             endif
             itry = itry + 1
          enddo
          if (ierr.ne.0) then
             print "(a)",' *** too many tries, aborting ***'
             ishape = ishape - 1
             cycle over_shapes
          endif
          poslabel = 'starting'
       end select
       if (itype.ne.7) then
          call prompt('enter '//trim(poslabel)//' x position (in '//trim(labelunits(iunits))//') ',shape(ishape)%xpos)
          call prompt('enter '//trim(poslabel)//' y position (in '//trim(labelunits(iunits))//') ',shape(ishape)%ypos)
       endif
       if (itype.eq.1 .or. itype.eq.2 .or. itype.eq.4) then
          call prompt('enter fill style (1=solid,2=outline,3=hatch,4=crosshatch) for '// &
                      trim(labelshapetype(itype)),shape(ishape)%ifillstyle,0,5)
       endif
       if (itype.ne.6) then
          call prompt('enter line style (1=solid,2=dash,3=dotdash,4=dot,5=dashdot) for '// &
                      trim(labelshapetype(itype)),shape(ishape)%linestyle,0,5)
       endif
       if (itype.ne.6) then
          call prompt('enter line width for '//trim(labelshapetype(itype)),shape(ishape)%linewidth,0)
       endif
       call prompt('enter '//trim(labelshapetype(itype))//' colour (0=background, 1=foreground, 2-16=pgplot colour indices)', &
                   shape(ishape)%icolour,0,16)

       print "(/,'  0 : plot on every panel ',/,"// &
              "' -1 : plot on first row only ',/,"// &
              "' -2 : plot on first column only ',/,"// &
              "'  n : plot on nth panel only ')"

       !--make sure the current setting falls within the allowed bounds
       if (shape(ishape)%iplotonpanel.lt.-2 .or. &
           shape(ishape)%iplotonpanel.gt.maxplot) shape(:)%iplotonpanel = 0

       call prompt('Enter selection ',shape(ishape)%iplotonpanel,-2,maxplot)
    endif
 enddo over_shapes
 nshapes = ishape - 1
 if (nshapes.gt.0) then
    print "(/,a,/,15('-'),10(/,i2,')',1x,a10,' (x,y) = (',1pe10.2,',',1pe10.2,') [',a,']'))",' SHAPES SET: ', &
           (ishape,labelshapetype(shape(ishape)%itype),shape(ishape)%ypos,shape(ishape)%ypos, &
           trim(labelunits(shape(ishape)%iunits)),ishape=1,nshapes)
 else
    print "(a)",' NO SHAPES SET '
 endif

 return
end subroutine submenu_shapes

subroutine plot_shapes(ipanel,irow,icolumn,itransx,itransy)
 use exactfunction, only:exact_function
 use transforms,    only:transform_inverse,transform
 use plotlib, only:plot_qci,plot_qls,plot_qlw,plot_qfs,plot_qwin,plot_sci,plot_sfs,plot_slw, &
      plot_sci,plot_rect,plot_sls,plot_line,plot_arro,plot_circ,plot_ptxt
 implicit none
 integer, intent(in) :: ipanel,irow,icolumn,itransx,itransy
 integer :: icolourprev,linestyleprev,linewidthprev,ifillstyle
 integer :: i,j,ierr,iplotonthispanel
 integer, parameter :: maxfuncpts = 1000
 real :: xmin,xmax,ymin,ymax,dxplot,dyplot
 real :: xpos,ypos,xlen,ylen,anglerad,dx,dy
 real, dimension(2) :: xline,yline
 real, dimension(maxfuncpts) :: xfunc,yfunc
!
!--store current settings
!
 call plot_qci(icolourprev)
 call plot_qls(linestyleprev)
 call plot_qlw(linewidthprev)
 call plot_qfs(ifillstyle)
 !
!--convert hpos and vpos to x, y to plot arrow
!
 call plot_qwin(xmin,xmax,ymin,ymax)
 dxplot = xmax - xmin
 dyplot = ymax - ymin
!
!--query window size in a variety of other units
!
 do i=1,nshapes

    iplotonthispanel = shape(i)%iplotonpanel
    if (iplotonthispanel.eq.0 &
       .or.(iplotonthispanel.gt.0  .and. ipanel.eq.iplotonthispanel) &
       .or.(iplotonthispanel.eq.-1 .and. irow.eq.1) &
       .or.(iplotonthispanel.eq.-2 .and. icolumn.eq.1)) then

       call plot_sci(shape(i)%icolour)
       call plot_sls(shape(i)%linestyle)
       call plot_slw(shape(i)%linewidth)
       call plot_sfs(shape(i)%ifillstyle)

       anglerad = shape(i)%angle*(pi/180.)

       call convert_units(shape(i),xpos,ypos,xlen,ylen, &
                          xmin,ymin,dxplot,dyplot,itransx,itransy)

       print "(a)",'> plotting shape: '//trim(labelshapetype(shape(i)%itype))
       select case(shape(i)%itype)
       case(1,2) ! square, rectangle
          if (xlen.gt.dxplot .or. ylen.gt.dyplot) then
             print "(2x,a)",'Error: shape size exceeds plot dimensions: not plotted'
          else
             call plot_rect(xpos-0.5*xlen,xpos+0.5*xlen,ypos-0.5*ylen,ypos + 0.5*ylen)
          endif
       case(3) ! arrow
          dx = xlen*cos(anglerad)
          dy = xlen*sin(anglerad)
          !--do not plot if length > size of plot
          if (dx.gt.dxplot .or. dy.gt.dyplot) then
             print "(2x,a)",'Error: arrow length exceeds plot dimensions: arrow not plotted'
          else
             call plot_arro(xpos-dx,ypos-dy,xpos,ypos)
          endif
       case(4) ! circle
          if (xlen.gt.dxplot .or. xlen.gt.dyplot) then
             print "(2x,a)",'Error: circle radius exceeds plot dimensions: circle not plotted'
          else
             call plot_circ(xpos,ypos,xlen)
          endif
       case(5) ! line
          xline(1) = xpos
          yline(1) = ypos
          xline(2) = xpos + xlen*cos(anglerad)
          yline(2) = ypos + xlen*sin(anglerad)
          call plot_line(2,xline,yline)
       case(6) ! text
          call plot_ptxt(xpos,ypos,shape(i)%angle,shape(i)%fjust,trim(shape(i)%text))
       case(7) ! arbitrary function
          !--set x to be evenly spaced in transformed (plot) coordinates
          dx = (xmax-xmin)/real(maxfuncpts - 1)
          do j=1,maxfuncpts
             xfunc(j) = xmin + (j-1)*dx
          enddo
          !--transform x array back to untransformed space to evaluate f(x)
          if (itransx.gt.0) call transform_inverse(xfunc,itransx)
          call exact_function(shape(i)%text,xfunc,yfunc,0.,ierr)
          if (ierr.eq.0) then
             !--reset x values
             do j=1,maxfuncpts
                xfunc(j) = xmin + (j-1)*dx
             enddo
             !--transform y if necessary
             if (itransy.gt.0) call transform(yfunc,itransy)
             !--plot the line
             call plot_line(maxfuncpts,xfunc,yfunc)
          endif
       end select
    endif
 enddo

 call plot_sci(icolourprev)
 call plot_sls(linestyleprev)
 call plot_slw(linewidthprev)
 call plot_sfs(ifillstyle)

end subroutine plot_shapes

integer function inshape(xpt,ypt,itransx,itransy)
 use plotlib, only:plot_qwin,plot_qtxt
 implicit none
 real, intent(in) :: xpt,ypt
 integer, intent(in) :: itransx,itransy
 integer :: i
 real :: xpos,ypos,xlen,ylen
 real :: xmin,ymin,xmax,ymax,dxplot,dyplot
 real, dimension(4) :: xbox,ybox

 call plot_qwin(xmin,xmax,ymin,ymax)
 dxplot = xmax - xmin
 dyplot = ymax - ymin

 inshape = 0
 do i=1,nshapes

    call convert_units(shape(i),xpos,ypos,xlen,ylen, &
                       xmin,ymin,dxplot,dyplot,itransx,itransy)

    select case(shape(i)%itype)
    case(1,2)  ! square, rectangle

    case(3) ! arrow

    case(4) ! circle

    case(5) ! line

    case(6) ! text
       call plot_qtxt(xpos,ypos,shape(i)%angle,shape(i)%fjust,trim(shape(i)%text),xbox,ybox)
       if (xpt.gt.minval(xbox) .and. xpt.le.maxval(xbox) &
          .and. ypt.gt.minval(ybox) .and. ypt.le.maxval(ybox)) then
          inshape = i
       endif
    end select
 enddo

end function inshape

subroutine edit_shape(i,xpt,ypt,itransx,itransy)
 use plotlib, only:plot_qwin
 implicit none
 integer, intent(in) :: i,itransx,itransy
 real, intent(in)    :: xpt,ypt
 real :: xmin,xmax,ymin,ymax,dxplot,dyplot,xlen,ylen
 real :: xpos,ypos

 call plot_qwin(xmin,xmax,ymin,ymax)
 dxplot = xmax - xmin
 dyplot = ymax - ymin

 call convert_units(shape(i),xpos,ypos,xlen,ylen, &
                    xmin,ymin,dxplot,dyplot,itransx,itransy)

 select case(shape(i)%itype)
 case(6)
    call edit_textbox(xpos,ypos,shape(i)%angle,shape(i)%text)
 case default

 end select

end subroutine edit_shape

subroutine delete_shape(ishape)
 implicit none
 integer, intent(in) :: ishape
 integer :: i

 if (ishape.gt.0) then
    do i=ishape+1,nshapes
       shape(i-1) = shape(i)
    enddo
    print "(a)",'> deleted shape: '//trim(labelshapetype(shape(ishape)%itype))
    nshapes = nshapes - 1
 endif

end subroutine delete_shape

subroutine add_textshape(xpt,ypt,itransx,itransy,ipanel)
 use plotlib, only:plot_qwin
 implicit none
 real, intent(in)    :: xpt,ypt
 integer, intent(in) :: itransx,itransy,ipanel
 integer :: i
 real :: xmin,xmax,ymin,ymax,xposi,yposi

 nshapes = nshapes + 1
 if (nshapes.gt.maxshapes) then
    print*,' *** cannot add shape: array limits reached, delete some shapes first ***'
    nshapes = maxshapes
    return
 endif
 i = nshapes
 shape(i)%itype = 6
 print*,' adding shape '//trim(labelshapetype(shape(i)%itype))
 shape(i)%icolour = 1
 shape(i)%linestyle = 1
 shape(i)%linewidth = 1
 shape(i)%ifillstyle = 2
 shape(i)%iplotonpanel = ipanel
!
!--position text relative to viewport
!
 shape(i)%iunits = 2
 call plot_qwin(xmin,xmax,ymin,ymax)
 xposi = (xpt - xmin)/(xmax-xmin)
 yposi = (ypt - ymin)/(ymax-ymin)
 shape(i)%xpos = xposi
 shape(i)%ypos = yposi

 shape(i)%xlen = 1.
 shape(i)%ylen = 1.
 shape(i)%angle = 0.
 shape(i)%text = 'click to edit'
 shape(i)%fjust = 0.
 call edit_shape(i,xposi,yposi,itransx,itransy)

end subroutine add_textshape


subroutine convert_units(shape,xpos,ypos,xlen,ylen,xmin,ymin,dxplot,dyplot,itransx,itransy)
 use transforms, only:transform
 implicit none
 type(shapedef), intent(in) :: shape
 real, intent(out)   :: xpos,ypos,xlen,ylen
 real, intent(in)    :: xmin,ymin,dxplot,dyplot
 integer, intent(in) :: itransx,itransy

 xpos = shape%xpos
 ypos = shape%ypos
 xlen = shape%xlen
 ylen = shape%ylen

 select case(shape%iunits)
 case(2) ! translate from viewport coordinates into plot coordinates
    xpos = xmin + xpos*dxplot
    ypos = ymin + ypos*dyplot
    xlen = xlen*dxplot
    ylen = ylen*dyplot
 case(1)
    if (itransx.gt.0) then
       call transform(xpos,itransx)
       call transform(xlen,itransx)
    endif
    if (itransy.gt.0) then
       call transform(ypos,itransy)
       call transform(ylen,itransy)
    endif
    ! do nothing here
 case default ! should never happen
    print "(a)",' INTERNAL ERROR: unknown units whilst plotting shape'
 end select

end subroutine convert_units


subroutine edit_textbox(xpt,ypt,angle,string)
use plotlib, only:plot_stbg,plot_ptxt,plot_curs
 implicit none
 real, intent(in) :: xpt,ypt,angle
 character(len=1) :: mychar
 real :: xpt2,ypt2
 character(len=*), intent(inout) :: string
 character(len=len(string)) :: oldstring
 integer :: i,ierr

 print*,'editing text box, esc or ctrl-c to quit'
 call plot_stbg(0)
 mychar = ' '
 oldstring = string
 i = max(len_trim(string)+1,1)
 call plot_ptxt(xpt,ypt,angle,0.,string(1:i)//'_')

 xpt2 = xpt
 ypt2 = ypt
 ierr = plot_curs(xpt2,ypt2,mychar)
 do while (mychar.ne.achar(13) &   ! carriage return
     .and. mychar.ne.achar(27) &   ! ctrl-c
     .and. mychar.ne.achar(3))     ! esc
    if (mychar.eq.achar(8)) then   ! backspace
       i = max(i - 1,1)
       string(i:i) = '_'
       call plot_ptxt(xpt,ypt,angle,0.,string(1:i))
       string(i:i) = ' '
    else
       if (trim(string).eq.'click to edit') then
          !print*,'erasing string'
          string = ' '
          i = 1
       endif
       string(i:i) = mychar
       call plot_ptxt(xpt,ypt,angle,0.,string(1:i))
       i = min(i + 1,len(string))
       if (i.eq.len(string)) print*,' reached end of string'
    endif
    ierr = plot_curs(xpt2,ypt2,mychar)
 enddo

 !--if ctrl-c or esc, restore original string
 if (mychar.eq.achar(3) .or. mychar.eq.achar(27)) then
    string = oldstring
    print*,'cancelled'
 else
    print*,'done: text = "'//trim(string)//'"'
 endif
 call plot_stbg(-1)

end subroutine edit_textbox


end module shapes
