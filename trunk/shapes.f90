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
   integer :: iunits
   integer :: iplotonpanel
   real :: xpos
   real :: ypos
   real :: xlen,ylen
   real :: angle,fjust
   character(len=40) :: text
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
 use params, only:maxplot
 use prompting, only:prompt
 implicit none
 integer :: i,ishape,itype,indexi,iunits
 character(len=8) :: poslabel
 character(len=80) :: string
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
       print "(a,i1,a,/)",'shape ',ishape,': type = '//trim(labelshapetype(itype))
       
       !--choose units
       do i=1,maxunits
          print "(i1,')',1x,a)",i,trim(labelunits(i))
       enddo
       call prompt('enter units to use for plotting shape',shape(ishape)%iunits,0,maxunits)
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
       end select
       call prompt('enter '//trim(poslabel)//' x position (in '//trim(labelunits(iunits))//') ',shape(ishape)%xpos)
       call prompt('enter '//trim(poslabel)//' y position (in '//trim(labelunits(iunits))//') ',shape(ishape)%ypos)
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

       print "(/,'  0 : plot on every panel ',/,"// &
              "' -1 : plot on first row only ',/,"// &
              "' -2 : plot on first column only ',/,"// &
              "'  n : plot on nth panel only ')"

       !--make sure the current setting falls within the allowed bounds
       if (shape(ishape)%iplotonpanel.lt.-2 .or. &
           shape(ishape)%iplotonpanel.gt.maxplot) shape(:)%iplotonpanel = 0

       call prompt('Enter selection ',shape(ishape)%iplotonpanel,-2,maxplot)
    endif
 enddo
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
 implicit none
 integer, intent(in) :: ipanel,irow,icolumn,itransx,itransy
 integer :: icolourprev,linestyleprev,linewidthprev,ifillstyle
 integer :: i,iplotonthispanel
 real :: xmin,xmax,ymin,ymax,dxplot,dyplot
 real :: xpos,ypos,xlen,ylen,anglerad,dx,dy
 real, dimension(2) :: xline,yline
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
!
!--query window size in a variety of other units
!
 do i=1,nshapes

    iplotonthispanel = shape(i)%iplotonpanel
    if (iplotonthispanel.eq.0 &
       .or.(iplotonthispanel.gt.0  .and. ipanel.eq.iplotonthispanel) &
       .or.(iplotonthispanel.eq.-1 .and. irow.eq.1) &
       .or.(iplotonthispanel.eq.-2 .and. icolumn.eq.1)) then

       call pgsci(shape(i)%icolour)
       call pgsls(shape(i)%linestyle)
       call pgslw(shape(i)%linewidth)
       call pgsfs(shape(i)%ifillstyle)

       anglerad = shape(i)%angle*(pi/180.)

       call convert_units(shape(i),xpos,ypos,xlen,ylen, &
                          xmin,ymin,dxplot,dyplot,itransx,itransy)

       print "(a)",'> plotting shape: '//trim(labelshapetype(shape(i)%itype))
       select case(shape(i)%itype)
       case(1,2) ! square, rectangle
          if (xlen.gt.dxplot .or. ylen.gt.dyplot) then
             print "(2x,a)",'Error: shape size exceeds plot dimensions: not plotted'
          else
             call pgrect(xpos-0.5*xlen,xpos+0.5*xlen,ypos-0.5*ylen,ypos + 0.5*ylen)   
          endif
       case(3) ! arrow
          dx = xlen*cos(anglerad)
          dy = xlen*sin(anglerad)
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
          xline(2) = xpos + xlen*cos(anglerad)
          yline(2) = ypos + xlen*sin(anglerad)
          call pgline(2,xline,yline)
       case(6) ! text
          call pgptext(xpos,ypos,shape(i)%angle,shape(i)%fjust,trim(shape(i)%text))
       end select
    endif
 enddo
 
 call pgsci(icolourprev)
 call pgsls(linestyleprev)
 call pgslw(linewidthprev)
 call pgsfs(ifillstyle)
 
end subroutine plot_shapes

integer function inshape(xpt,ypt,itransx,itransy)
 implicit none
 real, intent(in) :: xpt,ypt
 integer, intent(in) :: itransx,itransy
 integer :: i
 real :: xpos,ypos,xlen,ylen
 real :: xmin,ymin,xmax,ymax,dxplot,dyplot
 real, dimension(4) :: xbox,ybox

 call pgqwin(xmin,xmax,ymin,ymax)
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
       call pgqtxt(xpos,ypos,shape(i)%angle,shape(i)%fjust,trim(shape(i)%text),xbox,ybox)
       if (xpt.gt.minval(xbox) .and. xpt.le.maxval(xbox) &
          .and. ypt.gt.minval(ybox) .and. ypt.le.maxval(ybox)) then
          inshape = i
       endif
    end select
 enddo
 
end function inshape

subroutine edit_shape(i,xpt,ypt,itransx,itransy)
 implicit none
 integer, intent(in) :: i,itransx,itransy
 real, intent(in)    :: xpt,ypt
 real :: xmin,xmax,ymin,ymax,dxplot,dyplot,xlen,ylen
 real :: xpos,ypos

 call pgqwin(xmin,xmax,ymin,ymax)
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
 implicit none
 real, intent(in) :: xpt,ypt,angle
 character(len=1) :: mychar
 real :: xpt2,ypt2
 character(len=*), intent(inout) :: string
 character(len=len(string)) :: oldstring
 integer :: i
 
 print*,'editing text box, esc or ctrl-c to quit'
 call pgstbg(0)
 mychar = ' '
 oldstring = string
 i = max(len_trim(string)+1,1)
 call pgptxt(xpt,ypt,angle,0.,string(1:i)//'_')
 xpt2 = xpt
 ypt2 = ypt
 call pgcurs(xpt2,ypt2,mychar)
 do while (mychar.ne.achar(13) &   ! carriage return
     .and. mychar.ne.achar(27) &   ! ctrl-c
     .and. mychar.ne.achar(3))     ! esc
    if (mychar.eq.achar(8)) then   ! backspace
       i = max(i - 1,1)
       string(i:i) = '_'
       call pgptxt(xpt,ypt,angle,0.,string(1:i))
       string(i:i) = ' '
    else
       string(i:i) = mychar
       call pgptxt(xpt,ypt,angle,0.,string(1:i))
       i = min(i + 1,len(string))
       if (i.eq.len(string)) print*,' reached end of string'
    endif
    call pgcurs(xpt2,ypt2,mychar)
 enddo
 
 !--if ctrl-c or esc, restore original string
 if (mychar.eq.achar(3) .or. mychar.eq.achar(27)) then
    string = oldstring
    print*,'cancelled'
 else
    print*,'done: text = "'//trim(string)//'"'
 endif
 call pgstbg(-1)
 
end subroutine edit_textbox


end module shapes
