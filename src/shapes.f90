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
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! module which handles plotting of arbitrary shapes
! written by Daniel Price 2008
! as part of the SPLASH SPH visualisation package
!-----------------------------------------------------------------
module shapes
 implicit none
 integer, parameter, private :: maxshapes = 32
 integer, parameter, private :: maxshapetype = 7
 integer :: nshapes
 integer, parameter, private :: lentext = 120

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
   real :: opacity
   character(len=lentext) :: text
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

 integer, parameter, private :: maxunits = 2
 character(len=20), dimension(maxunits), &
    parameter, private :: labelunits = &
    (/'units of plot       ', &
      'viewport coordinates'/)
!      'inches             ', &
!      'millimeters        ', &
!      'pixels             '/)

 procedure(check_shapes), pointer, private :: checkshapes => null()
 procedure(add_shape), pointer, private :: addshape => null()
 procedure(delete_shape), pointer, private :: delshape => null()

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
 shape(:)%opacity = 1.

 return
end subroutine defaults_set_shapes

!-----------------------------------------------------------------
! shape submenu
!-----------------------------------------------------------------
subroutine submenu_shapes()
 use promptlist, only:prompt_list
 implicit none
 
 checkshapes => check_shapes
 addshape => add_shape
 delshape => delete_shape
 call prompt_list(nshapes,maxshapes,'shape',checkshapes,addshape,delshape)
 
end subroutine submenu_shapes

!-----------------------------------
! print the current list of shapes
!-----------------------------------
subroutine check_shapes(nshape)
 implicit none
 integer, intent(in) :: nshape
 integer :: ishape

 print "(/,a)", ' Current list of plot annotations:'
 if (nshape.gt.0) then
    do ishape=1,nshape
       call print_shapeinfo(ishape,shape(ishape)%itype,shape(ishape))
    enddo
 else
    print "(a)",' (none)'
 endif

end subroutine check_shapes

!----------------------------------------
! pretty-print information about a shape
!----------------------------------------
subroutine print_shapeinfo(inum,itype,shapein)
 implicit none
 integer, intent(in) :: inum,itype
 type(shapedef), intent(in), optional :: shapein
 character(len=20), parameter :: fmtstring = "('Shape ',i2,': ',a)"

 select case(itype)
 case(1)
    print "(10x,a)",'    --'
    write(*,fmtstring,advance='no') inum,'   |  |   '//labelshapetype(itype)
    if (present(shapein)) then
       print "(5x,es10.2,' x ',es10.2)",shapein%xlen,shapein%ylen
    else
       print*
    endif
    print "(10x,a)",'    --'
 case(2)
    print "(10x,a)",'  -----'
    write(*,fmtstring,advance='no') inum,' |     |  '//labelshapetype(itype)
    if (present(shapein)) then
       print "(5x,es10.2,' x ',es10.2)",shapein%xlen,shapein%ylen
    else
       print*
    endif
    print "(10x,a)",'  -----'
 case(3)
    print "(10x,a)"
    write(*,fmtstring,advance='no') inum,' ------>  '//labelshapetype(itype)
    if (present(shapein)) then
       print "(6x,'length = ',es10.2,', angle = ',f5.1,' deg.')",&
             shapein%xlen,shapein%angle
    else
       print*
    endif
    print "(10x,a)"
 case(4)
    print "(10x,a)",'   ___ '
    print "(10x,a)",'  /   \ '
    write(*,fmtstring,advance='no') inum,' (     )  '//labelshapetype(itype)
    if (present(shapein)) then
       print "(6x,'radius = ',es10.2)",&
             shapein%xlen
    else
       print*
    endif
    print "(10x,a)",'  \___/ '
 case(5)
    print "(10x,a)"
    write(*,fmtstring,advance='no') inum,' -------- '//labelshapetype(itype)
    if (present(shapein)) then
       print "(6x,'length = ',es10.2)",shapein%xlen
    else
       print*
    endif
    print "(10x,a)"
 case(6)
    print "(10x,a)",                  '        '
    write(*,fmtstring,advance='no') inum,'   TEXT   '
    if (present(shapein)) then
       print "('""',a,'""')",trim(shapein%text)
    else
       print*
    endif
    print "(10x,a)",                  '        '
 case(7)
    print "(10x,a)",'  _     / '
    write(*,fmtstring,advance='no') inum,' / \   /  '//trim(labelshapetype(itype))
    if (present(shapein)) then
       print "(' = ',a)",trim(shapein%text)
    else
       print*
    endif
    print "(10x,a)",'/   \_/     '
 case default
    print "(a)"
    write(*,fmtstring,advance='no') inum,'   '//labelshapetype(itype)
    print "(a)"
 end select

end subroutine print_shapeinfo

!------------------------------------------
! utility routine to add new shape object
!------------------------------------------
subroutine add_shape(istart,iend,nshape)
 use params,        only:maxplot
 use prompting,     only:prompt
 use exactfunction, only:check_function
 use plotlib,       only:plotlib_maxlinestyle,plotlib_maxlinecolour,plotlib_maxfillstyle,plotlib_supports_alpha
 implicit none
 integer, intent(in) :: istart,iend
 integer, intent(inout) :: nshape
 integer            :: i,ishape,itype,indexi,iunits,ierr,itry
 character(len=10)  :: poslabel
 character(len=80)  :: string

 itype   = 1
 ishape  = istart + 1
 if (ishape.gt.maxshapes) then
    print "(/,a,i2,a)",' *** Error, maximum number of shapes (',maxshapes,') reached, cannot add any more.'
    print "(a)",       ' *** If you hit this limit, *please email me* so I can change the default limits!'
    print "(a)",       ' *** (and then edit shapes.f90, changing the parameter "maxshapes" to something higher...)'
    return
 endif
 !
 !--fill prompt string with list of shapes
 !
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
 
 over_shapes: do while(ishape.le.iend .and. i.le.maxshapes)
    if (istart.eq.0 .or. shape(ishape)%itype.le.0 .or. shape(ishape)%itype.gt.maxshapetype) then
       call prompt('choose an object type (0=none) ',shape(ishape)%itype,0,maxshapetype)
    endif
    itype = shape(ishape)%itype
    if (itype.eq.0) then
       call delete_shape(ishape,nshape)
       exit over_shapes
    else
       call print_shapeinfo(ishape,itype)

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
          poslabel = ' centre'
       case(2) ! rectangle
          call prompt('enter x length of side (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen,0.)
          call prompt('enter y length of side (in '//trim(labelunits(iunits))//')',shape(ishape)%ylen,0.)
          poslabel = ' centre'
       case(3) ! arrow
          call prompt('enter arrow length (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen,0.)
          call prompt('enter angle in degrees (0 = horizontal) ',shape(ishape)%angle)
          call prompt('enter justification factor (0.0=tail at x,y 1.0=head at x,y)',shape(ishape)%fjust)
          poslabel = ''
       case(4) ! circle
          call prompt('enter radius (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen,0.)
          poslabel = ' centre'
       case(5) ! line
          call prompt('enter line length (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen,0.)
          call prompt('enter angle of line in degrees (0.0 = horizontal) ',shape(ishape)%angle)
          poslabel = ' starting'
       case(6) ! text
          call prompt('enter text string ',shape(ishape)%text)
          call prompt('enter angle for text in degrees (0 = horizontal) ',shape(ishape)%angle)
          call prompt('enter justification factor (0.0=left 1.0=right)',shape(ishape)%fjust)
          poslabel = ' starting'
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
          poslabel = ' starting'
       end select
       if (itype.ne.7) then
          call prompt('enter'//trim(poslabel)//' x position (in '//trim(labelunits(iunits))//') ',shape(ishape)%xpos)
          call prompt('enter'//trim(poslabel)//' y position (in '//trim(labelunits(iunits))//') ',shape(ishape)%ypos)
       endif
       if (itype.eq.1 .or. itype.eq.2 .or. itype.eq.4) then
          call prompt('enter fill style (1=solid,2=outline,3=hatch,4=crosshatch) for '// &
                      trim(labelshapetype(itype)),shape(ishape)%ifillstyle,0,plotlib_maxfillstyle)
       endif
       if (itype.ne.6) then
          call prompt('enter line style (1=solid,2=dash,3=dotdash,4=dot,5=dashdot) for '// &
                      trim(labelshapetype(itype)),shape(ishape)%linestyle,0,plotlib_maxlinestyle)
       endif
       if (itype.ne.6) then
          call prompt('enter line width for '//trim(labelshapetype(itype)),shape(ishape)%linewidth,0)
       endif
       call prompt('enter '//trim(labelshapetype(itype))//' colour (0=background, 1=foreground, 2-16=plot lib colour indices)', &
                   shape(ishape)%icolour,0,plotlib_maxlinecolour)
       if (plotlib_supports_alpha) then
          call prompt('enter '//trim(labelshapetype(itype))//' opacity (0.0-1.0)', &
                   shape(ishape)%opacity,0.,1.)
       endif

       print "(/,'  0 : plot on every panel ',/,"// &
              "' -1 : plot on first row only ',/,"// &
              "' -2 : plot on first column only ',/,"// &
              "'  n : plot on nth panel only ')"

       !--make sure the current setting falls within the allowed bounds
       if (shape(ishape)%iplotonpanel.lt.-2 .or. &
           shape(ishape)%iplotonpanel.gt.maxplot) shape(:)%iplotonpanel = 0

       call prompt('Enter selection ',shape(ishape)%iplotonpanel,-2,maxplot)
       if (ishape.gt.nshape) nshape = ishape
       ishape = ishape + 1
    endif
 enddo over_shapes

end subroutine add_shape

!------------------------------------------
! utility routine to delete a shape object
!------------------------------------------
subroutine delete_shape(ishape,nshape)
 implicit none
 integer, intent(in)    :: ishape
 integer, intent(inout) :: nshape
 integer :: i

 if (ishape.gt.0 .and. nshape.gt.0 .and. ishape.le.maxshapes) then
    do i=ishape+1,nshape
       shape(i-1) = shape(i)
    enddo
    print "(a)",'> deleted shape: '//trim(labelshapetype(shape(ishape)%itype))
    !--restore defaults
    shape(nshape)%itype = 0
    shape(nshape)%icolour = 1
    shape(nshape)%linestyle = 1
    shape(nshape)%linewidth = 1
    shape(nshape)%ifillstyle = 2
    shape(nshape)%iunits = 1
    shape(nshape)%iplotonpanel = 0
    shape(nshape)%xpos = 0.5
    shape(nshape)%ypos = 0.5
    shape(nshape)%xlen = 1.
    shape(nshape)%ylen = 1.
    shape(nshape)%angle = 0.
    shape(nshape)%text = ' '
    shape(nshape)%fjust = 0.
    nshape = nshape - 1
 endif

end subroutine delete_shape

!------------------------------------------------------------
! actual routine that implements plotting of various shapes
!------------------------------------------------------------
subroutine plot_shapes(ipanel,irow,icolumn,itransx,itransy,time)
 use exactfunction, only:exact_function
 use transforms,    only:transform_inverse,transform
 use asciiutils,    only:string_replace
 use plotlib, only:plot_qci,plot_qls,plot_qlw,plot_qfs,plot_qwin,plot_sci,plot_sfs,plot_slw, &
      plot_sci,plot_rect,plot_sls,plot_line,plot_arro,plot_circ,plot_ptxt,plot_numb,&
      plotlib_supports_alpha,plot_set_opacity
 implicit none
 integer, intent(in) :: ipanel,irow,icolumn,itransx,itransy
 real,    intent(in) :: time
 integer :: icolourprev,linestyleprev,linewidthprev,ifillstyle
 integer :: i,j,ierr,iplotonthispanel,ndec,nc
 integer, parameter :: maxfuncpts = 1000
 real :: xmin,xmax,ymin,ymax,dxplot,dyplot
 real :: xpos,ypos,xlen,ylen,anglerad,dx,dy,fjust
 real, dimension(2) :: xline,yline
 real, dimension(maxfuncpts) :: xfunc,yfunc
 character(len=lentext) :: text
 character(len=30)      :: string
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
       if (plotlib_supports_alpha) call plot_set_opacity(shape(i)%opacity)

       anglerad = shape(i)%angle*(pi/180.)

       call convert_units(shape(i),xpos,ypos,xlen,ylen, &
                          xmin,ymin,dxplot,dyplot,itransx,itransy)

       !call print_shapeinfo(i,shape(i)%itype,shape(i))
       !print "(a)",'> plotting shape: '//trim(labelshapetype(shape(i)%itype))
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
             fjust = shape(i)%fjust
             call plot_arro(xpos-fjust*dx,ypos-fjust*dy,xpos+(1.-fjust)*dx,ypos+(1.-fjust)*dy)
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
          text = trim(shape(i)%text)
          !--handle special characters in text strings (e.g. replace %t with time)
          if (index(text,'%t').ne.0) then
             ndec = 3
             call plot_numb(nint(time/10.**(int(log10(time)-ndec))),int(log10(time)-ndec),1,string,nc)
             call string_replace(text,'%t',string(1:nc))
          endif
          call plot_ptxt(xpos,ypos,shape(i)%angle,shape(i)%fjust,trim(text))
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
 if (plotlib_supports_alpha) call plot_set_opacity(1.0)

end subroutine plot_shapes

!------------------------------------------------------------
! query function asking whether or not a point falls within
! a shape object
!------------------------------------------------------------
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

!---------------------------------------
! routine to edit shapes interactively
!---------------------------------------
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

!--------------------------------------------------------
! utility routine to add a new text shape interactively
!--------------------------------------------------------
subroutine add_textshape(xpt,ypt,itransx,itransy,ipanel,ierr)
 use plotlib, only:plot_qwin
 implicit none
 real, intent(in)     :: xpt,ypt
 integer, intent(in)  :: itransx,itransy,ipanel
 integer, intent(out) :: ierr
 integer :: i
 real :: xmin,xmax,ymin,ymax,xposi,yposi

 ierr = 0
 nshapes = nshapes + 1
 if (nshapes.gt.maxshapes) then
    print*,' *** cannot add shape: array limits reached, delete some shapes first ***'
    nshapes = maxshapes
    ierr = 1
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

!-----------------------------------------------------------------
! utility routine to convert between units used in shape plotting
!-----------------------------------------------------------------
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

!--------------------------------------------------------
! utility routine to edit a text object interactively
!--------------------------------------------------------
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
