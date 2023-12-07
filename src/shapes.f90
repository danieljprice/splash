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
!  Copyright (C) 2005-2019 Daniel Price. All rights reserved.
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
 integer, parameter, private :: maxshapetype = 8
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
      'f(x)     ', &
      'marker   '/)
 namelist /shapeopts/ nshapes,shape

 integer, parameter, public :: &
     ishape_square    = 1, &
     ishape_rectangle = 2, &
     ishape_arrow     = 3, &
     ishape_circle    = 4, &
     ishape_line      = 5, &
     ishape_text      = 6, &
     ishape_function  = 7, &
     ishape_marker    = 8

 integer, parameter, private :: maxunits = 2
 character(len=20), dimension(maxunits), &
    parameter, private :: labelunits = &
    (/'units of plot       ', &
      'viewport coordinates'/)
!      'inches             ', &
!      'millimeters        ', &
!      'pixels             '/)

 procedure(check_shapes), pointer, private :: checkshapes => null()
 procedure(add_shape),    pointer, private :: addshape => null()
 procedure(delete_shape), pointer, private :: delshape => null()

 real, parameter, private :: pi = 4.*atan(1.)

contains
!-----------------------------------------------------------------
! shape default settings
!-----------------------------------------------------------------
subroutine defaults_set_shapes

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

end subroutine defaults_set_shapes

!-----------------------------------------------------------------
! shape submenu
!-----------------------------------------------------------------
subroutine submenu_shapes()
 use promptlist, only:prompt_list

 checkshapes => check_shapes
 addshape => add_shape
 delshape => delete_shape
 call prompt_list(nshapes,maxshapes,'shape',checkshapes,addshape,delshape)

end subroutine submenu_shapes

!-----------------------------------
! print the current list of shapes
!-----------------------------------
subroutine check_shapes(nshape)
 integer, intent(in) :: nshape
 integer :: ishape

 print "(/,a)", ' Current list of plot annotations:'
 if (nshape > 0) then
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
 integer, intent(in) :: inum,itype
 type(shapedef), intent(in), optional :: shapein
 character(len=20), parameter :: fmtstring = "('Shape ',i2,': ',a)"

 select case(itype)
 case(ishape_square)
    print "(10x,a)",'    --'
    write(*,fmtstring,advance='no') inum,'   |  |   '//labelshapetype(itype)
    if (present(shapein)) then
       print "(5x,es10.2,' x ',es10.2)",shapein%xlen,shapein%ylen
    else
       print*
    endif
    print "(10x,a)",'    --'
 case(ishape_rectangle)
    print "(10x,a)",'  -----'
    write(*,fmtstring,advance='no') inum,' |     |  '//labelshapetype(itype)
    if (present(shapein)) then
       print "(5x,es10.2,' x ',es10.2)",shapein%xlen,shapein%ylen
    else
       print*
    endif
    print "(10x,a)",'  -----'
 case(ishape_arrow)
    print "(10x,a)"
    write(*,fmtstring,advance='no') inum,' ------>  '//labelshapetype(itype)
    if (present(shapein)) then
       print "(2x,'(x,y) = (',es10.2,es10.2,') to (x,y) = (',es10.2,es10.2,')')",&
             shapein%xpos,shapein%ypos,shapein%xlen,shapein%ylen
    else
       print*
    endif
    print "(10x,a)"
 case(ishape_circle)
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
 case(ishape_line)
    print "(10x,a)"
    write(*,fmtstring,advance='no') inum,' -------- '//labelshapetype(itype)
    if (present(shapein)) then
       print "(6x,'length = ',es10.2)",shapein%xlen
    else
       print*
    endif
    print "(10x,a)"
 case(ishape_text)
    print "(10x,a)",                  '        '
    write(*,fmtstring,advance='no') inum,'   TEXT   '
    if (present(shapein)) then
       print "('""',a,'""')",trim(shapein%text)
    else
       print*
    endif
    print "(10x,a)",                  '        '
 case(ishape_function)
    print "(10x,a)",'  _     / '
    write(*,fmtstring,advance='no') inum,' / \   /  '//trim(labelshapetype(itype))
    if (present(shapein)) then
       print "(' = ',a)",trim(shapein%text)
    else
       print*
    endif
    print "(10x,a)",'/   \_/     '
 case(ishape_marker)
    write(*,fmtstring,advance='no') inum,' (x)  '//labelshapetype(itype)
    print "(a)"
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
 integer, intent(in) :: istart,iend
 integer, intent(inout) :: nshape
 integer            :: i,ishape,itype,indexi,iunits,ierr,itry
 character(len=10)  :: poslabel
 character(len=80)  :: string

 itype   = 1
 ishape  = istart + 1
 if (ishape > maxshapes) then
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
    if (i < 10) then
       write(string(indexi:),"(1x,i1,') ',a)") i,trim(labelshapetype(i))
    else
       write(string(indexi:),"(1x,i2,') ',a)") i,trim(labelshapetype(i))
    endif
    indexi = len_trim(string) + 1
 enddo
 print "(/,a)",trim(string)
 !print "(i2,a)",(i,') '//trim(labelshapetype(i)),i=1,maxshapetype)

 over_shapes: do while(ishape <= iend .and. i <= maxshapes)
    if (istart==0 .or. shape(ishape)%itype <= 0 .or. shape(ishape)%itype > maxshapetype) then
       call prompt('choose an object type (0=none) ',shape(ishape)%itype,0,maxshapetype)
    endif
    itype = shape(ishape)%itype
    if (itype==0) then
       call delete_shape(ishape,nshape)
       exit over_shapes
    else
       call print_shapeinfo(ishape,itype)

       if (itype==ishape_function) then
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
       case(ishape_square) ! square
          call prompt('enter length of side (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen,0.)
          shape(ishape)%ylen = shape(ishape)%xlen
          poslabel = ' centre'
       case(ishape_rectangle) ! rectangle
          call prompt('enter x position of left edge (in '//trim(labelunits(iunits))//')',shape(ishape)%xpos)
          call prompt('enter x position of right edge (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen)
          call prompt('enter y position of bottom edge (in '//trim(labelunits(iunits))//')',shape(ishape)%ypos)
          call prompt('enter y position of top edge (in '//trim(labelunits(iunits))//')',shape(ishape)%ylen)
          poslabel = ''
       case(ishape_arrow) ! arrow
          call prompt('enter starting x position (in '//trim(labelunits(iunits))//') ',shape(ishape)%xpos)
          call prompt('enter starting y position (in '//trim(labelunits(iunits))//') ',shape(ishape)%ypos)
          call prompt('enter finishing x position (in '//trim(labelunits(iunits))//') ',shape(ishape)%xlen)
          call prompt('enter finishing y position (in '//trim(labelunits(iunits))//') ',shape(ishape)%ylen)
          call prompt('enter (optional) text string ',shape(ishape)%text)
          if (len_trim(shape(ishape)%text) > 0) then
             call prompt('enter justification factor (0.0=left 1.0=right)',shape(ishape)%fjust)
          endif
          poslabel = ''
       case(ishape_circle) ! circle
          call prompt('enter radius (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen,0.)
          poslabel = ' centre'
       case(ishape_line) ! line
          call prompt('enter line length (in '//trim(labelunits(iunits))//')',shape(ishape)%xlen,0.)
          call prompt('enter angle of line in degrees (0.0 = horizontal) ',shape(ishape)%angle)
          poslabel = ' starting'
       case(ishape_text) ! text
          call prompt('enter text string ',shape(ishape)%text)
          call prompt('enter angle for text in degrees (0 = horizontal) ',shape(ishape)%angle)
          call prompt('enter justification factor (0.0=left 1.0=right)',shape(ishape)%fjust)
          poslabel = ' starting'
       case(ishape_function) ! arbitrary function
          ierr = 1
          itry = 1
          do while(ierr /= 0 .and. itry <= 3)
             if (itry > 1) print "(a,i1,a)",'attempt ',itry,' of 3:'
             print "(a,6(/,11x,a),/)",' Examples: sin(2*pi*x)','sqrt(0.5*x)','x^2', &
             'exp(-2*x**2)','log10(x/2)','exp(p),p=sin(pi*x)','cos(z/d),z=acos(d),d=x^2'
             call prompt('enter function f(x) to plot ',shape(ishape)%text)
             call check_function(shape(ishape)%text,ierr)
             if (ierr /= 0 .and. len(shape(ishape)%text)==len_trim(shape(ishape)%text)) then
                print "(a,i3,a)",' (errors are probably because string is too long, max length = ',len(shape(ishape)%text),')'
             endif
             itry = itry + 1
          enddo
          if (ierr /= 0) then
             print "(a)",' *** too many tries, aborting ***'
             ishape = ishape - 1
             cycle over_shapes
          endif
          poslabel = ' starting'
       case(ishape_marker)
          print "(/,' Marker options (for all from -8->31, see plot library userguide):',12(/,i2,') ',a))", &
                      0,'square',1,'.',2,'+',3,'*',4,'o',5,'x',12,'5-pointed star',17,'bold circle',-8,'large bold circle'
          call prompt('Enter marker type ',shape(ishape)%ifillstyle)
          poslabel = ''
       end select
       if (itype /= ishape_function .and. itype /= ishape_rectangle .and. itype /= ishape_arrow) then
          call prompt('enter'//trim(poslabel)//' x position (in '//trim(labelunits(iunits))//') ',shape(ishape)%xpos)
          call prompt('enter'//trim(poslabel)//' y position (in '//trim(labelunits(iunits))//') ',shape(ishape)%ypos)
       elseif (itype==7) then
          call prompt('enter xmin for line segment (0=ignore)',shape(ishape)%xpos)
          call prompt('enter xmax for line segment (0=ignore)',shape(ishape)%xlen)
       endif
       if (itype==ishape_square .or. itype==ishape_rectangle .or. itype==ishape_circle) then
          call prompt('enter fill style (1=solid,2=outline,3=hatch,4=crosshatch) for '// &
                      trim(labelshapetype(itype)),shape(ishape)%ifillstyle,0,plotlib_maxfillstyle)
       endif
       if (itype /= ishape_text .and. itype /= ishape_marker) then
          call prompt('enter line style (1=solid,2=dash,3=dotdash,4=dot,5=dashdot) for '// &
                      trim(labelshapetype(itype)),shape(ishape)%linestyle,0,plotlib_maxlinestyle)
       endif
       if (itype==ishape_text .or. itype==ishape_marker) then
          call prompt('enter character height for '//trim(labelshapetype(itype)),shape(ishape)%xlen,0.,10.)
       else
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
       if (shape(ishape)%iplotonpanel < -2 .or. &
           shape(ishape)%iplotonpanel > maxplot) shape(:)%iplotonpanel = 0

       call prompt('Enter selection ',shape(ishape)%iplotonpanel,-2,maxplot)
       if (ishape > nshape) nshape = ishape
       ishape = ishape + 1
    endif
 enddo over_shapes

end subroutine add_shape

!------------------------------------------
! utility routine to delete a shape object
!------------------------------------------
subroutine delete_shape(ishape,nshape)
 integer, intent(in)    :: ishape
 integer, intent(inout) :: nshape
 integer :: i

 if (ishape > 0 .and. nshape > 0 .and. ishape <= maxshapes) then
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
subroutine plot_shapes(ipanel,irow,icolumn,itransx,itransy,time,nvar,allvars,tags)
 use exactfunction, only:exact_function
 use transforms,    only:transform_inverse,transform
 use asciiutils,    only:string_replace
 use plotlib, only:plot_qci,plot_qls,plot_qlw,plot_qfs,plot_qwin,plot_sci,plot_sfs,plot_slw, &
      plot_sci,plot_rect,plot_sls,plot_line,plot_arro,plot_circ,plot_ptxt,plot_numb,&
      plotlib_supports_alpha,plot_set_opacity,plot_pt1,plot_sch,plot_qch,plot_qcs
 use parsetext,     only:parse_text,rn
 use params,        only:ltag
 integer, intent(in) :: ipanel,irow,icolumn,itransx,itransy
 real,    intent(in) :: time
 integer, intent(in) :: nvar
 real, intent(in), dimension(nvar) :: allvars
 character(len=*), dimension(nvar), intent(in) :: tags
 integer :: icolourprev,linestyleprev,linewidthprev,ifillstyle
 integer :: i,j,ierr,iplotonthispanel
 integer, parameter :: maxfuncpts = 1000
 real :: xmin,xmax,ymin,ymax,dxplot,dyplot,charheightprev,ysign
 real :: xpos,ypos,xlen,ylen,anglerad,dx,xmini,xmaxi,xch,ych
 real, dimension(2) :: xline,yline
 real, dimension(maxfuncpts) :: xfunc,yfunc
 character(len=lentext) :: text
 real(kind=rn),       dimension(nvar+1) :: vals
 character(len=ltag), dimension(nvar+1) :: vars
!
!--store current settings
!
 call plot_qci(icolourprev)
 call plot_qls(linestyleprev)
 call plot_qlw(linewidthprev)
 call plot_qfs(ifillstyle)
 call plot_qch(charheightprev)
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
    if (iplotonthispanel==0 &
       .or.(iplotonthispanel > 0  .and. ipanel==iplotonthispanel) &
       .or.(iplotonthispanel==-1 .and. irow==1) &
       .or.(iplotonthispanel==-2 .and. icolumn==1)) then

       call plot_sci(shape(i)%icolour)
       call plot_sls(shape(i)%linestyle)
       if (shape(i)%itype==ishape_text .or. shape(i)%itype==ishape_marker) then
          call plot_sch(shape(i)%xlen)
       endif
       call plot_slw(shape(i)%linewidth)
       call plot_sfs(shape(i)%ifillstyle)
       if (plotlib_supports_alpha) call plot_set_opacity(shape(i)%opacity)

       anglerad = shape(i)%angle*(pi/180.)

       call convert_units(shape(i),xpos,ypos,xlen,ylen, &
                          xmin,ymin,dxplot,dyplot,itransx,itransy)

       !call print_shapeinfo(i,shape(i)%itype,shape(i))
       !print "(a)",'> plotting shape: '//trim(labelshapetype(shape(i)%itype))
       select case(shape(i)%itype)
       case(ishape_square,ishape_rectangle) ! square, rectangle
          if (xlen > dxplot .or. ylen > dyplot) then
             print "(2x,a)",'Error: shape size exceeds plot dimensions: not plotted'
          else
             if (shape(i)%itype==1) then
                call plot_rect(xpos-0.5*xlen,xpos+0.5*xlen,ypos-0.5*ylen,ypos + 0.5*ylen)
             else
                call plot_rect(xpos,xlen,ypos,ylen)
             endif
          endif
       case(ishape_arrow) ! arrow
          !print "(1x,2(a,es10.2,es10.2))",trim(shape(i)%text),xpos,ypos,'->',xlen,ylen
          call plot_arro(xpos,ypos,xlen,ylen)
          !--plot text label under the arrow, if present
          call plot_qcs(4,xch,ych)
          text = trim(shape(i)%text)
          if (len_trim(text) > 0) then
             ysign = max(sign(1.,ylen-ypos),-0.15)
             call plot_ptxt(xpos,ypos-ysign*ych,shape(i)%angle,shape(i)%fjust,trim(text))
          endif
       case(ishape_circle) ! circle
          if (xlen > dxplot .or. xlen > dyplot) then
             print "(2x,a)",'Error: circle radius exceeds plot dimensions: circle not plotted'
          else
             call plot_circ(xpos,ypos,xlen)
             if (shape(i)%ifillstyle > 2) then
                call plot_sfs(2)  ! also plot outline if fill style is hatched
                call plot_circ(xpos,ypos,xlen)
             endif
          endif
       case(ishape_line) ! line
          xline(1) = xpos
          yline(1) = ypos
          xline(2) = xpos + xlen*cos(anglerad)
          yline(2) = ypos + xlen*sin(anglerad)
          call plot_line(2,xline,yline)
       case(ishape_text) ! text
          text = trim(shape(i)%text)
          !--handle special characters in text strings (e.g. replace %t with time)
          if (index(text,'%') /= 0) then
             vars = (/'t'/)
             vals(1) = real(time,kind=rn)
             if (nvar > 0) then
                vars(2:1+nvar) = tags(1:nvar)
                vals(2:1+nvar) = allvars(1:nvar)
             endif
             call parse_text(text,vars,vals)
          endif
          !--plot text string
          call plot_ptxt(xpos,ypos,shape(i)%angle,shape(i)%fjust,trim(text))
       case(ishape_function) ! arbitrary function
          !--set x to be evenly spaced in transformed (plot) coordinates
          xmini = xmin
          xmaxi = xmax
          if (abs(shape(i)%xpos) > 0.) xmini = shape(i)%xpos
          if (abs(shape(i)%xlen) > 0.) xmaxi = shape(i)%xlen
          dx = (xmaxi-xmini)/real(maxfuncpts - 1)
          do j=1,maxfuncpts
             xfunc(j) = xmini + (j-1)*dx
          enddo
          !--transform x array back to untransformed space to evaluate f(x)
          if (itransx > 0) call transform_inverse(xfunc,itransx)
          call exact_function(shape(i)%text,xfunc,yfunc,0.,ierr)
          if (ierr==0) then
             !--reset x values
             do j=1,maxfuncpts
                xfunc(j) = xmini + (j-1)*dx
             enddo
             !--transform y if necessary
             if (itransy > 0) call transform(yfunc,itransy)
             !--plot the line
             call plot_line(maxfuncpts,xfunc,yfunc)
          endif
       case(ishape_marker) ! marker
          call plot_pt1(xpos,ypos,shape(i)%ifillstyle)
       end select
    endif
 enddo

 call plot_sci(icolourprev)
 call plot_sls(linestyleprev)
 call plot_slw(linewidthprev)
 call plot_sfs(ifillstyle)
 call plot_sch(charheightprev)
 if (plotlib_supports_alpha) call plot_set_opacity(1.0)

end subroutine plot_shapes

!------------------------------------------------------------
! query function asking whether or not a point falls within
! a shape object
!------------------------------------------------------------
integer function inshape(xpt,ypt,itransx,itransy,xmin,xmax,ymin,ymax)
 use plotlib, only:plot_qwin,plot_qtxt,plot_qcs
 real, intent(in) :: xpt,ypt,xmin,ymin,xmax,ymax
 integer, intent(in) :: itransx,itransy
 integer :: i
 real :: xpos,ypos,xlen,ylen,ysign,xch,ych
 real :: dxplot,dyplot
 real, dimension(4) :: xbox,ybox

 dxplot = xmax - xmin
 dyplot = ymax - ymin

 inshape = 0
 do i=1,nshapes

    call convert_units(shape(i),xpos,ypos,xlen,ylen, &
                       xmin,ymin,dxplot,dyplot,itransx,itransy)

    select case(shape(i)%itype)
    case(ishape_square,ishape_rectangle)  ! square, rectangle

    case(ishape_arrow) ! arrow
       if (near_pt(xpos,ypos,xpt,ypt,dxplot,dyplot) .or. &
           near_pt(xlen,ylen,xpt,ypt,dxplot,dyplot)) then
          inshape = i
       endif
       call plot_qtxt(xpos,ypos,shape(i)%angle,shape(i)%fjust,trim(shape(i)%text),xbox,ybox)
       call plot_qcs(4,xch,ych)
       ysign = max(sign(1.,ylen-ypos),0.)
       call plot_qtxt(xpos,ypos-ysign*ych,shape(i)%angle,shape(i)%fjust,trim(shape(i)%text),xbox,ybox)
       if (xpt > minval(xbox) .and. xpt <= maxval(xbox) &
           .and. ypt > minval(ybox) .and. ypt <= maxval(ybox)) then
          inshape = i
       endif
    case(ishape_circle) ! circle

    case(ishape_line) ! line

    case(ishape_text) ! text
       call plot_qtxt(xpos,ypos,shape(i)%angle,shape(i)%fjust,trim(shape(i)%text),xbox,ybox)
       if (xpt > minval(xbox) .and. xpt <= maxval(xbox) &
          .and. ypt > minval(ybox) .and. ypt <= maxval(ybox)) then
          inshape = i
       endif
    end select
 enddo

end function inshape

!--------------------------------------------
! work out whether a mouse click is "close"
! to a plotted point
!--------------------------------------------
logical function near_pt(x1,y1,xc,yc,dx,dy)
 real, intent(in) :: x1,y1,xc,yc,dx,dy
 real :: xtol,ytol

 xtol = dx/128.  ! as in 1/128th of the page width
 ytol = dy/128.  ! as in 1/128th of the page height
 near_pt = ((x1-xc)**2 < xtol*xtol) .and. &
           ((y1-yc)**2 < ytol*ytol)

end function near_pt

!---------------------------------------
! routine to edit shapes interactively
!---------------------------------------
subroutine edit_shape(i,xpt,ypt,itransx,itransy,first)
 use plotlib, only:plot_qwin,plot_qtxt,plot_qcs
 integer, intent(in) :: i,itransx,itransy
 real, intent(in)    :: xpt,ypt
 logical, intent(in) :: first
 real :: xmin,xmax,ymin,ymax,dxplot,dyplot,xlen,ylen
 real :: xpos,ypos,xbox(4),ybox(4),xch,ych,ysign

 call plot_qwin(xmin,xmax,ymin,ymax)
 dxplot = xmax - xmin
 dyplot = ymax - ymin

 call convert_units(shape(i),xpos,ypos,xlen,ylen, &
                    xmin,ymin,dxplot,dyplot,itransx,itransy)

 select case(shape(i)%itype)
 case(ishape_text)
    call edit_textbox(xpos,ypos,shape(i)%angle,shape(i)%fjust,shape(i)%text)
 case(ishape_arrow)
    call plot_qtxt(xpos,ypos,shape(i)%angle,shape(i)%fjust,trim(shape(i)%text),xbox,ybox)
    call plot_qcs(4,xch,ych)
    ysign = max(sign(1.,ylen-ypos),0.0)
    if (near_pt(xlen,ylen,xpt,ypt,dxplot,dyplot)) then
       ! if click on head of arrow, move the tail
       call edit_arrow(xlen,ylen,shape(i)%xpos,shape(i)%ypos)
    elseif (xpt > minval(xbox) .and. xpt <= maxval(xbox) &
      .and. ypt-0.5*ysign*ych > minval(ybox) .and. ypt-0.5*ysign*ych <= maxval(ybox)) then
       ! if click on text box, edit text
       call edit_textbox(xpos,ypos,shape(i)%angle,shape(i)%fjust,shape(i)%text)
    else
       ! if click on tail of arrow, move the head
       call edit_arrow(xpos,ypos,shape(i)%xlen,shape(i)%ylen)
       if (first) call edit_textbox(xpos,ypos,shape(i)%angle,shape(i)%fjust,shape(i)%text)
    endif
 end select

end subroutine edit_shape

!--------------------------------------------------------
! utility routine to add a new text shape interactively
!--------------------------------------------------------
subroutine add_shape_interactive(xpt,ypt,itransx,itransy,ipanel,ierr,shape_type)
 use plotlib, only:plot_qwin,plot_qch
 real, intent(in)     :: xpt,ypt
 integer, intent(in)  :: itransx,itransy,ipanel
 integer, intent(out) :: ierr
 integer, intent(in), optional :: shape_type
 integer :: i,itype
 real :: xmin,xmax,ymin,ymax,xposi,yposi,charheight

 itype = ishape_text ! text shape by default
 if (present(shape_type)) itype = shape_type

 ierr = 0
 nshapes = nshapes + 1
 if (nshapes > maxshapes) then
    print*,' *** cannot add shape: array limits reached, delete some shapes first ***'
    nshapes = maxshapes
    ierr = 1
    return
 endif
 i = nshapes
 shape(i)%itype = itype
 print*,' adding shape '//trim(labelshapetype(shape(i)%itype))
 shape(i)%icolour = 1
 shape(i)%linestyle = 1
 shape(i)%linewidth = 1
 shape(i)%ifillstyle = 2
 shape(i)%iplotonpanel = ipanel
!
!--position shapes relative to viewport
!
 shape(i)%iunits = 2
 call plot_qwin(xmin,xmax,ymin,ymax)
 xposi = (xpt - xmin)/(xmax-xmin)
 yposi = (ypt - ymin)/(ymax-ymin)
 shape(i)%xpos = xposi
 shape(i)%ypos = yposi

 call plot_qch(charheight)
 shape(i)%xlen = charheight
 shape(i)%ylen = 1.
 shape(i)%angle = 0.
 shape(i)%text = 'click to edit'
 shape(i)%fjust = 0.
 if (itype==ishape_arrow) shape(i)%fjust = 0.5
 call edit_shape(i,xposi,yposi,itransx,itransy,first=.true.)

end subroutine add_shape_interactive

!------------------------------------------------------------
! utility routine to add a new text shape non-interactively
!-----------------------------------------------------------
subroutine add_text(xpos,ypos,string)
 use plotlib, only:plot_qch
 real, intent(in) :: xpos,ypos
 character(len=*), intent(in) :: string
 integer :: i,ierr
 real :: charheight

 ierr = 0
 ! do nothing if string is blank
 if (len_trim(string) <= 0) return

 call delete_text(string) ! delete if already exists
 nshapes = nshapes + 1
 if (nshapes > maxshapes) then
    print*,' *** cannot add shape: array limits reached, delete some shapes first ***'
    nshapes = maxshapes
    ierr = 1
    return
 endif
 i = nshapes
 shape(i)%itype = 6
 !print*,' adding shape '//trim(labelshapetype(shape(i)%itype))
 shape(i)%icolour = 1
 shape(i)%linestyle = 1
 shape(i)%linewidth = 1
 shape(i)%ifillstyle = 2
 shape(i)%iplotonpanel = 1
 shape(i)%iunits = 2
 shape(i)%xpos = xpos
 shape(i)%ypos = ypos

 call plot_qch(charheight)
 shape(i)%xlen = charheight
 shape(i)%ylen = 1.
 shape(i)%angle = 0.
 shape(i)%text = string
 shape(i)%fjust = 0.

end subroutine add_text

!------------------------------------------------------------
! delete matching text shape
!-----------------------------------------------------------
subroutine delete_text(string)
 character(len=*), intent(in) :: string
 integer :: i,nshapesold

 nshapesold = nshapes
 do i=1,nshapesold
    if (shape(i)%itype == 6 .and. shape(i)%text == trim(string)) then
       call delete_shape(i,nshapes)
    endif
 enddo

end subroutine delete_text

!-----------------------------------------------------------------
! utility routine to convert between units used in shape plotting
!-----------------------------------------------------------------
subroutine convert_units(shape,xpos,ypos,xlen,ylen,xmin,ymin,dxplot,dyplot,itransx,itransy)
 use transforms, only:transform
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
    if (shape%itype == ishape_circle) then
       xlen = xlen*dxplot
       ylen = ylen*dyplot
    else
       xlen = xmin + xlen*dxplot
       ylen = ymin + ylen*dyplot
    endif
 case(1)
    if (itransx > 0) then
       call transform(xpos,itransx)
       call transform(xlen,itransx)
    endif
    if (itransy > 0) then
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
subroutine edit_textbox(xpt,ypt,angle,fjust,string)
 use plotlib, only:plot_stbg,plot_ptxt,plot_curs
 real, intent(in) :: xpt,ypt,angle,fjust
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
 call plot_ptxt(xpt,ypt,angle,fjust,string(1:i)//'_')

 xpt2 = xpt
 ypt2 = ypt
 ierr = plot_curs(xpt2,ypt2,mychar)
 do while (mychar /= achar(13) &   ! carriage return
     .and. mychar /= achar(27) &   ! ctrl-c
     .and. mychar /= achar(3))     ! esc
    if (mychar==achar(8)) then   ! backspace
       i = max(i - 1,1)
       string(i:i) = '_'
       call plot_ptxt(xpt,ypt,angle,fjust,string(1:i))
       string(i:i) = ' '
    else
       if (trim(string)=='click to edit') then
          !print*,'erasing string'
          string = ' '
          i = 1
       endif
       string(i:i) = mychar
       call plot_ptxt(xpt,ypt,angle,fjust,string(1:i))
       i = min(i + 1,len(string))
       if (i==len(string)) print*,' reached end of string'
    endif
    ierr = plot_curs(xpt2,ypt2,mychar)
 enddo

 !--if ctrl-c or esc, restore original string
 if (mychar==achar(3) .or. mychar==achar(27)) then
    string = oldstring
    print*,'cancelled'
 else
    print*,'done: text = "'//trim(string)//'"'
 endif
 call plot_stbg(-1)

end subroutine edit_textbox

!--------------------------------------------------------
! utility routine to edit a text object interactively
!--------------------------------------------------------
subroutine edit_arrow(xpt,ypt,x2,y2)
 use plotlib, only:plot_band,plot_qwin
 real, intent(in)  :: xpt,ypt
 real, intent(out) :: x2,y2
 character(len=1) :: char2
 integer :: ierr
 real :: xmin,xmax,ymin,ymax

 ! get point in plot units
 ierr = plot_band(1,1,xpt,ypt,x2,y2,char2)

 ! transform point to viewport coordinates
 call plot_qwin(xmin,xmax,ymin,ymax)
 x2 = (x2 - xmin)/(xmax-xmin)
 y2 = (y2 - ymin)/(ymax-ymin)

end subroutine edit_arrow

end module shapes
