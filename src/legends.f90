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
!     module containing routines for plotting legends in PGPLOT
!     subroutines:
!      legend         : plots time on plots
!      legend_vec     : plots legend with vector arrow
!      legend_markers : plots different particle types
!      legend_scale   : plots a scale on co-ordinate plots
!-----------------------------------------------------------------
module legends
 implicit none
 
 public :: legend, legend_vec, legend_markers, legend_scale
 public :: prompt_panelselect, ipanelselect
 
 private

contains

!-----------------------------------------------------------------
!     plots time on plot
!     arguments:
!           t : current time
!        hpos : horizontal position as fraction of viewport
!        vpos : vertical position in character heights from top
!-----------------------------------------------------------------

subroutine legend(legendtext,t,unitslabel,hpos,vpos,fjust,usebox)
 use plotlib,    only:plot_annotate
 use asciiutils, only:string_replace
 use parsetext,  only:parse_text,rn
 implicit none
 real, intent(in) :: t,hpos,vpos,fjust
 character(len=*), intent(in) :: legendtext,unitslabel
 logical,          intent(in) :: usebox
 character(len=len(legendtext)+len(unitslabel)+20) :: label
 integer, parameter :: nvars = 1
 real(kind=rn),    dimension(nvars) :: vals
 character(len=1), dimension(nvars) :: vars
 
 label = trim(legendtext)

 !
 !  if string does not contain any formatting
 !  append the time variable to it
 !
 if (index(label,'%') <= 0) then
    if (t < 1.) then
       label = trim(label)//'%t.2'
    else
       label = trim(label)//'%t.3'    
    endif
 endif
 !
 !  parse string for functions of time and formatting
 !  i.e. %t.5
 !
 vars = (/'t'/)
 vals(1) = real(t,kind=rn)
 call parse_text(label,vars,vals)

 if (index(label,'%ut').gt.0) then
    call string_replace(label,'%ut',trim(unitslabel))
 else
    label = trim(label)//trim(unitslabel)
 endif
 
 if (usebox) call plot_box_around_text(trim(label),hpos,vpos,fjust)
 call plot_annotate('T',-vpos,hpos,fjust,trim(label))

 return
end subroutine legend

!-----------------------------------------------------------------
!  utility routine for plotting translucent box in legends
!-----------------------------------------------------------------
subroutine plot_box_around_text(string,hpos,vpos,fjust)
 use plotlib, only:plot_qwin,plot_qcs,plot_qtxt,plot_qci,plot_sci,plot_sfs, &
              plot_set_opacity, plot_rect
 implicit none
 character(len=*), intent(in) :: string
 real, intent(in) :: hpos,vpos,fjust
 real :: xmin,xmax,ymin,ymax,xpos,ypos
 real :: xbuf,ybuf,dx,dy,xch,ych,x1,x2,y1,y2
 real, dimension(4) :: xbox,ybox
 integer :: ic
!
!--convert hpos and vpos to x, y to plot arrow
!
 call plot_qwin(xmin,xmax,ymin,ymax)
 xpos = xmin + hpos*(xmax-xmin)
 call plot_qcs(4,xch,ych)
 ypos = ymax - (vpos + 1.)*ych
!
!--enquire bounding box of string
!
 call plot_qtxt(xpos,ypos,0.0,0.0,trim(string),xbox,ybox)

 xbuf = 0.25*xch
 ybuf = 0.5*ych
 dx = xbox(3) - xbox(1)
 dy = ybox(3) - ybox(1) + 0.25*ych
 x1 = xpos - fjust*dx - xbuf
 x2 = x1 + dx + 2.*xbuf
 y1 = ypos
 y2 = y1 + dy + ybuf
!
!--draw box around the string
!
 call plot_qci(ic) ! query colour index
 call plot_sci(0)  ! background colour
 call plot_sfs(1)  ! solid fill style
 call plot_set_opacity(0.5)
 call plot_rect(x1,x2,y1,y2,0.2*ych) ! draw a (rounded) rectangle
 call plot_set_opacity(1.0)
 call plot_sci(ic) ! restore colour index

end subroutine plot_box_around_text

!-----------------------------------------------------------------
!     plots vector plot legend
!     arguments:
!           t : current time
!        hpos : horizontal position as fraction of viewport
!        vpos : vertical position in character heights from top
!  charheight : this is the text character height
!               (legend_vec is called directly after plotting the arrows,
!                which may use a different character size - so we plot
!                the arrow here in the same way, but then revert to
!                the text character height to write the text)
!-----------------------------------------------------------------

subroutine legend_vec(label,unitslabel,vecmax,dx,hpos,vpos,charheight)
 use plotlib, only:plot_qwin,plot_qch,plot_sch,plot_qcs,plot_numb,plot_qtxt, &
                   plot_qci,plot_sci,plot_sfs,plot_rect,plot_sci,plot_text, &
                   plot_qvp,plot_svp,plot_swin,plot_arro,plot_set_opacity
 implicit none
 real, intent(in) :: vecmax,dx,hpos,vpos,charheight
 character(len=*), intent(in) :: label,unitslabel
 real :: xmin,xmax,ymin,ymax
 real :: xch,ych,charheightarrow,adjustlength,vecmaxnew
 real :: xpos,ypos,xbox(4),ybox(4),dxlabel,dxstring
 real :: dxbuffer,dybuffer,dxbox,dybox
 real :: xminnew,xmaxnew,yminnew,ymaxnew,x1,x2,y1,y2
 integer :: icolindex,mm,pp,nc,ndec
 character(len=len(label)+20) :: string

!
!--convert hpos and vpos to x, y to plot arrow
!
 call plot_qwin(xmin,xmax,ymin,ymax)
 call plot_qch(charheightarrow)
 call plot_sch(charheight)
 xpos = xmin + hpos*(xmax-xmin)
 call plot_qcs(4,xch,ych)
 ypos = ymax - (vpos + 1.)*ych
!
!--format string containing numerical value
!  vecmax corresponds to arrow of length dx
!  we will draw an arrow of length sqrt(dx^2 + ych^2)
!  so adjust vecmax accordingly
!
 adjustlength = sqrt(0.5*dx**2 + ych**2)/dx
 vecmaxnew = adjustlength*vecmax
 ndec = 2
 if (vecmaxnew.lt.tiny(vecmaxnew)) then
    string = '0'
    nc = 1
 else
    mm=int(vecmaxnew/10.**(int(log10(vecmaxnew)-ndec)))
    pp=int(log10(vecmaxnew)-ndec)
    call plot_numb(mm,pp,0,string,nc)
 endif
 string = '='//trim(string)
! write(string,"('=',1pe7.1)") vecmax
!
!--enquire size of label
!
 call plot_qtxt(xpos,ypos,0.0,0.0,trim(label),xbox,ybox)
 dxlabel = xbox(3) - xbox(2) + 0.5*xch
!
!--enquire size of string
!
 call plot_qtxt(xpos,ypos,0.0,0.0,trim(string),xbox,ybox)
 dxstring = xbox(3) - xbox(2)
!
!--set size of box in x direction
!
 dxbuffer = 0.25*xch ! these are size of margins (x and y)
 dybuffer = 0.25*ych
 dxbox = dxlabel + dxstring + 1.1*dx/sqrt(2.) + dxbuffer
 dybox = ych + 0.5*dybuffer
!
!--draw box around all of the legend
!
 call plot_qci(icolindex)
! draw a (rounded) rectangle in the background colour with solid fill style
 call plot_sci(0)
 call plot_sfs(1)
 call plot_set_opacity(0.66)
 call plot_rect(xpos-dxbuffer,xpos+dxbox,ypos-dybuffer,ypos + dybox,0.33*ych)
 call plot_set_opacity(1.0)
! change to foreground colour index
 call plot_sci(1)
! draw an outline around the box
! call pgsfs(2)
! call pgrect(xpos-dxbuffer,xpos+dxbox,ypos-dybuffer,ypos + dybox)
! call pgsfs(1)
!
!--write label
!
 call plot_text(xpos,ypos,trim(label))
 xpos = xpos + dxlabel
!
!--Draw arrow. Here we have to perform tricks to get the arrow
!  to appear even if outside the usual plotting area
!
!--save viewport settings
 call plot_qvp(0,x1,x2,y1,y2)
!--now allow the whole screen to be the viewport...
 call plot_svp(0.0,1.0,0.0,1.0)
!  ...but correspondingly adjust window so that x and y positions
!  are the same as in the old viewport
 xminnew = xmin - x1*(xmax-xmin)/(x2-x1)
 xmaxnew = xmax + (1.-x2)*(xmax-xmin)/(x2-x1)
 yminnew = ymin - y1*(ymax-ymin)/(y2-y1)
 ymaxnew = ymax + (1.-y2)*(ymax-ymin)/(y2-y1)
 call plot_swin(xminnew,xmaxnew,yminnew,ymaxnew)
!--use character height original arrows were drawn with
!  (this is to get the arrow head size right)
 call plot_sch(charheightarrow)
!--draw arrow
 call plot_arro(xpos,ypos,xpos + dx/sqrt(2.),ypos + ych)
!--restore viewport settings
 call plot_svp(x1,x2,y1,y2)
 call plot_swin(xmin,xmax,ymin,ymax)
 xpos = xpos + 1.1*dx/sqrt(2.)
!
!--write numerical value and units label
!
 call plot_sch(charheight)
!! call pgmtext('t',-vpos,hpos+0.02,0.0,trim(string))
 call plot_text(xpos,ypos,trim(string)//trim(unitslabel))
!
!--restore colour index
 call plot_sci(icolindex)

 return
end subroutine legend_vec

!-------------------------------------------------------------------------
!  draw a legend for different line/marker styles
!  uses current line style and colour
!  plots this below the time legend
!-------------------------------------------------------------------------
subroutine legend_markers(icall,icolour,imarkerstyle,ilinestyle, &
           iplotpts,iplotline,text,hposlegend,vposlegend,alphalegend)
  use plotlib, only:plot_qwin,plot_qcs,plot_qci,plot_qls,plot_sci,plot_sls, &
                    plot_line,plot_pt,plot_text,plot_stbg,plot_slc,plot_qlc,plot_set_opacity
  implicit none
  integer, intent(in) :: icall,icolour,imarkerstyle,ilinestyle
  logical, intent(in) :: iplotpts,iplotline
  character(len=*), intent(in) :: text
  real, intent(in) :: hposlegend,vposlegend,alphalegend
  integer :: icolourprev, ilinestyleprev,ilinecapprev
  real, dimension(3) :: xline,yline
  real :: xch, ych, xmin, xmax, ymin, ymax
  real :: vspace, vpos
!
!--do not plot anything if string is blank
!
  if (len_trim(text).le.0) return
  !call pgstbg(0)           ! opaque text to overwrite previous
!
!--set horizontal and vertical position and spacing
!  in units of the character height
!
  vspace = 1.5  ! (in units of character heights)
  vpos = vposlegend + icall*vspace + 0.5 ! distance from top, in units of char height

  call plot_qwin(xmin,xmax,ymin,ymax) ! query xmax, ymax
  call plot_qcs(4,xch,ych) ! query character height in x and y units
  call plot_qci(icolourprev)     ! save current colour index
  call plot_qls(ilinestyleprev)  ! save current line style
  call plot_qlc(ilinecapprev)    ! save the current line cap


  yline(:) = ymax - ((vpos - 0.5)*ych)
  xline(1) = xmin + hposlegend*(xmax-xmin)
  xline(2) = xline(1) + 1.5*xch
  xline(3) = xline(1) + 3.*xch

  call plot_sci(icolour)
  call plot_set_opacity(alphalegend)
  call plot_sls(ilinestyle)
!
!--set round caps
!
  !call plot_slc(1)
!
!--draw a small line segment
!
  if (iplotline) call plot_line(3,xline,yline)

  call plot_slc(ilinecapprev)
  call plot_sls(ilinestyleprev)
!
!--draw points, only two if line is also plotted so that you can see the line
!               three otherwise
!
  if (iplotpts .and. iplotline) then
     xline(2) = xline(3)
     call plot_pt(2,xline(1:2),yline(1:2),imarkerstyle)
  elseif (iplotpts) then
     call plot_pt(3,xline,yline,imarkerstyle)
  endif
!
!--add text
!
  if (iplotline .or. iplotpts .and. len_trim(text).gt.0) then
     call plot_text(xline(3) + 0.75*xch,yline(1)-0.25*ych,trim(text))
  endif

  call plot_sci(icolourprev)    ! reset colour index
  call plot_set_opacity(1.0)
  call plot_stbg(-1) ! reset text background to transparent

end subroutine legend_markers

!-------------------------------------------------------------------
!     plots labelled scale (horizontal error bar of a given length)
!     can be used on co-ordinate plots to give a length scale
!
!     e.g. would produce something like:
!
!                   |----|
!                   10 AU
!
!     arguments:
!        dxscale : length of scale in current x units
!        hpos : horizontal position as fraction of viewport
!        vpos : vertical position in character heights from top
!        text : label to print above scale
!-----------------------------------------------------------------
subroutine legend_scale(dxscale,hpos,vpos,text)
  use plotlib, only:plot_qwin,plot_qcs,plot_err1,plot_annotate
  implicit none
  real, intent(in) :: dxscale,hpos,vpos
  character(len=*), intent(in) :: text
  real :: xmin,xmax,ymin,ymax,xch,ych,xpos,ypos

  call plot_qwin(xmin,xmax,ymin,ymax)
  if (dxscale.gt.(xmax-xmin)) then
     print "(a)",'Error: scale size exceeds x dimensions: scale not plotted'
  else
     call plot_qcs(4,xch,ych)
     !--draw horizontal "error bar" above text
     ypos = ymin + (vpos+1.25)*ych
     xpos = xmin + hpos*(xmax-xmin)
     call plot_err1(5,xpos,ypos,0.5*dxscale,1.0)

     !--write text at the position specified
     call plot_annotate('B',-vpos,hpos,0.5,trim(text))
  endif

end subroutine legend_scale

!-------------------------------------------------------------------
!  The following subroutines handle the plotting of annotation
!  and legends only on particular panels
!-------------------------------------------------------------------
subroutine prompt_panelselect(string,iselect)
 use prompting, only:prompt
 implicit none
 character(len=*), intent(in) :: string
 integer,       intent(inout) :: iselect
 
 print "(4(/,a))", &
   '  0 : plot '//trim(string)//' on every panel ', &
   '  n : plot '//trim(string)//' on nth panel only ', &
   ' -1 : plot '//trim(string)//' on first row only ', &
   ' -2 : plot '//trim(string)//' on first column only '
 call prompt('Enter selection ',iselect,-2)

end subroutine prompt_panelselect

!-------------------------------------------------------------------
!  Function that evaluates the logic required to determine
!  whether the annotation should be plotted on the current panel
!  as per the prompts in prompt_panelselect
!-------------------------------------------------------------------
logical function ipanelselect(iselect,ipanel,irow,icolumn)
 implicit none
 integer,       intent(in) :: iselect,ipanel,irow,icolumn

 ipanelselect = ((iselect.gt.0 .and. ipanel.eq.iselect) &
              .or.(iselect.eq.-1 .and. irow.eq.1) &
              .or.(iselect.eq.-2 .and. icolumn.eq.1) &
              .or.(iselect.eq.0))

end function ipanelselect

end module legends
