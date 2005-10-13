!-----------------------------------------------------------------
!     module containing routines for plotting legends
!     subroutines:
!      legend      : plots time on plots
!      legend_vec  : plots legend with vector arrow
!      legend_part : plots different particle types
!-----------------------------------------------------------------
module legends
 implicit none

contains

!-----------------------------------------------------------------
!     plots time on plot
!     arguments:
!           t : current time
!        hpos : horizontal position as fraction of viewport
!        vpos : vertical position in character heights from top
!-----------------------------------------------------------------

subroutine legend(t)
 use settings_page, only:hposlegend, vposlegend, legendtext
 implicit none
 real, intent(in) :: t    
 integer :: mm,pp,nc,ndecimal,ndec
 real :: tplot
 character(len=30) :: string

 ndecimal = 2        ! number of decimal places to display
 ndec = 10**ndecimal
 if (t.eq.0.0) then
    tplot = 1e-6
 else
    tplot = t    !/(2.*3.1415926536)
 endif
 mm=nint(tplot*ndec)
 pp=nint(log10(tplot)-log10(tplot*ndec))
 call pgnumb(mm,pp,1,string,nc)
 call pgmtext('T',-vposlegend,hposlegend,0.0,trim(legendtext)//string(1:nc))

 return
end subroutine legend

!-----------------------------------------------------------------
!     plots vector plot legend
!     arguments:
!           t : current time
!        hpos : horizontal position as fraction of viewport
!        vpos : vertical position in character heights from top
!-----------------------------------------------------------------

subroutine legend_vec(label,vecmax,dx,hpos,vpos,charheight)
 implicit none
 real, intent(in) :: vecmax,dx,hpos,vpos,charheight
 character(len=*), intent(in) :: label
 real :: xmin,xmax,ymin,ymax
 real :: xch,ych,charheightarrow
 real :: xpos,ypos,xbox(4),ybox(4),dxlabel,dxstring
 integer :: icolindex
 character(len=len(label)+10) :: string

 write(string,"('=',1pe7.1)") vecmax
!
!--convert hpos and vpos to x, y to plot arrow
!
 call pgqwin(xmin,xmax,ymin,ymax)
 call pgqch(charheightarrow)
 call pgsch(charheight)
 xpos = xmin + hpos*(xmax-xmin)
 call pgqcs(4,xch,ych) 
 ypos = ymax - (vpos + 1.)*ych
!
!--enquire size of label
!
 call pgqtxt(xpos,ypos,0.0,0.0,trim(label),xbox,ybox)
 dxlabel = xbox(3) - xbox(2) + 0.5*xch
!
!--enquire size of string
!
 call pgqtxt(xpos,ypos,0.0,0.0,trim(string),xbox,ybox)
 dxstring = xbox(3) - xbox(2)
!
!--draw box around all of the legend
!
 call pgqci(icolindex)
 call pgsci(0)
 call pgrect(xpos-0.25*dx,xpos+dxlabel+dxstring+dx,ypos-0.5*dx,ypos + dx)
 call pgsci(icolindex)
 call pgsfs(2)
 call pgrect(xpos-0.25*dx,xpos+dxlabel+dxstring+dx,ypos-0.5*dx,ypos + dx)
 call pgsfs(1)
!
!--write label
!
 call pgtext(xpos,ypos,trim(label))
 xpos = xpos + dxlabel
!
!--draw arrow
!
 call pgsch(charheightarrow)
 call pgarro(xpos,ypos,xpos+dx/sqrt(2.),ypos+dx/sqrt(2.))
 xpos = xpos + dx/sqrt(2.)
!
!--write numerical value
!
 call pgsch(charheight)
!! call pgmtext('t',-vpos,hpos+0.02,0.0,trim(string))
 call pgtext(xpos,ypos,trim(string))

 return
end subroutine legend_vec

!-------------------------------------------------------------------------
!  draw a legend for different line/marker styles
!  uses current line style and colour
!  plots this below the time legend
!-------------------------------------------------------------------------
subroutine legend_markers(icall,icolour,imarkerstyle,ilinestyle,iplotpts,iplotline,text)
  use settings_page, only:hposlegend, vposlegend
  implicit none
  integer, intent(in) :: icall,icolour,imarkerstyle,ilinestyle
  logical, intent(in) :: iplotpts,iplotline
  character(len=*), intent(in) :: text
  integer :: icolourprev, ilinestyleprev
  real, dimension(3) :: xline,yline
  real :: xch, ych, xmin, xmax, ymin, ymax
  real :: vspace, vpos

  call pgstbg(0)           ! opaque text to overwrite previous
!
!--set horizontal and vertical position and spacing
!  in units of the character height
!
  vspace = 1.5  ! (in units of character heights)
  vpos = vposlegend + icall*vspace + 0.5 ! distance from top, in units of char height

  call pgqwin(xmin,xmax,ymin,ymax) ! query xmax, ymax
  call pgqcs(4,xch,ych) ! query character height in x and y units 
  call pgqci(icolourprev)     ! save current colour index
  call pgqls(ilinestyleprev)  ! save current line style

  yline(:) = ymax - ((vpos - 0.5)*ych)
  xline(1) = xmin + hposlegend*(xmax-xmin)
  xline(2) = xline(1) + xch
  xline(3) = xline(1) + 2.*xch

  call pgsci(icolour)
  call pgsls(ilinestyle)
!
!--draw a small line segment
!
  if (iplotline) call pgline(3,xline,yline)
!
!--draw points, only two if line is also plotted so that you can see the line
!               three otherwise
!
  if (iplotpts .and. iplotline) then
     xline(2) = xline(3)
     call pgpt(2,xline(1:2),yline(1:2),imarkerstyle)
  elseif (iplotpts) then
     call pgpt(3,xline,yline,imarkerstyle)
  endif
!
!--add text
!
  if (iplotline .or. iplotpts .and. len_trim(text).gt.0) then
     call pgtext(xline(3) + 0.75*xch,yline(1)-0.25*ych,trim(text))
  endif

  call pgsci(icolourprev)    ! reset colour index
  call pgsls(ilinestyleprev) ! reset line style
  call pgstbg(-1) ! reset text background to transparent

end subroutine legend_markers

end module legends
