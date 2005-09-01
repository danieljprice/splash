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

subroutine legend_vec(vecmax,scale,label,charheight)
 use settings_vecplot, only: hposlegendvec, vposlegendvec
 implicit none
 real, intent(in) :: vecmax,scale,charheight
 character(len=*), intent(in) :: label
 integer, parameter :: npixx=1, npixy=5
 integer :: j
 real, dimension(npixx,npixy) :: vecx,vecy
 real :: xmin,xmax,ymin,ymax
 real :: xch,ych,dx
 real :: xpos,ypos !!,xpos2,ypos2,rarrow
 character(len=len(label)+10) :: string
 real :: trans(6)
 
!
!--convert hpos and vpos to x, y to plot arrow
!
 call pgqwin(xmin,xmax,ymin,ymax)
!! print*,'lims = ',xmin,xmax,ymin,ymax
 xpos = xmin + hposlegendvec*(xmax-xmin)
 call pgqcs(0,xch,ych) 
 ypos = ymin - vposlegendvec*ych

!! print*,'legendvec: xpos, ypos = ',xpos,ypos
 
 dx = ych
 
 do j=1,npixy
    vecx(npixx,j) = (npixy-j)*vecmax
    vecy(npixx,j) = vecx(npixx,j)
 enddo

!set up grid for rendering 

 trans(1) = xpos - 0.5*dx                ! this is for the pgvect call
 trans(2) = dx
 trans(3) = 0.0
 trans(4) = ypos - 0.5*dx
 trans(5) = 0.0
 trans(6) = dx
!
!--this should match the call in render_vec
!
! call pgvect(vecx(:,:),vecy(:,:),npixx,npixy,1,npixx,1,npixy, &
!      scale,0,trans,-1000.0)
!
!--plot a box around the legend
! 
!! call pgbox(xpos,ypos,xpos2,ypos2)
!! write(string,"(a,'=',f6.2)") trim(label),scale
 

! call pgarro(xpos,ypos,xpos2,ypos2)
! call pgmtext('t',-vposlegendvec,hposlegendvec+0.02,0.0,trim(string))

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
