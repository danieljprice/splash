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
 use settings_page, only:hposlegend, vposlegend
 implicit none
 real, intent(in) :: t    
 integer :: mm,pp,nc,ndecimal,ndec
 real :: tplot
 character string*15

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
 call pgmtext('t',-vposlegend,hposlegend,0.0,'t='//string(1:nc))

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
 real :: xpos,ypos,xpos2,ypos2,rarrow
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

end module legends
