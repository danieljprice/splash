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

subroutine legend(t,hpos,vpos)
 implicit none    
 integer mm,pp,nc,ndecimal,ndec
 real hpos,vpos,t,tplot
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
 call pgmtext('t',-vpos,hpos,0.0,'t='//string(1:nc))

 return
end subroutine legend

!-----------------------------------------------------------------
!     plots vector plot legend
!     arguments:
!           t : current time
!        hpos : horizontal position as fraction of viewport
!        vpos : vertical position in character heights from top
!-----------------------------------------------------------------

subroutine legend_vec(scale,label)
 use settings_vecplot, only: hposlegendvec, vposlegendvec
 implicit none
 real, intent(in) :: scale
 character(len=*), intent(in) :: label
 real :: xmin,xmax,ymin,ymax
 real :: xch,ych
 real :: xpos,ypos,xpos2,ypos2,rarrow
 character(len=len(label)+10) :: string
!
!--convert hpos and vpos to x, y to plot arrow
!
 call pgqwin(xmin,xmax,ymin,ymax)
 print*,'lims = ',xmin,xmax,ymin,ymax
 xpos = xmin + hposlegendvec*(xmax-xmin)
 call pgqcs(0,xch,ych) 
 ypos = ymax - vposlegendvec*ych
!
!--plot arrow diagonally
! 
 rarrow = scale
 xpos2 = xpos + sqrt(0.5)*rarrow
 ypos2 = ypos + sqrt(0.5)*rarrow
 
 write(string,"(a,'=',f6.2)") trim(label),scale
 
 print*,'legendvec: xpos, ypos = ',xpos,ypos
 call pgarro(xpos,ypos,xpos2,ypos2)
 call pgmtext('t',-vposlegendvec,hposlegendvec+0.02,0.0,trim(string))

 return
end subroutine legend_vec

end module legends
