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

end module legends
