!-----------------------------------------------------------------
!     plots time on plot
!     arguments:
!           t : current time
!        hpos : horizontal position as fraction of viewport
!        vpos : vertical position in character heights from top
!-----------------------------------------------------------------

      SUBROUTINE legend(t,hpos,vpos)
      IMPLICIT NONE    
      INTEGER MM,PP,NC,ndecimal,ndec
      REAL hpos,vpos,t,tplot
      REAL xch, ych
      CHARACTER string*15
      
      !!CALL PGQVP(0,x1,x2,y1,y2)
      CALL PGQCS(0,xch,ych) ! get character height in normalised device coords
      
!      vspace=1.2
      ndecimal = 2	! number of decimal places to display
      ndec = 10**ndecimal
      IF (t.eq.0.0) THEN
         tplot = 1E-6
      ELSE
         tplot = t    !/(2.*3.1415926536)
      ENDIF
      MM=nint(tplot*ndec)
      PP=nint(log10(tplot)-log10(tplot*ndec))
      CALL PGNUMB(MM,PP,1,STRING,NC)
      CALL PGMTEXT('T',-vpos,hpos,0.0,'t='//STRING(1:NC))
   !!   CALL PGMTEXT('T',vpos,hpos,0.5,
   !!  &             STRING(1:NC)//' Rotational periods')

      RETURN
      END
