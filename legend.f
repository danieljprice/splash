!---------------------------------------------------
!     plots time on plot
!---------------------------------------------------

      SUBROUTINE legend(t)
      IMPLICIT NONE    
      INTEGER MM,PP,NC,ndecimal,ndec
      REAL hpos,vpos,t,tplot
      CHARACTER string*15
      
      !hpos=0.5	! either 0.1 or 0.75 for 1D shocks
      hpos=0.85
!      vspace=1.2
      !vpos=0.5	    
      vpos=-2.0
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
      CALL PGMTEXT('T',vpos,hpos,0.5,'t='//STRING(1:NC))
   !!   CALL PGMTEXT('T',vpos,hpos,0.5,
   !!  &             STRING(1:NC)//' Rotational periods')

      RETURN
      END
