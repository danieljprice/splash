!---------------------------------------------------
!     plots time on plot
!---------------------------------------------------

      SUBROUTINE legend(t)
      IMPLICIT NONE    
      INTEGER MM,PP,NC,ndecimal,ndec
      REAL hpos,vpos,t
      CHARACTER string*15
      
      hpos=0.75	! either 0.1 or 0.75 for 1D shocks
!      vspace=1.2
      vpos=-2.5	    
      ndecimal = 2	! number of decimal places to display
      ndec = 10**ndecimal
      IF (t.eq.0.0) t=1E-6
      MM=nint(t*ndec)
      PP=nint(log10(t)-log10(t*ndec))
      CALL PGNUMB(MM,PP,1,STRING,NC)
      CALL PGMTEXT('T',vpos,hpos,0.0,'t='//STRING(1:NC))

      RETURN
      END
