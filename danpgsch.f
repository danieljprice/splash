c
c  simple subroutine to set the PGPLOT character height
c  in a variety of units, independent of the page size
c
c Arguments:
c  SIZE   (input)  : character height in the appropriate units
c  UNITS  (input)  : Used to specify the units of the output value:
c                    UNITS = 0 : normalized device coordinates
c                    UNITS = 1 : inches
c                    UNITS = 2 : millimeters
c                    UNITS = 3 : pixels
c                    UNITS = 4 : world coordinates
c                    Other values give an error message, and are
c                    treated as 0.
c
c  Daniel Price, Institute of Astronomy, Cambridge, 2004.
c
      SUBROUTINE DANPGSCH(SIZE,UNITS)
      IMPLICIT NONE
      INTEGER UNITS
      REAL SIZE
      REAL XCH, YCH, SCALEFAC, CHARHEIGHT
c
c query current character height in the appropriate units
c
      CALL PGQCH(CHARHEIGHT)
      CALL PGQCS(UNITS,XCH,YCH)
c
c scale the character height appropriately to the desired value
c
      SCALEFAC = SIZE/YCH
c
c reset character height
c
      CALL PGSCH(SCALEFAC*CHARHEIGHT)
      CALL PGQCH(SCALEFAC)
      
      RETURN      
      END SUBROUTINE
