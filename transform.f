!------------------------------------------------------------------------
!
! module containing global parameter for transform subroutines
! this is the maximum number of transformations you can do in a row      
!
!------------------------------------------------------------------------
      MODULE transform_params
         INTEGER, PARAMETER :: nmax = 10
      END MODULE transform_params

!------------------------------------------------------------------------
!
!  subroutine returns log, 1/x of a given array
!
!  * can specify up to 9 individual operations to perform
!  * combinations of transformations are done when the 
!    input number is > 10 (e.g. 321 means 1/x, then abs, then log10)
!
!------------------------------------------------------------------------
      SUBROUTINE transform(array,arrayout,itrans,isize)
      USE transform_params
      IMPLICIT NONE
      INTEGER :: i,ndigits,itransmulti
      INTEGER, INTENT(IN) :: isize,itrans
      REAL, DIMENSION(isize), INTENT(IN) :: array
      REAL, DIMENSION(isize), INTENT(OUT) :: arrayout
      REAL, DIMENSION(isize) :: arraytemp
      INTEGER, DIMENSION(nmax) :: digit     
!
!--extract the digits from the input number 
!            
      IF (itrans.GT.0) THEN      
         CALL get_digits(itrans,digit,ndigits,nmax)         
!
!--do a transformation for each digit     
!
         arraytemp = array
	 
	 DO i=1,ndigits
	    itransmulti = digit(i)
!
!--perform transformation appropriate to this digit     
!
            SELECT CASE(itransmulti)
	       CASE(1)
                  WHERE (arraytemp > 0)
	             arraytemp = LOG10(arraytemp)
	          ELSEWHERE
	             arraytemp = 0.
     	          END WHERE
	       CASE(2)
                  arraytemp = ABS(arraytemp)    
	       CASE(3)
	          WHERE (arraytemp .NE. 0)
	             arraytemp = 1./arraytemp
	          ELSEWHERE
	             arraytemp = 0.
  	          END WHERE
	       CASE(4) 
	          WHERE (arraytemp .GT. 0)
	             arraytemp = SQRT(arraytemp)
		  ELSEWHERE
		     arraytemp = 0.
		  END WHERE
	       CASE(5)
	          arraytemp = arraytemp**2	     
            END SELECT
         ENDDO
	 
	 arrayout = arraytemp
	 
      ELSE
         arrayout = array
      ENDIF 
      
      END SUBROUTINE transform
      
!------------------------------------------------------------------------
!
!  same as transform but for a two dimensional array
!
!------------------------------------------------------------------------
      SUBROUTINE transform2(array,arrayout,itrans,isizex,isizey)
      USE transform_params
      IMPLICIT NONE
      INTEGER :: i,ndigits,itransmulti
      INTEGER, INTENT(IN) :: itrans,isizex,isizey
      REAL, DIMENSION(isizex,isizey), INTENT(IN) :: array
      REAL, DIMENSION(isizex,isizey), INTENT(OUT) :: arrayout
      REAL, DIMENSION(isizex,isizey) :: arraytemp
      INTEGER, DIMENSION(nmax) :: digit     
!
!--extract the digits from the input number 
!            
      IF (itrans.GT.0) THEN      
         CALL get_digits(itrans,digit,ndigits,nmax)         
!
!--do a transformation for each digit     
!
         arraytemp = array
	 
	 DO i=1,ndigits
	    itransmulti = digit(i)
!
!--perform transformation appropriate to this digit     
!
            SELECT CASE(itransmulti)
	       CASE(1)
                  WHERE (arraytemp > 0)
	             arraytemp = LOG10(arraytemp)
	          ELSEWHERE
	             arraytemp = 0.
	          END WHERE
	      CASE(2)
	         arraytemp = ABS(arraytemp)    
	      CASE(3)
	         WHERE (arraytemp .NE. 0)
	            arraytemp = 1./arraytemp
	         ELSEWHERE
	            arraytemp = 0.
	         END WHERE
	      CASE(4) 
	         WHERE (arraytemp .GT. 0)
	            arraytemp = SQRT(arraytemp)
	         ELSEWHERE
	            arraytemp = 0.
	         END WHERE
	      CASE(5)
	         arraytemp = arraytemp**2	 
           END SELECT
         ENDDO
	
	 arrayout = arraytemp
	 
      ELSE
         arrayout = array      
      ENDIF
      
      END SUBROUTINE transform2
      
      
!------------------------------------------------------------------------
!
!  function to adjust the label of a plot if log, 1/x etc
!
!  Note: *cannot* put print or write statements into this function
!        as it is used in the middle of write or print statements
!      
!------------------------------------------------------------------------
      FUNCTION transform_label(label,itrans)
      USE transform_params
      IMPLICIT NONE
      INTEGER :: itrans,itransmulti,i,itransprev,ndigits
      INTEGER, DIMENSION(nmax) :: digit
      CHARACTER(LEN=*) :: label, transform_label
      CHARACTER(LEN=120) :: temp_label      
!
!--extract the digits from the input number 
!            
      IF (itrans.GT.0) THEN      
         CALL get_digits(itrans,digit,ndigits,nmax)         
         temp_label = label      
!
!--do a transformation for each digit     
!
         DO i=1,ndigits
	    itransmulti = digit(i)
!
!--perform transformation appropriate to this digit     
!
            SELECT CASE(itransmulti)
	       CASE(1)
	          temp_label = 'log\d10\u'//TRIM(temp_label)
	       CASE(2)
	          temp_label = '|'//TRIM(temp_label)//'|'
	       CASE(3)
                  temp_label = '1/'//TRIM(temp_label)
	       CASE(4)
	          temp_label = 'SQRT('//TRIM(temp_label)//')'
	       CASE(5)
	          temp_label = TRIM(temp_label)//'\u2\d'
	       CASE DEFAULT
                  temp_label = TRIM(temp_label)
            END SELECT
         ENDDO 
             
         transform_label = temp_label
      ELSE
         transform_label = label     
      ENDIF     
      
      END FUNCTION transform_label

!------------------------------------------------------------------------
!     get_digits: for an integer i returns number of digits it contains
!     and a list of these
!
!     i            : integer to split into digits
!     nmax	   : dimensions of digits array
!     digits(nmax) : array of digits
!     ndigits      : number of digits in i
!------------------------------------------------------------------------
      
      SUBROUTINE get_digits(i,digits,ndigits,nmax)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i,nmax
      INTEGER, INTENT(OUT) :: ndigits
      INTEGER, INTENT(OUT), DIMENSION(nmax) :: digits
      INTEGER :: j,isubtract,idigit
      
      ndigits = 0
      isubtract = 0      

      DO j=nmax,0,-1
         IF (i.GE.10**j) THEN 
	    ndigits = ndigits + 1
	    idigit = (i - isubtract)/10**j
	    digits(ndigits) = idigit
	    isubtract = isubtract + digits(ndigits)*10**j
	 ENDIF   
      ENDDO
      
      END SUBROUTINE get_digits
