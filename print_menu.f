!!
!!  Prints supersphplot menu, returns ipicky,ipickx
!!
!!
      SUBROUTINE print_menu(ipicky,ipickx)
      USE labels
      USE multiplot
      USE settings
      USE prompting
      IMPLICIT NONE
      INTEGER :: i,ihalf,iadjust
      INTEGER, INTENT(OUT) :: ipickx,ipicky
      CHARACTER(LEN=25) :: transform_label      
!
!--print labels of data columns
!        
      PRINT*,' You may choose from a delectable sample of plots '
      PRINT 12
      ihalf = numplot/2		! print in two columns
      iadjust = MOD(numplot,2)
      PRINT 11, (i,transform_label(label(i),itrans(i)),
     &             ihalf + i + iadjust,
     &             transform_label(label(ihalf + i +
     &              iadjust),itrans(ihalf+i+iadjust)),i=1,ihalf)
      IF (iadjust.NE.0) THEN
         PRINT 13, ihalf + iadjust,transform_label(label(ihalf +
     &	 iadjust),itrans(ihalf+iadjust))
      ENDIF
      PRINT 12
      PRINT 18,numplot+1,
     &        'multiplot ',
     &         multiploty(1:nyplotmulti)
!
!--supersphplot options 
!     
      menuitems = 23	! this is the number of options in the menu (set here)

      IF (ishowopts) THEN
     
      PRINT 12
      PRINT 13,numplot+2,'read dat'
      PRINT 14,numplot+3,
     &        'Change number of timesteps read  ',(n_end-nstart+1)/nfreq 
      PRINT 15,numplot+4,
     &        'toggle animate                   ',animate
      PRINT 15,numplot+5,
     &        'toggle adaptive/fixed limits     ',iadapt
      PRINT 15,numplot+6,
     &        'toggle mag field                 ',magfield
      PRINT 15,numplot+7,
     &        'toggle cross section/projection  ',xsec_nomulti
      PRINT 15,numplot+8,
     &        'toggle circles of interaction    ',plotcirc
      PRINT 16,numplot+9,
     &        'Change graph markers             ',imark,imarkg
      PRINT 16,numplot+10,
     &        'Change plots per page            ',nacross,ndown
      PRINT 14,numplot+11,
     &        'Set multiplot                    ',nyplotmulti
      PRINT 15,numplot+12,
     &        'toggle page change               ',ipagechange
      PRINT 21,numplot+13,
     &        'change paper size                ',papersizex,aspectratio
      PRINT 19,numplot+14,
     &        'toggle plot line (all, init)     ',iplotline,iplotlinein 
      PRINT 14,numplot+15,
     &        'toggle exact solution            ',iexact
      PRINT 20,numplot+16,
     &        'toggle plot average line         ',iplotav,nbins
      PRINT 15,numplot+17,
     &        'toggle label particles           ',ilabelpart
      PRINT 19,numplot+18,
     &        'toggle plot ghosts/sinks         ',iplotghost,iplotsink
      PRINT 16,numplot+19,
     &        'set rendering/vector plots       ',
     &         irender,ivecplot_nomulti
      PRINT 13,numplot+20,
     &        'apply transformations (log10,1/x)'
      IF (iadapt) THEN
         PRINT 17,numplot+21,
     &        'Zoom out/in                      ',scalemax      
      ELSE
         PRINT 17,numplot+21,
     &        'Zoom out/in/set manual limits    ',zoom
      ENDIF
      
      ENDIF	! show/hide opts
      
      PRINT 12
      IF (ishowopts) THEN
         PRINT 13,numplot+22,'hide options'      
      ELSE
         PRINT 13,numplot+22,'supersphplot options'
      ENDIF
      PRINT 13,numplot+23,'Save defaults'
      PRINT 13,numplot+24,'Exit supersphplot'
      PRINT 12
     
11    FORMAT(1x,i2,')',1x,a20,1x,i2,')',1x,a)
12    FORMAT(1x,45('-'))
13    FORMAT(1x,i2,')',1x,a)
14    FORMAT(1x,i2,')',1x,a,'( ',i2, ' )')
15    FORMAT(1x,i2,')',1x,a,'( ',L1,' )')
16    FORMAT(1x,i2,')',1x,a,'( ',i2,',',i2,' )')
17    FORMAT(1x,i2,')',1x,a,'( ',f4.2,' )')
18    FORMAT(1x,i2,')',1x,a,'( ',10(i2,1x),' )')
19    FORMAT(1x,i2,')',1x,a,'( ',L1,',',1x,L1,' )')
20    FORMAT(1x,i2,')',1x,a,'( ',L1,',',i2,' )')
21    FORMAT(1x,i2,')',1x,a,'( ',f5.2,',',1x,f5.2,' )')
!
!--prompt for selection
!

9901  CONTINUE
      PRINT*,'Please enter your selection now (y axis or option): '
      READ (*,*,ERR=9901) ipicky
!
!--if needed prompt for x axis selection
!
      IF ((ipicky.le.(numplot-nextra)).and.(ipicky.gt.0)) THEN
         CALL prompt(' (x axis) ',ipickx)
	 IF (ipickx.GT.numplot .OR. ipickx.LE.0) GOTO 9901
      ELSEIF (ipicky.GT.numplot+menuitems .OR. ipicky.LE.0) THEN
         STOP
      ENDIF   

      RETURN
      END SUBROUTINE print_menu
