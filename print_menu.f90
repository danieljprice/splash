!!
!!  Prints supersphplot menu, returns ipicky,ipickx
!!
!!
subroutine print_menu(ipicky,ipickx,irender)
  use labels
  use multiplot
  use settings
  use prompting
  implicit none
  integer :: i,ihalf,iadjust
  integer, intent(INOUT) :: ipickx,ipicky,irender
  character(LEN=25) :: transform_label      
  !
  !--print labels of data columns
  !        
  print*,' You may choose from a delectable sample of plots '
  print 12
  ihalf = numplot/2		! print in two columns
  iadjust = mod(numplot,2)
  print 11, (i,transform_label(label(i),itrans(i)), &
       ihalf + i + iadjust, transform_label(label(ihalf + i + iadjust), &
       itrans(ihalf+i+iadjust)),i=1,ihalf)
  if (iadjust.ne.0) then
     print 13, ihalf + iadjust,transform_label(label(ihalf + iadjust), &
          itrans(ihalf+iadjust))
  endif
  print 12
  print 18,numplot+1,'multiplot ',multiploty(1:nyplotmulti)
  !
  !--supersphplot options 
  !     
  menuitems = 23	! this is the number of options in the menu (set here)

  if (ishowopts) then

     print 12
     print 13,numplot+2,'read dat'
     print 14,numplot+3, &
          'Change number of timesteps read  ',(n_end-nstart+1)/nfreq 
     print 15,numplot+4, &
          'toggle interactive mode          ',interactive
     print 15,numplot+5, &
          'toggle axes                      ',axes
     print 15,numplot+6, &
          'toggle adaptive/fixed limits     ',iadapt
     print 15,numplot+7, &
          'toggle cross section/projection  ',xsec_nomulti
     print 15,numplot+8, &
          'toggle circles of interaction    ',plotcirc
     print 16,numplot+9, &
          'Change graph markers             ',imark,imarkg
     print 16,numplot+10, &
          'Change plots per page            ',nacross,ndown
     print 14,numplot+11, &
          'Set multiplot                    ',nyplotmulti
     print 15,numplot+12, &
          'toggle page change               ',ipagechange
     print 21,numplot+13, &
          'change paper size                ',papersizex,aspectratio
     print 19,numplot+14, &
          'toggle plot line (all, init)     ',iplotline,iplotlinein 
     print 14,numplot+15, &
          'toggle exact solution            ',iexact
     print 20,numplot+16, &
          'toggle plot average line         ',iplotav,nbins
     print 15,numplot+17, &
          'toggle label particles           ',ilabelpart
     print 19,numplot+18, &
          'toggle plot ghosts/sinks         ',iplotghost,iplotsink
     print 14,numplot+19, &
          'rendering/vector plot options    ',ivecplot_nomulti
     print 13,numplot+20, &
          'apply transformations (log10,1/x)'
     if (iadapt) then
        print 17,numplot+21, &
             'Zoom out/in                      ',scalemax      
     else
        print 17,numplot+21, &
             'Zoom out/in/set manual limits    ',zoom
     endif

  endif	! show/hide opts

  print 12
  if (ishowopts) then
     print 13,numplot+22,'hide options'      
  else
     print 13,numplot+22,'supersphplot options'
  endif
  print 13,numplot+23,'Save defaults'
  print 13,numplot+24,'Exit supersphplot'
  print 12

11 format(1x,i2,')',1x,a20,1x,i2,')',1x,a)
12 format(1x,45('-'))
13 format(1x,i2,')',1x,a)
14 format(1x,i2,')',1x,a,'( ',i2, ' )')
15 format(1x,i2,')',1x,a,'( ',L1,' )')
16 format(1x,i2,')',1x,a,'( ',i2,',',i2,' )')
17 format(1x,i2,')',1x,a,'( ',f4.2,' )')
18 format(1x,i2,')',1x,a,'( ',10(i2,1x),' )')
19 format(1x,i2,')',1x,a,'( ',L1,',',1x,L1,' )')
20 format(1x,i2,')',1x,a,'( ',L1,',',i2,' )')
21 format(1x,i2,')',1x,a,'( ',f5.2,',',1x,f5.2,' )')
  !
  !--prompt for selection
  !

9901 continue
  print*,'Please enter your selection now (y axis or option): '
  read (*,*,ERR=9901) ipicky
  !
  !--if needed prompt for x axis selection
  !
  if ((ipicky.le.(numplot-nextra)).and.(ipicky.gt.0)) then
     call prompt(' (x axis) ',ipickx)
     if (ipickx.gt.numplot .or. ipickx.le.0) then
        goto 9901
     elseif (ipicky.le.ndim .and. ipickx.le.ndim) then
        call prompt('(render) (0=none)',irender,0,numplot)
     else
        irender = 0
     endif

  elseif (ipicky.gt.numplot+menuitems .or. ipicky.le.0) then
     stop
  endif

  return
end subroutine print_menu
