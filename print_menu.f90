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
  menuitems = 14	! this is the number of options in the menu (set here)

  if (ishowopts) then

     print 12
     print 13,numplot+2,'read dat'
     print 14,numplot+3, &
          'change number of timesteps read  ',(n_end-nstart+1)/nfreq 
     print 15,numplot+4, &
          'toggle interactive mode          ',interactive
     print 15,numplot+5, &
          'toggle animate                   ',animate
     print 14,numplot+6, &
          'set multiplot                    ',nyplotmulti
     print 15,numplot+7, &
          'toggle cross section/projection  ',xsec_nomulti 
     print 14,numplot+8, &
          'toggle exact solution            ',iexact
     print 13,numplot+9, &
          'page options'
     print 13,numplot+10,'particle plot options'
     print 14,numplot+11, &
          'rendering/vector plot options    ',ivecplot_nomulti
     print 13,numplot+12,'plot limits'

  endif	! show/hide opts

  print 12
  if (ishowopts) then
     print 13,numplot+13,'hide options'      
  else
     print 13,numplot+13,'supersphplot options'
  endif
  print 13,numplot+14,'Save defaults'
  print 13,numplot+15,'Exit supersphplot'
  print 12

11 format(1x,i2,')',1x,a20,1x,i2,')',1x,a)
12 format(1x,45('-'))
13 format(1x,i2,')',1x,a)
14 format(1x,i2,')',1x,a,'( ',i2, ' )')
15 format(1x,i2,')',1x,a,'( ',L1,' )')
18 format(1x,i2,')',1x,a,'( ',10(i2,1x),' )')
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
