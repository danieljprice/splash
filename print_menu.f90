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
!
!--multiplot
!  
  print 12
  print 18,numplot+1,'multiplot ',nyplotmulti,numplot+2,'set multiplot '
18 format(1x,i2,')',1x,a,'[ ',i2,' ]',5x,i2,') ',a)
  !
  !--supersphplot options 
  !     
  menuitems = 10	! this is the number of options in the menu (set here)

!  if (ishowopts) then

!     print 14,numplot+2, &
!          'set multiplot                    ',nyplotmulti
     print 12
     print 13,numplot+3,'read dat'
     print 14,numplot+4, &
          'change number of timesteps read  ',(n_end-nstart+1)/nfreq 
     print 15,numplot+5, &
          'toggle interactive mode          ',interactive
     print 13,numplot+6, &
          'page options'
     print 13,numplot+7,'particle plot options'
     print 14,numplot+8, &
          'rendering/vector plot options    ',ivecplot_nomulti
     print 13,numplot+9,'change plot limits'

!  endif	! show/hide opts

  print 12
!  if (ishowopts) then
!     print 13,numplot+10,'hide options'      
!  else
!     print 13,numplot+10,'supersphplot options'
!  endif
  print 11,numplot+10,'save defaults       ',numplot+11,'exit supersphplot'
  print 12

11 format(1x,i2,')',1x,a20,1x,i2,')',1x,a)
12 format(1x,50('-'))
13 format(1x,i2,')',1x,a)
14 format(1x,i2,')',1x,a,'( ',i2, ' )')
15 format(1x,i2,')',1x,a,'( ',L1,' )')
  !
  !--prompt for selection
  !

9901 continue
  if (ipicky.eq.0) ipicky = ndim+1
  write(*,"(a)",ADVANCE='NO') 'Please enter your selection now (y axis or option):'
  read(*,*,ERR=9901) ipicky
  !call prompt('Please enter your selection now (y axis or option): ',ipicky)
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
