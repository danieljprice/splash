! ------------------------------------------------------------------------
!      defines colour schemes for rendering
!      ** add your own here **
! ------------------------------------------------------------------------     
subroutine colour_set(icolours)
  implicit none
  integer, intent(in) :: icolours
  integer :: i,icolourmin,icolmin,icolmax,ncolmax
  integer :: ncolours
  real :: red,green,blue,hue,light,sat
  real :: dhue,dlight,dsat,dred,dgreen,dblue
  logical :: rgb
  character(len=30) :: scheme
      
  red = 0.0
  green = 0.0
  blue = 0.0
  rgb = .false.
!
!--set number of colours
!      
  ncolours = 256
!
!--set first colour index (warning: colours 1-16 have presets, so
!  overwriting these means that line graphs that use colour will come
!  out funny). Best to leave 0 and 1 alone as these are black and white.
!
  icolourmin = 2
!
!--inquire as to colour range available on current device
!  adjust ncolours if necessary
!      
  call PGQCOL(icolmin,icolmax)
  if (icolourmin.lt.icolmin) icolourmin = icolmin
  ncolmax = icolmax - icolourmin
  if (ncolours.gt.ncolmax) then  
     ncolours = ncolmax
     print*,'Warning: Device allows only ',ncolours+1,' colours'
  endif
!
!--starting values
!      
  select case(icolours)
  case(2)  
     scheme = 'red'
     rgb = .false.
     hue=100
     light=0.0
     sat=1.0      
     dlight = 1.0/FLOAT(ncolours)
     dsat = 0.0
     dhue = 80.0/FLOAT(ncolours)	 
  case(3)	! blue
     scheme = 'ice blue'
     rgb = .false.
     hue=330
     light=0.0
     sat=1.0
     dlight = 1.0/FLOAT(ncolours)
     dsat = 0.0
     dhue = 40.0/FLOAT(ncolours)
  case(4) 
     scheme = 'rainbow'
     rgb = .false.         
     hue=100
     light=0.5
     sat=1.0
     dlight = 0.0	!/FLOAT(ncolours)
     dsat = 0.0
     dhue = 320.0/FLOAT(ncolours)
  case(5)	! green
     scheme = 'frog monster'
     rgb = .false.         
     hue=200
     light=0.0
     sat=1.0
     dlight = 1.0/FLOAT(ncolours)
     dsat = 0.0
     dhue = 80.0/FLOAT(ncolours)
  case(6)
     scheme = 'some other crap'
     rgb = .false.
     hue=100
     light=0.5
     sat= 1.    
     dlight = 0.  !!!1.0/FLOAT(ncolours)
     dsat = 0.	!1.0/FLOAT(ncolours)
     dhue = 500.0/FLOAT(ncolours)	 
  case(10)	! rgb attempts
     scheme = 'really useless'
     rgb = .true.
     red = 0.3333
     green = 0.0
     blue = 0.0      
     dred = 0.0	!0.333/FLOAT(ncolours)
     dblue = 0.3333/FLOAT(ncolours)
     dgreen = 0.0	!0.333/FLOAT(ncolours)	 
  case default
     scheme = 'none'
     return
  end select
!
!--set colour indexes from icolourmin-> icolourmin+ncolours
!  increment values in steps
!
  do i=icolourmin,icolourmin+ncolours
     if (rgb) then
        red = red + dred
        blue = blue + dblue
        green = green + dgreen
        call PGSCR(i,red,green,blue)	 
     else   
        hue = hue + dhue
        sat = sat + dsat
        light = light + dlight
        call PGSHLS(i,hue,light,sat)
     endif
  enddo
  !
  !--set this as the range of colour indices to use
  !
  call PGSCIR(icolourmin,icolourmin+ncolours)  
  print*,'using colour scheme ',trim(scheme)
  return
  
end subroutine colour_set
