!
!--interactive tools on particle plots
!  (experimental at this stage)
!
subroutine interactive_part(npart,xcoords,ycoords,iadvance)
  implicit none
  integer, intent(in) :: npart
  integer, intent(out) :: iadvance
  real, dimension(npart), intent(in) :: xcoords, ycoords
  integer :: i,iclosest,nc,ipts,iadvance
  real :: xpt,ypt,rmin,rr,gradient,yint,xlength
  real, dimension(4) :: xline,yline
  character(len=1) :: char
  character(len=20) :: string
  logical :: iexit

  print*,'entering interactive mode...press h for help'
  char = 'A'
  ipts = 0
  xline = 0.
  yline = 0.
  iexit = .false.

  do while (.not.iexit)
     call pgcurs(xpt,ypt,char)
     
     print*,'location: x, y = ',xpt,ypt,' function = ',char
     !
     !--find closest particle
     !  
     rmin = 1.e6
     do i=1,npart
        rr = (xcoords(i)-xpt)**2 + (ycoords(i)-ypt)**2
        if (rr.lt.rmin) then
           iclosest = i
           rmin = rr
        endif
     enddo
     
     select case(char)
     case('p')
        print*,' closest particle = ',iclosest,'x = ',xcoords(iclosest),' y =',ycoords(iclosest)
        call pgnumb(iclosest,0,1,string,nc)
        call pgsch(2.0)
        call pgtext(xcoords(iclosest),ycoords(iclosest),string(1:nc))
        call pgsch(1.0)
     case('c','C')
        print*,'plotting circle of interaction (not implemented)'
     case('l','L')   ! draw a line between two points
        ipts = ipts + 1
        xline(2) = xline(3)
	yline(2) = yline(3)     
        xline(3) = xcoords(iclosest)
	yline(3) = ycoords(iclosest)
        call pgpt(1,xline(3),yline(3),4)
	if (ipts.gt.1) then
	   gradient = (yline(3)-yline(2))/(xline(3)-xline(2))
	   yint = yline(3) - gradient*xline(3)
	   xlength = sqrt((xline(3)-xline(2))**2 + (yline(3)-yline(2))**2) 
	   print*,' gradient = ',gradient,' y intercept = ',yint, 'length = ',xlength
	   if (xline(2).lt.xline(3)) then 
	      xline(1) = minval(xcoords)
	      xline(4) = maxval(xcoords)
	   else
	      xline(1) = maxval(xcoords)
	      xline(4) = minval(xcoords)
	   endif
	   yline(1) = gradient*xline(1) + yint
	   yline(4) = gradient*xline(4) + yint
	   call pgline(4,xline,yline) 
	else
	   print*,' select another point and press l again'
	endif
     case('h')
        print*,'-------------- interactive mode commands --------------'
	print*,' advance forward   : space, left click'
	print*,' go backwards      : right click'
	print*,' (r)eplot current plot        : r'
        print*,' label closest (p)article     : p'
	print*,' plot connecting (l)ine       : l, L'
	print*,' plot (c)ircle of interaction : c, C'
	print*,' (h)elp                       : h'
        print*,' (q)uit plotting              : q, Q'     	
        print*,'-------------------------------------------------------'
     case('q','Q')
        iadvance = 666666666
	print*,'quitting...'
	iexit = .true.
     case('X') ! right click -> go back
        iadvance = -abs(iadvance)
	iexit = .true.
     case('r','R') ! replot
        iadvance = 0
	iexit = .true.
     case(' ','A') ! space, left click
        iexit = .true.
     end select

  enddo

  return
end subroutine interactive_part
