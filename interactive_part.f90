!
!--interactive tools on particle plots
!  (experimental at this stage)
!
subroutine interactive_part(npart,iplotx,iploty,irender,xcoords,ycoords,iadvance)
  use limits
  implicit none
  integer, intent(in) :: npart,iplotx,iploty,irender
  integer, intent(out) :: iadvance
  real, dimension(npart), intent(in) :: xcoords, ycoords
  integer :: i,iclosest,nc,ipts,iadvance,int_from_string
  real :: xpt,ypt,xpt2,ypt2,rmin,rr,gradient,yint,xlength
  real :: xlength, ylength
  real, dimension(4) :: xline,yline
  character(len=1) :: char,char2
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
     !
     !--exit if the device is not interactive
     !
     if (char.eq.achar(0)) return
  
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
     case('g','G')   ! draw a line between two points
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
	print*,' next timestep/plot   : space, left click'
	print*,' previous timestep    : right click'
	print*,' jump by n timesteps  : 1,2,3..9 then left or right click'
	print*,' (r)eplot current plot        : r'
        print*,' label closest (p)article     : p'
	print*,' plot a line and find its g)radient : l, L'
	print*,' plot (c)ircle of interaction : c, C'
	print*,' (h)elp                       : h'
        print*,' (q)uit plotting              : q, Q'     	
        print*,'-------------------------------------------------------'
     case('A') ! left click
        !
	!--draw rectangle from the point and reset the limits
	!
	print*,'please select area to zoom in on'
	call pgband(2,1,xpt,ypt,xpt2,ypt2,char2)
	print*,xpt,ypt,xpt2,ypt2,char2
	if (char2.eq.'A') then   ! zoom if another left click
	   lim(iplotx,1) = min(xpt,xpt2)
	   lim(iplotx,2) = max(xpt,xpt2)
	   lim(iploty,1) = min(ypt,ypt2)
	   lim(iploty,2) = max(ypt,ypt2)
	   iadvance = 0
	   iexit = .true.
	endif
     case('-') ! zoom out by 10%
	print*,'zooming out'
	xlength = lim(iplotx,2)-lim(iplotx,1)
	ylength = lim(iploty,2)-lim(iploty,1)
	xlength = 1.1*xlength
	ylength = 1.1*ylength
	lim(iplotx,1) = xpt - 0.5*xlength
	lim(iplotx,2) = xpt + 0.5*xlength
	lim(iploty,1) = ypt - 0.5*ylength
	lim(iploty,2) = ypt + 0.5*ylength
	iadvance = 0
	iexit = .true.
     case('+') ! zoom in by 10%
	print*,'zooming in'
	xlength = lim(iplotx,2)-lim(iplotx,1)
	ylength = lim(iploty,2)-lim(iploty,1)
	xlength = 0.9*xlength
	ylength = 0.9*ylength
	lim(iplotx,1) = xpt - 0.5*xlength
	lim(iplotx,2) = xpt + 0.5*xlength
	lim(iploty,1) = ypt - 0.5*ylength
	lim(iploty,2) = ypt + 0.5*ylength
	iadvance = 0
	iexit = .true.
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
     case(' ') ! space
        iexit = .true.
     case('1','2','3','4','5','6','7','8','9')
        iadvance = int_from_string(char)
	print*,' setting timestep jump = ',iadvance
     end select

  enddo

  return
end subroutine interactive_part
