!
!--interactive tools on particle plots
!  allows user to change settings interactively
!
!  Arguments:
!
!   npart   : number of particles plotted
!   iplotx  : quantity plotted as x axis
!   iploty  : quantity plotted as y axis 
!   irender : quantity rendered
!   xcoords(npart) : x coordinates of particles
!   ycoords(npart) : y coordinates of particles
!   xmin, xmax, ymin, ymax : current plot limits
!   iadvance : integer telling the loop how to advance the timestep
!
subroutine interactive_part(npart,iplotx,iploty,irender,xcoords,ycoords, &
  xmin,xmax,ymin,ymax,anglerot,angletilt,iadvance,isave)
  implicit none
  integer, intent(in) :: npart,iplotx,iploty,irender
  integer, intent(out) :: iadvance
  real, dimension(npart), intent(in) :: xcoords, ycoords
  real, intent(inout) :: xmin,xmax,ymin,ymax
  real, intent(inout) :: anglerot,angletilt
  logical, intent(out) :: isave
  integer :: i,iclosest,nc,ipts,int_from_string
  real :: xpt,ypt,xpt2,ypt2,rmin,rr,gradient,yint
  real :: xlength, ylength
  real, dimension(4) :: xline,yline
  character(len=1) :: char,char2
  character(len=20) :: string
  logical :: iexit

  print*,'entering interactive mode...press h in plot window for help'
  char = 'A'
  ipts = 0
  xline = 0.
  yline = 0.
  iexit = .false.
  isave = .false.
  
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
        print*,' select area and zoom : left click (or A)'
        print*,' (r)eplot current plot        : r'
        print*,' label closest (p)article     : p'
        print*,' plot a line and find its g)radient : g'
        print*,' rotate about z axis by +(-) 15 degrees : , (.)'
        print*,' rotate about x axis by +(-) 15 degrees : / ('')' 
        print*,' rotate about z axis by +(-) 30 degrees : < (>)'
        print*,' rotate about z axis by +(-) 30 degrees : ? (")'        
        print*,' next timestep/plot   : space, n'
        print*,' previous timestep    : right click (or X), b'
        print*,' jump by n timesteps  : 0,1,2,3..9 then left or right click'
        print*,' (h)elp                       : h'
        print*,' (s)ave current settings for all steps : s'
        print*,' (q)uit plotting              : q, Q'             
        print*,'-------------------------------------------------------'
     case('s','S')
        isave = .not.isave
        print*,'save settings on exit = ',isave
     !
     !--zoom
     !
     case('A') ! left click
        !
        !--draw rectangle from the point and reset the limits
        !
        print*,'please select area to zoom in on'
        call pgband(2,1,xpt,ypt,xpt2,ypt2,char2)
        print*,xpt,ypt,xpt2,ypt2,char2
        if (char2.eq.'A') then   ! zoom if another left click
           xmin = min(xpt,xpt2)
           xmax = max(xpt,xpt2)
           ymin = min(ypt,ypt2)
           ymax = max(ypt,ypt2)
           iadvance = 0
           iexit = .true.
        endif     
     case('-') ! zoom out by 10%
        print*,'zooming out'
        xlength = xmax - xmin
        ylength = ymax - ymin
        xlength = 1.1*xlength
        ylength = 1.1*ylength
        xmin = xpt - 0.5*xlength
        xmax = xpt + 0.5*xlength
        ymin = ypt - 0.5*ylength
        ymax = ypt + 0.5*ylength
        iadvance = 0
        iexit = .true.
     case('+') ! zoom in by 10%
        print*,'zooming in'
        xlength = xmax - xmin
        ylength = ymax - ymin
        xlength = 0.9*xlength
        ylength = 0.9*ylength
        xmin = xpt - 0.5*xlength
        xmax = xpt + 0.5*xlength
        ymin = ypt - 0.5*ylength
        ymax = ypt + 0.5*ylength
        iadvance = 0
        iexit = .true.
     !
     !--rotation
     !
     case(',') ! rotate left by 15 degrees
        print*,'rotating about z by -15 degrees...'
        anglerot = anglerot - 15.
        iadvance = 0
        iexit = .true.
     case('<') ! rotate left by 30 degrees
        print*,'rotating about z by -30 degrees...'
        anglerot = anglerot - 30.
        iadvance = 0
        iexit = .true.
     case('.') ! rotate right by 15 degrees
        print*,'rotating about z by 15 degrees...'
        anglerot = anglerot + 15.
        iadvance = 0
        iexit = .true.
     case('>') ! rotate left by 30 degrees
        print*,'rotating about z by 30 degrees...'
        anglerot = anglerot + 30.
        iadvance = 0
        iexit = .true.
     case('''') ! rotate up by 15 degrees
        print*,'rotating up (about x) by 15 degrees...'
        angletilt = angletilt - 15.
        iadvance = 0
        iexit = .true.
     case('"') ! rotate up by 30 degrees
        print*,'rotating up (about x) by 30 degrees...'
        angletilt = angletilt - 30.
        iadvance = 0
        iexit = .true.
     case('/') ! rotate down by 15 degrees
        print*,'rotating down (about x) by 15 degrees...'
        angletilt = angletilt + 15.
        iadvance = 0
        iexit = .true.
     case('?') ! rotate down by 30 degrees
        print*,'rotating down (about x) by 30 degrees...'
        angletilt = angletilt + 30.
        iadvance = 0
        iexit = .true.
     !
     !--timestepping
     !
     case('q','Q')
        iadvance = 666666666
        print*,'quitting...'
        iexit = .true.
     case('X','b','B') ! right click -> go back
        iadvance = -abs(iadvance)
        iexit = .true.
     case('r','R') ! replot
        iadvance = 0
        iexit = .true.
     case(' ','n','N') ! space
        iexit = .true.
     case('0','1','2','3','4','5','6','7','8','9')
        iadvance = int_from_string(char)
        print*,' setting timestep jump = ',iadvance
     end select

  enddo
  return
end subroutine interactive_part
