module interactive_routines
 implicit none
 public :: interactive_part,interactive_step
 private :: mvlegend,mvtitle,save_limits,save_rotation
 
contains
!
!--interactive tools on particle plots
!  allows user to change settings interactively
!
!  Arguments:
!
!  INPUT:
!   npart   : number of particles plotted
!   iplotx  : quantity plotted as x axis
!   iploty  : quantity plotted as y axis 
!   irender : quantity rendered
!   xcoords(npart) : x coordinates of particles
!   ycoords(npart) : y coordinates of particles
!   hi(npart)      : smoothing lengths of particles
! CHANGEABLE:
!   icolourpart(npart) : flag indicating colour of particles
!   xmin, xmax, ymin, ymax : current plot limits
! OUTPUT:
!   iadvance : integer telling the loop how to advance the timestep
!   isave    : integer telling the loop to save the settings
!
subroutine interactive_part(npart,iplotx,iploty,irender,xcoords,ycoords,hi, &
  icolourpart,xmin,xmax,ymin,ymax,rendermin,rendermax, &
  anglex,angley,anglez,ndim,iadvance,isave)
  implicit none
  integer, intent(in) :: npart,iplotx,iploty,irender,ndim
  integer, intent(out) :: iadvance
  integer, dimension(npart), intent(inout) :: icolourpart
  real, dimension(npart), intent(in) :: xcoords, ycoords, hi
  real, intent(inout) :: xmin,xmax,ymin,ymax,rendermin,rendermax
  real, intent(inout) :: anglex,angley,anglez
  logical, intent(out) :: isave
  integer :: i,iclosest,nc,ipts,ierr
  integer :: nmarked
  real :: xpt,ypt,xpt2,ypt2,xptmin,xptmax,yptmin,yptmax
  real :: rmin,rr,gradient,yint
  real :: xlength, ylength, drender
  real, dimension(4) :: xline,yline
  character(len=1) :: char,char2
  character(len=20) :: string
  logical :: iexit, rotation

  call pgqinf('CURSOR',string,nc)
  if (string(1:nc).eq.'YES') then
     print*,'entering interactive mode...press h in plot window for help'
  else
     print*,'cannot enter interactive mode: device has no cursor'
     return
  endif
  char = 'A'
  ipts = 0
  xline = 0.
  yline = 0.
  xpt = 0.
  ypt = 0.
  xpt2 = 0.
  ypt2 = 0.
  nc = 0
  iexit = .false.
  isave = .false.
  rotation = .false.
  if (iplotx.le.ndim .and. iploty.le.ndim .and. ndim.ge.2) rotation = .true.
  
  do while (.not.iexit)
     call pgcurs(xpt,ypt,char)
     !
     !--exit if the device is not interactive
     !
     if (char.eq.achar(0)) return
  
     print*,' x, y = ',xpt,ypt,' function = ',char
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
     !
     !--particle plot stuff
     !
     case('p')
        print*,' closest particle = ',iclosest,'x = ',xcoords(iclosest),' y =',ycoords(iclosest)
        call pgnumb(iclosest,0,1,string,nc)
        call pgsch(2.0)
        call pgtext(xcoords(iclosest),ycoords(iclosest),string(1:nc))
        call pgsch(1.0)
     case('c','C')
        print*,'plotting circle of interaction on particle ',iclosest
        call pgsfs(2)
        call pgcirc(xcoords(iclosest),ycoords(iclosest),2.*hi(iclosest))
     case('g')   ! draw a line between two points
        xline(2) = xpt
        yline(2) = ypt
        !--mark first point
        call pgpt1(xpt,ypt,4)
        !--select second point
        print*,' select another point (using left click or g) to plot line '
        call pgband(1,1,xline(2),yline(2),xline(3),yline(3),char2)
        !--draw line if left click or g
        select case(char2)
        case('A','g')
           !--mark second point
           call pgpt1(xline(3),yline(3),4)
           xlength = xline(3)-xline(2)
           if (abs(xlength).lt.tiny(xlength)) then
              xline(1) = xline(2)
              xline(4) = xline(2)
              yline(1) = ymin
              yline(4) = ymax
              print*,' error: gradient = infinite'
           elseif (xline(2).lt.xline(3)) then 
              xline(1) = xmin
              xline(4) = xmax
           else
              xline(1) = xmax
              xline(4) = xmin
           endif           
           if (abs(xlength).gt.tiny(xlength)) then
              ylength = yline(3)-yline(2)
              gradient = ylength/xlength
              yint = yline(3) - gradient*xline(3)
              print*,' dx = ',xlength,' dy = ',ylength
              print*,' gradient = ',gradient,' y intercept = ',yint
              yline(1) = gradient*xline(1) + yint
              yline(4) = gradient*xline(4) + yint
           endif
           !--plot line joining the two points
           call pgline(4,xline,yline)           
        end select
     !
     !--help
     !
     case('h')
        print*,'-------------- interactive mode commands --------------'
        print*,' select region : left click (or A)'
        print*,': left click again to zoom on selection'
        if (irender.ne.0) then
           print*,': or select colour bar to change rendering limits'
        else
           print*,': or press 1-9 to mark selected particles with colour 1-9'
        endif
        print*,' zoom in by 10%       : +'
        print*,' zoom out by 10(20)%      : - (_)'
        print*,' (a)djust/reset plot limits to fit '
        print*,' (r)eplot current plot        : r'
        print*,' label closest (p)article     : p'
        print*,' plot a line and find its g)radient : g'
        print*,' G : move legend to current position'
        print*,' T : move title to current position'
        if (rotation) then
           print*,' rotate about z axis by +(-) 15 degrees : , (.)'
           print*,' rotate about z axis by +(-) 30 degrees : < (>)'
           if (ndim.ge.3) then
              print*,' rotate about x axis by +(-) 15 degrees : / ('')'
              print*,' rotate about x axis by +(-) 30 degrees : ? (")'
              print*,' rotate about y axis by +(-) 15 degrees : l (;)'
              print*,' rotate about y axis by +(-) 30 degrees : L (:)'
           endif
        endif
        print*,' next timestep/plot   : space, n'
        print*,' previous timestep    : right click (or X), b'
        print*,' jump forward (back) by n timesteps  : 0,1,2,3..9 then left (right) click'
        print*,' (h)elp                       : h'
        print*,' (s)ave current settings for all steps : s'
        print*,' (q)uit plotting              : q, Q'             
        print*,'-------------------------------------------------------'
     case('s','S')
        isave = .not.isave
        print*,'saving plot settings...',isave
        call save_limits(iplotx,xmin,xmax)
        call save_limits(iploty,ymin,ymax)
        if (irender.gt.0) call save_limits(irender,rendermin,rendermax)
        if (rotation) call save_rotation(ndim,anglex,angley,anglez)
     !
     !--zoom
     !
     case('A') ! left click
        !
        !--draw rectangle from the point and reset the limits
        !
        print*,'select area: '
        print*,'left click : zoom'
        if (irender.le.0) then
           print*,'1-9 = mark selected particles with colours 1-9'
        else
           print*,'select colour bar to change rendering limits'
        endif
        call pgband(2,1,xpt,ypt,xpt2,ypt2,char2)
        print*,xpt,ypt,xpt2,ypt2,char2
        select case (char2)
        case('A')   ! zoom if another left click and not around colour bar
           call pgsfs(2)
           call pgrect(xpt,xpt2,ypt,ypt2)
           if (xpt.gt.xmax .and. xpt2.gt.xmax .and. irender.gt.0) then
              drender = (rendermax-rendermin)/(ymax-ymin)
              rendermax = rendermin + (max(ypt,ypt2)-ymin)*drender
              rendermin = rendermin + (min(ypt,ypt2)-ymin)*drender
              print*,'setting render min = ',rendermin
              print*,'setting render max = ',rendermax
           else
              xmin = min(xpt,xpt2)
              xmax = max(xpt,xpt2)
              ymin = min(ypt,ypt2)
              ymax = max(ypt,ypt2)
           endif
           iadvance = 0
           iexit = .true.
        case('1','2','3','4','5','6','7','8','9')
           if (irender.le.0) then
              xptmin = min(xpt,xpt2)
              xptmax = max(xpt,xpt2)
              yptmin = min(ypt,ypt2)
              yptmax = max(ypt,ypt2)
           
              nmarked = 0
              do i=1,npart
                 if ((xcoords(i).ge.xptmin .and. xcoords(i).le.xptmax) &
                 .and.(ycoords(i).ge.yptmin .and. ycoords(i).le.yptmax)) then
                     read(char2,*,iostat=ierr) icolourpart(i)
                     if (ierr /=0) then
                        print*,'*** error marking particle' 
                        icolourpart(i) = 1
                     endif
                     nmarked = nmarked + 1
                 endif
              enddo
              print*,'marked ',nmarked,' particles in selected region'
              iadvance = 0
              iexit = .true.
           endif
        end select    
     case('-','_') ! zoom out by 10 or 20%
        print*,'zooming out'
        xlength = xmax - xmin
        ylength = ymax - ymin
        select case(char)
        case('-')
           xlength = 1.1*xlength
           ylength = 1.1*ylength
        case('_')
           xlength = 1.2*xlength
           ylength = 1.2*ylength
        end select
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
     case('a') ! reset plot limits
        print*,'resetting plot limits...'
        xmin = minval(xcoords)
        xmax = maxval(xcoords)
        ymin = minval(ycoords)
        ymax = maxval(ycoords)
     !
     !--rotation
     !
     case(',')
        if (rotation) then
           print*,'changing z rotation angle by -15 degrees...'
           anglez = anglez - 15.
           iadvance = 0
           iexit = .true.
        endif
     case('<')
        if (rotation) then
           print*,'changing z rotation angle by -30 degrees...'
           anglez = anglez - 30.
           iadvance = 0
           iexit = .true.
        endif
     case('.')
        if (rotation) then
           print*,'changing z rotation angle by 15 degrees...'
           anglez = anglez + 15.
           iadvance = 0
           iexit = .true.
        endif
     case('>')
        if (rotation) then
           print*,'changing z rotation angle by 30 degrees...'
           anglez = anglez + 30.
           iadvance = 0
           iexit = .true.
        endif
     case('l')
        if (rotation .and. ndim.ge.3) then
           print*,'changing y rotation angle by -15 degrees...'
           angley = angley - 15.
           iadvance = 0
           iexit = .true.
        endif
     case('L')
        if (rotation .and. ndim.ge.3) then
           print*,'changing y rotation angle by -30 degrees...'
           angley = angley - 30.
           iadvance = 0
           iexit = .true.
        endif
     case(';')
        if (rotation .and. ndim.ge.3) then
           print*,'changing y rotation angle by 15 degrees...'
           angley = angley + 15.
           iadvance = 0
           iexit = .true.
        endif
     case(':')
        if (rotation .and. ndim.ge.3) then
           print*,'changing y rotation angle by 30 degrees...'
           angley = angley + 30.
           iadvance = 0
           iexit = .true.
        endif
     case('''')
        if (rotation .and. ndim.ge.3) then
           print*,'changing x rotation angle by -15 degrees...'
           anglex = anglex - 15.
           iadvance = 0
           iexit = .true.
        endif
     case('"')
        if (rotation .and. ndim.ge.3) then
           print*,'changing x rotation angle by -30 degrees...'
           anglex = anglex - 30.
           iadvance = 0
           iexit = .true.
        endif
     case('/')
        if (rotation .and. ndim.ge.3) then
           print*,'changing x rotation angle by 15 degrees...'
           anglex = anglex + 15.
           iadvance = 0
           iexit = .true.
        endif
     case('?')
        if (rotation .and. ndim.ge.3) then
           print*,'changing x rotation angle by 30 degrees...'
           anglex = anglex + 30.
           iadvance = 0
           iexit = .true.
        endif
     !
     !--general plot stuff
     !
     case('G') ! move legend here
        print*,'setting legend position to current location...'
        call mvlegend(xpt,ypt,xmin,xmax,ymax)
     case('T') ! move title here
        print*,'setting title position to current location...'
        call mvtitle(xpt,ypt,xmin,xmax,ymax)
     !
     !--timestepping
     !
     case('q','Q')
        iadvance = -666
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
        read(char,*,iostat=ierr) iadvance
        if (ierr /=0) then
           print*,'*** error setting timestep jump' 
           iadvance = 1
        endif
        print*,' setting timestep jump = ',iadvance
     case(')')
        iadvance = 10
        print*,' setting timestep jump = ',iadvance
     end select

  enddo
  return
end subroutine interactive_part

!
! cut down version of interactive mode -> controls timestepping only
! used in powerspectrum / extra plots
!
subroutine interactive_step(iadvance,xmin,xmax,ymin,ymax)
 implicit none
 integer, intent(inout) :: iadvance
 real, intent(inout) :: xmin,xmax,ymin,ymax
 integer :: nc,ierr
 real :: xpt,ypt,xpt2,ypt2
 real :: xlength, ylength
 character(len=1) :: char,char2
 character(len=5) :: string
 logical :: iexit
 
  call pgqinf('CURSOR',string,nc)
  if (string(1:nc).eq.'YES') then
     print*,'entering interactive mode...press h in plot window for help'
  else
     print*,'cannot enter interactive mode: device has no cursor'
     return
  endif
  char = 'A'
  xpt = 0.
  ypt = 0.
  iexit = .false.
  
  do while (.not.iexit)
     call pgcurs(xpt,ypt,char)
     !
     !--exit if the device is not interactive
     !
     if (char.eq.achar(0)) return
  
     print*,'x, y = ',xpt,ypt,' function = ',char
     
     select case(char)
     case('h')
        print*,'-------------- interactive mode commands --------------'
        print*,' select area and zoom : left click (or A)'
        print*,' zoom in by 10%       : +'
        print*,' zoom out by 10(20)%      : - (_)'
        print*,' (r)eplot current plot        : r'
        print*,' next timestep/plot   : space, n'
        print*,' previous timestep    : right click (or X), b'
        print*,' jump forward (back) by n timesteps  : 0,1,2,3..9 then left (right) click'
        print*,' G : move legend to current position'
        print*,' T : move title to current position'
        print*,' (h)elp                       : h'
        print*,' (q)uit plotting              : q, Q'             
        print*,'-------------------------------------------------------'
     !
     !--zoom
     !
     case('A') ! left click
        !
        !--draw rectangle from the point and reset the limits
        !
        print*,'select area: '
        print*,'left click : zoom'
        call pgband(2,1,xpt,ypt,xpt2,ypt2,char2)
        print*,xpt,ypt,xpt2,ypt2,char2
        select case (char2)
        case('A')   ! zoom if another left click
           call pgrect(xpt,xpt2,ypt,ypt2)
           xmin = min(xpt,xpt2)
           xmax = max(xpt,xpt2)
           ymin = min(ypt,ypt2)
           ymax = max(ypt,ypt2)
           iadvance = 0
           iexit = .true.
        end select    
     case('-','_') ! zoom out by 10 or 20%
        print*,'zooming out'
        xlength = xmax - xmin
        ylength = ymax - ymin
        select case(char)
        case('-')
           xlength = 1.1*xlength
           ylength = 1.1*ylength
        case('_')
           xlength = 1.2*xlength
           ylength = 1.2*ylength
        end select
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
     !--general plot stuff
     !
     case('G') ! move legend here
        print*,'setting legend position to current location...'
        call mvlegend(xpt,ypt,xmin,xmax,ymax)
     case('T') ! move title here
        print*,'setting title position to current location...'
        call mvtitle(xpt,ypt,xmin,xmax,ymax)
     !
     !--timestepping
     !
     case('q','Q')
        iadvance = -666
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
        read(char,*,iostat=ierr) iadvance
        if (ierr /=0) then
           print*,'*** error setting timestep jump' 
           iadvance = 1
        endif
        print*,' setting timestep jump = ',iadvance
     case(')')
        iadvance = 10
        print*,' setting timestep jump = ',iadvance
     end select

  enddo
  return
end subroutine interactive_step

!
!--this subroutine moves the legend to the current position
!
subroutine mvlegend(xi,yi,xmin,xmax,ymax)
 use settings_page, only:hposlegend,vposlegend
 implicit none
 real, intent(in) :: xi,yi,xmin,xmax,ymax
 real :: xch,ych
 
 hposlegend = (xi - xmin)/(xmax-xmin)
 !--query character height in world coordinates
 call pgqcs(4,xch,ych)
 vposlegend = (ymax - yi)/ych
 print*,'hpos = ',hposlegend,' vpos = ',vposlegend
 
 return
end subroutine mvlegend

!
!--this subroutine moves the title to the current position
!
subroutine mvtitle(xi,yi,xmin,xmax,ymax)
 use settings_page, only:hpostitle,vpostitle
 implicit none
 real, intent(in) :: xi,yi,xmin,xmax,ymax
 real :: xch,ych
 
 hpostitle = (xi - xmin)/(xmax-xmin)
 !--query character height in world coordinates
 call pgqcs(4,xch,ych)
 vpostitle = (yi - ymax)/ych
 print*,'hpos = ',hpostitle,' vpos = ',vpostitle
 
 return
end subroutine mvtitle

!
!--subroutines to save current plot limits
!
subroutine save_limits(iplot,xmin,xmax)
 use limits, only:lim
 implicit none
 integer, intent(in) :: iplot
 real, intent(in) :: xmin,xmax
 
 lim(iplot,1) = xmin
 lim(iplot,2) = xmax
 
end subroutine save_limits

!
!--subroutines to save current plot limits
!
subroutine save_rotation(ndim,anglexi,angleyi,anglezi)
 use settings_xsecrot, only:anglex,angley,anglez
 implicit none
 integer, intent(in) :: ndim
 real, intent(in) :: anglexi,angleyi,anglezi
 
 anglez = anglezi
 if (ndim.ge.3) then
    anglex = anglexi
    angley = angleyi
 endif
 
end subroutine save_rotation

end module interactive_routines
