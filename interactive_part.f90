module interactive_routines
 implicit none
contains
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
subroutine interactive_part(npart,iplotx,iploty,irender,xcoords,ycoords,hi, &
  icolourpart,xmin,xmax,ymin,ymax,anglex,angley,anglez,ndim,iadvance,isave)
  implicit none
  integer, intent(in) :: npart,iplotx,iploty,irender,ndim
  integer, intent(out) :: iadvance
  integer, dimension(npart), intent(inout) :: icolourpart
  real, dimension(npart), intent(in) :: xcoords, ycoords, hi
  real, intent(inout) :: xmin,xmax,ymin,ymax
  real, intent(inout) :: anglex,angley,anglez
  logical, intent(out) :: isave
  integer :: i,iclosest,nc,ipts,ierr
  integer :: nmarked
  real :: xpt,ypt,xpt2,ypt2,xptmin,xptmax,yptmin,yptmax
  real :: rmin,rr,gradient,yint
  real :: xlength, ylength
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
        print*,'plotting circle of interaction on particle ',iclosest
        call pgsfs(2)
        call pgcirc(xcoords(iclosest),ycoords(iclosest),2.*hi(iclosest))
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
        print*,' zoom in by 10%       : +'
        print*,' zoom out by 10(20)%      : - (_)'
        print*,' (a)djust/reset plot limits to fit '
        print*,' (r)eplot current plot        : r'
        print*,' label closest (p)article     : p'
        print*,' plot a line and find its g)radient : g'
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
        print*,'save settings on exit = ',isave
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
        endif
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
subroutine interactive_step(iadvance)
 implicit none
 integer, intent(inout) :: iadvance
 integer :: nc,ierr
 real :: xpt,ypt
 character(len=1) :: char
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
  
     !print*,'location: x, y = ',xpt,ypt,' function = ',char
     
     select case(char)
     case('h')
        print*,'-------------- interactive mode commands --------------'
        print*,' (r)eplot current plot        : r'
        print*,' next timestep/plot   : space, n'
        print*,' previous timestep    : right click (or X), b'
        print*,' jump forward (back) by n timesteps  : 0,1,2,3..9 then left (right) click'
        print*,' (h)elp                       : h'
        print*,' (q)uit plotting              : q, Q'             
        print*,'-------------------------------------------------------'

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

subroutine mvlegend(hpos,vpos)
 use settings_page, only:hposlegend,vposlegend
 implicit none
 real, intent(in) :: hpos,vpos
 
 hposlegend = hpos
 vposlegend = vpos
 
 return
end subroutine mvlegend

end module interactive_routines
