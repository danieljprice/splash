!
!--interactive tools on particle plots
!  (experimental at this stage)
!
subroutine interactive_part(npart,xcoords,ycoords)
  implicit none
  integer, intent(in) :: npart
  real, dimension(npart), intent(in) :: xcoords, ycoords
  integer :: i,iclosest,nc
  real :: xpt,ypt,rmin,rr
  character(len=1) :: char
  character(len=20) :: string

  print*,'entering interactive mode... type q to quit'
  char = 'A'

  do while (char.ne.'q' .and. char.ne.'Q' .and. char.ne.' ')
     call pgcurs(xpt,ypt,char)
     
     print*,'location: x, y = ',xpt,ypt,' function = ',char
     
     rmin = 1.e6
     do i=1,npart
        rr = (xcoords(i)-xpt)**2 + (ycoords(i)-ypt)**2
        if (rr.lt.rmin) then
           iclosest = i
           rmin = rr
        endif
     enddo
     
     print*,' closest particle = ',iclosest,'x = ',xcoords(iclosest),' y =',ycoords(iclosest)
     if (char.eq.'A') then
        call pgnumb(iclosest,0,1,string,nc)
        call pgsch(2.0)
        call pgtext(xcoords(iclosest),ycoords(iclosest),string(1:nc))
        call pgsch(1.0)
     endif
     
     if (char.eq.'h') then
        print*,'plotting circle of interaction (not implemented)'
     endif

  enddo

  return
end subroutine interactive_part
