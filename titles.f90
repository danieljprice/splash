module titles
 implicit none
 integer, parameter, private :: maxtitles = 50
 integer, parameter, private :: maxsteplegend = 100
 character(len=60), dimension(maxtitles), public :: pagetitles
 character(len=60), dimension(maxsteplegend), public :: steplegend
 public :: read_titles, read_steplegend
 
 private
 
contains
!
!--reads a list of titles (one per line), to be used to label each plot on page
!
subroutine read_titles(ntitles)
 use asciiutils, only:read_asciifile
 implicit none
 integer, intent(out) :: ntitles
 integer :: i
 character(len=50) :: titlefile

 titlefile = 'titlelist'
 ntitles = 0

 print*,'reading plot titles from file ',trim(titlefile)  
 call read_asciifile(titlefile,ntitles,pagetitles)

 return
end subroutine read_titles
!
!--reads a list of labels (one per line) to be used in the timestep legend
! (ie. for multiple timesteps on same page)
!
subroutine read_steplegend(nsteptitles)
 use asciiutils, only:read_asciifile
 implicit none
 integer, intent(out) :: nsteptitles
 integer :: i
 character(len=50) :: legendfile

 legendfile = 'legend'
 nsteptitles = 0

 print*,'reading legend text from file ''',trim(legendfile),''''
 call read_asciifile(legendfile,nsteptitles,steplegend)

 return
end subroutine read_steplegend

end module titles
