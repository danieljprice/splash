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
 logical :: iexist

 titlefile = 'splash.titles'
 !--also allow obsolete title filename
 inquire(file=titlefile,exist=iexist)
 if (.not.iexist) then
    inquire(file='titlelist',exist=iexist)
    if (iexist) titlefile='titlelist'
 endif
 ntitles = 0

 call read_asciifile(titlefile,ntitles,pagetitles)
 if (ntitles.gt.0) print "(a)",'read plot titles from file '//trim(titlefile)  

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
 logical :: iexist

 legendfile = 'splash.legend'
 !--also allow obsolete legend filename
 inquire(file=legendfile,exist=iexist)
 if (.not.iexist) then
    inquire(file='legend',exist=iexist)
    if (iexist) legendfile='legend'
 endif
 nsteptitles = 0

 call read_asciifile(legendfile,nsteptitles,steplegend)
 if (nsteptitles.gt.0) print "(a)"//'read legend text from file '''//trim(legendfile)//''''

 return
end subroutine read_steplegend

end module titles
