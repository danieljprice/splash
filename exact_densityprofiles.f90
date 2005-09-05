! ----------------------------------------------------------------------
! compute exact solution for a linear wave
! plots a sine function with a given amplitude, period and wavelength
! ----------------------------------------------------------------------
module densityprofiles
  implicit none

contains

subroutine exact_densityprofiles(iprofile,Msphere,rsoft,xplot,yplot,ierr)
  implicit none
  real, parameter :: pi = 3.1415926536
  integer, intent(in) :: iprofile
  real, intent(in) :: Msphere,rsoft 
  real, intent(in), dimension(:) :: xplot
  real, intent(out), dimension(size(xplot)) :: yplot
  integer, intent(out) :: ierr
  integer :: i
!
! check for errors
!
  ierr = 0
  if (Msphere.le.0.) then
     print*,'error: mass <= 0 in exact_densityprofile'
     ierr = 2
     return
  endif
  if (rsoft.lt.0.) then
     print*,'error: rsoft < 0 in exact_densityprofile'
     ierr = 3
     return
  endif

  select case(iprofile)
  case(1)
!
!--Plummer sphere
!
     do i=1,size(xplot)
        yplot(i) = 3.*Msphere*rsoft**2/(4.*pi*(rsoft**2 + xplot(i)**2)**2.5)
     enddo
  case(2)
!
!--Hernquist model (use tiny to prevent divergences in cusp)
!
     do i=1,size(xplot)
        yplot(i) = Msphere*rsoft/ &
                   ((2.*pi*xplot(i)*(rsoft + xplot(i))**3) + tiny(rsoft))
     enddo
  case default
     ierr = 1  
  end select
    
  return
end subroutine exact_densityprofiles

end module densityprofiles
