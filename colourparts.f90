module colourparts
 implicit none

contains

subroutine colour_particles(dat,datmin,datmax,icolour,npart)
 implicit none
 integer, intent(in) :: npart
 real, dimension(npart), intent(in) :: dat
 real, intent(in) :: datmin, datmax
 integer, dimension(npart), intent(out) :: icolour
 integer :: i,icolourmin,icolourmax,icolourtemp
 real :: dx
 
 call pgqcir(icolourmin,icolourmax)
 
 dx = (datmax - datmin)/real(icolourmax - icolourmin)
 
 do i=1,npart
    icolourtemp = int((dat(i) - datmin)/dx) + icolourmin
    if (icolourtemp.gt.icolourmax) icolourtemp = icolourmax
    if (icolourtemp.lt.icolourmin) icolourtemp = icolourmin
    icolour(i) = icolourtemp
 enddo

 return
end subroutine colour_particles

end module colourparts
