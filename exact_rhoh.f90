!----------------------------------------------------------
! plots the relation between smoothing length and density
!
! ie. h = h_fact*(pmass/rho)^(1/ndim)
!
!----------------------------------------------------------
module rhoh
 implicit none
 public :: exact_rhoh

contains

subroutine exact_rhoh(ndim,hfact,pmassval,xplot,yplot,ierr)
 implicit none
 integer, intent(in) :: ndim
 integer, intent(out) :: ierr
 real, intent(in) :: hfact,pmassval
 real, dimension(:), intent(in) :: xplot
 real, dimension(size(xplot)), intent(out) :: yplot
       
 if (hfact.gt.0.01) then
    ierr = 0

    where (xplot > 0)    
      yplot(:) = hfact*(pmassval/xplot(:))**(1./FLOAT(ndim))
    elsewhere
      yplot(:) = huge(yplot)
    end where
    write(*,"(a,f5.2,a,1pe8.2,a,i1,a)") ' plotting h = ',hfact, &
                               '*(',pmassval,'/rho)**(1/',ndim,')'
 else
    print "(a)",'error: hfact = 0: can''t plot h vs rho exact solution'
    ierr = 1
 endif  
 
 return
end subroutine exact_rhoh

end module rhoh
