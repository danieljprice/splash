!----------------------------------------------------------------
!
! module for probability distribution function calculation
! and/or plotting on the particles
!
!----------------------------------------------------------------
module pdfs
 implicit none
 public :: pdfcalc

contains

subroutine pdfcalc(npart,xpart,xmin,xmax,nbins,xbin,pdf,pdfmin,pdfmax,itransx,itransy)
 use transforms, only:transform,transform_limits
 implicit none
 integer, intent(in) :: npart,nbins
 real, dimension(:), intent(in) :: xpart
 real, intent(in) :: xmin,xmax
 real, intent(out), dimension(nbins) :: xbin,pdf
 real, intent(out) :: pdfmin,pdfmax
 integer, intent(in) :: itransx,itransy
 integer :: ibin,i
 real :: dx,totprob,fi,fprev
 
 print "(a,i3,a)",' calculating PDF using ',nbins,' bins'
 !
 !--note that xmin, xmax and xpart will already 
 !  have been transformed prior to input
 !
 dx = (xmax - xmin)/real(nbins - 1)
 do ibin=1,nbins
    xbin(ibin) = xmin + (ibin-0.5)*dx
 enddo
! print*,'xmin,max = ',xmin,xmax
!
!--now calculate probability of finding a particle at each x
!
 pdf(:) = 0.
 do i=1,npart
    ibin = int((xpart(i) - xmin)/dx) + 1
    if (ibin.lt.1) ibin = 1
    if (ibin.gt.nbins) ibin = nbins
    pdf(ibin) = pdf(ibin) + 1.
 enddo
!
!--get total area under pdf by trapezoidal rule
!
! totprob = 0.
! fprev = 0.
! do ibin=1,nbins
!    fi = pdf(ibin)
!   ! if (itransx.eq.1) then
!   ! ! \int pdf dx = ln(10)*\int x*pdf d(log_10 x)
!   !    totprob = totprob + 0.5*log(10.)*dx*xbin(ibin)*(fi + fprev)    
!   ! else
!       totprob = totprob + 0.5*dx*(fi + fprev)
!   ! endif
!    fprev = fi
! enddo
!
!--normalise pdf so total area is unity
!
! print*,'normalisation factor = ',totprob,npart*dx
 pdf(1:nbins) = pdf(1:nbins)/(npart*dx)
 
 pdfmin = minval(pdf(1:nbins),mask=(pdf(1:nbins).gt.0.))
 pdfmax = maxval(pdf(1:nbins))
 if (itransy.gt.0) then
    call transform(pdf,itransy)
    call transform_limits(pdfmin,pdfmax,itransy)
 endif
 
end subroutine pdfcalc

subroutine write_pdf(iunit,basename,variable,nbins,xbin,pb)
 implicit none
 integer, intent(in) :: iunit
 character(len=*), intent(in) :: basename,variable
 integer, intent(in) :: nbins
 real, intent(in), dimension(nbins) :: xbin,pb
 integer :: i
 
 print "(a)",' writing to '//trim(basename)//'_'//trim(variable)//'.pdf'
 open(unit=iunit,file=trim(basename)//'_'//trim(variable)//'.pdf', &
      form='formatted',status='replace')
 write(iunit,*) nbins,trim(variable)
 do i=1,nbins
    write(iunit,*) xbin(i),pb(i)
 enddo
 close(iunit)
 
end subroutine write_pdf

end module pdfs
