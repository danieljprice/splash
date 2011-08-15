!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2009 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
!
!-----------------------------------------------------------------

!----------------------------------------------------------------
!
! module for probability density function calculation
! and/or plotting on the particles
!
!----------------------------------------------------------------
module pdfs
 implicit none
 public :: pdf_calc,pdf_write,mean_variance

contains
!-----------------------------------------------------------------
!
! subroutine bins particles into x, works out number in each bin,
! calculates normalisation for PDF.
!
!-----------------------------------------------------------------
subroutine pdf_calc(npart,xpart,xminplot,xmaxplot,nbins,xbin,pdf,pdfmin,pdfmax,&
                    usefixedbins,volweighted,ierr,icolours,rhopart,pmass)
 !use transforms, only:transform,transform_inverse,transform_limits,convert_to_ln_fac
 implicit none
 integer, intent(in)                 :: npart,nbins
 real, dimension(:), intent(in)      :: xpart
 real, intent(in)                    :: xminplot,xmaxplot
 real, intent(out), dimension(nbins) :: xbin,pdf
 real, intent(out)                   :: pdfmin,pdfmax
! integer, intent(in)                 :: itransx
 logical, intent(in)                 :: usefixedbins
 logical, intent(out)                :: volweighted
 integer, intent(out)                :: ierr
 integer, intent(in), dimension(:), optional :: icolours
 real, intent(in), dimension(:),    optional :: rhopart,pmass
 integer :: ibin,i
 real    :: dx,totprob,fi!,xbinprev !xbini,dxprev
 real    :: xmin,xmax,xminpart,xmaxpart,weighti,totvol
 logical :: use_part

 ierr = 0
 volweighted = .false.

 if (present(rhopart) .and. present(pmass)) then
    print "(a,i3,a)",' calculating (volume weighted) PDF using ',nbins,' bins'
    volweighted = .true.
 else
    print "(a,i3,a)",' calculating (mass weighted) PDF using ',nbins,' bins'
 endif
 !
 !--set bins in PDF: must always use all the particles
 !  note that xpart will already have been transformed
 !  so these are min and max in transformed space
 !
 xminpart = minval(xpart(1:npart))
 xmaxpart = maxval(xpart(1:npart))
 if (usefixedbins) then
    xmin = xminplot
    xmax = xmaxplot
    print "(a,1pe10.3,a,1pe10.3)",' PDF bins are fixed between the current x limits, min = ',xminplot,' max = ',xmaxplot
    if (xminpart.lt.xmin) print "(a)",' WARNING: particles fall outside of (fixed) bin range, will pile up on first bin'
    if (xmaxpart.gt.xmax) print "(a)",' WARNING: particles fall outside of (fixed) bin range, will pile up on last bin'
    dx = (xmax - xmin)/real(nbins)
    print "(a,1pe10.3)",' bin width = ',dx
    do ibin=1,nbins
       xbin(ibin) = xmin + (ibin-1)*dx
    enddo
 else
    xmin = xminpart
    xmax = xmaxpart
    dx = (xmax - xmin)/real(nbins-1)
    print "(a,1pe10.3)",' bin width = ',dx
    do ibin=1,nbins
       xbin(ibin) = xmin + (ibin-0.5)*dx
    enddo
 endif

!
!--now calculate probability of finding a particle at each x
!
 pdf(:) = 0.
 totvol = 0.
 do i=1,npart
    if (present(icolours)) then
       use_part = (icolours(i).ge.0)
    else
       use_part = .true.
    endif
    !--do not use hidden particles
    if (use_part) then
       ibin = int((xpart(i) - xmin)/dx) + 1
       if (ibin.lt.1) ibin = 1
       if (ibin.gt.nbins) ibin = nbins

       if (volweighted) then
          if (rhopart(i).gt.0.) then
             weighti = pmass(i)/rhopart(i)
          else
             weighti = 0.
          endif
       else
          weighti = 1.
       endif
       totvol = totvol + weighti

       !!--take the PDF of ln(x) if quantity is logged
       !if (itransx.gt.0) then
       !    weighti = weighti*convert_to_ln_fac(itransx)
       !endif
       pdf(ibin) = pdf(ibin) + weighti
    endif
 enddo
 print*,' sum of weights = ',totvol
!
!--get total area under pdf by trapezoidal rule
!
 totprob = 0.
 do ibin=1,nbins
    fi = pdf(ibin)
    totprob = totprob + dx*fi !!0.5*dx*(fi + fprev)
 enddo
!
!--normalise pdf so total area is unity
!
 print*,'normalisation factor = ',totprob,totvol*dx ! =npart*dx for mass-weighted, totvol*dx for volume weighted
 !totprob = totvol*dx
 !totprob = dx

 if (totprob.le.0.) then
    ierr = 1
    print "(a)",' error in normalisation factor: returning non-normalised PDF'
 else
    pdf(1:nbins) = pdf(1:nbins)/totprob

    !call pdf_write(nbins,xbin,pdf,labelx,itransx,volweighted)
   !
   !--return min and max for adaptive plot limit setting
   !  (exclude zero as min)
   !
    pdfmin = minval(pdf(1:nbins),mask=(pdf(1:nbins).gt.0.))
    pdfmax = maxval(pdf(1:nbins))
 endif

end subroutine pdf_calc

!-----------------------------------------------------------------
! interface which controls plotting of PDF
! (so can easily change properties of PDF plotting,
!  e.g. histogram vs. line)
!-----------------------------------------------------------------
!subroutine pdf_plot(nbins,xbin,pb)
! use plotutils, only:plotline !,plotbins
! implicit none
! integer, intent(in) :: nbins
! real, dimension(:), intent(in) :: xbin,pb

!
!--plot as line segment, with blanking at zero
!
! call plotline(nbins,xbin,pb,blank=0.)
!
!--plot as histogram, with blanking of zero
!
! call plotbins(nbins,xbin,pb,blank=0.)

!end subroutine pdf_plot

!-----------------------------------------------------------------
! routine to write pdf to file
!-----------------------------------------------------------------
subroutine pdf_write(nbins,xbin,pb,labelx,volweighted,rootname,tagline)
 use asciiutils, only:safename
 implicit none
 character(len=*), intent(in)       :: labelx,rootname,tagline
 integer, intent(in)                :: nbins !,itransx
 real, intent(in), dimension(nbins) :: xbin,pb
 logical, intent(in)                :: volweighted
 integer                :: i,ierr
 integer, parameter     :: iunit = 86
 logical                :: warned

 print "(a)",' writing to '//trim(rootname)//'_pdf_'//trim(safename(labelx))//'.dat'
 open(unit=iunit,file=trim(rootname)//'_pdf_'//trim(safename(labelx))//'.dat', &
      form='formatted',status='replace',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",'ERROR: could not open file: no output'
    return
 endif

 if (volweighted) then
    write(iunit,"(a)",iostat=ierr) '# volume weighted PDF, calculated using '//trim(tagline)
 else
    write(iunit,"(a)",iostat=ierr) '# density weighted PDF, calculated using '//trim(tagline)
 endif
 if (ierr /= 0) print "(a)",' ERROR writing header line'
 write(iunit,"(a,i5,a)",iostat=ierr) '# ',nbins,' bins evenly spaced in '//trim(labelx)

 warned = .false.
 !--dump bins to file
 do i=1,nbins
    write(iunit,*,iostat=ierr) xbin(i),pb(i)
    if (ierr /= 0 .and. .not.warned) then
       print "(a)",' ERRORS during write'
       warned = .true.
    endif
 enddo
 close(iunit)

 return
end subroutine pdf_write

!-------------------------------------------------
! Subroutine to calculate the mean and variance
! of a set of data points
! Mean is trivial but variance uses a special
! formula to reduce round-off error
! see Press et al Numerical Recipes, section 14.2
! this is similar to their subroutine avevar
!-------------------------------------------------
subroutine mean_variance(x,npts,xmean,xvariance)
 implicit none
 integer, intent(in) :: npts
 real, intent(in), dimension(npts) :: x
 real, intent(out) :: xmean, xvariance
 real :: roundoff, delta
 integer :: i
!
!--calculate average
!
 xmean = 0.
 do i=1,npts
    xmean = xmean + x(i)
 enddo
 xmean = xmean/real(npts)
!
!--calculate variance using the corrected two-pass formula
!
!    var = 1/(n-1)*( sum (x-\bar{x}) - 1/n * (sum(x-\bar{x}) )^2 )
!
!  where the last term corrects for the roundoff error
!  in the first term
!
 xvariance = 0.
 roundoff = 0.

 do i=1,npts
    delta = x(i) - xmean
    roundoff = roundoff + delta
    xvariance = xvariance + delta*delta
 enddo
 xvariance = (xvariance - roundoff**2/npts)/real(npts-1)

 return
end subroutine mean_variance

end module pdfs
