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

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! the data is stored in the global array dat
!
! THIS VERSION FOR DAN'S SPMHD CODE (BINARY DUMPS)
! PRE NOV 2005 FORMAT
! -> Now automatically handles single/double precision
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! maxplot,maxpart,maxstep      : dimensions of main data array
! dat(maxplot,maxpart,maxstep) : main data array
!
! npartoftype(maxstep) : number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!-------------------------------------------------------------------------

subroutine read_data(rootname,indexstart,ipos,nstepsread)
  use exact, only:hfact
  use particle_data
  use params
  use labels
  use filenames, only:nfiles
  use settings_data, only:ndim,ndimV,ncolumns,ncalc,icoords,iformat, &
                          buffer_data
  use mem_allocation
  use geometry, only:labelcoordsys
  implicit none
  integer, intent(in) :: indexstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  character(len=len(rootname)+4) :: datfile
  integer :: i,j,icol,ierr,iunit
  integer :: ncol_max,ndim_max,npart_max,ndimV_max,nstep_max
  integer :: npartin,ntotin,ncolstep,nparti,ntoti
  integer, dimension(3) :: ibound
  logical :: reallocate, singleprecision

  real :: timein, gammain, hfactin
  real, dimension(3) :: xmin, xmax
  real(doub_prec) :: timeind,gammaind,hfactind
  real(doub_prec), dimension(3) :: xmind, xmaxd

  iunit = 11 ! file unit number
  ndim_max = 1
  ndimV_max = 1
  nstepsread = 0
  if (rootname(1:1).ne.' ') then
     datfile = trim(rootname)
     !print*,'rootname = ',rootname
  else
     print*,' **** no data read **** '
     return
  endif

  print "(1x,a)",'reading old ndspmhd format'
  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  !--open data file and read data
  !
  open(unit=iunit,iostat=ierr,file=datfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print*,'*** Error opening '//trim(datfile)//' ***'
     return
  endif
!
!--read first header line
!
  singleprecision = .false.
  read(iunit,iostat=ierr,end=80) timeind,npartin,ntotin,gammaind, &
       hfactind,ndim_max,ndimV_max,ncol_max,icoords
  !!print*,'time = ',timeind,' hfact = ',hfactind,' ndim=',ndim_max,'ncol=',ncol_max
  !!print*,'npart = ',npartin,ntotin
  if (ierr /= 0 .or. ndim_max.le.0 .or. ndim_max.gt.3 &
     .or. ndimV_max.le.0 .or. ndimV_max.gt.3 &
     .or. ncol_max.le.0 .or. ncol_max.gt.100 &
     .or. npartin.le.0 .or. npartin.gt.1e7 .or. ntotin.le.0 .or. ntotin.gt.1e7 &
     .or. icoords.le.0 .or. icoords.gt.10) then
     !
     !--try single precision
     !
     rewind(iunit)
     read(iunit,iostat=ierr,end=80) timein,npartin,ntotin,gammain, &
         hfactin,ndim_max,ndimV_max,ncol_max,icoords
     singleprecision = .true.

     if (ierr /= 0) then
        print "(a)",'*** Error reading first header ***'
        print*,' time = ',timein,' hfact = ',hfactin,' ndim=',ndim_max,'ncol=',ncol_max
        close(iunit)
        return
     endif
  endif
!
!--allocate memory for data arrays
!
  if (buffer_data) then
     nstep_max = max(nfiles,maxstep,indexstart)
  else
     nstep_max = max(1,maxstep,indexstart)
  endif
  npart_max = max(int(1.5*ntotin),maxpart)
  if (.not.allocated(dat) .or. ntotin.gt.maxpart  &
       .or. nstep_max.gt.maxstep .or. ncol_max.gt.maxcol) then
     call alloc(npart_max,nstep_max,ncol_max+ncalc)
  endif
!
!--rewind file
!
  rewind(iunit)

  i = indexstart
  nstepsread = 0

  overstepsinfile: do while (i <= maxstep)
     !!print*,' reading step ',i
     reallocate = .false.
     npart_max = maxpart
     nstep_max = maxstep
     !
     !--read header line for this timestep
     !
     if (singleprecision) then
        print "(a)",'single precision dump'
        read(iunit,iostat=ierr,end=67) timein,nparti,ntoti,gammain, &
          hfactin,ndim,ndimV,ncolstep,icoords,iformat,ibound(1:ndim), &
          xmin(1:ndim),xmax(1:ndim)
     else
        print "(a)",'double precision dump'
        read(iunit,iostat=ierr,end=67) timeind,nparti,ntoti,gammaind, &
          hfactind,ndim,ndimV,ncolstep,icoords,iformat,ibound(1:ndim), &
          xmind(1:ndim),xmaxd(1:ndim)
        timein = real(timeind)
        gammain = real(gammaind)
        hfactin = real(hfactind)
        xmin = real(xmind)
        xmax = real(xmaxd)
     endif
     if (ierr /= 0) then
        print*,'*** error reading timestep header ***'
        close(iunit)
        return
     else ! count this as a successfully read timestep, even if data is partial
        nstepsread = nstepsread + 1
     endif

     time(i) = timein
     gamma(i) = gammain
     hfact = hfactin
     npartoftype(1,i) = nparti
     npartoftype(2,i) = ntoti - nparti
     print "(/a14,':',f8.4,a8,':',i8,a8,':',i8)",' time',time(i),'npart',nparti,'ntotal',ntoti
     print "(a14,':',i8,a8,':',f8.4,a8,':',f8.4)",' ncolumns',ncolstep,'gamma',gamma(i),'hfact',hfact
     print "(a14,':',i8,a8,':',i8)",'ndim',ndim,'ndimV',ndimV
     if (icoords.gt.1) print "(a14,':',2x,a)",' geometry',labelcoordsys(icoords)
     if (any(ibound(1:ndim).ne.0)) then
        print "(a14,':',a15,' =',3(f8.4))",'boundaries','xmin',xmin(1:ndim)
        print "(15x,a15,' =',3(f8.4))",'xmax',xmax(1:ndim)
     endif
     !
     !--check for errors in timestep header
     !
     if (ndim.gt.3 .or. ndimV.gt.3) then
        print*,'*** error in header: ndim or ndimV in file> 3'
        nstepsread = nstepsread - 1
        ndim = ndim_max
        ndimV = ndimV_max
        close(iunit)
        return
     endif
     if (ndim.gt.ndim_max) ndim_max = ndim
     if (ndimV.gt.ndimV_max) ndimV_max = ndimV

     if (ncolstep.ne.ncol_max) then
        print*,'*** Warning number of columns not equal for timesteps'
        ncolumns = ncolstep
        print*,'ncolumns = ',ncolumns,ncol_max
        if (ncolumns.gt.ncol_max) ncol_max = ncolumns
     endif
     if (ncolstep.gt.maxcol) then
        reallocate = .true.
        ncolumns = ncolstep
        ncol_max = ncolumns
     else
        ncolumns = ncolstep
     endif

     if (ntoti.gt.maxpart) then
        !print*, 'ntot greater than array limits!!'
        reallocate = .true.
        npart_max = int(1.5*ntoti)
     endif
     if (i.gt.maxstep) then
        nstep_max = i + max(10,INT(0.1*nstep_max))
        reallocate = .true.
     endif
     !
     !--reallocate memory for main data array
     !
     if (reallocate) then
        call alloc(npart_max,nstep_max,ncol_max+ncalc)
     endif


     if (ntoti.gt.0) then
        !
        !--read position vector
        !
        icol = 1
        call readvec(dat(1:ntoti,1:ndim,i),ntoti,ndim,singleprecision,ierr)
        if (ierr /= 0) then
           print "(a)",'*** error reading particle positions ***'
           exit overstepsinfile
        endif
        icol = icol + ndim
        !
        !--read velocity vector
        !
        call readvec(dat(1:ntoti,icol:icol+ndimV-1,i),ntoti,ndimV,singleprecision,ierr)
        if (ierr /= 0) then
           print "(a)",'*** error reading velocities ***'
           exit overstepsinfile
        endif
        icol = icol + ndimV
        !
        !--read scalar variables
        !
        do j=1,4
           call readcol(dat(1:ntoti,icol,i),ntoti,singleprecision,ierr)
           if (ierr /= 0) print "(a)",'*** error reading column data ***'
           icol = icol + 1
        enddo
        !
        !--non-MHD output
        !
          if (iformat.ne.2) then
        !
        !--read alpha, alphau
        !
           call readvec(dat(1:ntoti,icol:icol+1,i),ntoti,2,singleprecision,ierr)
           if (ierr /= 0) then
              print "(a)",'*** error reading alphas ***'
              exit overstepsinfile
           endif
           icol = icol + 2

        !
        !--pr, div v, gradh
        !
           do j=icol,ncolstep
              call readcol(dat(1:ntoti,j,i),ntoti,singleprecision,ierr)
              if (ierr /= 0) print "(a)",'*** error reading column data ***'
           enddo

        else
        !
        !--MHD output
        !
        !
        !--read alpha, alphau, alphaB
        !
           call readvec(dat(1:ntoti,icol:icol+2,i),ntoti,3,singleprecision,ierr)
           if (ierr /= 0) then
              print "(a)",'*** error reading alphas ***'
              exit overstepsinfile
           endif
           icol = icol + 3
        !
        !--Bfield
        !
           call readvec(dat(1:ntoti,icol:icol+ndimV-1,i),ntoti,ndimV,singleprecision,ierr)
           if (ierr /= 0) then
              print "(a)",'*** error reading B ***'
              exit overstepsinfile
           endif
           icol = icol + ndimV
        !
        !--psi, pr, div v, div B
        !
           do j = 1,4
              call readcol(dat(1:ntoti,icol,i),ntoti,singleprecision,ierr)
              if (ierr /= 0) print "(a)",'*** error reading column data ***'
              icol = icol + 1
           enddo
        !
        !--curl B
        !
           call readvec(dat(1:ntoti,icol:icol+ndimV-1,i),ntoti,ndimV,singleprecision,ierr)
           if (ierr /= 0) then
              print "(a)",'*** error reading curl B ***'
              exit overstepsinfile
           endif
           icol = icol + ndimV

        endif
        !!print*,'columns read = ',icol,' should be = ',ncolumns

     else
        npartoftype(1,i) = 1
        npartoftype(2,i) = 0
        dat(:,:,i) = 0.
     endif

     i = i + 1
  enddo overstepsinfile

67 continue
   !!!print "(a)",' > end of file <'
  !
  !--close data file and return
  !
close(unit=11)

ncolumns = ncol_max
ndim = ndim_max
ndimV = ndimV_max

print*,'> Read steps ',indexstart,'->',indexstart + nstepsread - 1, &
       ' last step ntot = ',sum(npartoftype(:,indexstart+nstepsread-1))
return
!
!--errors
!
80 continue
print*,' *** data file empty, no steps read ***'
return

contains

 subroutine readvec(datin,ntotal,ndims,singleprec,ierr)
   implicit none
   integer, intent(in) :: ndims,ntotal
   integer, intent(out) :: ierr
   real, intent(out), dimension(ntotal,ndims) :: datin
   logical, intent(in) :: singleprec

   integer :: ipos
   real, dimension(ndims,ntotal) :: datvec
   real(doub_prec), dimension(ndims,ntotal) :: datvecd
!
!--read a vector quantity and restructure into columns
!
   if (singleprec) then
      read (iunit,iostat=ierr) datvec(1:ndims,1:ntotal)
      do ipos = 1,ndims
         datin(1:ntotal,ipos) = datvec(ipos,1:ntotal)
      enddo
   else
      read (iunit,iostat=ierr) datvecd(1:ndims,1:ntotal)
      do ipos = 1,ndims
         datin(1:ntotal,ipos) = real(datvecd(ipos,1:ntotal))
      enddo
   endif

 end subroutine readvec
!
!--read scalar quantity and convert to single precision
!
 subroutine readcol(datin,ntotal,singleprec,ierr)
   implicit none
   integer, intent(in) :: ntotal
   integer, intent(out) :: ierr
   real, intent(out), dimension(ntotal) :: datin
   logical, intent(in) :: singleprec

   real(doub_prec), dimension(ntotal) :: dattempd
!
!--read several scalar quantities
!
   if (singleprec) then
      read (iunit,iostat=ierr) datin(1:ntotal)
   else
      read (iunit,iostat=ierr) dattempd(1:ntotal)
      datin(1:ntotal) = real(dattempd(1:ntotal))
   endif

 end subroutine readcol

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
 use labels
 use params
 use settings_data, only:ndim,ndimV,ncolumns,iformat,ntypes, &
                    UseTypeInRenderings
 use geometry, only:labelcoord
 implicit none
 integer :: i

 if (ndim.le.0 .or. ndim.gt.3) then
    print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
    return
 endif
 if (ndimV.le.0 .or. ndimV.gt.3) then
    print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
    return
 endif

 do i=1,ndim
    ix(i) = i
 enddo
 ivx = ndim + 1
 ih = ndim + ndimV + 1        !  smoothing length
 irho = ndim + ndimV + 2      ! location of rho in data array
 iutherm = ndim + ndimV + 3   !  thermal energy
 ipmass = ndim + ndimV + 4    !  particle mass

 label(ix(1:ndim)) = labelcoord(1:ndim,1)
 !
 !--label vector quantities (e.g. velocity) appropriately
 !
 iamvec(ivx:ivx+ndimV-1) = ivx
 labelvec(ivx:ivx+ndimV-1) = 'v'
 do i=1,ndimV
    label(ivx+i-1) = trim(labelvec(ivx+i-1))//'\d'//labelcoord(i,1)
 enddo

 label(irho) = '\gr'
 label(iutherm) = 'u'
 label(ih) = 'h       '
 label(ipmass) = 'particle mass'
 label(ndim + ndimV+5) = '\ga'
 label(ndim + ndimV+6) = '\ga\du'
 if (iformat.eq.2) then

    !
    !--mag field (vector)
    !
    label(ndim + ndimV+7) = '\ga\dB'
    iBfirst = ndim + ndimV+7+1        ! location of Bx
    iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
    labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
    do i=1,ndimV
       label(iBfirst+i-1) = trim(labelvec(iBfirst))//'\d'//labelcoord(i,1) !' (x10\u-3\d)' !//'/rho'
    enddo
    !
    !--more scalars
    !
    label(ndim+2*ndimV+8) = 'psi'

    ipr = ndim + 2*ndimV + 9 !  pressure
    label(ipr) = 'P'
    label(ndim+2*ndimV+10) = 'div v'

    idivB = ndim+2*ndimV+11
    label(idivB) = 'div B'
    !
    !--current density (vector)
    !
    iJfirst = ndim+2*ndimV+11+1
    iamvec(iJfirst:iJfirst+ndimV-1) = iJfirst
    labelvec(iJfirst:iJfirst+ndimV-1) = 'J'
    do i=1,ndimV
       label(iJfirst+i-1) = trim(labelvec(iJfirst))//labelcoord(i,1)
    enddo
 else
    ipr = ndim + ndimV + 7 !  pressure
    label(ipr) = 'P'
    label(ndim+ndimV+8) = 'grad h'
    label(ndim+ndimV+9) = 'grad soft'
    label(ndim+ndimV+10) = 'phi'
    label(ndim+ndimV+11) = 'f_grav'
!    label(ndim+ndimV+8) = 'div v'
!    label(ndim+ndimV+9) = 'grad h'
    if (iformat.eq.3) then
       !!!irho = ndim+ndimV+9
       label(ndim+ndimV+9) = 'rho*'
       label(ndim+ndimV+10) = 'sqrt g'
       iamvec(ndim+ndimV+11:ndim+ndimV+10+ndimV) = ndim+ndimV+11
       labelvec(ndim+ndimV+11:ndim+ndimV+10+ndimV) = 'pmom'
       do i=1,ndimV
          label(ndim+ndimV+10+i) = labelvec(ndim+ndimV+11)//labelcoord(i,1)
       enddo
    endif
    iBfirst = 0
 endif

 if (ncolumns.gt.ndim+3*ndimV+11) then
    label(ndim+3*ndimV+12) = 'f_visc_x'
    label(ndim+3*ndimV+13) = 'f_visc_y'
    label(ndim+3*ndimV+14) = 'f_x'
    label(ndim+3*ndimV+15) = 'f_y'
 endif
!
!--these are here for backwards compatibility -- could be removed
!  if (ncolumns.gt.ndim+3*ndimV+7) then
!     label(ndim + 3*ndimV+8) = 'v_parallel'
!     label(ndim + 3*ndimV+9) = 'v_perp'
!     label(ndim + 3*ndimV+10) = 'B_parallel'
!     label(ndim + 3*ndimV+11) = 'B_perp'
!  endif

!
!--set labels for each type of particles
!
 ntypes = 2
 labeltype(1) = 'gas'
 labeltype(2) = 'ghost'
 UseTypeInRenderings(1) = .true.
 UseTypeInRenderings(2) = .true.

!-----------------------------------------------------------

 return
end subroutine set_labels
