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
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! the data is stored in the global array dat
!
! THIS VERSION FOR DAN'S SPMHD CODE (BINARY DUMPS)
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
  use exact,          only:hfact
  use particle_data,  only:npartoftype,time,gamma,dat,maxpart,maxstep,maxcol,iamtype
  use params
  use filenames,      only:nfiles
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc,icoords,iformat, &
                          buffer_data,iverbose,debugmode
  use mem_allocation, only:alloc
  use geometry,       only:labelcoordsys
  use system_utils,   only:lenvironment
  use labels,         only:labeltype,print_types
  implicit none
  integer,          intent(in)  :: indexstart,ipos
  integer,          intent(out) :: nstepsread
  character(len=*), intent(in)  :: rootname
  character(len=len(rootname)+4) :: datfile
  integer :: i,icol,ierr,iunit,ilen,j,ilast
  integer :: ncol_max,ndim_max,npart_max,ndimV_max,nstep_max
  integer :: npartin,ntotin,ncolstep,nparti,ntoti
  integer, dimension(3) :: ibound
  logical :: reallocate, singleprecision

  real :: timein, gammain, hfactin
  real, dimension(3) :: xmin, xmax
  real(doub_prec) :: timeind,gammaind,hfactind
  real(doub_prec), dimension(3) :: xmind, xmaxd
  real(doub_prec), dimension(:), allocatable :: dattempd
  integer, dimension(:), allocatable :: itype
  character(len=20) :: geomfile

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

  if (iverbose.ge.1) print "(1x,a)",'reading ndspmhd format'
  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  !--open data file and read data
  !
  open(unit=iunit,iostat=ierr,file=datfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print*,' *** Error opening '//trim(datfile)//' ***'
     return
  endif
!
!--read first header line
!
  singleprecision = .false.
  read(iunit,iostat=ierr,end=80) timeind,npartin,ntotin,gammaind, &
       hfactind,ndim_max,ndimV_max,ncol_max,iformat
!  print*,'time = ',timeind,' hfact = ',hfactind,' ndim=',ndim_max,'ncol=',ncol_max
!  print*,'npart = ',npartin,ntotin,geomfile
  if (ierr /= 0 .or. ndim_max.le.0 .or. ndim_max.gt.3 &
     .or. ndimV_max.le.0 .or. ndimV_max.gt.3 &
     .or. ncol_max.le.0 .or. ncol_max.gt.100 &
     .or. npartin.le.0 .or. npartin.gt.1e7 .or. ntotin.le.0 .or. ntotin.gt.1e7 &
     .or. iformat.lt.0 .or. iformat.gt.10) then
     !
     !--try single precision
     !
     rewind(iunit)
     read(iunit,iostat=ierr,end=80) timein,npartin,ntotin,gammain, &
         hfactin,ndim_max,ndimV_max,ncol_max,iformat
     singleprecision = .true.
     if (ierr /= 0 .or. ndim_max.le.0 .or. ndim_max.gt.3 &
        .or. ndimV_max.le.0 .or. ndimV_max.gt.3 &
        .or. ncol_max.le.0 .or. ncol_max.gt.100 &
        .or. npartin.le.0 .or. npartin.gt.1e7 .or. ntotin.le.0 .or. ntotin.gt.1e7 &
        .or. iformat.lt.0 .or. iformat.gt.10) then

        print "(a)",' *** Error reading first header ***'
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
     call alloc(npart_max,nstep_max,ncol_max+ncalc,mixedtypes=.true.)
  endif
!
!--rewind file
!
  rewind(iunit)

  i = indexstart
  nstepsread = 0

     reallocate = .false.
     npart_max = maxpart
     nstep_max = maxstep
     geomfile = ' '
     !
     !--read header line for this timestep
     !
     if (singleprecision) then
        if (debugmode) print "(a)",'DEBUG: single precision dump'
        read(iunit,iostat=ierr) timein,nparti,ntoti,gammain, &
          hfactin,ndim,ndimV,ncolstep,iformat,ibound(1:ndim), &
          xmin(1:ndim),xmax(1:ndim),ilen,geomfile(1:ilen)
     else
        if (debugmode) print "(a)",'DEBUG: double precision dump'
        read(iunit,iostat=ierr) timeind,nparti,ntoti,gammaind, &
          hfactind,ndim,ndimV,ncolstep,iformat,ibound(1:ndim), &
          xmind(1:ndim),xmaxd(1:ndim),ilen,geomfile(1:ilen)
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
     npartoftype(3,i) = ntoti - nparti
     if (iverbose.ge.1) then
        print "(a14,':',es10.3,a6,':',i8,a8,':',i8)",' time',time(i),'npart',nparti,'ntotal',ntoti
        print "(a14,':',i8,a8,':',f8.4,a8,':',f8.4)",' ncolumns',ncolstep,'gamma',gamma(i),'hfact',hfact
        print "(a14,':',i8,a8,':',i8)",'ndim',ndim,'ndimV',ndimV
     else
        print "(1x,a,':',es10.3,a8,':',i8,a8,':',i8)",'time',time(i),'npart',nparti,'ntotal',ntoti
     endif
     select case(geomfile(1:6))
     case('cylrpz')
        icoords = 2
     case('sphrpt')
        icoords = 3
     case default
        icoords = 1
     end select
     if (icoords.ne.1) print "(a14,a)",' geometry: ',trim(geomfile)//' ('//trim(labelcoordsys(icoords))//')'
     if (iverbose.ge.1 .and. any(ibound(1:ndim).ne.0)) then
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
        if (iverbose.ge.1) print*,'ncolumns = ',ncolumns,ncol_max
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
        call alloc(npart_max,nstep_max,ncol_max+ncalc,mixedtypes=.true.)
     endif


     if (ntoti.gt.0) then
        if (.not.singleprecision) allocate(dattempd(ntoti))
        do icol=1,ncolstep
           if (singleprecision) then
              read (iunit,iostat=ierr,end=67) dat(1:ntoti,icol,i)
           else
              read (iunit,iostat=ierr,end=67) dattempd(1:ntoti)
              dat(1:ntoti,icol,i) = real(dattempd(1:ntoti))
           endif
           if (ierr /= 0) print "(a,i2,a)",'*** error reading column ',icol,' ***'
        enddo
        if (allocated(dattempd)) deallocate(dattempd)

        allocate(itype(ntoti))
        read(iunit,iostat=ierr) itype(1:ntoti)
        if (ierr.ne.0) then
           if (debugmode) print "(a)",'DEBUG: itype not found in dump file'
           iamtype(1:nparti,i) = 1
           iamtype(nparti+1:ntoti,i) = 3
        else
           !
           !--assign SPLASH types from ndspmhd types
           !
           npartoftype(:,i) = 0
           do j=1,ntoti
              if (j.gt.nparti) then
                 iamtype(j,i) = 3
                 npartoftype(3,i) = npartoftype(3,i) + 1
              elseif (itype(j).eq.2) then
                 iamtype(j,i) = 2
                 npartoftype(2,i) = npartoftype(2,i) + 1
              else
                 iamtype(j,i) = 1
                 npartoftype(1,i) = npartoftype(1,i) + 1
              endif
           enddo
        endif
        if (allocated(itype)) deallocate(itype)
     else
        npartoftype(1,i) = 1
        npartoftype(2:3,i) = 0
        dat(:,:,i) = 0.
     endif

   goto 68
67 continue
   print "(a)",' > end of file reached <'
68 continue
  !
  !--close data file and return
  !
close(unit=11)

ilast = i
!
!--ONE FLUID DUST: FAKE IT AS IF IT IS TWO FLUIDS
!  (copy the particles, then copy gas properties onto first lot, then dust properties onto second lot)
!
ncolumns = ncol_max
ndim = ndim_max
ndimV = ndimV_max

call set_labels

!if (iformat.eq.5 .and. .not.lenvironment('NSPLASH_BARYCENTRIC')) then
!   call fake_twofluids
!   iformat = 1
!endif

if (any(npartoftype(2:,ilast) > 0)) call print_types(npartoftype(:,ilast),labeltype)

if (debugmode) print*,'DEBUG> Read steps ',indexstart,'->',indexstart + nstepsread - 1, &
       ' last step ntot = ',sum(npartoftype(:,indexstart+nstepsread-1))
return

80 continue
print*,' *** data file empty : no timesteps ***'
return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
 use labels, only:ix,ivx,ih,irho,iutherm,ipmass,ipr,iBfirst, &
             idivB,iJfirst,iamvec,labelvec,label,labeltype, &
             irhorestframe,idustfrac,ideltav
 use params
 use settings_data, only:ndim,ndimV,iformat,ntypes, &
                    UseTypeInRenderings,ncolumns
 use geometry, only:labelcoord
 implicit none
 integer :: i,icol

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

 label(irho) = '\gr'
 label(iutherm) = 'u'
 label(ih) = 'h       '
 label(ipmass) = 'particle mass'
 label(ndim + ndimV+5) = '\ga'
 label(ndim + ndimV+6) = '\ga\du'
 icol = ndim+ndimV + 7
 if (iformat.eq.2 .or. iformat.eq.4) then
    !
    !--mag field (vector)
    !
    label(icol) = '\ga\dB'
    iBfirst = icol+1        ! location of Bx
    iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
    labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
    icol = icol + ndimV
    !
    !--more scalars
    !
    icol = icol + 1
    label(icol) = 'psi'

    icol = icol + 1
    ipr = icol !  pressure
    label(ipr) = 'P'

    icol = icol + 1
    label(icol) = 'div v'

    icol = icol + 1
    idivB = icol
    label(idivB) = 'div B'
    !
    !--current density (vector)
    !
    iJfirst = icol + 1
    iamvec(icol+1:icol+ndimV) = icol + 1
    labelvec(icol+1:icol+ndimV) = 'J'
    icol = icol + ndimV

    icol = icol + 1
    label(icol) = 'grad h'

    iamvec(icol+1:icol+ndimV) = icol + 1
    labelvec(icol+1:icol+ndimV) = 'force'
    icol = icol + ndimV

    iamvec(icol+1:icol+ndimV) = icol + 1
    labelvec(icol+1:icol+ndimV) = 'A'
    icol = icol + ndimV

 else
    ipr = icol !  pressure
    label(ipr) = 'P'
    icol = icol + 1
    label(icol) = 'div v'
    icol = icol + 1
    label(icol) = 'grad h'

    iamvec(icol+1:icol+ndimV) = icol + 1
    labelvec(icol+1:icol+ndimV) = 'force'
    icol = icol + ndimV
    iBfirst = 0

    if (ncolumns > 20) then
       iamvec(icol+1:icol+ndimV) = icol + 1
       labelvec(icol+1:icol+ndimV) = 'del^{2} v'
       icol = icol + ndimV

       iamvec(icol+1:icol+ndimV) = icol + 1
       labelvec(icol+1:icol+ndimV) = 'grad (div {\bf v})'
       icol = icol + ndimV
    endif
 endif
 if (iformat.eq.5) then
    icol = icol + 1
    label(icol) = 'Dust fraction'
    idustfrac = icol
    iamvec(icol+1:icol+ndimV) = icol + 1
    labelvec(icol+1:icol+ndimV) = '\Deltav'
    ideltav = icol + 1
    icol = icol + ndimV
 elseif (iformat.gt.2) then
    irhorestframe = irho
    icol = icol + 1
    irho = icol
    label(icol) = 'rho*'
    !irho = icol
    icol = icol + 1
    label(icol) = 'sqrt g'
    iamvec(icol+1:icol+ndimV) = icol + 1
    labelvec(icol+1:icol+ndimV) = 'pmom'
    icol = icol + ndimV
 endif
!
!--set labels for each type of particles
!
 ntypes = 3
 labeltype(1) = 'gas'
 labeltype(2) = 'dust'
 labeltype(3) = 'ghost'
 UseTypeInRenderings(1) = .true.
 UseTypeInRenderings(2) = .true.
 UseTypeInRenderings(3) = .true.

!-----------------------------------------------------------

 return
end subroutine set_labels
