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
!  Copyright (C) 2005-2012 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

module mem_allocation
 implicit none

contains
!----------------------------------------------------------------------------
!
!  memory allocation/reallocation for main data arrays
!
!  the global parameters maxpart, maxstep and maxcol are set to
!  the dimensions allocated
!
!----------------------------------------------------------------------------
subroutine alloc(npartin,nstep,ncolumnsin,mixedtypes)
  use particle_data
  implicit none
  integer, intent(in) :: npartin,nstep,ncolumnsin
  logical, intent(in), optional :: mixedtypes
  integer :: maxpartold,maxstepold,maxcolold
  integer :: ierr,ncolumns
  logical :: reallocate,reallocate_part,reallocate_step,reallocate_itype
  integer, dimension(:), allocatable :: icolourmetemp
  integer(kind=int1), dimension(:,:), allocatable :: iamtypetemp
  integer, dimension(:,:), allocatable :: npartoftypetemp
  real, dimension(:,:), allocatable :: masstypetemp
  real, dimension(:), allocatable :: timetemp, gammatemp
  real, dimension(:,:,:), allocatable :: dattemp
!
!--check for errors in input
!
  if (npartin.le.0) then
     print*,'allocate: error in input, npartin = ',npartin
     return
  endif
  if (nstep.le.0) then
     print*,'allocate: error in input, nstep = ',nstep
     return
  endif
  if (ncolumnsin.lt.0) then
     print*,'allocate: error in input, ncolumns = ',ncolumnsin
     return
  elseif (ncolumnsin.eq.0) then
     print*,'WARNING: allocate: ncolumns = 0 in input'
  endif
  !--do nothing if array sizes are the same
  if (npartin.eq.maxpart .and. ncolumnsin.eq.maxcol .and. nstep.eq.maxstep) then
     return
  endif

!
!--save array sizes
!
  if (npartin.lt.maxpart) print "(a)",' WARNING: # particles < previous in allocate'
  if (nstep.lt.maxstep) print "(a)",' WARNING: # steps < previous in allocate'
!--at the moment ncolumns cannot be decreased (due to calc_quantities)
  if (ncolumnsin.lt.maxcol) then
     ncolumns = maxcol
     !!print "(a)",' WARNING: # columns < previous in allocate'
  else
     ncolumns = ncolumnsin
  endif

  maxpartold = min(maxpart,npartin)
  maxstepold = min(maxstep,nstep)
  maxcolold = min(maxcol,ncolumns)

  reallocate = .false.
  reallocate_part = .false.
  reallocate_step = .false.
  reallocate_itype = .false.
!
!--if re-allocating, copy arrays to temporary versions
!
  ierr = 0
  if (allocated(dat)) then
     reallocate = .true.
     if (maxpart.ne.npartin) reallocate_part = .true.
     if (maxstep.ne.nstep) reallocate_step = .true.

     print 10,'> reallocating memory:',npartin,nstep,ncolumns
10   format (a,' parts = ',i10,' steps = ',i6,' cols = ',i4)
     allocate(dattemp(maxpartold,maxcolold,maxstepold), stat=ierr)
     if (ierr /= 0) stop 'error allocating memory (dattemp)'

     if (reallocate_part) then
        allocate(icolourmetemp(maxpartold),stat=ierr)
        if (ierr /= 0) stop 'error allocating memory (icolourmetemp)'
        icolourmetemp(1:maxpartold) = icolourme(1:maxpartold)
        deallocate(icolourme)
     endif

     dattemp = dat
     deallocate(dat)

     if (allocated(iamtype)) then ! should always be true
        !--if iamtype has meaningful contents and reallocation is necessary
        reallocate_itype = (reallocate_part .or. reallocate_step) .and. (size(iamtype(:,1)).eq.maxpartold)
        if (reallocate_itype) then
           allocate(iamtypetemp(maxpartold,maxstepold), stat=ierr)
           if (ierr /= 0) stop 'error allocating memory (iamtypetemp)'
           iamtypetemp(1:maxpartold,1:maxstepold) = iamtype(1:maxpartold,1:maxstepold)
           deallocate(iamtype)

        elseif (present(mixedtypes)) then
        !--if iamtype has size 1 or 0 but should be allocated here,
        !  deallocate so we can give it correct size
           if (mixedtypes .and. size(iamtype(:,1)).lt.maxpart) deallocate(iamtype)
        endif
     endif

     if (reallocate_step) then
        allocate(npartoftypetemp(maxparttypes,maxstep),stat=ierr)
        if (ierr /= 0) stop 'error allocating memory (npartoftypetemp)'
        npartoftypetemp = npartoftype
        deallocate(npartoftype)

        allocate(masstypetemp(maxparttypes,maxstep),stat=ierr)
        if (ierr /= 0) stop 'error allocating memory (npartoftypetemp)'
        masstypetemp = masstype
        deallocate(masstype)

        allocate(timetemp(maxstep),gammatemp(maxstep),stat=ierr)
        if (ierr /= 0) stop 'error allocating memory (timetemp,gammatemp)'
        timetemp = time
        gammatemp = gamma
        deallocate(time,gamma)
     endif

  else
     print 10,'> allocating memory:',npartin,nstep,ncolumns
     maxpart = npartin
     maxstep = nstep
     maxcol = ncolumns
  endif


  maxpart = npartin
  maxstep = nstep
  maxcol = ncolumns

!
!--main data array
!
  allocate(dat(maxpart,maxcol,maxstep), stat=ierr)
  if (ierr /= 0) then
     print*,' parts = ',maxpart,' columns = ',maxcol,' steps = ',maxstep
     stop 'error allocating memory for dat array'
  endif

  if (reallocate) then
     dat(1:maxpartold,1:maxcolold,1:maxstepold) = dattemp(1:maxpartold,1:maxcolold,1:maxstepold)
     deallocate(dattemp)
  else
     dat = 0.
  endif
!
!--type array if necessary
!
  if (present(mixedtypes)) then
     if (mixedtypes .and. .not.allocated(iamtype)) then
        allocate(iamtype(maxpart,maxstep), stat=ierr)
        if (ierr /= 0) stop 'error allocating memory for type array'
        iamtype = 1
        !--copy contents if reallocating
        if (reallocate_itype) then
           iamtype(1:maxpartold,1:maxstepold) = iamtypetemp(1:maxpartold,1:maxstepold)
           deallocate(iamtypetemp)
        endif
     elseif (.not.mixedtypes) then
        !--if called with mixedtypes explictly false, deallocate itype array
        if (allocated(iamtype)) deallocate(iamtype)
     endif
  elseif (reallocate_itype) then
     !--if called without mixedtypes, preserve contents of itype array
     allocate(iamtype(maxpart,maxstep), stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for type array'
     iamtype = 1

     iamtype(1:maxpartold,1:maxstepold) = iamtypetemp(1:maxpartold,1:maxstepold)
     deallocate(iamtypetemp)
  endif
  !--make sure iamtype is always allocated for safety, just with size=1 if not used
  if (.not.allocated(iamtype)) then
     allocate(iamtype(1,maxstep),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for type array (1)'
  endif
!
!--particle arrays
!
  if (.not.allocated(icolourme) .or. reallocate_part) then
     allocate(icolourme(maxpart),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for icolourme array'
     icolourme = 1

     if (reallocate_part) then
        icolourme(1:maxpartold) = icolourmetemp(1:maxpartold)
        deallocate(icolourmetemp)
     endif
  endif
!
!--other arrays
!
  if (.not.allocated(npartoftype)) then
     allocate(npartoftype(maxparttypes,maxstep),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for header arrays'

     allocate(masstype(maxparttypes,maxstep),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for header arrays'

     allocate(time(maxstep),gamma(maxstep),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for header arrays'

     npartoftype = 0
     masstype = 0.
     time = time_not_read_val ! initialise like this so we know if has not been read
     gamma = 0.

     if (reallocate_step) then
        npartoftype(:,1:maxstepold) = npartoftypetemp(:,1:maxstepold)
        masstype(:,1:maxstepold) = masstypetemp(:,1:maxstepold)
        time(1:maxstepold) = timetemp(1:maxstepold)
        gamma(1:maxstepold) = gammatemp(1:maxstepold)
        deallocate(npartoftypetemp,masstypetemp)
        deallocate(timetemp,gammatemp)
     endif
  endif

  return
end subroutine alloc


!-----------------------------------------
!
!  deallocation of remaining memory
!  (for tidiness - not strictly necessary)
!
!-----------------------------------------
subroutine deallocate_all
 use particle_data, only:dat,icolourme,iamtype,npartoftype,masstype,time,gamma
 implicit none

 if (allocated(dat)) deallocate(dat)
 if (allocated(icolourme)) deallocate(icolourme)
 if (allocated(iamtype)) deallocate(iamtype)
 if (allocated(npartoftype)) deallocate(npartoftype)
 if (allocated(masstype)) deallocate(masstype)
 if (allocated(time)) deallocate(time)
 if (allocated(gamma)) deallocate(gamma)

 return
end subroutine deallocate_all

end module mem_allocation
