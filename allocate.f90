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
subroutine alloc(npartin,nstep,ncolumnsin)
  use particle_data
  implicit none
  integer, intent(in) :: npartin,nstep,ncolumnsin
  integer :: maxpartold,maxstepold,maxcolold
  integer :: ierr,ncolumns
  logical :: reallocate,reallocate_part,reallocate_step
  integer, dimension(:), allocatable :: icolourmetemp
  integer, dimension(:,:), allocatable :: npartoftypetemp
  real, dimension(:,:), allocatable :: massoftypetemp
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
  if (ncolumnsin.le.0) then
     print*,'allocate: error in input, ncolumns = ',ncolumnsin
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
     endif

     dattemp = dat
     deallocate(dat)

     if (reallocate_part) then
        icolourmetemp(1:maxpartold) = icolourme(1:maxpartold)
        deallocate(icolourme)
     endif

     if (reallocate_step) then
        allocate(npartoftypetemp(maxparttypes,maxstep),stat=ierr)
        if (ierr /= 0) stop 'error allocating memory (npartoftypetemp)'
        npartoftypetemp = npartoftype
        deallocate(npartoftype)

        allocate(massoftypetemp(maxparttypes,maxstep),stat=ierr)
        if (ierr /= 0) stop 'error allocating memory (npartoftypetemp)'
        massoftypetemp = massoftype
        deallocate(massoftype)

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

     allocate(massoftype(maxparttypes,maxstep),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for header arrays'

     allocate(time(maxstep),gamma(maxstep),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for header arrays'

     npartoftype = 0
     massoftype = 0.
     time = -huge(time) ! initialise like this so we know if has not been read
     gamma = 0.

     if (reallocate_step) then
        npartoftype(:,1:maxstepold) = npartoftypetemp(:,1:maxstepold)
        massoftype(:,1:maxstepold) = massoftypetemp(:,1:maxstepold)
        time(1:maxstepold) = timetemp(1:maxstepold)
        gamma(1:maxstepold) = gammatemp(1:maxstepold)
        deallocate(npartoftypetemp,massoftypetemp)
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
 use particle_data
 implicit none
 
 if (allocated(dat)) deallocate(dat)
 if (allocated(icolourme)) deallocate(icolourme)
 if (allocated(npartoftype)) deallocate(npartoftype)
 if (allocated(massoftype)) deallocate(massoftype)
 if (allocated(time)) deallocate(time)
 if (allocated(gamma)) deallocate(gamma)
 
 return
end subroutine deallocate_all

end module mem_allocation
