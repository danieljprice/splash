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
  use settings_part, only:icircpart
  implicit none
  integer, intent(in) :: npartin,nstep,ncolumnsin
  integer :: maxpartold,maxstepold,maxcolold
  integer :: ierr,ncolumns
  logical :: reallocate,reallocate_part,reallocate_step
  integer, dimension(:), allocatable :: ntottemp,icolourmetemp,icircparttemp
  integer, dimension(:,:), allocatable :: iamtemp, npartoftypetemp
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
     
     print 10,'> reallocating memory: ',npartin,nstep,ncolumns
10   format (a,' ntot = ',i10,' nstep = ',i6,' ncol = ',i4)
     allocate(dattemp(maxpartold,maxcolold,maxstepold), stat=ierr)
     if (ierr /= 0) stop 'error allocating memory (dattemp)'
     allocate(iamtemp(maxpartold,maxstepold),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory (iamtemp)'
     
     if (reallocate_part) then
        allocate(icolourmetemp(maxpartold),stat=ierr)
        if (ierr /= 0) stop 'error allocating memory (icolourmetemp)'
        allocate(icircparttemp(maxpartold),stat=ierr)
        if (ierr /= 0) stop 'error allocating memory (icircparttemp)'
     endif

     dattemp = dat
     iamtemp = iam
     deallocate(dat,iam)

     if (reallocate_part) then
        icolourmetemp(1:maxpartold) = icolourme(1:maxpartold)
        icircparttemp(1:maxpartold) = icircpart(1:maxpartold)
        deallocate(icolourme,icircpart)
     endif

     if (reallocate_step) then
        allocate(ntottemp(maxstep),npartoftypetemp(maxparttypes,maxstep),stat=ierr)
        if (ierr /= 0) stop 'error allocating memory (ntottemp)'
        ntottemp = ntot
        npartoftypetemp = npartoftype
        deallocate(ntot,npartoftype)

        allocate(timetemp(maxstep),gammatemp(maxstep),stat=ierr)
        if (ierr /= 0) stop 'error allocating memory (timetemp,gammatemp)'
        timetemp = time
        gammatemp = gamma
        deallocate(time,gamma)
     endif

  else
     print 10,'> allocating memory: ',npartin,nstep,ncolumns
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
  if (ierr /= 0) stop 'error allocating memory for dat array'

  allocate(iam(maxpart,maxstep), stat=ierr)
  if (ierr /= 0) stop 'error allocating memory for iam array'
  iam = 0
    
  if (reallocate) then
     dat(1:maxpartold,1:maxcolold,1:maxstepold) = dattemp(1:maxpartold,1:maxcolold,1:maxstepold)
     iam(1:maxpartold,1:maxstepold) = iamtemp(1:maxpartold,1:maxstepold)
  endif
!
!--particle arrays
!
  if (.not.allocated(icolourme) .or. reallocate_part) then
     allocate(icolourme(maxpart),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for icolourme array'
     icolourme = 1

     allocate(icircpart(maxpart),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for icolourme array'
     icircpart = 0
     
     if (reallocate_part) then
        icolourme(1:maxpartold) = icolourmetemp(1:maxpartold)
        icircpart(1:maxpartold) = icircparttemp(1:maxpartold)
        deallocate(icolourmetemp,icircparttemp)
     endif
  endif
!
!--other arrays
!
  if (.not.allocated(npartoftype)) then
     allocate(ntot(maxstep),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for header arrays'

     allocate(npartoftype(maxparttypes,maxstep),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for header arrays'

     allocate(time(maxstep),gamma(maxstep),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory for header arrays'

     ntot = 0
     npartoftype = 0
     time = 0.
     gamma = 0.

     if (reallocate_step) then
        ntot(1:maxstepold) = ntottemp(1:maxstepold)
        npartoftype(:,1:maxstepold) = npartoftypetemp(:,1:maxstepold)
        time(1:maxstepold) = timetemp(1:maxstepold)
        gamma(1:maxstepold) = gammatemp(1:maxstepold)
        deallocate(ntottemp,npartoftypetemp)
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
 use settings_part, only:icircpart
 implicit none
 
 if (allocated(dat)) deallocate(dat)
 if (allocated(iam)) deallocate(iam)
 if (allocated(icolourme)) deallocate(icolourme)
 if (allocated(icircpart)) deallocate(icircpart)
 if (allocated(ntot)) deallocate(ntot)
 if (allocated(npartoftype)) deallocate(npartoftype)
 if (allocated(time)) deallocate(time)
 if (allocated(gamma)) deallocate(gamma)
 
 return
end subroutine deallocate_all

end module mem_allocation
