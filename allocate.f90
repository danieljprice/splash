!----------------------------------------------------------------------------
!
!  memory allocation/reallocation for main data arrays
!
!  the global parameters maxpart, maxstep and maxcol are set to 
!  the dimensions allocated
!
!----------------------------------------------------------------------------

subroutine alloc(npartin,nstep,ncolumns)
  use particle_data
  implicit none
  integer, intent(in) :: npartin,nstep,ncolumns
  integer :: maxpartold,maxstepold,maxcolold
  integer :: ierr
  logical :: reallocate
  integer, dimension(:), allocatable :: ntottemp, icolourmetemp
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
  if (ncolumns.le.0) then
     print*,'allocate: error in input, ncolumns = ',ncolumns
     return
  endif

  reallocate = .false.
!
!--if re-allocating, copy arrays to temporary versions
!
  ierr = 0
  if (allocated(dat)) then
     reallocate = .true.
     print 10,'> reallocating memory: ',npartin,nstep,ncolumns
10   format (a,' ntot = ',i10,' nstep = ',i6,' ncol = ',i4)
     allocate(dattemp(maxpart,maxcol,maxstep), stat=ierr)
     if (ierr /= 0) stop 'error allocating memory (dattemp)'
     allocate(iamtemp(maxpart,maxstep),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory (iamtemp)'
     allocate(icolourmetemp(maxpart),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory (icolourmetemp)'

     dattemp = dat
     iamtemp = iam
     icolourmetemp = icolourme
     deallocate(dat,iam,icolourme)

     allocate(ntottemp(maxstep),npartoftypetemp(maxparttypes,maxstep),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory'
     ntottemp = ntot
     npartoftypetemp = npartoftype
     deallocate(ntot,npartoftype)

     allocate(timetemp(maxstep),gammatemp(maxstep),stat=ierr)
     if (ierr /= 0) stop 'error allocating memory'
     timetemp = time
     gammatemp = gamma
     deallocate(time,gamma)

  else
     print 10,'> allocating memory: ',npartin,nstep,ncolumns
     maxpart = npartin
     maxstep = nstep
     maxcol = ncolumns
  endif
!
!--save array sizes
!
  if (npartin.lt.maxpart) print "(a)",' WARNING: # particles < previous in allocate'
  if (nstep.lt.maxstep) print "(a)",' WARNING: # steps < previous in allocate'
  if (ncolumns.lt.maxcol) print "(a)",' WARNING: # columns < previous in allocate'
  maxpartold = min(maxpart,npartin)
  maxstepold = min(maxstep,nstep)
  maxcolold = min(maxcol,ncolumns)

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
  
  allocate(icolourme(maxpart),stat=ierr)
  if (ierr /= 0) stop 'error allocating memory for icolourme array'
  icolourme = 1
  
  if (reallocate) then
     dat(1:maxpartold,1:maxcolold,1:maxstepold) = dattemp(1:maxpartold,1:maxcolold,1:maxstepold)
     iam(1:maxpartold,1:maxstepold) = iamtemp(1:maxpartold,1:maxstepold)
     icolourme(1:maxpartold) = icolourmetemp(1:maxpartold)
  endif
!
!--other arrays
!
  allocate(ntot(maxstep),stat=ierr)
  if (ierr /= 0) stop 'error allocating memory for header arrays'

  allocate(npartoftype(maxparttypes,maxstep),stat=ierr)
  if (ierr /= 0) stop 'error allocating memory for header arrays'

  allocate(time(maxstep),gamma(maxstep),stat=ierr)
  if (ierr /= 0) stop 'error allocating memory for header arrays'

  if (reallocate) then
     ntot(1:maxstepold) = ntottemp(1:maxstepold)
     npartoftype(:,1:maxstepold) = npartoftypetemp(:,1:maxstepold)
     time(1:maxstepold) = timetemp(1:maxstepold)
     gamma(1:maxstepold) = gammatemp(1:maxstepold)
     deallocate(ntottemp,npartoftypetemp)
     deallocate(timetemp,gammatemp)
  endif

  return 
end subroutine alloc

