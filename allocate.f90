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
  integer, dimension(:), allocatable :: nparttemp,ntottemp,nghosttemp
  integer, dimension(:), allocatable :: ntotplottemp
  integer, dimension(:,:), allocatable :: iamtemp
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
!
!--if re-allocating, copy arrays to temporary versions
!
  ierr = 0
  if (allocated(dat)) then
     reallocate = .true.
     print*,' reallocating memory: ntot = ', &
          npartin,' nstep = ',nstep,' ncol = ',ncolumns
     allocate(dattemp(maxcol,maxpart,maxstep), stat=ierr)
     allocate(iamtemp(maxpart,maxstep))
     if (ierr.ne.0) print*,'error allocating memory'
     dattemp = dat
     iamtemp = iam
     deallocate(dat,iam)

     allocate(nparttemp(maxstep),ntottemp(maxstep),nghosttemp(maxstep))
     allocate(ntotplottemp(maxstep))
     nparttemp = npart
     ntottemp = ntot
     nghosttemp = nghost
     ntotplottemp = ntotplot
     deallocate(npart,ntot,nghost,ntotplot)

     allocate(timetemp(maxstep),gammatemp(maxstep))
     timetemp = time
     gammatemp = gamma
     deallocate(time,gamma)

  else
     print*,' allocating memory: ntot = ', &
          npartin,' nstep = ',nstep,' ncol = ',ncolumns
  endif
!
!--save array sizes
!
  if (npartin.lt.maxpart) print*,' WARNING: # particles < previous in allocate'
  if (nstep.lt.maxstep) print*,' WARNING: # steps < previous in allocate'
  if (ncolumns.lt.maxcol) print*,' WARNING: # columns < previous in allocate'
  maxpartold = min(maxpart,npartin)
  maxstepold = min(maxstep,nstep)
  maxcolold = min(maxcol,ncolumns)

  maxpart = npartin
  maxstep = nstep
  maxcol = ncolumns

!
!--main data array
!
  allocate(dat(maxcol,maxpart,maxstep), stat=ierr)
  allocate(iam(maxpart,maxstep))
  if (ierr.ne.0) stop 'error allocating memory for dat array'
  if (reallocate) then
     dat(1:maxcolold,1:maxpartold,1:maxstepold) = dattemp(1:maxcolold,1:maxpartold,1:maxstepold)
     iam(1:maxpartold,1:maxstepold) = iamtemp(1:maxpartold,1:maxstepold)
  endif
!
!--other arrays
!
  allocate(npart(maxstep),ntot(maxstep),nghost(maxstep),stat=ierr)
  allocate(ntotplot(maxstep),stat=ierr)
  allocate(time(maxstep),gamma(maxstep),stat=ierr)
  if (ierr.ne.0) stop 'error allocating memory for header arrays'
  if (reallocate) then
     npart(1:maxstepold) = nparttemp(1:maxstepold)
     ntot(1:maxstepold) = ntottemp(1:maxstepold)
     nghost(1:maxstepold) = nghosttemp(1:maxstepold)
     ntotplot(1:maxstepold) = ntotplottemp(1:maxstepold)
     time(1:maxstepold) = timetemp(1:maxstepold)
     gamma(1:maxstepold) = gammatemp(1:maxstepold)
     deallocate(nparttemp,ntottemp,nghosttemp,ntotplottemp)
     deallocate(timetemp,gammatemp)
  endif

  return 
end subroutine alloc

