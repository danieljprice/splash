!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! the data is stored in the global array dat
!
! THIS VERSION FOR DAN'S SPMHD CODE (BINARY DUMPS)
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
! ntot(maxstep)       : total number of particles in each timestep
! iam(maxpart,maxstep): integer identification of particle type
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step 
!
! most of these values are stored in global arrays 
! in the module 'particle_data'
!-------------------------------------------------------------------------

subroutine read_data(rootname,indexstart,nstepsread)
  use exact, only:hfact
  use particle_data
  use params
  use labels
  use settings_data
  implicit none
  integer, intent(IN) :: indexstart
  integer, intent(OUT) :: nstepsread
  character(LEN=*), intent(IN) :: rootname
  character(LEN=LEN(rootname)+4) :: datfile
  integer :: i,j,k,ifile,icol,ipos,ierr
  integer :: ncol_max,ndim_max,npart_max,ndimV_max,nstep_max
  integer :: npartin,ntotin,ncolstep,nparti,ntoti
  logical :: reallocate
  real(doub_prec) :: timein,gammain,hfactin
  real(doub_prec), dimension(:), allocatable :: dattemp
  real(doub_prec), dimension(:,:), allocatable :: dattempvec

  ndim_max = 1
  ndimV_max = 1
  nstepsread = 0
  if (rootname(1:1).ne.' ') then
     !
     !--if rootname does not contain .dat, make it end in .dat
     !
     if (index(rootname,'.dat').eq.0) then
        datfile = trim(rootname)//'.dat'
     else
        datfile = trim(rootname)  
     endif
     ifile = 1
     !print*,'rootname = ',rootname
  else
     print*,' **** no data read **** '
     return
  endif

  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  k=1      

  !
  !--open data file and read data
  !
  open(unit=11,iostat=ierr,file=datfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print*,'*** Error opening '//trim(datfile)//' ***'
     return
  endif
!
!--read first header line
!
  read(11,iostat=ierr,end=80) timein,npartin,ntotin,gammain, &
       hfactin,ndim_max,ndimV_max,ncol_max
  if (ierr /= 0) then
     print "(a,/,a,/,a,/,a)",'*** Error reading first header ***', &
           '*** This is possibly because the dump is in single precision', &
           '*** whilst this subroutine assumes double precision. ', &
           '*** Edit read_data subroutine and recompile '
     close(11)
     return
  endif
!  print*,'reading time = ',timein,npartin,ntotin,gammain, &
!       ndim_max,ndimV_max,ncol_max     
!
!--allocate memory for data arrays (initially for 11 timesteps)
!
  if (ntotin.lt.5500) then
     nstep_max = max(111,maxstep)
     elseif (ntotin.lt.111111) then
     nstep_max = max(11,maxstep)
  else 
     nstep_max = max(5,maxstep)
  endif
  npart_max = max(ntotin,maxpart)
  if (.not.allocated(dat) .or. ntotin.gt.maxpart  &
       .or. nstep_max.gt.maxstep .or. ncol_max.gt.maxcol) then
     call alloc(npart_max,nstep_max,ncol_max)
  endif
!
!--rewind file
!
  rewind(11)

  i = indexstart - 1
  nstepsread = 0
  
  do while (i <= maxstep)
     i = i + 1
     !!print*,' reading step ',i
     reallocate = .false.
     npart_max = maxpart
     nstep_max = maxstep
     !
     !--read header line for this timestep
     !
     read(11,iostat=ierr,end=67) timein,nparti,ntoti,gammain, &
          hfactin,ndim,ndimV,ncolstep
     if (ierr /= 0) then
        print*,'*** error reading timestep header ***'
        close(11)     
        return
     else ! count this as a successfully read timestep, even if data is partial
        nstepsread = nstepsread + 1
     endif
          
     time(i) = real(timein)
     gamma(i) = real(gammain)
     hfact = real(hfactin)
     npartoftype(1,i) = nparti
     npartoftype(2,i) = ntoti - nparti
     ntot(i) = ntoti
     print*,'reading time = ',time(i),nparti,ntoti,gamma(i), &
          hfact,ndim,ndimV,ncolstep
     !
     !--check for errors in timestep header
     !
     if (ndim.gt.3 .or. ndimV.gt.3) then
        print*,'*** error in header: ndim or ndimV in file> 3'
        nstepsread = nstepsread - 1
        ndim = ndim_max
        ndimV = ndimV_max
        close(11)
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

     if (ntot(i).gt.maxpart) then
        !print*, 'ntot greater than array limits!!'    
        reallocate = .true.
        npart_max = int(1.1*ntot(i))
     endif
     if (i.eq.maxstep) then
        nstep_max = i + max(10,INT(0.1*nstep_max))
        reallocate = .true.
     endif
     !
     !--reallocate memory for main data array
     !
     if (reallocate) then
        call alloc(npart_max,nstep_max,ncol_max)
     endif

  
     if (ntot(i).gt.0) then
        !
        !--read position vector
        !
        icol = 1
        if (allocated(dattempvec)) deallocate(dattempvec)
        allocate(dattempvec(3,ntot(i)))
        read (11,end=66,ERR=67) dattempvec(1:ndim,1:ntot(i))
        do ipos = 1,ndim
           dat(1:ntot(i),ipos,i) = real(dattempvec(ipos,1:ntot(i)))
        enddo
        icol = icol + ndim
        !
        !--read velocity vector
        !
        read (11,end=66,ERR=67) dattempvec(1:ndimV,1:ntot(i))        
        do ipos = icol,icol+ndimV-1
           dat(1:ntot(i),ipos,i) = real(dattempvec(ipos-icol+1,1:ntot(i)))
        enddo
        icol = icol + ndimV
        !
        !--read scalar variables
        !
        if (allocated(dattemp)) deallocate(dattemp)
        allocate(dattemp(ntot(i)))
        do j = 1,5
           read (11,end=66,ERR=67) dattemp(1:ntot(i))
           dat(1:ntot(i),icol,i) = real(dattemp(1:ntot(i)))
           icol = icol + 1
        enddo
        
        !
        !--non-MHD output
        !
        if (ncolumns.le.ndim+7+ndimV) then
        !
        !--read alpha, alphau
        !
           read (11,end=66,ERR=67) dattempvec(1:2,1:ntot(i))
           do ipos = icol,icol+1
              dat(1:ntot(i),ipos,i) = real(dattempvec(ipos-icol+1,1:ntot(i)))
           enddo
           icol = icol + 2
        
        else
        !
        !--MHD output
        !
        
        !
        !--read alpha, alphau, alphaB
        !
           read (11,end=66,ERR=67) dattempvec(1:3,1:ntot(i))
           do ipos = icol,icol+2
              dat(1:ntot(i),ipos,i) = real(dattempvec(ipos-icol+1,1:ntot(i)))
           enddo
           icol = icol + 3
        !
        !--Bfield
        !
           read (11,end=66,ERR=67) dattempvec(1:ndimV,1:ntot(i))
           do ipos = icol, icol+ndimV-1
              dat(1:ntot(i),ipos,i) = real(dattempvec(ipos-icol+1,1:ntot(i)))
           enddo
           icol = icol + ndimV
        !
        !--curl B
        !
           read (11,end=66,ERR=67) dattempvec(1:ndimV,1:ntot(i))
           do ipos = icol, icol+ndimV-1
              dat(1:ntot(i),ipos,i) = real(dattempvec(ipos-icol+1,1:ntot(i)))
           enddo
           icol = icol + ndimV
        !
        !--div B
        !
           read (11,end=66,ERR=67) dattemp(1:ntot(i))
           dat(1:ntot(i),icol,i) = real(dattemp(1:ntot(i)))
           icol = icol + 1
        !
        !--psi
        !
           read (11,end=66,ERR=67) dattemp(1:ntot(i))
           dat(1:ntot(i),icol,i) = real(dattemp(1:ntot(i)))
           icol = icol + 1        
        endif
        !!print*,'columns read = ',icol,' should be = ',ncolumns

        deallocate(dattemp,dattempvec)
     else
        ntot(i) = 1
        npartoftype(1,i) = 1
        npartoftype(2,i) = 0
        dat(:,:,i) = 0.
     endif
     iam(:,i) = 0
  enddo

  print*,' REACHED ARRAY LIMITS IN READFILE'

  ! this is if reached array limits
  ntot(i-1) = j-1
  goto 68

66 continue
   print "(a)",'*** WARNING: incomplete data on last timestep'
nstepsread = i - indexstart + 1                ! timestep there but data incomplete
ntot(i) = j-1
goto 68

67 continue
   print "(a)",' > end of file <'
68 continue
  !
  !--close data file and return
  !                    
close(unit=11)

ncolumns = ncol_max
ndim = ndim_max
ndimV = ndimV_max
print*,'ncolumns = ',ncolumns

print*,'>> READ steps ',indexstart,'->',indexstart + nstepsread - 1, &
       ' last step ntot = ',ntot(indexstart+nstepsread-1)
return    
!
!--errors
!
80 continue
print*,' *** data file empty, no steps read ***'
return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
 use labels
 use params
 use settings_data
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
 irho = ndim + ndimV + 1      ! location of rho in data array
 ipr = ndim + ndimV + 2       !  pressure 
 iutherm = ndim + ndimV + 3   !  thermal energy
 ih = ndim + ndimV + 4        !  smoothing length
 ipmass = ndim + ndimV + 5    !  particle mass      

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
 label(ipr) = 'P      '
 label(iutherm) = 'u'
 label(ih) = 'h       '
 label(ipmass) = 'particle mass'
 label(ndim + ndimV+6) = '\ga'
 label(ndim + ndimV+7) = '\ga\du'
 if (ncolumns.gt.ndim+ndimV+7) then
    !
    !--mag field (vector)
    !
    label(ndim + ndimV+8) = '\ga\dB'
    iBfirst = ndim + ndimV+8+1        ! location of Bx
    iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
    labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
    do i=1,ndimV
       label(ndim + ndimV+8+i) = trim(labelvec(iBfirst))//'\d'//labelcoord(i,1) !' (x10\u-3\d)' !//'/rho'
    enddo
    !
    !--current density (vector)
    !
    iJfirst = ndim+2*ndimV+8+1
    iamvec(iJfirst:iJfirst+ndimV-1) = iJfirst
    labelvec(iJfirst:iJfirst+ndimV-1) = 'J'
    do i=1,ndimV
       label(ndim + ndimV+ndimV+8 + i) = labelvec(iJfirst)//labelcoord(i,1)
    enddo
    idivB = ndim+3*ndimV+9 
    label(idivB) = 'div B'
 else
    iBfirst = 0
 endif
 if (ncolumns.gt.ndim+3*ndimV+9) then
    label(ndim+3*ndimV+10) = 'psi'
 endif
 if (ncolumns.gt.ndim+3*ndimV+10) then
    label(ndim+3*ndimV+11) = 'f_visc_x'
    label(ndim+3*ndimV+12) = 'f_visc_y'
    label(ndim+3*ndimV+13) = 'f_x'
    label(ndim+3*ndimV+14) = 'f_y'
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
 
!-----------------------------------------------------------

 return 
end subroutine set_labels
