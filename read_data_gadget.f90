!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR OUTPUT FROM THE GADGET CODE
!
! NOTE THAT THIS ONLY "OFFICIALLY" WORKS WITH THE PARALLEL CODE AS WE
! REQUIRE KNOWLEDGE OF THE PARTICLE SMOOTHING LENGTHS
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! dat(maxpart,maxplot,maxstep) : main data array
!
! npartoftype(maxstep): number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step 
!                      (used in calc_quantities for calculating the pressure)
!
! most of these values are stored in global arrays 
! in the module 'particle_data'
!
! Partial data read implemented Nov 2006 means that columns with
! the 'required' flag set to false are not read (read is therefore much faster)
!-------------------------------------------------------------------------

subroutine read_data(rootname,istepstart,nstepsread)
  use particle_data, only:dat,npartoftype,time,gamma,maxpart,maxcol,maxstep
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc,iformat,required,ipartialread
  use mem_allocation, only:alloc
  use labels, only:ih,irho
  use system_utils, only:renvironment
  implicit none
  integer, intent(in) :: istepstart
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  character(len=len(rootname)+10) :: datfile
  integer, dimension(maxparttypes) :: npartoftypei,Nall
  integer, dimension(:), allocatable :: iamtemp
  integer :: i,j,itype,icol,ierr
  integer :: index1,index2,indexstart,indexend,Nmassesdumped
  integer :: ncolstep,npart_max,nstep_max,ntoti
  integer :: iFlagSfr,iFlagFeedback,iFlagCool,nfiles
  logical :: iexist,reallocate
  real(doub_prec) :: timetemp,ztemp
  real(doub_prec), dimension(6) :: Massoftype
  real, dimension(:), allocatable :: dattemp1
  real :: hsoft

  nstepsread = 0
  if (maxparttypes.lt.6) then
     print*,' *** ERROR: not enough particle types for GADGET data read ***'
     print*,' *** you need to edit splash parameters and recompile ***'
     stop
  endif
  
  if (len_trim(rootname).gt.0) then
     datfile = trim(rootname)
  else
     print*,' **** no data read **** ' 
     return
  endif
!
!--check if first data file exists
!
  inquire(file=datfile,exist=iexist)
  if (.not.iexist) then
     print "(a)",' *** error: '//trim(datfile)//': file not found ***'    
     return
  endif
!
!--set parameters which do not vary between timesteps
!
  ndim = 3
  ndimV = 3
!
!--read data from snapshots
!  
  i = istepstart

  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  !--open data file and read data
  !
  open(11,iostat=ierr,file=datfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print "(a)", '*** ERROR OPENING FILE ***'
     return
  endif
  !
  !--read header for this timestep
  !
  read(11,iostat=ierr) npartoftypei(1:6),Massoftype,timetemp,ztemp, &
      iFlagSfr,iFlagFeedback,Nall(1:6),iFlagCool,nfiles
  if (ierr /= 0) then
     print "(a)", '*** ERROR READING TIMESTEP HEADER ***'
     return
  endif

  iformat = 0
  if (iFlagCool.gt.0) then
     iformat = 1
     ncolstep = 12 ! 3 x pos, 3 x vel, pmass, utherm, rho, Ne, Nh, h
     ncolumns = ncolstep
  else
     iformat = 0
     ncolstep = 10 ! 3 x pos, 3 x vel, pmass, utherm, rho, h
     ncolumns = ncolstep  
  endif
  irho = 9
  ih = ncolstep
  
  ntoti = int(sum(npartoftypei(1:6)))
  print*,'time             : ',timetemp
  print*,'Npart (by type)  : ',npartoftypei
  print*,'Mass  (by type)  : ',Massoftype
  print*,'N_gas            : ',npartoftypei(1)
  print*,'N_total          : ',ntoti
  print*,'N data columns   : ',ncolstep

  if (nfiles.gt.1) then
     print*,' nfiles = ',nfiles
     print*,'*** ERROR: read from > 1 files not implemented'
     return
  endif
  !
  !--if successfully read header, increment the nstepsread counter
  !
  nstepsread = nstepsread + 1
  !
  !--now read data
  !
  reallocate = .false.
  npart_max = maxpart
  nstep_max = max(maxstep,1)

  if (ntoti.gt.maxpart) then
     reallocate = .true.
     if (maxpart.gt.0) then
        ! if we are reallocating, try not to do it again
        npart_max = int(1.1*ntoti)
     else
        ! if first time, save on memory
        npart_max = int(ntoti)
     endif
  endif
  if (i.ge.maxstep .and. i.ne.1) then
     nstep_max = i + max(10,INT(0.1*nstep_max))
     reallocate = .true.
  endif
  !
  !--reallocate memory for main data array
  !
  if (reallocate .or. .not.(allocated(dat))) then
     call alloc(npart_max,nstep_max,max(ncolstep+ncalc,maxcol))
  endif
  !
  !--copy header into header arrays
  !
  npartoftype(:,i) = npartoftypei
!  time(i) = real(timetemp)
!--use this line for redshift
  time(i) = real(ztemp)
  
  !
  !--read particle data
  !
  if (ntoti.gt.0) then
     !
     !--read positions of all particles
     !
     if (any(required(1:3))) then
        print*,'positions ',ntoti
        read (11, iostat=ierr) (dat(j,1:3,i),j=1,ntoti)
        if (ierr /= 0) then
           print "(a)",'error encountered whilst reading positions '
           return
        endif
     else
        read(11, iostat=ierr)
        if (ierr /= 0) then
           print "(a)",'error skipping positions '
           return
        endif
     endif
     !
     !--same for velocities
     !
     if (any(required(4:6))) then
        print*,'velocities ',ntoti
        read (11, iostat=ierr) (dat(j,4:6,i),j=1,ntoti)
        if (ierr /= 0) then
           print "(a)",'error encountered whilst reading velocities'
        endif
     else
        read(11, iostat=ierr)
        if (ierr /= 0) then
           print "(a)",'error skipping velocities '
           return
        endif
     endif
     !
     !--skip read of particle ID (only required if we sort the particles
     !  back into their correct order, which is not implemented at present)
     !
     !print*,'particle ID ',ntoti
     !if (allocated(iamtemp)) deallocate(iamtemp)
     !allocate(iamtemp(npart_max))
     read (11,iostat=ierr) ! iamtemp(1:ntoti)
     !deallocate(iamtemp)
     if (ierr /= 0) then
        print "(a)",'error encountered whilst reading particle ID'
     endif
     !
     !--read particle masses
     !
     !--work out total number of masses dumped 
     Nmassesdumped = 0
     do itype = 1,6
        if (abs(Massoftype(itype)).lt.tiny(Massoftype)) then
           Nmassesdumped = Nmassesdumped + Npartoftype(itype,i)
        endif
     enddo
     
     if (required(7)) then
        print*,'particle masses ',Nmassesdumped
        !--read this number of entries
        if (allocated(dattemp1)) deallocate(dattemp1)
        allocate(dattemp1(Nmassesdumped))
        if (Nmassesdumped.gt.0) then
           read(11,iostat=ierr) dattemp1(1:Nmassesdumped)
        endif
        if (ierr /= 0) then
           print "(a)",'error reading particle masses'
        endif
        !--now copy to the appropriate sections of the .dat array
        indexstart = 1
        index1 = 1

        do itype=1,6
           if (Npartoftype(itype,i).ne.0) then
              index2 = index1 + Npartoftype(itype,i) -1
              if (abs(Massoftype(itype)).lt.1.e-8) then ! masses dumped
                 indexend = indexstart + Npartoftype(itype,i) - 1
                 print*,'read ',Npartoftype(itype,i),' masses for type ', &
                        itype,index1,'->',index2,indexstart,'->',indexend
                 dat(index1:index2,7,i) = dattemp1(indexstart:indexend)
                 indexstart = indexend + 1
              else  ! masses not dumped
                 print*,'setting masses for type ',itype,' = ', &
                        real(Massoftype(itype)),index1,'->',index2
                 dat(index1:index2,7,i) = real(Massoftype(itype))
              endif
              index1 = index2 + 1
           endif
        enddo
        deallocate(dattemp1)
     elseif (Nmassesdumped.gt.0) then
        read(11,iostat=ierr)
        if (ierr /= 0) then
           print "(a)",'error reading particle masses'
        endif
     endif
     !
     !--read other quantities for rest of particles
     !
     print*,'gas properties ',npartoftype(1,i)
     do icol=8,ncolstep
        !!print*,icol
        if (npartoftype(1,i).gt.0) then
           if (required(icol)) then
              read (11,iostat=ierr) dat(1:npartoftype(1,i),icol,i)
           else
              read (11,iostat=ierr)
           endif
           if (ierr /= 0) then
              print "(a,i3)",'error reading particle data from column ',icol
           endif
        !
        !--for some reason the smoothing length output by GADGET is
        !  twice the usual SPH smoothing length
        !
           if (icol.eq.ncolstep .and. required(icol)) then
              dat(1:npartoftype(1,i),icol,i) = 0.5*dat(1:npartoftype(1,i),icol,i)
           endif
        endif
     enddo
     !
     !--if a value for the dark matter smoothing length is set
     !  via the environment variable GSPLASH_DARKMATTER_HSOFT,
     !  give dark matter particles this smoothing length
     !  and a density of 1 (so column density plots work)
     !
     hsoft = renvironment('GSPLASH_DARKMATTER_HSOFT')
     if (hsoft.gt.tiny(hsoft)) then
        if (required(ih)) then
           print "(a,1pe10.3,a)",' Assigning smoothing length of h = ',hsoft, &
                                 ' to dark matter particles'
           dat(npartoftype(1,i)+1:npartoftype(1,i)+npartoftype(2,i),ih,i) = hsoft
        endif
        if (required(irho)) then
           dat(npartoftype(1,i)+1:npartoftype(1,i)+npartoftype(2,i),irho,i) = 1.0
        endif
     endif

  else
     ntoti = 1
     npartoftype(1,i) = 1
     dat(:,:,i) = 0.
  endif
!
!--now memory has been allocated, set arrays which are constant for all time
!
  gamma = 5./3.
!
!--set flag to indicate that only part of this file has been read 
!
  if (.not.all(required(1:ncolstep))) ipartialread = .true.
!
!--close data file and return
!                    
  close(unit=11)

  if (nstepsread.gt.0) then
     print*,'>> last step ntot =',sum(npartoftype(:,istepstart+nstepsread-1))
  endif
  return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels, only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass,ih,irho,ipr,iutherm
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ntypes,UseTypeInRenderings
  use geometry, only:labelcoord
  use system_utils, only:renvironment
  implicit none
  integer :: i
  real :: hsoft

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
  ivx = 4
  ipmass = 7
  irho = 9        ! location of rho in data array
  ipr = 0
  iutherm = 8     !  thermal energy
  !
  !--set labels of the quantities read in
  !
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  label(irho) = 'density'
  label(iutherm) = 'u'
  label(ipmass) = 'particle mass'
  
  if (ncolumns.gt.10) then
     label(10) = 'Ne'
     label(11) = 'Nh'
     ih = 12        !  smoothing length
  else
     ih = 10
  endif
  label(ih) = 'h'
  !
  !--set labels for vector quantities
  !
  iamvec(ivx:ivx+ndimV-1) = ivx
  labelvec(ivx:ivx+ndimV-1) = 'v'
  do i=1,ndimV
     label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
  enddo
  
  !--set labels for each particle type
  !
  ntypes = 5
  labeltype(1) = 'gas'
  labeltype(2) = 'dark matter'
  labeltype(5) = 'star'
  UseTypeInRenderings(1) = .true.
  !
  !--dark matter particles are of non-SPH type (ie. cannot be used in renderings)
  !  unless they have had a smoothing length defined
  !
  hsoft = renvironment('GSPLASH_DARKMATTER_HSOFT')
  if (hsoft.gt.tiny(hsoft)) then
     UseTypeInRenderings(2) = .true.
  else
     UseTypeInRenderings(2) = .false.  
  endif
  UseTypeInRenderings(3:5) = .false.

!-----------------------------------------------------------
  return
end subroutine set_labels
