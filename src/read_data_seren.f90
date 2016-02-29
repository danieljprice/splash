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
!  Copyright (C) 2005-2013 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR OUTPUT FROM THE SEREN CODE
! HANDLES BOTH ASCII AND BINARY FILES
!
! THE FOLLOWING ENVIRONMENT VARIABLES AFFECT THIS FORMAT:
!
! DSPLASH_EXTRACOLS : set to number of extra columns to read (after itype)
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
! Partial data read implemented means that columns with
! the 'required' flag set to false are not read (read is therefore much faster)
!-------------------------------------------------------------------------

module seren_data_store
  implicit none
  integer           :: seren_maxparttypes      ! Number of types we are using

  character(len=20) :: format_id               ! File format (for verification)
  integer           :: nunits                  ! Number of units
  integer           :: ndata                   ! Number of data entries
  integer           :: ptot, stot              ! Number of particles/sinks
  integer           :: pboundary, picm, pgas   ! Number of each type of particle
  integer           :: pcdm, pdust, pion       ! Number of each type of particle
  integer           :: PR, NDIMtemp, VDIMtemp, BDIMtemp ! Important parameters
  integer           :: dmdt_range              ! DMDT_RANGE
  character(len=20) :: data_id(1:500)          ! Char ids of arrays written
  character(len=20) :: unit_data(1:500)        ! Unit data
  real              :: unit_coeff(1:500)       ! Unit multiplier
  integer :: typedata(1:5,1:500)               ! type data header array
  integer           :: itemp, iporig           ! Since SPLASH does not have one
  integer           :: iunknown(1:1000)        ! For unknown data types
  character(len=20) :: r_unit                  ! length unit
  character(len=20) :: m_unit                  ! mass unit
  character(len=20) :: rho_unit                ! density unit
  character(len=20) :: h_unit                  ! smoothing length unit
  integer, parameter :: DP = selected_real_kind(p=15) ! double precision
  integer, parameter :: SP = selected_real_kind(p=6)  ! single precision
  integer, parameter :: ILP = selected_int_kind(r=15)  ! Integer long precision

  ! Length units in metres
  real(kind=DP),parameter :: r_pc    = 3.08568E16_DP     ! parsec
  real(kind=DP),parameter :: r_au    = 1.49597870E11_DP  ! astronomical unit
  real(kind=DP),parameter :: r_sun   = 6.96E8_DP         ! solar radius
  real(kind=DP),parameter :: r_earth = 6.371E6_DP        ! Earth radius

  ! Mass units in kilograms
  real(kind=DP),parameter :: m_sun   = 1.98892E30_DP     ! solar mass
  real(kind=DP),parameter :: m_jup   = 1.8986E27_DP      ! Jupiter mass
  real(kind=DP),parameter :: m_earth = 5.9736E24_DP      ! Earth mass

  ! Time units in seconds
  real(kind=DP),parameter :: myr = 3.1556952E13_DP       ! megayear
  real(kind=DP),parameter :: yr  = 3.1556952E7_DP        ! year
  real(kind=DP),parameter :: day = 8.64E4_DP             ! day

end module seren_data_store

subroutine read_data(rootname,istepstart,ipos,nstepsread)
  use particle_data, only:dat,iamtype,npartoftype,time,gamma,maxpart,maxcol,maxstep
  use params
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc,ipartialread,ntypes
  use settings_units, only:unitzintegration, unit_interp
  use mem_allocation, only:alloc
  use labels,         only:labeltype,labelzintegration
  use system_utils,   only:ienvironment
  use seren_data_store
  implicit none
  integer, intent(in) :: istepstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  character(len=len(rootname)+10) :: datfile
  integer, parameter :: iunit = 16
  integer :: i,step,ierr,iambinaryfile,itype
  integer :: npart_max,nstep_max
  logical :: iexist,reallocate,doubleprec
  character(len=50) :: string

  integer            :: idata(1:50)
  integer (kind=ILP) :: ilpdata(1:50)
  real               :: rdata(1:50)
  real(doub_prec)    :: rdata_dp(1:50)
  real(doub_prec)    :: dpdata(1:50)
  real :: timetemp,gammatemp

  unit_coeff = 1. ! not yet used
  m_unit = ""
  rho_unit = ""
  h_unit = ""
  !iRescale = .TRUE.
  ipartialread = .false. ! we always read full data file

  seren_maxparttypes = min(maxparttypes,7)

  nstepsread = 0

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
!--read data from snapshots
!
  step = istepstart

  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  !--open data file and read data
  !
  !
  !--determine whether file is binary or ascii, open it and read the header
  !
  ! Try binary first, then ascii
  open(unit=iunit,file=datfile,status='old',form='unformatted',iostat=ierr)

  if (ierr /= 0) then
     print "(a)",'*** ERROR OPENING '//trim(datfile)//' ***'
     return
  endif
  !
  !--read the file header
  !  try binary format first, and if unsuccessful try ascii
  !
  read (unit=iunit,iostat=ierr) format_id

  if (ierr /= 0 .OR. trim(adjustl(format_id)) /= "SERENBINARYDUMPV2") then
     ! Ascii format
     iambinaryfile = 0
     close (unit=iunit)
     open(unit=iunit,file=datfile,status='old',form='formatted',iostat=ierr)
     if (ierr /= 0) then
        print "(a)",'*** ERROR OPENING '//trim(datfile)//' AS ASCII ***'
        return
     endif
     rewind(unit=iunit)
     read (unit=iunit,fmt=*,iostat=ierr) format_id
     if (ierr /= 0) then
        print "(a)",'*** ERROR OPENING '//trim(datfile)//' - UNKNOWN FILE FORMAT '//format_id//' ***'
        return
     endif
     if (trim(adjustl(format_id)) /= "SERENASCIIDUMPV2") then
        print "(a)",'*** ERROR OPENING '//trim(datfile)//' AS ASCII - WRONG FILE FORMAT ***'
        return
     end if
  else
     iambinaryfile = 1
  end if

  if (iambinaryfile==1) then
     print "(a)",' reading binary seren v2 format '
     read (iunit) PR
     read (iunit) NDIMtemp
     read (iunit) VDIMtemp
     read (iunit) BDIMtemp
  else
     print "(a)",' reading ascii seren v2 format '
     read (iunit,*) PR
     read (iunit,*) NDIMtemp
     read (iunit,*) VDIMtemp
     read (iunit,*) BDIMtemp
  end if

  if (iambinaryfile==0) then
     ! Don't care about precision
     doubleprec = .FALSE.
  else if (PR == 8 .OR. PR == 2) then
     ! Double precision file
     print "(a)",' Double precision file'
     doubleprec = .TRUE.
  else if (PR == 4 .OR. PR == 1) then
     ! Single precision file
     print "(a)",' Single precision file'
     doubleprec = .FALSE.
  else
     print "(a)",'*** WARNING OPENING '//trim(datfile)//' - ASSUMING SINGLE PRECISION ***'
     doubleprec = .FALSE.
  end if

  typedata = 0

  if (iambinaryfile.eq.1) then
     call read_serenheader_binary(iunit)
  else if (iambinaryfile.eq.0) then
     call read_serenheader_ascii(iunit)
  end if
  !
  !--get values of quantities from the header
  !
  ndim = NDIMtemp
  ndimv = VDIMtemp

  ptot = idata(1)
  stot = idata(2)
  pboundary = idata(3)
  picm = idata(4)
  pgas = idata(5)
  pcdm = idata(6)
  pdust = idata(7)
  pion = idata(8)
  dmdt_range = idata(30)

  !--check for errors in integer header (either from corrupt file or wrong endian)
  if (ptot+stot.le.0 .or. ptot+stot.gt.1.e10) then
     if (iambinaryfile.eq.1) then
        print "(a)",' ERROR reading binary file header: wrong endian? '
     else
        print "(a)",' ERROR reading ascii file header '
     endif
     close(unit=iunit)
     return
  endif

  timetemp = real(dpdata(1))
  gammatemp = 1.  ! Not saved in file

  print*,'time             : ',timetemp
  print*,'gamma            : ',gammatemp
  print*,'n_total          : ',ptot

  call set_labels
    !
  !--if successfully read header, increment the nstepsread counter
  !
  nstepsread = nstepsread + 1
  !
  !-- now work out dimensionless weight unit and z integration unit
  !
  call find_weights(unit_interp,unitzintegration,labelzintegration)
  !
  !--now read data
  !
  reallocate = .false.
  npart_max = maxpart
  nstep_max = max(maxstep,1)

  if ((ptot+stot).gt.maxpart) then
     reallocate = .true.
     if (maxpart.gt.0) then
        ! if we are reallocating, try not to do it again
        npart_max = int(1.1*(ptot+stot))
     else
        ! if first time, save on memory
        npart_max = int(ptot+stot)
     endif
  endif
  if (step.ge.maxstep .and. step.ne.1) then
     nstep_max = step + max(10,INT(0.1*nstep_max))
     reallocate = .true.
  endif
  !
  !--reallocate memory for main data array
  !
  if (reallocate .or. .not.(allocated(dat))) then
     call alloc(npart_max,nstep_max,max(ncolumns+ncalc,maxcol),mixedtypes=.true.)
  endif
  !
  !--copy header into header arrays
  !
  npartoftype(:,step) = 0
!   npartoftype(1,step) = ptot
  time(step) = timetemp
  gamma(step) = gammatemp
  !
  !--read particle data
  !
  if (ptot.gt.0) then
!      if (iambinaryfile.eq.1) then
!         call read_dragonbody_binary(iunit,ierr)
!      else
!         call read_dragonbody_ascii(iunit,ierr)
!      endif
     call read_serenbody(iunit,ierr)
  else
     ptot = 0
!      npartoftype(1,step) = 0
!      npartoftype(:,step) = 0
     dat(:,:,step) = 0.
  endif

!   if (allocated(iamtype)) then
!      !--relabel particle types
  call set_types(iamtype(:,step),ptot+stot,npartoftype(:,step))
!   endif
  if (any(npartoftype(2:,step).ne.0)) then
     do itype=1,ntypes
        if (npartoftype(itype,step).gt.0) then
           string = ' '
           write(string,"(a)") 'n_'//trim(labeltype(itype))
           write(string(18:len(string)),"(a)") ':'
           print*,trim(string),' ',npartoftype(itype,step)
        endif
     enddo
  endif
! !
! !--set flag to indicate that only part of this file has been read
! !
!   if (.not.all(required(1:ncolumns))) ipartialread = .true.
!
!--close data file and return
!
  close(unit=iunit)

  return

contains

!----------------------------------------------------
! binary header read
!----------------------------------------------------
subroutine read_serenheader_binary(iunitb)
   implicit none
   integer, intent(in) :: iunitb

   read (iunitb,end=55) idata
   read (iunitb,end=55) ilpdata
   if (doubleprec) then
      read (iunitb,end=55) rdata_dp
   else
      read (iunitb,end=55) rdata
   endif
   read (iunitb,end=55) dpdata

   ndata = idata(21)
   nunits = idata(20)

   if (nunits > 0) read (iunitb) unit_data(1:nunits)
   if (ndata > 0)  read (iunitb) data_id(1:ndata)
   if (ndata > 0)  read (iunitb) typedata(1:5,1:ndata)

   return

55 continue
   print "(a)",' ERROR: end of file in binary header read'
   stop
   !return

end subroutine read_serenheader_binary

!----------------------------------------------------
! ascii header read
!----------------------------------------------------
subroutine read_serenheader_ascii(iunita)
   implicit none
   integer, intent(in) :: iunita

   do i=1,size(idata)
      read (iunita,*,end=55) idata(i)
   end do

   do i=1,size(ilpdata)
      read (iunita,*,end=55) ilpdata(i)
   end do

   do i=1,size(rdata)
      read (iunita,*,end=55) rdata(i)
   end do

   do i=1,size(dpdata)
      read (iunita,*,end=55) dpdata(i)
   end do

   ndata = idata(21)
   nunits = idata(20)

   do i=1,nunits
      read (iunita,'(A)') unit_data(i)
   end do
   do i=1,ndata
      read (iunita,'(A)') data_id(i)
   end do
   do i=1,ndata
      read (iunita,*) typedata(1:5,i)
   end do

   return

55 continue
   print "(a)",' ERROR: end of file in ascii header read'
   stop
   !return

end subroutine read_serenheader_ascii

!----------------------------------------------------
! body read
!----------------------------------------------------

  subroutine read_serenbody(iunit,ierr_out)
     use seren_data_store
     use labels, only:ix,ivx,ipmass,ih,irho,iBfirst,iutherm
     implicit none
     integer, intent(in)  :: iunit
     integer, intent(out) :: ierr_out
     integer              :: ierr, ierr1
     integer              :: j, k
     integer, dimension(:), allocatable         :: dummy_int
     integer, dimension(:,:), allocatable       :: dummy_int_2D
     real(kind=SP), dimension(:,:), allocatable :: dummy
     real(kind=DP), dimension(:,:), allocatable :: dummy_dp
     real(kind=SP), dimension(:), allocatable :: dummy_scalar
     real(kind=DP), dimension(:), allocatable :: dummy_dp_scalar
     integer(kind=ILP), dimension(:), allocatable :: dummy_ilp_int
     logical, dimension(:), allocatable           :: dummy_logical
     integer(kind=ILP), dimension(:,:), allocatable :: dummy_ilp_int_2D
     logical, dimension(:,:), allocatable           :: dummy_logical_2D
     logical       :: ldummy(1:2)
     integer       :: idummy2(1:2)
     real(kind=SP), allocatable :: raux(:)
     real(kind=DP), allocatable :: raux_dp(:)
     real, allocatable :: sink_dat_array(:,:,:)
     integer       :: s
     integer       :: unknown
     integer :: pfirst
     integer :: plast
     integer :: width
     integer :: typecode
     integer :: unit_id
     integer :: sink_data_length
     character(len=30)  :: sink_format_string
     character(len=50)  :: format_string
     integer :: nl,ni,nli,npr,ndp,nchar
     logical, allocatable           :: l_data_st(:)
     integer, allocatable           :: i_data_st(:)
     integer(kind=ILP), allocatable :: ilp_data_st(:)
     real(kind=SP), allocatable     :: sp_data_st(:)
     real(kind=DP), allocatable     :: dp_data_st(:)
     real(kind=DP), allocatable     :: dp_data_st2(:)

     unknown = 0

     if (stot>0) then
        sink_data_length = 11+NDIMtemp+VDIMtemp+2*dmdt_range
        allocate(raux(1:sink_data_length))
        if (doubleprec) allocate(raux_dp(1:sink_data_length))
        allocate(sink_dat_array(1:stot,1:ncolumns,1))
        sink_dat_array = 0.
        write (sink_format_string,'(A,I0,A)') "(", sink_data_length, "E18.10)"
     end if

     ierr_out = 0

     allocate(dummy_int(1:ptot),stat=ierr)
     if (ierr /= 0) then
        print *,' ERROR allocating memory'
        goto 56
     endif
     if (doubleprec) allocate(dummy_dp_scalar(1:ptot),stat=ierr)
     if (ierr /= 0) then
        print *,' ERROR allocating memory'
        goto 56
     endif
     allocate(dummy_scalar(1:ptot),stat=ierr)
     if (ierr /= 0) then
        print *,' ERROR allocating memory'
        goto 56
     endif

     write (6,'(A,A)',ADVANCE="NO") "Loading: "
     do i=1,ndata
        if (i==1) then
           write (6,'(A)',ADVANCE="NO") trim(data_id(i))
        else
           write (6,'(A,A)',ADVANCE="NO") ", ", trim(data_id(i))
        end if
        select case (trim(data_id(i)))
           case ("porig")
              ! Original particle number
              ! Read through porig numbers
              pfirst = typedata(2,i); plast = typedata(3,i)
              if (iambinaryfile==1) then
                 read(iunit,end=55,iostat=ierr) dummy_int(pfirst:plast)
              else
                 do k=pfirst,plast
                    read(iunit,fmt=*,end=55,iostat=ierr1) dummy_int(k)
                    if (ierr1 /= 0) ierr = ierr1
                 end do
              end if
              if (ierr /= 0) then
                 print*,' WARNING: errors reading through porig '
                 ierr_out = -1
              end if
              do k=pfirst,plast
                 dat(k,iporig,step) = real(dummy_int(k))
              enddo
              if (ierr /= 0) then
                 print*,' WARNING: errors reading unknown data type ', trim(data_id(i))
                 ierr_out = -1
              end if
           case ("r")
              ! Position
              if (doubleprec) allocate(dummy_dp(1:NDIMtemp,1:ptot),stat=ierr)
              allocate(dummy(1:NDIMtemp,1:ptot),stat=ierr)
              pfirst = typedata(2,i); plast = typedata(3,i)
              if (iambinaryfile==1) then
                 if (doubleprec) then
                    read(iunit,end=55,iostat=ierr) dummy_dp(1:NDIMtemp,pfirst:plast)
                    dummy(1:NDIMtemp,pfirst:plast) = real(dummy_dp(1:NDIMtemp,pfirst:plast))
                 else
                    read(iunit,end=55,iostat=ierr) dummy(1:NDIMtemp,pfirst:plast)
                 end if
              else
                 do k=pfirst,plast
                    read(iunit,*,end=55,iostat=ierr1) dummy(1:NDIMtemp,k)
                    if (ierr1 /= 0) ierr = ierr1
                 end do
              end if
              do k=pfirst,plast
                 dat(k,ix(1):ix(1)+NDIMtemp-1,step) = dummy(1:NDIMtemp,k)
              enddo
              if (ierr /= 0) then
                 print*,' WARNING: errors reading positions '
                 ierr_out = -1
              end if
              deallocate(dummy)
              if (doubleprec) deallocate(dummy_dp)
           case ("h")
              ! Smoothing lengths
              pfirst = typedata(2,i); plast = typedata(3,i)
              if (iambinaryfile==1) then
                 if (doubleprec) then
                    read(iunit,end=55,iostat=ierr) dummy_dp_scalar(pfirst:plast)
                    dummy_scalar(pfirst:plast) = real(dummy_dp_scalar(pfirst:plast))
                 else
                    read(iunit,end=55,iostat=ierr) dummy_scalar(pfirst:plast)
                 end if
              else
                 do k=pfirst,plast
                    read(iunit,*,end=55,iostat=ierr1) dummy_scalar(k)
                    if (ierr1 /= 0) ierr = ierr1
                 end do
              end if
              do k=pfirst,plast
                 dat(k,ih,step) = dummy_scalar(k)
              enddo
              if (ierr /= 0) then
                 print*,' WARNING: errors reading smoothing lengths'
                 ierr_out = -1
              end if
           case ("m")
              ! Mass
              pfirst = typedata(2,i); plast = typedata(3,i)
              if (iambinaryfile==1) then
                 if (doubleprec) then
                    read(iunit,end=55,iostat=ierr) dummy_dp_scalar(pfirst:plast)
                    dummy_scalar(pfirst:plast) = real(dummy_dp_scalar(pfirst:plast))
                 else
                    read(iunit,end=55,iostat=ierr) dummy_scalar(pfirst:plast)
                 end if
              else
                 do k=pfirst,plast
                    read(iunit,*,end=55,iostat=ierr1) dummy_scalar(k)
                    if (ierr1 /= 0) ierr = ierr1
                 end do
              end if
              do k=pfirst,plast
                 dat(k,ipmass,step) = dummy_scalar(k)
              enddo
              if (ierr /= 0) then
                 print*,' WARNING: errors reading masses'
                 ierr_out = -1
              end if
           case ("v")
              ! Velocities
              if (doubleprec) allocate(dummy_dp(1:VDIMtemp,1:ptot),stat=ierr)
              allocate(dummy(1:VDIMtemp,1:ptot),stat=ierr)
              pfirst = typedata(2,i); plast = typedata(3,i)
              if (iambinaryfile==1) then
                 if (doubleprec) then
                    read(iunit,end=55,iostat=ierr) dummy_dp(1:VDIMtemp,pfirst:plast)
                    dummy(1:VDIMtemp,pfirst:plast) = real(dummy_dp(1:VDIMtemp,pfirst:plast))
                 else
                    read(iunit,end=55,iostat=ierr) dummy(1:VDIMtemp,pfirst:plast)
                 end if
              else
                 do k=pfirst,plast
                    read(iunit,*,end=55,iostat=ierr1) dummy(1:VDIMtemp,k)
                    if (ierr1 /= 0) ierr = ierr1
                 end do
              end if
              do k=pfirst,plast
                 dat(k,ivx:ivx+VDIMtemp-1,step) = dummy(1:VDIMtemp,k)
              enddo
              if (ierr /= 0) then
                 print*,' WARNING: errors reading velocities '
                 ierr_out = -1
              end if
              deallocate(dummy)
              if (doubleprec) deallocate(dummy_dp)
           case ("rho")
              ! Densities
              pfirst = typedata(2,i); plast = typedata(3,i)
              if (iambinaryfile==1) then
                 if (doubleprec) then
                    read(iunit,end=55,iostat=ierr) dummy_dp_scalar(pfirst:plast)
                    dummy_scalar(pfirst:plast) = real(dummy_dp_scalar(pfirst:plast))
                 else
                    read(iunit,end=55,iostat=ierr) dummy_scalar(pfirst:plast)
                 end if
              else
                 do k=pfirst,plast
                    read(iunit,*,end=55,iostat=ierr1) dummy_scalar(k)
                    if (ierr1 /= 0) ierr = ierr1
                 end do
              end if
              do k=pfirst,plast
                 dat(k,irho,step) = dummy_scalar(k)
              enddo
              if (ierr /= 0) then
                 print*,' WARNING: errors reading densities'
                 ierr_out = -1
              end if
           case ("temp")
              ! Temperatures
              pfirst = typedata(2,i); plast = typedata(3,i)
              if (iambinaryfile==1) then
                 if (doubleprec) then
                    read(iunit,end=55,iostat=ierr) dummy_dp_scalar(pfirst:plast)
                    dummy_scalar(pfirst:plast) = real(dummy_dp_scalar(pfirst:plast))
                 else
                    read(iunit,end=55,iostat=ierr) dummy_scalar(pfirst:plast)
                 end if
              else
                 do k=pfirst,plast
                    read(iunit,*,end=55,iostat=ierr1) dummy_scalar(k)
                    if (ierr1 /= 0) ierr = ierr1
                 end do
              end if
              do k=pfirst,plast
                 dat(k,itemp,step) = dummy_scalar(k)
              enddo
              if (ierr /= 0) then
                 print*,' WARNING: errors reading temperatures'
                 ierr_out = -1
              end if
           case ("u")
              ! Internal energy
              pfirst = typedata(2,i); plast = typedata(3,i)
              if (iambinaryfile==1) then
                 if (doubleprec) then
                    read(iunit,end=55,iostat=ierr) dummy_dp_scalar(pfirst:plast)
                    dummy_scalar(pfirst:plast) = real(dummy_dp_scalar(pfirst:plast))
                 else
                    read(iunit,end=55,iostat=ierr) dummy_scalar(pfirst:plast)
                 end if
              else
                 do k=pfirst,plast
                    read(iunit,*,end=55,iostat=ierr1) dummy_scalar(k)
                    if (ierr1 /= 0) ierr = ierr1
                 end do
              end if
              do k=pfirst,plast
                 dat(k,iutherm,step) = dummy_scalar(k)
              enddo
              if (ierr /= 0) then
                 print*,' WARNING: errors reading internal energy'
                 ierr_out = -1
              end if
           case ("B")
              ! Magnetic fields
              if (doubleprec) allocate(dummy_dp(1:BDIMtemp,1:ptot),stat=ierr)
              allocate(dummy(1:BDIMtemp,1:ptot),stat=ierr)
              pfirst = typedata(2,i); plast = typedata(3,i)
              if (iambinaryfile==1) then
                 if (doubleprec) then
                    read(iunit,end=55,iostat=ierr) dummy_dp(1:BDIMtemp,pfirst:plast)
                    dummy(1:BDIMtemp,pfirst:plast) = real(dummy_dp(1:BDIMtemp,pfirst:plast))
                 else
                    read(iunit,end=55,iostat=ierr) dummy(1:BDIMtemp,pfirst:plast)
                 end if
              else
                 do k=pfirst,plast
                    read(iunit,*,end=55,iostat=ierr1) dummy(1:BDIMtemp,k)
                    if (ierr1 /= 0) ierr = ierr1
                 end do
              end if
              do k=pfirst,plast
                 dat(k,iBfirst:iBfirst+BDIMtemp-1,step) = dummy(1:BDIMtemp,k)
              enddo
              if (ierr /= 0) then
                 print*,' WARNING: errors reading magnetic fields'
                 ierr_out = -1
              end if
              deallocate(dummy)
              if (doubleprec) deallocate(dummy_dp)
           case ("sink_v1")
              ! Load sinks in sink data storage, will add them later
              pfirst = typedata(2,i); plast = typedata(3,i)
              if (iambinaryfile==1) then
                 read(iunit,end=55,iostat=ierr) nl,ni,nli,npr,ndp,nchar
              else
                 read(iunit,fmt=*,end=55,iostat=ierr) nl,ni,nli,npr,ndp,nchar
              end if
              do s=pfirst,plast
                 if (iambinaryfile==1) then
                    read(iunit,end=55,iostat=ierr) ldummy
                    read(iunit,end=55,iostat=ierr) idummy2
                    if (doubleprec) then
                       read(iunit,end=55,iostat=ierr) raux_dp
                       raux = real(raux_dp)
                    else
                       read(iunit,end=55,iostat=ierr) raux
                    end if
                 else
                    read(iunit,'(2L1)',end=55,iostat=ierr) ldummy
                    read(iunit,fmt=*,end=55,iostat=ierr) idummy2
                    read(iunit,sink_format_string) raux(1:sink_data_length)
                 end if
                 if (ix(1)/=0) sink_dat_array(s,ix(1):ix(NDIMtemp),1) = raux(2:NDIMtemp+1)
                 if (ivx/=0) sink_dat_array(s,ivx:ivx+VDIMtemp-1,1) = raux(NDIMtemp+2:NDIMtemp+VDIMtemp+1)
                 if (ipmass/=0) sink_dat_array(s,ipmass,1) = raux(NDIMtemp+VDIMtemp+2)
                 if (ih/=0) sink_dat_array(s,ih,1) = raux(NDIMtemp+VDIMtemp+3)
                 if (itemp/=0) sink_dat_array(s,itemp,1) = raux(NDIMtemp+VDIMtemp+11)
              end do
           case default
              !print*,' WARNING: unknown data type ', trim(data_id(i))
              ! Assume this is an unknown data type
              !ierr_out = -4
              width = typedata(1,i)
              pfirst = typedata(2,i); plast = typedata(3,i)
              typecode = typedata(4,i); unit_id = typedata(5,i)
              unknown = unknown + 1
              if (typecode == 7) then
                 ! Special data structure we don't understand; read and skip
                 if (iambinaryfile==1) then
                    read(iunit,end=55,iostat=ierr) nl,ni,nli,npr,ndp,nchar
                 else
                    read(iunit,fmt=*,end=55,iostat=ierr) nl,ni,nli,npr,ndp,nchar
                 end if
                 if (nchar > 0) stop "Fail! character data :("
                 if (nl > 0) allocate(l_data_st(1:nl))
                 if (ni > 0) allocate(i_data_st(1:ni))
                 if (nli > 0) allocate(ilp_data_st(1:nli))
                 if (npr > 0) then
                    allocate(sp_data_st(1:npr))
                    if (doubleprec) allocate(dp_data_st(1:npr))
                 end if
                 if (ndp > 0) allocate(dp_data_st2(1:ndp))
                 if (iambinaryfile==1) then
                    do j=pfirst, plast
                       if (nl > 0) read(iunit,end=55,iostat=ierr) l_data_st
                       if (ni > 0) read(iunit,end=55,iostat=ierr) i_data_st
                       if (nli > 0) read(iunit,end=55,iostat=ierr) ilp_data_st
                       if (npr > 0) then
                          if (doubleprec) then
                             read(iunit,end=55,iostat=ierr) dp_data_st
                          else
                             read(iunit,end=55,iostat=ierr) sp_data_st
                          end if
                       end if
                       if (ndp > 0) read(iunit,end=55,iostat=ierr) dp_data_st2
                    end do
                 else
                    do j=pfirst, plast
                       if (nl > 0) then
                          write (format_string,'(A,I0,A)') "(",nl,"L1)"
                          read(iunit,end=55,iostat=ierr,fmt=format_string) l_data_st
                       end if
                       if (ni > 0) then
                          read(iunit,end=55,iostat=ierr,fmt=*) i_data_st
                       end if
                       if (nli > 0) then
                          read(iunit,end=55,iostat=ierr,fmt=*) ilp_data_st
                       end if
                       if (npr > 0) then
                          if (doubleprec) then
                             read(iunit,end=55,iostat=ierr,fmt=*) dp_data_st
                          else
                             read(iunit,end=55,iostat=ierr,fmt=*) sp_data_st
                          end if
                       end if
                       if (ndp > 0) then
                          read(iunit,end=55,iostat=ierr,fmt=*) dp_data_st2
                       end if
                    end do
                 end if
                 if (ierr /= 0) then
                    print*,' WARNING: errors reading unknown data structure ', trim(data_id(i))
                    ierr_out = -1
                 end if
                 if (allocated(i_data_st)) deallocate(i_data_st)
                 if (allocated(ilp_data_st)) deallocate(ilp_data_st)
                 if (allocated(sp_data_st)) deallocate(sp_data_st)
                 if (allocated(dp_data_st)) deallocate(dp_data_st)
                 if (allocated(dp_data_st2)) deallocate(dp_data_st2)
              else if (typecode == 6) then
                 ! Character data; read and skip
                 stop "Fail! character data :("
                 ! I have realised this was a silly idea
                 ! If we work out how this should work, I can put it in
              else if (typecode >= 1 .AND. typecode <= 5) then
                 ! 1 = logical, 2 = integer, 3 = long integer, 4 = PR, 5 = DP
                 ! Normal data set; either scalar or vector
                 if (width == 1) then
                    ! Scalar data
                    if (typecode==1) then
                       ! Logical data array
                       allocate(dummy_logical(1:ptot))
                       dummy_logical = .FALSE.
                       if (iambinaryfile==1) then
                          read(iunit,end=55,iostat=ierr) dummy_logical(pfirst:plast)
                       else
                          do k=pfirst,plast
                             read(iunit,'(L1)',end=55,iostat=ierr1) dummy_logical(k)
                             if (ierr1 /= 0) ierr = ierr1
                          end do
                       end if
                       where (dummy_logical)
                          dummy_scalar=1.d0
                       elsewhere
                          dummy_scalar=0.d0
                       end where
                       deallocate(dummy_logical)
                    else if (typecode==2) then
                       ! Integer data array
                       dummy_int = 0
                       if (iambinaryfile==1) then
                          read(iunit,end=55,iostat=ierr) dummy_int(pfirst:plast)
                       else
                          do k=pfirst,plast
                             read(iunit,*,end=55,iostat=ierr1) dummy_int(k)
                             if (ierr1 /= 0) ierr = ierr1
                          end do
                       end if
                       dummy_scalar(pfirst:plast) = real(dummy_int(pfirst:plast))
                    else if (typecode==3) then
                       ! Long integer data array
                       allocate(dummy_ilp_int(1:ptot))
                       dummy_ilp_int = 0
                       if (iambinaryfile==1) then
                          read(iunit,end=55,iostat=ierr) dummy_ilp_int(pfirst:plast)
                       else
                          do k=pfirst,plast
                             read(iunit,*,end=55,iostat=ierr1) dummy_ilp_int(k)
                             if (ierr1 /= 0) ierr = ierr1
                          end do
                       end if
                       dummy_scalar(pfirst:plast) = real(dummy_ilp_int(pfirst:plast))
                       deallocate(dummy_ilp_int)
                    else if (typecode==4) then
                       ! PR data array
                       if (iambinaryfile==1) then
                          if (doubleprec) then
                             read(iunit,end=55,iostat=ierr) dummy_dp_scalar(pfirst:plast)
                             dummy_scalar(pfirst:plast) = real(dummy_dp_scalar(pfirst:plast))
                          else
                             read(iunit,end=55,iostat=ierr) dummy_scalar(pfirst:plast)
                          end if
                       else
                          do k=pfirst,plast
                             read(iunit,*,end=55,iostat=ierr1) dummy_scalar(k)
                             if (ierr1 /= 0) ierr = ierr1
                          end do
                       end if
                    else if (typecode==5) then
                       ! DP data array
                       if (iambinaryfile==1) then
                          if (.NOT.doubleprec) allocate(dummy_dp_scalar(pfirst:plast))
                          read(iunit,end=55,iostat=ierr) dummy_dp_scalar(pfirst:plast)
                          dummy_scalar(pfirst:plast) = real(dummy_dp_scalar(pfirst:plast))
                          if (.NOT.doubleprec) deallocate(dummy_dp_scalar)
                       else
                          do k=pfirst,plast
                             read(iunit,*,end=55,iostat=ierr1) dummy_scalar(k)
                             if (ierr1 /= 0) ierr = ierr1
                          end do
                       end if
                    end if
                    do k=pfirst,plast
                       dat(k,iunknown(unknown),step) = dummy_scalar(k)
                    enddo
                    if (ierr /= 0) then
                       print*,' WARNING: errors reading unknown data type ', trim(data_id(i))
                       ierr_out = -1
                    end if
                 else
                    ! Vector data
                    allocate(dummy(1:width,1:ptot),stat=ierr)
                    if (typecode == 1) then
                       ! Logical data array
                       allocate(dummy_logical_2D(1:width,1:ptot))
                       dummy_logical_2D = .FALSE.
                       if (iambinaryfile==1) then
                          read(iunit,end=55,iostat=ierr) dummy_logical_2D(1:width,pfirst:plast)
                       else
                          write (format_string,'(A,I0,A)') "(",width,"L1)"
                          do k=pfirst,plast
                             read(iunit,fmt=format_string,end=55,iostat=ierr1) dummy_logical_2D(1:width,k)
                             if (ierr1 /= 0) ierr = ierr1
                          end do
                       end if
                       where (dummy_logical_2D)
                          dummy=1.d0
                       elsewhere
                          dummy=0.d0
                       end where
                       deallocate(dummy_logical_2D)
                    else if (typecode == 2) then
                       ! Integer data array
                       allocate(dummy_int_2D(1:width,1:ptot))
                       dummy_int_2D = 0
                       if (iambinaryfile==1) then
                          read(iunit,end=55,iostat=ierr) dummy_int_2D(1:width,pfirst:plast)
                       else
                          do k=pfirst,plast
                             read(iunit,*,end=55,iostat=ierr1) dummy_int_2D(1:width,k)
                             if (ierr1 /= 0) ierr = ierr1
                          end do
                       end if
                       dummy(1:width,pfirst:plast) = real(dummy_int_2D(1:width,pfirst:plast))
                       deallocate(dummy_int_2D)
                    else if (typecode == 3) then
                       ! Long integer data array
                       allocate(dummy_ilp_int_2D(1:width,1:ptot))
                       dummy_ilp_int_2D = 0
                       if (iambinaryfile==1) then
                          read(iunit,end=55,iostat=ierr) dummy_ilp_int_2D(1:width,pfirst:plast)
                       else
                          do k=pfirst,plast
                             read(iunit,*,end=55,iostat=ierr1) dummy_ilp_int_2D(1:width,k)
                             if (ierr1 /= 0) ierr = ierr1
                          end do
                       end if
                       dummy(1:width,pfirst:plast) = real(dummy_ilp_int_2D(1:width,pfirst:plast))
                       deallocate(dummy_ilp_int_2D)
                    else if (typecode == 4) then
                       ! PR data array
                       if (doubleprec) allocate(dummy_dp(1:width,1:ptot),stat=ierr)
                       if (iambinaryfile==1) then
                          if (doubleprec) then
                             read(iunit,end=55,iostat=ierr) dummy_dp(1:width,pfirst:plast)
                             dummy(1:width,pfirst:plast) = real(dummy_dp(1:width,pfirst:plast))
                          else
                             read(iunit,end=55,iostat=ierr) dummy(1:width,pfirst:plast)
                          end if
                       else
                          do k=pfirst,plast
                             read(iunit,*,end=55,iostat=ierr1) dummy(1:width,k)
                             if (ierr1 /= 0) ierr = ierr1
                          end do
                       end if
                       if (doubleprec) deallocate(dummy_dp)
                    else if (typecode == 5) then
                       ! DP data array
                       if (iambinaryfile==1) then
                          allocate(dummy_dp(1:width,1:ptot),stat=ierr)
                          read(iunit,end=55,iostat=ierr) dummy_dp(1:width,pfirst:plast)
                          dummy(1:width,pfirst:plast) = real(dummy_dp(1:width,pfirst:plast))
                          deallocate(dummy_dp)
                       else
                          do k=pfirst,plast
                             read(iunit,*,end=55,iostat=ierr1) dummy(1:width,k)
                             if (ierr1 /= 0) ierr = ierr1
                          end do
                       end if
                    end if
                    do k=pfirst,plast
                       dat(k,iunknown(unknown):iunknown(unknown)+width-1,step) = dummy(1:width,k)
                    enddo
                    if (ierr /= 0) then
                       print*,' WARNING: errors reading unknown data type ', trim(data_id(i))
                       ierr_out = -1
                    end if
                    deallocate(dummy)
                 end if
              end if

        end select
     end do
     write (6,*)

     if (stot>0) then
        ! Load sink stuff into end of dat array
        dat(ptot+1:ptot+stot,1:ncolumns,step) = sink_dat_array(1:stot,1:ncolumns,1)
     end if

     return

55   continue
     if (iambinaryfile==1) print "(a)",' ERROR: end of file in binary read'
     if (iambinaryfile==0) print "(a)",' ERROR: end of file in ascii read'
     ierr_out = -3
     return

56   continue
     ierr_out = -2

     return
  end subroutine

!----------------------------------------------------
! translate types into order (for old dragon read)
!----------------------------------------------------
subroutine set_types(itypei,ntotal,noftype)
 implicit none
 integer(kind=int1), dimension(:), intent(inout) :: itypei
 integer, intent(in) :: ntotal
 integer, dimension(:), intent(out) :: noftype
 integer             :: noftype_temp(1:7)

 noftype = 0

 noftype_temp(1) = pgas
 noftype_temp(2) = pboundary
 noftype_temp(3) = stot
 noftype_temp(4) = picm
 noftype_temp(5) = pcdm
 noftype_temp(6) = pdust
 noftype_temp(7) = pion

 if (sum(noftype_temp(1:7)).ne.ntotal) then
    print "(a)",' INTERNAL ERROR setting number in each type in dragon read'
 endif

 do i=1,7
    if (i > seren_maxparttypes .AND. noftype_temp(i) > 0) then
       print*,' *** ERROR: not enough particle types for SEREN data read ***'
       print*,' *** you need to edit splash parameters and recompile ***'
       stop
    end if
 end do

 noftype(1:seren_maxparttypes) = noftype_temp(1:seren_maxparttypes)

 if (pboundary>0) itypei(1:pboundary) = 2
 if (picm>0) itypei(pboundary+1:pboundary+picm) = 4
 if (pgas>0) itypei(pboundary+picm+1:pboundary+picm+pgas) = 1
 if (pcdm>0) itypei(pboundary+picm+pgas+1:pboundary+picm+pgas+pcdm) = 5
 if (pdust>0) itypei(pboundary+picm+pgas+pcdm+1:pboundary+picm+pgas+pcdm+pdust) = 6
 if (pion>0) itypei(pboundary+picm+pgas+pcdm+pdust+1:ptot) = 7
 if (stot>0) itypei(ptot+1:ptot+stot) = 3

 return
end subroutine set_types

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels, only:label,iamvec,labelvec,labeltype,unitslabel,&
     &ix,ivx,ipmass,ih,irho,iBfirst,iutherm,lenlabel,lenunitslabel
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ntypes,UseTypeInRenderings
  use geometry, only:labelcoord
  use settings_units, only:units
  use seren_data_store
  implicit none
  integer :: i, j, width, unit_no
  integer :: nunknown                ! Number of unknown data types
  character(len=lenunitslabel) :: unit_base, unit_string
  character(len=lenlabel)      :: type_names(1:7)
  logical                      :: type_use_render(1:7)

  if (ndim.le.0 .or. ndim.gt.3) then
     print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
     return
  endif
  if (ndimV.le.0 .or. ndimV.gt.3) then
     print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
     return
  endif

  ! Calculate number of columns to read
  iporig = 0
  ncolumns = 0
  nunknown = 0
  do i=1,ndata
     unit_no = typedata(5,i)
     !write (6,*) "i = ", i, "; data_id(i) = ", data_id(i)
     unit_base = ""
     unit_string = ""
     if (unit_no < 0 .OR. unit_no > nunits) then
        print*,'*** ERROR: unit_no = ',unit_no,' in set_labels ***'
     else if (unit_no /= 0) then
        unit_base = trim(adjustl(unit_data(unit_no)))
        unit_string = unit_base
        call translate_unit_names(unit_string)
        unit_string = ' ['//trim(adjustl(unit_string))//']'
     end if
     !write (6,*) "unit_base, unit_string = ", unit_base, unit_string
     select case (trim(data_id(i)))
        case ("porig")
           ! Original particle number
           iporig = -1 ! We always want this last
        case ("r")
           ! Position
           do j=1,NDIMtemp
              ix(j) = ncolumns + j
              units(ix(j)) = unit_coeff(ix(j))
              unitslabel(ix(j)) = unit_string
           end do
!            unitzintegration = units(ix(1))
!            labelzintegration = ' ['//trim(adjustl(unit_string))//']'
           label(ix(1:ndim)) = labelcoord(1:ndim,1)
           ncolumns = ncolumns + NDIMtemp
           r_unit = trim(adjustl(unit_base))
        case ("h")
           ! Smoothing lengths
           ih = ncolumns + 1
           units(ih) = unit_coeff(ih)
           unitslabel(ih) = unit_string
           label(ih) = 'h'
           ncolumns = ncolumns + 1
           h_unit = trim(adjustl(unit_base))
        case ("m")
           ! Mass
           ipmass = ncolumns + 1
           units(ipmass) = unit_coeff(ipmass)
           unitslabel(ipmass) = unit_string
           label(ipmass) = 'particle mass'
           ncolumns = ncolumns + 1
           m_unit = trim(adjustl(unit_base))
        case ("v")
           ! Velocities
           ivx = ncolumns + 1
           iamvec(ivx:ivx+VDIMtemp-1) = ivx
           labelvec(ivx:ivx+VDIMtemp-1) = 'v'
           do j=1,VDIMtemp
              label(ivx+j-1) = 'v\d'//labelcoord(j,1)
              units(ivx+j-1) = unit_coeff(ivx+j-1)
              unitslabel(ivx+j-1) = unit_string
           enddo
           ncolumns = ncolumns + VDIMtemp
        case ("rho")
           ! Densities
           irho = ncolumns + 1
           units(irho) = unit_coeff(irho)
           unitslabel(irho) = unit_string
           label(irho) = 'density'
           ncolumns = ncolumns + 1
           rho_unit = trim(adjustl(unit_base))
        case ("temp")
           ! Temperatures
           itemp = ncolumns + 1 ! NOT A PROPER SPLASH i_quantity
           units(itemp) = unit_coeff(itemp)
           unitslabel(itemp) = unit_string
           label(ncolumns + 1) = 'temperature'
           ncolumns = ncolumns + 1
        case ("u")
           ! Internal energy
           iutherm = ncolumns + 1
           units(iutherm) = unit_coeff(iutherm)
           unitslabel(iutherm) = unit_string
           label(ncolumns + 1) = 'internal energy'
           ncolumns = ncolumns + 1
        case ("B")
           ! Magnetic fields
           iBfirst = ncolumns + 1
           iamvec(iBfirst:iBfirst+BDIMtemp-1) = iBfirst
           labelvec(iBfirst:iBfirst+BDIMtemp-1) = 'B'
           do j=1,BDIMtemp
              label(iBfirst+j-1) = 'B\d'//labelcoord(j,1)
              units(iBfirst+j-1) = unit_coeff(iBfirst+j-1)
              unitslabel(iBfirst+j-1) = unit_string
           enddo
           ncolumns = ncolumns + BDIMtemp
        case ("sink_v0")
           ! Do nothing yet, sinks are stored separately
        case ("sink_v1")
           ! Do nothing yet, sinks are stored separately
        case default
           print*,' WARNING reading file: unknown data type ', trim(data_id(i))
           if (typedata(4,i) == 7) cycle ! Special data module we don't understand; ignore
           if (typedata(4,i) == 6) cycle ! Not a lot we can do with character data here!
           width = typedata(1,i)
           nunknown = nunknown + 1
           iunknown(nunknown) = ncolumns + 1 ! NOT A PROPER SPLASH i_quantity
           if (width == 1) then
              label(iunknown(nunknown)) = trim(data_id(i))
              units(iunknown(nunknown)) = unit_coeff(iunknown(nunknown))
              unitslabel(iunknown(nunknown)) = unit_string
           else if (width <= NDIMtemp) then
              do j=1,width
                 label(iunknown(nunknown)+j-1) = trim(data_id(i))//'\d'//labelcoord(j,1)
                 units(iunknown(nunknown)+j-1) = unit_coeff(iunknown(nunknown))
                 unitslabel(iunknown(nunknown)+j-1) = unit_string
              end do
           else
              do j=1,width
                 write(label(iunknown(nunknown)+j-1),'(A,A,I0)') trim(data_id(i)),'\d',j
                 units(iunknown(nunknown)+j-1) = unit_coeff(iunknown(nunknown))
                 unitslabel(iunknown(nunknown)+j-1) = unit_string
              end do
           end if
           ncolumns = ncolumns + width ! Add width of data to ncolumns
     end select
  end do
  if (iporig == -1) then
     ! If there is porig, add as last column
     iporig = ncolumns + 1 ! NOT A PROPER SPLASH i_quantity
     label(ncolumns + 1) = 'particle id'
     ncolumns = ncolumns + 1
  end if

  !--set labels for each particle type
  !
  ntypes = seren_maxparttypes
  type_names = (/'gas     ','boundary','sink    ','icm     ','cdm     ','dust    ','ion     '/)
  type_use_render = (/.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.TRUE./)

  labeltype(1:ntypes) = type_names(1:ntypes)
  UseTypeInRenderings(1:ntypes) = type_use_render(1:ntypes)

!-----------------------------------------------------------
  return
end subroutine set_labels

subroutine find_weights(out_unit_interp,out_unitzintegration,out_labelzintegration)
   use labels, only:lenunitslabel
   use params
   use seren_data_store
   implicit none

   real(doub_prec), intent(out)      :: out_unit_interp
   real, intent(out)                 :: out_unitzintegration
   character(len=lenunitslabel), intent(out) :: out_labelzintegration
   real(doub_prec)                   :: dm_unit, dh_unit, drho_unit, dr_unit
   logical                           :: do_dimweight, do_zintegration
   character(len=lenunitslabel)      :: rho_length_label
   real(doub_prec)                   :: rho_length

   out_unit_interp = 1.0
   out_unitzintegration = 1.0
   out_labelzintegration = ""

   do_dimweight = .TRUE.
   do_zintegration = .TRUE.

   if (m_unit=="") then
      print*,'No masses or no mass units!'
      print*,'Cannot create dimensionless weight (unnormalised rendered plots may be incorrect)'
      do_dimweight = .FALSE.
   end if
   if (h_unit=="") then
      print*,'No smoothing lengths or no smoothing length units!'
      print*,'Cannot create dimensionless weight (unnormalised rendered plots may be incorrect)'
      do_dimweight = .FALSE.
   end if
   if (rho_unit=="") then
      print*,'No densities or no density units!'
      print*,'Cannot create dimensionless weight (unnormalised rendered plots may be incorrect)'
      do_dimweight = .FALSE.
      print*,'Cannot set unitzintegration (column density plots may be incorrect)'
      do_zintegration = .FALSE.
   end if
   if (r_unit=="") then
      print*,'No positions or no position units!'
      print*,'Cannot set unitzintegration (column density plots may be incorrect)'
      do_zintegration = .FALSE.
   end if


! Length unit in S.I. units (m)
   if (r_unit=="pc") then
      dr_unit = r_pc
   else if (r_unit=="au") then
      dr_unit = r_au
   else if (r_unit=="r_sun") then
      dr_unit = r_sun
   else if (r_unit=="r_earth") then
      dr_unit = r_earth
   else if (r_unit=="km") then
      dr_unit = 1000.0_DP
   else if (r_unit=="m") then
      dr_unit = 1.0_DP
   else if (r_unit=="cm") then
      dr_unit = 0.01_DP
   else
      print*,'Unknown position unit ', r_unit, '!'
      print*,'Cannot set unitzintegration (column density plots may be incorrect)'
      do_zintegration = .FALSE.
      dr_unit = 1.0_DP
   end if

! Length unit in S.I. units (m)
   if (h_unit=="pc") then
      dh_unit = r_pc
   else if (h_unit=="au") then
      dh_unit = r_au
   else if (h_unit=="r_sun") then
      dh_unit = r_sun
   else if (h_unit=="r_earth") then
      dh_unit = r_earth
   else if (h_unit=="km") then
      dh_unit = 1000.0_DP
   else if (h_unit=="m") then
      dh_unit = 1.0_DP
   else if (h_unit=="cm") then
      dh_unit = 0.01_DP
   else
      print*,'Unknown smoothing length unit ', h_unit, '!'
      print*,'Cannot create dimensionless weight (unnormalised rendered plots may be incorrect)'
      do_dimweight = .FALSE.
      dh_unit = 1.0_DP
   end if

! Mass units in S.I. units (kg)
   if (m_unit=="m_sun") then
      dm_unit = m_sun
   else if (m_unit=="m_jup") then
      dm_unit = m_jup
   else if (m_unit=="m_earth") then
      dm_unit = m_earth
   else if (m_unit=="kg") then
      dm_unit = 1._DP
   else if (m_unit=="g") then
      dm_unit = 1.0E-3_DP
   else
      print*,'Unknown mass unit ', m_unit, '!'
      print*,'Cannot create dimensionless weight (unnormalised rendered plots may be incorrect)'
      do_dimweight = .FALSE.
      dm_unit = 1._DP
   end if

 ! Density units in S.I. units (i.e. kg/m^3)
   if (rho_unit=="m_sun_pc3") then
      drho_unit = m_sun / (r_pc**3)
      rho_length = r_pc
      rho_length_label = "pc"
   else if (rho_unit=="m_sun_pc2") then
      drho_unit = m_sun / (r_pc**2)
      rho_length = r_pc
      rho_length_label = "pc"
   else if (rho_unit=="kg_m3") then
      drho_unit = 1.0_DP
      rho_length = 1.0_DP
      rho_length_label = "m"
   else if (rho_unit=="kg_m2") then
      drho_unit = 1.0_DP
      rho_length = 1.0_DP
      rho_length_label = "m"
   else if (rho_unit=="g_cm3") then
      drho_unit = 1.0E3_DP
      rho_length = 0.01_DP
      rho_length_label = "cm"
   else if (rho_unit=="g_cm2") then
      drho_unit = 10.0_DP
      rho_length = 0.01_DP
      rho_length_label = "cm"
   else
      print*,'Unknown density unit ', rho_unit, '!'
      print*,'Cannot create dimensionless weight (unnormalised rendered plots may be incorrect)'
      do_dimweight = .FALSE.
      print*,'Cannot set unitzintegration (column density plots may be incorrect)'
      do_zintegration = .FALSE.
      rho_length = 1.0_DP
   end if

   if (do_dimweight) then
      out_unit_interp = dm_unit/(drho_unit*dh_unit**NDIMtemp)
   end if

   if (do_zintegration) then
      out_unitzintegration = dr_unit / rho_length
      out_labelzintegration = rho_length_label
   end if

   return
end subroutine find_weights

subroutine translate_unit_names(unit_name)
   implicit none
   character(len=*), intent(inout) :: unit_name

   select case (trim(unit_name))
      case ("r_sun")
         unit_name = "r\dSun\u"
      case ("au")
         unit_name = "AU"
      case ("r_earth")
         unit_name = "r\dEarth\u"
      case ("m_sun")
         unit_name = "M\d\(2281)\u"
      case ("m_jup")
         unit_name = "M\dJupiter\u"
      case ("m_earth")
         unit_name = "M\dEarth\u"
      case ("myr")
         unit_name = "Myrs"
      case ("km_s")
         unit_name = "km s\u-1\d"
      case ("au_yr")
         unit_name = "AU / yr"
      case ("m_s")
         unit_name = "m s\u-1\d"
      case ("cm_s")
         unit_name = "cm s\u-1\d"
      case ("km_s2")
         unit_name = "km \u-2\d"
      case ("au_yr2")
         unit_name = "AU yr\u-2\d"
      case ("m_s2")
         unit_name = "m s\u-2\d"
      case ("cm_s2")
         unit_name = "cm s\u-2\d"
      case ("m_sun_pc3")
         unit_name = "M\d\(2281)\u pc\u-3\d"
      case ("kg_m_3")
         unit_name = "kg m\u-3\d"
      case ("g_cm_3")
         unit_name = "g cm\u-3\d"
      case ("m_sun_pc2")
         unit_name = "M\d\(2281)\u pc\u-2\d"
      case ("kg_m_2")
         unit_name = "kg m\u-2\d"
      case ("g_cm_2")
         unit_name = "g cm\u-2\d"
      case ("g_cms2")
         unit_name = "g cm\u-2\d"
      case ("10^40erg")
         unit_name = "\x 10\u40\d ergs"
      case ("m_sunkm_s")
         unit_name = "M\d\(2281)\u km s\u-1\d"
      case ("m_sunau_yr")
         unit_name = "M\d\(2281)\u AU yr\u-1\d"
      case ("kgm_s")
         unit_name = "kg m s\u-1\d"
      case ("gcm_s")
         unit_name = "g cm s\u-1\d"
      case ("m_sunkm2_s")
         unit_name = "M\d\(2281)\u km\u2\d s\u-1\d"
      case ("m_sunau2_yr")
         unit_name = "M\d\(2281)\u AU\u2\d yr\u-1\d"
      case ("kgm2_s")
         unit_name = "kg m\u2\d s\u-1\d"
      case ("gcm2_s")
         unit_name = "g cm\u2\d s\u-1\d"
      case ("rad_s")
         unit_name = "radians s\u-1\d"
      case ("m_sun_myr")
         unit_name = "M\d\(2281)\u Myr\u-1\d"
      case ("m_sun_yr")
         unit_name = "M\d\(2281)\u yr\u-1\d"
      case ("kg_s")
         unit_name = "kg s\u-1\d"
      case ("g_s")
         unit_name = "g s\u-1\d"
      case ("L_sun")
         unit_name = "Ld\(2281)\u"
      case ("J_s")
         unit_name = "J s\u-1\d"
      case ("ergs_s")
         unit_name = "ergs s\u-1\d"
      case ("m2_kg")
         unit_name = "m\u2\d kg\u-1\d"
      case ("cm2_g")
         unit_name = "cm\u2\d g\u-1\d"
      case ("tesla")
         unit_name = "Tesla"
      case ("gauss")
         unit_name = "Gauss"
      case ("C_s_m2")
         unit_name = "C s\u-1\d m\u-2\d"
      case ("J_kg")
         unit_name = "J kg\u-1\d"
      case ("erg_g")
         unit_name = "ergs g\u-1\d"
      case default
         unit_name = unit_name ! this could be improved
   end select

   return
end subroutine translate_unit_names

