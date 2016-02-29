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
!  Copyright (C) 2005-2015 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM
! THE NEXT GENERATION SPH CODE (sphNG)
!
! (also my Phantom SPH code which uses a similar format)
!
! *** CONVERTS TO SINGLE PRECISION ***
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
! SSPLASH_RESET_CM if 'YES' then centre of mass is reset to origin
! SSPLASH_OMEGA if non-zero subtracts corotating velocities with omega as set
! SSPLASH_OMEGAT if non-zero subtracts corotating positions and velocities with omega as set
! SSPLASH_TIMEUNITS sets default time units, either 's','min','hrs','yrs' or 'tfreefall'
!
! the data is stored in the global array dat
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
! npartoftype(1:6,maxstep) : number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!
! Partial data read implemented Nov 2006 means that columns with
! the 'required' flag set to false are not read (read is therefore much faster)
!-------------------------------------------------------------------------
module sphNGread
 use params
 implicit none
 real(doub_prec) :: udist,umass,utime,umagfd
 real :: tfreefall
 integer :: istartmhd,istartrt,nmhd,idivvcol,icurlvxcol,icurlvycol,icurlvzcol
 integer :: nhydroreal4,istart_extra_real4
 integer :: nhydroarrays,nmhdarrays
 logical :: phantomdump,smalldump,mhddump,rtdump,usingvecp,igotmass,h2chem,rt_in_header
 logical :: usingeulr,cleaning
 logical :: batcode,tagged,debug
 integer, parameter :: maxarrsizes = 10
 integer, parameter :: maxinblock = 128 ! max allowed in each block
 integer, parameter :: lentag = 16
 character(len=lentag) :: tagarr(maxplot)
 integer, parameter :: itypemap_sink_phantom = 3
 integer, parameter :: itypemap_dust_phantom = 2
 integer, parameter :: itypemap_unknown_phantom = 8

 !------------------------------------------
 ! generic interface to utilities for tagged
 ! dump format
 !------------------------------------------
 interface extract
  module procedure extract_int, extract_real4, extract_real8, &
         extract_intarr, extract_real4arr, extract_real8arr
 end interface extract
 
contains

 !-------------------------------------------------------------------
 ! function mapping iphase setting in sphNG to splash particle types
 !-------------------------------------------------------------------
 elemental integer function itypemap_sphNG(iphase)
  integer*1, intent(in) :: iphase

  select case(int(iphase))
  case(0)
    itypemap_sphNG = 1
  case(1:9)
    itypemap_sphNG = 3
  case(10:)
    itypemap_sphNG = 4
  case default
    itypemap_sphNG = 5
  end select  
  
 end function itypemap_sphNG

 !---------------------------------------------------------------------
 ! function mapping iphase setting in Phantom to splash particle types
 !---------------------------------------------------------------------
 elemental integer function itypemap_phantom(iphase)
  integer*1, intent(in) :: iphase
  
  select case(int(iphase))
  case(1:2)
    itypemap_phantom = iphase
  case(3:6) ! put sinks as type 3, everything else shifted by one
    itypemap_phantom = iphase + 1
  case(-3) ! sink particles, either from external_binary or read from dump
    itypemap_phantom = itypemap_sink_phantom
  case default
    itypemap_phantom = itypemap_unknown_phantom
  end select
  
 end function itypemap_phantom

 !------------------------------------------
 ! extraction of single integer variables
 !------------------------------------------
 subroutine extract_int(tag,ival,intarr,tags,ntags,ierr)
  character(len=*),      intent(in)  :: tag
  integer,               intent(out) :: ival
  integer,               intent(in)  :: ntags,intarr(:)
  character(len=lentag), intent(in)  :: tags(:)
  integer,               intent(out) :: ierr
  logical :: matched
  integer :: i

  ierr = 1
  matched = .false.
  ival = 0 ! default if not found
  over_tags: do i=1,min(ntags,size(tags))
     if (trim(tags(i))==trim(adjustl(tag))) then
        if (size(intarr) >= i) then
           ival = intarr(i)
           matched = .true.
        endif
        exit over_tags  ! only match first occurrence
     endif
  enddo over_tags
  if (matched) ierr = 0
  if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
  
 end subroutine extract_int

 !------------------------------------------
 ! extraction of single real*8 variables
 !------------------------------------------
 subroutine extract_real8(tag,rval,r8arr,tags,ntags,ierr)
  character(len=*),      intent(in)  :: tag
  real*8,                intent(out) :: rval
  real*8,                intent(in)  :: r8arr(:)
  character(len=lentag), intent(in)  :: tags(:)
  integer,               intent(in)  :: ntags
  integer,               intent(out) :: ierr
  logical :: matched
  integer :: i

  ierr = 1
  matched = .false.
  rval = 0.d0 ! default if not found
  over_tags: do i=1,min(ntags,size(tags))
     if (trim(tags(i))==trim(adjustl(tag))) then
        if (size(r8arr) >= i) then
           rval = r8arr(i)
           matched = .true.
        endif
        exit over_tags  ! only match first occurrence
     endif
  enddo over_tags
  if (matched) ierr = 0
  if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
  
 end subroutine extract_real8

 !------------------------------------------
 ! extraction of single real*4 variables
 !------------------------------------------
 subroutine extract_real4(tag,rval,r4arr,tags,ntags,ierr)
  character(len=*),      intent(in)  :: tag
  real*4,                intent(out) :: rval
  real*4,                intent(in)  :: r4arr(:)
  character(len=lentag), intent(in)  :: tags(:)
  integer,               intent(in)  :: ntags
  integer,               intent(out) :: ierr
  logical :: matched
  integer :: i

  ierr = 1
  matched = .false.
  rval = 0. ! default if not found
  over_tags: do i=1,min(ntags,size(tags))
     if (trim(tags(i))==trim(adjustl(tag))) then
        if (size(r4arr) >= i) then
           rval = r4arr(i)
           matched = .true.
        endif
        exit over_tags  ! only match first occurrence
     endif
  enddo over_tags
  if (matched) ierr = 0
  if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
  
 end subroutine extract_real4

 !------------------------------------------
 ! extraction of integer arrays
 !------------------------------------------
 subroutine extract_intarr(tag,ival,intarr,tags,ntags,ierr)
  character(len=*),      intent(in)  :: tag
  integer,               intent(out) :: ival(:)
  integer,               intent(in)  :: ntags,intarr(:)
  character(len=lentag), intent(in)  :: tags(:)
  integer,               intent(out) :: ierr
  integer :: i,nmatched

  ierr = 1
  nmatched = 0
  ival(:) = 0 ! default if not found
  over_tags: do i=1,min(ntags,size(tags))
     if (trim(tags(i))==trim(adjustl(tag))) then
        if (size(intarr) >= i .and. size(ival) > nmatched) then
           nmatched = nmatched + 1
           ival(nmatched) = intarr(i)
        endif
     endif
  enddo over_tags
  if (nmatched==size(ival)) ierr = 0
  if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 
 end subroutine extract_intarr

 !------------------------------------------
 ! extraction of real*8 arrays
 !------------------------------------------
 subroutine extract_real8arr(tag,rval,r8arr,tags,ntags,ierr)
  character(len=*),      intent(in)  :: tag
  real*8,                intent(out) :: rval(:)
  real*8,                intent(in)  :: r8arr(:)
  character(len=lentag), intent(in)  :: tags(:)
  integer,               intent(in)  :: ntags
  integer,               intent(out) :: ierr
  integer :: i,nmatched

  ierr = 1
  nmatched = 0
  rval = 0.d0 ! default if not found
  over_tags: do i=1,min(ntags,size(tags))
     if (trim(tags(i))==trim(adjustl(tag))) then
        if (size(r8arr) >= i .and. size(rval) > nmatched) then
           nmatched = nmatched + 1
           rval(nmatched) = r8arr(i)
        endif
     endif
  enddo over_tags
  if (nmatched==size(rval)) ierr = 0
  if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
  
 end subroutine extract_real8arr

 !------------------------------------------
 ! extraction of real*4 arrays
 !------------------------------------------
 subroutine extract_real4arr(tag,rval,r4arr,tags,ntags,ierr)
  character(len=*),      intent(in)  :: tag
  real*4,                intent(out) :: rval(:)
  real*4,                intent(in)  :: r4arr(:)
  character(len=lentag), intent(in)  :: tags(:)
  integer,               intent(in)  :: ntags
  integer,               intent(out) :: ierr
  integer :: i,nmatched

  ierr = 1
  nmatched = 0
  rval = 0. ! default if not found
  over_tags: do i=1,min(ntags,size(tags))
     if (trim(tags(i))==trim(adjustl(tag))) then
        if (size(r4arr) >= i .and. size(rval) > nmatched) then
           nmatched = nmatched + 1
           rval(nmatched) = r4arr(i)
        endif
     endif
  enddo over_tags
  if (nmatched==size(rval)) ierr = 0
  if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
  
 end subroutine extract_real4arr

 !----------------------------------------------------------------------
 ! Extract various options from the fileident string
 !----------------------------------------------------------------------
 subroutine get_options_from_fileident(fileident,smalldump,tagged,phantomdump,&
                       usingvecp,usingeulr,cleaning,h2chem,rt_in_header,batcode)
  character(len=*), intent(in) :: fileident
  logical,          intent(out) :: smalldump,tagged,phantomdump,batcode
  logical,          intent(out) :: usingvecp,usingeulr,cleaning,h2chem,rt_in_header

  smalldump = .false.
  phantomdump = .false.
  usingvecp = .false.
  usingeulr = .false.
  cleaning = .false.
  h2chem = .false.
  rt_in_header = .false.
  batcode = .false.
  tagged = .false.
  if (fileident(1:1).eq.'S') then
     smalldump = .true.
  endif
  if (fileident(2:2).eq.'T') then
     tagged = .true.
  endif
  if (index(fileident,'Phantom').ne.0) then
     phantomdump = .true.
  else
     phantomdump = .false.
  endif
  if (index(fileident,'vecp').ne.0) then
     usingvecp = .true.
  endif
  if (index(fileident,'eulr').ne.0) then
     usingeulr = .true.
  endif
  if (index(fileident,'clean').ne.0) then
     cleaning = .true.
  endif
  if (index(fileident,'H2chem').ne.0) then
     h2chem = .true.
  endif
  if (index(fileident,'RT=on').ne.0) then
     rt_in_header = .true.
  endif
  if (index(fileident,'This is a test').ne.0) then
     batcode = .true.
  endif

 end subroutine get_options_from_fileident

 !----------------------------------------------------------------------
 ! Routine to read the header of sphNG dump files and extract relevant
 ! information
 !----------------------------------------------------------------------
 subroutine read_header(iunit,iverbose,debug,doubleprec,&
                        npart,npartoftypei,n1,ntypes,nblocks,&
                        narrsizes,realarr,tagsreal,nreals,ierr)
  integer, intent(in)  :: iunit,iverbose
  logical, intent(in)  :: debug,doubleprec
  integer, intent(out) :: npart,npartoftypei(:),n1,ntypes,nblocks,narrsizes,nreals,ierr
  real,    intent(out) :: realarr(maxinblock)
  character(len=lentag), intent(out) :: tagsreal(maxinblock)
  character(len=lentag) :: tags(maxinblock)
  integer               :: intarr(maxinblock)
  real(doub_prec)       :: real8arr(maxinblock)
  real(sing_prec)       :: real4arr(maxinblock)
  integer               :: i,ierr1,ierr2,ierrs(4)
  integer               :: nints,ninttypes,nreal4s,nreal8s,n2,nreassign,naccrete,nkill
  real(doub_prec), allocatable :: dattemp(:)
  real(sing_prec), allocatable :: dattempsingle(:)

  ! initialise empty tag array
  tags(:) = ''
  intarr(:) = 0

  nblocks = 1 ! number of MPI blocks
  npartoftypei(:) = 0
  read(iunit,iostat=ierr) nints
  if (ierr /=0) then
     print "(a)",'error reading nints'
     return
  else
     if (tagged) then
        if (nints > maxinblock) then
           print*,'WARNING: number of ints in header exceeds splash array limit, ignoring some'
           nints = maxinblock
        endif
        read(iunit,iostat=ierr1) tags(1:nints)
        read(iunit,iostat=ierr2) intarr(1:nints)
        if (ierr1 /= 0 .or. ierr2 /= 0) then
           print "(a)",'error reading integer header'
           ierr = 1
           return
        endif
        if (debug) print*,'DEBUG: got tags = ',tags(1:nints)
        call extract('nblocks',nblocks,intarr,tags,nints,ierr)
        if (ierr /= 0) return
        call extract('nparttot',npart,intarr,tags,nints,ierr)
        if (ierr /= 0) return
        if (phantomdump) then
           call extract('ntypes',ntypes,intarr,tags,nints,ierr)
           if (ierr /= 0) return
           call extract('npartoftype',npartoftypei(1:ntypes),intarr,tags,nints,ierr)
           if (ierr /= 0) return
        endif
        if (phantomdump .and. nints < 7) ntypes = nints - 1
        if (iverbose.ge.1) print *,'npart = ',npart,' MPI blocks = ',nblocks
        if (phantomdump) then
           n1 = npartoftypei(1)
        else
           call extract('n1',n1,intarr,tags,nints,ierr)
        endif
     else
        if (nints.lt.3) then
           if (.not.phantomdump) print "(a)",'WARNING: npart,n1,n2 NOT IN HEADER??'
           read(iunit,iostat=ierr) npart
           npartoftypei(1) = npart
        elseif (phantomdump) then
           if (nints.lt.7) then
              ntypes = nints - 1
              read(iunit,iostat=ierr) npart,npartoftypei(1:ntypes)
           else
              ntypes = 5
              read(iunit,iostat=ierr) npart,npartoftypei(1:5),nblocks
           endif
           if (debug) then
              print*,'DEBUG: ntypes = ',ntypes,' npartoftype = ',npartoftypei(:)
           endif
           n1 = npartoftypei(1)
           n2 = 0
        elseif (nints.ge.7) then
           read(iunit,iostat=ierr) npart,n1,n2,nreassign,naccrete,nkill,nblocks
        else
           print "(a)",'warning: nblocks not read from file (assuming non-MPI dump)'
           read(iunit,iostat=ierr) npart,n1,n2
        endif
        if (ierr /=0) then
           print "(a)",'error reading npart,n1,n2 and/or number of MPI blocks'
           return
        elseif (nblocks.gt.2000) then
           print *,'npart = ',npart,' MPI blocks = ',nblocks
           nblocks = 1
           print*,' corrupt number of MPI blocks, assuming 1 '
        else
           if (iverbose.ge.1) print *,'npart = ',npart,' MPI blocks = ',nblocks
        endif
     endif
  endif

  if (ntypes > maxparttypes) then
     print "(a,i2)",' WARNING: number of particle types exceeds array limits: ignoring types > ',maxparttypes
     ntypes = maxparttypes
  endif

!--int*1, int*2, int*4, int*8
  ierr1 = 0
  ierr2 = 0
  do i=1,4
     read(iunit,iostat=ierr) ninttypes
     if (ninttypes > 0) then
        if (tagged) read(iunit,iostat=ierr1)
        read(iunit,iostat=ierr2)
     endif
     if (ierr /= 0 .or. ierr1 /= 0 .or. ierr2 /= 0) then
        print "(a)",'error skipping int types'
        return
     endif
  enddo
!--default reals
  read(iunit,iostat=ierr) nreals
  if (ierr /=0) then
     print "(a)",'error reading default reals'
     return
  else
     if (nreals > maxinblock) then
        print*,'WARNING: number of reals in header exceeds splash array limit, ignoring some'
        nreals = maxinblock
     endif
     if (tagged) then
        read(iunit,iostat=ierr) tagsreal(1:nreals)
     else
     !
     !--set the tags manually for older formats
     !  (but only the ones we care about)
     !
        if (phantomdump) then
           tagsreal(1) = 'time'
           tagsreal(3) = 'gamma'
           tagsreal(4) = 'rhozero'
           tagsreal(6) = 'hfact'
           tagsreal(7) = 'tolh'
           tagsreal(15:19) = 'massoftype'
        elseif (batcode) then
           tagsreal(1) = 'time'
           tagsreal(3) = 'gamma'
           tagsreal(4) = 'radL1'
           tagsreal(5) = 'PhiL1'
           tagsreal(15) = 'Er'
        else
           tagsreal(1) = 'gt'
           tagsreal(2) = 'dtmax'
           tagsreal(3) = 'gamma'
           tagsreal(4) = 'rhozero'
           tagsreal(5) = 'RK2'
           if (smalldump) then ! sphNG small dump
              if (nreals.eq.15) then
                 tagsreal(15) = 'pmassinitial'
              else
                 tagsreal(23) = 'pmassinitial'
              endif
           endif
        endif
     endif
     if (doubleprec) then
        read(iunit,iostat=ierr) real8arr(1:nreals)
        realarr(1:nreals) = real(real8arr(1:nreals))
     else
        read(iunit,iostat=ierr) real4arr(1:nreals)
        realarr(1:nreals) = real(real4arr(1:nreals))
     endif
  endif
!--real*4, real*8
  read(iunit,iostat=ierr) nreal4s
  if (nreal4s > 0) then
     if (tagged) read(iunit,iostat=ierr1)
     read(iunit,iostat=ierr2)
  endif
  if (ierr /= 0 .or. ierr1 /= 0 .or. ierr2 /= 0) then
     print "(a)",'error skipping real*4''s in header'
     return
  endif

  read(iunit,iostat=ierr) nreal8s
  if (ierr /= 0 .or. nreal8s < 0) then
     print "(a)",'error reading nreal8s'
     return
  endif
!   print "(a,i3)",' ndoubles = ',nreal8s
  if (iverbose.ge.1) print "(4(a,i3),a)",' header contains ',nints,' ints, ',&
                           nreals,' reals,',nreal4s,' real4s, ',nreal8s,' doubles'
  if (tagged) then
     if (nreal8s > maxinblock) then
        print*,'WARNING: number of real8''s in header exceeds splash array limit, ignoring some'
        nreal8s = maxinblock
     endif
     read(iunit,iostat=ierr) tags(1:nreal8s)
     read(iunit,iostat=ierr) real8arr(1:nreal8s)
     call extract('udist',udist,real8arr,tags,nreal8s,ierrs(1))
     call extract('umass',umass,real8arr,tags,nreal8s,ierrs(2))
     call extract('utime',utime,real8arr,tags,nreal8s,ierrs(3))
     call extract('umagfd',umagfd,real8arr,tags,nreal8s,ierrs(4))
     if (any(ierrs /= 0)) then
        print "(a)",' *** error reading units'
     endif
  else
     if (nreal8s.ge.4) then
        read(iunit,iostat=ierr) udist,umass,utime,umagfd
     elseif (nreal8s.ge.3) then
        read(iunit,iostat=ierr) udist,umass,utime
        umagfd = 1.0
     else
        print "(a)",'*** WARNING: units not found in file'
        udist = 1.0
        umass = 1.0
        utime = 1.0
        umagfd = 1.0
     endif
  endif
  if (ierr /= 0) then
     print "(a)",'*** error reading units'
  endif
!
!--Total number of array blocks in the file
!
  read(iunit,iostat=ierr) narrsizes
  if (ierr /= 0) return
  if (debug) print*,' nblocks(total)=',narrsizes
  narrsizes = narrsizes/nblocks
  if (ierr /= 0) then
     print "(a)",'*** error reading number of array sizes ***'
     close(iunit)
     return
  elseif (narrsizes.gt.maxarrsizes) then
     narrsizes = maxarrsizes
     print "(a,i2)",'WARNING: too many array sizes: reading only ',narrsizes
  endif
  if (narrsizes.ge.4 .and. nreal8s.lt.4) then
     print "(a)",' WARNING: could not read magnetic units from dump file'
  endif
  if (debug) print*,' number of array sizes = ',narrsizes

 end subroutine read_header

 !----------------------------------------------------------------------
 ! Read the header to each array block
 !----------------------------------------------------------------------
 subroutine read_block_header(iunit,iblock,iarr,iverbose,debug,&
                              isize,nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8,&
                              ntotblock,npart,ntotal,nptmasstot,ncolstep,ierr)
  integer,   intent(in)    :: iunit,iblock,iarr,iverbose
  logical,   intent(in)    :: debug
  integer*8, intent(out)   :: isize(:)
  integer,   intent(out)   :: nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8,ierr
  integer,   intent(inout) :: ntotblock,npart,ntotal,nptmasstot,ncolstep
  
  read(iunit,iostat=ierr) isize(iarr),nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8
  if (iarr.eq.1) then
     ntotblock = isize(iarr)
     if (npart.le.0) npart = ntotblock
     ntotal = ntotal + ntotblock
  elseif (iarr.eq.2) then
     nptmasstot = nptmasstot + isize(iarr)
  endif
  if (debug) print*,'DEBUG: array size ',iarr,' size = ',isize(iarr)
  if (isize(iarr).gt.0 .and. iblock.eq.1) then
     if (iverbose.ge.1) print "(1x,a,i1,a,i12,a,5(i2,1x),a,3(i2,1x))", &
        'block ',iarr,' dim = ',isize(iarr),' nint =',nint,nint1,nint2,nint4,nint8,&
        'nreal =',nreal,nreal4,nreal8
  endif
!--we are going to read all real arrays but need to convert them all to default real
  if (iarr.ne.2 .and. isize(iarr).eq.isize(1) .and. iblock.eq.1) then
     ncolstep = ncolstep + nreal + nreal4 + nreal8
  endif

 end subroutine read_block_header

 !----------------------------------------------------------------------
 ! Extract and print relevant variables from the header block
 !----------------------------------------------------------------------
 subroutine extract_variables_from_header(tags,realarr,nreals,iverbose,debug,&
            gotbinary,nblocks,nptmasstot,npartoftypei,ntypes,&
            time,gamma,hfact,npart,ntotal,npartoftype,massoftype,dat,ix,ih,ipmass,ivx)
  character(len=lentag), intent(in) :: tags(maxinblock)
  real, intent(in) :: realarr(maxinblock)
  integer, intent(in) :: nreals,iverbose,nblocks,nptmasstot,npartoftypei(:),ntypes
  integer, intent(in) :: ix(3),ih,ipmass,ivx
  real, intent(out) :: time,gamma,hfact,massoftype(:)
  real, intent(inout) :: dat(:,:)
  integer, intent(out) :: npartoftype(:)
  integer, intent(inout) :: npart,ntotal
  logical, intent(in)  :: debug
  logical, intent(out) :: gotbinary
  real :: rhozero,tfreefall,tff,radL1,PhiL1,Er,RK2,dtmax,tolh
  real :: massoftypei(ntypes)
  integer :: i,ierrs(10),ipos
  integer :: itype
  integer, parameter :: ilocbinary = 24
  real,    parameter :: pi=3.141592653589
  
  if (phantomdump) then
     call extract('time',time,realarr,tags,nreals,ierrs(1))
  else
     call extract('gt',time,realarr,tags,nreals,ierrs(1))  
  endif
  call extract('gamma',gamma,realarr,tags,nreals,ierrs(2))
  call extract('rhozero',rhozero,realarr,tags,nreals,ierrs(3))

!--extract required information from the first block header
  if (rhozero.gt.0.) then
     tfreefall = SQRT((3. * pi) / (32. * rhozero))
     tff = time/tfreefall
  else
     tfreefall = 0.
     tff = 0.
  endif
  if (phantomdump) then
     call extract('massoftype',massoftypei(1:ntypes),realarr,tags,nreals,ierrs(4))
     npartoftype(:) = 0
     do i=1,ntypes !--map from phantom types to splash types
        itype = itypemap_phantom(int(i,kind=1))
        if (debug) print*,'DEBUG: npart of type ',itype,' += ',npartoftypei(i)
        npartoftype(itype) = npartoftype(itype) + npartoftypei(i)
        massoftype(itype)  = massoftypei(i)
     enddo
     npartoftype(itypemap_sink_phantom) = nptmasstot  ! sink particles
     if (nblocks.gt.1) then
        print "(a)",' setting ngas=npart for MPI code '
        npartoftype(1)  = npart
        npartoftype(2:) = 0
     endif
     !
     !--if Phantom calculation uses the binary potential
     !  then read this as two point mass particles
     !
     if (nreals.ge.ilocbinary + 14) then
        if (nreals.ge.ilocbinary + 15) then
           ipos = ilocbinary
        else
           print*,'*** WARNING: obsolete header format for external binary information ***'
           ipos = ilocbinary + 1
        endif
        if (debug) print*,'DEBUG: reading binary information from header ',ilocbinary
        if (any(realarr(ilocbinary:ilocbinary+14).ne.0.)) then
           gotbinary = .true.
           npartoftype(itypemap_sink_phantom) = npartoftype(itypemap_sink_phantom) + 2
           ntotal = ntotal + 2
           dat(npart+1,ix(1)) = realarr(ipos)
           dat(npart+1,ix(2)) = realarr(ipos+1)
           dat(npart+1,ix(3)) = realarr(ipos+2)
           if (debug) print *,npart+1,npart+2
           if (iverbose.ge.1) print *,'binary position:   primary: ',realarr(ipos:ipos+2)
           if (nreals.ge.ilocbinary+15) then
              if (ipmass.gt.0) dat(npart+1,ipmass) = realarr(ipos+3)
              dat(npart+1,ih)     = realarr(ipos+4)
              dat(npart+2,ix(1))  = realarr(ipos+5)
              dat(npart+2,ix(2))  = realarr(ipos+6)
              dat(npart+2,ix(3))  = realarr(ipos+7)
              if (ipmass.gt.0) dat(npart+2,ipmass) = realarr(ipos+8)
              dat(npart+2,ih)     = realarr(ipos+9)
              if (iverbose.ge.1) then
                 print *,'                 secondary: ',realarr(ipos+5:ipos+7)
                 print *,' m1: ',realarr(ipos+3),' m2:',realarr(ipos+8),&
                         ' h1: ',realarr(ipos+4),' h2:',realarr(ipos+9)
              endif
              ipos = ipos + 10
           else
              dat(npart+1,ih)    = realarr(ipos+3)
              dat(npart+2,ix(1)) = realarr(ipos+4)
              dat(npart+2,ix(2)) = realarr(ipos+5)
              dat(npart+2,ix(3)) = realarr(ipos+6)
              dat(npart+2,ih)    = realarr(ipos+7)
              print *,'                 secondary: ', realarr(ipos+4:ipos+6)
              ipos = ipos + 8
           endif
           if (ivx.gt.0) then
              dat(npart+1,ivx)   = realarr(ipos)
              dat(npart+1,ivx+1) = realarr(ipos+1)
              dat(npart+1,ivx+2) = realarr(ipos+2)
              dat(npart+2,ivx)   = realarr(ipos+3)
              dat(npart+2,ivx+1) = realarr(ipos+4)
              dat(npart+2,ivx+2) = realarr(ipos+5)
           endif
           npart  = npart  + 2
        endif
     endif
  else
     npartoftype(:) = 0
     npartoftype(1) = npart
     npartoftype(2) = max(ntotal - npart,0)
  endif
  hfact = 1.2
  if (phantomdump) then
     call extract('hfact',hfact,realarr,tags,nreals,ierrs(1))
     call extract('tolh',tolh,realarr,tags,nreals,ierrs(2))
     print "(a,es12.4,a,f6.3,a,f5.2,a,es8.1)", &
           ' time = ',time,' gamma = ',gamma, &
           ' hfact = ',hfact,' tolh = ',tolh
  elseif (batcode) then
     call extract('radL1',radL1,realarr,tags,nreals,ierrs(1))
     call extract('PhiL1',PhiL1,realarr,tags,nreals,ierrs(2))
     call extract('Er',Er,realarr,tags,nreals,ierrs(3))
     print "(a,es12.4,a,f9.5,a,f8.4,/,a,es12.4,a,es9.2,a,es10.2)", &
           '   time: ',time,  '   gamma: ',gamma, '   tsph: ',realarr(2), &
           '  radL1: ',radL1,'   PhiL1: ',PhiL1,'     Er: ',Er
  else
     call extract('RK2',RK2,realarr,tags,nreals,ierrs(1))
     call extract('dtmax',dtmax,realarr,tags,nreals,ierrs(2))
     print "(a,es12.4,a,f9.5,a,f8.4,/,a,es12.4,a,es9.2,a,es10.2)", &
           '   time: ',time,  '   gamma: ',gamma, '   RK2: ',RK2, &
           ' t/t_ff: ',tff,' rhozero: ',rhozero,' dtmax: ',dtmax
  endif
 end subroutine extract_variables_from_header
 
!---------------------------------------------------------------
! old subroutine for guessing labels in non-tagged sphNG format
!---------------------------------------------------------------
 subroutine guess_labels(ncolumns,iamvec,label,labelvec,istartmhd,istart_extra_real4,nmhd,nhydroreal4, &
                         ndimV,irho,iBfirst,ivx,iutherm,idivB,iJfirst,iradenergy,icv,&
                         udist,utime,units,unitslabel)
  use geometry, only:labelcoord
  integer, intent(in) :: ncolumns,istartmhd,istart_extra_real4,nmhd,nhydroreal4,ndimV,irho
  integer, intent(out) :: iBfirst,ivx,iutherm,idivB,iJfirst,iradenergy,icv
  integer, intent(inout) :: iamvec(:)
  character(len=*), intent(inout) :: label(:),labelvec(:),unitslabel(:)
  real(doub_prec),  intent(in)    :: udist,utime
  real,    intent(inout) :: units(:)
  integer :: i
  real(doub_prec) :: uergg
  
!--the following only for mhd small dumps or full dumps
  if (ncolumns.ge.7) then
     if (mhddump) then
        iBfirst = irho+1
        if (.not.smalldump) then
           ivx = iBfirst+ndimV
           iutherm = ivx+ndimV

           if (phantomdump) then
              !--phantom MHD full dumps
              if (nmhd.ge.4) then
                 iamvec(istartmhd:istartmhd+ndimV-1) = istartmhd
                 labelvec(istartmhd:istartmhd+ndimV-1) = 'A'
                 do i=1,ndimV
                    label(istartmhd+i-1) = trim(labelvec(istartmhd))//'\d'//labelcoord(i,1)
                 enddo
                 if (nmhd.ge.7) then
                    label(istartmhd+3) = 'Euler beta\dx'
                    label(istartmhd+4) = 'Euler beta\dy'
                    label(istartmhd+5) = 'Euler beta\dz'
                    idivB = istartmhd+2*ndimV
                 else
                    idivB = istartmhd+ndimV
                 endif
              elseif (nmhd.ge.3) then
                 label(istartmhd) = 'Euler alpha'
                 label(istartmhd+1) = 'Euler beta'
                 idivB = istartmhd + 2
              elseif (nmhd.ge.2) then
                 label(istartmhd) = 'Psi'
                 idivB = istartmhd + 1
              elseif (nmhd.ge.1) then
                 idivB = istartmhd
              endif
              iJfirst = 0
              if (ncolumns.ge.idivB+1) then
                 label(idivB+1) = 'alpha\dB\u'
              endif

           else
              !--sphNG MHD full dumps
              label(iutherm+1) = 'grad h'
              label(iutherm+2) = 'grad soft'
              label(iutherm+3) = 'alpha'
              if (nmhd.ge.7 .and. usingvecp) then
                 iamvec(istartmhd:istartmhd+ndimV-1) = istartmhd
                 labelvec(istartmhd:istartmhd+ndimV-1) = 'A'
                 do i=1,ndimV
                    label(istartmhd+i-1) = trim(labelvec(16))//'\d'//labelcoord(i,1)
                 enddo
                 idivB = istartmhd+ndimV
              elseif (nmhd.ge.6 .and. usingeulr) then
                 label(istartmhd) = 'Euler alpha'
                 label(istartmhd+1) = 'Euler beta'
                 idivB = istartmhd + 2
              elseif (nmhd.ge.6) then
                 label(istartmhd) = 'psi'
                 idivB = istartmhd + 1
                 if (nmhd.ge.8) then
                    label(istartmhd+2+ndimV+1) = '\eta_{real}'
                    label(istartmhd+2+ndimV+2) = '\eta_{art}'
                    units(istartmhd+2+ndimV+1:istartmhd+2+ndimV+2) = udist*udist/utime
                    unitslabel(istartmhd+2+ndimV+1:istartmhd+2+ndimV+2) = ' [cm\u2\d/s]'
                 endif
                 if (nmhd.ge.14) then
                    label(istartmhd+2+ndimV+3) = 'fsym\dx'
                    label(istartmhd+2+ndimV+4) = 'fsym\dy'
                    label(istartmhd+2+ndimV+5) = 'fsym\dz'
                    labelvec(istartmhd+ndimV+5:istartmhd+ndimV+7) = 'fsym'
                    iamvec(istartmhd+ndimV+5:istartmhd+ndimV+7) = istartmhd+ndimV+5
                    label(istartmhd+2+ndimV+6) = 'faniso\dx'
                    label(istartmhd+2+ndimV+7) = 'faniso\dy'
                    label(istartmhd+2+ndimV+8) = 'faniso\dz'
                    labelvec(istartmhd+ndimV+8:istartmhd+ndimV+10) = 'faniso'
                    iamvec(istartmhd+ndimV+8:istartmhd+ndimV+10) = istartmhd+ndimV+8
                 endif
              elseif (nmhd.ge.1) then
                 idivB = istartmhd
              endif
              iJfirst = idivB + 1
              if (ncolumns.ge.iJfirst+ndimV) then
                 label(iJfirst+ndimV) = 'alpha\dB\u'
              endif
           endif
        else ! mhd small dump
           if (nhydroreal4.ge.3) iutherm = iBfirst+ndimV
        endif
     elseif (.not.smalldump) then
        ! pure hydro full dump
        ivx = irho+1
        iutherm = ivx + ndimV
        if (phantomdump) then
           if (istart_extra_real4.gt.0 .and. istart_extra_real4.lt.100) then
              label(istart_extra_real4) = 'alpha'
              label(istart_extra_real4+1) = 'alphau'
           endif
        else
           if (istart_extra_real4.gt.0 .and. istart_extra_real4.lt.100) then
              label(istart_extra_real4) = 'grad h'
              label(istart_extra_real4+1) = 'grad soft'
              label(istart_extra_real4+2) = 'alpha'
           endif
        endif
     endif

     if (phantomdump .and. h2chem) then
        if (smalldump) then
           label(nhydroarrays+nmhdarrays+1) = 'H_2 ratio'
        elseif (.not.smalldump .and. iutherm.gt.0) then
           label(iutherm+1) = 'H_2 ratio'
           label(iutherm+2) = 'HI abundance'
           label(iutherm+3) = 'proton abundance'
           label(iutherm+4) = 'e^- abundance'
           label(iutherm+5) = 'CO abundance'
        endif
     endif
     if (istartrt.gt.0 .and. istartrt.le.ncolumns .and. rtdump) then ! radiative transfer dump
        iradenergy = istartrt
        label(iradenergy) = 'radiation energy'
        uergg = (udist/utime)**2
        units(iradenergy) = uergg
        if (smalldump) then
           icv = istartrt+1
        else
           label(istartrt+1) = 'opacity'
           units(istartrt+1) = udist**2/umass
           icv = istartrt+2
           label(istartrt+3) = 'lambda'
           units(istartrt+3) = 1.0

           label(istartrt+4) = 'eddington factor'
           units(istartrt+4) = 1.0
        endif
       if (icv.gt.0) then
          label(icv) = 'u/T'
          units(icv) = uergg
       endif
    else
       iradenergy = 0
       icv = 0
    endif
  endif
 end subroutine guess_labels

 integer function assign_column(tag,iarr,ipos,ikind,imaxcolumnread) result(icolumn)
  use labels, only:ih,irho,ix,ipmass
  character(len=lentag), intent(in) :: tag
  integer,               intent(in) :: iarr,ipos,ikind
  integer,               intent(inout) :: imaxcolumnread

  if (tagged .and. len_trim(tag) > 0) then
     !
     ! use the tags to put certain arrays in an assigned place
     ! no matter what type is used for the variable in the file
     ! and no matter what order they appear in the dump file
     !
     select case(trim(tag))
     case('x')
        icolumn = ix(1)
     case('y')
        icolumn = ix(2)
     case('z')
        icolumn = ix(3)
     case('m')
        icolumn = ipmass
     case('h')
        icolumn = ih
     case('rho')
        icolumn = irho
     case('Bx')
        icolumn = nhydroarrays + 1
     case('By')
        icolumn = nhydroarrays + 2
     case('Bz')
        icolumn = nhydroarrays + 3
     case default
        icolumn = max(nhydroarrays + nmhdarrays + 1,imaxcolumnread + 1)
        if (iarr==1) then
           if (ikind==4) then  ! real*4 array
              istart_extra_real4 = min(istart_extra_real4,icolumn)
              if (debug) print*,' istart_extra_real4 = ',istart_extra_real4
           endif
        endif
     end select
  else
     !
     ! this is old code handling the non-tagged format where
     ! particular arrays are assumed to be in particular places
     !
     if (ikind==6) then   ! default reals
        if (iarr.eq.1.and.((phantomdump.and.ipos.eq.4) &
           .or.(.not.phantomdump.and.ipos.eq.6))) then
           ! read x,y,z,m,h and then place arrays after always-present ones
           ! (for phantom read x,y,z only)
           icolumn = nhydroarrays+nmhdarrays + 1
        elseif (.not.phantomdump .and. (iarr.eq.4 .and. ipos.le.3)) then
           icolumn = nhydroarrays + ipos
        else
           icolumn = imaxcolumnread + 1
        endif
     elseif (ikind==4) then  ! real*4s
        if (phantomdump) then
           if (iarr.eq.1 .and. ipos.eq.1) then
              icolumn = ih ! h is always first real4 in phantom dumps
              !!--density depends on h being read
              !required(ih) = .true.
           elseif (iarr.eq.4 .and. ipos.le.3) then
              icolumn = nhydroarrays + ipos
           else
              icolumn = max(nhydroarrays+nmhdarrays + 1,imaxcolumnread + 1)
              if (iarr.eq.1) then
                 istart_extra_real4 = min(istart_extra_real4,icolumn)
                 if (debug) print*,' istart_extra_real4 = ',istart_extra_real4
              endif
           endif
        else
           if (iarr.eq.1 .and. ipos.eq.1) then
              icolumn = irho ! density
           elseif (iarr.eq.1 .and. smalldump .and. ipos.eq.2) then
              icolumn = ih ! h which is real4 in small dumps
           !--this was a bug for sphNG files...
           !elseif (iarr.eq.4 .and. i.le.3) then
           !   icolumn = nhydroarrays + ipos
           else
              icolumn = max(nhydroarrays+nmhdarrays + 1,imaxcolumnread + 1)
              if (iarr.eq.1) then
                 istart_extra_real4 = min(istart_extra_real4,icolumn)
                 if (debug) print*,' istart_extra_real4 = ',istart_extra_real4
              endif
           endif
        endif
     else ! used for untagged format with real*8's
        icolumn = imaxcolumnread + 1
     endif
  endif
  imaxcolumnread = max(imaxcolumnread,icolumn)

 end function assign_column

end module sphNGread

!----------------------------------------------------------------------
!  Main read_data routine for splash
!----------------------------------------------------------------------
subroutine read_data(rootname,indexstart,iposn,nstepsread)
  use particle_data,  only:dat,gamma,time,iamtype,npartoftype,maxpart,maxstep,maxcol,masstype
  !use params,         only:int1,int8
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc,required,ipartialread,&
                      lowmemorymode,ntypes,iverbose
  use mem_allocation, only:alloc
  use system_utils,   only:lenvironment,renvironment
  use labels,         only:ipmass,irho,ih,ix,ivx,labeltype,print_types
  use calcquantities, only:calc_quantities
  use sphNGread
  implicit none
  integer, intent(in)  :: indexstart,iposn
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: i,j,k,ierr,iunit
  integer :: intg1,int2,int3,ilocvx,iversion
  integer :: i1,iarr,i2,iptmass1,iptmass2
  integer :: npart_max,nstep_max,ncolstep,icolumn,nptmasstot
  integer :: narrsizes
  integer :: nskip,ntotal,npart,n1,ngas,nreals
  integer :: iblock,nblocks,ntotblock,ncolcopy
  integer :: ipos,nptmass,nptmassi,nstar,nunknown,ilastrequired
  integer :: imaxcolumnread,nhydroarraysinfile,nremoved
  integer :: itype,iphaseminthistype,iphasemaxthistype,nthistype,iloc
  integer, dimension(maxparttypes) :: npartoftypei
  real,    dimension(maxparttypes) :: massoftypei
  real    :: pmassi,hi,rhoi,hrlim,rad2d
  logical :: iexist, doubleprec,imadepmasscolumn,gotbinary,gotiphase

  character(len=len(rootname)+10) :: dumpfile
  character(len=100) :: fileident

  integer*8, dimension(maxarrsizes) :: isize
  integer, dimension(maxarrsizes) :: nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8
  integer*1, dimension(:), allocatable :: iphase
  integer, dimension(:), allocatable :: listpm
  real(doub_prec), dimension(:), allocatable :: dattemp
  real*4, dimension(:), allocatable :: dattempsingle
  real(doub_prec) :: r8
  real(sing_prec) :: r4
  real, dimension(:,:), allocatable :: dattemp2
  real, dimension(maxinblock) :: dummyreal
  real :: hfact,omega
  logical :: skip_corrupted_block_3
  character(len=lentag) :: tagsreal(maxinblock), tagtmp
  
  integer, parameter :: splash_max_iversion = 1

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  npart = 0
  iunit = 15
  ipmass = 4
  idivvcol = 0
  icurlvxcol = 0
  icurlvycol = 0
  icurlvzcol = 0
  nhydroreal4 = 0
  umass = 1.d0
  utime = 1.d0
  udist = 1.d0
  umagfd = 1.d0
  istartmhd = 0
  istartrt  = 0
  istart_extra_real4 = 100
  nmhd      = 0
  igotmass    = .false.
  tfreefall   = 1.d0
  gotbinary   = .false.
  gotiphase   = .false.
  skip_corrupted_block_3 = .false.

  dumpfile = trim(rootname)
  !
  !--check if data file exists
  !
  inquire(file=dumpfile,exist=iexist)
  if (.not.iexist) then
     print "(a)",' *** error: '//trim(dumpfile)//': file not found ***'
     return
  endif
  !
  !--fix number of spatial dimensions
  !
  ndim = 3
  ndimV = 3

  j = indexstart
  nstepsread = 0
  doubleprec = .true.
  ilastrequired = 0
  do i=1,size(required)-1
     if (required(i)) ilastrequired = i
  enddo

  if (iverbose.ge.1) print "(1x,a)",'reading sphNG format'
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)

  debug = lenvironment('SSPLASH_DEBUG')
  if (debug) iverbose = 1
!
!--open the (unformatted) binary file
!
   open(unit=iunit,iostat=ierr,file=dumpfile,status='old',form='unformatted')
   if (ierr /= 0) then
      print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
      return
   else
      !
      !--read header key to work out precision
      !
      doubleprec = .true.
      read(iunit,iostat=ierr) intg1,r8,int2,iversion,int3
      if (intg1.ne.690706 .and. intg1.ne.060769) then
         print "(a)",'*** ERROR READING HEADER: corrupt file/zero size/wrong endian?'
         close(iunit)
         return
      endif
      if (int2.ne.780806 .and. int2.ne.060878) then
         print "(a)",' single precision dump'
         rewind(iunit)
         read(iunit,iostat=ierr) intg1,r4,int2,iversion,int3
         if (int2.ne.780806 .and. int2.ne.060878) then
            print "(a)",'ERROR determining single/double precision in file header'
         endif
         doubleprec = .false.
      elseif (int3.ne.690706) then
          print*,' got ',intg1,r4,int2,iversion,int3
         print "(a)",'*** WARNING: default int appears to be int*8: not implemented'
      else
         if (debug) print "(a)",' double precision dump' ! no need to print this
      endif
      if (iversion==690706) then ! handle old-format files (without version number) gracefully
         iversion = 0
      endif
   endif
   if (iversion > splash_max_iversion) then
      print "(/a,i2,/,a,i2)",&
        ' *** WARNING: this copy of splash can only read version ',splash_max_iversion, &
        '              but the file format version is ',iversion
      if (.not.lenvironment('SSPLASH_IGNORE_IVERSION')) then
         print "(2(/,a))",'   ** press any key to bravely proceed anyway ** ', &
                          '   (set SSPLASH_IGNORE_IVERSION=yes to silence this warning)'
         read*
      endif
   endif
!
!--read file ID
!
   read(iunit,iostat=ierr) fileident
   if (ierr /=0) then
      print "(a)",'*** ERROR READING FILE ID ***'
      close(iunit)
      return
   else
      print "(a)",' File ID: '//trim(fileident)
   endif
   mhddump = .false.
   rtdump = .false.
   call get_options_from_fileident(fileident,smalldump,tagged,phantomdump,&
                                   usingvecp,usingeulr,cleaning,h2chem,rt_in_header,batcode)
   if (tagged .and. iversion < 1) print "(a)",'ERROR: got tagged format but iversion is ',iversion
!
!--read variables from header
!
   call read_header(iunit,iverbose,debug,doubleprec, &
                    npart,npartoftypei,n1,ntypes,nblocks,narrsizes,dummyreal,tagsreal,nreals,ierr)
   if (ierr /= 0) then
      print "(a)",' *** ERROR READING HEADER ***'
      close(iunit)
      return
   endif
!
!--Attempt to read all MPI blocks
!
   ntotal = 0
   ntotblock = 0
   nptmasstot = 0
   i2 = 0
   iptmass2 = 0
   igotmass = .true.
   imadepmasscolumn = .false.
   massoftypei(:) = 0.

   over_MPIblocks: do iblock=1,nblocks
!
!--read array header from this block
!
   if (iblock.eq.1) ncolstep = 0
   do iarr=1,narrsizes
      call read_block_header(iunit,iblock,iarr,iverbose,debug, &
           isize,nint(iarr),nint1(iarr),nint2(iarr),nint4(iarr),nint8(iarr),&
           nreal(iarr),nreal4(iarr),nreal8(iarr),&
           ntotblock,npart,ntotal,nptmasstot,ncolstep,ierr)
      if (ierr /= 0) then
         print "(a)",' *** ERROR READING ARRAY SIZES ***'
         close(iunit)
         return
      endif
   enddo
   if (debug) print*,'DEBUG: ncolstep=',ncolstep,' from file header, also nptmasstot = ',nptmasstot
!
!--this is a bug fix for a corrupt version of wdump outputting bad
!  small dump files
!
   if (smalldump .and. nreal(1).eq.5 .and. iblock.eq.1 .and. lenvironment('SSPLASH_FIX_CORRUPT')) then
      print*,'FIXING CORRUPT HEADER ON SMALL DUMPS: assuming nreal=3 not 5'
      nreal(1) = 3
      ncolstep = ncolstep - 2
   endif

   npart_max = maxval(isize(1:narrsizes))
   npart_max = max(npart_max,npart+nptmasstot,ntotal)
!
!--work out from array header how many columns we are going to read
!  in order to allocate memory
!
   if (iblock.eq.1) then
      igotmass = .true.
      if (smalldump .or. phantomdump) then
         if (phantomdump) then
            if (tagged) then
               call extract('massoftype',massoftypei(1:ntypes),dummyreal,tagsreal,nreals,ierr)
            else
            ! old phantom dumps had only 5 types
               call extract('massoftype',massoftypei(1:5),dummyreal,tagsreal,nreals,ierr)
            endif
         else
            call extract('pmassinitial',massoftypei(1),dummyreal,tagsreal,nreals,ierr)
            if (ierr /= 0) then
               print "(a)",' error extracting particle mass from small dump file'
               massoftypei(1) = 0.
               igotmass = .false.
            endif
         endif
         if (debug) print*,'DEBUG: got massoftype(gas) = ',massoftypei(1)
         if (any(massoftypei(1:ntypes).gt.tiny(0.)) .and. .not.lowmemorymode) then
            ncolstep = ncolstep + 1  ! make an extra column to contain particle mass
            imadepmasscolumn = .true.
         elseif (lowmemorymode) then
            igotmass = .false.
         else
            igotmass = .false.
         endif
         if (all(abs(massoftypei(1:ntypes)).lt.tiny(0.)) .and. nreal(1).lt.4) then
            print "(a)",' error: particle masses not present in small dump file'
            igotmass = .false.
         endif
      endif
      if (debug) print*,'DEBUG: gotmass = ',igotmass, ' ncolstep = ',ncolstep
!
!--   to handle both small and full dumps, we need to place the quantities dumped
!     in both small and full dumps at the start of the dat array
!     quantities only in the full dump then come after
!     also means that hydro/MHD are "semi-compatible" in the sense that x,y,z,m,h and rho
!     are in the same place for both types of dump
!
      ix(1) = 1
      ix(2) = 2
      ix(3) = 3
      if (igotmass) then
         ipmass = 4
         ih = 5
         irho = 6
         nhydroarrays = 6 ! x,y,z,m,h,rho
      else
         ipmass = 0
         ih = 4
         irho = 5
         nhydroarrays = 5 ! x,y,z,h,rho
      endif
      nhydroarraysinfile = nreal(1) + nreal4(1) + nreal8(1)
      nhydroreal4 = nreal4(1)
      if (imadepmasscolumn) nhydroarraysinfile = nhydroarraysinfile + 1
      if (nhydroarraysinfile .lt.nhydroarrays .and. .not.phantomdump) then
         print "(a)",' ERROR: one of x,y,z,m,h or rho missing in small dump read'
         nhydroarrays = nreal(1)+nreal4(1)+nreal8(1)
      elseif (phantomdump .and. (nreal(1).lt.3 .or. nreal4(1).lt.1)) then
         print "(a)",' ERROR: x,y,z or h missing in phantom read'
      endif
      if (narrsizes.ge.4) then
         nmhdarrays = 3 ! Bx,By,Bz
         nmhd = nreal(4) + nreal4(4) + nreal8(4) - nmhdarrays ! how many "extra" mhd arrays
         if (debug) print*,'DEBUG: ',nmhd,' extra MHD arrays'
      else
         nmhdarrays = 0
      endif

      !--radiative transfer dump?
      if (narrsizes.ge.3 .and. isize(3).eq.isize(1)) rtdump = .true.
      !--mhd dump?
      if (narrsizes.ge.4) mhddump = .true.

      if (.not.(mhddump.or.smalldump)) then
         ivx = nhydroarrays+1
      elseif (mhddump .and. .not.smalldump) then
         ivx = nhydroarrays+nmhdarrays+1
      else
         ivx = 0
      endif
      !--need to force read of velocities e.g. for corotating frame subtraction
      if (any(required(ivx:ivx+ndimV-1))) required(ivx:ivx+ndimV-1) = .true.

      !--for phantom dumps, also make a column for density
      !  and divv, if a .divv file exists
      if (phantomdump) then
         ncolstep = ncolstep + 1
         inquire(file=trim(dumpfile)//'.divv',exist=iexist)
         if (iexist) then
            idivvcol   = ncolstep + 1
            icurlvxcol = ncolstep + 2
            icurlvycol = ncolstep + 3
            icurlvzcol = ncolstep + 4
            ncolstep   = ncolstep + 4
         endif
      endif
   endif
!
!--allocate memory now that we know the number of columns
!
   if (iblock.eq.1) then
      ncolumns = ncolstep + ncalc
      if (ncolumns.gt.maxplot) then
         print*,'ERROR with ncolumns = ',ncolumns,' in data read'
         return
      endif
      ilastrequired = 0
      do i=1,ncolumns
         if (required(i)) ilastrequired = i
      enddo
   endif

   if (npart_max.gt.maxpart .or. j.gt.maxstep .or. ncolumns.gt.maxcol) then
      if (lowmemorymode) then
         call alloc(max(npart_max+2,maxpart),j,ilastrequired)
      else
         call alloc(max(npart_max+2,maxpart),j,ncolumns,mixedtypes=.true.)
      endif
   endif
!
!--now that memory has been allocated, copy info from the header into
!  the relevant arrays
!
   if (iblock.eq.1) then
      call extract_variables_from_header(tagsreal,dummyreal,nreals,iverbose,debug, &
           gotbinary,nblocks,nptmasstot,npartoftypei,ntypes,&
           time(j),gamma(j),hfact,npart,ntotal,npartoftype(:,j),masstype(:,j), &
           dat(:,:,j),ix,ih,ipmass,ivx)

      nstepsread = nstepsread + 1
      !
      !--stop reading file here if no columns required
      !
      if (ilastrequired.eq.0) exit over_MPIblocks

      if (allocated(iphase)) deallocate(iphase)
      allocate(iphase(npart_max+2))
      if (phantomdump) then
         iphase(:) = 1
      else
         iphase(:) = 0
      endif

      if (gotbinary) then
         iphase(npart-1) = -3
         iphase(npart)   = -3
      endif
   endif
!
!--Arrays
!
   imaxcolumnread = 0
   icolumn = 0
   istartmhd = 0
   istartrt = 0
   i1 = i2 + 1
   i2 = i1 + isize(1) - 1
   if (debug) then
      print "(1x,a10,i4,3(a,i12))",'MPI block ',iblock,':  particles: ',i1,' to ',i2,' of ',npart
   elseif (nblocks.gt.1) then
      if (iblock.eq.1) write(*,"(a,i1,a)",ADVANCE="no") ' reading MPI blocks: .'
      write(*,"('.')",ADVANCE="no")
   endif
   iptmass1 = iptmass2 + 1
   iptmass2 = iptmass1 + isize(2) - 1
   nptmass = nptmasstot
   if (nptmass.gt.0 .and. debug) print "(15x,3(a,i12))",'  pt. masses: ',iptmass1,' to ',iptmass2,' of ',nptmass

   do iarr=1,narrsizes
      if (nreal(iarr) + nreal4(iarr) + nreal8(iarr).gt.0) then
         if (iarr.eq.4) then
            istartmhd = imaxcolumnread + 1
            if (debug) print*,' istartmhd = ',istartmhd
         elseif (iarr.eq.3 .and. rtdump) then
            istartrt = max(nhydroarrays+nmhdarrays+1,imaxcolumnread + 1)
            if (debug) print*,' istartrt = ',istartrt
         endif
      endif
!--read iphase from array block 1
      if (iarr.eq.1) then
         !--skip default int
         nskip = nint(iarr)
         do i=1,nskip
            if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
            read(iunit,end=33,iostat=ierr)
         enddo
         if (nint1(iarr).lt.1) then
            if (.not.phantomdump .or. any(npartoftypei(2:).gt.0)) then
               print "(a)",' WARNING: can''t locate iphase in dump'
            elseif (phantomdump) then
               print "(a)",' WARNING: can''t locate iphase in dump'
            endif
            gotiphase = .false.
            !--skip remaining integer arrays
            nskip = nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
         else
            gotiphase = .true.
            if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
            read(iunit,end=33,iostat=ierr) iphase(i1:i2)
            !--skip remaining integer arrays
            nskip = nint1(iarr) - 1 + nint2(iarr) + nint4(iarr) + nint8(iarr)
         endif
      elseif (smalldump .and. iarr.eq.2 .and. isize(iarr).gt.0 .and. .not.phantomdump) then
!--read listpm from array block 2 for small dumps (needed here to extract sink masses)
         if (allocated(listpm)) deallocate(listpm)
         allocate(listpm(isize(iarr)))
         if (nint(iarr).lt.1) then
            print "(a)",'ERROR: can''t locate listpm in dump'
            nskip = nint(iarr) + nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
         else
            if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
            read(iunit,end=33,iostat=ierr) listpm(1:isize(iarr))
            nskip = nint(iarr) - 1 + nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
         endif
      else
!--otherwise skip all integer arrays (not needed for plotting)
         nskip = nint(iarr) + nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
      endif

      if (iarr.eq.3 .and. lenvironment('SSPLASH_BEN_HACKED')) then
         nskip = nskip - 1
         print*,' FIXING HACKED DUMP FILE'
      endif
      !print*,'skipping ',nskip
      do i=1,nskip
         if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
         read(iunit,end=33,iostat=ierr)
      enddo
!
!--real arrays
!
      if (iarr.eq.2) then
!--read sink particles from phantom dumps
         if (phantomdump .and. iarr.eq.2 .and. isize(iarr).gt.0) then
            if (nreal(iarr).lt.5) then
               print "(a)",'ERROR: not enough arrays written for sink particles in phantom dump'
               nskip = nreal(iarr)
            else
               iphase(npart+1:npart+isize(iarr)) = -3
               ilocvx = nreal(iarr)-2 ! velocity is always last 3 numbers for phantom sinks
               if (doubleprec) then
                  !--convert default real to single precision where necessary
                  if (debug) print*,'DEBUG: reading sink data, converting from double precision ',isize(iarr)
                  if (allocated(dattemp)) deallocate(dattemp)
                  allocate(dattemp(isize(iarr)),stat=ierr)
                  if (ierr /= 0) then
                     print "(a)",'ERROR in memory allocation'
                     return
                  endif
                  tagtmp = ''
                  do k=1,nreal(iarr)
                     if (tagged) read(iunit,end=33,iostat=ierr) tagtmp
                     if (debug) print*,'DEBUG: reading sink array ',k,isize(iarr),' tag = ',trim(tagtmp)
                     read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
                     if (ierr /= 0) print*,' ERROR during read of sink particle data, array ',k
                     
                     select case(k)
                     case(1:3)
                        iloc = ix(k)
                     case(4)
                        iloc = ipmass
                     case(5)
                        iloc = ih
                     case default
                        if (k >= ilocvx .and. k < ilocvx+3 .and. ivx > 0) then
                           iloc = ivx + k-ilocvx ! put velocity into correct arrays
                        else
                           iloc = 0
                        endif
                     end select
                     if (iloc.gt.size(dat(1,:,j))) then; print*,' error iloc = ',iloc,ivx; stop; endif
                     if (iloc.gt.0) then
                        do i=1,isize(iarr)
                           dat(npart+i,iloc,j) = real(dattemp(i))
                        enddo
                     else
                        if (debug) print*,'DEBUG: skipping sink particle array ',k
                     endif
                  enddo
               else
                  if (debug) print*,'DEBUG: reading sink data, converting from single precision ',isize(iarr)
                  if (allocated(dattempsingle)) deallocate(dattempsingle)
                  allocate(dattempsingle(isize(iarr)),stat=ierr)
                  if (ierr /= 0) then
                     print "(a)",'ERROR in memory allocation'
                     return
                  endif
                  do k=1,nreal(iarr)
                     select case(k)
                     case(1:3)
                        iloc = ix(k)
                     case(4)
                        iloc = ipmass
                     case(5)
                        iloc = ih
                     case default
                        if (k >= ilocvx .and. k < ilocvx+3 .and. ivx > 0) then
                           iloc = ivx + k-ilocvx ! put velocity into correct arrays
                        else
                           iloc = 0
                        endif
                     end select
                     if (iloc.gt.0) then
                        if (debug) print*,'DEBUG: reading sinks into ',npart+1,'->',npart+isize(iarr),iloc
                        if (tagged) read(iunit,end=33,iostat=ierr) !tagarr(iloc)
                        read(iunit,end=33,iostat=ierr) dattempsingle(1:isize(iarr))
                        do i=1,isize(iarr)
                           dat(npart+i,iloc,j) = real(dattempsingle(i))
                        enddo
                        if (ierr /= 0) print*,' ERROR during read of sink particle data, array ',k
                     else
                        if (debug) print*,'DEBUG: skipping sink particle array ',k
                        if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
                        read(iunit,end=33,iostat=ierr)
                     endif
                  enddo
               endif
               npart  = npart + isize(iarr)
            endif
         elseif (smalldump .and. iarr.eq.2 .and. allocated(listpm)) then
!--for sphNG, read sink particle masses from block 2 for small dumps
            if (nreal(iarr).lt.1) then
               if (isize(iarr).gt.0) print "(a)",'ERROR: sink masses not present in small dump'
               nskip = nreal(iarr) + nreal4(iarr) + nreal8(iarr)
            else
               if (doubleprec) then
                  !--convert default real to single precision where necessary
                  if (allocated(dattemp)) deallocate(dattemp)
                  allocate(dattemp(isize(iarr)),stat=ierr)
                  if (ierr /=0) print "(a)",'ERROR in memory allocation'
                  if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
                  read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
                  if (nptmass.ne.isize(iarr)) print "(a)",'ERROR: nptmass.ne.block size'
                  if (ipmass.gt.0) then
                     do i=1,isize(iarr)
                        dat(listpm(iptmass1+i-1),ipmass,j) = real(dattemp(i))
                     enddo
                  else
                     print*,'WARNING: sink particle masses not read because no mass array allocated'
                  endif
               else
                  !--convert default real to double precision where necessary
                  if (allocated(dattempsingle)) deallocate(dattempsingle)
                  allocate(dattempsingle(isize(iarr)),stat=ierr)
                  if (ierr /=0) print "(a)",'ERROR in memory allocation'
                  if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
                  read(iunit,end=33,iostat=ierr) dattempsingle(1:isize(iarr))
                  if (nptmass.ne.isize(iarr)) print "(a)",'ERROR: nptmass.ne.block size'
                  if (ipmass.gt.0) then
                     do i=1,isize(iarr)
                        dat(listpm(iptmass1+i-1),ipmass,j) = real(dattempsingle(i))
                     enddo
                  else
                     print*,'WARNING: sink particle masses not read because no mass array allocated'
                  endif
               endif
               nskip = nreal(iarr) - 1 + nreal4(iarr) + nreal8(iarr)
            endif
         else
!--for other blocks, skip real arrays if size different
            nskip = nreal(iarr) + nreal4(iarr) + nreal8(iarr)
         endif
         do i=1,nskip
            if (tagged) read(iunit,end=33,iostat=ierr) ! skip tags
            read(iunit,end=33,iostat=ierr)
         enddo
         ! deallocate dattempsingle
         if (allocated(dattempsingle)) deallocate(dattempsingle)

      elseif (isize(iarr).eq.isize(1)) then
!
!--read all real arrays defined on all the particles (same size arrays as block 1)
!
         if ((doubleprec.and.nreal(iarr).gt.0).or.nreal8(iarr).gt.0) then
            if (allocated(dattemp)) deallocate(dattemp)
            allocate(dattemp(isize(iarr)),stat=ierr)
            if (ierr /=0) print "(a)",'ERROR in memory allocation (read_data_sphNG: dattemp)'
         elseif (nreal(iarr).gt.0 .or. nreal8(iarr).gt.0) then
            if (allocated(dattempsingle)) deallocate(dattempsingle)
            allocate(dattempsingle(isize(iarr)),stat=ierr)
            if (ierr /=0) print "(a)",'ERROR in memory allocation (read_data_sphNG: dattempsingle)'
         endif
!        default reals may need converting
         do i=1,nreal(iarr)
            tagtmp = ''
            if (tagged) read(iunit,end=33,iostat=ierr) tagtmp
            icolumn = assign_column(tagtmp,iarr,i,6,imaxcolumnread)
            if (tagged) tagarr(icolumn) = tagtmp
            if (debug)  print*,' reading real ',icolumn,' tag = ',trim(tagtmp)
            if (required(icolumn)) then
               if (doubleprec) then
                  read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
                  dat(i1:i2,icolumn,j) = real(dattemp(1:isize(iarr)))
               else
                  read(iunit,end=33,iostat=ierr) dattempsingle(1:isize(iarr))
                  dat(i1:i2,icolumn,j) = real(dattempsingle(1:isize(iarr)))
               endif
            else
               read(iunit,end=33,iostat=ierr)
            endif
         enddo
!
!        set masses for equal mass particles (not dumped in small dump or in phantom)
!
         if (((smalldump.and.nreal(1).lt.ipmass).or.phantomdump).and. iarr.eq.1) then
            if (abs(masstype(1,j)).gt.tiny(masstype)) then
               icolumn = ipmass
               if (required(ipmass) .and. ipmass.gt.0) then
                  if (phantomdump) then
                     dat(i1:i2,ipmass,j) = masstype(itypemap_phantom(iphase(i1:i2)),j)
                  else
                     where (iphase(i1:i2).eq.0) dat(i1:i2,icolumn,j) = masstype(1,j)
                  endif
               endif
               !--dust mass for phantom particles
               if (phantomdump .and. npartoftypei(itypemap_dust_phantom).gt.0 .and. ipmass.gt.0) then
                  print*,'dust particle mass = ',masstype(itypemap_dust_phantom,j),&
                         ' ratio m_dust/m_gas = ',masstype(itypemap_dust_phantom,j)/masstype(1,j)
               endif
               if (debug) print*,'mass ',icolumn
            elseif (phantomdump .and. npartoftypei(1).gt.0) then
               print*,' ERROR: particle mass zero in Phantom dump file!'
            endif
         endif
!
!        real4 arrays (may need converting if splash is compiled in double precision)
! 
         if (nreal4(iarr).gt.0 .and. kind(dat).eq.doub_prec) then
            if (allocated(dattempsingle)) deallocate(dattempsingle)
            allocate(dattempsingle(isize(iarr)),stat=ierr)
            if (ierr /=0) print "(a)",'ERROR in memory allocation (read_data_sphNG: dattempsingle)'
         endif

         if (debug) print*,'DEBUG: SIZE of dattempsingle',size(dattempsingle)
!        real4s may need converting
         imaxcolumnread = max(imaxcolumnread,icolumn)
         if ((nreal(iarr)+nreal4(iarr)).gt.6) imaxcolumnread = max(imaxcolumnread,6)

         do i=1,nreal4(iarr)
            tagtmp = ''
            if (tagged) read(iunit,end=33,iostat=ierr) tagtmp
            icolumn = assign_column(tagtmp,iarr,i,4,imaxcolumnread)
            if (debug) print*,'reading real4 ',icolumn,' tag = ',trim(tagtmp)
            if (tagged) tagarr(icolumn) = tagtmp

            if (phantomdump .and. icolumn==ih) required(ih) = .true. ! h always required for density

            if (required(icolumn)) then
               if (allocated(dattempsingle)) then
                  read(iunit,end=33,iostat=ierr) dattempsingle(1:isize(iarr))
                  dat(i1:i2,icolumn,j) = real(dattempsingle(1:isize(iarr)))
               else
                  read(iunit,end=33,iostat=ierr) dat(i1:i2,icolumn,j)
               endif
            else
               read(iunit,end=33,iostat=ierr)
            endif
            !--construct density for phantom dumps based on h, hfact and particle mass
            if (phantomdump .and. icolumn.eq.ih) then
               icolumn = irho ! density
               !
               !--dead particles have -ve smoothing lengths in phantom
               !  so use abs(h) for these particles and hide them
               !
               if (any(npartoftypei(2:).gt.0)) then
                  if (.not.required(ih)) print*,'ERROR: need to read h, but required=F'
                  !--need masses for each type if not all gas
                  if (debug) print*,'DEBUG: phantom: setting h for multiple types ',i1,i2
                  if (debug) print*,'DEBUG: massoftype = ',masstype(:,j)
                  do k=i1,i2
                     itype = itypemap_phantom(iphase(k))
                     pmassi = masstype(itype,j)
                     hi = dat(k,ih,j)
                     if (hi > 0.) then
                        if (required(irho)) dat(k,irho,j) = pmassi*(hfact/hi)**3
                     elseif (hi < 0.) then
                        npartoftype(itype,j) = npartoftype(itype,j) - 1
                        npartoftype(itypemap_unknown_phantom,j) = npartoftype(itypemap_unknown_phantom,j) + 1
                        if (required(irho)) dat(k,irho,j) = pmassi*(hfact/abs(hi))**3
                     else
                        if (required(irho)) dat(k,irho,j) = 0.
                     endif
                  enddo
               else
                  if (.not.required(ih)) print*,'ERROR: need to read h, but required=F'
                  if (debug) print*,'debug: phantom: setting rho for all types'
                  !--assume all particles are gas particles
                  do k=i1,i2
                     hi = dat(k,ih,j)
                     if (hi.gt.0.) then
                        rhoi = massoftypei(1)*(hfact/hi)**3
                     elseif (hi.lt.0.) then
                        rhoi = massoftypei(1)*(hfact/abs(hi))**3
                        iphase(k) = -1
                     else ! if h = 0.
                        rhoi = 0.
                        iphase(k) = -2
                     endif
                     if (required(irho)) dat(k,irho,j) = rhoi
                  enddo
               endif

               if (debug) print*,'debug: making density ',icolumn
            endif
         enddo
!        real 8's need converting
         do i=1,nreal8(iarr)
            tagtmp = ''
            if (tagged) read(iunit,end=33,iostat=ierr) tagtmp
            icolumn = assign_column(tagtmp,iarr,i,8,imaxcolumnread)
            if (debug) print*,'reading real8 ',icolumn,' tag = ',trim(tagtmp)
            if (tagged) tagarr(icolumn) = tagtmp
            if (required(icolumn)) then
               read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
               dat(i1:i2,icolumn,j) = real(dattemp(1:isize(iarr)))
            else
               read(iunit,end=33,iostat=ierr)
            endif
         enddo
      endif
   enddo ! over array sizes
   enddo over_MPIblocks
!
!--reached end of file (during data read)
!
   goto 34
33 continue
   print "(/,1x,a,/)",'*** WARNING: END OF FILE DURING READ ***'
   print*,'Press any key to continue (but there is likely something wrong with the file...)'
   read*
34 continue
 !
 !--read .divv file for phantom dumps
 !
    if (phantomdump .and. idivvcol.ne.0 .and. any(required(idivvcol:icurlvzcol))) then
       print "(a)",' reading divv from '//trim(dumpfile)//'.divv'
       open(unit=66,file=trim(dumpfile)//'.divv',form='unformatted',status='old',iostat=ierr)
       if (ierr /= 0) then
          print "(a)",' ERROR opening '//trim(dumpfile)//'.divv'
       else
          read(66,iostat=ierr) dat(1:ntotal,idivvcol,j)
          if (ierr /= 0) print "(a)",' WARNING: ERRORS reading divv from file'
          if (any(required(icurlvxcol:icurlvzcol))) then
             read(66,iostat=ierr) dat(1:ntotal,icurlvxcol,j)
             read(66,iostat=ierr) dat(1:ntotal,icurlvycol,j)
             read(66,iostat=ierr) dat(1:ntotal,icurlvzcol,j)
          endif
          if (ierr /= 0) print "(a)",' WARNING: ERRORS reading curlv from file'
          close(66)
       endif
    endif
 !
 !--reset centre of mass to zero if environment variable "SSPLASH_RESET_CM" is set
 !
    if (allocated(dat) .and. n1.GT.0 .and. n1 <= size(dat(:,1,1)) &
       .and. lenvironment('SSPLASH_RESET_CM') .and. allocated(iphase)) then
       call reset_centre_of_mass(dat(1:n1,1:3,j),dat(1:n1,4,j),iphase(1:n1),n1)
    endif
 ! 
 !--remove particles at large H/R is "SSPLASH_REMOVE_LARGE_HR" is set
 !
    if (lenvironment('SSPLASH_REMOVE_LARGE_HR')) then
       hrlim = renvironment('SSPLASH_HR_LIMIT')
       print "(a)", 'SSPLASH_REMOVE_LARGE_HR set:'
       print "(a)", 'Removing particles at large H/R values'
       print "(a,F7.4)", 'H/R limit set to ',hrlim
       
       nremoved = 0
       do i = 1,npart
          if (int(iphase(i)) == 0) then
             rad2d = sqrt(dat(i,1,j)**2 + dat(i,2,j)**2)
             if (abs(dat(i,3,j) / rad2d) >= hrlim) then
                iphase(i) = -1
                nremoved  = nremoved + 1
             endif
          endif
       enddo
       print "(I5,a)", nremoved, ' particles removed at large H/R'
    endif
 !
 !--reset corotating frame velocities if environment variable "SSPLASH_OMEGA" is set
 !
    if (allocated(dat) .and. n1.GT.0 .and. all(required(1:2))) then
       omega = renvironment('SSPLASH_OMEGAT')
       if (abs(omega).gt.tiny(omega) .and. ndim.ge.2) then
          call reset_corotating_positions(n1,dat(1:n1,1:2,j),omega,time(j))
       endif

       if (.not. smalldump) then
          if (abs(omega).lt.tiny(omega)) omega = renvironment('SSPLASH_OMEGA')
          if (abs(omega).gt.tiny(omega) .and. ivx.gt.0) then
             if (.not.all(required(1:2)) .or. .not.all(required(ivx:ivx+1))) then
                print*,' ERROR subtracting corotating frame with partial data read'
             else
                call reset_corotating_velocities(n1,dat(1:n1,1:2,j),dat(1:n1,ivx:ivx+1,j),omega)
             endif
          endif
       endif
    endif

    !--set flag to indicate that only part of this file has been read
    if (.not.all(required(1:ncolstep))) ipartialread = .true.


    nptmassi = 0
    nunknown = 0
    ngas = 0
    nstar = 0
    !--can only do this loop if we have read the iphase array
    iphasealloc: if (allocated(iphase)) then
!
!--translate iphase into particle types (mixed type storage)
!
    if (size(iamtype(:,j)).gt.1) then
       if (phantomdump) then
       !
       !--phantom: translate iphase to splash types
       !
          do i=1,npart
             itype = itypemap_phantom(iphase(i))
             iamtype(i,j) = itype
             select case(itype)
             case(1,2,4) ! remove accreted particles
                if (ih.gt.0 .and. required(ih)) then
                   if (dat(i,ih,j) <= 0.) then
                      iamtype(i,j) = itypemap_unknown_phantom
                   endif
                endif
             case(itypemap_unknown_phantom)
                nunknown = nunknown + 1
             end select
          enddo
       else
       !
       !--sphNG: translate iphase to splash types
       !
          do i=1,npart
             itype = itypemap_sphNG(iphase(i))
             iamtype(i,j) = itype
             select case(itype)
             case(1)
               ngas = ngas + 1
             case(3)
               nptmassi = nptmassi + 1
             case(4)
               nstar = nstar + 1
             case default
               nunknown = nunknown + 1
             end select
          enddo
          do i=npart+1,ntotal
             iamtype(i,j) = 2
          enddo
       endif
       !print*,'mixed types: ngas = ',ngas,nptmassi,nunknown

    elseif (any(iphase(1:ntotal).ne.0)) then
       if (phantomdump) then
          print*,'ERROR: low memory mode will not work correctly with phantom + multiple types'
          print*,'press any key to ignore this and continue anyway (at your own risk...)'
          read*
       endif
!
!--place point masses after normal particles
!  if not storing the iamtype array
!
       print "(a)",' sorting particles by type...'
       nunknown = 0
       do i=1,npart
          if (iphase(i).ne.0) nunknown = nunknown + 1
       enddo
       ncolcopy = min(ncolstep,maxcol)
       allocate(dattemp2(nunknown,ncolcopy))

       iphaseminthistype = 0  ! to avoid compiler warnings
       iphasemaxthistype = 0
       do itype=1,3
          nthistype = 0
          ipos = 0
          select case(itype)
          case(1) ! ptmass
             iphaseminthistype = 1
             iphasemaxthistype = 9
          case(2) ! star
             iphaseminthistype = 10
             iphasemaxthistype = huge(iphasemaxthistype)
          case(3) ! unknown
             iphaseminthistype = -huge(iphaseminthistype)
             iphasemaxthistype = -1
          end select

          do i=1,ntotal
             ipos = ipos + 1
             if (iphase(i).ge.iphaseminthistype .and. iphase(i).le.iphasemaxthistype) then
                nthistype = nthistype + 1
                !--save point mass information in temporary array
                if (nptmassi.gt.size(dattemp2(:,1))) stop 'error: ptmass array bounds exceeded in data read'
                dattemp2(nthistype,1:ncolcopy) = dat(i,1:ncolcopy,j)
   !             print*,i,' removed', dat(i,1:3,j)
                ipos = ipos - 1
             endif
            !--shuffle dat array
             if (ipos.ne.i .and. i.lt.ntotal) then
     !           print*,'copying ',i+1,'->',ipos+1
                dat(ipos+1,1:ncolcopy,j) = dat(i+1,1:ncolcopy,j)
                !--must also shuffle iphase (to be correct for other types)
                iphase(ipos+1) = iphase(i+1)
             endif
          enddo

          !--append this type to end of dat array
          do i=1,nthistype
             ipos = ipos + 1
   !          print*,ipos,' appended', dattemp2(i,1:3)
             dat(ipos,1:ncolcopy,j) = dattemp2(i,1:ncolcopy)
             !--we make iphase = 1 for point masses (could save iphase and copy across but no reason to)
             iphase(ipos) = iphaseminthistype
          enddo

          select case(itype)
          case(1)
             nptmassi = nthistype
             if (nptmassi.ne.nptmass) print *,'WARNING: nptmass from iphase =',nptmassi,'not equal to nptmass =',nptmass
          case(2)
             nstar = nthistype
          case(3)
             nunknown = nthistype
          end select
       enddo

     endif

     endif iphasealloc

     if (allocated(dattemp)) deallocate(dattemp)
     if (allocated(dattempsingle)) deallocate(dattempsingle)
     if (allocated(dattemp2)) deallocate(dattemp2)
     if (allocated(iphase)) deallocate(iphase)
     if (allocated(listpm)) deallocate(listpm)

     call set_labels
     if (.not.phantomdump) then
        npartoftype(:,j) = 0
        npartoftype(1,j) = npart - nptmassi - nstar - nunknown
        npartoftype(2,j) = ntotal - npart
        npartoftype(3,j) = nptmassi
        npartoftype(4,j) = nstar
        npartoftype(5,j) = nunknown
     else
        npartoftype(1,j) = npartoftype(1,j) - nunknown
        npartoftype(itypemap_unknown_phantom,j) = npartoftype(itypemap_unknown_phantom,j) + nunknown
     endif

     call print_types(npartoftype(:,j),labeltype)

     close(15)
     if (debug) print*,' finished data read, npart = ',npart, ntotal, npartoftype(1:ntypes,j)

     return

contains

!
!--reset centre of mass to zero
!
 subroutine reset_centre_of_mass(xyz,pmass,iphase,np)
  implicit none
  integer, intent(in) :: np
  real, dimension(np,3), intent(inout) :: xyz
  real, dimension(np), intent(in) :: pmass
  integer(kind=int1), dimension(np), intent(in) :: iphase
  real :: masstot,pmassi
  real, dimension(3) :: xcm
  integer :: i

  !
  !--get centre of mass
  !
  xcm(:) = 0.
  masstot = 0.
  do i=1,np
     if (iphase(i).ge.0) then
        pmassi = pmass(i)
        masstot = masstot + pmass(i)
        where (required(1:3)) xcm(:) = xcm(:) + pmassi*xyz(i,:)
     endif
  enddo
  xcm(:) = xcm(:)/masstot
  print*,'RESETTING CENTRE OF MASS (',pack(xcm,required(1:3)),') TO ZERO '

  if (required(1)) xyz(1:np,1) = xyz(1:np,1) - xcm(1)
  if (required(2)) xyz(1:np,2) = xyz(1:np,2) - xcm(2)
  if (required(3)) xyz(1:np,3) = xyz(1:np,3) - xcm(3)

  return
 end subroutine reset_centre_of_mass

 subroutine reset_corotating_velocities(np,xy,velxy,omeg)
  implicit none
  integer, intent(in) :: np
  real, dimension(np,2), intent(in) :: xy
  real, dimension(np,2), intent(inout) :: velxy
  real, intent(in) :: omeg
  integer :: ip

  print*,'SUBTRACTING COROTATING VELOCITIES, OMEGA = ',omeg
  do ip=1,np
     velxy(ip,1) = velxy(ip,1) + xy(ip,2)*omeg
  enddo
  do ip=1,np
     velxy(ip,2) = velxy(ip,2) - xy(ip,1)*omeg
  enddo

  return
 end subroutine reset_corotating_velocities

 subroutine reset_corotating_positions(np,xy,omeg,t)
  implicit none
  integer, intent(in) :: np
  real, dimension(np,2), intent(inout) :: xy
  real, intent(in) :: omeg,t
  real :: phii,phinew,r
  integer :: ip

  print*,'SUBTRACTING COROTATING POSITIONS, OMEGA = ',omeg,' t = ',t
!$omp parallel default(none) &
!$omp shared(xy,np) &
!$omp firstprivate(omeg,t) &
!$omp private(ip,r,phii,phinew)
!$omp do
  do ip=1,np
     r = sqrt(xy(ip,1)**2 + xy(ip,2)**2)
     phii = atan2(xy(ip,2),xy(ip,1))
     phinew = phii + omeg*t
     xy(ip,1) = r*COS(phinew)
     xy(ip,2) = r*SIN(phinew)
  enddo
!$omp end do
!$omp end parallel

  return
 end subroutine reset_corotating_positions

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels, only:label,unitslabel,labelzintegration,labeltype,labelvec,iamvec, &
              ix,ipmass,irho,ih,iutherm,ivx,iBfirst,idivB,iJfirst,icv,iradenergy,idustfrac
  use params
  use settings_data,   only:ndim,ndimV,ntypes,ncolumns,UseTypeInRenderings,debugmode
  use geometry,        only:labelcoord
  use settings_units,  only:units,unitzintegration
  use sphNGread
  use asciiutils,      only:lcase
  use system_commands, only:get_environment
  use system_utils,    only:lenvironment
  implicit none
  integer :: i
  real(doub_prec)   :: uergg
  character(len=20) :: string

  if (ndim.le.0 .or. ndim.gt.3) then
     print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
     return
  endif
  if (ndimV.le.0 .or. ndimV.gt.3) then
     print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
     return
  endif
!--all formats read the following columns
  do i=1,ndim
     ix(i) = i
  enddo
  if (igotmass) then
     ipmass = 4   !  particle mass
     ih = 5       !  smoothing length
  else
     ipmass = 0
     ih = 4       !  smoothing length
  endif
  irho = ih + 1     !  density
  if (smalldump .and. nhydroreal4.ge.3) iutherm = irho+1

  !
  !--translate array tags into column labels, where necessary
  !
  if (tagged) then
     do i=1,ncolumns
        label(i) = tagarr(i)
        select case(trim(tagarr(i)))
        case('m')
           ipmass = i
        case('h')
           ih = i
        case('rho')
           irho = i
        case('vx')
           ivx = i
        case('u')
           iutherm = i
        case('divv')
           idivvcol = i
        case('curlvx')
           icurlvxcol = i
        case('curlvy')
           icurlvycol = i
        case('curlvz')
           icurlvzcol = i
        case('Bx')
           iBfirst = i
        case('divB')
           idivB = i
        case('curlBx')
           iJfirst = i
        case('psi')
           label(i) = '\psi'
        case('dustfrac')
           idustfrac = i
        case('alpha')
           label(i) = '\alpha'
        case('alphaB')
           label(i) = '\alpha_B'
        case('EulerAlpha')
           label(i) = 'Euler \alpha'
        case('EulerBeta')
           label(i) = 'Euler \beta'
        case('EtaReal')
           label(i) = '\eta_{real}'
        case('EtaArtificial')
           label(i) = '\eta_{art}'
        case('Erad')
           iradenergy = i
           label(i) = 'radiation energy'
           units(iradenergy) = (udist/utime)**2
        case('opacity')
           label(i) = 'opacity'
           units(i) = udist**2/umass
        case('EddingtonFactor')
           label(i) = 'Eddington Factor'
        case('Cv')
           label(i) = 'u/T'
           icv = i
           units(icv) = (udist/utime)**2
        case('h2ratio')
           label(i) = 'H_2 ratio'
        case('abH1q','abHIq')
           label(i) = 'HI abundance'
        case('abhpq')
           label(i) = 'proton abundance'
        case('abeq')
           label(i) = 'e^- abundance'
        case('abco')
           label(i) = 'CO abundance'
        case default
           if (debugmode) print "(a,i2)",' DEBUG: Unknown label '''//trim(tagarr(i))//''' in column ',i
           label(i) = tagarr(i)
        end select
     enddo
  else
     call guess_labels(ncolumns,iamvec,label,labelvec,istartmhd,istart_extra_real4,&
          nmhd,nhydroreal4,ndimV,irho,iBfirst,ivx,iutherm,idivB,iJfirst,&
          iradenergy,icv,udist,utime,units,unitslabel)
  endif

  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  if (irho.gt.0) label(irho) = 'density'
  if (iutherm.gt.0) label(iutherm) = 'u'
  if (ih.gt.0) label(ih) = 'h       '
  if (ipmass.gt.0) label(ipmass) = 'particle mass'
  if (idivB.gt.0) label(idivB) = 'div B'
  if (idivvcol.gt.0) label(idivvcol) = 'div v'
  if (icurlvxcol.gt.0) label(icurlvxcol) = 'curl v_x'
  if (icurlvycol.gt.0) label(icurlvycol) = 'curl v_y'
  if (icurlvzcol.gt.0) label(icurlvzcol) = 'curl v_z'
  if (icurlvxcol.gt.0 .and. icurlvycol.gt.0 .and. icurlvzcol.gt.0) then
     iamvec(icurlvxcol:icurlvzcol) = icurlvxcol
     labelvec(icurlvxcol:icurlvzcol) = 'curl v'
  endif
  !
  !--set labels for vector quantities
  !
  if (ivx.gt.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
        label(ivx+i-1) = trim(labelvec(ivx))//'_'//labelcoord(i,1)
     enddo
  endif
  if (iBfirst.gt.0) then
     iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
     labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
     do i=1,ndimV
        label(iBfirst+i-1) = trim(labelvec(iBfirst))//'_'//labelcoord(i,1)
     enddo
  endif
  if (iJfirst.gt.0) then
     iamvec(iJfirst:iJfirst+ndimV-1) = iJfirst
     labelvec(iJfirst:iJfirst+ndimV-1) = 'J'
     do i=1,ndimV
        label(iJfirst+i-1) = trim(labelvec(iJfirst))//'_'//labelcoord(i,1)
     enddo
  endif
  !
  !--set units for plot data
  !
!   npower = int(log10(udist))
!   udist = udist/10.**npower
!   udistAU = udist/1.495979e13
   if (ndim.ge.3) then
      units(1:3) = udist
      unitslabel(1:3) = ' [cm]'
   endif
!   do i=1,3
!      write(unitslabel(i),"('[ 10\u',i2,'\d cm]')") npower
!   enddo
   if (ipmass.gt.0) then
      units(ipmass) = umass
      unitslabel(ipmass) = ' [g]'
   endif
   units(ih) = udist
   unitslabel(ih) = ' [cm]'
   if (ivx.gt.0) then
      units(ivx:ivx+ndimV-1) = udist/utime
      unitslabel(ivx:ivx+ndimV-1) = ' [cm/s]'
   endif
   if (iutherm.gt.0) then
      units(iutherm) = (udist/utime)**2
      unitslabel(iutherm) = ' [erg/g]'
   endif
   units(irho) = umass/udist**3
   unitslabel(irho) = ' [g/cm\u3\d]'
   if (iBfirst.gt.0) then
      units(iBfirst:iBfirst+ndimV-1) = umagfd
      unitslabel(iBfirst:iBfirst+ndimV-1) = ' [G]'
   endif

   !--use the following two lines for time in years
   call get_environment('SSPLASH_TIMEUNITS',string)
   select case(trim(lcase(adjustl(string))))
   case('s','seconds')
      units(0) = utime
      unitslabel(0) = trim(string)
   case('min','minutes','mins')
      units(0) = utime/60.d0
      unitslabel(0) = trim(string)
   case('h','hr','hrs','hours','hour')
      units(0) = utime/3600.d0
      unitslabel(0) = trim(string)
   case('y','yr','yrs','years','year')
      units(0) = utime/3.1536d7
      unitslabel(0) = trim(string)
   case('d','day','days')
      units(0) = utime/(3600.d0*24.d0)
      unitslabel(0) = trim(string)
   case('tff','freefall','tfreefall')
   !--or use these two lines for time in free-fall times
      units(0) = 1./tfreefall
      unitslabel(0) = ' '
   case default
      units(0) = utime/3.1536d7
      unitslabel(0) = ' yrs'
   end select
   !--or use these two lines for time in free-fall times
   !units(0) = 1./tfreefall
   !unitslabel(0) = ' '

  unitzintegration = udist
  labelzintegration = ' [cm]'
  !
  !--set labels for each particle type
  !
  if (phantomdump) then  ! phantom
     ntypes = itypemap_unknown_phantom
     labeltype(1) = 'gas'
     labeltype(2) = 'dust'
     labeltype(3) = 'sink'
     labeltype(4) = 'ghost'
     labeltype(5) = 'star'
     labeltype(6) = 'dark matter'
     labeltype(7) = 'bulge'
     labeltype(8) = 'unknown/dead'
     UseTypeInRenderings(:) = .true.
     UseTypeInRenderings(3) = .false.
     if (lenvironment('SSPLASH_PLOT_DUST')) then
        UseTypeInRenderings(2) = .false.
     endif
     if (lenvironment('SSPLASH_PLOT_STARS')) then
        UseTypeInRenderings(5) = .false.
     endif
     if (lenvironment('SSPLASH_PLOT_DM')) then
        UseTypeInRenderings(6) = .false.
     endif
  else
     ntypes = 5
     labeltype(1) = 'gas'
     labeltype(2) = 'ghost'
     labeltype(3) = 'sink'
     labeltype(4) = 'star'
     labeltype(5) = 'unknown/dead'
     UseTypeInRenderings(1) = .true.
     UseTypeInRenderings(2) = .true.
     UseTypeInRenderings(3) = .false.
     UseTypeInRenderings(4) = .true.
     UseTypeInRenderings(5) = .true.  ! only applies if turned on
  endif

!-----------------------------------------------------------

  return
end subroutine set_labels
