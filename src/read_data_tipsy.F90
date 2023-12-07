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
!  Copyright (C) 2005-2023 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING TIPSY FILES
!
! => HANDLES BOTH BINARY AND ASCII TIPSY FORMATS
!    (DETECTS WHICH ONE AUTOMATICALLY)
!
!  BINARY FORMAT READING REQUIRES F2003 STREAM I/O
!  WHICH MAY NOT BE IMPLEMENTED ON OLDER COMPILERS
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
!-------------------------------------------------------------------------

module readdata_tipsy
 use settings_data, only:debugmode
 use byteswap,      only:bs
 implicit none

 public :: read_data_tipsy, set_labels_tipsy, file_format_is_tipsy

 private

contains

subroutine read_data_tipsy(rootname,indexstart,ipos,nstepsread)
 use particle_data,  only:dat,time,npartoftype,gamma,maxpart
 use params
 use settings_data,  only:ndim,ndimV,ncolumns
 use mem_allocation, only:alloc
 use labels,         only:label,ih,ipmass,irho,ivx
 integer, intent(in) :: indexstart,ipos
 integer, intent(out) :: nstepsread
 character(len=*), intent(in) :: rootname
 integer, parameter :: iunit = 16
 integer :: j,ierr
 integer :: nprint,ngas,ndark,nptmass,npart_max,nstep_max
 integer :: ncol,nread,iambinaryfile
 logical :: iexist,do_byteswap
 character(len=len(rootname)) :: dumpfile
 character(len=11) :: fmt
 real :: timei, hfact

 nstepsread = 0
 nstep_max = 0
 npart_max = maxpart
 dumpfile = trim(rootname)
 !
 !--check if first data file exists
 !
 inquire(file=dumpfile,exist=iexist)
 if (.not.iexist) then
    print "(a)",' *** error: '//trim(dumpfile)//': file not found ***'
    return
 endif

 nstep_max = max(nstep_max,indexstart,1)
 j = indexstart
 nstepsread = 0

 write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
 !
 !--determine whether file is binary or ascii and open it
 !
 inquire(file=dumpfile,form=fmt)

 select case(trim(adjustl(fmt)))
 case('UNFORMATTED')
    iambinaryfile = 1
    open(unit=iunit,file=dumpfile,status='old',form='unformatted',access='stream',iostat=ierr)
 case('FORMATTED')
    iambinaryfile = 0
    open(unit=iunit,file=dumpfile,status='old',form='formatted',iostat=ierr)
 case default
    !--if compiler cannot distinguish the two, try ascii first, then binary
    iambinaryfile = -1
    open(unit=iunit,file=dumpfile,status='old',form='formatted',iostat=ierr)
 end select

 if (ierr /= 0) then
    print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
    return
 endif
 !
 !--read the file header
 !  try ascii format first, and if unsuccessful try binary
 !
 if (iambinaryfile==1) then
    write(*,"(a)",advance='no') ' reading binary tipsy format:'
    call read_tipsyheader_binary(iunit,do_byteswap,ierr)
 else
    if (iambinaryfile==0) print "(a)",' reading ascii tipsy format '
    call read_tipsyheader_ascii(iunit,ierr,iambinaryfile)
    if (iambinaryfile < 0) then
       if (ierr==0) then
          !--if successful ascii header read, file is ascii
          iambinaryfile = 0
          print "(a)",' reading ascii tipsy format '
       else
          !--otherwise, close ascii file, and assume file is binary
          close(unit=iunit)
          iambinaryfile = 1
          open(unit=iunit,file=dumpfile,status='old',form='unformatted',&
               access='stream',iostat=ierr)
           write(*,"(a)",advance='no') ' reading binary tipsy format:'
          call read_tipsyheader_binary(iunit,do_byteswap,ierr)
       endif
    endif
 endif
 if (ierr /= 0) then
    print*
    ndim = 0
    ncolumns = 0
    close(unit=iunit)
    return
 endif

 print "(a,f10.2,1(a,i1))",' time: ',timei,' ndim: ',ndim
 print "(4(a,i10))",' ntot: ',nprint,' ngas: ',ngas,' ndark: ',ndark,' nstar: ',nptmass

 ndimV = ndim
 ncol = 2*ndim + 4
 ncolumns = ncol
 !
 !--allocate memory
 !
 if (.not.allocated(dat) .or. nprint > npart_max) then
    npart_max = max(npart_max,nprint)
    call alloc(npart_max,nstep_max,ncolumns)
 endif
 !
 !--now read the timestep data in the dumpfile
 !
 dat(:,:,j) = 0.
 time(j) = timei

 nread = 0
 call set_labels_tipsy

 if (iambinaryfile==1) then
    call read_tipsybody_binary(iunit,do_byteswap,ierr,nread)
 else
    call read_tipsybody_ascii(iunit,ierr,nread)
 endif
 close(unit=iunit)

 if (nread < ncol) then
    print "(a,i2)",' WARNING: END OF FILE: READ TO COLUMN ',nread
    ncolumns = nread
 endif
 !
 !--often tipsy dumps contain only a (fixed) gravitational softening length
 ! for sph particles. In this case we need to create a sensible smoothing length
 ! (and warn about the evils of using fixed softening lengths for sph particles)
 !
 if (ngas >= 0 .and. nread >= irho .and. all(abs(dat(1:ngas,ih,j)-dat(1,ih,j)) <= tiny(dat))) then
    hfact=1.2
    print "(a)",' WARNING: fixed softening lengths detected: simulation may contain artificial fragmentation!'
    print "(a,f5.2,a,i1,a)",'        : creating SPH smoothing lengths using h = ',hfact,'*(m/rho)**(1/',ndim,')'
    dat(1:ngas,ih,j) = hfact*(dat(1:ngas,ipmass,j)/(dat(1:ngas,irho,j) + tiny(dat)))**(1./ndim)
 endif

 nstepsread = nstepsread + 1
 npartoftype(1,j) = ngas
 npartoftype(2,j) = ndark
 npartoftype(3,j) = nptmass
 gamma(j) = 1.666666666667
 j = j + 1

 if (allocated(npartoftype)) then
    print*,'>> end of dump file: nsteps =',j-1,'ntot = ',sum(npartoftype(:,j-1))
 endif

 return

contains

!----------------------------------------------------
! ascii header read
!----------------------------------------------------
subroutine read_tipsyheader_ascii(iunit,ierr,iwarn)
 integer, intent(in) :: iunit,iwarn
 integer, intent(out) :: ierr

 read(iunit,*,end=55,iostat=ierr) nprint,ngas,nptmass
 read(iunit,*,end=55,iostat=ierr) ndim
 read(iunit,*,end=55,iostat=ierr) timei
 ndark = nprint - ngas - nptmass
 !--errors in header read
 if (nprint <= 0 .or. nprint > 1e10 .or. ndim <= 0 .or. ndim > 3 .or. ndark < 0) then
    if (iwarn >= 0) print "(a)",' ERROR reading ascii file header '
    ierr = 2
    return
 endif

 return

55 continue
 if (iwarn >= 0) print "(a)",' ERROR: end of file in ascii header read '
 ierr = -1
 return

end subroutine read_tipsyheader_ascii

!----------------------------------------------------
! binary header read
!----------------------------------------------------
subroutine read_tipsyheader_binary(iunitb,do_byteswap,ierr)
 integer, intent(in)  :: iunitb
 logical, intent(out) :: do_byteswap
 integer, intent(out) :: ierr
 real(doub_prec) :: timedb
 integer :: ipad

 ierr = 0
 do_byteswap = .false.
 read(iunitb,iostat=ierr,end=55) timedb,nprint,ndim,ngas,ndark,nptmass,ipad
 if (debugmode) print*,'header = ',timedb,nprint,ndim,ngas,ndark,nptmass

 !--check for wrong endianness and byte-swap if necessary
 if (ierr /= 0 .or. timedb < 0. .or. bad_header(ndim,nprint,ngas,ndark,nptmass)) then
    timedb = bs(timedb); ndim = bs(ndim); nprint = bs(nprint); ngas = bs(ngas)
    ndark = bs(ndark); nptmass = bs(nptmass)
    if (ierr /= 0 .or. timedb < 0. .or. &
        bad_header(ndim,nprint,ngas,ndark,nptmass)) then
       print "(a)",' ERROR reading binary file header'
       ierr = 2
    else
       do_byteswap = .true.
       write(*,"(a)",advance='no') ' big endian: '
    endif
 endif
 timei = real(timedb)

 if (ndim==0) ndim = 3

 return

55 continue
 print "(a)",' ERROR: end of file in binary header read'
 ierr = -1

end subroutine read_tipsyheader_binary

!----------------------------------------------------
! ascii body read
!----------------------------------------------------
subroutine read_tipsybody_ascii(iunit,ierr,nread)
 integer, intent(in) :: iunit
 integer, intent(out) :: ierr, nread
 integer :: i,ic,icol,nerr

 !--pmass,x,y,z,vx,vy,vz
 do ic=1,2*ndim+1
    nerr = 0
    nread = nread + 1
    if (ic==1) then ! pmass
       icol = ndim + 1
    elseif (ic >= 2 .and. ic <= ndim+1) then ! x, y, z
       icol = ic - 1
    else ! everything after
       icol = ic
    endif
    !print "(1x,a)",trim(label(icol))
    nerr = 0
    do i=1,nprint
       read(iunit,*,end=44,iostat=ierr) dat(i,icol,j)
       if (ierr /= 0) nerr = nerr + 1
    enddo
    if (nerr > 0) print *,'*** WARNING: ERRORS READING '//trim(label(icol))//' ON ',nerr,' LINES'
 enddo
 !--h dark matter
 if (ndark > 0) then
    nerr = 0
    do i=ngas+1,ngas+ndark-1
       read(iunit,*,end=44,iostat=ierr) dat(i,ih,j)
       if (ierr /= 0) nerr = nerr + 1
    enddo
    if (nerr > 0) print *,'*** WARNING: ERRORS READING DARK MATTER H ON ',nerr,' LINES'
 endif
 !--h star particles
 if (nptmass > 0) then
    nerr = 0
    do i=ngas+ndark+1,ngas+ndark+nptmass
       read(iunit,*,end=44,iostat=ierr) dat(i,ih,j)
       if (ierr /= 0) nerr = nerr + 1
    enddo
    if (nerr > 0) print *,'*** WARNING: ERRORS READING PTMASS H ON ',nerr,' LINES'
 endif
 !--density, temperature, sph smoothing length
 do icol=2*ndim+2,ncol
    nread = nread + 1
    !print "(1x,a)",trim(label(icol))
    do i=1,ngas
       read(iunit,*,end=44,iostat=ierr) dat(i,icol,j)
       if (ierr /= 0) nerr = nerr + 1
    enddo
    if (nerr > 0) print *,'*** WARNING: ERRORS READING '//trim(label(icol))//' ON ',nerr,' LINES'
 enddo

 ierr = 0
 return

44 continue
 ierr = -1

end subroutine read_tipsybody_ascii

!----------------------------------------------------
! binary body read
!----------------------------------------------------
subroutine read_tipsybody_binary(iunitb,do_byteswap,ierr,nread)
 use settings_data, only:debugmode
 integer, intent(in)  :: iunitb
 logical, intent(in)  :: do_byteswap
 integer, intent(out) :: ierr,nread
 integer :: i,nerr
 real(kind=4) :: pmass,xyz(3),vxyz(3),rho,temp,h,dummy

 !--gas particles
 nerr = 0
 if (debugmode) print*,'DEBUG: reading ',ngas,' gas particles'
 do i=1,ngas
    !--pmass,x,y,z,vx,vy,vz,rho,temp,h
    read(iunitb,end=44,iostat=ierr) pmass,xyz(1:ndim),vxyz(1:ndim),&
                                    rho,temp,h,dummy,dummy
    if (do_byteswap) then
       pmass = bs(pmass); xyz = bs(xyz); vxyz = bs(vxyz); rho = bs(rho)
       temp  = bs(temp); h = bs(h)
    endif
    dat(i,ipmass,j)         = pmass
    dat(i,1:ndim,j)         = xyz(1:ndim)
    dat(i,ivx:ivx+ndim-1,j) = vxyz(1:ndim)
    dat(i,irho,j)           = rho
    dat(i,irho+1,j)         = temp
    dat(i,ih,j)             = h
    if (ierr /= 0) nerr = nerr + 1
 enddo
 nread = ncolumns
 if (nerr > 0) print *,'*** WARNING: ERRORS READING GAS PARTICLES ON ',nerr,' LINES'

 !--dark matter
 if (ndark > 0) then
    nerr = 0
    do i=ngas+1,ngas+ndark
       !--only read as far as velocities, then eps as smoothing length
       read(iunitb,end=44,iostat=ierr) pmass,xyz,vxyz,h,dummy
       if (do_byteswap) then
          pmass = bs(pmass); xyz = bs(xyz); vxyz = bs(vxyz); h = bs(h)
       endif
       dat(i,ipmass,j)         = pmass
       dat(i,1:ndim,j)         = xyz(1:ndim)
       dat(i,ivx:ivx+ndim-1,j) = vxyz(1:ndim)
       dat(i,ih,j)             = h
       if (ierr /= 0) nerr = nerr + 1
    enddo
    if (nerr > 0) print *,'*** WARNING: ERRORS READING DARK MATTER PARTICLES ON ',nerr,' LINES'
 endif

 !--star particles
 if (nptmass > 0) then
    nerr = 0
    do i=ngas+ndark+1,ngas+ndark+nptmass
       !--only read as far as velocities, then eps as smoothing length
       read(iunitb,end=44,iostat=ierr) pmass,xyz(1:ndim),vxyz(1:ndim),&
                                       dummy,dummy,h,dummy
       if (do_byteswap) then
          pmass = bs(pmass); xyz = bs(xyz); vxyz = bs(vxyz); h = bs(h)
       endif
       dat(i,ipmass,j)         = pmass
       dat(i,1:ndim,j)         = xyz(1:ndim)
       dat(i,ivx:ivx+ndim-1,j) = vxyz(1:ndim)
       dat(i,ih,j)             = h
       if (ierr /= 0) nerr = nerr + 1
    enddo
    if (nerr > 0) print *,'*** WARNING: ERRORS READING STAR PARTICLES ON ',nerr,' LINES'
 endif

 ierr = 0
 return

44 continue
 ierr = -1

end subroutine read_tipsybody_binary

end subroutine read_data_tipsy

!--------------------------------------------------------
! function to check if values read from
! the tipsy header are sensible
!--------------------------------------------------------
logical function bad_header(ndim,nprint,ngas,ndark,nptmass)
 integer, intent(in) :: ndim,nprint,ngas,ndark,nptmass

 bad_header = (ndim < 0 .or. ndim > 3 &
    .or. nprint <= 0 .or. ngas < 0 .or. ndark < 0 .or. nptmass < 0 &
    .or. nprint > 1e10 .or. ngas > 1.e10 &
    .or. ndark > 1.e10 .or. nptmass > 1.e8 )

end function bad_header

!-----------------------------------------------------------
!
! check if a file is in phantom/sphNG format
!
!-----------------------------------------------------------
logical function file_format_is_tipsy(filename) result(is_tipsy)
 character(len=*), intent(in) :: filename
 integer :: iunit,nprint,ndim,ngas,ndark,nptmass,ierr
 real(kind=8) :: timedb

 is_tipsy = .false.
 !
 ! open file and read the first line
 !
 open(newunit=iunit,iostat=ierr,file=filename,status='old',form='unformatted',access='stream')
 if (ierr /= 0) return
 !
 ! check the first line to determine tipsy format
 !
 read(iunit,iostat=ierr) timedb,nprint,ndim,ngas,ndark,nptmass
 if (.not.(ierr /= 0 .or. timedb < 0. .or. &
           bad_header(ndim,nprint,ngas,ndark,nptmass))) then
    is_tipsy = .true.
 elseif (.not. ierr /= 0) then
    ! try byte-swapping
    timedb = bs(timedb); ndim = bs(ndim); nprint = bs(nprint); ngas = bs(ngas)
    ndark = bs(ndark); nptmass = bs(nptmass)
    if (.not. (timedb < 0. .or. bad_header(ndim,nprint,ngas,ndark,nptmass))) then
       is_tipsy =  .true.
    endif
 endif
 close(iunit)    ! close the file

end function file_format_is_tipsy

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------
subroutine set_labels_tipsy
 use labels,        only:label,labelvec,labeltype,iamvec,&
                         ix,ivx,ih,irho,ipmass,make_vector_label
 use settings_data, only:ndim,ndimV,ntypes,UseTypeInRenderings
 use geometry,      only:labelcoord
 integer :: i

 if (ndim <= 0 .or. ndim > 3) then
    print*,'*** ERROR: ndim = ',ndim,' in set_labels_tipsy ***'
    return
 endif
 if (ndimV <= 0 .or. ndimV > 3) then
    print*,'*** ERROR: ndimV = ',ndimV,' in set_labels_tipsy ***'
    return
 endif

 do i=1,ndim
    ix(i) = i
 enddo
 ipmass = ndim + 1
 ivx = ndim + 2
 irho = ivx + ndim
 !iutherm = irho + 1
 label(irho+1) = 'temperature'
 ih = irho + 2

 label(ix(1:ndim)) = labelcoord(1:ndim,1)
 label(ih) = 'h'
 !if (iutherm > 0) label(iutherm) = 'temperature'
 label(ipmass) = 'particle mass'
 label(irho) = 'density'

 call make_vector_label('v',ivx,ndimV,iamvec,labelvec,label,labelcoord(:,1))
 !
 !--set labels for each particle type
 !
 ntypes = 3
 labeltype(1) = 'gas'
 labeltype(2) = 'dark matter'
 labeltype(3) = 'star'
 UseTypeInRenderings(1) = .true.
 UseTypeInRenderings(2) = .false.
 UseTypeInRenderings(3) = .false.

end subroutine set_labels_tipsy

end module readdata_tipsy
