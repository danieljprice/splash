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
!  Copyright (C) 2005-2019 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Module implementing "splash to phantom" operation, writing
! a binary dump file suitable for input to the PHANTOM code
!-----------------------------------------------------------------
module write_data_phantom
 use iso_c_binding, only:c_float,c_double
 implicit none
 integer, parameter :: int8 = selected_int_kind(10)
 integer, parameter :: sing_prec = c_float
 integer, parameter :: doub_prec = c_double
 character(len=10), parameter, public :: formatname='phantom'
 integer, parameter :: lentag = 16

 public :: write_sphdata_phantom
 private

contains
!--------------------------------------------------------------------
!+
!  write header tag for tagged format
!+
!--------------------------------------------------------------------
function tag(string)
 character(len=lentag) :: tag
 character(len=*), intent(in) :: string

 tag = adjustl(string)

end function tag

!--------------------------------------------------------------------
!+
!  write header tag for tagged format
!+
!--------------------------------------------------------------------
function convert_tag_to_phantom(string) result(mytag)
 character(len=lentag) :: mytag
 character(len=*), intent(in) :: string

 mytag = adjustl(string)
 select case(trim(mytag))
 case('Potential')
    mytag = 'poten'
 end select

end function convert_tag_to_phantom

!--------------------------------------------------------------------
!+
!  write output data file that is readable by phantom
!+
!--------------------------------------------------------------------
subroutine write_sphdata_phantom(time,gamma,dat,ndim,ntotal,ntypes,npartoftype, &
                                  masstype,ncolumns,udist,umass,utime,umagfd,labeltype,&
                                  label_dat,ix,ih,ivx,iBfirst,ipmass,irho,iutherm,filename,hsoft_sink)

 integer, intent(in)          :: ndim,ntotal,ntypes,ncolumns
 integer, intent(in)          :: npartoftype(:)

 real, intent(in)             :: time !! real8 vs real4 error

 real, intent(in)             :: gamma
 real, intent(in)             :: dat(ntotal,ncolumns)
 real, intent(in)             :: masstype(:)
 real(doub_prec), intent(in)  :: udist,umass,utime,umagfd
 character(len=*), intent(in) :: labeltype(ntypes),label_dat(ncolumns)
 integer,          intent(in) :: ix(3),ivx,ih,iBfirst,ipmass,irho,iutherm
 character(len=*), intent(in) :: filename
 real,             intent(in) :: hsoft_sink

 integer, parameter    :: i_int   = 1, &
                          i_int1  = 2, &
                          i_int2  = 3, &
                          i_int4  = 4, &
                          i_int8  = 5, &
                          i_real  = 6, &
                          i_real4 = 7, &
                          i_real8 = 8
 integer, parameter    :: idump = 83
 character(len=len(filename)+10) :: outfile

 integer, parameter :: intval1=690706,intval2=780806
 integer, parameter :: int1o=690706 !,int2o=780806
 integer, parameter :: idimhead = 22
 integer(kind=int8) :: nparttot,npartoftypetot(ntypes),number8
 integer            :: nums(8),idot
 integer            :: narraylengths,nblocks,nblockarrays,ntypesi
 integer            :: i,j,ierr,i1,index1,number,nptmass,iversion,np,maxrhead
 real               :: rheader(idimhead)
 character(len=lentag) :: rheader_tags(idimhead)
 real               :: r1,hfact,macc,spinx,spiny,spinz
 character(len=2), parameter :: vlabel(3) = (/'vx','vy','vz'/)
 character(len=2), parameter :: Blabel(3) = (/'Bx','By','Bz'/)
 logical            :: mhd
 logical, allocatable :: mask(:)
!
! sink particle locations in dat array
!
 integer, allocatable :: ilocsink(:)
!
!--define output file name
!
 idot = index(filename,'.')-1
 if (idot <= 0) idot = len_trim(filename)
 outfile=filename(1:idot)//'.tmp'

 narraylengths = 2
 nblocks = 1          ! not parallel dump
 hfact   = 1.2        ! must be specified in phantom dumps
 nptmass = 0          ! work this out later
!
!--check if we have enough data to write a PHANTOM dump
!
 allocate(mask(ncolumns))
 mask = .false.
 if (ndim < 3) then
    print "(a)",' ERROR: ndim < 3 but must be 3 for PHANTOM data -- cannot write PHANTOM dump, skipping...'
    return
 endif
 if (any(ix(:) <= 0)) then
    print "(a)",' ERROR: position labels not set -- cannot write PHANTOM dump, skipping...'
    return
 endif
 if (ivx <= 0) then
    print "(a)",' ERROR: velocity not found in data -- cannot write PHANTOM dump, skipping...'
    return
 endif
 if (ih <= 0) then
    print "(a)",' ERROR: smoothing length not found in data -- cannot write PHANTOM dump, skipping...'
    return
 endif
 mhd = .false.
 if (iBfirst > 0) then
    mhd = .true.
    narraylengths = 4
    mask(iBfirst:iBfirst+ndim-1) = .true.
 endif
 ! do not write known quantities twice
 mask(ix) = .true.
 mask(ivx:ivx+ndim-1) = .true.
 mask(ih) = .true.
 if (iutherm > 0) mask(iutherm) = .true.
 if (ipmass > 0) mask(ipmass) = .true.
 if (irho > 0) mask(irho) = .true.

!
!--figure out whether we have sink particles
!
 call extract_sink_particles_from_data(ntypes,npartoftype,labeltype,np,&
      npartoftypetot,nptmass,ntypesi,ilocsink)

 nparttot = np

!--fill rheader and check that we have equal mass particles
 rheader_tags = ' '
 rheader(:) = 0.
 rheader(1) = time
 rheader(2) = gamma
 rheader(3) = hfact
 rheader_tags(1) = 'time'
 rheader_tags(2) = 'gamma'
 rheader_tags(3) = 'hfact'
 if (ipmass > 0) then
    index1 = 1
    do i=1,ntypesi
       rheader(3+i) = dat(index1,ipmass)
       rheader_tags(3+i) = 'massoftype'
       if (npartoftype(i) > 0) then
          if (any(dat(index1:index1+npartoftype(i)-1,ipmass) /= dat(index1,ipmass))) then
             print "(a)",' WARNING: unequal mass particles detected but PHANTOM only accepts equal mass...'
          endif
          index1 = index1 + npartoftype(i) - 1
       endif
    enddo
 else
    do i=1,ntypesi
       rheader(3+i) = masstype(i)
       rheader_tags(3+i) = 'massoftype'
    enddo
 endif
 maxrhead = 3+ntypesi

 write(*,"(/,/,'-------->   TIME = ',f10.4,"// &
              "': full dump written to file ',a,' on unit ',i2,'   <--------',/)") &
       time,trim(outfile),idump

 open(unit=idump,file=outfile,status='replace',form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    write(*,*) 'error: can''t create new dumpfile '//trim(outfile)
    return
 endif
!
!--write full dump Phantom/sphNG file
!
 i1 = intval1
 r1 = real(intval2)
 iversion = 1 ! file version to write
 write (idump, err=100) intval1,r1,intval2,iversion,int1o
 write (idump, err=100) fileident('F','Phantom',mhd=mhd)
!
!--single values
!
!--default int
 number = 4+ntypesi
 write (idump, err=100) number
 write (idump, err=100) tag('nparttot'),tag('ntypes'),(tag('npartoftype'),i=1,ntypesi),tag('nblocks'),tag('nptmass')
 write (idump, err=100) int(nparttot),ntypesi,(int(npartoftypetot(i)),i=1,ntypesi),nblocks,nptmass

!--int*1, int*2, int*4
 number = 0
 do i = 1, 3
    write (idump, err=100) number
 enddo
!--int*8
 number = 2 + ntypesi
 write (idump, err=100) number
 write (idump, err=100) tag('nparttot'),tag('ntypes'),(tag('npartoftype'),i=1,ntypesi)
 write (idump, err=100) nparttot,int(ntypesi,kind=int8),npartoftypetot(1:ntypesi)

!--default real
 write (idump, err=100) maxrhead
 write (idump, err=100) rheader_tags(1:maxrhead)
 write (idump, err=100) rheader(1:maxrhead)

!--real*4
 number = 0
 write (idump, err=100) number
!--real*8
 if (umagfd > 0.) then
    number = 4
    write (idump, err=100) number
    write (idump, err=100) tag('udist'),tag('umass'),tag('utime'),tag('umagfd')
    write (idump, err=100) udist, umass, utime, umagfd
 else
    number = 3
    write (idump, err=100) number
    write (idump, err=100) tag('udist'),tag('umass'),tag('utime')
    write (idump, err=100) udist, umass, utime
 endif

 nblockarrays = narraylengths*nblocks
 write (idump, err=100) nblockarrays
!
!--array length 1 header
!
 number8 = np
 nums(:) = 0
 if (iutherm > 0) then
    nums(i_real) = 7
 else
    nums(i_real) = 6
 endif
 ! write any array not already counted
 nums(i_real) = nums(i_real) + count(.not.mask)

 nums(i_real4) = 1
 write (idump, err=100) number8, (nums(i), i=1,8)
!
!--array length 2 header
!
 number8 = nptmass
 nums(:) = 0
 if (nptmass > 0) nums(i_real) = 13
 write (idump, err=100) number8, (nums(i), i=1,8)
!
!--array length 3 header
!
 if (narraylengths >= 3) then
    number8 = 0
    nums(1:8) = 0
    write (idump, err=100) number8, (nums(i), i=1,8)
 endif
!
!--array length 4 header
!
 if (narraylengths >= 4) then
    if (mhd) then
       number8 = np
    else
       number8 = 0
    endif
    nums(:) = 0
    if (mhd) nums(i_real4) = 3
    write (idump, err=100) number8, (nums(i), i=1,8)
 endif

!
!--array length 1 arrays
!
!--default int
!--int*1
!--int*2
!--int*4
!--int*8
!--default real
 do j = 1, 3
    write (idump, err=100) tag(label_dat(ix(j)))
    write (idump, err=100) (dat(i,ix(j)), i=1, np)
 enddo

 do j = 1, 3
    write (idump, err=100) tag(vlabel(j))
    write (idump, err=100) (dat(i,ivx+j-1), i=1, np)
 enddo

 if (iutherm > 0) then
    write (idump, err=100) tag(label_dat(iutherm))
    write (idump, err=100) (dat(i,iutherm), i=1, np)
 endif

 do j=1,ncolumns
    if (.not.mask(j)) then
       print*,tag(label_dat(j)),'->',convert_tag_to_phantom(label_dat(j))
       write (idump, err=100) convert_tag_to_phantom(label_dat(j))
       write (idump, err=100) (dat(i,j), i=1, np)
    endif
 enddo

!--real*4
!   dump smoothing length as a real*4 to save space
 write (idump, err=100) tag('h')
 write (idump, err=100) (real(dat(i,ih),kind=sing_prec), i=1, np)
!
!--sink particle arrays
!
 if (nptmass > 0 .and. allocated(ilocsink)) then
    do j = 1, 3
       write (idump, err=100) tag(label_dat(ix(j)))
       write (idump, err=100) (dat(ilocsink(i),ix(j)),i=1,nptmass)
    enddo
    write (idump, err=100) tag(label_dat(ipmass))
    write (idump, err=100) (dat(ilocsink(i),ipmass),i=1,nptmass)
    write (idump, err=100) tag('h')
    write (idump, err=100) (dat(ilocsink(i),ih),i=1,nptmass)

    ! extra sink information
    macc = 0.
    spinx = 0.
    spiny = 0.
    spinz = 0.
    write (idump, err=100) tag('hsoft')
    write (idump, err=100) (hsoft_sink,i=1,nptmass)
    write (idump, err=100) tag('maccreted')
    write (idump, err=100) (macc,i=1,nptmass)
    write (idump, err=100) tag('spinx')
    write (idump, err=100) (spinx,i=1,nptmass)
    write (idump, err=100) tag('spiny')
    write (idump, err=100) (spiny,i=1,nptmass)
    write (idump, err=100) tag('spinz')
    write (idump, err=100) (spinz,i=1,nptmass)
    do j = 1, 3
       write (idump, err=100) tag(label_dat(ivx+j-1))
       write (idump, err=100) (dat(ilocsink(i),ivx+j-1),i=1,nptmass)
    enddo
 endif

 if (allocated(ilocsink)) deallocate(ilocsink)

 if (mhd) then
    do j=1,3
       write (idump,err=100) tag(Blabel(j))
       write (idump,err=100) (real(dat(i,iBfirst+j-1),kind=sing_prec),i=1, np)
    enddo
 endif

 close(unit=idump)
 return

100 continue
 write(*,*) 'error whilst writing dumpfile '//trim(outfile)
 close(unit=idump)

end subroutine write_sphdata_phantom

!--------------------------------------------------------------------
!+
!  contruct header string based on compile-time options
!  these are for information only (ie. not important for restarting)
!+
!--------------------------------------------------------------------
character(len=100) function fileident(firstchar,codestring,mhd)
 character(len=1), intent(in) :: firstchar
 character(len=*), intent(in), optional :: codestring
 logical,          intent(in), optional :: mhd
 character(len=10) :: datestring, timestring, string
 logical :: gotmhd
!
!--print date and time stamp in file header
!
 call date_and_time(datestring,timestring)
 datestring = datestring(7:8)//'/'//datestring(5:6)//'/'//datestring(1:4)
 timestring = timestring(1:2)//':'//timestring(3:4)//':'//timestring(5:)

 string = ' '

 if (present(codestring)) then
    fileident = firstchar//'T:'//trim(codestring)
 else
    fileident = firstchar//'T:Phantom'
 endif

 gotmhd = .false.
 if (present(mhd)) gotmhd = mhd
 if (gotmhd) then
    fileident = trim(fileident)//' (mhd'//trim(string)//')  : '//trim(datestring)//' '//trim(timestring)
 else
    fileident = trim(fileident)//' (hydro'//trim(string)//'): ' &
                //trim(datestring)//' '//trim(timestring)
 endif

end function fileident

!--------------------------------------------------------------------
!+
!  extract sink particle information from dat array, as these arrays
!  are written separately in phantom data files
!+
!--------------------------------------------------------------------
subroutine extract_sink_particles_from_data(ntypes,npartoftype,labeltype,np,noftype,nptmass,ntypesi,ilocsink)
 integer,              intent(in)  :: ntypes,npartoftype(ntypes)
 character(len=*),     intent(in)  :: labeltype(ntypes)
 integer,              intent(out) :: np,nptmass,ntypesi
 integer(kind=int8),   intent(out) :: noftype(ntypes)
 integer, allocatable, intent(out) :: ilocsink(:)
 integer :: i,j

 np = 0
 ntypesi = ntypes
 noftype(:) = int(npartoftype(:),kind=8)
 over_types: do i=1,ntypes
    if (trim(labeltype(i))=='sink' .or. (trim(labeltype(i))=='dark matter' .and. npartoftype(i) <= 1000)) then
       nptmass = npartoftype(i)
       noftype(i) = 0
       allocate(ilocsink(nptmass))
       do j=1,nptmass
          ilocsink(j) = np+j
       enddo
       ntypesi = ntypes - 1  ! do not write types after this
       exit over_types
    else
       np = np + npartoftype(i)
    endif
 enddo over_types
 if (nptmass > 0) print "(/,a,i2,a)",' WRITING ',nptmass,' SINK PARTICLES'

end subroutine extract_sink_particles_from_data

end module write_data_phantom
