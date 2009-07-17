!-----------------------------------------------------------------
! Module implementing "splash to phantom" operation, writing
! a binary dump file suitable for input to the PHANTOM code
!-----------------------------------------------------------------
module write_data_phantom
 implicit none
 character(len=10), parameter, public :: formatname='phantom'

 public :: write_sphdata_phantom
 private

contains

subroutine write_sphdata_phantom(time,gamma,dat,ntotal,ntypes,npartoftype, &
                                 masstype,ncolumns,filename)
 use labels,         only:labeltype,label,ih,ivx,iBfirst,ipmass,ix,iutherm
 use settings_units, only:units
 use settings_data,  only:ndim,UseTypeInRenderings
 use params,         only:int1,int8,doub_prec
 implicit none
 integer, intent(in)                          :: ntotal,ntypes,ncolumns
 integer, intent(in), dimension(:)            :: npartoftype
 real, intent(in)                             :: time,gamma
 real, intent(in), dimension(ntotal,ncolumns) :: dat
 real, intent(in), dimension(:)               :: masstype
 character(len=*), intent(in)                 :: filename

 integer, parameter    :: idump = 83
 character(len=len(filename)+10) :: outfile

 integer, parameter    :: intval1=690706,intval2=780806
 integer, parameter    :: idimhead = 22
 integer(kind=int8)    :: nparttot,npartoftypetot(5),number8
 integer, dimension(8) :: nums
 integer               :: narraylengths,nblocks,nblockarrays
 integer               :: i,j,ierr,i1,index1,number,npart
 real, dimension(idimhead) :: rheader
 real(doub_prec)           :: udist,umass,utime,umagfd
 real                      :: r1,hfact
!
!--statement function defining conversion to a 4 byte real
!
 real*4 :: real4
 real   :: dum
 real4(dum) = dum
!
!--define output file name
!
 outfile=trim(filename)//'.init'
 narraylengths = 2
 nblocks = 1          ! not parallel dump
 hfact   = 1.2        ! must be specified in phantom dumps

!
!--check if we have enough data to write a PHANTOM dump
!
 if (ndim.lt.3) then
    print "(a)",' ERROR: ndim < 3 but must be 3 for PHANTOM data -- cannot write PHANTOM dump, skipping...'
    return
 endif
 if (any(ix(:).le.0)) then
    print "(a)",' ERROR: position labels not set -- cannot write PHANTOM dump, skipping...'
    return
 endif
 if (ivx.le.0) then
    print "(a)",' ERROR: velocity not found in data -- cannot write PHANTOM dump, skipping...'
    return
 endif
 if (ih.le.0) then
    print "(a)",' ERROR: smoothing length not found in data -- cannot write PHANTOM dump, skipping...'
    return
 endif
!--fill rheader and check that we have equal mass particles
 rheader(:) = 0.
 rheader(1) = time
 rheader(3) = gamma
 rheader(6) = hfact
 if (ipmass.gt.0) then
    index1 = 1
    do i=1,ntypes
       rheader(14+i) = dat(index1,ipmass)
       if (any(dat(index1:index1+npartoftype(i),ipmass).ne.dat(index1,ipmass))) then
          print*,' ERROR: unequal mass particles detected but PHANTOM only accepts equal mass, skipping...'
          return
       endif
       index1 = index1 + npartoftype(i)
    enddo
 else
    do i=1,ntypes
       rheader(14+i) = masstype(i)
    enddo
 endif

 write(*,"(/,/,'-------->   TIME = ',f10.4,"// &
              "': full dump written to file ',a,' on unit ',i2,'   <--------',/)") &
       time,trim(outfile),idump

 open(unit=idump,file=outfile,status='new',form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    write(*,*) 'error: can''t create new dumpfile '//trim(outfile)
    return
 endif
!
!--write full dump Phantom/sphNG file
!
 i1 = intval1
 r1 = real(intval2)
 write (idump, err=100) intval1,r1,intval2,i1,intval1
 write (idump, err=100) fileident('F','Phantom')
 
 npart = npartoftype(1)
 npartoftypetot(:) = 0
 do i=2,ntypes
    if (all(UseTypeInRenderings(1:i))) then
       npart = npart + npartoftype(i)
       if (npartoftype(i).gt.0) print "(a)",' WARNING: assuming '// &
          trim(labeltype(i))//' particles are same as gas particles'
       if (rheader(15).le.0.) then
          rheader(15) = masstype(i)
          rheader(15+i) = 0.
       elseif (masstype(i).ne.rheader(15)) then
          print*,' WARNING! WARNING! mass of '//trim(labeltype(i))// &
                ' particles differs from '//trim(labeltype(1))//' particles'
          print*,' Assuming all particles have '//trim(labeltype(1))//' particle mass'
       endif          
    endif
 enddo
 npartoftypetot(1) = npart
 nparttot = npart
!
!--single values
!
!--default int
 number = 7
 write (idump, err=100) number
 write (idump, err=100) int(nparttot),(int(npartoftypetot(i)),i=1,5),nblocks
!--int*1, int*2, int*4
 number = 0
 do i = 1, 3
    write (idump, err=100) number
 end do
!--int*8 
 number = 1 + ntypes
 write (idump, err=100) number
 write (idump, err=100) nparttot,npartoftypetot(1:ntypes)

!--default real

 write (idump, err=100) idimhead
 write (idump, err=100) rheader(1:idimhead)

!--real*4
 number = 0
 write (idump, err=100) number
!--real*8
 udist = units(ix(1))
 utime = units(0)
 if (ipmass.gt.0) then
    umass = units(ipmass)
 else
    print "(a)",' WARNING: units for mass unknown, written as 1.0'
    umass = 1.0d0
 endif
 if (iBfirst.gt.0) then
    umagfd = units(iBfirst)
    number = 4
    write (idump, err=100) number
    write (idump, err=100) udist, umass, utime, umagfd
 else
    number = 3
    write (idump, err=100) number
    write (idump, err=100) udist, umass, utime
 endif

 nblockarrays = narraylengths*nblocks
 write (idump, err=100) nblockarrays
!
!--array length 1 header
!
 number8 = npart
 nums(:) = 0
 if (iutherm.gt.0) then
    nums(6) = 7
 else
    nums(6) = 6
 endif
 nums(7) = 1
 write (idump, err=100) number8, (nums(i), i=1,8)
!
!--array length 2 header
!
 number8 = 0
 nums(:) = 0
 write (idump, err=100) number8, (nums(i), i=1,8)
!
!--array length 3 header
!
 if (narraylengths.ge.3) then
    number8 = 0
    nums(1:8) = 0
    write (idump, err=100) number8, (nums(i), i=1,8)
 endif
!
!--array length 4 header
!
 if (narraylengths.ge.4) then
    !if (imhd.eq.idim) then
    !   number8 = npart
    !else
       number8 = 0
    !endif
    nums(:) = 0
!    if (imhd.eq.idim) then
!       if (ivecp.eq.idim) then
!          nums(6) = 0
!          nums(7) = 3 + iBevol + 1
!       else
!          nums(6) = 0
!          nums(7) = 4
!       endif
!    endif
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
    write (idump, err=100) (dat(i,ix(j)), i=1, npart)
 end do
 
 do j = 1, 3
    write (idump, err=100) (dat(i,ivx+j), i=1, npart)
 end do

 if (iutherm.gt.0) then
    write (idump, err=100) (dat(i,iutherm), i=1, npart)
 endif

!--real*4
!   dump smoothing length as a real*4 to save space
  write (idump, err=100) (real4(dat(i,ih)), i=1, npart)

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
character(len=100) function fileident(firstchar,codestring)
 implicit none
 character(len=1), intent(in) :: firstchar
 character(len=*), intent(in), optional :: codestring
 character(len=10) :: datestring, timestring, string
!
!--print date and time stamp in file header
!
 call date_and_time(datestring,timestring)
 datestring = datestring(7:8)//'/'//datestring(5:6)//'/'//datestring(1:4)
 timestring = timestring(1:2)//':'//timestring(3:4)//':'//timestring(5:)

 string = ' '

 if (present(codestring)) then
    fileident = firstchar//':'//trim(codestring)
 else
    fileident = firstchar//':Phantom'
 endif

 !if (imhd.eq.idim) then
 !   fileident = trim(fileident)//' (mhd'//trim(string)//')  : '//trim(datestring)//' '//trim(timestring)
 !else
    fileident = trim(fileident)//' (hydro'//trim(string)//'): ' &
                //trim(datestring)//' '//trim(timestring)
 !endif

end function fileident

end module write_data_phantom
