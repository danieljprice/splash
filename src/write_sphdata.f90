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

!-----------------------------------------------------------------
!     module containing output routines for writing SPH data
!     as read by SPLASH to output file in various formats
!-----------------------------------------------------------------
module write_sphdata
 implicit none
 public :: issphformat,write_sphdump
 private

contains

!-----------------------------------------------------------------
! utility to check if a format selection is valid
!-----------------------------------------------------------------
logical function issphformat(string)
 use asciiutils, only:lcase
 character(len=*), intent(in) :: string
 logical :: verbose

 issphformat = .false.
 verbose = .true.
 select case(trim(lcase(string)))
 case('ascii','csv')
    issphformat = .true.
 case('binary','binarystream','stream')
    issphformat = .true.
 case('rsph')
    issphformat = .true.
 case('phantom')
    issphformat = .true.
 case('gadget')
    issphformat = .true.
 case('none')
    verbose = .false.
 end select

 if (.not.issphformat .and. verbose) then
    print "(a)",' convert mode ("splash to X dumpfiles"): '
    print "(a,/)",' splash to ascii   : convert SPH data to ascii file dumpfile.ascii'
    print "(a)",  '        to csv     : convert SPH data to csv file dumpfile.csv'
    print "(a)",  '        to binary  : convert SPH data to simple unformatted binary dumpfile.binary: '
    print "(a)",  '                      write(1) time,np,ncols'
    print "(a)",  '                      do i=1,np'
    print "(a)",  '                         write(1) dat(1:ncols),itype'
    print "(a)",  '                      enddo'
    print "(a)",  '        to stream  : convert SPH data to streamed binary format dumpfile.binarystream: '
    if (kind(1.)==8) then
       print "(a)",  '                      t,np,ncols,np*(dat(1:ncols),itype) [8,4,4,np*(ncols*8,4)] bytes'
    else
       print "(a)",  '                      t,np,ncols,np*(dat(1:ncols),itype) [4,4,4,np*(ncols*4,4)] bytes'
    endif
    print "(a)",  '        to phantom : convert SPH data to binary dump file for PHANTOM'
    print "(a)",  '        to gadget  : convert SPH data to default GADGET snapshot file format'
 elseif (.not.issphformat) then
    print "(a)",'Convert mode: '
    print "(a)",'   splash to ascii     : convert to ascii; type "splash to" for other formats'
 endif

 return
end function issphformat

subroutine write_sphdump(time,gamma,dat,npart,ntypes,npartoftype,masstype,itype,ncolumns,filename,outformat,listofcolumns)
 use labels,         only:labeltype,label,unitslabel,irho,ipmass,ix,iBfirst,ih,iutherm,ivx
 use settings_units, only:units
 use settings_data,  only:ndim,icoords,icoordsnew,xorigin
 use params,         only:int1,maxplot,doub_prec
 use write_data_phantom, only:write_sphdata_phantom
 use write_data_gadget,  only:write_sphdata_gadget
 use filenames,      only:tagline
 use geomutils,      only:change_coords
 use asciiutils,     only:lcase
 integer, intent(in)                          :: npart,ntypes,ncolumns
 integer, intent(in), dimension(:)            :: npartoftype
 integer(kind=int1), intent(in), dimension(:) :: itype
 real, intent(in)                             :: time,gamma
 real, intent(in), dimension(npart,ncolumns)  :: dat
 real, intent(in), dimension(:)               :: masstype
 character(len=*), intent(in)                 :: filename,outformat
 integer, intent(in), dimension(:), optional  :: listofcolumns
 integer, parameter :: iunit = 83
 integer, parameter :: maxline = 1000
 integer            :: ierr,i,idim,i1,i2,ncols
 integer, dimension(ncolumns) :: iorder
 character(len=40)  :: fmtstring,fmtstring2,fmtstringlab,outfile
 character(len=6)   :: ext
 real(kind=doub_prec), dimension(maxplot) :: vals
 real(kind=doub_prec) :: udist,umass,utime,umagfd
 real, dimension(3) :: x0,v0
 logical :: change_coordsys,use_csv
 !
 !--allow for the possibility of writing columns in a different order
 !
 if (present(listofcolumns)) then
    call get_order_from_list(listofcolumns,iorder,ncols)
 else
    ncols = ncolumns
    do i=1,ncols
       iorder(i) = i
    enddo
 endif

 select case(trim(lcase(outformat)))
 case ('ascii','csv')
    ext = '.ascii'
    use_csv = (trim(lcase(outformat))=='csv')
    if (use_csv) ext = '.csv'

    print "(/,5('-'),'>',a,i2,a,1x,'<',5('-'),/)",' WRITING TO FILE '//trim(filename)//trim(ext)//' WITH ',ncols,' COLUMNS'

    !--format the header lines to go in the ascii file
    if (kind(1.)==8) then
       if (use_csv) then
          write(fmtstring,"(i10,a)") ncols,"(es23.15,',')"
          write(fmtstringlab,"(i10,a)") ncols,"(a23,','),a"
       else
          write(fmtstring,"(i10,a)") ncols,'(es23.15,1x)'
          write(fmtstringlab,"(i10,a)") ncols,'(a23,1x),a'
       endif
    else
       if (use_csv) then
          write(fmtstring,"(i10,a)") ncols,"(es15.7,',')"
          write(fmtstringlab,"(i10,a)") ncols,"(a15,','),a"
       else
          write(fmtstring,"(i10,a)") ncols,'(es15.7,1x)'
          write(fmtstringlab,"(i10,a)") ncols,'(a15,1x),a'
       endif
    endif
    fmtstring2 = '('//trim(adjustl(fmtstring))//',i0)'
    fmtstring = '('//trim(adjustl(fmtstring))//')'
    if (use_csv) then
       fmtstringlab = '('//trim(adjustl(fmtstringlab))//')'
    else
       fmtstringlab = '(''#'',1x,'//trim(adjustl(fmtstringlab))//')'
    endif

    open(unit=iunit,file=trim(filename)//trim(ext),status='replace',form='formatted',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR OPENING FILE FOR WRITING'
       return
    endif
    if (.not.use_csv) then
       write(iunit,"(a)",iostat=ierr) '# '//trim(filename)//trim(ext)//', created by '//trim(tagline)
       write(iunit,"('#')",iostat=ierr)
       write(iunit,"('#',1x,'time:',13x,'time unit (',a,')')",iostat=ierr) trim(unitslabel(0))
       write(iunit,"('#',2(1x,1pe15.7))",iostat=ierr) time,units(0)
       write(iunit,"('#')",iostat=ierr)
       write(iunit,"('#',1x,'npart:',30(1x,a12))",iostat=ierr) (trim(labeltype(i)),i=1,ntypes)
       write(iunit,"('#',7x,30(1x,i12))",iostat=ierr) npartoftype(1:ntypes)
       write(iunit,"('# units:')",iostat=ierr)
       write(iunit,"('#'"//fmtstring(2:),iostat=ierr) units(iorder(1:ncols))
       write(iunit,fmtstringlab,iostat=ierr) unitslabel(iorder(1:ncols))
       write(iunit,"('#')",iostat=ierr)
    endif
    if (ierr /= 0) print "(a)",' ERROR WRITING ASCII HEADER'
    !
    !--write body
    !
    change_coordsys = (icoordsnew /= icoords .and. ndim > 0 .and. all(ix(1:ndim) > 0))
    x0 = xorigin(:)  ! note that it is not currently possible to do splash to ascii
    v0 = 0.          ! with coords set relative to a tracked particle, so just use xorigin

    if (size(itype) > 1) then
       write(iunit,fmtstringlab,iostat=ierr) label(iorder(1:ncols)),'itype'
       do i=1,npart
          vals(1:ncolumns) = dat(i,1:ncolumns)
          if (change_coordsys) call change_coords(vals,ncolumns,ndim,icoords,icoordsnew,x0,v0)
          write(iunit,fmtstring2,iostat=ierr) vals(iorder(1:ncols)),itype(i)
       enddo
    else
       write(iunit,fmtstringlab,iostat=ierr) label(iorder(1:ncols))
       if (change_coordsys) then
          do i=1,npart
             vals(1:ncolumns) = dat(i,1:ncolumns)
             call change_coords(vals,ncolumns,ndim,icoords,icoordsnew,x0,v0)
             write(iunit,fmtstring,iostat=ierr) vals(iorder(1:ncols))
          enddo
       else
          do i=1,npart
             write(iunit,fmtstring,iostat=ierr) dat(i,iorder(1:ncols))
          enddo
       endif
    endif
    close(iunit)
    if (ierr /= 0) then
       print "(a)",' ERROR WRITING ASCII FILE'
    endif
    return

 case ('binary','binarystream','stream')
!
!--This is the most basic binary (ie. unformatted) file format I could think of,
!  as an alternative to ascii for large files.
!
    select case(trim(outformat))
    case('binarystream','BINARYSTREAM','stream')
       print "(/,5('-'),'>',a,i2,a,1x,'<',5('-'),/)",' WRITING TO FILE '//trim(filename)//'.binarystream WITH ',ncolumns,' COLUMNS'
       open(unit=iunit,file=trim(filename)//'.binarystream',status='replace',form='unformatted',access='stream',iostat=ierr)
    case default
       print "(/,5('-'),'>',a,i2,a,1x,'<',5('-'),/)",' WRITING TO FILE '//trim(filename)//'.binary WITH ',ncolumns,' COLUMNS'
       open(unit=iunit,file=trim(filename)//'.binary',status='replace',form='unformatted',iostat=ierr)
    end select
    if (ierr /= 0) then
       print "(a)",' ERROR OPENING FILE FOR WRITING'
       return
    endif
    write(iunit,iostat=ierr) time,npart,ncols
    if (ierr /= 0) print "(a)",' ERROR WRITING HEADER LINE TO BINARY FILE '
    !
    !--write body
    !
    if (size(itype) > 1) then
       do i=1,npart
          write(iunit,iostat=ierr) dat(i,iorder(1:ncols)),int(itype(i))
       enddo
    else
       do i=1,npart
          write(iunit,iostat=ierr) dat(i,iorder(1:ncols))
       enddo
    endif
    close(iunit)
    if (ierr /= 0) print "(a)",' ERROR WRITING BINARY FILE'
    return

 case ('rsph')
!
!--Files for Steinar Borve's RSPH format
!
    if (ndim < 2) then
       print "(a)",' ERROR: cannot write RSPH format for < 2D'
       return
    endif
    outfile = 'rsph2D_pos.dat'
    print "(a)",' writing to '//trim(outfile)
    open(unit=iunit,file=outfile,status='replace',form='formatted',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR OPENING '//trim(outfile)//' FOR WRITING'
       return
    endif
    write(iunit,"(i1)") ndim
    write(iunit,"(a)") 'position'
    write(iunit,"(i4)") maxline
    write(iunit,*) (minval(dat(1:npart,ix(i))),i=1,ndim)
    write(iunit,*) (maxval(dat(1:npart,ix(i))),i=1,ndim)
    write(iunit,*) time
    write(iunit,*) npart
    do idim=1,2
       i1 = 1
       i2 = 0
       do while (i2 < npart)
          i2 = min(i2 + maxline,npart)
          write(iunit,*) dat(i1:i2,ix(idim))
          i1 = i2 + 1
       enddo
    enddo
    close(unit=iunit)

    outfile = 'rsph2D_rho.dat'
    print "(a)",' writing to '//trim(outfile)
    open(unit=iunit,file=outfile,status='replace',form='formatted',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR OPENING FILE FOR WRITING'
       return
    endif
    write(iunit,"(i1)") ndim
    write(iunit,"(a)") 'density'
    write(iunit,"(i4)") maxline
    write(iunit,*) (minval(dat(1:npart,ix(i))),i=1,ndim)
    write(iunit,*) (maxval(dat(1:npart,ix(i))),i=1,ndim)
    i1 = 1
    i2 = 0
    do while (i2 < npart)
       i2 =  min(i2 + maxline,npart)
       write(iunit,*) dat(i1:i2,irho)
       i1 = i2 + 1
    enddo
    close(unit=iunit)

    outfile = 'rsph2D_siz.dat'
    print "(a)",' writing to '//trim(outfile)
    open(unit=iunit,file=outfile,status='replace',form='formatted',iostat=ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR OPENING FILE FOR WRITING'
       return
    endif
    write(iunit,"(i1)") ndim
    write(iunit,"(a)") 'size'
    write(iunit,"(i4)") maxline
    write(iunit,*) (minval(dat(1:npart,ix(i))),i=1,ndim)
    write(iunit,*) (maxval(dat(1:npart,ix(i))),i=1,ndim)
    i1 = 1
    i2 = 0
    do while (i2 < npart)
       i2 =  min(i2 + maxline,npart)
       write(iunit,*) ((dat(i,ipmass)/dat(i,irho))**(1./ndim),i=i1,i2)
       i1 = i2 + 1
    enddo
    close(unit=iunit)

 case('phantom')
    udist = 0.d0
    umagfd = 0.d0
    if (ix(1) > 0) udist = units(ix(1))
    utime = units(0)
    if (ipmass > 0) then
       umass = units(ipmass)
    else
       print "(a)",' WARNING: units for mass unknown, written as 1.0'
       umass = 1.0d0
    endif
    if (iBfirst > 0) umagfd = units(iBfirst)
    call write_sphdata_phantom(time,gamma,dat,ndim,npart,ntypes,npartoftype,&
                               masstype,ncolumns,udist,utime,umass,umagfd,&
                               labeltype,label,ix,ih,ivx,iBfirst,&
                               ipmass,irho,iutherm,filename,0.)
 case('gadget')
    call write_sphdata_gadget(time,dat,itype,npart,ntypes,npartoftype,&
                               masstype,ncolumns,filename)
 case default
    print "(a)",' ERROR: unknown output format '''//trim(outformat)//''' in write_sphdump'
    return
 end select

end subroutine write_sphdump

subroutine get_order_from_list(list_in,list_out,n)
 integer, intent(in)  :: list_in(:)
 integer, intent(out) :: list_out(:)
 integer, intent(out) :: n
 integer :: i
 logical :: done
 !
 ! set list_out = list_in, but only while list_in
 ! is non-zero and numbers are in range
 !
 i = 1
 done = .false.
 list_out = 0
 do while(.not.done .and. i <= size(list_in) .and. i <= size(list_out))
    if (list_in(i) > 0 .and. list_in(i) <= size(list_out)) then
       list_out(i) = list_in(i)
    else
       done = .true.
    endif
    i = i + 1
 enddo
 n = count(list_out > 0)
 ! if no valid entries, set list_out to 1..n
 if (n==0) then
    n = min(size(list_out),size(list_in))
    do i=1,n
       list_out(i) = i
    enddo
 endif

end subroutine get_order_from_list

end module write_sphdata
