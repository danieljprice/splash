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
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Module implementing "splash to gadget" operation, writing
! a binary dump file suitable for input to the GADGET code
!-----------------------------------------------------------------
module write_data_gadget
 implicit none
 character(len=10), parameter, public :: formatname='gadget'

 public :: write_sphdata_gadget
 private

contains

subroutine write_sphdata_gadget(time,dat,iamtype,ntotal,ntypes,npartoftype, &
                                 masstype,ncolumns,filename)
 use labels,         only:ih,ivx,ix,iutherm,irho,ipmass
 use settings_data,  only:ndim
 use limits,         only:lim
 use params,         only:int1
 implicit none
 integer, intent(in)                          :: ntotal,ntypes,ncolumns
 integer, intent(in), dimension(:)            :: npartoftype
 real, intent(in)                             :: time
 real, intent(in), dimension(ntotal,ncolumns) :: dat
 integer(kind=int1), intent(in), dimension(:) :: iamtype
 real, intent(in), dimension(:)               :: masstype
 character(len=*), intent(in)                 :: filename

 integer, parameter    :: idump = 83
 character(len=len(filename)+10) :: outfile
 integer, dimension(6)                 :: nall,ncrap,noftype
 real(kind=8), dimension(6)            :: massoftype
 real(kind=8)                          :: boxsize,dtime
 real(kind=8), parameter               :: dumz = 0.d0
 real(kind=4), dimension(15)           :: unused
 integer, parameter :: iflagsfr = 0, iflagfeedback = 0, iflagcool = 0
 integer, parameter :: nfiles = 1
 integer            :: ierr,i,j,nmasses,ngas,itype
 integer, dimension(:), allocatable :: iorder
!
!--define output file name
!
 outfile=trim(filename)//'.gadget'
!
!--check if we have enough data to write a GADGET dump
!
 if (ndim.lt.3) then
    print "(a)",' ERROR: ndim < 3 but must be 3 for GADGET data -- cannot write PHANTOM dump, skipping...'
    return
 endif
 if (any(ix(:).le.0)) then
    print "(a)",' ERROR: position labels not set -- cannot write GADGET dump, skipping...'
    return
 endif
 if (ivx.le.0) then
    print "(a)",' ERROR: velocity not found in data -- cannot write GADGET dump, skipping...'
    return
 endif
 if (iutherm.le.0) then
    print "(a)",' ERROR: thermal energy not found in data -- cannot write GADGET dump, skipping...'
    return
 endif
 if (irho.le.0) then
    print "(a)",' ERROR: density not found in data -- cannot write GADGET dump, skipping...'
    return
 endif
 if (ih.le.0) then
    print "(a)",' ERROR: smoothing length not found in data -- cannot write GADGET dump, skipping...'
    return
 endif
!
!--open dumpfile
!
 write(*,"(/,/,'-------->   TIME = ',f10.4,"// &
         "': writing GADGET snapshot file ',a,'   <--------',/)")  time,trim(outfile)

 print "(a)",  ' WARNING: conversion to GADGET format is LIMITED in scope...'
 print "(a)",  '          Currently converts only basic hydro quantities (x,v,m,u,rho,h) '
 print "(a)",  '          and header quantities may be guessed/fudged ...your mileage may vary'
 print "(a,/)",'          (if you use this functionality and want it improved, just let me know)'

 print "(a,i2)",' writing to unit ',idump
 open(unit=idump,file=outfile,status='replace',form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    write(*,*) 'error: can''t create new dumpfile ',trim(outfile)
    return
 endif
 
 massoftype(:)  = 0.
 nall(:)        = 0
 noftype(:)     = 0
 massoftype(1:ntypes) = masstype(1:ntypes)
 nall(1:ntypes)       = npartoftype(1:ntypes)
 noftype(1:ntypes)    = npartoftype(1:ntypes)
 ncrap(:)  = 0
 boxsize   = abs(lim(ix(1),2) - lim(ix(1),1))
 unused(:) = 0
 dtime     = time
 
 write(idump,err=100) noftype(1:6),massoftype(1:6),dtime,dumz, &
                      iflagsfr,iflagfeedback,nall(1:6),iflagcool,nfiles,boxsize, &
                      dumz,dumz,dumz,iflagsfr,iflagsfr,ncrap(1:6),iflagsfr,unused(:)
 !
 !--work out how many particle masses to write
 !
 nmasses = 0
 do j=1,6
    if (massoftype(j).le.0.) then
       nmasses = nmasses + noftype(j)
    endif
 enddo
 print*,'nmasses = ',nmasses
 ngas = npartoftype(1)
 print*,'ngas = ',ngas
 
 if (ntotal > ngas .and. (size(iamtype).gt.1)) then
 !--must print the particles ordered by type
    allocate(iorder(ntotal),stat=ierr)
    j = 0
    do itype=1,min(ntypes,6)
       if (npartoftype(itype).gt.0) then
          do i=1,ntotal
             if (iamtype(i).eq.itype) then
                j = j + 1
                iorder(j) = i
             endif
          enddo
       endif
    enddo
    if (j.lt.ntotal) then
       print*,' ERROR: too many particle types in conversion to gadget format'
       do i=j+1,ntotal
          iorder(i) = i     
       enddo
    endif
    write(idump,err=100) ((dat(iorder(i),ix(j)),j=1,3),i=1,ntotal)
    write(idump,err=100) ((dat(iorder(i),j),j=ivx,ivx+2),i=1,ntotal)
    write(idump,err=100) (iorder(i),i=1,ntotal)  ! particle id
    write(idump,err=100) (dat(iorder(i),ipmass), i=1,nmasses)
    write(idump,err=100) (dat(iorder(i),iutherm),i=1,ngas)
    write(idump,err=100) (dat(iorder(i),irho),   i=1,ngas)
    write(idump,err=100) (2.*dat(iorder(i),ih),  i=1,ngas)
    deallocate(iorder)
 else
    write(idump,err=100) ((dat(i,ix(j)),j=1,3),i=1,ntotal)
    write(idump,err=100) ((dat(i,j),j=ivx,ivx+2),i=1,ntotal)
    write(idump,err=100) (i,i=1,ntotal)  ! particle id
    write(idump,err=100) (dat(i,ipmass), i=1,nmasses)
    write(idump,err=100) (dat(i,iutherm),i=1,ngas)
    write(idump,err=100) (dat(i,irho),   i=1,ngas)
    write(idump,err=100) (2.*dat(i,ih),  i=1,ngas)
 endif
 
 print*,'finished writing file -- OK'

 close(unit=idump)
 return

100 continue
 write(*,*) 'error whilst writing dumpfile '//trim(outfile)
 close(unit=idump)

end subroutine write_sphdata_gadget

end module write_data_gadget
