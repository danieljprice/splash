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
!  Copyright (C) 2005-2026 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Module implementing "splash to ndspmhd" operation, writing
! a binary dump file suitable for input to the ndspmhd code
!-----------------------------------------------------------------
module write_data_ndspmhd
 implicit none
 character(len=10), parameter, public :: formatname='ndspmhd'

 public :: write_sphdata_ndspmhd
 private

contains

!-----------------------------------------------------------------
!+
!  write output dump readable by ndspmhd (and splash -f ndspmhd)
!+
!-----------------------------------------------------------------
subroutine write_sphdata_ndspmhd(time,gamma,dat,iamtype,ntotal,ntypes,npartoftype, &
                                 masstype,ncolumns,filename)
 use labels,         only:ih,ivx,ix,iutherm,irho,ipmass,iBfirst,headertags,&
                          find_column
 use settings_data,  only:ndim,ndimV
 use params,         only:int1
 use particle_data,  only:headervals
 use asciiutils,     only:get_value
 integer, intent(in)                          :: ntotal,ntypes,ncolumns
 integer, intent(in), dimension(:)            :: npartoftype
 real, intent(in)                             :: time,gamma
 real, intent(in), dimension(ntotal,ncolumns) :: dat
 integer(kind=int1), intent(in), dimension(:) :: iamtype
 real, intent(in), dimension(:)               :: masstype
 character(len=*), intent(in)                 :: filename

 integer, parameter :: idump = 83
 character(len=len(filename)+10) :: outfile
 character(len=12) :: geom
 integer :: ierr,i,ngas,iformat,ncol_out,ialpha,ialphau,idivv
 integer :: ibound(3)
 integer, allocatable :: igas(:)
 real, allocatable :: coldat(:)
 real :: hfact,xmin(3),xmax(3),pmass_default

 outfile = trim(filename)//'.dat'
 geom = 'cartesian'
 ibound = 0
 xmin = 0.
 xmax = 0.
!
!--require positions and enough spatial dimensions
!
 if (ndim <= 0 .or. ndim > 3) then
    print "(a)",' ERROR: ndim not set -- cannot write ndspmhd dump, skipping...'
    return
 endif
 if (any(ix(1:ndim) <= 0)) then
    print "(a)",' ERROR: position labels not set -- cannot write ndspmhd dump, skipping...'
    return
 endif
 if (ndimV <= 0) then
    print "(a)",' ERROR: ndimV not set -- cannot write ndspmhd dump, skipping...'
    return
 endif
!
!--gas particles only
!
 call select_gas_particles(ntotal,npartoftype,iamtype,igas,ngas)
 if (ngas <= 0) then
    print "(a)",' ERROR: no gas particles found -- cannot write ndspmhd dump, skipping...'
    return
 endif
!
!--header quantities from splash header tags when present
!
 hfact = 1.2
 if (allocated(headervals)) then
    hfact = get_value('hfact',headertags,headervals(:,1),default=1.2)
    do i=1,ndim
       xmin(i) = get_header_bound(headertags,headervals(:,1),'xmin',i)
       xmax(i) = get_header_bound(headertags,headervals(:,1),'xmax',i)
    enddo
 endif
!
!--format and column count (match ndspmhd write_dump hydro / MHD)
!
 if (iBfirst > 0) then
    iformat = 2
    ncol_out = ndim + 2*ndimV + 4 + 8 + 2*ndimV
 else
    iformat = 1
    ncol_out = ndim + 2*ndimV + 4 + 5
 endif

 write(*,"(/,/,'-------->   TIME = ',es12.4,"// &
         "': writing ndspmhd dump file ',a,'   <--------',/)") time,trim(outfile)
 print "(a)",' writing gas particles only (sinks/other types dropped)'
 print "(a,i8,a,i2,a,i2)",' ngas = ',ngas,'  iformat = ',iformat,'  ncolumns = ',ncol_out

 open(unit=idump,file=outfile,status='replace',form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR: cannot create '//trim(outfile)
    deallocate(igas)
    return
 endif

 write(idump,iostat=ierr) time,ngas,ngas,gamma,hfact,ndim,ndimV, &
      ncol_out,iformat,ibound(1:ndim),xmin(1:ndim),xmax(1:ndim),len(geom),geom
 if (ierr /= 0) then
    print "(a)",' ERROR writing ndspmhd header'
    close(idump)
    deallocate(igas)
    return
 endif

 allocate(coldat(ngas),stat=ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR allocating temporary array for ndspmhd write'
    close(idump)
    deallocate(igas)
    return
 endif

 pmass_default = 0.
 if (ntypes >= 1) pmass_default = masstype(1)
 ialpha  = find_column('alpha',ncolumns)
 ialphau = find_column('alphau',ncolumns)
 idivv   = find_column('div v',ncolumns)
!
!--essential columns: x, v, h, rho, u, pmass
!  (stop on first write failure so later iostats cannot mask it)
!
 ierr = 0
 do i=1,ndim
    call fill_column(coldat,dat,igas,ngas,ix(i),0.)
    call write_coldat(idump,coldat,ierr)
 enddo
 do i=1,ndimV
    if (ivx > 0) then
       call fill_column(coldat,dat,igas,ngas,ivx+i-1,0.)
    else
       coldat = 0.
    endif
    call write_coldat(idump,coldat,ierr)
 enddo
 call fill_column(coldat,dat,igas,ngas,ih,0.)
 call write_coldat(idump,coldat,ierr)
 call fill_column(coldat,dat,igas,ngas,irho,0.)
 call write_coldat(idump,coldat,ierr)
 call fill_column(coldat,dat,igas,ngas,iutherm,0.)
 call write_coldat(idump,coldat,ierr)
 if (ipmass > 0) then
    call fill_column(coldat,dat,igas,ngas,ipmass,pmass_default)
 else
    coldat = pmass_default
 endif
 call write_coldat(idump,coldat,ierr)

 if (iformat==2) then
    !
    !--MHD: alpha(3), B, psi, then info columns
    !
    call fill_column(coldat,dat,igas,ngas,ialpha,0.)
    call write_coldat(idump,coldat,ierr)
    call fill_column(coldat,dat,igas,ngas,ialphau,0.)
    call write_coldat(idump,coldat,ierr)
    coldat = 0.  ! alphaB
    call write_coldat(idump,coldat,ierr)
    do i=1,ndimV
       call fill_column(coldat,dat,igas,ngas,iBfirst+i-1,0.)
       call write_coldat(idump,coldat,ierr)
    enddo
    coldat = 0.  ! psi
    call write_coldat(idump,coldat,ierr)
    call write_pressure(idump,coldat,dat,igas,ngas,gamma,ierr)
    call fill_column(coldat,dat,igas,ngas,idivv,0.)
    call write_coldat(idump,coldat,ierr)
    coldat = 0.  ! div B
    call write_coldat(idump,coldat,ierr)
    do i=1,ndimV
       coldat = 0.  ! J
       call write_coldat(idump,coldat,ierr)
    enddo
    coldat = 0.  ! grad h
    call write_coldat(idump,coldat,ierr)
    do i=1,ndimV
       coldat = 0.  ! force
       call write_coldat(idump,coldat,ierr)
    enddo
 else
    !
    !--hydro: alpha(2), P, div v, grad h, force
    !
    call fill_column(coldat,dat,igas,ngas,ialpha,0.)
    call write_coldat(idump,coldat,ierr)
    call fill_column(coldat,dat,igas,ngas,ialphau,0.)
    call write_coldat(idump,coldat,ierr)
    call write_pressure(idump,coldat,dat,igas,ngas,gamma,ierr)
    call fill_column(coldat,dat,igas,ngas,idivv,0.)
    call write_coldat(idump,coldat,ierr)
    coldat = 0.  ! grad h
    call write_coldat(idump,coldat,ierr)
    do i=1,ndimV
       coldat = 0.  ! force
       call write_coldat(idump,coldat,ierr)
    enddo
 endif

 if (ierr /= 0) then
    print "(a)",' ERROR writing ndspmhd particle data'
    close(idump)
    deallocate(coldat,igas)
    return
 endif
!
!--particle types (all gas)
!
 write(idump,iostat=ierr) (0,i=1,ngas)
 if (ierr /= 0) then
    print "(a)",' ERROR writing ndspmhd particle types'
 else
    print "(a)",' finished writing file -- OK'
 endif

 close(idump)
 deallocate(coldat,igas)

end subroutine write_sphdata_ndspmhd

!-----------------------------------------------------------------
!+
!  select gas particle indices (type 1 only)
!+
!-----------------------------------------------------------------
subroutine select_gas_particles(ntotal,npartoftype,iamtype,igas,ngas)
 use params, only:int1
 integer, intent(in) :: ntotal
 integer, intent(in) :: npartoftype(:)
 integer(kind=int1), intent(in) :: iamtype(:)
 integer, allocatable, intent(out) :: igas(:)
 integer, intent(out) :: ngas
 integer :: i,ierr

 ngas = 0
 if (size(iamtype) > 1) then
    do i=1,ntotal
       if (iamtype(i)==1) ngas = ngas + 1
    enddo
 else
    ngas = 0
    if (size(npartoftype) >= 1) ngas = npartoftype(1)
    if (ngas <= 0) ngas = ntotal
 endif

 allocate(igas(max(ngas,1)),stat=ierr)
 if (ierr /= 0) then
    ngas = 0
    return
 endif

 if (size(iamtype) > 1) then
    ngas = 0
    do i=1,ntotal
       if (iamtype(i)==1) then
          ngas = ngas + 1
          igas(ngas) = i
       endif
    enddo
 else
    do i=1,ngas
       igas(i) = i
    enddo
 endif

end subroutine select_gas_particles

!-----------------------------------------------------------------
!+
!  write one column array; skip if a previous write already failed
!+
!-----------------------------------------------------------------
subroutine write_coldat(iunit,coldat,ierr)
 integer, intent(in) :: iunit
 real, intent(in) :: coldat(:)
 integer, intent(inout) :: ierr
 integer :: ierr_write

 if (ierr /= 0) return
 write(iunit,iostat=ierr_write) coldat
 ierr = ierr_write

end subroutine write_coldat

!-----------------------------------------------------------------
!+
!  copy one data column for the selected gas particles
!+
!-----------------------------------------------------------------
subroutine fill_column(coldat,dat,igas,ngas,icol,default)
 real, intent(out) :: coldat(:)
 real, intent(in)  :: dat(:,:)
 integer, intent(in) :: igas(:),ngas,icol
 real, intent(in) :: default
 integer :: i

 if (icol > 0 .and. icol <= size(dat,2)) then
    do i=1,ngas
       coldat(i) = dat(igas(i),icol)
    enddo
 else
    coldat(1:ngas) = default
 endif

end subroutine fill_column

!-----------------------------------------------------------------
!+
!  write pressure column (use ipr, else (gamma-1)*rho*u, else 0)
!+
!-----------------------------------------------------------------
subroutine write_pressure(iunit,coldat,dat,igas,ngas,gamma,ierr)
 use labels, only:ipr,irho,iutherm
 integer, intent(in) :: iunit,igas(:),ngas
 real, intent(inout) :: coldat(:)
 real, intent(in) :: dat(:,:),gamma
 integer, intent(inout) :: ierr
 integer :: i
 real :: rhoi,ui

 if (ierr /= 0) return
 if (ipr > 0) then
    call fill_column(coldat,dat,igas,ngas,ipr,0.)
 elseif (irho > 0 .and. iutherm > 0) then
    do i=1,ngas
       rhoi = dat(igas(i),irho)
       ui   = dat(igas(i),iutherm)
       coldat(i) = (gamma - 1.)*rhoi*ui
    enddo
 else
    coldat(1:ngas) = 0.
 endif
 call write_coldat(iunit,coldat(1:ngas),ierr)

end subroutine write_pressure

!-----------------------------------------------------------------
!+
!  look up xmin/xmax from header tags (xmin, xmax, or xmin1..)
!+
!-----------------------------------------------------------------
real function get_header_bound(tags,vals,base,idim)
 use asciiutils, only:match_tag
 character(len=*), intent(in) :: tags(:)
 real, intent(in) :: vals(:)
 character(len=*), intent(in) :: base
 integer, intent(in) :: idim
 character(len=16) :: tag
 integer :: itag
 character(len=1), parameter :: xyz(3) = (/'x','y','z'/)

 get_header_bound = 0.
 if (idim < 1 .or. idim > 3) return

 ! plain xmin / xmax applies to first dimension only
 if (idim==1) then
    itag = match_tag(tags,trim(base))
    if (itag > 0 .and. itag <= size(vals)) then
       get_header_bound = vals(itag)
       return
    endif
 endif

 ! xmin1 / xmax1 style
 write(tag,"(a,i1)") trim(base),idim
 itag = match_tag(tags,trim(tag))
 if (itag > 0 .and. itag <= size(vals)) then
    get_header_bound = vals(itag)
    return
 endif

 ! x/y/z min/max
 if (trim(base)=='xmin') then
    tag = xyz(idim)//'min'
 elseif (trim(base)=='xmax') then
    tag = xyz(idim)//'max'
 else
    return
 endif
 itag = match_tag(tags,trim(tag))
 if (itag > 0 .and. itag <= size(vals)) get_header_bound = vals(itag)

end function get_header_bound

end module write_data_ndspmhd
