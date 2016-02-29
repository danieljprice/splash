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
! the data is stored in the global array dat
!
! THIS VERSION FOR ABDEL REHEAM's SPH CODE
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
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!-------------------------------------------------------------------------
!--local module to store header information so we can later set the labels
module alydataread
 use params
 use labels, only:lenlabel
 implicit none
 character(len=16), dimension(maxplot) :: compName

end module alydataread

subroutine read_data(rootname,indexstart,ipos,nstepsread)
  use particle_data,  only:npartoftype,masstype,time,gamma,dat,maxpart,maxstep,maxcol,iamtype
  use params
  use filenames,      only:nfiles
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc, &
                          buffer_data,iverbose,debugmode,ntypes
  use mem_allocation, only:alloc
  use labels,         only:ipr,ivx,ih,irho,labeltype
  use alydataread,    only:compName
  implicit none
  integer,          intent(in)  :: indexstart,ipos
  integer,          intent(out) :: nstepsread
  character(len=*), intent(in)  :: rootname
  character(len=len(rootname)+4) :: datfile
  integer :: i,ierr,iunit,j,iblock
  integer :: npart_max,nstep_max

  integer, dimension(:), allocatable :: itype
  character(len=20) :: geomfile
  character(len=7)  :: keyword
  character(len=70) :: title
  character(len=1)  :: dumchar
  character(len=16) :: unitsys
  integer :: istep,jtype,np,ione,kk,idum,nblock
  real(kind=sing_prec) :: version,timesingle,dum
  real(kind=doub_prec) :: versiond,timedbl,dumd
  real :: timein,dx,dy
  logical :: singleprecision

  iunit = 11 ! file unit number
  nstepsread = 0
  if (rootname(1:1).ne.' ') then
     datfile = trim(rootname)
  else
     print*,' **** no data read **** '
     return
  endif

  if (iverbose.ge.1) print "(1x,a)",'reading Aly format'
  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)

  ndim = 2
  ndimV = 2
  !
  !--open data file and read data
  !
  open(unit=iunit,iostat=ierr,file=datfile,status='old',form='unformatted',access='stream')
  if (ierr /= 0) then
     print*,' *** Error opening '//trim(datfile)//' ***'
     return
  endif
!
!--read first header line
!
  read(iunit,iostat=ierr,end=80) keyword
  read(iunit,iostat=ierr,end=80) version
  read(iunit,iostat=ierr,end=80) title
  read(iunit,iostat=ierr,end=80) istep
  read(iunit,iostat=ierr,end=80) timesingle
  timein = timesingle
  read(iunit,iostat=ierr,end=80) np
  read(iunit,iostat=ierr,end=80) ione
  if (ierr /= 0 .or. np <= 0 .or. timesingle.lt.0.) then
     !
     !--try single precision
     !
     rewind(iunit)
     read(iunit,iostat=ierr,end=80) keyword
     read(iunit,iostat=ierr,end=80) versiond
     version = versiond
     read(iunit,iostat=ierr,end=80) title
     read(iunit,iostat=ierr,end=80) istep
     read(iunit,iostat=ierr,end=80) timedbl
     timein = timedbl
     read(iunit,iostat=ierr,end=80) np
     read(iunit,iostat=ierr,end=80) ione
     singleprecision = .false.
     if (ierr /= 0 .or. np < 0 .or. timedbl < 0.) then
        print "(a)",' *** Error reading first header ***'
        close(iunit)
        return
     endif
  endif
  print "(a,f4.2)",' keyword = '//trim(keyword)//' version = ',version
  !print "(a)",' title = '//trim(title)
  print "(a,i6,a,es10.3,a,i6)",' step = ',istep,' time = ',timein,' np = ',np
!
!--allocate memory for data arrays
!
  if (buffer_data) then
     nstep_max = max(nfiles,maxstep,indexstart)
  else
     nstep_max = max(1,maxstep,indexstart)
  endif
  npart_max = max(int(1.1*np),maxpart)
  ncolumns = 7
  if (.not.allocated(dat) .or. npart_max.gt.maxpart  &
       .or. nstep_max.gt.maxstep .or. ncolumns+ncalc.gt.maxcol) then
     call alloc(npart_max,nstep_max,ncolumns+ncalc,mixedtypes=.true.)
  endif

  i = indexstart
  nstepsread = 0
  time(i)  = timein
  gamma(i) = 5./3.
  npartoftype(1,i) = np
  
  do j=1,np
     read(iunit,iostat=ierr,end=67) idum,(dat(j,kk,i),kk=1,2),dum
  enddo
  read(iunit,iostat=ierr,end=67) nblock
  read(iunit,iostat=ierr,end=67) (idum,j=1,nblock)

  ntypes = 4
  allocate(itype(np))
  read(iunit,iostat=ierr,end=67) (itype(j),j=1,nblock)
  npartoftype(:,i) = 0
  do j=1,np
     jtype = itype(j)
     !--map types from code types to splash types
     select case(jtype)
     case(1)
        jtype = 2 ! boundary
     case(2)
        jtype = 1 ! water
     case default ! unknown
        jtype = 4
     end select
     iamtype(j,i) = jtype
     npartoftype(jtype,i) = npartoftype(jtype,i) + 1
  enddo
  deallocate(itype)
  read(iunit,iostat=ierr,end=67) (dumchar,j=1,nblock)
  read(iunit,iostat=ierr,end=67) (idum,j=1,nblock)

  call set_labels
  do iblock=1,1
  
     read(iunit,iostat=ierr,end=67) nblock
     read(iunit,iostat=ierr,end=67) idum

     do j=1,nblock
        read(iunit,iostat=ierr,end=67) compName(j),unitsys,idum,idum,dum
        !print*,trim(compName(j))
     enddo

     do j=1,np
        read(iunit,iostat=ierr,end=67) dum,dat(j,ipr,i),dat(j,ivx,i),dat(j,ivx+1,i),dat(j,1,i),dat(j,2,i)
        if (abs(dum).lt.tiny(0.)) then
           npartoftype(iamtype(j,i),i) = npartoftype(iamtype(j,i),i) - 1 ! remove from previous type
           iamtype(j,i) = 3
           npartoftype(3,i) = npartoftype(3,i) + 1 ! add to "box" type
        endif
     enddo
  enddo
  !
  !--fake other properties: density, mass, smoothing length etc.
  !
  masstype(1,i) = 1./npartoftype(1,i)
  !
  !--assume smoothing length to be the max dimension divided by the number of particles^(1/ndim)
  !
  dx = maxval(dat(1:np,1,i)) - minval(dat(1:np,1,i))
  dy = maxval(dat(1:np,2,i)) - minval(dat(1:np,2,i))
  dat(:,ih,i)   = dx/(npartoftype(1,i))**(1./ndim)*(dy/dx)
  print*,' WARNING: ASSUMING SMOOTHING LENGTH = ',dat(1,ih,i),' AND ARBITRARY PARTICLE MASSES'
  dat(:,irho,i) = 1.
  nstepsread = 1
  do jtype=1,ntypes
     write(*,"(' n(',a,') = ',i6)",advance="no") trim(labeltype(jtype)),npartoftype(jtype,i)
  enddo
  write(*,*)
  !read*


   goto 68
67 continue
   print "(a)",' > end of file reached <'
68 continue
  !
  !--close data file and return
  !
close(unit=11)

if (debugmode) print*,'DEBUG> Read steps ',indexstart,'->',indexstart + nstepsread - 1, &
       ' last step ntot = ',sum(npartoftype(:,indexstart+nstepsread-1))
return

80 continue
print*,' *** data file empty : no timesteps ***'
return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
 use labels, only:ix,ivx,ih,irho,ipr,&
             iamvec,labelvec,label,labeltype
 use params
 use settings_data, only:ndim,ndimV,UseTypeInRenderings
 use geometry,      only:labelcoord
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
 ipr = ndim + ndimV + 1
 irho = ipr + 1
 ih  = irho + 1
 label(ipr) = 'pressure'
 label(irho) = 'density'
 label(ih) = 'h'

 label(ix(1:ndim)) = labelcoord(1:ndim,1)
 !
 !--label vector quantities (e.g. velocity) appropriately
 !
 iamvec(ivx:ivx+ndimV-1) = ivx
 labelvec(ivx:ivx+ndimV-1) = 'v'

!
!--set labels for each type of particles
!
 labeltype(1) = 'water'
 labeltype(2) = 'boundary'
 labeltype(3) = 'box'
 labeltype(4) = 'unknown'
 UseTypeInRenderings(1) = .true.
 UseTypeInRenderings(2) = .false.
 UseTypeInRenderings(3) = .false.
 UseTypeInRenderings(4) = .false.

!-----------------------------------------------------------

 return
end subroutine set_labels
