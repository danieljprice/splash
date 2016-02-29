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
!  Copyright (C) 2005-2016 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! the data is stored in the global array dat
!
! THIS VERSION FOR RSPH BINARY DUMPS
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

subroutine read_data(rootname,indexstart,ipos,nstepsread)
  use exact, only:hfact
  use particle_data, only:npartoftype,time,gamma,dat,maxpart,maxstep,maxcol
  use params
!  use labels
  use filenames, only:nfiles
  use settings_data, only:ndim,ndimV,ncolumns,ncalc,ntypes, &
                          buffer_data
  use mem_allocation, only:alloc
  use geometry, only:labelcoordsys
  implicit none
  integer, intent(in) :: indexstart,ipos
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  character(len=len(rootname)+4) :: datfile
  integer :: i,icol,ierr,iunit,isizeheader
  integer :: npart_max,nstep_max
  integer :: ntoti
  logical :: singleprecision

  integer*2, dimension(10) :: sheader
  integer, dimension(10) :: iheader
  integer*8, dimension(10):: lheader
  real, dimension(10) :: rheader
  real(doub_prec), dimension(10) :: dheader
  character(len=100) :: headerstring
  character(len=10), dimension(maxplot) :: cheader
  real, dimension(:,:), allocatable :: dattemp
  real(doub_prec), dimension(:,:), allocatable :: dattempd
  common /chead/ cheader

  iunit = 11 ! file unit number
  nstepsread = 0
  if (rootname(1:1).ne.' ') then
     datfile = trim(rootname)
     !print*,'rootname = ',rootname
  else
     print*,' **** no data read **** '
     return
  endif

  print "(1x,a)",'reading RSPH format'
  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  !--open data file and read data
  !
  open(unit=iunit,iostat=ierr,file=datfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print*,'*** Error opening '//trim(datfile)//' ***'
     return
  endif

  i = indexstart
!
!--read first header line
!
  read(iunit,iostat=ierr,end=80) headerstring
  print "(1x,a)",'header string="'//trim(headerstring)//'"'
!
!--read other header lines (short/normal/long ints, reals, doubles)
!
  read(iunit,iostat=ierr,end=80) sheader(1:10)
  if (ierr /= 0) then
     print*,'WARNING: errors during sheader read'
  endif
  print "(1x,a,10(1x,i2))",'sheader = ',sheader(1:10)

  isizeheader = sheader(1)
  ndim = sheader(2)
  ndimV = sheader(3)
  ntypes = sheader(4)
  ncolumns = sheader(5) + ndim  ! ncolumns in sheader + positions
  !
  !--check for errors in sheader
  !
  if (ndim.gt.3 .or. ndimV.gt.3 .or. ndim.le.0 .or. ndimV.le.0 .or. &
      ncolumns.le.0 ) then
     print*,'*** ERROR: header corrupted: ndim = ',ndim,' ndimV = ', ndimV
     ndim = 0
     ndimV = 0
     ntypes = 0
     ncolumns = 0
     close(iunit)
     return
  elseif (ncolumns.gt.maxplot) then
     print "(1x,a)",'*** WARNING: too many columns for array limits'
     ncolumns = maxplot
     print "(1x,a,i2,a)",'    reading only first ',ncolumns,' columns'
  endif

  read(iunit,iostat=ierr,end=80) iheader(1:isizeheader)
  if (ierr /= 0) then
     print*,'WARNING: errors during iheader read'
  endif
  print "(1x,a,20(1x,i6))",'iheader = ',iheader(1:isizeheader)

  ntoti = sum(iheader(1:ntypes))

  read(iunit,iostat=ierr,end=80) lheader(1:isizeheader)
  if (ierr /= 0) then
     print*,'WARNING: errors during lheader read'
  endif
  print "(1x,a,20(1x,i6))",'lheader = ',lheader(1:isizeheader)

  read(iunit,iostat=ierr,end=80) rheader(1:isizeheader)
  if (ierr /= 0) then
     print*,'WARNING: errors during rheader read'
  endif
  print "(1x,a,20(1x,f8.2))",'rheader = ',rheader
  read(iunit,iostat=ierr,end=80) dheader(1:isizeheader)
  if (ierr /= 0) then
     print*,'WARNING: errors during dheader read'
  endif
  print "(1x,a,20(1x,1pe8.2))",'dheader = ',dheader(1:isizeheader)
  do icol=1,ncolumns-ndim
     read(iunit,iostat=ierr,end=80) cheader(icol)(1:isizeheader)
  enddo
!  print "(a)",(trim(cheader(icol)),icol=1,ncolumns-ndim)
!
!--allocate/reallocate memory for data arrays
!
  if (buffer_data) then
     nstep_max = max(nfiles,maxstep,indexstart)
  else
     nstep_max = max(1,maxstep,indexstart)
  endif
  npart_max = max(int(2.0*ntoti),maxpart)
  if (.not.allocated(dat) .or. ntoti.gt.maxpart  &
       .or. nstep_max.gt.maxstep .or. ncolumns.gt.maxcol) then
     call alloc(npart_max,nstep_max,ncolumns+ncalc)
  endif
!
!--count this as a successful read
!
  nstepsread = nstepsread + 1
!
!--copy header data into now-allocated arrays
!
  npartoftype(1,i) = iheader(1)
  npartoftype(2,i) = iheader(2)
  npartoftype(3,i) = iheader(3)
  ntoti = sum(npartoftype(1:ntypes,i))
!
!--determine whether dump is single or double precision
!
  singleprecision = .true.
  if (all(abs(rheader(1:isizeheader)).lt.tiny(rheader))) then
     singleprecision = .false.
     print "(1x,a)",'double precision dump'
     time(i) = real(dheader(1))
     print*,'particles per smoothing length = ',dheader(2)
     print*,'kernel range = ',dheader(3)
     gamma(i) = real(dheader(4))
  else
     print "(1x,a)",'single precision dump'
     time(i) = rheader(1)
     print*,'particles per smoothing length = ',rheader(2)
     print*,'kernel range = ',rheader(3)
     gamma(i) = rheader(4)
  endif
  hfact = 0.
  print "(/a14,':',f8.4,a8,':',i8,a8,':',i8)",' time',time(i),'npart',npartoftype(1,i),'ntotal',ntoti
  print "(a14,':',i8,a8,':',f8.4,a8,':',f8.4)",' ncolumns',ncolumns,'gamma',gamma(i),'hfact',hfact
  print "(a14,':',i8,a8,':',i8)",'ndim',ndim,'ndimV',ndimV
!
!--read data arrays
!
  if (singleprecision) then
     if (allocated(dattemp)) deallocate(dattemp)
     allocate(dattemp(ndim,ntoti))
     !
     !--read positions of all particles
     !
     read(iunit,iostat=ierr) dattemp(1:ndim,1:ntoti)
     if (ierr /=0 ) then
        print "(a)",'error reading particle positions'
     else
        do icol=1,ndim
           dat(1:ntoti,icol,i) = dattemp(icol,1:ntoti)
        enddo
     endif
     !
     !--read rest of data columns
     !
     do icol=ndim+1,ncolumns
        read(iunit,iostat=ierr) dat(1:ntoti,icol,i)
        if (ierr /= 0) print "(a,i2)", 'error reading column ',icol
     enddo
  else
     if (allocated(dattempd)) deallocate(dattempd)
     allocate(dattempd(ndim,ntoti))
     !
     !--read positions of all particles
     !
     read(iunit,iostat=ierr) dattempd(1:ndim,1:ntoti)
     if (ierr /=0 ) then
        print "(a)",'error reading particle positions'
     else
        do icol=1,ndim
           dat(1:ntoti,icol,i) = real(dattempd(icol,1:ntoti))
        enddo
     endif
     !
     !--read rest of data columns
     !
     do icol=ndim+1,ncolumns
        read(iunit,iostat=ierr) dattempd(1,1:ntoti)
        if (ierr /= 0) print "(a,i2)", 'error reading column ',icol
        !--convert to single precision
        dat(1:ntoti,icol,i) = real(dattempd(1,1:ntoti))
     enddo

  endif

  if (allocated(dattemp)) deallocate(dattemp)
  if (allocated(dattempd)) deallocate(dattempd)
  !
  !--close data file and return
  !
close(unit=iunit)

print "(a)",' finished data read '

return

80 continue
print*,' *** data file empty : no timesteps ***'
return

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
 use labels, only:ix,ivx,ih,irho,iutherm,ipmass,ipr,iBfirst, &
             iamvec,labelvec,label,labeltype
 use params
 use settings_data, only:ndim,ndimV,ncolumns,ntypes, &
                    UseTypeInRenderings
 use geometry, only:labelcoord
 implicit none
 integer :: i,j
 character(len=10), dimension(maxplot) :: cheader
 common /chead/ cheader

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

 do i=ndim+1,ncolumns
    label(i) = cheader(i-ndim)
    !--blank characters in c are ascii zero - correct these to spaces
    do j=1,len(label(i))
       if (iachar(label(i)(j:j)).eq.0) label(i)(j:j) = ' '
    enddo
    !--set positions of various quantities depending on labels
    if (label(i)(1:1)=='m' .or. label(i)(1:4)=='mass') then
       ipmass = i
    elseif (label(i)(1:3)=='rho' .or. label(i)(1:4)=='dens') then
       irho = i
    elseif (label(i)(1:1)=='h' .or. label(i)(1:6)=='smooth') then
       ih = i
    elseif (label(i)(1:2)=='u ' .or. label(i)(1:1)=='e') then
       iutherm = i
    elseif (label(i)(1:2)=='pr' .or. trim(label(i))=='P') then
       ipr = i
    elseif (label(i)(1:1)=='v') then
       if (ivx.eq.0 .or. i.lt.ivx) ivx = i
    elseif (label(i)(1:1)=='B') then
       if (iBfirst.eq.0 .or. i.lt.iBfirst) iBfirst = i
    endif
 enddo

 label(ix(1:ndim)) = labelcoord(1:ndim,1)
 !
 !--label vector quantities (e.g. velocity) appropriately
 !
 if (ivx.gt.0) then
    iamvec(ivx:ivx+ndimV-1) = ivx
    labelvec(ivx:ivx+ndimV-1) = 'v'
    do i=1,ndimV
       label(ivx+i-1) = trim(labelvec(ivx+i-1))//'\d'//labelcoord(i,1)
    enddo
 endif
 if (iBfirst.gt.0) then
    iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
    labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
 endif
!
!--set labels for each type of particles
!
 labeltype(1) = 'gas'
 UseTypeInRenderings(1) = .true.
 if (ntypes.ge.2) then
    labeltype(2) = 'auxiliary'
    UseTypeInRenderings(2) = .true.
 endif
 if (ntypes.ge.3) then
    labeltype(3) = 'mirror'
    UseTypeInRenderings(3) = .true.
 endif

!-----------------------------------------------------------

 return
end subroutine set_labels
