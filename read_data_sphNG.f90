!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM
! THE NEXT GENERATION SPH CODE
!
! *** CONVERTS TO SINGLE PRECISION ***
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

subroutine read_data(rootname,indexstart,nstepsread)
  use particle_data
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc
  use mem_allocation
  implicit none
  integer, intent(in) :: indexstart
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer, parameter :: maxarrsizes = 10, maxreal = 50
  real, parameter :: pi=3.141592653589
  integer :: i,j,ifile,ierr,iunit,int1,int2,int3,i1,iarr
  integer :: npart_max,nstep_max,ncolstep,icolumn
  integer :: narrsizes,nints,nreals,nreal4s,nreal8s
  integer :: nskip,nprint,itype,ntypes
  logical :: iexist, doubleprec
    
  character(len=len(rootname)+10) :: dumpfile
  character(len=100) :: fileident
  
  integer*8, dimension(maxarrsizes) :: isize
  integer, dimension(maxarrsizes) :: nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8
  real(doub_prec), dimension(:), allocatable :: dattemp
  real(doub_prec) :: udist, utime, umass, r8
  real, dimension(maxreal) :: dummyreal


  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  ifile = 1
  iunit = 15

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
  call alloc(max(1,maxpart),max(ncolumns,1),max(indexstart,maxstep))
  
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
!
!--open the (unformatted) binary file
!
   open(unit=15,iostat=ierr,file=dumpfile,status='old',form='unformatted')
   if (ierr /= 0) then
      print "(a)",'*** ERROR OPENING '//trim(dumpfile)//' ***'
   else
      !
      !--read header key to work out precision
      !
      doubleprec = .true.
      read(iunit,end=55,iostat=ierr) int1,r8,int2,i1,int3
      if (int1.ne.690706) then
         print "(a)",'*** ERROR READING HEADER: wrong endian?'
         return
      endif
      if (int2.ne.780806) then
         print "(a)",'single precision dump'
         doubleprec = .false.
      elseif (int3.ne.690706) then
         print "(a)",'*** WARNING: default int appears to be int*8: not implemented'
      else
         print "(a)",'double precision dump'
      endif
   endif
!
!--read file ID
!     
   read(iunit,end=55,iostat=ierr) fileident
   if (ierr /=0) then
      print "(a)",'*** ERROR READING FILE ID ***'
      return
   else
      print "(a)",'File ID: '//trim(fileident)
   endif
!
!--read number of default ints
!   
   read(iunit,end=55,iostat=ierr) nints
   if (ierr /=0) then
      print "(a)",'error reading nints'
      return
   else
      read(iunit,end=55,iostat=ierr) nprint
      print*,'nprint = ',nprint
   endif
!--int*1, int*2, int*4, int*8
   do i=1,4
      read(iunit,end=55,iostat=ierr) ntypes
      do itype=1,ntypes
         read(iunit,end=55,iostat=ierr)
      enddo
      if (ierr /=0) print "(a)",'error skipping int types'
   enddo
!--default reals
   read(iunit,end=55,iostat=ierr) nreals
   if (ierr /=0) then
      print "(a)",'error reading default reals'
      return
   else
      print*,'nreals = ',nreals
      if (nreals.gt.maxreal) then
         print*,'WARNING: nreal> array size'
         nreals = maxreal
      endif
      if (doubleprec) then
         if (allocated(dattemp)) deallocate(dattemp)
         allocate(dattemp(nreals),stat=ierr)
         if (ierr /=0) print*,'ERROR in memory allocation'
         read(iunit,end=55,iostat=ierr) dattemp(1:nreals)
         dummyreal(1:nreals) = real(dattemp(1:nreals))
      else
         read(iunit,end=55,iostat=ierr) dummyreal(1:nreals)
      endif
!--extract required information
      time(j) = dummyreal(1)
      gamma(j) = dummyreal(3)
      print "(a,1pe12.4,a,0pf8.4)",' time = ',time(j),' gamma = ',gamma(j)
   endif
!--real*4, real*8
   read(iunit,end=55,iostat=ierr) nreal4s
   print "(a,i3)",' nreal4s = ',nreal4s
   if (nreal4s.gt.0) read(iunit,end=55,iostat=ierr) 

   read(iunit,end=55,iostat=ierr) nreal8s
   print "(a,i3)",' ndoubles = ',nreal8s
   read(iunit,end=55,iostat=ierr) udist,umass,utime
   if (ierr /= 0) then
      print "(a)",'*** error reading units'
   endif

!
!--Array headers
!
   read(iunit,end=55,iostat=ierr) narrsizes
   if (ierr /= 0) then 
      print "(a)",'*** error reading number of array sizes ***'
      return
   elseif (narrsizes.gt.maxarrsizes) then
      narrsizes = maxarrsizes
      print "(a,i2)",'WARNING: too many array sizes: reading only ',narrsizes
   endif
   ncolstep = 0
   do iarr=1,narrsizes
      read(iunit,end=55,iostat=ierr) isize(iarr),nint(iarr),nint1(iarr),nint2(iarr), &
                 nint4(iarr),nint8(iarr),nreal(iarr),nreal4(iarr),nreal8(iarr)
      if (iarr.eq.1) npartoftype(1,j) = isize(iarr)
      print *,'size ',iarr,' dim = ',isize(iarr),'nint=',nint(iarr),nint1(iarr), &
            nint2(iarr),nint4(iarr),nint8(iarr),'nreal =',nreal(iarr),nreal4(iarr),nreal8(iarr)
!--we are going to read all real arrays but need to convert them all to default real
      if (isize(iarr).eq.isize(1)) then
         ncolstep = ncolstep + nreal(iarr) + nreal4(iarr) + nreal8(iarr)
      endif
   enddo
   
   npart_max = maxval(isize(1:narrsizes))
!
!--allocate memory for all columns
!
   call alloc(npart_max,j,ncolstep)
   nstepsread = nstepsread + 1
   ncolumns = ncolstep + ncalc
   icolumn = 0
!
!--Arrays
!
   do iarr=1,narrsizes
!--skip integer arrays (not needed for plotting)
      nskip = nint(iarr) + nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
      !!print*,'skipping ',nskip,' isize = ',isize(iarr)
      do i=1,nskip
         read(iunit,end=55,iostat=ierr)
      enddo
!--skip real arrays if size different      
      if (isize(iarr).ne.isize(1)) then
         nskip = nreal(iarr) + nreal4(iarr) + nreal8(iarr)
         do i=1,nskip
            read(iunit,end=55,iostat=ierr)
         enddo
      else
!--otherwise read them      
         if (allocated(dattemp)) deallocate(dattemp)
         allocate(dattemp(isize(iarr)),stat=ierr)
         if (ierr /=0) print "(a)",'ERROR in memory allocation'
!        default reals may need converting
         do i=1,nreal(iarr)
            if (doubleprec) then
               read(iunit,end=55,iostat=ierr) dattemp(1:isize(iarr))
               icolumn = icolumn + 1
               dat(1:isize(iarr),icolumn,j) = real(dattemp(1:isize(iarr)))
            else
               icolumn = icolumn + 1
               read(iunit,end=55,iostat=ierr) dat(1:isize(iarr),icolumn,j)
            endif
         enddo
!        real4's go straight into dat
         do i=1,nreal4(iarr)
            icolumn = icolumn + 1
            read(iunit,end=55,iostat=ierr) dat(1:isize(iarr),icolumn,j) 
         enddo
!        real 8's need converting
         do i=1,nreal8(iarr)
            icolumn = icolumn + 1
            read(iunit,end=55,iostat=ierr) dattemp(1:isize(iarr))
            dat(1:isize(iarr),icolumn,j) = real(dattemp(1:isize(iarr)))
         enddo
      endif
   enddo

55 continue
!
!--reached end of file
!
close(15)
   
return
                    
end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels
  use params
  use settings_data
  use geometry, only:labelcoord
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
  ipmass = 4   !  particle mass
  ih = 5       !  smoothing length
  ivx = 6
  iutherm = 9  !  thermal energy
  irho = 10     ! location of rho in data array
  if (ncolumns.gt.10) then
     label(11) = 'dgrav'
  endif
  
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  do i=1,ndimV
     label(ivx+i-1) = 'v\d'//labelcoord(i,1)
  enddo
  label(irho) = 'density (g/cm\u3\d)'
  label(iutherm) = 'u'
  label(ih) = 'h       '
  label(ipmass) = 'particle mass'     

    !
  !--set labels for vector quantities
  !
  iamvec(ivx:ivx+ndimV-1) = ivx
  labelvec(ivx:ivx+ndimV-1) = 'v'
  do i=1,ndimV
     label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
  enddo
  !
  !--set labels for each particle type
  !
  ntypes = 2  !!maxparttypes
  labeltype(1) = 'gas'
  !!labeltype(2) = 'ghost'
  labeltype(2) = 'sink'
 
!-----------------------------------------------------------

  return 
end subroutine set_labels
