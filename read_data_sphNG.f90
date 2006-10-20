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
  use settings_data, only:ndim,ndimV,ncolumns,ncalc,units,unitslabel,iformat
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
  integer :: nskip,ntotal,npart,itype,ntypes
  integer :: ipos,nptmass,nptmassi,nunknown
  integer :: nhydroarrays,nmhdarrays,imaxcolumnread,nhydroarraysinfile
  logical :: iexist, doubleprec, smalldump,imadepmasscolumn
    
  character(len=len(rootname)+10) :: dumpfile
  character(len=100) :: fileident
  
  integer*8, dimension(maxarrsizes) :: isize
  integer, dimension(maxarrsizes) :: nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8
  integer*1, dimension(:), allocatable :: iphase
  real(doub_prec), dimension(:), allocatable :: dattemp
  real(doub_prec) :: udist, utime, umass, umagfd, r8
  real, dimension(maxreal) :: dummyreal
  real, dimension(:,:), allocatable :: dattemp2
  real :: pmassinitial


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
 !!! call alloc(max(1,maxpart),max(ncolumns,1),max(indexstart,maxstep))
  
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)
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
      read(iunit,end=55,iostat=ierr) int1,r8,int2,i1,int3
      if (ierr /= 0) then
         print "(a)",'*** ERROR READING FILE : corrupted/zero size/wrong endian?'
         close(iunit)
         return
      elseif (int1.ne.690706) then
         print "(a)",'*** ERROR READING HEADER: wrong endian?'
         close(iunit)
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
   read(iunit,iostat=ierr) fileident
   if (ierr /=0) then
      print "(a)",'*** ERROR READING FILE ID ***'
      close(iunit)
      return
   else
      print "(a)",'File ID: '//trim(fileident)
   endif
   iformat = 0
   smalldump = .false.
   imadepmasscolumn = .false.
   if (fileident(1:1).eq.'S') then
      smalldump = .true.
   endif
!
!--read number of default ints
!   
   read(iunit,iostat=ierr) nints
   if (ierr /=0) then
      print "(a)",'error reading nints'
      close(iunit)
      return
   else
      read(iunit,iostat=ierr) npart
      if (ierr /=0) then
         print "(a)",'error reading npart'
         close(iunit)
         return
      else
         print*,'npart = ',npart
      endif
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
      close(iunit)
      return
   else
!      print*,'nreals = ',nreals
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
   endif
!--real*4, real*8
   read(iunit,end=55,iostat=ierr) nreal4s
!   print "(a,i3)",' nreal4s = ',nreal4s
   if (nreal4s.gt.0) read(iunit,end=55,iostat=ierr) 

   read(iunit,end=55,iostat=ierr) nreal8s
!   print "(a,i3)",' ndoubles = ',nreal8s
   print "(4(a,i3))",' header contains ',nreals,' reals,',nreal4s,' real4s, ',nreal8s,' doubles'
   if (nreal8s.ge.4) then
      read(iunit,end=55,iostat=ierr) udist,umass,utime,umagfd   
   elseif (nreal8s.ge.3) then
      read(iunit,end=55,iostat=ierr) udist,umass,utime
   else
      print "(a)",'*** error: units not found in file'
      udist = 1.0
      umass = 1.0
      utime = 1.0
      umagfd = 1.0
      read(iunit,end=55,iostat=ierr)
   endif
   if (ierr /= 0) then
      print "(a)",'*** error reading units'
   endif
!   npower = int(log10(udist))
!   udist = udist/10.**npower
!   udistAU = udist/1.495979e13
   units(1:3) = udist
   unitslabel(1:3) = ' [cm]'
!   do i=1,3
!      write(unitslabel(i),"('[ 10\u',i2,'\d cm]')") npower
!   enddo
   units(4) = umass
   unitslabel(4) = ' [g]'
   units(5) = udist
   unitslabel(5) = ' [cm]'
   units(6:8) = udist/utime
   unitslabel(6:8) = ' [cm/s]'
   units(9) = (udist/utime)**2
   unitslabel(9) = ' [erg/g]'
   units(10) = umass/udist**3
   unitslabel(10) = ' [g/cm\u3\d]'
   units(0) = utime/3.1536e7
   unitslabel(0) = ' yrs'
   
!
!--Array headers
!
   read(iunit,end=55,iostat=ierr) narrsizes
   if (ierr /= 0) then 
      print "(a)",'*** error reading number of array sizes ***'
      close(iunit)
      return
   elseif (narrsizes.gt.maxarrsizes) then
      narrsizes = maxarrsizes
      print "(a,i2)",'WARNING: too many array sizes: reading only ',narrsizes
   endif
   ncolstep = 0
   do iarr=1,narrsizes
      read(iunit,end=55,iostat=ierr) isize(iarr),nint(iarr),nint1(iarr),nint2(iarr), &
                 nint4(iarr),nint8(iarr),nreal(iarr),nreal4(iarr),nreal8(iarr)
      if (iarr.eq.1) ntotal = isize(iarr)
      if (isize(iarr).gt.0) then
         print *,'block ',iarr,' dim = ',isize(iarr),'nint=',nint(iarr),nint1(iarr), &
            nint2(iarr),nint4(iarr),nint8(iarr),'nreal =',nreal(iarr),nreal4(iarr),nreal8(iarr)
      endif
!--we are going to read all real arrays but need to convert them all to default real
      if (isize(iarr).eq.isize(1)) then
         ncolstep = ncolstep + nreal(iarr) + nreal4(iarr) + nreal8(iarr)
      endif
   enddo
   
   npart_max = maxval(isize(1:narrsizes))
   if (smalldump) then
      if (nreals.ge.15) then
         pmassinitial = dummyreal(15)
         if (pmassinitial.gt.tiny(0.)) then
            ncolstep = ncolstep + 1  ! make an extra column to contain particle mass
            imadepmasscolumn = .true.
         endif
      else
         print "(a)",' error extracting pmassinitial from small dump file'
         pmassinitial = 0.
      endif
      if (abs(pmassinitial).lt.tiny(0.) .and. nreal(1).lt.4) then
         print "(a)",' error: particle masses not present in small dump file'   
      endif
   endif
!
!--allocate memory for all columns
!
   if (npart_max.gt.maxpart .or. j.gt.maxstep .or. (ncolstep+ncalc).gt.maxcol) then
      call alloc(npart_max,j,ncolstep+ncalc)
   endif
!--extract required information
   time(j) = dummyreal(1)
   gamma(j) = dummyreal(3)
   npartoftype(1,j) = npart
   npartoftype(2,j) = ntotal - npart
   nptmass = isize(2)
   print "(a,1pe12.4,a,0pf8.4)",' time = ',time(j),' gamma = ',gamma(j)
   nstepsread = nstepsread + 1
   ncolumns = ncolstep + ncalc
   icolumn = 0
   if (allocated(iphase)) deallocate(iphase)
   allocate(iphase(npart_max))
   iphase(:) = 0
!
!--to handle both small and full dumps, we need to place the quantities dumped
!  in both small and full dumps at the start of the dat array
!  quantities only in the full dump then come after
!  also means that hydro/MHD are "semi-compatible" in the sense that x,y,z,m,h and rho
!  are in the same place for both types of dump
!
   nhydroarrays = 6 ! x,y,z,m,h,rho
   nhydroarraysinfile = nreal(1) + nreal4(1) + nreal8(1)
   if (imadepmasscolumn) nhydroarraysinfile = nhydroarraysinfile + 1
   if (nhydroarraysinfile .lt.nhydroarrays) then
      print "(a)",' ERROR: one of x,y,z,m,h or rho missing in small dump read'
      nhydroarrays = nreal(1)+nreal4(1)+nreal8(1)
   endif
   if (narrsizes.ge.4) then
      nmhdarrays = 3 ! Bx,By,Bz
   else
      nmhdarrays = 0
   endif
   imaxcolumnread = 0
   iformat = 0 ! hydro full dump
   if (smalldump) iformat = 1 ! hydro small dump
   if (narrsizes.ge.4) then
      iformat = 2 ! mhd full dump
      if (smalldump) iformat = 3 ! mhd small dump
   endif
!
!--Arrays
!
   do iarr=1,narrsizes
!--read iphase from array block 1 if sinks are present
      if (iarr.eq.1 .and. narrsizes.gt.1 .and. isize(2).gt.0) then
         !--skip default int
         nskip = nint(iarr)
         do i=1,nskip
            read(iunit,end=33,iostat=ierr)
         enddo
         if (nint1(iarr).lt.1) then
            print "(a)",'ERROR: sinks present but can''t locate iphase in dump'
            !--skip remaining integer arrays
            nskip = nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
         else
            read(iunit,end=33,iostat=ierr) iphase(1:isize(iarr))
            !--skip remaining integer arrays
            nskip = nint1(iarr) - 1 + nint2(iarr) + nint4(iarr) + nint8(iarr)
         endif
      else
!--otherwise skip all integer arrays (not needed for plotting)
         nskip = nint(iarr) + nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
      endif
      !!print*,'skipping ',nskip,' isize = ',isize(iarr)
      do i=1,nskip
         read(iunit,end=33,iostat=ierr)
      enddo
!--skip real arrays if size different      
      if (isize(iarr).ne.isize(1)) then
         nskip = nreal(iarr) + nreal4(iarr) + nreal8(iarr)
         do i=1,nskip
            read(iunit,end=33,iostat=ierr)
         enddo
      else
!--otherwise read them      
         if (allocated(dattemp)) deallocate(dattemp)
         allocate(dattemp(isize(iarr)),stat=ierr)
         if (ierr /=0) print "(a)",'ERROR in memory allocation'

!        default reals may need converting
         do i=1,nreal(iarr)
            if (iarr.eq.1 .and. i.eq.6) then
               ! read x,y,z,m,h and then place arrays after always-present ones
               icolumn = nhydroarrays+nmhdarrays + 1
               imaxcolumnread = icolumn
            elseif (iarr.eq.4 .and. i.le.3) then
               icolumn = 6 + i
            else
               icolumn = imaxcolumnread + 1
               imaxcolumnread = icolumn
            endif
            if (doubleprec) then
               read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
               dat(1:isize(iarr),icolumn,j) = real(dattemp(1:isize(iarr)))
            else
               read(iunit,end=33,iostat=ierr) dat(1:isize(iarr),icolumn,j)
            endif
         enddo
!        set masses for equal mass particles (not dumped in small dump)
         if (smalldump .and. iarr.eq.1 .and. nreal(1).lt.4) then
            if (abs(pmassinitial).gt.tiny(pmassinitial)) then
               icolumn = 4
               dat(1:isize(iarr),icolumn,j) = pmassinitial
            endif
         endif
!        real4's go straight into dat
         imaxcolumnread = max(icolumn,6)
         do i=1,nreal4(iarr)
            if (iarr.eq.1 .and. i.eq.1) then
               icolumn = 6 ! density
            elseif (iarr.eq.1 .and. smalldump .and. i.eq.2) then
               icolumn = 5 ! h which is real4 in small dumps
            else
               icolumn = imaxcolumnread + 1
               imaxcolumnread = icolumn
            endif
            read(iunit,end=33,iostat=ierr) dat(1:isize(iarr),icolumn,j) 
         enddo
         icolumn = imaxcolumnread
!        real 8's need converting
         do i=1,nreal8(iarr)
            icolumn = icolumn + 1 !!nextcolumn(icolumn,iarr,nhydroarrays,nmhdarrays,imaxcolumnread) 
            read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
            dat(1:isize(iarr),icolumn,j) = real(dattemp(1:isize(iarr)))
         enddo
      endif
   enddo
!
!--reached end of file (during data read)
!
33 continue

    nptmassi = 0
    nunknown = 0
!
!--place point masses after normal particles
!     
    if (any(iphase(1:ntotal).ne.0)) then
       nunknown = 0
       do i=1,npart
          if (iphase(i).ne.0) nunknown = nunknown + 1
       enddo
       allocate(dattemp2(nunknown,ncolstep))

       nptmassi = 0
       ipos = 0
       do i=1,ntotal
          ipos = ipos + 1
          if (iphase(i).ge.1) then
             nptmassi = nptmassi + 1
             !--save point mass information in temporary array
             if (nptmassi.gt.size(dattemp2(:,1))) stop 'error: ptmass array bounds exceeded in data read'
             dattemp2(nptmassi,1:ncolstep) = dat(i,1:ncolstep,j)
!             print*,i,' removed', dat(i,1:3,j)
             ipos = ipos - 1
          endif
         !--shuffle dat array
          if (ipos.ne.i .and. i.lt.ntotal) then
  !           print*,'copying ',i+1,'->',ipos+1
             dat(ipos+1,1:ncolstep,j) = dat(i+1,1:ncolstep,j)
             !--must also shuffle iphase (to be correct for other types)
             iphase(ipos+1) = iphase(i+1)
          endif
       enddo
       if (nptmassi.ne.nptmass) print *,'WARNING: nptmass from iphase =',nptmassi,'not equal to nptmass =',nptmass
       !--append ptmasses to end of dat array
       do i=1,nptmassi
          ipos = ipos + 1             
!          print*,ipos,' appended', dattemp2(i,1:3)
          dat(ipos,1:ncolstep,j) = dattemp2(i,1:ncolstep)
          !--we make iphase = 1 for point masses (could save iphase and copy across but no reason to)
          iphase(ipos) = 1
       enddo
!
!--do the same with unknown/dead particles
!     
       nunknown = 0
       ipos = 0
       do i=1,ntotal
          ipos = ipos + 1
          if (iphase(i).lt.0) then
             nunknown = nunknown + 1
             !--save information in temporary array
             if (nunknown.gt.size(dattemp2(:,1))) stop 'error: array bounds for dead particles exceeded in data read'
             dattemp2(nunknown,1:ncolstep) = dat(i,1:ncolstep,j)
  !           print*,i,' removed'
             ipos = ipos - 1
          endif
          !--shuffle dat array
          if (ipos.ne.i .and. i.lt.ntotal) then
             !print*,'copying ',i+1,'->',ipos+1
             dat(ipos+1,1:ncolstep,j) = dat(i+1,1:ncolstep,j)
             iphase(ipos+1) = iphase(i+1) ! for completeness (ie. if more types used in future)
          endif
       enddo
       !--append dead particles to end of dat array
       do i=1,nunknown
          ipos = ipos + 1
          dat(ipos,1:ncolstep,j) = dattemp2(i,1:ncolstep)
       enddo

     endif

     if (allocated(dattemp)) deallocate(dattemp)
     if (allocated(dattemp2)) deallocate(dattemp2)
     if (allocated(iphase)) deallocate(iphase)

     npartoftype(1,j) = npart - nptmassi - nunknown
     npartoftype(2,j) = ntotal - npart
     npartoftype(3,j) = nptmassi
     npartoftype(4,j) = nunknown
     if (npartoftype(2,j).ne.0) then
        print*,' n(gas) = ',npartoftype(1,j),' n(ghost) = ',npartoftype(2,j),' n(sinks) = ',nptmassi, ' n(unknown) = ',nunknown
     else
        print*,' n(gas) = ',npartoftype(1,j),' n(sinks) = ',nptmassi, ' n(unknown) = ',nunknown
     endif
     
     close(15)
     
     return

55 continue
   print "(a)", ' *** ERROR: end of file during header read ***'

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
!--all formats read the following columns    
  do i=1,ndim
     ix(i) = i
  enddo
  ipmass = 4   !  particle mass
  ih = 5       !  smoothing length
  irho = 6     !  density
!--the following only for mhd small dumps or full dumps
  if (ncolumns.ge.7) then
  select case(iformat)
     case(0) ! hydro full dump
        ivx = 7
        iutherm = 10
        label(11) = 'grad h'
        label(12) = 'grad soft'
     case(2) ! mhd full dump
        iBfirst = 7
        ivx = 10
        iutherm = 13
        label(14) = 'grad h'
        label(15) = 'grad soft'
        if (ncolumns.ge.21) then
           label(16) = 'Euler alpha'
           label(17) = 'Euler beta'
           idivB = 18
           iJfirst = 19
        elseif (ncolumns.ge.19) then
           idivB = 16
           label(idivB) = 'div B'
           iJfirst = 17
        endif
        if (ncolumns.ge.iJfirst+ndimV) then
           label(iJfirst+ndimV) = 'alpha\dB\u'
        endif
     case(3) ! mhd small dump
        iBfirst = 7
     end select
  endif
  
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  if (irho.gt.0) label(irho) = 'density'
  if (iutherm.gt.0) label(iutherm) = 'u'
  if (ih.gt.0) label(ih) = 'h       '
  if (ipmass.gt.0) label(ipmass) = 'particle mass'     
  if (idivB.gt.0) label(idivB) = 'div B'

  !
  !--set labels for vector quantities
  !
  if (ivx.gt.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
        label(ivx+i-1) = trim(labelvec(ivx))//'\d'//labelcoord(i,1)
     enddo
  endif
  if (iBfirst.gt.0) then
     iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
     labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
  endif
  if (iJfirst.gt.0) then
     iamvec(iJfirst:iJfirst+ndimV-1) = iJfirst
     labelvec(iJfirst:iJfirst+ndimV-1) = 'J'  
  endif
  !
  !--set labels for each particle type
  !
  ntypes = 4  !!maxparttypes
  labeltype(1) = 'gas'
  labeltype(2) = 'ghost'
  labeltype(3) = 'sink'
  labeltype(4) = 'unknown/dead'
  UseTypeInRenderings(1) = .true.
  UseTypeInRenderings(2) = .true.
  UseTypeInRenderings(3) = .false.
  UseTypeInRenderings(4) = .true.  ! only applies if turned on

!-----------------------------------------------------------

  return 
end subroutine set_labels
