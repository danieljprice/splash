!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR READING UNFORMATTED OUTPUT FROM
! THE NEXT GENERATION SPH CODE (sphNG)
!
! (also my Phantom SPH code which uses a similar format)
!
! *** CONVERTS TO SINGLE PRECISION ***
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
! SSPLASH_RESET_CM if 'YES' then centre of mass is reset to origin
! SSPLASH_OMEGA if non-zero subtracts corotating velocities with omega as set
! SSPLASH_OMEGAT if non-zero subtracts corotating positions and velocities with omega as set
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
!
! Partial data read implemented Nov 2006 means that columns with
! the 'required' flag set to false are not read (read is therefore much faster)
!-------------------------------------------------------------------------
module sphNGread
 use params
 implicit none
 real(doub_prec) :: udist,umass,utime,umagfd
 real :: tfreefall
 integer :: istartmhd,istartrt,nmhd,idivvcol,nhydroreal4
 logical :: phantomdump,smalldump,mhddump,rtdump,usingvecp,igotmass
 
end module sphNGread

subroutine read_data(rootname,indexstart,nstepsread)
  use particle_data, only:dat,gamma,time,iamtype,npartoftype,maxpart,maxstep,maxcol,icolourme,masstype
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc,required,ipartialread,&
                     lowmemorymode,ntypes
  use mem_allocation, only:alloc
  use system_utils, only:lenvironment,renvironment
  use labels, only:ipmass,irho,ih,ix,ivx
  use calcquantities, only:calc_quantities
  use sphNGread
  implicit none
  integer, intent(in) :: indexstart
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  integer, parameter :: maxarrsizes = 10, maxreal = 50
  real, parameter :: pi=3.141592653589
  integer :: i,j,ierr,iunit
  integer :: intg1,int2,int3
  integer :: i1,iarr,i2,iptmass1,iptmass2,ilocpmassinitial
  integer :: npart_max,nstep_max,ncolstep,icolumn,nptmasstot
  integer :: narrsizes,nints,nreals,nreal4s,nreal8s
  integer :: nskip,ntotal,npart,n1,n2,ninttypes,ngas
  integer :: nreassign,naccrete,nkill,iblock,nblocks,ntotblock,ncolcopy
  integer :: ipos,nptmass,nptmassi,nstar,nunknown,isink,ilastrequired
  integer :: nhydroarrays,nmhdarrays,imaxcolumnread,nhydroarraysinfile
  integer :: itype,iphaseminthistype,iphasemaxthistype,nthistype
  integer, dimension(maxparttypes) :: npartoftypei
  real, dimension(maxparttypes) :: massoftypei
  logical :: iexist, doubleprec,imadepmasscolumn
  logical :: debug

  character(len=len(rootname)+10) :: dumpfile
  character(len=100) :: fileident
  
  integer*8, dimension(maxarrsizes) :: isize
  integer, dimension(maxarrsizes) :: nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8
  integer*1, dimension(:), allocatable :: iphase
  integer, dimension(:), allocatable :: listpm
  real(doub_prec), dimension(:), allocatable :: dattemp
  real(doub_prec) :: r8
  real, dimension(maxreal) :: dummyreal
  real, dimension(:,:), allocatable :: dattemp2
  real, dimension(3) :: xyzsink
  real :: rhozero,hfact,omega,r4,tff

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart
  npart = 0
  iunit = 15
  ipmass = 4
  idivvcol = 0
  nhydroreal4 = 0

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
  doubleprec = .true.
  
  write(*,"(26('>'),1x,a,1x,26('<'))") trim(dumpfile)

  debug = lenvironment('SSPLASH_DEBUG')
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
      read(iunit,iostat=ierr) intg1,r8,int2,i1,int3
      if (intg1.ne.690706) then
         print "(a)",'*** ERROR READING HEADER: corrupt file/zero size/wrong endian?'
         close(iunit)
         return
      endif
      if (int2.ne.780806) then
         print "(a)",'single precision dump'
         rewind(iunit)
         read(iunit,iostat=ierr) intg1,r4,int2,i1,int3
         if (int2.ne.780806) then
            print "(a)",'ERROR determining single/double precision in file header'
         endif
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
   smalldump = .false.
   mhddump = .false.
   usingvecp = .false.
   rtdump = .false.
   imadepmasscolumn = .false.
   if (fileident(1:1).eq.'S') then
      smalldump = .true.
   endif
   if (index(fileident,'Phantom').ne.0) then
      phantomdump = .true.
   else
      phantomdump = .false.
   endif
   if (index(fileident,'vecp').ne.0) then
      usingvecp = .true.
   endif
!
!--read global dump header
!   
   nblocks = 1 ! number of MPI blocks
   npartoftypei(:) = 0
   read(iunit,iostat=ierr) nints
   if (ierr /=0) then
      print "(a)",'error reading nints'
      close(iunit)
      return
   else
      if (nints.lt.3) then
         if (.not.phantomdump) print "(a)",'WARNING: npart,n1,n2 NOT IN HEADER??'
         read(iunit,iostat=ierr) npart
         npartoftypei(1) = npart
      elseif (phantomdump) then
         if (nints.lt.7) then
            ntypes = nints - 1
            read(iunit,iostat=ierr) npart,npartoftypei(1:ntypes)
         else
            ntypes = 1
            read(iunit,iostat=ierr) npart,npartoftypei(1:5),nblocks
         endif
         n1 = npartoftypei(1)
         n2 = 0
      elseif (nints.ge.7) then
         read(iunit,iostat=ierr) npart,n1,n2,nreassign,naccrete,nkill,nblocks
      else
         print "(a)",'warning: nblocks not read from file (assuming non-MPI dump)'         
         read(iunit,iostat=ierr) npart,n1,n2
      endif
      if (ierr /=0) then
         print "(a)",'error reading npart,n1,n2 and/or number of MPI blocks'
         close(iunit)
         return
      elseif (nblocks.gt.2000) then
         print *,'npart = ',npart,' MPI blocks = ',nblocks
         nblocks = 1
         print*,' corrupt number of MPI blocks, assuming 1 '
      else
         print *,'npart = ',npart,' MPI blocks = ',nblocks
      endif
   endif
!--int*1, int*2, int*4, int*8
   do i=1,4
      read(iunit,end=55,iostat=ierr) ninttypes
      if (ninttypes.gt.0) read(iunit,end=55,iostat=ierr)
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
      umagfd = 1.0
   else
      print "(a)",'*** WARNING: units not found in file'
      udist = 1.0
      umass = 1.0
      utime = 1.0
      umagfd = 1.0
   endif
   if (ierr /= 0) then
      print "(a)",'*** error reading units'
   endif
!
!--Total number of array blocks in the file
!
   read(iunit,end=55,iostat=ierr) narrsizes
   narrsizes = narrsizes/nblocks
   if (ierr /= 0) then 
      print "(a)",'*** error reading number of array sizes ***'
      close(iunit)
      return
   elseif (narrsizes.gt.maxarrsizes) then
      narrsizes = maxarrsizes
      print "(a,i2)",'WARNING: too many array sizes: reading only ',narrsizes
   endif
   if (narrsizes.ge.4 .and. nreal8s.lt.4) then
      print "(a)",' WARNING: could not read magnetic units from dump file'
   endif
   print*,' number of array sizes = ',narrsizes
!
!--Attempt to read all MPI blocks
!
   ntotal = 0
   ntotblock = 0
   nptmasstot = 0
   i2 = 0
   iptmass2 = 0
   igotmass = .true.
   massoftypei(:) = 0.

   over_MPIblocks: do iblock=1,nblocks
      
      !if (nblocks.gt.1) print "(10('-'),' MPI block ',i4,1x,10('-'))",iblock
!
!--read array header from this block
!  
   if (iblock.eq.1) ncolstep = 0
   do iarr=1,narrsizes
      read(iunit,end=55,iostat=ierr) isize(iarr),nint(iarr),nint1(iarr),nint2(iarr), &
                 nint4(iarr),nint8(iarr),nreal(iarr),nreal4(iarr),nreal8(iarr)
      if (iarr.eq.1) then
         ntotblock = isize(iarr)
         if (npart.le.0) npart = ntotblock
         ntotal = ntotal + ntotblock
      elseif (iarr.eq.2) then
         nptmasstot = nptmasstot + isize(iarr)
      endif
      
      if (isize(iarr).gt.0 .and. iblock.eq.1) then
         print "(1x,a,i1,a,i12,a,5(i2,1x),a,3(i2,1x))", &
            'block ',iarr,' dim = ',isize(iarr),' nint =',nint(iarr),nint1(iarr), &
            nint2(iarr),nint4(iarr),nint8(iarr),'nreal =',nreal(iarr),nreal4(iarr),nreal8(iarr)
      endif
!--we are going to read all real arrays but need to convert them all to default real
      if (isize(iarr).eq.isize(1) .and. iblock.eq.1) then
         ncolstep = ncolstep + nreal(iarr) + nreal4(iarr) + nreal8(iarr)
      endif
   enddo 
!
!--this is a bug fix for a corrupt version of wdump outputting bad
!  small dump files
!
   if (smalldump .and. nreal(1).eq.5 .and. iblock.eq.1) then
      print*,'FIXING CORRUPT HEADER ON SMALL DUMPS: assuming nreal=3 not 5'
      nreal(1) = 3
      ncolstep = ncolstep - 2
   endif
   
   npart_max = maxval(isize(1:narrsizes))
   npart_max = max(npart_max,npart,ntotal)
!
!--work out from array header what sort of dump this is and where things should lie
!
   if (iblock.eq.1) then
      igotmass = .true.
      if (smalldump .or. phantomdump) then
         if (phantomdump) then
            ilocpmassinitial = 15
         else
            ilocpmassinitial = 23
         endif
         if (nreals.ge.ilocpmassinitial) then
            massoftypei(1) = dummyreal(ilocpmassinitial)
            if (massoftypei(1).gt.tiny(0.) .and. .not.lowmemorymode) then
               ncolstep = ncolstep + 1  ! make an extra column to contain particle mass
               imadepmasscolumn = .true.
            elseif (lowmemorymode) then
               igotmass = .false.
            endif
            !--read dust mass from phantom dumps
            if (phantomdump .and. nreals.ge.ilocpmassinitial+1) then
               massoftypei(2) = dummyreal(ilocpmassinitial+1)
            else
               massoftypei(2) = 0.
            endif
         else
            print "(a)",' error extracting particle mass from small dump file'
            massoftypei(1) = 0.
            igotmass = .false.
         endif
         if (abs(massoftypei(1)).lt.tiny(0.) .and. nreal(1).lt.4) then
            print "(a)",' error: particle masses not present in small dump file'
            igotmass = .false.
         endif
      endif
!
!--   to handle both small and full dumps, we need to place the quantities dumped
!     in both small and full dumps at the start of the dat array
!     quantities only in the full dump then come after
!     also means that hydro/MHD are "semi-compatible" in the sense that x,y,z,m,h and rho
!     are in the same place for both types of dump
!
      ix(1) = 1
      ix(2) = 2
      ix(3) = 3
      if (igotmass) then  
         ipmass = 4
         ih = 5
         irho = 6
         nhydroarrays = 6 ! x,y,z,m,h,rho
      else
         ipmass = 0
         ih = 4
         irho = 5
         nhydroarrays = 5 ! x,y,z,h,rho      
      endif
      nhydroarraysinfile = nreal(1) + nreal4(1) + nreal8(1)
      nhydroreal4 = nreal4(1)
      if (imadepmasscolumn) nhydroarraysinfile = nhydroarraysinfile + 1
      if (nhydroarraysinfile .lt.nhydroarrays .and. .not.phantomdump) then
         print "(a)",' ERROR: one of x,y,z,m,h or rho missing in small dump read'
         nhydroarrays = nreal(1)+nreal4(1)+nreal8(1)
      elseif (phantomdump .and. (nreal(1).lt.3 .or. nreal4(1).lt.1)) then
         print "(a)",' ERROR: x,y,z or h missing in phantom read'
      endif
      if (narrsizes.ge.4) then
         nmhdarrays = 3 ! Bx,By,Bz
         nmhd = nreal(4) + nreal4(4) + nreal8(4) - nmhdarrays ! how many "extra" mhd arrays
      else
         nmhdarrays = 0
      endif
 
      !--radiative transfer dump?
      if (narrsizes.ge.3 .and. isize(3).eq.isize(1)) rtdump = .true.
      !--mhd dump?
      if (narrsizes.ge.4) mhddump = .true.

      if (.not.(mhddump.or.smalldump)) then
         ivx = nhydroarrays+1
      elseif (mhddump .and. .not.smalldump) then
         ivx = nhydroarrays+nmhdarrays+1
      else
         ivx = 0
      endif
      !--need to force read of velocities e.g. for corotating frame subtraction
      if (any(required(ivx:ivx+ndimV-1))) required(ivx:ivx+ndimV-1) = .true.
      
      !--for phantom dumps, also make a column for density
      !  and divv, if a .divv file exists
      if (phantomdump) then
         ncolstep = ncolstep + 1
         inquire(file=trim(dumpfile)//'.divv',exist=iexist)
         if (iexist) then
            ncolstep = ncolstep + 1
            idivvcol = ncolstep
         endif
      endif
   endif
!
!--allocate memory now that we know the number of columns
!
   if (iblock.eq.1) ncolumns = ncolstep + ncalc
   if (npart_max.gt.maxpart .or. j.gt.maxstep .or. ncolumns.gt.maxcol) then
      if (lowmemorymode) then
         ilastrequired = 0
         do i=1,ncolumns
            if (required(i)) ilastrequired = i
         enddo
         call alloc(max(npart_max,maxpart),j,ilastrequired)
      else
         call alloc(max(npart_max,maxpart),j,ncolumns,mixedtypes=.true.)
      endif
   endif
   
   if (iblock.eq.1) then
!--extract required information from the first block header
      time(j) = dummyreal(1)
      gamma(j) = dummyreal(3)
      rhozero = dummyreal(4)
      masstype(:,j) = massoftypei(:)

      if (rhozero.gt.0.) then
         tfreefall = SQRT((3. * pi) / (32. * rhozero))
         tff = time(j)/tfreefall
      else
         tfreefall = 0.
         tff = 0.
      endif
      if (phantomdump) then
         npartoftype(1:ntypes,j) = npartoftypei(1:ntypes)
         if (nblocks.gt.1) then
            print "(a)",' setting ngas=npart for MPI code '
            npartoftype(1,j) = npart
            npartoftype(2:,j) = 0
         endif
      else
         npartoftype(1,j) = npart
         npartoftype(2,j) = max(ntotal - npart,0)
      endif
      hfact = 1.2
      if (phantomdump) then
         if (nreals.lt.6) then
            print "(a)",' error: hfact not present in phantom dump'
         else
            hfact = dummyreal(6)
         endif
         print "(a,1pe12.4,a,0pf6.3,a,0pf5.2,a,1pe7.1)", &
               ' time = ',time(j),' gamma = ',gamma(j), &
               ' hfact = ',hfact,' tolh = ',dummyreal(7)
      else
         print "(a,1pe12.4,a,0pf9.5,a,f8.4,/,a,1pe12.4,a,1pe9.2,a,1pe10.2)", &
               '   time: ',time(j),  '   gamma: ',gamma(j), '   RK2: ',dummyreal(5), &
               ' t/t_ff: ',tff,' rhozero: ',rhozero,' dtmax: ',dummyreal(2)
      endif
      nstepsread = nstepsread + 1

      if (allocated(iphase)) deallocate(iphase)
      allocate(iphase(npart_max))
      iphase(:) = 0
   endif ! iblock = 1
!
!--Arrays
!
   imaxcolumnread = 0
   icolumn = 0
   istartmhd = 0
   istartrt = 0
   i1 = i2 + 1
   i2 = i1 + isize(1) - 1
   print "(1x,a10,i4,3(a,i12))",'MPI block ',iblock,':  particles: ',i1,' to ',i2,' of ',npart
   iptmass1 = iptmass2 + 1
   iptmass2 = iptmass1 + isize(2) - 1
   nptmass = nptmasstot
   if (nptmass.gt.0) print "(15x,3(a,i12))",'  pt. masses: ',iptmass1,' to ',iptmass2,' of ',nptmass

   do iarr=1,narrsizes
      if (nreal(iarr) + nreal4(iarr) + nreal8(iarr).gt.0) then
         if (iarr.eq.4) then
            istartmhd = imaxcolumnread + 1
            if (debug) print*,' istartmhd = ',istartmhd
         elseif (iarr.eq.3) then
            istartrt = max(nhydroarrays+nmhdarrays+1,imaxcolumnread + 1)
            if (debug) print*,' istartrt = ',istartrt
         endif
      endif 
!--read iphase from array block 1
      if (iarr.eq.1) then
         !--skip default int
         nskip = nint(iarr)
         do i=1,nskip
            read(iunit,end=33,iostat=ierr)
         enddo
         if (nint1(iarr).lt.1) then
            if (.not.phantomdump) print "(a)",'ERROR: can''t locate iphase in dump'
            !--skip remaining integer arrays
            nskip = nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
         else
            read(iunit,end=33,iostat=ierr) iphase(i1:i2)
            !--skip remaining integer arrays
            nskip = nint1(iarr) - 1 + nint2(iarr) + nint4(iarr) + nint8(iarr)
         endif
      elseif (smalldump .and. iarr.eq.2 .and. isize(iarr).gt.0) then
!--read listpm from array block 2 for small dumps (needed here to extract sink masses)
         if (allocated(listpm)) deallocate(listpm)
         allocate(listpm(isize(iarr)))
         if (nint(iarr).lt.1) then      
            print "(a)",'ERROR: can''t locate listpm in dump'
            nskip = nint(iarr) + nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
         else
            read(iunit,end=33,iostat=ierr) listpm(1:isize(iarr))
            nskip = nint(iarr) - 1 + nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
         endif
      else
!--otherwise skip all integer arrays (not needed for plotting)
         nskip = nint(iarr) + nint1(iarr) + nint2(iarr) + nint4(iarr) + nint8(iarr)
      endif
      
      if (iarr.eq.3 .and. lenvironment('SSPLASH_BEN_HACKED')) then
         nskip = nskip - 1
         print*,' FIXING HACKED DUMP FILE'
      endif
      !print*,'skipping ',nskip
      do i=1,nskip
         read(iunit,end=33,iostat=ierr)
      enddo
!      
!--real arrays
!
      if (isize(iarr).ne.isize(1)) then
         if (smalldump .and. iarr.eq.2 .and. allocated(listpm)) then
!--read sink particle masses from block 2 for small dumps
            if (nreal(iarr).lt.1) then
               if (isize(iarr).gt.0) print "(a)",'ERROR: sink masses not present in small dump'
               nskip = nreal(iarr) + nreal4(iarr) + nreal8(iarr)
            else
               if (doubleprec) then
                  !--convert default real to single precision where necessart
                  if (allocated(dattemp)) deallocate(dattemp)
                  allocate(dattemp(isize(iarr)),stat=ierr)
                  if (ierr /=0) print "(a)",'ERROR in memory allocation'
                  read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
                  if (nptmass.ne.isize(iarr)) print "(a)",'ERROR: nptmass.ne.block size'
                  if (ipmass.gt.0) then
                     do i=1,isize(iarr)
                        dat(listpm(iptmass1+i-1),ipmass,j) = real(dattemp(i))
                     enddo
                  else
                     print*,'WARNING: sink particle masses not read because no mass array allocated' 
                  endif
               else
                  read(iunit,end=33,iostat=ierr) (dat(listpm(i),ipmass,j),i=iptmass1,iptmass2)
               endif 
               nskip = nreal(iarr) - 1 + nreal4(iarr) + nreal8(iarr)
            endif
         else
!--for other blocks, skip real arrays if size different
            nskip = nreal(iarr) + nreal4(iarr) + nreal8(iarr)
         endif
         do i=1,nskip
            read(iunit,end=33,iostat=ierr)
         enddo
      else
!--otherwise read them      
         if ((doubleprec.and.nreal(iarr).gt.0).or.nreal8(iarr).gt.0) then
            if (allocated(dattemp)) deallocate(dattemp)
            allocate(dattemp(isize(iarr)),stat=ierr)
            if (ierr /=0) print "(a)",'ERROR in memory allocation (read_data_sphNG: dattemp)'
         endif

!        default reals may need converting
         do i=1,nreal(iarr)
            if (iarr.eq.1.and.((phantomdump.and.i.eq.4) &
               .or.(.not.phantomdump.and.i.eq.6))) then
               ! read x,y,z,m,h and then place arrays after always-present ones
               ! (for phantom read x,y,z only)
               icolumn = nhydroarrays+nmhdarrays + 1
            elseif (.not.phantomdump .and. (iarr.eq.4 .and. i.le.3)) then
               icolumn = nhydroarrays + i
            else
               icolumn = imaxcolumnread + 1
            endif
            imaxcolumnread = max(imaxcolumnread,icolumn)
            if (debug) print*,' reading real ',icolumn
            if (required(icolumn)) then
               if (doubleprec) then
                  read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
                  dat(i1:i2,icolumn,j) = real(dattemp(1:isize(iarr)))
               else
                  read(iunit,end=33,iostat=ierr) dat(i1:i2,icolumn,j)
               endif
            else
               read(iunit,end=33,iostat=ierr)
            endif
         enddo
!        set masses for equal mass particles (not dumped in small dump)
         if (((smalldump.and.nreal(1).lt.ipmass).or.phantomdump).and. iarr.eq.1) then
            if (abs(massoftypei(1)).gt.tiny(massoftypei)) then
               icolumn = ipmass
               if (required(ipmass) .and. ipmass.gt.0) then
                  where (iphase(i1:i2).eq.0) dat(i1:i2,icolumn,j) = massoftypei(1)
               endif
               !--dust mass for phantom particles
               if (phantomdump .and. npartoftypei(2).gt.0 .and. ipmass.gt.0) then
                  print *,' dust particle mass = ',massoftypei(2),' ratio dust/gas = ',massoftypei(2)/massoftypei(1)
                  dat(npartoftypei(1)+1:npartoftypei(1)+npartoftypei(2),icolumn,j) = massoftypei(2)
               endif
               if (debug) print*,'mass ',icolumn
            elseif (phantomdump) then
               print*,' ERROR: particle mass zero in Phantom dump file!'
            endif
         endif
!        real4's go straight into dat
         imaxcolumnread = max(imaxcolumnread,icolumn,6)
         do i=1,nreal4(iarr)
            if (phantomdump) then
               if (iarr.eq.1 .and. i.eq.1) then
                  icolumn = ih ! h is always first real4 in phantom dumps
                  !--density depends on h being read
                  if (required(irho)) required(ih) = .true.
               elseif (iarr.eq.4 .and. i.le.3) then
                  icolumn = nhydroarrays + i
               else
                  icolumn = max(nhydroarrays+nmhdarrays + 1,imaxcolumnread + 1)
               endif
            else
               if (iarr.eq.1 .and. i.eq.1) then
                  icolumn = irho ! density
               elseif (iarr.eq.1 .and. smalldump .and. i.eq.2) then
                  icolumn = ih ! h which is real4 in small dumps
               elseif (iarr.eq.4 .and. i.le.3) then
                  icolumn = nhydroarrays + i
               else
                  icolumn = max(nhydroarrays+nmhdarrays + 1,imaxcolumnread + 1)
               endif
            endif
            imaxcolumnread = max(imaxcolumnread,icolumn)
            if (debug) print*,'reading real4 ',icolumn
            if (required(icolumn)) then
               read(iunit,end=33,iostat=ierr) dat(i1:i2,icolumn,j)
            else
               read(iunit,end=33,iostat=ierr)
            endif
            !--construct density for phantom dumps based on h, hfact and particle mass
            if (phantomdump .and. icolumn.eq.ih) then
               icolumn = irho ! density
               !
               !--dead particles have -ve smoothing lengths in phantom
               !  so use abs(h) for these particles and hide them
               !
               if (required(irho)) then
                  where (dat(i1:i2,ih,j).gt.0.)
                     dat(i1:i2,irho,j) = &
                        massoftypei(1)*(hfact/dat(i1:i2,ih,j))**3
                  elsewhere (dat(i1:i2,ih,j).lt.0.)
                     dat(i1:i2,irho,j) = &
                        massoftypei(1)*(hfact/abs(dat(i1:i2,ih,j)))**3
                     icolourme(i1:i2) = -1
                  elsewhere ! if h = 0.
                     dat(i1:i2,irho,j) = 0.
                     icolourme(i1:i2) = -1
                  end where
               endif
                   
               if (debug) print*,'making density ',icolumn
            endif
         enddo
         icolumn = imaxcolumnread
!        real 8's need converting
         do i=1,nreal8(iarr)
            icolumn = icolumn + 1 !!nextcolumn(icolumn,iarr,nhydroarrays,nmhdarrays,imaxcolumnread) 
            if (required(icolumn)) then
               read(iunit,end=33,iostat=ierr) dattemp(1:isize(iarr))
               dat(i1:i2,icolumn,j) = real(dattemp(1:isize(iarr)))
            else
               read(iunit,end=33,iostat=ierr)
            endif
         enddo
      endif
   enddo ! over array sizes
   enddo over_MPIblocks
!
!--reached end of file (during data read)
!
   goto 34
33 continue
   print "(1x,a)",'WARNING: end of file during read'
34 continue
 !
 !--read .divv file for phantom dumps
 !
    if (phantomdump .and. idivvcol.ne.0 .and. required(idivvcol)) then
       print "(a)",' reading divv from '//trim(dumpfile)//'.divv'
       open(unit=66,file=trim(dumpfile)//'.divv',form='unformatted',status='old',iostat=ierr)
       if (ierr /= 0) then
          print "(a)",' ERROR opening '//trim(dumpfile)//'.divv'
       else
          read(66,iostat=ierr) dat(1:ntotal,idivvcol,j)
          if (ierr /= 0) print "(a)",' WARNING: ERRORS reading file'
          close(66)
       endif
    endif
 !
 !--reset centre of mass to zero if environment variable "SSPLASH_RESET_CM" is set
 !
    if (allocated(dat) .and. n1.GT.0 .and. lenvironment('SSPLASH_RESET_CM')) then
       call reset_centre_of_mass(dat(1:n1,1:3,j),dat(1:n1,4,j),iphase(1:n1),n1)
    endif
 !
 !--centre on sink if "SSPLASH_CENTRE_ON_SINK" is set
 !
    if (lenvironment('SSPLASH_CENTRE_ON_SINK')) then
       if (nptmass.EQ.1) then
          isink = 0
          xyzsink = 0.
          do i=1,ntotal
             if (iphase(i).GE.1) then
                isink = i
                xyzsink(1:3) = dat(isink,1:3,j)
             endif
          enddo
          if (isink.EQ.0) then
             print "(a)",'WARNING: SSPLASH_CENTRE_ON_SINK set but cannot find sink'
          else
             print "(a,3(1pe10.3,1x))",' CENTREING ON SINK PARTICLE ',PACK(xyzsink(1:3),required(1:3))
             do i=1,3
                dat(1:ntotal,i,j) = dat(1:ntotal,i,j) - xyzsink(i)
             enddo
          endif
       endif
    endif
 !
 !--reset corotating frame velocities if environment variable "SSPLASH_OMEGA" is set
 !
    if (allocated(dat) .and. n1.GT.0 .and. all(required(1:2))) then
       omega = renvironment('SSPLASH_OMEGAT')
       if (abs(omega).gt.tiny(omega) .and. ndim.ge.2) then
          call reset_corotating_positions(n1,dat(1:n1,1:2,j),omega,time(j))
       endif

       if (.not. smalldump) then
          if (abs(omega).lt.tiny(omega)) omega = renvironment('SSPLASH_OMEGA')
          if (abs(omega).gt.tiny(omega) .and. ivx.gt.0) then
             if (.not.all(required(1:2)) .or. .not.all(required(ivx:ivx+1))) then
                print*,' ERROR subtracting corotating frame with partial data read'
             else
                call reset_corotating_velocities(n1,dat(1:n1,1:2,j),dat(1:n1,ivx:ivx+1,j),omega)
             endif
          endif
       endif
    endif

    !--set flag to indicate that only part of this file has been read 
    if (.not.all(required(1:ncolstep))) ipartialread = .true.
    
    
    nptmassi = 0
    nunknown = 0
    ngas = 0
    nstar = 0
!
!--translate iphase into particle types (mixed type storage)
!
    if (size(iamtype(:,j)).gt.1) then
       do i=1,npart
          select case(int(iphase(i)))
          case(0)
            iamtype(i,j) = 1
            ngas = ngas + 1
          case(1:9)
            iamtype(i,j) = 3
            nptmassi = nptmassi + 1
          case(10:)
            iamtype(i,j) = 4
            nstar = nstar + 1
          case default
            iamtype(i,j) = 5
            nunknown = nunknown + 1
          end select
       enddo
       do i=npart+1,ntotal
          iamtype(i,j) = 2
       enddo 
       !print*,'mixed types: ngas = ',ngas,nptmassi,nunknown

    elseif (any(iphase(1:ntotal).ne.0)) then
!
!--place point masses after normal particles
!  if not storing the iamtype array
!     
       print "(a)",' sorting particles by type...'
       nunknown = 0
       do i=1,npart
          if (iphase(i).ne.0) nunknown = nunknown + 1
       enddo
       ncolcopy = min(ncolstep,maxcol)
       allocate(dattemp2(nunknown,ncolcopy))

       do itype=1,3       
          nthistype = 0
          ipos = 0
          select case(itype)
          case(1) ! ptmass
             iphaseminthistype = 1
             iphasemaxthistype = 9
          case(2) ! star
             iphaseminthistype = 10
             iphasemaxthistype = huge(iphasemaxthistype)
          case(3) ! unknown
             iphaseminthistype = -huge(iphaseminthistype)
             iphasemaxthistype = -1
          end select

          do i=1,ntotal
             ipos = ipos + 1
             if (iphase(i).ge.iphaseminthistype .and. iphase(i).le.iphasemaxthistype) then
                nthistype = nthistype + 1
                !--save point mass information in temporary array
                if (nptmassi.gt.size(dattemp2(:,1))) stop 'error: ptmass array bounds exceeded in data read'
                dattemp2(nthistype,1:ncolcopy) = dat(i,1:ncolcopy,j)
   !             print*,i,' removed', dat(i,1:3,j)
                ipos = ipos - 1
             endif
            !--shuffle dat array
             if (ipos.ne.i .and. i.lt.ntotal) then
     !           print*,'copying ',i+1,'->',ipos+1
                dat(ipos+1,1:ncolcopy,j) = dat(i+1,1:ncolcopy,j)
                !--must also shuffle iphase (to be correct for other types)
                iphase(ipos+1) = iphase(i+1)
             endif
          enddo

          !--append this type to end of dat array
          do i=1,nthistype
             ipos = ipos + 1             
   !          print*,ipos,' appended', dattemp2(i,1:3)
             dat(ipos,1:ncolcopy,j) = dattemp2(i,1:ncolcopy)
             !--we make iphase = 1 for point masses (could save iphase and copy across but no reason to)
             iphase(ipos) = iphaseminthistype
          enddo
          
          select case(itype)
          case(1)
             nptmassi = nthistype
             if (nptmassi.ne.nptmass) print *,'WARNING: nptmass from iphase =',nptmassi,'not equal to nptmass =',nptmass
          case(2)
             nstar = nthistype
          case(3)
             nunknown = nthistype
          end select          
       enddo

     endif

     if (allocated(dattemp)) deallocate(dattemp)
     if (allocated(dattemp2)) deallocate(dattemp2)
     if (allocated(iphase)) deallocate(iphase)
     if (allocated(listpm)) deallocate(listpm)

     if (.not.phantomdump) then
        npartoftype(1,j) = npart - nptmassi - nstar - nunknown
        npartoftype(2,j) = ntotal - npart
        npartoftype(3,j) = nptmassi
        npartoftype(4,j) = nstar
        npartoftype(5,j) = nunknown
     endif
     
     if (phantomdump) then
        print*,' n(gas) = ',npartoftype(1,j),' n(dust) = ',npartoftype(2,j)
     else
        if (npartoftype(2,j).ne.0) then
           print "(5(a,i10))",' n(gas) = ',npartoftype(1,j),' n(ghost) = ',npartoftype(2,j), &
                  ' n(sinks) = ',nptmassi,' n(stars) = ',nstar,' n(unknown) = ',nunknown        
        else
           print "(5(a,i10))",' n(gas) = ',npartoftype(1,j),' n(sinks) = ',nptmassi, &
                              ' n(stars) = ',nstar,' n(unknown) = ',nunknown
        endif
     endif
     
     close(15)
     
     return

55 continue
   print "(a)", ' *** ERROR: end of file during header read ***'

close(15)
   
return

contains

!
!--reset centre of mass to zero
!
 subroutine reset_centre_of_mass(xyz,pmass,iphase,np)
  implicit none
  integer, intent(in) :: np
  real, dimension(np,3), intent(inout) :: xyz
  real, dimension(np), intent(in) :: pmass
  integer*1, dimension(np), intent(in) :: iphase
  real :: masstot,pmassi
  real, dimension(3) :: xcm
  integer :: i
  
  !
  !--get centre of mass
  !
  xcm(:) = 0.
  masstot = 0.
  do i=1,np
     if (iphase(i).ge.0) then
        pmassi = pmass(i)
        masstot = masstot + pmass(i)
        where (required(1:3)) xcm(:) = xcm(:) + pmassi*xyz(i,:)
     endif
  enddo
  xcm(:) = xcm(:)/masstot
  print*,'RESETTING CENTRE OF MASS (',pack(xcm,required(1:3)),') TO ZERO '
  
  if (required(1)) xyz(1:np,1) = xyz(1:np,1) - xcm(1)
  if (required(2)) xyz(1:np,2) = xyz(1:np,2) - xcm(2)
  if (required(3)) xyz(1:np,3) = xyz(1:np,3) - xcm(3)
  
  return
 end subroutine reset_centre_of_mass

 subroutine reset_corotating_velocities(np,xy,velxy,omeg)
  implicit none
  integer, intent(in) :: np
  real, dimension(np,2), intent(in) :: xy
  real, dimension(np,2), intent(inout) :: velxy
  real, intent(in) :: omeg
  integer :: ip
  
  print*,'SUBTRACTING COROTATING VELOCITIES, OMEGA = ',omeg
  do ip=1,np
     velxy(ip,1) = velxy(ip,1) + xy(ip,2)*omeg
  enddo
  do ip=1,np
     velxy(ip,2) = velxy(ip,2) - xy(ip,1)*omeg
  enddo
  
  return
 end subroutine reset_corotating_velocities

 subroutine reset_corotating_positions(np,xy,omeg,t)
  implicit none
  integer, intent(in) :: np
  real, dimension(np,2), intent(inout) :: xy
  real, intent(in) :: omeg,t
  real :: phii,phinew,r
  integer :: ip
  
  print*,'SUBTRACTING COROTATING POSITIONS, OMEGA = ',omeg,' t = ',t
!$omp parallel default(none) &
!$omp shared(xy,np) &
!$omp firstprivate(omeg,t) &
!$omp private(ip,r,phii,phinew)
!$omp do
  do ip=1,np
     r = sqrt(xy(ip,1)**2 + xy(ip,2)**2)
     phii = atan2(xy(ip,2),xy(ip,1))
     phinew = phii + omeg*t
     xy(ip,1) = r*COS(phinew)
     xy(ip,2) = r*SIN(phinew)
  enddo
!$omp end do
!$omp end parallel

  return
 end subroutine reset_corotating_positions

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels, only:label,labeltype,labelvec,iamvec, &
              ix,ipmass,irho,ih,iutherm,ivx,iBfirst,idivB,iJfirst,icv,iradenergy
  use params
  use settings_data, only:ndim,ndimV,ntypes,ncolumns,UseTypeInRenderings
  use geometry, only:labelcoord
  use settings_units, only:units,unitslabel,unitzintegration,labelzintegration
  use sphNGread
  implicit none
  integer :: i
  real(doub_prec) :: uergg
  
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
  if (igotmass) then
     ipmass = 4   !  particle mass
     ih = 5       !  smoothing length
  else
     ipmass = 0
     ih = 4       !  smoothing length
  endif
  irho = ih + 1     !  density
  if (smalldump .and. nhydroreal4.ge.3) iutherm = irho+1
  
!--the following only for mhd small dumps or full dumps
  if (ncolumns.ge.7) then
     if (mhddump) then
        iBfirst = irho+1
        if (.not.smalldump) then
           ivx = iBfirst+ndimV
           iutherm = ivx+ndimV

           if (phantomdump) then
              !--phantom MHD full dumps
              if (nmhd.ge.4) then
                 iamvec(istartmhd:istartmhd+ndimV-1) = istartmhd
                 labelvec(istartmhd:istartmhd+ndimV-1) = 'A'
                 do i=1,ndimV
                    label(istartmhd+i-1) = trim(labelvec(istartmhd))//'\d'//labelcoord(i,1)
                 enddo
                 idivB = istartmhd+ndimV
              elseif (nmhd.ge.3) then
                 label(istartmhd) = 'Euler alpha'
                 label(istartmhd+1) = 'Euler beta'
                 idivB = istartmhd + 2
              elseif (nmhd.ge.1) then
                 idivB = istartmhd
              endif
              iJfirst = 0
              if (ncolumns.ge.idivB+1) then
                 label(idivB+1) = 'alpha\dB\u'
              endif

           else
              !--sphNG MHD full dumps
              label(iutherm+1) = 'grad h'
              label(iutherm+2) = 'grad soft'
              label(iutherm+3) = 'alpha'
              if (nmhd.ge.7 .and. usingvecp) then
                 iamvec(istartmhd:istartmhd+ndimV-1) = istartmhd
                 labelvec(istartmhd:istartmhd+ndimV-1) = 'A'
                 do i=1,ndimV
                    label(istartmhd+i-1) = trim(labelvec(16))//'\d'//labelcoord(i,1)
                 enddo
                 idivB = istartmhd+ndimV
              elseif (nmhd.ge.6) then
                 label(istartmhd) = 'Euler alpha'
                 label(istartmhd+1) = 'Euler beta'
                 idivB = istartmhd + 2
              elseif (nmhd.ge.1) then
                 idivB = istartmhd
              endif
              iJfirst = idivB + 1
              if (ncolumns.ge.iJfirst+ndimV) then
                 label(iJfirst+ndimV) = 'alpha\dB\u'
              endif
           endif
        else ! mhd small dump
           if (nhydroreal4.ge.3) iutherm = iBfirst+ndimV
        endif
     elseif (.not.smalldump) then
        ! pure hydro full dump
        ivx = irho+1
        iutherm = ivx + ndimV
        if (phantomdump) then
           label(iutherm+1) = 'alpha'
           label(iutherm+2) = 'alphau'
        else
           label(iutherm+1) = 'grad h'
           label(iutherm+2) = 'grad soft'
           label(iutherm+3) = 'alpha'
        endif
     endif 

     if (istartrt.gt.0 .and. istartrt.le.ncolumns) then ! radiative transfer dump
        iradenergy = istartrt
        label(iradenergy) = 'radiation energy'
        uergg = (udist/utime)**2
        units(iradenergy) = uergg

        if (smalldump) then
           icv = istartrt+1
           !--the following lines refer to a format
           !  which was hopefully never used
           !iutherm = istartrt + 1
           !label(iutherm) = 'u'
           !icv = istartrt+2
        else
           label(istartrt+1) = 'opacity'
           units(istartrt+1) = udist**2/umass

           icv = istartrt+2

           label(istartrt+3) = 'lambda'
           units(istartrt+3) = 1.0

           label(istartrt+4) = 'eddington factor'
           units(istartrt+4) = 1.0
        endif

       if (icv.gt.0) then
          label(icv) = 'u/T'
          units(icv) = uergg
       endif
    else
       iradenergy = 0
       icv = 0
    endif
  endif

  
  label(ix(1:ndim)) = labelcoord(1:ndim,1)
  if (irho.gt.0) label(irho) = 'density'
  if (iutherm.gt.0) label(iutherm) = 'u'
  if (ih.gt.0) label(ih) = 'h       '
  if (ipmass.gt.0) label(ipmass) = 'particle mass'     
  if (idivB.gt.0) label(idivB) = 'div B'
  if (idivvcol.gt.0) label(idivvcol) = 'div v'

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
  !--set units for plot data
  !
!   npower = int(log10(udist))
!   udist = udist/10.**npower
!   udistAU = udist/1.495979e13
   units(1:3) = udist
   unitslabel(1:3) = ' [cm]'
!   do i=1,3
!      write(unitslabel(i),"('[ 10\u',i2,'\d cm]')") npower
!   enddo
   if (ipmass.gt.0) then
      units(ipmass) = umass
      unitslabel(ipmass) = ' [g]'
   endif
   units(ih) = udist
   unitslabel(ih) = ' [cm]'
   if (ivx.gt.0) then
      units(ivx:ivx+ndimV-1) = udist/utime
      unitslabel(ivx:ivx+ndimV-1) = ' [cm/s]'
   endif
   if (iutherm.gt.0) then
      units(iutherm) = (udist/utime)**2
      unitslabel(iutherm) = ' [erg/g]'
   endif
   units(irho) = umass/udist**3
   unitslabel(irho) = ' [g/cm\u3\d]'
   if (iBfirst.gt.0) then
      units(iBfirst:iBfirst+ndimV-1) = umagfd
      unitslabel(iBfirst:iBfirst+ndimV-1) = ' [G]'
   endif
   
   !--use the following two lines for time in years
   units(0) = utime/3.1536e7
   unitslabel(0) = ' yrs'
   !--or use these two lines for time in free-fall times
   !units(0) = 1./tfreefall
   !unitslabel(0) = ' '
  
  unitzintegration = udist
  labelzintegration = ' [cm]'
  !
  !--set labels for each particle type
  !
  if (ntypes.eq.2) then  ! phantom
     ntypes = 2
     labeltype(1) = 'gas'
     labeltype(2) = 'dust'
     UseTypeInRenderings(1) = .true.
     UseTypeInRenderings(2) = .false.  
  else
     ntypes = 5
     labeltype(1) = 'gas'
     labeltype(2) = 'ghost'
     labeltype(3) = 'sink'
     labeltype(4) = 'star'
     labeltype(5) = 'unknown/dead'
     UseTypeInRenderings(1) = .true.
     UseTypeInRenderings(2) = .true.
     UseTypeInRenderings(3) = .false.
     UseTypeInRenderings(4) = .true.
     UseTypeInRenderings(5) = .true.  ! only applies if turned on
  endif

!-----------------------------------------------------------

  return 
end subroutine set_labels
