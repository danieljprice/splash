!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR OUTPUT FROM THE GADGET CODE
!
! NOTE THAT THIS ONLY "OFFICIALLY" WORKS WITH THE PARALLEL CODE AS WE
! REQUIRE KNOWLEDGE OF THE PARTICLE SMOOTHING LENGTHS
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
! GSPLASH_FORMAT if = 2 then reads the block-labelled GADGET format
!   rather than the default format.
! GSPLASH_USE_Z if 'YES' uses redshift in the legend instead of time
! GSPLASH_DARKMATTER_HSOFT if given a value > 0.0 will assign a
!  smoothing length to dark matter particles which can then be
!  used in the rendering
! GSPLASH_EXTRACOLS if set to a comma separated list of column labels,
!  will attempt to read additional columns containing gas particle 
!  properties beyond the end of file
! GSPLASH_STARPARTCOLS if set to a comma separated list of column labels,
!  will attempt to read additional columns containing star particle 
!  properties beyond the end of file
! GSPLASH_CHECKIDS if 'YES','yes','TRUE' or 'true' then reads and 
!  checks particle IDS for negative values and flags these as accreted particles
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! dat(maxpart,maxplot,maxstep) : main data array
!
! npartoftype(maxstep): number of particles of each type in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step 
!                      (used in calc_quantities for calculating the pressure)
!
! most of these values are stored in global arrays 
! in the module 'particle_data'
!
! Partial data read implemented Nov 2006 means that columns with
! the 'required' flag set to false are not read (read is therefore much faster)
!-------------------------------------------------------------------------
module gadgetread
 use params, only:maxplot
 implicit none
 real :: hsoft
 character(len=4), dimension(maxplot) :: blocklabelgas
 
end module gadgetread

subroutine read_data(rootname,istepstart,nstepsread)
  use particle_data, only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol,maxstep
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ncalc,iformat,required,ipartialread
  use settings_page, only:legendtext
  use mem_allocation, only:alloc
  use labels, only:ih,irho,ipmass
  use system_utils, only:renvironment,lenvironment,ienvironment,envlist
  use gadgetread, only:hsoft,blocklabelgas
  implicit none
  integer, intent(in) :: istepstart
  integer, intent(out) :: nstepsread
  character(len=*), intent(in) :: rootname
  character(len=len(rootname)+10) :: datfile,densfile,hfile
  character(len=4) :: blocklabel
  integer, dimension(maxparttypes) :: npartoftypei,Nall
  integer, dimension(:), allocatable :: iamtemp
  integer :: i,j,k,itype,icol,ierr,ierrh,ierrrho,nhset,nvec
  integer :: index1,index2,indexstart,indexend,Nmassesdumped
  integer :: ncolstep,npart_max,nstep_max,ntoti,nacc
  integer :: iFlagSfr,iFlagFeedback,iFlagCool,nfiles,istart
  integer :: nextracols,nstarcols,i1,i2,i3,i4,lenblock,idumpformat
  integer, parameter :: iunit = 11, iunitd = 102, iunith = 103
  logical :: iexist,reallocate,mixedtypes,checkids
  real(doub_prec) :: timetemp,ztemp
  real(doub_prec), dimension(6) :: massoftypei
  real, dimension(:), allocatable :: dattemp1
  real :: hfact,dmdensi

  nstepsread = 0
  if (maxparttypes.lt.6) then
     print*,' *** ERROR: not enough particle types for GADGET data read ***'
     print*,' *** you need to edit splash parameters and recompile ***'
     stop
  endif
  
  if (len_trim(rootname).gt.0) then
     datfile = trim(rootname)
  else
     print*,' **** no data read **** ' 
     return
  endif
!
!--check if first data file exists
!
  inquire(file=datfile,exist=iexist)
  if (.not.iexist) then
     print "(a)",' *** error: '//trim(datfile)//': file not found ***'    
     return
  endif
!
!--set parameters which do not vary between timesteps
!
  ndim = 3
  ndimV = 3
  idumpformat = 0
  idumpformat = ienvironment('GSPLASH_FORMAT')
  checkids = lenvironment('GSPLASH_CHECKIDS')
!
!--read data from snapshots
!  
  i = istepstart

  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  !
  !--open data file and read data
  !
  open(iunit,iostat=ierr,file=datfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print "(a)", '*** ERROR OPENING FILE ***'
     return
  endif
  !
  !--read header for this timestep
  !
  if (idumpformat.eq.2) then
     print "(a)",' >> reading block labelled Gadget format <<'
     read(iunit,iostat=ierr) blocklabel,lenblock
     !print*,ierr,blocklabel,lenblock
     if (ierr /= 0 .or. lenblock.ne.264) then
        print "(a)",'*** ERROR READING HEADER: wrong endian? ***'
        close(iunit)
        return
     endif
  else
     print "(a)",' >> reading default Gadget format <<'
  endif

  read(iunit,iostat=ierr) npartoftypei(1:6),massoftypei,timetemp,ztemp, &
      iFlagSfr,iFlagFeedback,Nall(1:6),iFlagCool,nfiles

  ntoti = int(sum(npartoftypei(1:6)))
  if (ierr /= 0 .or. ntoti.le.0 .or. any(npartoftypei.lt.0)) then
     print "(/,a)", '*** ERROR READING TIMESTEP HEADER: wrong endian? ***'
     print "(/,a)", '   (see splash userguide for compiler-dependent'
     print "(a)", '    ways to change endianness on the command line)'
     print "(/,a)", '   (set environment variable GSPLASH_FORMAT to 2 '
     print "(a,/)", '    if you are using the block-labelled Gadget format)'
     close(iunit)
     return
  endif

  if (idumpformat.eq.2) then
     ncolstep = 1
     do while (ierr.eq.0)
        call read_blockheader(idumpformat,iunit,0,index2,blocklabelgas(ncolstep),lenblock,nvec)
        read(iunit,iostat=ierr)
        if ((ierr.eq.0 .and. index2.gt.0) .and. (index2.eq.ntoti &
            .or. index2.eq.npartoftypei(1) &
            .or. index2.eq.npartoftypei(2) &
            .or. index2.eq.npartoftypei(5) &
            .or. index2.eq.(npartoftypei(1)+npartoftypei(5)) & 
            .or. index2.eq.(npartoftypei(1)+npartoftypei(2)))) then
           select case(blocklabelgas(ncolstep))
           case('ID  ')
              ! not a column
           case default
              ncolstep = ncolstep + nvec
           end select
        endif
     enddo
     ncolstep = ncolstep - 1
     rewind(iunit)
     read(iunit,iostat=ierr)
     read(iunit,iostat=ierr)
     iformat = 2
     nextracols = 0
     nstarcols = 0

  else

     iformat = 0
     if (iFlagCool.gt.0 .and. .not.lenvironment('GSPLASH_IGNORE_IFLAGCOOL')) then
        iformat = 1
        ncolstep = 12 ! 3 x pos, 3 x vel, pmass, utherm, rho, Ne, Nh, h
        print "(a)",' cooling flag on  : assuming Ne, Nh dumped before h'
     else
        iformat = 0
        ncolstep = 10 ! 3 x pos, 3 x vel, pmass, utherm, rho, h
     endif
     if (iFlagSfr.gt.0) then
        print "(a)",' star formation flag on: assuming star formation rate dumped '
        ncolstep = ncolstep + 1
        iformat = iformat + 10
     endif

     call envlist('GSPLASH_EXTRACOLS',nextracols)
     if (nextracols.gt.0) then
        print "(a,i2,a)",' READING ',nextracols,' EXTRA COLUMNS '
        ncolstep = ncolstep + nextracols
     endif
     call envlist('GSPLASH_STARPARTCOLS',nstarcols)
     if (nstarcols.gt.0) then
        print "(a,i2,a)",' READING ',nstarcols,' STAR PARTICLE COLUMN(S) '
        ncolstep = ncolstep + nstarcols
     endif  
     !call envlist('GSPLASH_EXTRAVECCOLS',nextraveccols)
     !if (nextraveccols.gt.0) then
     !   print "(a,i2,a)",' READING ',nextraveccols,' EXTRA COLUMNS '
     !   ncolstep = ncolstep + nextraveccols
     !endif
  endif
  
  ncolumns = ncolstep
  !
  !--call set labels to get ih, ipmass, irho for use in the read routine
  !
  call set_labels
  
  print*,'time             : ',timetemp
  print*,'Npart (by type)  : ',npartoftypei
  print*,'Mass  (by type)  : ',massoftypei
  print*,'N_gas            : ',npartoftypei(1)
  print*,'N_total          : ',ntoti
  print*,'N data columns   : ',ncolstep

  if (nfiles.gt.1) then
     print*,' nfiles = ',nfiles
     print*,'*** ERROR: read from > 1 files not implemented'
     return
  endif
  !
  !--if successfully read header, increment the nstepsread counter
  !
  nstepsread = nstepsread + 1
  !
  !--now read data
  !
  reallocate = .false.
  npart_max = maxpart
  nstep_max = max(maxstep,1)

  if (ntoti.gt.maxpart) then
     reallocate = .true.
     if (maxpart.gt.0) then
        ! if we are reallocating, try not to do it again
        npart_max = int(1.1*ntoti)
     else
        ! if first time, save on memory
        npart_max = int(ntoti)
     endif
  endif
  if (i.ge.maxstep .and. i.ne.1) then
     nstep_max = i + max(10,INT(0.1*nstep_max))
     reallocate = .true.
  endif
  !
  !--reallocate memory for main data array
  !
  if (reallocate .or. .not.(allocated(dat))) then
     call alloc(npart_max,nstep_max,max(ncolstep+ncalc,maxcol))
  endif
  !
  !--copy header into header arrays
  !
  npartoftype(:,i) = npartoftypei
  !
  !--set time to be used in the legend
  !
  if (lenvironment('GSPLASH_USE_Z')) then
     !--use this line for redshift
     legendtext = 'z='
     time(i) = real(ztemp)
  else
     !--use this line for code time
     time(i) = real(timetemp) 
  endif  
  !
  !--read particle data
  !
  if (ntoti.gt.0) then
     !
     !--read positions of all particles
     !
     call read_blockheader(idumpformat,iunit,ntoti,index2,blocklabel,lenblock,nvec)
     if (iformat.eq.2 .and. blocklabel.ne.'POS ')  then
        print "(a)",' WARNING: expecting positions, got '//blocklabel//' in data read'
     endif
     if (any(required(1:3))) then
        print*,'positions ',ntoti
        read (iunit, iostat=ierr) (dat(j,1:3,i),j=1,index2)
        if (ierr /= 0) then
           print "(a)",'error encountered whilst reading positions '
           return
        endif
     else
        read(iunit, iostat=ierr)
        if (ierr /= 0) then
           print "(a)",'error skipping positions '
           return
        endif
     endif
     !
     !--same for velocities
     !
     call read_blockheader(idumpformat,iunit,ntoti,index2,blocklabel,lenblock,nvec)
     if (iformat.eq.2 .and. blocklabel.ne.'VEL ')  then
        print "(a)",' WARNING: expecting velocity, got '//blocklabel//' in data read'
     endif
     if (any(required(4:6))) then
        print*,'velocities ',index2
        read (iunit, iostat=ierr) (dat(j,4:6,i),j=1,index2)
        if (ierr /= 0) then
           print "(a)",'error encountered whilst reading velocities'
        endif
     else
        read(iunit, iostat=ierr)
        if (ierr /= 0) then
           print "(a)",'error skipping velocities '
           return
        endif
     endif
     !
     !--skip read of particle ID (only required if we sort the particles
     !  back into their correct order, which is not implemented at present)
     !  OR if using particle ID to flag dead particles
     !
     if (checkids) then
        print*,'particle ID ',ntoti
        if (allocated(iamtemp)) deallocate(iamtemp)
        allocate(iamtemp(npart_max))
        call read_blockheader(idumpformat,iunit,ntoti,index2,blocklabel,lenblock,nvec)
        if (iformat.eq.2 .and. blocklabel.ne.'ID  ') then
           print "(a)",' WARNING: expecting particle ID, got '//blocklabel//' in data read'
        endif
        if (index2.gt.0) read (iunit,iostat=ierr) iamtemp(1:index2)
     else
        call read_blockheader(idumpformat,iunit,ntoti,index2,blocklabel,lenblock,nvec)
        if (iformat.eq.2 .and. blocklabel.ne.'ID  ') then
           print "(a)",' WARNING: expecting particle ID, got '//blocklabel//' in data read'
        endif
        if (index2.gt.0) read (iunit,iostat=ierr) ! skip this line
     endif
     if (ierr /= 0) then
        print "(a)",'error encountered whilst reading particle ID'
     endif
     !
     !--read particle masses
     !
     !--work out total number of masses dumped 
     Nmassesdumped = 0
     do itype = 1,6
        if (abs(massoftypei(itype)).lt.tiny(massoftypei)) then
           Nmassesdumped = Nmassesdumped + Npartoftype(itype,i)
        endif
     enddo
     
     if (ipmass.eq.0) then
        masstype(1:6,i) = real(massoftypei(1:6))
     else
        if (required(ipmass)) then
           print*,'particle masses ',Nmassesdumped
           !--read this number of entries
           if (Nmassesdumped.gt.0) then
              if (allocated(dattemp1)) deallocate(dattemp1)
              allocate(dattemp1(Nmassesdumped))
              call read_blockheader(idumpformat,iunit,Nmassesdumped,index2,blocklabel,lenblock,nvec)
              if (iformat.eq.2 .and. blocklabel.ne.'MASS')  then
                 print "(a)",' WARNING: expecting particle masses, got '//blocklabel//' in data read'
              endif

           else
              index2 = 0
           endif

           if (index2.gt.0) then
              read(iunit,iostat=ierr) dattemp1(1:index2)
           endif
           if (ierr /= 0) then
              print "(a)",'error reading particle masses'
           endif
           !--now copy to the appropriate sections of the .dat array
           indexstart = 1
           index1 = 1

           do itype=1,6
              if (Npartoftype(itype,i).ne.0) then
                 index2 = index1 + Npartoftype(itype,i) -1
                 if (abs(massoftypei(itype)).lt.tiny(massoftypei)) then ! masses dumped
                    indexend = indexstart + Npartoftype(itype,i) - 1
                    print*,'read ',Npartoftype(itype,i),' masses for type ', &
                           itype,index1,'->',index2,indexstart,'->',indexend
                    dat(index1:index2,ipmass,i) = dattemp1(indexstart:indexend)
                    indexstart = indexend + 1
                 else  ! masses not dumped
                    print*,'setting masses for type ',itype,' = ', &
                           real(massoftypei(itype)),index1,'->',index2
                    dat(index1:index2,ipmass,i) = real(massoftypei(itype))
                 endif
                 index1 = index2 + 1
              endif
           enddo
           if (allocated(dattemp1)) deallocate(dattemp1)
        elseif (Nmassesdumped.gt.0) then
           read(iunit,iostat=ierr)
           if (ierr /= 0) then
              print "(a)",'error reading particle masses'
           endif
        endif
     endif
     !
     !--read other quantities for rest of particles
     !
     print*,'gas properties ',npartoftypei(1)
     if (ipmass.eq.0) then
        istart = 7
     else
        istart = 8
     endif
     icol = istart-1
     do while (icol.lt.ncolstep) !icol=istart,ncolstep !-nextraveccols
        !!print*,icol
        i3 = 0
        i4 = 0
        if (idumpformat.eq.2) then
           if (icol+1.le.ih) then
              call read_blockheader(idumpformat,iunit,npartoftypei(1),index2,blocklabel,lenblock,nvec)           
           else
              call read_blockheader(idumpformat,iunit,0,index2,blocklabel,lenblock,nvec)
           endif
           icol = icol + nvec
           
           if (index2.eq.ntoti) then
              i1 = 1
              i2 = i1 + index2 - 1
              print*,blocklabel//' (',index2,': all particles)'
           elseif (index2.eq.npartoftypei(1)) then
              i1 = 1
              i2 = i1 + index2 - 1
              print*,blocklabel//' (',index2,': gas particles only)'
           elseif (index2.eq.npartoftypei(2)) then
              i1 = npartoftypei(1) + 1
              i2 = i1 + index2 - 1
              print*,blocklabel//' (',index2,': dark matter particles only)'
           elseif (index2.eq.npartoftypei(1)+npartoftypei(2)) then
              i1 = 1
              i2 = i1 + index2 - 1
              print*,blocklabel//' (',index2,': gas+dark matter particles only)'
           elseif (index2.eq.npartoftypei(5)) then
              i1 = sum(npartoftypei(1:4)) + 1
              i2 = i1 + index2 - 1
              print*,blocklabel//' (',index2,': star particles only)'
           elseif (index2.eq.npartoftypei(1)+npartoftypei(5)) then
              i1 = 1
              i2 = npartoftypei(1)
              i3 = sum(npartoftypei(1:4)) + 1
              i4 = i3 + npartoftypei(5) - 1
              print*,blocklabel//' (',index2,': gas+star particles only)'
           else
              print*,blocklabel//': ERROR in block length/quantity defined on unknown mix of types n = (',index2,')'
              i1 = 1
              i2 = index2
           endif
        else
           nvec = 1
           icol = icol + nvec
           if (icol.gt.ncolstep-nstarcols) then
              i1 = sum(npartoftypei(1:4)) + 1
              i2 = i1 + npartoftypei(5) - 1
              print*,'star particle properties ',icol,i1,i2
           else
              i1 = 1
              i2 = npartoftypei(1)
           endif
        endif
        

        if (npartoftype(1,i).gt.0) then
           if (required(icol)) then
              if (i3.gt.0) then
                 read (iunit,iostat=ierr) dat(i1:i2,icol,i),dat(i3:i4,icol,i)
              else
                 if (nvec.gt.1) then
                    read (iunit, iostat=ierr) ((dat(k,j,i),j=icol-nvec+1,icol),k=i1,i2)
                 else
                    read (iunit,iostat=ierr) dat(i1:i2,icol,i)
                 endif
              endif
           else
              read (iunit,iostat=ierr)
           endif
           if (ierr /= 0) then
              print "(1x,a,i3)",'ERROR READING PARTICLE DATA from column ',icol
           endif
        !
        !--for some reason the smoothing length output by GADGET is
        !  twice the usual SPH smoothing length
        !
           if (icol.eq.ih .and. required(icol)) then
              dat(1:npartoftype(1,i),icol,i) = 0.5*dat(1:npartoftype(1,i),icol,i)
           endif
        endif
     enddo
     
     !if (nextraveccols.gt.0) then
     !   print*,'chemical species ',index2
     !   read (iunit, iostat=ierr) (dat(j,4:6,i),j=1,index2)
     !   if (ierr /= 0) then
     !      print "(a)",'error encountered whilst reading velocities'
     !   endif
     !endif
     
     !
     !--try to read dark matter and star particle smoothing lengths and/or density from a separate
     !  one column ascii file. If only density, use this to compute smoothing lengths.
     !
     densfile = trim(rootname)//'.dens'
     hfile = trim(rootname)//'.hsml'
     hfact = 1.2 ! related to the analytic neighbour number (hfact=1.2 gives 58 neighbours in 3D)
     open(unit=iunitd,file=densfile,iostat=ierrrho,status='old',form='formatted')
     open(unit=iunith,file=hfile,iostat=ierrh,status='old',form='formatted')
     
     if (ierrh.eq.0 .or. ierrrho.eq.0) then
        if (ierrh.eq.0) then
           print "(a)",' READING DARK MATTER SMOOTHING LENGTHS from '//trim(densfile)
           j = 1
           ierr = 0
           do while (ierr.eq.0 .and. j.le.sum(npartoftype(2:,i)))
              read(iunith,*,iostat=ierr) dat(npartoftype(1,i)+j,ih,i)
              j = j + 1
           enddo
           print "(a,i10,a)",' SMOOTHING LENGTHS READ for ',j-1,' dark matter / star particles '
           hsoft = 1.0 ! just so dark matter rendering is allowed in set_labels routine
           close(unit=iunith)
        endif
        
        if (ierrrho.eq.0) then
           print "(a)",' READING DARK MATTER DENSITIES FROM '//trim(densfile)
           j = 1
           ierr = 0
           nhset = 0
           do while (ierr.eq.0 .and. j.le.sum(npartoftype(2:,i)))
              read(iunitd,*,iostat=ierr) dmdensi
              if (ierr.eq.0) then
                 index1 = npartoftype(1,i) + j
                 dat(index1,irho,i) = dmdensi
                 if (ierrh.ne.0 .and. dmdensi.gt.tiny(dmdensi)) then
                    dat(index1,ih,i) = hfact*(dat(index1,ipmass,i)/dmdensi)**(1./3.)
                    nhset = nhset + 1
                 endif
              endif
              j = j + 1
           enddo
           print "(a,i10,a)",' DENSITIES READ FOR ',j-1,' dark matter / star particles '
           if (ierrh.ne.0) print "(a,i10,a,f5.2,a)", &
              ' SMOOTHING LENGTHS SET for ',nhset,' DM/star particles using h = ',hfact,'*(m/rho)**(1/3)'
           hsoft = 1.0 ! just so dark matter rendering is allowed in set_labels routine
           close(iunitd)
        endif
     else
     !
     !--if a value for the dark matter smoothing length is set
     !  via the environment variable GSPLASH_DARKMATTER_HSOFT,
     !  give dark matter particles this smoothing length
     !  and a density of 1 (so column density plots work)
     !
        hsoft = renvironment('GSPLASH_DARKMATTER_HSOFT')
        if (hsoft.gt.tiny(hsoft)) then
           if (required(ih)) then
              print "(a,1pe10.3,a)",' ASSIGNING SMOOTHING LENGTH of h = ',hsoft, &
                                    ' to dark matter particles'
              dat(npartoftype(1,i)+1:npartoftype(1,i)+npartoftype(2,i),ih,i) = hsoft
           endif
           if (required(irho)) then
              dat(npartoftype(1,i)+1:npartoftype(1,i)+npartoftype(2,i),irho,i) = 1.0
           endif
        else
           if (npartoftype(1,i).le.0 .and. sum(npartoftype(:,i)).gt.0) then
              print "(6(/,a),/)",' NOTE!! For gadget data using dark matter only, column density ',&
                                 ' plots can be produced by setting the GSPLASH_DARKMATTER_HSOFT ',&
                                 ' environment variable to give the dark matter smoothing length', &
                                 ' (for a fixed smoothing length)', &
                                 ' or by creating a one-column ascii file called '//trim(hfile), &
                                 ' containing the smoothing lengths for the dark matter particles'
           endif
        endif
     endif
     
     !
     !--DEAL WITH ACCRETED PARTICLES
     !  if particle ID is less than zero, treat this as an accreted particle
     !  (give it a negative smoothing length)
     !
     if (checkids) then
        nacc = 0
        do j=1,index2
           if (iamtemp(j) < 0) then
              if (required(ih)) then
                 dat(j,ih,i) = -abs(dat(j,ih,i))
              endif
              nacc = nacc + 1
           endif
        enddo
        if (nacc.gt.0) then
           print "(a,i10,a,/,a)",' marking ',nacc,' particles with negative ID as accreted/dead', &
             ' (giving them a negative smoothing length so they will be ignored in renderings)'        
        else
           print "(a)",' no particles with negative ID (ie. accreted particles) found'
        endif
        deallocate(iamtemp)
     endif
  else
     ntoti = 1
     npartoftype(1,i) = 1
     dat(:,:,i) = 0.
  endif
!
!--now memory has been allocated, set arrays which are constant for all time
!
  gamma = 5./3.
!
!--set flag to indicate that only part of this file has been read 
!
  if (.not.all(required(1:ncolstep))) ipartialread = .true.
!
!--close data file and return
!                    
  close(unit=iunit)

  if (nstepsread.gt.0) then
     print*,'>> last step ntot =',sum(npartoftype(:,istepstart+nstepsread-1))
  endif
  return
  
contains

!!-----------------------------------------------------------------
!! small utility to transparently handle block labelled data read
!!-----------------------------------------------------------------
 subroutine read_blockheader(idumpfmt,lun,nexpected,ndumped,blklabel,lenblk,nvec)
  implicit none
  integer, intent(in) :: idumpfmt,lun,nexpected
  integer, intent(out) :: ndumped
  character(len=4), intent(out) :: blklabel
  integer, intent(out) :: lenblk
  integer, intent(out) :: nvec
  
  blklabel = '    '
  if (idumpformat.eq.2) then
     read(lun, iostat=ierr) blklabel,lenblk
     if (ierr /= 0) then
        ndumped = 0
        return
     endif
     if (blklabel.eq.'POS ' .OR. blklabel.eq.'VEL ' .OR. blklabel.eq.'BFLD') then
        ndumped = (lenblk-8)/12
        nvec = 3
     else
        ndumped = (lenblk-8)/4
        nvec = 1
     endif
     !print*,blklabel,lenblk,ndumped
     !if (nexpected.gt.0) then
     !   if (ndumped.ne.nexpected) then
     !      !print*,'warning: number of '//blklabel//' dumped (',ndumped,') /= expected (',nexpected,')'
     !   endif
     !endif
  else
     ndumped = nexpected
  endif
  
  return
 end subroutine read_blockheader

end subroutine read_data

!!------------------------------------------------------------
!! set labels for each column of data
!!------------------------------------------------------------

subroutine set_labels
  use labels, only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass, &
              ih,irho,ipr,iutherm,iBfirst,idivB
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ntypes,UseTypeInRenderings,iformat
  use geometry, only:labelcoord
  use system_utils, only:envlist,ienvironment
  use gadgetread, only:hsoft,blocklabelgas
  use prompting, only:lcase
  implicit none
  integer :: i,nextracols,nstarcols,icol,ihset
  character(len=30), dimension(10) :: labelextra

  if (ndim.le.0 .or. ndim.gt.3) then
     print*,'*** ERROR: ndim = ',ndim,' in set_labels ***'
     return
  endif
  if (ndimV.le.0 .or. ndimV.gt.3) then
     print*,'*** ERROR: ndimV = ',ndimV,' in set_labels ***'
     return
  endif

  if (iformat.eq.2) then
     icol = 0
     do i=1,size(blocklabelgas)
        icol = icol + 1
        select case(blocklabelgas(i))
        case('POS ')
           ix(1) = icol
           ix(2) = icol+1
           ix(3) = icol+2
        case('VEL ')
           ivx = icol
        case('BFLD')
           iBfirst = icol
        case('MASS')
           ipmass = icol
        case('U   ')
           iutherm = icol
        case('RHO ')
           irho = icol
        case('NE  ')
           label(icol) = 'N\de\u'
        case('NH  ')
           label(icol) = 'N\dH\u'
        case('HSML')
           ih = icol
        case('NHP ')
           label(icol) = 'N\dH+\u'
        case('NHE ')
           label(icol) = 'N\dHe\u'
        case('NHEP')
           label(icol) = 'N\dHe+\u'
        case('elec')
           label(icol) = 'N\de\u'
        case('HI  ')
           label(icol) = 'HI'
        case('HII ')
           label(icol) = 'HII'
        case('HeI ')
           label(icol) = 'HeI'
        case('HeII')
           label(icol) = 'HeII'
        case('H2I ')
           label(icol) = 'H\d2\uI'
        case('H2II')
           label(icol) = 'H\d2\uII'
        case('HM  ')
           label(icol) = 'HM'
        case('SFR ')
           label(icol) = 'Star formation rate'
        case('TEMP')
           label(icol) = 'temperature'
        case('POT ')
           label(icol) = 'potential'
        case('AGE ')
           label(icol) = 'Stellar formation time'
        case('Z   ')
           label(icol) = 'Metallicity'
        case('ACCE')
           label(icol) = 'Acceleration'
        case('ENDT')
           label(icol) = 'd(Entropy)/dt'
        case('STRD')
           label(icol) = 'Stress (diagonal)'
        case('STRO')
           label(icol) = 'Stress (off-diagonal)'
        case('STRB')
           label(icol) = 'Stress (bulk)'
        case('SHCO')
           label(icol) = 'Shear coefficient'
        case('TSTP')
           label(icol) = 'Time step'
        case('DBDT')
           label(icol) = 'dB/dt'
        case('DIVB')
           label(icol) = 'div B'
           idivB = icol
        case('ABVC')
           label(icol) = 'alpha\dvisc\u'
        case('AMDC')
           label(icol) = 'alpha\dresist\u'
        case('PHI ')
           label(icol) = 'div B cleaning function'
        case('COOR')
           label(icol) = 'Cooling Rate'
        case('CONR')
           label(icol) = 'Conduction Rate'
        case('BFSM')
           label(icol) = 'B\dsmooth\u'
        case('DENN')
           label(icol) = 'Denn'
        case('CRC0')
           label(icol) = 'Cosmic Ray C0'
        case('CRP0')
           label(icol) = 'Cosmic Ray P0'
        case('CRE0')
           label(icol) = 'Cosmic Ray E0'
        case('CRn0')
           label(icol) = 'Cosmic Ray n0'
        case('CRco')
           label(icol) = 'Cosmic Ray Thermalization Time'
        case('CRdi')
           label(icol) = 'Cosmic Ray Dissipation Time'
        case('BHMA')
           label(icol) = 'Black hole mass'
        case('BHMD')
           label(icol) = 'black hole mass accretion rate'
        case('MACH')
           label(icol) = 'Mach number'
        case('DTEG')
           label(icol) = 'dt (energy)'
        case('PSDE')
           label(icol) = 'Pre-shock density'
        case('PSEN')
           label(icol) = 'Pre-shock energy'
        case('PSXC')
           label(icol) = 'Pre-shock X\d\u'
        case('DJMP')
           label(icol) = 'Density jump'
        case('EJMP')
           label(icol) = 'Energy jump'
        case('CRDE')
           label(icol) = 'Cosmic Ray injection'        
        case('ID  ')
           icol = icol - 1
        case default
           label(icol) = trim(lcase(blocklabelgas(i)))
        end select
     enddo
  else
     do i=1,ndim
        ix(i) = i
     enddo
     ivx = 4
     ipmass = 7
     irho = 9        ! location of rho in data array
     ipr = 0
     iutherm = 8     !  thermal energy
     if (iformat.eq.1 .or. iformat.eq.11 .and. ncolumns.gt.10) then
        label(10) = 'Ne'
        label(11) = 'Nh'
        ih = 12        !  smoothing length
        if (iformat.eq.11) label(13) = 'Star formation rate'
     else
        ih = 10
        if (iformat.eq.1) label(11) = 'Star formation rate'
     endif
     ihset = ienvironment('GSPLASH_HSML_COLUMN',errval=-1)
     if (ihset.gt.0) ih = ihset
     !
     !--deal with extra columns
     !
     if (ncolumns.gt.ih) then
        call envlist('GSPLASH_EXTRACOLS',nextracols,labelextra)
        do i=ih+1,ih+nextracols
           label(i) = trim(labelextra(i-ih))
        enddo
        call envlist('GSPLASH_STARPARTCOLS',nstarcols,labelextra)
        do i=ih+nextracols+1,ih+nextracols+nstarcols
           label(i) = trim(labelextra(i-ih-nextracols))
        enddo
     endif
  endif
  !
  !--set labels of the quantities read in
  !
  if (ix(1).gt.0) label(ix(1:ndim)) = labelcoord(1:ndim,1)
  if (irho.gt.0) label(irho) = 'density'
  if (iutherm.gt.0) label(iutherm) = 'u'
  if (ipmass.gt.0) label(ipmass) = 'particle mass'
  if (ih.gt.0) label(ih) = 'h'
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
     do i=1,ndimV
        label(iBfirst+i-1) = trim(labelvec(iBfirst))//'\d'//labelcoord(i,1)
     enddo
  endif
  
  !--set labels for each particle type
  !
  ntypes = 5
  labeltype(1) = 'gas'
  labeltype(2) = 'dark matter'
  labeltype(5) = 'star'
  UseTypeInRenderings(1) = .true.
  !
  !--dark matter particles are of non-SPH type (ie. cannot be used in renderings)
  !  unless they have had a smoothing length defined
  !
  if (hsoft.gt.tiny(hsoft)) then
     UseTypeInRenderings(2) = .true.
  else
     UseTypeInRenderings(2) = .false.  
  endif
  UseTypeInRenderings(3:5) = .false.

!-----------------------------------------------------------
  return
end subroutine set_labels
