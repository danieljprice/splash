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

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR OUTPUT FROM THE GADGET CODE
! (works with GADGET v1.0, v2.0 and v3.0)
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
! GSPLASH_CHECKIDS if 'YES','yes','TRUE' or 'true' then reads and checks
!  particle IDs for negative values and flags these as accreted particles
! GSPLASH_HSML_COLUMN if set to a positive integer, specifies the location
!  of the smoothing length in the columns, overriding any default settings.
! GSPLASH_IGNORE_IFLAGCOOL if set to 'YES' or `TRUE', does not assume that
!  extra columns are present even if the cooling flag is set in the header.
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
 logical :: havewarned = .false.

end module gadgetread

subroutine read_data(rootname,istepstart,ipos,nstepsread)
  use particle_data,  only:dat,npartoftype,masstype,time,gamma,maxpart,maxcol,maxstep
  use params,         only:doub_prec,sing_prec,maxparttypes
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc,iformat,required,ipartialread, &
                           ntypes,debugmode,iverbose
  use settings_page,  only:legendtext
  use mem_allocation, only:alloc
  use labels,         only:ih,irho,ipmass,labeltype
  use system_utils,   only:renvironment,lenvironment,ienvironment,envlist
  use gadgetread,     only:hsoft,blocklabelgas,havewarned
  implicit none
  integer, intent(in)                :: istepstart,ipos
  integer, intent(out)               :: nstepsread
  character(len=*), intent(in)       :: rootname
  character(len=len(rootname)+10)    :: datfile,densfile,hfile
  character(len=4)                   :: blocklabel
  character(len=20)                  :: string
  integer, dimension(maxparttypes)   :: npartoftypei,Nall
  integer, dimension(:), allocatable :: iamtemp
  integer               :: i,j,k,n,itype,icol,ierr,ierrh,ierrrho,nhset,nvec,ifile
  integer               :: index1,index2,indexstart,indexend,nmassesdumped,ntypesused
  integer               :: ncolstep,npart_max,nstep_max,ntoti,nacc,ntotall,idot
  integer               :: iFlagSfr,iFlagFeedback,iFlagCool,nfiles,istart,nhfac
  integer               :: nextracols,nstarcols,i1,i2,i3,i4,lenblock,idumpformat
  integer, dimension(6) :: i0,i1all,i2all
  integer, parameter    :: iunit = 11, iunitd = 102, iunith = 103
  logical               :: iexist,reallocate,checkids,usez,goterrors
  logical, dimension(6) :: ireadtype
  real(doub_prec)                    :: timetemp,ztemp
  real(doub_prec), dimension(6)      :: massoftypei
  real(sing_prec), dimension(:),   allocatable :: dattemp1
  real(sing_prec), dimension(:,:), allocatable :: dattemp
  real :: hfact,hfactmean
  real, parameter :: pi = 3.1415926536

  nstepsread = 0
  goterrors  = .false.
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
     !
     !--look for a file with .0 on the end for multiple-file reads
     !
     datfile=trim(rootname)//'.0'
     inquire(file=datfile,exist=iexist)
     if (.not.iexist) then
        print "(a)",' *** error: '//trim(rootname)//': file not found ***'
        return
     endif
  endif
!
!--set parameters which do not vary between timesteps
!
  ndim  = 3
  ndimV = 3
  idumpformat = 0
  idumpformat = ienvironment('GSPLASH_FORMAT')
  checkids    = lenvironment('GSPLASH_CHECKIDS')
  usez        = lenvironment('GSPLASH_USE_Z')
!
!--read data from snapshots
!
  i = istepstart
!
!--i0 is the offset used to read the data into the arrays
!  (non-zero for read from multiple files)
!  The offset is different for each particle type, somewhat
!  complicating the data read -- we shuffle the particles from
!  multiple files so that they are in type order.
!
  i0(:) = 0
!
!--loop over the number of files
!
  ifile = 0
  ntotall = 0
  over_files: do while(iexist)

  write(*,"(23('-'),1x,a,1x,23('-'))") trim(datfile)
  ifile = ifile + 1

  !
  !--open data file and read data
  !
  open(iunit,iostat=ierr,file=datfile,status='old',form='unformatted')
  if (ierr /= 0) then
     print "(a)", '*** ERROR OPENING FILE ***'
     return
  endif

  !if (any(i0.gt.0)) print*,'starting position for each type in data array: ',i0(:)
  !
  !--read header for this timestep
  !
  if (idumpformat.eq.2) then
     print "(a)",' >> reading block labelled Gadget format <<'
     read(iunit,iostat=ierr) blocklabel,lenblock
     !print*,ierr,blocklabel,lenblock
     if (ierr /= 0 .or. lenblock.ne.264) then
        print "(/,a,/)",'*** ERROR READING HEADER: wrong endian? or wrong format? ***'
        close(iunit)
        if (ifile.eq.1) then
           return
        else
           exit over_files
        endif
     endif
  else
     if (ifile.eq.1) print "(a)",' >> reading default Gadget format <<'
  endif

  npartoftypei(:) = 0
  Nall(:) = 0
  massoftypei(:) = 0.
  iFlagCool = 0
  nfiles = 0
  read(iunit,iostat=ierr) npartoftypei(1:6),massoftypei(1:6),timetemp,ztemp, &
      iFlagSfr,iFlagFeedback,Nall(1:6),iFlagCool,nfiles

  ntoti = int(sum(npartoftypei(1:6)))  ! int here is unnecessary, but avoids compiler warnings

  if (nfiles.gt.1) then
     ntotall = int(sum(Nall(1:6)))
  else
     ntotall = ntoti
  endif
  if (debugmode) then
     print*,'DEBUG: ierr = ',ierr
     print*,'DEBUG: ntoti = ',ntoti,' ntotall = ',ntotall,' nfiles = ',nfiles
     print*,'DEBUG: npartoftype = ',npartoftypei(1:6),' Nall = ',Nall(1:6)
     print*,'DEBUG: iFlagSfr = ',iFlagSfr,' iFlagFeedback = ',iFlagFeedback,' iFlagCool = ',iFlagCool
     print*,'DEBUG: time = ',timetemp,' z = ',ztemp
  endif

  if (ierr /= 0 .or. ntoti.le.0 .or. ntotall.le.0 .or. any(npartoftypei.lt.0) .or. nfiles.lt.0 &
      .or. nfiles.gt.1e6) then
     print "(/,a)", '*** ERROR READING TIMESTEP HEADER: wrong endian? ***'
     print "(/,a)", '   (see splash userguide for compiler-dependent'
     print "(a)",   '    ways to change endianness on the command line)'
     print "(/,a)", '   (set environment variable GSPLASH_FORMAT to 2 '
     print "(a,/)", '    if you are using the block-labelled Gadget format)'
     close(iunit)
     if (ifile.eq.1) then
        return
     else
        exit over_files
     endif
  endif

  !
  !--if we are reading from multiple files,
  !  check that the sequence starts from the correct file
  !
  if (nfiles.gt.1) then
     idot = len_trim(datfile)-1
     if (ifile.eq.1 .and. datfile(idot:idot+1).ne.'.0') then
        if (nfiles.lt.100) then
           string = "(/,a,i2,a,/,a,/)"
        else
           string = "(/,a,i7,a,/,a,/)"
        endif
        if (nfiles.gt.10000) then
           !--this is the most likely scenario here
           print "(a)",'*** ERROR reading timestep header: wrong endian? ***'
        else
           print string,' ERROR: read is from multiple files (nfiles = ',nfiles,')',&
                     '        but this is not the first file (does not end in .0): skipping...'
        endif
        close(iunit)
        return
     endif
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
     if (iFlagCool.eq.1 .and. .not.lenvironment('GSPLASH_IGNORE_IFLAGCOOL')) then
        iformat = 1
        ncolstep = 12 ! 3 x pos, 3 x vel, pmass, utherm, rho, Ne, Nh, h
        if (ifile.eq.1) print "(a)",' cooling flag on : assuming Ne, Nh dumped before h'
     else
        iformat = 0
        ncolstep = 10 ! 3 x pos, 3 x vel, pmass, utherm, rho, h
     endif
     if (iFlagSfr.eq.1) then
        if (ifile.eq.1) print "(a)",' star formation flag on: assuming star formation rate dumped '
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

  if (ifile.eq.1) then
     ncolumns = ncolstep
  !
  !--call set labels to get ih, ipmass, irho for use in the read routine
  !
     hsoft = 0. ! to avoid unset variable
     call set_labels
  endif

  if (ifile.eq.1) then
     print*,'time            : ',timetemp
     if (usez) then
        print "(1x,a,f8.2,a)",'z (redshift)    : ',ztemp,' (using in legend from GSPLASH_USE_Z setting)'
     else
        print "(1x,a,f8.2,a)",'z (redshift)    : ',ztemp,' (set GSPLASH_USE_Z=yes to use in legend)'
     endif
  endif
  print*,'Npart (by type) : ',npartoftypei(1:6)
  if (ifile.eq.1) print*,'Mass  (by type) : ',massoftypei(1:6)
!  print "(10x,'|',6(1x,a12,'|'))",   (labeltype(itype),itype=1,ntypes)
!  print "(a10,'|',6(i11,2x,'|'))",   'Npart  : ',npartoftypei
!  print "(a10,'|',6(es11.3,2x,'|'))",'Mass   : ',massoftypei
  print*,'N_gas           : ',npartoftypei(1)
  print*,'N_total         : ',ntoti
  if (ifile.eq.1) print*,'N data columns  : ',ncolstep
  if (nfiles.gt.1 .and. ifile.eq.1) then
     print*,'Nall            : ',Nall(1:6)
  endif

  if (nfiles.gt.1) then
     if (ifile.eq.1) print "(a,i4,a)",' reading from ',nfiles,' files'
  elseif (nfiles.lt.0) then
     print*,'*** ERROR: nfiles = ',nfiles,' in file header: aborting'
     return
  endif

  if (ifile.eq.1) then
     !--Softening lengths for Dark Matter Particles...
     hsoft = renvironment('GSPLASH_DARKMATTER_HSOFT')
     !
     !--try to read dark matter and star particle smoothing lengths and/or density from a separate
     !  one column ascii file. If only density, use this to compute smoothing lengths.
     !
     densfile = trim(rootname)//'.dens'
     hfile = trim(rootname)//'.hsml'
     hfact = 1.2 ! related to the analytic neighbour number (hfact=1.2 gives 58 neighbours in 3D)
     open(unit=iunitd,file=densfile,iostat=ierrrho,status='old',form='formatted')
     open(unit=iunith,file=hfile,iostat=ierrh,status='old',form='formatted')

     if (idumpformat.eq.2) then
        if (ih.eq.0 .and. (hsoft.gt.tiny(hsoft) .or. ierrrho.eq.0 .or. ierrh.eq.0)) then
           ncolumns = ncolumns + 1
           blocklabelgas(ncolumns) = 'HSML'
           ih = ncolumns
           call set_labels
        endif
        if (irho.eq.0 .and. (hsoft.gt.tiny(hsoft) .or. ierrrho.eq.0 .or. ierrh.eq.0)) then
           ncolumns = ncolumns + 1
           blocklabelgas(ncolumns) = 'RHO '
           irho = ncolumns
           call set_labels
        endif
     endif

  !
  !--if successfully read header, increment the nstepsread counter
  !
     nstepsread = nstepsread + 1
  endif
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
        npart_max = int(1.1*ntotall)
     else
        ! if first time, save on memory
        npart_max = int(ntotall)
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
     call alloc(npart_max,nstep_max,max(ncolumns+ncalc,maxcol))
  endif

  !
  !--copy npartoftypei into allocated header arrays
  !  and set the offset position of particle types in the main data arrays
  !
  if (nfiles.eq.1 .or. ifile.eq.1) then
     i0(1) = 0
     do itype=2,ntypes
        if (nfiles.eq.1) then
           i0(itype) = sum(npartoftypei(1:itype-1)) ! this is avoid depending on Nall at all for single file read
        else
           i0(itype) = sum(Nall(1:itype-1))
        endif
     enddo
     npartoftype(:,i) = npartoftypei
  else
     i0(1) = npartoftype(1,i)
     do itype=2,ntypes
        i0(itype) = sum(Nall(1:itype-1)) + npartoftype(itype,i)
     enddo
     npartoftype(:,i) = npartoftype(:,i) + npartoftypei
  endif
  if (debugmode) print*,'DEBUG: starting position for each type in data array: ',i0(:)
  !
  !--set time to be used in the legend
  !
  if (ifile.eq.1) then
     if (usez) then
        !--use this line for redshift
        legendtext = 'z='
        time(i) = real(ztemp)
     else
        !--use this line for code time
        time(i) = real(timetemp)
     endif
  else
     if (usez) then
        if (abs(real(ztemp)-time(i)).gt.tiny(0.)) print*,'ERROR: redshift different between files in multiple-file read'
     else
        if (abs(real(timetemp)-time(i)).gt.tiny(0.)) print*,'ERROR: time different between files in multiple-file read'
     endif
     if (sum(Nall).ne.ntotall) then
        print*,' ERROR: Nall differs between files'
        goterrors = .true.
     endif
  endif
  !
  !--read particle data
  !
  got_particles: if (ntoti.gt.0) then
     !
     !--read positions of all particles
     !  (note that errors on position read are fatal)
     !
     call read_blockheader(idumpformat,iunit,ntoti,index2,blocklabel,lenblock,nvec)
     if (iformat.eq.2 .and. blocklabel.ne.'POS ')  then
        print "(a)",' WARNING: expecting positions, got '//blocklabel//' in data read'
     endif
     if (any(required(1:3))) then
        print*,'positions ',index2
        if (allocated(dattemp)) deallocate(dattemp)
        allocate(dattemp(3,ntoti))
        read(iunit,iostat=ierr) (dattemp(:,j),j=1,index2)
        if (nfiles.gt.1) then
           !
           !--read data into type order if multiple files are present:
           !  this means the offset position is different for each type
           !
           if (sum(npartoftypei).ne.index2) print*,' ERROR: number of positions .ne. sum of types'
           n = 0
           do itype=1,ntypes
              do j=i0(itype)+1,i0(itype)+npartoftypei(itype)
                 n = n + 1
                 dat(j,1:3,i) = dattemp(1:3,n)
              enddo
           enddo
           !read (iunit, iostat=ierr) ((dat(j,1:3,i),j=i0(itype)+1,i0(itype)+npartoftypei(itype)),itype=1,ntypes)
        else
           do j=1,index2
              dat(j,1:3,i) = dattemp(1:3,j)
           enddo
!           read (iunit, iostat=ierr) (dat(j,1:3,i),j=1,index2)
        endif
        if (ierr /= 0) then
           print "(a)",'error encountered whilst reading positions '
           deallocate(dattemp)
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
        if (.not.allocated(dattemp)) allocate(dattemp(3,ntoti))

        read (iunit, iostat=ierr) (dattemp(:,j),j=1,index2)
        if (nfiles.gt.1) then
           !--see above re: type order
           if (sum(npartoftypei).ne.index2) print*,' ERROR: number of velocities .ne. sum of types'
           n = 0
           do itype=1,ntypes
              do j=i0(itype)+1,i0(itype)+npartoftypei(itype)
                 n = n + 1
                 dat(j,4:6,i) = dattemp(1:3,n)
              enddo
           enddo
           !read (iunit, iostat=ierr) ((dat(j,4:6,i),j=i0(itype)+1,i0(itype)+npartoftypei(itype)),itype=1,ntypes)
        else
           do j=1,index2
              dat(j,4:6,i) = dattemp(1:3,j)
           enddo
           !read (iunit, iostat=ierr) (dat(j,4:6,i),j=1,index2)
        endif
        if (ierr /= 0) then
           print "(a)",'error encountered whilst reading velocities'
           goterrors = .true.
        endif
     else
        read(iunit, iostat=ierr)
        if (ierr /= 0) then
           print "(a)",'error skipping velocities '
           if (allocated(dattemp)) deallocate(dattemp)
           return
        endif
     endif
     if (allocated(dattemp)) deallocate(dattemp)
     !
     !--skip read of particle ID (only required if we sort the particles
     !  back into their correct order, which is not implemented at present)
     !  OR if using particle ID to flag dead particles
     !
     !  For multiple files we only allocate and read the IDs for one file
     !
     if (checkids) then
        print*,'particle ID ',ntoti
        if (allocated(iamtemp)) deallocate(iamtemp)
        allocate(iamtemp(ntoti))
     endif

     call read_blockheader(idumpformat,iunit,ntoti,index2,blocklabel,lenblock,nvec)
     if (iformat.eq.2 .and. blocklabel.ne.'ID  ') then
        print "(a)",' WARNING: expecting particle ID, got '//blocklabel//' in data read'
     endif

     if (index2.gt.0) then
        if (checkids .and. required(ih)) then
           !--particle IDs are currently only used to set h -ve for accreted particles
           !  so do not read if h not required
           read (iunit,iostat=ierr) iamtemp(1:index2)
        else
           read (iunit,iostat=ierr) ! skip this line
        endif
        if (ierr /= 0) then
           print "(a)",'error encountered whilst reading particle ID'
           goterrors = .true.
        endif
     endif
     !
     !--read particle masses
     !
     !--work out total number of masses dumped
     nmassesdumped = 0
     do itype = 1,6
        if (abs(massoftypei(itype)).lt.tiny(massoftypei)) then
           nmassesdumped = nmassesdumped + npartoftypei(itype)
        endif
     enddo

     if (ipmass.eq.0) then
        masstype(1:6,i) = real(massoftypei(1:6))
     else
        if (required(ipmass)) then
           print*,'particle masses ',nmassesdumped
           !--read this number of entries
           if (nmassesdumped.gt.0) then
              if (allocated(dattemp1)) deallocate(dattemp1)
              allocate(dattemp1(nmassesdumped))
              call read_blockheader(idumpformat,iunit,nmassesdumped,index2,blocklabel,lenblock,nvec)
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
              goterrors = .true.
           endif
           !--now copy to the appropriate sections of the dat array
           indexstart = 1
           !index1 = 1

           do itype=1,6
              if (npartoftypei(itype).ne.0) then
                 !--work out the appropriate section of the dat array for this particle type
                 index1 = i0(itype) + 1
                 index2 = i0(itype) + npartoftypei(itype)

                 if (abs(massoftypei(itype)).lt.tiny(massoftypei)) then ! masses dumped
                    indexend = indexstart + npartoftypei(itype) - 1
                    if (debugmode) &
                       print*,' read ',npartoftypei(itype),' masses for '//trim(labeltype(itype))// &
                              ' particles',index1,'->',index2,indexstart,'->',indexend

                    dat(index1:index2,ipmass,i) = dattemp1(indexstart:indexend)
                    indexstart = indexend + 1
                 else  ! masses not dumped
                    if (debugmode) print "(a,es10.3,i10,a,i10)",&
                      ' setting masses for '//trim(labeltype(itype))//' particles = ', &
                       real(massoftypei(itype)),index1,'->',index2

                    dat(index1:index2,ipmass,i) = real(massoftypei(itype))
                 endif
                 !index1 = index2 + 1
              endif
           enddo
           if (allocated(dattemp1)) deallocate(dattemp1)
        elseif (nmassesdumped.gt.0) then
           read(iunit,iostat=ierr)
           if (ierr /= 0) then
              print "(a)",'error reading particle masses'
              goterrors = .true.
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
     gas_properties: do while (icol.lt.ncolstep) !icol=istart,ncolstep !-nextraveccols
        !!print*,icol
        i3 = 0
        i4 = 0
        ireadtype(:) = .false.

        if (idumpformat.eq.2) then
           if (icol+1.le.ih) then
              call read_blockheader(idumpformat,iunit,npartoftypei(1),index2,blocklabel,lenblock,nvec)
           else
              call read_blockheader(idumpformat,iunit,0,index2,blocklabel,lenblock,nvec)
           endif
           icol = icol + nvec
           !
           !--work out from the number of entries what mix of particle types
           !  the quantity is defined on
           !
           if (index2.eq.ntoti) then
              i1 = i0(1) + 1
              i2 = i1 + ntoti - 1
              print*,blocklabel//' (',index2,': all particles)'
              ireadtype(:) = .true.
           elseif (index2.eq.npartoftypei(1)) then
              i1 = i0(1) + 1
              i2 = i1 + index2 - 1
              print*,blocklabel//' (',index2,': gas particles only)'
              ireadtype(1) = .true.
           elseif (index2.eq.npartoftypei(2)) then
              i1 = i0(2) + 1
              i2 = i1 + index2 - 1
              print*,blocklabel//' (',index2,': dark matter particles only)'
              ireadtype(2) = .true.
           elseif (index2.eq.npartoftypei(1)+npartoftypei(2)) then
              i1 = i0(1) + 1
              i2 = i1 + index2 - 1
              print*,blocklabel//' (',index2,': gas+dark matter particles only)'
              ireadtype(1:2) = .true.
           elseif (index2.eq.npartoftypei(5)) then
              i1 = i0(5) + 1
              i2 = i1 + index2 - 1
              print*,blocklabel//' (',index2,': star particles only)'
              ireadtype(5) = .true.
           elseif (index2.eq.npartoftypei(1)+npartoftypei(5)) then
              i1 = i0(1) + 1
              i2 = i1 + npartoftypei(1) - 1
              i3 = i0(5) + 1
              i4 = i3 + npartoftypei(5) - 1
              print*,blocklabel//' (',index2,': gas+star particles only)'
              ireadtype(1) = .true.
              ireadtype(5) = .true.
           else
              print*,blocklabel//': ERROR in block length/quantity defined on unknown mix of types n = (',index2,')'
              i1 = i0(1)+1
              i2 = i0(1)+index2
           endif
        else
           nvec = 1
           icol = icol + nvec
           if (icol.gt.ncolstep-nstarcols) then
              i1 = i0(5) + 1
              i2 = i1 + npartoftypei(5) - 1
              print*,'star particle properties ',icol,i1,i2
              ireadtype(5) = .true.
           else
              !--default is a quantity defined only on gas particles
              i1 = i0(1) + 1
              i2 = i1 + npartoftypei(1) - 1
              ireadtype(1) = .true.
           endif
        endif

        !
        !--construct the array offsets required when reading from multiple files
        !
        ntypesused = 0
        do itype=1,6
           if (ireadtype(itype) .and. npartoftypei(itype).gt.0) then
              ntypesused = ntypesused + 1
              i1all(ntypesused) = i0(itype) + 1
              i2all(ntypesused) = i0(itype) + npartoftypei(itype)
           endif
        enddo

        if (npartoftypei(1).gt.0) then
           if (required(icol)) then
              if (i3.gt.0) then
                 if (nfiles.gt.1) then
                    read (iunit,iostat=ierr) (dat(i1all(itype):i2all(itype),icol,i),itype=1,ntypesused)
                 else
                    read (iunit,iostat=ierr) dat(i1:i2,icol,i),dat(i3:i4,icol,i)
                 endif
              else
                 if (nvec.gt.1) then
                    if (nfiles.gt.1) then
                       read (iunit,iostat=ierr) &
                        (((dat(k,j,i),j=icol-nvec+1,icol),k=i1all(itype),i2all(itype)),itype=1,ntypesused)
                    else
                       read (iunit,iostat=ierr) ((dat(k,j,i),j=icol-nvec+1,icol),k=i1,i2)
                    endif
                 else
                    if (nfiles.gt.1) then
                       read (iunit,iostat=ierr) (dat(i1all(itype):i2all(itype),icol,i),itype=1,ntypesused)
                    else
                       read (iunit,iostat=ierr) dat(i1:i2,icol,i)
                    endif
                 endif
              endif
           else
              read (iunit,iostat=ierr)
           endif
           if (ierr /= 0) then
              print "(1x,a,i3)",'ERROR READING PARTICLE DATA from column ',icol
              goterrors = .true.
           endif
        endif
     enddo gas_properties

     !if (nextraveccols.gt.0) then
     !   print*,'chemical species ',index2
     !   read (iunit, iostat=ierr) (dat(j,4:6,i),j=1,index2)
     !   if (ierr /= 0) then
     !      print "(a)",'error encountered whilst reading velocities'
     !   endif
     !endif
!
!--close data file now that we have finished reading data
!
     close(unit=iunit)

     !
     !--DEAL WITH ACCRETED PARTICLES (in this file only)
     !  if particle ID is less than zero, treat this as an accreted particle
     !  (give it a negative smoothing length)
     !
     if (checkids) then
        nacc = 0
        !--only do this if the smoothing length is required in the data read
        if (required(ih)) then
           n = 0
           !do itype=1,ntypes
           itype = 1
           do j=1,npartoftypei(itype)
              n = n + 1
              if (iamtemp(n) < 0) then
                 !if (itype.gt.1) print*,' id -ve on non-gas particle ',itype,j
                 dat(i0(itype)+j,ih,i) = -abs(dat(i0(itype)+j,ih,i))
                 nacc = nacc + 1
              endif
           enddo
           !enddo
           if (nacc.gt.0) then
              print "(a,i10,a,/,a)",' marking ',nacc,' '//trim(labeltype(1))// &
                ' particles with negative ID as accreted/dead', &
                ' (giving them a negative smoothing length so they will be ignored in renderings)'
           else
              print "(a)",' no particles with negative ID (i.e. accreted particles) found'
           endif
        endif
        if (allocated(iamtemp)) deallocate(iamtemp)
     endif

  endif got_particles
!
!--now memory has been allocated, set arrays which are constant for all time
!
  gamma = 5./3.
!
!--set flag to indicate that only part of this file has been read
!
  if (.not.all(required(1:ncolstep))) ipartialread = .true.
!
!--for read from multiple files, work out the next file in the sequence
!
  iexist = .false.
  if (nfiles.gt.1 .and. ifile.lt.nfiles) then
     !--see if the next file exists
     idot = index(datfile,'.',back=.true.)
     if (idot.le.0) then
        print "(a)",' ERROR: read from multiple files but could not determine next file in sequence'
        goterrors = .true.
     else
        write(string,*) ifile
        write(datfile,"(a,i1)") trim(datfile(1:idot))//trim(adjustl(string))
        iexist = .false.
        inquire(file=datfile,exist=iexist)
        if (.not.iexist) then
           print "(a)",' ERROR: read from multiple files '// &
           'but could not find '//trim(datfile)//': next in sequence'
           goterrors = .true.
        endif
     endif
  endif

  enddo over_files
  !
  !--for some reason the smoothing length output by GADGET is
  !  twice the usual SPH smoothing length
  !  (do this after we have read data from all of the files)
  !
  if (required(ih) .and. ih.gt.0 .and. size(dat(1,:,:)).ge.ih .and. npartoftype(1,i).gt.0) then
     print "(a)",' converting GADGET smoothing length on gas particles to usual SPH definition (x 0.5)'
     dat(1:npartoftype(1,i),ih,i) = 0.5*dat(1:npartoftype(1,i),ih,i)
  endif

  if (nfiles.gt.1. .and. any(npartoftype(:,i).ne.Nall(:))) then
     print*,'ERROR: sum of Npart across multiple files .ne. Nall in data read '
     print*,'Npart = ',npartoftype(:,i)
     print*,'Nall  = ',Nall(:)
     goterrors = .true.
  endif
  !
  !--look for dark matter smoothing length/density files
  !
  if (ierrh.eq.0 .or. ierrrho.eq.0) then
     if (ierrh.eq.0) then
        print "(a)",' READING DARK MATTER SMOOTHING LENGTHS from '//trim(hfile)
        ierr = 0
        index1 = npartoftype(1,i)+1
        index2 = npartoftype(1,i)+sum(npartoftype(2:,i))
        read(iunith,*,iostat=ierr) (dat(j,ih,i),j=index1,index2)
        close(unit=iunith)
        if (ierr.lt.0) then
           nhset = 0
           do j=index1,index2
              if (dat(j,ih,i).gt.0.) nhset = nhset + 1
           enddo
           print "(a,i10,a,/)",' *** END-OF-FILE: GOT ',nhset,' SMOOTHING LENGTHS ***'
        elseif (ierr.gt.0) then
           print "(a)", ' *** ERROR reading smoothing lengths from file'
           goterrors = .true.
        else
           print "(a,i10,a)",' SMOOTHING LENGTHS READ OK for ',index2-index1+1,' dark matter / star particles '
        endif
        hsoft = 1.0 ! just so dark matter rendering is allowed in set_labels routine
     endif

     if (ierrrho.eq.0) then
        print "(a)",' READING DARK MATTER DENSITIES FROM '//trim(densfile)
        ierr = 0
        index1 = npartoftype(1,i)+1
        index2 = npartoftype(1,i)+sum(npartoftype(2:,i))
        read(iunitd,*,iostat=ierr) (dat(j,irho,i),j=index1,index2)
        close(iunitd)
        if (ierr.lt.0) then
           nhset = 0
           do j=index1,index2
              if (dat(j,irho,i).gt.0.) nhset = nhset + 1
           enddo
           print "(a,i10,a,/)",' *** END-OF-FILE: GOT ',nhset,' DENSITIES ***'
        elseif (ierr.gt.0) then
           print "(a)", ' *** ERROR reading dark matter densities from file'
           goterrors = .true.
        else
           print "(a,i10,a)",' DENSITY READ OK for ',index2-index1+1,' dark matter / star particles '
        endif
        if (ierrh.ne.0 .and. ipmass.gt.0) then
           where(dat(:,irho,i) > tiny(dat))
              dat(:,ih,i) = hfact*(dat(:,ipmass,i)/dat(:,irho,i))**(1./3.)
           elsewhere
              dat(:,ih,i) = 0.
           end where
           print "(a,i10,a,f5.2,a)", &
            ' SMOOTHING LENGTHS SET for ',j-1-index1,' DM/star particles using h = ',hfact,'*(m/rho)**(1/3)'
        endif
        hsoft = 1.0 ! just so dark matter rendering is allowed in set_labels routine
     endif
  else
  !
  !--if a value for the dark matter smoothing length is set
  !  via the environment variable GSPLASH_DARKMATTER_HSOFT,
  !  give dark matter particles this smoothing length
  !  and a density of 1 (so column density plots work)
  !
     if (hsoft.gt.tiny(hsoft)) then
        if (required(ih)) then
           print "(a,1pe10.3,a)",' ASSIGNING SMOOTHING LENGTH of h = ',hsoft, &
                                 ' to dark matter particles'
           !print*,'ih = ',ih,' npartoftype = ',npartoftype(1:2,i), shape(dat)
           if (ih.gt.0) then
              dat(npartoftype(1,i)+1:npartoftype(1,i)+npartoftype(2,i),ih,i) = hsoft
           else
              print*,' ERROR: smoothing length not found in data arrays'
              goterrors = .true.
           endif
        endif
        if (required(irho)) then
           if (irho.gt.0) then
              dat(npartoftype(1,i)+1:npartoftype(1,i)+npartoftype(2,i),irho,i) = 1.0
           else
              print*,' ERROR: place for density not found in data arrays'
              goterrors = .true.
           endif
        endif
     else
        if (npartoftype(1,i).le.0 .and. sum(npartoftype(:,i)).gt.0) then
           print "(66('*'),4(/,a),/)",'* NOTE!! For GADGET data using dark matter only, column density ',&
                              '* plots can be produced by setting the GSPLASH_DARKMATTER_HSOFT ',&
                              '* environment variable to give the dark matter smoothing length', &
                              '* (for a fixed smoothing length)'
           hsoft = (maxval(dat(:,1,i)) - minval(dat(:,1,i)))/sum(npartoftype(2:,i))**(1./3.)
           print*,' suggested value for GSPLASH_DARKMATTER_HSOFT = ',hsoft
           hsoft = 0.

           print "(7(/,a),/)",'* Alternatively, and for best results, calculate a number density', &
                                        '* on dark matter particles, set individual smoothing lengths from', &
                                        '* this using h = hfact*(n)**(-1/3), with hfact=1.2 and either ', &
                                        '* dump the results back into the HSML array in the original dump ', &
                                        '* file (if using the block-labelled format), or create an ascii ',&
                                        '* file called '//trim(hfile)//' containing the smoothing length ',&
                                        '* values for the dark matter particles.'
           print "(2(/,a),/,66('*'),/)",  '* Also make sure normalised interpolations are OFF when plotting ',&
                                        '* dark matter density '
        endif
     endif
  endif
!
!--pause with fatal errors
!
  if (goterrors .and. .not.lenvironment('GSPLASH_IGNORE_ERRORS')) then
     print "(/,a)",'*** ERRORS detected during data read: data will be corrupted'
     print "(a,/)",'    Please REPORT this and/or fix your file ***'
     print "(a)",'     (set GSPLASH_IGNORE_ERRORS=yes to skip this message)'
     if (iverbose.ge.1) then
        print "(a)",'    > Press any key to bravely proceed anyway  <'
        read*
     endif
  endif
!
!--give a friendly warning about using too few or too many neighbours
!  (only works with equal mass particles because otherwise we need the number density estimate)
!
  if (ih.gt.0 .and. required(ih) .and. ipmass.gt.0 .and. required(ipmass) &
      .and. abs(massoftypei(1)).lt.tiny(0.) .and. ndim.eq.3 .and. .not.havewarned) then
     nhfac = 100
     if (npartoftype(1,i).gt.nhfac) then
        hfactmean = 0.
        do j=1,nhfac
           hfact = dat(j,ih,i)*(dat(j,irho,i)/(dat(j,ipmass,i)))**(1./ndim)
           hfactmean = hfactmean + hfact
        enddo
        hfact = hfactmean/real(nhfac)
        havewarned = .true.
        if (hfact.lt.1.125 .or. hfact.gt.1.45) then
           print "(/,a)",'** FRIENDLY NEIGHBOUR WARNING! **'
           print "(3x,a,f5.1,a,/,3x,a,f5.2,a,i1,a)", &
                 'It looks like you are using around ',4./3.*pi*(2.*hfact)**3,' neighbours,', &
                 'corresponding to h = ',hfact,'*(m/rho)^(1/',ndim,') in 3D:'

           if (hfact.lt.1.15) then
              print "(4(/,3x,a))",'This is a quite a low number of neighbours for the cubic spline and ', &
                                  'may result in increased noise and inaccurate wave propagation speeds', &
                                  '(a cubic lattice is also an unstable initial configuration for the ',&
                                  ' particles in this regime -- see Morris 1996, Borve et al. 2004).'
           elseif (hfact.gt.1.45) then
              print "(4(/,3x,a))",'Using h >~ 1.5*(m/rho)^(1/3) with the cubic spline results in the', &
                                'particle pairing instability due to the first neighbour being placed under', &
                                'the hump in the kernel gradient. Whilst not fatal, it results in a', &
                                'loss of resolution so is a bit of a waste of cpu time.'
              print "(4(/,3x,a))",'If you are attempting to perform a "resolution study" by increasing the', &
                                'neighbour number, this is a *bad idea*, as you are also increasing h.',      &
                                '(a better way is to increase the smoothness of the integrals without changing h',   &
                                ' by adopting a smoother kernel such as the M5 Quintic that goes to 3h).'
           endif
           print "(/,3x,a,/,3x,a,/)", &
              'A good default range is h = 1.2-1.3 (m/rho)^1/ndim ', &
              'corresponding to around 58-75 neighbours in 3D.'
        else
           print "(/,1x,a,f5.1,a,/,1x,a,f5.2,a,i1,a,/)", &
                'Simulations employ ',4./3.*pi*(2.*hfact)**3,' neighbours,', &
                'corresponding to h = ',hfact,'*(m/rho)^(1/',ndim,') in 3D'
        endif
     endif
  else
     !print*,'not true'
  endif
!
!--cover the special case where no particles have been read
!
  if (ntotall.le.0) then
     npartoftype(1,i) = 1
     dat(:,:,i) = 0.
  endif

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
  if (idumpfmt.eq.2) then
     read(lun, iostat=ierr) blklabel,lenblk
     if (ierr /= 0) then
        ndumped = 0
        return
     endif
     if (blklabel.eq.'POS ' .OR. blklabel.eq.'VEL ' .OR. blklabel.eq.'ACCE' .OR. blklabel.eq.'BFLD' .OR. &
         blklabel.eq.'BPOL' .OR. blklabel.eq.'BTOR') then
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
  use labels,        only:label,iamvec,labelvec,labeltype,ix,ivx,ipmass, &
                          ih,irho,ipr,iutherm,iBfirst,iBpol,iBtor,idivB,iax
  use params
  use settings_data, only:ndim,ndimV,ncolumns,ntypes,UseTypeInRenderings,iformat
  use geometry,      only:labelcoord
  use system_utils,  only:envlist,ienvironment
  use gadgetread,    only:hsoft,blocklabelgas
  use asciiutils,    only:lcase
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
        case('ACCE')
           iax = icol
        case('BFLD')
           iBfirst = icol
        case('BPOL')
           iBpol = icol
        case('BTOR')
           iBtor = icol
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
        case('PRES')
           label(icol) = 'pressure'
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
  if (ix(1).gt.0)   label(ix(1:ndim)) = labelcoord(1:ndim,1)
  if (irho.gt.0)    label(irho)       = 'density'
  if (iutherm.gt.0) label(iutherm)    = 'u'
  if (ipmass.gt.0)  label(ipmass)     = 'particle mass'
  if (ih.gt.0)      label(ih)         = 'h'
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

  if (iax.gt.0) then
     iamvec(iax:iax+ndimV-1) = iax
     labelvec(iax:iax+ndimV-1) = 'a'
     do i=1,ndimV
        label(iax+i-1) = trim(labelvec(iax))//'\d'//labelcoord(i,1)
     enddo
   endif

  if (iBfirst.gt.0) then
     iamvec(iBfirst:iBfirst+ndimV-1) = iBfirst
     labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
     do i=1,ndimV
        label(iBfirst+i-1) = trim(labelvec(iBfirst))//'\d'//labelcoord(i,1)
     enddo
  endif

  if (iBpol.gt.0) then
     iamvec(iBpol:iBpol+ndimV-1) = iBpol
     labelvec(iBpol:iBpol+ndimV-1) = 'B\dpol'
     do i=1,ndimV
        label(iBpol+i-1) = trim(labelvec(iBpol))//'\d'//labelcoord(i,1)
     enddo
   endif

  if (iBtor.gt.0) then
     iamvec(iBtor:iBtor+ndimV-1) = iBtor
     labelvec(iBtor:iBtor+ndimV-1) = 'B\dtor'
     do i=1,ndimV
        label(iBtor+i-1) = trim(labelvec(iBtor))//'\d'//labelcoord(i,1)
     enddo
   endif

  !--set labels for each particle type
  !
  ntypes = 6
  labeltype(1) = 'gas'
  labeltype(2) = 'dark matter'
  labeltype(3) = 'boundary 1'
  labeltype(4) = 'boundary 2'
  labeltype(5) = 'star'
  labeltype(6) = 'sink / black hole'
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
  UseTypeInRenderings(3:6) = .false.

!-----------------------------------------------------------
  return
end subroutine set_labels
