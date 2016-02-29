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
!  Copyright (C) 2005-2012 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-------------------------------------------------------------------------
! this subroutine reads from the data file(s)
! change this to change the format of data input
!
! THIS VERSION IS FOR DATA FORMATS WRITTEN WITH THE H5PART LIBRARY
!
! SOME CHOICES FOR THIS FORMAT CAN BE SET USING THE FOLLOWING
!  ENVIRONMENT VARIABLES:
!
!  H5SPLASH_NDIM=2     : number of spatial dimensions (overrides value inferred from data)
!  H5SPLASH_HFAC=1.2   : factor to use in h= hfac*(m/rho)**(1/ndim) if h not present in data
!  H5SPLASH_HSML=1.0   : value for global smoothing length if h not present in data
!  H5SPLASH_TYPEID='MatID'  : name of dataset containing the particle types
!
! the data is stored in the global array dat
!
! >> this subroutine must return values for the following: <<
!
! ncolumns    : number of data columns
! ndim, ndimV : number of spatial, velocity dimensions
! nstepsread  : number of steps read from this file
!
! dat(maxplot,maxpart,maxstep) : main data array
!
! npartoftype(1:6,maxstep) : number of particles of each type in each timestep
! ntot(maxstep)       : total number of particles in each timestep
!
! time(maxstep)       : time at each step
! gamma(maxstep)      : gamma at each step
!
! most of these values are stored in global arrays
! in the module 'particle_data'
!-------------------------------------------------------------------------

!--local module to store header information so we can later set the labels
module h5partdataread
 use params
 use labels, only:lenlabel
 implicit none
 character(len=lenlabel), dimension(maxplot) :: datasetnames
 logical :: warn_labels = .true.

end module h5partdataread

subroutine read_data(rootname,indexstart,ipos,nstepsread)
  use particle_data,   only:dat,iamtype,npartoftype,time,gamma,maxpart,maxcol,maxstep
  use params
  use settings_data,   only:ndim,ndimV,ncolumns,ncalc,debugmode,ntypes,iverbose
  use mem_allocation,  only:alloc
  use iso_c_binding,   only:c_double,c_int64_t
  use asciiutils,      only:lcase
  use system_utils,    only:renvironment
  use system_commands, only:get_environment
  use labels,          only:ih,ipmass,irho,ix
  use h5part
  use h5partattrib,    only:h5pt_readstepattrib,h5pt_getnstepattribs,h5pt_getstepattribinfo, &
                            h5pt_readstepattrib_r8
  use h5partdataread
  implicit none
  integer, intent(in)               :: indexstart,ipos
  integer, intent(out)              :: nstepsread
  character(len=*), intent(in)      :: rootname
  integer                           :: i,j,ncolstep,nsteps,ncolsfile,icol,itypeidcol
  integer(kind=c_int64_t)           :: nattrib,iattrib,nelem,idatasettype,k,icolsfile,ierr
  integer                           :: nprint,npart_max,nstep_max,maxcolsfile,itype
  integer, dimension(maxplot)       :: iorder
  integer, dimension(maxparttypes)  :: itypemap
  integer, dimension(:),allocatable :: itypefile
  logical                           :: iexist,typeiddefault
  integer(kind=c_int64_t)           :: ifile,istep
  character(len=len(rootname)+5)    :: dumpfile
  character(len=64)                 :: attribname
  character(len=lenlabel)           :: datasetname,type_datasetname
  real(kind=doub_prec),dimension(1) :: dtime
  real                              :: hsmooth,hfac,dndim

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart

  dumpfile = trim(rootname)

  if (iverbose.ge.1) print "(1x,a)",'reading h5part format'
  print "(26('>'),1x,a,1x,26('<'))",trim(dumpfile)
  !
  !--check if first data file exists
  !
  inquire(file=dumpfile,exist=iexist)
  if (.not.iexist) then
     print "(a)",' *** error: '//trim(dumpfile)//': file not found ***'
     return
  endif
  !
  !--fix number of spatial dimensions (0 means no particle coords)
  !
  ndim = 0
  ndimV = 0
  j = indexstart
  nstepsread = 0

  !
  !--open the file and read the number of particles
  !
  if (debugmode) print*,'DEBUG: opening '//trim(dumpfile)
  !ierr = h5pt_set_verbosity_level(6_8)

  ifile = h5pt_openr(trim(dumpfile))
  if (ifile.le.0) then
     print "(a)",'*** ERROR opening '//trim(dumpfile)//' ***'
     return
  endif
  if (debugmode) print*,'DEBUG: file opened ok'

  !
  !--get number of steps and particles in the file
  !
  nsteps = int(h5pt_getnsteps(ifile))
  if (debugmode) print*,'DEBUG: nsteps = ',nsteps

  !
  !--read environment variable giving the name of the dataset
  !  containing the particle type ID
  !  give default value if this is not set
  !
  call get_environment('H5SPLASH_TYPEID',type_datasetname)
  if (len_trim(type_datasetname).le.0) then
     typeiddefault = .true.
     type_datasetname = 'MatID'
  else
     typeiddefault = .false.
  endif
  !
  !--read "header" information from all steps in file:
  !  get maximum number of particles for all steps in file
  !  and maximum number of columns (datasets) in file
  !
  npart_max = 0
  ncolstep  = 0
  maxcolsfile = 0
  itypeidcol = 0
  do istep=1,nsteps
     if (debugmode) print "(a,i2)",'DEBUG: setting step ',istep
     ierr = h5pt_setstep(ifile,istep)
     if (ierr.eq.0) then
        npart_max = max(npart_max,int(h5pt_getnpoints(ifile)))
        ncolsfile = int(h5pt_getndatasets(ifile))
        icol = 0
        if ((istep.eq.nsteps .and. ncolsfile.gt.0)) then
           do icolsfile=0,ncolsfile-1
              !ierr = h5pt_getdatasetname(ifile,icol,datasetnames(icol+1))
              if (debugmode) print*,'DEBUG: getting datasetinfo'
              ierr = h5pt_getdatasetinfo(ifile,icolsfile,datasetname,idatasettype,nelem)
              if (ierr.ne.0) then
                 print "(a,i3)",' ERROR reading dataset name for column ',icolsfile+1
              else
                 !
                 !--read only columns which contain real or double precision data
                 !
                 select case(idatasettype)
                 case(H5PART_FLOAT32,H5PART_FLOAT64)
                    icol = icol + 1
                    datasetnames(icol) = trim(adjustl(datasetname))
                    !print*,' data set info = ',icol,trim(datasetnames(icol)),idatasettype,nelem
                 case(H5PART_INT32,H5PART_INT64)
                    !
                    !--try to recognise the dataset giving the particle types
                    !
                    if (trim(type_datasetname)==trim(datasetname) .or. trim(datasetname)=='Phase') then
                       type_datasetname = trim(datasetname)
                       if (itypeidcol.le.0) itypeidcol = int(icolsfile) + 1
                       if (iverbose.ge.1) print "(a)",' getting particle types from data set '//trim(type_datasetname)
                    else
                       if (iverbose.ge.1) print "(a)",' skipping data set '//trim(datasetname)// &
                                                      ' of type '//h5part_type(int(idatasettype))
                    endif
                 case default
                    if (iverbose.ge.1) print "(a)",' skipping data set '//trim(datasetname)//&
                                                   ' of type '//h5part_type(int(idatasettype))
                 end select
              endif
           enddo
        elseif (ncolsfile.le.0) then
           print*,'ERROR: number of datasets in step ',istep,' = ',ncolsfile
        endif
        ncolstep = max(ncolstep,icol)
        maxcolsfile = max(ncolsfile,maxcolsfile)
     else
        print "(a,i3)",' ERROR, could not choose step ',istep
        return
     endif
 enddo
  nprint = npart_max

  !
  !--warn if no particle type data has been read
  !
  if (itypeidcol.le.0) then
     if (typeiddefault) then
        if (iverbose.ge.1) print "(a)",' Particle type dataset not found in file: Use H5SPLASH_TYPEID to give dataset name'
     else
        print "(a)",' WARNING: Particle type dataset '//trim(type_datasetname)//' (from H5SPLASH_TYPEID) not found in file '
     endif
  endif

  !
  !--call the set_labels routine to get the initial location of coords, smoothing length etc. given dataset labels
  !
  warn_labels = .false.
  call set_labels()
  warn_labels = .true.

  !
  !--set default ordering of columns
  !
  do i=1,size(iorder)
     iorder(i) = i
  enddo
  !
  !--if coordinates are not in the first 3 columns, shift data so that they are
  !

  if (ndim.gt.0 .and. ix(1).ne.1) then
     do i=1,ndim
        iorder(ix(i)) = i
     enddo
     !--preserve the order of things after the coordinates
     icol = ndim
     do i=ix(1)+1,ncolstep
        if (.not.any(ix(1:ndim).eq.i)) then
           icol = icol + 1
           iorder(i) = icol
        endif
     enddo
     !--shuffle things before the coordinates to the end
     do i=1,ix(1)-1
        if (.not.any(ix(1:ndim).eq.i)) then
           icol = icol + 1
           iorder(i) = icol
        endif
     enddo
  endif
  if (debugmode) print*,'DEBUG: iorder = ',iorder

  !
  !--if smoothing length has not been set, look for an environment variable
  !  giving the smoothing length value
  !
  hsmooth = -1.
  if (ih.eq.0) then
     hsmooth = renvironment('H5SPLASH_HSML',errval=-1.)
     if (hsmooth.ge.0.) then
        ncolstep = ncolstep + 1
     elseif (ipmass.gt.0. .and. irho.gt.0 .and. ndim.gt.0) then
        hfac = renvironment('H5SPLASH_HFAC',errval=-1.)
        if (hfac.gt.0.) then
           if (iverbose.ge.1) print "(/,a,f6.2,a,/)",' Setting smoothing length using h = ',hfac,&
                                  '*(m/rho)**(1/ndim) (from H5SPLASH_HFAC setting)'
        else
           hfac = 1.2
           if (iverbose.ge.1) then
              print "(/,a)",' WARNING: Smoothing length not found in data: using h = hfac*(m/rho)**(1/ndim)'
              print "(a)",  '          (hfac = 1.2 by default, set H5SPLASH_HFAC to change this)'
              print "(a,/)",'          (set constant h with H5SPLASH_HSML to give a global value)'
           endif
        endif
        ncolstep = ncolstep + 1
     else
        if (iverbose.ge.1) print "(/,a,/)",' WARNING: Smoothing length not found in data: Set H5SPLASH_HSML to give a global value'
     endif
  endif
  ncolumns = ncolstep

  !
  !--allocate memory for all data in the file
  !
  nstep_max = max(nsteps,indexstart,1,maxstep)
  npart_max = max(maxpart,npart_max)
  if (.not.allocated(dat) .or. (nprint.gt.maxpart) .or. (ncolstep+ncalc).gt.maxcol) then
     if (itypeidcol.gt.0) then
        call alloc(npart_max,nstep_max,ncolstep+ncalc,mixedtypes=.true.)
     else
        call alloc(npart_max,nstep_max,ncolstep+ncalc)
     endif
  endif

!
!--now read the timestep data in the dumpfile (for all steps)
!
  istep = 0
  do j=indexstart,indexstart+nsteps-1
     istep = istep + 1
     print "(a,i4,a,i10)",' step ',istep,': ntotal = ',nprint
     ierr = h5pt_setstep(ifile,istep)
     nprint = int(h5pt_getnpoints(ifile))  ! use int() to avoid compiler warning about type conversion
!
!--get the time from the step attributes
!
     nattrib = h5pt_getnstepattribs(ifile)
     if (nattrib.gt.0) then
        do iattrib=0,nattrib-1   ! yes, it's written in C
           ierr = h5pt_getstepattribinfo(ifile,iattrib,attribname,nelem)
           !print*,' step attribute '//trim(attribname),' nelem = ',nelem
           if (ierr.eq.0) then
              !
              !--match anything that looks vaguely like the time
              !
              if (nelem.eq.1 .and. (index(lcase(attribname),'time').ne.0  &
                 .or. index(lcase(attribname),'t ').ne.0)) then
                 ierr = h5pt_readstepattrib_r8(ifile,attribname,dtime)
                 if (ierr.eq.0) then
                    time(j) = real(dtime(1))
                    print "(12x,a,es10.3,a)",'time   = ',time(j),' (from '//trim(attribname)//')'
                 else
                    print "(a,i2,a)",' ERROR could not read time from step ',istep,' (from '//trim(attribname)//')'
                 endif
              !
              !--match gamma if possible
              !
              elseif (nelem.eq.1 .and. (index(lcase(attribname),'gamma').ne.0  &
                 .or. index(lcase(attribname),'gam ').ne.0)) then
                 ierr = h5pt_readstepattrib_r8(ifile,attribname,dtime)

                 if (ierr.eq.0) then
                    gamma(j) = real(dtime(1))
                    print "(12x,a,es10.3,a)",'gamma  = ',gamma(j),' (from '//trim(attribname)//')'
                 else
                    print "(a,i2,a)",' ERROR could not read gamma from step ',istep,' (from '//trim(attribname)//')'
                 endif
              else
                 print "(a)",' unknown attribute '//trim(attribname)
              endif
           else
              print "(a,i3,a,i2)",' ERROR reading attribute info for step ',istep,', attribute #',iattrib
           endif
        enddo
     endif
!
!--now read the data for this step
!
    icol = 0
    do k=0,maxcolsfile-1
       ierr = h5pt_getdatasetinfo(ifile,k,datasetname,idatasettype,nelem)
       select case(idatasettype)
       case(H5PART_FLOAT32,H5PART_FLOAT64)
          icol = icol + 1
          datasetnames(iorder(icol)) = trim(datasetname)
          if (debugmode) print "(a,i3,a,i3)",'DEBUG: reading data set ',icol,&
                                             ': '//trim(datasetnames(iorder(icol)))//' into column ',iorder(icol)
          ierr = h5pt_readdata(ifile,datasetnames(iorder(icol)),dat(:,iorder(icol),j))
          if (ierr.ne.0) print "(a)",' ERROR reading dataset '//trim(datasetnames(iorder(icol)))
       case default
          if (debugmode) print "(a)",' skipping data set '//trim(datasetname)//' of type '//lcase(h5part_type(int(idatasettype)))
       end select
    enddo
!
!--read the particle types for this step from the typeid dataset (specified from the type_datasetname setting)
!
    npartoftype(:,j) = 0
    if (itypeidcol.gt.0 .and. size(iamtype(:,1)).gt.1) then
       if (debugmode) print "(a)",'DEBUG: reading particle types from '//trim(type_datasetname)
       !
       !--allocate temporary memory
       !
       if (allocated(itypefile)) deallocate(itypefile)
       allocate(itypefile(nprint),stat=ierr)
       if (ierr.ne.0) stop 'ERROR allocating temporary memory for particle types'
       !
       !--read type array from file
       !
       ierr = h5pt_readdata(ifile,trim(type_datasetname),itypefile(:))
       if (ierr.ne.0) then
          print "(a)",' ERROR reading dataset '//trim(type_datasetname)
       else
       !
       !--work out the number of unique particle types
       !  and map these into SPLASH particle types (1->maxtypes)
       !
          if (j.eq.1 .and. istep.eq.1) then
             ntypes      = 1
             itypemap(1) = minval(itypefile)
          endif
          do i=1,nprint
             !--increase the number of particle types if a particle of new type is found
             if (.not.any(itypemap(1:ntypes).eq.itypefile(i))) then
                ntypes = ntypes + 1
                if (ntypes.le.size(itypemap)) then
                   itypemap(ntypes) = itypefile(i)
                   npartoftype(ntypes,j) = npartoftype(ntypes,j) + 1
                   iamtype(i,j) = ntypes
                endif
             else
                do itype=1,ntypes
                   if (itypefile(i).eq.itypemap(itype)) then
                      npartoftype(itype,j) = npartoftype(itype,j) + 1
                      iamtype(i,j) = itype
                   endif
                enddo
             endif
          enddo
          if (nprint.lt.1e6) then
             print "(12x,a,10(i5,1x))",'npart (by type) = ',npartoftype(1:ntypes,j)
          else
             print "(12x,a,10(i10,1x))",'npart (by type) = ',npartoftype(1:ntypes,j)
          endif
          !
          !--warn if the number of types exceeds the current limit
          !
          if (ntypes.gt.maxparttypes) &
             print "(/,2(a,i2),a/)", &
              ' WARNING: too many particle types in dataset '//trim(type_datasetname)// &
              ' (got ',ntypes,': maximum is currently ',maxparttypes,')'
       endif
       !
       !--clean up
       !
       if (allocated(itypefile)) deallocate(itypefile)
    else
!--only one particle type
       ntypes = 1
       npartoftype(1,j) = nprint
    endif
!
!--reset the labels now that the columns have been read in the correct order
!
    if (j.eq.indexstart) then ! set labels based on the first step read from the file
       warn_labels = .false.
       call set_labels()
       warn_labels = .true.
    endif
!
!--if smoothing length set via environment variable, fill the extra column with the smoothing length value
!
    if (ih.eq.0) then
       if (hsmooth.ge.0.) then
          datasetnames(ncolstep) = 'h'
          ih = ncolstep
          dat(:,ih,j) = hsmooth
       elseif (ipmass.gt.0 .and. irho.gt.0 .and. ndim.gt.0) then
          ih = ncolstep
          datasetnames(ncolstep) = 'h'
          dndim = 1./ndim
          where (dat(:,irho,j).gt.tiny(0.))
            dat(:,ih,j) = hfac*(dat(:,ipmass,j)/dat(:,irho,j))**dndim
          elsewhere
            dat(:,ih,j) = 0.
          end where
       endif
    endif
!    read(iunit,*,iostat=ierr) (dat(i,icol,j),icol = 1,ncolstep)
    nstepsread = nstepsread + 1
  enddo

  ierr = h5pt_close(ifile)

return 
end subroutine read_data

!!-------------------------------------------------------------------
!! set labels for each column of data
!!
!! read these from a file called 'columns' in the current directory
!! then take sensible guesses as to which quantities are which
!! from the column labels
!!
!!-------------------------------------------------------------------

subroutine set_labels()
  use asciiutils,      only:lcase
  use labels,          only:label,ix,irho,ipmass,ih,iutherm, &
                            ipr,ivx,iBfirst,iamvec,labelvec,lenlabel !,labeltype
  !use params,          only:maxparttypes
  use settings_data,   only:ndim,ndimV,UseTypeInRenderings,iverbose
  use geometry,        only:labelcoord
  use system_utils,    only:ienvironment
  use h5partdataread
  implicit none
  integer                 :: i,ndimset,ndim_max
  character(len=lenlabel) :: labeli

  ndim = 0
  ndimV = 0
  ndimset = ienvironment('H5SPLASH_NDIM',errval=-1)
  ndim_max = 3
  if (ndimset.ge.0) ndim_max = ndimset
  irho = 0
  ih = 0
  ipmass = 0

  do i=1,size(datasetnames)
     if (len_trim(datasetnames(i)).gt.0) then
        label(i) = trim(datasetnames(i))
     else
        label(i) = ' '
     endif
!--now try to recognise the column based on the dataset name
!  compare all strings in lower case, trimmed and with no preceding spaces
!
     labeli = trim(adjustl(lcase(label(i))))
     if (index(labeli,'coords').ne.0 .and. index(labeli,'_').ne.0 .or. labeli(1:1).eq.'x') then
        if (ndim.lt.ndim_max) then
           ndim = ndim + 1
           ix(ndim) = i
           label(ix(ndim)) = labelcoord(ndim,1)
        endif
     elseif (index(labeli,'vel_').ne.0 .and. (ivx.eq.0 .or. i.le.ivx+ndim)) then
        if (ndimV.lt.3) ndimV = ndimV + 1
        if (index(labeli,'_0').ne.0) ivx = i
     elseif (index(labeli,'dens').ne.0 .and. irho.eq.0) then
        irho = i
     elseif (index(labeli,'mass').ne.0 .and. ipmass.eq.0) then
        ipmass = i
     elseif (ih.eq.0 .and. (index(labeli,'smoothing').ne.0 .or. labeli(1:1).eq.'h')) then
        ih = i
     elseif (labeli(1:1).eq.'u') then
        iutherm = i
     !--identify vector quantities based on _0, _1, _2 labelling
     elseif (index(labeli,'_0').ne.0) then
        !print*,'labelling ',labeli(1:index(labeli,'_0')-1),' as vector, column ',i
        iamvec(i) = i
        labelvec(i) = labeli(1:index(labeli,'_0')-1)
     elseif (index(labeli,'_1').ne.0 .and. i.gt.1 .and. ndim.ge.2) then
        if (iamvec(i-1).gt.0) then
           iamvec(i) = i-1
           labelvec(i) = labelvec(i-1)
        endif
     elseif (index(labeli,'_2').ne.0 .and. i.gt.2 .and. ndim.ge.3) then
        if (iamvec(i-2).gt.0) then
           iamvec(i) = i-2
           labelvec(i) = labelvec(i-2)
        endif
     endif
  enddo

  if (ndim.lt.1) ndimV = 0
  if (ndimV.gt.ndim) ndimV = ndim

  if (warn_labels .and. iverbose.ge.1) then
     if (ndimset.gt.0) then
        if (ndim.ne.ndimset) then
           print "(2(a,i1))",' WARNING: ndim = ',ndimset, &
                             ' from H5SPLASH_NDIM setting but coords not found in data: using ndim = ',ndim
        else
           print "(a,i1,a)",' Assuming number of dimensions = ',ndim,' from H5SPLASH_NDIM setting'
        endif
     else
        if (ndim.gt.0) print "(a,i1,a)",' Assuming number of dimensions = ',ndim,' (set H5SPLASH_NDIM to override)'
     endif

     if (ndimV.gt.0) print "(a,i1)",' Assuming vectors have dimension = ',ndimV
     if (irho.gt.0) print "(a,i2)",' Assuming density in column ',irho
     if (ipmass.gt.0) print "(a,i2)",' Assuming particle mass in column ',ipmass
     if (ih.gt.0) print "(a,i2)",' Assuming smoothing length in column ',ih
     if (iutherm.gt.0) print "(a,i2)",' Assuming thermal energy in column ',iutherm
     if (ipr.gt.0) print "(a,i2)",' Assuming pressure in column ',ipr
     if (ivx.gt.0) then
        if (ndimV.gt.1) then
           print "(a,i2,a,i2)",' Assuming velocity in columns ',ivx,' to ',ivx+ndimV-1
        else
           print "(a,i2)",' Assuming velocity in column ',ivx
        endif
     endif
     if (ndim.eq.0 .or. irho.eq.0 .or. ipmass.eq.0 .or. ih.eq.0) then
        print "(4(/,a))",' NOTE: Rendering capabilities cannot be enabled', &
                    '  until positions of density, smoothing length and particle', &
                    '  mass are known (for the h5part read this means labelling ', &
                    '  the dataset appropriately)'
     endif

     !
     !--assign vectors (don't do this on the first call otherwise it will remain assigned to the wrong columns)
     !
     if (ivx.gt.0) then
        iamvec(ivx:ivx+ndimV-1) = ivx
        labelvec(ivx:ivx+ndimV-1) = 'v'
     endif
     if (iBfirst.gt.0) then
        iamvec(iBfirst:iBfirst+ndimV-1) = ivx
        labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
     endif

     !
     !--set labels for vector quantities
     !
     do i=1,size(datasetnames)
        if (iamvec(i).ne.0) then
           label(i) = trim(labelvec(iamvec(i)))//'_'//trim(labelcoord(i-iamvec(i)+1,1))
        endif
     enddo
  endif
  !
  !--set labels for each particle type
  !  (for h5part this is done in the read_data routine)
  !

  !ntypes = 1 !!maxparttypes
!  labeltype(1) = 'gas'
!  labeltype(2) = 'gas'
!  labeltype(3) = 'gas'
!  labeltype(4) = 'gas'
  UseTypeInRenderings(:) = .true.


!-----------------------------------------------------------

  return
end subroutine set_labels
