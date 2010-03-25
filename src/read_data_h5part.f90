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
!  Copyright (C) 2005-2010 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
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
!  H5SPLASH_HSML=1.0   : value for global smoothing length if h not present in data
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

end module h5partdataread

subroutine read_data(rootname,indexstart,nstepsread)
  use particle_data,  only:dat,npartoftype,time,gamma,maxpart,maxcol,maxstep
  use params
  use settings_data,  only:ndim,ndimV,ncolumns,ncalc,debugmode
  use mem_allocation, only:alloc
  use iso_c_binding,  only:c_double,c_int64_t
  use asciiutils,     only:lcase
  use system_utils,   only:renvironment
  use labels,         only:ih,ipmass,irho
  use h5part
  use h5partdataread
  implicit none
  integer, intent(in)          :: indexstart
  integer, intent(out)         :: nstepsread
  character(len=*), intent(in) :: rootname
  integer :: j,ierr,ncolstep,nsteps,ncolsfile,nrealcols
  integer(kind=c_int64_t) :: nattrib,iattrib,nelem,icol,idatasettype,k,icolsfile
  integer :: nprint,npart_max,nstep_max,maxcolsfile
  logical :: iexist
  integer(kind=c_int64_t) :: ifile,istep
  character(len=len(rootname)+5) :: dumpfile
  character(len=64)              :: attribname
  character(len=lenlabel)        :: datasetname
  real(kind=doub_prec),dimension(1) :: dtime
  real :: hsmooth,hfac,dndim

  nstepsread = 0
  nstep_max = 0
  npart_max = maxpart

  dumpfile = trim(rootname)

  print "(1x,a)",'reading h5part format'
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
  ifile = h5pt_openr(trim(dumpfile))

  !if (h5pt_isvalid(ifile).ne.1) then
  if (ifile.le.0) then
     print "(a)",'*** ERROR opening '//trim(dumpfile)//' ***'
     return
  endif
  !print*,' ifile = ',ifile
  
  !
  !--get number of steps and particles in the file
  !
  nsteps = h5pt_getnsteps(ifile)
  if (debugmode) print*,'DEBUG: nsteps = ',nsteps

  !
  !--read "header" information from all steps in file:
  !  get maximum number of particles for all steps in file
  !  and maximum number of columns (datasets) in file
  !
  npart_max = 0
  ncolstep  = 0
  maxcolsfile = 0
  do istep=1,nsteps
     ierr = h5pt_setstep(ifile,istep)
     if (ierr.eq.0) then
        npart_max = max(npart_max,h5pt_getnpoints(ifile))
        ncolsfile = h5pt_getndatasets(ifile)
        icol = 0
        if ((istep.eq.nsteps .and. ncolsfile.gt.0)) then
           do icolsfile=0,ncolsfile-1
!              ierr = h5pt_getdatasetname(ifile,icol,datasetnames(icol+1))
!              print*,' data set name = ',icol,trim(datasetnames(icol+1))
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
                 case default
                    print "(a)",' skipping data set '//trim(datasetname)//' of type '//h5pt_type(idatasettype)
                 end select
              endif
           enddo
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
  !--call the set_labels routine to get ix, ih etc. from labels
  !
  call set_labels
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
           print "(/,a,f6.2,a,/)",' Setting smoothing length using h = ',hfac,&
                                  '*(m/rho)**(1/ndim) (from H5SPLASH_HFAC setting)'       
        else
           hfac = 1.2
           print "(/,a)",' WARNING: Smoothing length not found in data: using h = hfac*(m/rho)**(1/ndim)'
           print "(a)",  '          (hfac = 1.2 by default, set H5SPLASH_HFAC to change this)'
           print "(a,/)",'          (set constant h with H5SPLASH_HSML to give a global value)'
        endif
        ncolstep = ncolstep + 1
     else
        print "(/,a,/)",' WARNING: Smoothing length not found in data: Set H5SPLASH_HSML to give a global value'
     endif
  endif
  ncolumns = ncolstep

  !
  !--allocate memory for all data in the file
  !
  nstep_max = max(nsteps,indexstart,1)
  npart_max = max(maxpart,npart_max)
  if (.not.allocated(dat) .or. (nprint.gt.maxpart) .or. (ncolstep+ncalc).gt.maxcol) then
     call alloc(npart_max,nstep_max,ncolstep+ncalc)
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
              !--match anything that looks vaguely like the time
              if (nelem.eq.1 .and. (index(lcase(attribname),'time').ne.0  &
                 .or. index(lcase(attribname),'t ').ne.0)) then
                 ierr =h5pt_readstepattrib_r8(ifile,attribname,dtime)
                 if (ierr.eq.0) then
                    time(j) = real(dtime(1))
                    print "(12x,a,1pe10.3,a)",'time   = ',time(j),' (from '//trim(attribname)//')'
                 else
                    print "(a,i2,a)",' ERROR could not read time from step ',istep,' (from '//trim(attribname)//')'
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
    nrealcols = 0
    icol = 0
    do k=0,maxcolsfile-1
       ierr = h5pt_getdatasetinfo(ifile,k,datasetname,idatasettype,nelem)
       select case(idatasettype)
       case(H5PART_FLOAT32,H5PART_FLOAT64)
          icol = icol + 1
          datasetnames(icol) = trim(datasetname)
          !print "(a,i3,a)",' reading data set ',icol,': '//trim(datasetnames(icol))
          ierr = h5pt_readdata_r4(ifile,datasetnames(icol),dat(:,icol,j))
       case default
          !print "(a)",' skipping data set '//trim(datasetname)//' of type '//lcase(h5pt_type(idatasettype))
       end select
    enddo
!
!--if smoothing length set via environment variable, fill the extra column with the smoothing length value
!
    if (hsmooth.ge.0.) then
       dat(:,ih,j) = hsmooth
       datasetnames(ncolstep) = 'h'
    elseif (ipmass.gt.0 .and. irho.gt.0 .and. ndim.gt.0) then
       dndim = 1./ndim
       where (dat(:,irho,j).gt.tiny(0.))
         dat(:,ih,j) = hfac*(dat(:,ipmass,j)/dat(:,irho,j))**dndim
       elsewhere
         dat(:,ih,j) = 0.
       end where
       datasetnames(ncolstep) = 'h'
    endif
!    read(iunit,*,iostat=ierr) (dat(i,icol,j),icol = 1,ncolstep)
    nstepsread = nstepsread + 1

    npartoftype(:,j) = 0
    npartoftype(1,j) = nprint
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

subroutine set_labels
  use asciiutils,      only:lcase
  use labels,          only:label,labeltype,ix,irho,ipmass,ih,iutherm, &
                            ipr,ivx,iBfirst,iamvec,labelvec,lenlabel
  !use params,          only:maxparttypes
  use settings_data,   only:ncolumns,ntypes,ndim,ndimV,UseTypeInRenderings
  use geometry,        only:labelcoord
  use system_utils,    only:ienvironment
  use h5partdataread
  implicit none
  integer                 :: i,ierr,ndimVtemp,ndimset,ndim_max
  character(len=120)      :: columnfile
  character(len=lenlabel) :: labeli
  logical                 :: iexist
  
  ndim = 0
  ndimV = 0
  ndimset = ienvironment('H5SPLASH_NDIM',errval=-1)
  ndim_max = 3
  if (ndimset.ge.0) ndim_max = ndimset
  
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
     elseif ((index(labeli,'vel').ne.0 .or. labeli(1:1).eq.'v') .and. index(labeli,'_').ne.0) then
        if (ndimV.le.3) ndimV = ndimV + 1
        if (index(labeli,'_0').ne.0) ivx = i
     elseif (index(labeli,'dens').ne.0) then
        irho = i
     elseif (index(labeli,'mass').ne.0) then
        ipmass = i
     elseif (index(labeli,'smoothing').ne.0 .or. labeli(1:1).eq.'h') then
        ih = i
     elseif (labeli(1:1).eq.'u') then
        iutherm = i
     endif
  enddo

  if (ndim.lt.1) ndimV = 0
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

  if (ivx.gt.0) then
     iamvec(ivx:ivx+ndimV-1) = ivx
     labelvec(ivx:ivx+ndimV-1) = 'v'
     do i=1,ndimV
       label(ivx+i-1) = 'v\d'//labelcoord(i,1)
     enddo
  endif
  if (iBfirst.gt.0) then
     iamvec(iBfirst:iBfirst+ndimV-1) = ivx
     labelvec(iBfirst:iBfirst+ndimV-1) = 'B'
     do i=1,ndimV
       label(iBfirst+i-1) = 'B\d'//labelcoord(i,1)
     enddo
  endif
  !
  !--set labels for each particle type
  !
  ntypes = 1 !!maxparttypes
  labeltype(1) = 'gas'
  UseTypeInRenderings(1) = .true.
  
 
!-----------------------------------------------------------

  return 
end subroutine set_labels
