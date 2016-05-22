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

!-----------------------------------------------------------------
!  Module containing routines for converting 3D SPH dump
!  files to 3D gridded data.
!-----------------------------------------------------------------
module convert_grid
 private
 public :: convert_to_grid

contains

!-----------------------------------------------------------------
! interpolate 3D SPH data to grid and interface to grid
! data output routines
!-----------------------------------------------------------------
subroutine convert_to_grid(time,dat,ntypes,npartoftype,masstype,itype,ncolumns,filename,&
                           outformat,interpolateall)
 use labels,             only:label,labelvec,irho,ih,ipmass,ix,ivx,iBfirst
 use limits,             only:lim,get_particle_subset
 use settings_units,     only:units,unit_interp
 use settings_data,      only:ndim,ndimV,UseTypeInRenderings,iRescale,required,lowmemorymode,debugmode
 use settings_part,      only:iplotpartoftype
 use settings_render,    only:npix,inormalise_interpolations,idensityweightedinterpolation
 use params,             only:int1
 use interpolation,      only:set_interpolation_weights
 use interpolations3D,   only:interpolate3D,interpolate3D_vec
 use interpolations2D,   only:interpolate2D,interpolate2D_vec
 use system_utils,       only:lenvironment,renvironment,envlist,lenvstring,ienvstring
 use readwrite_griddata, only:open_gridfile_w,write_grid,write_gridlimits
 use particle_data,      only:icolourme
 use params,             only:int8
 implicit none
 integer, intent(in)                          :: ntypes,ncolumns
 integer, intent(in), dimension(:)            :: npartoftype
 integer(kind=int1), intent(in), dimension(:) :: itype
 real, intent(in)                             :: time
 real, intent(in), dimension(:,:)             :: dat
 real, intent(in), dimension(:)               :: masstype
 character(len=*), intent(in)                 :: filename,outformat
 logical, intent(in)                          :: interpolateall
 integer, parameter :: iunit = 89
 integer            :: ierr,i,k,ncolsgrid,ivec,nvec,iloc,j,nzero
 integer            :: npixx,ntoti,ninterp,nstring
 character(len=40)  :: fmtstring
 character(len=64)  :: fmtstring1

 real, dimension(:,:,:), allocatable   :: datgrid
 real, dimension(:,:),   allocatable   :: datgrid2D
 real, dimension(:,:,:,:), allocatable :: datgridvec
 real, dimension(:,:,:),   allocatable :: datgridvec2D
 real, dimension(:), allocatable       :: weight
 real, dimension(3)    :: xmin,xmax
 real, dimension(3)    :: partmin,partmax,partmean
 real, dimension(3)    :: datmin,datmax,datmean
 integer, dimension(3) :: npixels
 integer(kind=int8), dimension(3) :: npixels8
 integer, dimension(12) :: icoltogrid
 integer :: ncolstogrid,icol
 real    :: hmin,pixwidth,rhominset,rhomin,gridmin,gridmax,gridmean
 logical :: inormalise,lowmem
 logical, dimension(3) :: isperiodic
 character(len=30), dimension(12) :: strings
 character(len=1), dimension(3), parameter :: xlab = (/'x','y','z'/)

 !
 !--check for errors in input settings
 !
 if (ndim.lt.2 .or. ndim.gt.3) then
    print "(/,a,i2,a,/)",' ERROR: SPH data has ',ndim,' spatial dimensions: cannot convert to 3D grid'
    return
 endif

 print "(/,'----->',1x,a,i1,a,/)",'CONVERTING SPH DATA -> ',ndim,'D GRID'
 xmin(1:ndim) = lim(ix(1:ndim),1)
 xmax(1:ndim) = lim(ix(1:ndim),2)
 !
 !--print limits information
 !
 call write_gridlimits(ndim,xmin,xmax,label(ix(1:ndim)))

 !
 !--SPLASH_TO_GRID can be set to comma separated list of columns
 !  in order to select particular quantities for interpolation to grid
 !
 ncolstogrid   = 0
 icoltogrid(:) = 0
 call envlist('SPLASH_TO_GRID',nstring,strings)
 if (nstring.gt.0) then
    do i=1,nstring
       icol = ienvstring(strings(i))
       if (ienvstring(strings(i)).gt.0) then
          ncolstogrid = ncolstogrid + 1
          icoltogrid(ncolstogrid) = icol
       endif
    enddo
 endif
 !
 !--for backwards compatibility, support the SPLASH_TO_GRID_DENSITY_ONLY option
 !  but only if SPLASH_TO_GRID is not set
 !
 if (ncolstogrid.eq.0 .and. lenvironment('SPLASH_TO_GRID_DENSITY_ONLY')) then
    ncolstogrid = 1
    icoltogrid(1) = irho
 endif
 
 !
 !--whether or not to wrap particle contributions across boundaries
 !
 isperiodic(:) = .false.
 call envlist('SPLASH_TO_GRID_PERIODIC',nstring,strings)
 if (nstring.gt.ndim) then
    print "(a)",' ERROR in SPLASH_TO_GRID_PERIODIC setting'
    nstring = ndim
 endif
 do i=1,nstring
    isperiodic(i) = lenvstring(strings(i))
 enddo
 if (nstring.eq.1) isperiodic(2:ndim) = isperiodic(1)

 if (all(isperiodic(1:ndim))) then
    print "(/,a)",' using PERIODIC boundaries (from SPLASH_TO_GRID_PERIODIC setting)'
 elseif (isperiodic(1) .or. isperiodic(2) .or. isperiodic(3)) then
    print*
    do i=1,ndim
       if (isperiodic(i)) then
          print "(a)",' using PERIODIC boundaries in '//xlab(i)//' (from SPLASH_TO_GRID_PERIODIC setting)'
       else
          print "(a)",' using NON-PERIODIC bounds in '//xlab(i)//' (from SPLASH_TO_GRID_PERIODIC setting)'    
       endif
    enddo
 else
    print "(/,a)",' using NON-PERIODIC boundaries'
    print "(a)",' (set SPLASH_TO_GRID_PERIODIC=yes for periodic'
    if (ndim.eq.3) then
       print "(a)",'   or SPLASH_TO_GRID_PERIODIC=yes,no,yes for mixed)'    
    else
       print "(a)",'   or SPLASH_TO_GRID_PERIODIC=yes,no for mixed)'
    endif
 endif

 ierr = 0
 do i=1,ndim
    if ((xmax(i)-xmin(i)).lt.tiny(0.)) then
       print "(a)",' ERROR: min = max in '//trim(label(ix(i)))//&
                   ' coordinate: cannot interpolate to zero-sized grid!'
       ierr = 1
    endif
 enddo

 if (irho.le.0 .or. irho.gt.ncolumns) then
    print "(a)",' ERROR: density not found in data read.'
    ierr = 2
 endif
 if (ih.le.0 .or. ih.gt.ncolumns) then
    print "(a)",' ERROR: smoothing length not found in data read.'
    ierr = 3
 endif
 if (ipmass.le.0 .or. ipmass.gt.ncolumns) then
    if (all(masstype(:).lt.tiny(0.))) then
       print "(a)",' ERROR: particle masses not read as column, and mass per type not set.'
       ierr = 4
    endif
 endif

 if (ierr /= 0) then
    print "(a,i1,a)",' cannot perform SPH interpolation to ',ndim,'D grid, skipping file...'
    return
 endif
 ierr = 0

 !
 !--set number of particles to use in the interpolation routines
 !  (by default, only the gas particles)
 !
 ntoti = sum(npartoftype)
 ninterp = npartoftype(1)
 if (any(UseTypeInRenderings(2:ntypes).and.iplotpartoftype(2:ntypes)) &
     .or. size(itype).gt.1) ninterp = ntoti

 allocate(weight(ninterp),stat=ierr)
 if (ierr /= 0) then
    print*,' ERROR allocating memory for interpolation weights, aborting...'
    return
 endif
 !
 !--set interpolation weights (w = m/(rho*h^ndim)
 !
 inormalise = inormalise_interpolations
 call set_interpolation_weights(weight,dat,itype,(iplotpartoftype .and. UseTypeInRenderings),&
      ninterp,npartoftype,masstype,ntypes,ncolumns,irho,ipmass,ih,ndim,iRescale,&
      idensityweightedinterpolation,inormalise,units,unit_interp,required,.false.)
 !
 !--set colours (just in case)
 !
 icolourme(:) = 1
 !
 !--apply range restrictions to data
 !
 call get_particle_subset(icolourme,dat,ncolumns)

 !
 !--work out how many pixels to use
 !
 npixx = npix
 if (npixx.le.0) then
    print "(/,a)",' WARNING: number of pixels = 0, using automatic pixel numbers'
    hmin = 0.
    call minmaxmean_part(dat(:,ih:ih),weight,ninterp,partmin,partmax,partmean,nonzero=.true.)
    hmin = partmin(1)
    if (hmin.gt.0.) then
       print*,'based on the minimum smoothing length of hmin = ',hmin
       npixels8(1:ndim) = int((xmax(1:ndim) - xmin(1:ndim))/hmin,kind=int8) + 1
       if (ndim.eq.3) then
          print "(a,i6,2(' x',i6),a)",' requires ',npixels8(1:ndim),' pixels to capture the full resolution'       
          if (product(npixels8(1:ndim)).gt.512**3 .or. product(npixels8(1:ndim)).le.0) then
             npixx = 512
             print "(a,i4)",' but this is ridiculous, so instead we choose ',npixx
          else
             npixx = npixels8(1)
          endif
       else
          print "(a,i6,1(' x',i6),a)",' requires ',npixels8(1:ndim),' pixels to capture the full resolution'
          if (product(npixels8(1:ndim)).gt.1024**ndim .or. product(npixels8(1:ndim)).le.0) then
             npixx = 1024
             print "(a,i4)",' but this is very large, so instead we choose ',npixx
          else
             npixx = npixels8(1)
          endif
       endif
    else
       npixx = 512
       print "(a)",' ...but cannot get auto pixel number because hmin = 0'
       print "(a)",'    so instead we choose npixels = ',npixx
    endif
 endif
 print*
 pixwidth        = (xmax(1)-xmin(1))/npixx
 npixels(1:ndim) = int((xmax(1:ndim)-xmin(1:ndim) - 0.5*pixwidth)/pixwidth) + 1

 !
 !--work out how many columns will be written to file
 !
 nvec = 0
 if (ncolstogrid.gt.0) then
    ncolsgrid = ncolstogrid
 elseif (interpolateall) then
    ncolsgrid = 0
    do i=1,ncolumns
       if (.not.any(ix(1:ndim).eq.i) .and. i.ne.ih .and. i.ne.ipmass) then
          ncolsgrid = ncolsgrid + 1
       endif
    enddo
 else
    if (ndimV.eq.ndim) then
       if (ivx.gt.0 .and. ivx+ndimV-1.le.ncolumns) nvec = nvec + 1
       if (iBfirst.gt.0 .and. iBfirst+ndimV-1.le.ncolumns) nvec = nvec + 1
    endif
    ncolsgrid = 1 + ndimV*nvec
 endif

 !
 !--use low memory mode for large grids
 !
 if (trim(outformat)=='gridascii2') then
    lowmem = .false.
 else
    lowmem = .true.
 endif

! if ((ndim.eq.3 .and. product(npixels).gt.256**3)) then
!    lowmem = .true.
! else
!    lowmem = lowmemorymode
! endif
 if (lowmem .and. nvec.gt.0) &
    print "(a,/)",' [doing velocity field components separately (low memory mode)]'
 !
 !--allocate memory for the grid
 !
 if (allocated(datgrid))   deallocate(datgrid)
 if (allocated(datgrid2D)) deallocate(datgrid2D)

 if (ndim.eq.3) then
    write(*,"(a,i5,2(' x',i5),a)",advance='no') ' >>> allocating memory for ',npixels(1:ndim),' grid ...'
    allocate(datgrid(npixels(1),npixels(2),npixels(3)),stat=ierr)
 elseif (ndim.eq.2) then
    write(*,"(a,i5,1(' x',i5),a)",advance='no') ' >>> allocating memory for ',npixels(1:ndim),' grid ...'
    allocate(datgrid2D(npixels(1),npixels(2)),stat=ierr)
 endif
 if (ierr /= 0) then
    write(*,*) 'FAILED: NOT ENOUGH MEMORY'
    if (allocated(weight)) deallocate(weight)
    return
 else
    write(*,*) 'OK'
 endif

 !
 !--open grid file for output (also checks format is OK)
 !
 call open_gridfile_w(iunit,filename,outformat,ndim,ncolsgrid,npixels(1:ndim),time,ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR: could not open grid file for output, skipping...'
    if (allocated(datgrid)) deallocate(datgrid)
    if (allocated(datgrid)) deallocate(datgrid2D)
    if (allocated(weight)) deallocate(weight)
    return
 endif

 fmtstring1 = "(12x,a20,1x,'    min    ',1x,'    max    ',1x,'    mean    ')"
 fmtstring = "(22x,a10,1x,3(es10.2,1x))"
 !
 !--interpolate density to the 3D grid
 !
 print "(/,a,i1,a)",' interpolating density to ',ndim,'D grid...'
 if (debugmode) print*,'DEBUG: density in column ',irho,' vals = ',dat(1:10,irho)

 call minmaxmean_part(dat(1:ninterp,irho:irho),weight,ninterp,partmin,partmax,partmean)
 print fmtstring1,trim(label(irho))
 print fmtstring,' on parts:',partmin(1),partmax(1),partmean(1)

 if (ndim.eq.3) then
    call interpolate3D(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)),&
         dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,irho),icolourme,ninterp,&
         xmin(1),xmin(2),xmin(3),datgrid,npixels(1),npixels(2),npixels(3),&
         pixwidth,pixwidth,inormalise,isperiodic(1),isperiodic(2),isperiodic(3))
    !
    !--set minimum density on the grid
    !
    call minmaxmean_grid(datgrid,npixels,gridmin,gridmax,gridmean,nonzero=.true.)
 else
    call interpolate2D(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),&
         dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,irho),icolourme,ninterp,&
         xmin(1),xmin(2),datgrid2D,npixels(1),npixels(2),&
         pixwidth,pixwidth,inormalise,isperiodic(1),isperiodic(2))
    !
    !--set minimum density on the grid
    !
    call minmaxmean_grid2D(datgrid2D,npixels,gridmin,gridmax,gridmean,nonzero=.true.) 
 endif

 print fmtstring1,trim(label(irho))
 print fmtstring,' on grid :',gridmin,gridmax,gridmean

 rhomin = gridmin
 rhominset = renvironment('SPLASH_TO_GRID_RHOMIN',errval=-1.)

 print*
 if (rhominset.ge.0.) then
    rhomin = rhominset
    print*,'enforcing minimum density on grid = ',rhomin
    print*,'(based on SPLASH_TO_GRID_RHOMIN setting)'
 elseif (rhomin.gt.0.) then
    print*,'enforcing minimum density on grid = ',rhomin
    print*,'set SPLASH_TO_GRID_RHOMIN=minval to manually set this (e.g. to zero)'
 endif

 if (rhomin.gt.0.) then
    nzero = 0
    if (ndim.eq.3) then
       !$omp parallel do private(k,j,i) reduction(+:nzero) schedule(static)
       do k=1,npixels(3)
          do j=1,npixels(2)
             do i=1,npixels(1)
                if (datgrid(i,j,k).le.tiny(datgrid)) then
                   datgrid(i,j,k) = rhomin
                   nzero = nzero + 1
                endif
             enddo
          enddo
       enddo
    else
       !$omp parallel do private(j,i) reduction(+:nzero) schedule(static)
       do j=1,npixels(2)
          do i=1,npixels(1)
             if (datgrid2D(i,j).le.tiny(datgrid2D)) then
                datgrid2D(i,j) = rhomin
                nzero = nzero + 1
             endif
          enddo
       enddo
    endif
    print "(a,i8,a)",' minimum density enforced on ',nzero,' grid cells'
 else
    print*,'minimum density NOT enforced'
 endif

 !
 !--write density to grid data file
 !
 print*
 if (ndim.eq.3) then
    if (lowmem .or. interpolateall .or. ncolstogrid.gt.0) then
       call write_grid(iunit,filename,outformat,ndim,1,npixels,trim(label(irho)),&
                       time,pixwidth,xmin,ierr,dat=datgrid)
    endif
 else
    call write_grid(iunit,filename,outformat,ndim,1,npixels,trim(label(irho)),&
                    time,pixwidth,xmin,ierr,dat2D=datgrid2D)
 endif
 !
 !--interpolate remaining quantities to the 3D grid
 !
 if (interpolateall .or. ncolstogrid.gt.0) then
    if (ncolstogrid.gt.0) then
       print "(/,a,i2,a)",' Interpolating ',ncolstogrid,' columns to grid from SPLASH_TO_GRID setting:'
       print "(' got SPLASH_TO_GRID=',10(i2,1x))",icoltogrid(1:ncolstogrid)
    endif

    do i=1,ncolumns
       if ((ncolstogrid.gt.0 .and. any(icoltogrid.eq.i) .and. i.ne.irho) .or.  &
           (interpolateall .and. &
            .not.any(ix(:).eq.i) .and. i.ne.ih .and. i.ne.ipmass .and. i.ne.irho)) then

          print "(/,a)",' interpolating '//trim(label(i))
          print fmtstring1,trim(label(i))
          call minmaxmean_part(dat(1:ninterp,i:i),weight,ninterp,partmin,partmax,partmean)
          print fmtstring,' on parts:',partmin(1),partmax(1),partmean(1)

          if (iszero(partmin,partmax,1)) then
             datgrid = 0.
          else
             if (ndim.eq.3) then
                call interpolate3D(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)),&
                     dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,i),icolourme,ninterp,&
                     xmin(1),xmin(2),xmin(3),datgrid,npixels(1),npixels(2),npixels(3),&
                     pixwidth,pixwidth,.true.,isperiodic(1),isperiodic(2),isperiodic(3))             
             else
                call interpolate2D(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),&
                     dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,i),icolourme,ninterp,&
                     xmin(1),xmin(2),datgrid2D,npixels(1),npixels(2),&
                     pixwidth,pixwidth,.true.,isperiodic(1),isperiodic(2))
             endif
          endif

          !
          !--write gridded data to file
          !
          if (ndim.eq.3) then
             call minmaxmean_grid(datgrid,npixels,gridmin,gridmax,gridmean,.false.)
             print fmtstring,' on grid :',gridmin,gridmax,gridmean
             call write_grid(iunit,filename,outformat,ndim,1,npixels,trim(label(i)),&
                  time,pixwidth,xmin,ierr,dat=datgrid)
          else
             call minmaxmean_grid2D(datgrid2D,npixels,gridmin,gridmax,gridmean,.false.)          
             print fmtstring,' on grid :',gridmin,gridmax,gridmean
             call write_grid(iunit,filename,outformat,ndim,1,npixels,trim(label(i)),&
                  time,pixwidth,xmin,ierr,dat2D=datgrid2D)
          endif
       endif
    enddo

 else

    if (nvec.gt.0) then

       print "(/,a,i2,a)",' set SPLASH_TO_GRID=',irho,' to interpolate density ONLY and skip remaining columns'
       print "(a,i2,a)",  '     SPLASH_TO_GRID=6,8,10 to select particular columns'

       if (.not.lowmem) then
          if (ndim.eq.3) then
             write(*,"(/,a,i5,2(' x',i5),a)",advance='no') ' >>> allocating memory for ',npixels(1:ndim),' x 3 grid ...'
             allocate(datgridvec(3,npixels(1),npixels(2),npixels(3)),stat=ierr)
          else
             write(*,"(/,a,i5,1(' x',i5),a)",advance='no') ' >>> allocating memory for ',npixels(1:ndim),' x 3 grid ...'
             allocate(datgridvec2D(2,npixels(1),npixels(2)),stat=ierr)
          endif
          if (ierr /= 0) then
             write(*,*) 'FAILED: NOT ENOUGH MEMORY'
             return
          else
             write(*,*) 'OK'
          endif
       endif

       !
       !--interpolate velocity field and magnetic fields and write to file
       !
       over_vec: do ivec=1,nvec
          select case(ivec)
          case(1)
             iloc = ivx
             print "(a,i1,a)",' interpolating velocity field to ',ndim,'D grid...'
          case(2)
             iloc = iBfirst
             print "(a,i1,a)",' interpolating magnetic field to ',ndim,'D grid...'
          case default
             iloc = 0
             exit over_vec
          end select
          if (iloc.le.0 .or. iloc.ge.ncolumns) cycle over_vec

          if (lowmem) then

             do i=iloc,iloc+ndimV-1
                print "(/,a)",' interpolating '//trim(label(i))
                print fmtstring1,trim(label(i))
                call minmaxmean_part(dat(1:ninterp,i:i),weight,ninterp,partmin,partmax,partmean)
                print fmtstring,' on parts:',partmin(1),partmax(1),partmean(1)

                if (iszero(partmin,partmax,1)) then
                   datgrid = 0.
                else
                   if (ndim.eq.3) then
                      call interpolate3D(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)),&
                           dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,i),icolourme,ninterp,&
                           xmin(1),xmin(2),xmin(3),datgrid,npixels(1),npixels(2),npixels(3),&
                           pixwidth,pixwidth,.true.,isperiodic(1),isperiodic(2),isperiodic(3))
                   else
                      call interpolate2D(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),&
                           dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,i),icolourme,ninterp,&
                           xmin(1),xmin(2),datgrid2D,npixels(1),npixels(2),&
                           pixwidth,pixwidth,.true.,isperiodic(1),isperiodic(2))
                   endif
                endif

                !
                !--write gridded data to file
                !
                if (ndim.eq.3) then
                   call minmaxmean_grid(datgrid,npixels,gridmin,gridmax,gridmean,.false.)
                   print fmtstring,' on grid :',gridmin,gridmax,gridmean
                   call write_grid(iunit,filename,outformat,ndim,1,npixels,trim(label(i)),&
                        time,pixwidth,xmin,ierr,dat=datgrid)
                else
                   call minmaxmean_grid2D(datgrid2D,npixels,gridmin,gridmax,gridmean,.false.)                
                   print fmtstring,' on grid :',gridmin,gridmax,gridmean
                   call write_grid(iunit,filename,outformat,ndim,1,npixels,trim(label(i)),&
                        time,pixwidth,xmin,ierr,dat2D=datgrid2D)
                endif
             enddo
          else

             print fmtstring1,trim(labelvec(iloc))
             call minmaxmean_part(dat(1:ninterp,iloc:iloc+ndimV-1),weight,ninterp,partmin,partmax,partmean)
             do i=1,ndimV
                print fmtstring,' on parts:',partmin(i),partmax(i),partmean(i)
             enddo

             if (iszero(partmin,partmax,ndimV)) then
                datgridvec = 0.
             else
                if (ndim.eq.3) then
                   call interpolate3D_vec(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)),&
                        dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,iloc:iloc+ndimV-1),icolourme,ninterp,&
                        xmin(1),xmin(2),xmin(3),datgridvec,npixels(1),npixels(2),npixels(3),&
                        pixwidth,pixwidth,.true.,isperiodic(1),isperiodic(2),isperiodic(3))                
                else
                   call interpolate2D_vec(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),&
                        dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,iloc),dat(1:ninterp,iloc+1), &
                        icolourme,ninterp,xmin(1),xmin(2),datgridvec2D(1,:,:),datgridvec2D(2,:,:), &
                        npixels(1),npixels(2),pixwidth,pixwidth,.true.,isperiodic(1),isperiodic(2))
                endif
             endif

             if (ndim.eq.3) then
                call minmaxmean_gridvec(datgridvec,npixels,ndimV,datmin,datmax,datmean)
             else
                call minmaxmean_gridvec2D(datgridvec2D,npixels,ndimV,datmin,datmax,datmean)             
             endif
             do i=1,ndimV
                print fmtstring,' on grid :',datmin(i),datmax(i),datmean(i)
             enddo
             !
             !--write result to grid file
             !
             if (ndim.eq.3) then
                call write_grid(iunit,filename,outformat,ndim,ndimV,npixels,&
                                label(irho),time,pixwidth,xmin,ierr,&
                                dat=datgrid,dat3D=datgridvec,label3D=label(iloc:iloc+ndimV))                
             else
                do i=1,ndimV
                   call write_grid(iunit,filename,outformat,ndim,ndimV,npixels,&
                                   label(iloc+i-1),time,pixwidth,xmin,ierr,dat2D=datgridvec2D(i,:,:))
                enddo
             endif
          endif
          print*
       enddo over_vec
    endif
! else
!    print "(/,a)",' skipping remaining quantities (from SPLASH_TO_GRID_DENSITY_ONLY setting)'
 endif

 close(iunit)

 if (allocated(datgrid))      deallocate(datgrid)
 if (allocated(datgrid2D))    deallocate(datgrid2D)
 if (allocated(datgridvec))   deallocate(datgridvec)
 if (allocated(datgridvec2D)) deallocate(datgridvec2D)
 if (allocated(weight)) deallocate(weight)

 return
end subroutine convert_to_grid

!-----------------------------------------------
! calculate max and min and mean values on grid
!-----------------------------------------------
subroutine minmaxmean_grid(datgrid,npixels,gridmin,gridmax,gridmean,nonzero)
 implicit none
 real, dimension(:,:,:), intent(in) :: datgrid
 integer, dimension(3), intent(in)  :: npixels
 real, intent(out)                  :: gridmin,gridmax,gridmean
 logical, intent(in)                :: nonzero
 real    :: dati
 integer :: i,j,k

 gridmax = -huge(gridmax)
 gridmin = huge(gridmin)
 gridmean = 0.
 !$omp parallel do schedule(static) &
 !$omp reduction(min:gridmin) &
 !$omp reduction(max:gridmax) reduction(+:gridmean) &
 !$omp private(k,j,i,dati)
 do k=1,npixels(3)
    do j=1,npixels(2)
       do i=1,npixels(1)
          dati = datgrid(i,j,k)
          gridmax  = max(gridmax,dati)
          if (nonzero) then
             if (dati.gt.tiny(0.)) gridmin  = min(gridmin,dati)
          else
             gridmin  = min(gridmin,dati)
          endif
          gridmean = gridmean + dati
       enddo
    enddo
 enddo
 gridmean = gridmean/product(npixels(1:3))

 return
end subroutine minmaxmean_grid

!----------------------------------------------------
! calculate max and min and mean values on grid (2D)
!----------------------------------------------------
subroutine minmaxmean_grid2D(datgrid,npixels,gridmin,gridmax,gridmean,nonzero)
 implicit none
 real, dimension(:,:), intent(in)   :: datgrid
 integer, dimension(2), intent(in)  :: npixels
 real, intent(out)                  :: gridmin,gridmax,gridmean
 logical, intent(in)                :: nonzero
 real    :: dati
 integer :: i,j

 gridmax = -huge(gridmax)
 gridmin = huge(gridmin)
 gridmean = 0.
 !$omp parallel do schedule(static) &
 !$omp reduction(min:gridmin) &
 !$omp reduction(max:gridmax) reduction(+:gridmean) &
 !$omp private(j,i,dati)
 do j=1,npixels(2)
    do i=1,npixels(1)
       dati = datgrid(i,j)
       gridmax  = max(gridmax,dati)
       if (nonzero) then
          if (dati.gt.tiny(0.)) gridmin  = min(gridmin,dati)
       else
          gridmin  = min(gridmin,dati)
       endif
       gridmean = gridmean + dati
    enddo
 enddo
 gridmean = gridmean/product(npixels(1:2))

 return
end subroutine minmaxmean_grid2D

!-----------------------------------------------
! calculate max and min and mean values on grid
! (for vector quantities)
!-----------------------------------------------
subroutine minmaxmean_gridvec(datgridvec,npixels,jlen,gridmin,gridmax,gridmean)
 implicit none
 real, dimension(:,:,:,:), intent(in) :: datgridvec
 integer, dimension(3), intent(in)    :: npixels
 integer, intent(in)                  :: jlen
 real, dimension(jlen), intent(out)   :: gridmin,gridmax,gridmean
 real    :: dati
 integer :: ivec,i,j,k

 gridmax(:)  = -huge(gridmax)
 gridmin(:)  = huge(gridmin)
 gridmean(:) = 0.
 !$omp parallel do schedule(static) &
 !$omp reduction(min:gridmin) &
 !$omp reduction(max:gridmax) reduction(+:gridmean) &
 !$omp private(k,j,i,dati)
 do k=1,npixels(3)
    do j=1,npixels(2)
       do i=1,npixels(1)
          do ivec=1,jlen
             dati = datgridvec(ivec,i,j,k)
             gridmax(ivec)  = max(gridmax(ivec),dati)
             gridmin(ivec)  = min(gridmin(ivec),dati)
             gridmean(ivec) = gridmean(ivec) + dati
          enddo
       enddo
    enddo
 enddo
 gridmean(1:jlen) = gridmean(1:jlen)/real(product(npixels(1:3)))

 return
end subroutine minmaxmean_gridvec

!-----------------------------------------------
! calculate max and min and mean values on grid
! (for vector quantities)
!-----------------------------------------------
subroutine minmaxmean_gridvec2D(datgridvec,npixels,jlen,gridmin,gridmax,gridmean)
 implicit none
 real, dimension(:,:,:), intent(in) :: datgridvec
 integer, dimension(2), intent(in)  :: npixels
 integer, intent(in)                :: jlen
 real, dimension(jlen), intent(out) :: gridmin,gridmax,gridmean
 real    :: dati
 integer :: ivec,i,j

 gridmax(:)  = -huge(gridmax)
 gridmin(:)  = huge(gridmin)
 gridmean(:) = 0.
 !$omp parallel do schedule(static) &
 !$omp reduction(min:gridmin) &
 !$omp reduction(max:gridmax) reduction(+:gridmean) &
 !$omp private(j,i,dati)
 do j=1,npixels(2)
    do i=1,npixels(1)
       do ivec=1,jlen
          dati = datgridvec(ivec,i,j)
          gridmax(ivec)  = max(gridmax(ivec),dati)
          gridmin(ivec)  = min(gridmin(ivec),dati)
          gridmean(ivec) = gridmean(ivec) + dati
       enddo
    enddo
 enddo
 gridmean(1:jlen) = gridmean(1:jlen)/real(product(npixels(1:2)))

 return
end subroutine minmaxmean_gridvec2D

!----------------------------------------------------
! calculate max and min and mean values on particles
!----------------------------------------------------
subroutine minmaxmean_part(dat,weight,npart,partmin,partmax,partmean,nonzero)
 implicit none
 real, dimension(:,:), intent(in)   :: dat
 real, dimension(:), intent(in)     :: weight
 integer, intent(in)                :: npart
 real, dimension(3), intent(out)    :: partmin,partmax,partmean
 logical, intent(in), optional      :: nonzero
 real    :: partval
 integer :: np,jlen,i,j
 logical :: usenonzero

 usenonzero = .false.
 if (present(nonzero)) usenonzero = nonzero

 partmax(:)  = -huge(partmax)
 partmin(:)  = huge(partmin)
 partmean(:) = 0.
 np          = 0
 jlen        = min(size(dat(1,:)),3)

 !--could do this in parallel but reduction on arrays
 !  does not seem to work in ifort
 !!$omp parallel do default(none) schedule(static) &
 !!$omp shared(dat,weight,jlen,npart,usenonzero) &
 !!$omp reduction(min:partmin) &
 !!$omp reduction(max:partmax) &
 !!$omp reduction(+:partmean,np) &
 !!$omp private(i,j,partval)
 do i=1,npart
    !--only count particles used in the rendering
    if (weight(i).gt.tiny(0.)) then
       np = np + 1
       do j=1,jlen
          partval = dat(i,j)
          if (usenonzero) then
             if (partval.gt.tiny(0.)) partmin(j) = min(partmin(j),partval)
          else
             partmin(j) = min(partmin(j),partval)
          endif
          partmax(j) = max(partmax(j),partval)
          partmean(j) = partmean(j) + partval
       enddo
    endif
 enddo
 !!$omp end parallel do

 if (np.gt.0) then
    partmean(:) = partmean(:)/real(np)
 endif

 return
end subroutine minmaxmean_part

!----------------------------------------------------
! calculate max and min and mean values on particles
!----------------------------------------------------
logical function iszero(partmin,partmax,ndim)
 implicit none
 real, dimension(:), intent(in) :: partmin,partmax
 integer, intent(in)            :: ndim

 if (all(abs(partmin(1:ndim)).lt.tiny(0.)) .and. &
     all(abs(partmax(1:ndim)).lt.tiny(0.))) then
    iszero = .true.
    print "(a)",' min=max=0 on particles: skipping pointless interpolation and setting dat = 0.'
 else
    iszero = .false.
 endif

end function iszero

end module convert_grid
