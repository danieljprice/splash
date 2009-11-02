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
!  Copyright (C) 2005-2009 Daniel Price. All rights reserved.
!  Contact: daniel.price@sci.monash.edu.au
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
subroutine convert_to_grid(time,dat,npart,ntypes,npartoftype,masstype,itype,ncolumns,filename,outformat,&
                           interpolateall)
 use labels,          only:label,labelvec,irho,ih,ipmass,ix,ivx,iBfirst
 use limits,          only:lim
 use settings_units,  only:units
 use settings_data,   only:ndim,ndimV,ndataplots,UseTypeInRenderings,iRescale,required
 use settings_part,   only:iplotpartoftype
 use settings_render, only:npix,inormalise_interpolations,idensityweightedinterpolation
 use params,          only:int1
 use interpolation,   only:set_interpolation_weights
 use interpolations3D,only:interpolate3D,interpolate3D_vec
 use system_utils,    only:lenvironment,renvironment
 use write_griddata,  only:open_gridfile,write_grid
 use particle_data,   only:icolourme
 implicit none
 integer, intent(in)                          :: npart,ntypes,ncolumns
 integer, intent(in), dimension(:)            :: npartoftype
 integer(kind=int1), intent(in), dimension(:) :: itype
 real, intent(in)                             :: time
 real, intent(in), dimension(npart,ncolumns)  :: dat
 real, intent(in), dimension(:)               :: masstype
 character(len=*), intent(in)                 :: filename,outformat
 logical, intent(in)                          :: interpolateall
 integer, parameter :: iunit = 89
 integer            :: ierr,i,k,ncolsgrid,ivec,nvec,iloc
 integer            :: npixx,ntoti,ninterp
 character(len=40)  :: fmtstring
 character(len=64)  :: fmtstring1

 real, dimension(:,:,:), allocatable   :: datgrid
 real, dimension(:,:,:,:), allocatable :: datgridvec
 real, dimension(:), allocatable       :: weight
 real, dimension(3)    :: xmin,xmax
 real, dimension(3)    :: partmin,partmax,partmean
 real, dimension(3)    :: datmin,datmax,datmean
 integer, dimension(3) :: npixels
 real    :: hmin,pixwidth,rhominset,rhomin,gridmin,gridmax,gridmean
 logical :: isperiodic,inormalise
 
 !
 !--check for errors in input settings
 !
 if (ndim.ne.3) then
    print "(/,a,/)",' ERROR: SPH data has < 3 spatial dimensions: cannot convert to 3D grid'
    return
 endif

 print "(/,'----->',1x,a,/)",'CONVERTING SPH DATA -> 3D GRID'
 xmin(:) = lim(ix(:),1)
 xmax(:) = lim(ix(:),2)
 !
 !--print limits information
 !
 print "(a)",' grid dimensions:'
 do i=1,3
    if (maxval(abs(xmax)).lt.1.e7) then
       print "(1x,a,': ',f8.2,' -> ',f8.2)",trim(label(ix(i))),xmin(i),xmax(i)
    else
       print "(1x,a,': ',es10.2,' -> ',es10.2)",trim(label(ix(i))),xmin(i),xmax(i)
    endif
 enddo
 !
 !--whether or not to wrap particle contributions across boundaries
 !
 isperiodic = lenvironment('SPLASH_TO_GRID_PERIODIC')
 if (isperiodic) then
    print "(/,a)",' using PERIODIC boundaries (from SPLASH_TO_GRID_PERIODIC setting)'
 else
    print "(/,a)",' using NON-PERIODIC boundaries'
    print "(a)",' (set SPLASH_TO_GRID_PERIODIC=yes for periodic)'
 endif

 ierr = 0
 do i=1,ndim
    if ((xmax(i)-xmin(i)).lt.tiny(0.)) then
       print "(a)",' ERROR: min = max in '//trim(label(ix(i)))//&
                   ' coordinate: cannot interpolate to zero-sized grid!'
       ierr = 1
    endif
 enddo
 
 if (irho.le.0 .or. irho.gt.ndataplots) then
    print "(a)",' ERROR: density not found in data read.'
    ierr = 2
 endif
 if (ih.le.0 .or. ih.gt.ndataplots) then
    print "(a)",' ERROR: smoothing length not found in data read.'
    ierr = 3
 endif
 if (ipmass.le.0 .or. ipmass.gt.ndataplots) then
    if (all(masstype(:).lt.tiny(0.))) then
       print "(a)",' ERROR: particle masses not read as column, and mass per type not set.'
       ierr = 4
    endif
 endif
 
 if (ierr /= 0) then
    print "(a)",' cannot perform SPH interpolation to 3D grid, skipping file...'
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
      ninterp,npartoftype,masstype,ntypes,ndataplots,irho,ipmass,ih,ndim,iRescale,&
      idensityweightedinterpolation,inormalise,units,required)
 !
 !--set colours (just in case)
 !
 icolourme(:) = 1

 !
 !--work out how many pixels to use
 !
 npixx = npix
 if (npixx.le.0) then
    print "(a)",' WARNING: number of pixels = 0, using automatic pixel numbers'
    hmin = 0.
    call minmaxmean_part(dat(:,ih:ih),weight,ninterp,partmin,partmax,partmean,nonzero=.true.)
    hmin = partmin(1)
    if (hmin.gt.0.) then
       print*,'based on the minimum smoothing length of hmin = ',hmin
       npixels(:) = int(xmax(:) - xmin(:))/hmin + 1
       print "(a,i6,2(' x',i6),a)",' requires ',npixels(:),' pixels to capture the full resolution'
       if (product(npixels(:)).gt.huge(0) .or. product(npixels(:)).le.0) then
          npixx = 512
          print "(a,i4)",' but this is ridiculous, so instead we choose ',npixx
       else
          npixx = npixels(1)
       endif
    else
       npixx = 512
       print "(a)",' ...but cannot get auto pixel number because hmin = 0'
       print "(a)",'    so instead we choose npixels = ',npixx
    endif
 endif
 pixwidth    = (xmax(1)-xmin(1))/npixx
 npixels(:)  = int((xmax(:)-xmin(:) - epsilon(0.))/pixwidth) + 1

 !
 !--work out how many columns will be written to file
 !
 if (interpolateall) then
    ncolsgrid = 0
    do i=1,ndataplots
       if (.not.any(ix(:).eq.i) .and. i.ne.ih .and. i.ne.ipmass) then
          ncolsgrid = ncolsgrid + 1
       endif
    enddo
    nvec = 0
 else
    nvec = 0
    if (ndimV.eq.3) then
       if (ivx.gt.0 .and. ivx+ndimV-1.le.ndataplots) nvec = nvec + 1
       if (iBfirst.gt.0 .and. iBfirst+ndimV-1.le.ndataplots) nvec = nvec + 1
    endif
    ncolsgrid = 1 + ndimV*nvec
 endif
 !
 !--allocate memory for the grid
 ! 
 if (allocated(datgrid)) deallocate(datgrid)

 write(*,"(a,i5,2(' x',i5),a)",advance='no') ' >>> allocating memory for ',npixels(:),' grid ...'
 allocate(datgrid(npixels(1),npixels(2),npixels(3)),stat=ierr)
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
 call open_gridfile(iunit,filename,outformat,npixels,ncolsgrid,time,ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR: could not open grid file for output, skipping...'
    if (allocated(datgrid)) deallocate(datgrid)
    if (allocated(weight)) deallocate(weight)
    return
 endif

 fmtstring1 = "(12x,a20,1x,'    min    ',1x,'    max    ',1x,'    mean    ')"
 fmtstring = "(22x,a10,1x,3(es10.2,1x))"
 !
 !--interpolate density to the 3D grid
 !
 print "(/,a)",' interpolating density to 3D grid...'

 call minmaxmean_part(dat(1:ninterp,irho:irho),weight,ninterp,partmin,partmax,partmean)
 print fmtstring1,trim(label(irho))
 print fmtstring,' on parts:',partmin(1),partmax(1),partmean(1)

 call interpolate3D(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)),&
      dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,irho),icolourme,ninterp,&
      xmin(1),xmin(2),xmin(3),datgrid,npixels(1),npixels(2),npixels(3),&
      pixwidth,pixwidth,inormalise,isperiodic)
 !
 !--set minimum density on the grid
 !
 call minmaxmean_grid(datgrid,npixels,gridmin,gridmax,gridmean,nonzero=.true.)
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
    !$omp parallel do private(k) schedule(static)
    do k=1,npixels(3)
       where (datgrid(:,:,k).le.tiny(datgrid))
          datgrid(:,:,k) = rhomin
       end where
    enddo
 else
    print*,'minimum density NOT enforced'
 endif

 !
 !--write density to grid data file
 !
 print*
 call write_grid(iunit,filename,outformat,datgrid,npixels,trim(label(irho)),time,ierr)

 !
 !--interpolate remaining quantities to the 3D grid
 !
 if (interpolateall) then
    do i=1,ndataplots
       if (.not.any(ix(:).eq.i) .and. i.ne.ih .and. i.ne.ipmass .and. i.ne.irho) then
 
          print "(/,a)",' interpolating '//trim(label(i))
          print fmtstring1,trim(label(i))
          call minmaxmean_part(dat(1:ninterp,i:i),weight,ninterp,partmin,partmax,partmean)
          print fmtstring,' on parts:',partmin(1),partmax(1),partmean(1)

          call interpolate3D(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)),&
               dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,i),icolourme,ninterp,&
               xmin(1),xmin(2),xmin(3),datgrid,npixels(1),npixels(2),npixels(3),&
               pixwidth,pixwidth,.true.,isperiodic)

          call minmaxmean_grid(datgrid,npixels,gridmin,gridmax,gridmean,.false.)
          print fmtstring,' on grid :',gridmin,gridmax,gridmean
          !
          !--write gridded data to file
          !
          call write_grid(iunit,filename,outformat,datgrid,npixels,trim(label(i)),time,ierr)

       endif
    enddo

 else
    if (allocated(datgrid)) deallocate(datgrid)
        
    if (nvec.gt.0) then

       write(*,"(/,a,i5,2(' x',i5),a)",advance='no') ' >>> allocating memory for ',npixels(:),' x 3 grid ...'
       allocate(datgridvec(3,npixels(1),npixels(2),npixels(3)),stat=ierr)
       if (ierr /= 0) then
          write(*,*) 'FAILED: NOT ENOUGH MEMORY'
          return
       else
          write(*,*) 'OK'
       endif

       !
       !--interpolate velocity field and magnetic fields and write to file
       !
       over_vec: do ivec=1,nvec
          select case(ivec)
          case(1)
             iloc = ivx
             print "(a)",' interpolating velocity field to 3D grid...'
          case(2)
             iloc = iBfirst
             print "(a)",' interpolating magnetic field to 3D grid...'
          case default
             iloc = 0
             exit over_vec
          end select
          
          print fmtstring1,trim(labelvec(iloc))
          call minmaxmean_part(dat(1:ninterp,iloc:iloc+ndimV-1),weight,ninterp,partmin,partmax,partmean)
          do i=1,ndimV
             print fmtstring,' on parts:',partmin(i),partmax(i),partmean(i)
          enddo

          call interpolate3D_vec(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)),&
               dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,iloc:iloc+ndimV-1),icolourme,ninterp,&
               xmin(1),xmin(2),xmin(3),datgridvec,npixels(1),npixels(2),npixels(3),&
               pixwidth,pixwidth,.true.,isperiodic)

          call minmaxmean_gridvec(datgridvec,npixels,ndimV,datmin,datmax,datmean)
          do i=1,ndimV
             print fmtstring,' on grid :',datmin(i),datmax(i),datmean(i)
          enddo
          !
          !--write result to grid file
          !
          do i=1,ndimV
             call write_grid(iunit,filename,outformat,datgridvec(i,:,:,:),npixels,&
                             trim(label(iloc+i-1)),time,ierr)
          enddo
          print*
       enddo over_vec
    endif

 endif

 close(iunit)

 if (allocated(datgrid)) deallocate(datgrid)
 if (allocated(datgridvec)) deallocate(datgridvec)
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

 !$omp parallel do reduction(min:partmin) &
 !$omp reduction(max:partmax) reduction(+:partmean,np) &
 !$omp private(j,partval)
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
 if (np.gt.0) then
    partmean(:) = partmean(:)/real(np)
 endif

 return
end subroutine minmaxmean_part

end module convert_grid
