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
!  Copyright (C) 2005-2023 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!  Module containing routines for converting 3D SPH dump
!  files to 3D gridded data.
!-----------------------------------------------------------------
module convert_grid
 use params, only:doub_prec
 implicit none
 private
 public :: convert_to_grid

contains

!-----------------------------------------------------------------
! interpolate 3D SPH data to grid and interface to grid
! data output routines
!-----------------------------------------------------------------
subroutine convert_to_grid(time,dat,ntypes,npartoftype,masstype,itype,ncolumns,filename,&
                           outformat,interpolateall,icols,rhogrid,dat3D)
 use labels,               only:label,labelvec,irho,ih,ipmass,ix,ivx,iBfirst,get_sink_type
 use limits,               only:lim,get_particle_subset
 use settings_units,       only:units,unit_interp
 use settings_data,        only:ndim,ndimV,UseTypeInRenderings,iRescale,required,debugmode,icoordsnew,xorigin,iverbose
 use settings_part,        only:iplotpartoftype
 use settings_render,      only:npix,inormalise=>inormalise_interpolations,&
                                idensityweightedinterpolation,exact_rendering
 use settings_xsecrot,     only:anglex,angley,anglez
 use rotation,             only:rotate3D
 use params,               only:int1
 use interpolation,        only:set_interpolation_weights,get_n_interp
 use interpolations3D,     only:interpolate3D,interpolate3D_vec
 use interpolations3Dgeom, only:interpolate3Dgeom,interpolate3Dgeom_vec
 use interpolations2D,     only:interpolate2D,interpolate2D_vec
 use system_utils,         only:renvironment,ienvlist
 use readwrite_griddata,   only:open_gridfile_w,write_grid,write_gridlimits
 use particle_data,        only:icolourme
 use params,               only:int8
 use geometry,             only:coord_is_length,igeom_cartesian,labelcoord,labelcoordsys
 use asciiutils,           only:strip
 use timing,               only:wall_time,print_time
 use filenames,            only:tagline
 integer, intent(in)                          :: ntypes,ncolumns
 integer, intent(in), dimension(:)            :: npartoftype
 integer(kind=int1), intent(in), dimension(:) :: itype
 real, intent(in)                             :: time
 real, intent(inout), dimension(:,:)             :: dat
 real, intent(in), dimension(:)               :: masstype
 character(len=*), intent(in)                 :: filename,outformat
 logical, intent(in)                          :: interpolateall
 integer, intent(in),  optional               :: icols(:)
 real,    intent(out), allocatable, optional  :: rhogrid(:,:,:),dat3D(:,:,:)
 integer, parameter :: iunit = 89
 integer            :: ierr,i,k,ncolsgrid,ivec,nvec,iloc,j,nzero
 integer            :: npixx,ninterp,isinktype
 character(len=40)  :: fmtstring
 character(len=64)  :: fmtstring1

 real(doub_prec), dimension(:,:,:), allocatable   :: datgrid,rhogrid_tmp
 real, dimension(:,:),   allocatable   :: datgrid2D
 real(doub_prec), dimension(:,:,:,:), allocatable :: datgridvec
 real, dimension(:,:,:),   allocatable :: datgridvec2D
 real, dimension(:), allocatable       :: weight,x,y,z
 logical, dimension(:), allocatable    :: mask
 real, dimension(3)    :: xmin,xmax
 real, dimension(3)    :: partmin,partmax,partmean,part_integral
 real, dimension(3)    :: datmin,datmax,datmean
 real(doub_prec)       :: dtime
 real(doub_prec), dimension(3) :: xmind,xmaxd
 integer, dimension(3) :: npixels
 integer(kind=int8), dimension(3) :: npixels8
 integer, dimension(12) :: icoltogrid
 integer :: ncolstogrid,igeom
 real    :: hmin,pixwidth,pixwidthx(3),rhominset,rhomin
 real    :: gridmin,gridmax,gridmean,grid_integral
 real    :: mtot,mtotgrid,err,t2,t1,xi(3),ax,ay,az
 real, parameter :: pi=4.0*atan(1.0)
 logical :: lowmem,do_output
 logical, dimension(3) :: isperiodic
 character(len=len(labelcoord)), dimension(3) :: xlab
 character(len=120) :: origin

 dtime = real(time,kind=doub_prec)
 do_output = .not.(trim(filename)=='none' .or. trim(outformat)=='none')
 !
 !--check for errors in input settings
 !
 if (ndim < 2 .or. ndim > 3) then
    print "(/,a,i2,a,/)",' ERROR: SPH data has ',ndim,' spatial dimensions: cannot convert to 3D grid'
    return
 endif

 print "(/,'----->',1x,a,i1,a,/)",'CONVERTING SPH DATA -> ',ndim,'D GRID'
 xmin(1:ndim) = lim(ix(1:ndim),1)
 xmax(1:ndim) = lim(ix(1:ndim),2)
 xlab(:) = (/'x','y','z'/)
 igeom = max(icoordsnew,1)  ! ensure it is not zero
 if (igeom /= igeom_cartesian) xlab = strip(labelcoord(:,igeom),'\')
 !
 !--print limits information
 !
 call write_gridlimits(ndim,xmin,xmax,label(ix(1:ndim)))
 !
 !--get environment variable options
 !
 if (present(icols)) then
    icoltogrid(:) = 0
    icoltogrid(1:size(icols)) = icols(:)
    ncolstogrid = count(icols > 0)
 else
    call get_splash2grid_options(ndim,ncolstogrid,icoltogrid,isperiodic,xlab)
 endif
 !
 !--check for errors
 !
 ierr = 0
 do i=1,ndim
    if ((xmax(i)-xmin(i)) < tiny(0.)) then
       print "(a)",' ERROR: min = max in '//trim(label(ix(i)))//&
                   ' coordinate: cannot interpolate to zero-sized grid!'
       ierr = 1
    endif
 enddo

 if (irho <= 0 .or. irho > ncolumns) then
    print "(a)",' ERROR: density not found in data read.'
    ierr = 2
 endif
 if (ih <= 0 .or. ih > ncolumns) then
    print "(a)",' ERROR: smoothing length not found in data read.'
    ierr = 3
 endif
 if (ipmass <= 0 .or. ipmass > ncolumns) then
    if (all(masstype(:) < tiny(0.))) then
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
 !  and allocate memory for weights
 !
 ninterp = get_n_interp(ntypes,npartoftype,UseTypeInRenderings,iplotpartoftype,size(itype),.false.)
 allocate(weight(ninterp),stat=ierr)
 allocate(x(ninterp),y(ninterp),z(ninterp),mask(ninterp),stat=ierr)
 if (ierr /= 0) then
    print*,' ERROR allocating memory for interpolation weights, aborting...'
    return
 endif
 !
 !--set interpolation weights (w = m/(rho*h^ndim)
 !
 isinktype = get_sink_type(ntypes)
 call set_interpolation_weights(weight,dat,itype,(iplotpartoftype .and. UseTypeInRenderings),&
      ninterp,npartoftype,masstype,ntypes,ncolumns,irho,ipmass,ih,ndim,iRescale,&
      idensityweightedinterpolation,inormalise,units,unit_interp,required,.false.,isinktype)
 !
 !--set default mask and apply range restrictions to data
 !
 icolourme(:) = 1
 call get_particle_subset(icolourme,dat,ncolumns)
 mask = (icolourme > 0 .and. weight > 0.)

 !
 !--work out how many pixels to use
 !
 npixx = npix
 npixels = ienvlist('SPLASH_TO_GRID_NPIX',3)
 if (product(npixels) > 0) then
    print "(a,2(i5,','),i5)",' Using --npix=',npixels
 else
    if (npixx <= 0) then
       print "(/,a)",' WARNING: number of pixels = 0, using automatic pixel numbers'
       hmin = 0.
       call minmaxmean_part(dat(:,ih:ih),mask,ninterp,partmin,partmax,partmean,nonzero=.true.)
       hmin = partmin(1)
       if (hmin > 0. .and. igeom==igeom_cartesian) then
          print*,'based on the minimum smoothing length of hmin = ',hmin
          npixels8(1:ndim) = int((xmax(1:ndim) - xmin(1:ndim))/hmin,kind=int8) + 1
          if (ndim==3) then
             print "(a,i6,2(' x',i6),a)",' requires ',npixels8(1:ndim),' pixels to capture the full resolution'
             if (product(npixels8(1:ndim)) > 512**3 .or. product(npixels8(1:ndim)) <= 0) then
                npixx = 512
                print "(a,i4)",' but this is ridiculous, so instead we choose ',npixx
             else
                npixx = int(npixels8(1))
             endif
          else
             print "(a,i6,1(' x',i6),a)",' requires ',npixels8(1:ndim),' pixels to capture the full resolution'
             if (product(npixels8(1:ndim)) > 1024**ndim .or. product(npixels8(1:ndim)) <= 0) then
                npixx = 1024
                print "(a,i4)",' but this is very large, so instead we choose ',npixx
             else
                npixx = int(npixels8(1))
             endif
          endif
       else
          npixx = 512
          if (hmin <= 0.) then
             print "(a)",' ...but cannot get auto pixel number because hmin = 0'
          else
             print "(a)",' ...but cannot get auto pixel number because of non-cartesian geometry'
          endif
          print "(a)",'    so instead we choose npixels = ',npixx
       endif
    endif
    print "(a,/)",' Set this manually using --npix=100,100,100'
    pixwidth = (xmax(1)-xmin(1))/npixx
    do i=1,ndim
       if (coord_is_length(i,igeom)) then
          npixels(i) = int((xmax(i)-xmin(i) - 0.5*pixwidth)/pixwidth) + 1
       else
          npixels(i) = npixx
       endif
    enddo
 endif
 pixwidthx(:) = (xmax(:) - xmin(:))/npixels(:)
 !
 !--work out how many columns will be written to file
 !
 nvec = 0
 if (ncolstogrid > 0) then
    ncolsgrid = ncolstogrid
 elseif (interpolateall) then
    ncolsgrid = 0
    do i=1,ncolumns
       if (.not.any(ix(1:ndim)==i) .and. i /= ih .and. i /= ipmass) then
          ncolsgrid = ncolsgrid + 1
       endif
    enddo
 else
    if (ndimV==ndim) then
       if (ivx > 0 .and. ivx+ndimV-1 <= ncolumns) nvec = nvec + 1
       if (iBfirst > 0 .and. iBfirst+ndimV-1 <= ncolumns) nvec = nvec + 1
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

 if (lowmem .and. nvec > 0) &
    print "(a,/)",' [doing velocity field components separately (low memory mode)]'
 !
 !--allocate memory for the grid
 !
 if (allocated(datgrid))   deallocate(datgrid)
 if (allocated(datgrid2D)) deallocate(datgrid2D)

 if (ndim==3) then
    write(*,"(a,i5,2(' x',i5),a)",advance='no') ' >>> allocating memory for ',npixels(1:ndim),' grid ...'
    allocate(datgrid(npixels(1),npixels(2),npixels(3)),stat=ierr)
 elseif (ndim==2) then
    write(*,"(a,i5,1(' x',i5),a)",advance='no') ' >>> allocating memory for ',npixels(1:ndim),' grid ...'
    allocate(datgrid2D(npixels(1),npixels(2)),stat=ierr)
 endif
 if (ierr /= 0) then
    write(*,*) 'FAILED: NOT ENOUGH MEMORY'
    call deallocate_memory()
    return
 else
    write(*,*) 'OK'
 endif

 !
 !--open grid file for output (also checks format is OK)
 !
 if (do_output) then
    origin='"splash to '//trim(outformat)//'" on file '//trim(filename)
    xmind = xmin; xmaxd = xmax  ! convert to double precision
    call open_gridfile_w(iunit,filename,outformat,ndim,ncolsgrid,npixels(1:ndim),&
                         xmind(1:ndim),xmaxd(1:ndim),dtime,ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR: could not open grid file for output, skipping...'
       call deallocate_memory()
       return
    endif
 endif

 fmtstring1 = "(12x,a20,1x,'    min    ',1x,'    max    ',1x,'    mean    ')"
 fmtstring = "(22x,a10,1x,3(es10.2,1x))"
 !
 !--interpolate density to the 3D grid
 !
 print "(/,a,i1,a)",' interpolating density to ',ndim,'D grid...'
 if (debugmode) print*,'DEBUG: density in column ',irho,' vals = ',dat(1:10,irho)

 call minmaxmean_part(dat(1:ninterp,irho:irho),mask,ninterp,partmin,partmax,partmean)
 print fmtstring1,trim(label(irho))
 print fmtstring,' on parts:',partmin(1),partmax(1),partmean(1)

 if (ipmass > 0.) then
    mtot = sum(dat(1:ninterp,ipmass),mask=(icolourme(1:ninterp) > 0 .and. weight(1:ninterp) > 0.))
    print "(9x,a23,1x,es10.4)",'total mass on parts:',mtot
 endif

 if (ndim==3) then
    call wall_time(t1)
    if (igeom /= igeom_cartesian) then
       call interpolate3Dgeom(igeom,dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)),&
            dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,irho),icolourme,ninterp,&
            xmin,datgrid,npixels,pixwidthx,xorigin,inormalise,isperiodic)
    else
       x = dat(1:ninterp,ix(1))
       y = dat(1:ninterp,ix(2))
       z = dat(1:ninterp,ix(3))
       if (abs(anglez)>0. .or. abs(angley)>0. .or. abs(anglex)>0.) then
          print*, 'Rotating particles around (z,y,x) by',anglez,angley,anglex
          print*, 'WARNING: This does not rotate vector components'
          ax = anglex*pi/180.0 ! convert degrees to radians to pass into rotate
          ay = angley*pi/180.0
          az = anglez*pi/180.0
          do i=1,ninterp
             xi = (/x(i),y(i),z(i)/)
             call rotate3D(xi,ax,ay,az,0.,0.)
             x(i) = xi(1)
             y(i) = xi(2)
             z(i) = xi(3)
          enddo
       endif
       mask = (mask .and. (x >= xmin(1) .and. x <= xmax(1)) .and. &
                          (y >= xmin(2) .and. y <= xmax(2)) .and. &
                          (z >= xmin(3) .and. z <= xmax(3)))
       if (ipmass > 0.) then
          mtot = sum(dat(1:ninterp,ipmass),mask=mask)
          print "(9x,a23,1x,es10.4,/)",'mass on parts in box:',mtot
       endif

       call interpolate3D(x,y,z,&
            dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,irho),icolourme,ninterp,&
            xmin(1),xmin(2),xmin(3),datgrid,npixels(1),npixels(2),npixels(3),&
            pixwidthx(1),pixwidthx(2),pixwidthx(3),inormalise,&
            isperiodic(1),isperiodic(2),isperiodic(3))
       mtotgrid = sum(datgrid)*product(pixwidthx)
    endif
    if (present(rhogrid)) rhogrid = datgrid
    rhogrid_tmp = datgrid
    call wall_time(t2)
    if (t2 - t1 > 1.) call print_time(t2-t1)
    !
    !--set minimum density on the grid
    !
    call minmaxmean_grid(datgrid,npixels,gridmin,gridmax,gridmean,nonzero=.true.)
 else
    call interpolate2D(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),&
         dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,irho),icolourme,ninterp,&
         xmin(1),xmin(2),datgrid2D,npixels(1),npixels(2),&
         pixwidth,pixwidth,inormalise,exact_rendering,isperiodic(1),isperiodic(2),iverbose)
    !
    !--set minimum density on the grid
    !
    call minmaxmean_grid2D(datgrid2D,npixels,gridmin,gridmax,gridmean,nonzero=.true.)
    mtot = sum(datgrid2D)*pixwidth**2*npixels(1)*npixels(2)
 endif

 print fmtstring1,trim(label(irho))
 print fmtstring,' on grid :',gridmin,gridmax,gridmean
 !
 !--print error if mass is not conserved
 !
 if (mtotgrid > 0.) then
    print "(9x,a23,1x,es10.4,/)",'total mass on grid:',mtotgrid
    err = 100.*(mtotgrid - mtot)/mtot
    if (abs(err) > 1) print "(/,a,1pg8.1,a,/)",' WARNING! MASS NOT CONSERVED BY ',err,&
    '% BY INTERPOLATION'
 endif

 rhomin = gridmin
 rhominset = renvironment('SPLASH_TO_GRID_RHOMIN',errval=-1.)

 print*
 if (rhominset >= 0.) then
    rhomin = rhominset
    print*,'enforcing minimum density on grid = ',rhomin
    print*,'(based on --rhomin setting)'
 elseif (rhomin > 0.) then
    print*,'enforcing minimum density on grid = ',rhomin
    print "(a)",' ** set --rhomin=minval to manually set this (e.g. to zero) **'
 endif

 if (rhomin > 0.) then
    nzero = 0
    if (ndim==3) then
       !$omp parallel do private(k,j,i) reduction(+:nzero) schedule(static)
       do k=1,npixels(3)
          do j=1,npixels(2)
             do i=1,npixels(1)
                if (datgrid(i,j,k) <= tiny(datgrid)) then
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
             if (datgrid2D(i,j) <= tiny(datgrid2D)) then
                datgrid2D(i,j) = rhomin
                nzero = nzero + 1
             endif
          enddo
       enddo
    endif
    print "(a,i0,a)",' minimum density enforced on ',nzero,' grid cells'
 else
    print*,'minimum density NOT enforced'
 endif

 !
 !--write density to grid data file
 !
 print*
 if (ndim==3 .and. do_output) then
    if (lowmem .or. interpolateall .or. ncolstogrid > 0) then
       call write_grid(iunit,filename,outformat,ndim,1,npixels,trim(label(irho)),&
                       labelcoordsys(igeom),xlab,dtime,pixwidthx,xmin,xmax,ierr,dat=datgrid,&
                       tagline=tagline,origin=origin)
    endif
 elseif (do_output) then
    call write_grid(iunit,filename,outformat,ndim,1,npixels,trim(label(irho)),&
                    labelcoordsys(igeom),xlab,dtime,pixwidthx,xmin,xmax,ierr,dat2D=datgrid2D,&
                    tagline=tagline,origin=origin)
 endif
 !
 !--use density weighted interpolation for remaining quantities
 !
 idensityweightedinterpolation = .true.
 call set_interpolation_weights(weight,dat,itype,(iplotpartoftype .and. UseTypeInRenderings),&
     ninterp,npartoftype,masstype,ntypes,ncolumns,irho,ipmass,ih,ndim,iRescale,&
     idensityweightedinterpolation,inormalise,units,unit_interp,required,.false.,isinktype)

 !
 !--interpolate remaining quantities to the 3D grid
 !
 if (interpolateall .or. ncolstogrid > 0) then
    if (ncolstogrid > 0 .and. .not.present(icols)) then
       print "(/,a,i2,a)",' Interpolating ',ncolstogrid,' columns to grid from --grid setting:'
       print "(' got --grid=',10(i2,1x))",icoltogrid(1:ncolstogrid)
    endif

    do i=1,ncolumns
       if ((ncolstogrid > 0 .and. any(icoltogrid==i) .and. i /= irho) .or.  &
           (interpolateall .and. &
            .not.any(ix(:)==i) .and. i /= ih .and. i /= ipmass .and. i /= irho)) then

          print "(/,a)",' interpolating '//trim(label(i))
          print fmtstring1,trim(label(i))
          call minmaxmean_part(dat(1:ninterp,i:i),mask,ninterp,partmin,partmax,partmean)
          print fmtstring,' on parts:',partmin(1),partmax(1),partmean(1)

          if (iszero(partmin,partmax,1)) then
             datgrid = 0.
          else
             if (ndim==3) then
                if (igeom /= igeom_cartesian) then
                   call interpolate3Dgeom(igeom,dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)),&
                        dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,i),icolourme,ninterp,&
                        xmin,datgrid,npixels,pixwidthx,xorigin,.true.,isperiodic)
                else
                   call interpolate3D(x,y,z,&
                        dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,i),icolourme,ninterp,&
                        xmin(1),xmin(2),xmin(3),datgrid,npixels(1),npixels(2),npixels(3),&
                        pixwidthx(1),pixwidthx(2),pixwidthx(3),.true.,isperiodic(1),isperiodic(2),isperiodic(3))
                endif
                if (present(dat3D)) dat3D = datgrid
             else
                call interpolate2D(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),&
                     dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,i),icolourme,ninterp,&
                     xmin(1),xmin(2),datgrid2D,npixels(1),npixels(2),&
                     pixwidth,pixwidth,.true.,exact_rendering,isperiodic(1),isperiodic(2),iverbose)
             endif
          endif

          !
          !--write gridded data to file
          !
          if (ndim==3) then
             call minmaxmean_grid(datgrid,npixels,gridmin,gridmax,gridmean,.false.)
             print fmtstring,' on grid :',gridmin,gridmax,gridmean
             if (do_output) call write_grid(iunit,filename,outformat,ndim,1,npixels,trim(label(i)),&
                  labelcoordsys(igeom),xlab,dtime,pixwidthx,xmin,xmax,ierr,dat=datgrid)
          else
             call minmaxmean_grid2D(datgrid2D,npixels,gridmin,gridmax,gridmean,.false.)
             print fmtstring,' on grid :',gridmin,gridmax,gridmean
             if (do_output) call write_grid(iunit,filename,outformat,ndim,1,npixels,trim(label(i)),&
                  labelcoordsys(igeom),xlab,dtime,pixwidthx,xmin,xmax,ierr,dat2D=datgrid2D)
          endif
       endif
    enddo

 else

    if (nvec > 0) then

       print "(/,a,i2,a)",' set --grid=',irho,' to interpolate density ONLY and skip remaining columns'
       print "(a,i2,a)",  '     --grid=6,8,10 to select particular columns'

       if (.not.lowmem) then
          if (ndim==3) then
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
          if (iloc <= 0 .or. iloc >= ncolumns) cycle over_vec

          if (lowmem) then

             do i=iloc,iloc+ndimV-1
                print "(/,a)",' interpolating '//trim(label(i))
                print fmtstring1,trim(label(i))
                if (ipmass > 0) then
                   call minmaxmean_part(dat(1:ninterp,i:i),mask,ninterp,partmin,partmax,&
                        partmean,mass=dat(1:ninterp,ipmass),part_integral=part_integral)
                   print fmtstring,' on parts:',partmin(1),partmax(1),partmean(1)
                   print "(9x,a23,1x,es10.2,/)",' sum(m*'//trim(label(i))//'):',part_integral(1)
                else
                   call minmaxmean_part(dat(1:ninterp,i:i),mask,ninterp,partmin,partmax,partmean)
                   print fmtstring,' on parts:',partmin(1),partmax(1),partmean(1)
                endif

                if (iszero(partmin,partmax,1)) then
                   datgrid = 0.
                else
                   if (ndim==3) then
                      if (igeom /= igeom_cartesian) then
                         call interpolate3Dgeom(igeom,dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)),&
                           dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,i),icolourme,ninterp,&
                           xmin,datgrid,npixels,pixwidthx,xorigin,.true.,isperiodic)
                      else
                         call interpolate3D(x,y,z,&
                              dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,i),icolourme,ninterp,&
                              xmin(1),xmin(2),xmin(3),datgrid,npixels(1),npixels(2),npixels(3),&
                              pixwidthx(1),pixwidthx(2),pixwidthx(3),.true.,isperiodic(1),isperiodic(2),isperiodic(3))
                      endif
                   else
                      call interpolate2D(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),&
                           dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,i),icolourme,ninterp,&
                           xmin(1),xmin(2),datgrid2D,npixels(1),npixels(2),&
                           pixwidth,pixwidth,.true.,exact_rendering,isperiodic(1),isperiodic(2),iverbose)
                   endif
                endif

                !
                !--write gridded data to file
                !
                if (ndim==3) then
                   call minmaxmean_grid(datgrid,npixels,gridmin,gridmax,gridmean,.false.,&
                                        rho=rhogrid_tmp,grid_integral=grid_integral)
                   print fmtstring,' on grid :',gridmin,gridmax,gridmean
                   grid_integral = grid_integral*product(pixwidthx)
                   print "(9x,a23,1x,es10.2,/)",' int(rho*'//trim(label(i))//')dV:',grid_integral

                   call write_grid(iunit,filename,outformat,ndim,1,npixels,trim(label(i)),&
                        labelcoordsys(igeom),xlab,dtime,pixwidthx,xmin,xmax,ierr,&
                        dat=datgrid,tagline=tagline,origin=origin)
                else
                   call minmaxmean_grid2D(datgrid2D,npixels,gridmin,gridmax,gridmean,.false.)
                   print fmtstring,' on grid :',gridmin,gridmax,gridmean
                   call write_grid(iunit,filename,outformat,ndim,1,npixels,trim(label(i)),&
                        labelcoordsys(igeom),xlab,dtime,pixwidthx,xmin,xmax,ierr,&
                        dat2D=datgrid2D,tagline=tagline,origin=origin)
                endif
             enddo
          else

             print fmtstring1,trim(labelvec(iloc))
             call minmaxmean_part(dat(1:ninterp,iloc:iloc+ndimV-1),mask,ninterp,partmin,partmax,partmean)
             do i=1,ndimV
                print fmtstring,' on parts:',partmin(i),partmax(i),partmean(i)
             enddo

             if (iszero(partmin,partmax,ndimV)) then
                datgridvec = 0.
             else
                if (ndim==3) then
                   if (igeom /= igeom_cartesian) then
                      call interpolate3Dgeom_vec(igeom,dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),dat(1:ninterp,ix(3)),&
                           dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,iloc:iloc+ndimV-1),icolourme,ninterp,&
                           xmin,datgridvec,npixels,pixwidthx,xorigin,.true.,isperiodic)
                   else
                      call interpolate3D_vec(x,y,z,&
                        dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,iloc:iloc+ndimV-1),icolourme,ninterp,&
                        xmin(1),xmin(2),xmin(3),datgridvec,npixels(1),npixels(2),npixels(3),&
                        pixwidthx(1),pixwidthx(2),pixwidthx(3),.true.,isperiodic(1),isperiodic(2),isperiodic(3))
                   endif
                else
                   call interpolate2D_vec(dat(1:ninterp,ix(1)),dat(1:ninterp,ix(2)),&
                        dat(1:ninterp,ih),weight(1:ninterp),dat(1:ninterp,iloc),dat(1:ninterp,iloc+1), &
                        icolourme,ninterp,xmin(1),xmin(2),datgridvec2D(1,:,:),datgridvec2D(2,:,:), &
                        npixels(1),npixels(2),pixwidth,pixwidth,.true.,exact_rendering,isperiodic(1),isperiodic(2))
                endif
             endif

             if (ndim==3) then
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
             if (ndim==3 .and. do_output) then
                call write_grid(iunit,filename,outformat,ndim,ndimV,npixels,&
                                label(irho),labelcoordsys(igeom),xlab,dtime,pixwidthx,xmin,xmax,ierr,&
                                dat=datgrid,dat3D=datgridvec,label3D=label(iloc:iloc+ndimV),tagline=tagline,origin=origin)
             elseif (do_output) then
                do i=1,ndimV
                   call write_grid(iunit,filename,outformat,ndim,ndimV,npixels,&
                                   label(iloc+i-1),labelcoordsys(igeom),xlab,dtime,pixwidthx,&
                                   xmin,xmax,ierr,dat2D=datgridvec2D(i,:,:),tagline=tagline,origin=origin)
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

 call deallocate_memory()
 return

contains
!------------------------------------------------------------
! manual mop-up of allocatable memory
!------------------------------------------------------------
 subroutine deallocate_memory()

  if (allocated(datgrid))      deallocate(datgrid)
  if (allocated(datgrid2D))    deallocate(datgrid2D)
  if (allocated(datgridvec))   deallocate(datgridvec)
  if (allocated(datgridvec2D)) deallocate(datgridvec2D)
  if (allocated(weight)) deallocate(weight)
  if (allocated(x)) deallocate(x)
  if (allocated(y)) deallocate(y)
  if (allocated(z)) deallocate(z)

 end subroutine deallocate_memory

end subroutine convert_to_grid

!------------------------------------------------------------
! get options for splash to grid from environment variables
!------------------------------------------------------------
subroutine get_splash2grid_options(ndim,ncolstogrid,icoltogrid,isperiodic,xlab)
 use system_utils, only:lenvironment,renvironment,envlist,lenvstring,ienvstring
 use labels,       only:irho
 integer, intent(in)  :: ndim
 integer, intent(out) :: ncolstogrid,icoltogrid(:)
 logical, intent(out) :: isperiodic(3)
 character(len=1), intent(in) :: xlab(3)
 integer :: nstring,i,icol
 character(len=30), dimension(12) :: strings
 !
 !--SPLASH_TO_GRID can be set to comma separated list of columns
 !  in order to select particular quantities for interpolation to grid
 !
 ncolstogrid   = 0
 icoltogrid(:) = 0
 call envlist('SPLASH_TO_GRID',nstring,strings)
 if (nstring > 0) then
    do i=1,nstring
       icol = ienvstring(strings(i))
       if (ienvstring(strings(i)) > 0) then
          ncolstogrid = ncolstogrid + 1
          icoltogrid(ncolstogrid) = icol
       endif
    enddo
 endif
 !
 !--for backwards compatibility, support the SPLASH_TO_GRID_DENSITY_ONLY option
 !  but only if SPLASH_TO_GRID is not set
 !
 if (ncolstogrid==0 .and. lenvironment('SPLASH_TO_GRID_DENSITY_ONLY')) then
    ncolstogrid = 1
    icoltogrid(1) = irho
 endif
 !
 !--whether or not to wrap particle contributions across boundaries
 !
 isperiodic(:) = .false.
 call envlist('SPLASH_TO_GRID_PERIODIC',nstring,strings)
 if (nstring > ndim) then
    print "(a)",' ERROR in --periodic setting'
    nstring = ndim
 endif
 do i=1,nstring
    isperiodic(i) = lenvstring(strings(i))
 enddo
 if (nstring==1) isperiodic(2:ndim) = isperiodic(1)

 if (all(isperiodic(1:ndim))) then
    print "(/,a)",' using PERIODIC boundaries (from --periodic setting)'
 elseif (isperiodic(1) .or. isperiodic(2) .or. isperiodic(3)) then
    print*
    do i=1,ndim
       if (isperiodic(i)) then
          print "(a)",' using PERIODIC boundaries in '//xlab(i)//' (from --periodic flag)'
       else
          print "(a)",' using NON-PERIODIC bounds in '//xlab(i)//' (from --periodic flag)'
       endif
    enddo
 else
    print "(/,a)",' using NON-PERIODIC boundaries: use --periodic=yes for periodic'
    if (ndim==3) then
       print "(a,/)",'                                 or --periodic=yes,no,yes for mixed'
    else
       print "(a,/)",'                                 or --periodic=yes,no for mixed'
    endif
 endif

end subroutine

!-----------------------------------------------
! calculate max and min and mean values on grid
!-----------------------------------------------
subroutine minmaxmean_grid(datgrid,npixels,gridmin,gridmax,gridmean,nonzero,rho,grid_integral)
 real(doub_prec), dimension(:,:,:), intent(in) :: datgrid
 integer, dimension(3), intent(in)  :: npixels
 real, intent(out)                  :: gridmin,gridmax,gridmean
 logical, intent(in)                :: nonzero
 real(doub_prec), dimension(:,:,:), intent(in), optional :: rho
 real, intent(out), optional        :: grid_integral

 real    :: dati,gridint
 integer :: i,j,k

 gridmax = -huge(gridmax)
 gridmin = huge(gridmin)
 gridmean = 0.
 gridint = 0.
 !$omp parallel do schedule(static) default(shared) &
 !$omp reduction(min:gridmin) &
 !$omp reduction(max:gridmax) reduction(+:gridmean,gridint) &
 !$omp private(k,j,i,dati)
 do k=1,npixels(3)
    do j=1,npixels(2)
       do i=1,npixels(1)
          dati = datgrid(i,j,k)
          gridmax  = max(gridmax,dati)
          if (nonzero) then
             if (dati > tiny(0.)) gridmin  = min(gridmin,dati)
          else
             gridmin  = min(gridmin,dati)
          endif
          gridmean = gridmean + dati
          if (present(grid_integral) .and. present(rho)) then
             gridint = gridint + dati*rho(i,j,k)
          endif
       enddo
    enddo
 enddo
 if (gridmin >= huge(gridmin)) gridmin = 0.
 gridmean = gridmean/product(npixels(1:3))
 if (present(grid_integral)) grid_integral = gridint

end subroutine minmaxmean_grid

!----------------------------------------------------
! calculate max and min and mean values on grid (2D)
!----------------------------------------------------
subroutine minmaxmean_grid2D(datgrid,npixels,gridmin,gridmax,gridmean,nonzero)
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
          if (dati > tiny(0.)) gridmin  = min(gridmin,dati)
       else
          gridmin  = min(gridmin,dati)
       endif
       gridmean = gridmean + dati
    enddo
 enddo
 if (gridmin >= huge(gridmin)) gridmin = 0.
 gridmean = gridmean/product(npixels(1:2))

end subroutine minmaxmean_grid2D

!-----------------------------------------------
! calculate max and min and mean values on grid
! (for vector quantities)
!-----------------------------------------------
subroutine minmaxmean_gridvec(datgridvec,npixels,jlen,gridmin,gridmax,gridmean)
 real(doub_prec), dimension(:,:,:,:), intent(in) :: datgridvec
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
 !$omp end parallel do
 where (gridmin >= huge(gridmin)) gridmin = 0.
 gridmean(1:jlen) = gridmean(1:jlen)/real(product(npixels(1:3)))

end subroutine minmaxmean_gridvec

!-----------------------------------------------
! calculate max and min and mean values on grid
! (for vector quantities)
!-----------------------------------------------
subroutine minmaxmean_gridvec2D(datgridvec,npixels,jlen,gridmin,gridmax,gridmean)
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
 !$omp end parallel do
 gridmean(1:jlen) = gridmean(1:jlen)/real(product(npixels(1:2)))

end subroutine minmaxmean_gridvec2D

!----------------------------------------------------
! calculate max and min and mean values on particles
!----------------------------------------------------
subroutine minmaxmean_part(dat,mask,npart,partmin,partmax,partmean,nonzero,mass,part_integral)
 real, dimension(:,:), intent(in)   :: dat
 logical, dimension(:), intent(in)  :: mask
 integer, intent(in)                :: npart
 real, dimension(3), intent(out)    :: partmin,partmax,partmean
 logical, intent(in), optional      :: nonzero
 real, dimension(:), intent(in),  optional :: mass
 real, dimension(3), intent(out), optional :: part_integral
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
 if (present(part_integral)) part_integral(:) = 0.

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
    if (mask(i)) then
       np = np + 1
       do j=1,jlen
          partval = dat(i,j)
          if (usenonzero) then
             if (partval > tiny(0.)) partmin(j) = min(partmin(j),partval)
          else
             partmin(j) = min(partmin(j),partval)
          endif
          partmax(j) = max(partmax(j),partval)
          partmean(j) = partmean(j) + partval
          if (present(mass) .and. present(part_integral)) then
             part_integral(j) = part_integral(j) + mass(i)*partval
          endif
       enddo
    endif
 enddo
 !!$omp end parallel do

 if (np > 0) then
    partmean(:) = partmean(:)/real(np)
 endif

end subroutine minmaxmean_part

!----------------------------------------------------
! calculate max and min and mean values on particles
!----------------------------------------------------
logical function iszero(partmin,partmax,ndim)
 real, dimension(:), intent(in) :: partmin,partmax
 integer, intent(in)            :: ndim

 if (all(abs(partmin(1:ndim)) < tiny(0.)) .and. &
     all(abs(partmax(1:ndim)) < tiny(0.))) then
    iszero = .true.
    print "(a)",' min=max=0 on particles: skipping pointless interpolation and setting dat = 0.'
 else
    iszero = .false.
 endif

end function iszero

end module convert_grid
