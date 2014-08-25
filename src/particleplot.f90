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

module particleplots
 implicit none
 public :: particleplot,plot_errorbarsx,plot_errorbarsy
 public :: plot_kernel_gr

contains
!
!  Drives raw particle plots
!  Handles different particle types, particle cross-sections, particle labelling
!  fast-plotting added 12.10.06 (excludes particles in crowded fields)
!
!  Arguments:
!
!
subroutine particleplot(x,y,z,h,ntot,iplotx,iploty,icolourpart,iamtype,noftype,iplot_type, &
                        use_zrange,zmin,zmax,labelz,xmin,xmax,ymin,ymax, &
                        fast,datpix,nx,ny,dval,brightness)
  use params,           only:int1
  use labels,           only:labeltype, maxparttypes,is_coord
  use settings_data,    only:ndim,icoords,ntypes
  use settings_part,    only:imarktype,ncircpart,icoordsnew,icircpart,itypeorder, &
                             ilabelpart,iplotline,linestylethisstep,linecolourthisstep
  use interpolations2D, only:interpolate_part,interpolate_part1
  use transforms,       only:transform
  use part_utils,       only:igettype
  use sort,             only:indexx
  use plotlib,          only:plot_qci,plot_bbuf,plot_ebuf,plot_sci,plot_sfs,plot_circ, &
                             plot_pt,plot_numb,plot_text,plot_pt1,plot_qls,plot_sls, &
                             plot_line,plot_qlw,plot_slw,plot_errb,plotlib_maxlinestyle
  implicit none
  integer,            intent(in) :: ntot,iplotx, iploty
  integer(kind=int1), intent(in) :: iamtype(:)
  integer,            intent(in) :: icolourpart(:)
  integer,            intent(in) :: noftype(maxparttypes)
  real,               intent(in) :: x(:), y(:), z(:), h(:)
  real,               intent(in) :: zmin,zmax,xmin,xmax,ymin,ymax
  logical,            intent(in) :: use_zrange,fast
  logical,            intent(in) :: iplot_type(maxparttypes)
  character(len=*),   intent(in) :: labelz

  integer,            intent(in),    optional :: nx,ny
  real,               intent(inout), optional :: datpix(:,:),brightness(:,:)
  real,               intent(in),    optional :: dval

  integer :: j,n,itype,linewidth,icolourindex,nplotted,oldlinestyle,ierr
  integer :: lenstring,index1,index2,ntotplot,icolourstart,nlooptypes,ilooptype
  integer             :: nplottedtype(maxparttypes)
  character(len=20)   :: string
  integer, parameter  :: ncellx = 500, ncelly = 500 ! for crowded field reduction
  integer(kind=int1)  :: nincell(ncellx,ncelly,maxparttypes)
  integer             :: icellx,icelly,maxz
  real                :: dx1,dy1,dxpix
  logical             :: mixedtypes
  real, allocatable   :: xerrb(:), yerrb(:), herr(:)

  !--query current character height and colour
  call plot_qci(icolourstart)
  !
  !--check for errors in input
  !
  ntotplot = sum(noftype(1:ntypes))
  if (ntot.lt.ntotplot) then
     print "(a)",' ERROR: number of particles input < number of each type '
     print*,ntot,noftype(1:ntypes)
     return
  elseif (ntot.ne.ntotplot) then
     print "(a)",' WARNING: particleplot: total not equal to sum of types on input'
     print*,' ntotal = ',ntot,' sum of types = ',ntotplot
  endif
  maxz = size(z)
  if (maxz > ntot) maxz = ntot
  if (use_zrange .and. maxz.lt.ntot) then
     print "(a)",' WARNING: particleplot: slice plot but z array too small - excluding particles > z array size'
  endif
  dxpix = 0.
  if (present(datpix)) then
     if (.not.(present(nx).and.present(ny).and.present(dval))) then
        print "(a)",' INTERNAL ERROR in call to particleplot: optional args not present'
        return
     else
        dxpix = (xmax - xmin)/real(nx)
     endif
  endif

  !
  !--loop over all particle types
  !
  index1 = 1
  nplottedtype = 0
  nlooptypes = ntypes
  mixedtypes = size(iamtype).gt.1
  if (mixedtypes .or. use_zrange) nlooptypes = 1
  dx1 = (ncellx - 1)/(xmax-xmin + tiny(xmin))
  dy1 = (ncelly - 1)/(ymax-ymin + tiny(ymin))
  nincell = 0

  over_types: do ilooptype=1,nlooptypes
     call plot_bbuf !--buffer plot output until each particle type finished
     if (mixedtypes .or. use_zrange) then
        index1 = 1
        index2 = ntot
        itype = 0
     else
        itype = itypeorder(ilooptype)
        if (itype.eq.1) then
           index1 = 1
        else
           index1 = sum(noftype(1:itype-1))+1
        endif
        index2 = index1 + noftype(itype) - 1
        if (.not.iplot_type(itype)) then
           call plot_ebuf
           cycle over_types
        endif
     endif
     if (index2.gt.ntot) then
        index2 = ntot
        print "(a)",' WARNING: incomplete data'
     endif
     if (index2.lt.index1) then
        call plot_ebuf
        cycle over_types
     endif

     if (use_zrange) then
        !
        !--if particle cross section, plot particles only in a defined (z) coordinate range
        !
        nplotted = 0        
        overj: do j=1,ntot
           if (mixedtypes) then
              itype = min(max(int(iamtype(j)),1),maxparttypes)
           else
              itype = igettype(j,noftype)
           endif
           if (.not. iplot_type(itype)) cycle overj
           if (j.le.maxz) then
              if (z(j) > zmin .and. z(j) < zmax) then
                 if (icolourpart(j).ge.0) then
                    nplotted = nplotted + 1
                    nplottedtype(itype) = nplottedtype(itype) + 1

                    if (fast .and. noftype(itype) > 100) then
                       if (in_cell(icellx,icelly,x(j),y(j),xmin,ymin,dx1,dy1,ncellx,ncelly)) then
                          if (nincell(icellx,icelly,itype).eq.0) then
                             nincell(icellx,icelly,itype) = nincell(icellx,icelly,itype) + 1_int1  ! this +1 of type int*1
                             call plot_sci(icolourpart(j))
                             call plot_particle(imarktype(itype),x(j),y(j),h(j))
                          endif
                       endif
                    else
                       call plot_sci(icolourpart(j))
                       call plot_particle(imarktype(itype),x(j),y(j),h(j))
                    endif

                    if (present(datpix)) then
                       if (present(brightness)) then
                          call interpolate_part1(x(j),y(j),h(j),xmin,ymin,datpix,nx,ny,dxpix,dval,brightness)
                       else
                          call interpolate_part1(x(j),y(j),h(j),xmin,ymin,datpix,nx,ny,dxpix,dval)
                       endif
                    endif
                 endif
                 !--plot circle of interaction if gas particle
                 if (itype.eq.1 .and. ncircpart.gt.0 .and. ANY(icircpart(1:ncircpart).eq.j)) then
                    call plot_circ(x(j),y(j),2*h(j))
                 endif
                 !!--plot particle label
                 if (ilabelpart) then
                    call plot_numb(j,0,1,string,lenstring)
                    call plot_text(x(j),y(j),string(1:lenstring))
                 endif
              endif
           endif
        enddo overj

        do itype=1,ntypes
           if (iplot_type(itype) .and. nplottedtype(itype).gt.0) then
              if (zmin < -0.1*huge(zmin)) then
                 print*,'plotted ',nplottedtype(itype),' of ',noftype(itype), &
                  trim(labeltype(itype))//' particles with ', trim(labelz),' < ',zmax
              else
                 print*,'plotted ',nplottedtype(itype),' of ',noftype(itype), &
                  trim(labeltype(itype))//' particles in range ', trim(labelz),' = ',zmin,' -> ',zmax                 
              endif
           endif
        enddo
     else
        !
        !--otherwise plot all particles of this type using appropriate marker and colour
        !
        call plot_qci(icolourindex)
        !
        !--all particles in range have same colour and type
        !
        if (.not.mixedtypes .and. all(icolourpart(index1:index2).eq.icolourpart(index1)) &
            .and. icolourpart(index1).ge.0) then
           call plot_sci(icolourpart(index1))
           if (fast .and. (index2-index1).gt.100) then
              !--fast-plotting only allows one particle per "grid cell" - avoids crowded fields
              write(*,"(a,i8,1x,a)") ' fast-plotting ',index2-index1+1,trim(labeltype(itype))//' particles'
              nincell = 0
              do j=index1,index2
                 if (in_cell(icellx,icelly,x(j),y(j),xmin,ymin,dx1,dy1,ncellx,ncelly)) then
                    if (nincell(icellx,icelly,itype).eq.0) then
                       nincell(icellx,icelly,itype) = nincell(icellx,icelly,itype) + 1_int1  ! this +1 of type int*1

                       call plot_particle(imarktype(itype),x(j),y(j),h(j))

                       if (present(datpix)) then
                          if (present(brightness)) then
                             call interpolate_part1(x(j),y(j),h(j),xmin,ymin,datpix,nx,ny,dxpix,dval,brightness)
                          else
                             call interpolate_part1(x(j),y(j),h(j),xmin,ymin,datpix,nx,ny,dxpix,dval)
                          endif
                       endif
                    endif
                 endif
              enddo
           else
              !--plot all particles of this type
              print "(a,i8,1x,a)",' plotting ',index2-index1+1,trim(labeltype(itype))//' particles'
              select case(imarktype(itype))
              case(32:35)
                 do j=1,noftype(itype)
                    call plot_particle(imarktype(itype),x(j),y(j),h(j))
                 enddo
                 call plot_sfs(1)
              case default
                 call plot_pt(noftype(itype),x(index1:index2),y(index1:index2),imarktype(itype))
              end select
              if (present(datpix)) then
                 if (present(brightness)) then
                    call interpolate_part(x(index1:index2),y(index1:index2),h(index1:index2), &
                                          noftype(itype),xmin,ymin,datpix,nx,ny,dxpix,dval,brightness)
                 else
                    call interpolate_part(x(index1:index2),y(index1:index2),h(index1:index2), &
                                          noftype(itype),xmin,ymin,datpix,nx,ny,dxpix,dval)
                 endif
              endif
           endif
        else
        !
        !--mixed colours and/or mixed types
        !
           nplotted = 0
           nplottedtype = 0

           overj2: do j=index1,index2
              if (icolourpart(j).ge.0) then
                 if (mixedtypes) then
                    itype = int(iamtype(j))
                    if (.not.iplot_type(itype)) cycle overj2
                    nplottedtype(itype) = nplottedtype(itype) + 1
                 endif
                 nplotted = nplotted + 1
                 if (fast .and. noftype(itype) > 100) then
                    if (in_cell(icellx,icelly,x(j),y(j),xmin,ymin,dx1,dy1,ncellx,ncelly)) then
                    !--exclude particles if there are more than 2 particles per cell
                    !  (two here because particles can have different colours)
                       if (nincell(icellx,icelly,itype).le.0) then
                          nincell(icellx,icelly,itype) = nincell(icellx,icelly,itype) + 1_int1  ! this +1 of type int*1
                          
                          call plot_sci(icolourpart(j))
                          call plot_particle(imarktype(itype),x(j),y(j),h(j))

                          if (present(datpix)) then
                             if (present(brightness)) then
                                call interpolate_part1(x(j),y(j),h(j),xmin,ymin,datpix,nx,ny,dxpix,dval,brightness)
                             else
                                call interpolate_part1(x(j),y(j),h(j),xmin,ymin,datpix,nx,ny,dxpix,dval)
                             endif
                          endif
                       endif
                    endif
                 else
                    call plot_sci(icolourpart(j))
                    call plot_particle(imarktype(itype),x(j),y(j),h(j))

                    if (present(datpix)) then
                       if (present(brightness)) then
                          call interpolate_part1(x(j),y(j),h(j),xmin,ymin,datpix,nx,ny,dxpix,dval,brightness)
                       else
                          call interpolate_part1(x(j),y(j),h(j),xmin,ymin,datpix,nx,ny,dxpix,dval)
                       endif
                    endif
                 endif
              endif
           enddo overj2
           if (mixedtypes) then
              do itype=1,ntypes
                 if (iplot_type(itype)) then
                    if (fast .and. noftype(itype).gt.100) then
                       print*,' fast-plotted ',nplottedtype(itype),' of ',noftype(itype),trim(labeltype(itype))//' particles'
                    elseif (noftype(itype).gt.0) then
                       print*,' plotted ',nplottedtype(itype),' of ',noftype(itype),trim(labeltype(itype))//' particles'
                    endif
                 endif
              enddo
           else
              if (fast .and. noftype(itype).gt.100) then
                 print*,' fast-plotted ',nplotted,' of ',index2-index1+1,trim(labeltype(itype))//' particles'
              else
                 print*,' plotted ',nplotted,' of ',index2-index1+1,trim(labeltype(itype))//' particles'
              endif
           endif
        endif
        call plot_sci(icolourindex)

        if (ilabelpart) then
           !!--plot particle labels
           print*,'plotting particle labels ',index1,':',index2
           do j=index1,index2
              call plot_numb(j,0,1,string,lenstring)
              call plot_text(x(j),y(j),string(1:lenstring))
           enddo
        endif
     endif
     index1 = index2 + 1
     call plot_ebuf !--flush PGPLOT buffer at end of each type
  enddo over_types

  !
  !--plot lines joining particles if relevant
  !
  call plot_qci(icolourindex)
  call plot_sci(linecolourthisstep)
  !  i.e., don't plot a line for cross section plots (would plot all particles)
  !        but do if there is 3D perspective --> in which case zmin = -huge(x)
  if (iplotline .and. .not.(use_zrange .and. abs(zmax-zmin).lt.0.5*huge(0.))) then
     call plot_qls(oldlinestyle)
     call plot_sls(linestylethisstep)
     call plot_line(noftype(1),x(1:noftype(1)),y(1:noftype(1)))

     if (noftype(2).gt.0 .and. iplot_type(2)) then
        call plot_sls(mod(linestylethisstep+1,plotlib_maxlinestyle) + 1)
        call plot_line(noftype(2),x(noftype(1)+1:sum(noftype(1:2))),y(noftype(1)+1:sum(noftype(1:2))))
     endif
     call plot_sls(oldlinestyle)! reset
  endif

  call plot_sci(icolourindex)
  !
  !--plot circles of interaction (ie a circle of radius 2h)
  !  around all or selected particles. For plots with only one coordinate axis,
  !  these are plotted as error bars in the coordinate direction.
  !
  !--this bit is also used for error bar plotting on x or y axis.
  !
  if (ncircpart.gt.0) then
     !
     !--set fill area style and line width
     !
     call plot_qlw(linewidth)
     call plot_slw(2)
     call plot_qci(icolourindex)
     call plot_sci(2)
     call plot_sfs(2)

     if (ncircpart.gt.0) then

        if (is_coord(iplotx,ndim) .and. is_coord(iploty,ndim) .and. ncircpart.gt.0) then
           print*,'plotting ',ncircpart,' circles of interaction'
           do n = 1,ncircpart
              if (icircpart(n).gt.ntot) then
                 print*,'error: particle index > number of particles'
              else
                 if (icoordsnew.ne.icoords) then
                    call plot_kernel_gr(icoordsnew,icoords,x(icircpart(n)),y(icircpart(n)),2*h(icircpart(n)))
                 else
                    call plot_circ(x(icircpart(n)),y(icircpart(n)),2*h(icircpart(n)))
                 endif
              endif
           enddo
        else
           if (.not.allocated(herr)) then
              allocate(xerrb(ncircpart),yerrb(ncircpart),herr(ncircpart),stat=ierr)
              if (ierr /= 0) &
                 stop ' Error allocating memory in particleplot for circles of interaction'
           endif
           !!--only on specified particles
           do n=1,ncircpart
              if (icircpart(n).gt.ntot) then
                 print*,'error: particle index > number of particles'
                 xerrb(n) = 0.
                 yerrb(n) = 0.
                 herr(n) = 0.
              else
                 xerrb(n) = x(icircpart(n))
                 yerrb(n) = y(icircpart(n))
                 herr(n) = 2.*h(icircpart(n))
              endif
           enddo
           if (is_coord(iplotx,ndim)) then
              print*,'plotting ',ncircpart,' error bars x axis '
              call plot_errb(5,ncircpart,xerrb(1:ncircpart),yerrb(1:ncircpart),herr(1:ncircpart),1.0)
           elseif (is_coord(iploty,ndim)) then
              print*,'plotting ',ncircpart,' error bars y axis'
              call plot_errb(6,ncircpart,xerrb(1:ncircpart),yerrb(1:ncircpart),herr(1:ncircpart),1.0)
           endif
           if (allocated(herr)) deallocate(herr)
           if (allocated(xerrb)) deallocate(xerrb)
           if (allocated(yerrb)) deallocate(yerrb)
        endif
     endif

     call plot_slw(linewidth)
     call plot_sci(icolourindex)

  endif

!
!--reset colour
!
  call plot_sci(icolourstart)

  return

end subroutine particleplot
!--------------------------------------------------------------------------------
!
! subroutine implementing scalable markers
! default case is just an interface to usual particle plotting routine
!
!--------------------------------------------------------------------------------
subroutine plot_particle(imarktype,x,y,h)
 use plotlib,       only:plot_circ,plot_sfs,plot_sci,plot_pt1
 use settings_part, only:hfacmarkers
 implicit none
 integer, intent(in) :: imarktype
 real,    intent(in) :: x,y,h
 integer :: imarker
 real    :: size

 select case(imarktype)
 case(32:35)
    imarker = imarktype - 31
    size = hfacmarkers*h
    if (imarker.le.2) then
       call plot_sfs(imarker)
       call plot_circ(x,y,size)
       call plot_sfs(1)
    elseif (imarker.eq.3) then
       call plot_sfs(1)
       call plot_circ(x,y,size)
       call plot_sfs(2)
       call plot_sci(0)
       call plot_circ(x,y,size)
       call plot_sfs(1)
    elseif (imarker.eq.4) then
       call plot_sfs(1)
       call plot_circ(x,y,size)
       call plot_sfs(2)
       call plot_sci(1)
       call plot_circ(x,y,size)
       call plot_sfs(1)
    else
       call plot_circ(x,y,size)
    endif
 case default
    call plot_pt1(x,y,imarktype)
 end select

end subroutine plot_particle
!------------------------------------------------------------
!
! function used to determine which cell a particle lies in
! returns TRUE if within allowed limits, FALSE if not
!
!------------------------------------------------------------
logical function in_cell(ix,iy,x,y,xmin,ymin,dx1,dy1,nx,ny)
 integer, intent(out) :: ix,iy
 real,    intent(in)  :: x,y,xmin,ymin,dx1,dy1
 integer, intent(in)  :: nx,ny

 ix = int((x - xmin)*dx1) + 1
 iy = int((y - ymin)*dy1) + 1
 !--exclude particles if there are more than 2 particles per cell
 !  (two here because particles can have different colours)
 in_cell = (ix.gt.0 .and. ix.le.nx .and. iy.gt.0 .and. iy.le.ny)

end function in_cell

!--------------------------------------------------------------------------------
!
! subroutine to plot the circle of interaction for a given particle
! in general coordinate systems (e.g. cylindrical coordinates)
!
! input:  igeom : coordinate system (0,1=cartesian, 2=cylindrical, 3=spherical)
!         x,y   : particle location in cartesian space
!         h     : size of smoothing sphere
!                 (assumed isotropic in coordinate space)
!
! PGPLOT page must already be set up - this just draws the "circle"
!
!--------------------------------------------------------------------------------
subroutine plot_kernel_gr(igeom,igeomold,x,y,h)
  use geometry, only:coord_transform,maxcoordsys,labelcoordsys
  use plotlib,  only:plot_line
  implicit none
  integer, intent(in) :: igeom, igeomold
  real, intent(in) :: x,y,h

  integer, parameter :: npts = 100 ! big enough to give a smooth circle
  real, parameter :: pi = 3.1415926536
  integer :: i
  real, dimension(2) :: xtemp
  real, dimension(2,npts) :: xpts
  real :: angle, dangle, xi, yi

  if (igeom.gt.1 .and. igeom.le.maxcoordsys) then
     print 10,labelcoordsys(igeom)
  else
     print 10,labelcoordsys(1)
  endif
10 format('coordinate system = ',a)

  xtemp(1) = x
  xtemp(2) = y
  !--e.g. from cylindricals TO cartesians
  call coord_transform(xtemp,2,igeom,xpts(:,1),2,igeomold)
  xi = xpts(1,1)
  yi = xpts(2,1)
!
!--step around a circle in co-ordinate space of radius h and store the
!  location of the points in cartesian space in the 2D array xpts
!
  dangle = 2.*pi/REAL(npts-1)
  do i=1,npts
     angle = (i-1)*dangle
     xtemp(1) = xi + h*COS(angle)
     xtemp(2) = yi + h*SIN(angle)
!
!--translate back to actual coordinate system plotted
!
     call coord_transform(xtemp,2,igeomold,xpts(:,i),2,igeom)
  enddo
!
!--now plot the circle using pgline
!
  call plot_line(npts,xpts(1,:),xpts(2,:))

  return
end subroutine plot_kernel_gr

!--------------------------------------------------------------------------------
!
!  Plot y-axis error bars, handling the case where the axes are transformed
!
!  input x,y are in transformed space (i.e., already logged)
!  input err is not transformed, (i.e., not logged)
!
!--------------------------------------------------------------------------------
subroutine plot_errorbarsy(npts,x,y,err,itrans)
 use plotlib,    only:plot_bbuf,plot_ebuf,plot_err1,plot_errb
 use transforms, only:transform,transform_inverse,islogged
 use settings_part, only:ErrorBarType
 use settings_data, only:iverbose
 implicit none
 integer, intent(in) :: npts,itrans
 real, intent(in), dimension(:) :: x,y,err
 real :: yval,errval
 real, dimension(2) :: val
 real, dimension(npts) :: errp,errm
 integer :: i

 if (iverbose >= 1) then
    if (npts < 10000) then
       print "(a,i4,a)",' plotting ',npts,' error bars y axis' 
    else
       print "(a,i10,a)",' plotting ',npts,' error bars y axis'
    endif
 endif
 if (itrans.ne.0) then
    if (islogged(itrans)) then
       errval = 0. !-300.
    else
       errval = 0.
    endif
    !call plot_bbuf
    do i=1,npts
       yval = y(i)
       call transform_inverse(yval,itrans)
       val(1) = yval + err(i)
       val(2) = yval - err(i)
       call transform(val,itrans,errval=errval)
       errp(i) = val(1) - y(i)
       errm(i) = y(i) - val(2)
       val(1) = val(1) - y(i)
       val(2) = y(i) - val(2)
       if (ErrorBarType /= 1) then
          call plot_err1(2,x(i),y(i),val(1),1.0)
          call plot_err1(4,x(i),y(i),val(2),1.0)
       endif
    enddo
    if (ErrorBarType==1) then
       call plot_errb(7,npts,x,y,errp,1.0)
       call plot_errb(8,npts,x,y,errm,1.0)
    endif
    !call plot_ebuf
 else
    if (ErrorBarType==1) then
       call plot_errb(9,npts,x,y,err,1.0)
    else
       call plot_errb(6,npts,x,y,err,1.0)
    endif
 endif

end subroutine plot_errorbarsy

!--------------------------------------------------------------------------------
!
!  Plot x-axis error bars, handling the case where the axes are transformed
!
!  input x,y are in transformed space (i.e., already logged)
!  input err is not transformed, (i.e., not logged)
!
!--------------------------------------------------------------------------------
subroutine plot_errorbarsx(npts,x,y,err,itrans)
 use transforms, only:transform,transform_inverse,islogged
 use plotlib,    only:plot_bbuf,plot_ebuf,plot_err1,plot_errb
 implicit none
 integer, intent(in) :: npts,itrans
 real, intent(in), dimension(:) :: x,y,err
 real :: xval,errval
 real, dimension(2) :: val
 integer :: i

 print*,'plotting ',npts,' error bars x axis '
 if (itrans.ne.0) then
    if (islogged(itrans)) then
       errval = -300.
    else
       errval = 0.
    endif
    call plot_bbuf
    do i=1,npts
       xval = x(i)
       call transform_inverse(xval,itrans)
       val(1) = xval + err(i)
       val(2) = xval - err(i)
       call transform(val,itrans,errval=errval)
       val(1) = val(1) - x(i)
       val(2) = x(i) - val(2)
       call plot_err1(1,x(i),y(i),val(1),1.0)
       call plot_err1(3,x(i),y(i),val(2),1.0)
    enddo
    call plot_ebuf
 else
    call plot_errb(5,npts,x,y,err,1.0)
 endif

end subroutine plot_errorbarsx

end module particleplots
