module particleplots
 implicit none
 public :: particleplot
 private :: plot_kernel_gr
 
contains
!
!  Drives raw particle plots
!  Handles different particle types, particle cross-sections, particle labelling
!  fast-plotting added 12.10.06 (excludes particles in crowded fields)
!
!  Arguments:
!
!
subroutine particleplot(xplot,yplot,zplot,h,ntot,iplotx,iploty, &
                        icolourpart,iamtype,npartoftype,iplotpartoftype, &
                        use_zrange,zmin,zmax,labelz,xmin,xmax,ymin,ymax, &
                        fastparticleplot,datpix,npixx,npixy,dval,brightness)
  use params, only:int1
  use labels, only:labeltype, maxparttypes
  use settings_data, only:ndim,icoords,ntypes
  use settings_part, only:imarktype,ncircpart,icoordsnew,icircpart,itypeorder, &
                          ilabelpart,iplotline,linestylethisstep,linecolourthisstep
  use interpolations2D, only:interpolate_part,interpolate_part1
  implicit none
  integer, intent(in) :: ntot, iplotx, iploty
  integer(kind=int1), dimension(:), intent(in) :: iamtype
  integer, intent(in), dimension(ntot) :: icolourpart
  integer, dimension(maxparttypes), intent(in) :: npartoftype
  real, dimension(ntot), intent(in) :: xplot, yplot, zplot, h
  real, dimension(:), allocatable :: xerrb, yerrb, herr
  real, intent(in) :: zmin,zmax,xmin,xmax,ymin,ymax
  logical, intent(in) :: use_zrange,fastparticleplot
  logical, dimension(maxparttypes), intent(in) :: iplotpartoftype
  character(len=*), intent(in) :: labelz
  
  integer, intent(in), optional :: npixx,npixy
  real, dimension(:,:), intent(inout), optional :: datpix,brightness
  real, intent(in), optional :: dval
  
  integer :: j,n,itype,linewidth,icolourindex,nplotted,oldlinestyle
  integer :: lenstring,index1,index2,ntotplot,icolourstart,nlooptypes,ilooptype
  integer, dimension(maxparttypes) :: nplottedtype
!  real :: charheight
  character(len=20) :: string
  integer, parameter :: ncellx = 500, ncelly = 500 ! for crowded field reduction
  integer(kind=int1), dimension(ncellx,ncelly) :: nincell
  integer :: icellx,icelly !,notplotted
  real :: dxcell1,dycell1,dxpix
  logical :: mixedtypes
  
  !--query current character height and colour
!  call pgqch(charheight)
  call pgqci(icolourstart)
  !!print "(a,i8)",' entering particle plot, total particles = ',ntot
  !
  !--check for errors in input
  !
  ntotplot = sum(npartoftype(1:ntypes))
  if (ntot.lt.ntotplot) then
     print "(a)",' ERROR: number of particles input < number of each type '
     print*,ntot,npartoftype(1:ntypes)
     return
  elseif (ntot.ne.ntotplot) then
     print "(a)",' WARNING: particleplot: total not equal to sum of types on input'
     print*,' ntotal = ',ntot,' sum of types = ',ntotplot
  endif
  dxpix = 0.
  if (present(datpix)) then
     if (.not.(present(npixx).and.present(npixy).and.present(dval))) then
        print "(a)",' INTERNAL ERROR in call to particleplot: optional args not present'
        return
     else
        dxpix = (xmax - xmin)/real(npixx)
     endif
  endif
  
  !
  !--loop over all particle types
  !
  index1 = 1
  nplottedtype = 0
  nlooptypes = ntypes
  mixedtypes = size(iamtype).gt.1
  if (mixedtypes) nlooptypes = 1
  
  over_types: do ilooptype=1,nlooptypes
     call pgbbuf !--buffer PGPLOT output until each particle type finished
     if (mixedtypes) then
        index1 = 1
        index2 = ntot
        itype = 0
     else
        itype = itypeorder(ilooptype)
        if (itype.eq.1) then
           index1 = 1
        else
           index1 = sum(npartoftype(1:itype-1))+1
        endif
        index2 = index1 + npartoftype(itype) - 1
        if (.not.iplotpartoftype(itype)) then
           call pgebuf
           cycle over_types
        endif
     endif
     if (index2.gt.ntot) then 
        index2 = ntot
        print "(a)",' WARNING: incomplete data'
     endif
     if (index2.lt.index1) then
        call pgebuf
        cycle over_types
     endif

     if (use_zrange) then
        !
        !--if particle cross section, plot particles only in a defined (z) coordinate range
        !
        nplotted = 0
        overj: do j=index1,index2
           if (mixedtypes) then
              itype = min(max(int(iamtype(j)),1),maxparttypes)
              if (.not.iplotpartoftype(itype)) cycle overj
           endif
           if (zplot(j).lt.zmax .and. zplot(j).gt.zmin) then
              if (icolourpart(j).ge.0) then
                 nplotted = nplotted + 1
                 nplottedtype(itype) = nplottedtype(itype) + 1
                 call pgsci(icolourpart(j))
                 call pgpt(1,xplot(j),yplot(j),imarktype(itype))
                 if (present(datpix)) then
                    if (present(brightness)) then
                       call interpolate_part1(xplot(j),yplot(j),h(j),xmin,ymin,datpix, &
                                              npixx,npixy,dxpix,dval,brightness)                    
                    else
                       call interpolate_part1(xplot(j),yplot(j),h(j),xmin,ymin,datpix, &
                                              npixx,npixy,dxpix,dval)
                    endif
                 endif
              endif
              !--plot circle of interaction if gas particle
              if (itype.eq.1 .and. ncircpart.gt.0 .and. ANY(icircpart(1:ncircpart).eq.j)) then
                 call pgcirc(xplot(j),yplot(j),2*h(j))
              endif
              !!--plot particle label
              if (ilabelpart) then
                 call pgnumb(j,0,1,string,lenstring)
                 call pgtext(xplot(j),yplot(j),string(1:lenstring))
              endif
           endif
        enddo overj
        if (mixedtypes) then
           do itype=1,ntypes
              if (iplotpartoftype(itype) .and. nplottedtype(itype).gt.0) then
                 print*,' plotted ',nplottedtype(itype),' of ', &
                   npartoftype(itype),trim(labeltype(itype))//' particles in range ', &
                   trim(labelz),' = ',zmin,' -> ',zmax
              endif
           enddo          
        else
           print*,' plotted ',nplotted,' of ', &
             index2-index1+1,trim(labeltype(itype))//' particles in range ', &
             trim(labelz),' = ',zmin,' -> ',zmax
        endif
     else
        !
        !--otherwise plot all particles of this type using appropriate marker and colour
        !
        call pgqci(icolourindex)
        dxcell1 = (ncellx - 1)/(xmax-xmin + tiny(xmin))
        dycell1 = (ncelly - 1)/(ymax-ymin + tiny(ymin))
        !
        !--all particles in range have same colour and type
        !
        if (.not.mixedtypes .and. all(icolourpart(index1:index2).eq.icolourpart(index1)) &
            .and. icolourpart(index1).ge.0) then
           call pgsci(icolourpart(index1))
           if (fastparticleplot .and. (index2-index1).gt.100) then
              !--fast-plotting only allows one particle per "grid cell" - avoids crowded fields
              write(*,"(a,i8,1x,a)") &
                   ' fast-plotting ',index2-index1+1,trim(labeltype(itype))//' particles'
              nincell(1:ncellx,1:ncelly) = 0
!                 notplotted = 0
              do j=index1,index2
                 icellx = int((xplot(j) - xmin)*dxcell1) + 1
                 icelly = int((yplot(j) - ymin)*dycell1) + 1
                 !--exclude particles if there are more than one particle per cell
                 if (icellx.gt.0 .and. icellx.le.ncellx &
                    .and. icelly.gt.0 .and. icelly.le.ncelly) then
                    if (nincell(icellx,icelly).eq.0) then
                       nincell(icellx,icelly) = nincell(icellx,icelly) + 1_int1  ! this +1 of type int*1
                       call pgpt1(xplot(j),yplot(j),imarktype(itype))
                       if (present(datpix)) then
                          if (present(brightness)) then
                             call interpolate_part1(xplot(j),yplot(j),h(j),xmin,ymin,datpix, &
                                                    npixx,npixy,dxpix,dval,brightness)                          
                          else
                             call interpolate_part1(xplot(j),yplot(j),h(j),xmin,ymin,datpix, &
                                                    npixx,npixy,dxpix,dval)
                          endif
                       endif
!                       else
!                         notplotted = notplotted + 1
                    endif
                 endif
              enddo
!                 write(*,"(a,i7,a)") ' (minus ',notplotted,' in crowded fields)'
           else
              !--plot all particles of this type
              print "(a,i8,1x,a)",' plotting ',index2-index1+1,trim(labeltype(itype))//' particles'
              call pgpt(npartoftype(itype),xplot(index1:index2),yplot(index1:index2),imarktype(itype))
              if (present(datpix)) then
                 if (present(brightness)) then
                    call interpolate_part(xplot(index1:index2),yplot(index1:index2),h(index1:index2), &
                                          npartoftype(itype),xmin,ymin,datpix,npixx,npixy,dxpix,dval,brightness)                 
                 else
                    call interpolate_part(xplot(index1:index2),yplot(index1:index2),h(index1:index2), &
                                          npartoftype(itype),xmin,ymin,datpix,npixx,npixy,dxpix,dval)
                 endif
              endif
           endif
        else
        !
        !--mixed colours and/or mixed types
        !
           nplotted = 0
           nplottedtype = 0
           nincell(1:ncellx,1:ncelly) = 0

           overj2: do j=index1,index2
              if (icolourpart(j).ge.0) then
                 if (mixedtypes) then
                    itype = int(iamtype(j))
                    if (.not.iplotpartoftype(itype)) cycle overj2
                    nplottedtype(itype) = nplottedtype(itype) + 1
                 endif
                 nplotted = nplotted + 1
                 if (fastparticleplot .and. npartoftype(itype).gt.100) then
                    icellx = int((xplot(j) - xmin)*dxcell1) + 1
                    icelly = int((yplot(j) - ymin)*dycell1) + 1
                    !--exclude particles if there are more than 2 particles per cell
                    !  (two here because particles can have different colours)
                    if (icellx.gt.0 .and. icellx.le.ncellx &
                       .and. icelly.gt.0 .and. icelly.le.ncelly) then
                       if (nincell(icellx,icelly).le.0) then
                          nincell(icellx,icelly) = nincell(icellx,icelly) + 1_int1  ! this +1 of type int*1
                          call pgsci(icolourpart(j))
                          call pgpt1(xplot(j),yplot(j),imarktype(itype))
                          if (present(datpix)) then
                             if (present(brightness)) then
                                call interpolate_part1(xplot(j),yplot(j),h(j),xmin,ymin,datpix, &
                                                       npixx,npixy,dxpix,dval,brightness)                             
                             else
                                call interpolate_part1(xplot(j),yplot(j),h(j),xmin,ymin,datpix, &
                                                       npixx,npixy,dxpix,dval)
                             endif
                          endif
!                       else
!                         notplotted = notplotted + 1
                       endif
                    endif
                 else
                    call pgsci(icolourpart(j))
                    call pgpt1(xplot(j),yplot(j),imarktype(itype))
                    if (present(datpix)) then
                       if (present(brightness)) then
                          call interpolate_part1(xplot(j),yplot(j),h(j),xmin,ymin,datpix, &
                                                 npixx,npixy,dxpix,dval,brightness)                       
                       else
                          call interpolate_part1(xplot(j),yplot(j),h(j),xmin,ymin,datpix, &
                                                 npixx,npixy,dxpix,dval)
                       endif
                    endif
                 endif
              endif
           enddo overj2
           if (mixedtypes) then
              do itype=1,ntypes
                 if (iplotpartoftype(itype)) then
                    if (fastparticleplot .and. npartoftype(itype).gt.100) then
                       print*,' fast-plotted ',nplottedtype(itype),' of ',npartoftype(itype),trim(labeltype(itype))//' particles'
                    elseif (npartoftype(itype).gt.0) then
                       print*,' plotted ',nplottedtype(itype),' of ',npartoftype(itype),trim(labeltype(itype))//' particles'
                    endif
                 endif
              enddo     
           else
              if (fastparticleplot .and. npartoftype(itype).gt.100) then
                 print*,' fast-plotted ',nplotted,' of ',index2-index1+1,trim(labeltype(itype))//' particles'
              else
                 print*,' plotted ',nplotted,' of ',index2-index1+1,trim(labeltype(itype))//' particles'
              endif
           endif
        endif
        call pgsci(icolourindex)

        if (ilabelpart) then
           !!--plot particle labels
           print*,'plotting particle labels ',index1,':',index2
           do j=index1,index2
              call pgnumb(j,0,1,string,lenstring)
              call pgtext(xplot(j),yplot(j),string(1:lenstring))
           enddo
        endif
     endif
     index1 = index2 + 1
     call pgebuf !--flush PGPLOT buffer at end of each type
  enddo over_types

  !
  !--plot lines joining particles if relevant
  !
  if (iplotline .and. .not.use_zrange) then
     call pgqls(oldlinestyle)
     call pgqci(icolourindex)
     call pgsls(linestylethisstep)
     call pgsci(linecolourthisstep)

     call pgline(npartoftype(1),xplot(1:npartoftype(1)), &
                 yplot(1:npartoftype(1)))

     call pgsls(oldlinestyle)! reset 
     call pgsci(icolourindex)
  endif
  !
  !--plot circles of interaction (ie a circle of radius 2h)
  !  around all or selected particles. For plots with only one coordinate axis, 
  !  these are plotted as error bars in the coordinate direction.
  !
  if (ncircpart.gt.0) then
     !
     !--set fill area style and line width
     !
     call pgqlw(linewidth)
     call pgslw(2)
     call pgqci(icolourindex)
     call pgsci(2)
     call pgsfs(2)

     if (iplotx.le.ndim .and. iploty.le.ndim) then
        print*,'plotting ',ncircpart,' circles of interaction'
        do n = 1,ncircpart
           if (icircpart(n).gt.ntot) then 
              print*,'error: particle index > number of particles'
           else
              if (icoordsnew.ne.icoords) then   
                 call plot_kernel_gr(icoordsnew,icoords,xplot(icircpart(n)),  &
                      yplot(icircpart(n)),2*h(icircpart(n)))
              else
                 call pgcirc(xplot(icircpart(n)),  &
                      yplot(icircpart(n)),2*h(icircpart(n)))
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
              xerrb(n) = xplot(icircpart(n))
              yerrb(n) = yplot(icircpart(n))
              herr(n) = 2.*h(icircpart(n))
           endif
        enddo         
        if (iplotx.le.ndim) then
           print*,'plotting ',ncircpart,' error bars x axis '
           call pgerrb(5,ncircpart,xerrb(1:ncircpart), &
                yerrb(1:ncircpart),herr(1:ncircpart),1.0)
        elseif (iploty.le.ndim) then
           print*,'plotting ',ncircpart,' error bars y axis'
           call pgerrb(6,ncircpart,xerrb(1:ncircpart), &
                yplot(1:ncircpart),herr(1:ncircpart),1.0)      
        endif
        if (allocated(herr)) deallocate(herr)
        if (allocated(xerrb)) deallocate(xerrb)
        if (allocated(yerrb)) deallocate(yerrb)
     endif
     
     call pgslw(linewidth)
     call pgsci(icolourindex)
     
  endif

!
!--reset colour
!
  call pgsci(icolourstart)

  return
     
end subroutine particleplot

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
subroutine plot_kernel_gr(igeom,igeomold,x,y,h)
  use geometry, only:coord_transform,maxcoordsys,labelcoordsys
  implicit none
  integer, intent(in) :: igeom, igeomold
  real, intent(in) :: x,y,h
  
  integer, parameter :: npts = 100 ! big enough to give a smooth circle
  real, parameter :: pi = 3.1415926536
  integer :: i
  real, dimension(2) :: xtemp
  real, dimension(2,npts) :: xpts
  real :: angle, dangle  

  if (igeom.gt.1 .and. igeom.le.maxcoordsys) then
     print 10,labelcoordsys(igeom) 
  else
     print 10,labelcoordsys(1)
  endif
10 format('coordinate system = ',a)
  
!
!--step around a circle in co-ordinate space of radius h and store the 
!  location of the points in cartesian space in the 2D array xpts
!
  dangle = 2.*pi/REAL(npts-1)
  do i=1,npts
     angle = (i-1)*dangle
     xtemp(1) = x + h*COS(angle) 
     xtemp(2) = y + h*SIN(angle)
!
!--translate back to actual coordinate system plotted
!
     call coord_transform(xtemp,2,igeomold,xpts(:,i),2,igeom)
  enddo
!
!--now plot the circle using pgline
!
  call pgline(npts,xpts(1,:),xpts(2,:))

  return
end subroutine plot_kernel_gr

end module particleplots
