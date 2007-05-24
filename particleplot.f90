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
                        icolourpart,npartoftype,iplotpartoftype, &
                        use_zrange,zmin,zmax,labelz,xmin,xmax,ymin,ymax, &
                        fastparticleplot)
  use labels, only:labeltype, maxparttypes
  use settings_data, only:ndim,icoords,ntypes
  use settings_part, only:imarktype,ncircpart,icoordsnew,icircpart, &
                          ilabelpart,iplotline,linestylethisstep,linecolourthisstep
  implicit none
  integer, intent(in) :: ntot, iplotx, iploty
  integer, intent(in), dimension(ntot) :: icolourpart
  integer, dimension(maxparttypes), intent(in) :: npartoftype
  real, dimension(ntot), intent(in) :: xplot, yplot, zplot, h
  real, dimension(ntot) :: xerrb, yerrb, herr
  real, intent(in) :: zmin,zmax,xmin,xmax,ymin,ymax
  logical, intent(in) :: use_zrange,fastparticleplot
  logical, dimension(maxparttypes), intent(in) :: iplotpartoftype
  character(len=*), intent(in) :: labelz
  integer :: j,n,itype,linewidth,icolourindex,nplotted,oldlinestyle
  integer :: lenstring,index1,index2,ntotplot,icolourstart
!  real :: charheight
  character(len=20) :: string
  integer, parameter :: ncellx = 500, ncelly = 500 ! for crowded field reduction
  integer(kind=1), dimension(ncellx,ncelly) :: nincell
  integer :: icellx,icelly !,notplotted
  real :: dxcell1,dycell1
  
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
  !
  !--loop over all particle types
  !
  index1 = 1
  over_types: do itype=1,ntypes
     call pgbbuf !--buffer PGPLOT output until each particle type finished
     index2 = index1 + npartoftype(itype) - 1
     if (index2.gt.ntot) then 
        index2 = ntot
        print "(a)",' WARNING: incomplete data'
     endif
     if (index2.lt.index1) then
        call pgebuf
        cycle over_types
     endif

     if (iplotpartoftype(itype) .and. npartoftype(itype).gt.0) then
        if (use_zrange) then
           !
           !--if particle cross section, plot particles only in a defined (z) coordinate range
           !
           nplotted = 0
           do j=index1,index2
              if (zplot(j).lt.zmax .and. zplot(j).gt.zmin) then
                 if (icolourpart(j).ge.0) then
                    nplotted = nplotted + 1
                    call pgsci(icolourpart(j))
                    call pgpt(1,xplot(j),yplot(j),imarktype(itype))
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
           enddo
           print*,' plotted ',nplotted,' of ', &
             index2-index1+1,trim(labeltype(itype))//' particles in range ', &
             trim(labelz),' = ',zmin,' -> ',zmax
        else
           !
           !--otherwise plot all particles of this type using appropriate marker and colour
           !
           call pgqci(icolourindex)
           if (all(icolourpart(index1:index2).eq.icolourpart(index1) &
               .and. icolourpart(index1).ge.0)) then
              call pgsci(icolourpart(index1))
              if (fastparticleplot .and. (index2-index1).gt.100) then
                 !--fast-plotting only allows one particle per "grid cell" - avoids crowded fields
                 write(*,"(a,i8,1x,a)") &
                      ' fast-plotting ',index2-index1+1,trim(labeltype(itype))//' particles'
                 dxcell1 = (ncellx - 1)/(xmax-xmin + tiny(xmin))
                 dycell1 = (ncelly - 1)/(ymax-ymin + tiny(ymin))
                 nincell(1:ncellx,1:ncelly) = 0
!                 notplotted = 0
                 do j=index1,index2
                    icellx = int((xplot(j) - xmin)*dxcell1) + 1
                    icelly = int((yplot(j) - ymin)*dycell1) + 1
                    !--exclude particles if there are more than one particle per cell
                    if (icellx.gt.0 .and. icellx.le.ncellx &
                       .and. icelly.gt.0 .and. icelly.le.ncelly) then
                       if (nincell(icellx,icelly).eq.0) then
                          nincell(icellx,icelly) = nincell(icellx,icelly) + 1_1  ! this +1 of type int*1
                          call pgpt1(xplot(j),yplot(j),imarktype(itype))
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
              endif
           else
              nplotted = 0
              do j=index1,index2
                 if (icolourpart(j).ge.0) then
                    nplotted = nplotted + 1
                    call pgsci(icolourpart(j))
                    call pgpt1(xplot(j),yplot(j),imarktype(itype))
                 endif
              enddo
              print*,' plotted ',nplotted,' of ',index2-index1+1,trim(labeltype(itype))//' particles'
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
