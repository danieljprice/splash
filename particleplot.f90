!
!  Drives raw particle plots
!  Handles different particle types, particle cross-sections, particle labelling
!
!  Arguments:
!
!
subroutine particleplot(xplot,yplot,zplot,h,ntot,iplotx,iploty,npartoftype,x_sec,xsecmin,xsecmax)
  use labels
  use settings_data ! ndim and icoords
  use settings_part
  implicit none
  integer, intent(in) :: ntot, iplotx, iploty
  integer, dimension(maxparttypes), intent(in) :: npartoftype
  real, dimension(ntot), intent(in) :: xplot, yplot, zplot, h
  real, dimension(ntot) :: xerrb, yerrb, herr
  real, intent(in) :: xsecmin,xsecmax
  logical, intent(in) :: x_sec
  integer :: j,n,itype,linewidth,icolourindex
  integer :: lenstring,index1,index2
  real :: charheight
  character(len=20) :: string
  
  !--query current character height
  call pgqch(charheight)
  print "(a,i8)",' entering particle plot, total particles = ',ntot
  !
  !--loop over all particle types
  !
  index1 = 1
  over_types: do itype=1,ntypes
     index2 = index1 + npartoftype(itype) - 1
     if (iplotpartoftype(itype) .and. npartoftype(itype).gt.0) then
        if (x_sec .and. ndim.eq.3) then
           !
           !--if particle cross section, plot particles only in a defined (z) coordinate range
           !
           print "(a,i8,a,f7.2,a,f7.2)",' plotting ',npartoftype(itype),trim(labeltype(itype))//' particles in range '// &
                ' z = ',xsecmin,' -> ',xsecmax
           do j=index1,index2
              if (zplot(j).lt.xsecmax .and. zplot(j).gt.xsecmin) then
                 call pgpt(1,xplot(j),yplot(j),imarktype(itype))
                 !--plot circle of interaction if gas particle
                 if (itype.eq.1 .and. ncircpart.gt.0 .and. ANY(icircpart(1:ncircpart).eq.j)) then
                    call pgcirc(xplot(j),yplot(j),2*h(j))
                 endif
                 !!--plot particle label
                 if (ilabelpart) then
                    call pgnumb(j,0,1,string,lenstring)
                    call danpgsch(4.0,2)
                    call pgtext(xplot(j),yplot(j),string(1:lenstring))
                    call pgsch(charheight)
                 endif
              endif
           enddo
        else
           !
           !--otherwise plot all particle of this type using appropriate marker and colour
           !
           print "(a,i8,1x,a)",' plotting ',npartoftype(itype),trim(labeltype(itype))//' particles'
           call pgpt(npartoftype(itype),xplot(index1:index2),yplot(index1:index2),imarktype(itype))
           if (ilabelpart) then
              !!--plot particle labels
              print*,'plotting particle labels ',index1,':',index2
              do j=index1,index2
                 call pgnumb(j,0,1,string,lenstring)
                 call danpgsch(4.0,2)
                 call pgtext(xplot(j),yplot(j),string(1:lenstring))
                 call pgsch(charheight)
              enddo
           endif
        endif
     endif
     index1 = index2 + 1
  enddo over_types

  !
  !--plot circles of interaction (ie a circle of radius 2h)
  !  around all or selected particles. For plots with only one coordinate axis, 
  !  these are plotted as error bars in the coordinate direction.
  !
  if (ncircpart.gt.0) then
     !
     !--set fill area style and line width
     !
     call pgsfs(2)
     call pgqlw(linewidth)
     call pgslw(2)
     call pgqci(icolourindex)
     call pgsci(2)
     
     if (iplotx.le.ndim .and. iploty.le.ndim) then
        print*,'plotting circles of interaction',ncircpart
        do n = 1,ncircpart
           if (icircpart(n).gt.ntot) then 
              print*,'error: particle index > number of particles'
           else
              if (icoords.gt.1) then   
                 call plot_kernel_gr(icoords,xplot(icircpart(n)),  &
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
           call pgerrb(6,ncircpart,xerrb(icircpart), &
                yplot(1:ncircpart),herr(1:ncircpart),1.0)      
        endif
     endif
     
     call pgslw(linewidth)
     call pgsci(icolourindex)
     
  endif

  return
     
end subroutine particleplot
