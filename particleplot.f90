!
!  Drives raw particle plots
!  Handles different particle types, particle cross-sections, particle labelling
!
!  Arguments:
!
!
subroutine particleplot(xplot,yplot,zplot,h,ntot,iplotx,iploty,npartoftype,x_sec,xsecmin,xsecmax)
  use labels
  use settings
  implicit none
  integer, intent(in) :: ntot, iplotx, iploty
  integer, dimension(maxparttypes), intent(in) :: npartoftype
  real, dimension(ntot), intent(in) :: xplot, yplot, zplot, h
  real, intent(in) :: xsecmin,xsecmax
  logical, intent(in) :: x_sec
  integer :: i,j,n,itype
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
                 if (itype.eq.1 .and. plotcirc) then
                    call pgcirc(xplot(j),yplot(j),2*h(j))
                 endif
                 !!--plot particle label
                 if (ilabelpart) then
                    call pgnumb(j,0,1,string,lenstring)
                    call pgsch(0.5*charheight)
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
                 call pgsch(0.5*charheight)
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
  if (plotcirc) then
     if (iplotx.le.ndim .and. iploty.le.ndim) then
        if (plotcircall) then
           print*,'plotting circles of interaction'
           do j=1,ntot
              call pgcirc(xplot(j),yplot(j),2*h(j))
           enddo
        else
           print*,'plotting circles of interaction',ncircpart
           do n = 1,ncircpart   
              if (icoords.gt.1) then   
                 call plot_kernel_gr(icoords,xplot(icircpart(n)),  &
                      yplot(icircpart(n)),2*h(icircpart(n)))
              else
                 call pgcirc(xplot(icircpart(n)),  &
                      yplot(icircpart(n)),2*h(icircpart(n)))
              endif
           enddo
        endif

     else
        if (plotcircall) then
           !!--on all particles
           if (iplotx.le.ndim) then
              print*,'plotting error bars x axis',ntot
              call pgerrb(5,ntot,xplot(1:ntot), &
                   yplot(1:ntot),2.*h(1:ntot),1.0)
           elseif (iploty.le.ndim) then
              print*,'plotting error bars y axis',ntot 
              call pgerrb(6,ntot,xplot(1:ntot), &
                   yplot(1:ntot),2.*h(1:ntot),1.0)
           endif
        else 
           !!--only on specified particles
           do n=1,ncircpart
              if (iplotx.le.ndim) then
                 print*,'plotting error bar x axis',icircpart(n)
                 call pgerrb(5,1,xplot(icircpart(n)), &
                      yplot(icircpart(n)), &
                      2.*h(icircpart(n)),1.0)
              elseif (iploty.le.ndim) then
                 print*,'plotting error bar y axis',icircpart(n)
                 call pgerrb(6,1,xplot(icircpart(n)),yplot(icircpart(n)), &
                      2.*h(icircpart(n)),1.0)      
              endif
           enddo
        endif
     endif
  endif

  return
     
end subroutine particleplot
