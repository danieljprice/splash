!----------------------------------------------------------------------
! sets options relating to cross sectioning / rotation
!----------------------------------------------------------------------
subroutine options_xsecrotate
 use labels
 use limits, only:lim
 use prompting
 use settings_data ! for ndim
 use settings_xsecrot
 implicit none
 integer :: ians,i
 character(len=1) :: char
 logical :: interact
 
 if (ndim.eq.1) print*,' WARNING: none of these options have any effect in 1D'
 ians = 0
 interact = .true.
 print 10,xsec_nomulti,xsecpos_nomulti,irotate
10  format(' 0) exit ',/,                 &
           ' 1) toggle cross section/projection           (',L1,' )',/, &
           ' 2) set cross section position                (',f5.2,' )',/, &
           ' 3) toggle rotation                           (',L1,' )',/, &
           ' 4) change rotation options')
 call prompt('enter option',ians,0,4)
!
!--options
!
 select case(ians)
!------------------------------------------------------------------------
 case(1)
    xsec_nomulti = .not.xsec_nomulti 
    print *,' Cross section = ',xsec_nomulti
!------------------------------------------------------------------------
 case(2)
    flythru = .false.
    if (ndim.eq.3) then
       call prompt('Do you want a fly-through',flythru)
       if (.not.flythru) then
          call prompt('enter co-ordinate location of cross section slice', &
               xsecpos_nomulti)
       endif
    elseif (ndim.eq.2) then
       call prompt('set cross section position interactively?',interact)
       
       if (interact) then
       !
       !--set cross section position interactively
       !
          call pgbegin(0,'/xw',1,1)
          call pgenv(lim(1,1),lim(1,2),lim(2,1),lim(2,2),1,0)
          call pgcurs(xseclineX1,xseclineY1,char)
          print*,'please select cross section line'
          call pgband(1,1,xseclineX1,xseclineY1,xseclineX2,xseclineY2,char)
          print*,'cross section line: xmin = ',xseclineX1,' xmax = ',xseclineX2
          print*,'                    ymin = ',xseclineY1,' ymax = ',xseclineY2
          call pgend       
       else
       !
       !--set position manually
       !
          if (abs(xseclineX2-xseclineX1).lt.1.e-5 .and. &
              abs(xseclineY2-xseclineY1).lt.1.e-5) then
          !--if not already set (ie. if all = 0.0)
          !  then set default line to diagonal across the domain
             xseclineX1 = lim(1,1)
             xseclineX2 = lim(1,2)
             xseclineY1 = lim(2,1)
             xseclineY2 = lim(2,2)
          endif
          print*,'please set position of cross section through 2D data:'
          call prompt('enter xmin of cross section line',xseclineX1)
          call prompt('enter xmax of cross section line',xseclineX2)
          call prompt('enter ymin of cross section line',xseclineY1)
          call prompt('enter ymax of cross section line',xseclineY2)
       endif
    else
       print*,'WARNING: this option has no effect in 1D'
    endif
!------------------------------------------------------------------------
 case(3,4)
    if (ians.eq.3) irotate = .not.irotate
    print*,'rotate = ',irotate
    if (irotate .or. ians.eq.4) then
       call prompt('enter rotation angle about z axis (deg)',anglez,0.,360.)
       if (ndim.eq.3) then
          call prompt('enter rotation angle about y axis (deg)',angley,0.,360.)
          call prompt('enter rotation angle about z axis (deg)',anglex,0.,360.)
       endif
       !xorigin(1:ndim) = 0.5*(lim(1:ndim,1) + lim(1:ndim,2))
       do i=1,ndim
          call prompt('enter location of origin '//trim(label(ix(i))),xorigin(i))
       enddo
    endif
 end select

 return
end subroutine options_xsecrotate
