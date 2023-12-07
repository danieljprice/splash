module initialise
 implicit none
!!
!! initialise some variables
!! this should only be called at code start
!!
contains

subroutine defaults_set_initial
 use filenames, only:rootname
 use labels, only:label,labeltype,iamvec,labelvec,set_default_labels,reset_columnids
 use limits, only:lim,range
 use particle_data, only:maxpart,maxstep,maxcol
 use settings_data, only:UseTypeInRenderings,device
 integer :: i
!
!--limits (could set them to anything but min & max must be different
!          to enable them to be reset interactively if not set elsewhere)
!
 lim(:,1) = 0.
 lim(:,2) = 1.
 range(:,:) = 0.
 call reset_columnids()
 !
 !--filenames
 !
 rootname = ' '
 !
 !--data array sizes
 !
 maxpart = 0
 maxcol = 0
 maxstep = 0
 !
 !--labels
 !
 call set_default_labels(label)

 !  particle types
 labeltype(1) = 'gas'
 do i=2,size(labeltype)
    if (i > 9) then
       write(labeltype(i),"(a,1x,i2)") 'type',i
    else
       write(labeltype(i),"(a,1x,i1)") 'type',i
    endif
 enddo
 UseTypeInRenderings(:) = .false.
 UseTypeInRenderings(1) = .true.

 !  vector labels
 iamvec(:) = 0
 labelvec = ' '

 !  device from command line
 device = ' '

end subroutine defaults_set_initial

end module initialise
