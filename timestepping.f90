module timestepping
 implicit none
 public :: timestep_loop
 private :: get_nextstep
 
contains
!
! This subroutine drives the main plotting loop
!
subroutine timestep_loop(ipicky,ipickx,irender,ivecplot)
  use particle_data, only:ntot,npartoftype,time,gamma,dat
  use settings_data, only:nstart,n_end,nfreq,numplot,buffer_data
  use settings_page, only:interactive
  use timestep_plotting, only:initialise_plotting,plotstep
  implicit none
  integer, intent(in) :: ipicky,ipickx,irender,ivecplot
  integer :: i, ifile, iadvance
  
  call initialise_plotting(ipicky,ipickx,irender)

  !------------------------------------------------------------------------      
  ! loop over timesteps (flexible to allow going forwards/backwards in
  !                      interactive mode)
  !------------------------------------------------------------------------            
  i = nstart
  iadvance = nfreq   ! amount to increment timestep by (changed in interactive)
  ifile = 1

  over_timesteps: do while (i.le.n_end)

     if (.not.buffer_data) then    
        !
        !--make sure we have data for this timestep
        !
        call get_nextstep(i,ifile)
        if (i.eq.-666) exit over_timesteps
     endif
     !
     !--check timestepping
     !
     if (i.lt.1) then
        print*,'reached first step: can''t go back'
        i = 1
     endif
     if (i.lt.nstart) then
        print*,'warning: i < nstart'
     endif

     print 33, time(i),i
33   format (5('-'),' t = ',f9.4,', dump #',i5,1x,18('-'))

     call plotstep(i,irender,ivecplot,npartoftype(:,i), &
                   dat(:,:,i),time(i),gamma(i),iadvance)
!
!--increment timestep
!
     if (iadvance.eq.-666) exit over_timesteps
     i = i + iadvance
     if (interactive .and. i.gt.n_end) i = n_end

  enddo over_timesteps

  if (.not.interactive) then
     !!call pgebuf
     print*,'press return to finish'
     read*
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call pgend

  return

end subroutine timestep_loop

!-------------------------------------------------------------------
! works out whether or not we need to read another dump into memory
!-------------------------------------------------------------------
subroutine get_nextstep(i,ifile)
 use filenames
 use settings_page, only:interactive
 implicit none
 integer, intent(inout) :: i,ifile
 integer :: iskipfiles

 !
 !--if data is not stored in memory, read next step from file
 !  skip files if necessary. At the moment assumes number of steps in
 !  each file are the same
 !
 !  request is for timestep i from file ifile
 !  if i > number of steps in this file, tries to skip
 !  appropriate number of files to get to timestep requested
 !
 if (i.gt.nstepsinfile(ifile)) then
    if (nstepsinfile(i).ge.1) then
       iskipfiles = (i-nstepsinfile(ifile))/nstepsinfile(ifile)
    else
       print*,'*** error in timestepping: file contains zero timesteps'
       iskipfiles = 1
    endif
    if (iskipfiles.gt.1) then
       print*,'skipping ',iskipfiles,' files '
    elseif (iskipfiles.le.0) then
       print*,'error with iskipfiles = ',iskipfiles
       iskipfiles = 1
    endif
    ifile = ifile+iskipfiles
 elseif (i.lt.1) then
    print*,' stepping back'
    ifile = ifile-1
    if (ifile.ge.1) then
       iskipfiles = (i-1)/nstepsinfile(ifile)
       if (abs(iskipfiles).gt.0) print*,'skipping back ',abs(iskipfiles),' files'
       ifile = ifile + iskipfiles + 1
       if (ifile.lt.1) ifile = 1
    else
       ifile = 1
    endif
    i = 1
 endif

 if (ifile.gt.nfiles) then
    if (interactive) then  ! freeze on last step
       ifile = nfiles
       i = nstepsinfile(ifile)
    else ! exit timestepping loop
       ifile = nfiles
       i = -666
    endif
 elseif (ifile.lt.1) then
    print*,'*** get_nextstep: error: request for file < 1'
 elseif (ifile.ne.ifileopen) then
    call get_data(ifile,.true.)
    i = MOD(i-1,nstepsinfile(ifile)) + 1
    if (i.ne.1) then
       print*,'starting at step ',i
    endif
 endif

 return
end subroutine get_nextstep     

end module timestepping
