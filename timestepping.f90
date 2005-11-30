module timestepping
 implicit none
 public :: timestep_loop
 private :: get_nextstep
 
contains
!
! This subroutine drives the main plotting loop
!
subroutine timestep_loop(ipicky,ipickx,irender,ivecplot)
  use particle_data, only:npartoftype,time,gamma,dat
  use settings_data, only:nstart,n_end,nfreq,DataIsBuffered,iUsesteplist,isteplist
  use settings_page, only:interactive,nstepsperpage,iColourEachStep,iChangeStyles
  use timestep_plotting, only:initialise_plotting,plotstep
  implicit none
  integer, intent(in) :: ipicky,ipickx,irender,ivecplot
  integer :: i, icount, istep, ifile, iadvance, istepsonpage
  logical :: ipagechange
  
  call initialise_plotting(ipicky,ipickx,irender)

  !----------------------------------------------------------------------------  
  ! loop over timesteps (flexible to allow going forwards/backwards in
  !                      interactive mode)
  !
  ! bookkeeping is as follows:
  !            i : current step number or position in steplist array
  !        istep : current or requested location in the dat array
  !                (requested location is sent to get_nextstep which skips files
  !                or steps appropriately and sets actual location)
  !       icount : current step number
  !                (as if all steps were in memory in sequential order)
  ! istepsonpage : number of steps which have been plotted on current page
  !                this is used because steps are coloured/marked differently
  !                from this routine (call to colour_timestep)
  !----------------------------------------------------------------------------         
  i = nstart
  iadvance = nfreq   ! amount to increment timestep by (changed in interactive)
  ifile = 1
  istepsonpage = 0
  icount = nstart

  over_timesteps: do while (i.le.n_end .and. icount.le.n_end)

     if (iUseStepList) then
        istep = isteplist(i)
        icount = isteplist(i)
     else
        istep = i
     endif

     if (.not.DataIsBuffered) then    
        !
        !--make sure we have data for this timestep
        !
        call get_nextstep(istep,ifile)
        if (.not.iUseStepList) then ! istep can be reset to last in file
           i = istep
        endif
        if (istep.eq.-666) exit over_timesteps
     endif
     !
     !--check timestepping
     !
     if (istep.lt.1) then
        print*,'reached first step: can''t go back'
        istep = 1
        i = 1
        icount = 1
     endif
     if (istep.lt.nstart) then
        print*,'warning: i < nstart'
     endif

     print 33, time(istep),icount
33   format (5('-'),' t = ',f9.4,', dump #',i5,1x,18('-'))

     istepsonpage = istepsonpage + 1
     if (nstepsperpage.gt.1 .and. istepsonpage.le.nstepsperpage) then
        ipagechange = .false.
     else
        istepsonpage = 1
        ipagechange = .true.
     endif
     !--colour the timestep if appropriate
     if (nstepsperpage.gt.1 .and. (iColourEachStep .or. iChangeStyles)) then
        call colour_timestep(istepsonpage,iColourEachStep,iChangeStyles)
     endif

     call plotstep(icount,istepsonpage,irender,ivecplot,npartoftype(:,istep), &
                   dat(:,:,istep),time(istep),gamma(istep),ipagechange,iadvance)
!
!--increment timestep
!
     if (iadvance.eq.-666) exit over_timesteps
     i = i + iadvance
     icount = icount + iadvance

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
 use filenames, only:nstepsinfile,nfiles,ifileopen
 use settings_page, only:interactive
 use getdata, only:get_data
 implicit none
 integer, intent(inout) :: i,ifile
 integer :: iskipfiles,ifileprev

 !
 !--if data is not stored in memory, read next step from file
 !  skip files if necessary. At the moment assumes number of steps in
 !  each file are the same
 !
 !  request is for timestep i from file ifile
 !  if i > number of steps in this file, tries to skip
 !  appropriate number of files to get to timestep requested
 !
 ifileprev = ifile
 !!print*,'request step ',i,' from file #',ifile,' nstepsinfile= ',nstepsinfile(ifile)
 
 if (i.gt.nstepsinfile(ifile)) then
    if (nstepsinfile(ifile).ge.1) then
       iskipfiles = (i-nstepsinfile(ifile)-1)/nstepsinfile(ifile)
    else
       print*,'*** error in timestepping: file contains zero timesteps'
       iskipfiles = 0
    endif
    if (iskipfiles.ge.1) then
       print*,'skipping ',iskipfiles,' files '
    elseif (iskipfiles.lt.0) then
       print*,'error with iskipfiles = ',iskipfiles
       iskipfiles = 0
    endif
    ifile = ifile+iskipfiles+1
 elseif (i.lt.1) then
    ifile = ifile-1
    if (ifile.ge.1) then
       iskipfiles = (i)/nstepsinfile(ifile)
       if (abs(iskipfiles).gt.0) print*,'skipping back ',abs(iskipfiles),' files'
       ifile = ifile + iskipfiles
       if (ifile.lt.1) then
          ifile = 1
          print*,'can''t skip back that far, starting at file ',ifile
       endif
    else
       ifile = 1
    endif
 endif

 if (ifile.gt.nfiles) then
    if (interactive) then  ! freeze on last step
       ifile = nfiles
       i = nstepsinfile(ifile)
    else ! if non-interactive, exit timestepping loop
       ifile = nfiles
       i = -666
       return
    endif
 elseif (ifile.lt.1) then
    print*,'*** get_nextstep: error: request for file < 1'
 elseif (ifile.ne.ifileopen) then
    call get_data(ifile,.true.)
    if (i.gt.nstepsinfile(ifileprev)) then
       i = MOD(i-nstepsinfile(ifileprev)-1,nstepsinfile(ifileprev)) + 1
    elseif (i.lt.1) then
       i = nstepsinfile(ifile) + MOD(i,nstepsinfile(ifileprev))
    endif
    if (i.ne.1) then
       print*,'starting at step ',i
    endif
    !!print*,'getting file ',ifile,' step ',i
 endif
 
 return

 return
end subroutine get_nextstep     

!-------------------------------------------------------------
! colours all the particles a given colour for this timestep
! and/or changes the marker type for type 1 particles
!-------------------------------------------------------------
subroutine colour_timestep(istep,iChangeColours,iChangeStyles)
  use particle_data, only:icolourme
  use settings_part, only:linecolourthisstep,linestylethisstep,imarktype
  implicit none
  integer, intent(in) :: istep
  logical, intent(in) :: iChangeColours, iChangeStyles
  integer :: icolour,imarkernumber
  
  if (iChangeColours) then
     if (allocated(icolourme)) then
        icolour = istep + 1
        if (icolour.gt.16) then
           print "(a)",'warning: step colour > 16: re-using colours'
           icolour = mod(icolour-1,16) + 1
        endif
        icolourme = icolour
        linecolourthisstep = icolour
     else
        print "(a)",'***error: array not allocated in colour_timestep***'
     endif
  endif
  if (iChangeStyles) then
     linestylethisstep = mod(istep-1,5) + 1
     imarkernumber = istep
     select case(imarkernumber)
     case(1)
        imarktype(1) = 4
     case(2)
        imarktype(1) = 17
     case(3)
        imarktype(1) = 2
     case(4)
        imarktype(1) = 3
     case(5:16)
        imarktype(1) = imarkernumber
     case(17:)
        imarktype(1) = imarkernumber + 1
     end select
  endif

  return
end subroutine colour_timestep

end module timestepping
