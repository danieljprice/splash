subroutine set_options_from_dataread
 use exact,          only:read_exactparams
 use settings_part,  only:iplotpartoftype,UseTypeInRenderings
 use particle_data,  only:npartoftype
 use settings_data,  only:ndim,ncolumns,nsteps,ifileopen,ntypes,iexact
 use labels,         only:labeltype,irho,ih
 use asciiutils,     only:ucase
 use settings_render, only:icolour_particles

 integer   :: itype, nplot

 !
 !--check for errors in data read / print warnings
 !
  if (ndim /= 0 .and. ncolumns > 0 .and. nsteps > 0 .and. iverbose==1) then
     if (sum(npartoftype(:,1)) > 0 .and. npartoftype(1,1)==0 .and. .not.any(iplotpartoftype(2:))) then
        print "(/,a)",' WARNING! DATA APPEARS TO CONTAIN NO '//trim(ucase(labeltype(1)))//' PARTICLES'
        itype = 0
        nplot = 0
        do while (nplot==0 .and. itype < ntypes)
           itype = itype + 1
           if (npartoftype(itype,1) > 0) then
              iplotpartoftype(itype) = .true.
              if (UseTypeInRenderings(itype)) nplot = nplot + npartoftype(itype,1)
              print "(a)",' (plotting of '//trim(labeltype(itype))//' particles turned ON)'
           endif
        enddo
        print*
     endif
        ! Set particle colouring if it can't render
        if (irho==0 .or. ih==0) then
           icolour_particles = .true.
        endif
  endif



   !
   !--read exact solution parameters from files if present
   !
   if (iexact /= 0) then
      if (ifileopen > 0) then
         call read_exactparams(iexact,rootname(ifileopen),ierr)
      endif
   endif

end subroutine set_options_from_dataread
