module calcquantities
 implicit none
 public :: calc_quantities
 private
 
contains

!!
!!   calculates various additional quantities from the input data
!!
subroutine calc_quantities(ifromstep,itostep,dontcalculate)
  use labels, only:label,labelvec,iamvec,ix,irho,ih,ipmass,iutherm,ipr,ivx,ike, &
                   irad,iBfirst,idivB,icv,iradenergy
  use particle_data, only:dat,npartoftype,gamma,maxpart,maxstep,maxcol
  use settings_data, only:ndim,ndimV,ncolumns,ncalc,icoords,iRescale,xorigin,itrackpart
  use settings_part, only:iexact
  use mem_allocation, only:alloc
  use settings_units, only:unitslabel,units
  implicit none
  integer, intent(in) :: ifromstep, itostep
  logical, intent(in), optional :: dontcalculate
  integer :: i,j,ncolsnew
  integer :: ientrop,idhdrho,ivalfven,imach,ideltarho,ivol
  integer :: ipmag,ibeta,itotpr,idivBerr,icrosshel,ithermal
  integer :: irad2,ivpar,ivperp,iBpar,iBperp,ntoti
  integer :: iamvecprev,ivec,nveclist,ivecstart,inewcol
  integer :: imri,ipk
  integer :: itempgas,itemprad,idudtrad
  integer, dimension(ncolumns) :: iveclist,ivecmagcol
  logical :: skip
  real :: Bmag, veltemp, spsound, gmw
  real, parameter :: mhonkb = 1.6733e-24/1.38e-16
  real, parameter :: pi = 3.1415926536
  real, parameter :: Omega0 = 1.e-3 ! for MRI delta v
  real :: angledeg,anglexy,runit(3)  ! to plot r at some angle
  real, parameter :: radconst = 7.5646e-15
  real, parameter :: lightspeed = 3.e10   ! in cm/s (cgs)
  
  !
  !--allow dummy call to set labels without actually calculating stuff
  !
  if (present(dontcalculate)) then
     skip = dontcalculate
  else
     skip = .false.
  endif  
  !
  !--initialise extra quantities to zero
  !
  ike = 0
  ientrop = 0
  idhdrho = 0
  ipmag = 0
  ibeta = 0
  itotpr = 0
  idivBerr = 0
  icrosshel = 0
  irad2 = 0
  ivpar = 0
  ivperp = 0
  iBpar = 0
  iBperp = 0
  ivalfven = 0
  imach = 0
  ideltarho = 0
  ivol = 0
  ithermal = 0
  imri = 0
  ipk = 0
  itempgas = 0
  itemprad = 0
  idudtrad = 0
  !
  !--specify which of the possible quantities you would like to calculate
  !  (0 = not calculated)
  !
  !--specify hydro quantities
  !
  ncalc = 0
  
  !--radius
  if (ndim.gt.0) call addcolumn(irad,'radius ')
  if (irad.gt.0) unitslabel(irad) = unitslabel(ix(1))
  !--thermal energy per unit volume
  if (irho.ne.0 .and. iutherm.ne.0) call addcolumn(ithermal,'\gr u')
  !--entropy
  if (irho.ne.0 .and. iutherm.ne.0) call addcolumn(ientrop,'entropy')
  !--pressure
  if ((ipr.eq.0 .or. ipr.gt.ncolumns) .and. iutherm.ne.0.and.irho.ne.0) then
     call addcolumn(ipr,'pressure')
  endif
  !--mach number
  if (ipr.ne.0 .and. irho.ne.0 .and. ivx.ne.0) call addcolumn(imach,'mach number')
  !--deltarho for toy star
  if (irho.ne.0 .and. irad.ne.0 .and. iexact.eq.4) call addcolumn(ideltarho,'\gd \gr')
  !--mean particle spacing (m/rho)**(1/ndim)
  !if (ipmass.ne.0 .and. irho.ne.0 .and. ndim.ge.1) call addcolumn(ivol,'m/rho^(1/ndim)')
  !--dh/drho
  !if (ih.ne.0 .and. irho.ne.0) call addcolumn(idhdrho,'dh/d\gr')
  !--pk
  !call addcolumn(ipk,'P(k) k\u2')

  !
  !--radiative transfer stuff
  !
  if (ndim.gt.0 .and. iutherm.gt.0) call addcolumn(itempgas,'gas temperature')
  if (ndim.gt.0 .and. irho.gt.0 .and. iradenergy.gt.0) call addcolumn(itemprad,'radiation temperature')
  if (ndim.gt.0 .and. irho.gt.0 .and. iradenergy.gt.0 .and. icv.gt.0 .and. iutherm.gt.0 .and. iRescale) &
     call addcolumn(idudtrad,'du/dt\drad\u')
  !
  !--specify MHD quantities
  !
  if (iBfirst.ne.0) then
     call addcolumn(ipmag,'1/2 B\u2\d')
     if (ipr.ne.0 .and. ipmag.ne.0) then
        call addcolumn(ibeta,'plasma \gb')
!        call addcolumn(itotpr,'P_gas + P_mag')
     endif
     if (idivB.ne.0 .and. ih.ne.0) call addcolumn(idivBerr,'h |div B| / |B|')
!     if (ivx.ne.0) call addcolumn(icrosshel,'B dot v')
     if (ipmag.ne.0 .and. (irho.ne.0)) call addcolumn(ivalfven,'v\dalfven\u')
     !--MRI perturbed velocity (delta v)
     if (ivx.ne.0 .and. ndim.ge.2 .and. ndimV.ge.3) call addcolumn(imri,'\gd v\d3\u')
  else
!     call addcolumn(ivpar,'v\d\(0737)')
!     call addcolumn(ivperp,'v\d\(0738)')
  endif

!  if (ndim.eq.2 .and. iBfirst.ne.0 .and. ivx.ne.0) then
!     call addcolumn(irad2,'r\d\(0737)')
!     call addcolumn(ivpar,'v\d\(0737)')
!     call addcolumn(ivperp,'v\d\(0738)')
!     call addcolumn(iBpar,'B\d\(0737)')
!     call addcolumn(iBperp,'B\d\(0738)')
!  endif
!
!--magnitudes of all vector quantities (cartesian only)
!
  iamvecprev = 0
  nveclist = 0
  iveclist(:) = 0
  if (icoords.eq.1) then
     do i=1,ncolumns
        if (iamvec(i).gt.0 .and. iamvec(i).le.ncolumns .and. iamvec(i).ne.iamvecprev) then
           nveclist = nveclist + 1
           iveclist(nveclist) = iamvec(i)
           call addcolumn(ivecmagcol(nveclist),'|'//trim(labelvec(iamvec(i)))//'|')           
           !--set units label for vector magnitudes to be the same as the vectors
           unitslabel(ivecmagcol(nveclist)) = unitslabel(iamvec(i))
           iamvecprev = iamvec(i)
        endif
     enddo
  endif
  
  if (.not.skip) print*,'calculating ',ncalc,' additional quantities...'
  ncolsnew = ncolumns + ncalc
  if (ncolsnew.gt.maxcol) call alloc(maxpart,maxstep,ncolsnew) 
!
!--reset iamvec to zero for calculated columns
!
  iamvec(ncolumns+1:ncolsnew) = 0
  labelvec(ncolumns+1:ncolsnew) = ' '

  if (.not.skip .and. ncalc.gt.0) then
   do i=ifromstep,itostep
      ntoti = SUM(npartoftype(:,i))
      !!--pressure if not in data array
      if ((ipr.gt.ncolumns).and.(irho.ne.0).and.(iutherm.ne.0)) then
         if (gamma(i).gt.1.00001) then
            dat(1:ntoti,ipr,i) = dat(1:ntoti,irho,i)*dat(1:ntoti,iutherm,i)*(gamma(i)-1.)
         else
         !--for isothermal it depends what utherm is set to. This is the case in sphNG:
            dat(1:ntoti,ipr,i) = dat(1:ntoti,irho,i)*dat(1:ntoti,iutherm,i)*2./3.
         endif
      endif
      !!--entropy
      if (ientrop.ne.0 .and. ipr.ne.0) then
         where (dat(1:ntoti,irho,i).gt.tiny(0.)) 
            dat(1:ntoti,ientrop,i) = dat(1:ntoti,ipr,i)/dat(1:ntoti,irho,i)**gamma(i)
         elsewhere
            dat(1:ntoti,ientrop,i) = 0.
         endwhere
      endif
      !!--mach number
      if (imach.ne.0 .and. ivx.ne.0 .and. irho.ne.0 .and. ipr.ne.0) then
         do j=1,ntoti
            veltemp = dot_product(dat(j,ivx:ivx+ndimV-1,i), &
                                  dat(j,ivx:ivx+ndimV-1,i))
            if (dat(j,irho,i).gt.tiny(0.)) then
               spsound = gamma(i)*dat(j,ipr,i)/dat(j,irho,i)
               if (spsound.gt.tiny(spsound)) then
                  dat(j,imach,i) = sqrt(veltemp/spsound)
               else
                  dat(j,imach,i) = 0.
               endif
            else
               dat(j,imach,i) = 0.
            endif
         enddo
      endif
      !!--dh/drho
      if (idhdrho.ne.0) then
         where (dat(1:ntoti,irho,i).gt.tiny(0.)) 
            dat(1:ntoti,idhdrho,i) = &
               -dat(1:ntoti,ih,i)/(ndim*(dat(1:ntoti,irho,i)))
         elsewhere
            dat(1:ntoti,idhdrho,i) = 0.
         endwhere
      endif
      !!--radius         
      if (irad.ne.0) then
         if (icoords.gt.1) then
            dat(1:ntoti,irad,i) = dat(1:ntoti,ix(1),i)
         else
            !--make radius relative to tracked particle if particle tracking is set
            if (itrackpart.gt.0 .and. itrackpart.le.ntoti) then
               print "(a,i10)",' radius relative to particle ',itrackpart
               do j=1,ntoti
                  dat(j,irad,i) = sqrt(dot_product( &
                                  dat(j,ix(1:ndim),i)-dat(itrackpart,ix(1:ndim),i), &
                                  dat(j,ix(1:ndim),i)-dat(itrackpart,ix(1:ndim),i)))
               enddo
            else
            !--calculate radius using origin settings for rotation
               if (any(abs(xorigin(1:ndim)).gt.tiny(xorigin))) then
                  !--only print origin setting if non-zero
                  print*,'radius calculated relative to origin = ',xorigin(1:ndim)
               endif
               do j=1,ntoti
                  dat(j,irad,i) = sqrt(dot_product(dat(j,ix(1:ndim),i)-xorigin(1:ndim),   &
                    dat(j,ix(1:ndim),i)-xorigin(1:ndim)))
               enddo
            endif    
         endif
      endif
      !!--specific KE
      if ((ike.ne.0).and.(ivx.ne.0)) then
         do j=1,ntoti
            dat(j,ike,i) = 0.5*dot_product(dat(j,ivx:ivx+ndimV-1,i), &
                                           dat(j,ivx:ivx+ndimV-1,i))
         enddo
      endif
      if (ithermal.ne.0 .and. iutherm.ne.0 .and. irho.ne.0) then
         dat(1:ntoti,ithermal,i) = dat(1:ntoti,irho,i)*dat(1:ntoti,iutherm,i)
      endif
      !!--volume - (m/rho)**(1/ndim)
      if (ivol.ne.0 .and. ipmass.ne.0 .and. irho.ne.0 .and. ndim.gt.0) then
         where (dat(1:ntoti,irho,i).gt.tiny(0.))
            dat(1:ntoti,ivol,i) = (dat(1:ntoti,ipmass,i)/dat(1:ntoti,irho,i))**(1./real(ndim))    
         elsewhere
            dat(1:ntoti,ivol,i) = 0.
         end where
      endif
      !!--delta rho for toy star
      if ((ideltarho.ne.0).and.(irho.ne.0).and.(irad.ne.0)) then
         do j=1,ntoti
            dat(j,ideltarho,i) = dat(j,irho,i) - (1.-dat(j,irad,i)**2)
         enddo
      endif
      !!--P(k)*k-2
      if ((ipk.ne.0)) then
         do j=1,ntoti
            dat(j,ipk,i) = dat(j,2,i)*dat(j,1,i)**2
         enddo
      endif
      !!--distance along the line with angle 30deg w.r.t. x axis
      if (irad2.ne.0) then
         angledeg = 30.
         anglexy = angledeg*pi/180.
         runit(1) = cos(anglexy)
         if (ndim.ge.2) runit(2) = sin(anglexy)
         do j=1,ntoti
            dat(j,irad2,i) = dot_product(dat(j,ix(1:ndim),i),runit(1:ndim))
         enddo
         if (ivpar.ne.0) then
            dat(1:ntoti,ivpar,i) = dat(1:ntoti,ivx+1,i)*SIN(anglexy) + dat(1:ntoti,ivx,i)*COS(anglexy)
         endif
         if (ivperp.ne.0) then
            dat(1:ntoti,ivperp,i) = dat(1:ntoti,ivx+1,i)*COS(anglexy) - dat(1:ntoti,ivx,i)*SIN(anglexy)           
         endif
         if (iBpar.ne.0) then
            dat(1:ntoti,iBpar,i) = dat(1:ntoti,iBfirst+1,i)*SIN(anglexy) + dat(1:ntoti,iBfirst,i)*COS(anglexy)
         endif
         if (iBperp.ne.0) then
            dat(1:ntoti,iBperp,i) = dat(1:ntoti,iBfirst+1,i)*COS(anglexy) - dat(1:ntoti,iBfirst,i)*SIN(anglexy)
         endif
      endif
      !
      !--radiative transfer quantities
      !
      !!--gas temperature
      if (itempgas.gt.0 .and. ndim.gt.0 .and. iutherm.gt.0) then
         if (icv.gt.0) then
            print*,' TEMP USES CV'
            where(abs(dat(1:ntoti,icv,i)).gt.tiny(0.))
               dat(1:ntoti,itempgas,i) = dat(1:ntoti,iutherm,i)/dat(1:ntoti,icv,i)
            elsewhere  
               dat(1:ntoti,itempgas,i) = 0.
            endwhere
         else
            gmw = 4.0/(2.*0.7 + 0.28)
            print*,' CALCULATING GAS TEMPERATURE USING ASSUMED MU=',gmw
            print*,' PHYSICAL UNITS MUST BE ON TO GET A RESULT IN K'
            dat(1:ntoti,itempgas,i) = 2./3.*gmw*mhonkb*dat(1:ntoti,iutherm,i)
         endif
      endif
      !!--radiation temperature
      if (itemprad.gt.0 .and. ndim.gt.0 .and. irho.gt.0 .and. iradenergy.gt.0) then
         if (iRescale) then
            dat(1:ntoti,itemprad,i) = abs(dat(1:ntoti,irho,i)*dat(1:ntoti,iradenergy,i)/radconst)**0.25
         else ! if not using physical units, still give radiation temperature in physical units
            dat(1:ntoti,itemprad,i) = abs(dat(1:ntoti,irho,i)*units(irho)*dat(1:ntoti,iradenergy,i) &
                                      *units(iradenergy)/radconst)**0.25         
         endif
         if (idudtrad.gt.0) then
            if (iRescale) then
               dat(1:ntoti,idudtrad,i) = lightspeed*dat(1:ntoti,iradenergy+1,i)* &
                                         (abs(dat(1:ntoti,irho,i))*dat(1:ntoti,iradenergy,i) &
                                          - radconst*(dat(1:ntoti,iutherm,i)/dat(1:ntoti,icv,i))**4)
            else
               dat(1:ntoti,idudtrad,i) = 0.
            endif
         endif
      endif
      !
      !--magnetic quantities
      !
      if (iBfirst.ne.0) then
         !!--magnetic pressure
         if (ipmag.ne.0) then
            do j=1,ntoti
               dat(j,ipmag,i) = 0.5*dot_product(dat(j,iBfirst:iBfirst+ndimV-1,i), &
                                                dat(j,iBfirst:iBfirst+ndimV-1,i))
            enddo
            !!--plasma beta
            if (ibeta.ne.0) then
               where(abs(dat(1:ntoti,ipmag,i)).gt.tiny(0.))
                  dat(1:ntoti,ibeta,i) = dat(1:ntoti,ipr,i)/dat(1:ntoti,ipmag,i)
               elsewhere  
                  dat(1:ntoti,ibeta,i) = 0.
               endwhere
            endif
         endif

         !!--total pressure (gas + magnetic)     
         if (itotpr.ne.0) then
            dat(1:ntoti,itotpr,i) = dat(1:ntoti,ipr,i) + dat(1:ntoti,ipmag,i)
         endif
         !!--div B error        (h*divB / abs(B))
         if (idivBerr.ne.0) then
            do j=1,ntoti
               Bmag = sqrt(dot_product(dat(j,iBfirst:iBfirst+ndimV-1,i), &
                                       dat(j,iBfirst:iBfirst+ndimV-1,i)))
               if (Bmag.gt.tiny(Bmag)) then
                  dat(j,idivBerr,i) = abs(dat(j,idivB,i))*dat(j,ih,i)/Bmag
               else
                  dat(j,idivBerr,i) = 0.
               endif
            enddo
         endif
         if (icrosshel.ne.0) then
            do j=1,ntoti
               dat(j,icrosshel,i) = dot_product(dat(j,iBfirst:iBfirst+ndimV-1,i),  &
                    dat(j,ivx:ivx+ndimV-1,i))
            enddo
         endif
         if (ivalfven.ne.0) then
            where (dat(:,irho,i).gt.tiny(0.))
               dat(:,ivalfven,i) = sqrt(dat(:,ipmag,i)/dat(:,irho,i))
            elsewhere
               dat(:,ivalfven,i) = 0.
            end where
         endif
         if (imri.gt.0 .and. ivx.gt.0 .and. ndim.ge.2 .and. ndimV.ge.3) then
            dat(:,imri,i) = dat(:,ivx+2,i) + 1.5*Omega0*dat(:,ix(1),i)
         endif
      endif
      !
      !--magnitudes of all vector quantities
      !
      do ivec=1,nveclist
         inewcol = ivecmagcol(ivec)
         ivecstart = iveclist(ivec)
         do j=1,ntoti
            dat(j,inewcol,i) = sqrt(dot_product(dat(j,ivecstart:ivecstart+ndimV-1,i), &
                                                dat(j,ivecstart:ivecstart+ndimV-1,i)))
         enddo
      enddo

   enddo
   
  endif  ! skip
  !
  !--override units of calculated quantities if necessary
  !
  if (iRescale .and. any(abs(units(ncolumns+1:ncolumns+ncalc)-1.0).gt.tiny(0.)) &
      .and. .not.skip) then
     write(*,"(/a)") ' rescaling data...'
     do i=ncolumns+1,ncolumns+ncalc
        if (abs(units(i)-1.0).gt.tiny(0.) .and. units(i).gt.tiny(0.)) then
           dat(:,i,ifromstep:itostep) = dat(:,i,ifromstep:itostep)*units(i)
           if (index(label(i),trim(unitslabel(i))).eq.0) label(i) = trim(label(i))//trim(unitslabel(i))
        endif
     enddo
  elseif (iRescale) then
     do i=ncolumns+1,ncolumns+ncalc
        if (index(label(i),trim(unitslabel(i))).eq.0) label(i) = trim(label(i))//trim(unitslabel(i))
     enddo
  endif
  
  return
end subroutine calc_quantities

subroutine addcolumn(inewcolumn,labelin)
 use labels, only:label
 use settings_data, only:ncolumns,ncalc 
 implicit none
 integer, intent(out) :: inewcolumn
 character(len=*), intent(in) :: labelin

 ncalc = ncalc + 1
 inewcolumn = ncolumns + ncalc
 if (inewcolumn.le.size(label)) then
    label(inewcolumn) = trim(labelin)
 else
    print*,' WARNING!!! too many columns for array dimensions'
    print*,' => change parameter ''maxplot'' in globaldata.f90 and recompile'
 endif

 return
end subroutine addcolumn

end module calcquantities
