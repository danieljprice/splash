module calcquantities
 implicit none
 public :: calc_quantities
 private
 
contains

!!
!!   calculates various additional quantities from the input data
!!
subroutine calc_quantities(ifromstep,itostep)
  use labels
  use particle_data, only:dat,npartoftype,gamma,maxpart,maxstep,maxcol
  use settings_data, only:ndim,ndimV,ncolumns,ncalc,icoords,unitslabel
  use settings_part, only:iexact
  use mem_allocation, only:alloc
  implicit none
  integer, intent(in) :: ifromstep, itostep
  integer :: i,j,nstartfromcolumn,ncolsnew
  integer :: ientrop,idhdrho,ivalfven,imach,ideltarho,ivol
  integer :: ipmag,ibeta,itotpr,idivBerr,icrosshel
  integer :: irad2,ivpar,ivperp,iBpar,iBperp,ntoti
  integer :: iamvecprev,ivec,nveclist,ivecstart,inewcol
  integer, dimension(ncolumns) :: iveclist,ivecmagcol
  real :: Bmag, veltemp, spsound
  real, parameter :: pi = 3.1415926536
  real :: angledeg,anglexy,runit(3)  ! to plot r at some angle
  !
  !--initialise extra quantities to zero
  !
  ientrop = 0
  ike = 0
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
  !
  !--specify which of the possible quantities you would like to calculate
  !  (0 = not calculated)
  !
  !--specify hydro quantities
  !
  ncalc = 0
  !--radius
  if (ndim.gt.0) then
     nstartfromcolumn = ncolumns + ncalc
     ncalc = ncalc + 1
     irad = nstartfromcolumn + 1
  endif
  !--entropy
  if (irho.ne.0 .and. iutherm.ne.0) then
     nstartfromcolumn = ncolumns + ncalc
     ncalc = ncalc + 1
     ientrop = nstartfromcolumn + 1
  endif
  !--pressure
  if ((ipr.eq.0 .or. ipr.gt.ncolumns) .and. iutherm.ne.0.and.irho.ne.0) then
     nstartfromcolumn = ncolumns + ncalc
     ncalc = ncalc + 1
     ipr = nstartfromcolumn + 1
  endif
  !--mach number
  if (ivx.ne.0 .and. ipr.ne.0 .and. irho.ne.0) then
     nstartfromcolumn = ncolumns + ncalc
     ncalc = ncalc + 1
     imach = nstartfromcolumn + 1
  endif
  !--deltarho for toy star
  if (irho.ne.0 .and. irad.ne.0 .and. iexact.eq.4) then
     nstartfromcolumn = ncolumns + ncalc
     ncalc = ncalc + 1
     ideltarho = nstartfromcolumn + 1
  endif
  !--mean particle spacing (m/rho)**(1/ndim)
  !if (ipmass.ne.0 .and. irho.ne.0 .and. ndim.ge.1) then
  !   nstartfromcolumn = ncolumns + ncalc
  !   ncalc = ncalc + 1
  !   ivol = nstartfromcolumn + 1
  !endif
  !--dh/drho
  !if (ih.ne.0 .and. irho.ne.0) then
  !  nstartfromcolumn = ncolumns + ncalc
  !   ncalc = ncalc + 1
  !   idhdrho = nstartfromcolumn + 1
  !endif
  !
  !--specify MHD quantities
  !
  if (iBfirst.ne.0) then
!     nstartfromcolumn = ncolumns + ncalc
!     ncalc = ncalc + 1
!     ipmag = nstartfromcolumn + 1
     if (ipr.ne.0 .and. ipmag.ne.0) then
        nstartfromcolumn = ncolumns + ncalc
        ncalc = ncalc + 2
        ibeta = nstartfromcolumn + 1
        itotpr = nstartfromcolumn + 2
     endif
     if (idivB.ne.0 .and. ih.ne.0) then
        nstartfromcolumn = ncolumns + ncalc
        ncalc = ncalc + 1
        idivBerr = nstartfromcolumn + 1
     endif
!     if (ivx.ne.0) then
!        nstartfromcolumn = ncolumns + ncalc
!        ncalc = ncalc + 1
!        icrosshel = nstartfromcolumn + 1
!     endif
     if (ipmag.ne.0 .and. (irho.ne.0)) then
        nstartfromcolumn = ncolumns + ncalc
        ncalc = ncalc + 1
        ivalfven = nstartfromcolumn + 1
     endif
!     if (iJfirst.ne.0 .and.irho.ne.0) then
!        nstartfromcolumn = ncolumns + ncalc
!        ncalc = ncalc + 1
!        icurr = nstartfromcolumn + 1
!     endif
  else
!     nstartfromcolumn = ncolumns + ncalc
!     ncalc = ncalc + ndimV
!     ivpar = nstartfromcolumn + 1
!     ivperp= nstartfromcolumn + 2
  endif

!  if (ndim.eq.2 .and. iBfirst.ne.0 .and. ivx.ne.0) then
!     nstartfromcolumn = ncolumns + ncalc
!     ncalc = ncalc + 5
!     irad2 = nstartfromcolumn + 1
!     ivpar = nstartfromcolumn + 2
!     ivperp = nstartfromcolumn + 3
!     iBpar = nstartfromcolumn + 4
!     iBperp = nstartfromcolumn + 5
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
           ncalc = ncalc + 1
           nveclist = nveclist + 1
           iveclist(nveclist) = iamvec(i)
           ivecmagcol(nveclist) = ncolumns + ncalc
           iamvecprev = iamvec(i)
        endif
     enddo
  endif
  
  print*,'calculating ',ncalc,' additional quantities...'
  ncolsnew = ncolumns + ncalc
  if (ncolsnew.gt.maxcol) call alloc(maxpart,maxstep,ncolsnew) 

  do i=ifromstep,itostep
     ntoti = SUM(npartoftype(:,i))
     !!--pressure if not in data array
     if ((ipr.gt.ncolumns).and.(irho.ne.0).and.(iutherm.ne.0)) then
        dat(1:ntoti,ipr,i) = dat(1:ntoti,irho,i)*dat(1:ntoti,iutherm,i)*(gamma(i)-1.)
     endif
     !!--entropy
     if (ientrop.ne.0 .and. ipr.ne.0) then
        where (dat(1:ntoti,irho,i).gt.tiny(dat)) 
           dat(1:ntoti,ientrop,i) = dat(1:ntoti,ipr,i)/dat(1:ntoti,irho,i)**gamma(i)
        elsewhere
           dat(1:ntoti,ientrop,i) = 0.
        endwhere
     endif
     !!--mach number
     if (imach.ne.0) then
        do j=1,ntoti
           veltemp = dot_product(dat(j,ivx:ivx+ndimV-1,i), &
                                 dat(j,ivx:ivx+ndimV-1,i))
           if (dat(j,irho,i).gt.tiny(dat)) then
              spsound = gamma(i)*dat(j,ipr,i)/dat(j,irho,i)
              dat(j,imach,i) = sqrt(veltemp/spsound)
           else
              dat(j,imach,i) = 0.
           endif
        enddo
     endif
     !!--dh/drho
     if (idhdrho.ne.0) then
        where (dat(1:ntoti,irho,i).gt.tiny(dat)) 
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
           do j=1,ntoti
              dat(j,irad,i) = sqrt(dot_product(dat(j,ix(1:ndim),i),   &
                dat(j,ix(1:ndim),i)))
           enddo        
        endif
     endif
     !!--specific KE
     if ((ike.ne.0).and.(ivx.ne.0)) then
        do j=1,ntoti
           dat(j,ike,i) = 0.5*dot_product(dat(j,ivx:ivx+ndimV-1,i), &
                                          dat(j,ivx:ivx+ndimV-1,i))
        enddo
     endif
     !!--volume - (m/rho)**(1/ndim)
     if (ivol.ne.0) then
        where (dat(1:ntoti,irho,i).gt.tiny(dat))
           dat(1:ntoti,ivol,i) = (dat(1:ntoti,ipmass,i)/dat(1:ntoti,irho,i))**(1./real(ndim))    
        else where
           dat(1:ntoti,ivol,i) = 0.
        end where
     endif
     !!--delta rho for toy star
     if ((ideltarho.ne.0).and.(irho.ne.0).and.(irad.ne.0)) then
        do j=1,ntoti
           dat(j,ideltarho,i) = dat(j,irho,i) - (1.-dat(j,irad,i)**2)
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
              where(abs(dat(1:ntoti,ipmag,i)).gt.tiny(dat))
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
              if (Bmag.gt.0.) then
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
           where (dat(:,irho,i).gt.tiny(dat))
              dat(:,ivalfven,i) = sqrt(dat(:,ipmag,i)/dat(:,irho,i))
           elsewhere
              dat(:,ivalfven,i) = 0.
           end where
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
  !
  !--set labels for calculated quantities
  !  also units label where dimensions have not changed
  !
  if (ientrop.ne.0) label(ientrop) = 'entropy'
  if (idhdrho.ne.0) label(idhdrho) = 'dh/d\gr'
  if (irad.ne.0) then
     label(irad) = 'radius '
     unitslabel(irad) = unitslabel(ix(1))
  endif
  if (irad2.ne.0) label(irad2) = 'r\d\(0737)'    !!!parallel'
  if (ike.ne.0) label(ike) = 'v\u2\d/2'
  if (ipr.gt.ncolumns) label(ipr) = 'P_gas'
  if (imach.ne.0) label(imach) = '|v|/c\ds'
  if (ideltarho.ne.0) label(ideltarho) = '\gd \gr'
  if (ivol.ne.0) label(ivol) = '(m/rho)^1/ndim'
  
  if (ipmag.ne.0) label(ipmag) = 'B\u2\d/2'
  if (itotpr.ne.0) label(itotpr) = 'P_gas + P_mag'
  if (ibeta.ne.0) label(ibeta) = 'plasma \gb'
  if (idivberr.ne.0) label(idivberr) = 'h |div B| / |B|'
  if (icrosshel.ne.0) label(icrosshel) = 'B dot v'
  if (ivalfven.ne.0) label(ivalfven) = 'v\dalfven'
  !
  !--magnitudes of all vector quantities
  !
  do ivec=1,nveclist
     label(ivecmagcol(ivec)) = '|'//trim(labelvec(iveclist(ivec)))//'|'
  enddo
  
  if (ivpar.ne.0) label(ivpar) = 'v\d\(0737)'  !!!_parallel'
  if (ivperp.ne.0) label(ivperp) = 'v\d\(0738)' !!_perp'
  if (iBpar.ne.0) label(iBpar) = 'B\d\(0737)'  !!_parallel'
  if (iBperp.ne.0) label(iBperp) = 'B\d\(0738)' !!_perp'
  
  return
end subroutine calc_quantities

end module calcquantities
