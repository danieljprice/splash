!!
!!   calculates various additional quantities from the input data
!!
subroutine calc_quantities(ifromstep,itostep)
  use labels
  use particle_data
  use settings_data
  use mem_allocation
  implicit none
  integer, intent(in) :: ifromstep, itostep
  integer :: i,j,nstartfromcolumn
  integer :: ientrop
  integer :: ipmag,ibeta,itotpr,idivBerr,icurr,icrosshel
  integer :: irad2,ivpar,ivperp,iBpar,iBperp
  real :: Bmag, Jmag
  real, parameter :: pi = 3.1415926536
  real :: angledeg,anglexy,runit(3)  ! to plot r at some angle
  !
  !--initialise extra quantities to zero
  !
  ientrop = 0
  ike = 0
  ipmag = 0
  ibeta = 0
  itotpr = 0
  idivBerr = 0
  icurr = 0
  icrosshel = 0
  irad2 = 0
  ivpar = 0
  ivperp = 0
  iBpar = 0
  iBperp = 0
  !
  !--specify which of the possible quantities you would like to calculate
  !  (0 = not calculated)
  !
  !--specify hydro quantities
  !
  ncalc = 3
  ientrop = ncolumns + 1      
  irad = ncolumns + 2
  ike = ncolumns + 3
  if (ipr.eq.0 .or. ipr.gt.ncolumns) then
     nstartfromcolumn = ncolumns + ncalc
     ncalc = ncalc + 1
     ipr = nstartfromcolumn + 1
  endif
  !
  !--specify MHD quantities
  !
  if (iBfirst.ne.0) then
     nstartfromcolumn = ncolumns + ncalc
     ncalc = ncalc + 5
     ipmag = nstartfromcolumn + 1
     ibeta = nstartfromcolumn + 2
     itotpr = nstartfromcolumn + 3      
     idivBerr = nstartfromcolumn + 4
     icrosshel = nstartfromcolumn + 5
     !!icurr = ncolumns + 8
     icurr = 0
  else
     nstartfromcolumn = ncolumns + ncalc
     ncalc = ncalc + ndimV
     ivpar = nstartfromcolumn + 1
     ivperp= nstartfromcolumn + 2
     ipmag = 0
     ibeta = 0
     itotpr = 0
     idivBerr = 0
  endif

  if (ndim.eq.2 .and. iBfirst.ne.0 .and. ivx.ne.0) then
     nstartfromcolumn = ncolumns + ncalc
     ncalc = ncalc + 5
     irad2 = nstartfromcolumn + 1
     ivpar = nstartfromcolumn + 2
     ivperp = nstartfromcolumn + 3
     iBpar = nstartfromcolumn + 4
     iBperp = nstartfromcolumn + 5
!  else
!     irad2 = 0
!     ivpar = 0
!     ivperp = 0
!     iBpar = 0
!     iBperp = 0         
  endif
  
  print*,'calculating ',ncalc,' additional quantities...',ifromstep,itostep
  numplot = ncolumns + ncalc
  if (numplot.gt.maxcol) call alloc(maxpart,maxstep,numplot) 

  do i=ifromstep,itostep
     !!--pressure if not in data array
     if ((ipr.gt.ncolumns).and.(irho.ne.0).and.(iutherm.ne.0)) then
        dat(1:ntot(i),ipr,i) = dat(1:ntot(i),irho,i)*dat(1:ntot(i),iutherm,i)*(gamma(i)-1.)
     endif
     !!--entropy
     if (ientrop.ne.0 .and. ipr.ne.0) then
        where (dat(1:ntot(i),irho,i).gt.1.e-10) 
           dat(1:ntot(i),ientrop,i) = dat(1:ntot(i),ipr,i)/dat(1:ntot(i),irho,i)**gamma(i)
        elsewhere
           dat(1:ntot(i),ientrop,i) = 0.
        endwhere
     endif
     !!--radius         
     if (irad.ne.0) then
        if (icoords.gt.1) then  
           dat(1:ntot(i),irad,i) = dat(1:ntot(i),ix(1),i)
        else
           do j=1,ntot(i)
              dat(j,irad,i) = sqrt(dot_product(dat(j,ix(1:ndim),i),   &
                dat(j,ix(1:ndim),i)))
           enddo        
        endif
     endif
     !!--specific KE
     if ((ike.ne.0).and.(ivx.ne.0)) then
        do j=1,ntot(i)
           dat(j,ike,i) = 0.5*dot_product(dat(j,ivx:ivx+ndimV-1,i), &
                                          dat(j,ivx:ivx+ndimV-1,i))
        enddo
     endif

     !!--distance along the line with angle 30deg w.r.t. x axis
     if (irad2.ne.0) then
        angledeg = 30.
        anglexy = angledeg*pi/180.
        runit(1) = cos(anglexy)
        if (ndim.ge.2) runit(2) = sin(anglexy)
        do j=1,ntot(i)
           dat(j,irad2,i) = dot_product(dat(j,ix(1:ndim),i),runit(1:ndim))
        enddo
        if (ivpar.ne.0) then
           dat(1:ntot(i),ivpar,i) = dat(1:ntot(i),ivx+1,i)*SIN(anglexy) + dat(1:ntot(i),ivx,i)*COS(anglexy)
        endif
        if (ivperp.ne.0) then
           dat(1:ntot(i),ivperp,i) = dat(1:ntot(i),ivx+1,i)*COS(anglexy) - dat(1:ntot(i),ivx,i)*SIN(anglexy)           
        endif
        if (iBpar.ne.0) then
           dat(1:ntot(i),iBpar,i) = dat(1:ntot(i),iBfirst+1,i)*SIN(anglexy) + dat(1:ntot(i),iBfirst,i)*COS(anglexy)
        endif
        if (iBperp.ne.0) then
           dat(1:ntot(i),iBperp,i) = dat(1:ntot(i),iBfirst+1,i)*COS(anglexy) - dat(1:ntot(i),iBfirst,i)*SIN(anglexy)
        endif
     endif
     !
     !--magnetic quantities
     !
     if (iBfirst.ne.0) then
        !!--magnetic pressure
        if (ipmag.ne.0) then
           do j=1,ntot(i)
              dat(j,ipmag,i) = 0.5*dot_product(dat(j,iBfirst:iBfirst+ndimV-1,i), &
                                               dat(j,iBfirst:iBfirst+ndimV-1,i))
           enddo
           !!--plasma beta
           if (ibeta.ne.0 .and. ipr.ne.0) then
              where(abs(dat(1:ntot(i),ipmag,i)).gt.1.e-10)
                 dat(1:ntot(i),ibeta,i) = dat(1:ntot(i),ipr,i)/dat(1:ntot(i),ipmag,i)
              elsewhere  
                 dat(1:ntot(i),ibeta,i) = 0.
              endwhere
           endif
        endif
           
        !!--total pressure (gas + magnetic)     
        if ((itotpr.ne.0).and.(ipr.ne.0).and.(ipmag.ne.0)) then
           dat(1:ntot(i),itotpr,i) = dat(1:ntot(i),ipr,i) + dat(1:ntot(i),ipmag,i)
        endif
        !!--div B error        (h*divB / abs(B))
        if ((idivBerr.ne.0).and.(idivB.ne.0).and.(ih.ne.0)) then
           do j=1,ntot(i)
              Bmag = sqrt(dot_product(dat(j,iBfirst:iBfirst+ndimV-1,i), &
                                      dat(j,iBfirst:iBfirst+ndimV-1,i)))
              if (Bmag.gt.0.) then
                 dat(j,idivBerr,i) = abs(dat(j,idivB,i))*dat(j,ih,i)/Bmag
              else
                 dat(j,idivBerr,i) = 0.
              endif
           enddo
        endif
        if ((icurr.ne.0).and.(iJfirst.ne.0).and.irho.ne.0) then
           do j=1,ntot(i)
              Jmag = sqrt(dot_product(dat(j,iJfirst:iJfirst+ndimV-1,i), &
                   dat(j,iJfirst:iJfirst+ndimV-1,i)))
              if (dat(j,irho,i).ne.0) then
                 dat(j,icurr,i) = Jmag !!/sqrt(dat(irho,j,i))
              endif
           enddo
        endif
        if ((icrosshel.ne.0).and.(iBfirst.ne.0).and.ivx.ne.0) then
           do j=1,ntot(i)
              dat(j,icrosshel,i) = dot_product(dat(j,iBfirst:iBfirst+ndimV-1,i),  &
                   dat(j,ivx:ivx+ndimV-1,i))
           enddo
        endif
     endif
  enddo
  !
  !--set labels for calculated quantities
  !
  if (ientrop.ne.0) label(ientrop) = 'entropy'
  if (irad.ne.0) label(irad) = 'radius '
  if (irad2.ne.0) label(irad2) = 'r\d\(0737)'    !!!parallel'
  if (ike.ne.0) label(ike) = 'specific KE'
  if (ipr.gt.ncolumns) label(ipr) = 'P_gas'
  if (ipmag.ne.0) label(ipmag) = 'P_mag'
  if (itotpr.ne.0) label(itotpr) = 'P_gas + P_mag'
  if (ibeta.ne.0) label(ibeta) = 'plasma \gb'
  if (idivberr.ne.0) label(idivberr) = 'h |div B| / |B|'
  if (icurr.ne.0) label(icurr) = '|J|'
  if (icrosshel.ne.0) label(icrosshel) = 'B dot v'
  
  !
  !--calculate the vector quantities in the new co-ordinate basis
  !
  select case(icoords)
  case(2)
     if (ivpar.ne.0) label(ivpar) = 'v_r'
     if (ivperp.ne.0) label(ivperp) = 'v_phi'
     if (iBpar.ne.0) label(iBpar) = 'B_r'
     if (iBperp.ne.0) label(iBperp) = 'B_phi'  
  case default
     if (ivpar.ne.0) label(ivpar) = 'v\d\(0737)'  !!!_parallel'
     if (ivperp.ne.0) label(ivperp) = 'v\d\(0738)' !!_perp'
     if (iBpar.ne.0) label(iBpar) = 'B\d\(0737)'  !!_parallel'
     if (iBperp.ne.0) label(iBperp) = 'B\d\(0738)' !!_perp'
  end select
  
  return
end subroutine calc_quantities
