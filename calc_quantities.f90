!!
!!   calculates various additional quantities from the input data
!!
subroutine calc_quantities
  use labels
  use particle_data
  use settings
  implicit none
  integer :: i,j,icurr,icrosshel
  real :: Bmag, Jmag
  real, parameter :: pi = 3.1415926536
  real :: angledeg,anglexy,runit(3)  ! to plot r at some angle
  !
  !--specify which of the possible quantities you would like to calculate
  !  (0 = not calculated)
  ncalc = 8	! specify number to calculate
  ientrop = ncolumns + 1      
  irad = ncolumns + 2
  ike = ncolumns + 3
  if (iBfirst.ne.0) then
     ipmag = ncolumns + 4
     ibeta = ncolumns + 5
     itotpr = ncolumns + 6      
     idivBerr = ncolumns + 7
     icrosshel = ncolumns + 8
     !!icurr = ncolumns + 8
     icurr = 0
  else
     ncalc = 3 + ndimV
     ivpar = ncolumns + 4
     ivperp= ncolumns + 5
     ipmag = 0
     ibeta = 0
     itotpr = 0
     idivBerr = 0
  endif
  itimestep = 0
  if (ndim.eq.2 .and. iBfirst.ne.0 .and. ivx.ne.0) then
     ncalc = ncalc + 5
     irad2 = ncolumns + 9
     ivpar = ncolumns + 10
     ivperp = ncolumns + 11
     iBpar = ncolumns + 12
     iBperp = ncolumns + 13
!  else
!     irad2 = 0
!     ivpar = 0
!     ivperp = 0
!     iBpar = 0
!     iBperp = 0	 
  endif
  
  print*,'calculating ',ncalc,' additional quantities...'
  numplot = ncolumns + ncalc
  if (numplot.gt.maxcol) call alloc(maxpart,maxstep,numplot) 

  do i=nstart,n_end
     do j=1,ntot(i)
        !!--pressure if not in data array
        if ((ipr.gt.ncolumns)                            &
             .and.(irho.ne.0).and.(iutherm.ne.0))         &
             dat(ipr,j,i) = dat(irho,j,i)*dat(iutherm,j,i)*(gamma(i)-1.)
        !!--entropy
        if (ientrop.ne.0) then
           if (dat(irho,j,i).gt.1.e-10) then
              dat(ientrop,j,i) = dat(ipr,j,i)/dat(irho,j,i)**gamma(i)
           else  
              dat(ientrop,j,i) = 0.
           endif
        endif
        !!--radius	 
        if (irad.ne.0) &    
             dat(irad,j,i) = sqrt(dot_product(dat(ix(1:ndim),j,i),   &
             dat(ix(1:ndim),j,i)))
        !!--distance along the line with angle 30deg w.r.t. x axis
        if (irad2.ne.0) then
           angledeg = 30.
           anglexy = angledeg*pi/180.
           runit(1) = cos(anglexy)
           if (ndim.ge.2) runit(2) = sin(anglexy)
           dat(irad2,j,i) = dot_product(dat(ix(1:ndim),j,i),runit(1:ndim))
           if (ivpar.ne.0) then
              dat(ivpar,j,i) = dat(ivx+1,j,i)*SIN(anglexy) + dat(ivx,j,i)*COS(anglexy)
           endif
           if (ivperp.ne.0) then
              dat(ivperp,j,i) = dat(ivx+1,j,i)*COS(anglexy) - dat(ivx,j,i)*SIN(anglexy)           
           endif
           if (iBpar.ne.0) then
              dat(iBpar,j,i) = dat(iBfirst+1,j,i)*SIN(anglexy) + dat(iBfirst,j,i)*COS(anglexy)
           endif
           if (iBperp.ne.0) then
              dat(iBperp,j,i) = dat(iBfirst+1,j,i)*COS(anglexy) - dat(iBfirst,j,i)*SIN(anglexy)           
           endif
        endif
        !!--specific KE
        if ((ike.ne.0).and.(ivx.ne.0)) &
             dat(ike,j,i) = 0.5*dot_product(dat(ivx:ivlast,j,i),dat(ivx:ivlast,j,i))
!
!--magnetic quantities
!
        if (iBfirst.ne.0) then
           !!--magnetic pressure
           if (ipmag.ne.0) then
              dat(ipmag,j,i) = 0.5*dot_product(dat(iBfirst:iBlast,j,i), &
                   dat(iBfirst:iBlast,j,i))
              !!--plasma beta
              if (ibeta.ne.0) then
                 if (abs(dat(ipmag,j,i)).gt.1e-10) then
                    dat(ibeta,j,i) = dat(ipr,j,i)/dat(ipmag,j,i)
                 else  
                    dat(ibeta,j,i) = 0.
                 endif
              endif
           endif
           
           !!--total pressure (gas + magnetic)     
           if ((itotpr.ne.0).and.(ipr.ne.0).and.(ipmag.ne.0)) then
              dat(itotpr,j,i) = dat(ipr,j,i) + dat(ipmag,j,i)
           endif
           !!--div B error	(h*divB / abs(B))	
           if ((idivBerr.ne.0).and.(idivB.ne.0).and.(ih.ne.0)) then
              Bmag = sqrt(dot_product(dat(iBfirst:iBlast,j,i),dat(iBfirst:iBlast,j,i)))
              if (Bmag.gt.0.) then
                 dat(idivBerr,j,i) = abs(dat(idivB,j,i))*dat(ih,j,i)/Bmag
              else
                 dat(idivBerr,j,i) = 0.
              endif
           endif
           !!--h*SQRT(rho)/abs(B)	(timestep)
           if ((itimestep.ne.0).and.(ih.ne.0).and.(irho.ne.0)) then
              Bmag = sqrt(dot_product(dat(iBfirst:iBlast,j,i),  &
                   dat(iBfirst:iBlast,j,i)))
              if (Bmag.ne.0.) then
                 dat(itimestep,j,i) = dat(ih,j,i)*sqrt(dat(irho,j,i))/Bmag
              else
                 dat(itimestep,j,i) = 0.
              endif
           endif
	   if ((icurr.ne.0).and.(iJfirst.ne.0).and.irho.ne.0) then
	      Jmag = sqrt(dot_product(dat(iJfirst:iJfirst+ndimV,j,i), &
	                              dat(iJfirst:iJfirst+ndimV,j,i)))
	      if (dat(irho,j,i).ne.0) then
	         dat(icurr,j,i) = Jmag     !!/sqrt(dat(irho,j,i))
              endif
	   endif
	   if ((icrosshel.ne.0).and.(iBfirst.ne.0).and.ivx.ne.0) then
	      dat(icrosshel,j,i) = dot_product(dat(iBfirst:iBlast,j,i),  &
                   dat(ivx:ivlast,j,i))
	   endif
	endif
     enddo
  enddo
  !
  !--set labels for calculated quantities
  !
  if (ientrop.ne.0) label(ientrop) = 'entropy'
  if (irad.ne.0) label(irad) = 'radius '
  if (irad2.ne.0) label(irad2) = 'r\d\(0737)'    !!!parallel'
  if (ike.ne.0) label(ike) = 'specific KE'
  if (ipr.ne.0) label(ipr) = 'P'   !'p_gas '
  if (ipmag.ne.0) label(ipmag) = 'P_mag'
  if (itotpr.ne.0) label(itotpr) = 'P_gas + P_mag'
  if (ibeta.ne.0) label(ibeta) = 'plasma \gb'
  if (idivberr.ne.0) label(idivberr) = 'h |div B| / |B|'
  if (itimestep.ne.0) label(itimestep) = 'h sqrt(\gr) / |B|'
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
