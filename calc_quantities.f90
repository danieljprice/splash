!!
!!   calculates various additional quantities from the input data
!!
subroutine calc_quantities
  use labels
  use particle_data
  use settings
  implicit none
  integer :: i,j
  real :: Bmag
  real, parameter :: pi = 3.1415926536
  real :: angledeg,anglexy,runit(ndimmax)	! to plot r at some angle
  
  print*,'calculating ',ncalc,' additional quantities...'
  numplot = ncolumns + ncalc
  do i=nstart,n_end
     do j=1,npart(i)
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
        !!--magnetic pressure
        if ((ipmag.ne.0).and.(iBfirst.ne.0))   &
             dat(ipmag,j,i) = 0.5*dot_product(dat(iBfirst:iBlast,j,i), &
             dat(iBfirst:iBlast,j,i))
        !!--plasma beta
        if ((ibeta.ne.0).and.(ipmag.ne.0)) then
           if (abs(dat(ipmag,j,i)).gt.1e-10) then
              dat(ibeta,j,i) = dat(ipr,j,i)/dat(ipmag,j,i)
           else  
              dat(ibeta,j,i) = 0.
           endif
        endif
        !!--total pressure (gas + magnetic)     
        if ((itotpr.ne.0).and.(ipr.ne.0).and.(ipmag.ne.0)) then
           dat(itotpr,j,i) = dat(ipr,j,i) + dat(ipmag,j,i)
        endif
        !!--div B error	(h*divB / abs(B))	
        if ((idivBerr.ne.0).and.    &
             (idivB.ne.0).and.(ih.ne.0).and.(iBfirst.ne.0)) then
           Bmag = sqrt(dot_product(dat(iBfirst:iBlast,j,i),dat(iBfirst:iBlast,j,i)))
           if (Bmag.gt.0.) then
              dat(idivBerr,j,i) = abs(dat(idivB,j,i))*dat(ih,j,i)/Bmag
           else
              dat(idivBerr,j,i) = 0.
           endif
        endif
        !!--h*SQRT(rho)/abs(B)	(timestep)
        if ((itimestep.ne.0).and.   &
             (ih.ne.0).and.(irho.ne.0).and.(iBfirst.ne.0)) then
           Bmag = sqrt(dot_product(dat(iBfirst:iBlast,j,i),  &
                dat(iBfirst:iBlast,j,i)))
           if (Bmag.ne.0.) then
              dat(itimestep,j,i) = dat(ih,j,i)*sqrt(dat(irho,j,i))/Bmag
           else
              dat(itimestep,j,i) = 0.
           endif
        endif
        
     enddo
  enddo
  !
  !--set labels for calculated quantities
  !
  if (ientrop.ne.0) label(ientrop) = 'entropy'
  if (irad.ne.0) label(irad) = 'radius '
  if (irad2.ne.0) label(irad2) = 'r_parallel'
  if (ike.ne.0) label(ike) = 'specific KE'
  if (ipr.ne.0) label(ipr) = 'P'	!'p_gas '
  if (ipmag.ne.0) label(ipmag) = 'P_mag'
  if (itotpr.ne.0) label(itotpr) = 'P_gas + P_mag'
  if (ibeta.ne.0) label(ibeta) = 'plasma \gb'
  if (idivberr.ne.0) label(idivberr) = 'h |div B| / |B|'
  if (itimestep.ne.0) label(itimestep) = 'h sqrt(\gr) / |B|'
  if (ivpar.ne.0) label(ivpar) = 'v_parellel'
  if (ivperp.ne.0) label(ivperp) = 'v_perp'
  if (iBpar.ne.0) label(iBpar) = 'B_parallel'
  if (iBperp.ne.0) label(iBperp) = 'B_perp'
  
  return
end subroutine calc_quantities
