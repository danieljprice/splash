!
!--unit test for interpolation routines
!
program test_interpolation
 use projections3D
 use xsections3D
 implicit none
 integer, parameter :: idimx = 100
 integer, parameter :: idim = idimx**3
 integer, parameter :: ipixx = 1000, ipixy = 1000
 integer :: npart,npartx,nparty,npartz
 integer :: npixx, npixy,i
 real, parameter :: errtol = 1.e-7
 real, dimension(idim) :: x,y,z,pmass,h,rho
 real, dimension(idim) :: dat,weight
 integer, dimension(idim) :: itype
 real, dimension(ipixx,ipixy) :: datpix
 real, dimension(0:maxcoltable) :: q,w
 real :: xmin,xmax,ymin,ymax,zmin,zmax
 real :: columndens,dxpix,err,dens,datmax
 real :: trans(6)
 logical :: ifastrender,normalise
 
 xmin = -0.5
 xmax = 0.5
 ymin = -0.5
 ymax = 0.5
 zmin = -0.5
 zmax = 0.5
 itype = 0
 ifastrender = .true.
 print*,'accelerated rendering = ',ifastrender
 call pgopen('?')
! call pgenv(xmin,xmax,ymin,ymax,0,0)
! call pglabel('x','y',' ')
! call pgpt(npart,x,y,1)
! call pgenv(xmin,xmax,ymin,ymax,0,0)
! call pglabel('y','z',' ')
! call pgpt(npart,y,z,1)
! call pgenv(xmin,xmax,ymin,ymax,0,0)
! call pglabel('x','z',' ')
! call pgpt(npart,x,z,1)
!
!--setup integrated kernel table
!
 call setup_integratedkernel 
!
!--check value of the integration at q=zero (can do this analytically)
!
 if (abs(coltable(0)-1.5/3.1415926536).lt.errtol) then
    print*,'CENTRAL KERNEL TABLE OK'
 else
    print*,'coltable(0) = ',coltable(0),' should be ',2.*0.75/3.1415926536
    print*,'error = ',abs(coltable(0)-1.5/3.1415926536)
    print*,'ERROR: CENTRAL INTEGRATED KERNEL VALUE WRONG'
 endif
!
!--plot integrated kernel
!
 call pgenv(0.,2.,0.,coltable(0)*1.1,0,0)
 call pglabel('r','int W',' ')
 print*,'radkernel = ',radkernel
 
 do i=0,50
    q(i) = i*radkernel/50.
    !print*,q(i)
    w(i) = wfromtable(q(i)*q(i))
 enddo
 call pgsci(3) 
 call pgline(50,q,w)

! do i=1,maxcoltable
!    q(i) = sqrt((i-1)*radkernel*radkernel*dmaxcoltable)
! enddo
! call pgsci(2)
! call pgline(maxcoltable,q,coltable)
 call pgsci(1)
 !call pgpage
!
!--setup one particle
!
 print*,'SINGLE PARTICLE TEST'
 npart = 1
 npixx = 10
 npixy = 10
 x(1) = 0.5*(xmin + xmax)
 y(1) = 0.5*(ymin + ymax)
 z(1) = 0.5*(zmin + zmax)
 rho(1) = 1.0
 pmass(1) = 2.0
 h(1) = 0.35*xmax
 weight(1) = 1./1.5**3
 dat(1) = rho(1)
 dxpix = (xmax-xmin)/real(npixx)
 normalise = .false.
 datpix = 0.
 call interpolate3D_projection(x(1:npart),y(1:npart),z(1:npart),h(1:npart), &
                               weight(1:npart),dat(1:npart),itype(1:npart),npart,xmin,ymin, &
                               datpix(1:npixx,1:npixy),npixx,npixy,dxpix,dxpix,normalise,0.,0.,.false.)
 call pgenv(xmin,xmax,ymin,ymax,1,0)
 trans = 0.
 trans(1) = xmin - 0.5*dxpix
 trans(2) = dxpix
 trans(4) = ymin - 0.5*dxpix
 trans(6) = dxpix
 datmax = maxval(datpix(1:npixx,1:npixy))
 print*,'max datpix = ',datmax,maxloc(datpix(1:npixx,1:npixy))
 if (abs(datmax-0.035367373).gt.errtol) then
    print*,'FAILED: central maximum wrong, error = ',abs(datmax-0.035367373)
 else
    print*,'OK: central maximum seems fine'
 endif
 call pgimag(datpix,ipixx,ipixy,1,npixx,1,npixy,0.0,datmax,trans)
 
 print*,'TEST WITH ACCELERATION'
 call interpolate3D_projection(x(1:npart),y(1:npart),z(1:npart),h(1:npart), &
                               weight(1:npart),dat(1:npart),itype(1:npart),npart,xmin,ymin, &
                               datpix(1:npixx,1:npixy),npixx,npixy,dxpix,dxpix,normalise,0.,0.,.true.)
 call pgenv(xmin,xmax,ymin,ymax,1,0)
 trans = 0.
 trans(1) = xmin - 0.5*dxpix
 trans(2) = dxpix
 trans(4) = ymin - 0.5*dxpix
 trans(6) = dxpix
 datmax = maxval(datpix(1:npixx,1:npixy))
 print*,'max datpix = ',datmax,maxloc(datpix(1:npixx,1:npixy))
 if (abs(datmax-0.035367373).gt.errtol) then
    print*,'FAILED: central maximum wrong, error = ',abs(datmax-0.035367373)
 else
    print*,'OK: central maximum seems fine'
 endif
 call pgimag(datpix,ipixx,ipixy,1,npixx,1,npixy,0.0,datmax,trans)
!
!--setup two overlapping particles
!
 print*,'TWO PARTICLE TEST'
 npart = 2
 npixx = 1000
 npixy = 1000
 x(1) = -0.25
 x(2) = 0.25
 y(2) = 0.5*(ymin + ymax)
 z(2) = 0.5*(zmin + zmax)
 rho(2) = 1.0
 pmass(2) = 2.0
 h(1:2) = 0.5*xmax
 weight(2) = 1./1.5**3
 dat(2) = rho(2)
 dxpix = (xmax-xmin)/real(npixx)
 call interpolate3D_projection(x(1:npart),y(1:npart),z(1:npart),h(1:npart), &
                               weight(1:npart),dat(1:npart),itype(1:npart),npart,xmin,ymin, &
                               datpix(1:npixx,1:npixy),npixx,npixy,dxpix,dxpix,normalise,0.,0.,ifastrender)
 call pgenv(xmin,xmax,ymin,ymax,1,0)
 trans = 0.
 trans(1) = xmin - 0.5*dxpix
 trans(2) = dxpix
 trans(4) = ymin - 0.5*dxpix
 trans(6) = dxpix
 datmax = maxval(datpix(1:npixx,1:npixy))
 print*,'max datpix = ',datmax,maxloc(datpix(1:npixx,1:npixy))

 !!print*,'datpix = ',datpix(1:npixx,1:npixy)
 call pgimag(datpix,ipixx,ipixy,1,npixx,1,npixy,0.0,datmax,trans)
!
!--set up a cubic lattice of particles
! 
 print*,'NORMAL LATTICE TEST'
 npartx = 50
 nparty = 50
 npartz = 50
 npart = npartx*nparty*npartz
 npixx = 500
 npixy = 500
 dxpix = (xmax-xmin)/real(npixx)
 call setgrid(npartx,nparty,npartz,x,y,z,pmass,rho,h,weight,xmin,xmax,ymin,ymax,zmin,zmax) 
!
!--now call interpolation routine to pixels
!
 call interpolate3D_projection(x(1:npart),y(1:npart),z(1:npart),h(1:npart), &
                               weight(1:npart),dat(1:npart),itype(1:npart),npart,xmin,ymin, &
                               datpix(1:npixx,1:npixy),npixx,npixy,dxpix,dxpix,normalise,0.,0.,ifastrender)
!
!--check output
!
 dens = rho(1)
 columndens = dens*(zmax-zmin)
 call geterr(datpix(1:npixx,1:npixy),npixx,npixy,columndens,err)
 print "(70('-'))"
 print*,'average error in column density interpolation = ',err
 if (err.gt.0.05) then
    print*,'FAILED: average error > usual'
 else
    print*,'OK: average error same as usual'
 endif

 call pgenv(xmin,xmax,ymin,ymax,0,0)
! call pgpixl(datpix,ipixx,ipixy,1,npixx,1,npixy,xmin,xmax,ymin,ymax)
 trans = 0.
 trans(1) = xmin - 0.5*dxpix
 trans(2) = dxpix
 trans(4) = ymin - 0.5*dxpix
 trans(6) = dxpix
 call pgimag(datpix,ipixx,ipixy,1,npixx,1,npixy,0.0,1.0,trans)
!
!--NORMALISED VERSION OF ABOVE
!
 normalise = .true.
 call interpolate3D_projection(x(1:npart),y(1:npart),z(1:npart),h(1:npart), &
                               weight(1:npart),dat(1:npart),itype(1:npart),npart,xmin,ymin, &
                               datpix(1:npixx,1:npixy),npixx,npixy,dxpix,dxpix,normalise,0.,0.,ifastrender)
!
!--check output
!
 dens = rho(1)
 call geterr(datpix(1:npixx,1:npixy),npixx,npixy,dens,err)
 print "(70('-'))"
 print*,'average error in <density> interpolation = ',err
 print*,' dens = ',dens,' datpix = ',datpix(1:10,1:10)
 if (err.gt.0.05) then
    print*,'FAILED: average error > usual'
 else
    print*,'OK: average error same as usual'
 endif
 call pgenv(xmin,xmax,ymin,ymax,0,0)
! call pgpixl(datpix,ipixx,ipixy,1,npixx,1,npixy,xmin,xmax,ymin,ymax)
 trans = 0.
 trans(1) = xmin - 0.5*dxpix
 trans(2) = dxpix
 trans(4) = ymin - 0.5*dxpix
 trans(6) = dxpix
 call pgimag(datpix,ipixx,ipixy,1,npixx,1,npixy,0.0,1.0,trans)
!
!--take cross section at midplane and check density
!
 print "(70('-'))"
 call interpolate3D_fastxsec(x(1:npart),y(1:npart),z(1:npart), &
      h(1:npart),weight(1:npart),dat(1:npart),itype(1:npart),npart,&
      xmin,ymin,0.0,datpix(1:npixx,1:npixy),npixx,npixy,dxpix,.false.)     

 call geterr(datpix(1:npixx,1:npixy),npixx,npixy,dens,err)
 print*,'average error in non-normalised xsec interpolation = ',err
 print "(70('-'))"

 
 call pgenv(xmin,xmax,ymin,ymax,0,0)
! call pgpixl(datpix,ipixx,ipixy,1,npixx,1,npixy,xmin,xmax,ymin,ymax)
! trans = 0.
! trans(1) = xmin - 0.5*dxpix
! trans(2) = dxpix
! trans(4) = ymin - 0.5*dxpix
! trans(6) = dxpix
 call pgimag(datpix,ipixx,ipixy,1,npixx,1,npixy,0.0,1.0,trans)
! call pgend
!
!--take normalised cross section at midplane and check density
!
 call interpolate3D_fastxsec(x(1:npart),y(1:npart),z(1:npart), &
      h(1:npart),weight(1:npart),dat(1:npart),itype(1:npart),npart,&
      xmin,ymin,0.0,datpix(1:npixx,1:npixy),npixx,npixy,dxpix,.true.)     

 call geterr(datpix(1:npixx,1:npixy),npixx,npixy,dens,err)
 print*,'average error in normalised xsec interpolation = ',err
 call pgenv(xmin,xmax,ymin,ymax,0,0)
 call pgimag(datpix,ipixx,ipixy,1,npixx,1,npixy,0.0,1.0,trans)

 print*,'closing PGPLOT'
 call pgend

 print "(70('-'))"

 print*,'SPEED CHECKS...'
 
 normalise = .true.
 npixx = 1
 npixy = 1
 npartx = idimx
 nparty = idimx
 npartz = idimx
 npart = npartx*nparty*npartz
 call setgrid(npartx,nparty,npartz,x,y,z,pmass,rho,h,weight,xmin,xmax,ymin,ymax,zmin,zmax)
 
 dxpix = (xmax-xmin)/real(npixx)
 call interpolate3D_projection(x(1:npart),y(1:npart),z(1:npart),h(1:npart), &
                               weight(1:npart),dat(1:npart),itype(1:npart),npart,xmin,ymin, &
                               datpix(1:npixx,1:npixy),npixx,npixy,dxpix,dxpix,normalise,0.,0.,ifastrender)
 
 call geterr(datpix(1:npixx,1:npixy),npixx,npixy,columndens,err)
 print*,'average error in projection = ',err

 npixx = 1000
 npixy = 1000
 npartx = 2
 nparty = 2
 npartz = 2
 npart = npartx*nparty*npartz
 call setgrid(npartx,nparty,npartz,x,y,z,pmass,rho,h,weight,xmin,xmax,ymin,ymax,zmin,zmax)
 
 dxpix = (xmax-xmin)/real(npixx)
 call interpolate3D_projection(x(1:npart),y(1:npart),z(1:npart),h(1:npart), &
                               weight(1:npart),dat(1:npart),itype(1:npart),npart,xmin,ymin, &
                               datpix(1:npixx,1:npixy),npixx,npixy,dxpix,dxpix,normalise,0.,0.,ifastrender)
 
 call geterr(datpix(1:npixx,1:npixy),npixx,npixy,columndens,err)
 print*,'average error in projection = ',err


 
contains

subroutine setgrid(npartx,nparty,npartz,x,y,z,pmass,rho,h,weight,xmin,xmax,ymin,ymax,zmin,zmax)
 implicit none
 integer , intent(in) :: npartx,nparty,npartz 
 real, dimension(:), intent(out) :: x,y,z,pmass,rho,h,weight
 real, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
 integer :: ipart,k,j,i
 real :: dx,dy,dz,ypos,zpos
 real :: totmass,massp,vol,dens,h0
 
 dz = (zmax-zmin)/real(npartz - 1)
 dy = (ymax-ymin)/real(nparty - 1)
 dx = (xmax-xmin)/real(npartx - 1)
 ipart = 0
 
 do k=1,npartz
    zpos = zmin + (k-1)*dz
    do j=1,nparty
       ypos = ymin + (j-1)*dy
       do i=1,npartx
          ipart = ipart + 1
          x(ipart) = xmin + (i-1)*dx
          y(ipart) = ypos
          z(ipart) = zpos
!          print*,ipart,'x,y,z=',x(ipart),y(ipart),z(ipart)
       enddo
    enddo
 enddo
 
 npart = npartx*nparty*npartz
!
!--set other properties
!
 totmass = 3.1415926536
 massp = totmass/real(npart)
 vol = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
 dens = totmass/vol
 h0 = 1.5*(massp/dens)**(1./3)
 print*,' testing ',npart,' particles in a cube configuration'
 print*,' dx = ',dx,' dy = ',dy,' dz = ',dz
 print*,' mass = ',massp,' dens = ',dens,' h = ',h0
 print*,' approx density = ',massp/(dx*dy*dz)
 do i = 1,npart
    pmass(i) = massp
    rho(i) = dens
    h(i) = h0
    dat(i) = rho(i)
    weight(i) = pmass(i)/(rho(i)*h(i)**3)
 enddo
  
end subroutine setgrid

subroutine geterr(datpix,npixx,npixy,datexact,err)
 implicit none
 integer, intent(in) :: npixx,npixy
 real, dimension(:,:), intent(in) :: datpix
 real, intent(in) :: datexact
 real, intent(out) :: err
 integer :: icalc,j,i
 real :: erri
 
 err = 0.
 icalc = 0
 do j=2,npixy-1
    do i=2,npixx-1
       icalc = icalc + 1
       erri = abs(datpix(i,j)-datexact)/datexact
       err = err + erri
       !if (erri.gt.0.05) print*,i,j,' xsec dens = ',datpix(i,j),' should be ',dens
    enddo
 enddo
 if (icalc.le.0) then
    print*,'cannot calculate error => npix too small'
    err = -1.0
 else
    err = err/real(icalc)
 endif
 
end subroutine geterr
 
end program test_interpolation
