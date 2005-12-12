!
!--unit test for interpolation routines
!
program test_interpolation
 use projections3D
 use xsections3D
 implicit none
 integer, parameter :: npartx = 50
 integer, parameter :: nparty = npartx, npartz = npartx
 integer, parameter :: npart = npartx*nparty*npartz
 integer, parameter :: npixx = 10, npixy = 10
 integer :: i,j,k,ipart,icalc
 real, dimension(npart) :: x,y,z,pmass,h,rho
 real, dimension(npart) :: dat
 real, dimension(npixx,npixy) :: datpix
 real :: xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz,ypos,zpos
 real :: totmass,massp,vol,dens,h0,columndens,dxpix
 real :: trans(6),err,erri
!
!--set up a cubic lattice of particles
!
 xmin = -0.5
 xmax = 0.5
 ymin = -0.5
 ymax = 0.5
 zmin = -0.5
 zmax = 0.5
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
! call pgopen('/xw')
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
!--set other properties
!
 totmass = 1.0
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
 enddo
!
!--setup integrated kernel table
!
 call setup_integratedkernel 
!
!--check value of the integration at q=zero (can do this analytically)
!
 print*,'coltable(0) = ',coltable(1),' should be ',2.*0.75/3.1415926536
!
!--now call interpolation routine to 10x10x10 pixels
!
 dxpix = (xmax-xmin)/real(npixx)
 call interpolate3D_projection(x,y,z,pmass,rho,h,dat,npart,xmin,ymin, &
                               datpix,npixx,npixy,dxpix,0.,0.)
!
!--check output
!

 columndens = dens*(zmax-zmin)

 err = 0.
 icalc = 0
 do j=2,npixy-1
    do i=2,npixx-1
       icalc = icalc + 1
       erri = abs(datpix(i,j)-columndens)/columndens
       err = err + erri
       !if (erri.gt.0.05) print*,i,j,' column dens = ',datpix(i,j),' should be ',columndens
    enddo
 enddo
 print*,'average error in column density interpolation = ',err/real(icalc)

! call pgenv(xmin,xmax,ymin,ymax,0,0)
! call pgpixl(datpix,npixx,npixy,1,npixx,1,npixy,xmin,xmax,ymin,ymax)
! trans = 0.
! trans(1) = xmin - 0.5*dxpix
! trans(2) = dxpix
! trans(4) = ymin - 0.5*dxpix
! trans(6) = dxpix
! call pgimag(datpix,npixx,npixy,1,npixx,1,npixy,0.0,1.0,trans)

!
!--take cross section at midplane and check density
!
 call interpolate3D_fastxsec(x,y,z,pmass,rho,h,dat,npart,&
     xmin,ymin,0.0,datpix,npixx,npixy,dxpix)
     
 print*,'cross section at z=0 :'
 err = 0.
 icalc = 0
 do j=2,npixy-1
    do i=2,npixx-1
       icalc = icalc + 1
       erri = abs(datpix(i,j)-dens)/dens
       err = err + erri
       !if (erri.gt.0.05) print*,i,j,' xsec dens = ',datpix(i,j),' should be ',dens
    enddo
 enddo
 print*,'average error in xsec interpolation = ',err/real(icalc)

 
! call pgenv(xmin,xmax,ymin,ymax,0,0)
! call pgpixl(datpix,npixx,npixy,1,npixx,1,npixy,xmin,xmax,ymin,ymax)
! trans = 0.
! trans(1) = xmin - 0.5*dxpix
! trans(2) = dxpix
! trans(4) = ymin - 0.5*dxpix
! trans(6) = dxpix
! call pgimag(datpix,npixx,npixy,1,npixx,1,npixy,0.0,1.0,trans)
! call pgend

end program test_interpolation
