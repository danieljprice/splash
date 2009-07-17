!
! subroutine to plot isosurface of 3D SPH data 
! by reflecting light off a surface
!
! Daniel Price, Institute of Astronomy, Cambridge
! June 2004
!
! inputs:
!    npart        : number of particles
!    x(3,npart)   : particle coordinates
!    hh(npart)    : smoothing length of particles
!    pmass(npart) : particle masses
!    rho(npart)   : density
!    dat(npart)   : data to take isosurface of
!    datlevel     : value of data at isosurface
!
subroutine isosurface(npart,x,hh,pmass,rho,dat,datlevel)
implicit none
 integer, intent(in) :: npart
 real, intent(in), dimension(3,npart) :: x
 real, intent(in), dimension(npart) :: hh,pmass,rho
 real, intent(in) :: datlevel
 integer :: i, iheadoflist, npartinsurface
 integer, dimension(npart) :: linklist
 real, dimension(3,npart) :: svec
 real, dimension(3) :: lightpos, dx
 real :: hi,hav,q2,wab,rij2
!
!--set position of light source
!
 lightpos(1) = 1.0  ! position in x
 lightpos(2) = 1.0  !             y
 lightpos(3) = 1.0  !             z
!
!--tag particles with dat < datlevel as being within the surface
!  store the particles in the surface using a linked list
!
 linklist(i) = 666
 iheadoflist = -1
 npartinsurface = 0

 do i=1,npart
    if (dat(i).le.datlevel) then
       linklist(i) = iheadoflist
       iheadoflist = i
       npartinsurface = npartinsurface + 1
    endif
 enddo

 if (npartinsurface.le.0) then
    print*,'isosurface: no particles in surface: doing nothing...'
    return
 else
    print*,'particles in surface = ',npartinsurface
 endif
!
!--initialise svec
!
 svec = 0.
!
!--loop over all the particles in the surface and calculate the surface vector (svec)
!
 ipart = iheadoflist

 do i=1,npartinsurface
    rho1i = 1./rho(ipart)
    hi = hh(ipart)
   
    jpart = iheadoflist
    do j=i,npartinsurface
       dx(:) = x(:,ipart) - x(:,jpart)    
       rij2 = DOT_PRODUCT(dx,dx)
       ! use average h
       hav = 0.5*(hi + hh(jpart))
       q2 = rij2/hav**2
       if (q2.lt.2.) then
          svec(:,ipart) = svec(:,ipart) + pmass(jpart)*dx(:)*wab/rho(jpart)
          svec(:,jpart) = svec(:,jpart) - pmass(ipart)*dx(:)*wab*rho1i
       enddo
       jpart = ll(jpart)
    enddo
    ipart = ll(ipart)
 enddo
!
!--now interpolate to an array of pixels by bouncing light off the particles
!
! bounce light off surface
 
end subroutine isosurface
