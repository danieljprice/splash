!
!--unit test for field line plotting routine
!
program test_fieldlines
 use fieldlines, only:streamlines
 use render, only:render_pix
 implicit none
 integer, parameter :: ipixx = 1000, ipixy = 1000, nc = 100
 integer :: npixx, npixy,i,j
 real, parameter :: errtol = 1.e-7
 real, dimension(ipixx,ipixy) :: datpix,vecpixx,vecpixy,datpix1,datpix2,vecpixx1,vecpixy1
 real :: xmin,xmax,ymin,ymax,zmin,zmax,dxpix
 real :: xi,yi,datmin,datmax,err,dcont
 real, parameter :: pi=3.1415926536
 real, dimension(6) :: trans
 real, dimension(nc) :: levels
 
 xmin = -0.5
 xmax = 0.5
 ymin = -0.5
 ymax = 0.5

 call pgopen('/xw')
 
 print "(70('-'))"
 print*,'ORSZAG-TANG TEST'

 npixx = 400
 npixy = 400
 dxpix = (xmax-xmin)/real(npixx)
 do j = 1,npixy
    yi = ymin + (j-0.5)*dxpix
    if (j.le.10) print*,'y = ',yi
    do i = 1,npixx
       xi = xmin + (i-0.5)*dxpix
       vecpixx(i,j) = func_vecx(xi,yi)
       vecpixy(i,j) = func_vecy(xi,yi)
       datpix1(i,j) = func_stream(xi,yi)
    enddo
 enddo

 call streamlines(vecpixx(1:npixx,1:npixy),vecpixy(1:npixx,1:npixy), &
      datpix(1:npixx,1:npixy),npixx,npixy,dxpix)

 do j = 1,npixy,2
    yi = ymin + (j-0.5)*2.*dxpix
    do i = 1,npixx,2
       xi = xmin + (i-0.5)*2.*dxpix
       vecpixx1((i+1)/2,(j+1)/2) = func_vecx(xi,yi) !vecpixx(i,j)
       vecpixy1((i+1)/2,(j+1)/2) = func_vecy(xi,yi) !vecpixy(i,j)
    enddo
 enddo
 npixx = npixx/2
 npixy = npixy/2
 dxpix = (xmax-xmin)/real(npixx)
 call streamlines(vecpixx1(1:npixx,1:npixy),vecpixy1(1:npixx,1:npixy), &
      datpix2(1:npixx,1:npixy),npixx,npixy,dxpix)
 npixx = npixx*2
 npixy = npixy*2 
 dxpix = (xmax-xmin)/real(npixx)
 !do j=1,npixy,2
 !   do i=1,npixx,2
 !      datpix(i,j) = (4.*datpix(i,j) - datpix2((i+1)/2,(j+1)/2))/3.
 !   enddo
 !enddo

 call pgenv(xmin,xmax,ymin,ymax,1,0)
 trans = 0.
 trans(1) = xmin - 0.5*dxpix
 trans(2) = dxpix
 trans(4) = ymin - 0.5*dxpix
 trans(6) = dxpix
 datmax = maxval(datpix(1:npixx,1:npixy))
 datmin = minval(datpix(1:npixx,1:npixy))
 print*,'min,max datpix = ',datmin,datmax
 datpix(1:npixx,1:npixy) = datpix(1:npixx,1:npixy) - 0.5*(datmax + datmin)
 datmax = maxval(datpix(1:npixx,1:npixy))
 datmin = minval(datpix(1:npixx,1:npixy))
 print*,'min,max datpix = ',datmin,datmax

 print*,'plotting integrated array...'
! call render_pix(datpix(1:npixx,1:npixy),datmin,datmax,'crap',npixx,npixy, &
!                 xmin,ymin,dxpix,0,.true.,.false.,30,.false.)
 call pgimag(datpix,ipixx,ipixy,1,npixx,1,npixy,datmin,datmax,trans)    
!
!--set contour levels
! 
  dcont = (datmax-datmin)/real(nc+1)   ! even contour levels
  do i=1,nc
     levels(i) = datmin + real(i)*dcont
  enddo
!
!--plot contours (use pgcont if pgcons causes trouble)
!
 call pgcons(datpix(1:npixx,1:npixy),npixx,npixy,1,npixx,1,npixy,levels(1:nc),nc,trans)

! call pgenv(xmin,xmax,ymin,ymax,1,0)

 datmax = maxval(datpix1(1:npixx,1:npixy))
 datmin = minval(datpix1(1:npixx,1:npixy)) 
 print*,'min,max datpix1 = ',datmin,datmax

 !print*,' plotting exact solution ...'
! call pgimag(datpix1,ipixx,ipixy,1,npixx,1,npixy,datmin,datmax,trans)
 !call pgcons(datpix1(1:npixx,1:npixy),npixx,npixy,1,npixx,1,npixy,levels(1:nc),nc,trans)

 call geterr(datpix(1:npixx,1:npixy),npixx,npixy,datpix1(1:npixx,1:npixy),err)
 print*,'average error in stream line calculation  = ',err
 if (npixx.eq.400 .and. npixy.eq.400) then
    if (err < 0.00019) then
       print*,'PASSED: error within limits'
    else
       print*,'FAILED: error too large!'    
    endif
 else
    print*,'setup different to usual one'
 endif
 
!---------------------
!  dipole field test
!----------------------

 print "(70('-'))"
 print*,'DIPOLE TEST'
 npixx = 400
 npixy = 400
 dxpix = (xmax-xmin)/real(npixx)
 do j = 1,npixy
    yi = ymin + (j-0.5)*dxpix
    do i = 1,npixx
       xi = xmin + (i-0.5)*dxpix
       vecpixx(i,j) = func_vecx_dipole(xi,yi)
       vecpixy(i,j) = func_vecy_dipole(xi,yi)
       datpix1(i,j) = func_stream_dipole(xi,yi)
    enddo
 enddo

 call streamlines(vecpixx(1:npixx,1:npixy),vecpixy(1:npixx,1:npixy), &
      datpix(1:npixx,1:npixy),npixx,npixy,dxpix)

 call pgenv(xmin,xmax,ymin,ymax,1,0)
 trans = 0.
 trans(1) = xmin - 0.5*dxpix
 trans(2) = dxpix
 trans(4) = ymin - 0.5*dxpix
 trans(6) = dxpix
 datmax = maxval(datpix(1:npixx,1:npixy))
 datmin = minval(datpix(1:npixx,1:npixy))
 print*,'min,max datpix = ',datmin,datmax
 datpix(1:npixx,1:npixy) = datpix(1:npixx,1:npixy) - 0.5*(datmax + datmin)
 datmax = maxval(datpix(1:npixx,1:npixy))
 datmin = minval(datpix(1:npixx,1:npixy))
 print*,'min,max datpix = ',datmin,datmax

! call pgimag(datpix,ipixx,ipixy,1,npixx,1,npixy,datmin,datmax,trans)    
!
!--set contour levels
! 
  dcont = (datmax-datmin)/real(nc+1)   ! even contour levels
  do i=1,nc
     levels(i) = datmin + real(i)*dcont
  enddo
!
!--plot contours (use pgcont if pgcons causes trouble)
!
 call pgcons(datpix(1:npixx,1:npixy),npixx,npixy,1,npixx,1,npixy,levels(1:nc),nc,trans)

! call pgenv(xmin,xmax,ymin,ymax,1,0)

 datmax = maxval(datpix1(1:npixx,1:npixy))
 datmin = minval(datpix1(1:npixx,1:npixy)) 
 print*,'min,max datpix1 = ',datmin,datmax

! call pgimag(datpix1,ipixx,ipixy,1,npixx,1,npixy,datmin,datmax,trans)
! call pgcons(datpix1(1:npixx,1:npixy),npixx,npixy,1,npixx,1,npixy,levels(1:nc),nc,trans)

 call geterr(datpix(1:npixx,1:npixy),npixx,npixy,datpix1(1:npixx,1:npixy),err)
 print*,'average error in stream line calculation  = ',err
 if (npixx.eq.400 .and. npixy.eq.400) then
    if (err < 0.00019) then
       print*,'PASSED: error within limits'
    else
       print*,'FAILED: error too large!'    
    endif
 else
    print*,'setup different to usual one'
 endif
 
 
 call pgend

 print "(70('-'))"
 
contains

real function func_vecx(xi,yi)
 implicit none
 real :: xi,yi
 
 func_vecx = -sin(2.*pi*yi)
 
end function func_vecx

real function func_vecy(xi,yi)
 implicit none
 real :: xi,yi

 func_vecy = sin(4.*pi*xi)
 
end function func_vecy

real function func_stream(xi,yi)
 implicit none
 real :: xi,yi
 
 func_stream = 0.5/pi*(cos(2.*pi*yi) + 0.5*cos(4.*pi*xi))
 
end function func_stream

!----------------
! dipole 
!----------------
real function func_vecx_dipole(xi,yi)
 implicit none
 real :: xi,yi
 real, parameter :: Bdipole = 1.0, Rdipole = 0.3, eps=0.3

 func_vecx_dipole = 3.*Bdipole*Rdipole**3*(Rdipole + xi)*yi/ &
        sqrt((Rdipole + xi)**2 + yi**2 + (eps*Rdipole)**2)**5
 
end function func_vecx_dipole

real function func_vecy_dipole(xi,yi)
 implicit none
 real :: xi,yi
 real, parameter :: Bdipole = 1.0, Rdipole = 0.3, eps=0.3
 
 func_vecy_dipole = Bdipole*Rdipole**3/ &
        sqrt((Rdipole + xi)**2 + yi**2 + (eps*Rdipole)**2)**3
 
end function func_vecy_dipole

real function func_stream_dipole(xi,yi)
 implicit none
 real :: xi,yi
 real, parameter :: Bdipole = 1.0, Rdipole = 0.3, eps=0.3
 
 func_stream_dipole = -Bdipole*Rdipole**3*(Rdipole + xi)/ &
        sqrt((Rdipole + xi)**2 + yi**2 + (eps*Rdipole)**2)**3
  
end function func_stream_dipole


subroutine geterr(datpix,npixx,npixy,datexact,err)
 implicit none
 integer, intent(in) :: npixx,npixy
 real, dimension(:,:), intent(in) :: datpix
 real, dimension(:,:), intent(in) :: datexact
 real, intent(out) :: err
 integer :: icalc,j,i
 real :: errij
 
 err = 0.
 icalc = 0
 do j=1,npixy
    do i=1,npixx
       icalc = icalc + 1
       errij = abs(datpix(i,j)-datexact(i,j))
       err = err + errij
    enddo
 enddo
 if (icalc.le.0) then
    print*,'cannot calculate error => npix too small'
    err = -1.0
 else
    err = err/(icalc*maxval(datexact))
 endif
 
end subroutine geterr
 
end program test_fieldlines
