!
! This module contains subroutines and variables for setting the 
! colour schemes for rendered plots
!
module colours
 implicit none
 integer, parameter :: ncolourmax = 256
 integer, parameter :: ncolourschemes = 14
 character(len=17), dimension(ncolourschemes), parameter :: schemename = &
    (/'greyscale        ', &
      'red              ', &
      'ice blue         ', &
      'heat             ', &
      'rainbow          ', &
      'prism            ', &
      'red-blue-yellow  ', &
      'blue-yellow-red  ', &
      'purple-blue-green', &
      'gamma            ', &
      'red-green-blue   ', &
      'blue-green-red   ', &
      'rainbow II       ', &
      'rainbow III      '/)
!
!--rgb colours of the colour table are stored in the array below
!  this is used for colour blending (opacity rendering)
!
  integer :: ifirstcolour, ncolours
  real, dimension(3,ncolourmax) :: rgbtable
  
contains

! ------------------------------------------------------------------------
!      defines colour schemes for rendering
!      ** add your own here **
! ------------------------------------------------------------------------     
subroutine colour_set(icolourscheme)
  implicit none
  integer, intent(in) :: icolourscheme
  integer :: i,icolmin,icolmax,ncolmax,nset,index
  real :: hue,light,sat,brightness,red,green,blue
  real :: dhue,dlight,dsat,dred,dgreen,dblue,contrast
  real, dimension(18) :: lumarr,redarr,greenarr,bluearr
  logical :: rgb
      
  red = 0.0
  green = 0.0
  blue = 0.0
  rgb = .false.
  ncolours = ncolourmax
  nset = 5
!
!--set first colour index (warning: colours 1-16 have presets, so
!  overwriting these means that line graphs that use colour will come
!  out funny). Best to leave 0 and 1 alone as these are black and white.
!
  ifirstcolour = 2
!
!--inquire as to colour range available on current device
!  adjust ncolours if necessary
!      
  call PGQCOL(icolmin,icolmax)
!  print*,' from device = ',icolmin,icolmax
  call PGQCIR(icolmin,icolmax)
!  print*,' other = ',icolmin,icolmax
  if (ifirstcolour.lt.icolmin) ifirstcolour = icolmin
  ncolmax = icolmax - ifirstcolour
  if (ncolours.gt.ncolmax) then  
     ncolours = ncolmax
     print*,'Warning: Device allows only ',ncolours+1,' colours'
  endif
  !
  !--set this as the range of colour indices to use
  !
  call PGSCIR(ifirstcolour,ifirstcolour+ncolours)  

!
!--starting values
!      
  select case(icolourscheme)
!  case(2)  
!    red
!     rgb = .false.
!     hue=100
!     light=0.0
!     sat=1.0      
!     dlight = 1.0/FLOAT(ncolours)
!     dsat = 0.0
!     dhue = 80.0/FLOAT(ncolours)         
  case(2)        ! green
!    universe
     rgb = .true.         
     red=0.
     green=0.
     blue=1./6.
     dred =1.0/FLOAT(ncolours)
     dgreen = 1.0/FLOAT(ncolours)
     dblue = (1./6.)/FLOAT(ncolours)
!  case(3)        ! blue
!    ice blue
!     rgb = .false.
!     hue=330
!     light=0.0
!     sat=1.0
!     dlight = 1.0/FLOAT(ncolours)
!     dsat = 0.0
!     dhue = 40.0/FLOAT(ncolours)
  case(3) 
!    rainbow
     rgb = .false.         
     hue=100
     light=0.5
     sat=1.0
     dlight = 0.0        !/FLOAT(ncolours)
     dsat = 0.0
     dhue = 320.0/FLOAT(ncolours)

!  case default
!    default is greyscale
!!     return
  end select
!
!--set colour indexes from ifirstcolour-> ifirstcolour+ncolours
!  increment values in steps
!
  if (abs(icolourscheme) .lt. 2) then
     do i=ifirstcolour,ifirstcolour+ncolours
        if (rgb) then
           red = red + dred
           blue = blue + dblue
           green = green + dgreen
           red = min(red,1.0)
           blue = min(blue,1.0)
           green = min(green,1.0)
           red = max(red,0.)
           blue = max(blue,0.)
           green = max(green,0.)
           call PGSCR(i,red,green,blue)
        else   
           hue = hue + dhue
           sat = sat + dsat
           light = light + dlight
           call PGSHLS(i,hue,light,sat)
        endif
     enddo
  elseif (abs(icolourscheme).lt.15) then
     brightness = 0.5
     contrast = 1.0
     !--invert colour table for negative values
     if (icolourscheme.lt.0) contrast = -1.0
     
     select case(abs(icolourscheme))
     case(2)
     !--red temperature (IDL red-temperature)
     nset = 5
     lumarr(1:nset) =  (/0.0,0.69,0.75,1.0/)
     redarr(1:nset) =  (/0.0,1.00,1.00,1.0/)
     greenarr(1:nset)= (/0.0,0.41,0.52,1.0/)
     bluearr(1:nset) = (/0.0,0.00,0.00,1.0/)
     case(3)
     !--ice blue (IDL blue-white)
     nset =  5
     lumarr(1:nset)  = (/0.000,0.376,0.737,0.753,1.000/)
     redarr(1:nset)  = (/0.000,0.000,0.000,0.000,1.000/)
     greenarr(1:nset)= (/0.000,0.000,0.580,0.604,1.000/)
     bluearr(1:nset) = (/0.000,0.510,1.000,1.000,1.000/)
     case(4)
     !--heat
     nset = 5
     lumarr(1:nset) =  (/0.,0.25,0.5,0.75,1.0/)
     redarr(1:nset) =  (/0.0,0.0,0.0,1.0,1.0/)
     greenarr(1:nset)= (/0.0,1.0,1.0,1.0,0.0/)
     bluearr(1:nset) = (/1.0,1.0,0.0,0.0,0.0/)
     case(5)
     !--rainbow
     nset = 8
     lumarr(1:nset) =  (/0.0,0.125,0.225,0.25,0.425,0.625,0.8125,1.0/)
     redarr(1:nset) =  (/0.0,0.341,0.100,0.00,0.000,0.000,1.0000,1.0/)
     greenarr(1:nset)= (/0.0,0.000,0.000,0.00,1.000,1.000,1.0000,0.0/)
     bluearr(1:nset) = (/0.0,0.569,1.000,1.00,1.000,0.000,0.0000,0.0/)
     case(6)
     !--prism (IDL prism)
     nset =  8
     lumarr(1:nset)  = (/0.000,0.251,0.263,0.494,0.502,0.749,0.753,1.000/)
     redarr(1:nset)  = (/0.000,0.953,1.000,0.035,0.000,0.000,0.000,0.000/)
     greenarr(1:nset)= (/0.000,0.000,0.043,0.969,1.000,0.000,0.000,0.000/)
     bluearr(1:nset) = (/0.000,0.000,0.000,0.000,0.027,0.984,1.000,0.000/)
!     case(6)
!     !--universe (this is my attempt at ripping off a spiral galaxy colour table)
!     nset = 6
!     lumarr(1:nset) =  (/0.0,0.20,0.4,0.60,0.95,1.0/)
!     redarr(1:nset) =  (/0.0,0.10,0.2,0.40,1.00,1.0/)
!     greenarr(1:nset)= (/0.0,0.15,0.2,0.42,1.00,1.0/)
!     bluearr(1:nset) = (/0.0,0.30,0.4,0.50,0.95,1.0/)  
     case(7)
     !--red-blue-yellow (IDL stern special)
     nset =  7
     lumarr(1:nset)  = (/0.000,0.055,0.247,0.251,0.502,0.737,1.000/)
     redarr(1:nset)  = (/0.000,0.996,0.000,0.251,0.502,0.737,1.000/)
     greenarr(1:nset)= (/0.000,0.055,0.247,0.251,0.502,0.737,1.000/)
     bluearr(1:nset) = (/0.000,0.106,0.490,0.498,1.000,0.000,1.000/)
!    this was another attempt at the same thing
!     nset = 7
!     lumarr(1:nset) =   (/0.0,0.2,0.35,0.4, 0.7, 0.9,1.0/)
!     redarr(1:nset) =   (/0.0,0.9,0.01,0.25,0.5,0.75,1.0/)
!     greenarr(1:nset) = (/0.0,0.1,0.24,0.25,0.5,0.75,1.0/)
!     bluearr(1:nset) =  (/0.0,0.2,0.48,0.50,1.0,0.05,1.0/)
     case(8)
     !--blue-yellow-red (IDL blue-red 0.34)
     nset = 10
     lumarr(1:nset)  = (/0.000,0.004,0.125,0.129,0.380,0.384,0.635,0.886,0.996,1.000/)
     redarr(1:nset)  = (/0.000,0.000,0.000,0.000,0.000,0.000,1.000,1.000,0.514,0.514/)
     greenarr(1:nset)= (/0.000,0.000,0.000,0.000,1.000,1.000,1.000,0.000,0.000,0.000/)
     bluearr(1:nset) = (/0.514,0.514,1.000,1.000,1.000,1.000,0.000,0.000,0.000,0.000/)
!     nset = 6
!     lumarr(1:nset) =  (/0.0,0.2,0.4,0.6,0.8,1.0/)
!     redarr(1:nset) =  (/0.0,0.0,0.5,1.0,1.0,0.5/)
!     bluearr(1:nset) = (/1.0,1.0,0.5,0.0,0.0,0.0/)
!     greenarr(1:nset)= (/0.0,1.0,0.5,1.0,0.0,0.0/)
     case(9)
     !--purple-blue-green
     nset = 6
     lumarr(1:nset) =  (/0.0,0.1,0.2,0.5,0.8,1.0/)
     redarr(1:nset) =  (/0.0,0.1,0.5,0.02,0.0,0.0/)
     bluearr(1:nset) = (/0.0,0.2,0.5,0.98,0.0,0.0/)
     greenarr(1:nset)= (/0.0,0.0,0.0,0.0,0.62,0.98/)
     case(10)
     !--gamma (IDL stdgamma-ii)
     nset = 18
     lumarr(1:nset)  =(/0.,0.184,0.192,0.251,0.31,0.376,0.427,0.431,0.443,0.502,0.569,0.624,0.635,0.682,0.69,0.749,0.753,1./)
     redarr(1:nset)  =(/0.,0.000,0.035,0.318,0.31,0.643,0.914,1.000,1.000,1.000,1.000,1.000,1.000,0.639,0.678,0.976,1.00,1./)
     greenarr(1:nset)=(/0.,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.000,0.318,0.639,0.639,0.639,0.639,0.639,1.000,1.00,1./)
     bluearr(1:nset) =(/0.,0.957,1.000,0.682,0.365,0.00,0.000,0.000,0.000,0.000,0.322,0.000,0.000,0.000,0.000,0.188,0.20,1./)
     case(11)
     nset = 3
     !--red-greeny-blue
     lumarr(1:nset) =   (/0.0,0.5,1.0/)
     redarr(1:nset) =   (/1.0,0.66,0.0/)
     greenarr(1:nset) = (/0.0,0.66,0.0/)
     bluearr(1:nset) =  (/0.0,0.66,1.0/)
     case(12)
     nset = 3
     !--blue-green-red ("highlight")
     lumarr(1:nset) =   (/0.0,0.5,1.0/)
     redarr(1:nset) =   (/0.0,0.5,1.0/)
     greenarr(1:nset) = (/0.0,1.0,0.0/)
     bluearr(1:nset) =  (/1.0,0.5,0.0/)
!     nset = 5
     !--dolag other
!     lumarr(1:nset) =   (/0.0,0.33,0.5,0.66,1.0/)
!     redarr(1:nset) =   (/0.0,1.00,0.5,0.00,1.0/)
!     greenarr(1:nset) = (/0.0,0.66,1.0,0.66,0.0/)
!     bluearr(1:nset) =  (/1.0,0.66,0.5,0.33,0.0/)
     case(13)
     !--rainbow II (as used in NS merger I)
     nset = 10
     lumarr(1:nset)  = (/0.000,0.153,0.157,0.310,0.314,0.467,0.471,0.624,0.627,1.000/)
     redarr(1:nset)  = (/1.000,1.000,0.996,0.016,0.000,0.000,0.000,0.000,0.020,1.000/)
     greenarr(1:nset)= (/0.000,0.980,1.000,1.000,1.000,1.000,0.984,0.004,0.000,0.000/)
     bluearr(1:nset) = (/0.000,0.000,0.000,0.000,0.012,0.988,1.000,1.000,1.000,1.000/)
     case(14)
     !--rainbow III
     nset = 13
     lumarr(1:nset)  = (/0.000,0.004,0.110,0.114,0.333,0.557,0.561,0.565,0.569,0.776,0.780,0.996,1.000/)
     redarr(1:nset)  = (/0.486,0.486,0.012,0.000,0.000,0.004,0.020,0.051,0.055,0.992,1.000,1.000,1.000/)
     greenarr(1:nset)= (/0.000,0.000,0.000,0.008,0.996,1.000,1.000,1.000,1.000,1.000,0.988,0.020,0.020/)
     bluearr(1:nset) = (/1.000,1.000,1.000,1.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/)
     end select

     call PGCTAB(lumarr(1:nset),redarr(1:nset),greenarr(1:nset),bluearr(1:nset), &
                 nset,contrast,brightness)
  endif
!
!--if icolourscheme = ncolourschemes+1 set the PGPLOT colour indices 
!  from the contents of the rgbtable array
!
  if (abs(icolourscheme).eq.ncolourschemes+1) then
     call PGSCIR(ifirstcolour,ifirstcolour+ncolourmax)
     do i=1,ncolourmax
        index = ifirstcolour + (i-1)
        call pgscr(index,rgbtable(1,i),rgbtable(2,i),rgbtable(3,i))
     enddo
     print "(1x,a)",'using colour scheme other'

  elseif (abs(icolourscheme).le.ncolourschemes) then
!
!--also store the colour table as a list of r,g,b values
!
     do i=1,ncolours+1
        index = ifirstcolour + (i-1)
        call pgqcr(index,rgbtable(1,i),rgbtable(2,i),rgbtable(3,i))
     enddo

     if (icolourscheme.lt.0) then
        print "(1x,a)",'using colour scheme inverse '//trim(schemename(abs(icolourscheme)))
     else
        print "(1x,a)",'using colour scheme '//trim(schemename(icolourscheme))
     endif
  else
     print "(1x,a)",'warning: unknown colour scheme - uses default greyscale'

  endif
  
  return
  
end subroutine colour_set

!------------------------------------------------
! demonstration plot of all the colour schemes
!------------------------------------------------
subroutine colour_demo
  implicit none
  integer :: i,j,nc
  !
  !--npixx should be >= ncolours in setcolours.f
  !      
  integer, parameter :: npixx = ncolourmax
  integer, parameter :: npixy = npixx/10
  real, dimension(npixx,npixy) :: sample
  real :: xmin,xmax,ymin,ymax,dx,dy,trans(6)
  character(len=10) :: string

  call pgbegin(0,'/xw',1,ncolourschemes)
!!  call pgpaper(6.0,8.0) !!!0.25/sqrt(2.))

  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 0.1
  dx = (xmax-xmin)/float(npixx)
  dy = (ymax-ymin)/float(npixy)
  trans(1) = xmin - 0.5*dx
  trans(2) = dx
  trans(3) = 0.0
  trans(4) = ymin - 0.5*dy
  trans(5) = 0.0
  trans(6) = dy

  do j=1,npixy
     do i=1,npixx
        sample(i,j) = (i-1)*dx
     enddo
  enddo
  call pgsch(2.0)
  call pgenv(xmin,xmax,ymin,ymax,0,-1)
  call pgsch(1.0)
  call pggray(sample,npixx,npixy,1,npixx,1,npixy, &
              minval(sample),maxval(sample),trans)
  call pgnumb(1,0,0,string,nc)
  call pgsch(7.0)
  call pgmtxt('t',0.5,0.5,0.5,string(1:nc)//': '//trim(schemename(1)))

  do i=2,ncolourschemes
     call pgsch(2.0)      
     call pgenv(xmin,xmax,ymin,ymax,0,-1) 
     call pgsch(7.0)
     call pgnumb(i,0,0,string,nc)
     call pgmtxt('t',0.5,0.5,0.5,string(1:nc)//': '//trim(schemename(i)))     
     call colour_set(i)
     call pgimag(sample,npixx,npixy,1,npixx,1,npixy, &
                 minval(sample),maxval(sample),trans)
  enddo

  call pgsch(1.0)
  call pgend 

end subroutine colour_demo


end module colours
