 !---------------------------------------
 ! numerically integrate a polytrope
 ! with the density = 1 at the centre
 ! This uses scaled variables (no sigma in the equation)
 ! then the radius is scaled to unity and the central
 ! density is calculated to get unit mass.
 ! Finally the polytropic K is calculated
 !---------------------------------------

       subroutine exact_polytrope(gamma)
       implicit none
       integer :: i,j
       real, parameter :: pi = 3.1415926536
       real, intent(in) :: gamma
       real, dimension(1000) :: v,r,den,rplot,denplot
       real :: dr,an,rhs,radv,sigma,totmassf,totmass,akf
              
       dr=0.005
       an = 1./(gam-1.)
       v(1) = 0.0
       v(2) = dr*(1.0 - dr*dr/6. )
       
       i = 2
       r(1) = 0.
       r(2) = dr
      
       do while (v(i).ge.0.)         
         r(i) = (i - 1)*dr         
         rhs = - r(i)*(v(i)/r(i))**an
         v(i+1) = 2*v(i) - v(i-1) + dr*dr*rhs         
         i = i + 1                   
       enddo
       !----------------------------------
       ! now scale the radius
       !----------------------------------
       radv = r(i-1)
       print *,' rad v ',r(i-1),v(i-1)
       
       sigma = radv/2.	! divide by two as a kernel
       print *,' sigma ',sigma
       
       !--------------------------------------
       ! calculate the mass out to radius r
       ! using the density without the central
       ! density multiplier- call this totmassf
       ! the true scaled totmass = 1.
       !----------------------------------------
       
       den(1) = 1.0
       totmassf = 0.
       do j = 2,(i-1)
         den(j) = (v(j)/r(j))**an
         totmassf = totmassf + 4.*pi*r(j)*r(j)*den(j)*dr
       enddo
       
       !-------------------------------------
       ! calculate sigma to give unit radius
       ! calculate the central density
       !-------------------------------------
       
       rhocentre = sigma**3/totmassf
       totmass = rhocentre*totmassf/sigma**3
       print *,' rhocentre totmass ',rhocentre,totmass
       
       akf = 4.*pi*rhocentre**(1.-1./an)/((an+1)*sigma**2)
       print *,' akf ',akf
       
       print *,'   den   r'
       do j = 2,(i-1)
         if(r(j).gt.0.0)then
           realden = rhocentre*den(j)
	   realrad = r(j)/sigma
	   denplot(j-1) = realden
	   rplot(j-1) = realrad
         endif.
       
       !---------------------------------------
       ! plot results using PGPLOT       
       !---------------------------------------
       
       print*,' plotting results ',j-2,rplot(j-2),denplot(j-2)
       print*,rplot(1),denplot(1)
       CALL PGLINE(j-2,rplot(1:j-2),denplot(1:j-2))
              
       return
       end subroutine exact_polytrope
       
       
          
       
       
       
       
       
