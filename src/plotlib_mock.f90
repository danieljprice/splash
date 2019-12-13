!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2014 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

!---------------------------------------------------------------------------
!  The mocking plotlib module to fake a real plotting module behaviour
!  in external libraries (libexact)
!---------------------------------------------------------------------------
module plotlib
 implicit none
 public

contains

subroutine plot_pt1(R1,R2,I)
 real :: R1, R2
 integer :: I
end subroutine

subroutine plot_swin(R1,R2,R3,R4)
 real :: R1, R2, R3, R4
end subroutine

subroutine plot_box(S1,R1,I1,S2,R2,I2)
 real :: R1, R2
 integer :: I1, I2
 character(len=*) :: S1, S2
end subroutine

subroutine plot_funx(F1,I1,R1,R2,I2)
 real :: R1, R2, F1
 integer :: I1, I2
 external F1
end subroutine

subroutine plot_label (S1,S2,S3)
 character(len=*) :: S1, S2, S3
end subroutine

subroutine plot_line(I,R1,R2)
 real :: R1(:), R2(:)
 integer :: I
end subroutine

end module plotlib
