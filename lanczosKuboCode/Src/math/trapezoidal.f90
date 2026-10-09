!---------------------------------------------------------!
!                 MATH LIBRARY                            !
!                                                         !
!     Copyright (C) 2012 Rafael Martinez-Gordillo         !
!                                                         !
!   This file is distributed under the terms of the       !
!   GNU General Public License. See the file `LICENSE'    !
!   in the root directory of the present distribution,    !
!   or http://www.gnu.org/copyleft/gpl.txt                !
!                                                         !
!---------------------------------------------------------!


!****** Function: TrapezoidalInt **********************************************
!******************************************************************************
!
!  Calculate the integral of a function
!  References:
!   [1] Johnson, W. R.; 'Atomic Structure Theory: Lectures on Atomic Physics',
!       Springer:Berlin, 2007; Chapter 2, pp. 35-45
!
! INPUT -----------------------------------------------------------------------
!
! real f(n)                       : Function
! integer n                       : Size of f
! real h                          : Step size
! integer k                       : Order of the integral
!
! RESUTL ----------------------------------------------------------------------
!
! real s                          : Integral of f
!
!******************************************************************************
function TrapezoidalInt(f,n,h,k) result(s)

   implicit none

   include 'parameters.h'

   integer, intent(in) :: n, k
   real(dp), intent(in) :: f(n), h
   real(dp) :: s

   integer :: i

   if (k==1 .or. k==2) then
      s = 0.5_dp*(f(1)+f(n))
      s = s+sum(f(2:n-1))
   elseif (k==3 .or. k==4) then
      s = (9.0_dp*(f(1)+f(n))+28.0_dp*(f(2)+f(n-1))+23.0_dp*(f(3)+f(n-2)))/24.0_dp
      s = s+sum(f(4:n-3))
   elseif (k==5 .or. k==6) then
      s = (47.5_dp*(f(1)+f(n))+190.2_dp*(f(2)+f(n-1))+110.4_dp*(f(3)+f(n-2))+ &
        158.6_dp*(f(4)+f(n-3))+141.3_dp*(f(5)+f(n-4)))/144.0_dp
      s = s+sum(f(6:n-5))
   elseif (k==7 .or. k==8) then
      s = (367.99_dp*(f(1)+f(n))+1766.48_dp*(f(2)+f(n-1))+548.51_dp*(f(3)+f(n-2))+ &
        1779.84_dp*(f(4)+f(n-3))+894.37_dp*(f(5)+f(n-4))+1309.36_dp*(f(6)+f(n-5))+ &
        1195.85_dp*(f(7)+f(n-6)))/1209.60_dp
      s = s+sum(f(8:n-7))
   end if
   s = s*h

end function TrapezoidalInt
!****** End function: TrapezoidalInt ******************************************
!******************************************************************************
