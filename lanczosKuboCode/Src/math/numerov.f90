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

!****** Subroutine: NumerovIntIn **********************************************
!******************************************************************************
!
!  Integral with the Numerov method inwards
! Reference:
!   [1] Kenhere, D. G.; PHYSICS THROUGH COMPUTATION – II, 'Bound State of One
!       Dimensional Potential by Numerov Method'
!       (2006).
!
! INPUT -----------------------------------------------------------------------
!
! real f(*)                       : Function
! integer k                       : Order of the polynomial
!
! OUTPUT ----------------------------------------------------------------------
!
! real df(*)                      : Differential of f
!
!******************************************************************************
function NumerovIn(f,n,i,h) result(y)

   implicit none

   include 'parameters.h'

   integer, intent(in) :: n, i
   real(dp), intent(in) :: f(n)
   real(dp), intent(in) :: h
   real(dp) :: y

   y = (f(i-1)+10.0_dp*f(i)+f(i+1))*(h**2)/12.0_dp

end function NumerovIn
!****** End subroutine: NumerovIntIn ******************************************
!******************************************************************************

