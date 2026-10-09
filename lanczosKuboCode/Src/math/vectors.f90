!---------------------------------------------------------!
!                  ADVENTURE                              !
!                                                         !
!     Copyright (C) 2012 Rafael Martinez-Gordillo         !
!                                                         !
!   This file is distributed under the terms of the       !
!   GNU General Public License. See the file `LICENSE'    !
!   in the root directory of the present distribution,    !
!   or http://www.gnu.org/copyleft/gpl.txt                !
!                                                         !
!---------------------------------------------------------!

!****** Function: CrossProd ***************************************************
!******************************************************************************
!
!  Calculate the cross product of two vectors
!
! INPUT -----------------------------------------------------------------------
!
! real a, b                       : Vectors
!
! RESUTL ----------------------------------------------------------------------
!
! real axb                        : Cross product of a and b
!
!******************************************************************************
function CrossProd(a,b) result(axb)

    implicit none

    include 'parameters.h'

    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: axb(3)

    axb(1) = a(2)*b(3) - a(3)*b(2)
    axb(2) = a(3)*b(1) - a(1)*b(3)
    axb(3) = a(1)*b(2) - a(2)*b(1)

end function CrossProd
!****** End function: CrossProd ***********************************************
!******************************************************************************

!****** Function: norm ********************************************************
!******************************************************************************
!
!  Calculate the norm of a vector
!
! INPUT -----------------------------------------------------------------------
!
! real a                          : Vector
!
! RESUTL ----------------------------------------------------------------------
!
! real norm                       : The norm of the vector
!
!******************************************************************************
function norm(a)

    implicit none

    include 'parameters.h'

    real(dp), intent(in) :: a(3)
    real(dp) :: norm

    norm = sqrt(dot_product(a,a))

end function norm
!****** End function: norm ****************************************************
!******************************************************************************
