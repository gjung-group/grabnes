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

module qsort

   implicit none

   PRIVATE

   public :: QuickSort

contains

!****** Subroutine: QuickSort *************************************************
!******************************************************************************
!
!  Sort an array by numerical ascending order
! Reference:
!   [1] http://www.fortran.com/qsort_c.f95
!
! INPUT -----------------------------------------------------------------------
!
! real array                      : Array to be sorted
!
!******************************************************************************
recursive subroutine QuickSort(array)

   implicit none

   include 'parameters.h'

   real(dp), intent(inout) :: array(:)

end subroutine QuickSort
!****** End subroutine: QuickSort *********************************************
!******************************************************************************

end module qsort
