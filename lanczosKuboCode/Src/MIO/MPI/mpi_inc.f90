!---------------------------------------------------------!
!             MEMORY, INPUT & OUTPUT LIBRARY              !
!                                                         !
!     Copyright (C) 2012 Rafael Martinez-Gordillo         !
!                                                         !
!   This file is distributed under the terms of the       !
!   GNU General Public License. See the file `LICENSE'    !
!   in the root directory of the present distribution,    !
!   or http://www.gnu.org/copyleft/gpl.txt                !
!                                                         !
!---------------------------------------------------------!

module mpi_inc

   implicit none

   include 'mpif.h'

   integer, parameter :: defPNode = 0 ! Default principal node

   integer, public, save :: nProc=1, Node=0 ! Number of processes and rank of node
   logical, public, save :: IsPNode = .true. ! True if it's the principal node
   integer, public, save :: rankPNode = defPNode  ! Rank of principal node

end module mpi_inc
