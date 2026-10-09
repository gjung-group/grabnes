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

module mpi

   use mpitime
   use mpi_inc

   use mpi_c
   use mpi_z
   use mpi_r
   use mpi_d
   use mpi_i
   use mpi_l
   use mpi_s

   implicit none

   integer, private :: ierr

   public :: MPIIni
   public :: MPIFin
   public :: MPIGetNodes
   public :: MPIIsPNode
   public :: MPIBarrier

contains

!****** Subroutine: MPIIni ****************************************************
!******************************************************************************
!
!  Initialize MPI
!
!******************************************************************************
subroutine MPIIni()

   call MPITimer(1)
   call MPI_Init(ierr)
   call MPI_Comm_Size(MPI_COMM_WORLD,nProc,ierr)
   call MPI_Comm_Rank(MPI_COMM_WORLD,Node,ierr)
   if (Node /= rankPNode) IsPNode=.false.
   call MPITimer(0)

end subroutine MPIIni
!****** End subroutine: MPIIni ************************************************
!******************************************************************************


!****** Subroutine: MPIFin ****************************************************
!******************************************************************************
!
!  Terminate MPI
!
!******************************************************************************
subroutine MPIFin()

   call MPI_Finalize(ierr)

end subroutine MPIFin
!****** End subroutine: MPIFin ************************************************
!******************************************************************************


!****** Subroutine: MPIGetNodes ***********************************************
!******************************************************************************
!
!  Return the rank of the process and the total number of processes
!
!******************************************************************************
subroutine MPIGetNodes(n,tn)

   integer, intent(out) :: n, tn

   n = Node
   tn = nProc

end subroutine MPIGetNodes
!****** End subroutine: MPIGetNodes *******************************************
!******************************************************************************


!****** Function: MPIIsPNode **************************************************
!******************************************************************************
!
!  True if is the principal node
!
!******************************************************************************
function MPIIsPNode() result(is)

   logical :: is

   is = IsPNode

end function MPIIsPNode
!****** End function: MPIIsPNode **********************************************
!******************************************************************************


subroutine MPIBarrier()

   call MPITimer(1)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   call MPITimer(0)

end subroutine MPIBarrier

end module mpi
