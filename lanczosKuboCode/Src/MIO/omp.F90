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

module omp

   use omp_lib

   implicit none

   PRIVATE

#ifndef MPI
   integer, public, save :: Node=0, nProc=1
#endif
   integer, public, save :: numThreads=1, nThread=0, minCy
   !$OMP THREADPRIVATE(nThread)

   public :: OmpInit

contains

subroutine OmpInit()

   !$OMP PARALLEL
   !$ nThread = OMP_Get_Thread_Num()
   !$OMP END PARALLEL
   !$ numThreads = OMP_Get_Max_Threads()
   minCy = numThreads*10

end subroutine

end module omp
