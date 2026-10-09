module parallel

   use mio

   implicit none

   PRIVATE

   integer, public, save :: xDiv, yDiv, nDiv, xSubdiv, ySubdiv
   integer, public, save :: procID ! From 1 to nDiv
   !$OMP ThreadPrivate(procID)
   integer, public, pointer :: nTh(:)

   public :: ParallelDiv

contains

subroutine ParallelDiv()

   integer, parameter :: primes(17) = (/1,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53/)
   integer, parameter :: szpr = size(primes)

   integer ::sq, i, j

#ifdef MPI
   call MPIAllRedSum(numThreads,nDiv,1,MPI_INTEGER)
   call MIO_Allocate(nTh,nProc,'nTh','parallel')
   call MPIAllGather(numThreads,1,nTh,1,MPI_INTEGER)
   nDiv = 1
   nTh = 1
   ! With nDiv forced to 1 the first process owns every atom and the others
   ! own none, while most routines still loop over all atoms on every
   ! process. Stop here instead of overrunning the per-process arrays.
   if (nProc > 1) then
      call MIO_Kill('This version of GRABNES must be run with a single MPI '// &
        'process: the MPI domain decomposition is disabled. Use OpenMP '// &
        'threads for parallel execution.','parallel','ParallelDiv')
   end if
#else
   nDiv = numThreads
   call MIO_Allocate(nTh,1,'nTh','parallel')
   nTh = numThreads
#endif /* MPI */
   sq = ceiling(sqrt(real(nDiv)))
   if (mod(nDiv,sq)==0) then
      xDiv = sq
      yDiv = nDiv/xDiv
   else
      do i=1,szpr-1
         if (primes(i)>sq) exit
      end do
      do j=i,1,-1
         if (mod(nDiv,primes(j))==0) then
            xDiv = primes(j)
            yDiv = nDiv/xDiv
            exit
         end if
      end do
   end if
   if (nDiv/=1) then
      call MIO_Print('Cell divided in '//trim(num2str(xDiv))//'x'// &
        trim(num2str(yDiv))//' subcells','parallel')
      call MIO_Print('')
   end if
#ifdef MPI
   !$OMP Parallel
   procID = 1
   !$OMP End Parallel
#else
   !$OMP Parallel
   procID = nThread + 1
   !$OMP End Parallel
#endif /* MPI */

end subroutine ParallelDiv

end module parallel
