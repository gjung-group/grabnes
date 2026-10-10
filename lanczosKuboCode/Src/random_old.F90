module random

   use mio

   implicit none

#ifdef MPI
#define USE_MPI 1
#include <sprng_f.h>
#endif /* MPI */

   PRIVATE

#ifdef MPI
   SPRNG_POINTER, save :: stream
   !$OMP THREADPRIVATE(stream)
#endif /* MPI */

   public :: RandInit
   public :: RandGen

contains

subroutine RandInit()

   use parallel,             only : nDiv, procID
   use name,                 only : prefix

   integer :: seed
   integer, pointer :: seedS(:)
   integer :: nk, clock, ik
   logical :: test
   integer :: ntest, i, u
   type(cl_file) :: file
   real(dp) :: r, s

#ifdef DEBUG
   call MIO_Debug('RandInit',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('random')
#endif /* TIMER */

   ! Initialize random seed
#ifdef MPI
   if (nDiv>1) then
      seed = make_sprng_seed()
      !$OMP PARALLEL
      stream = init_sprng(4,procID-1,nDiv,seed,SPRNG_DEFAULT)
      !$OMP END PARALLEL
   else
#endif /* MPI */
      call random_seed(size=nk)
      allocate(seedS(nk))
      call system_clock(clock)
      seedS = clock + 37*(/(ik,ik=0,nk-1)/)
      call random_seed(put=seedS)
      deallocate(seedS)
#ifdef MPI
   end if
#endif /* MPI */
   call MIO_InputParameter('Random.Test',test,.false.)
   if (test) then
      call MIO_InputParameter('Random.TestNumber',ntest,100)
      call file%Open(name=trim(prefix)//'.RNDM',serial=.true.)
      u = file%GetUnit()
      s = 0.0_dp
      do i=1,ntest
         r = RandGen()
         write(u,*) r
         s = s + r
      end do
      s = s/real(ntest,8)
      call MIO_Print('Average: '//trim(num2str(s,12)),'random')
      call MIO_Print('')
      call file%Close()
   end if

#ifdef TIMER
   call MIO_TimerStop('random')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('RandInit',1)
#endif /* DEBUG */

end subroutine RandInit


function RandGen() result(rand)

   use parallel,             only : nDiv

   real(dp) :: rand

#ifdef MPI
   if (nDiv>1) then
      rand = sprng(stream)
   else
#endif /* MPI */
      call random_number(rand)
#ifdef MPI
   end if
#endif /* MPI */

end function RandGen

end module random
