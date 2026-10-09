module random
! http://jblevins.org/log/openmp
!https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
   use mio

   implicit none

   PRIVATE

   !public :: RandSeed
   !public :: RandNum
   public :: RandTest
   public :: RandSeedFromInput

   public :: rand_t

   integer, parameter :: ns = 4
   integer, parameter :: default_seed(ns) = [521288629, 362436069, 16163801, 1131199299]
   integer, parameter :: prime(8) = [99991, 287291, 299977, 301237, 388673, 456623, 499819, 501173]

   type :: rand_t
      integer :: state(ns) = default_seed
   end type rand_t

contains

!> @brief Seed the intrinsic random-number generator for the disorder models.
!! @details With setSeed .true. the seed is seedValue, so that a disordered
!!          Hamiltonian can be reproduced; otherwise it comes from the system
!!          clock, as before.
subroutine RandSeedFromInput()

   integer :: n, clock, i, seedValue
   integer, allocatable :: seed(:)
   logical :: setSeed

   call random_seed(size = n)
   allocate(seed(n))
   call MIO_InputParameter('setSeed',setSeed,.false.)
   call MIO_InputParameter('seedValue',seedValue,123456)
   if (setSeed) then
      seed = seedValue
   else
      call system_clock(COUNT=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   end if
   call random_seed(PUT = seed)
   deallocate(seed)

end subroutine RandSeedFromInput

function RandNum(self) result(rand)

   type(rand_t), intent(inout) :: self
   real(dp) :: rand

   integer :: imz

   imz = self%state(1) - self%state(3)

   if (imz < 0) imz = imz + 2147483579

   self%state(1) = self%state(2)
   self%state(2) = self%state(3)
   self%state(3) = imz
   self%state(4) = 69069 * self%state(4) + 1013904243
   imz = imz + self%state(4)
   rand = 0.5_dp + 0.23283064d-9 * imz

end function RandNum

subroutine RandTest(rng,n)

   use name,                 only : prefix

   type(rand_t), intent(inout) :: rng
   integer, intent(in) :: n

   character(60) :: filename
   integer :: u, i
   real(dp) :: s, r
   type(cl_file) :: file

#ifdef DEBUG
   call MIO_Debug('RandTest',0)
#endif /* DEBUG */

   filename = trim(prefix)//'.'//trim(num2str(nThread))//'.RNDM'
   u = 200+nThread
   open(u,FILE=filename,STATUS='replace')
   s = 0.0_dp
   do i=1,n
      r = RandNum(rng)
      write(u,*) r
      s = s + r
   end do
   s = s/real(n,8)
   write(u,*)
   write(u,*) 'Av:', s
   call MIO_Print('Average: '//trim(num2str(s,12)),'random')
   call MIO_Print('')
   close(u)

#ifdef DEBUG
   call MIO_Debug('RandTest',1)
#endif /* DEBUG */

end subroutine RandTest

end module random
