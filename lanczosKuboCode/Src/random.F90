module random
! http://jblevins.org/log/openmp
!https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
   use mio

   implicit none

   PRIVATE

   !public :: RandSeed
   !public :: RandNum
   public :: RandTest

   public :: rand_t

   integer, parameter :: ns = 4
   integer, parameter :: default_seed(ns) = [521288629, 362436069, 16163801, 1131199299]
   integer, parameter :: prime(8) = [99991, 287291, 299977, 301237, 388673, 456623, 499819, 501173]

   type :: rand_t
      integer :: state(ns) = default_seed
   end type rand_t

contains

subroutine RandSeed(self,thread)

   type(rand_t), intent(inout) :: self
   integer, intent(in) :: thread

   integer :: clock, seed

   logical :: randomSeed

   call MIO_InputParameter('randomSeed',randomSeed,.true.)
   !call system_clock(clock)
   if (randomSeed) then
      seed = 932117 + thread
   else
      seed = 932117
   end if
   !seed = 99991 + clock + thread
   !seed = prime(mod(thread,8)+1) + clock + thread
   self%state(1) = seed
   self%state(2:ns) = default_seed(2:ns)

end subroutine RandSeed


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

   !call MIO_InputParameter('RandomTest',test,.false.)
   !call MIO_InputParameter('RandomTestNumber',ntest,100)
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

end subroutine RandTest

end module random
