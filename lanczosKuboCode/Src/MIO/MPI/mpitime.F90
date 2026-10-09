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

module mpitime

   use mpikinds,             only : dp

   implicit none

   PRIVATE

   real(dp), save :: ctime=0.0_dp, itime=0.0_dp
   integer, save :: ncalls=0
   logical, save :: active=.false.

   public :: MPITimer
   public :: MPITimerGet

contains


!****** Subroutine: MPITimer **************************************************
!******************************************************************************
!
!  Count the time for MPI modules
!
! INPUT -----------------------------------------------------------------------
!
! integer action                  : if action = 1, count
!                                   else, stop
!
!******************************************************************************
subroutine MPITimer(action)
   integer, intent(in) :: action
#ifdef TIMER

   real(dp) :: t

   call cpu_time(t)
   if (action==1) then
      if (active) then
         !stop 'MPITimer: the counter is already active'
      else
         itime = t
      end if
      ncalls = ncalls+1
      active = .true.
   else
      if (active) then
         ctime = ctime + t - itime
      else
         !stop 'MPITimer: the counter is not active'
      end if
      active = .false.
   end if

#endif /* TIMER */
end subroutine MPITimer
!****** End subroutine: MPITimer **********************************************
!******************************************************************************


!****** Subroutine: MPITimerGet ***********************************************
!******************************************************************************
!
!  Get the time
!
!******************************************************************************
subroutine MPITimerGet(t,nc)

   real(dp), intent(out) :: t
   integer, intent(out) :: nc

   t = ctime
   nc = ncalls

end subroutine MPITimerGet
!****** End subroutine: MPITimerGet *******************************************
!******************************************************************************

end module mpitime
