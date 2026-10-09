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

module timer
#ifdef TIMER

   use precision,            only : dp

   implicit none

   PRIVATE

   integer, parameter :: tmax = 100, printmax = 8, lmod = 30

   type :: time_array
      integer :: ncalls              ! Number of calls to this module
      real(dp) :: itime              ! Initial time
      real(dp) :: ctime              ! Current time
      logical :: active              ! Flag indicating if the timer is active
      character(len=lmod) :: module    ! Name of module
   end type time_array

   integer, save :: ntimers=0, ncalls=0
#ifdef MPI
   type(time_array), save :: timers(tmax+2)
#else
   type(time_array), save :: timers(tmax+1)
#endif /* MPI */

   real(dp), save :: timetime=0.0_dp

   public :: TimerCount
   public :: TimerStop
   public :: TimerPrint

contains

!****** Subroutine: TimerCount ************************************************
!******************************************************************************
!
!  Count the time for a specified module
!
! INPUT -----------------------------------------------------------------------
!
! character(*) module            : Name of the module
!
!******************************************************************************
subroutine TimerCount(module)

   use sys,                  only : SysKill

   implicit none

   character(len=*), intent(in) :: module
   integer :: it
   real(dp) :: t, t1

   call cpu_time(t)              ! Get current time

   !$OMP CRITICAL(timer_critical)
   ncalls = ncalls+1
   !$OMP END CRITICAL(timer_critical)

   ! Look for existing timers for this module
   do it=ntimers,1,-1
      if(timers(it)%module == module) then
         if (timers(it)%active) then
            call SysKill('Timer '//trim(module)//' already running','timer','TimerCount')
         end if
         !$OMP CRITICAL(timer_critical)
         timers(it)%itime = t
         timers(it)%ncalls = timers(it)%ncalls+1
         timers(it)%active = .true.
         call cpu_time(t1)
         timetime = timetime + t1-t
         !$OMP END CRITICAL(timer_critical)
         return
      end if
   end do

   ! New timer initialization
   !$OMP CRITICAL(timer_critical)
   ntimers = ntimers+1
   if (ntimers > tmax) &
     call SysKill('Not enough timers, parameter tmax too small','timer',&
     'TimerCount')
   timers(ntimers)%module = module
   timers(ntimers)%itime = t
   timers(ntimers)%ctime = 0.0_dp
   timers(ntimers)%ncalls = 1
   timers(ntimers)%active = .true.

   call cpu_time(t1)
   timetime = timetime + t1-t
   !$OMP END CRITICAL(timer_critical)

end subroutine TimerCount
!****** Subroutine: TimerCount ************************************************
!******************************************************************************


!****** Subroutine: TimerStop *************************************************
!******************************************************************************
!
!  Stop the timer for a specified module
!
! INPUT -----------------------------------------------------------------------
!
! character(*) module            : Name of the module
!
!******************************************************************************
subroutine TimerStop(module)

   use sys,                  only : SysKill

   implicit none

   character(len=*), intent(in) :: module
   integer :: it
   real(dp) :: t, t1

   call cpu_time(t)              ! Get current time

   !$OMP CRITICAL(timer_critical)
   ncalls=ncalls+1
   !$OMP END CRITICAL(timer_critical)

   ! Look for existing timers
   do it=ntimers,1,-1
      if(timers(it)%module == module) then
         if (.not. timers(it)%active) then
            call SysKill('Timer '//trim(module)//' not running','timer','TimerStop')
         end if
         !$OMP CRITICAL(timer_critical)
         timers(it)%ctime = timers(it)%ctime + t - timers(it)%itime
         timers(it)%active = .false.
         call cpu_time(t1)
         timetime = timetime + t1-t
         !$OMP END CRITICAL(timer_critical)
         return
      end if
   end do

   ! Timer not found
   call SysKill('Timer '//trim(module)//' not found','timer','TimerStop')

end subroutine TimerStop
!****** End subroutine: TimerStop *********************************************
!******************************************************************************


!****** Subroutine: TimerPrint ************************************************
!******************************************************************************
!
!  Print timers
!
!******************************************************************************
subroutine TimerPrint()

#ifdef MPI
   use mpi
#endif /* MPI */
   use sys,                  only : SysWarning, SysPrint
   use mem,                  only : MemAlloc, MemDealloc

   implicit none

   type(time_array) :: tprint(printmax)
   integer :: it, imt, nt, i
   real(dp) :: t, tot, t1
   character(len=50) :: message
#ifdef MPI
   character(len=lmod), pointer :: modules(:) => NULL()
   real(dp),pointer :: tArr(:) => NULL(), tArr0(:) => NULL()
   integer, pointer :: callsArr(:) => NULL(), callsArr0(:) => NULL()
#endif /* MPI */

   call cpu_time(t)              ! Get current time
#ifdef MPI
   call MemAlloc(tArr,ntimers+2,'tArr','TimerPrint')
   call MemAlloc(tArr0,ntimers+2,'tArr0','TimerPrint')
   call MemAlloc(callsArr,ntimers+2,'callsArr','TimerPrint')
   call MemAlloc(callsArr0,ntimers+2,'callsArr0','TimerPrint')
   call MemAlloc(modules,ntimers+2,'modules','TimerPrint')
#endif /* MPI */
   do it=2,ntimers
      if (timers(it)%active) then
         call SysWarning('Timer not closed for module '// &
           trim(timers(it)%module),'timer','TimerPrint')
         timers(it)%active = .false.
         timers(it)%ctime = timers(it)%ctime + t - timers(it)%itime
      end if
   end do
#ifdef MPI
   if (IsPNode) then
      do it=1,ntimers
         modules(it) = timers(it)%module
      end do
   endif
   call MPIBCast(modules,lmod*ntimers,MPI_CHARACTER)
   do it=1,ntimers
      imt = it
      if (timers(imt)%module==modules(it)) then
         tArr(it) = timers(imt)%ctime
         callsArr(it) = timers(imt)%ncalls
      else
         tArr(it) = 0.0_dp
         callsArr(it) = 0
         do imt=2,ntimers
            if (timers(imt)%module==modules(it)) then
               tArr(it) = timers(imt)%ctime
               callsArr(it) = timers(imt)%ncalls
               exit
            end if
         end do
      end if
   end do
   ntimers = ntimers+1
   call MPITimerGet(tArr(ntimers),callsArr(ntimers))
   timers(ntimers)%module = 'mpi'
#endif /* MPI */
   !$OMP CRITICAL(timer_critical)
   ncalls=ncalls+1
   ntimers = ntimers+1
   call cpu_time(t1)
   timetime = timetime + t1-t
   timers(ntimers)%module = 'timer'
   timers(ntimers)%ncalls = ncalls
   !$OMP END CRITICAL(timer_critical)
   if (timers(1)%active) then
      call TimerStop(timers(1)%module)
   end if
#ifdef MPI
   tArr(1) = timers(1)%ctime
   callsArr(1) = timers(1)%ncalls
   tArr(ntimers) = timetime
   callsArr(ntimers) = ncalls
   call MPIRedSum(callsArr,callsArr0,ntimers,MPI_INTEGER)
   call MPIRedSum(tArr,tArr0,ntimers,MPI_DOUBLE_PRECISION)
   if (IsPNode) then
      do it=1,ntimers
         timers(it)%ncalls = callsArr0(it)
         timers(it)%ctime = tArr0(it)
      end do
#endif /* MPI */
      ! Find timers with greater times
      nt = 0
      do it=2,ntimers
         if (nt < printmax) then
            nt = nt+1
            tprint(nt)%ncalls = timers(it)%ncalls
            tprint(nt)%itime = timers(it)%itime
            tprint(nt)%ctime = timers(it)%ctime
            tprint(nt)%active = timers(it)%active
            tprint(nt)%module = timers(it)%module
         end if
         do imt=1,nt
            if (timers(it)%ctime > tprint(imt)%ctime) then
               if (imt/=printmax) then
                  do i=printmax,imt+1,-1
                     tprint(i)%ncalls = tprint(i-1)%ncalls
                     tprint(i)%itime = tprint(i-1)%itime
                     tprint(i)%ctime = tprint(i-1)%ctime
                     tprint(i)%active = tprint(i-1)%active
                     tprint(i)%module = tprint(i-1)%module
                  end do
               end if
               tprint(imt)%ncalls = timers(it)%ncalls
               tprint(imt)%itime = timers(it)%itime
               tprint(imt)%ctime = timers(it)%ctime
               tprint(imt)%active = timers(it)%active
               tprint(imt)%module = timers(it)%module
               exit
            end if
         end do
      end do
#ifdef MPI
   end if
#endif /* MPI */
   call cpu_time(t1)
   tot = timers(1)%ctime + t1-t
   if (tot==0.0_dp) tot=epsilon(1.0_dp)
   call SysPrint('')
   call SysPrint('    Module     Calls      Time (s)       % Total','timer')
   call SysPrint('  -----------------------------------------------','timer')
#ifdef MPI
   if (IsPNode) then
#endif /* MPI */
      do it=1,nt
         write(message,'(a,i8,f14.2,f14.2)') '   '//tprint(it)%module(1:9),&
           tprint(it)%ncalls,tprint(it)%ctime, 100.0_dp*tprint(it)%ctime/tot
         call SysPrint(trim(message),'timer')
      end do
      write(message,'(a,i8,f14.2,f14.2)') '   '//timers(1)%module(1:9),&
        timers(1)%ncalls,tot, 100.0_dp
      call SysPrint(message,'timer')
#ifdef MPI
   end if
   call MemDealloc(modules,'modules','TimerPrint')
   call MemDealloc(callsArr0,'callsArr0','TimerPrint')
   call MemDealloc(callsArr,'callsArr','TimerPrint')
   call MemDealloc(tArr0,'tArr0','TimerPrint')
   call MemDealloc(tArr,'tArr','TimerPrint')
#endif /* MPI */

end subroutine TimerPrint
!****** End subroutine: TimerPrint ********************************************
!******************************************************************************


#endif /*TIMER */
end module timer
