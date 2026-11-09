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

module sys

   implicit none

   PRIVATE

   integer, save :: nerr = 0, nwarn = 0

   ! Debug output (Output.Debug / Output.DebugModules, read in MIO_Initialize).
   ! Off by default: SysPrintDebug prints nothing unless it is switched on for
   ! every module or for the calling module.
   integer, parameter :: maxdbgmod = 32
   logical, save :: debugAll = .false.
   integer, save :: ndbgmod = 0
   character(len=16), save :: dbgmod(maxdbgmod)

   ! Width of a printed line (module prefix included) before SysPrint wraps it.
   integer, parameter :: printWidth = 110
   ! Width of a section header.
   integer, parameter :: headerWidth = 78

   public :: SysKill
   public :: SysWarning
   public :: SysIOErr
#ifdef MPI
   public :: SysIOErrMPI
#endif /* MPI */
   public :: SysPrint
   public :: SysPrintErrors
   public :: SysKillFlush
   public :: SysSection
   public :: SysSubsection
   public :: SysPrintDebug
   public :: SysDebugging
   public :: SysSetDebug

   interface SysPrint
      module procedure SysPrint_mod, SysPrint__
   end interface


contains

!****** Subroutine: SysKill ***************************************************
!******************************************************************************
!
!  Display an error notification and terminates the program. This function
! shows the module from where it was called, and, in DEBUG mode, the function
! inside that module. With 'hold' variable is possible to wait for a
! SysKillFlush call to terminate the program.
!
! INPUT -----------------------------------------------------------------------
!
! character(*) message            : Error message
! character(*) module             : Module where the error occurred
! character(*) func               : Function where the error occurred
! logical hold                    : If true, the program continues executing
!
!******************************************************************************
subroutine SysKill(message,module,func,hold)

#ifdef MPI
   use mpi
#endif /* MPI */

   implicit none

   character(len=*), intent(in) :: message, module
   character(len=*), intent(in), optional :: func
   logical, intent(in), optional :: hold

   integer :: ierr
#ifdef MPI
   integer(8) :: c0, c1, crate
#endif /* MPI */

   nerr = nerr+1
   call SysPrint('#################################################################')
   call SysPrint('ERROR: '//trim(message),module)
#ifdef DEBUG
   if (present(func)) then
      call SysPrint('*** Error from function '// &
        trim(func)//' in module '//trim(module))
   end if
#endif /* DEBUG */
   if (present(hold)) then
     if (hold) return
   end if
   call SysPrintErrors()
   call SysPrint('*** PROGRAM ABORTED')
   ! MPI_Abort may end the process without flushing buffered output; make
   ! sure the error message above is written first.
   flush(6)
#ifdef MPI
   ! Only the first process prints. When every process stops at the same
   ! check, give it time to write the message before another process aborts
   ! the run.
   if (.not. IsPNode) then
      call system_clock(c0, crate)
      c1 = c0
      do while (real(c1 - c0, kind(1.0d0)) < 2.0d0*crate)
         call system_clock(c1)
      end do
   end if
   call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#endif /* MPI */
   stop

end subroutine SysKill
!****** End subroutine: SysKill ***********************************************
!******************************************************************************

!****** Subroutine: SysWarning ************************************************
!******************************************************************************
!
!  Display a warning notification
!
! INPUT -----------------------------------------------------------------------
!
! character(*) message            : Warning message
! character(*) module             : Module that called this subroutine
! character(*) func               : Function that called this subroutine
!
!******************************************************************************
subroutine SysWarning(message,module,func)

#ifdef MPI
   use mpi
#endif /* MPI */

   implicit none

   character(len=*), intent(in) :: message, module
   character(len=*), intent(in), optional :: func

   nwarn = nwarn+1
   call SysPrint('')
   call SysPrint('WARNING: '//trim(message),module)
#ifdef DEBUG
   if (present(func)) then
      call SysPrint('*** Warning from function '// &
        trim(func)//' in module '//trim(module))
   end if
#endif /* DEBUG */

end subroutine SysWarning
!****** End subroutine: SysWarning ********************************************
!******************************************************************************



!****** Subroutine: SysPrint_mod **********************************************
!******************************************************************************
!
!  Display a message to standard output indicating the module who made the
! call. The message will be split into as many lines as needed so the number of
! characters in each line is no more than 'maxl'.
!
! INPUT -----------------------------------------------------------------------
!
! character(*) message            : General message
! character(*) module             : Module that called this subroutine
!
!******************************************************************************
subroutine SysPrint_mod(message,module)

   use io,                   only : stdout
#ifdef MPI
   use mpi
#endif /* MPI */

   implicit none

   character(len=*), intent(in) :: message, module
   integer :: istr, lstr, indx, iostat, i, maxl
   logical :: first

   lstr = len_trim(message)
   ! room left on a line after "module: "; continuation lines are indented
   ! by two more spaces so that a wrapped message reads as one item
   maxl = max(printWidth - len_trim(module) - 4, 20)
#ifdef MPI
   iostat = 0
   if (IsPNode) then
#endif /* MPI */
      if (lstr > maxl) then
         i = 1
         first = .true.
         do
            indx = min(maxl+i-1,lstr)
            if (lstr-i >= maxl) then
               do istr=indx+1,i+1,-1
                  if (message(istr:istr) == ' ') then
                     indx = istr
                     exit
                  end if
               end do
            end if
            if (first) then
               write(stdout,'(a)',IOSTAT=iostat) trim(module)//': ' &
                 //trim(message(i:indx))
            else
               write(stdout,'(a)',IOSTAT=iostat) trim(module)//':   ' &
                 //trim(adjustl(message(i:indx)))
            end if
            first = .false.
            i = indx+1
            if (i > lstr) exit
            if (len_trim(message(i:))==0 .or. iostat/=0) exit
         end do
      else
         write(stdout,'(a)',IOSTAT=iostat) trim(module)//': '//trim(message)
      end if
#ifdef MPI
   end if
#endif /* MPI */
   call SysIOErr(iostat)

end subroutine SysPrint_mod
!****** End subroutine: SysPrint_mod ******************************************
!******************************************************************************


!****** Subroutine: SysPrint__ ************************************************
!******************************************************************************
!
!  Display a message
!
! INPUT -----------------------------------------------------------------------
!
! character(*) message            : General message
!
!******************************************************************************
subroutine SysPrint__(message)

   use io,                 only : stdout
#ifdef MPI
   use mpi
#endif /* MPI */

   implicit none

   character(len=*), intent(in) :: message
   integer :: iostat

#ifdef MPI
   iostat = 0
   if (IsPNode) then
#endif /* MPI */
      write(stdout,'(a)',IOSTAT=iostat) message
#ifdef MPI
   end if
#endif /* MPI */
   call SysIOErr(iostat)

end subroutine SysPrint__
!****** End subroutine: SysPrint__ ********************************************
!******************************************************************************


!****** Subroutine: SysPrintErrors ********************************************
!******************************************************************************
!
!  Print the count of errors and warnings.
!
! INPUT -----------------------------------------------------------------------
!
! logical final                   : If is not the final count a message is
!                                   displayed ('Current count')
!
!******************************************************************************
subroutine SysPrintErrors(final)

   use format,               only : num2string

   implicit none

   logical, intent(in), optional :: final

   character(len=60) :: mess

   call SysPrint('')
   if (present(final)) then
      if (.not. final) then
        call SysPrint('Current count:')
      end if
   end if
   mess = num2string(nerr)
   if (nerr==1) then
      mess = trim(mess)//' error, '//num2string(nwarn)
   else
      mess = trim(mess)//' errors, '//num2string(nwarn)
   end if
   if (nwarn==1) then
      mess = trim(mess)//' warning'
   else
      mess = trim(mess)//' warnings'
   end if
   call SysPrint(trim(mess),'mio')

end subroutine SysPrintErrors
!****** End subroutine: SysPrintErrors ****************************************
!******************************************************************************


!****** Subroutine: SysKillFlush **********************************************
!******************************************************************************
!
!  If there are errors, the program stops
!
!******************************************************************************
subroutine SysKillFlush()

#ifdef MPI
   use mpi
#endif /* MPI */

   implicit none

   integer :: ierr

   if (nerr==0) return
   call SysPrintErrors()
   call SysPrint('*** PROGRAM ABORTED')
   ! MPI_Abort may end the process without flushing buffered output; make
   ! sure the error message above is written first.
   flush(6)
#ifdef MPI
   call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#endif /* MPI */
   stop

end subroutine SysKillFlush
!****** End subroutine: SysKillFlush ******************************************
!******************************************************************************

!****** Subroutine: SysSection ************************************************
!******************************************************************************
!
!  Print a section header: a blank line, then one line of '=' that carries
! the title and the module. One header per stage of a run (system set-up, neighbour search,
! Hamiltonian, band structure, ...), so that a reader can find the stage in
! the output at a glance.
!
! INPUT -----------------------------------------------------------------------
!
! character(*) title              : Title of the section
! character(*) module             : Module that called this subroutine
!
!******************************************************************************
subroutine SysSection(title,module)

   implicit none

   character(len=*), intent(in) :: title, module
   character(len=headerWidth) :: line
   integer :: lt, lm

   ! "===== Title ========== [module]", one line, headerWidth wide
   call SysPrint('')
   line = repeat('=',headerWidth)
   lt = min(len_trim(title), headerWidth-10)
   lm = len_trim(module)
   line(7:7+lt+1) = ' '//title(1:lt)//' '
   if (8 + lt + lm + 4 <= headerWidth) line(headerWidth-lm-2:headerWidth) = ' ['//trim(module)//']'
   call SysPrint(trim(line))

end subroutine SysSection
!****** End subroutine: SysSection ********************************************
!******************************************************************************


!****** Subroutine: SysSubsection *********************************************
!******************************************************************************
!
!  Print a sub-section title inside a section: "module: --- title ---".
!
!******************************************************************************
subroutine SysSubsection(title,module)

   implicit none

   character(len=*), intent(in) :: title, module

   call SysPrint('--- '//trim(title)//' ---',module)

end subroutine SysSubsection
!****** End subroutine: SysSubsection *****************************************
!******************************************************************************


!****** Subroutine: SysSetDebug ***********************************************
!******************************************************************************
!
!  Switch debug output on or off, for every module (no module given) or for
! one module.
!
!******************************************************************************
subroutine SysSetDebug(on,module)

   implicit none

   logical, intent(in) :: on
   character(len=*), intent(in), optional :: module

   integer :: i

   if (.not. present(module)) then
      debugAll = on
      return
   end if
   do i=1,ndbgmod
      if (trim(dbgmod(i)) == trim(module)) then
         if (.not. on) then
            dbgmod(i) = dbgmod(ndbgmod)
            ndbgmod = ndbgmod - 1
         end if
         return
      end if
   end do
   if (on .and. ndbgmod < maxdbgmod) then
      ndbgmod = ndbgmod + 1
      dbgmod(ndbgmod) = module
   end if

end subroutine SysSetDebug
!****** End subroutine: SysSetDebug *******************************************
!******************************************************************************


!****** Function: SysDebugging ************************************************
!******************************************************************************
!
!  True if debug output is on for this module (or for every module). Use it
! to skip the computation of a diagnostic that is only printed.
!
!******************************************************************************
logical function SysDebugging(module)

   implicit none

   character(len=*), intent(in) :: module

   integer :: i

   SysDebugging = debugAll
   if (SysDebugging) return
   do i=1,ndbgmod
      if (trim(dbgmod(i)) == trim(module)) then
         SysDebugging = .true.
         return
      end if
   end do

end function SysDebugging
!****** End function: SysDebugging ********************************************
!******************************************************************************


!****** Subroutine: SysPrintDebug *********************************************
!******************************************************************************
!
!  Print a debug line, "module: DEBUG: message", only when debug output is on
! for the module (Output.Debug, Output.DebugModules) or when 'force' is true
! (a module's own debug key, e.g. Diag.TAPWDebug).
!
!******************************************************************************
subroutine SysPrintDebug(message,module,force)

   implicit none

   character(len=*), intent(in) :: message, module
   logical, intent(in), optional :: force

   logical :: on

   on = SysDebugging(module)
   if (present(force)) on = on .or. force
   if (on) call SysPrint('DEBUG: '//trim(message),module)

end subroutine SysPrintDebug
!****** End subroutine: SysPrintDebug *****************************************
!******************************************************************************

#include "ioerror.f90"

end module sys
