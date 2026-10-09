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

   public :: SysKill
   public :: SysWarning
   public :: SysIOErr
#ifdef MPI
   public :: SysIOErrMPI
#endif /* MPI */
   public :: SysPrint
   public :: SysPrintErrors
   public :: SysKillFlush

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

   integer, parameter :: maxl = 66

   character(len=*), intent(in) :: message, module
   integer :: istr, lstr, indx, iostat, nl, il, i

   lstr = len_trim(message)
#ifdef MPI
   iostat = 0
   if (IsPNode) then
#endif /* MPI */
      if (lstr > maxl) then
         nl = lstr/maxl+2
         i = 1
         do il = 1,nl
            indx = min(maxl+i-1,lstr)
            if (lstr-i >= maxl) then
               do istr=indx+1,i,-1
                  if (message(istr:istr) == ' ') then
                     indx = istr
                     exit
                  end if
               end do
               if (indx==i) then
                  indx = min(maxl+i-1,lstr)
                  write(*,*) 'MIN', maxl+i-1,lstr,indx
               end if
            end if
            write(stdout,'(a)',IOSTAT=iostat) module//': ' &
              //message(i:indx)
            i = indx+1
            if (len_trim(message(i:))==0 .or. iostat/=0) exit
         end do
      else
         write(stdout,'(a)',IOSTAT=iostat) module//': '//trim(message)
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
   call SysPrint(mess)

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

#include "ioerror.f90"

end module sys
