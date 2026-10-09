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

module debug

   use file,                 only : cl_file
   use io,                   only : udebug

   implicit none

   PRIVATE

   type(cl_file), save :: debugFile

   public :: DebugInit
   public :: DebugPrint
   public :: DebugClose

contains


!****** Subroutine: DebugInit *************************************************
!******************************************************************************
!
!  Initialize debug file 'debug.log'
!
!******************************************************************************
subroutine DebugInit()
#ifdef DEBUG

   call debugFile%Open(name='debug.log',unit=udebug)

#endif /* DEBUG */
end subroutine DebugInit
!****** End subroutine: DebugInit *********************************************
!******************************************************************************


!****** Subroutine: DebugPrint ************************************************
!******************************************************************************
!
!  Write debug information. All debug information is written to file debug.log.
! This routine prints a report on the routines that are being used to follow
! track of the execution of the program.
!
! INPUT -----------------------------------------------------------------------
!
! character(*) subroutine         : Calling subroutine
! integer i                       : 0 - Entering routine
!                                   otherwise - Routine finished
!
!******************************************************************************
subroutine DebugPrint(subroutine, i)

#ifdef MPI
   use mpi,                  only : IsPNode
#endif /* MPI */

   character(len=*), intent(in) :: subroutine
   integer, intent(in) :: i

#ifdef DEBUG
#ifdef MPI
   if (IsPNode) then
#endif /* MPI */
      if (i==0) then
         write(udebug,'(a)') 'DEBUG: Entering subroutine '//trim(subroutine)
      else
         write(udebug,'(a)') 'DEBUG: Finished subroutine '//trim(subroutine)
      end if
      call flush(udebug)
#ifdef MPI
   end if
#endif /* MPI */

#endif /* DEBUG */
end subroutine DebugPrint
!****** End subroutine: DebugPrint ********************************************
!******************************************************************************


!****** Subroutine: DebugClose ************************************************
!******************************************************************************
!
!  Close debug file
!
!******************************************************************************
subroutine DebugClose()
#ifdef DEBUG

   call debugFile%Close()

#endif /* DEBUG */
end subroutine DebugClose
!****** End subroutine: DebugClose ********************************************
!******************************************************************************

end module debug
