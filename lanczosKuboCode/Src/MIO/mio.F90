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

module mio

   use debug
   use file
   use format
   use input,      MIO_InputParameter => InputParameter
   use input,      MIO_InputParameterVariable => InputParameterVariable
   use input,      MIO_InputFindBlock =>  InputFindBlock, &
                   MIO_InputBlock => InputBlock, &
                   MIO_InputSearchLabel => InputSearchLabel
   use io
   use mem,        MIO_Allocate => MemAlloc, MIO_Deallocate => MemDealloc, &
                   MIO_MemGetName  => MemGetName, MIO_MemCount => MemCount
   use omp
   use parser,     MIO_ParserFindValue => ParserFindValue
   use precision
   use string,     MIO_StringComp => StringComp
   use sys,        MIO_Warning => SysWarning, MIO_Kill => SysKill, &
                   MIO_KillFlush => SysKillFlush, MIO_Print => SysPrint, &
                   MIO_IOErr => SysIOErr
   use timer

#ifdef MPI
   use mpi
#endif

   implicit none

   !debug module
   public :: MIO_Debug
   interface MIO_Debug
      module procedure DebugPrint
   end interface

   ! file module
   public :: cl_file

   public :: MIO_FileCloseAll
   interface MIO_FileCloseAll
      module procedure CloseAll
   end interface

   ! format module
   public :: num2str
   interface num2str
      module procedure num2stri, num2strd, num2strz, num2strr, num2strc
   end interface

   public :: num2month
   interface num2month
      module procedure num2month
   end interface

   public :: num2bytes
   interface num2bytes
      module procedure num2bytes
   end interface

   ! input module
   public :: MIO_InputParameter
   public :: MIO_InputParameterVariable

   ! mem module
   public :: MIO_Allocate
   public :: MIO_Deallocate

   ! mio module
   public :: MIO_Initialize
   public :: MIO_Finalize

   ! sys module
   public :: MIO_Print

   public :: MIO_Warning
   public :: MIO_Kill, MIO_KillFlush

   ! timer module
   public :: MIO_TimerCount
   interface MIO_TimerCount
      module procedure TimerCount
   end interface

   public :: MIO_TimerStop
   interface MIO_TimerStop
      module procedure TimerStop
   end interface

contains

!****** Subroutine: MIO_Initialize ********************************************
!******************************************************************************
!
!  Initialize MIO library
!
! INPUT -----------------------------------------------------------------------
!
! character(*) program            : Name of the program
!
!******************************************************************************
subroutine MIO_Initialize(program,success)

   use input,                only : InputInit
   use options,              only : OptionsProcess

   character(len=*), intent(in) :: program
   logical, intent(out) :: success

   character(len=10) :: date, time, zone
   character(len=65) :: mess
   integer :: values(8)

#ifdef TIMER
   call TimerCount(program)
#endif /* TIMER */
#ifdef MPI
   call MPIIni()
#endif /* MPI */
#ifdef DEBUG
   call DebugInit()
#endif /* DEBUG */
   call date_and_time(date,time,zone,values)
   mess = '* Program started on '//trim(num2month(values(2)))//' '&
     //date(7:8)//' '//date(1:4)//' at '//time(1:2)//':'//time(3:4)// &
     ':'//time(5:6)
   call MIO_Print(mess)
   call OmpInit()
   call OptionsProcess()
   call InputInit(success)

end subroutine MIO_Initialize
!****** End subroutine: MIO_Initialize ****************************************
!******************************************************************************


!****** Subroutine: MIO_Finalize **********************************************
!******************************************************************************
!
!  Finalize MIO library
!
!******************************************************************************
subroutine MIO_Finalize()

   use input,                only : InputClose

   character(len=10) :: date, time, zone
   character(len=55) :: mess
   integer :: values(8)

   call InputClose()
#ifdef DEBUG
   call DebugClose()
#endif /* DEBUG */
   call MIO_Print('.......................................................'&
     //'..................')
#ifdef TIMER
   call TimerPrint()
#endif /* TIMER */
   call MemPrint()
   call CloseAll()
   call SysPrintErrors()
   call date_and_time(date,time,zone,values)
   mess = '* Program finished on '//trim(num2month(values(2)))//' '&
     //date(7:8)//' '//date(1:4)//' at '//time(1:2)//':'//time(3:4)// &
     ':'//time(5:6)
   call MIO_Print(mess)
#ifdef MPI
   call MPIFin()
#endif /* MPI */

end subroutine MIO_Finalize
!****** End subroutine: MIO_Finalize ******************************************
!******************************************************************************

end module mio
