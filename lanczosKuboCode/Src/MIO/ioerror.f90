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


!****** Subroutine: SysIOErr **************************************************
!******************************************************************************
!
!  Error handler
!
! INPUT -----------------------------------------------------------------------
!
! integer iostat                  : Status ID
! char module                     : Module that called this function
! char function                   : The function inside the module
!
!******************************************************************************
subroutine SysIOErr(iostat,module,function)

   use format,               only : num2string

   implicit none

   integer, intent(in) :: iostat
   character(len=*), intent(in), optional :: module, function

   if (iostat==0) return
   if (present(module)) then
      call SysPrint('')
#ifdef DEBUG
      if (present(function)) then
         call SysPrint('I/O error in function '//trim(function),module)
      else
#endif /* DEBUG */
         call SysPrint('I/O error',module)
#ifdef DEBUG
      endif
#endif /* DEBUG */
   end if
   select case(iostat)
   case(-3)
      call SysKill('runtime error: Bad logical expression', &
        'ioerror')
   case(-2)
      call SysKill('runtime error: Bad integer number', &
        'ioerror')
   case(-1)
      call SysKill('severe (-1): End-of-record during read', &
        'ioerror')
   case(1)
      call SysKill('severe (1): Not a Fortran-specific error', &
        'ioerror')
   case(8)
      call SysKill('severe (8): Internal consistency check failure', &
        'ioerror')
   case(9)
      call SysKill('severe (9): Permission to access file denied', &
        'ioerror')
   case(10)
      call SysKill('severe (10): Cannot overwrite existing file', &
        'ioerror')
   case(11)
      call SysKill('info (11): Unit not connected', &
        'ioerror')
   case(12)
      call SysKill('syntax error: Uncompensated parentheses', &
        'ioerror')
   case(13)
      call SysKill('syntax error: Empty field', &
        'ioerror')
   case(14)
      call SysKill('syntax error: Unrecognized statement', &
        'ioerror')
   case(15)
      call SysKill('runtime error: Bad real number', &
        'ioerror')
   case(16)
      call SysKill('syntax error: Wrong expression','ioerror')
   case(17)
      call SysKill('severe (17): Syntax error in NAMELIST input', &
        'ioerror')
   case(18)
      call SysKill('severe (18): Too many values for NAMELIST variable', &
        'ioerror')
   case(19)
      call SysKill('severe (19): Invalid reference to variable in NAMELIST' &
        //' input', 'ioerror')
   case(20)
      call SysKill('severe (20): REWIND error', &
        'ioerror')
   case(21)
      call SysKill('severe (21): Duplicate file specifications', &
        'ioerror')
   case(22)
      call SysKill('severe (22): Input record too long', &
        'ioerror')
   case(23)
      call SysKill('severe (23): BACKSPACE error', &
        'ioerror')
   case(24)
      call SysKill('severe (24): End-of-file during read', &
        'ioerror')
   case(25)
      call SysKill('severe (25): Record number outside range', &
        'ioerror')
   case(26)
      call SysKill('severe (26): OPEN or DEFINE FILE required', &
        'ioerror')
   case(27)
      call SysKill('severe (27): Too many records in I/O statement', &
        'ioerror')
   case(28)
      call SysKill('severe (28): CLOSE error', &
        'ioerror')
   case(29)
      call SysKill('severe (29): File not found', &
        'ioerror')
   case(30)
      call SysKill('severe (30): Open failure', &
        'ioerror')
   case(31)
      call SysKill('severe (31): Mixed file access modes', &
        'ioerror')
   case(32)
      call SysKill('severe (32): Invalid logical unit number', &
        'ioerror')
   case(33)
      call SysKill('severe (33): ENDFILE error', &
        'ioerror')
   case(34)
      call SysKill('severe (34): Unit already open', &
        'ioerror')
   case(35)
      call SysKill('severe (35): Segmented record format error', &
        'ioerror')
   case(36)
      call SysKill('severe (36): Attempt to access non-existent record', &
        'ioerror')
   case(37)
      call SysKill('severe (37): Inconsistent record length', &
        'ioerror')
   case(38)
      call SysKill('severe (38): Error during write', &
        'ioerror')
   case(39)
      call SysKill('severe (39): Error during read', &
        'ioerror')
   case(40)
      call SysKill('severe (40): Recursive I/O operation', &
        'ioerror')
   case(41)
      call SysKill('severe (41): Insufficient virtual memory', &
        'ioerror')

   case(45)
      call SysKill('severe (45): Keyword value error in OPEN statement', &
        'ioerror')
   case(47)
      call SysKill('severe (47): Write to READONLY file', &
        'ioerror')
   case(59)
      call SysKill('severe (59): list-directed I/O syntax error', &
        'ioerror')
   case(64)
      call SysKill('severe (64): Input conversion error', &
        'ioerror')
   case(66)
      call SysKill('severe (66): Output statement overflows record', &
        'ioerror')
   case default
      call SysKill('Unrecognized error ('//&
        trim(num2string(iostat))//')','ioerror')
   end select

end subroutine SysIOErr
!****** End subroutine: SysIOErr **********************************************
!******************************************************************************

#ifdef MPI
!****** Subroutine: SysIOErrMPI ***********************************************
!******************************************************************************
!
!  Error handler
!
! INPUT -----------------------------------------------------------------------
!
! integer iostat                  : Status ID
! char module                     : Module that called this function
! char function                   : The function inside the module
!
!******************************************************************************
subroutine SysIOErrMPI(iostat,module,function)

   use format,               only : num2string
   use mpi

   implicit none

   integer, intent(in) :: iostat
   character(len=*), intent(in), optional :: module, function

   if (iostat==0) return
   if (present(module)) then
      call SysPrint('')
#ifdef DEBUG
      if (present(function)) then
         call SysPrint('I/O error in function '//trim(function),module)
      else
#endif /* DEBUG */
         call SysPrint('I/O error',module)
#ifdef DEBUG
      endif
#endif /* DEBUG */
   end if
   select case(iostat)
      case (MPI_ERR_FILE)
         call SysKill('MPI I/O error: Invalid file handle', &
           'ioerror')
      case (MPI_ERR_NOT_SAME)
         call SysKill('Collective argument not identical on all processes,'//&
           'or collective routines called in a different order by different'//&
           'processes', 'ioerror')
      case (MPI_ERR_AMODE)
         call SysKill('Error related to the amode passed to MPI_FILE_OPEN', &
           'ioerror')
      case (MPI_ERR_UNSUPPORTED_DATAREP)
         call SysKill('Unsupported datarep passed to MPI_FILE_SET_VIEW', &
           'ioerror')
      case (MPI_ERR_UNSUPPORTED_OPERATION)
         call SysKill('Unsupported operation, such as seeking on a file '//&
           'which supports sequential access only', 'ioerror')
      case (MPI_ERR_NO_SUCH_FILE)
         call SysKill('File does not exist', 'ioerror')
      case (MPI_ERR_FILE_EXISTS)
         call SysKill('File exists', 'ioerror')
      case (MPI_ERR_BAD_FILE)
         call SysKill('Invalid file name', 'ioerror')
      case (MPI_ERR_ACCESS)
         call SysKill('Permission denied', 'ioerror')
      case (MPI_ERR_NO_SPACE)
         call SysKill('Not enough space', 'ioerror')
      case (MPI_ERR_QUOTA)
         call SysKill('Quota exceeded', 'ioerror')
      case (MPI_ERR_READ_ONLY)
         call SysKill('Read-only file or file system', 'ioerror')
      case (MPI_ERR_FILE_IN_USE)
         call SysKill('File operation could not be completed, as the file'//&
           ' is currently open by some process', 'ioerror')
      case (MPI_ERR_DUP_DATAREP)
         call SysKill('Conversion functions could not be registered because'//&
           ' a data representation identifier that was already defined was '//&
           'passed to MPI_REGISTER_DATAREP', 'ioerror')
      case (MPI_ERR_CONVERSION)
         call SysKill('An error occurred in a user supplied data conversion'//&
           ' function', 'ioerror')
      case default
         call SysKill('Unrecognized error ('//&
           trim(num2string(iostat))//')','ioerror')
   end select

end subroutine SysIOErrMPI
!****** End subroutine: SysIOErrMPI *******************************************
!******************************************************************************
#endif /* MPI */
