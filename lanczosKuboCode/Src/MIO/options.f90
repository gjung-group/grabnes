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

module options

   use io,                   only : maxfnlen

   implicit none

   PRIVATE

   character(len=maxfnlen) :: input

   public :: OptionsProcess
   public :: OptionsInputName

contains

!****** Subroutine: OptionsProcess ********************************************
!******************************************************************************
!
!  Process input arguments to the program
!
!******************************************************************************
subroutine OptionsProcess()

   implicit none

   integer :: narg, ia
   character(len=maxfnlen) :: arg

   narg = command_argument_count()
   ia = 0
   input = ''
   do while(ia/=narg)
      ia = ia+1
      call get_command_argument(ia,arg)
      if (ia==1 .or. ia==narg) then
         input = arg
      end if
   end do

end subroutine OptionsProcess
!****** End subroutine: OptionsProcess ****************************************
!******************************************************************************


!****** Subroutine: OptionsInputName ******************************************
!******************************************************************************
!
!  Returns the input file name found in input options
!
!******************************************************************************
function OptionsInputName() result(inname)

   implicit none

   character(len=maxfnlen) :: inname

   inname = input

end function OptionsInputName
!****** End subroutine: OptionsInputName **************************************
!******************************************************************************

end module options
