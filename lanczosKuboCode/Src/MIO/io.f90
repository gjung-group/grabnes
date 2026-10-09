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

module io

   implicit none

   PRIVATE

   integer, public, save :: stdin = 5, stdout = 6, stderr = 0
   integer, public, parameter :: unit0 = 11, udebug = 10
   integer, public, parameter :: maxlinel=80, maxfnlen=60

end module io
