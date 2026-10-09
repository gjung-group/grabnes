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

module format

   implicit none

   PRIVATE

   integer, parameter :: strlen = 25

   public :: num2string
   public :: num2month
   public :: num2bytes
   public :: lowercase
   public :: uppercase

   public :: num2stri
   public :: num2strr
   public :: num2strc
   public :: num2strd
   public :: num2strz

   interface num2string
      module procedure num2stri, num2strd, num2strz, num2strr, num2strc
   end interface

contains

!****** Function: num2string **************************************************
!******************************************************************************
!
!  Converts a number to a character string. Integer, real and complex version
! are available (in single and double precision for real and complex). The
! number of decimal digits can be chosen for real and complex conversions, the
! precision should be between 1 and 15.
!
! INPUT -----------------------------------------------------------------------
!
! * num                           : The number to be converted
! integer dec                     : Decimal precision for real and complex
!                                     numbers
!
!******************************************************************************
elemental function num2stri(num) result(string)

   character(len=strlen) :: string
   integer, intent(in) :: num

   write(string,'(i0)') num

end function num2stri
!******************************************************************************
elemental function num2strr(num,dec) result(string)

   use precision,            only : sp

   character(len=strlen) :: string
   real(sp), intent(in) :: num
   integer, intent(in), optional :: dec

   character(len=strlen) :: s
   integer :: i, prec
   integer*8 :: n

   if (present(dec)) then
      prec = dec
   else
      prec = 4
   end if
   if (prec < 1) then
      prec = 1
   elseif (prec > 15) then
      prec = 15
   end if
   n = nint(num*10**prec)
   write(s,'(i0)') n
   i = len_trim(s)
   string = s(1:i-prec)//'.'//s(i-prec+1:i)

end function num2strr
!******************************************************************************
elemental function num2strc(num,dec) result(string)

   use precision,            only : sp

   character(len=strlen*2) :: string
   complex(sp), intent(in) :: num
   integer, intent(in), optional :: dec

   string = trim(num2string(real(num,sp),dec))//'+i'// &
     trim(num2string(aimag(num),dec))

end function num2strc
!******************************************************************************
elemental function num2strd(num,dec) result(string)

   use precision,            only : dp

   character(len=strlen) :: string
   real(dp), intent(in) :: num
   integer, intent(in), optional :: dec

   character(len=strlen) :: s
   integer :: i, prec
   real(dp) :: n

   if (present(dec)) then
      prec = dec
   else
      prec = 4
   end if
   if (prec < 1) then
      prec = 1
   elseif (prec > 15) then
      prec = 15
   end if
   write(s,'(f25.16)') num
   s = adjustl(s)
   i = index(s,'.')
   i = i+prec
   string = s(1:i)
   if (ichar(s(i+1:i+1))>53) then
      n = num + 5*10**(-prec-1)
      write(s,'(f25.16)') n
      s = adjustl(s)
      i = index(s,'.')
      i = i+prec
      string = s(1:i)
   end if

end function num2strd
!******************************************************************************
elemental function num2strz(num,dec) result(string)

   use precision,            only : dp


   character(len=strlen*2) :: string
   complex(dp), intent(in) :: num
   integer, intent(in), optional :: dec

   string = trim(num2string(real(num,dp),dec))//'+i'// &
     trim(num2string(aimag(num),dec))

end function num2strz
!****** End function: num2string **********************************************
!******************************************************************************


!****** Function: num2month ***************************************************
!******************************************************************************
!
!  Converts a number to a string with the corresponding month
!
! INPUT -----------------------------------------------------------------------
!
! integer num                     : The number of month
!
!******************************************************************************
elemental function num2month(num) result(string)

   character(9), parameter :: month(12) = &
                              (/'January  ','February ','March    ', &
                                'April    ','May      ', 'June     ', &
                                'July     ','August   ','September', &
                                'October  ','November ','December '/)
   integer, intent(in) :: num
   character(len=9) :: string

   if (num > 0 .and. num < 13) then
      string = month(num)
   else
      string = '*'
   end if

end function num2month
!****** End function: num2month ***********************************************
!******************************************************************************


!****** Function: num2bytes ***************************************************
!******************************************************************************
!
!  Converts a number of bytes to a string with the proper representation
!
! INPUT -----------------------------------------------------------------------
!
! integer num                     : Number of bytes
!
!******************************************************************************
elemental function num2bytes(num) result(string)

   use precision,            only : sp

   integer, intent(in) :: num
   character(len=30) :: string

   character(len=*), parameter :: u(5) = (/'bytes','kB   ','MB   ','GB   ','TB   '/)
   real(sp) :: q
   integer :: i

   q = real(num)
   do i=1,5
      if (q < 1000.0 .or. i>=5) then
         if (i==1) then
            string = num2stri(int(q))
         else
            string = num2strr(q,2)
         end if
         string = trim(string)//' '//u(i)
         exit
      end if
      q = q/1024.0
   end do

end function num2bytes
!****** End function: num2month ***********************************************
!******************************************************************************


!****** Function: lowercase ***************************************************
!******************************************************************************
!
!  Converts a string to lowercase only
!
! INPUT -----------------------------------------------------------------------
!
! character(*) str                : String to convert
!
! OUTPUT ----------------------------------------------------------------------
!
! character(*) lcstr              : Converted string
!
!******************************************************************************
elemental function lowercase(str) result(lcstr)

   use io,                   only : maxlinel

   integer, parameter :: idiff=ichar('Z')-ichar('z')

   character(*), intent(in) :: str
   character(maxlinel) :: lcstr

   integer :: l, i
   character :: s

   l=len_trim(str)
   lcstr = ''
   do i=1,l
      s = str(i:i)
      if (lge(s,'A') .and. lle(s,'Z')) s=char(ichar(s)-idiff)
      lcstr(i:i) = s
   end do

end function lowercase
!****** End function: lowercase ***********************************************
!******************************************************************************


!****** Function: uppercase ***************************************************
!******************************************************************************
!
!  Converts a string to uppercase only
!
! INPUT -----------------------------------------------------------------------
!
! character(*) str                : String to convert
!
! OUTPUT ----------------------------------------------------------------------
!
! character(*) lcstr              : Converted string
!
!******************************************************************************
elemental function uppercase(str) result(lcstr)

   use io,                   only : maxlinel

   integer, parameter :: idiff=ichar('Z')-ichar('z')

   character(*), intent(in) :: str
   character(maxlinel) :: lcstr

   integer :: l, i
   character :: s

   l=len_trim(str)
   lcstr = ''
   do i=1,l
      s = str(i:i)
      if (lge(s,'a') .and. lle(s,'z')) s=char(ichar(s)+idiff)
      lcstr(i:i) = s
   end do

end function uppercase
!****** End function: uppercase ***********************************************
!******************************************************************************

end module format
