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

module string

   use io,                   only : maxlinel

   implicit none

   PRIVATE

   public :: StringComp

   interface StringComp
      module procedure StringComp_s
      module procedure StringComp_v
   end interface

contains

!****** Function: StringComp **************************************************
!******************************************************************************
!
!  Compare two strings and check whether they are equal
!
! INPUT -----------------------------------------------------------------------
!
! character string1, string2     : Strings to be compared
!
! OUTPUT ----------------------------------------------------------------------
!
! logical equal                  : True if the strings are equal
!
!******************************************************************************
function StringComp_s(string1,string2) result(equal)

   implicit none

   character(*), intent(in) :: string1, string2
   logical :: equal

   integer :: len1, len2
   integer :: idiff, il
   character :: s1, s2
   character(len=maxlinel) :: str1, str2
   logical :: ignore

   str1 = string1
   ignore = .true.
   do il=1,len(string1)
      if (str1(il:il) == ' ') then
         if (str1(1:1)=='&' .and. ignore) then
            ignore = .false.
         else
            str1 = str1(1:il-1)
            exit
         end if
      else if (str1(il:il) == '=' .or. str1(il:il) == ':') then
         str1 = str1(1:il-1)
         exit
      else if (str1(il:il)=='_' .or. str1(il:il)=='-') then
         str1 = str1(1:il-1)//str1(il+1:)
      end if
   end do
   str2 = string2
   ignore = .true.
   do il=1,len(string2)
      if (str2(il:il) == ' ') then
         if (str2(1:1)=='&' .and. ignore) then
            ignore = .false.
         else
            str2 = str2(1:il-1)
            exit
         end if
      else if (str2(il:il) == '=' .or. str2(il:il) == ':') then
         str2 = str2(1:il-1)
         exit
      else if (str2(il:il)=='_' .or. str2(il:il)=='-') then
         str2 = str2(1:il-1)//str2(il+1:)
      end if
   end do
   len1 = len_trim(str1)
   len2 = len_trim(str2)
   equal = .false.
   if (len1 /= len2) then
      return
   end if
   idiff = ichar('Z')-ichar('z')
   do il=1,len1
      s1 = str1(il:il)
      s2 = str2(il:il)
      if (lge(s1,'A') .and. lle(s1,'Z')) s1=char(ichar(s1)-idiff)
      if (lge(s2,'A') .and. lle(s2,'Z')) s2=char(ichar(s2)-idiff)
      if (s1 .ne. s2) return
   end do
   equal = .true.

 end function StringComp_s
!******************************************************************************
function StringComp_v(string1,string2) result(equal)

   implicit none

   character(*), intent(in) :: string1(:), string2
   logical :: equal

   integer :: len1, len2
   integer :: idiff, il, istr
   character :: s1, s2
   character(len=maxlinel) :: str1, str2
   logical :: ignore

   str2 = string2
   ignore = .true.
   do il=1,len(string2)
      if (str2(il:il) == ' ') then
         if (str2(1:1)=='&' .and. ignore) then
            ignore = .false.
         else
            str2 = str2(1:il-1)
            exit
         end if
      else if (str2(il:il) == '=' .or. str2(il:il) == ':') then
         str2 = str2(1:il-1)
         exit
      else if (str2(il:il)=='_' .or. str2(il:il)=='-') then
         str2 = str2(1:il-1)//str2(il+1:)
      end if
   end do
   len2 = len_trim(str2)
   equal = .false.
cy:do istr=1,size(string1)
      str1 = string1(istr)
      ignore = .true.
      do il=1,len_trim(string1(istr))
         if (str1(il:il) == ' ') then
            if (str1(1:1)=='&' .and. ignore) then
               ignore = .false.
            else
               str1 = str1(1:il-1)
               exit
            end if
         else if (str1(il:il) == '=' .or. str1(il:il) == ':') then
            str1 = str1(1:il-1)
            exit
         else if (str1(il:il)=='_' .or. str1(il:il)=='-') then
            str1 = str1(1:il-1)//str1(il+1:)
         end if
      end do
      len1 = len_trim(str1)
      if (len1 /= len2) then
         cycle
      end if
      idiff = ichar('Z')-ichar('z')
        do il=1,len1
         s1 = str1(il:il)
         s2 = str2(il:il)
         if (lge(s1,'A') .and. lle(s1,'Z')) s1=char(ichar(s1)-idiff)
         if (lge(s2,'A') .and. lle(s2,'Z')) s2=char(ichar(s2)-idiff)
         if (s1 .ne. s2) cycle cy
      end do
      equal = .true.
      return
   end do cy

end function StringComp_v
!****** End function: StringComp **********************************************
!******************************************************************************


end module string
