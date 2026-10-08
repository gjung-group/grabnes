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

module parser

   implicit none

   PRIVATE

   public :: ParserCheck
   public :: ParserFindValue

   interface ParserCheck
      module procedure ParserCheck_i, ParserCheck_r, ParserCheck_d
      module procedure ParserCheck_s, ParserCheck_l
   end interface

   interface ParserFindValue
      module procedure ParserFindValue_i, ParserFindValue_d
   end interface

contains

!****** Subroutine: ParserCheck ***********************************************
!******************************************************************************
!
!  Check if some value within a string is correct for its type
!
! INPUT -----------------------------------------------------------------------
!
! character string               : String where to look for
! integer enum                   : Number of the element inside string
!
! OUTPUT ----------------------------------------------------------------------
!
! * value                        : Value obtained from the string
! integer err                    : Number of error
!
!******************************************************************************
subroutine ParserCheck_i(value,string,enum,err)

   use sys,                  only : SysIOErr
   use io,                   only : maxlinel

   implicit none

   integer, intent(out) :: value
   character(*), intent(in) :: string
   integer, intent(in), optional :: enum
   integer, intent(out), optional :: err

   integer :: en, e
   character(len=maxlinel) :: str

   e = 0
   if (present(enum)) then
      en = enum
   else
      en = 1
   end if
   if (ParserGetElement(string,str,en)) then
      call ParserInt(str,value,e)
   else
      e = -1
   end if

   if (present(err)) then
      err = e
   else
      call SysIOErr(e)
   end if

end subroutine ParserCheck_i
!******************************************************************************
subroutine ParserCheck_r(value,string,enum,err)

   use sys,                  only : SysIOErr
   use precision,            only : sp


   implicit none

   real(sp), intent(out) :: value
   character(*), intent(in) :: string
   integer, intent(in), optional :: enum
   integer, intent(out), optional :: err

   integer :: en, e
   character(len=50) :: str

   e = 0
   if (present(enum)) then
      en = enum
   else
      en = 1
   end if
   if (ParserGetElement(string,str,en)) then
      call ParserR(str,value,e)
   else
      e = -1
   end if

   if (present(err)) then
      err = e
   else
      call SysIOErr(e)
   end if

end subroutine ParserCheck_r
!******************************************************************************
subroutine ParserCheck_d(value,string,enum,err)

   use sys,                  only : SysIOErr
   use precision,            only : dp
   use io,                   only : maxlinel

   implicit none

   real(dp), intent(out) :: value
   character(*), intent(in) :: string
   integer, intent(in), optional :: enum
   integer, intent(out), optional :: err

   integer :: en, e
   character(len=maxlinel) :: str

   e = 0
   if (present(enum)) then
      en = enum
   else
      en = 1
   end if
   if (ParserGetElement(string,str,en)) then
      call ParserD(str,value,e)
   else
      e = -1
   end if

   if (present(err)) then
      err = e
   else
      call SysIOErr(e)
   end if

end subroutine ParserCheck_d
!******************************************************************************
subroutine ParserCheck_s(value,string,enum,err,multiword)

   use sys,                  only : SysIOErr
   use io,                   only : maxlinel

   implicit none

   character(*), intent(out) :: value
   character(*), intent(in) :: string
   integer, intent(in), optional :: enum
   integer, intent(out), optional :: err
   logical, intent(in), optional :: multiword

   integer :: en, e
   character(len=maxlinel) :: str

   e = 0
   if (present(enum)) then
      en = enum
   else
      en = 1
   end if
   if (ParserGetElement(string,str,en,multiword)) then
      value = str
   else
      e = -1
   end if

   if (present(err)) then
      err = e
   else
      call SysIOErr(e)
   end if

end subroutine ParserCheck_s
!******************************************************************************
subroutine ParserCheck_l(value,string,enum,err)

   use sys,                  only : SysIOErr
   use io,                   only : maxlinel

   implicit none

   logical, intent(out) :: value
   character(*), intent(in) :: string
   integer, intent(in), optional :: enum
   integer, intent(out), optional :: err

   integer :: en, e
   character(len=maxlinel) :: str

   e = 0
   if (present(enum)) then
      en = enum
   else
      en = 1
   end if
   if (ParserGetElement(string,str,en)) then
      call ParserL(str,value,e)
   else
      e = -1
   end if

   if (present(err)) then
      err = e
   else
      call SysIOErr(e)
   end if

end subroutine ParserCheck_l
!****** End subroutine: ParserCheck *******************************************
!******************************************************************************

!****** Subroutine: ParserGetElement ******************************************
!******************************************************************************
!
!  Given a string, get the n-th element. Each element is delimited by spaces or
!   an equal sign.
!
! INPUT -----------------------------------------------------------------------
!
! character string               : String where to look for
! integer n                      : Number of the element inside string
! logical multiword              : If true, this function returns in element
!                                  a string including the n-th element up to
!                                  the end of string
!
! OUTPUT ----------------------------------------------------------------------
!
! character element              : String of that element
! logical success                : If true, the n-th element of string was
!                                  found
!
!******************************************************************************
 function ParserGetElement(string,element,n,multiword) result(success)

   implicit none

   character(*), intent(in) :: string
   character(*), intent(out) :: element
   integer, intent(in) :: n
   logical, intent(in), optional :: multiword
   logical :: success

   integer :: i1, i2, il, ie, lstr
   character :: s
   logical :: inStr, inSp, mw

   if (present(multiword)) then
      mw = multiword
   else
      mw = .false.
   end if
   success = .true.
   i1 = 1
   i2 = 1
   ie = 0
   inStr = .false.
   inSp = .true.
   lstr = len_trim(string)
   do il = 1, lstr
      if (inStr) then
         if (string(il:il) == s) then
            inStr = .false.
         end if
         cycle
      else if (string(il:il) == "'") then
         inStr = .true.
         s = "'"
      else if (string(il:il) == '"') then
         inStr = .true.
         s = '"'
      else if (string(il:il) == ' ' .and. .not. inSp) then
         ie = ie+1
         inSp = .true.
         if (ie == n) then
            exit
         else
            i1 = il
         end if
      else if (il==lstr) then
         ie = ie+1
      else
         inSp = .false.
      end if
      i2 = il
   end do
   if (mw) then
      if (ie < n-1) then
         success = .false.
      else
         element = adjustl(string(i1:))
      end if
   else
      if (ie < n) then
         success = .false.
      else
         element = adjustl(string(i1:i2))
      end if
   end if

end function ParserGetElement
!****** End subroutine: ParserGetElement **************************************
!******************************************************************************

!****** Subroutine: ParserInt *************************************************
!******************************************************************************
!
!  Traduce an integer expression to a number
!
! INPUT -----------------------------------------------------------------------
!
! character string               : String with the integer expression
!
! OUTPUT ----------------------------------------------------------------------
!
! integer value                  : Value obtained from the string
! integer error                  : Value of the error
!
!******************************************************************************
subroutine ParserInt(string,value,error)

   use sys,                  only : SysIOErr
   use precision,            only : dp
   use format,               only : lowercase

   implicit none

   integer, intent(out) :: value
   character(*), intent(in) :: string
   integer, intent(out), optional :: error

   integer :: err, l
   real(dp) :: v

   err = 0
   l = len_trim(string)
   if (l==0) then
      err = 13
   else
      call ParserProcess(lowercase(string),v,err)
      if (mod(v,1.0_dp) /= 0.0_dp .and. err==0) err = -2
      value = int(v)
   end if
   if (present(error)) then
      error = err
   else
      call SysIOErr(err,'parser','ParserInt')
   end if

end subroutine ParserInt
!****** End subroutine: ParserInt *********************************************
!******************************************************************************


!****** Subroutine: ParserD ***************************************************
!******************************************************************************
!
!  Traduce a double expression to a number
!
! INPUT -----------------------------------------------------------------------
!
! character string               : String with the integer expression
!
! OUTPUT ----------------------------------------------------------------------
!
! real value                     : Value obtained from the string
! integer error                  : Value of the error
!
!******************************************************************************
subroutine ParserR(string,value,error)

   use sys,                  only : SysIOErr
   use precision,            only : sp, dp
   use format,               only : lowercase

   implicit none

   real(sp), intent(out) :: value
   character(*), intent(in) :: string
   integer, intent(out), optional :: error

   integer :: err, l
   real(dp) :: v

   err = 0
   l = len_trim(string)
   if (l==0) then
      err = 13
   else
      call ParserProcess(lowercase(string),v,err)
      value = v
   end if
   if (present(error)) then
      error = err
   else
      call SysIOErr(err,'parser','ParserInt')
   end if

end subroutine ParserR
!******************************************************************************
subroutine ParserD(string,value,error)

   use sys,                  only : SysIOErr
   use precision,            only : dp
   use format,               only : lowercase

   implicit none

   real(dp), intent(out) :: value
   character(*), intent(in) :: string
   integer, intent(out), optional :: error

   integer :: err, l

   err = 0
   l = len_trim(string)
   if (l==0) then
      err = 13
   else
      call ParserProcess(lowercase(string),value,err)
   end if
   if (present(error)) then
      error = err
   else
      call SysIOErr(err,'parser','ParserInt')
   end if

end subroutine ParserD
!****** End subroutine: ParserD ***********************************************
!******************************************************************************

!****** Subroutine: ParserL ***************************************************
!******************************************************************************
!
!  Traduce an integer expression to a number
!
! INPUT -----------------------------------------------------------------------
!
! character string               : String with the integer expression
!
! OUTPUT ----------------------------------------------------------------------
!
! logical value                  : Value obtained from the string
! integer error                  : Value of the error
!
!******************************************************************************
subroutine ParserL(string,value,error)

   use sys,                  only : SysIOErr
   use precision,            only : dp
   use format,               only : lowercase
   use io,                   only : maxlinel

   implicit none

   logical, intent(out) :: value
   character(*), intent(in) :: string
   integer, intent(out), optional :: error

   integer :: err, l
   character(len=maxlinel) :: str

   err = 0
   l = len_trim(string)
   if (l==0) then
      err = 13
   else
      str = lowercase(string)
      if (str=='.true.' .or. str=='true' .or. str=='t' .or. str=='yes' &
        .or. str=='y') then
         value = .true.
      else if (str=='.false.' .or. str=='false' .or. str=='f' .or. str=='no' &
        .or. str=='n') then
         value = .false.
      else
         err = -3
      end if
   end if
   if (present(error)) then
      error = err
   else
      call SysIOErr(err,'parser','ParserInt')
   end if

end subroutine ParserL
!****** End subroutine: ParserL ***********************************************
!******************************************************************************

!****** Subroutine: ParserProcess *********************************************
!******************************************************************************
!
!  Process an expression an returns a valid number
!
! INPUT -----------------------------------------------------------------------
!
! character string               : String with the expression
!
! OUTPUT ----------------------------------------------------------------------
!
! real value                     : Value obtained from the string
! integer error                  : Value of the error
!
!******************************************************************************
recursive subroutine ParserProcess(string,value,error,reset)

   use precision,            only : dp
   use sys,                  only : SysIOErr, SysKill
   use io,                   only : maxlinel

   implicit none

   integer, parameter :: nvalues = 64
   character(*), parameter :: op(8) = (/'acos','asin','atan','cos ','sin ','tan '&
                                    ,'sqrt','exp '/)
   character(*), parameter :: arit(5) = (/'**','/ ','* ','- ','+ '/)
   character(*), parameter :: dec(2) = (/'e','d'/)
   integer, parameter :: nop = size(op), ndec = size(dec), narit = size(arit)

   real(dp), intent(out) :: value
   character(*), intent(in) :: string
   integer, intent(out), optional :: error
   logical, intent(in), optional :: reset

   integer :: err
   integer :: indx1, indx2, indx, iop, lop, i
   character(len=maxlinel) :: str
   real(dp), save :: v(nvalues)
   integer, save :: iv = 0
   character(len=32) :: s
   integer :: iop2, indx0, ii
   real(dp) :: val1, val2

   if (present(reset)) then
      if (reset) then
         iv = 0
      end if
   else
      iv = 0
   end if
   err = 0
   str = string
   ! Process expressions with functions
cy:do iop=1,nop
      indx = 1
      lop = len_trim(op(iop))
      do
         indx1 = index(str(indx:),trim(op(iop)))
         if (indx1/=0) then
            indx1 = indx1 + indx + lop - 1
            call ParserParentheses(str,indx1,indx2)
            call ParserProcess(str(indx1+1:indx2-1),val1,err,.false.)
            if (err/=0) exit cy
            iv = iv+1
            if (iv > nvalues) then
                call SysKill("Cannot parse expression. Try increasing "// &
                  "'nvalues' in parser.f90 and recompile",'parser','ParserProcess')
            end if
            v(iv) = ParserOperation(val1,iop)
            write(s,'(a,i0,a)') '@',iv,'@'
            indx1 = indx1 - lop
            str = str(:indx1-1)//trim(s)//str(indx2+1:)
            indx = indx1 + len_trim(s)
         else
            exit
         end if
      end do
   end do cy

   if (err/=0) then
      if (present(error)) then
         error = err
         return
      else
         call SysIOErr(err,'parser','ParserProcess')
      end if
   end if

   ! Process expressions in parenthesis
   indx = 1
   do
      indx1 = index(str(indx:),'(')
      if (indx1 /= 0) then
         indx1 = indx1 + indx - 1
         call ParserParentheses(str,indx1,indx2)
         call ParserProcess(str(indx1+1:indx2-1),val1,err,.false.)
         if (err/=0) exit
         iv = iv+1
         if (iv > nvalues) then
            call SysKill("Cannot parse expression. Try increasing "// &
              "'nvalues' in parser.f90 and recompile",'parser','ParserProcess')
         end if
         v(iv) = val1
         write(s,'(a,i0,a)') '@',iv,'@'
         str = str(:indx1-1)//trim(s)//str(indx2+1:)
         indx = indx1 + len_trim(s)
      else
         exit
      end if
   end do

   if (err/=0) then
      if (present(error)) then
         error = err
         return
      else
         call SysIOErr(err,'parser','ParserProcess')
      end if
   end if

   ! Process numbers with exponential notation
l1:do iop=1,ndec
      indx = 1
      do
         indx1 = index(str(indx:),dec(iop))
         if (indx1/=0) then
            indx = indx1 + indx - 1
            do i=indx-1,1,-1
               if (lge(str(i:i),'0') .and. lle(str(i:i),'9')) then
                  indx1 = i
               else if (str(i:i)=='.') then
                  indx1 = i
               else
                  exit
               end if
            end do
            if (str(indx+1:indx+1)=='+' .or. str(indx+1:indx+1)=='-') indx = indx + 1
            indx2 = len_trim(str)+1
            do i=indx+1,len_trim(str)
               if (lge(str(i:i),'0') .and. lle(str(i:i),'9')) then
                  indx2 = i
               else
                  exit
               end if
            end do
            iv = iv+1
            if (iv > nvalues) then
                call SysKill("Cannot parse expression. Try increasing "// &
                "'nvalues' in parser.f90 and recompile",'parser','ParserProcess')
            end if
            read(str(indx1:indx2),*,IOSTAT=err) v(iv)
            if (err/=0) then
               err=15
               exit l1
            end if
            write(s,'(a,i0,a)') '@',iv,'@'
            str = str(:indx1-1)//trim(s)//str(indx2+1:)
            indx = indx1 + len_trim(s)
         else
            exit
         end if
      end do
   end do l1

   if (err/=0) then
      if (present(error)) then
         error = err
         return
      else
         call SysIOErr(err,'parser','ParserProcess')
      end if
   end if

   do i=1,len_trim(str)
      if (llt(str(i:i),'*') .or. lgt(str(i:i),'9')) then
         if (str(i:i)/='@') then
            err = 16
            if (present(error)) then
               error = err
               return
            else
               call SysIOErr(err,'parser','ParserProcess')
            end if
         end if
      end if
   end do

l2:do iop = 1,narit
      indx = 1
      do
         indx1 = index(str(indx:),trim(arit(iop)))
         if (indx1/=0) then
            indx1 = indx1 + indx - 1
            indx0 = 0
            do iop2=1,narit
               indx2 = index(str(1:indx1-1),trim(arit(iop2)),back=.true.)
               indx0 = max(indx0,indx2)
            end do
            s = str(indx0+1:indx1-1)
            ii = indx0
            if (s(1:1)=='@') then
               indx0 = index(s(2:),'@')
               if (indx0==0) then
                  call SysKill('Wrong format for saved values','parser', &
                    'ParserProcess')
               end if
               read(s(2:indx0),*) indx2
               val1 = v(indx2)
            else
               read(s,*,IOSTAT=err) val1
               if (err/=0) then
                  err=15
                  exit l2
               end if
            end if
            indx0 = len_trim(str)+1
            do iop2=1,narit
               indx2 = index(str(indx1+len_trim(arit(iop)):),trim(arit(iop2)))
               if (indx2/=0) indx0 = min(indx0,indx2 + indx1 + len_trim(arit(iop)) - 1)
            end do
            s = str(indx1+len_trim(arit(iop)):indx0-1)
            if (s(1:1)=='@') then
               indx1 = index(s(2:),'@')
               if (indx1==0) then
                  call SysKill('Wrong format for saved values','parser', &
                    'ParserProcess')
               end if
               read(s(2:indx1),*) indx2
               val2 = v(indx2)
            else
               read(s,*,IOSTAT=err) val2
               if (err/=0) then
                  err=15
                  exit l2
               end if
            end if
            iv = iv+1
            if (iv > nvalues) then
               call SysKill("Cannot parse expression. Try increasing "// &
                 "'nvalues' in parser.f90 and recompile",'parser','ParserProcess')
            end if
            select case(iop)
            case(1)
               v(iv) = val1**val2
            case(2)
               v(iv) = val1/val2
            case(3)
               v(iv) = val1*val2
            case(4)
               v(iv) = val1-val2
            case(5)
               v(iv) = val1+val2
            end select
            write(s,'(a,i0,a)') '@',iv,'@'
            str = str(:ii)//trim(s)//str(indx0:)
            indx = ii + len_trim(s) + 1
         else
            exit
         end if
      end do
   end do l2

   if (err/=0) then
      if (present(error)) then
         error = err
         return
      else
         call SysIOErr(err,'parser','ParserProcess')
      end if
   end if

   if (str(1:1)=='@') then
      indx0 = index(str(2:),'@')
      if (indx0==0) then
         call SysKill('Wrong format for saved values','parser','ParserProcess')
      end if
      s = str(2:indx0)
      read(s,*) indx
      value = v(indx)
   else
      read(str,*,IOSTAT=err) value
      if (err/=0) then
         err=15
      end if
   end if

   if (err/=0) then
      if (present(error)) then
         error = err
         return
      else
         call SysIOErr(err,'parser','ParserProcess')
      end if
   end if

   ! error is intent(out): it must also be defined on success
   if (present(error)) error = 0

end subroutine ParserProcess
!****** End subroutine: ParserProcess *****************************************
!******************************************************************************


!****** Subroutine: ParserParentheses *****************************************
!******************************************************************************
!
!  Finds closing parentheses that correspond to the one in position i1 inside
! string
!
! INPUT -----------------------------------------------------------------------
!
! character string               : String where to find parentheses
! integer i1                     : Position of opening parentheses
!
! OUTPUT ----------------------------------------------------------------------
!
! integer i2                     : Position of closing parentheses
! integer error                  : Value of the error
!
!******************************************************************************
subroutine ParserParentheses(string,i1,i2,error)

   use sys,                  only : SysIOErr

   implicit none

   character(*), intent(in) :: string
   integer, intent(inout) :: i1
   integer, intent(out) :: i2
   integer, intent(out), optional :: error

   integer :: l, i, np
   integer :: err

   l = len_trim(string)
   np = 1
   err = 0
   i2 = 0
   if (string(i1:i1)/='(') then
      i1 = l+1
      err = 14
   end if
   do i=i1+1,l
      if (string(i:i)==')') then
         np = np - 1
         if (np==0) then
            i2 = i
            exit
         end if
      end if
      if (string(i:i)=='(') then
         np = np + 1
      end if
   end do
   if (i2==0 .and. err==0) err = 12
   if (present(error)) then
      error = err
   else
      call SysIOErr(err,'parser','ParserParentheses')
   end if

end subroutine ParserParentheses
!****** End subroutine: ParserParentheses *************************************
!******************************************************************************


!****** Subroutine: ParserOperation *******************************************
!******************************************************************************
!
!  Performs the specified operation to value
!
! INPUT -----------------------------------------------------------------------
!
! real value                     : Value to be inserted in the function
! integer n                      : Number of operation
!
! OUTPUT ----------------------------------------------------------------------
!
! real res                      : Result
!
!******************************************************************************
function ParserOperation(value,n) result(res)

   use precision,            only : dp
   use sys,                  only : SysKill

   implicit none

   real(dp), intent(in) :: value
   integer, intent(in) :: n
   real(dp) :: res

   if (n==1) then
      res = acos(value)
   else if (n==2) then
      res = asin(value)
   else if (n==3) then
      res = atan(value)
   else if (n==4) then
      res = cos(value)
   else if (n==5) then
      res = sin(value)
   else if (n==6) then
      res = tan(value)
   else if (n==7) then
      res = sqrt(value)
   else if (n==8) then
      res = exp(value)
   else
      call SysKill('Undefined parser operation','parser','ParserOperation')
   end if

end function ParserOperation
!****** End subroutine: ParserOperation ***************************************
!******************************************************************************

!****** Subroutine: ParserFindValue *******************************************
!******************************************************************************
!
!  Find a value in a string. The expression is of the form: 'name = value'
!
! INPUT -----------------------------------------------------------------------
!
! character name                 : Name of the value to find
! character string               : String where to look for the value
!
! OUTPUT ----------------------------------------------------------------------
!
! * value                       : The value found
! logical found                 : True if the value is found
!
!******************************************************************************
function ParserFindValue_i(name,strng,value,err) result(found)

   use sys,                  only : SysKill
   use io,                   only : maxlinel
   use string,               only : StringComp

   implicit none

   character(*), intent(in) :: strng, name
   integer, intent(out) :: value, err
   logical :: found

   integer :: i, sz
   character(len=maxlinel) :: str

   err = 0
   i = 0
   do
      i = i+1
      if (ParserGetElement(strng,str,i,.true.)) then
         if (StringComp(str,name)) then
            sz = len_trim(name)
            str = adjustl(str(sz+1:))
            found = .true.
            if (str(1:1)=='=') then
               str = adjustl(str(2:))
               call ParserCheck(value,str,1,err)
            else
               err = -4
            end if
            exit
         end if
      else
         found = .false.
         exit
      end if
   end do

end function ParserFindValue_i
!******************************************************************************
function ParserFindValue_d(name,strng,value,err) result(found)

   use sys,                  only : SysKill
   use io,                   only : maxlinel
   use string,               only : StringComp
   use precision,            only : dp

   implicit none

   character(*), intent(in) :: strng, name
   integer, intent(out) :: err
   real(dp), intent(out) :: value
   logical :: found

   integer :: i, sz
   character(len=maxlinel) :: str

   err = 0
   i = 0
   do
      i = i+1
      if (ParserGetElement(strng,str,i,.true.)) then
         if (StringComp(str,name)) then
            sz = len_trim(name)
            str = adjustl(str(sz+1:))
            found = .true.
            if (str(1:1)=='=') then
               str = adjustl(str(2:))
               call ParserCheck(value,str,1,err)
            else
               err = -4
            end if
            exit
         end if
      else
         found = .false.
         exit
      end if
   end do

end function ParserFindValue_d
!****** End subroutine: ParserFindValue ***************************************
!******************************************************************************

end module parser
