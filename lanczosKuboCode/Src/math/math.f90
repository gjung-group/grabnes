!---------------------------------------------------------!
!                  ADVENTURE                              !
!                                                         !
!     Copyright (C) 2012 Rafael Martinez-Gordillo         !
!                                                         !
!   This file is distributed under the terms of the       !
!   GNU General Public License. See the file `LICENSE'    !
!   in the root directory of the present distribution,    !
!   or http://www.gnu.org/copyleft/gpl.txt                !
!                                                         !
!---------------------------------------------------------!

module math

   implicit none

interface
   subroutine AdamsIn(f,df,d2f,c1,n,mn,h,k)
   include 'parameters.h'
   real(dp), intent(inout) :: f(n), df(n), d2f(n)
   real(dp), intent(in) :: h,c1(n)
   integer, intent(in) :: n, mn, k
   end subroutine
end interface

interface
   subroutine AdamsOut(f,df,d2f,c1,n,mn,h,k,nodes)
   include 'parameters.h'
   real(dp), intent(inout) :: f(n), df(n), d2f(n)
   real(dp), intent(in) :: h, c1(n)
   integer, intent(in) :: n, mn, k
   integer, intent(out) :: nodes
   end subroutine
end interface

interface
   function AdamsOE(df,n,i,k) result(f)
   include 'parameters.h'
   integer, intent(in) :: n, i, k
   real(dp), intent(in) :: df(n)
   real(dp) :: f
   end function
end interface

interface
   function AdamsOI(df,n,i,k) result(f)
   include 'parameters.h'
   integer, intent(in) :: n, i, k
   real(dp), intent(in) :: df(n)
   real(dp) :: f
   end function
end interface

interface
   function AdamsIE(df,n,i,k) result(f)
   include 'parameters.h'
   integer, intent(in) :: n, i, k
   real(dp), intent(in) :: df(n)
   real(dp) :: f
   end function
end interface

interface
   function AdamsII(df,n,i,k) result(f)
   include 'parameters.h'
   integer, intent(in) :: n, i, k
   real(dp), intent(in) :: df(n)
   real(dp) :: f
   end function
end interface

interface
   subroutine AdamsOutX(f,df,d2f,c1,c2,n,mn,h,k,nodes)
   include 'parameters.h'
   real(dp), intent(inout) :: f(n), df(n), d2f(n)
   real(dp), intent(in) :: h,c1(n), c2(n)
   integer, intent(in) :: n, mn, k
   integer, intent(out) :: nodes
   end subroutine
end interface

interface
   function CrossProd(a,b) result(axb)
   include 'parameters.h'
   real(dp), intent(in) :: a(3), b(3)
   real(dp) :: axb(3)
   end function
end interface

interface
   subroutine LagrangeDiff(f,df,k)
   include 'parameters.h'
   real(dp), intent(in) :: f(k+1)
   real(dp), intent(out) :: df(k+1)
   integer, intent(in) :: k
   end subroutine
end interface

interface
   subroutine LagrangeDiffX(f,df,E,V,r,dr,Z,l,h,k)
   include 'parameters.h'
   real(dp), intent(in) :: V(0:k), r(0:k), dr(0:k), E, h
   real(dp), intent(out) :: f(0:k), df(0:k)
   integer, intent(in) :: k, l, Z
   end subroutine
end interface

interface
   function norm(a)
   include 'parameters.h'
   real(dp), intent(in) :: a(*)
   real(dp) :: norm
   end function
end interface

interface
   function NumerovIn(f,n,i,h) result(y)
   include 'parameters.h'
   integer, intent(in) :: n, i
   real(dp), intent(in) :: f(n)
   real(dp), intent(in) :: h
   real(dp) :: y
   end function
end interface

interface
   function TrapezoidalInt(f,n,h,k) result(s)
   include 'parameters.h'
   integer, intent(in) :: n, k
   real(dp), intent(in) :: f(n), h
   real(dp) :: s
   end function
end interface

end module math
