!---------------------------------------------------------!
!                 MATH LIBRARY                            !
!                                                         !
!     Copyright (C) 2012 Rafael Martinez-Gordillo         !
!                                                         !
!   This file is distributed under the terms of the       !
!   GNU General Public License. See the file `LICENSE'    !
!   in the root directory of the present distribution,    !
!   or http://www.gnu.org/copyleft/gpl.txt                !
!                                                         !
!---------------------------------------------------------!

!
! This routines are used for the integration of the atomic wavefunctions.
! References:
!   [1] Johnson, W. R.; 'Atomic Structure Theory: Lectures on Atomic Physics',
!       Springer:Berlin, 2007; Chapter 2, pp. 35-45
!   [2] Math Library-Source Code. 2004. RLH. 29 Jan. 2012.
!       <http://www.mymathlib.com/>.
!

!****** Subroutine: AdamsOut **************************************************
!******************************************************************************
!
!  Integrate the radial wavefunction outwards. The second derivative is
! calculated as:
!            d2f(i) = c(i)*f(i)
!
! INPUT -----------------------------------------------------------------------
!
! real f(n), df(n), d2f(n)        : The function and its derivatives
! real c(n)                       : Coefficients to calculate d2f
! integer n                       : Size of the function
! intger mn                       : Point up to which integrate
! real h                          : Mesh spacing
! integer k                       : The order of the polynomial
!
! OUTPUT ----------------------------------------------------------------------
!
! real f(n), df(n), d2f(n)        : The function and its derivates
! integer nodes                   : The number of nodes in the function
!
!******************************************************************************
subroutine AdamsOut(f,df,d2f,c1,n,mn,h,k,nodes)

   implicit none

   include 'parameters.h'

   real(dp), intent(inout) :: f(n), df(n), d2f(n)
   real(dp), intent(in) :: h, c1(n)
   integer, intent(in) :: n, mn, k
   integer, intent(out) :: nodes

   integer :: m, i, ii
   real(dp) :: AdamsOE, AdamsOI

   if (mn<n) then
      m = mn
   else
      m = n
   end if
   nodes = 0
   do i=k+1,m-1
      f(i+1) = f(i) + AdamsOE(df,n,i,k)
      df(i+1) = df(i) + AdamsOE(d2f,n,i,k)
      do ii=1,2
         d2f(i+1) = c1(i+1)*f(i+1) + df(i+1)*h
         df(i+1) = df(i) + AdamsOI(d2f,n,i,k)
         f(i+1) = f(i) + AdamsOI(df,n,i,k)
      end do
      if (f(i)*f(i+1) <= 0.0_dp) then
         nodes = nodes+1
      end if
   end do

end subroutine AdamsOut
!****** End subroutine: AdamsOut **********************************************
!******************************************************************************


!****** Subroutine: AdamsIn ***************************************************
!******************************************************************************
!
!  Integrate the radial wavefunction inwards. The second derivative is
! calculated as:
!            d2f(i) = c(i)*f(i)
!
! INPUT -----------------------------------------------------------------------
!
! real f(n), df(n), d2f(n)        : The function and its derivatives
! real c(n)                       : Coefficients to calculate d2f
! integer n                       : Size of the function
! intger mn                       : Point up to which integrate
! real h                          : Mesh spacing
! integer k                       : The order of the polynomial
!
! OUTPUT ----------------------------------------------------------------------
!
! real f(n), df(n), d2f(n)        : The function and its derivates
!
!******************************************************************************
subroutine AdamsIn(f,df,d2f,c1,n,mn,h,k)

   implicit none

   include 'parameters.h'

   real(dp), intent(inout) :: f(n), df(n), d2f(n)
   real(dp), intent(in) :: h,c1(n)
   integer, intent(in) :: n, mn, k

   integer :: i, ii
   real(dp) :: AdamsIE, AdamsII

   do i=n-k,mn+1,-1
      f(i-1) = f(i) + AdamsIE(df,n,i,k)
      df(i-1) = df(i) + AdamsIE(d2f,n,i,k)
      do ii=1,2
         d2f(i-1) = c1(i-1)*f(i-1) + df(i-1)*h
         df(i-1) = df(i) + AdamsII(d2f,n,i,k)
         f(i-1) = f(i) + AdamsII(df,n,i,k)
      end do
   end do

end subroutine AdamsIn
!****** End subroutine: AdamsIn ***********************************************
!******************************************************************************


!****** Function: Adams* ******************************************************
!******************************************************************************
!
!  Functions for outward/inwards extrapolation/interpolation.
!
! INPUT -----------------------------------------------------------------------
!
! real df(n)                      : Derivative of the function
! integer n                       : Size of df
! integer i                       : Index of the point to be calculated
! integer k                       : The order of the polynomial
!
! RESULT ----------------------------------------------------------------------
!
! real f                          : The function at point i
!
!******************************************************************************
function AdamsOE(df,n,i,k) result(f)

   implicit none

   include 'parameters.h'

   integer, intent(in) :: n, i, k
   real(dp), intent(in) :: df(n)
   real(dp) :: f

   if (k==3) then
      f = (55.0_dp*df(i)-59.0_dp*df(i-1)+37.0_dp*df(i-2)-9.0_dp*df(i-3))/24.0_dp
   elseif (k==4) then
      f = (190.1_dp*df(i)-277.4_dp*df(i-1)+261.6_dp*df(i-2)-127.4_dp*df(i-3)+25.1_dp*df(i-4))/72.0_dp
   elseif (k==5) then
      f = (427.7_dp*df(i)-792.3_dp*df(i-1)+998.2_dp*df(i-2)-729.8_dp*df(i-3)+287.7_dp*df(i-4)-47.5*df(i-5))/144.0_dp
   elseif (k==6) then
      f = (1987.21_dp*df(i)-4472.88_dp*df(i-1)+7055.49_dp*df(i-2)-6882.56_dp*df(i-3)+4071.39_dp*df(i-4)- &
        1344.72_dp*df(i-5)+190.87_dp*df(i-6))/604.8_dp
   end if

end function AdamsOE
!******************************************************************************
function AdamsOI(df,n,i,k) result(f)

   implicit none

   include 'parameters.h'

   integer, intent(in) :: n, i, k
   real(dp), intent(in) :: df(n)
   real(dp) :: f

   if (k==3) then
      f = (9.0_dp*df(i+1)+19.0_dp*df(i)-5.0_dp*df(i-1)+df(i-2))/24.0_dp
   elseif (k==4) then
      f = (25.1_dp*df(i+1)+64.6_dp*df(i)-26.4_dp*df(i-1)+10.6_dp*df(i-2)-1.9_dp*df(i-3))/72.0_dp
   elseif (k==5) then
      f = (47.5_dp*df(i+1)+142.7_dp*df(i)-79.8_dp*df(i-1)+48.2_dp*df(i-2)-17.3_dp*df(i-3)+2.7*df(i-4))/144.0_dp
   elseif (k==6) then
      f = (190.87_dp*df(i+1)+651.12_dp*df(i)-464.61_dp*df(i-1)+375.04_dp*df(i-2)-202.11_dp*df(i-3)+ &
        63.12_dp*df(i-4)-8.63_dp*df(i-5))/604.8_dp
   end if

end function AdamsOI
!******************************************************************************
function AdamsIE(df,n,i,k) result(f)

   implicit none

   include 'parameters.h'

   integer, intent(in) :: n, i, k
   real(dp), intent(in) :: df(n)
   real(dp) :: f

   if (k==3) then
      f = -(55.0_dp*df(i)-59.0_dp*df(i+1)+37.0_dp*df(i+2)-9.0_dp*df(i+3))/24.0_dp
   elseif (k==4) then
      f = -(190.1_dp*df(i)-277.4_dp*df(i+1)+261.6_dp*df(i+2)-127.4_dp*df(i+3)+25.1_dp*df(i+4))/72.0_dp
   elseif (k==5) then
      f = -(427.7_dp*df(i)-792.3_dp*df(i+1)+998.2_dp*df(i+2)-729.8_dp*df(i+3)+287.7_dp*df(i+4)-47.5*df(i+5))/144.0_dp
   elseif (k==6) then
      f = -(1987.21_dp*df(i)-4472.88_dp*df(i+1)+7055.49_dp*df(i+2)-6882.56_dp*df(i+3)+4071.39_dp*df(i+4)- &
        1344.72_dp*df(i+5)+190.87_dp*df(i+6))/604.8_dp
   end if

end function AdamsIE
!******************************************************************************
function AdamsII(df,n,i,k) result(f)

   implicit none

   include 'parameters.h'

   integer, intent(in) :: n, i, k
   real(dp), intent(in) :: df(n)
   real(dp) :: f

   if (k==3) then
      f = -(9.0_dp*df(i-1)+19.0_dp*df(i)-5.0_dp*df(i+1)+df(i+2))/24.0_dp
   elseif (k==4) then
      f = -(25.1_dp*df(i-1)+64.6_dp*df(i)-26.4_dp*df(i+1)+10.6_dp*df(i+2)-1.9_dp*df(i+3))/72.0_dp
   elseif (k==5) then
      f = -(47.5_dp*df(i-1)+142.7_dp*df(i)-79.8_dp*df(i+1)+48.2_dp*df(i+2)-17.3_dp*df(i+3)+2.7*df(i+4))/144.0_dp
   elseif (k==6) then
      f = -(190.87_dp*df(i-1)+651.12_dp*df(i)-464.61_dp*df(i+1)+375.04_dp*df(i+2)-202.11_dp*df(i+3)+ &
        63.12_dp*df(i+4)-8.63_dp*df(i+5))/604.8_dp
   end if

end function AdamsII
!****** End function: Adams* **************************************************
!******************************************************************************

!****** Subroutine: AdamsOutX *************************************************
!******************************************************************************
!
!  Integrate the radial wavefunction outwards. In this subroutine the formula
! from Ref. 1 (2.57 page 39)
!
! INPUT -----------------------------------------------------------------------
!
! real f(n), df(n), d2f(n)        : The function and its derivatives
! real c1(n), c2(2)               : Coefficients of the derivatives
! integer n                       : Size of the function
! intger mn                       : Point up to which integrate
! real h                          : Mesh spacing
! integer k                       : The order of the polynomial
!
! OUTPUT ----------------------------------------------------------------------
!
! real f(n), df(n), d2f(n)        : The function and its derivates
! integer nodes                   : The number of nodes in the function
!
!******************************************************************************
subroutine AdamsOutX(f,df,d2f,c1,c2,n,mn,h,k,nodes)

   implicit none

   include 'parameters.h'

   real(dp), intent(inout) :: f(n), df(n), d2f(n)
   real(dp), intent(in) :: h,c1(n), c2(n)
   integer, intent(in) :: n, mn, k
   integer, intent(out) :: nodes

   real(dp), parameter :: fD(2:8) = (/12.0_dp,24.0_dp,720.0_dp,1440.0_dp, &
                                      60480.0_dp,120960.0_dp,3628800.0_dp/)
   real(dp), parameter :: a(1:9,2:8) = reshape((/&
          -1.0_dp,8.0_dp,5.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
          1.0_dp,-5.0_dp,19.0_dp,9.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
          -19.0_dp,106.0_dp,-264.0_dp,646.0_dp,251.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,&
          27.0_dp,-173.0_dp,482.0_dp,-798.0_dp,1427.0_dp,475.0_dp,0.0_dp,0.0_dp,0.0_dp,&
          -863.0_dp,6312.0_dp,-20211.0_dp,37504.0_dp,-46461.0_dp,65112.0_dp,19087.0_dp,0.0_dp,0.0_dp,&
          1375.0_dp,-11351.0_dp,41499.0_dp,-88547.0_dp,123133.0_dp,-121797.0_dp,139849.0_dp,36799.0_dp,0.0_dp, &
          -33953.0_dp,312874.0_dp,-1291214.0_dp,3146338.0_dp,-5033120.0_dp,5595358.0_dp,-4604594.0_dp,4467094.0_dp, &
          1070017.0_dp/),(/9,7/))

   integer :: m, i, j
   real(dp) :: D, lambda, s1, s2

   if (mn<n) then
      m = mn
   else
      m = n
   end if
   nodes = 0
   lambda = h*a(k+1,k)/fD(k)
   do i=k+1,m-1
      D = 1.0_dp-lambda**2*c1(i+1)*c2(i+1)
      s1 = 0.0_dp
      s2 = 0.0_dp
      do j=1,k
         s1 = s1 + a(j,k)*df(i-k+j)
         s2 = s2 + a(j,k)*d2f(i-k+j)
      end do
      f(i+1) = (f(i)+s1*h/fD(k) + lambda*c1(i+1)*(df(i)+s2*h/fD(k)))/D
      df(i+1) = (lambda*c2(i+1)*(f(i)+s1*h/fD(k)) + df(i)+s2*h/fD(k))/D
      d2f(i+1) = c2(i+1)*f(i+1)
      if (f(i)*f(i+1) <= 0.0_dp) then
         nodes = nodes+1
      end if
   end do

end subroutine AdamsOutX
!****** End subroutine: AdamsOutX *********************************************
!******************************************************************************
