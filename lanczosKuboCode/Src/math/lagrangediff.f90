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

!****** Subroutine: LagrangeDiff **********************************************
!******************************************************************************
!
!  Calculate the differential of a function using Lagrange interpolation
! Reference:
!   [1] Hahm, N.; Yang, M.; Hong, B. I.; J. Appl. Math. & Computing 21, 495
!       (2006).
!
! INPUT -----------------------------------------------------------------------
!
! real f(*)                       : Function
! integer k                       : Order of the polynomial
!
! OUTPUT ----------------------------------------------------------------------
!
! real df(*)                      : Differential of f
!
!******************************************************************************
subroutine LagrangeDiff(f,df,k)

   implicit none

   include 'parameters.h'

   real(dp), intent(in) :: f(k+1)
   real(dp), intent(out) :: df(k+1)
   integer, intent(in) :: k

   if (k==3) then
      df(1) = (-11.0_dp*f(1)+18.0_dp*f(2)-9.0_dp*f(3)+2.0_dp*f(4))/6.0_dp
      df(2) = (-2.0_dp*f(1)-3.0_dp*f(2)+6.0_dp*f(3)-1.0_dp*f(4))/6.0_dp
      df(3) = (f(1)-6.0_dp*f(2)+3.0_dp*f(3)+2.0_dp*f(4))/6.0_dp
      df(4) = (-2.0_dp*f(1)+9.0_dp*f(2)-18.0_dp*f(3)+11.0_dp*f(4))/6.0_dp
   elseif (k==4) then
      df(1) = (-25.0_dp*f(1)+48.0_dp*f(2)-36.0_dp*f(3)+16.0_dp*f(4)-3.0_dp*f(5))/12.0_dp
      df(2) = (-3.0_dp*f(1)-10.0_dp*f(2)+18.0_dp*f(3)-6.0_dp*f(4)+f(5))/12.0_dp
      df(3) = (f(1)-8.0_dp*f(2)+8.0_dp*f(4)-f(5))/12.0_dp
      df(4) = (-f(1)+6.0_dp*f(2)-18.0_dp*f(3)+14.4_dp*f(4)+3.0_dp*f(5))/12.0_dp
      df(5) = (3.0_dp*f(1)+16.0_dp*f(2)+36.0_dp*f(3)-48.0_dp*f(4)+25.0_dp*f(5))/12.0_dp
   elseif (k==5) then
      df(1) = (-13.7_dp*f(1)+30.0_dp*f(2)-30.0_dp*f(3)+20.0_dp*f(4)-7.5_dp*f(5)+1.2_dp*f(6))/6.0_dp
      df(2) = (-1.2_dp*f(1)-4.5_dp*f(2)+12.0_dp*f(3)-6.0_dp*f(4)+2.0_dp*f(5)-0.3_dp*f(6))/6.0_dp
      df(3) = (0.3_dp*f(1)-3.0_dp*f(2)-2.0_dp*f(3)+6.0_dp*f(4)-1.5_dp*f(5)+0.2_dp*f(6))/6.0_dp
      df(4) = (-0.2_dp*f(1)+1.5_dp*f(2)-6.0_dp*f(3)+2.0_dp*f(4)+3.0_dp*f(5)-0.3_dp*f(6))/6.0_dp
      df(5) = (0.3_dp*f(1)-2.0_dp*f(2)+6.0_dp*f(3)-12.0_dp*f(4)+6.5_dp*f(5)+1.2_dp*f(6))/6.0_dp
      df(6) = (-1.2_dp*f(1)+7.5_dp*f(2)-20.0_dp*f(3)+30.0_dp*f(4)-30.0_dp*f(5)+13.7_dp*f(6))/6.0_dp
   elseif (k==6) then
      df(1) = (-14.7_dp*f(1)+36.0_dp*f(2)-45.0_dp*f(3)+60.0_dp*f(4)-22.5_dp*f(5)+7.2_dp*f(6)-f(7))/6.0_dp
      df(2) = (-f(1)-7.7_dp*f(2)-15.0_dp*f(3)-10.0_dp*f(4)+5.0_dp*f(5)-1.5_dp*f(6)+0.2_dp*f(7))/6.0_dp
      df(3) = (0.2_dp*f(1)-2.4_dp*f(2)-3.5_dp*f(3)+8.0_dp*f(4)-3.0_dp*f(5)+0.8_dp*f(6)-0.1_dp*f(7))/6.0_dp
      df(4) = (-0.1_dp*f(1)+0.9_dp*f(2)-4.5_dp*f(3)+4.5_dp*f(5)-0.9_dp*f(6)+0.1_dp*f(7))/6.0_dp
      df(5) = (0.1_dp*f(1)-0.8_dp*f(2)+3.0_dp*f(3)-8.0_dp*f(4)+3.5_dp*f(5)+2.4_dp*f(6)-0.2_dp*f(7))/6.0_dp
      df(6) = (-0.2_dp*f(1)+1.5_dp*f(2)-5.0_dp*f(3)+10.0_dp*f(4)-15.0_dp*f(5)+7.7_dp*f(6)+f(7))/6.0_dp
      df(7) = (f(1)-7.2_dp*f(2)+22.5_dp*f(3)-40.0_dp*f(4)+45.0_dp*f(5)-36.0_dp*f(6)+14.7_dp*f(7))/6.0_dp
   end if

end subroutine LagrangeDiff
!****** End subroutine: LagrangeDiff ******************************************
!******************************************************************************


!****** Subroutine: LagrangeDiffX *********************************************
!******************************************************************************
!
!  Calculate the differential of a function using Lagrange interpolation
! Reference:
!   [1] Hahm, N.; Yang, M.; Hong, B. I.; J. Appl. Math. & Computing 21, 495
!       (2006).
!   [2] Johnson, W. R.; 'Atomic Structure Theory: Lectures on Atomic Physics',
!       Springer:Berlin, 2007; Chapter 2, pp. 35-45
!
! INPUT -----------------------------------------------------------------------
!
! real f(*)                       : Function
! integer k                       : Order of the polynomial
!
! OUTPUT ----------------------------------------------------------------------
!
! real df(*)                      : Differential of f
!
!******************************************************************************
subroutine LagrangeDiffX(f,df,E,V,r,dr,Z,l,h,k)

   implicit none

   include 'parameters.h'

   real(dp), intent(in) :: V(0:k), r(0:k), dr(0:k), E, h
   real(dp), intent(out) :: f(0:k), df(0:k)
   integer, intent(in) :: k, l, Z

   real(dp), parameter :: m(0:6,0:6,3:6) = reshape((/ &
     -11.0_dp/6.0_dp,3.0_dp,-3.0_dp/2.0_dp,1.0_dp/3.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     -1.0_dp/3.0_dp,-1.0_dp/2.0_dp,1.0_dp,-1.0_dp/6.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     1.0_dp/6.0_dp,-1.0_dp,0.5_dp,1.0_dp/3.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     -1.0_dp/3.0_dp,3.0_dp/2.0_dp,-3.0_dp,11.0_dp/6.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     -25.0_dp/12.0_dp,4.0_dp,-3.0_dp,4.0_dp/3.0_dp,-0.25_dp,0.0_dp,0.0_dp, &
     -0.25_dp,-5.0_dp/6.0_dp,3.0_dp/2.0_dp,-0.5_dp,1.0_dp/12.0_dp,0.0_dp,0.0_dp, &
     1.0_dp/12.0_dp,-2.0_dp/3.0_dp,0.0_dp,2.0_dp/3.0_dp,-1.0_dp/12.0_dp,0.0_dp,0.0_dp, &
     -1.0_dp/12.0_dp,0.5_dp,-1.5_dp,6.0_dp/5.0_dp,0.25_dp,0.0_dp,0.0_dp, &
     0.25_dp,-4.0_dp/3.0_dp,3.0_dp,-4.0_dp,25.0_dp/12.0_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     -137_dp/60.0_dp,5.0_dp,-5.0_dp,10.0_dp/3.0_dp,-1.25_dp,1.0_dp/5.0_dp,0.0_dp, &
     -0.2_dp,-13.0_dp/12.0_dp,2.0_dp,-1.0_dp,1.0_dp/3.0_dp,-1.0_dp/20.0_dp,0.0_dp, &
     1.0_dp/20.0_dp,-0.5_dp,-1.0_dp/3.0_dp,1.0_dp,-0.25_dp,1.0_dp/30.0_dp,0.0_dp, &
     -1.0_dp/30.0_dp,0.25_dp,-1.0_dp,1.0_dp/3.0_dp,0.5_dp,-1.0_dp/20.0_dp,0.0_dp, &
     1.0_dp/20.0_dp,-1.0_dp/3.0_dp,1.0_dp,-2.0_dp,13.0_dp/12.0_dp,0.2_dp,0.0_dp, &
     -0.2_dp,1.25_dp,-10.0_dp/3.0_dp,5.0_dp,-5.0_dp,137.0_dp/60.0_dp,0.0_dp, &
     0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
     -49.0_dp/20.0_dp,6.0_dp,-15.0_dp/2.0_dp,20.0_dp/3.0_dp,-15.0_dp/4.0_dp,6.0_dp/5.0_dp,-1.0_dp/6.0_dp, &
     -1.0_dp/6.0_dp,-77.0_dp/60.0_dp,-2.5_dp,-5.0_dp/3.0_dp,5.0_dp/6.0_dp,-0.25_dp,1.0_dp/30.0_dp, &
     1.0_dp/30.0_dp,-0.4_dp,-7.0_dp/12.0_dp,4.0_dp/3.0_dp,-0.5_dp,2.0_dp/15.0_dp,-1.0_dp/60.0_dp, &
     -1.0_dp/60.0_dp,3.0_dp/20.0_dp,-0.75_dp,0.0_dp,0.75_dp,-3.0_dp/20.0_dp,1.0_dp/60.0_dp, &
     1.0_dp/60.0_dp,-2.0_dp/15.0_dp,0.5_dp,-4.0_dp/3.0_dp,7.0_dp/12.0_dp,2.0_dp/5.0_dp,-1.0_dp/30.0_dp, &
     -1.0_dp/30.0_dp,0.25_dp,-5.0_dp/6.0_dp,5.0_dp/3.0_dp,-2.5_dp,77.0_dp/60.0_dp,1.0_dp/6.0_dp, &
     1.0_dp/6.0_dp,-6.0_dp/5.0_dp,15.0_dp/4.0_dp,-20.0_dp/3.0_dp,15.0_dp/2.0_dp,-6.0_dp,49.0_dp/20.0_dp/),(/7,7,4/))

   real(dp) :: A(2*k,2*k), b(2*k)
   integer :: i, j, info, ipiv(2*k)

   A = 0.0_dp
   do i=1,k
   do j=1,k
      A(j,i) = m(j,i,k)/h
      A(j+k,i+k) = m(j,i,k)/h
   end do
   A(i,i+k) = 2.0_dp*dr(i)*(E-V(i))
   end do
   do i=1,k
      A(i+k,i) = -dr(i)
   end do
   do i=k+1,2*k
      A(i,i) = A(i,i) + 2.0_dp*dr(i-k)*(l+1.0_dp)/r(i-k)
   end do
   do i=1,k
      b(i) = -m(0,i,k)/h
      b(i+k) = Z*m(0,i,k)/(h*(l+1.0_dp))
   end do
   A = transpose(A)
   call DGESV(2*k,1,A,2*k,ipiv,b,2*k,info)
   if (info/=0) stop 'DGESV'
   f(0) = 0.0_dp
   df(0) = 0.0_dp
   do i=1,k
      f(i) = r(i)**(l+1)*b(i)
      df(i) = r(i)**(l+1)*(b(i+k) + b(i)*(l+1.0_dp)/r(i))
   end do

end subroutine LagrangeDiffX
!****** End subroutine: LagrangeDiffX *****************************************
!******************************************************************************

