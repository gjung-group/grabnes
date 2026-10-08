module tbpar

   use mio

   implicit none

   PRIVATE

   real(dp), public, save :: g0 ! Hopping parameter
   real(dp), public, save :: e0_C, e0_B, e0_N ! On-site energy for carbon, boron and nitrogen
   real(dp), public, save :: e0_C1, e0_C2, e0_C1_LB, e0_C2_LB, e0_C1_LT, e0_C2_LT
   real(dp), public, save :: gIntLay, U(4)
   integer, public, save :: tbnn
   real(dp), public, pointer :: gn(:,:,:)

   public :: TBInit

contains

subroutine TBInit()

   use cell,                 only : aG, aBN

   integer :: i, i1, i2
   character(len=80) :: line
   integer :: id
   real(dp) :: g, gdef
   real(dp) :: t2

   call MIO_InputParameter('TB.NeighLevels',tbnn,1)
   if (tbnn < 1) then
      call MIO_Kill('TB.NeighLevels must be at least 1 (requested: '//trim(num2str(tbnn))//').','tbpar','TBInit')
   end if
   call MIO_Allocate(gn,[4,4,tbnn],'gn','tbpar')
   if (MIO_InputSearchLabel('&begin TB.Hoppings.1',line,id)) then
      call MIO_InputBlock('TB.Hoppings.1',gn(:,:,1))
      g0 = gn(2,1,1)
   else
      ! The default must not be the variable that receives the value: with
      ! both arguments aliased the result is undefined (zero in optimised builds).
      gdef = -3.72_dp*aG + 12.14_dp
      call MIO_InputParameter('TB.Hopping',g0,gdef)
      call MIO_Print('Hopping g0: '//trim(num2str(g0)),'tbpar')
      gn(2,1,1) = g0
      gn(1,1,1) = g0
      gn(2,2,1) = g0
      gdef = -3.11_dp*aBN + 10.68_dp
      call MIO_InputParameter('TB.BNHopping',g,gdef)
      gn(4,3,1) = g
      call MIO_InputParameter('TB.CBHopping',g,2.68_dp)
      gn(3,1,1) = g
      gn(3,2,1) = g
      call MIO_InputParameter('TB.CNHopping',g,2.79_dp)
      gn(4,1,1) = g
      gn(4,2,1) = g
   end if
   call MIO_InputParameter('TB.OnSiteC',e0_C,0.0_dp)
   e0_C = e0_C/g0
   !call MIO_InputParameter('TB.OnSiteCA',e0_CA,0.0_dp)
   !e0_CA = e0_CA/g0
   !call MIO_InputParameter('TB.OnSiteCB',e0_CB,0.0_dp)
   !e0_CB = e0_CB/g0
   call MIO_InputParameter('TB.OnSiteC1',e0_C1,0.0_dp)
   e0_C1 = e0_C1/g0
   call MIO_InputParameter('TB.OnSiteC2',e0_C2,0.0_dp)
   e0_C2 = e0_C2/g0
   call MIO_InputParameter('TB.OnSiteC1_LB',e0_C1_LB,0.0_dp)
   e0_C1_LB = e0_C1_LB/g0
   call MIO_InputParameter('TB.OnSiteC2_LB',e0_C2_LB,0.0_dp)
   e0_C2_LB = e0_C2_LB/g0
   call MIO_InputParameter('TB.OnSiteC1_LT',e0_C1_LT,0.0_dp)
   e0_C1_LT = e0_C1_LT/g0
   call MIO_InputParameter('TB.OnSiteC2_LT',e0_C2_LT,0.0_dp)
   e0_C2_LT = e0_C2_LT/g0
   call MIO_InputParameter('TB.OnSiteB',e0_B,3.09_dp)
   e0_B = e0_B/g0
   call MIO_InputParameter('TB.OnSiteN',e0_N,-1.89_dp)
   e0_N = e0_N/g0
   call MIO_InputParameter('TB.InterLayerHopping',gIntLay,0.39_dp)
   gIntLay = gIntLay/g0
   if (tbnn>=2) then
      if (MIO_InputSearchLabel('&begin TB.Hoppings.1',line,id)) then
          do i=2,tbnn
             call MIO_InputBlock('TB.Hoppings.'//trim(num2str(i)),gn(:,:,i))
          end do
      else
          do i=2,tbnn
             call MIO_InputParameter('HaldaneT2',t2,0.0_dp)
             gn(2,1,i) = 0.0_dp !t2  ! We add t2 in the hams.f90 file for the Haldane model
             gn(1,1,i) = 0.0_dp !t2
             gn(2,2,i) = 0.0_dp !t2
          end do
          print*, "HaldaneT2 = ", t2
      end if
   end if
   do i=1,tbnn; do i1=1,3; do i2=i1+1,4
      gn(i1,i2,i) = gn(i2,i1,i)
   end do; end do; end do
   gn = gn/g0
   call MIO_InputParameter('TypeOfSystem',line,'Graphene')
   if (MIO_StringComp(line,'Ribbons')  .or. MIO_StringComp(line,'Hybrid') .or. MIO_StringComp(line,'ReadXYZ')) then
      call MIO_Print('Hopping C-C: '//trim(num2str(g0,3)),'tbpar')
      call MIO_Print('Hopping B-N: '//trim(num2str(gn(4,3,1)*g0,3)),'tbpar')
      call MIO_Print('')
   end if
   call MIO_InputParameter('HubbardU_C',U(1),3.0_dp)
   U(2) = U(1)
   call MIO_InputParameter('HubbardU_B',U(3),3.0_dp)
   call MIO_InputParameter('HubbardU_N',U(4),3.0_dp)
   U = U/g0

end subroutine TBInit

end module tbpar
