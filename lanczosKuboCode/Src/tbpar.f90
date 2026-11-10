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

   call MIO_Section('Tight-binding parameters','tbpar')
   call MIO_InputParameter('TB.NeighLevels',tbnn,1)
     !! - `TB.NeighLevels`: Number of neighbour shells of the hopping table (at least 1).
   if (tbnn < 1) then
      call MIO_Kill('TB.NeighLevels must be at least 1 (requested: '//trim(num2str(tbnn))//')','tbpar','TBInit')
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
        !! - `TB.Hopping`: Nearest-neighbour hopping g0 [eV], the internal energy unit: every energy is divided by it. Default
        !!   -3.72 aG + 12.14 (2.9887 eV at aG 2.46).
      call MIO_Print('Hopping g0: '//trim(num2str(g0,4))//' eV','tbpar')
      gn(2,1,1) = g0
      gn(1,1,1) = g0
      gn(2,2,1) = g0
      gdef = -3.11_dp*aBN + 10.68_dp
      call MIO_InputParameter('TB.BNHopping',g,gdef)
        !! - `TB.BNHopping`: Nearest-neighbour B-N hopping of hBN [eV] (default: TB.Hopping).
      gn(4,3,1) = g
      call MIO_InputParameter('TB.CBHopping',g,2.68_dp)
        !! - `TB.CBHopping`: Nearest-neighbour C-B hopping [eV].
      gn(3,1,1) = g
      gn(3,2,1) = g
      call MIO_InputParameter('TB.CNHopping',g,2.79_dp)
        !! - `TB.CNHopping`: Nearest-neighbour C-N hopping [eV].
      gn(4,1,1) = g
      gn(4,2,1) = g
   end if
   call MIO_InputParameter('TB.OnSiteC',e0_C,0.0_dp)
     !! - `TB.OnSiteC`: On-site energy of carbon [eV].
   e0_C = e0_C/g0
   call MIO_InputParameter('TB.OnSiteC1',e0_C1,0.0_dp)
     !! - `TB.OnSiteC1`: On-site energy of carbon sublattice 1 (Species 1) [eV].
   e0_C1 = e0_C1/g0
   call MIO_InputParameter('TB.OnSiteC2',e0_C2,0.0_dp)
     !! - `TB.OnSiteC2`: On-site energy of carbon sublattice 2 (Species 2) [eV].
   e0_C2 = e0_C2/g0
   call MIO_InputParameter('TB.OnSiteC1_LB',e0_C1_LB,0.0_dp)
     !! - `TB.OnSiteC1_LB`: On-site energy of carbon sublattice 1 of the bottom layer [eV], for the stack models that use it.
   e0_C1_LB = e0_C1_LB/g0
   call MIO_InputParameter('TB.OnSiteC2_LB',e0_C2_LB,0.0_dp)
     !! - `TB.OnSiteC2_LB`: On-site energy of carbon sublattice 2 of the bottom layer [eV], for the stack models that use it.
   e0_C2_LB = e0_C2_LB/g0
   call MIO_InputParameter('TB.OnSiteC1_LT',e0_C1_LT,0.0_dp)
     !! - `TB.OnSiteC1_LT`: On-site energy of carbon sublattice 1 of the top layer [eV], for the stack models that use it.
   e0_C1_LT = e0_C1_LT/g0
   call MIO_InputParameter('TB.OnSiteC2_LT',e0_C2_LT,0.0_dp)
     !! - `TB.OnSiteC2_LT`: On-site energy of carbon sublattice 2 of the top layer [eV], for the stack models that use it.
   e0_C2_LT = e0_C2_LT/g0
   call MIO_InputParameter('TB.OnSiteB',e0_B,3.09_dp)
     !! - `TB.OnSiteB`: On-site energy of boron [eV].
   e0_B = e0_B/g0
   call MIO_InputParameter('TB.OnSiteN',e0_N,-1.89_dp)
     !! - `TB.OnSiteN`: On-site energy of nitrogen [eV].
   e0_N = e0_N/g0
   call MIO_InputParameter('TB.InterLayerHopping',gIntLay,0.39_dp)
     !! - `TB.InterLayerHopping`: Interlayer hopping t exp(-d/dIntLay) of the simplest interlayer branch of HamHopping, d =
     !!   in-plane distance [eV].
   gIntLay = gIntLay/g0
   if (tbnn>=2) then
      if (MIO_InputSearchLabel('&begin TB.Hoppings.1',line,id)) then
          do i=2,tbnn
             call MIO_InputBlock('TB.Hoppings.'//trim(num2str(i)),gn(:,:,i))
          end do
      else
          do i=2,tbnn
             call MIO_InputParameter('Haldane.T2',t2,0.0_dp)
               !! - `Haldane.T2`: Here: read to leave the second-neighbour table entries at zero; the Haldane term itself is added
               !!   in ham.F90.
             gn(2,1,i) = 0.0_dp !t2  ! We add t2 in the hams.f90 file for the Haldane model
             gn(1,1,i) = 0.0_dp !t2
             gn(2,2,i) = 0.0_dp !t2
          end do
      end if
   end if
   do i=1,tbnn; do i1=1,3; do i2=i1+1,4
      gn(i1,i2,i) = gn(i2,i1,i)
   end do; end do; end do
   gn = gn/g0
   call MIO_InputParameter('Run.TypeOfSystem',line,'Graphene')
     !! - `Run.TypeOfSystem`: Here: Ribbons, Hybrid and ReadXYZ structures get the hopping table of mixed C/B/N systems.
   if (MIO_StringComp(line,'Ribbons')  .or. MIO_StringComp(line,'Hybrid') .or. MIO_StringComp(line,'ReadXYZ')) then
      call MIO_Print('Hopping C-C: '//trim(num2str(g0,4))//' eV','tbpar')
      call MIO_Print('Hopping B-N: '//trim(num2str(gn(4,3,1)*g0,4))//' eV','tbpar')
   end if
   call MIO_InputParameter('Hubbard.U_C',U(1),3.0_dp)
     !! - `Hubbard.U_C`: Hubbard U of carbon [eV] (two-spin runs with Hubbard.EnableSCF).
   U(2) = U(1)
   call MIO_InputParameter('Hubbard.U_B',U(3),3.0_dp)
     !! - `Hubbard.U_B`: Hubbard U of boron [eV] (two-spin runs with Hubbard.EnableSCF).
   call MIO_InputParameter('Hubbard.U_N',U(4),3.0_dp)
     !! - `Hubbard.U_N`: Hubbard U of nitrogen [eV] (two-spin runs with Hubbard.EnableSCF).
   U = U/g0

end subroutine TBInit

end module tbpar
