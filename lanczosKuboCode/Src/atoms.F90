module atoms

   use mio
   use math

   implicit none

   PRIVATE

   real(dp), public, pointer :: Rat(:,:), RatC(:,:), RatBN(:,:), RatC1(:,:), RatC2(:,:), RatC3(:,:)
   real(dp), public, pointer :: RatTemp(:,:), RatInit(:,:)
   integer, public, save :: nAt, nAtC, nAtBN, nAtC1, nAtC2, nAtC3
   integer, public, pointer :: indxDiv(:), procAt(:), indxNode(:), secAt(:)
   integer, public, pointer :: Species(:), SpeciesTemp(:)
   integer, public, pointer :: layerIndex(:)
   real(dp), public, pointer :: interlayerDistances(:)=>NULL(), displacements(:,:)=>NULL()
   real(dp), public, pointer :: displacements_b(:,:), displacements_t(:,:)
   real(dp), public :: dIntLay, nEl
   real(dp), public :: phiForEffectiveModel
   ! For species:
   ! 1 - Carbon site A
   ! 2 - Carbon site B
   ! 3 - Boron
   ! 4 - Nitrogen
   integer, public, save :: inode1, inode2, in1, in2
   logical, public, save :: frac=.true.
   !$OMP THREADPRIVATE(in1,in2)

   public :: AtomsPos
   public :: AtomsSetFrac
   public :: AtomsSetCart
   public :: AtomsRotate

contains

! counterclockwise if angle is positive
subroutine AtomsRotate(angle)

   use constants,            only : pi

   integer i
   real(dp), intent(in) :: angle
   real(dp) :: a, tempRat1

   a = angle*pi/180.0_dp
   !$OMP PARALLEL DO
   do i=1,nAt
        tempRat1 = Rat(1,i)
        Rat(1,i) = Rat(1,i) * cos(a) - Rat(2,i) * sin(a)
        Rat(2,i) = tempRat1 * sin(a) + Rat(2,i) * cos(a)
   end do
   !$OMP END PARALLEL DO

end subroutine AtomsRotate

subroutine AtomsPos()

   use math
   use parallel,             only : nDiv, nTh, procID
   use cell,                 only : sCell, sCell2, aG, ucell, rcell, volume, area, aBN
   use constants,            only : pi, twopi
   use hybrid,               only : HybridGen
   !use gauss,                only : GaussHeight
   use name,                 only : prefix

   character(len=80) :: str, str2, str3, str4, str5
   integer :: m(2), m2(2), i, j, n(4), sh, nC, nBN
   real(dp) :: g, delta, phi, h, yShift, BLAngle
   real(dp) :: vn(3)
   integer :: shiftFactor, ncell
   character(len=3) :: spc

   logical :: l, ll
   logical :: randomStrain
   logical :: addShift, Bernal, bridge

   real(dp) :: limit0, limit1, limit2, limit3, DBShiftRatio, boundaryWidth
   logical :: createBLDomainBoundary, AtomsOrderDeactivated
   logical :: boundaryTypeArmAA, boundaryTypeArmSP
   integer :: numberOfShifts

   real(dp) :: ggg, ORZ

   logical :: twoLayers, oneLayer, GBNtwoLayers, GBNuseDisplacementFile, encapsulatedFourLayers, tBGuseDisplacementFile
   logical :: t2GBN
   logical :: t2BG
   logical :: t3BG
  logical :: fourLayers
  logical :: fourLayersSandwiched
  logical :: helicalTwistedMBM
  logical :: fiveLayersSandwiched
   logical :: sixLayersSandwiched
   logical :: sevenLayersSandwiched
   logical :: eightLayersSandwiched
   logical :: tenLayersSandwiched
   logical :: twentyLayersSandwiched
   logical :: threeLayers

   logical :: useSublatticeFile, readInterlayerDistances, distanceDependentEffectiveModel
   logical :: readRigidXYZ, invertDisplacements, readLayerIndex

   real(dp) :: twoLayersZ1
   real(dp) :: twoLayersZ2
   real(dp) :: fourLayersZ1
   real(dp) :: fourLayersZ2
   real(dp) :: fourLayersZ3
   real(dp) :: fourLayersZ4
   real(dp) :: threeLayersZ1
   real(dp) :: threeLayersZ2
   real(dp) :: threeLayersZ3
   real(dp) :: fiveLayersZ1
   real(dp) :: fiveLayersZ2
   real(dp) :: fiveLayersZ3
   real(dp) :: fiveLayersZ4
   real(dp) :: fiveLayersZ5
   real(dp) :: sixLayersZ1
   real(dp) :: sixLayersZ2
   real(dp) :: sixLayersZ3
   real(dp) :: sixLayersZ4
   real(dp) :: sixLayersZ5
   real(dp) :: sixLayersZ6
   real(dp) :: sevenLayersZ1
   real(dp) :: sevenLayersZ2
   real(dp) :: sevenLayersZ3
   real(dp) :: sevenLayersZ4
   real(dp) :: sevenLayersZ5
   real(dp) :: sevenLayersZ6
   real(dp) :: sevenLayersZ7
   real(dp) :: eightLayersZ1
   real(dp) :: eightLayersZ2
   real(dp) :: eightLayersZ3
   real(dp) :: eightLayersZ4
   real(dp) :: eightLayersZ5
   real(dp) :: eightLayersZ6
   real(dp) :: eightLayersZ7
   real(dp) :: eightLayersZ8
   real(dp) :: tenLayersZ1
   real(dp) :: tenLayersZ2
   real(dp) :: tenLayersZ3
   real(dp) :: tenLayersZ4
   real(dp) :: tenLayersZ5
   real(dp) :: tenLayersZ6
   real(dp) :: tenLayersZ7
   real(dp) :: tenLayersZ8
   real(dp) :: tenLayersZ9
   real(dp) :: tenLayersZ10
   real(dp) :: twentyLayersZ1
   real(dp) :: twentyLayersZ2
   real(dp) :: twentyLayersZ3
   real(dp) :: twentyLayersZ4
   real(dp) :: twentyLayersZ5
   real(dp) :: twentyLayersZ6
   real(dp) :: twentyLayersZ7
   real(dp) :: twentyLayersZ8
   real(dp) :: twentyLayersZ9
   real(dp) :: twentyLayersZ10
   real(dp) :: twentyLayersZ11
   real(dp) :: twentyLayersZ12
   real(dp) :: twentyLayersZ13
   real(dp) :: twentyLayersZ14
   real(dp) :: twentyLayersZ15
   real(dp) :: twentyLayersZ16
   real(dp) :: twentyLayersZ17
   real(dp) :: twentyLayersZ18
   real(dp) :: twentyLayersZ19
   real(dp) :: twentyLayersZ20

   integer :: dontcare
   integer :: ii, jj, kk
   integer :: SuperCellX, SuperCellY
   logical :: encapsulatedThreeLayers

   !real(dp) :: phiForEffectiveModel

   !real(dp), pointer :: RatTemp(:,:)
   !integer, pointer :: SpeciesTemp(:)

#ifdef DEBUG
   call MIO_Debug('AtomsPos',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('atoms')
#endif /* TIMER */

   call MIO_InputParameter('encapsulatedThreeLayers',encapsulatedThreeLayers,.false.)
   call MIO_InputParameter('TypeOfSystem',str,'Graphene')
   frac = .true.
   if (MIO_StringComp(str,'Graphene') .or. MIO_StringComp(str,'BoronNitride') .or. &
     MIO_StringComp(str,'Hybrid') .or. MIO_StringComp(str,'BLtoSLYoungju')) then
      call MIO_InputParameter('basedOnMoireCellParameters',ll,.false.)
      if (ll) then
          call MIO_InputParameter('MoireCellParameters',n,[0,0,0,0])
          call MIO_InputParameter('CellHeight',h,40.0_dp)
          call MIO_InputParameter('InterlayerDistance',dIntLay,3.22_dp) ! Ref. PRB 76, 73103
          g = n(1)**2 + n(2)**2 + n(1)*n(2)
          delta = sqrt(real(n(3)**2 + n(4)**2 + n(3)*n(4))/g)
          phi = acos((2.0_dp*n(1)*n(3)+2.0_dp*n(2)*n(4) + n(1)*n(4) + n(2)*n(3))/(2.0_dp*delta*g))
          if (delta<1.0_dp) then
             delta = 1.0_dp/delta
          else
             i = n(1)
             n(1) = n(3)
             n(3) = i
             i = n(2)
             n(2) = n(4)
             n(4) = i
          end if
          n = n*sCell
          g = dIntLay/(2.0_dp*h)
          call AtomsConstruct(RatC,n(1:2),0.0_dp,'graphene',g, 0.0_dp)
          nAtC = (n(1)**2 + n(1)*n(2) + n(2)**2)*2
          call MIO_Print('Atoms of C:  '//trim(num2str(nAtC)),'atoms')
          phiForEffectiveModel = -phi*180.0_dp/pi
          !call AtomsConstruct(RatBN,n(3:4),-phi*180.0_dp/pi,'BN',-g, 0.0_dp)
          !nAtBN = (n(3)**2 + n(3)*n(4) + n(4)**2)*2
          !call MIO_Print('Atoms of BN: '//trim(num2str(nAtBN)),'atoms')
          nAt = nAtC
          call MIO_Allocate(Rat,[3,nAt],'Rat','atoms')
          Rat(:,1:nAtC) = RatC
          !Rat(:,nAtC+1:nAt) = RatBN
          call MIO_Print('Total number of atoms: '//trim(num2str(nAt)),'atoms')
          call MIO_Allocate(Species,nAt,'Species','atoms')
          call MIO_Allocate(layerIndex,nAt,'layerIndex','atoms')
          !$OMP PARALLEL DO
          do i=1,nAtC
             Species(i) = mod(i+1,2) + 1
             layerIndex(i) = 1
          end do
          !$OMP END PARALLEL DO

          !!$OMP PARALLEL DO
          !do i=nAtC+1,nAt
          !   Species(i) = mod(i+1,2) + 3
          !end do
          !!$OMP END PARALLEL DO
      else
          call MIO_InputParameter('CellSize',m(1),50)
          call MIO_InputParameter('SuperCellAsymmetric',l,.false.)
          if (l) then
             m2(1) = m(1)*sCell2
             m2(2) = 0
             m(1) = m(1)*sCell
             m(2) = 0
          else
             m(1) = m(1)*sCell
             m(2) = 0
          end if
          if (MIO_StringComp(str,'Graphene') .or. MIO_StringComp(str,'Hybrid') &
               .or. MIO_StringComp(str,'BLtoSLYoungju')) then
             if (l) then
                call AtomsConstructAsymm(Rat,m,m2,0.0_dp,'graphene',0.0_dp, 0.0_dp)
             else
                call AtomsConstruct(Rat,m,0.0_dp,'graphene',0.0_dp, 0.0_dp)
             end if
             sh = 1
          else if (MIO_StringComp(str,'BoronNitride')) then
             call AtomsConstruct(Rat,m,0.0_dp,'BN',0.0_dp, 0.0_dp)
             sh = 3
          end if
          nAt = size(Rat,2)
          call MIO_Print('Number of atoms: '//trim(num2str(nAt)),'atoms')
          call MIO_Allocate(Species,nAt,'Species','atoms')
          call MIO_Allocate(layerIndex,nAt,'layerIndex','atoms')
          !$OMP PARALLEL DO &
          !$OMP& SHARED (Species, sh, nAt, layerIndex), &
          !$OMP& DEFAULT(NONE)
          do i=1,nAt
             Species(i) = mod(i+1,2) + sh
             layerIndex(i) = 1
             !print*, i, Species(i)
          end do
          !$OMP END PARALLEL DO
          call MIO_InputParameter('distanceDependentEffectiveModel',distanceDependentEffectiveModel,.false.)
          if (distanceDependentEffectiveModel) then
             call MIO_Allocate(interlayerDistances,nAt*sCell*sCell,'interlayerDistances','atoms')
             do i=1,nAt
                interlayerDistances(i) = 3.35_dp
             end do
          end if
          !print*, interlayerDistances
      end if
   else if (MIO_StringComp(str,'Graphene_Over_BN')) then
      call MIO_InputParameter('MoireCellParameters',n,[0,0,0,0])
      call MIO_InputParameter('CellHeight',h,40.0_dp)
      call MIO_InputParameter('InterlayerDistance',dIntLay,3.22_dp) ! Ref. PRB 76, 73103
      g = n(1)**2 + n(2)**2 + n(1)*n(2)
      delta = sqrt(real(n(3)**2 + n(4)**2 + n(3)*n(4))/g)
      phi = acos((2.0_dp*n(1)*n(3)+2.0_dp*n(2)*n(4) + n(1)*n(4) + n(2)*n(3))/(2.0_dp*delta*g))
      if (delta<1.0_dp) then
         delta = 1.0_dp/delta
      else
         i = n(1)
         n(1) = n(3)
         n(3) = i
         i = n(2)
         n(2) = n(4)
         n(4) = i
      end if
      n = n*sCell
      g = dIntLay/(2.0_dp*h)
      call AtomsConstruct(RatC,n(1:2),0.0_dp,'graphene',g, 0.0_dp)
      nAtC = (n(1)**2 + n(1)*n(2) + n(2)**2)*2
      call MIO_Print('Atoms of C:  '//trim(num2str(nAtC)),'atoms')
      call AtomsConstruct(RatBN,n(3:4),-phi*180.0_dp/pi,'BN',-g, 0.0_dp)
      nAtBN = (n(3)**2 + n(3)*n(4) + n(4)**2)*2
      call MIO_Print('Atoms of BN: '//trim(num2str(nAtBN)),'atoms')
      nAt = nAtC + nAtBN
      call MIO_Allocate(Rat,[3,nAt],'Rat','atoms')
      Rat(:,1:nAtC) = RatC
      Rat(:,nAtC+1:nAt) = RatBN
      call MIO_Print('Total number of atoms: '//trim(num2str(nAt)),'atoms')
      call MIO_Allocate(Species,nAt,'Species','atoms')
      !$OMP PARALLEL DO
      do i=1,nAtC
         Species(i) = mod(i+1,2) + 1
      end do
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO
      do i=nAtC+1,nAt
         Species(i) = mod(i+1,2) + 3
      end do
      !$OMP END PARALLEL DO
   else if (MIO_StringComp(str,'BilayerGraphene_Over_BN')) then
      call MIO_InputParameter('MoireCellParameters',n,[0,0,0,0])
      call MIO_InputParameter('CellHeight',h,40.0_dp)
      call MIO_InputParameter('InterlayerDistance',dIntLay,3.22_dp) ! Ref. PRB 76, 73103
      g = n(1)**2 + n(2)**2 + n(1)*n(2)
      delta = sqrt(real(n(3)**2 + n(4)**2 + n(3)*n(4))/g)
      phi = acos((2.0_dp*n(1)*n(3)+2.0_dp*n(2)*n(4) + n(1)*n(4) + n(2)*n(3))/(2.0_dp*delta*g))
      if (delta<1.0_dp) then
         delta = 1.0_dp/delta
      else
         i = n(1)
         n(1) = n(3)
         n(3) = i
         i = n(2)
         n(2) = n(4)
         n(4) = i
      end if
      n = n*sCell
      g = dIntLay/(2.0_dp*h)
      call AtomsConstruct(RatC,n(1:2),0.0_dp,'graphene',g, 0.0_dp)
      nAtC = (n(1)**2 + n(1)*n(2) + n(2)**2)*2
      nAtC2 = nAtC
      call MIO_Print('Atoms of C:  '//trim(num2str(nAtC)),'atoms')
      call AtomsConstruct(RatBN,n(3:4),-phi*180.0_dp/pi,'BN',-g, 0.0_dp)
      nAtBN = (n(3)**2 + n(3)*n(4) + n(4)**2)*2
      call MIO_Print('Atoms of BN: '//trim(num2str(nAtBN)),'atoms')
      call MIO_InputParameter('twistedBLAddShift',addShift,.false.)
      call MIO_InputParameter('BernalShift',Bernal,.false.)
      call MIO_InputParameter('bridgeShift',bridge,.false.)
      call MIO_InputParameter('BilayerShiftFactor',shiftFactor,0)
      !print*, "shiftinfo", n(1), shiftFactor
      if (addShift) then
         if (Bernal) then
            yShift = (aG/3.0_dp)/(n(1)*aG)+shiftFactor!*aG !
            !print*, "yShift= ", yShift, aG/3.0_dp
         else if (bridge) then
            yShift = (aG/3.0_dp/2.0_dp)/(n(1)*aG)+shiftFactor!*aG !
            !print*, "yShift= ", yShift, aG/3.0_dp
         end if
      else
         yShift = 0.0_dp
         !print*, "yShift= ", yShift, aG/3.0_dp
      end if
      call AtomsConstruct(RatC2,n(1:2),0.0_dp,'graphene',3.0_dp*g, yShift)
      nAt = nAtC + nAtBN + nAtC2
      call MIO_Allocate(Rat,[3,nAt],'Rat','atoms')
      Rat(:,1:nAtC) = RatC
      Rat(:,nAtC+1:nAtC+nAtBN) = RatBN
      Rat(:,nAtC+nAtBN+1:nAt) = RatC2
      call MIO_Print('Total number of atoms: '//trim(num2str(nAt)),'atoms')
      call MIO_Allocate(Species,nAt,'Species','atoms')
      !$OMP PARALLEL DO
      do i=1,nAtC
         Species(i) = mod(i+1,2) + 1
      end do
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO
      do i=nAtC+1,nAt
         Species(i) = mod(i+1,2) + 3
      end do
      !$OMP END PARALLEL DO
   else if (MIO_StringComp(str,'TwistedBilayerBasedOnMoireCell')) then ! Other twisted bilayers seem to have a wrong shift
      call MIO_InputParameter('CellSize',m(1),50)
      call MIO_InputParameter('MoireCellParameters',n,[0,0,0,0])
      call MIO_InputParameter('CellHeight',h,40.0_dp)
      call MIO_InputParameter('InterlayerDistance',dIntLay,3.22_dp) ! Ref. PRB 76, 73103
      g = n(1)**2 + n(2)**2 + n(1)*n(2)
      delta = sqrt(real(n(3)**2 + n(4)**2 + n(3)*n(4))/g)
      phi = acos((2.0_dp*n(1)*n(3)+2.0_dp*n(2)*n(4) + n(1)*n(4) + n(2)*n(3))/(2.0_dp*delta*g))
      !print*, "phi = ", phi*180_dp/pi
      if (delta<1.0_dp) then
         delta = 1.0_dp/delta
      else
         i = n(1)
         n(1) = n(3)
         n(3) = i
         i = n(2)
         n(2) = n(4)
         n(4) = i
      end if
      n = n*sCell
      g = dIntLay/(2.0_dp*h)
      call AtomsConstruct(RatC,n(1:2),0.0_dp,'graphene',g, 0.0_dp)
      nAtC = (n(1)**2 + n(1)*n(2) + n(2)**2)*2
      call MIO_Print('Atoms of C:  '//trim(num2str(nAtC)),'atoms')
      call MIO_InputParameter('twistedBLAddShift',addShift,.false.)
      call MIO_InputParameter('BernalShift',Bernal,.false.)
      call MIO_InputParameter('bridgeShift',bridge,.false.)
      call MIO_InputParameter('BilayerShiftFactor',shiftFactor,0)
      !print*, "shiftinfo", n(1), shiftFactor
      if (addShift) then
         if (Bernal) then
            yShift = (aG/3.0_dp)/(n(1)*aG)+shiftFactor!*aG !
            !print*, "yShift= ", yShift, aG/3.0_dp
         else if (bridge) then
            yShift = (aG/3.0_dp/2.0_dp)/(n(1)*aG)+shiftFactor!*aG !
            !print*, "yShift= ", yShift, aG/3.0_dp
         end if
      else
         yShift = 0.0_dp
         !print*, "yShift= ", yShift, aG/3.0_dp
      end if
      call AtomsConstruct(RatBN,n(3:4),-phi*180.0_dp/pi,'graphene',-g, yShift)
      nAtBN = (n(3)**2 + n(3)*n(4) + n(4)**2)*2
      call MIO_Print('Atoms of BN: '//trim(num2str(nAtBN)),'atoms')
      nAt = nAtC + nAtBN
      call MIO_Allocate(Rat,[3,nAt],'Rat','atoms')
      Rat(:,1:nAtC) = RatC
      Rat(:,nAtC+1:nAt) = RatBN
      call MIO_InputParameter('createBLDomainBoundary',createBLDomainBoundary,.false.)
      m(1) = m(1)*sCell
      if (createBLDomainBoundary) then
          if (frac) then
              call AtomsSetCart()
              !print*, "Using cartesian coordinates"
          end if
          call MIO_InputParameter('BLDomainBoundaryWidth',boundaryWidth,60.0_dp)
          call MIO_InputParameter('BLDomainBoundaryTypeArmAA',boundaryTypeArmAA,.false.)
          call MIO_InputParameter('BLDomainBoundaryTypeArmSP',boundaryTypeArmSP,.false.)
          call MIO_InputParameter('BLDomainBoundaryNRShifts',numberOfShifts,1)
          limit0 = 0.0_dp
          limit1 = (m(1)*aG)*1.0_dp/2.0_dp
          limit2 = limit1 + boundaryWidth
          !print*, "BL Domain Boundary limits: ", limit1, limit2
             !do i=nAtC+1,nAt
             !   if (Rat(1,i).gt.limit2) then
             !       Rat(2,i) = Rat(2,i) - 1.42_dp * numberOfShifts
             !   else if (Rat(1,i).gt.limit1 .and. Rat(1,i).lt.limit2) then
             !       DBShiftRatio = (Rat(1,i)-limit1)/(limit2-limit1)
             !       Rat(2,i) = Rat(2,i) - DBShiftRatio * 1.42_dp * numberOfShifts
             !   end if
             !end do
          if (boundaryTypeArmAA) then ! check PRB 88 115409
             do i=1,nAtC
                if (Rat(1,i).gt.limit2) then
                    Rat(2,i) = Rat(2,i) - (aG/sqrt(3.0_dp) + (numberOfShifts-1)*aG/sqrt(3.0_dp))
                else if (Rat(1,i).gt.limit1 .and. Rat(1,i).lt.limit2) then
                    DBShiftRatio = (Rat(1,i)-limit1)/(limit2-limit1)
                    Rat(2,i) = Rat(2,i) - DBShiftRatio * (aG/sqrt(3.0_dp) + (numberOfShifts-1)*aG/sqrt(3.0_dp))
                end if
             end do
             do i=nAtC+1,nAt
                if (Rat(1,i).gt.limit2) then
                    Rat(2,i) = Rat(2,i) + (aG/sqrt(3.0_dp) + (numberOfShifts-1)*aG/sqrt(3.0_dp))
                else if (Rat(1,i).gt.limit1 .and. Rat(1,i).lt.limit2) then
                    DBShiftRatio = (Rat(1,i)-limit1)/(limit2-limit1)
                    Rat(2,i) = Rat(2,i) + DBShiftRatio * (aG/sqrt(3.0_dp) + (numberOfShifts-1)*aG/sqrt(3.0_dp))
                end if
             end do
          else if (boundaryTypeArmSP) then ! move top and bottom layer
             do i=1,nAtC
                if (Rat(1,i).gt.limit2) then
                    Rat(2,i) = Rat(2,i) + (0.71_dp + (numberOfShifts-1)*aG/sqrt(3.0_dp))
                else if (Rat(1,i).gt.limit1 .and. Rat(1,i).lt.limit2) then
                    DBShiftRatio = (Rat(1,i)-limit1)/(limit2-limit1)
                    Rat(2,i) = Rat(2,i) + DBShiftRatio * (0.71_dp + (numberOfShifts-1)*aG/sqrt(3.0_dp))
                end if
             end do
             do i=nAtC+1,nAt
                if (Rat(1,i).gt.limit2) then
                    Rat(2,i) = Rat(2,i) - (0.71_dp + (numberOfShifts-1)*aG/sqrt(3.0_dp))
                else if (Rat(1,i).gt.limit1 .and. Rat(1,i).lt.limit2) then
                    DBShiftRatio = (Rat(1,i)-limit1)/(limit2-limit1)
                    Rat(2,i) = Rat(2,i) - DBShiftRatio * (0.71_dp + (numberOfShifts-1)*aG/sqrt(3.0_dp))
                end if
             end do
          end if
          ggg = norm(ucell(:,1))
          ORZ = sqrt(3.0)/2.0 * ggg
          !print*, "ORZ = ", ORZ
          do i=1,nAt
            if (Rat(2,i).gt.ORZ) then
                Rat(2,i) = Rat(2,i) - ORZ
            else if (Rat(2,i).lt.0) then
                Rat(2,i) = Rat(2,i) + ORZ
            end if
          end do
      end if
      call MIO_Print('Total number of atoms: '//trim(num2str(nAt)),'atoms')
      call MIO_Allocate(Species,nAt,'Species','atoms')
      !$OMP PARALLEL DO
      do i=1,nAtC
         Species(i) = mod(i+1,2) + 1
      end do
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO
      do i=nAtC+1,nAt
         !Species(i) = mod(i+1,2) + 3
         Species(i) = mod(i+1,2) + 1
      end do
      !$OMP END PARALLEL DO
   else if (MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell')) then
      call MIO_InputParameter('CellSize',m(1),50)
      call MIO_InputParameter('MoireCellParameters',n,[0,0,0,0])
      call MIO_InputParameter('CellHeight',h,40.0_dp)
      call MIO_InputParameter('InterlayerDistance',dIntLay,3.22_dp) ! Ref. PRB 76, 73103
      g = n(1)**2 + n(2)**2 + n(1)*n(2)
      delta = sqrt(real(n(3)**2 + n(4)**2 + n(3)*n(4))/g)
      phi = acos((2.0_dp*n(1)*n(3)+2.0_dp*n(2)*n(4) + n(1)*n(4) + n(2)*n(3))/(2.0_dp*delta*g))
      !print*, "phi = ", phi*180_dp/pi
      if (delta<1.0_dp) then
         delta = 1.0_dp/delta
      else ! invert the indices of top and bottom layers
         i = n(1)
         n(1) = n(3)
         n(3) = i
         i = n(2)
         n(2) = n(4)
         n(4) = i
      end if
      n = n*sCell
      g = dIntLay/(2.0_dp*h)
      call AtomsConstruct(RatC,n(1:2),0.0_dp,'graphene',g, 0.0_dp)
      nAtC = (n(1)**2 + n(1)*n(2) + n(2)**2)*2
      call MIO_Print('Atoms of C:  '//trim(num2str(nAtC)),'atoms')
      call MIO_InputParameter('twistedBLAddShift',addShift,.false.)
      call MIO_InputParameter('BernalShift',Bernal,.false.)
      call MIO_InputParameter('bridgeShift',bridge,.false.)
      call MIO_InputParameter('BilayerShiftFactor',shiftFactor,0)
      !print*, "shiftinfo", n(1), shiftFactor
      if (addShift) then
         if (Bernal) then
            yShift = (aG/3.0_dp)/(n(1)*aG)+shiftFactor!*aG !
            !print*, "yShift= ", yShift, aG/3.0_dp
         else if (bridge) then
            yShift = (aG/3.0_dp/2.0_dp)/(n(1)*aG)+shiftFactor!*aG !
            !print*, "yShift= ", yShift, aG/3.0_dp
         end if
      else
         yShift = 0.0_dp
         !print*, "yShift= ", yShift, aG/3.0_dp
      end if
      call AtomsConstruct(RatBN,n(3:4),-phi*180.0_dp/pi,'graphene',-g, yShift)
      nAtBN = (n(3)**2 + n(3)*n(4) + n(4)**2)*2
      call MIO_Print('Atoms of BN: '//trim(num2str(nAtBN)),'atoms')
      nAt = nAtC + nAtBN
      call MIO_Allocate(Rat,[3,nAt],'Rat','atoms')
      Rat(:,1:nAtC) = RatC
      Rat(:,nAtC+1:nAt) = RatBN
      call MIO_Print('Total number of atoms: '//trim(num2str(nAt)),'atoms')
      call MIO_Allocate(Species,nAt,'Species','atoms')
      !$OMP PARALLEL DO
      do i=1,nAtC
         Species(i) = mod(i+1,2) + 1
      end do
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO
      do i=nAtC+1,nAt
         !Species(i) = mod(i+1,2) + 3
         Species(i) = mod(i+1,2) + 1
      end do
      !$OMP END PARALLEL DO
   else if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer')) then
      call MIO_InputParameter('CellSize',m(1),50)
      call MIO_InputParameter('CellHeight',h,40.0_dp)
      call MIO_InputParameter('InterlayerDistance',dIntLay,3.22_dp) ! Ref. PRB 76, 73103
      call MIO_InputParameter('twistedBilayerAngle',BLAngle,0.0_dp) ! Ref. PRB 76, 73103
      !print*, "BLAngle ", BLAngle
      call MIO_InputParameter('BilayerShiftFactor',shiftFactor,0) ! Ref. PRB 76, 73103
      m(1) = m(1)*sCell
      m(2) = 0
      !print*, m(1)
      !print*, m(2)
      g = dIntLay/(2.0_dp*h)
      !print*, "g= ", g
      call AtomsConstruct(RatC1,m,0.0_dp,'graphene',g, 0.0_dp)
      !print*, "RatC1", RatC1(:,10)
      nAtC1 = size(RatC1,2)
      call MIO_Print('Number of atoms in layer1: '//trim(num2str(nAtC1)),'atoms')
      !yShift = (aG/3.0_dp+shiftFactor*aG)/(m(1)*aG)
      if (MIO_StringComp(str,'MoireEncapsulatedBilayer')) then
          !yShift = (aG/3.0_dp)/(m(1)*aG)+shiftFactor*aG ! You can put shiftFactor to 1. This shifts the whole unit cell outside the first supercell. Consequently, when using the modulus and Rat()-1 routine, you shift everything back into this first supercell. If you don't proceed like this, one part is inside, the other is outside. When shifting the part outside back inside, you screw up the neighbor routine.
          yShift = (aG/3.0_dp)/(n(1)*aG)+shiftFactor
          !print*, "yShift= ", yShift, aG/3.0_dp
      end if
      if (MIO_StringComp(str,'TwistedBilayer')) then
          !yShift = shiftFactor*aG ! Remove the shift from normal bilayer graphene
          yShift = (aG/3.0_dp)/(m(1)*aG)+shiftFactor*aG
          !print*, "yShift= ", yShift, aG/3.0_dp
      end if
      !print*, "yShift= ", yShift, aG/3.0_dp
      !print*, "shiftFactor= ", shiftFactor
      !call AtomsConstruct(RatC2,m,BL Angle,'graphene',-g, yShift)
      ncell = m(1)**2 + m(1)*m(2) + m(2)**2
      !call AtomsConstructBasedOnFastNN(RatC2,m,BLAngle,'graphene',-g, yShift,ncell*2)
      call AtomsConstruct(RatC2,m,BLAngle,'graphene',-g, yShift)
      !print*, "RatC2", RatC2(:,10)
      nAtC2 = size(RatC2,2)
      call MIO_Print('Number of atoms in layer2: '//trim(num2str(nAtC2)),'atoms')
      sh = 1
      nAt = nAtC1 + nAtC2
      call MIO_Allocate(Rat,[3,nAt],'Rat','atoms')
      Rat(:,1:nAtC1) = RatC1
      Rat(:,nAtC1+1:nAt) = RatC2
      call MIO_Print('Total number of atoms: '//trim(num2str(nAt)),'atoms')
      call MIO_Allocate(Species,nAt,'Species','atoms')
      if (frac) then
          call AtomsSetCart()
          !print*, "hi1"
      end if
      !$OMP PARALLEL DO
      do i=1,nAt
         Species(i) = mod(i+1,2) + sh
         !if (i.gt.nAtC1) Rat(2,i) = Rat(2,i) + aG/sqrt(3.0_dp)
      end do
      !$OMP END PARALLEL DO
      if (.not. frac) call AtomsSetFrac()
      !!$OMP PARALLEL DO
      !do i=1,nAt
      !    if (Rat(2,i)>1.0_dp) print*, Rat(2,i), "hey"
      !    if (Rat(2,i)>1.0_dp) Rat(2,i) = mod(Rat(2,i) - 1.0_dp,1.0_dp)
      !end do
      !!$OMP END PARALLEL DO
   else if (MIO_StringComp(str,'TrilayerBasedOnMoireCell')) then ! Other twisted bilayers seem to have a wrong shift
      call MIO_InputParameter('CellSize',m(1),50)
      call MIO_InputParameter('MoireCellParameters',n,[0,0,0,0])
      call MIO_InputParameter('CellHeight',h,40.0_dp)
      call MIO_InputParameter('InterlayerDistance',dIntLay,3.22_dp) ! Ref. PRB 76, 73103
      g = n(1)**2 + n(2)**2 + n(1)*n(2)
      delta = sqrt(real(n(3)**2 + n(4)**2 + n(3)*n(4))/g)
      phi = acos((2.0_dp*n(1)*n(3)+2.0_dp*n(2)*n(4) + n(1)*n(4) + n(2)*n(3))/(2.0_dp*delta*g))
      !print*, "phi = ", phi*180_dp/pi
      if (delta<1.0_dp) then
         delta = 1.0_dp/delta
      else
         i = n(1)
         n(1) = n(3)
         n(3) = i
         i = n(2)
         n(2) = n(4)
         n(4) = i
      end if
      n = n*sCell
      g = dIntLay/(2.0_dp*h)
      call AtomsConstruct(RatC,n(1:2),0.0_dp,'graphene',g, 0.0_dp)
      nAtC = (n(1)**2 + n(1)*n(2) + n(2)**2)*2
      call MIO_Print('Atoms of C:  '//trim(num2str(nAtC)),'atoms')
      call MIO_InputParameter('TrilayerAddShift',addShift,.false.)
      call MIO_InputParameter('BernalShift',Bernal,.false.)
      call MIO_InputParameter('bridgeShift',bridge,.false.)
      call MIO_InputParameter('TrilayerShiftFactor',shiftFactor,0)
      !print*, "shiftinfo", n(1), shiftFactor
      if (addShift) then
         if (Bernal) then
            yShift = (aG/3.0_dp)/(n(1)*aG)+shiftFactor!*aG !
            !print*, "yShift= ", yShift, aG/3.0_dp
         else if (bridge) then
            yShift = (aG/3.0_dp/2.0_dp)/(n(1)*aG)+shiftFactor!*aG !
            !print*, "yShift= ", yShift, aG/3.0_dp
         end if
      else
         yShift = 0.0_dp
         !print*, "yShift= ", yShift, aG/3.0_dp
      end if
      call AtomsConstruct(RatC2,n(3:4),-phi*180.0_dp/pi,'graphene',-g, yShift)
      nAtC2 = (n(3)**2 + n(3)*n(4) + n(4)**2)*2
      call MIO_Print('Atoms of C2: '//trim(num2str(nAtC2)),'atoms')

      call AtomsConstruct(RatC3,n(3:4),-phi*180.0_dp/pi,'graphene',-g*3.0, 2.0*yShift)
      nAtC3 = (n(3)**2 + n(3)*n(4) + n(4)**2)*2
      call MIO_Print('Atoms of C3: '//trim(num2str(nAtC3)),'atoms')

      nAt = nAtC + nAtC2 + nAtC3
      call MIO_Allocate(Rat,[3,nAt],'Rat','atoms')
      Rat(:,1:nAtC) = RatC
      Rat(:,nAtC+1:nAtC+nAtC2) = RatC2
      Rat(:,nAtC+nAtC2+1:nAt) = RatC3
      call MIO_Print('Total number of atoms: '//trim(num2str(nAt)),'atoms')
      call MIO_Allocate(Species,nAt,'Species','atoms')
      !$OMP PARALLEL DO
      do i=1,nAtC
         Species(i) = mod(i+1,2) + 1
      end do
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO
      do i=nAtC+1,nAt
         !Species(i) = mod(i+1,2) + 3
         Species(i) = mod(i+1,2) + 1
      end do
      !$OMP END PARALLEL DO
   else if (MIO_StringComp(str,'Ribbons')) then
      call MIO_InputParameter('RibbonType',str,'Zigzag')
      call MIO_InputParameter('GrapheneWidth',nC,9)
      call MIO_InputParameter('BNWidth',nBN,9)
      nAt = (nC+nBN)*2
      call MIO_Allocate(Rat,[3,nAt],'Rat','atoms')
      call MIO_Allocate(Species,nAt,'Species','atoms')
      if (MIO_StringComp(str,'Zigzag')) then
         call AtomsRibbons('Z',nC,nBN,nAt,Rat,Species)
      else if (MIO_StringComp(str,'Armchair')) then
         call AtomsRibbons('A',nC,nBN,nAt,Rat,Species)
      else
         call MIO_Kill("Type of system not recognized. Check 'RibbonType'", &
           'atoms','AtomsPos')
      end if
   else if (MIO_StringComp(str,'ReadXYZ')) then
      call MIO_InputParameter('XYZFile',str,trim(prefix)//'.xyz')
      open(1,FILE=str,STATUS='old')
      call MIO_InputParameter('useSublatticeFile',useSublatticeFile,.true.)
      if (useSublatticeFile) then
         call MIO_Print('Let us read in the sublattice file, must be sorted after LAMMPS calculation','atoms')
         call MIO_InputParameter('SublatticeFile',str2,'sublatticesSorted.dat')
         open(2,FILE=str2,STATUS='old')
      end if
      call MIO_InputParameter('readInterlayerDistances',readInterlayerDistances,.false.)
      if (readInterlayerDistances) then
         call MIO_Print('Let us read in the interlayerDistance and displacements files, must be sorted after LAMMPS calculation','atoms')
         call MIO_InputParameter('interlayerDistanceFile',str3,'interlayerDistances.dat')
         !call MIO_InputParameter('displacementsFile',str4,'displacements.txt')
         open(3,FILE=str3,STATUS='old')
         !open(4,FILE=str4,STATUS='old')
      end if
      call MIO_InputParameter('GBNuseDisplacementFile',GBNuseDisplacementFile,.false.)
      call MIO_InputParameter('tBGuseDisplacementFile',tBGuseDisplacementFile,.false.)
      if (GBNuseDisplacementFile .or. tBGuseDisplacementFile) then ! Carr format here
         call MIO_InputParameter('displacementsFile',str4,'displacements.txt')
         open(4,FILE=str4,STATUS='old')
      end if
      do i=1,3
         read(1,*)
      end do
      !print*, "inside atoms.f90, ucell: ", ucell
      read(1,*) nAt
      call MIO_Allocate(RatTemp,[3,nAt],'RatTemp','atoms')
      call MIO_InputParameter('SuperCellX',SuperCellX,1)
      call MIO_InputParameter('SuperCellY',SuperCellY,1)
      if (SuperCellX .ne. SuperCellY) then
         call MIO_Allocate(Rat,[3,nAt*SuperCellX*SuperCellY],'Rat','atoms')
         call MIO_Allocate(Species,nAt*SuperCellX*SuperCellY,'Species','atoms')
      else if (sCell .ne. 1) then
         call MIO_Allocate(Rat,[3,nAt*sCell*sCell],'Rat','atoms')
         call MIO_Allocate(Species,nAt*sCell*sCell,'Species','atoms')
      else
         call MIO_Allocate(Rat,[3,nAt],'Rat','atoms')
         call MIO_Allocate(Species,nAt,'Species','atoms')
      end if
      call MIO_Allocate(SpeciesTemp,nAt,'SpeciesTemp','atoms')
      call MIO_Allocate(interlayerDistances,nAt*sCell*sCell,'interlayerDistances','atoms')
      call MIO_Allocate(displacements,[3,nAt*sCell*sCell],'displacements','atoms')
      call MIO_Allocate(displacements_b,[3,nAt*sCell*sCell],'displacements_b','atoms')
      call MIO_Allocate(displacements_t,[3,nAt*sCell*sCell],'displacements_t','atoms')
      call MIO_Allocate(layerIndex,nAt*sCell*sCell,'layerIndex','atoms')
      if (readInterlayerDistances) then  ! Konda, read in the displacement and interlayer distance file
         do i=1,nAt
            read(3,*) interlayerDistances(i)
            !read(4,*) (displacements(j,i),j=1,3) ! First is in x, second in y and third is the distance
         end do
      end if
      do i=1,nAt
         read(1,*) spc, (RatTemp(j,i),j=1,3)
         if (useSublatticeFile) then
            read(2,*) dontcare, SpeciesTemp(i)
         else
            if (spc=='C') then
               SpeciesTemp(i) = 1
            else if (spc=='B') then
               SpeciesTemp(i) = 3
            else if (spc=='N') then
               SpeciesTemp(i) = 4
            else
               call MIO_Kill("Wrong atomic number in file '"//trim(str)//"'", &
                 'atoms','AtomsPos')
            end if
         end if
      end do
      if (GBNuseDisplacementFile .or. tBGuseDisplacementFile) then ! Carr format here
         do i=1,nAt
            if (encapsulatedThreeLayers .and. (SpeciesTemp(i).eq.1 .or. SpeciesTemp(i).eq.2)) then
               read(4,*) displacements_b(1,i), displacements_b(2,i), displacements_b(3,i), displacements_t(1,i), displacements_t(2,i), displacements_t(3,i)
               !print*, "assigning bottom and top displacements"
            else
               read(4,*) (displacements(j,i),j=1,3) ! First is in x, second in y and third is the distance
            end if
         end do
      end if
      call MIO_InputParameter('invertDisplacements',invertDisplacements,.false.)
      if (invertDisplacements) then
         do i=1,nAt
             displacements(1,i) = -displacements(1,i)
             displacements(2,i) = -displacements(2,i)
         end do
      end if
      call MIO_InputParameter('readRigidXYZ',readRigidXYZ,.false.)
      if (readRigidXYZ) then
         call MIO_InputParameter('rigidPositions',str5,'generateInit.xyz')
         call MIO_Allocate(RatInit,[3,nAt*sCell*sCell],'RatInit','atoms')
         open(9,FILE=str5,STATUS='old')
         do i=1,4
            read(9,*)
         end do
         do i=1,nAt
            read(9,*) spc, (RatInit(j,i),j=1,3)
         end do
      end if
      close(1)
      close(9)
      !call MIO_InputParameter('SuperCellX',SuperCellX,1)
      !call MIO_InputParameter('SuperCellY',SuperCellY,1)
      !if (SuperCellX .ne. SuperCellY) then
      !   print*, "adjusting lattice vectors"
      !   ucell(:,1) = ucell(:,1)*SuperCellX
      !   ucell(:,2) = ucell(:,2)*SuperCellY
      !else if (sCell.ne.1) then
      !   print*, "adjusting lattice vectors"
      !   ucell(:,:2) = ucell(:,:2)*sCell
      !end if
      if (SuperCellX .ne. SuperCellY) then
          call MIO_Print('Different number of supercells in X and Y direction','atoms')
          kk = 0
          frac = .false.
          call AtomsSetFracTemp()
          !print*, "blube", RatTemp(2,:)
          call MIO_InputParameter('readInterlayerDistances',readInterlayerDistances,.false.)
          do ii=1,SuperCellX
             do jj=1,SuperCellY
                do i=1,nAt
                   Rat(1,kk*nAt+i) = RatTemp(1,i)+(ii-1)*1.0
                   Rat(2,kk*nAt+i) = RatTemp(2,i)+(jj-1)*1.0
                   Rat(3,kk*nAt+i) = RatTemp(3,i)
                   Species(kk*nAt+i) = SpeciesTemp(i)
                   if (readInterlayerDistances) then
                        interlayerDistances(kk*nAt+i) = interlayerDistances(i)
                   end if
                end do
                kk = kk + 1
             end do
          end do
          frac = .true.
          nAt = nAt*SuperCellX*SuperCellY
          call AtomsSetCart()
          ucell(:,1) = ucell(:,1)*SuperCellX
          ucell(:,2) = ucell(:,2)*SuperCellY
          vn = CrossProd(ucell(:,1),ucell(:,2))
          volume = dot_product(ucell(:,3),vn)
          area = norm(vn)
          rcell(:,1) = twopi*CrossProd(ucell(:,2),ucell(:,3))/volume
          rcell(:,2) = twopi*CrossProd(ucell(:,3),ucell(:,1))/volume
          rcell(:,3) = twopi*CrossProd(ucell(:,1),ucell(:,2))/volume
      else if (sCell.eq.1) then
          !print*, "only one supercell for READXYZ"
          !call MIO_Print('Different number of supercells in X and Y direction','atoms')
          Rat = RatTemp
          !print*, Rat
          Species = SpeciesTemp
          call MIO_InputParameter('readLayerIndex',readLayerIndex,.false.)
          if (readLayerIndex) then
             open(987,FILE='layerIndex.dat',STATUS='old')
             do i=1,nAt
                read(987,*) layerIndex(i)
             end do
             close(987)
          end if
      else
        call MIO_Print('We repeat the ReadXYZ several times using the supercell flag','atoms')
        kk = 0
        frac = .false.
        call AtomsSetFracTemp()
        !print*, "blube", RatTemp(2,:)
        call MIO_InputParameter('readInterlayerDistances',readInterlayerDistances,.false.)
        call MIO_InputParameter('readLayerIndex',readLayerIndex,.false.)
        call MIO_InputParameter('GBNuseDisplacementFile',GBNuseDisplacementFile,.false.)
        call MIO_InputParameter('tBGuseDisplacementFile',tBGuseDisplacementFile,.false.)
        if (readLayerIndex) then
           open(987,FILE='layerIndex.dat',STATUS='old')
           do i=1,nAt
              read(987,*) layerIndex(i)
           end do
           close(987)
        end if
        do ii=1,sCell
           do jj=1,sCell
              do i=1,nAt
                 Rat(1,kk*nAt+i) = RatTemp(1,i)+(ii-1)*1.0
                 Rat(2,kk*nAt+i) = RatTemp(2,i)+(jj-1)*1.0
                 Rat(3,kk*nAt+i) = RatTemp(3,i)
                 Species(kk*nAt+i) = SpeciesTemp(i)
                 if (readInterlayerDistances) then
                      interlayerDistances(kk*nAt+i) = interlayerDistances(i)
                 end if
                 if (readLayerIndex) then
                      layerIndex(kk*nAt+i) = layerIndex(i)
                 end if
                 if (GBNuseDisplacementFile .or. tBGuseDisplacementFile) then ! Carr format here
                      displacements(:,kk*nAt+i) = displacements(:,i)
                 end if
              end do
              kk = kk + 1
           end do
        end do
        !print*, "blabe", RatTemp(2,:)
        !print*, "blebe", Rat(2,:)
        !print*, "blobo", Rat(2,:)
        frac = .true.
        nAt = nAt*sCell*sCell
        call AtomsSetCart()
        ucell = ucell*sCell
        vn = CrossProd(ucell(:,1),ucell(:,2))
        volume = dot_product(ucell(:,3),vn)
        area = norm(vn)
        rcell(:,1) = twopi*CrossProd(ucell(:,2),ucell(:,3))/volume
        rcell(:,2) = twopi*CrossProd(ucell(:,3),ucell(:,1))/volume
        rcell(:,3) = twopi*CrossProd(ucell(:,1),ucell(:,2))/volume
      end if
      frac = .false.
   else
      call MIO_Kill('Type of system not recognized','atoms','AtomsPos')
   end if
   !call MIO_InputParameter('RandomStrain',randomStrain,.false.)
   !if (randomStrain) then
   !   if (frac) call AtomsSetCart()
   !   call GaussHeight()
   !end if
   call MIO_Allocate(indxDiv,nDiv+1,'indxDiv','atoms')
   call MIO_Allocate(procAt,nAt,'procAt','atoms')
   call MIO_Allocate(secAt,nAt,'secAt','atoms')
   nEl = 0
   do i=1,nAt
      if (Species(i)==1 .or. Species(i)==2) then
         nEl = nEl + 1.0_dp
      else if (Species(i)==4) then
         nEl = nEl + 2.0_dp
      end if
   end do
   call MIO_Print('Number of electrons: '//trim(num2str(nEl,1)),'atoms')
   !call MIO_InputParameter('AtomsOrderDeactivated',AtomsOrderDeactivated,.false.)
   !if (AtomsOrderDeactivated) then
   !    nDiv = 1
   !end if
   !print*, "nDiv= ", nDiv
   !if (nDiv /=1) then
      call AtomsOrder()
   !else
   !   indxDiv = (/1,nAt+1/)
   !end if
!#ifdef MPI
   call MIO_Allocate(indxNode,nProc+1,'indxNode','atoms')
   if (nProc>1) then
      do i=1,nProc
         indxNode(i) = indxDiv(sum(nTh(:i-1))+1)
      end do
   else
      indxNode = [1,nAt+1]
   end if
   indxNode(nProc+1) = nAt+1
   inode1 = indxNode(Node+1)
   inode2 = indxNode(Node+2)-1
   if (nDiv == 1) then
      ! in1 and in2 are THREADPRIVATE. Without a thread decomposition the first
      ! thread owns every atom; the others must get an EMPTY range, otherwise
      ! their unset bounds (0,0) make every "do i=in1,in2" loop inside a parallel
      ! region touch element 0.
      !$OMP PARALLEL
      in1 = 1
      in2 = 0
      !$OMP END PARALLEL
      in1 = 1
      in2 = nAt
   else
      !$OMP PARALLEL
      in1 = indxDiv(procID)
      in2 = indxDiv(procID+1) - 1
      !$OMP END PARALLEL
   end if
   !in1 = 1
   !in2 = nAt
   !inode1 = 1
   !inode2 = nAt
!#else
!   call MIO_Allocate(indxNode,2,'indxNode','atoms')
!   indxNode = [1,nAt+1]
!   inode1 = 1
!   inode2 = nAt
!   in1 = 1
!   in2 = nAt
!#endif /* MPI */
   call MIO_InputParameter('TypeOfSystem',str,'Graphene')
   if (MIO_StringComp(str,'Hybrid')) then
      call AtomsSetCart()
      call HybridGen(nAt,Rat,Species,in1,in2)
   end if
   call MIO_Allocate(layerIndex,nAt,'layerIndex','atoms')
   call MIO_Allocate(interlayerDistances,nAt,'interlayerDistance','atoms')
   if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer') &
        .or. MIO_StringComp(str,'TwistedBilayerBasedOnMoireCell') .or. MIO_StringComp(str,'MoireEncapsulatedBilayer') &
        .or. MIO_StringComp(str,'ReadXYZ') &
        .or. MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell')) then
      call AtomsSetCart()
      !print*, "lets assign the indices for the layers"
      call MIO_InputParameter('fourLayers',fourLayers,.false.)
      call MIO_InputParameter('fourLayersSandwiched',fourLayersSandwiched,.false.)
      call MIO_InputParameter('helicalTwistedMBM',helicalTwistedMBM,.false.)
      call MIO_InputParameter('fiveLayersSandwiched',fiveLayersSandwiched,.false.)
      call MIO_InputParameter('sixLayersSandwiched',sixLayersSandwiched,.false.)
      call MIO_InputParameter('sevenLayersSandwiched',sevenLayersSandwiched,.false.)
      call MIO_InputParameter('eightLayersSandwiched',eightLayersSandwiched,.false.)
      call MIO_InputParameter('tenLayersSandwiched',tenLayersSandwiched,.false.)
      call MIO_InputParameter('twentyLayersSandwiched',twentyLayersSandwiched,.false.)
      call MIO_InputParameter('threeLayers',threeLayers,.false.)
      call MIO_InputParameter('oneLayer',oneLayer,.false.)
      call MIO_InputParameter('twoLayers',twoLayers,.false.)
      call MIO_InputParameter('twoLayersZ1',twoLayersZ1,-3.3_dp)
      call MIO_InputParameter('twoLayersZ2',twoLayersZ2,3.3_dp)
      call MIO_InputParameter('threeLayersZ1',threeLayersZ1,-3.3_dp)
      call MIO_InputParameter('threeLayersZ2',threeLayersZ2,0.0_dp)
      call MIO_InputParameter('threeLayersZ3',threeLayersZ3,3.3_dp)
      call MIO_InputParameter('fourLayersZ1',fourLayersZ1,-3.3_dp)
      call MIO_InputParameter('fourLayersZ2',fourLayersZ2,-0.0_dp)
      call MIO_InputParameter('fourLayersZ3',fourLayersZ3,3.3_dp)
      call MIO_InputParameter('fourLayersZ4',fourLayersZ4,6.6_dp)
      call MIO_InputParameter('fiveLayersZ1',fiveLayersZ1,-3.3_dp)
      call MIO_InputParameter('fiveLayersZ2',fiveLayersZ2,-0.0_dp)
      call MIO_InputParameter('fiveLayersZ3',fiveLayersZ3,3.3_dp)
      call MIO_InputParameter('fiveLayersZ4',fiveLayersZ4,6.6_dp)
      call MIO_InputParameter('fiveLayersZ5',fiveLayersZ5,6.6_dp)
      call MIO_InputParameter('sixLayersZ1',sixLayersZ1,-3.3_dp)
      call MIO_InputParameter('sixLayersZ2',sixLayersZ2,-0.0_dp)
      call MIO_InputParameter('sixLayersZ3',sixLayersZ3,3.3_dp)
      call MIO_InputParameter('sixLayersZ4',sixLayersZ4,6.6_dp)
      call MIO_InputParameter('sixLayersZ5',sixLayersZ5,6.6_dp)
      call MIO_InputParameter('sixLayersZ6',sixLayersZ6,6.6_dp)
      call MIO_InputParameter('sevenLayersZ1',sevenLayersZ1,-3.3_dp)
      call MIO_InputParameter('sevenLayersZ2',sevenLayersZ2,-0.0_dp)
      call MIO_InputParameter('sevenLayersZ3',sevenLayersZ3,3.3_dp)
      call MIO_InputParameter('sevenLayersZ4',sevenLayersZ4,6.6_dp)
      call MIO_InputParameter('sevenLayersZ5',sevenLayersZ5,6.6_dp)
      call MIO_InputParameter('sevenLayersZ6',sevenLayersZ6,6.6_dp)
      call MIO_InputParameter('sevenLayersZ7',sevenLayersZ7,6.6_dp)
      call MIO_InputParameter('eightLayersZ1',eightLayersZ1,-3.3_dp)
      call MIO_InputParameter('eightLayersZ2',eightLayersZ2,-0.0_dp)
      call MIO_InputParameter('eightLayersZ3',eightLayersZ3,3.3_dp)
      call MIO_InputParameter('eightLayersZ4',eightLayersZ4,6.6_dp)
      call MIO_InputParameter('eightLayersZ5',eightLayersZ5,6.6_dp)
      call MIO_InputParameter('eightLayersZ6',eightLayersZ6,6.6_dp)
      call MIO_InputParameter('eightLayersZ7',eightLayersZ7,6.6_dp)
      call MIO_InputParameter('eightLayersZ8',eightLayersZ8,6.6_dp)
      call MIO_InputParameter('tenLayersZ1',tenLayersZ1,-3.3_dp)
      call MIO_InputParameter('tenLayersZ2',tenLayersZ2,-0.0_dp)
      call MIO_InputParameter('tenLayersZ3',tenLayersZ3,3.3_dp)
      call MIO_InputParameter('tenLayersZ4',tenLayersZ4,6.6_dp)
      call MIO_InputParameter('tenLayersZ5',tenLayersZ5,6.6_dp)
      call MIO_InputParameter('tenLayersZ6',tenLayersZ6,6.6_dp)
      call MIO_InputParameter('tenLayersZ7',tenLayersZ7,6.6_dp)
      call MIO_InputParameter('tenLayersZ8',tenLayersZ8,6.6_dp)
      call MIO_InputParameter('tenLayersZ9',tenLayersZ9,6.6_dp)
      call MIO_InputParameter('tenLayersZ10',tenLayersZ10,6.6_dp)
      call MIO_InputParameter('twentyLayersZ1',twentyLayersZ1,-3.3_dp)
      call MIO_InputParameter('twentyLayersZ2',twentyLayersZ2,-0.0_dp)
      call MIO_InputParameter('twentyLayersZ3',twentyLayersZ3,3.3_dp)
      call MIO_InputParameter('twentyLayersZ4',twentyLayersZ4,6.6_dp)
      call MIO_InputParameter('twentyLayersZ5',twentyLayersZ5,6.6_dp)
      call MIO_InputParameter('twentyLayersZ6',twentyLayersZ6,6.6_dp)
      call MIO_InputParameter('twentyLayersZ7',twentyLayersZ7,6.6_dp)
      call MIO_InputParameter('twentyLayersZ8',twentyLayersZ8,6.6_dp)
      call MIO_InputParameter('twentyLayersZ9',twentyLayersZ9,6.6_dp)
      call MIO_InputParameter('twentyLayersZ10',twentyLayersZ10,6.6_dp)
      call MIO_InputParameter('twentyLayersZ11',twentyLayersZ11,-3.3_dp)
      call MIO_InputParameter('twentyLayersZ12',twentyLayersZ12,-0.0_dp)
      call MIO_InputParameter('twentyLayersZ13',twentyLayersZ13,3.3_dp)
      call MIO_InputParameter('twentyLayersZ14',twentyLayersZ14,6.6_dp)
      call MIO_InputParameter('twentyLayersZ15',twentyLayersZ15,6.6_dp)
      call MIO_InputParameter('twentyLayersZ16',twentyLayersZ16,6.6_dp)
      call MIO_InputParameter('twentyLayersZ17',twentyLayersZ17,6.6_dp)
      call MIO_InputParameter('twentyLayersZ18',twentyLayersZ18,6.6_dp)
      call MIO_InputParameter('twentyLayersZ19',twentyLayersZ19,6.6_dp)
      call MIO_InputParameter('twentyLayersZ20',twentyLayersZ20,6.6_dp)
      call MIO_InputParameter('GBNtwoLayers',GBNtwoLayers,.false.)
      call MIO_InputParameter('t2GBN',t2GBN,.false.)
      call MIO_InputParameter('t2BG',t2BG,.false.)
      call MIO_InputParameter('t3BG',t3BG,.false.)
      call MIO_InputParameter('encapsulatedFourLayers',encapsulatedFourLayers,.false.)
      call MIO_InputParameter('readLayerIndex',readLayerIndex,.false.)
      !if (GBNtwoLayers) then
      !    if (Species(i).eq.1 .or. Species(i).eq.2) then
      !       layerIndex(i) = 1
      !    else if (Species(i).eq.3 .or. Species(i).eq.4) then
      !       layerIndex(i) = 2
      !    end if
      if (readLayerIndex) then
          call MIO_Print('We read the layer indices from an external file, make sure to provide it','atoms')
      else if (fourLayers .or. fourLayersSandwiched .or. helicalTwistedMBM .or. encapsulatedFourLayers) then
        do i=1,nAt
          if (Rat(3,i).lt.((fourLayersZ1+fourLayersZ2)/2.0_dp)) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          elseif (Rat(3,i).gt.((fourLayersZ1+fourLayersZ2)/2.0_dp) .and. Rat(3,i) .lt. (fourLayersZ2+fourLayersZ3)/2.0_dp) then
             layerIndex(i) = 2
          elseif (Rat(3,i).gt.((fourLayersZ2+fourLayersZ3)/2.0_dp) .and. Rat(3,i) .lt. (fourLayersZ3+fourLayersZ4)/2.0_dp) then
             layerIndex(i) = 3
          elseif (Rat(3,i).gt. (fourLayersZ3+fourLayersZ4)/2.0_dp) then
             layerIndex(i) = 4
          end if
        end do
      else if (fiveLayersSandwiched) then
        do i=1,nAt
          if (Rat(3,i).lt.((fiveLayersZ1+fiveLayersZ2)/2.0_dp)) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          elseif (Rat(3,i).gt.((fiveLayersZ1+fiveLayersZ2)/2.0_dp) .and. Rat(3,i) .lt. (fiveLayersZ2+fiveLayersZ3)/2.0_dp) then
             layerIndex(i) = 2
          elseif (Rat(3,i).gt.((fiveLayersZ2+fiveLayersZ3)/2.0_dp) .and. Rat(3,i) .lt. (fiveLayersZ3+fiveLayersZ4)/2.0_dp) then
             layerIndex(i) = 3
          elseif (Rat(3,i).gt.((fiveLayersZ3+fiveLayersZ4)/2.0_dp) .and. Rat(3,i) .lt. (fiveLayersZ4+fiveLayersZ5)/2.0_dp) then
             layerIndex(i) = 4
          elseif (Rat(3,i).gt. (fiveLayersZ4+fiveLayersZ5)/2.0_dp) then
             layerIndex(i) = 5
          end if
        end do
      else if (sixLayersSandwiched) then
        do i=1,nAt
          if (Rat(3,i).lt.((sixLayersZ1+sixLayersZ2)/2.0_dp)) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          elseif (Rat(3,i).gt.((sixLayersZ1+sixLayersZ2)/2.0_dp) .and. Rat(3,i) .lt. (sixLayersZ2+sixLayersZ3)/2.0_dp) then
             layerIndex(i) = 2
          elseif (Rat(3,i).gt.((sixLayersZ2+sixLayersZ3)/2.0_dp) .and. Rat(3,i) .lt. (sixLayersZ3+sixLayersZ4)/2.0_dp) then
             layerIndex(i) = 3
          elseif (Rat(3,i).gt.((sixLayersZ3+sixLayersZ4)/2.0_dp) .and. Rat(3,i) .lt. (sixLayersZ4+sixLayersZ5)/2.0_dp) then
             layerIndex(i) = 4
          elseif (Rat(3,i).gt.((sixLayersZ4+sixLayersZ5)/2.0_dp) .and. Rat(3,i) .lt. (sixLayersZ5+sixLayersZ6)/2.0_dp) then
             layerIndex(i) = 5
          elseif (Rat(3,i).gt. (sixLayersZ5+sixLayersZ6)/2.0_dp) then
             layerIndex(i) = 6
          end if
        end do
      else if (sevenLayersSandwiched) then
        do i=1,nAt
          if (Rat(3,i).lt.((sevenLayersZ1+sevenLayersZ2)/2.0_dp)) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          elseif (Rat(3,i).gt.((sevenLayersZ1+sevenLayersZ2)/2.0_dp) .and. Rat(3,i) .lt. (sevenLayersZ2+sevenLayersZ3)/2.0_dp) then
             layerIndex(i) = 2
          elseif (Rat(3,i).gt.((sevenLayersZ2+sevenLayersZ3)/2.0_dp) .and. Rat(3,i) .lt. (sevenLayersZ3+sevenLayersZ4)/2.0_dp) then
             layerIndex(i) = 3
          elseif (Rat(3,i).gt.((sevenLayersZ3+sevenLayersZ4)/2.0_dp) .and. Rat(3,i) .lt. (sevenLayersZ4+sevenLayersZ5)/2.0_dp) then
             layerIndex(i) = 4
          elseif (Rat(3,i).gt.((sevenLayersZ4+sevenLayersZ5)/2.0_dp) .and. Rat(3,i) .lt. (sevenLayersZ5+sevenLayersZ6)/2.0_dp) then
             layerIndex(i) = 5
          elseif (Rat(3,i).gt.((sevenLayersZ5+sevenLayersZ6)/2.0_dp) .and. Rat(3,i) .lt. (sevenLayersZ6+sevenLayersZ7)/2.0_dp) then
             layerIndex(i) = 6
          elseif (Rat(3,i).gt. (sevenLayersZ6+sevenLayersZ7)/2.0_dp) then
             layerIndex(i) = 7
          end if
        end do
      else if (eightLayersSandwiched) then
        do i=1,nAt
          if (Rat(3,i).lt.((eightLayersZ1+eightLayersZ2)/2.0_dp)) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          elseif (Rat(3,i).gt.((eightLayersZ1+eightLayersZ2)/2.0_dp) .and. Rat(3,i) .lt. (eightLayersZ2+eightLayersZ3)/2.0_dp) then
             layerIndex(i) = 2
          elseif (Rat(3,i).gt.((eightLayersZ2+eightLayersZ3)/2.0_dp) .and. Rat(3,i) .lt. (eightLayersZ3+eightLayersZ4)/2.0_dp) then
             layerIndex(i) = 3
          elseif (Rat(3,i).gt.((eightLayersZ3+eightLayersZ4)/2.0_dp) .and. Rat(3,i) .lt. (eightLayersZ4+eightLayersZ5)/2.0_dp) then
             layerIndex(i) = 4
          elseif (Rat(3,i).gt.((eightLayersZ4+eightLayersZ5)/2.0_dp) .and. Rat(3,i) .lt. (eightLayersZ5+eightLayersZ6)/2.0_dp) then
             layerIndex(i) = 5
          elseif (Rat(3,i).gt.((eightLayersZ5+eightLayersZ6)/2.0_dp) .and. Rat(3,i) .lt. (eightLayersZ6+eightLayersZ7)/2.0_dp) then
             layerIndex(i) = 6
          elseif (Rat(3,i).gt.((eightLayersZ6+eightLayersZ7)/2.0_dp) .and. Rat(3,i) .lt. (eightLayersZ7+eightLayersZ8)/2.0_dp) then
             layerIndex(i) = 7
          elseif (Rat(3,i).gt. (eightLayersZ7+eightLayersZ8)/2.0_dp) then
             layerIndex(i) = 8
          end if
        end do
      else if (tenLayersSandwiched) then
        do i=1,nAt
          if (Rat(3,i).lt.((tenLayersZ1+tenLayersZ2)/2.0_dp)) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          elseif (Rat(3,i).gt.((tenLayersZ1+tenLayersZ2)/2.0_dp) .and. Rat(3,i) .lt. (tenLayersZ2+tenLayersZ3)/2.0_dp) then
             layerIndex(i) = 2
          elseif (Rat(3,i).gt.((tenLayersZ2+tenLayersZ3)/2.0_dp) .and. Rat(3,i) .lt. (tenLayersZ3+tenLayersZ4)/2.0_dp) then
             layerIndex(i) = 3
          elseif (Rat(3,i).gt.((tenLayersZ3+tenLayersZ4)/2.0_dp) .and. Rat(3,i) .lt. (tenLayersZ4+tenLayersZ5)/2.0_dp) then
             layerIndex(i) = 4
          elseif (Rat(3,i).gt.((tenLayersZ4+tenLayersZ5)/2.0_dp) .and. Rat(3,i) .lt. (tenLayersZ5+tenLayersZ6)/2.0_dp) then
             layerIndex(i) = 5
          elseif (Rat(3,i).gt.((tenLayersZ5+tenLayersZ6)/2.0_dp) .and. Rat(3,i) .lt. (tenLayersZ6+tenLayersZ7)/2.0_dp) then
             layerIndex(i) = 6
          elseif (Rat(3,i).gt.((tenLayersZ6+tenLayersZ7)/2.0_dp) .and. Rat(3,i) .lt. (tenLayersZ7+tenLayersZ8)/2.0_dp) then
             layerIndex(i) = 7
          elseif (Rat(3,i).gt.((tenLayersZ7+tenLayersZ8)/2.0_dp) .and. Rat(3,i) .lt. (tenLayersZ8+tenLayersZ9)/2.0_dp) then
             layerIndex(i) = 8
          elseif (Rat(3,i).gt.((tenLayersZ8+tenLayersZ9)/2.0_dp) .and. Rat(3,i) .lt. (tenLayersZ9+tenLayersZ10)/2.0_dp) then
             layerIndex(i) = 9
          elseif (Rat(3,i).gt. (tenLayersZ9+tenLayersZ10)/2.0_dp) then
             layerIndex(i) = 10
          end if
        end do
      else if (twentyLayersSandwiched) then
        do i=1,nAt
          if (Rat(3,i).lt.((twentyLayersZ1+twentyLayersZ2)/2.0_dp)) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          elseif (Rat(3,i).gt.((twentyLayersZ1+twentyLayersZ2)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ2+twentyLayersZ3)/2.0_dp) then
             layerIndex(i) = 2
          elseif (Rat(3,i).gt.((twentyLayersZ2+twentyLayersZ3)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ3+twentyLayersZ4)/2.0_dp) then
             layerIndex(i) = 3
          elseif (Rat(3,i).gt.((twentyLayersZ3+twentyLayersZ4)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ4+twentyLayersZ5)/2.0_dp) then
             layerIndex(i) = 4
          elseif (Rat(3,i).gt.((twentyLayersZ4+twentyLayersZ5)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ5+twentyLayersZ6)/2.0_dp) then
             layerIndex(i) = 5
          elseif (Rat(3,i).gt.((twentyLayersZ5+twentyLayersZ6)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ6+twentyLayersZ7)/2.0_dp) then
             layerIndex(i) = 6
          elseif (Rat(3,i).gt.((twentyLayersZ6+twentyLayersZ7)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ7+twentyLayersZ8)/2.0_dp) then
             layerIndex(i) = 7
          elseif (Rat(3,i).gt.((twentyLayersZ7+twentyLayersZ8)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ8+twentyLayersZ9)/2.0_dp) then
             layerIndex(i) = 8
          elseif (Rat(3,i).gt.((twentyLayersZ8+twentyLayersZ9)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ9+twentyLayersZ10)/2.0_dp) then
             layerIndex(i) = 9
          elseif (Rat(3,i).gt.((twentyLayersZ9+twentyLayersZ10)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ10+twentyLayersZ11)/2.0_dp) then
             layerIndex(i) = 10
          elseif (Rat(3,i).gt.((twentyLayersZ10+twentyLayersZ11)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ11+twentyLayersZ12)/2.0_dp) then
             layerIndex(i) = 11
          elseif (Rat(3,i).gt.((twentyLayersZ11+twentyLayersZ12)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ12+twentyLayersZ13)/2.0_dp) then
             layerIndex(i) = 12
          elseif (Rat(3,i).gt.((twentyLayersZ12+twentyLayersZ13)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ13+twentyLayersZ14)/2.0_dp) then
             layerIndex(i) = 13
          elseif (Rat(3,i).gt.((twentyLayersZ13+twentyLayersZ14)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ14+twentyLayersZ15)/2.0_dp) then
             layerIndex(i) = 14
          elseif (Rat(3,i).gt.((twentyLayersZ14+twentyLayersZ15)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ15+twentyLayersZ16)/2.0_dp) then
             layerIndex(i) = 15
          elseif (Rat(3,i).gt.((twentyLayersZ15+twentyLayersZ16)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ16+twentyLayersZ17)/2.0_dp) then
             layerIndex(i) = 16
          elseif (Rat(3,i).gt.((twentyLayersZ16+twentyLayersZ17)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ17+twentyLayersZ18)/2.0_dp) then
             layerIndex(i) = 17
          elseif (Rat(3,i).gt.((twentyLayersZ17+twentyLayersZ18)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ18+twentyLayersZ19)/2.0_dp) then
             layerIndex(i) = 18
          elseif (Rat(3,i).gt.((twentyLayersZ18+twentyLayersZ19)/2.0_dp) .and. Rat(3,i) .lt. (twentyLayersZ19+twentyLayersZ20)/2.0_dp) then
             layerIndex(i) = 19
          elseif (Rat(3,i).gt. (twentyLayersZ19+twentyLayersZ20)/2.0_dp) then
             layerIndex(i) = 20
          end if
        end do
      else if (threeLayers) then
        do i=1,nAt
          if (Rat(3,i).lt.((threeLayersZ1+threeLayersZ2)/2.0_dp)) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          elseif (Rat(3,i).gt.((threeLayersZ1+threeLayersZ2)/2.0_dp) .and. Rat(3,i) .lt. ((threeLayersZ2+threeLayersZ3)/2.0_dp)) then
             layerIndex(i) = 2
          elseif ( Rat(3,i) .gt. ((threeLayersZ2+threeLayersZ3)/2.0_dp)) then
             layerIndex(i) = 3
          end if
        end do
      else if (twoLayers) then
        do i=1,nAt
          if (Rat(3,i).lt.((twoLayersZ1+twoLayersZ2)/2.0_dp)) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          else
             layerIndex(i) = 2
          end if
        end do
      else if (GBNtwoLayers) then
        do i=1,nAt
          if (Species(i).eq.1 .or. Species(i).eq.2) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          else
             layerIndex(i) = 2
          end if
        end do
      else if (t2GBN) then
        do i=1,nAt
          if (Rat(3,i).lt.((threeLayersZ1+threeLayersZ2)/2.0_dp)) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          elseif (Rat(3,i).gt.((threeLayersZ1+threeLayersZ2)/2.0_dp) .and. Rat(3,i) .lt. ((threeLayersZ2+threeLayersZ3)/2.0_dp)) then
             layerIndex(i) = 2
          elseif ( Rat(3,i) .gt. ((threeLayersZ2+threeLayersZ3)/2.0_dp)) then
             layerIndex(i) = 3
          end if
        end do
      else if (t3BG) then
        do i=1,nAt
          if (Rat(3,i).lt.((sixLayersZ1+sixLayersZ2)/2.0_dp)) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          elseif (Rat(3,i).gt.((sixLayersZ1+sixLayersZ2)/2.0_dp) .and. Rat(3,i) .lt. (sixLayersZ2+sixLayersZ3)/2.0_dp) then
             layerIndex(i) = 2
          elseif (Rat(3,i).gt.((sixLayersZ2+sixLayersZ3)/2.0_dp) .and. Rat(3,i) .lt. (sixLayersZ3+sixLayersZ4)/2.0_dp) then
             layerIndex(i) = 3
          elseif (Rat(3,i).gt.((sixLayersZ3+sixLayersZ4)/2.0_dp) .and. Rat(3,i) .lt. (sixLayersZ4+sixLayersZ5)/2.0_dp) then
             layerIndex(i) = 4
          elseif (Rat(3,i).gt.((sixLayersZ4+sixLayersZ5)/2.0_dp) .and. Rat(3,i) .lt. (sixLayersZ5+sixLayersZ6)/2.0_dp) then
             layerIndex(i) = 5
          elseif (Rat(3,i).gt. (sixLayersZ5+sixLayersZ6)/2.0_dp) then
             layerIndex(i) = 6
          end if
        end do
      else if (oneLayer) then
        do i=1,nAt
          layerIndex(i) = 1
        end do
      else
        do i=1,nAt
          if (Rat(3,i).lt.20.0_dp) then ! we first start by considering only the atoms in the bottom layer
             layerIndex(i) = 1
          else
             layerIndex(i) = 2
          end if
        end do
      end if
   else if (MIO_StringComp(str,'TrilayerBasedOnMoireCell')) then
      call AtomsSetCart()
      do i=1,nAt
        if (Rat(3,i).gt.20.0_dp) then ! we first start by considering only the atoms in the bottom layer
           layerIndex(i) = 3
        else if (Rat(3,i).gt.18.0_dp) then
           layerIndex(i) = 2
        else
           layerIndex(i) = 1
        end if
      end do
   !else
   !   do i=1,nAt
   !        layerIndex(i) = 1
   !   end do
   end if


#ifdef TIMER
   call MIO_TimerStop('atoms')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('AtomsPos',1)
#endif /* DEBUG */

end subroutine AtomsPos

subroutine AtomsConstruct(X,mm,angle,label,d,yShift)

   use constants,            only : pi

   real(dp), intent(out), pointer :: X(:,:)
   integer, intent(in) :: mm(2)
   real(dp), intent(in) :: angle,d, yShift
   character(*), intent(in) :: label

   real(dp), parameter :: bss(2,2) = reshape((/1.0_dp/3.0_dp,1.0_dp/3.0_dp, &
                          2.0_dp/3.0_dp,2.0_dp/3.0_dp/),(/2,2/))

   integer :: ncell, i1, i2, j, n1, n2, m(2,2)
   real(dp) :: ucell(3,2), Rcell(3,2), a
   real(dp) :: f(3)
   real(dp) :: basis(3,2)
   character(len=80) :: str

   if (mm(1)==0 .and. mm(2)==0) then
      call MIO_Kill('The supercell indices are zero: MoireCellParameters must be given for this type of system', &
        'atoms','AtomsLattice')
   end if
   ncell = mm(1)**2 + mm(1)*mm(2) + mm(2)**2
   call MIO_Allocate(X,[3,ncell*2],'X','atoms')
   a = angle*pi/180.0_dp
   !print*, "Angle= ", a, angle
   ucell(:,1) = [cos(a),sin(a),0.0_dp]
   ucell(:,2) = [cos(a+pi/3.0_dp),sin(a+pi/3.0_dp),0.0_dp]
   !print*, "ucell ", ucell
   m(:,1) = [mm(1),mm(2)]
   m(:,2) = [-mm(2),mm(1)+mm(2)]
   Rcell = 0.0_dp
   Rcell(:2,:2) = matmul(ucell(:2,:2),m)
   !print*, "Rcell ", Rcell
   basis = 0.0_dp
   basis(:2,:) = matmul(ucell(:2,:2),bss)
   !print*, "basis ", basis
   n1 = gcd(m(1,1),m(2,1))
   n2 = ncell/n1
   if (n2 > n1) then
      i1 = n1
      n1 = n2
      n2 = i1
   end if
   write(str,'(a,i0,a,i0)') 'Size of '//trim(label)//' layer:  ',n1,' x ',n2
   call MIO_Print(str,'atoms')
   !$OMP PARALLEL DO
   do i2=0,n2-1; do i1=0,n1-1
      do j=1,2
         X(:,2*(n1*i2+i1)+j) = i1*ucell(:,1) + i2*ucell(:,2) + basis(:,j)
      end do
   end do; end do
   !$OMP END PARALLEL DO
   !print*, "X", X(:,10)
   !yShift = yShift/2.0_dp
   ! We need to define a non rotated Rcell to be used in the loop below. Get the reduced coordinates of the rotated cell in the unrotated one
   !a = 0.0_dp
   !print*, "Angle= ", a, angle
   !ucell(:,1) = [cos(a),sin(a),0.0_dp]
   !ucell(:,2) = [cos(a+pi/3.0_dp),sin(a+pi/3.0_dp),0.0_dp]
   !print*, "ucell ", ucell
   !m(:,1) = [mm(1),mm(2)]
   !m(:,2) = [-mm(2),mm(1)+mm(2)]
   !Rcell = 0.0_dp
   !Rcell(:2,:2) = matmul(ucell(:2,:2),m)
   !$OMP PARALLEL DO PRIVATE(f)
   do i1=1,ncell*2
      !if (i1==10) print*, i1
      !if (i1==10) print*, "Xdiff", X(2,i1)-X(1,i1)
      f(2) = Rcell(1,1)*(X(2,i1)-X(1,i1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
      f(1) = (X(1,i1)-f(2)*Rcell(1,2))/Rcell(1,1)
      !if (i1==10) print*, "f", f
      !X(:,i1) = [mod(f(1)+yShift,1.0_dp), mod(f(2)+yShift,1.0_dp), 0.5_dp+d] ! first attemp, didn't work, if i remember well
      !X(:,i1) = [mod(f(1),1.0_dp)+yShift, mod(f(2),1.0_dp)+yShift, 0.5_dp+d]
      !if (X(1,i1)<0.0_dp) X(1,i1) = mod(X(1,i1) + 1.0_dp,1.0_dp) ! this is if f(1) was giving a negative value in previous line
      !if (X(2,i1)<0.0_dp) X(2,i1) = mod(X(2,i1) + 1.0_dp,1.0_dp)
      !X(1,i1) = X(1,i1)-yShift
      !X(2,i1) = X(2,i1)-yShift
      !if (X(1,i1)<0.0_dp) X(1,i1) = mod(X(1,i1) + 1.0_dp,1.0_dp)
      !if (X(1,i1)>1.0_dp) X(1,i1) = mod(X(1,i1) - 1.0_dp,1.0_dp)
      !if (X(2,i1)<0.0_dp) X(2,i1) = mod(X(2,i1) + 1.0_dp,1.0_dp)
      !if (X(2,i1)>1.0_dp) X(2,i1) = mod(X(2,i1) - 1.0_dp,1.0_dp)
      !if (i1==10) print*, "Xin0", X(:,i1)
      !print*, "yShiftShift", yShift
      X(:,i1) = [mod(f(1),1.0_dp)+yShift, mod(f(2),1.0_dp)+yShift, 0.5_dp+d]
      !if (i1==10) print*, "Xin1", X(:,i1)
      if (X(1,i1)<0.0_dp) X(1,i1) = mod(X(1,i1) + 1.0_dp,1.0_dp)
      !if (i1==10) print*, "Xin2", X(:,i1)
      if (X(2,i1)<0.0_dp) X(2,i1) = mod(X(2,i1) + 1.0_dp,1.0_dp)
      !if (i1==10) print*, "Xin3", X(:,i1)
      !X(1,i1) = X(1,i1)+yShift
      !X(2,i1) = X(2,i1)+yShift
      if (X(1,i1)>1.0_dp) X(1,i1) = mod(X(1,i1) - 1.0_dp,1.0_dp)
      !if (i1==10) print*, "Xin4", X(:,i1)
      if (X(2,i1)>1.0_dp) X(2,i1) = mod(X(2,i1) - 1.0_dp,1.0_dp)
      !if (i1==10) print*, "Xin5", X(:,i1)
      !if (X(1,i1)<0.0_dp) X(1,i1) = mod(X(1,i1) + 1.0_dp,1.0_dp)
      !if (X(2,i1)<0.0_dp) X(2,i1) = mod(X(2,i1) + 1.0_dp,1.0_dp)
   end do
   !$OMP END PARALLEL DO

end subroutine AtomsConstruct

subroutine AtomsConstructBasedOnFastNN(X,mm,angle,label,d,yShift, natoms)

   !use tbpar,                only : tbnn
   !use cell,                 only : aG, aBN,ucell
   use math
   use constants,            only : pi

   implicit none
   real(dp), intent(out), pointer :: X(:,:)
   integer, intent(in) :: mm(2)
   real(dp), intent(in) :: angle,d, yShift
   character(*), intent(in) :: label

   real(dp), parameter :: bss(2,2) = reshape((/1.0_dp/3.0_dp,1.0_dp/3.0_dp, &
                          2.0_dp/3.0_dp,2.0_dp/3.0_dp/),(/2,2/))

   integer :: ncell, i1, i2, j, n1, n2, m(2,2)
   real(dp) :: ucell(3,2), Rcell(3,2), a
   real(dp) :: f(3)
   real(dp) :: basis(3,2)
   character(len=80) :: str

   integer :: natoms
   !real(dp) :: x(natoms),y(natoms),z(natoms)
   !integer :: nn(natoms,maxnn)
   real(dp) ::  A1(3),A2(3)
   integer :: num0,num1,num2,num3,num4,num5
   integer :: i,k,l,n
   !real(dp) :: dx(natoms,maxnn),dy(natoms,maxnn),dz(natoms,maxnn)
   !real(dp), intent(out) :: dr(natoms,maxnn)

   real(dp) :: dx_t,dy_t,dz_t,d2_t
   real(dp) :: xmin,xmax,ymin,ymax,xcell,ycell,rho

   real(dp) :: xnew(9*natoms),ynew(9*natoms),znew(9*natoms)
   integer :: ixs(9*natoms),iys(9*natoms),vecino

   integer :: Nx,Ny,ix,iy
   integer, allocatable :: cells(:,:,:)
   integer :: NNcount,addx,addy

   real(dp), parameter :: rad(3) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp)]
   !character(len=50) :: str

   !integer :: ncell(3,9)
   real(dp) :: v(3)
   real(dp) :: rmax

   integer, pointer :: countN(:,:), cnt(:,:)
   integer :: np

   integer :: safetycounter


   if (mm(1)==0 .and. mm(2)==0) then
      call MIO_Kill('The supercell indices are zero: MoireCellParameters must be given for this type of system', &
        'atoms','AtomsLattice')
   end if
   ncell = mm(1)**2 + mm(1)*mm(2) + mm(2)**2
   call MIO_Allocate(X,[3,ncell*2],'X','atoms')
   a = angle*pi/180.0_dp
   !print*, "Angle= ", a, angle
   ucell(:,1) = [cos(a),sin(a),0.0_dp]
   ucell(:,2) = [cos(a+pi/3.0_dp),sin(a+pi/3.0_dp),0.0_dp]
   !print*, "ucell ", ucell
   m(:,1) = [mm(1),mm(2)]
   m(:,2) = [-mm(2),mm(1)+mm(2)]
   Rcell = 0.0_dp
   Rcell(:2,:2) = matmul(ucell(:2,:2),m)
   !print*, "Rcell ", Rcell
   basis = 0.0_dp
   basis(:2,:) = matmul(ucell(:2,:2),bss)
   !print*, "basis ", basis
   n1 = gcd(m(1,1),m(2,1))
   n2 = ncell/n1
   if (n2 > n1) then
      i1 = n1
      n1 = n2
      n2 = i1
   end if
   write(str,'(a,i0,a,i0)') 'Size of '//trim(label)//' layer:  ',n1,' x ',n2
   call MIO_Print(str,'atoms')
   !$OMP PARALLEL DO
   do i2=0,n2-1; do i1=0,n1-1
      do j=1,2
         X(:,2*(n1*i2+i1)+j) = i1*ucell(:,1) + i2*ucell(:,2) + basis(:,j)
      end do
   end do; end do
   !$OMP END PARALLEL DO

   !f(2) = Rcell(1,1)*(ynew(n)-xnew(n)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
   !f(1) = (X(1,n)-f(2)*Rcell(1,2))/Rcell(1,1)
   !if (f(1)<0 .or. f(1)>1 .or. f(2)<0 .or. f(2) > 0) cycle
   !i1 = i1+1
   !X(:,i1) = [mod(f(1),1.0_dp)+yShift, mod(f(2),1.0_dp)+yShift, 0.5_dp+d]

   !x(i) = X(:,1)
   !y(i) = X(:,2)
   !z(i) = Z(:,2)

   A1 = ucell(:,1)*n1
   A2 = ucell(:,2)*n2
   !print*, natoms, size(X(1,:))

   do i=1,natoms
     xnew(i)=X(1,i); ynew(i)=X(2,i); znew(i)=X(3,i)
   end do
   !print*, "hola1"
   do i=1,natoms
     n=i+natoms
     xnew(n)=X(1,i)+A1(1)       ! + A2(1)
     ynew(n)=X(2,i)+A1(2)       ! + A2(2)
     znew(n)=X(3,i)
   end do
   !print*, "hola2"
   do i=1,natoms
     n=i+2*natoms
     xnew(n)=X(1,i)+A1(1)+A2(1) !
     ynew(n)=X(2,i)+A1(2)+A2(2)
     znew(n)=X(3,i)
   end do
   !print*, "hola3"
   do i=1,natoms
     n=i+3*natoms
     xnew(n)=X(1,i)+A2(1)       ! + A1(1)
     ynew(n)=X(2,i)+A2(2)       ! + A1(2)
     znew(n)=X(3,i)
   end do
   !print*, "hola4"
   do i=1,natoms
     n=i+4*natoms
     xnew(n)=X(1,i)-A1(1)+A2(1) !
     ynew(n)=X(2,i)-A1(2)+A2(2)
     znew(n)=X(3,i)
   end do
   !print*, "hola5"
   do i=1,natoms
     n=i+5*natoms
     xnew(n)=X(1,i)-A1(1)       ! + A2(1)
     ynew(n)=X(2,i)-A1(2)       ! + A2(2)
     znew(n)=X(3,i)
   end do
   do i=1,natoms
     n=i+6*natoms
     xnew(n)=X(1,i)-A1(1)-A2(1)
     ynew(n)=X(2,i)-A1(2)-A2(2)
     znew(n)=X(3,i)
   end do
   do i=1,natoms
     n=i+7*natoms
     xnew(n)=X(1,i)-A2(1)       ! + A1(1)
     ynew(n)=X(2,i)-A2(2)       ! + A1(2)
     znew(n)=X(3,i)
   end do
   do i=1,natoms
     n=i+8*natoms
     xnew(n)=X(1,i)+A1(1)-A2(1)
     ynew(n)=X(2,i)+A1(2)-A2(2)
     znew(n)=X(3,i)
   end do
   !============================================================
   xmin = minval(xnew); xmax = maxval(xnew)
   ymin = minval(ynew); ymax = maxval(ynew)

   a = 0.0_dp
   !print*, "Angle= ", a, angle
   ucell(:,1) = [cos(a),sin(a),0.0_dp]
   ucell(:,2) = [cos(a+pi/3.0_dp),sin(a+pi/3.0_dp),0.0_dp]
   !print*, "ucell ", ucell
   m(:,1) = [mm(1),mm(2)]
   m(:,2) = [-mm(2),mm(1)+mm(2)]
   Rcell = 0.0_dp
   Rcell(:2,:2) = matmul(ucell(:2,:2),m)

   i1 = 0
   !!$OMP PARALLEL DO PRIVATE(f)
   do n = 1,9*natoms
      f(2) = Rcell(1,1)*(ynew(n)-xnew(n)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
      f(1) = (xnew(n)-f(2)*Rcell(1,2))/Rcell(1,1)
      !if (n.eq.54450) print*, "f =", f
      if (f(1).le.0.0_dp .or. f(1).ge.1.0_dp .or. f(2).le.0.0_dp .or. f(2).ge.1.0_dp) then
          !if (n.eq.54450) print*, "fa =", f
          cycle
      else
         !if (n.eq.54450) print*, "fb =", f
         i1 = i1+1
         if (i1.gt.natoms) then
           !print*, "cycled for ", i1
           cycle
         end if
         X(:,i1) = [mod(f(1),1.0_dp)+yShift, mod(f(2),1.0_dp)+yShift, 0.5_dp+d]
         if (X(1,i1)<0.0_dp) X(1,i1) = mod(X(1,i1) + 1.0_dp,1.0_dp)
         if (X(2,i1)<0.0_dp) X(2,i1) = mod(X(2,i1) + 1.0_dp,1.0_dp)
         if (X(1,i1)>1.0_dp) X(1,i1) = mod(X(1,i1) - 1.0_dp,1.0_dp)
         if (X(2,i1)>1.0_dp) X(2,i1) = mod(X(2,i1) - 1.0_dp,1.0_dp)
      end if
   end do
   !!$OMP END PARALLEL DO
   !do i=1,natoms
   !      X(:,i) = X(1,i)*ucell(:,1)*n1 + X(2,i)*ucell(:,2)*n2 + X(3,i)*1.0_dp
   !end do
   !open(4,FILE='posInt')
   !do i=1,natoms
   !   write(4,*) "C    ", X(1,i),X(2,i),X(3,i)
   !end do
   !close(4)

   !print*, "i1 = ", i1, natoms


end subroutine AtomsConstructBasedOnFastNN

subroutine AtomsConstructAsymm(X,mm,mm2,angle,label,d,yShift)

   use constants,            only : pi

   real(dp), intent(out), pointer :: X(:,:)
   integer, intent(in) :: mm(2),mm2(2)
   real(dp), intent(in) :: angle,d, yShift
   character(*), intent(in) :: label

   real(dp), parameter :: bss(2,2) = reshape((/1.0_dp/3.0_dp,1.0_dp/3.0_dp, &
                          2.0_dp/3.0_dp,2.0_dp/3.0_dp/),(/2,2/))

   integer :: ncell, i1, i2, j, n1, n2, m(2,2)
   real(dp) :: ucell(3,2), Rcell(3,2), a
   real(dp) :: f(3)
   real(dp) :: basis(3,2)
   character(len=80) :: str

   ! Modified
   !ncell = mm(1)**2 + mm(1)*mm(2) + mm(2)**2
   ncell = mm(1)*mm2(1) + mm(1)*mm(2) + mm(2)**2
   ! Modified

   call MIO_Allocate(X,[3,ncell*2],'X','atoms')
   a = angle*pi/180.0_dp
   ucell(:,1) = [cos(a),sin(a),0.0_dp]
   ucell(:,2) = [cos(a+pi/3.0_dp),sin(a+pi/3.0_dp),0.0_dp]

   m(:,1) = [mm(1),mm(2)]
   m(:,2) = [-mm2(2),mm2(1)+mm2(2)] ! I think this takes care of y direction (second column) of m matrix

   !m(:,1) = [mm(1),mm(2)]
   !m(:,2) = [-mm(2),mm(1)+mm(2)]

   !! Added
   !m2(:,1) = [mm2(1),mm2(2)]
   !m2(:,2) = [-mm2(2),mm2(1)+mm2(2)]
   !! Added

   Rcell = 0.0_dp
   Rcell(:2,:2) = matmul(ucell(:2,:2),m)
   basis = 0.0_dp
   basis(:2,:) = matmul(ucell(:2,:2),bss)
   n1 = gcd(m(1,1),m(2,1))
   n2 = ncell/n1
   if (n2 > n1) then
      i1 = n1
      n1 = n2
      n2 = i1
   end if
   write(str,'(a,i0,a,i0)') 'Size of '//trim(label)//' layer:  ',n1,' x ',n2
   call MIO_Print(str,'atoms')
   !$OMP PARALLEL DO
   do i2=0,n2-1; do i1=0,n1-1
      do j=1,2
         X(:,2*(n1*i2+i1)+j) = i1*ucell(:,1) + i2*ucell(:,2) + basis(:,j)
      end do
   end do; end do
   !$OMP END PARALLEL DO
   !print*, "d =", d
   !yShift = yShift/2.0_dp
   !$OMP PARALLEL DO PRIVATE(f)
   do i1=1,ncell*2
      f(2) = Rcell(1,1)*(X(2,i1)-X(1,i1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
      f(1) = (X(1,i1)-f(2)*Rcell(1,2))/Rcell(1,1)
      X(:,i1) = [mod(f(1),1.0_dp)+yShift, mod(f(2),1.0_dp)+yShift, 0.5_dp+d]
      if (X(1,i1)<0.0_dp) X(1,i1) = mod(X(1,i1) + 1.0_dp,1.0_dp)
      if (X(1,i1)>1.0_dp) X(1,i1) = mod(X(1,i1) - 1.0_dp,1.0_dp)
      if (X(2,i1)<0.0_dp) X(2,i1) = mod(X(2,i1) + 1.0_dp,1.0_dp)
      if (X(2,i1)>1.0_dp) X(2,i1) = mod(X(2,i1) - 1.0_dp,1.0_dp)
   end do
   !$OMP END PARALLEL DO

end subroutine AtomsConstructAsymm

subroutine AtomsRibbons(type,nC,nBN,n,X,Sp)

   real(dp), parameter :: v0Z(3)=[0.0_dp,1.0_dp/(2.0_dp*sqrt(3.0_dp)),0.0_dp]
   real(dp), parameter :: v0A(3)=[0.0_dp,1.0_dp/4.0_dp,0.0_dp]
   real(dp), parameter :: bssZ(3,2)= reshape([1.0_dp/4.0_dp,0.0_dp,0.0_dp, &
                                  3.0_dp/4.0_dp,1.0_dp/(2.0_dp*sqrt(3.0_dp)),0.0_dp], [3,2])
   real(dp), parameter :: bssA(3,2)= reshape([1.0_dp/sqrt(3.0_dp),0.0_dp,0.0_dp, &
                                  2.0_dp/sqrt(3.0_dp),0.0_dp,0.0_dp], [3,2])
   real(dp), parameter :: displxZ(3,2) = reshape([0.0_dp,0.0_dp,0.0_dp, &
                                  1.0_dp/2.0_dp,0.0_dp,0.0_dp], [3,2])
   real(dp), parameter :: displxA(3,2) = reshape([0.0_dp,0.0_dp,0.0_dp, &
                                  sqrt(3.0_dp)/2.0_dp,0.0_dp,0.0_dp], [3,2])
   real(dp), parameter :: displyZ(3) = [0.0_dp,sqrt(3.0_dp)/2.0_dp,0.0_dp]
   real(dp), parameter :: displyA(3) = [0.0_dp,1.0_dp/2.0_dp,0.0_dp]

   character, intent(in) :: type
   integer, intent(in) :: nC, nBN, n
   real(dp), intent(out) :: X(3,n)
   integer, intent(out) :: Sp(n)

   integer :: i

   if (type=='Z' .or. type=='z') then
      do i=1,nC
         X(:,2*i-1) = v0Z + (i-1.0_dp)*displyZ + displxZ(:,mod(i-1,2)+1) + bssZ(:,1)
         Sp(2*i-1) = 1
         X(:,2*i) = v0Z + (i-1.0_dp)*displyZ + displxZ(:,mod(i-1,2)+1) + bssZ(:,2)
         Sp(2*i) = 2
      end do
      do i=1,nBN
         X(:,2*(i+nC)-1) = v0Z + (i+nC-1.0_dp)*displyZ + displxZ(:,mod(i+nC-1,2)+1) + bssZ(:,1)
         Sp(2*(i+nC)-1) = 3
         X(:,2*(i+nC)) = v0Z + (i+nC-1.0_dp)*displyZ + displxZ(:,mod(i+nC-1,2)+1) + bssZ(:,2)
         Sp(2*(i+nC)) = 4
      end do
      do i=1,nAt
         X(1,i) = mod(X(1,i),1.0_dp)
         X(2,i) = 2.0_dp*X(2,i)/((nC+nBN)*sqrt(3.0_dp))
         X(3,i) = 0.50_dp
      end do
   else if (type=='A' .or. type=='a') then
      do i=1,nC
         X(:,2*i-1) = v0A + (i-1.0_dp)*displyA + displxA(:,mod(i-1,2)+1) + bssA(:,1)
         Sp(2*i-1) = 1
         X(:,2*i) = v0A + (i-1.0_dp)*displyA + displxA(:,mod(i-1,2)+1) + bssA(:,2)
         Sp(2*i) = 2
      end do
      do i=1,nBN
         X(:,2*(i+nC)-1) = v0A + (i+nC-1.0_dp)*displyA + displxA(:,mod(i+nC-1,2)+1) + bssA(:,1)
         Sp(2*(i+nC)-1) = 3
         X(:,2*(i+nC)) = v0A + (i+nC-1.0_dp)*displyA + displxA(:,mod(i+nC-1,2)+1) + bssA(:,2)
         Sp(2*(i+nC)) = 4
      end do
      do i=1,nAt
         X(1,i) = mod(X(1,i)/sqrt(3.0_dp),1.0_dp)
         X(2,i) = 2.0_dp*X(2,i)/real(nC+nBN,8)
         X(3,i) = 0.50_dp
      end do
   end if

end subroutine AtomsRibbons

subroutine AtomsSetCart()

   use cell,                 only : ucell

   integer :: i

#ifdef DEBUG
   call MIO_Debug('AtomsSetCart',0)
#endif /* DEBUG */

   if (frac) then
      !$OMP PARALLEL DO
      do i=1,nAt
         Rat(:,i) = Rat(1,i)*ucell(:,1) + Rat(2,i)*ucell(:,2) + Rat(3,i)*ucell(:,3)
      end do
      !$OMP END PARALLEL DO
      frac = .false.
   end if

#ifdef DEBUG
   call MIO_Debug('AtomsSetCart',1)
#endif /* DEBUG */

end subroutine AtomsSetCart

subroutine AtomsSetFrac()

   use cell,                 only : ucell

   integer :: i

#ifdef DEBUG
   call MIO_Debug('AtomsSetFrac',0)
#endif /* DEBUG */

   if (.not. frac) then
      !$OMP PARALLEL DO
      do i=1,nAt
         Rat(2,i) = ucell(1,1)*(Rat(2,i)-Rat(1,i)*ucell(2,1)/ucell(1,1))/(ucell(2,2)*ucell(1,1)-ucell(1,2)*ucell(2,1))
         Rat(1,i) = (Rat(1,i)-Rat(2,i)*ucell(1,2))/ucell(1,1)
         Rat(3,i) = Rat(3,i)/ucell(3,3)
      end do
      !$OMP END PARALLEL DO
      frac = .true.
   end if

#ifdef DEBUG
   call MIO_Debug('AtomsSetFrac',1)
#endif /* DEBUG */


end subroutine AtomsSetFrac

subroutine AtomsSetCartTemp()

   use cell,                 only : ucell

   integer :: i

#ifdef DEBUG
   call MIO_Debug('AtomsSetCartTemp',0)
#endif /* DEBUG */

   if (frac) then
      !$OMP PARALLEL DO
      do i=1,nAt
         RatTemp(:,i) = RatTemp(1,i)*ucell(:,1) + RatTemp(2,i)*ucell(:,2) + RatTemp(3,i)*ucell(:,3)
      end do
      !$OMP END PARALLEL DO
      frac = .false.
   end if

#ifdef DEBUG
   call MIO_Debug('AtomsSetCartTemp',1)
#endif /* DEBUG */

end subroutine AtomsSetCartTemp

subroutine AtomsSetFracTemp()

   use cell,                 only : ucell

   integer :: i

#ifdef DEBUG
   call MIO_Debug('AtomsSetFracTemp',0)
#endif /* DEBUG */

   if (.not. frac) then
      !$OMP PARALLEL DO
      do i=1,nAt
         RatTemp(2,i) = ucell(1,1)*(RatTemp(2,i)-RatTemp(1,i)*ucell(2,1)/ucell(1,1))/(ucell(2,2)*ucell(1,1)-ucell(1,2)*ucell(2,1))
         RatTemp(1,i) = (RatTemp(1,i)-RatTemp(2,i)*ucell(1,2))/ucell(1,1)
         RatTemp(3,i) = RatTemp(3,i)/ucell(3,3)
      end do
      !$OMP END PARALLEL DO
      frac = .true.
   end if

#ifdef DEBUG
   call MIO_Debug('AtomsSetFracTemp',1)
#endif /* DEBUG */


end subroutine AtomsSetFracTemp

subroutine AtomsOrder()

   use parallel,             only : xDiv, yDiv, procID, nDiv, nTh, xSubdiv, ySubdiv
   use cell,                 only : ucell

   integer :: i, ix, iy
   real(dp) :: dx, dy, a, d
   character(len=80) :: str
   logical :: AtomsOrderDeactivated
#ifdef MPI
   integer, pointer :: displ(:)
#endif /* MPI */

#ifdef DEBUG
   call MIO_Debug('AtomsOrder',0)
#endif /* DEBUG */
   call MIO_InputParameter('LatticeParameter',a,2.46_dp)
   a = a*2.0_dp
   d = norm(ucell(:,1))
   xSubdiv = max(1,floor(d/(a*xDiv)))
   ySubdiv = max(1,floor(d/(a*yDiv)))
   dx = 1.0_dp/xDiv
   dy = 1.0_dp/yDiv
   !$OMP PARALLEL DO PRIVATE(ix,iy)
   do i=1,nAt
      ix = Rat(1,i)/dx
      iy = Rat(2,i)/dy
      procAt(i) = iy*xDiv + ix + 1
      ix = xSubdiv*(Rat(1,i)-ix*dx)/dx
      iy = ySubdiv*(Rat(2,i)-iy*dy)/dy
      secAt(i) = xSubdiv*ySubdiv*(procAt(i)-1) + ix + iy*xSubdiv + 1
   end do
   !$OMP END PARALLEL DO
#ifdef DEBUG
   call MIO_Debug('AtomsSort',0)
#endif /* DEBUG */
   call MIO_InputParameter('AtomsOrderDeactivated',AtomsOrderDeactivated,.true.)
   if (.NOT. AtomsOrderDeactivated)then
      call AtomsSort(secAt,procAt,Rat,Species)
   end if
#ifdef DEBUG
   call MIO_Debug('AtomsSort',1)
#endif /* DEBUG */
   !$OMP PARALLEL PRIVATE(ix)
   !procAt = secAt ! garbage from Nicolas, don't add it again, unless I know
   !what I'm doing
   if (procID==1) then
      indxDiv(1) = 1
      !print*, "we are in the procID one loop in atoms.f90"
   else
      !print*, "we are NOT in the procID one loop in atoms.f90"
      ix = (procID-1)*nAt/nDiv + 1
      do
         if (procAt(ix)>=procID) then
            ix = ix - 1
         else
            if (procAt(ix+1)==procID) then
               indxDiv(procID) = ix+1
               exit
            else
               ix = ix + 1
            end if
         end if
      end do
   end if
   !$OMP END PARALLEL
   indxDiv(nDiv+1) = nAt+1
#ifdef MPI
   allocate(displ(nProc))
   !call MIO_Allocate(displ,nProc,'displ','atoms')
   displ(1) = 0
   do i=2,nProc
      displ(i) = displ(i-1) + nTh(i-1)
   end do
   !ix = sum(nTh(:Node)) + 1
   !iy = ix + numThreads - 1
   call MPIAllGatherV(MPI_IN_PLACE,0,indxDiv,nTh,displ,MPI_INTEGER)
   !call MIO_Deallocate(displ,'displ','atoms')
   deallocate(displ)
#endif /* MPI */
    !print*, "hi7"
   call MIO_Print('')
   call MIO_Print('Atoms distribution:','atoms')
   call MIO_Print('ProcessID   Num. of atoms           Indx1           Indx2','atoms')
   do ix=1,nDiv
      write(str,'(i8,i16,2i16)') ix, indxDiv(ix+1)-indxDiv(ix), indxDiv(ix), indxDiv(ix+1)-1
      call MIO_Print(str,'atoms')
   end do
   call MIO_Print('')
#ifdef DEBUG
   call MIO_Debug('AtomsOrder',1)
#endif /* DEBUG */

end subroutine AtomsOrder

recursive subroutine AtomsSort(array1,array2,array3,array4)

   integer, intent(inout) :: array1(:), array2(:), array4(:)
   real(dp), intent(inout) :: array3(:,:)

   integer :: iq

   if (size(array1) > 1) then
      call AtomsPart(array1,array2,array3,array4,iq)
      call AtomsSort(array1(:iq-1),array2(:iq-1),array3(:,:iq-1),array4(:iq-1))
      call AtomsSort(array1(iq:),array2(iq:),array3(:,iq:),array4(iq:))
   end if

end subroutine AtomsSort

subroutine AtomsPart(array1,array2,array3,array4,mrk)

   integer, intent(inout) :: array1(:), array2(:), array4(:)
   real(dp), intent(inout) :: array3(:,:)
   integer, intent(out) :: mrk

   integer :: i, j
   integer :: temp
   integer :: x      ! pivot point
   real(dp) :: t(3)

   x = array1(1)
   i= 0
   j= size(array1) + 1
   do
      j = j-1
      do
         if (array1(j) <= x) exit
         j = j-1
      end do
      i = i+1
      do
         if (array1(i) >= x) exit
         i = i+1
      end do
      if (i < j) then
         ! exchange array(i) and array(j)
         temp = array1(i)
         array1(i) = array1(j)
         array1(j) = temp
         temp = array2(i)
         array2(i) = array2(j)
         array2(j) = temp
         t = array3(:,i)
         array3(:,i) = array3(:,j)
         array3(:,j) = t
         temp = array4(i)
         array4(i) = array4(j)
         array4(j) = temp
      elseif (i == j) then
         mrk = i+1
         return
      else
         mrk = i
         return
     end if
  end  do

end subroutine AtomsPart

elemental function gcd(i1,i2)

    integer :: gcd
    integer, intent(in) :: i1,i2
    integer :: a,b

    a = i1
    b = i2
    if (a > b) then
        gcd = a
        a = b
        b = gcd
    end if

    do
      gcd = mod(a, b)
      if (gcd == 0) exit
      a = b
      b = gcd
    end do

    gcd = b

end function gcd

end module atoms
