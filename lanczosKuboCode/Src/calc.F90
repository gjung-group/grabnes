module calc

   use mio

   implicit none

   PRIVATE

   public :: CalcSelect

contains

subroutine CalcSelect()

   use ham,                  only : HamInit, HamHopping, HamOnSite, hopp
   use tbpar,                only : TBInit
   use magf,                 only : MagfInit, mBi, mBf, MagfValue, mStep, HaldPhaseInit, mPhii, mPhif, HaldPhaseValue, mPhiStep, magfield, Bmag
   use moireBLShift,         only : mSi, mSf, mSStep, moireBLShiftInit, moireBLShiftValue
   use atoms,                only : nAt, Rat, frac, AtomsSetCart
   use neigh,                only : Nneigh, neighCell, NList, maxNeigh, neighD
   use constants
   use cell,                 only : ucell, aG
   use math

   logical :: calcK, calcD, u, readDataFiles, calcT
   logical :: Frank
   integer :: mB, mS, mPhi, i, j, l, n
   integer :: n1, n2

   real(dp) :: hoppI(maxNeigh,nAt), hoppR(maxNeigh,nAt)
   real(dp) :: d, flux, phase, v1(3), v2(3)
   real(dp) :: diffx, diffy, mmphi

   call MIO_InputParameter('moireBLShift.Calc',u,.false.)
   if (u) then
       call TBInit()
       call HamInit()
       call moireBLShiftInit()
       do mS=mSi,mSf,1
          call moireBLShiftValue(mS)
          call HamOnSite()
          call MagfInit()
          call MIO_InputParameter('Kubo.Calc',calcK,.true.)
          call MIO_InputParameter('Diag.Calc',calcD,.false.)

          do mB=mBi,mBf, mStep
             call MagfValue(mB)
             call HamHopping()
             if(calcK) call CalcKubo()
             if (calcD) call CalcDiag()
          end do
       end do
   else
       call TBInit()
       call HamInit()
       call MagfInit()
       call HaldPhaseInit()
       call MIO_InputParameter('Kubo.Calc',calcK,.true.)
       call MIO_InputParameter('Diag.Calc',calcD,.false.)
       call MIO_InputParameter('ReadDataFiles',readDataFiles,.false.)
       call MIO_InputParameter('Tunn.Calc',calcT,.false.)
       call MIO_InputParameter('FrankMagneticField',Frank,.false.)

       do mB=mBi,mBf, mStep
          call MagfValue(mB)
          do mPhi=mPhii,mPhif,mPhiStep
             call HaldPhaseValue(mPhi)
             if (readDataFiles) then
                 call MIO_Print('Reading in the (already g0-renormalized) hopping terms from generate.s','ham')
                 open(222,FILE="generate.s",STATUS='old')
                 do i=1,nAt
                    do j=1,Nneigh(i)
                       read(222,*) hoppR(j,i), hoppI(j,i)
                       hopp(j,i) = (hoppR(j,i) + cmplx_i * hoppI(j,i))
                       !print*, i,j,Nneigh(i),hopp(j,i)
                    end do
                    !print*, hopp(:,i)
                 end do
                 close(222)
                 if (magfield) then
                    !if (Frank) then
                    !   if (frac) call AtomsSetCart()
                    !   do i=1,nAt
                    !      do j=1,Nneigh(i)
                    !         diffx = neighD(1,j,i)
                    !         diffy = neighD(2,j,i)
                    !         if (abs(diffy).lt.0.1d0 .and. diffx .gt. 0.1_dp) then
                    !             l = 1
                    !             n = 3
                    !         else if (abs(diffy).lt.0.1d0 .and. diffx .lt. -0.1_dp) then
                    !             l = 3
                    !             n = 1
                    !         else if (diffy.gt.0.1d0 .and. diffx .lt. -0.1_dp) then
                    !             l = 1
                    !             n = 2
                    !         else if (diffy.lt.-0.1d0 .and. diffx .gt. 0.1_dp) then
                    !             l = 2
                    !             n = 1
                    !         else if (diffy.lt.-0.1d0 .and. diffx .lt. -0.1_dp) then
                    !             l = 1
                    !             n = 2
                    !         else if (diffy.lt.0.1d0 .and. diffx .gt. 0.1_dp) then
                    !             l = 2
                    !             n = 1
                    !         else if (diffy.gt.0.1d0 .and. diffx .gt. 0.1_dp) then
                    !             l = 3
                    !             n = 4
                    !         else if (diffy.lt.-0.1d0 .and. diffx .lt. -0.1_dp) then
                    !             l = 4
                    !             n = 3
                    !         else if (diffy.lt.-0.1d0 .and. diffx .gt. 0.1_dp) then
                    !             l = 3
                    !             n = 4
                    !         else if (diffy.lt.0.1d0 .and. diffx .lt. -0.1_dp) then
                    !             l = 3
                    !             n = 4
                    !         end if
                    !         n2 = norm(ucell(:,2))/sqrt(3.0_dp)*3.0/2.0/aG ! n2 = N_y
                    !         n1 = norm(ucell(:,1))/(3.0*aG/sqrt(3.0_dp))/2.0_dp  ! n1 = Nx/2
                    !         if (l .eq. 1 .and. n .eq. 3) mmphi = -1.d0 / n2 * (j-1)
                    !         if (l .eq. 1 .and. n .eq. 2 .and. diffy .ge. 0.d0) mmphi = -1.d0 / n1 * (i-1)
                    !         if (l .eq. 1 .and. n .eq. 2 .and. diffy .lt. 0.d0) mmphi = 0.d0
                    !         if (l .eq. 2 .and. n .eq. 1 .and. diffy .ge. 0.d0) mmphi = 0.d0
                    !         if (l .eq. 2 .and. n .eq. 1 .and. diffy .lt. 0.d0) mmphi = 1.d0 / n1 * (i-1)
                    !         if (l .eq. 2 .and. n .eq. 0) mmphi= 1.d0 / n2 * (j-1)
                    !         if (l .eq. 3 .and. n .eq. 0 .and. diffy .ge. 0.d0) mmphi= 0.d0
                    !         if (l .eq. 3 .and. n .eq. 0 .and. diffy .lt. 0.d0) mmphi= 1.d0 / n1 * (i-1)+1.d0/ n1 / 2.d0
                    !         if (l .eq. 3 .and. n .eq. 1) mmphi= 1.d0 / n2 * (j-1)
                    !         if (l .eq. 4 .and. n .eq. 2) mmphi= -1.d0 / n2 * (j-1)
                    !         if (l .eq. 4 .and. n .eq. 3 .and. diffy .ge. 0.d0) mmphi= -1.d0 / n1 * (i-1)-1.d0/ n1 / 2.d0
                    !         if (l .eq. 4 .and. n .eq. 3 .and. diffy .lt. 0.d0) mmphi= 0.d0
                    !         !flux = Bmag*pi/fluxq
                    !         !print*, "flux, B, pi, fluxq", flux, Bmag, pi, fluxq
                    !         print*, "Magnetic field using Franks approach: ", mB / n2 / n1 * 39471.80806616257_dp, "T"
                    !         phase = -mmphi*2.0_dp*pi*mB
                    !         !print*, "phase= ", phase, v1, v2
                    !         hopp(j,i) = hopp(j,i)*exp(cmplx_i*phase)
                    !      end do
                    !   end do
                    !else
                       if (frac) call AtomsSetCart()
                       flux = Bmag*pi/fluxq
                       print*, "flux, B, pi, fluxq", flux, Bmag, pi, fluxq
                       do i=1,nAt
                          do j=1,Nneigh(i)
                             v1 = Rat(:,i) + Rat(:,NList(j,i))
                             v2 = neighCell(1,j,i)*ucell(:,1) + neighCell(2,j,i)*ucell(:,2) ! neighCell(1...) contains n1-m1, neighCell(2...) contains n2-m2 (see Cresti's notes)
                             v2 = CrossProd(v1,v2)
                             v1 = CrossProd(Rat(:,i),Rat(:,NList(j,i))) + v2 ! Implementation of third expersion in Cresti's notes
                             phase = flux*v1(3)/1.0d20
                             !print*, "phase= ", phase, v1, v2
                             hopp(j,i) = hopp(j,i)*exp(cmplx_i*phase)
                          end do
                       end do
                    !end if
                 end if
                 !print*, hopp
             else
                call HamHopping()
             end if
             if (calcK) call CalcKubo()
             if (calcD) call CalcDiag()
             if (calcT) call CalcTunn()
          end do
       end do
    end if

end subroutine CalcSelect

subroutine CalcKubo()

   use atoms,                only : inode1, inode2, Rat, in1, in2, frac, nAt
   use atoms,                only : AtomsSetFrac, AtomsSetCart, Species, layerIndex
   use kuboarrays
   use kubo,                 only : KuboInitWF, KuboInitWFPDOS, KuboDOS, KuboTEvol, KuboInitWFLayerDOS, KuboInitWFLayerAndSpeciesDOS
   use kubosubs,             only : KuboInterval, KuboCn
   use neigh,                only : NList2, NeighD
   use ham,                  only : H0, hopp
   use cell,                 only : sCell, aG
   use tbpar,                only : g0
   use name,                 only : prefix, sysname
#ifdef MPI
   use neigh,                only : rcvList, sndSz
#endif /* MPI */

   integer :: nRecurs, nT, nPol, nEn, nWr
   real(dp) :: dT, eps, Emin, Emax, ac, bc
   integer :: sz
   logical :: pol, onlydos
   real(dp) :: PDOSxmin, PDOSxmax, PDOSymin, PDOSymax, PDOSMoireSuperMoireLength
   integer :: i, j, ii
   logical :: PDOS, PDOSDist, PDOSMoireSC1, PDOSAtomList, PDOSByNumber, PDOSPNP, PDOSIgnoreLayer1and4
   logical :: PDOSMoireSuperMoire, PDOSMoireSuperMoire2, PDOSMoireSuperMoire3
   logical :: PDOSLayer, encapsulatedFourLayers, encapsulatedSixLayers, t3BG, PDOSLayerAndSpecies
   integer :: PDOSNumberOfAtoms, PDOSNumber, PDOSMinValue
   integer :: PDOSList(6)
   character(len=60) :: coordinates
   type(cl_file) :: file
   integer u, numberOfLayers
   integer flag1, flag2, n, cellSize, numberOfBNAtoms1, numberOfBNAtoms2, numberOfCAtoms, numberOfAtomsInLayer, maxSpecies
   real(dp) :: delta, limit1, aCC

   call MIO_InputParameter('RecursionNumber',nRecurs,700)
   call MIO_InputParameter('NumberofTimeSteps',nT,500)
   call MIO_InputParameter('TimeStep',dT,5.0_dp)
   call MIO_InputParameter('NumberofPolynomials',nPol,100)
   call MIO_InputParameter('NumberofEnergyPoints',nEn,1000)
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('WriteNumberofSteps',nWr,10)
   call MIO_InputParameter('EnergyMin',Emin,-4.0_dp)
   call MIO_InputParameter('EnergyMax',Emax,4.0_dp)
   eps = eps/g0
   Emin = Emin/g0
   Emax = Emax/g0

   call MIO_Allocate(a,nRecurs,'a','calc')
   call MIO_Allocate(b,nRecurs,'b','calc')
#ifdef MPI
   sz = 0
   if (associated(rcvList)) sz = size(rcvList)
#else
   sz = 0
#endif /* MPI */
   call MIO_Allocate(Psi,[inode1],[inode2+sz],'Psi','calc')
   call MIO_Allocate(Psin,[inode1],[inode2+sz],'Psin','calc')
   call MIO_Allocate(XpnPsi,[inode1],[inode2+sz],'XpnPsi','calc')
   call MIO_Allocate(XpnPsim1,[inode1],[inode2+sz],'XpnPsim1','calc')
   call MIO_Allocate(Psinm1,[inode1],[inode2+sz],'Psinm1','calc')
   call MIO_Allocate(ZUPsi,[inode1],[inode2+sz],'ZUPsi','calc')
   call MIO_Allocate(H,[inode1],[inode2],'H','calc')
#ifdef MPI
   if (nProc>1) then
      call MIO_Allocate(tempD,sndSz,'tempD','calc')
      call MIO_Allocate(tempZ,sndSz,'tempZ','calc')
   end if
#endif /* MPI */
   call MIO_InputParameter('PDOS',PDOS,.false.)
   if (PDOS) then
        call MIO_InputParameter('PDOSxmin',PDOSxmin,0.0_dp)
        call MIO_InputParameter('PDOSxmax',PDOSxmax,5.0_dp)
        call MIO_InputParameter('PDOSymin',PDOSymin,0.0_dp)
        call MIO_InputParameter('PDOSymax',PDOSymax,5.0_dp)
        call MIO_InputParameter('PDOSDist',PDOSDist,.false.)
        call MIO_InputParameter('PDOSLayer',PDOSLayer,.false.)
        call MIO_InputParameter('PDOSLayerAndSpecies',PDOSLayerAndSpecies,.false.)
        call MIO_InputParameter('PDOSMoireSC1',PDOSMoireSC1,.false.)
        call MIO_InputParameter('PDOSMoireSuperMoire',PDOSMoireSuperMoire,.true.)
        call MIO_InputParameter('PDOSMoireSuperMoire2',PDOSMoireSuperMoire2,.true.)
        call MIO_InputParameter('PDOSMoireSuperMoire3',PDOSMoireSuperMoire3,.true.)
        call MIO_InputParameter('PDOSMoireSuperMoireLength',PDOSMoireSuperMoireLength,135d0)
        call MIO_InputParameter('PDOSMinValue',PDOSMinValue,0)
        PDOSMoireSuperMoireLength = PDOSMoireSuperMoireLength/aG
        call MIO_InputParameter('CellSize',cellSize,50)
        call MIO_InputParameter('PDOSAtomList',PDOSAtomList,.false.)
        call MIO_InputParameter('PDOSByNumber',PDOSByNumber,.false.)
        call MIO_InputParameter('PDOSIgnoreLayer1and4',PDOSIgnoreLayer1and4,.false.)
        call MIO_InputParameter('numberOfLayers',numberOfLayers,2)
        call MIO_InputParameter('numberOfBNAtoms1',numberOfBNAtoms1,5832)
        call MIO_InputParameter('numberOfBNAtoms2',numberOfBNAtoms2,5832)
        call MIO_InputParameter('numberOfCAtoms',numberOfCAtoms,6050)
        call MIO_InputParameter('encapsulatedFourLayers',encapsulatedFourLayers,.false.)
        call MIO_InputParameter('encapsulatedSixLayers',encapsulatedSixLayers,.false.)
        call MIO_InputParameter('t3BG',t3BG,.false.)
        call MIO_InputParameter('PDOSPNP',PDOSPNP,.false.)
        !if (PDOSInterpolation) then
        !  !call PDOSByInterpolation()
        !  call MIO_InputParameter('PDOS_A',PDOS_A,1)
        !  call MIO_InputParameter('PDOS_Cp',PDOS_Cp,1)
        !  call MIO_InputParameter('PDOS_B',PDOS_B,1)
        !  call MIO_InputParameter('PDOS_Ap',PDOS_Ap,1)
        !  call MIO_InputParameter('PDOS_C',PDOS_C,1)
        !  call MIO_InputParameter('PDOS_Bp',PDOS_Bp,1)
        !  ! Calculate each of the 6 high symmetry sites
        !  ! PDOS_A, corresponds to AA site
        !  write(prefix,'(a4,a5,I8.8)') 'PDOS','_atom', PDOS_A
        !  call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'PDOS')
        !  prefix = './'//trim(prefix)
        !  call system('mkdir '//trim(prefix))
        !  coordinates = trim(prefix)//'/'//'coords'
        !  call file%Open(name=coordinates,serial=.true.)
        !  u = file%GetUnit()
        !  write(u,*) i, Rat(1,i), Rat(2,i), Rat(3,i)
        !  call file%Close()
        !  prefix = trim(prefix)//'/'//trim(sysname)
        !  call KuboInitWFPDOS(Psi,i)
        !  call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
        !  call KuboReadDOS(i, DOSA)
        !  do i=1, nenergy
        !     phi1 =
        !     C1 =
        !     C1p =
        !     phi2 =
        !     C2 =
        !     C0
        !     do i=inode1,inode2

        !     end do
        !  end do
        !else if (PDOSLayer) then
        if (PDOSLayer) then
          if (frac) call AtomsSetCart()
          do i=1, numberOfLayers
             write(prefix,'(a4,a6,I2.2)') 'PDOS','_layer', i
             call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'PDOS')
             prefix = './'//trim(prefix)
             call system('mkdir '//trim(prefix))
             !coordinates = trim(prefix)//'/'//'coords'
             !call file%Open(name=coordinates,serial=.true.)
             !u = file%GetUnit()
             !write(u,*) i, Rat(1,i), Rat(2,i), Rat(3,i)
             !call file%Close()
             prefix = trim(prefix)//'/'//trim(sysname)
             if ((i.eq.1 .or. i.eq.2 .or. i.eq.3 .or. i.eq.4  ) .and. t3BG) then
                numberOfAtomsInLayer = numberOfCAtoms*sCell*sCell
             else if (i.eq.5 .and. t3BG) then
!                print*, " t3BGPDOS: ", t3BG
                numberOfAtomsInLayer = numberOfCAtoms*sCell*sCell
             else if (i.eq.6 .and. t3BG) then
                numberOfAtomsInLayer = numberOfCAtoms*sCell*sCell
             else if ((i.eq.1) .and. encapsulatedFourLayers) then
                numberOfAtomsInLayer = numberOfBNAtoms1*sCell*sCell
             else if ((i.eq.4) .and. encapsulatedFourLayers) then
                numberOfAtomsInLayer = numberOfBNAtoms2*sCell*sCell
             else if ((i.eq.2 .or. i.eq.3) .and. encapsulatedFourLayers) then
                numberOfAtomsInLayer = numberOfCAtoms*sCell*sCell
             else
                numberOfAtomsInLayer = nAt/numberOfLayers
             end if
             call KuboInitWFLayerDOS(Psi,i,numberOfLayers,numberOfAtomsInLayer)
             call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
          end do
        else if (PDOSLayerAndSpecies) then
          if (frac) call AtomsSetCart()
          ! Find maximum species number in the system
          maxSpecies = 1
          do j=1, nAt
             if (Species(j).gt.maxSpecies) maxSpecies = Species(j)
          end do
          ! Process all layers and all species
          do i=1, numberOfLayers
             do ii=1, maxSpecies
                ! Dynamically count atoms in this layer and species
                numberOfAtomsInLayer = 0
                do j=1, nAt
                   if (layerIndex(j).eq.i .and. Species(j).eq.ii) then
                      numberOfAtomsInLayer = numberOfAtomsInLayer + 1
                   end if
                end do
                if (numberOfAtomsInLayer.eq.0) then
                   call MIO_Print("Warning: No atoms found for layer "//trim(num2str(i))//" species "//trim(num2str(ii))//", skipping",'PDOS')
                   cycle
                end if
                write(prefix,'(a4,a6,I2.2,a4,I2.2)') 'PDOS','_layer', i, '_sub',ii
                call MIO_Print("Writing files in folder '"//trim(prefix)//"' for "//trim(num2str(numberOfAtomsInLayer))//" atoms",'PDOS')
                prefix = './'//trim(prefix)
                call system('mkdir '//trim(prefix))
                !coordinates = trim(prefix)//'/'//'coords'
                !call file%Open(name=coordinates,serial=.true.)
                !u = file%GetUnit()
                !write(u,*) i, Rat(1,i), Rat(2,i), Rat(3,i)
                !call file%Close()
                prefix = trim(prefix)//'/'//trim(sysname)
                call KuboInitWFLayerAndSpeciesDOS(Psi,i,ii,numberOfLayers,numberOfAtomsInLayer)
                call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
             end do
          end do
        else if (PDOSDist) then
          if (frac) call AtomsSetCart()
          do i=inode1,inode2
            if(Rat(1,i) .gt. PDOSxmin .and. Rat(1,i) .lt. PDOSxmax .and.                        &
              Rat(2,i) .gt. PDOSymin .and. Rat(2,i) .lt. PDOSymax .and. i.ge.PDOSMinValue) then
                write(prefix,'(a4,a5,I8.8)') 'PDOS','_atom', i
                call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'PDOS')
                prefix = './'//trim(prefix)
                call system('mkdir '//trim(prefix))
                coordinates = trim(prefix)//'/'//'coords'
                call file%Open(name=coordinates,serial=.true.)
                u = file%GetUnit()
                write(u,*) i, Rat(1,i), Rat(2,i), Rat(3,i)
                call file%Close()
                prefix = trim(prefix)//'/'//trim(sysname)
                call KuboInitWFPDOS(Psi,i)
                call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
            end if
          end do
        else if (PDOSMoireSC1) then
          do i=inode1,inode2
            if (.not. frac) call AtomsSetFrac()
            if ((Rat(1,i) .lt. 1.0_dp/sCell) .and. (Rat(2,i) .lt. 1.0_dp/sCell) .and. i.ge.PDOSMinValue) then
               if (PDOSIgnoreLayer1and4 .and. (layerIndex(i).eq.1 .or. layerIndex(i).eq.4)) then
                   cycle
               else
                   if (frac) call AtomsSetCart()
                   write(prefix,'(a4,a5,I8.8)') 'PDOS','_atom', i
                   call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'PDOS')
                   prefix = './'//trim(prefix)
                   call system('mkdir '//trim(prefix))
                   coordinates = trim(prefix)//'/'//'coords'
                   call file%Open(name=coordinates,serial=.true.)
                   u = file%GetUnit()
                   write(u,*) i, Rat(1,i), Rat(2,i), Species(i), layerIndex(i)
                   call file%Close()
                   prefix = trim(prefix)//'/'//trim(sysname)
                   call KuboInitWFPDOS(Psi,i)
                   call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
               end if
            end if
          end do
        else if (PDOSMoireSuperMoire) then
          do i=inode1,inode2
            if (.not. frac) call AtomsSetFrac()
            if ((Rat(1,i) .lt. 1.0_dp/sCell/cellSize*PDOSMoireSuperMoireLength) .and. (Rat(2,i) .lt. 1.0_dp/scell/cellSize*PDOSMoireSuperMoireLength .and. Species(i).eq.1 .and. i.ge.PDOSMinValue)) then
               if (frac) call AtomsSetCart()
               write(prefix,'(a4,a5,I8.8)') 'PDOS','_atom', i
               call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'PDOS')
               prefix = './'//trim(prefix)
               call system('mkdir '//trim(prefix))
               coordinates = trim(prefix)//'/'//'coords'
               call file%Open(name=coordinates,serial=.true.)
               u = file%GetUnit()
               write(u,*) i, Rat(1,i), Rat(2,i), Rat(3,i)
               call file%Close()
               prefix = trim(prefix)//'/'//trim(sysname)
               call KuboInitWFPDOS(Psi,i)
               call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
            end if
          end do
        else if (PDOSMoireSuperMoire2) then
          do i=inode1,inode2
            if (.not. frac) call AtomsSetFrac()
            if ((Rat(1,i) .ge. 1.0_dp/sCell/cellSize*PDOSMoireSuperMoireLength) .and. ((Rat(2,i) .ge. 1.0_dp/scell/cellSize*PDOSMoireSuperMoireLength) .and. (Rat(1,i) .lt. 2.0_dp/sCell/cellSize*PDOSMoireSuperMoireLength) .and. (Rat(2,i) .lt. 2.0_dp/scell/cellSize*PDOSMoireSuperMoireLength) .and. (Species(i).eq.1) .and. (i.ge.PDOSMinValue))) then
               if (frac) call AtomsSetCart()
               write(prefix,'(a4,a5,I8.8)') 'PDOS','_atom', i
               call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'PDOS')
               prefix = './'//trim(prefix)
               call system('mkdir '//trim(prefix))
               coordinates = trim(prefix)//'/'//'coords'
               call file%Open(name=coordinates,serial=.true.)
               u = file%GetUnit()
               write(u,*) i, Rat(1,i), Rat(2,i), Rat(3,i)
               call file%Close()
               prefix = trim(prefix)//'/'//trim(sysname)
               call KuboInitWFPDOS(Psi,i)
               call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
            end if
          end do
        else if (PDOSMoireSuperMoire3) then
          do i=inode1,inode2
            if (.not. frac) call AtomsSetFrac()
            if ((Rat(1,i) .ge. 2.0_dp/sCell/cellSize*PDOSMoireSuperMoireLength) .and. ((Rat(2,i) .ge. 2.0_dp/scell/cellSize*PDOSMoireSuperMoireLength) .and. (Rat(1,i) .lt. 3.0_dp/sCell/cellSize*PDOSMoireSuperMoireLength) .and. (Rat(2,i) .lt. 3.0_dp/scell/cellSize*PDOSMoireSuperMoireLength) .and. (Species(i).eq.1) .and. (i.ge.PDOSMinValue))) then
               if (frac) call AtomsSetCart()
               write(prefix,'(a4,a5,I8.8)') 'PDOS','_atom', i
               call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'PDOS')
               prefix = './'//trim(prefix)
               call system('mkdir '//trim(prefix))
               coordinates = trim(prefix)//'/'//'coords'
               call file%Open(name=coordinates,serial=.true.)
               u = file%GetUnit()
               write(u,*) i, Rat(1,i), Rat(2,i), Rat(3,i)
               call file%Close()
               prefix = trim(prefix)//'/'//trim(sysname)
               call KuboInitWFPDOS(Psi,i)
               call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
            end if
          end do
        else if (PDOSAtomList) then
          !call MIO_InputParameter('PDOSNumberOfAtoms',PDOSNumberOfAtoms,10)
          !call MIO_Allocate(PDOSList,PDOSNumberOfAtoms,'PDOSList','calc')
          call MIO_InputParameter('PDOSList',PDOSList,[1,2,3,4,5,6])
          do j=1,size(PDOSList)
            i = PDOSList(j)
            write(prefix,'(a4,a5,I8.8)') 'PDOS','_atom', i
            call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'PDOS')
            prefix = './'//trim(prefix)
            call system('mkdir '//trim(prefix))
            coordinates = trim(prefix)//'/'//'coords'
            call file%Open(name=coordinates,serial=.true.)
            u = file%GetUnit()
            write(u,*) i, Rat(1,i), Rat(2,i), Rat(3,i)
            call file%Close()
            prefix = trim(prefix)//'/'//trim(sysname)
            call KuboInitWFPDOS(Psi,i)
            call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
          end do
        else if (PDOSByNumber) then
          call MIO_InputParameter('PDOSNumber',PDOSNumber,1)
          call KuboInitWFPDOS(Psi,PDOSNumber)
          call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
        else if (PDOSPNP) then
          if (frac) call AtomsSetCart()
          call MIO_InputParameter('CellSize',n,50)
          aCC = aG/sqrt(3.0_dp)
          limit1 = (n*sCell*aG)*1.0_dp/4.0_dp
          call MIO_InputParameter('PNPDelta',delta,10.0_dp)
          flag1 = 0
          flag2 = 0
          do i=inode1,inode2
            if ((Rat(1,i) .gt. limit1) .and. (Rat(1,i) .lt. limit1+aCC) .and. flag1.eq.0) then
               flag1 = 1
               write(prefix,'(a4,a5,I8.8)') 'PDOS','_atom', i
               call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'PDOS')
               call MIO_Print("limit1 = "//trim(num2str(limit1))//" and chosen atom X = "//trim(num2str(Rat(1,i))),'PDOS')
               call MIO_Print("atom i = "//trim(num2str(i)),'PDOS')
               prefix = './'//trim(prefix)
               call system('mkdir '//trim(prefix))
               coordinates = trim(prefix)//'/'//'coords'
               call file%Open(name=coordinates,serial=.true.)
               u = file%GetUnit()
               write(u,*) i, Rat(1,i), Rat(2,i), Rat(3,i)
               call file%Close()
               prefix = trim(prefix)//'/'//trim(sysname)
               call KuboInitWFPDOS(Psi,i)
               call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
            end if
            if ((Rat(1,i) .gt. limit1+delta*2.0_dp) .and. (Rat(1,i) .lt. limit1+delta*2.0_dp+aCC) .and. flag2.eq.0) then
               flag2 = 1
               write(prefix,'(a4,a5,I8.8)') 'PDOS','_atom', i
               call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'PDOS')
               call MIO_Print("limit1 = "//trim(num2str(limit1))//" and chosen atom X = "//trim(num2str(Rat(1,i))),'PDOS')
               prefix = './'//trim(prefix)
               call system('mkdir '//trim(prefix))
               coordinates = trim(prefix)//'/'//'coords'
               call file%Open(name=coordinates,serial=.true.)
               u = file%GetUnit()
               write(u,*) i, Rat(1,i), Rat(2,i), Rat(3,i)
               call file%Close()
               prefix = trim(prefix)//'/'//trim(sysname)
               call KuboInitWFPDOS(Psi,i)
               call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
            end if
          end do

        end if
   else
        call KuboInitWF(Psi)
        !print*, " nEn, Emin, Emax: ", nEn, Emin, Emax
        call KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList2,nRecurs,eps,nEn,Emin,Emax)
   end if
   !print*, "Psi =", Psi
   call MIO_InputParameter('Calculate.Polynomials',pol,.false.)
   call MIO_InputParameter('Calculate.OnlyDOS',onlydos,.false.)
   if (pol .or. .not. onlydos) then
      call KuboInterval(nRecurs,a,b,ac,bc)
      call MIO_Allocate(c,nPol,'c','calc')
      call KuboCn(dT,ac,bc,nPol,c)
   end if
   if (.not. onlydos) then
      call KuboTEvol(Psi,ZUPsi,Psin,Psinm1,XpnPsi,XpnPsim1,c,dT,nT,a,b,H,H0,hopp,NList2,NeighD, &
                    ac,bc,nEn,Emin,Emax,eps,nRecurs,nPol,nWr)
   end if

end subroutine CalcKubo

!subroutine CalcTunn()
!
!   use constants,             only : cmplx_i
!   use atoms,                only : Species, frac, AtomsSetCart
!   use neigh,                only : Nneigh, NList, NeighD
!   use atoms,                only : nAt, Species, layerIndex
!   use cell,                 only : aG, rcell
!   use tbpar,                only : g0
!   use ham,                  only : hopp
!   use math
!   use name,                 only : prefix
!
!   integer :: i, jj, j, u1, u2, u3, u4, numberOfMoires
!   real(dp) :: TAAR, TAAI, TABR, TABI
!   real(dp) :: dx,dy,bigKVec(3),bigKVecX,bigKVecY,bigKVecZ, dVec(3)
!   complex(dp) :: TAA, TAB
!   type(cl_file) :: file1, file2, file3, file4
!
!   if (frac) call AtomsSetCart()
!   call MIO_InputParameter('bigKVecX',bigKVecX,0.0_dp) ! give the coordinates of K like for the k-path, one by one
!   call MIO_InputParameter('bigKVecY',bigKVecY,0.0_dp)
!   call MIO_InputParameter('bigKVecZ',bigKVecZ,0.0_dp)
!   print*, bigKVecX
!   print*, bigKVecY
!   print*, bigKVecZ
!   bigKVec = [bigKVecX,bigKVecY,bigKVecZ]
!   bigKVec = bigKVec(1)*rcell(:,1) + bigKVec(2)*rcell(:,2) + bigKVec(3)*rcell(:,3)
!   print*, "bigKVec: ", bigKVec
!   call MIO_InputParameter('numberOfMoires',numberOfMoires,1)
!   call file1%Open(name=trim(prefix)//'.'//'TunnAAR',serial=.true.)
!   u1 = file1%GetUnit()
!   call file2%Open(name=trim(prefix)//'.'//'TunnAAI',serial=.true.)
!   u2 = file2%GetUnit()
!   call file3%Open(name=trim(prefix)//'.'//'TunnABR',serial=.true.)
!   u3 = file3%GetUnit()
!   call file4%Open(name=trim(prefix)//'.'//'TunnABI',serial=.true.)
!   u4 = file4%GetUnit()
!   do dx=0.0,aG*255.682280897,aG*255.682280897/11.0_dp
!      do dy=0,3.0_dp/sqrt(3.0_dp)*aG*255.682280897,3.0_dp/sqrt(3.0_dp)*aG*255.682280897/18.0_dp
!         dVec = [dx, dy, 0.0_dp]
!         !print*, "dVec"
!         !print*, dVec
!         !print*, dVec(1:2)
!         TAA = 0.0_dp
!         TAB = 0.0_dp
!         do i=1,nAt
!            do j=1,Nneigh(i)
!               jj = NList(j,i)
!               !print*, "dVec", dVec, NeighD(1:2,j,i)
!               if (layerIndex(i).eq.1 .and.  layerIndex(jj).eq.2) then ! focus on 1 layer only
!                  if (Species(i).eq.1 .and. Species(jj).eq. 1) then
!                     TAA = TAA + hopp(j,i)*g0 * exp(cmplx_i*dot_product(bigKVec(1:2),dVec(1:2)+NeighD(1:2,j,i)))
!                  else if (Species(i).eq.1 .and. Species(jj).eq. 2) then
!                     TAB = TAB + hopp(j,i)*g0 * exp(cmplx_i*dot_product(bigKVec(1:2),dVec(1:2)+NeighD(1:2,j,i)))
!                  !else if (Species(i).eq.2 .and. Species(jj).eq. 1) then
!                  !   TAA = TAA + hopp(j,i) * exp(dot_product(cmplx_i*bigKVec,NeighD(:,j,i)))
!                  !else if (Species(i).eq.2 .and. Species(jj).eq. 2) then
!                  !   TAB = TAB + hopp(j,i) * exp(dot_product(cmplx_i*bigKVec,NeighD(:,j,i)))
!                  endif
!                  !print*, "temp"
!               end if
!            end do
!         end do
!         !print*, "TAA: ", TAA
!         !print*, "TAB: ", TAB
!         TAAR = real(TAA)/numberOfMoires/(nAt/4.0_dp)
!         TAAI = imag(TAA)/numberOfMoires/(nAt/4.0_dp)
!         TABR = real(TAB)/numberOfMoires/(nAt/4.0_dp)
!         TABI = imag(TAB)/numberOfMoires/(nAt/4.0_dp)
!         !print*, TAAI
!         write(u1,*) dx, dy, TAAR
!         write(u2,*) dx, dy, TAAI
!         write(u3,*) dx, dy, TABR
!         write(u4,*) dx, dy, TABI
!      end do
!   end do
!   call file1%Close()
!   call file2%Close()
!   call file3%Close()
!   call file4%Close()
!
!end subroutine CalcTunn

subroutine CalcTunn()

   use constants,             only : cmplx_i
   use atoms,                only : Species, frac, AtomsSetCart
   use neigh,                only : Nneigh, NList, NeighD
   use atoms,                only : nAt, Species, layerIndex, Rat
   use cell,                 only : aG, rcell
   use tbpar,                only : g0
   use ham,                  only : hopp
   use math
   use name,                 only : prefix

   integer :: i, jj, j, u1, u2, u3, u4, numberOfMoires
   real(dp) :: TAAR, TAAI, TABR, TABI
   real(dp) :: dx,dy,bigKVec(3),bigKVecX,bigKVecY,bigKVecZ, dVec(3)
   complex(dp) :: TAA, TAB
   type(cl_file) :: file1, file2, file3, file4

   if (frac) call AtomsSetCart()
   call MIO_InputParameter('bigKVecX',bigKVecX,0.0_dp) ! give the coordinates of K like for the k-path, one by one
   call MIO_InputParameter('bigKVecY',bigKVecY,0.0_dp)
   call MIO_InputParameter('bigKVecZ',bigKVecZ,0.0_dp)
   print*, bigKVecX
   print*, bigKVecY
   print*, bigKVecZ
   bigKVec = [bigKVecX,bigKVecY,bigKVecZ]
   bigKVec = bigKVec(1)*rcell(:,1) + bigKVec(2)*rcell(:,2) + bigKVec(3)*rcell(:,3)
   print*, "bigKVec: ", bigKVec
   call MIO_InputParameter('numberOfMoires',numberOfMoires,1)
   call file1%Open(name=trim(prefix)//'.'//'TunnAAR',serial=.true.)
   u1 = file1%GetUnit()
   call file2%Open(name=trim(prefix)//'.'//'TunnAAI',serial=.true.)
   u2 = file2%GetUnit()
   call file3%Open(name=trim(prefix)//'.'//'TunnABR',serial=.true.)
   u3 = file3%GetUnit()
   call file4%Open(name=trim(prefix)//'.'//'TunnABI',serial=.true.)
   u4 = file4%GetUnit()
   !do dx=0.0,aG*255.682280897,aG*255.682280897/11.0_dp
   !   do dy=0,3.0_dp/sqrt(3.0_dp)*aG*255.682280897,3.0_dp/sqrt(3.0_dp)*aG*255.682280897/18.0_dp
   !      dVec = [dx, dy, 0.0_dp]
         !print*, "dVec"
         !print*, dVec
         !print*, dVec(1:2)
         do i=1,nAt
            if (Species(i).eq.1 .and. layerIndex(i).eq.1) then
               TAA = 0.0_dp
               TAB = 0.0_dp
               do j=1,Nneigh(i)
                  jj = NList(j,i)
                  !print*, "dVec", dVec, NeighD(1:2,j,i)
                  if (layerIndex(i).eq.1 .and.  layerIndex(jj).eq.2) then ! focus on 1 layer only
                     if (Species(i).eq.1 .and. Species(jj).eq. 1) then
                        TAA = TAA - hopp(j,i)*g0 * exp(cmplx_i*dot_product(bigKVec(1:2),NeighD(1:2,j,i)))
                     else if (Species(i).eq.1 .and. Species(jj).eq. 2) then
                        TAB = TAB - hopp(j,i)*g0 * exp(cmplx_i*dot_product(bigKVec(1:2),NeighD(1:2,j,i)))
                     !else if (Species(i).eq.2 .and. Species(jj).eq. 1) then
                     !   TAA = TAA + hopp(j,i) * exp(dot_product(cmplx_i*bigKVec,NeighD(:,j,i)))
                     !else if (Species(i).eq.2 .and. Species(jj).eq. 2) then
                     !   TAB = TAB + hopp(j,i) * exp(dot_product(cmplx_i*bigKVec,NeighD(:,j,i)))
                     endif
                     !print*, "temp"
                  end if
               end do
               TAAR = real(TAA)/numberOfMoires!/(nAt/4.0_dp)
               TAAI = imag(TAA)/numberOfMoires!/(nAt/4.0_dp)
               TABR = real(TAB)/numberOfMoires!/(nAt/4.0_dp)
               TABI = imag(TAB)/numberOfMoires!/(nAt/4.0_dp)
               !print*, TAAI
               write(u1,*) Rat(1,i), Rat(2,i), TAAR
               write(u2,*) Rat(1,i), Rat(2,i), TAAI
               write(u3,*) Rat(1,i), Rat(2,i), TABR
               write(u4,*) Rat(1,i), Rat(2,i), TABI
            end if
         end do
         !print*, "TAA: ", TAA
         !print*, "TAB: ", TAB
   !   end do
   !end do
   call file1%Close()
   call file2%Close()
   call file3%Close()
   call file4%Close()

end subroutine CalcTunn

subroutine CalcDiag()

   use diag,                 only : DiagInit, DiagDOS, DiagPDOS, DiagBands, Diag3DBands, DiagBandsG, DiagBandsAroundK, DiagSpectralFunction, DiagSpectralFunctionKGrid, DiagSpectralFunctionKGridInequivalent, DiagSpectralFunctionKGridInequivalentEnergyCut, DiagBandsRashba, DiagChern
   use diag,                 only : DiagSpectralFunctionKGridInequivalentEnergyCutNickDale
   use diag,                 only : gWeightsTAPW, gw_b1, gw_b2, gw_ntop, gWeightsHam
   use diag,                 only : layerWeightsTAPW, lw_b1, lw_b2
   use diag,                 only : DiagSpectralFunctionKGridInequivalent_v2, DiagSpectralFunctionKGridInequivalentEnergyCut_v2, moireAngle, gGridRotationAngle, tapwNG, calculateChern, nk_chern_x, nk_chern_y, fermi_energy, useTriangularTruncation, checkTAPWUnitary, physicalTwistAngle, useKprimeValley, tapwDebug, useRigidPositions, tapwLowdin, tapwBothValleys, tapwValleyDecouple, calculate3DTAPWBands, nk_3D_x, nk_3D_y, gammaCentred3D, socDebug, forceBlockTAPW
   use atoms,                only : nAt
   use ham,                  only : nspin
   use scf,                  only : SCFGetCharge, SCFInit, charge
#ifdef SEMICL
#define SEMICL_SECTION 1
#include "calc_semicl.inc"
#undef SEMICL_SECTION
#endif

   logical :: dos, bands, bands3D, spectral, grapheneUnitCell, aroundGrapheneK, spectralEnergyCut, RashbaCalculation, PDOS
   logical :: spectralEnergyCutNickDale, DirectBandVal, chern, enableSCF
#ifdef SEMICL
#define SEMICL_SECTION 2
#include "calc_semicl.inc"
#undef SEMICL_SECTION
#endif

   call MIO_InputParameter('Calculate.DOS',dos,.false.)
   call MIO_InputParameter('Calculate.PDOS',PDOS,.false.)
   call MIO_InputParameter('Calculate.Bands',bands,.false.)
   call MIO_InputParameter('Calculate.Chern',chern,.false.)
   call MIO_InputParameter('Diag.RashbaCalculation',RashbaCalculation,.false.)
   call MIO_InputParameter('Diag.MoireAngle',moireAngle,-33.004491598883078_dp)
   call MIO_InputParameter('Diag.GGridRotationAngle',gGridRotationAngle,0.0_dp)
   call MIO_InputParameter('Diag.N_G',tapwNG,180)
   call MIO_InputParameter('Diag.ChernCalculation',calculateChern,.false.)
   call MIO_InputParameter('Diag.ChernGridX',nk_chern_x,10)
   call MIO_InputParameter('Diag.ChernGridY',nk_chern_y,10)
   call MIO_InputParameter('Diag.Chern_EF',fermi_energy,0.0_dp)
   call MIO_InputParameter('Diag.TriangularTruncation',useTriangularTruncation,.false.)
   call MIO_InputParameter('Diag.CheckTAPWUnitary',checkTAPWUnitary,.false.)
   call MIO_InputParameter('Diag.PhysicalTwistAngle',physicalTwistAngle,1.08_dp)
   call MIO_InputParameter('Diag.UseKprimeValley',useKprimeValley,.false.)
   call MIO_InputParameter('Diag.TAPWDebug',tapwDebug,.false.)
   call MIO_InputParameter('Diag.socDebug',socDebug,.false.)
   call MIO_InputParameter('Diag.forceBlockTAPW',forceBlockTAPW,.false.)
   call MIO_InputParameter('Diag.UseRigidPositions',useRigidPositions,.false.)
   call MIO_InputParameter('Diag.TAPWLowdin',tapwLowdin,.false.)
   ! Both graphene valleys in one TAPW basis: G-shells around K AND K' (valley
   ! mixing enters through the exact projection).  Default .false. = one valley.
   call MIO_InputParameter('Diag.TAPWBothValleys',tapwBothValleys,.false.)
   ! Diagnostic: drop the K-K' block again (must reproduce the two one-valley runs).
   call MIO_InputParameter('Diag.TAPWValleyDecouple',tapwValleyDecouple,.false.)
   if (tapwValleyDecouple .and. .not. tapwBothValleys) &
        call MIO_Print('Diag.TAPWValleyDecouple ignored: needs Diag.TAPWBothValleys','calc')
   call MIO_InputParameter('Diag.Calculate3DTAPWBands',calculate3DTAPWBands,.false.)
   ! Default .true.: Monkhorst-Pack sits half a grid step off every high-
   ! symmetry point (measured: 0.492 steps from moire K'), which turns a band
   ! extremum there into a sampling artefact -- the G/hBN CellSize 55 CNP gap
   ! reads 7.14 meV on MP against 1.67 meV when K' is actually sampled.
   ! Contour AREAS are unaffected (median 1.3e-5 of A_BZ between the two).
   ! Set .false. to recover the old Monkhorst-Pack k-set.
   call MIO_InputParameter('Diag.GammaCentred3DGrid',gammaCentred3D,.true.)
   call MIO_InputParameter('Diag.3DBandsGridX',nk_3D_x,10)
   call MIO_InputParameter('Diag.3DBandsGridY',nk_3D_y,10)
#ifdef SEMICL
#define SEMICL_SECTION 3
#include "calc_semicl.inc"
#undef SEMICL_SECTION
#endif
   ! Plane-wave (G) composition of TAPW eigenstates -> <prefix>.GWeights (opt-in diagnostic: which moire
   ! harmonic couples two sheets at an avoided crossing).  Default .false.: a run that does not ask for it
   ! behaves exactly as before.
   call MIO_InputParameter('Diag.GWeights',gWeightsTAPW,.false.)
   call MIO_InputParameter('Diag.GWeightsBandMin',gw_b1,1)
   call MIO_InputParameter('Diag.GWeightsBandMax',gw_b2,1)
   call MIO_InputParameter('Diag.GWeightsTop',gw_ntop,8)
   call MIO_InputParameter('Diag.GWeightsHam',gWeightsHam,.false.)
   ! Layer / sublattice composition of TAPW eigenstates -> <prefix>.LayerWeights (opt-in; BandMin/Max 0 0 = all
   ! bands, as needed for layer charges).  Default .false.: a run that does not ask for it behaves exactly as before.
   call MIO_InputParameter('Diag.LayerWeights',layerWeightsTAPW,.false.)
   call MIO_InputParameter('Diag.LayerWeightsBandMin',lw_b1,0)
   call MIO_InputParameter('Diag.LayerWeightsBandMax',lw_b2,0)
   if (gWeightsTAPW .and. gWeightsHam) call MIO_Print('Diag.GWeightsHam: writing the projected TAPW Hamiltonian '// &
        'to <prefix>.TAPWHam','calc')
   if (gWeightsTAPW) call MIO_Print('Diag.GWeights: writing the plane-wave composition of bands '// &
        trim(num2str(gw_b1))//'..'//trim(num2str(gw_b2))//' (top '//trim(num2str(gw_ntop))//' G) to <prefix>.GWeights','calc')
   call MIO_InputParameter('Calculate.3DBands',bands3D,.false.)
   call MIO_InputParameter('Calculate.Spectral',spectral,.false.)
   call MIO_InputParameter('Calculate.SpectralEnergyCut',spectralEnergyCut,.false.)
   call MIO_InputParameter('Spectral.DirectBandVal',DirectBandVal,.false.)
   call MIO_InputParameter('Calculate.SpectralEnergyCutNickDale',spectralEnergyCutNickDale,.false.)
   call MIO_InputParameter('Bands.GrapheneUnitCell',grapheneUnitCell,.false.)
   call MIO_InputParameter('Bands.AroundGrapheneK',aroundGrapheneK,.false.)
   call MIO_InputParameter('EnableSCF',enableSCF,.true.)
   if (spectralEnergyCut .or. spectralEnergyCutNickDale .or. spectral .or. bands .or. dos .or. bands3D) then
      call DiagInit(nAt)
      if (nspin==2 .and. enableSCF) then
         ! The self-consistent Hubbard module (scf.F90) was never completed or
         ! tested and ends with a segmentation fault. Stop with an explanation.
         call MIO_Kill('Spin-polarized calculations with the self-consistent Hubbard module (EnableSCF, on by '// &
           'default) are not supported at present: the module is unfinished and untested. Set EnableSCF .false. '// &
           'to run a two-spin calculation without self-consistency.','calc','CalcDiag')
         call SCFGetCharge()
      else if (nspin==2 .and. .not. enableSCF) then
         ! Initialize arrays needed for two-spin calculations without SCF
         call MIO_Print('Running two-spin calculation without SCF','calc')
         ! Call SCFInit to allocate workspace arrays needed for diagonalization
         call SCFInit(nAt)
         ! Allocate charge array that diag.F90 needs (even without SCF)
         call MIO_Allocate(charge,[2,nAt],'charge','calc')
         ! Initialize charge array to default values
         charge = 0.0_dp
      end if
   end if
   if (dos) then
      if (PDOS .eqv. .true.) then
         call DiagPDOS()
      else
         call DiagDOS()
      end if
   end if
#ifdef SEMICL
#define SEMICL_SECTION 4
#include "calc_semicl.inc"
#undef SEMICL_SECTION
#endif
   if (chern) then
      call DiagInit(nAt)
      call DiagChern()
   end if
   if (bands) then
      if (grapheneUnitCell) then
         call DiagBandsG()
      else if (aroundGrapheneK) then
         call DiagBandsAroundK()
      else
         if (RashbaCalculation) then
            call DiagBandsRashba()
         else
            call DiagBands()
         end if
      end if
   end if
   if (bands3D) then
      if (grapheneUnitCell) then
         call DiagBandsG()
      else if (aroundGrapheneK) then
         call DiagBandsAroundK()
      else
         call Diag3DBands()
      end if
   end if
   if (spectral) then
      if (DirectBandVal) then
         call DiagSpectralFunctionKGridInequivalent_v2()
      else
         call DiagSpectralFunctionKGridInequivalent()
      end if
   else if (spectralEnergyCut) then
      if (DirectBandVal) then
         call DiagSpectralFunctionKGridInequivalentEnergyCut_v2()
      else
         call DiagSpectralFunctionKGridInequivalentEnergyCut()
      end if
   else if (spectralEnergyCutNickDale) then
      call DiagSpectralFunctionKGridInequivalentEnergyCutNickDale()
   end if

end subroutine CalcDiag

!subroutine PDOSByInterpolation()
!
!end subroutine PDOSByInterpolation

end module calc
