!> @file ham.F90
!! @brief Tight-binding Hamiltonian module for twisted bilayer graphene and related 2D materials.
!! @details Provides routines to initialize on-site terms and compute hopping integrals
!! for tight-binding Hamiltonians, including moiré strain, Zeeman terms, and interlayer
!! coupling models. Parallelization via OpenMP/MPI is supported.
!!
!! @section public_api Public API
!! - HamInit(): Initialize Hamiltonian parameters and neighbor lists
!! - HamOnSite(): Compute and add on-site energy terms
!! - HamHopping(): Compute hopping integrals (intra/inter-layer)
!!
!! @section deps Dependencies
!! - mio: I/O and utilities
!! - neigh: neighbor lists and distances
!! - atoms: atomic positions and species
!! - cell: unit cell parameters
!! - tbpar: tight-binding parameters
!!
!! @authors
!! - Rafael Martinez-Gordillo (original)
!! - Nicolas (modifications)
!! - Jinwoo (current version)
!!
!! @version 1.0
!! @date current

module ham
   use mio

   implicit none

   PRIVATE

!#ifdef MPI
!#define USE_MPI 1
!#include <sprng_f.h>
!#endif /* MPI */

   ! Hamiltonian matrices (public interface)
   real(dp), public, pointer :: H0(:)=>NULL()        ! On-site energy matrix
   real(dp), public, pointer :: Ho(:)=>NULL()        ! Alternative on-site energy matrix
   real(dp), public, pointer :: HABreal(:)=>NULL()   ! Real part of inter-layer coupling
   real(dp), public, pointer :: HABimag(:)=>NULL()   ! Imaginary part of inter-layer coupling
   real(dp), public, pointer :: H0Bottom(:)=>NULL()  ! Bottom layer on-site energies
   real(dp), public, pointer :: H0Top(:)=>NULL()     ! Top layer on-site energies
   complex(dp), public, pointer :: hopp(:,:)=>NULL() ! Hopping integral matrix

   ! Public subroutines
   public :: HamInit, HamOnSite, HamHopping

   ! Public functions
   public :: SOCEnabledForLayer, HaldaneEnabledForLayer

   ! System configuration flags (module-level state)
   logical, save :: moireStrain    ! Enable moiré pattern strain effects
   logical, save, public :: Zterm          ! Enable Zeeman term (magnetic field)
   logical, save :: PZterm         ! Enable pseudo-Zeeman term
   logical, save :: randomStrain   ! Enable random strain effects

   ! SOC terms flags
   logical, save, public :: IntrinsicSOCterm   ! Enable intrinsic SOC (λI)
   logical, save, public :: IsingSOCterm       ! Enable Ising SOC (λIsing) - spin-dependent
   logical, save, public :: RashbaSOCterm      ! Enable Rashba SOC (λR) - requires block Hamiltonian
   logical, save, public :: PIASOCterm         ! Enable pseudo-inversion asymmetry (λPIA)

   ! SOC layer control
   logical, save, public :: SOCLayerControl     ! Enable layer-specific SOC control
   integer, save, allocatable, public :: SOCLayersArray(:)  ! Array of layer indices for SOC
   integer, parameter :: maxSOCLayers = 50  ! Maximum number of layers that can have SOC

   ! Haldane layer control
   logical, save, public :: HaldaneLayerControl     ! Enable layer-specific Haldane control
   integer, save, allocatable, public :: HaldaneLayersArray(:)  ! Array of layer indices for Haldane
   integer, parameter :: maxHaldaneLayers = 50  ! Maximum number of layers that can have Haldane

   ! Physical parameters (module-level state)
   real(dp), save :: strFactor     ! Strain factor for moiré patterns
   real(dp), save :: dCC           ! Carbon-carbon bond distance
   real(dp), save, public :: gZeeman       ! Zeeman coupling strength
   real(dp), save :: gPZeeman      ! Pseudo-Zeeman coupling strength
   real(dp), save :: hStr          ! Strain height amplitude

   ! SOC parameters
   real(dp), save, public :: lambdaI       ! Intrinsic SOC parameter (λI)
   real(dp), save, public :: lambdaIsing   ! Ising SOC parameter (λIsing) - spin-dependent
   real(dp), save, public :: lambdaR       ! Rashba SOC parameter (λR)
   real(dp), save, public :: lambdaPIA     ! Pseudo-inversion asymmetry parameter (λPIA)

   ! Spin and system parameters
   integer, save :: spin           ! Spin index (1 or -1)
   integer, save, public :: nspin ! Number of spin components

   ! Special system features
   integer :: indexOfBigBubbleCenter  ! Index of the center of large bubbles

contains

!> @brief Initialize Hamiltonian parameters and neighbor lists.
!! @details Reads input parameters (strain, Zeeman terms, models), prepares neighbor
!! lists and internal arrays, and sets module state. On-site terms are computed in
!! HamOnSite; hopping is assembled elsewhere.
!! @par Side effects
!! - Initializes module-level variables
!! - Builds neighbor lists
!! - Configures strain and magnetic field parameters
!! @see HamOnSite
subroutine HamInit()

   use neigh,                only : NeighList, maxNeigh, Nneigh, neighCell, NList,neighD,Nradii, fastNNnotsquare, fastNNnotsquareSmall, fastNN, fastNNnotsquareBulk, fastNNnotsquareBulkSmall
   use neigh,                only : fastNNnotsquareNotRectangle
   use neigh,                only : NeighShellCutoff2
   use atoms,                only : inode1, inode2, Species, nAt, Rat, AtomsSetCart, layerIndex, in1, in2, frac, interlayerDistances
   use atoms,                only : AtomsSetFrac, displacements
   use interface,            only : InterfacePot2
   use cell,                 only : aG, ucell, aBN, sCell
   use tbpar,                only : g0, tbnn
   use gauss,                only : GaussHeight
   use math

   ! System configuration flags
   logical :: printBubble           ! Add a bubble to the system
   logical :: prnt                  ! Print flag for debugging
   logical :: createBLDomainBoundary ! Create bilayer domain boundary
   logical :: cutAtNN3              ! Cut at third nearest neighbor
   logical :: F2G2Model             ! Use F2G2 model
   logical :: readDataFiles         ! Read data from files instead of computing
   logical :: bulk, bulksmall, small ! System size flags

   ! Physical parameters
   integer, parameter :: numN(5) = [3,6,3,6,6]  ! Number of neighbors for each shell
   real(dp) :: aCC                  ! Carbon-carbon bond length
   real(dp) :: cutoff2, cutoff, cutoff2bis  ! Cutoff distances for neighbor search
   real(dp) :: A1(3), A2(3), A3    ! Lattice vectors
   real(dp) :: distFact, d, dist   ! Distance factors and distances
   real(dp) :: w, sigma, dist2, dy, dx  ! Strain and displacement parameters
   real(dp) :: limit0, limit1, limit2, limit3, DBShiftRatio, boundaryWidth  ! Domain boundary parameters

   ! Neighbor counting
   integer :: inplaneNeigh, outplaneNeigh  ! In-plane and out-of-plane neighbor counts
   integer :: nShells, modelShells         ! Intralayer shells searched / reached by the model
   integer :: maxnn                        ! Maximum number of neighbors

   ! Loop variables
   integer :: i, j, n
   integer :: kkk, kk

   ! OpenMP variables
   integer :: TID, OMP_GET_THREAD_NUM

   ! Temporary variables (consider removing if unused)
   logical :: l, ll, lll
   character(len=80) :: str
   integer :: tempSOCLayers(maxSOCLayers), nSOCLayersRead
   character(len=500) :: SOCLayersString
   integer :: tempHaldaneLayers(maxHaldaneLayers), nHaldaneLayersRead
   character(len=500) :: HaldaneLayersString


   ! Get OpenMP thread ID for potential thread-specific operations
   !$OMP PARALLEL PRIVATE(TID)
   TID = OMP_GET_THREAD_NUM()
   !$OMP END PARALLEL
#ifdef DEBUG
   call MIO_Debug('HamInit',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('ham')
#endif /* TIMER */

   call MIO_InputParameter('ZeemanTerm',Zterm,.false.)
   call MIO_InputParameter('PseudoZeemanTerm',PZterm,.false.)
   if (Zterm) then
      call MIO_InputParameter('Spin',spin,1)
      if (spin==2) then
         spin=-1
      else if (spin/=1 .and. spin/=-1) then
         call MIO_Kill('Wrong spin index','ham','HamInit')
      end if
      call MIO_InputParameter('ZeemanFactor',gZeeman,0.00033_dp)
      gZeeman = gZeeman/g0
      call MIO_Print('Calculation for spin '//trim(num2str(spin))//'/2','ham')
   end if
   if (PZterm) then
      call MIO_InputParameter('PseudoZeemanFactor',gPZeeman,0.00016_dp)
      gPZeeman = gPZeeman/g0
   end if
   call MIO_InputParameter('SpinPolarized',l,.false.)
   if (l) then
      nspin = 2
   else
      nspin = 1
   end if

   ! Read SOC term flags and parameters
   call MIO_InputParameter('IntrinsicSOCterm',IntrinsicSOCterm,.false.)
   if (IntrinsicSOCterm) then
      call MIO_InputParameter('LambdaI',lambdaI,0.0_dp)
      lambdaI = lambdaI / g0  ! Normalize by g0 for consistency with other SOC terms
      call MIO_Print('Intrinsic SOC enabled with λI = '//trim(num2str(lambdaI)),'ham')
   end if

   call MIO_InputParameter('IsingSOCterm',IsingSOCterm,.false.)
   if (IsingSOCterm) then
      call MIO_InputParameter('LambdaIsing',lambdaIsing,0.0_dp)
      lambdaIsing = lambdaIsing / g0  ! Normalize by g0 for consistency with other SOC terms
      call MIO_Print('Ising SOC enabled with λIsing = '//trim(num2str(lambdaIsing)),'ham')
   end if

   call MIO_InputParameter('RashbaSOCterm',RashbaSOCterm,.false.)
   if (RashbaSOCterm) then
      call MIO_InputParameter('LambdaR',lambdaR,0.0_dp)
      lambdaR = lambdaR / g0  ! Normalize by g0 for consistency with other SOC terms
      call MIO_Print('Rashba SOC enabled with λR = '//trim(num2str(lambdaR)),'ham')
   end if

   call MIO_InputParameter('PIASOCterm',PIASOCterm,.false.)
   if (PIASOCterm) then
      call MIO_InputParameter('LambdaPIA',lambdaPIA,0.0_dp)
      lambdaPIA = lambdaPIA / g0  ! Normalize by g0 for consistency with other SOC terms
      call MIO_Print('Pseudo-inversion asymmetry enabled with λPIA = '//trim(num2str(lambdaPIA)),'ham')
   end if

   ! Layer-specific SOC control
   call MIO_InputParameter('SOCLayerControl',SOCLayerControl,.false.)
   if (SOCLayerControl) then
      ! Use new variable-length input parameter function
      call MIO_InputParameterVariable('SOCLayers', tempSOCLayers, nSOCLayersRead, maxSOCLayers)

      if (nSOCLayersRead > 0) then
         ! Allocate and store the actual layers
         allocate(SOCLayersArray(nSOCLayersRead))
         SOCLayersArray = tempSOCLayers(1:nSOCLayersRead)

         ! Print detailed output
         call MIO_Print('SOC layer control enabled. Successfully parsed layers for SOC corrections:', 'ham')
         SOCLayersString = ''
         do i = 1, nSOCLayersRead
            if (i == 1) then
               SOCLayersString = trim(num2str(SOCLayersArray(i)))
            else
               SOCLayersString = trim(SOCLayersString) // ', ' // trim(num2str(SOCLayersArray(i)))
            end if
         end do
         call MIO_Print('  Layers: '//trim(SOCLayersString), 'ham')
         call MIO_Print('SOC corrections will ONLY be applied to the above '//trim(num2str(nSOCLayersRead))//' layer(s).', 'ham')
      else
         ! No layers found - this could mean user forgot to specify or parsing failed
         call MIO_Print('WARNING: SOCLayerControl is enabled but no valid layers were found in SOCLayers parameter.', 'ham')
         call MIO_Print('Expected format: space-separated integers (e.g., "1 5" or "1, 5")', 'ham')
         allocate(SOCLayersArray(0))  ! Empty array
      end if
   else
      ! SOCLayerControl disabled - allocate empty array (means apply to all layers)
      allocate(SOCLayersArray(0))
      call MIO_Print('SOC layer control disabled. SOC will be applied to all layers.','ham')
   end if

   ! Layer-specific Haldane control
   call MIO_InputParameter('HaldaneLayerControl',HaldaneLayerControl,.false.)
   if (HaldaneLayerControl) then
      ! Use new variable-length input parameter function
      call MIO_InputParameterVariable('HaldaneLayers', tempHaldaneLayers, nHaldaneLayersRead, maxHaldaneLayers)

      if (nHaldaneLayersRead > 0) then
         ! Allocate and store the actual layers
         allocate(HaldaneLayersArray(nHaldaneLayersRead))
         HaldaneLayersArray = tempHaldaneLayers(1:nHaldaneLayersRead)

         ! Print detailed output
         call MIO_Print('Haldane layer control enabled. Successfully parsed layers for Haldane corrections:', 'ham')
         HaldaneLayersString = ''
         do i = 1, nHaldaneLayersRead
            if (i == 1) then
               HaldaneLayersString = trim(num2str(HaldaneLayersArray(i)))
            else
               HaldaneLayersString = trim(HaldaneLayersString) // ', ' // trim(num2str(HaldaneLayersArray(i)))
            end if
         end do
         call MIO_Print('  Layers: '//trim(HaldaneLayersString), 'ham')
         call MIO_Print('Haldane corrections will ONLY be applied to the above '//trim(num2str(nHaldaneLayersRead))//' layer(s).', 'ham')
      else
         ! No layers found - this could mean user forgot to specify or parsing failed
         call MIO_Print('WARNING: HaldaneLayerControl is enabled but no valid layers were found in HaldaneLayers parameter.', 'ham')
         call MIO_Print('Expected format: space-separated integers (e.g., "1 5" or "1, 5")', 'ham')
         allocate(HaldaneLayersArray(0))  ! Empty array
      end if
   else
      ! HaldaneLayerControl disabled - allocate empty array (means apply to all layers)
      allocate(HaldaneLayersArray(0))
      call MIO_Print('Haldane layer control disabled. Haldane will be applied to all layers.','ham')
   end if

   ! Start with finding all the neighbors. This routine will use one of the
   ! possible routines. The standard, best tested one, is Neigh.fastNNnotsquare.
   ! It is based on a binning approach.
   call MIO_InputParameter('Neigh.fastNN',l,.false.)
   call MIO_InputParameter('Neigh.fastNNnotsquare',ll,.true.)
   call MIO_Print('We use Neigh.fastNNnotsquare by default. This one works as long as the system size is large enough. Otherwise, one has to use one of the specialized routines, for instance when working on commensurate cells that contain less than 10 atoms.','ham')
   call MIO_InputParameter('Neigh.fastNNnotsquareNotRectangle',lll,.false.)
   call MIO_Print('aG ='//trim(num2str(aG)),'ham')
   if (l) then
      aCC = aG/sqrt(3.0_dp)
      cutoff = 1.2_dp
      A1 = ucell(:,1)
      A2 = ucell(:,2)
      maxnn = 20 !to be adjsuted for accuracy
      call fastNN(nAt, Rat(1,:),Rat(2,:),Rat(3,:), aCC,cutoff, A1,A2, maxnn) ! change cutoff compared to other one
      call MIO_Print('fastNN finished','ham')
   else if (ll) then
      aCC = aG/sqrt(3.0_dp)
      ! TB.NeighLevels is the number of intralayer neighbour shells. It sets
      ! the search radius through the shell table of the neigh module; the
      ! neighbour arrays are sized by the search itself.
      nShells = tbnn
      call MIO_InputParameter('Neigh.CutAtNN3',cutAtNN3,.false.)
      if (cutAtNN3) then
         ! Legacy flag: it limited the search to three shells when
         ! TB.NeighLevels was larger than 2.
         call MIO_Print('Neigh.CutAtNN3 is deprecated; use TB.NeighLevels 3 instead','ham')
         nShells = min(nShells,3)
      end if
      cutoff2 = NeighShellCutoff2(nShells,aCC)
      call MIO_Print('Intralayer neighbours: '//trim(num2str(nShells))//' shell(s), in-plane radius '// &
        trim(num2str(sqrt(cutoff2),4))//' Ang','ham')
      call MIO_InputParameter('F2G2Model',F2G2Model,.true.)
      modelShells = 8
      if (F2G2Model) modelShells = 5
      if (nShells < modelShells) then
         call MIO_Print('Intralayer hopping terms beyond shell '//trim(num2str(nShells))// &
           ' are not included (the intralayer models reach shell '//trim(num2str(modelShells))// &
           '); raise TB.NeighLevels to include them','ham')
      end if
      call MIO_InputParameter('Neigh.LayerNeighbors',outplaneNeigh,0)
      call MIO_InputParameter('Bulk',bulk,.false.)
      call MIO_InputParameter('BulkSmall',bulksmall,.false.)
      call MIO_InputParameter('nonBulkSmall',small,.false.)
      call MIO_InputParameter('InterlayerDistance',d,3.22_dp)
      if (outplaneNeigh/=0) then
         call MIO_Print('There are out of plane neighbors','ham')
         call MIO_InputParameter('Neigh.LayerDistFactor',distFact,1.0_dp) ! to increase the cutoff for outerlayer neighbors
         !cutoff2bis = (sqrt((aBN*1.1_dp)**2/3.0_dp + d**2) * distFact)**2.0_dp
         cutoff2bis = (aG/sqrt(3.0_dp)*1.1_dp*distFact)**2.0_dp
      else
         call MIO_Print('There are NO out of plane neighbors','ham')
         cutoff2bis = 1.0_dp
      end if
      !cutoff2 = aCC**2 * 1.2_dp**2
      A1 = ucell(:,1)
      A2 = ucell(:,2)
      A3 = ucell(3,3)
      !ALLOCATE(nn(natoms,maxnn))
      !ALLOCATE(near(natoms))
      if (frac) call AtomsSetCart()
      !print*, cutoff2, cutoff2bis
      if (bulk) then
         if (bulksmall) then
             call fastNNnotsquareBulkSmall(nAt,Rat(1,:),Rat(2,:),Rat(3,:),cutoff2,cutoff2bis,A1,A2,A3)
         else
             call fastNNnotsquareBulk(nAt,Rat(1,:),Rat(2,:),Rat(3,:),cutoff2,cutoff2bis,A1,A2,A3)
         end if
      else
         if (small) then
            call fastNNnotsquareBulkSmall(nAt,Rat(1,:),Rat(2,:),Rat(3,:),cutoff2,cutoff2bis,A1,A2,A3)
            !call fastNNnotsquareSmall(nAt,Rat(1,:),Rat(2,:),Rat(3,:),aCC,cutoff2,cutoff2bis,A1,A2,maxnn) # this one doesn't work
            !for now, let's use bulk, but add a large amount of free space to avoid interactions between periodic images in z
            !direction
            !call fastNNnotsquareBulkSmall(nAt,Rat(1,:),Rat(2,:),Rat(3,:),cutoff2,cutoff2bis,A1,A2,A3)
         else
            call fastNNnotsquare(nAt,Rat(1,:),Rat(2,:),Rat(3,:),cutoff2,cutoff2bis,A1,A2)
         end if
      end if
      call MIO_Print('fastNNnotsquare finished','ham')
   else if (lll) then
      !natoms = nAt
      !x = Rat(1,:)
      !y = Rat(2,:)
      !z = Rat(3,:)
      aCC = aG/sqrt(3.0_dp)
      if (tbnn==1) then
          cutoff2 = aCC**2 * 1.2_dp**2
      else if (tbnn==2) then
          cutoff2 = aG**2 * 1.2_dp**2
      else if (tbnn==3) then
          cutoff2 = (aG**2 + aCC**2) * 1.2_dp**2
      end if
      call MIO_InputParameter('Neigh.LayerNeighbors',outplaneNeigh,0)
      call MIO_InputParameter('InterlayerDistance',d,3.22_dp)
      if (outplaneNeigh/=0) then
         call MIO_Print('There are out of plane neighbors','ham')
         call MIO_InputParameter('Neigh.LayerDistFactor',distFact,1.0_dp) ! to increase the cutoff for outerlayer neighbors
         !cutoff2bis = (sqrt((aBN*1.1_dp)**2/3.0_dp + d**2) * distFact)**2.0_dp
         cutoff2bis = (aG/sqrt(3.0_dp)*1.1_dp*distFact)**2.0_dp
      else
         call MIO_Print('There are NO out of plane neighbors','ham')
         cutoff2bis = 1.0_dp
      end if
      !cutoff2 = aCC**2 * 1.2_dp**2
      A1 = ucell(:,1)
      A2 = ucell(:,2)
      !print*, "A1 and A2", A1, A2
      inplaneNeigh=sum(numN(1:tbnn))
      maxnn = inplaneNeigh + outplaneNeigh !to be adjusted for accuracy
      !ALLOCATE(nn(natoms,maxnn))
      !ALLOCATE(near(natoms))
      if (.not. frac) call AtomsSetFrac()
      call fastNNnotsquareNotRectangle(nAt,Rat(1,:),Rat(2,:),Rat(3,:),aCC,cutoff2,cutoff2bis,A1,A2,maxnn)
      if (frac) call AtomsSetCart()
      !print*, "fastNNnotsquareNotRectangle finished"
      call MIO_Print('fastNNnotsquareNotRectangle finished','ham')
   else
      call NeighList()
   end if
   call MIO_InputParameter('ReadDataFiles',readDataFiles,.false.)
   if (readDataFiles) then
      call MIO_Print('Reading in the (non-remormalized by g0) onsite energies from generate.e','ham')
      call MIO_Allocate(H0,[inode1],[inode2],'H0','ham')
      open(222,FILE="generate.e",STATUS='old')
      do i=1,nAt
         read(222,*) H0(i)
         H0(i) = H0(i)/g0
      end do
      close(222)
   else
       call HamOnSite()
   end if
   !print*, "inode1, maxNeigh, inode2", inode1, maxNeigh, inode2
   call MIO_Allocate(hopp,[1,inode1],[maxNeigh,inode2],'hopp','ham')
   hopp = 0.0_dp
   call MIO_InputParameter('TypeOfSystem',str,'Graphene')
   !if (MIO_StringComp(str,'Ribbons') .or. MIO_StringComp(str,'Hybrid') .or. MIO_StringComp(str,'ReadXYZ')) then
   !   call InterfacePot2(H0,Nneigh,NList,neighCell,neighD,Species,Nradii)
   !end if
   call MIO_InputParameter('MoireStrain',moireStrain,.false.)
   if (moireStrain) then
      call MIO_InputParameter('MoireStrainFactor',strFactor,3.37_dp)
      call MIO_InputParameter('MoireHeightAmp',hStr,0.1_dp)
      call MIO_Print('Hopping modified by moire strain','ham')
      call MIO_Print('Factor: '//trim(num2str(strFactor,2)),'ham')
      call MIO_Print('')
      dCC = aG/sqrt(3.0_dp)
   end if
   call MIO_InputParameter('RandomStrain',randomStrain,.false.)
   if (randomStrain) then
     if (frac) call AtomsSetCart()
     call GaussHeight()
     call MIO_InputParameter('WriteDataFiles',prnt,.false.)
     if (prnt) then
        open(1,FILE='z')
        do i=1,nAt
           write(1,*) Rat(3,i)
        end do
        close(1)
     end if
   end if
   call MIO_InputParameter('printBubble',printBubble,.false.)
   if (printBubble) then
     call MIO_InputParameter('bubbleSigmaR',sigma,1.0_dp)
     call MIO_InputParameter('bubbleHeight',w,1.0_dp)
     if (frac) call AtomsSetCart()
     do i=1,nAt
         dx = Rat(1,3328) - Rat(1,i)
         dy = Rat(2,3328) - Rat(2,i)
         dist2 = dx**2 + dy**2
         Rat(3,i) = Rat(3,i) + w*exp(-dist2**2/(2.0_dp*sigma**2))
     end do
     open(1,FILE='z')
     do i=1,nAt
        write(1,*) Rat(3,i)
     end do
     close(1)
   end if


#ifdef TIMER
   call MIO_TimerStop('ham')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('HamInit',1)
#endif /* DEBUG */

end subroutine HamInit

!===============================================================================
! SUBROUTINE: HamOnSite
!
! DESCRIPTION:
!   Controls the on-site energy terms of the Hamiltonian. The values are stored
!   in H0 and one must always make sure to update rather than reinitialize
!   this value to avoid losing the additive approach when different on-site
!   energy terms exist in the same Hamiltonian.
!
!> @brief Compute and add on-site energy terms to H0.
!! @details Controls the on-site contributions to the Hamiltonian. Values are added to
!! H0 (never reinitialized) to preserve composition of effects (e.g., strain + gates).
!! Handles multiple layer configurations and optional random/moire strain.
!! @note Always update H0 additively; do not reinitialize.
!! @par Side effects
!! - Updates H0(:) in place
!! - May switch between fractional and cartesian coordinates temporarily
!! - Applies optional Gaussian/interface potentials
!! @see HamInit
subroutine HamOnSite()
   use atoms,                only : indxNode, Rat, Species, indxDiv, AtomsSetFrac, AtomsSetCart, AtomsRotate
   use atoms,                only : inode1, inode2, in1, in2, nAt, nAtC1, frac, layerIndex
   use atoms,                only : RatInit
   use atoms,                only : phiForEffectiveModel, interlayerDistances, displacements, displacements_b, displacements_t
   use tbpar,                only : e0_C, g0, e0_B, e0_N, e0_C1, e0_C2, e0_C1_LB, e0_C2_LB, e0_C1_LT, e0_C2_LT
   use cell,                 only : sCell, aG, aBN, ucell
   use constants,            only : twopi, pi
   use math
   use parallel,             only : nDiv, procID
   use gauss,                only : GaussPot, GaussHeight
   use moireBLShift,          only : tauX1, tauY1, tauX2, tauY2
   use interface,            only : InterfacePot1
   !use random,               only : RandNum, RandSeed, rand_t
   use neigh,                only : NeighD

   real(dp), parameter :: phi1= 1.884_dp, phi2 = 1.531_dp
   logical :: threeLayers, twoLayers, fourLayers, middleTwist, fourLayersSandwiched, fiveLayersSandwiched, sixLayersSandwiched, sevenLayersSandwiched, eightLayersSandwiched, tenLayersSandwiched, twentyLayersSandwiched, helicalTwistedMBM, helicalTwistedMBM_CDW
   logical :: l, w, v, u, z, twisted, twisted2, zz, deactivateASubLattice, deactivateBSubLattice
   logical :: ll
   real(dp) :: A, B, C, rand, A1, A2, B1, B2, lambda1, lambda2, CDWAmplitude, CDWPeriod
   real(dp) :: AA, randomPot, randomPot2, randomPot3, tot, tot2, tot3, nImp
   real(dp) :: rand2
   real(dp) :: C0, Cz, Cab, Phi0, Phiz, Phiab, dx, dy, dxPrime, dyPrime, H0jj, Hzjj, Hkjj, Habjj, eps, eps2, ggg, Hjj, Hjj_b, Hjj_t
   real(dp) :: dxTemp, dyTemp
   real(dp) :: dxTemp_b, dyTemp_b
   real(dp) :: dxTemp_t, dyTemp_t
   real(dp) :: dx_b, dy_b
   real(dp) :: dx_t, dy_t
   real(dp) :: CAA, CBB, CApAp, CBpBp, PhiAA, PhiBB, PhiApAp, PhiBpBp
   real(dp) :: CAA0, CBB0, CApAp0, CBpBp0
   real(dp) :: CAA_global, CBB_global, CApAp_global, CBpBp_global, PhiAA_global, PhiBB_global, PhiApAp_global, PhiBpBp_global
   real(dp) :: CAA0_global, CBB0_global, CApAp0_global, CBpBp0_global
   integer :: i, j,k,n1, n2, m1, m2, n
   character(len=80) :: str
   !type(rand_t) :: rng
   integer, pointer :: seed(:)
   integer :: clock
   real(dp) :: onsiteShift, KekuleAngle, KekuleEpsFactor, KekuleAngleGrad,H0_temp
   real(dp) :: twistAngle, twistAngle2, twistAngleGrad, twistAngleGrad2,shift2, shift1, sign2, sign1, minZ, phi2M
   real(dp) :: shift2_x, shift2_y
   logical :: rotateFirst
   real(dp) :: MoireBilayerTopAngle, MoireBilayerTopAngleGrad, ElectricShift, massterm, onsiteEnergyShift
   real(dp) :: MoireBilayerBottomAngle, MoireBilayerBottomAngleGrad
   !real(dp) :: tauX1, tauX2, tauY1, tauY2
   real(dp) :: xValueBeforeRotate, xValueAfterRotate
   real(dp) :: yValueBeforeRotate, yValueAfterRotate
   integer :: checkerDivider, activatedCheckers
   real(dp) :: checkerDensity
   logical, allocatable :: checkerActivate(:,:)
   integer :: sinusNumberOfPeriod
!   real(dp) :: AmplitudeOfSquare, A
   real(dp) :: H, P1, P2, Amp2
   integer :: NOWH
   real(dp) :: Amp3
   real(dp) :: PNPAmp, PNPLeft, PNPRight, limit0, limit1, limit2, limit3, limit4, delta
   real(dp) :: pct

   logical :: randomStrain, prnt
   logical :: deactivateUpperLayer, deactivateUpperLayers
   logical :: GBNtwoLayers, encapsulatedFourLayers, t3GwithBN, tBGDiag, tBGDiagPRB, dontUseInplaneMoire, BNt2GBN, BNBNDiag
   logical :: encapsulatedThreeLayers, removeTopMoireInL2
   logical :: encapsulatedFiveLayers
   logical :: encapsulatedSixLayers
   logical :: encapsulatedSevenLayers
   logical :: useLayerSpecificOnsiteEnergyTerms
   logical :: t2GBN
   logical :: t2BG
   logical :: t3BG
   logical :: BNBNtwoLayers
   logical :: GBNtwoLayersF2G2s, GBNuseDisplacementFile, GBNuseHarmonicApprox
   logical :: tBGuseDisplacementFile, BNBNuseDisplacementFile
   real(dp) :: GlobalPhi,Globalang, rotationAngle, epstBG
   real(dp) :: GlobalPhi2,Globalang2, GlobalPhiL1, GlobalPhiL2, GlobalPhiL3, GlobalPhiL2a, GlobalPhiL2b

   logical :: BNBNtwoLayersF2G2s

   integer :: nnn(4)

   logical :: tBGSwitchDxDy


   logical :: FanZhang, sinusModulationUsingPeriod, sinusModulationAddMassTerm
   logical :: cosinusModulationUsingPeriod, cosinusModulationAddMassTerm, sinusModulationUsingPeriodY, cosinusModulationUsingPeriodY
   real(dp) :: tr1, tr2, tr3, tr4, trDelta, u1, u2, u3, sinusModulationPeriod, cosinusModulationPeriod

   real(dp) :: epsFactor
   real(dp) :: invertEposX, sinusFactor, sinusFactor2, cosinusFactor, cosinusFactor2
   integer :: Efactor,iii
   logical :: invertE, switchHzjj, writeData

   logical :: singleLayerXYZ, distanceDependentEffectiveModel, onlyBottomLayerMassTerm
   real (dp) :: C0d, Czd, BfactorC0, BfactorCz, z0, GBNAngle, moirePreFactor, tBGAngle

   logical :: forceBilayerF2G2Intralayer, sublatticeBasis, addDisplacements, readRigidXYZ
   real(dp) :: layerShift1, layerShift2, layerShift3, layerShift4


#ifdef DEBUG
   call MIO_Debug('HamOnSite',0)
#endif /* DEBUG */

   !print*, displacements(1,:)
   !print*, "-----"
   !print*, displacements(2,:)
   !call MIO_InputParameter('readRigidXYZ',readRigidXYZ,.false.)
   !if (readRigidXYZ) then
   !   print*, Rat(:,1), RatInit(:,1)
   !   print*, "we replace Rat by RatInit... be careful, only tested for the specific case of the effective distant-dependent continuum model"
   !   Rat = RatInit
   !   print*, Rat(:,1), RatInit(:,1)
   !end if

   !print*, interlayerDistances
   ! Removed restriction: Zterm and nspin==2 can now work together
   ! if (Zterm .and. nspin==2) call MIO_Kill('Spin calculation inconsistency','ham','HamOnSite')
   !call MIO_InputParameter('MoireEncapsulatedBilayer',w,.false.)
   !if (w) then
   !    call MIO_Allocate(H0,[inode1],[inode2*2],'H0','ham')
   !else
   call MIO_Allocate(H0,[inode1],[inode2],'H0','ham')
   call MIO_Allocate(H0Bottom,[inode1],[inode2],'H0Bottom','ham')
   call MIO_Allocate(H0Top,[inode1],[inode2],'H0Top','ham')
   ! Initialize H0 arrays to zero to ensure deterministic behavior
   H0 = 0.0_dp
   H0Bottom = 0.0_dp
   H0Top = 0.0_dp
   call MIO_InputParameter('GBNtwoLayersF2G2s',GBNtwoLayersF2G2s,.false.)
   call MIO_InputParameter('GBNtwoLayers',GBNtwoLayers,.false.)
   call MIO_InputParameter('BNBNtwoLayers',BNBNtwoLayers,.false.)
   call MIO_InputParameter('t2GBN',t2GBN,.false.)
   call MIO_InputParameter('t2BG',t2BG,.false.)
   call MIO_InputParameter('t3BG',t3BG,.false.)
   call MIO_InputParameter('encapsulatedThreeLayers',encapsulatedThreeLayers,.false.)
   call MIO_InputParameter('encapsulatedFourLayers',encapsulatedFourLayers,.false.)
   call MIO_InputParameter('encapsulatedFiveLayers',encapsulatedFiveLayers,.false.)
   call MIO_InputParameter('encapsulatedSixLayers',encapsulatedSixLayers,.false.)
   call MIO_InputParameter('encapsulatedSevenLayers',encapsulatedSevenLayers,.false.)
   call MIO_InputParameter('useLayerSpecificOnsiteEnergyTerms',useLayerSpecificOnsiteEnergyTerms,.false.)
   call MIO_InputParameter('t3GwithBN',t3GwithBN,.false.)
   call MIO_InputParameter('BNt2GBN',BNt2GBN,.false.)
  call MIO_InputParameter('dontUseInplaneMoire',dontUseInplaneMoire,.false.)
  call MIO_InputParameter('twoLayers',twoLayers,.false.)
  call MIO_InputParameter('threeLayers',threeLayers,.false.)
  call MIO_InputParameter('fourLayersSandwiched',fourLayersSandwiched,.false.)
  call MIO_InputParameter('helicalTwistedMBM',helicalTwistedMBM,.false.)
  call MIO_InputParameter('fiveLayersSandwiched',fiveLayersSandwiched,.false.)
   call MIO_InputParameter('sixLayersSandwiched',sixLayersSandwiched,.false.)
   call MIO_InputParameter('sevenLayersSandwiched',sevenLayersSandwiched,.false.)
   call MIO_InputParameter('eightLayersSandwiched',eightLayersSandwiched,.false.)
   call MIO_InputParameter('tenLayersSandwiched',tenLayersSandwiched,.false.)
   call MIO_InputParameter('twentyLayersSandwiched',twentyLayersSandwiched,.false.)
   call MIO_InputParameter('middleTwist',middleTwist,.false.)
   call MIO_InputParameter('fourLayers',fourLayers,.false.)
   call MIO_InputParameter('forceBilayerF2G2Intralayer',forceBilayerF2G2Intralayer,.false.)
   call MIO_InputParameter('readRigidXYZ',readRigidXYZ,.false.)
   !end if
   print*, "onsite, B, N, C1, C2", e0_B, e0_N, e0_C1, e0_C2
   !$OMP PARALLEL DO
   do i=1,nAt
      if (Species(i)==3) then
         H0(i) = e0_B
      else if (Species(i)==4) then
         H0(i) = e0_N
      else if (e0_C1 .ne. 0.0_dp .or. e0_C2 .ne. 0.0) then ! we dont enter this clause for most single layer F2G2 systems not sure why i am specifying it
         if (fourLayers) then
             if (Species(i)==1 .and. (layerIndex(i).eq.1 .or. layerIndex(i).eq.3)) then
                H0(i) = e0_C1
             else if (Species(i)==2 .and. (layerIndex(i).eq.1 .or. layerIndex(i).eq.3)) then
                H0(i) = e0_C2
             else if (Species(i)==1 .and. (layerIndex(i).eq.2 .or. layerIndex(i).eq.4)) then
                H0(i) = e0_C2
             else if (Species(i)==2 .and. (layerIndex(i).eq.2 .or. layerIndex(i).eq.4)) then
                H0(i) = e0_C1
             end if
         else if (threeLayers .or. fourLayersSandwiched .or. helicalTwistedMBM .or. fiveLayersSandwiched .or. sixLayersSandwiched .or. sevenLayersSandwiched .or. eightLayersSandwiched .or. tenLayersSandwiched .or. twentyLayersSandwiched) then
             if ((layerIndex(i).eq.1 .or. layerIndex(i).eq.2) .and. (.not. middleTwist)) then
                if (Species(i)==1 .and. (layerIndex(i).eq.1)) then
                   H0(i) = e0_C1
                else if (Species(i)==2 .and. (layerIndex(i).eq.1)) then
                   H0(i) = e0_C2
                else if (Species(i)==1 .and. (layerIndex(i).eq.2)) then
                   H0(i) = e0_C2
                else if (Species(i)==2 .and. (layerIndex(i).eq.2)) then
                   H0(i) = e0_C1
                end if
             else
                H0(i) = (e0_C1+e0_C2)/2.0_dp ! Is this what we want for the N>2-layer systems?
             end if
         else if (GBNtwoLayers) then
            !if (GBNtwoLayersF2G2s) then
              ! add F2G2 model onsite parameters for G, B and N atoms
              if (layerIndex(i) .eq. 1) then ! graphene layer
                if (Species(i).eq.1) then
                   H0(i) = e0_C1
                else if (Species(i).eq.2) then
                   H0(i) = e0_C2
                end if
              else if (layerIndex(i) .eq. 2) then ! BN layer
                cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
              end if
            !else ! add dx, dy-dependent values using first order harmonic approximation (see summary by Jiaqi)
            !   cycle
            !end if
         else if (t3BG) then
              if (Species(i)==1 .and. (layerIndex(i).eq.1 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.5)) then
!                      print*, "taken onsite t3BG =",e0_C1
                H0(i) = e0_C1
             else if (Species(i)==2 .and. (layerIndex(i).eq.1 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.5)) then
 !                     print*, "taken onsite t3BG =",e0_C2
                H0(i) = e0_C2
             else if (Species(i)==1 .and. (layerIndex(i).eq.2 .or. layerIndex(i).eq.4 .or. layerIndex(i).eq.6)) then
                H0(i) = e0_C2
             else if (Species(i)==2 .and. (layerIndex(i).eq.2 .or. layerIndex(i).eq.4 .or. layerIndex(i).eq.6)) then
                H0(i) = e0_C1
            end if
         else if (t2GBN) then
              if (layerIndex(i) .eq. 1) then ! bottom substrate hBN
                cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
              else if (layerIndex(i) .eq. 2) then ! middle graphene layer
                if (Species(i).eq.1) then
                   H0(i) = e0_C1
                else if (Species(i).eq.2) then
                   H0(i) = e0_C2
                end if
              else if (layerIndex(i) .eq. 3) then ! top graphene layer
                if (Species(i).eq.1) then
                   H0(i) = e0_C1
                else if (Species(i).eq.2) then
                   H0(i) = e0_C2
                end if
              end if
         else if (encapsulatedThreeLayers) then
            !if (GBNtwoLayersF2G2s) then
              ! add F2G2 model onsite parameters for G, B and N atoms
              if (layerIndex(i) .eq. 2) then ! graphene layer
                if (Species(i).eq.1) then
                   H0(i) = e0_C1
                else if (Species(i).eq.2) then
                   H0(i) = e0_C2
                end if
              else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 3) then ! BN layer
                cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
              end if
            !else ! add dx, dy-dependent values using first order harmonic approximation (see summary by Jiaqi)
            !   cycle
            !end if
         else if (encapsulatedFourLayers) then
             if (useLayerSpecificOnsiteEnergyTerms) then
                 if (layerIndex(i) .eq. 2) then
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1_LB
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2_LB
                   end if
                 else if (layerIndex(i) .eq. 3) then
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1_LT
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2_LT
                   end if
                 else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 4) then ! BN layer
                   cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
                 end if
             else
                 if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3) then ! graphene layer
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2
                   end if
                 else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 4) then ! BN layer
                   cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
                 end if
             end if
         else if (encapsulatedFiveLayers) then
              if (useLayerSpecificOnsiteEnergyTerms) then
                 if (layerIndex(i) .eq. 2) then
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1_LB
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2_LB
                   end if
                 else if (layerIndex(i) .eq. 4) then
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1_LT
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2_LT
                   end if
                 else if (layerIndex(i) .eq. 3) then ! graphene layer
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2
                   end if
                 else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 5) then ! BN layer
                   cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
                 end if
              else
                 if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3 .or. layerIndex(i) .eq. 4) then ! graphene layer
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2
                   end if
                 else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 5) then ! BN layer
                   cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
                 end if
              end if
         else if (encapsulatedSixLayers) then
              if (useLayerSpecificOnsiteEnergyTerms) then
                 if (layerIndex(i) .eq. 2) then
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1_LB
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2_LB
                   end if
                 else if (layerIndex(i) .eq. 5) then
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1_LT
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2_LT
                   end if
                 else if (layerIndex(i) .eq. 3 .or. layerIndex(i) .eq. 4) then ! graphene layer
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2
                   end if
                 else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 6) then ! BN layer
                   cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
                 end if
              else
                 if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3 .or. layerIndex(i) .eq. 4 .or. layerIndex(i) .eq. 5) then ! graphene layer
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2
                   end if
                 else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 6) then ! BN layer
                   cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
                 end if
              end if
         else if (encapsulatedSevenLayers) then
              if (useLayerSpecificOnsiteEnergyTerms) then
                 if (layerIndex(i) .eq. 2) then
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1_LB
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2_LB
                   end if
                 else if (layerIndex(i) .eq. 6) then
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1_LT
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2_LT
                   end if
                 else if (layerIndex(i) .eq. 3 .or. layerIndex(i) .eq. 4 .or. layerIndex(i) .eq. 5 ) then ! graphene layer
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2
                   end if
                 else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 7) then ! BN layer
                   cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
                 end if
              else
                 if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3 .or. layerIndex(i) .eq. 4 .or. layerIndex(i) .eq. 5 .or. layerIndex(i).eq. 6) then ! graphene layer
                   if (Species(i).eq.1) then
                      H0(i) = e0_C1
                   else if (Species(i).eq.2) then
                      H0(i) = e0_C2
                   end if
                 else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 7) then ! BN layer
                   cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
                 end if
              end if
         else if (t3GWithBN) then
            !if (GBNtwoLayersF2G2s) then
              ! add F2G2 model onsite parameters for G, B and N atoms
              if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3 .or. layerIndex(i) .eq. 4) then ! graphene layer
                if (Species(i).eq.1) then
                   H0(i) = e0_C1
                else if (Species(i).eq.2) then
                   H0(i) = e0_C2
                end if
              else if (layerIndex(i) .eq. 1) then ! BN layer
                cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
              end if
            !else ! add dx, dy-dependent values using first order harmonic approximation (see summary by Jiaqi)
            !   cycle
            !end if
         else if (BNt2GBN) then
            !if (GBNtwoLayersF2G2s) then
              ! add F2G2 model onsite parameters for G, B and N atoms
              if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3) then ! graphene layer
                if (Species(i).eq.1) then
                   H0(i) = e0_C1
                else if (Species(i).eq.2) then
                   H0(i) = e0_C2
                end if
              else if (layerIndex(i) .eq. 1 .or. layerIndex(i).eq. 4) then ! BN layer
                cycle ! already assigned 40 lines higher (search for e0_B and e0_N)
              end if
            !else ! add dx, dy-dependent values using first order harmonic approximation (see summary by Jiaqi)
            !   cycle
            !end if
         else if (BNBNtwoLayers) then ! checked
              ! add F2G2 model parameters for BN and BN
            cycle ! already assigned 40 lines higher
            !if (BNBNtwoLayersF2G2s) then
            !  if (Species(i) .eq. 3) then ! boron atoms
            !      H0(i) = e0_B
            !  else if (Species(i) .eq. 4) then ! nitrogen atoms
            !      H0(i) = e0_N
            !  end if
            !end if
         else if (forceBilayerF2G2Intralayer) then
            if (Species(i)==1 .and. (layerIndex(i).eq.1)) then
               H0(i) = e0_C1
            else if (Species(i)==2 .and. (layerIndex(i).eq.1)) then
               H0(i) = e0_C2
            else if (Species(i)==1 .and. (layerIndex(i).eq.2)) then
               H0(i) = e0_C2
            else if (Species(i)==2 .and. (layerIndex(i).eq.2)) then
               H0(i) = e0_C1
            end if
         else if (twoLayers) then
            if (Species(i).eq.1) then
               H0(i) = e0_C1
            else if (Species(i).eq.2) then
               H0(i) = e0_C2
            end if
         end if
      else
         H0(i) = e0_C
      end if
   end do
   !$OMP END PARALLEL DO
   if (nspin==2) then
      call MIO_Allocate(Ho,[inode1],[inode2],'Ho','ham')
      Ho = H0
   end if
   call MIO_InputParameter('MoirePotential',l,.false.)
   call MIO_InputParameter('MoirePreFactor',moirePreFactor,1.0_dp)
   call MIO_InputParameter('tBGDiag',tBGDiag,.false.)
   call MIO_InputParameter('tBGDiagPRB',tBGDiagPRB,.false.)
   call MIO_InputParameter('BNBNDiag',BNBNDiag,.false.)
   call MIO_InputParameter('Latticepercent',eps,-0.018181818181818_dp)
   !call MIO_Print('We are adding an effective moire potential in the onsite energies','ham')
   if (l) then
      call MIO_InputParameter('MoireSachs',l,.false.)
      if (l) then
         call MIO_Print('Moire potential model [Sachs et al. PRB 84, 195414 (2011)]','ham')
         call MIO_InputParameter('MoirePotA',A,0.0186_dp)
         call MIO_InputParameter('MoirePotB',B,0.042_dp)
         call MIO_InputParameter('MoirePotC',C,0.0_dp)
      end if
! --- New parameters for the moire pattern ---
      call MIO_InputParameter('MoireJeil',l,.false.)
      if (l) then
         call MIO_Print('Moire potential model [Jeil Jung, PRB 89, 205414 (2014)]','ham')
         call MIO_InputParameter('MoirePrefactor',moirePreFactor,1.0_dp)
         call MIO_InputParameter('MoirePotC0',C0,-0.01013_dp)
         C0 = C0/g0*moirePreFactor
         call MIO_InputParameter('MoirePotCz',Cz,-0.00901_dp)
         Cz = Cz/g0*moirePrefactor
         call MIO_InputParameter('MoirePotCab',Cab,0.01134_dp)
         Cab = Cab/g0*moirePrefactor
         call MIO_InputParameter('MoirePotPhi0',Phi0,1.510233401750693_dp)
         call MIO_InputParameter('MoirePotPhiz',Phiz,0.147131255943122_dp)
         call MIO_InputParameter('MoirePotPhiab',Phiab,0.342084533390889_dp)

         call MIO_InputParameter('distanceDependentEffectiveModel',distanceDependentEffectiveModel,.false.)
         call MIO_InputParameter('BfactorC0',BfactorC0,3.1_dp)
         call MIO_InputParameter('BfactorCz',BfactorCz,3.1_dp)
         call MIO_InputParameter('InterlayerDistance',z0,3.35_dp)
         call MIO_InputParameter('sublatticeBasis',sublatticeBasis,.false.)
         call MIO_InputParameter('addDisplacements',addDisplacements,.false.)
      end if
! --- ---

      if (.not. frac) call AtomsSetFrac()
      call MIO_InputParameter('MoireTrilayer',l,.false.)
      call MIO_InputParameter('GBNtwoLayers',GBNtwoLayers,.false.)
      call MIO_InputParameter('BNBNtwoLayers',BNBNtwoLayers,.false.)
      call MIO_InputParameter('t2GBN',t2GBN,.false.)
      call MIO_InputParameter('t2BG',t2BG,.false.)
      call MIO_InputParameter('t3BG',t3BG,.false.)
      call MIO_InputParameter('encapsulatedThreeLayers',encapsulatedThreeLayers,.false.)
      call MIO_InputParameter('encapsulatedFourLayers',encapsulatedFourLayers,.false.)
      call MIO_InputParameter('encapsulatedFiveLayers',encapsulatedFiveLayers,.false.)
      call MIO_InputParameter('encapsulatedSixLayers',encapsulatedSixLayers,.false.)
      call MIO_InputParameter('encapsulatedSevenLayers',encapsulatedSevenLayers,.false.)
      call MIO_InputParameter('removeTopMoireInL2',removeTopMoireInL2,.false.)
      call MIO_InputParameter('t3GWithBN',t3GWithBN,.false.)
      call MIO_InputParameter('BNt2GBN',BNt2GBN,.false.)
      call MIO_InputParameter('GBNtwoLayersF2G2s',GBNtwoLayersF2G2s,.false.)
      call MIO_InputParameter('GBNuseDisplacementFile',GBNuseDisplacementFile,.false.)
      call MIO_InputParameter('tBGuseDisplacementFile',tBGuseDisplacementFile,.false.)
      call MIO_InputParameter('BNBNuseDisplacementFile',BNBNuseDisplacementFile,.false.)
      call MIO_InputParameter('GBNuseHarmonicApprox',GBNuseHarmonicApprox,.false.)
      ! For GBNtwoLayers with GBNuseHarmonicApprox, automatically disable tBGDiag and BNBNDiag
      ! to ensure the harmonic approximation path is taken instead
      if (GBNtwoLayers .and. GBNuseHarmonicApprox) then
         if (tBGDiag) then
            call MIO_Print('WARNING: tBGDiag is true but GBNuseHarmonicApprox is also true for GBNtwoLayers. Disabling tBGDiag.','ham')
            tBGDiag = .false.
         end if
         if (BNBNDiag) then
            call MIO_Print('WARNING: BNBNDiag is true but GBNuseHarmonicApprox is also true for GBNtwoLayers. Disabling BNBNDiag.','ham')
            BNBNDiag = .false.
         end if
      end if
      call MIO_InputParameter('GBNAngle',GBNAngle,0.0_dp)
      call MIO_InputParameter('GlobalTwist',Globalang,0.0_dp)
      call MIO_InputParameter('GlobalTwist2',Globalang2,0.0_dp)
      call MIO_InputParameter('GlobalPhiL1',GlobalPhiL1,0.0_dp)
      call MIO_InputParameter('GlobalPhiL2',GlobalPhiL2,0.0_dp)
      call MIO_InputParameter('GlobalPhiL2a',GlobalPhiL2a,0.0_dp)
      call MIO_InputParameter('GlobalPhiL2b',GlobalPhiL2b,0.0_dp)
      call MIO_InputParameter('GlobalPhiL3',GlobalPhiL3,0.0_dp)
      call MIO_InputParameter('rotationAngle',rotationAngle,0.0_dp)
      call MIO_InputParameter('tBGSwitchDxDy',tBGSwitchDxDy,.false.)
      ! Initialize global variables to avoid using uninitialized values
      CAA_global = 0.0_dp
      CBB_global = 0.0_dp
      CApAp_global = 0.0_dp
      CBpBp_global = 0.0_dp
      CAA0_global = 0.0_dp
      CBB0_global = 0.0_dp
      CApAp0_global = 0.0_dp
      CBpBp0_global = 0.0_dp
      PhiAA_global = 0.0_dp
      PhiBB_global = 0.0_dp
      PhiApAp_global = 0.0_dp
      PhiBpBp_global = 0.0_dp
      if (tBGDiag) then
         if (tBGDiagPRB) then
             call MIO_InputParameter('CAA',CAA_global,0.00110_dp)
             CAA_global = CAA_global/g0
             call MIO_InputParameter('CBB',CBB_global,0.00110_dp)
             CBB_global = CBB_global/g0
             call MIO_InputParameter('CAA0',CAA0_global,0.0_dp)
             CAA0_global = CAA0_global/g0
             call MIO_InputParameter('CBB0',CBB0_global,0.0_dp)
             CBB0_global = CBB0_global/g0
             call MIO_InputParameter('CApAp',CApAp_global,0.00110_dp)
             CApAp_global = CApAp_global/g0
             call MIO_InputParameter('CBpBp',CBpBp_global,0.00110_dp)
             CBpBp_global = CBpBp_global/g0
             call MIO_InputParameter('CApAp0',CApAp0_global,0.0_dp)
             CApAp0_global = CApAp0_global/g0
             call MIO_InputParameter('CBpBp0',CBpBp0_global,0.0_dp)
             CBpBp0_global = CBpBp0_global/g0
             call MIO_InputParameter('PhiAA',PhiAA_global,82.54_dp)
             PhiAA_global = PhiAA_global*pi/180.0_dp
             call MIO_InputParameter('PhiBB',PhiBB_global,-82.54_dp)
             PhiBB_global = PhiBB_global*pi/180.0_dp
             call MIO_InputParameter('PhiApAp',PhiApAp_global,-82.54_dp)
             PhiApAp_global = PhiApAp_global*pi/180.0_dp
             call MIO_InputParameter('PhiBpBp',PhiBpBp_global,82.54_dp)
             PhiBpBp_global = PhiBpBp_global*pi/180.0_dp
         else
             call MIO_Print('we are going to calculate interlayer-distance-dependent tBG intralayer  moire parametrization inside the loop using Srivani parametrization','ham')
         end if
      else if (BNBNDiag) then
         !if (tBGDiagPRB) then
             call MIO_InputParameter('CAA',CAA_global,0.00110_dp)
             CAA_global = CAA_global/g0
             call MIO_InputParameter('CBB',CBB_global,0.00110_dp)
             CBB_global = CBB_global/g0
             call MIO_InputParameter('CAA0',CAA0_global,0.0_dp)
             CAA0_global = CAA0_global/g0
             call MIO_InputParameter('CBB0',CBB0_global,0.0_dp)
             CBB0_global = CBB0_global/g0
             call MIO_InputParameter('CApAp',CApAp_global,0.00110_dp)
             CApAp_global = CApAp_global/g0
             call MIO_InputParameter('CBpBp',CBpBp_global,0.00110_dp)
             CBpBp_global = CBpBp_global/g0
             call MIO_InputParameter('CApAp0',CApAp0_global,0.0_dp)
             CApAp0_global = CApAp0_global/g0
             call MIO_InputParameter('CBpBp0',CBpBp0_global,0.0_dp)
             CBpBp0_global = CBpBp0_global/g0
             call MIO_InputParameter('PhiAA',PhiAA_global,82.54_dp)
             PhiAA_global = PhiAA_global*pi/180.0_dp
             call MIO_InputParameter('PhiBB',PhiBB_global,-82.54_dp)
             PhiBB_global = PhiBB_global*pi/180.0_dp
             call MIO_InputParameter('PhiApAp',PhiApAp_global,-82.54_dp)
             PhiApAp_global = PhiApAp_global*pi/180.0_dp
             call MIO_InputParameter('PhiBpBp',PhiBpBp_global,82.54_dp)
             PhiBpBp_global = PhiBpBp_global*pi/180.0_dp
         !else
         !    call MIO_Print('we are going to calculate interlayer-distance-dependent tBG intralayer  moire parametrization inside the loop using Srivani parametrization','ham')
         !end if
      else
         call MIO_InputParameter('CAA',CAA,0.005733_dp)
         CAA = CAA/g0
         call MIO_InputParameter('CBB',CBB,0.004826_dp)
         CBB = CBB/g0
         call MIO_InputParameter('CAA0',CAA0,3.332_dp)
         CAA0 = CAA0/g0
         call MIO_InputParameter('CBB0',CBB0,-1.493_dp)
         CBB0 = CBB0/g0
         call MIO_InputParameter('CApAp',CApAp,-0.005703_dp)
         CApAp = CApAp/g0
         call MIO_InputParameter('CBpBp',CBpBp,-0.003596_dp)
         CBpBp = CBpBp/g0
         call MIO_InputParameter('CApAp0',CApAp0,0.0_dp)
         CApAp0 = CApAp0/g0
         call MIO_InputParameter('CBpBp0',CBpBp0,0.0_dp)
         CBpBp0 = CBpBp0/g0
         call MIO_InputParameter('PhiAA',PhiAA,90.0_dp)
         PhiAA = PhiAA*pi/180.0_dp
         call MIO_InputParameter('PhiBB',PhiBB,65.49_dp)
         PhiBB = PhiBB*pi/180.0_dp
         call MIO_InputParameter('PhiApAp',PhiApAp,87.51_dp)
         PhiApAp = PhiApAp*pi/180.0_dp
         call MIO_InputParameter('PhiBpBp',PhiBpBp,65.06_dp)
         PhiBpBp = PhiBpBp*pi/180.0_dp
      end if
      !if (tBGDiag) then
      !    rotationAngle = rotationAngle*pi/180.0_dp
      !    epstBG = 0.0_dp
      !    GlobalPhi = atan((1.0_dp+epstBG)*sin(rotationAngle)/((1.0_dp+epstBG)*cos(rotationAngle)-1.0_dp))
      !    GlobalPhi = -(GlobalPhi+120.0_dp*pi/180.0_dp)
      !else
      GBNAngle = GBNAngle*pi/180.0_dp
         GlobalPhi = Globalang*pi/180.0_dp
         GlobalPhi2 = Globalang2*pi/180.0_dp
      GlobalPhiL1 = GlobalPhiL1*pi/180.0_dp
      GlobalPhiL2 = GlobalPhiL2*pi/180.0_dp
      GlobalPhiL2a = GlobalPhiL2a*pi/180.0_dp
      GlobalPhiL2b = GlobalPhiL2b*pi/180.0_dp
      GlobalPhiL3 = GlobalPhiL3*pi/180.0_dp
      !end if
      call MIO_InputParameter('WriteDataFiles',writeData,.false.)
      if (writeData) then
         !print*, "we are going to print the diagonal intralayer moire terms"
         open(588,FILE='Hjj1')
         open(589,FILE='Hjj2')
      end if
      if (dontUseInplaneMoire) then
          call MIO_Print('We are deactivating the inplane moire terms','ham')
      else if (GBNtwoLayers .or. BNBNtwoLayers .or. encapsulatedThreeLayers .or.  encapsulatedFourLayers .or. encapsulatedFiveLayers .or. t3GwithBN .or. BNt2GBN .or. tBGDiag .or. t2GBN) then ! adding in plane moire
         !call MIO_Print('We are working on the GBNtwoLayers system','ham')
         if (frac) call AtomsSetCart()
         !if (GBNtwoLayersF2G2s) then
            !print*, "we are not adding any other onsite energies as we are using the F2G2s from GBN, please make sure the lattice constant for the BN layer is correct: aBN= ", aBN
            !print*, "we are adding additional displacement-dependent terms for the intralayer model here, also make sure the lattice constant for the BN layer is correct: aBN= ", aBN
         !else
            !$OMP PARALLEL DO PRIVATE(i,dx,dy,dx_b,dx_t,dy_b,dy_t,dxTemp,dyTemp,dxTemp_b,dxTemp_t,dyTemp_b,dyTemp_t, Hjj, Hjj_b, Hjj_t, A, B, C, CAA, CBB, CAA0, CBB0, CApAp, CBpBp, CApAp0, CBpBp0, PhiAA, PhiBB, PhiApAp, PhiBpBp)
            do i=1,nAt
              ! Initialize all variables to avoid using uninitialized values
              Hjj = 0.0_dp
              Hjj_b = 0.0_dp
              Hjj_t = 0.0_dp
              dx = 0.0_dp
              dy = 0.0_dp
              dx_b = 0.0_dp
              dx_t = 0.0_dp
              dy_b = 0.0_dp
              dy_t = 0.0_dp
              dxTemp = 0.0_dp
              dyTemp = 0.0_dp
              dxTemp_b = 0.0_dp
              dxTemp_t = 0.0_dp
              dyTemp_b = 0.0_dp
              dyTemp_t = 0.0_dp
              A = 0.0_dp
              B = 0.0_dp
              C = 0.0_dp
              CAA = 0.0_dp
              CBB = 0.0_dp
              CAA0 = 0.0_dp
              CBB0 = 0.0_dp
              CApAp = 0.0_dp
              CBpBp = 0.0_dp
              CApAp0 = 0.0_dp
              CBpBp0 = 0.0_dp
              PhiAA = 0.0_dp
              PhiBB = 0.0_dp
              PhiApAp = 0.0_dp
              PhiBpBp = 0.0_dp
              if (tBGDiag) then
                 ! Debug: Check if we're entering tBGDiag path
                 if (writeData .and. layerIndex(i).eq.1 .and. (Species(i).eq.1 .or. Species(i).eq.2) .and. i <= 5) then
                    write(588,*) '*** WARNING: Entered tBGDiag path instead of GBNuseHarmonicApprox!'
                 end if
                 if (tBGDiagPRB) then
                    CAA = CAA_global
                    CBB = CBB_global
                    CAA0 = CAA0_global
                    CBB0 = CBB0_global
                    CApAp = CApAp_global
                    CBpBp = CBpBp_global
                    CApAp0 = CApAp0_global
                    CBpBp0 = CBpBp0_global
                    PhiAA = PhiAA_global
                    PhiBB = PhiBB_global
                    PhiApAp = PhiApAp_global
                    PhiBpBp = PhiBpBp_global
                 else
                     call fitFunSrivani(0.0082_dp, -0.0588_dp, 0.1060_dp, interlayerDistances(i), CAA)
                     CAA = CAA/g0
                     CBB = CAA
                     CApAp = CAA
                     CBpBp = CAA
                     !call fitFunSrivani(0.0082_dp, -0.0588_dp, 0.1060_dp, interlayerDistances(i), CBB)
                     !call fitFunSrivani(-0.0051_dp, 0.04_dp, -0.0794_dp, interlayerDistances(i), CAA0)
                     !CAA0 = CAA0/g0
                     CAA0 = 0.0_dp
                     !call fitFunSrivani(-0.0051_dp, 0.04_dp, -0.0794_dp, interlayerDistances(i), CBB0)
                     !CBB0 = CBB0/g0
                     CBB0 = 0.0_dp
                     !call fitFunSrivani(0.0082_dp, -0.0588_dp, 0.1060_dp, interlayerDistances(i), CApAp)
                     !call fitFunSrivani(0.0082_dp, -0.0588_dp, 0.1060_dp, interlayerDistances(i), CBpBp)
                     !CBpBp = CBpBp/g0
                     !call fitFunSrivani(-0.0051_dp, 0.04_dp, -0.0794_dp, interlayerDistances(i), CApAp0)
                     !CApAp0 = CApAp0/g0
                     CApAp0 = 0.0_dp
                     !call fitFunSrivani(-0.0051_dp, 0.04_dp, -0.0794_dp, interlayerDistances(i), CBpBp0)
                     !CBpBp0 = CBpBp0/g0
                     CBpBp0 = 0.0_dp
                     call fitFunSrivani(41.150_dp, -288.4_dp, 570.9_dp, interlayerDistances(i), PhiAA)
                     PhiAA = PhiAA*pi/180.0_dp
                     PhiBpBp = PhiAA
                     PhiBB = -PhiAA
                     PhiApAp = -PhiAA
                     !call fitFunSrivani(-25.90_dp, 188.7_dp, -407.9_dp, interlayerDistances(i), PhiBB)
                     !PhiBB = PhiBB*pi/180.0_dp
                     !call fitFunSrivani(-25.90_dp, 188.7_dp, -407.9_dp, interlayerDistances(i), PhiApAp)
                     !PhiApAp = PhiApAp*pi/180.0_dp
                     !call fitFunSrivani(41.150_dp, -288.4_dp, 570.9_dp, interlayerDistances(i), PhiBpBp)
                     !PhiBpBp = PhiBpBp*pi/180.0_dp
                 end if
              else if (BNBNDiag) then
                 ! Debug: Check if we're entering BNBNDiag path
                 if (writeData .and. layerIndex(i).eq.1 .and. (Species(i).eq.1 .or. Species(i).eq.2) .and. i <= 5) then
                    write(588,*) '*** WARNING: Entered BNBNDiag path instead of GBNuseHarmonicApprox!'
                 end if
                 !if (tBGDiagPRB) then
                    CAA = CAA_global
                    CBB = CBB_global
                    CAA0 = CAA0_global
                    CBB0 = CBB0_global
                    CApAp = CApAp_global
                    CBpBp = CBpBp_global
                    CApAp0 = CApAp0_global
                    CBpBp0 = CBpBp0_global
                    PhiAA = PhiAA_global
                    PhiBB = PhiBB_global
                    PhiApAp = PhiApAp_global
                    PhiBpBp = PhiBpBp_global
                 !else
                 !    call fitFunSrivani(0.0082_dp, -0.0588_dp, 0.1060_dp, interlayerDistances(i), CAA)
                 !    CAA = CAA/g0
                 !    CBB = CAA
                 !    CApAp = CAA
                 !    CBpBp = CAA
                 !    !call fitFunSrivani(0.0082_dp, -0.0588_dp, 0.1060_dp, interlayerDistances(i), CBB)
                 !    !call fitFunSrivani(-0.0051_dp, 0.04_dp, -0.0794_dp, interlayerDistances(i), CAA0)
                 !    !CAA0 = CAA0/g0
                 !    CAA0 = 0.0_dp
                 !    !call fitFunSrivani(-0.0051_dp, 0.04_dp, -0.0794_dp, interlayerDistances(i), CBB0)
                 !    !CBB0 = CBB0/g0
                 !    CBB0 = 0.0_dp
                 !    !call fitFunSrivani(0.0082_dp, -0.0588_dp, 0.1060_dp, interlayerDistances(i), CApAp)
                 !    !call fitFunSrivani(0.0082_dp, -0.0588_dp, 0.1060_dp, interlayerDistances(i), CBpBp)
                 !    !CBpBp = CBpBp/g0
                 !    !call fitFunSrivani(-0.0051_dp, 0.04_dp, -0.0794_dp, interlayerDistances(i), CApAp0)
                 !    !CApAp0 = CApAp0/g0
                 !    CApAp0 = 0.0_dp
                 !    !call fitFunSrivani(-0.0051_dp, 0.04_dp, -0.0794_dp, interlayerDistances(i), CBpBp0)
                 !    !CBpBp0 = CBpBp0/g0
                 !    CBpBp0 = 0.0_dp
                 !    call fitFunSrivani(41.150_dp, -288.4_dp, 570.9_dp, interlayerDistances(i), PhiAA)
                 !    PhiAA = PhiAA*pi/180.0_dp
                 !    PhiBpBp = PhiAA
                 !    PhiBB = -PhiAA
                 !    PhiApAp = -PhiAA
                 !    !call fitFunSrivani(-25.90_dp, 188.7_dp, -407.9_dp, interlayerDistances(i), PhiBB)
                 !    !PhiBB = PhiBB*pi/180.0_dp
                 !    !call fitFunSrivani(-25.90_dp, 188.7_dp, -407.9_dp, interlayerDistances(i), PhiApAp)
                 !    !PhiApAp = PhiApAp*pi/180.0_dp
                 !    !call fitFunSrivani(41.150_dp, -288.4_dp, 570.9_dp, interlayerDistances(i), PhiBpBp)
                 !    !PhiBpBp = PhiBpBp*pi/180.0_dp
                 !end if
              end if
              if (tBGuseDisplacementFile) then
                 if (layerIndex(i).eq.1) then
                     dx = displacements(1,i)*cos(GlobalPhi2)+displacements(2,i)*sin(GlobalPhi2)
                     dy = displacements(2,i)*cos(GlobalPhi2)-displacements(1,i)*sin(GlobalPhi2)
                 else if (layerIndex(i).eq.2) then
                     dx = displacements(1,i)*cos(GlobalPhi)+displacements(2,i)*sin(GlobalPhi)
                     dy = displacements(2,i)*cos(GlobalPhi)-displacements(1,i)*sin(GlobalPhi)
                 end if
                 !dx = displacements(1,i)
                 !dy = displacements(2,i)
                 !print*, "dx, dy", dx, dy
              else if (GBNuseDisplacementFile) then
                 !print*, "starting1"
                 !print*, displacements(1,i)
                 !print*, GlobalPhi
                 !dx = displacements(1,i)*cos(GlobalPhi)+displacements(2,i)*sin(GlobalPhi)
                 !dy = displacements(2,i)*cos(GlobalPhi)-displacements(1,i)*sin(GlobalPhi)
                 if (encapsulatedThreeLayers .and. layerIndex(i).eq.2) then
                    dx_b = displacements_b(1,i)*cos(GlobalPhiL2a)+displacements_b(2,i)*sin(GlobalPhiL2a)
                    dy_b = displacements_b(2,i)*cos(GlobalPhiL2a)-displacements_b(1,i)*sin(GlobalPhiL2a)
                    dx_t = displacements_t(1,i)*cos(GlobalPhiL2b)+displacements_t(2,i)*sin(GlobalPhiL2b)
                    dy_t = displacements_t(2,i)*cos(GlobalPhiL2b)-displacements_t(1,i)*sin(GlobalPhiL2b)
                 else if (encapsulatedThreeLayers .and. layerIndex(i).eq.1) then
                    dx = displacements(1,i)*cos(GlobalPhiL1)+displacements(2,i)*sin(GlobalPhiL1)
                    dy = displacements(2,i)*cos(GlobalPhiL1)-displacements(1,i)*sin(GlobalPhiL1)
                 else if (encapsulatedThreeLayers .and. layerIndex(i).eq.3) then
                    dx = displacements(1,i)*cos(GlobalPhiL3)+displacements(2,i)*sin(GlobalPhiL3)
                    dy = displacements(2,i)*cos(GlobalPhiL3)-displacements(1,i)*sin(GlobalPhiL3)
                 else
                    dx = displacements(1,i)*cos(GlobalPhi)+displacements(2,i)*sin(GlobalPhi)
                    dy = displacements(2,i)*cos(GlobalPhi)-displacements(1,i)*sin(GlobalPhi)
                 end if
                 !dx = displacements(1,i)
                 !dy = displacements(2,i)
                 !print*, "dx, dy", dx, dy
              !else if (t2GBNuseDisplacementFile) then
              !   dx = displacements(1,i)*cos(GlobalPhi)+displacements(2,i)*sin(GlobalPhi)
              !   dy = displacements(2,i)*cos(GlobalPhi)-displacements(1,i)*sin(GlobalPhi)
              !   !dx = displacements(1,i)
              !   !dy = displacements(2,i)
              !   !print*, "dx, dy", dx, dy
              else
                 dx = ((1.0_dp+eps) * cos(GBNAngle) - 1.0_dp) * (Rat(1,i)) &
                           - ((1.0_dp+eps) * sin(GBNAngle) * (Rat(2,i)))
                 dy = ((1.0_dp+eps) * cos(GBNAngle) - 1.0_dp) * (Rat(2,i)) &
                           + ((1.0_dp+eps) * sin(GBNAngle) * (Rat(1,i)))
              end if
              if (tBGDiag) then
                 if (tBGuseDisplacementFile) then
                    if (layerIndex(i).eq.1) then
                       dxTemp = -dx/aG
                       dyTemp = -dy/aG
                    else if (layerIndex(i).eq.2) then
                       dxTemp = dx/aG
                       dyTemp = dy/aG
                       !dxTemp = dx/aG
                       !dyTemp = dy/aG
                    end if
                 else
                    if (layerIndex(i).eq.1) then
                       dxTemp = dx
                       dyTemp = dy
                    else if (layerIndex(i).eq.2) then
                       dxTemp = -dx
                       dyTemp = -dy
                    end if
                 end if
                 !end if
                 !if (Species(i).eq.1 .and. layerIndex(i).eq.2) then
                 !   call diagotBG(Hjj,dxTemp,dyTemp,CBB,PhiBB,CBB0)
                 !else if (Species(i).eq.2 .and. layerIndex(i).eq.1) then
                 !   call diagotBG(Hjj,dxTemp,dyTemp,CBB,PhiBB,CBB0)
                 !else
                 !   call diagotBG(Hjj,dxTemp,dyTemp,CAA,PhiAA,CAA0)
                 !end if
                 if (tBGSwitchDxDy) then
                    if (Species(i).eq.1 .and. layerIndex(i).eq.2) then
                       call diagotBG(Hjj,dyTemp,dxTemp,CApAp,PhiApAp,CApAp0)
                    else if (Species(i).eq.2 .and. layerIndex(i).eq.2) then
                       call diagotBG(Hjj,dyTemp,dxTemp,CBpBp,PhiBpBp,CBpBp0)
                    else if (Species(i).eq.1 .and. layerIndex(i).eq.1) then
                       call diagotBG(Hjj,dyTemp,dxTemp,CAA,PhiAA,CAA0)
                    else if (Species(i).eq.2 .and. layerIndex(i).eq.1) then
                       call diagotBG(Hjj,dyTemp,dxTemp,CBB,PhiBB,CBB0)
                    end if
                 else
                    if (Species(i).eq.1 .and. layerIndex(i).eq.2) then
                       call diagotBG(Hjj,dxTemp,dyTemp,CApAp,PhiApAp,CApAp0)
                    else if (Species(i).eq.2 .and. layerIndex(i).eq.2) then
                       call diagotBG(Hjj,dxTemp,dyTemp,CBpBp,PhiBpBp,CBpBp0)
                    else if (Species(i).eq.1 .and. layerIndex(i).eq.1) then
                       call diagotBG(Hjj,dxTemp,dyTemp,CAA,PhiAA,CAA0)
                    else if (Species(i).eq.2 .and. layerIndex(i).eq.1) then
                       call diagotBG(Hjj,dxTemp,dyTemp,CBB,PhiBB,CBB0)
                    end if
                 end if
                 if (layerIndex(i).eq.1) then
                      if (writeData) then
                          write(588,*) Rat(1,i), Rat(2,i), Hjj
                      end if
                 else if (layerIndex(i).eq.2) then
                      if (writeData) then
                          write(589,*) Rat(1,i), Rat(2,i), Hjj
                      end if
                 end if
                 !Hjj = Hjj/g0
                 !H0(i) = H0(i)+Hjj
                 !if (i.eq.1) then
                 !end if
                 !print*, "Hjj=", Hjj
              else if (BNBNDiag) then
                 if (BNBNuseDisplacementFile) then
                    if (layerIndex(i).eq.1) then
                       dxTemp = -dx/aBN
                       dyTemp = -dy/aBN
                    else if (layerIndex(i).eq.2) then
                       dxTemp = dx/aBN
                       dyTemp = dy/aBN
                       !dxTemp = dx/aG
                       !dyTemp = dy/aG
                    end if
                 else
                    if (layerIndex(i).eq.1) then
                       dxTemp = dx
                       dyTemp = dy
                    else if (layerIndex(i).eq.2) then
                       dxTemp = -dx
                       dyTemp = -dy
                    end if
                 end if
                 !end if
                 !if (Species(i).eq.1 .and. layerIndex(i).eq.2) then
                 !   call diagotBG(Hjj,dxTemp,dyTemp,CBB,PhiBB,CBB0)
                 !else if (Species(i).eq.2 .and. layerIndex(i).eq.1) then
                 !   call diagotBG(Hjj,dxTemp,dyTemp,CBB,PhiBB,CBB0)
                 !else
                 !   call diagotBG(Hjj,dxTemp,dyTemp,CAA,PhiAA,CAA0)
                 !end if
                 !if (tBGSwitchDxDy) then
                 !   if (Species(i).eq.1 .and. layerIndex(i).eq.2) then
                 !      call diagotBG(Hjj,dyTemp,dxTemp,CApAp,PhiApAp,CApAp0)
                 !   else if (Species(i).eq.2 .and. layerIndex(i).eq.2) then
                 !      call diagotBG(Hjj,dyTemp,dxTemp,CBpBp,PhiBpBp,CBpBp0)
                 !   else if (Species(i).eq.1 .and. layerIndex(i).eq.1) then
                 !      call diagotBG(Hjj,dyTemp,dxTemp,CAA,PhiAA,CAA0)
                 !   else if (Species(i).eq.2 .and. layerIndex(i).eq.1) then
                 !      call diagotBG(Hjj,dyTemp,dxTemp,CBB,PhiBB,CBB0)
                 !   end if
                 !else
                    if (Species(i).eq.3 .and. layerIndex(i).eq.2) then
                       call diagotBG(Hjj,dxTemp,dyTemp,CApAp,PhiApAp,CApAp0)
                    else if (Species(i).eq.4 .and. layerIndex(i).eq.2) then
                       call diagotBG(Hjj,dxTemp,dyTemp,CBpBp,PhiBpBp,CBpBp0)
                    else if (Species(i).eq.3 .and. layerIndex(i).eq.1) then
                       call diagotBG(Hjj,dxTemp,dyTemp,CAA,PhiAA,CAA0)
                    else if (Species(i).eq.4 .and. layerIndex(i).eq.1) then
                       call diagotBG(Hjj,dxTemp,dyTemp,CBB,PhiBB,CBB0)
                    end if
                 !end if
                 if (layerIndex(i).eq.1) then
                      if (writeData) then
                          write(588,*) Rat(1,i), Rat(2,i), Hjj
                      end if
                 else if (layerIndex(i).eq.2) then
                      if (writeData) then
                          write(589,*) Rat(1,i), Rat(2,i), Hjj
                      end if
                 end if
                 !Hjj = Hjj/g0
                 !H0(i) = H0(i)+Hjj
                 !if (i.eq.1) then
                 !end if
                 !print*, "Hjj=", Hjj
              else if (GBNuseHarmonicApprox) then
                  if (GBNtwoLayers) then
                      ! Initialize A, B, C to zero to avoid using uninitialized values
                      A = 0.0_dp
                      B = 0.0_dp
                      C = 0.0_dp
                      if (layerIndex(i).eq.1) then ! First layer, graphene
                          if (Species(i).eq.1) then
                              A = -0.37776d0
                              B = -0.39969d0
                              C = -0.4277d0
                          else if (Species(i).eq.2) then
                              A = -0.36739d0
                              B = -0.42011d0
                              C = -0.42048d0
                          end if
                          ! Force compiler to actually use A, B, C by doing an explicit check
                          ! This prevents aggressive optimization from removing the assignments
                          ! Store in a way that compiler cannot optimize away
                          if (abs(A) < 1.0e-10_dp .or. abs(B) < 1.0e-10_dp .or. abs(C) < 1.0e-10_dp) then
                             ! This should never happen for carbon atoms - if it does, we have a serious problem
                             ! Force an error condition that cannot be optimized away
                             if (writeData) then
                                write(588,*) 'CRITICAL ERROR: A,B,C not set for carbon atom', i, &
                                             'Species=', Species(i), 'A=', A, 'B=', B, 'C=', C
                             end if
                             ! Set to a sentinel value to make the error obvious
                             A = 1.0e20_dp
                             B = 1.0e20_dp
                             C = 1.0e20_dp
                          end if
                      else if (layerIndex(i).eq.2) then ! Second layer, BN
                          if (Species(i).eq.3) then
                              B = 2.3222d0
                              C = 2.3476d0
                              A = 2.2731d0
                          else if (Species(i).eq.4) then
                              B = -1.7795d0
                              C = -1.7778d0
                              A = -1.8324d0
                          end if
                      end if
                      !if(GBNtwoLayers) then
                      ! Calculate dxTemp and dyTemp based on layerIndex
                      if (layerIndex(i).eq.1) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.2) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                      ! Check if A, B, C are properly set before calling harmonicApprox
                      ! If they're still zero, the conditions didn't match (wrong Species for this layerIndex)
                      if (abs(A) < 1.0e-10_dp .and. abs(B) < 1.0e-10_dp .and. abs(C) < 1.0e-10_dp) then
                         ! A, B, C are still zero - conditions didn't match
                         ! This should not happen for carbon atoms (Species 1 or 2) in layer 1
                         ! Skip harmonicApprox to avoid non-deterministic behavior
                         Hjj = 0.0_dp
                         if (writeData) then
                            write(588,*) 'WARNING: A,B,C=0 for atom', i, 'Species=', Species(i), 'layerIndex=', layerIndex(i)
                         end if
                      else
                         ! A, B, C are set - safe to call harmonicApprox
                         ! Force compiler to use A, B, C by explicitly passing them
                         ! Store them in local variables that cannot be optimized away
                         ! This is a workaround for compiler optimization issues
                         call harmonicApprox(dxTemp,dyTemp,A,B,C, Hjj)
                         ! For carbon atoms, verify the result is reasonable
                         if (layerIndex(i).eq.1 .and. (Species(i).eq.1 .or. Species(i).eq.2)) then
                            if (abs(Hjj) < 1.0e-10_dp) then
                               ! This is unexpected - harmonicApprox should return non-zero for carbon atoms
                               if (writeData) then
                                  write(588,*) 'WARNING: harmonicApprox returned zero for carbon atom', i, &
                                               'Hjj=', Hjj, 'dxTemp=', dxTemp, 'dyTemp=', dyTemp, &
                                               'A=', A, 'B=', B, 'C=', C, 'dx=', dx, 'dy=', dy
                               end if
                            end if
                         end if
                      end if
                      !end if
                      if (layerIndex(i).eq.1) then
                           if (writeData) then
                               write(588,*) Rat(1,i), Rat(2,i), Hjj
                           end if
                      else if (layerIndex(i).eq.2) then
                           if (writeData) then
                               write(589,*) Rat(1,i), Rat(2,i), Hjj
                           end if
                      end if
                      Hjj = Hjj/g0
                      !if (i.eq.1) then
                      !end if
                      !print*, "Hjj=", Hjj
                  else if (encapsulatedThreeLayers) then
                      ! Initialize A, B, C to zero to avoid using uninitialized values
                      A = 0.0_dp
                      B = 0.0_dp
                      C = 0.0_dp
                      if (layerIndex(i).eq.2) then ! middle layers, graphene
                          if (Species(i).eq.1) then
                              B = -0.39969d0
                              C = -0.4277d0
                              A = -0.37776d0
                          else if (Species(i).eq.2) then
                              B = -0.42011d0
                              C = -0.42048d0
                              A = -0.36739d0
                          end if
                      else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3) then ! hBN
                          if (Species(i).eq.3) then
                              B = 2.3222d0
                              C = 2.3476d0
                              A = 2.2731d0
                          else if (Species(i).eq.4) then
                              B = -1.7795d0
                              C = -1.7778d0
                              A = -1.8324d0
                          end if
                      end if
                      if (GBNuseDisplacementFile) then
                          if (layerIndex(i).eq.2) then
                             dxTemp_b = dx_b/aG
                             dyTemp_b = dy_b/aG
                             dxTemp_t = dx_t/aG
                             dyTemp_t = dy_t/aG
                             !print*, "dxTemp", dxTemp_b, dyTemp_b, dxTemp_t, dyTemp_t
                             !print*, "dx", dx_b, dy_b, dx_t, dy_t
                          else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3) then
                             dxTemp = -dx/aBN
                             dyTemp = -dy/aBN
                          end if
                          !if (layerIndex(i).eq.2) then
                          !   dxTemp = dx/aG
                          !   dyTemp = dy/aG
                          !else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3) then
                          !   dxTemp = -dx/aBN
                          !   dyTemp = -dy/aBN
                          !end if
                      else
                          if (layerIndex(i).eq.2) then
                             dxTemp = dx/aG
                             dyTemp = dy/aG
                             ! If GBNuseDisplacementFile is false, we still need dxTemp_b, dxTemp_t for layer 2
                             ! Use dx, dy as fallback (they're calculated from moire pattern, not from file)
                             dxTemp_b = dx/aG
                             dyTemp_b = dy/aG
                             dxTemp_t = dx/aG
                             dyTemp_t = dy/aG
                          else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3) then
                             dxTemp = -dx/aBN
                             dyTemp = -dy/aBN
                          end if
                      end if
                      if (layerIndex(i).eq.2) then
                          call harmonicApprox(dxTemp_b,dyTemp_b,A,B,C, Hjj_b)
                          call harmonicApprox(dxTemp_t,dyTemp_t,A,B,C, Hjj_t)
                          !Hjj = Hjj_b + (Hjj_t - (A+B+C)/3.0_dp) ! avoid double counting of the average C0 value
                          if (removeTopMoireInL2) then
                              Hjj = (Hjj_b - (A+B+C)/3.0_dp) ! only include the corrections... should also work for single layer F2G2 then
                          else
                              Hjj = (Hjj_b - (A+B+C)/3.0_dp) + (Hjj_t - (A+B+C)/3.0_dp) ! only include the corrections... should also work for single layer F2G2 then
                          end if
                          Hjj = Hjj/g0
                          !print*, "Hjj", Hjj_b, Hjj_t, Hjj
                      else
                          call harmonicApprox(dxTemp,dyTemp,A,B,C, Hjj)
                          Hjj = Hjj - (A+B+C)/3.0_dp
                          Hjj = Hjj/g0
                      end if
                      !print*, "Hjj=", Hjj
                  else if (encapsulatedFourLayers) then
                      ! Initialize A, B, C to zero to avoid using uninitialized values
                      A = 0.0_dp
                      B = 0.0_dp
                      C = 0.0_dp
                      if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then ! middle layers, graphene
                          if (Species(i).eq.1) then
                              B = -0.39969d0
                              C = -0.4277d0
                              A = -0.37776d0
                          else if (Species(i).eq.2) then
                              B = -0.42011d0
                              C = -0.42048d0
                              A = -0.36739d0
                          end if
                      else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then ! hBN
                          if (Species(i).eq.3) then
                              B = 2.3222d0
                              C = 2.3476d0
                              A = 2.2731d0
                          else if (Species(i).eq.4) then
                              B = -1.7795d0
                              C = -1.7778d0
                              A = -1.8324d0
                          end if
                      end if
                      if (GBNuseDisplacementFile) then
                          if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                             dxTemp = dx/aG
                             dyTemp = dy/aG
                          else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                             dxTemp = -dx/aBN
                             dyTemp = -dy/aBN
                          end if
                      else
                          if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                             dxTemp = dx/aG
                             dyTemp = dy/aG
                          else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                             dxTemp = -dx/aBN
                             dyTemp = -dy/aBN
                          end if
                      end if
                      call harmonicApprox(dxTemp,dyTemp,A,B,C, Hjj)
                      Hjj = Hjj/g0
                      !print*, "Hjj=", Hjj
                  else if (encapsulatedFiveLayers) then
                      ! Initialize A, B, C to zero to avoid using uninitialized values
                      A = 0.0_dp
                      B = 0.0_dp
                      C = 0.0_dp
                      if (layerIndex(i).eq.2  .or. layerIndex(i).eq. 4) then ! middle layers, graphene
                          if (Species(i).eq.1) then
                              B = -0.39969d0
                              C = -0.4277d0
                              A = -0.37776d0
                          else if (Species(i).eq.2) then
                              B = -0.42011d0
                              C = -0.42048d0
                              A = -0.36739d0
                          end if
                      else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.5) then ! hBN
                          if (Species(i).eq.3) then
                              B = 2.3222d0
                              C = 2.3476d0
                              A = 2.2731d0
                          else if (Species(i).eq.4) then
                              B = -1.7795d0
                              C = -1.7778d0
                              A = -1.8324d0
                          end if
                      end if
                      if (GBNuseDisplacementFile) then
                          if (layerIndex(i).eq.2 .or. layerIndex(i).eq.4) then
                             dxTemp = dx/aG
                             dyTemp = dy/aG
                          else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.5) then
                             dxTemp = -dx/aBN
                             dyTemp = -dy/aBN
                          end if
                      else
                          if (layerIndex(i).eq.2 .or. layerIndex(i).eq.4) then
                             dxTemp = dx/aG
                             dyTemp = dy/aG
                          else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.5) then
                             dxTemp = -dx/aBN
                             dyTemp = -dy/aBN
                          end if
                      end if
                      call harmonicApprox(dxTemp,dyTemp,A,B,C, Hjj)
                      Hjj = Hjj/g0
                      !print*, "Hjj=", Hjj
                  else if (t3GwithBN) then
                      ! Initialize A, B, C to zero to avoid using uninitialized values
                      A = 0.0_dp
                      B = 0.0_dp
                      C = 0.0_dp
                      if (layerIndex(i).eq.4 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.2) then ! G layers
                          if (Species(i).eq.1) then
                              B = -0.39969d0
                              C = -0.4277d0
                              A = -0.37776d0
                          else if (Species(i).eq.2) then
                              B = -0.42011d0
                              C = -0.42048d0
                              A = -0.36739d0
                          end if
                      else if (layerIndex(i).eq.1) then ! hBN
                          if (Species(i).eq.3) then
                              B = 2.3222d0
                              C = 2.3476d0
                              A = 2.2731d0
                          else if (Species(i).eq.4) then
                              B = -1.7795d0
                              C = -1.7778d0
                              A = -1.8324d0
                          end if
                      end if
                      !if (GBNuseDisplacementFile) then
                          if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.4) then
                             dxTemp = dx/aG
                             dyTemp = dy/aG
                          else if (layerIndex(i).eq.1) then
                             dxTemp = -dx/aBN
                             dyTemp = -dy/aBN
                          end if
                      !else
                      !    if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                      !       dxTemp = dx/aG
                      !       dyTemp = dy/aG
                      !    else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                      !       dxTemp = -dx/aBN
                      !       dyTemp = -dy/aBN
                      !    end if
                      !end if
                           call harmonicApprox(dxTemp,dyTemp,A,B,C, Hjj)
                           Hjj = Hjj/g0
                      !print*, "Hjj=", Hjj
                  else if (BNt2GBN) then
                      ! Initialize A, B, C to zero to avoid using uninitialized values
                      A = 0.0_dp
                      B = 0.0_dp
                      C = 0.0_dp
                      if (layerIndex(i).eq.3 .or. layerIndex(i).eq.2) then ! G layers
                          if (Species(i).eq.1) then
                              B = -0.39969d0
                              C = -0.4277d0
                              A = -0.37776d0
                          else if (Species(i).eq.2) then
                              B = -0.42011d0
                              C = -0.42048d0
                              A = -0.36739d0
                          end if
                      else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then ! hBN
                          if (Species(i).eq.3) then
                              B = 2.3222d0
                              C = 2.3476d0
                              A = 2.2731d0
                          else if (Species(i).eq.4) then
                              B = -1.7795d0
                              C = -1.7778d0
                              A = -1.8324d0
                          end if
                      end if
                      !if (GBNuseDisplacementFile) then
                          if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                             dxTemp = dx/aG
                             dyTemp = dy/aG
                          else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                             dxTemp = -dx/aBN
                             dyTemp = -dy/aBN
                          end if
                      !else
                      !    if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                      !       dxTemp = dx/aG
                      !       dyTemp = dy/aG
                      !    else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                      !       dxTemp = -dx/aBN
                      !       dyTemp = -dy/aBN
                      !    end if
                      !end if
                           call harmonicApprox(dxTemp,dyTemp,A,B,C, Hjj)
                           Hjj = Hjj/g0
                      !print*, "Hjj=", Hjj
                  else if (t2GBN) then
                      !print*, "starting"
                      !print*, layerIndex(i)
                      !print*, Species(i)
                      !print*, dx
                      !print*, dy
                      ! Initialize A, B, C to zero to avoid using uninitialized values
                      A = 0.0_dp
                      B = 0.0_dp
                      C = 0.0_dp
                      if (layerIndex(i).eq.3 .or. layerIndex(i).eq.2) then ! G layers
                          if (Species(i).eq.1) then
                              B = -0.39969d0
                              C = -0.4277d0
                              A = -0.37776d0
                          else if (Species(i).eq.2) then
                              B = -0.42011d0
                              C = -0.42048d0
                              A = -0.36739d0
                          end if
                      else if (layerIndex(i).eq.1) then ! hBN
                          if (Species(i).eq.3) then
                              B = 2.3222d0
                              C = 2.3476d0
                              A = 2.2731d0
                          else if (Species(i).eq.4) then
                              B = -1.7795d0
                              C = -1.7778d0
                              A = -1.8324d0
                          end if
                      end if
                      !if (GBNuseDisplacementFile) then
                          if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                             dxTemp = dx/aG
                             dyTemp = dy/aG
                          else if (layerIndex(i).eq.1) then
                             dxTemp = -dx/aBN
                             dyTemp = -dy/aBN
                          else
                             print*, "we have atoms in the wrong layer"
                          end if
                      !else
                      !    if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                      !       dxTemp = dx/aG
                      !       dyTemp = dy/aG
                      !    else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                      !       dxTemp = -dx/aBN
                      !       dyTemp = -dy/aBN
                      !    end if
                      !end if
                           call harmonicApprox(dxTemp,dyTemp,A,B,C, Hjj)
                           Hjj = Hjj/g0
                      !print*, "Hjj=", Hjj
                  else
                      if (layerIndex(i).eq.1) then ! First layer, graphene
                          if (Species(i).eq.1) then
                              call diagoGBN(Hjj,dx,dy,CApAp,PhiApAp,CApAp0)
                          else if (Species(i).eq.2) then
                              call diagoGBN(Hjj,dx,dy,CBpBp,PhiBpBp,CBpBp0)
                          end if
                      else ! second layer, hBN
                          if (Species(i).eq.3) then ! B
                              call diagoGBN(Hjj,dx,dy,CAA,PhiAA,CAA0)
                          else if (Species(i).eq.4) then ! N
                              call diagoGBN(Hjj,dx,dy,CBB,PhiBB,CBB0)
                          end if
                      end if
                  end if
              end if
              H0(i) = H0(i) + Hjj
            end do
            !$OMP END PARALLEL DO
         !end if

      else if (l) then
         call MIO_Print('Moire Trilayer','ham')
         ! L1 = n1*aG = m1*lambda1
         ! L = n2*L1 = m2*lambda2
         call MIO_InputParameter('MoirePotA1',A1,A)
         call MIO_InputParameter('MoirePotB1',B1,B)
         call MIO_InputParameter('MoirePotA2',A2,A)
         call MIO_InputParameter('MoirePotB2',B2,B)
         call MIO_InputParameter('MoireTrin1',n1,0)
         call MIO_InputParameter('MoireTrim1',m1,0)
         call MIO_InputParameter('MoireTrin2',n2,0)
         call MIO_InputParameter('MoireTrim2',m2,0)
         call MIO_Print('  A1: '//trim(num2str(A1,4)),'ham')
         call MIO_Print('  B1: '//trim(num2str(B1,4)),'ham')
         call MIO_Print('  A2: '//trim(num2str(A2,4)),'ham')
         call MIO_Print('  B2: '//trim(num2str(B2,4)),'ham')
         call MIO_Print('  C : '//trim(num2str(C,4)),'ham')
         A1 = A1/g0
         B1 = B1/g0
         A2 = A2/g0
         B2 = B2/g0
         C = C/g0
         lambda1 = n1*aG/real(m1,dp)
         lambda2 = n2*n1*aG/real(m2,dp)
         call MIO_Print('lambda 1: '//trim(num2str(lambda1,5))//' Ang','ham')
         call MIO_Print('lambda 2: '//trim(num2str(lambda2,5))//' Ang','ham')
         !$OMP PARALLEL DO PRIVATE(i)
         do i=1,nAt
            H0(i) = H0(i) - (-1.0_dp)**Species(i)*(A1*sin(twopi*Rat(1,i)*n2*m1) + &
              B1*sin(twopi*Rat(2,i)*n2*m1) +  A2*sin(twopi*Rat(1,i)*m2) + &
              B2*sin(twopi*Rat(2,i)*m2) + C)/2.0_dp
         end do
         !$OMP END PARALLEL DO
      else
         call MIO_InputParameter('MoireSachs',l,.false.)
         if (l) then
           call MIO_Print('Actually adding the Sachs potential...','ham')
           call MIO_Print('  A: '//trim(num2str(A,4)),'ham')
           call MIO_Print('  B: '//trim(num2str(B,4)),'ham')
           call MIO_Print('  C: '//trim(num2str(C,4)),'ham')
           A = A/g0
           B = B/g0
           C = C/g0
           !$OMP PARALLEL DO PRIVATE(i)
           do i=1,nAt
              H0(i) = H0(i) - (-1.0_dp)**Species(i)*(A*sin(twopi*Rat(1,i)*sCell+phi1) + &
                B*sin(twopi*Rat(2,i)*sCell+phi2) + C)/2.0_dp
           end do
           !$OMP END PARALLEL DO
         end if
         call MIO_InputParameter('MoireJeil',l,.false.)
         call MIO_InputParameter('singleLayerXYZ',singleLayerXYZ,.false.)
         if (l) then
           call MIO_Print('Actually adding the Jeil potential...','ham')
           call MIO_InputParameter('TypeOfSystem',str,'Graphene')
           call MIO_InputParameter('TrilayerFanZhang',FanZhang,.false.)
           if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell') .or. (MIO_StringComp(str,'ReadXYZ') .and. .not.(singleLayerXYZ))) then
             call MIO_Print('Working on MoireEncapsulatedBilayer or MoireEncapsulatedBilayerBasedOnMoireCell or ReadXYZ but not a single layer','ham')
             call MIO_InputParameter('MoireOnlyH0',l,.false.)
             call MIO_InputParameter('MoireOnlyHZ',w,.false.)
             call MIO_InputParameter('MoireNoH0AndHZ',v,.false.)
             call MIO_InputParameter('MoireH0AndHZ',u,.false.)
             call MIO_InputParameter('MoireBilayerTopAngle',MoireBilayerTopAngle,0.0_dp)
             MoireBilayerTopAngleGrad = MoireBilayerTopAngle*pi/180.0_dp
             call MIO_InputParameter('MoireBilayerBottomAngle',MoireBilayerBottomAngle,0.0_dp)
             MoireBilayerBottomAngleGrad = MoireBilayerBottomAngle*pi/180.0_dp
             !print*, "quick check:"
             !print*, "cos(MoireBilayerTopAngleGrad) - 1)",cos(MoireBilayerTopAngleGrad) - 1
             !print*, "cos(MoireBilayerTopAngleGrad) - 1.0)",cos(MoireBilayerTopAngleGrad) - 1.0
             call MIO_InputParameter('MoireBLDeactivateUpperLayer',deactivateUpperLayer,.false.)
             call MIO_InputParameter('MoiretDBLDeactivateUpperLayers',deactivateUpperLayers,.false.)
             if (frac) call AtomsSetCart()
             !print*, "taus =", tauX1, tauY1, tauX2, tauY2
             !$OMP PARALLEL DO PRIVATE(i,dx,dy,Hzjj,H0jj,C0d,Czd)
             do i=1,nAt
               !if (i<=nAtC1) then
               !if (Rat(3,i).lt.20.0_dp) then
               if (layerIndex(i).eq.1) then ! First layer, no rotation
                   !dx = (Rat(1,i)+ tauX1)*eps
                   !dy = (Rat(2,i)+ tauY1)*eps
                   dx = ((1.0_dp+eps) * cos(MoireBilayerBottomAngleGrad) - 1) * (Rat(1,i)) &
                             - ((1.0_dp+eps) * sin(MoireBilayerBottomAngleGrad) * (Rat(2,i)))
                   dy = ((1.0_dp+eps) * cos(MoireBilayerBottomAngleGrad) - 1) * (Rat(2,i)) &
                             + ((1.0_dp+eps) * sin(MoireBilayerBottomAngleGrad) * (Rat(1,i)))
                   dx = dx + tauX1
                   dy = dy + tauY1
                   if (distanceDependentEffectiveModel) then
                      call distanceDependentC(C0d, C0, BfactorC0, interlayerDistances(i), z0) ! it should be ok to take the first neighbor as it is effective model with only 3 neighbors
                      call distanceDependentC(Czd, Cz, BfactorCz, interlayerDistances(i), z0)
                   else
                      C0d = C0
                      Czd = Cz
                   end if
                   if (addDisplacements) then
                       dx = dx + displacements(1,i)
                       dy = dy + displacements(2,i)
                   end if
                   if (l) then
                     Hzjj = 0.0_dp
                     if (sublatticeBasis) then
                        call diagoH0FromAABB(H0jj,dx,dy,interlayerDistances(i),moirePreFactor)
                     else
                        call diago(H0jj,dx,dy,C0d,Phi0)
                     end if
                   else if (w) then
                     if (sublatticeBasis) then
                        call diagoHzFromAABB(Hzjj,dx,dy,interlayerDistances(i),moirePreFactor)
                        !print*, "using the sublattice dependent basis with Hz"
                     else
                        call diago(Hzjj,dx,dy,Czd,Phiz)
                     end if
                     H0jj = 0.0_dp
                   else if (v) then
                     Hzjj = 0.0_dp
                     H0jj = 0.0_dp
                   else if (u) then
                     if (sublatticeBasis) then
                        call diagoH0FromAABB(H0jj,dx,dy,interlayerDistances(i),moirePreFactor)
                        call diagoHzFromAABB(Hzjj,dx,dy,interlayerDistances(i),moirePreFactor)
                     else
                        call diago(Hzjj,dx,dy,Czd,Phiz)
                        call diago(H0jj,dx,dy,C0d,Phi0)
                     end if
                   end if
                   H0(i) = H0(i) + (-1.0_dp)**Species(i)*Hzjj + H0jj
               else
                   !dx = ((1.0_dp+eps) * cos(MoireBilayerTopAngleGrad) - 1) * (Rat(1,i)+tauX2) &
                   !          - ((1.0_dp+eps) * sin(MoireBilayerTopAngleGrad) * (Rat(2,i)+tauY2))
                   !dy = ((1.0_dp+eps) * cos(MoireBilayerTopAngleGrad) - 1) * (Rat(2,i)+tauY2) &
                   !          + ((1.0_dp+eps) * sin(MoireBilayerTopAngleGrad) * (Rat(1,i)+tauX2))
                   if (rotateFirst) then
                      dx = ((1.0_dp+eps) * cos(MoireBilayerTopAngleGrad) - 1) * (Rat(1,i)) &
                                - ((1.0_dp+eps) * sin(MoireBilayerTopAngleGrad) * (Rat(2,i)))
                      dy = ((1.0_dp+eps) * cos(MoireBilayerTopAngleGrad) - 1) * (Rat(2,i)) &
                                + ((1.0_dp+eps) * sin(MoireBilayerTopAngleGrad) * (Rat(1,i)))
                      dx = dx + tauX2
                      dy = dy + tauY2 !- aG/sqrt(3.0_dp)
                   else
                      dx = ((1.0_dp+eps) * cos(MoireBilayerTopAngleGrad) - 1) * (Rat(1,i)) &
                                - ((1.0_dp+eps) * sin(MoireBilayerTopAngleGrad) * (Rat(2,i) - tauY2/eps))
                      dy = ((1.0_dp+eps) * cos(MoireBilayerTopAngleGrad) - 1) * (Rat(2,i) - tauY2/eps) &
                                + ((1.0_dp+eps) * sin(MoireBilayerTopAngleGrad) * (Rat(1,i)))
                   end if
                   !if (deactivateUpperLayer) then
                   !if (Rat(3,i) .gt. 20.0_dp .and. deactivateUpperLayer) then
                   if (layerIndex(i).eq.2 .and. deactivateUpperLayer) then
                     Hzjj = 0.0_dp
                     H0jj = 0.0_dp
                   else if (layerIndex(i).ne.1 .and. deactivateUpperLayers) then
                     Hzjj = 0.0_dp
                     H0jj = 0.0_dp
                   else
                     if (distanceDependentEffectiveModel) then
                      call distanceDependentC(C0d, C0, BfactorC0, interlayerDistances(i), z0) ! it should be ok to take the first neighbor as it is effective model with only 3 neighbors
                      call distanceDependentC(Czd, Cz, BfactorCz, interlayerDistances(i), z0)
                     else
                        C0d = C0
                        Czd = Cz
                     end if
                     if (addDisplacements) then
                         dx = dx + displacements(1,i)
                         dy = dy + displacements(2,i)
                     end if
                     if (l) then
                       Hzjj = 0.0_dp
                       if (sublatticeBasis) then
                          call diagoH0FromAABB(H0jj,dx,dy, interlayerDistances(i),moirePreFactor)
                       else
                          call diago(H0jj,dx,dy,C0d,Phi0)
                       end if
                     else if (w) then
                       if (sublatticeBasis) then
                          call diagoHzFromAABB(Hzjj,dx,dy, interlayerDistances(i),moirePreFactor)
                          !print*, "using the sublattice dependent basis with Hz"
                       else
                          call diago(Hzjj,dx,dy,Czd,Phiz)
                       end if
                       H0jj = 0.0_dp
                     else if (v) then
                       Hzjj = 0.0_dp
                       H0jj = 0.0_dp
                     else if (u) then
                       if (sublatticeBasis) then
                          call diagoH0FromAABB(H0jj,dx,dy, interlayerDistances(i),moirePreFactor)
                          call diagoHzFromAABB(Hzjj,dx,dy, interlayerDistances(i),moirePreFactor)
                       else
                          call diago(Hzjj,dx,dy,Czd,Phiz)
                          call diago(H0jj,dx,dy,C0d,Phi0)
                       end if
                     end if
                   end if
                   H0(i) = H0(i) + (-1.0_dp)**Species(i)*Hzjj + H0jj
               end if
             end do
             !$OMP END PARALLEL DO
           !else if (MIO_StringComp(str,'TwistedBilayer')) then ! I don't think it's at the right place
           !        ! add onsite energies
           else
             call MIO_Print('Working on ANY moire system not satisfying MoireEncapsulatedBilayer or MoireEncapsulatedBilayerBasedOnMoireCell or ReadXYZ but not a single layer','ham')
             call MIO_InputParameter('MoireOnlyH0',l,.false.)
             call MIO_InputParameter('MoireOnlyHZ',w,.false.)
             call MIO_InputParameter('MoireNoH0AndHZ',v,.false.)
             call MIO_InputParameter('MoireH0AndHZ',u,.false.)
             call MIO_InputParameter('MoireKekule',z,.false.)
             call MIO_InputParameter('MoireAddSecondMoire',zz,.false.)
             call MIO_InputParameter('LatticepercentFactor',epsFactor,1.0_dp)
             call MIO_InputParameter('basedOnMoireCellParameters',ll,.false.)
             call MIO_InputParameter('switchHzjj',switchHzjj,.false.)

             if (ll) then
                 !print*, "phiForEffectiveModel: ", phiForEffectiveModel
                 twistAngleGrad = -phiForEffectiveModel*pi/180.0_dp
                 twistAngleGrad2 = phiForEffectiveModel*pi/180.0_dp
                 call MIO_InputParameter('MoireCellParameters',nnn,[0,0,0,0])
                 !eps = 1.0/nnn(1)
                 !eps2 = 1.0/nnn(3)
                 !eps = -1.0/(nnn(1))
                 !eps2 = -1.0/(nnn(1))
                 !eps = 0.0-0.02
                 !eps2 = 1.0/(nnn(2))
                 !ggg = nnn(1)**2 + nnn(2)**2 + nnn(1)*nnn(2)
                 !print*, "ggg", ggg
                 !eps = dsqrt(real(nnn(3)**2 + nnn(4)**2 + nnn(3)*nnn(4))/ggg)
                 !print*, "calculatedEps", eps
                 !eps = -0.02463661_dp
                 !eps2 = eps
                 !print*, "calculatedPhi", acos((2.0d0*nnn(1)*nnn(3)+2.0d0*nnn(2)*nnn(4) + nnn(1)*nnn(4) + nnn(2)*nnn(3))/(2.0d0*eps*ggg))
             else
                 call MIO_InputParameter('MoireTwistAngle',twistAngle,0.0_dp)
                 call MIO_InputParameter('MoireTwistAngle2',twistAngle2,0.0_dp)
                 twistAngleGrad = twistAngle*pi/180.0_dp
                 twistAngleGrad2 = twistAngle2*pi/180.0_dp
             end if
             call MIO_InputParameter('MoireTwisted',twisted,.false.)
             call MIO_InputParameter('MoireTwisted2',twisted2,.false.)
             call MIO_InputParameter('deactivateASubLattice',deactivateASubLattice,.false.)
             call MIO_InputParameter('deactivateBSubLattice',deactivateBSubLattice,.false.)
             call MIO_InputParameter('MoireLayerShift1',shift1,0.0_dp)
             call MIO_InputParameter('MoireLayerShift2',shift2,0.0_dp)
             phi2M= atan((1.0_dp+eps)*sin(twistAngleGrad2)/((1.0_dp+eps)*cos(twistAngleGrad2)-1.0_dp))
             shift2_x = shift2*cos(phi2M)
             shift2_y = shift2*sin(phi2M)
             call MIO_InputParameter('MoireSecondMoireRotateFirst',rotateFirst,.true.)
             call MIO_InputParameter('MoireSecondMoireMassFactor',sign2,1.0_dp)
             call MIO_InputParameter('MoireFirstMoireMassFactor',sign1,1.0_dp)
             call MIO_InputParameter('minZ',minZ,1.0_dp)
             !print*, "changing layer 1 or layer from BN to NB", sign2, sign1
             if (frac) call AtomsSetCart()
             !$OMP PARALLEL DO PRIVATE(i,dx,dy,Hzjj,H0jj,C0d,Czd)
             do i=1,nAt ! Konda, modify onsite energies
                if (twisted) then
                    dx = Rat(1,i) * ((1.0_dp+eps) * cos(twistAngleGrad) - 1) &
                         - Rat(2,i) * (1.0_dp+eps) * sin(twistAngleGrad)
                    dy = Rat(2,i) * ((1.0_dp+eps) * cos(twistAngleGrad) - 1) &
                         + Rat(1,i) * (1.0_dp+eps) * sin(twistAngleGrad)
                else
                    if (readRigidXYZ) then
                       dx = RatInit(1,i)*eps
                       dy = RatInit(2,i)*eps
                    else
                       dx = Rat(1,i)*eps
                       dy = Rat(2,i)*eps
                    end if
                end if
                if (.not. shift1==0.0_dp) then
                     call MIO_Print('we add a shift to the bottom moire layer','ham')
                end if
                dy = dy + shift1
                !if (sign1.lt.0.0) then ! First go from BN to NB orientation
                !    dx = -dx
                !    dy = -dy
                !end if
                if (layerIndex(i).eq.1) then
                    if (sign1.lt.0.0) then ! Go from BN to NB orientation if necessary
                        dx = -dx
                        dy = -dy
                    end if
                else if (layerIndex(i).eq.2) then
                    if (sign2.lt.0.0) then ! Go from BN to NB orientation for second layer of bilayer graphene
                        dx = -dx
                        dy = -dy
                    end if
                end if
                if (addDisplacements) then ! Use the displacements from matlab or lammps
                    dx = dx + displacements(1,i)
                    dy = dy + displacements(2,i)
                end if
                if (distanceDependentEffectiveModel) then
                   call distanceDependentC(C0d, C0, BfactorC0, interlayerDistances(i), z0)
                   call distanceDependentC(Czd, Cz, BfactorCz, interlayerDistances(i), z0)
                else
                   C0d = C0
                   Czd = Cz
                end if
                if (l) then
                  Hzjj = 0.0_dp
                  if (sublatticeBasis) then
                     call diagoH0FromAABB(H0jj,dx,dy,interlayerDistances(i),moirePreFactor)
                  else
                     call diago(H0jj,dx,dy,C0d,Phi0)
                  end if
                  !end if
                else if (w) then
                  if (sublatticeBasis) then
                     call diagoHzFromAABB(Hzjj,dx,dy,interlayerDistances(i),moirePreFactor)
                     !print*, "using the sublattice dependent basis with Hz"
                  else
                     call diago(Hzjj,dx,dy,Czd,Phiz)
                  end if
                  !end if
                  H0jj = 0.0_dp
                else if (v) then
                  Hzjj = 0.0_dp
                  H0jj = 0.0_dp
                else if (u) then
                  if (sublatticeBasis) then
                     call diagoH0FromAABB(H0jj,dx,dy,interlayerDistances(i),moirePreFactor)
                     call diagoHzFromAABB(Hzjj,dx,dy,interlayerDistances(i),moirePreFactor)
                     !print*, "using the sublattice dependent basis with Hz"
                  else
                     call diago(Hzjj,dx,dy,Czd,Phiz)
                     call diago(H0jj,dx,dy,C0d,Phi0)
                  end if
                end if
                if (deactivateASubLattice .and. Species(i).eq.2) then
                    Hzjj = 0
                    H0jj = 0
                else if (deactivateBSubLattice .and. Species(i).eq.1) then
                    Hzjj = 0
                    H0jj = 0
                end if
                if (FanZhang) then
                    if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                       Hzjj = 0.0_dp
                       H0jj = 0.0_dp
                    end if
                end if
                    if (layerIndex(i).eq.1) then
                       H0(i) = H0(i) + (-1.0_dp)**Species(i)*Hzjj * sign1 + H0jj
                    else if (layerIndex(i).eq.2) then
                       H0(i) = H0(i) + (-1.0_dp)**Species(i)*Hzjj * sign2 + H0jj
                    end if
             end do
             !$OMP END PARALLEL DO
             H0Bottom = H0
             if (zz) then ! only used when adding second moire to the same carbon atoms
                !$OMP PARALLEL DO PRIVATE(i,dx,dy,Hzjj,H0jj,C0d,Czd)
                do i=1,nAt
                   if (ll) then
                       if (twisted2) then
                           dx = Rat(1,i) * ((1.0_dp+eps2) * cos(twistAngleGrad2) - 1.0_dp) &
                                  - (Rat(2,i)) * (1.0_dp+eps2) * sin(twistAngleGrad2)
                           dy = (Rat(2,i)) * ((1.0_dp+eps2) * cos(twistAngleGrad2) - 1.0_dp) &
                                  + Rat(1,i) * (1.0_dp+eps2) * sin(twistAngleGrad2)
                       else
                           dx = Rat(1,i)*eps2
                           dy = Rat(2,i)*eps2
                       end if
                   else
                       if (twisted2) then
                           if (rotateFirst) then
                               dx = Rat(1,i) * ((1.0_dp+eps*epsFactor) * cos(twistAngleGrad2) - 1.0_dp) &
                                      - (Rat(2,i)) * (1.0_dp+eps*epsFactor) * sin(twistAngleGrad2)
                               dy = (Rat(2,i)) * ((1.0_dp+eps*epsFactor) * cos(twistAngleGrad2) - 1.0_dp) &
                                      + Rat(1,i) * (1.0_dp+eps*epsFactor) * sin(twistAngleGrad2)
                           else
                               dx = Rat(1,i)  * ((1.0_dp+eps*epsFactor) * cos(twistAngleGrad2) - 1.0_dp) &
                                      - (Rat(2,i) - shift2/eps) * (1.0_dp+eps*epsFactor) * sin(twistAngleGrad2)
                               dy = (Rat(2,i) - shift2/eps) * ((1.0_dp+eps*epsFactor) * cos(twistAngleGrad2) - 1.0_dp) &
                                      + Rat(1,i) * (1.0_dp+eps*epsFactor) * sin(twistAngleGrad2)
                           end if
                       else
                           dx = Rat(1,i)*eps*epsFactor
                           dy = Rat(2,i)*eps*epsFactor
                       end if
                   end if
                   if (rotateFirst) then
                       dy = dy + shift2_y
                       dx = dx + shift2_x
                   end if
                   if (distanceDependentEffectiveModel) then
                      call distanceDependentC(C0d, C0, BfactorC0, interlayerDistances(i), z0)
                      call distanceDependentC(Czd, Cz, BfactorCz, interlayerDistances(i), z0)
                   else
                      C0d = C0
                      Czd = Cz
                   end if
                   if (addDisplacements) then
                       dx = dx + displacements(1,i)
                       dy = dy + displacements(2,i)
                   end if
                   if (sign2.lt.0.0) then
                       dx = -dx
                       dy = -dy
                   end if
                   if (l) then
                     Hzjj = 0.0_dp
                     if (sublatticeBasis) then
                        call diagoH0FromAABB(H0jj,dx,dy,interlayerDistances(i),moirePreFactor)
                     else
                         call diago(H0jj,dx,dy,C0d,Phi0)
                     end if
                   else if (w) then
                     if (sublatticeBasis) then
                        call diagoHzFromAABB(Hzjj,dx,dy,interlayerDistances(i),moirePreFactor)
                        !print*, "using the sublattice dependent basis with Hz"
                     else
                        call diago(Hzjj,dx,dy,Czd,Phiz)
                     end if
                     H0jj = 0.0_dp
                   else if (v) then
                     Hzjj = 0.0_dp
                     H0jj = 0.0_dp
                   else if (u) then
                     if (sublatticeBasis) then
                        call diagoHzFromAABB(Hzjj,dx,dy, interlayerDistances(i),moirePreFactor)
                        call diagoH0FromAABB(H0jj,dx,dy, interlayerDistances(i),moirePreFactor)
                     else
                        call diago(Hzjj,dx,dy,Czd,Phiz)
                        call diago(H0jj,dx,dy,C0d,Phi0)
                     end if
                   end if
                   H0(i) = H0(i) + (-1.0_dp)**Species(i)*Hzjj * sign2 + H0jj
                end do
                !$OMP END PARALLEL DO
                H0Top = H0-H0Bottom
             end if
             if (z) then
                if (frac) call AtomsSetCart()
                call MIO_InputParameter('MoirePotC',C,0.0_dp)
                call MIO_InputParameter('MoireKekuleAngle',KekuleAngle,30.0_dp)
                KekuleAngleGrad = KekuleAngle*pi/180.0_dp
                call MIO_InputParameter('MoireKekuleEpsFactor',KekuleEpsFactor,sqrt(3.0_dp))
                !$OMP PARALLEL DO PRIVATE(i,dx,dy,dxPrime,dyPrime,Hkjj,C0d,Czd)
                do i=1,nAt
                  if (distanceDependentEffectiveModel) then
                     call distanceDependentC(C0d, C0, BfactorC0, interlayerDistances(i), z0)
                     call distanceDependentC(Czd, Cz, BfactorCz, interlayerDistances(i), z0)
                  else
                     C0d = C0
                     Czd = Cz
                  end if
                  dx = Rat(1,i) * ((1.0_dp+eps/KekuleEpsFactor) * cos(KekuleAngleGrad) - 1) &
                         - Rat(2,i) * (1.0_dp+eps/KekuleEpsFactor) * sin(KekuleAngleGrad)
                  dy = Rat(2,i) * ((1.0_dp+eps/KekuleEpsFactor) * cos(KekuleAngleGrad) - 1) &
                         + Rat(1,i) * (1.0_dp+eps/KekuleEpsFactor) * sin(KekuleAngleGrad)
                  if (addDisplacements) then
                      dx = dx + displacements(1,i)
                      dy = dy + displacements(2,i)
                  end if
                  if (sublatticeBasis) then
                     call diagoH0FromAABB(Hkjj,dx,dy,interlayerDistances(i),moirePreFactor)
                  else
                     call diago(Hkjj,dx,dy,C0d,Phi0)
                  end if
                  H0(i) = H0(i) + Hkjj
                end do
                !$OMP END PARALLEL DO
             end if
           end if
           if (.not. frac) call AtomsSetFrac()
         end if
         call MIO_Print('')
      end if
      call MIO_InputParameter('MoireSymmetricPot',l,.false.)
      if (l) then
         call MIO_InputParameter('MoireSymmA',A,0.01_dp)
         call MIO_InputParameter('MoireSymmB',B,0.01_dp)
         call MIO_print('Symmetric part:')
         call MIO_Print('  A: '//trim(num2str(A,4)),'ham')
         call MIO_Print('  B: '//trim(num2str(B,4)),'ham')
         A = A/g0
         B = B/g0
         !$OMP PARALLEL DO PRIVATE(i)
         do i=1,nAt
            H0(i) = H0(i) + (A*sin(twopi*Rat(1,i)*sCell) + &
              B*sin(twopi*Rat(2,i)*sCell))/2.0_dp
         end do
         !$OMP END PARALLEL DO
      end if
   end if
   call MIO_InputParameter('MoireBilayerElectricField',u,.false.)
   call MIO_InputParameter('MoireBilayerElectricShift',ElectricShift,0.150_dp)
   !ElectricShift = ElectricShift/g0
   call MIO_InputParameter('MoireBilayerElectricFieldInvert',invertE,.false.)
   !call MIO_InputParameter('MoireBilayerInvertELimitX',invertEposX,0.5)
   call MIO_InputParameter('CellSize',n,55)
   invertEposX = (n*sCell*aG)/2.0
   !print*, "switching the electric field at ", invertEposX
   ElectricShift = ElectricShift/g0
  call MIO_InputParameter('fourLayers',fourLayers,.false.)
  call MIO_InputParameter('fourLayersSandwiched',fourLayersSandwiched,.false.)
  call MIO_InputParameter('helicalTwistedMBM',helicalTwistedMBM,.false.)
  call MIO_InputParameter('fiveLayersSandwiched',fiveLayersSandwiched,.false.)
   call MIO_InputParameter('sixLayersSandwiched',sixLayersSandwiched,.false.)
   call MIO_InputParameter('sevenLayersSandwiched',sevenLayersSandwiched,.false.)
   call MIO_InputParameter('eightLayersSandwiched',eightLayersSandwiched,.false.)
   call MIO_InputParameter('tenLayersSandwiched',tenLayersSandwiched,.false.)
   call MIO_InputParameter('twentyLayersSandwiched',twentyLayersSandwiched,.false.)
   if (u) then
      call MIO_Print('Adding a energy shift beteen both layers in the bilayer graphene','ham')
      if (frac) call AtomsSetCart()
      !$OMP PARALLEL DO PRIVATE(i, EFactor)
      do i=1,nAt
         if (invertE) then
            if (Rat(1,i).gt.invertEposX) then
                Efactor = -1
            else
                Efactor = 1
            end if
         else
            Efactor = 1
         end if
         !if (i<=nAtC1) then
         !if (Rat(3,i).lt.20.0_dp) then
         if (fourLayers .or. fourLayersSandwiched .or. helicalTwistedMBM .or. BNt2GBN) then
            if (layerIndex(i).eq.1) then
               H0(i) = H0(i) - ElectricShift*Efactor/2.0_dp
            else if (layerIndex(i).eq.2) then
               H0(i) = H0(i) - ElectricShift*Efactor/6.0_dp
            else if (layerIndex(i).eq.3) then
               H0(i) = H0(i) + ElectricShift*Efactor/6.0_dp
            else if (layerIndex(i).eq.4) then
               H0(i) = H0(i) + ElectricShift*Efactor/2.0_dp
            end if
         else if (threeLayers) then
            if (layerIndex(i).eq.1) then
               H0(i) = H0(i) - ElectricShift*Efactor/2.0_dp
            else if (layerIndex(i).eq.3) then
               H0(i) = H0(i) + ElectricShift*Efactor/2.0_dp
            else
               H0(i) = H0(i)
            end if
         else if (fiveLayersSandwiched) then
            if (layerIndex(i).eq.1) then
               H0(i) = H0(i) - ElectricShift*Efactor/2.0_dp
            else if (layerIndex(i).eq.2) then
               H0(i) = H0(i) - ElectricShift*Efactor/4.0_dp
            else if (layerIndex(i).eq.4) then
               H0(i) = H0(i) + ElectricShift*Efactor/4.0_dp
            else if (layerIndex(i).eq.5) then
               H0(i) = H0(i) + ElectricShift*Efactor/2.0_dp
            else
               H0(i) = H0(i)
            end if
         else if (sixLayersSandwiched) then
            if (layerIndex(i).eq.1) then
               H0(i) = H0(i) - ElectricShift*Efactor/2.0_dp
            else if (layerIndex(i).eq.2) then
               H0(i) = H0(i) - ElectricShift*Efactor*0.3_dp
            else if (layerIndex(i).eq.3) then
               H0(i) = H0(i) - ElectricShift*Efactor*0.1_dp
            else if (layerIndex(i).eq.4) then
               H0(i) = H0(i) + ElectricShift*Efactor*0.1_dp
            else if (layerIndex(i).eq.5) then
               H0(i) = H0(i) + ElectricShift*Efactor*0.3_dp
            else if (layerIndex(i).eq.6) then
               H0(i) = H0(i) + ElectricShift*Efactor/2.0_dp
            end if
         else if (t3BG) then
            if (layerIndex(i).eq.1) then
               H0(i) = H0(i) - ElectricShift*Efactor/2.0_dp
            else if (layerIndex(i).eq.2) then
               H0(i) = H0(i) - ElectricShift*Efactor/2.0_dp + ElectricShift*Efactor/5.0_dp
            else if (layerIndex(i).eq.3) then
               H0(i) = H0(i) - ElectricShift*Efactor/2.0_dp + ElectricShift*Efactor*2.0_dp/5.0_dp
            else if (layerIndex(i).eq.4) then
               H0(i) = H0(i) + ElectricShift*Efactor/2.0_dp - ElectricShift*Efactor*2.0_dp/5.0_dp
            else if (layerIndex(i).eq.5) then
               H0(i) = H0(i) + ElectricShift*Efactor/2.0_dp - ElectricShift*Efactor/5.0_dp
            else if (layerIndex(i).eq.6) then
               H0(i) = H0(i) + ElectricShift*Efactor/2.0_dp
            end if
         else if (twentyLayersSandwiched) then
            do iii=1,20
               if (layerIndex(i).eq.iii) then
                  H0(i) = H0(i) + (ElectricShift*Efactor*(-0.5 + (iii-1)*1.0_dp/(20.0_dp-1.0_dp)))
               end if
            end do
         else
            if (layerIndex(i).eq.1) then
               H0(i) = H0(i) - ElectricShift*Efactor/2.0_dp
            else
               H0(i) = H0(i) + ElectricShift*Efactor/2.0_dp
            end if
         end if
      end do
      !$OMP END PARALLEL DO
   end if
   call MIO_InputParameter('TrilayerFanZhang',FanZhang,.false.)
   if (FanZhang) then
      do i=1,nAt
       call MIO_InputParameter('Trilayeru1',u1,0.0_dp)
       call MIO_InputParameter('Trilayeru2',u2,0.0_dp)
       call MIO_InputParameter('Trilayeru3',u3,0.0_dp)
       call MIO_InputParameter('TrilayerDelta',trDelta,0.0_dp)
       u1 = u1/g0
       u2 = u2/g0
       u3 = u3/g0
       trDelta = trDelta/g0
       if (layerIndex(i).eq.1) then
          H0(i) = H0(i) + u1 + trDelta
       else if (layerIndex(i).eq.2) then
          H0(i) = H0(i) + u2
       else if (layerIndex(i).eq.3) then
          H0(i) = H0(i) + u3
       end if
       !tAB1 = -tAB1/g0 ! Add minus sign to compensate for intrinsic minus sign
      end do
   !else
   !    print*, "You didnt specify any trilayer parameters, is that correct?"
   end if

   call MIO_InputParameter('addSublatticeMassterm',u,.false.)
   call MIO_InputParameter('onlyBottomLayerMassTerm',onlyBottomLayerMassTerm,.false.)
   call MIO_InputParameter('sublatticeMassterm',massterm,0.150_dp)
   massterm = massterm/g0
   if (u) then
      if (onlyBottomLayerMassTerm) then
          call MIO_Print('Adding a massterm between sublattices, but only for the bottom layer, i.e. layerIndex == 1','ham')
      else
          call MIO_Print('Adding a massterm between sublattices','ham')
      end if
      !$OMP PARALLEL DO PRIVATE(i)
      do i=1,nAt
         if (onlyBottomLayerMassTerm) then
             if (layerIndex(i)==1) then
                 H0(i) = H0(i) + (-1.0_dp)**Species(i)*massterm/2.0_dp
             end if
         else
             H0(i) = H0(i) + (-1.0_dp)**Species(i)*massterm/2.0_dp
         end if
      end do
      !$OMP END PARALLEL DO
   end if

   call MIO_InputParameter('addOnsiteEnergyShift',u,.false.)
   call MIO_InputParameter('onsiteEnergyShift',onsiteEnergyShift,0.150_dp)
   onsiteEnergyShift = onsiteEnergyShift/g0
   if (u) then
      call MIO_Print('Adding an onsite energy shift','ham')
      !$OMP PARALLEL DO PRIVATE(i)
      do i=1,nAt
         H0(i) = H0(i) + onsiteEnergyShift
      end do
      !$OMP END PARALLEL DO
   end if




   call MIO_InputParameter('Bubbles',l,.false.)
   if (l) then
      call MIO_InputParameter('onsiteShift',onsiteShift,0.1_dp)
      call MIO_InputParameter('checkerDivider',checkerDivider,10)
      call MIO_InputParameter('checkerDensity',checkerDensity,0.1_dp)
      allocate(checkerActivate(checkerDivider,checkerDivider))
      activatedCheckers = 0
      call random_seed(size = n)
      allocate(seed(n))
      call system_clock(COUNT=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(PUT = seed)
      deallocate(seed)
   20 do j = 1,checkerDivider
         do k = 1,checkerDivider
            checkerActivate(j,k) = .false.
            call random_number(rand)
            !print*, rand, checkerDensity
            if(rand.le.checkerDensity) then
                checkerActivate(j,k) = .true.
                activatedCheckers = activatedCheckers + 1
                if (dble(activatedCheckers)/dble(checkerDivider*checkerDivider)+0.01_dp.ge.checkerDensity) then
                   GO TO 10
                end if
            end if
         end do
      end do
      IF (dble(activatedCheckers)/dble(checkerDivider*checkerDivider)-0.01_dp.lt.checkerDensity) GO TO 20
   10 CONTINUE
      !print*, "density of checkers = ", dble(activatedCheckers)/dble(checkerDivider*checkerDivider)
      if (.not. frac) call AtomsSetFrac()
      !if (frac) call AtomsSetCart()
      !$OMP PARALLEL DO PRIVATE(i)
      do i = in1,in2
         !if (mod(CEILING(Rat(1,i)*checkerDivider),2).ne.mod(CEILING(Rat(2,i)*checkerDivider),2)) then
         !      H0(i) = H0(i)
         !else
         !      H0(i) = H0(i) + onsiteShift
         !endif
         if (checkerActivate(CEILING(Rat(1,i)*checkerDivider),CEILING(Rat(2,i)*checkerDivider)).eqv.(.true.)) then
               H0(i) = H0(i) + onsiteShift
         else
               H0(i) = H0(i)
         end if
      end do
      !$OMP END PARALLEL DO
   end if


   ! PNP
   call MIO_InputParameter('PNP',l,.false.)
   if (l) then
      if (frac) call AtomsSetCart()
      call MIO_InputParameter('CellSize',n,50)
      call MIO_InputParameter('PNPAmp',PNPAmp,0.01_dp)
      !PNPLeft = (n*sCell*aG)/8.0 * 3.5_dp
      !PNPRight = (n*sCell*aG)/8.0 * 4.5_dp
      PNPLeft = (n*sCell*aG)/6.0_dp * 2.0_dp
      PNPRight = (n*sCell*aG)/6.0_dp * 4.0_dp
      !$OMP PARALLEL DO PRIVATE(i,H)
      do i = in1,in2
        if (Rat(1,i).lt.PNPLeft) then
           H0(i) = H0(i) - PNPAmp
        else if (Rat(1,i).ge.PNPLeft .and. Rat(1,i).le.PNPRight) then
           H0(i) = H0(i) + PNPAmp
        else if (Rat(1,i).ge.(PNPLeft+(n*Scell*aG)) .and. Rat(1,i).le.(PNPRight+(n*Scell*aG))) then
           H0(i) = H0(i) + PNPAmp
        else if (Rat(1,i).gt.PNPRight) then
           H0(i) = H0(i) - PNPAmp
        else
           call MIO_Kill('PNP Inconsistency','ham','HamOnSite')
        end if
      end do
      !$OMP END PARALLEL DO
   end if

   call MIO_InputParameter('PNPKink',l,.false.)
   if (l) then
      if (frac) call AtomsSetCart()
      call MIO_InputParameter('CellSize',n,50)
      call MIO_InputParameter('PNPAmp',PNPAmp,0.01_dp)
      call MIO_InputParameter('PNPDelta',delta,10.0_dp)
      limit0 = 0.0_dp
      limit1 = (n*sCell*aG)*1.0_dp/4.0_dp
      limit2 = (n*sCell*aG)*2.0_dp/4.0_dp
      limit3 = (n*sCell*aG)*3.0_dp/4.0_dp
      limit4 = (n*sCell*aG)
      !$OMP PARALLEL DO PRIVATE(i,H)
      do i = in1,in2
        if (Rat(1,i).lt.limit2) then
          H = PNPAmp * tanh((Rat(1,i)-limit1)/delta)
        else if ((Rat(1,i).gt.limit2).and.(Rat(1,i).lt.limit4)) then
          H = - PNPAmp * tanh((Rat(1,i)-limit3)/delta)
        else if ((Rat(1,i).gt.limit4)) then
          H = PNPAmp * tanh((Rat(1,i)-(limit4+limit1))/delta)
        end if
        H0(i) = H0(i) + H
        !call MIO_Kill('PNP Inconsistency','ham','HamOnSite')
      end do
      !$OMP END PARALLEL DO
   end if


   ! changes for sinus function
   call MIO_InputParameter('sinusModulation',l,.false.)
   if (l) then
      call MIO_InputParameter('sinusNumberOfPeriod',sinusNumberOfPeriod,1)
      call MIO_InputParameter('CellSize',n,50)
      if (frac) call AtomsSetCart()
      !$OMP PARALLEL DO PRIVATE(i)
      do i = in1,in2
        H0(i) = H0(i) + sin(sinusNumberOfPeriod*2.0_dp*pi*Rat(1,i)/(n*sCell*aG))
      end do
      !$OMP END PARALLEL DO
      !print*, sCell
   end if
   ! end changes for sinus function

   ! changes for sinus function

   call MIO_InputParameter('sinusModulationUsingPeriod',sinusModulationUsingPeriod,.false.)
   if (sinusModulationUsingPeriod) then
      call MIO_InputParameter('sinusModulationPeriod',sinusModulationPeriod,135.0_dp)
      call MIO_InputParameter('sinusModulationAddMassTerm',sinusModulationAddMassTerm,.false.)
      call MIO_InputParameter('sinusFactor',sinusFactor,0.01_dp)
      call MIO_InputParameter('sinusMassTermFactor',sinusFactor2,0.01_dp)
      if (frac) call AtomsSetCart()
      !$OMP PARALLEL DO PRIVATE(i)
      do i = 1,nAt
        if (sinusModulationAddMassTerm) then
            H0(i) = H0(i) + sinusFactor2*sin(2.0_dp*pi*Rat(1,i)/(sinusModulationPeriod))*(-1.0_dp)**Species(i)
        else
            H0(i) = H0(i) + sinusFactor*sin(2.0_dp*pi*Rat(1,i)/(sinusModulationPeriod))
        !if (sinusModulationAddMassTerm) then
        !    H0(i) = H0(i) + sinusFactor2*sin(2.0_dp*pi*Rat(1,i)/(sinusModulationPeriod))*(-1.0_dp)**Species(i)
        end if
      end do
      !$OMP END PARALLEL DO
      !print*, sCell
   end if
   ! end changes for sinus function

   ! changes for sinus function

   call MIO_InputParameter('sinusModulationUsingPeriodYDirection',sinusModulationUsingPeriodY,.false.)
   if (sinusModulationUsingPeriodY) then
      call MIO_InputParameter('sinusModulationPeriod',sinusModulationPeriod,135.0_dp)
      call MIO_InputParameter('sinusModulationAddMassTerm',sinusModulationAddMassTerm,.false.)
      call MIO_InputParameter('sinusFactor',sinusFactor,0.01_dp)
      call MIO_InputParameter('sinusMassTermFactor',sinusFactor2,0.01_dp)
      if (frac) call AtomsSetCart()
      !$OMP PARALLEL DO PRIVATE(i)
      do i = 1,nAt
        if (sinusModulationAddMassTerm) then
            H0(i) = H0(i) + sinusFactor2*sin(2.0_dp*pi*Rat(2,i)/(sinusModulationPeriod))*(-1.0_dp)**Species(i)
        else
          H0(i) = H0(i) + sinusFactor*sin(2.0_dp*pi*Rat(2,i)/(sinusModulationPeriod))
        end if
      end do
      !$OMP END PARALLEL DO
      !print*, sCell
   end if
   ! end changes for sinus function

   ! changes for cosinus function

   call MIO_InputParameter('cosinusModulationUsingPeriod',cosinusModulationUsingPeriod,.false.)
   if (cosinusModulationUsingPeriod) then
      call MIO_InputParameter('cosinusModulationPeriod',cosinusModulationPeriod,135.0_dp)
      call MIO_InputParameter('cosinusModulationAddMassTerm',cosinusModulationAddMassTerm,.false.)
      call MIO_InputParameter('cosinusFactor',cosinusFactor,0.01_dp)
      call MIO_InputParameter('cosinusMassTermFactor',cosinusFactor2,0.01_dp)
      if (frac) call AtomsSetCart()
      !$OMP PARALLEL DO PRIVATE(i)
      do i = 1,nAt
        if (cosinusModulationAddMassTerm) then
            H0(i) = H0(i) + cosinusFactor2*cos(2.0_dp*pi*Rat(1,i)/(cosinusModulationPeriod))*(-1.0_dp)**Species(i)
        else
          H0(i) = H0(i) + cosinusFactor*cos(2.0_dp*pi*Rat(1,i)/(cosinusModulationPeriod))
        end if
      end do
      !$OMP END PARALLEL DO
      !print*, sCell
   end if
   ! end changes for sinus function

   ! changes for cosinus function

   call MIO_InputParameter('cosinusModulationUsingPeriodYDirection',cosinusModulationUsingPeriodY,.false.)
   if (cosinusModulationUsingPeriodY) then
      call MIO_InputParameter('cosinusModulationPeriod',cosinusModulationPeriod,135.0_dp)
      call MIO_InputParameter('cosinusModulationAddMassTerm',cosinusModulationAddMassTerm,.false.)
      call MIO_InputParameter('cosinusFactor',cosinusFactor,0.01_dp)
      call MIO_InputParameter('cosinusMassTermFactor',cosinusFactor2,0.01_dp)
      if (frac) call AtomsSetCart()
      !$OMP PARALLEL DO PRIVATE(i)
      do i = 1,nAt
        if (cosinusModulationAddMassTerm) then
            H0(i) = H0(i) + cosinusFactor2*cos(2.0_dp*pi*Rat(2,i)/(cosinusModulationPeriod))*(-1.0_dp)**Species(i)
        else
          H0(i) = H0(i) + cosinusFactor*cos(2.0_dp*pi*Rat(2,i)/(cosinusModulationPeriod))
        end if
      end do
      !$OMP END PARALLEL DO
      !print*, sCell
   end if
   ! end changes for sinus function

   ! square function
   call MIO_InputParameter('SquareFunction',l,.false.)
   if (l) then
      call MIO_InputParameter('sinusNumberOfPeriod',sinusNumberOfPeriod,1)
      call MIO_InputParameter('CellSize',n,50)
      call MIO_InputParameter('AmplitudeOfSquare',A,1.0_dp)
      if (frac) call AtomsSetCart()

      call MIO_InputParameter('TwoDimensional',l,.false.)
      if (l) then

        call MIO_InputParameter('AddZTerm',l,.false.)
           if (l) then

      !$OMP PARALLEL DO PRIVATE(i,H)
      do i = in1,in2
        H = sin(sinusNumberOfPeriod*2.0_dp*pi*Rat(1,i)/(n*sCell*aG))*sin(sinusNumberOfPeriod*2.0_dp*pi*Rat(2,i)/(n*sCell*aG))
         if (H.GT.0) then
           H0(i) = H0(i) - (-1.0_dp)**Species(i)*A
         else
           H0(i) = H0(i) + (-1.0_dp)**Species(i)*A
         end if
      end do
      !$OMP END PARALLEL DO

      else

      !$OMP PARALLEL DO PRIVATE(i,H)
      do i = in1,in2
        H = sin(sinusNumberOfPeriod*2.0_dp*pi*Rat(1,i)/(n*sCell*aG))*sin(sinusNumberOfPeriod*2.0_dp*pi*Rat(2,i)/(n*sCell*aG))
         if (H.GT.0) then
           H0(i) = H0(i) + A
         else
           H0(i) = H0(i)
         end if
      end do
      !$OMP END PARALLEL DO
        end if

      else
      !$OMP PARALLEL DO PRIVATE(i,H)
      do i = in1,in2
        H = sin(sinusNumberOfPeriod*2.0_dp*pi*Rat(1,i)/(n*sCell*aG))
         if (H.GT.0) then
           H0(i) = H0(i) + A
         else
           H0(i) = H0(i)
         end if
      end do
      !$OMP END PARALLEL DO
      end if
   end if
   ! square function

   ! square function 2
   call MIO_InputParameter('SquareFunction2',l,.false.)
   if (l) then
      call MIO_InputParameter('NumberOfWidthHoneycomb',NOWH,1)
      call MIO_InputParameter('CellSize',n,50)
      call MIO_InputParameter('AmplitudeOfSquare2',Amp2,0.01_dp)
      call MIO_InputParameter('PhaseOfSquareX',P1,0.0_dp)
      call MIO_InputParameter('PhaseOfSquareY',P2,0.0_dp)
      if (frac) call AtomsSetCart()

      call MIO_InputParameter('TwoDimension',l,.false.)
      if (l) then
        call MIO_InputParameter('AddZTerm',l,.false.)
           if (l) then

      !$OMP PARALLEL DO PRIVATE(i,H)
      do i = in1,in2
        H = sin(2.0_dp*pi*Rat(1,i)/(2.0_dp*aG*NOWH)+P1)*sin(2.0_dp*pi*Rat(2,i)*sqrt(3.0_dp)/(4.0_dp*aG*NOWH)+P2)
         if (H.GT.0) then
           H0(i) = H0(i) - (-1.0_dp)**Species(i)*Amp2
         else
           H0(i) = H0(i) + (-1.0_dp)**Species(i)*Amp2
         end if
      end do
      !$OMP END PARALLEL DO

      else
      !$OMP PARALLEL DO PRIVATE(i,H)
      do i = in1,in2
        H = sin(2.0_dp*pi*Rat(1,i)/(2.0_dp*aG*NOWH)+P1)*sin(2.0_dp*pi*Rat(2,i)*sqrt(3.0_dp)/(4.0_dp*aG*NOWH)+P2)
         if (H.GT.0) then
           H0(i) = H0(i) + Amp2
         else
           H0(i) = H0(i)
         end if
      end do
      !$OMP END PARALLEL DO
        end if

      else
      call MIO_InputParameter('ArmChairShape',l,.false.)
      if (l) then
      !!$OMP& SHARED (Rat, Amp2, aG, NOWH, P1, H0), &
      !!$OMP& DEFAULT(NONE)
      !$OMP PARALLEL DO PRIVATE(i,H)
      do i = in1,in2
        H = sin(2.0_dp*pi*Rat(1,i)/(2.0_dp*aG*NOWH)+P1)
         if (H.GT.0) then
           H0(i) = H0(i) + Amp2
         else
           H0(i) = H0(i)
         end if
      end do
      !$OMP END PARALLEL DO

      else
      !$OMP PARALLEL DO PRIVATE(i,H)
      do i = in1,in2
        H = sin(2.0_dp*pi*Rat(2,i)*sqrt(3.0_dp)/(4.0_dp*aG*NOWH)+P2)
         if (H.GT.0) then
           H0(i) = H0(i) + Amp2
         else
           H0(i) = H0(i)
         end if
      end do
      !$OMP END PARALLEL DO
      end if
      end if
   end if
   ! square function 2

   ! zterm1d
   call MIO_InputParameter('Zterm1D',l,.false.)
   if (l) then
      if (frac) call AtomsSetCart()
      call MIO_InputParameter('sinusNumberOfPeriod',sinusNumberOfPeriod,1)
      call MIO_InputParameter('CellSize',n,50)
      call MIO_InputParameter('AmplitudeOfSquare3',Amp3,0.01_dp)
      !print*, sinusNumberOfPeriod, n, Amp3, aG, sCell

      !$OMP PARALLEL DO PRIVATE(i,H)
      do i = in1,in2
        H = sin(sinusNumberOfPeriod*2.0_dp*pi*Rat(1,i)/(n*sCell*aG))
         if (H.GT.0) then
           H0(i) = H0(i) - (-1.0_dp)**Species(i)*Amp3
         else
           H0(i) = H0(i) + (-1.0_dp)**Species(i)*Amp3
         end if
      end do
      !$OMP END PARALLEL DO
   end if

   call MIO_InputParameter('Zterm1DKink',l,.false.)   ! only works for 2 periods
   if (l) then
      if (frac) call AtomsSetCart()
      call MIO_InputParameter('CellSize',n,50)
      call MIO_InputParameter('Zterm1DAmp',Amp3,0.01_dp)
      call MIO_InputParameter('Zterm1DDelta',delta,10.0_dp)
      limit0 = 0.0_dp
      limit1 = (n*sCell*aG)*1.0_dp/4.0_dp
      limit2 = (n*sCell*aG)*2.0_dp/4.0_dp
      limit3 = (n*sCell*aG)*3.0_dp/4.0_dp
      limit4 = (n*sCell*aG)
      !$OMP PARALLEL DO PRIVATE(i,H)
      do i = in1,in2
        if (Rat(1,i).lt.limit2) then
          H = Amp3 * tanh((Rat(1,i)-limit1)/delta)
        else if ((Rat(1,i).gt.limit2).and.(Rat(1,i).lt.limit4)) then
          H = - Amp3 * tanh((Rat(1,i)-limit3)/delta)
        else if ((Rat(1,i).gt.limit4)) then
          H = Amp3 * tanh((Rat(1,i)-(limit4+limit1))/delta)
        end if
        H0(i) = H0(i) + (-1.0_dp)**Species(i)*H
      end do
      !$OMP END PARALLEL DO
   end if

   ! zterm1d
   ! squarechecker
 call MIO_InputParameter('SquareChecker2219',l,.false.)
      if (l) then

      call MIO_InputParameter('CellSize',n,50)
      call MIO_InputParameter('AmplitudeOfSquare',A,1.0_dp)
      if (frac) call AtomsSetCart()


      !$OMP PARALLEL DO PRIVATE(i,H)

      do i = in1,in2

        H = sin(22*2.0_dp*pi*Rat(1,i)/(n*sCell*aG))*sin(19*4.0_dp*pi*Rat(2,i)/(sqrt(3.0_dp)*n*sCell*aG))

         if (H.GT.0) then

           H0(i) = H0(i) - (-1.0_dp)**Species(i)*A

         else

           H0(i) = H0(i) + (-1.0_dp)**Species(i)*A

         end if

      end do

      !$OMP END PARALLEL DO

end if

   ! squarechecker
   call MIO_InputParameter('SublatticeDisorder',l,.false.)
   if (l) then
      call random_seed(size = n)
      allocate(seed(n))
      call system_clock(COUNT=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(PUT = seed)
      deallocate(seed)
      call MIO_InputParameter('SublattAmp',A,2.0_dp)
      call MIO_InputParameter('SublattPct',pct,0.1_dp)
      call MIO_Print('Sublattice disorder','ham')
      call MIO_Print('  w: '//trim(num2str(A,4)),'ham')
      call MIO_Print('')
      !!$OMP PARALLEL DO
      !call RandSeed(rng,nThread)
      !do i=1,nAt
      do i=inode1,inode2
         call random_number(rand)
         if(rand.le.(pct*2.0_dp) .and. Species(i).eq.2) then  ! factor 2 to compensate for the sublattice restriction
             H0(i) = H0(i) + A
         end if
      end do
      !!$OMP END PARALLEL DO
   end if

   call MIO_InputParameter('Anderson',l,.false.)
   if (l) then
      call random_seed(size = n)
      allocate(seed(n))
      call system_clock(COUNT=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(PUT = seed)
      deallocate(seed)
      call MIO_InputParameter('AndersonAmp',A,1.0_dp)
      call MIO_Print('Anderson disorder','ham')
      call MIO_Print('  w: '//trim(num2str(A,4)),'ham')
      call MIO_Print('')
      !!$OMP PARALLEL DO
      !call RandSeed(rng,nThread)
      !do i=1,nAt
      do i=inode1,inode2
         call random_number(rand)
         H0(i) = H0(i) + (rand-0.5_dp)*A
      end do
      !!$OMP END PARALLEL DO
   end if
   call MIO_InputParameter('deltaDisorder',l,.false.)
   if (l) then
      call random_seed(size = n)
      allocate(seed(n))
      call system_clock(COUNT=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(PUT = seed)
      deallocate(seed)
      call MIO_InputParameter('deltaAmp',A,1.0_dp)
      call MIO_Print('delta disorder','ham')
      call MIO_Print('  w: '//trim(num2str(A,4)),'ham')
      call MIO_Print('')
      call MIO_InputParameter('deltaSkewFactor',AA,1.0_dp)
      call MIO_Print('delta disorder','ham')
      call MIO_Print('  skewFactor: '//trim(num2str(AA,4)),'ham')
      call MIO_Print('')
      tot = 0.0_dp
      tot2 = 0.0_dp
      tot3 = 0.0_dp
      !!$OMP PARALLEL DO
      !call RandSeed(rng,nThread)
      !do i=1,nAt
      do i=inode1,inode2
         call random_number(rand)
         call randomInRange(rand,rand2,-1.0_dp,1.0_dp/AA)
         if (rand2<0) then
             call random_number(rand)
             call randomInRange(rand,rand2,-1.0_dp/AA,0.0_dp)
         else
             call random_number(rand)
             call randomInRange(rand,rand2,0.0_dp,1.0_dp)
         end if
         randomPot = rand2
         randomPot2 = rand2**2.0_dp
         randomPot3 = rand2**3.0_dp
         tot = tot + randomPot
         tot2 = tot2 + randomPot2
         tot3 = tot3 + randomPot3
         H0(i) = H0(i) + randomPot*A
      end do
      !!$OMP END PARALLEL DO
      call MIO_Print('delta disorder','ham')
      nImp = inode2
      call MIO_Print('  mean: '//trim(num2str((tot/nImp),9)),'ham')
      call MIO_Print('  variance: '//trim(num2str((tot2/nImp),9)),'ham')
      call MIO_Print('  skewness: '//trim(num2str((tot3/nImp),9)),'ham')
   end if
   call MIO_InputParameter('GaussDisorder',l,.false.)
   if (l) then
      call GaussPot(H0)
   end if
   call MIO_InputParameter('TypeOfSystem',str,'Graphene')
   !if (MIO_StringComp(str,'Ribbons') .or. MIO_StringComp(str,'Hybrid') .or. MIO_StringComp(str,'ReadXYZ')) then
   !   call InterfacePot1()
   !end if
   if (Zterm .or. PZterm) then
      call MIO_Allocate(Ho,[inode1],[inode2],'Ho','ham')
      Ho = H0
   end if

   call MIO_InputParameter('fourLayerOnsiteShifts',l,.false.)
   call MIO_InputParameter('fourLayerShift1',layerShift1,0.0_dp)
   call MIO_InputParameter('fourLayerShift2',layerShift2,0.0_dp)
   call MIO_InputParameter('fourLayerShift3',layerShift3,0.0_dp)
   call MIO_InputParameter('fourLayerShift4',layerShift4,0.0_dp)
   layerShift1 = layerShift1/g0
   layerShift2 = layerShift2/g0
   layerShift3 = layerShift3/g0
   layerShift4 = layerShift4/g0
   if (l) then
      !$OMP PARALLEL DO PRIVATE(i)
      do i = 1,nAt
           if (layerIndex(i).eq.1) then
              H0(i) = H0(i) + layershift1
           else if (layerIndex(i).eq.2) then
              H0(i) = H0(i) + layershift2
           else if (layerIndex(i).eq.3) then
              H0(i) = H0(i) + layershift3
           else if (layerIndex(i).eq.4) then
              H0(i) = H0(i) + layershift4
           end if
      end do
      !$OMP END PARALLEL DO
   end if

   ! CDW modulation for helical twisted MBM system (applied at the end, after all other onsite modifications)
   if (helicalTwistedMBM) then
      call MIO_InputParameter('helicalTwistedMBM_CDW',helicalTwistedMBM_CDW,.false.)
      if (helicalTwistedMBM_CDW) then
         call MIO_Print('Adding CDW modulation for helical twisted MBM system','ham')
         call MIO_InputParameter('CDWAmplitude',CDWAmplitude,0.01_dp)
         CDWAmplitude = CDWAmplitude/g0
         call MIO_InputParameter('CDWPeriod',CDWPeriod,1.0_dp)
         call MIO_InputParameter('CDWUseMassTerm',l,.false.)
         if (.not. frac) call AtomsSetFrac()
         call MIO_Print('  Using fractional coordinates for CDW modulation (Rat in fractional units)','ham')
         call MIO_Print('  CDW period: '//trim(num2str(CDWPeriod,5))//' (in fractional units, 1.0 = full simulation cell)','ham')
         call MIO_Print('  CDW amplitude: '//trim(num2str(CDWAmplitude*g0,5))//' eV','ham')
         if (l) then
            call MIO_Print('  CDW using sublattice-dependent mass term (opens spatially varying gap)','ham')
            !$OMP PARALLEL DO PRIVATE(i)
            do i=1,nAt
               !H0(i) = H0(i) + (-1.0_dp)**Species(i) * CDWAmplitude * cos(twopi*Rat(1,i)*CDWPeriod) * cos(twopi*Rat(2,i)*CDWPeriod)
               H0(i) = H0(i) + (-1.0_dp)**Species(i) * (1.0_dp/3.0_dp) * CDWAmplitude * &
        ( cos(twopi*sCell*Rat(1,i)) + &
          cos(twopi*sCell*Rat(2,i)) + &
          cos(twopi*sCell*(Rat(1,i)+Rat(2,i))) )
            end do
            !$OMP END PARALLEL DO
         else
            call MIO_Print('  CDW using uniform energy shift (no gap opening)','ham')
            !$OMP PARALLEL DO PRIVATE(i)
            do i=1,nAt
               !H0(i) = H0(i) + CDWAmplitude * cos(twopi*Rat(1,i)*CDWPeriod) * cos(twopi*Rat(2,i)*CDWPeriod)
               H0(i) = H0(i) + (1.0_dp/3.0_dp) * CDWAmplitude * &
        ( cos(twopi*sCell*Rat(1,i)) + &
          cos(twopi*sCell*Rat(2,i)) + &
          cos(twopi*sCell*(Rat(1,i)+Rat(2,i))) )
            end do
            !$OMP END PARALLEL DO
         end if
         ! Debug: print first few H0 values to verify modulation
         call MIO_Print('  CDW: Sample H0 values after modulation (first 5 atoms):','ham')
         do i=1,min(5,nAt)
            call MIO_Print('    Atom '//trim(num2str(real(i,dp),0))//': H0='//trim(num2str(H0(i)*g0,8))//' eV, x='//trim(num2str(Rat(1,i),5))//' y='//trim(num2str(Rat(2,i),5)),'ham')
         end do
      end if
   end if

#ifdef DEBUG
   call MIO_Debug('HamOnSite',1)
#endif /* DEBUG */

end subroutine HamOnSite

!> @brief Rescale a random value to be inside a given range.
!! @param[in]  rand   Input random value in [0,1]
!! @param[out] rand2  Rescaled random value in [low,high]
!! @param[in]  low    Lower bound of output range
!! @param[in]  high   Upper bound of output range
subroutine randomInRange(rand,rand2, low,high)

    !real(dp), intent(in) :: low, high, rand
    !real(dp), intent(out) :: rand2
    real(dp), intent(in) :: low, high, rand
    real(dp), intent(out) :: rand2

    rand2 = (rand*(high-low))+low

end subroutine randomInRange

!> @brief Fit function for distance-dependent coupling parameters.
!! @param[in]  p1,p2,p3  Fitting coefficients
!! @param[in]  dist       Interlayer distance
!! @param[out] Cjj        Fitted coupling parameter
subroutine fitFunSrivani(p1, p2, p3, dist, Cjj)

    real(dp) :: p1 ! parabolic function coefficient
    real(dp) :: p2 ! parabolic function coefficient
    real(dp) :: p3 ! parabolic function coefficient
    real(dp) :: dist ! distance separating the two layers
    real(dp) :: Cjj ! fitted value

    Cjj = p1*abs(dist)**2 + p2*abs(dist) + p3
    return
end subroutine fitFunSrivani


!> @brief Compute distance-dependent coupling using exponential decay.
!! @param[out] Cd      Distance-dependent coupling
!! @param[in]  C0      Base coupling at reference distance
!! @param[in]  Bfactor Exponential decay factor
!! @param[in]  z       Current interlayer distance
!! @param[in]  z0      Reference interlayer distance
subroutine distanceDependentC(Cd, C0, Bfactor, z, z0)

    real(dp), intent(in) :: C0, z, z0, Bfactor
    real(dp), intent(out) :: Cd
    logical :: writeData
    !call MIO_InputParameter('WriteDataFiles',writeData,.false.)

    Cd = C0 * exp(-Bfactor * (z-z0))
    !if (writeData .and. Bfactor.eq.3.2_dp) then
    !    write(587,*) Cd
    !end if
    return
end subroutine distanceDependentC

!> @brief Angular / orbital-orientation factor of the interlayer p_z-p_z
!!        two-center hopping. Single implementation, shared by every interlayer
!!        Vpppi/Vppsigma site so the expression is never duplicated.
!! @param[in] vpi         Vpppi(d)    from the existing radial model (unchanged)
!! @param[in] vsig        Vppsigma(d) from the existing radial model (unchanged)
!! @param[in] rvec        bond vector NeighD(:,j,i)  (points atom i -> atom j)
!! @param[in] dist        distance used by the existing radial model (legacy path only)
!! @param[in] ni, nj      local unit p_z axes of the two atoms (corrugated path only)
!! @param[in] corrugated  .false. -> legacy global-z expression, same operands as before
!!                        .true.  -> general Slater-Koster with local axes ni, nj
!!
!! corrugated = .false. :  vpi*(1 - (rvec(3)/dist)**2) + vsig*(rvec(3)/dist)**2
!! corrugated = .true.  :  vpi*dot(ni,nj) + (vsig-vpi)*dot(ni,dhat)*dot(nj,dhat)
!!                         with dhat = rvec/|rvec|.  The angular direction comes
!!                         purely from the geometric NeighD vector; the radial
!!                         model keeps its own dist / Vpppi / Vppsigma.
!! For ni = nj = z_hat the two branches are algebraically equal.
pure function interlayerTwoCenter(vpi, vsig, rvec, dist, ni, nj, corrugated) result(t)

   real(dp), intent(in) :: vpi, vsig, rvec(3), dist, ni(3), nj(3)
   logical,  intent(in) :: corrugated
   real(dp) :: t, dgeom, nid, njd

   if (corrugated) then
      dgeom = sqrt(rvec(1)**2.0_dp + rvec(2)**2.0_dp + rvec(3)**2.0_dp)
      nid   = dot_product(ni, rvec) / dgeom
      njd   = dot_product(nj, rvec) / dgeom
      t     = vpi * dot_product(ni, nj) + (vsig - vpi) * nid * njd
   else
      t     = vpi  * (1.0_dp - ((rvec(3))/(dist))**2.0_dp) &
            + vsig * ((rvec(3))/(dist))**2.0_dp
   end if

end function interlayerTwoCenter

!> @brief Build one local p_z orbital axis per locally-owned atom from its
!!        first-shell same-layer neighbour bond vectors. Used only when
!!        corrugatedInterlayerTwoCenter = .true.
!! @details Same layer  = equal layerIndex.
!!          First shell  = |NeighD| < Nradii(1,nlay)  (existing shell radius).
!!          Normal       = normalised sum of +z-oriented pairwise cross products
!!                         of the (unit) same-layer first-shell bonds.
!!          -> exactly (0,0,1) for a flat xy layer, smooth for a corrugated one,
!!             independent of neighbour ordering. Atoms with fewer than two
!!             usable neighbours fall back to global z_hat.
subroutine ComputePzNormals(pzn, nfallback)

   use neigh,  only : NeighD, NList, Nneigh, Nradii
   use atoms,  only : Species, layerIndex, inode1, inode2
   use math,   only : CrossProd, norm

   real(dp), intent(inout) :: pzn(:, inode1:)   ! filled for the owned range inode1:inode2
   integer,  intent(out)   :: nfallback         ! # owned atoms that fell back to z_hat

   integer  :: i, ja, jb, nlay
   real(dp) :: acc(3), c(3), ba(3), bb(3), la, lb, rcut, nc, na

   nfallback = 0
   do i = inode1, inode2
      nlay = (Species(i) - 1) / 2 + 1
      rcut = Nradii(1, nlay)
      acc  = 0.0_dp
      do ja = 1, Nneigh(i)
         if (layerIndex(NList(ja,i)) /= layerIndex(i)) cycle
         la = norm(NeighD(:,ja,i))
         if (la >= rcut .or. la < 1.0e-6_dp) cycle
         ba = NeighD(:,ja,i) / la
         do jb = ja + 1, Nneigh(i)
            if (layerIndex(NList(jb,i)) /= layerIndex(i)) cycle
            lb = norm(NeighD(:,jb,i))
            if (lb >= rcut .or. lb < 1.0e-6_dp) cycle
            bb = NeighD(:,jb,i) / lb
            c  = CrossProd(ba, bb)
            nc = norm(c)
            if (nc < 1.0e-10_dp) cycle
            if (c(3) < 0.0_dp) c = -c
            acc = acc + c / nc
         end do
      end do
      na = norm(acc)
      if (na > 1.0e-10_dp) then
         pzn(:,i) = acc / na
      else
         pzn(:,i) = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
         nfallback = nfallback + 1
      end if
   end do

end subroutine ComputePzNormals

#ifdef DEBUG
!> @brief DEBUG-only concise summary of the pzNormal field (per layer): normal
!!        component ranges, |n| range, max/RMS tilt vs +z, max ||n|-1|, the
!!        z_hat-fallback count and the most tilted atom. One block, no per-atom
!!        dump. Compiled only with -DDEBUG.
subroutine HamPrintNormalStats(pzn, nfallback)

   use atoms,      only : layerIndex, inode1, inode2
   use math,       only : norm
   use constants,  only : pi

   real(dp), intent(in) :: pzn(:, inode1:)
   integer,  intent(in) :: nfallback

   real(dp), parameter :: r2d = 180.0_dp / pi
   integer  :: i, lay, laymin, laymax, ntot, imaxtilt
   real(dp) :: nx0, nx1, ny0, ny1, nz0, nz1, m0, m1, mag, tilt, tmax, tsum, maxmagdev

   laymin = minval(layerIndex(inode1:inode2))
   laymax = maxval(layerIndex(inode1:inode2))
   maxmagdev = 0.0_dp; tmax = 0.0_dp; imaxtilt = inode1
   do i = inode1, inode2
      mag = norm(pzn(:,i))
      maxmagdev = max(maxmagdev, abs(mag - 1.0_dp))
      tilt = acos(max(-1.0_dp, min(1.0_dp, pzn(3,i))))
      if (tilt > tmax) then
         tmax = tilt; imaxtilt = i
      end if
   end do

   call MIO_Print('--- pzNormal field summary (corrugatedInterlayerTwoCenter) ---','ham')
   do lay = laymin, laymax
      ntot = count(layerIndex(inode1:inode2) == lay)
      if (ntot == 0) cycle
      nx0 =  huge(1.0_dp); nx1 = -huge(1.0_dp); ny0 = huge(1.0_dp); ny1 = -huge(1.0_dp)
      nz0 =  huge(1.0_dp); nz1 = -huge(1.0_dp); m0  = huge(1.0_dp); m1  = -huge(1.0_dp)
      tsum = 0.0_dp
      do i = inode1, inode2
         if (layerIndex(i) /= lay) cycle
         nx0 = min(nx0, pzn(1,i)); nx1 = max(nx1, pzn(1,i))
         ny0 = min(ny0, pzn(2,i)); ny1 = max(ny1, pzn(2,i))
         nz0 = min(nz0, pzn(3,i)); nz1 = max(nz1, pzn(3,i))
         mag = norm(pzn(:,i)); m0 = min(m0, mag); m1 = max(m1, mag)
         tilt = acos(max(-1.0_dp, min(1.0_dp, pzn(3,i)))); tsum = tsum + tilt*tilt
      end do
      call MIO_Print('  layer '//trim(num2str(lay))//'  atoms '//trim(num2str(ntot)),'ham')
      call MIO_Print('    nx ['//trim(num2str(nx0,5))//','//trim(num2str(nx1,5))//']   ny ['// &
                     trim(num2str(ny0,5))//','//trim(num2str(ny1,5))//']','ham')
      call MIO_Print('    nz ['//trim(num2str(nz0,7))//','//trim(num2str(nz1,7))//']   |n| ['// &
                     trim(num2str(m0,10))//','//trim(num2str(m1,10))//']','ham')
      call MIO_Print('    tilt acos(nz) deg:  max '//trim(num2str(acos(max(-1.0_dp,min(1.0_dp,nz0)))*r2d,4))// &
                     '   RMS '//trim(num2str(sqrt(tsum/real(ntot,dp))*r2d,4)),'ham')
   end do
   call MIO_Print('  global: z_hat fallbacks = '//trim(num2str(nfallback))// &
                  '   most-tilted atom = '//trim(num2str(imaxtilt))// &
                  '   max||n|-1| = '//trim(num2str(maxmagdev*1.0e15_dp,4))//' e-15','ham')
   call MIO_Print('  global: max tilt = '//trim(num2str(tmax*r2d,4))//' deg','ham')

end subroutine HamPrintNormalStats
#endif /* DEBUG */

!> @brief Fill the MPI halo region of pzn for remote neighbour atoms, reusing
!!        the neigh send/receive lists (same exchange pattern as
!!        kubosubs::KuboUpdate_d). No-op for serial runs / nProc < 2.
!! @details After this call, pzn(:, NList2(j,i)) is valid for every neighbour j
!!          of every locally-owned atom i.
subroutine HamUpdateNormals(pzn)

   use atoms,  only : inode1, inode2
#ifdef MPI
   use neigh,  only : rcvIndx, sndList, sndIndx, sndSz, nodeSndRcv, rcvList
#endif /* MPI */

   real(dp), intent(inout) :: pzn(:, inode1:)

#ifdef MPI
   integer :: i, in, ip
   real(dp), allocatable :: sbuf(:), rbuf(:)

   if (nProc < 2) return

   ! Pack/unpack through flat rank-1 buffers (3 reals per atom) so the send and
   ! receive slices are contiguous rank-1 sections, matching MPISendRecv's
   ! assumed-size dummies (same exchange pattern as kubosubs::KuboUpdate_d).
   allocate(sbuf(3*sndSz), rbuf(3*size(rcvList)))
   do i = 1, sndSz
      sbuf(3*i-2:3*i) = pzn(:, sndList(i))
   end do
   do i = 1, nProc - 1
      in = mod(Node + i, nProc)
      ip = mod(nProc + Node - i, nProc)
      call MPISendRecv(sbuf(3*(sndIndx(1,in+1)-1)+1:), 3*nodeSndRcv(Node+1,in+1), in, 62, &
                       rbuf(3*(rcvIndx(1,ip+1)-1)+1:), 3*nodeSndRcv(ip+1,Node+1), ip, 62, &
                       MPI_DOUBLE_PRECISION)
      call MPIBarrier()
   end do
   do i = 1, size(rcvList)
      pzn(:, inode2+i) = rbuf(3*i-2:3*i)
   end do
   deallocate(sbuf, rbuf)
#endif /* MPI */

end subroutine HamUpdateNormals

!> @brief Compute H0 onsite term using AABB sublattice basis.
!! @param[out] H0             Onsite energy contribution
!! @param[in]  dx,dy          In-plane displacements
!! @param[in]  z              Interlayer distance
!! @param[in]  moirePreFactor Moiré pattern prefactor
subroutine diagoH0FromAABB(H0,dx,dy, z, moirePreFactor)

    use constants
    use cell,                 only : aG
    use tbpar,                only : g0
    real(dp), intent(in):: z, moirePreFactor
    real(dp):: dx, dy, G1, CAA0, CBB0, phiAA, phiBB
    real(dp) :: HdAA, HdBB, CAAd, CBBd
    real(dp), intent(out) :: H0

    !% gbn
    !Caa0 = 0.0021326; Cbb0 = 0.0060425; phi0 = 180*pi/180;
    !Caa =  0.0154451  ;    phiaa = -126.7332*pi/180 ;
    !Cbb =  0.0141258  ;    phibb =  -48.6533*pi/180 ;
    !Cab =  0.0120337  ;    phiab1 = 99.3551*pi/180;
    !phiab2 = 20.6449*pi/180; phiab3 = 140.6449*pi/180;

    !Choose the desired parameter set
    !CAA0 = 0.0021326_dp/g0
    !CBB0 = 0.0060425_dp/g0
    !phiAA = -126.7332_dp/180.0_dp*pi
    !phiBB = -48.6533_dp/180.0_dp*pi
    CAA0 = -0.01488_dp/g0*moirePreFactor
    CBB0 = 0.01209_dp/g0*moirePreFactor
    phiAA = 50.19_dp/180.0_dp*pi
    phiBB = -46.64_dp/180.0_dp*pi
    !CAA0 =  0.0154451_dp/g0
    !CBB0 =  0.0141258_dp/g0
    !phiAA = -126.7332*pi/180_dp
    !phiBB =  -48.6533*pi/180_dp

    call distanceDependentC(CAAd, CAA0, 3.0_dp, z, 3.35_dp)
    call distanceDependentC(CBBd, CBB0, 3.2_dp, z, 3.35_dp)
    !
    call diago(HdAA, dx, dy, CAAd, phiAA) ! distance-dependent
    call diago(HdBB, dx, dy, CBBd, phiBB)
    !call diago(HdAA, dx, dy, CAA0, phiAA) ! distance-independent
    !call diago(HdBB, dx, dy, CBB0, phiBB)


    H0 = (HdBB + HdAA)/2.0_dp

    return
end subroutine diagoH0FromAABB

!> @brief Compute Hz onsite term using AABB sublattice basis.
!! @param[out] Hz             Hz energy contribution
!! @param[in]  dx,dy          In-plane displacements
!! @param[in]  z              Interlayer distance
!! @param[in]  moirePreFactor Moiré pattern prefactor
subroutine diagoHzFromAABB(Hz,dx,dy,z,moirePreFactor)

    use constants
    use cell,                 only : aG
    use tbpar,                only : g0
    real(dp), intent(in):: z, moirePreFactor
    real(dp):: dx, dy, G1, CAA0, CBB0, phiAA, phiBB
    real(dp) :: HdAA, HdBB, CAAd, CBBd
    real(dp), intent(out) :: Hz

    !CAA0 = 0.0021326_dp/g0
    !CBB0 = 0.0060425_dp/g0
    !phiAA = -126.7332_dp/180.0_dp*pi
    !phiBB = -48.6533_dp/180.0_dp*pi
    CAA0 = -0.01488_dp/g0*moirePreFactor
    CBB0 = 0.01209_dp/g0*moirePreFactor
    phiAA = 50.19_dp/180.0_dp*pi
    phiBB = -46.64_dp/180.0_dp*pi
    !CAA0 =  0.0154451_dp/g0
    !CBB0 =  0.0141258_dp/g0
    !phiAA = -126.7332*pi/180_dp
    !phiBB =  -48.6533*pi/180_dp

    call distanceDependentC(CAAd, CAA0, 3.0_dp, z, 3.35_dp)
    call distanceDependentC(CBBd, CBB0, 3.2_dp, z, 3.35_dp)
    !!print*, z, CAAd, CAA0
    !!print*, z, CBBd, CBB0
    !
    call diago(HdAA, dx, dy, CAAd, phiAA)
    call diago(HdBB, dx, dy, CBBd, phiBB)
    !call diago(HdAA, dx, dy, CAA0, phiAA)
    !call diago(HdBB, dx, dy, CBB0, phiBB)

    Hz = (HdAA - HdBB)/2.0_dp

    return
end subroutine diagoHzFromAABB


!> @brief Compute diagonal onsite energy using cosine modulation.
!! @param[out] Hdjj   Diagonal onsite energy
!! @param[in]  dx,dy  In-plane displacements
!! @param[in]  Cjj    Coupling amplitude
!! @param[in]  Phijj  Phase factor
subroutine diago(Hdjj,dx,dy,Cjj,Phijj)

    use constants
    use cell,                 only : aG
    real(dp), intent(in):: Cjj, Phijj
    real(dp):: dx, dy, G1
    real(dp), intent(out) :: Hdjj

    !print*, dx, dy
    G1 = 4.0_dp*pi/sqrt(3.0_dp)/(aG)

    Hdjj = 2.0_dp*Cjj*real( (exp(-cmplx_i*G1*dy) &
               + 2.0_dp*exp(cmplx_i*G1*dy/2.0_dp)*cos(sqrt(3.0_dp)*G1*dx/2.0_dp))*exp(cmplx_i*Phijj) )

    return
end subroutine diago

!> @brief Compute diagonal onsite energy for GBN systems.
!! @param[out] Hdjj   Diagonal onsite energy
!! @param[in]  dx,dy  In-plane displacements
!! @param[in]  Cjj    Coupling amplitude
!! @param[in]  Phijj  Phase factor
!! @param[in]  Cjj0   Constant offset
subroutine diagoGBN(Hdjj,dx,dy,Cjj,Phijj,Cjj0)

    use constants
    use cell,                 only : aG
    real(dp), intent(in):: Cjj, Phijj
    real(dp), intent(in):: Cjj0
    real(dp):: dx, dy, G1
    real(dp), intent(out) :: Hdjj

    !print*, dx, dy
    G1 = 4.0_dp*pi/sqrt(3.0_dp)/(aG)

    Hdjj = Cjj0 + 2.0_dp*Cjj*real( (exp(-cmplx_i*G1*dy) &
               + 2.0_dp*exp(cmplx_i*G1*dy/2.0_dp)*cos(sqrt(3.0_dp)*G1*dx/2.0_dp))*exp(cmplx_i*Phijj) )

    return
end subroutine diagoGBN

!> @brief Compute diagonal onsite energy with conjugate phase.
!! @param[out] Hdjj   Diagonal onsite energy
!! @param[in]  dx,dy  In-plane displacements
!! @param[in]  Cjj    Coupling amplitude
!! @param[in]  Phijj  Phase factor (conjugated)
subroutine diagoConj(Hdjj,dx,dy,Cjj,Phijj)

    use constants
    use cell,                 only : aG
    real(dp), intent(in):: Cjj, Phijj
    real(dp):: dx, dy, G1
    real(dp), intent(out) :: Hdjj

    !print*, dx, dy
    G1 = 4.0_dp*pi/sqrt(3.0_dp)/(aG)

    Hdjj = 2.0_dp*Cjj*real( (exp(cmplx_i*G1*dy) &
               + 2.0_dp*exp(-cmplx_i*G1*dy/2.0_dp)*cos(sqrt(3.0_dp)*G1*dx/2.0_dp))*exp(-cmplx_i*Phijj) )

    return
end subroutine diagoConj

!> @brief Compute diagonal onsite energy with constant offset.
!! @param[out] Hdjj   Diagonal onsite energy
!! @param[in]  dx,dy  In-plane displacements
!! @param[in]  Cjj    Coupling amplitude
!! @param[in]  Phijj  Phase factor
!! @param[in]  Cjj0   Constant offset
subroutine diago2(Hdjj,dx,dy,Cjj,Phijj,Cjj0)

    use constants
    use cell,                 only : aG
    real(dp), intent(in):: Cjj, Phijj,Cjj0
    real(dp):: dx, dy, G1
    real(dp), intent(out) :: Hdjj

    !print*, dx, dy
    G1 = 4.0_dp*pi/sqrt(3.0_dp)/(aG)

    Hdjj = Cjj0 + 2.0_dp*Cjj*real( (exp(-cmplx_i*G1*dy) &
           + 2.0_dp*exp(cmplx_i*G1*dy/2.0_dp)*cos(sqrt(3.0_dp)*G1*dx/2.0_dp))*exp(cmplx_i*Phijj) )

    return
end subroutine diago2



!> @brief Compute off-diagonal hopping using cosine modulation.
!! @param[out] Hodjj  Off-diagonal hopping element
!! @param[in]  dx,dy  In-plane displacements
!! @param[in]  Cjj    Coupling amplitude
!! @param[in]  Phijj  Phase factor
subroutine offdiago(Hodjj,dx,dy,Cjj,Phijj)

   use constants
   use cell,                 only : aG
   real(dp), intent(in)::  dx, dy, Cjj, Phijj
   real(dp):: G1
   complex*16, intent(out) :: Hodjj

   G1 = 4.0_dp*pi/sqrt(3.0_dp)/(aG)

   Hodjj = 2.0_dp*Cjj*cos(sqrt(3.0_dp)*G1*dx/2.0_dp)*( cos(G1*dy/2.0_dp - Phijj) + sin(G1*dy/2.0_dp - Phijj - pi/6.0_dp)  ) &
 + 2.0_dp*Cjj*sin( G1*dy + Phijj - pi/6.0_dp )  + &
 cmplx_i*2.0_dp*Cjj*sin( sqrt(3.0_dp)*G1*dx/2.0_dp )*  &
 ( cos( G1*dy/2.0_dp - Phijj )  - sin( G1*dy/2.0_dp - Phijj - pi/6.0_dp )  )


   return
end subroutine offdiago

!> @brief Compute off-diagonal hopping with constant offset.
!! @param[out] Hodjj  Off-diagonal hopping element
!! @param[in]  dx,dy  In-plane displacements
!! @param[in]  Cjj    Coupling amplitude
!! @param[in]  Phijj  Phase factor
subroutine offdiago2(Hodjj,dx,dy,Cjj,Phijj)
   use constants
   use cell,                 only : aG
   real(dp), intent(in)::  dx, dy, Cjj, Phijj
   real(dp):: G1
   complex*16, intent(out) :: Hodjj

   G1 = 4.0_dp*pi/sqrt(3.0_dp)/(aG)

   Hodjj = 2.0_dp*Cjj*cos(sqrt(3.0_dp)*G1*dx/2.0_dp)* cos(G1*dy/2.0_dp - Phijj)  &
 - 2.0_dp*Cjj*cos( G1*dy + Phijj )  - &
 cmplx_i*2.0_dp*sqrt(3.0_dp)*Cjj*sin( sqrt(3.0_dp)*G1*dx/2.0_dp )*  &
  sin( G1*dy/2.0_dp - Phijj )


   return
end subroutine offdiago2

!> @brief Compute off-diagonal hopping for GBN systems.
!! @param[out] Hodjj  Off-diagonal hopping element
!! @param[in]  dx,dy  In-plane displacements
!! @param[in]  Cjj    Coupling amplitude
!! @param[in]  Phijj  Phase factor
subroutine offdiagoGBN(Hodjj,dx,dy,Cjj,Phijj)
   use constants
   use cell,                 only : aG
   real(dp), intent(in)::  dx, dy, Cjj, Phijj
   real(dp):: G1
   complex*16, intent(out) :: Hodjj

   !G1 = 4.0_dp*pi/sqrt(3.0_dp)/(aG)
   G1 = 4.0_dp*pi/(sqrt(3.0_dp))

   Hodjj = 2.0_dp*Cjj*cos(sqrt(3.0_dp)*G1*dx/2.0_dp)*( cos(G1*dy/2.0_dp - Phijj) ) &
 - 2.0_dp*Cjj*cos( G1*dy + Phijj )  - &
 cmplx_i*2.0_dp*sqrt(3.0_dp)*Cjj*sin( sqrt(3.0_dp)*G1*dx/2.0_dp )*  &
 ( sin( G1*dy/2.0_dp - Phijj)  )


   return
end subroutine offdiagoGBN

!> @brief Compute off-diagonal hopping for twisted bilayer graphene.
!! @param[out] Hodjj  Off-diagonal hopping element
!! @param[in]  dx,dy  In-plane displacements
!! @param[in]  Cjj    Coupling amplitude
!! @param[in]  Phijj  Phase factor
subroutine offdiagotBG(Hodjj,dx,dy,Cjj,Phijj)
   use constants
   use cell,                 only : aG
   real(dp), intent(in)::  dx, dy, Cjj, Phijj
   real(dp):: G1
   complex*16, intent(out) :: Hodjj

   !G1 = 4.0_dp*pi/sqrt(3.0_dp)/(aG)
   G1 = 4.0_dp*pi/(sqrt(3.0_dp))

   Hodjj = 2.0_dp*Cjj*cos(sqrt(3.0_dp)*G1*dx/2.0_dp)*( cos(G1*dy/2.0_dp - Phijj) ) &
 - 2.0_dp*Cjj*cos( G1*dy + Phijj )  - &
 cmplx_i*2.0_dp*sqrt(3.0_dp)*Cjj*sin( sqrt(3.0_dp)*G1*dx/2.0_dp )*  &
 ( sin( G1*dy/2.0_dp - Phijj)  )


   return
end subroutine offdiagotBG

!> @brief Compute diagonal onsite energy for twisted bilayer graphene.
!! @param[out] Hodjj  Diagonal onsite energy
!! @param[in]  dx,dy  In-plane displacements
!! @param[in]  Cjj    Coupling amplitude
!! @param[in]  Phijj  Phase factor
!! @param[in]  Cjj0   Constant offset
subroutine diagotBG(Hodjj,dx,dy,Cjj,Phijj,Cjj0)
   use constants
   use cell,                 only : aG
   real(dp), intent(in)::  dx, dy, Cjj, Phijj, Cjj0
   real(dp):: G1
   real(dp), intent(out) :: Hodjj

   !G1 = 4.0_dp*pi/sqrt(3.0_dp)/(aG)
   G1 = 4.0_dp*pi/(sqrt(3.0_dp))

   Hodjj = Cjj0 + 2.0_dp*Cjj*real((exp(-cmplx_i * G1 * dy) + 2.0_dp*exp(cmplx_i * G1 * dy / 2.0_dp) * cos(sqrt(3.0_dp)/2.0_dp * G1 * dx))*exp(cmplx_i*Phijj))

   return
end subroutine diagotBG

!> @brief Harmonic approximation for in-plane onsite variation.
!! @param[in]  dx,dy  Local displacements
!! @param[in]  A,B,C  Harmonic parameters
!! @param[out] Hjj    Onsite energy contribution
subroutine harmonicApprox(dx,dy,A,B,C, Hjj)
    use constants

    real(dp), intent(in) :: dx, dy, A, B, C
    real(dp), intent(out) :: Hjj
    real(dp) :: x1, y1, x2, y2, G1, m1, m2, a1, a2, a3, a4, alpha, beta, gammma
    real(dp) :: delta, Ax, Ay, acc, f1, D, c0, c1, phi

    if (abs(B-C) < 0.0000001) then
        print*, "WARNING: we have a singularity"
        D = (A-B)/(10**(-16))
    else
        D = (A-B)/(B-C)
    end if

    x1 = 0.0_dp
    y1 = 1.0_dp/sqrt(3.0_dp)
    x2 = 0.0_dp
    y2 = 0.0_dp
    G1 = 4.0_dp*pi/(sqrt(3.0_dp))
    m1 = cos(sqrt(3.0_dp)*G1*x1/2.0_dp)
    m2 = cos(sqrt(3.0_dp)*G1*x2/2.0_dp)
    a1 = cos(G1*y1) - cos(G1*y2)
    a3 = 2.0_dp*cos(G1*y1/2.0_dp) * m1 - 2.0_dp*cos(G1*y2/2.0_dp)* m2
    alpha = a1 + a3
    a2 = sin(G1*y1) - sin(G1*y2)
    a4 = 2.0_dp*sin(G1*y1/2.0_dp) * m1 - 2.0_dp*sin(G1*y2/2.0_dp)* m2
    beta = a2 - a4
    x1 = 0.0_dp
    y1 = 0.0_dp
    x2 = 0.0_dp
    y2 = 2.0_dp/sqrt(3.0_dp)
    G1 = 4.0_dp*pi/(sqrt(3.0_dp))
    m1 = cos(sqrt(3.0_dp)*G1*x1/2.0_dp)
    m2 = cos(sqrt(3.0_dp)*G1*x2/2.0_dp)
    a1 = cos(G1*y1) - cos(G1*y2)
    a3 = 2.0_dp*cos(G1*y1/2.0_dp) * m1 - 2.0_dp*cos(G1*y2/2.0_dp)* m2
    gammma = a1 + a3
    a2 = sin(G1*y1) - sin(G1*y2)
    a4 = 2.0_dp*sin(G1*y1/2.0_dp) * m1 - 2.0_dp*sin(G1*y2/2.0_dp)* m2
    delta = a2 - a4
    phi = atan((1.0_dp/(delta/beta*D-1.0_dp)*((delta*alpha-beta*gammma)/(beta*delta)))-(gammma)/(delta))
    if (abs(B-C) < 0.0000001_dp) then
        c1 = (10**(-16))/(2.0_dp*(gammma*cos(phi)+delta*sin(phi)))
        print*, "WARNING: we have a singularity2"
    else
        c1= (B-C)/(2.0*(gammma*cos(phi)+delta*sin(phi)))
    end if
    Ax = 0.0_dp
    Ay = 1_dp/sqrt(3.0_dp)
    G1 = 4.0_dp*pi/(sqrt(3.0_dp))
    c0 = A - 2.0_dp * c1 * cos(phi - G1 * Ay) - 4.0_dp * c1 * cos(G1 * Ay / 2.0_dp + phi) * cos(sqrt(3.0) * G1 * Ax / 2.0)
    acc = 1.0_dp
    G1 = 4.0_dp*pi/(sqrt(3.0_dp)*acc)
    f1 = 2.0_dp*c1*cos(phi-G1*dy) + 4.0_dp*c1*cos(G1*dy/2.0_dp + phi)*cos(sqrt(3.0_dp)*G1*dx/2.0_dp)
    Hjj = c0 + f1

    return

end subroutine

!> @brief Compute interlayer BL coupling phase factor HBL.
!! @param[out] HBL   Complex coupling value
!! @param[in]  dx    Relative x displacement
!! @param[in]  dy    Relative y displacement
!! @param[in]  tAB   Base interlayer hopping amplitude
subroutine interlayerBL(HBL,dx,dy,tAB)
   use constants
   use cell,                 only : aG
   real(dp), intent(in)::  dx, dy, tAB
   real(dp):: G1
   complex*16, intent(out) :: HBL

   HBL = tAB

   !HBL = exp(cmplx_i * ky * aG/sqrt(3.0_dp)) + 2.0_dp * exp(−cmplx_i * ky * aG/(2.0_dp * sqrt(3.0_dp)) * cos (kx * aG/2.0_dp)

   return
end subroutine interlayerBL

!> @brief Interlayer coupling for AB stacking between layers.
!! @param[out] HAB      Complex AB coupling
!! @param[in]  dx,dy    Relative in-plane displacement
!! @param[in]  tbt      Base coupling amplitude
!! @param[in]  posOrNeg Sign selector (+1/-1) for valley/rotation
subroutine interlayerBLAB(HAB,dx,dy,tbt,posOrNeg)
   use constants
   use cell,                 only : aG
   real(dp), intent(in) ::  dx, dy, tbt, posOrNeg
   real(dp) :: Gplusx, Gplusy, Gminx, Gminy, phi
   complex(dp), intent(out) :: HAB

   Gplusx = 4.0_dp * pi / (sqrt(3.0_dp) * aG) * (-sqrt(3.0_dp)/2.0_dp)
   Gplusy = 4.0_dp * pi / (sqrt(3.0_dp) * aG) * (1.0_dp/2.0_dp)
   Gminx = Gplusx
   Gminy = -Gplusy

   phi = 2.0_dp * pi / 3.0_dp

   HAB = tbt * (1.0_dp + exp(-cmplx_i * posOrNeg * phi) * exp(-cmplx_i * posOrNeg * (Gplusx * dx + Gplusy * dy)) &
            + exp(cmplx_i * posOrNeg * phi) * exp(-cmplx_i * posOrNeg * (Gminx * dx + Gminy * dy)))


   return
end subroutine interlayerBLAB

!> @brief Interlayer coupling for BA stacking between layers.
!! @param[out] HBA      Complex BA coupling
!! @param[in]  dx,dy    Relative in-plane displacement
!! @param[in]  tbt      Base coupling amplitude
!! @param[in]  posOrNeg Sign selector (+1/-1) for valley/rotation
subroutine interlayerBLBA(HBA,dx,dy,tbt,posOrNeg)

   use constants
   use cell,                 only : aG
   real(dp), intent(in) ::  dx, dy, tbt, posOrNeg
   real(dp) :: Gplusx, Gplusy, Gminx, Gminy, phi
   complex(dp), intent(out) :: HBA

   Gplusx = 4.0_dp * pi / (sqrt(3.0_dp) * aG) * (-sqrt(3.0_dp)/2.0_dp)
   Gplusy = 4.0_dp * pi / (sqrt(3.0_dp) * aG) * (1.0_dp/2.0_dp)
   Gminx = Gplusx
   Gminy = -Gplusy

   phi = -2.0_dp * pi / 3.0_dp ! negative sign for phi

   HBA = tbt * (1.0_dp + exp(-cmplx_i * posOrNeg * phi) * exp(-cmplx_i * posOrNeg * (Gplusx * dx + Gplusy * dy)) &
             + exp(cmplx_i * posOrNeg *phi) * exp(-cmplx_i * posOrNeg * (Gminx * dx + Gminy * dy)))


   return
end subroutine interlayerBLBA

!> @brief Interlayer coupling for AA stacking between layers.
!! @param[out] HAA      Complex AA coupling
!! @param[in]  dx,dy    Relative in-plane displacement
!! @param[in]  tbt      Base coupling amplitude
!! @param[in]  posOrNeg Sign selector (+1/-1) for valley/rotation
subroutine interlayerBLAA(HAA,dx,dy,tbt,posOrNeg)

   use constants
   use cell,                 only : aG
   real(dp), intent(in) ::  dx, dy, tbt, posOrNeg
   real(dp) :: Gplusx, Gplusy, Gminx, Gminy
   complex(dp), intent(out) :: HAA

   Gplusx = 4.0_dp * pi / (sqrt(3.0_dp) * aG) * (-sqrt(3.0_dp)/2.0_dp)
   Gplusy = 4.0_dp * pi / (sqrt(3.0_dp) * aG) * (1.0_dp/2.0_dp)
   Gminx = Gplusx
   Gminy = -Gplusy

   HAA = tbt * (1.0_dp +  exp(-cmplx_i * posOrNeg * (Gplusx * dx + Gplusy * dy)) &
                   +  exp(-cmplx_i * posOrNeg * (Gminx * dx + Gminy * dy)))

   return
end subroutine interlayerBLAA

!===============================================================================
! SUBROUTINE: HamHopping
!
! DESCRIPTION:
!   Calculates the hopping integrals between atomic sites in the tight-binding
!   Hamiltonian. This includes both in-plane and inter-layer hopping terms,
!   with support for various bilayer models and strain effects.
!
!> @brief Assemble hopping integrals (intra/inter-layer) into hopp(:,:).
!! @details Computes hopping amplitudes using selected bilayer/single-layer models,
!! including twist angles, strain/bubble effects, and optional magnetic phases.
!! Supports OpenMP for performance.
!! @par Side effects
!! - Updates hopp(:,:) in place
!! - May access/transform geometry data
!! - Applies magnetic field phase factors if enabled
!! @note Works with parameters set in HamInit and structures from neigh/atoms.
!! @see HamInit
subroutine HamHopping

   use neigh,                only : maxNeigh
   use atoms,                only : in1, in2, nAt, nAtC1, layerIndex, interlayerDistances
   use atoms,                only : phiForEffectiveModel, displacements, displacements_b, displacements_t
   use atoms,                only : RatInit
   use neigh,                only : Nneigh, NeighD, NList, neighCell, Nradii, NList2
#ifdef MPI
   use neigh,                only : rcvList
#endif /* MPI */
   use atoms,                only : Species, dIntLay, Rat, frac, AtomsSetCart, AtomsSetFrac, inode1, inode2
   use tbpar,                only : gIntLay, gn, tbnn, g0
   use magf,                 only : magfield, Bmag, BmagZeeman, HaldanePhase
   use constants,            only : pi, fluxq, cmplx_i, twopi, uB
   use cell,                 only : ucell, sCell, rcell
   use cell,                 only : aG, aBN
   use gauss,                only : gaussPotDefinedPositions
   use moireBLShift,         only : tauX1, tauY1, tauX2, tauY2
   use math
   use name,                 only : prefix, sysname

   integer :: i, j, ilvl, nlay
   real(dp) :: d, flux, phase, v1(3), v2(3), h1, h2
   real(dp) :: uvec(3), vvec(3), axb(3), vvecR(3)
   logical :: w, BilayerOneParameter, BilayerThreeParameters, FanZhang
   real(dp) :: Cab, tAB3, tAB4, Phiab, dx, dy, eps, realH, imagH, del(3), delta, eps2, ggg
   real(dp) :: realH_b, imagH_b
   real(dp) :: realH_t, imagH_t
   real(dp) :: dxTemp, dyTemp
   real(dp) :: dxTemp_b, dyTemp_b
   real(dp) :: dxTemp_t, dyTemp_t
   real(dp) :: dx_b, dy_b
   real(dp) :: dx_t, dy_t
   real(dp) :: PhiabG, PhiabBN
   real(dp) :: tr1, tr2, tr3, tr4, trDelta
   complex*16 :: Habjj, HBL
   complex*16 :: Habjj_b, Habjj_t
   logical :: GBNOffdiag

   real(dp) :: dxi, dyi, dxj, dyj

   complex(dp) :: HAA, HAB, HBA

   real(dp) :: posOrNeg

   logical :: l,ll, lll, twisted, twisted2, zz, deactivateUpperLayer, deactivateUpperLayers, paperOrientation

   integer :: numberOfDel1, numberOfDel2, numberOfDel3, numberOfInterlayerHoppings
   integer :: numberOfHAA1, numberOfHAA2, numberOfHAB1, numberOfHAB2, numberOfHBA1, numberOfHBA2

   real(dp) :: MoireBilayerTopAngle, MoireBilayerTopAngleGrad, twistAngle, twistAngleGrad, twistAngle2, twistAngleGrad2, shift2, shift1
   real(dp) :: phi2M
   real(dp) :: shift2_x, shift2_y
   real(dp) :: MoireBilayerBottomAngle, MoireBilayerBottomAngleGrad
   real(dp) :: sign2
   real(dp) :: twistedBilayerAngle, twistedBilayerAngleGrad, tbt

   logical :: rotateFirst

   character(len=80) :: str, BilayerModel, SinglelayerModel

   real(dp) :: t2Complex
   real(dp) :: t2

   real(dp) :: n1minm1, n2minm2, dist2, distXY, distXYZ
   real(dp) :: aCC, refDist

   logical :: prnt

   real(dp) :: dlij
   real(dp), pointer :: epsxx(:,:)=>NULL()
   real(dp), pointer :: epsyy(:,:)=>NULL()
   real(dp), pointer :: epsxy(:,:)=>NULL()
   real(dp), pointer :: fprimex(:,:)=>NULL()
   real(dp), pointer :: fprimey(:,:)=>NULL()
   integer, pointer :: indxImp(:)=>NULL()
   integer, pointer :: bubbleCenterX(:)=>NULL()
   integer, pointer :: bubbleCenterY(:)=>NULL()
   logical, pointer :: def(:)=>NULL()
   integer, pointer :: seed(:)

   real(dp), pointer :: onsiteShift(:)=>NULL()
   real(dp) :: bubbleSigmaR, bubbleSigmaR2, bubbleShift, bubbleRadius
   real(dp) :: bubbleC, bubbleR, bubbleTheta, minDiffBubbleX, minDiffBubbleY, diffBubbleX, diffBubbleY
   integer :: ic1, ic2, clock

   logical :: bigBubble, manyBubbles, bubbleInPlaneStrain, realisticBubbles, bubbleGauss

   real(dp) :: dist, vpppi, vpppi0, vppsigma, vppsigma0, BLdelta, per, rand
   integer :: nImp, totImp, ii, kk

   real(dp) :: yShift, xShift
   logical :: addShift, Bernal, bridge

   real(dp) :: gAngle, deltaAngle
   integer :: n(4)
   integer :: firstNN(3)
   integer :: k

   real(dp) :: expFactor

   real(dp) :: maxDist
   integer :: jj,jj1, jj2, iii

   integer :: nPeriod, nn
   real(dp) :: LMoire, u0

   logical :: randomStrain, periodicStrain, strainedMoire
   real(dp) :: umax

   type(cl_file) :: file
   integer u

   integer :: nnn(4)

   real (dp) :: sign1, minDelta

   real(dp) :: limit0, limit1, limit2, limit3, DBShiftRatio, boundaryWidth
   logical :: createBLDomainBoundary
   integer :: nnnn

   logical :: sublatticeDependent, sublatticeIndependent, useThetaIJ, useTheta

   real(dp) :: t1K, t2K, t3K, t4K, t5K, t2KA, t5KA, t2KB, t5KB, t2KN, t6K, t7K, t8K
   real(dp) :: t2KSL, t3KSL, t4KSL, t5KSL, t6KSL, t7KSL, t8KSL
   real(dp) :: t2KSLBN_B, t2KSLBN_N, t3KSLBN, t4KSLBN, t5KSLBN_B, t5KSLBN_N, t6KBN, t7KBN, t8KBN
   real(dp) :: V0, V3, V6, theta
   real(dp) :: V0AA, V3AA, V6AA
   real(dp) :: V0AB, V3AB, V6AB
   real(dp) :: V0AAp, V3AAp, V6AAp
   real(dp) :: V0ABp, V3ABp, V6ABp
   real(dp) :: V0BAp, V3BAp, V6BAp
   real(dp) :: V0BBp, V3BBp, V6BBp
   real(dp) :: lambda0, lambda3, lambda6, lambda6b
   real(dp) :: epsilon0, epsilon3, epsilon6, epsilon6b
   real(dp) :: x0, x3, x6, x6b
   real(dp) :: xi0, xi3, xi6, xi6b
   real(dp) :: kappa0, kappa3, kappa6, kappa6b
   real(dp) :: lambda0AA, lambda3AA, lambda6AA
   real(dp) :: epsilon0AA, epsilon3AA, epsilon6AA
   real(dp) :: x0AA, x3AA, x6AA
   real(dp) :: xi0AA, xi3AA, xi6AA
   real(dp) :: kappa0AA, kappa3AA, kappa6AA
   real(dp) :: lambda0AB, lambda3AB, lambda6AB
   real(dp) :: epsilon0AB, epsilon3AB, epsilon6AB
   real(dp) :: x0AB, x3AB, x6AB
   real(dp) :: xi0AB, xi3AB, xi6AB
   real(dp) :: kappa0AB, kappa3AB, kappa6AB
   real(dp) :: rbar

   real(dp) :: p1a0AA
   real(dp) :: p1b0AA
   real(dp) :: p1c0AA
   real(dp) :: p1d0AA
   real(dp) :: p1h0AA
   real(dp) :: p1j0AA
   real(dp) :: p1k0AA
   real(dp) :: p2a0AA
   real(dp) :: p2b0AA
   real(dp) :: p2c0AA
   real(dp) :: p2d0AA
   real(dp) :: p2h0AA
   real(dp) :: p2j0AA
   real(dp) :: p2k0AA
   real(dp) :: p3a0AA
   real(dp) :: p3b0AA
   real(dp) :: p3c0AA
   real(dp) :: p3d0AA
   real(dp) :: p3h0AA
   real(dp) :: p3j0AA
   real(dp) :: p3k0AA
   real(dp) :: p1a3AA
   real(dp) :: p1b3AA
   real(dp) :: p1c3AA
   real(dp) :: p1d3AA
   real(dp) :: p2a3AA
   real(dp) :: p2b3AA
   real(dp) :: p2c3AA
   real(dp) :: p2d3AA
   real(dp) :: p3a3AA
   real(dp) :: p3b3AA
   real(dp) :: p3c3AA
   real(dp) :: p3d3AA
   real(dp) :: p1a6AA
   real(dp) :: p1b6AA
   real(dp) :: p1c6AA
   real(dp) :: p1d6AA
   real(dp) :: p2a6AA
   real(dp) :: p2b6AA
   real(dp) :: p2c6AA
   real(dp) :: p2d6AA
   real(dp) :: p3a6AA
   real(dp) :: p3b6AA
   real(dp) :: p3c6AA
   real(dp) :: p3d6AA
   real(dp) :: a0AA
   real(dp) :: b0AA
   real(dp) :: c0AA
   real(dp) :: d0AA
   real(dp) :: h0AA
   real(dp) :: j0AA
   real(dp) :: k0AA
   real(dp) :: a3AA
   real(dp) :: b3AA
   real(dp) :: c3AA
   real(dp) :: d3AA
   real(dp) :: a6AA
   real(dp) :: b6AA
   real(dp) :: c6AA
   real(dp) :: d6AA
   real(dp) :: p1a0AB
   real(dp) :: p1b0AB
   real(dp) :: p1c0AB
   real(dp) :: p1d0AB
   real(dp) :: p1h0AB
   real(dp) :: p1j0AB
   real(dp) :: p1k0AB
   real(dp) :: p2a0AB
   real(dp) :: p2b0AB
   real(dp) :: p2c0AB
   real(dp) :: p2d0AB
   real(dp) :: p2h0AB
   real(dp) :: p2j0AB
   real(dp) :: p2k0AB
   real(dp) :: p3a0AB
   real(dp) :: p3b0AB
   real(dp) :: p3c0AB
   real(dp) :: p3d0AB
   real(dp) :: p3h0AB
   real(dp) :: p3j0AB
   real(dp) :: p3k0AB
   real(dp) :: p1a3AB
   real(dp) :: p1b3AB
   real(dp) :: p1c3AB
   real(dp) :: p1d3AB
   real(dp) :: p2a3AB
   real(dp) :: p2b3AB
   real(dp) :: p2c3AB
   real(dp) :: p2d3AB
   real(dp) :: p3a3AB
   real(dp) :: p3b3AB
   real(dp) :: p3c3AB
   real(dp) :: p3d3AB
   real(dp) :: p1a6AB
   real(dp) :: p1b6AB
   real(dp) :: p1c6AB
   real(dp) :: p1d6AB
   real(dp) :: p2a6AB
   real(dp) :: p2b6AB
   real(dp) :: p2c6AB
   real(dp) :: p2d6AB
   real(dp) :: p3a6AB
   real(dp) :: p3b6AB
   real(dp) :: p3c6AB
   real(dp) :: p3d6AB
   real(dp) :: a0AB
   real(dp) :: b0AB
   real(dp) :: c0AB
   real(dp) :: d0AB
   real(dp) :: h0AB
   real(dp) :: j0AB
   real(dp) :: k0AB
   real(dp) :: a3AB
   real(dp) :: b3AB
   real(dp) :: c3AB
   real(dp) :: d3AB
   real(dp) :: a6AB
   real(dp) :: b6AB
   real(dp) :: c6AB
   real(dp) :: d6AB
   real(dp) :: p1a0BA
   real(dp) :: p1b0BA
   real(dp) :: p1c0BA
   real(dp) :: p1d0BA
   real(dp) :: p1h0BA
   real(dp) :: p1j0BA
   real(dp) :: p1k0BA
   real(dp) :: p2a0BA
   real(dp) :: p2b0BA
   real(dp) :: p2c0BA
   real(dp) :: p2d0BA
   real(dp) :: p2h0BA
   real(dp) :: p2j0BA
   real(dp) :: p2k0BA
   real(dp) :: p3a0BA
   real(dp) :: p3b0BA
   real(dp) :: p3c0BA
   real(dp) :: p3d0BA
   real(dp) :: p3h0BA
   real(dp) :: p3j0BA
   real(dp) :: p3k0BA
   real(dp) :: p1a3BA
   real(dp) :: p1b3BA
   real(dp) :: p1c3BA
   real(dp) :: p1d3BA
   real(dp) :: p2a3BA
   real(dp) :: p2b3BA
   real(dp) :: p2c3BA
   real(dp) :: p2d3BA
   real(dp) :: p3a3BA
   real(dp) :: p3b3BA
   real(dp) :: p3c3BA
   real(dp) :: p3d3BA
   real(dp) :: p1a6BA
   real(dp) :: p1b6BA
   real(dp) :: p1c6BA
   real(dp) :: p1d6BA
   real(dp) :: p2a6BA
   real(dp) :: p2b6BA
   real(dp) :: p2c6BA
   real(dp) :: p2d6BA
   real(dp) :: p3a6BA
   real(dp) :: p3b6BA
   real(dp) :: p3c6BA
   real(dp) :: p3d6BA
   real(dp) :: a0BA
   real(dp) :: b0BA
   real(dp) :: c0BA
   real(dp) :: d0BA
   real(dp) :: h0BA
   real(dp) :: j0BA
   real(dp) :: k0BA
   real(dp) :: a3BA
   real(dp) :: b3BA
   real(dp) :: c3BA
   real(dp) :: d3BA
   real(dp) :: a6BA
   real(dp) :: b6BA
   real(dp) :: c6BA
   real(dp) :: d6BA
   real(dp) :: p1a0BB
   real(dp) :: p1b0BB
   real(dp) :: p1c0BB
   real(dp) :: p1d0BB
   real(dp) :: p1h0BB
   real(dp) :: p1j0BB
   real(dp) :: p1k0BB
   real(dp) :: p2a0BB
   real(dp) :: p2b0BB
   real(dp) :: p2c0BB
   real(dp) :: p2d0BB
   real(dp) :: p2h0BB
   real(dp) :: p2j0BB
   real(dp) :: p2k0BB
   real(dp) :: p3a0BB
   real(dp) :: p3b0BB
   real(dp) :: p3c0BB
   real(dp) :: p3d0BB
   real(dp) :: p3h0BB
   real(dp) :: p3j0BB
   real(dp) :: p3k0BB
   real(dp) :: p1a3BB
   real(dp) :: p1b3BB
   real(dp) :: p1c3BB
   real(dp) :: p1d3BB
   real(dp) :: p2a3BB
   real(dp) :: p2b3BB
   real(dp) :: p2c3BB
   real(dp) :: p2d3BB
   real(dp) :: p3a3BB
   real(dp) :: p3b3BB
   real(dp) :: p3c3BB
   real(dp) :: p3d3BB
   real(dp) :: p1a6BB
   real(dp) :: p1b6BB
   real(dp) :: p1c6BB
   real(dp) :: p1d6BB
   real(dp) :: p2a6BB
   real(dp) :: p2b6BB
   real(dp) :: p2c6BB
   real(dp) :: p2d6BB
   real(dp) :: p3a6BB
   real(dp) :: p3b6BB
   real(dp) :: p3c6BB
   real(dp) :: p3d6BB
   real(dp) :: a0BB
   real(dp) :: b0BB
   real(dp) :: c0BB
   real(dp) :: d0BB
   real(dp) :: h0BB
   real(dp) :: j0BB
   real(dp) :: k0BB
   real(dp) :: a3BB
   real(dp) :: b3BB
   real(dp) :: c3BB
   real(dp) :: d3BB
   real(dp) :: a6BB
   real(dp) :: b6BB
   real(dp) :: c6BB
   real(dp) :: d6BB

   real(dp) ::  lambda0AAp
   real(dp) ::  xi0AAp
   real(dp) ::  kappa0AAp
   real(dp) ::  lambda0bAAp
   real(dp) ::  xi0bAAp
   real(dp) ::  kappa0bAAp
   real(dp) ::  lambda0ABp
   real(dp) ::  xi0ABp
   real(dp) ::  kappa0ABp
   real(dp) ::  lambda0bABp
   real(dp) ::  xi0bABp
   real(dp) ::  kappa0bABp
   real(dp) ::  lambda0BAp
   real(dp) ::  xi0BAp
   real(dp) ::  kappa0BAp
   real(dp) ::  lambda0bBAp
   real(dp) ::  xi0bBAp
   real(dp) ::  kappa0bBAp
   real(dp) ::  lambda3ABp
   real(dp) ::  xi3ABp
   real(dp) ::  x3ABp
   real(dp) ::  lambda3bABp
   real(dp) ::  xi3bABp
   real(dp) ::  x3bABp
   real(dp) ::  lambda3BAp
   real(dp) ::  xi3BAp
   real(dp) ::  x3BAp
   real(dp) ::  lambda3bBAp
   real(dp) ::  xi3bBAp
   real(dp) ::  x3bBAp

   real(dp) :: theta12, theta21, d122, d132, d232, d232F
   integer :: mmm, mm, kkk

   real(dp) :: KaxirasCutoff, KaxirasCutoff2, distFact
   real(dp) :: SrivaniCutoff

   logical :: checking, SrivaniAA, SrivaniAB, useBNGKaxiras, fourLayers, threeLayers, middleTwist, useBNGSrivani, twoLayers, fourLayersSandwiched, fiveLayersSandwiched, sixLayersSandwiched, sevenLayersSandwiched, eightLayersSandwiched, tenLayersSandwiched, twentyLayersSandwiched, helicalTwistedMBM
   logical :: GBNtwoLayers, BNBNtwoLayers, encapsulatedFourLayers, encapsulatedFourLayersF2G2, t3GwithBN, t2GBN, t2BG, t3BG, BNt2GBN
   logical :: encapsulatedThreeLayers, encapsulatedFiveLayers, removeTopMoireInL2
   logical :: encapsulatedSixLayers, encapsulatedSevenLayers
   logical :: GBNtwoLayersF2G2s, BNBNtwoLayersF2G2s, GBNuseDisplacementFile, tBGOffDiag, tBGOffDiagPRB
   logical :: oneLayer
   logical :: tBGuseDisplacementFile
   logical :: addPressureDependence
   logical :: F2G2Model, useOldGrapheneF2G2
   real(dp) :: GlobalPhi,Globalang
   real(dp) :: GlobalPhi2,Globalang2
   real(dp) :: GlobalPhiL1, GlobalPhiL2, GlobalPhiL3, GlobalPhiL2a, GlobalPhiL2b

   real(dp) :: interlayerdistance
   real(dp) :: epsKax

   real (dp) :: a1AA
   real (dp) :: b1AA
   real (dp) :: c1AA
   real (dp) :: d1AA
   real (dp) :: a2AA
   real (dp) :: b2AA
   real (dp) :: c2AA
   !real (dp) :: a3AA
   !real (dp) :: b3AA
   !real (dp) :: c3AA
   real (dp) :: a1AB
   real (dp) :: b1AB
   real (dp) :: c1AB
   real (dp) :: d1AB
   real (dp) :: a2AB
   real (dp) :: b2AB
   real (dp) :: c2AB
   !real (dp) :: a3AB
   !real (dp) :: b3AB
   !real (dp) :: c3AB

   real (dp) :: V0_CB, lambda0_CB, xi0_CB, kappa0_CB
   real (dp) :: V3_BC, lambda3_BC, xi3_BC, x3_BC
   real (dp) :: V3_CB, lambda3_CB, xi3_CB, x3_CB
   real (dp) :: V0_CN, lambda0_CN, xi0_CN, kappa0_CN
   real (dp) :: V3_NC, lambda3_NC, xi3_NC, x3_NC
   real (dp) :: V3_CN, lambda3_CN, xi3_CN, x3_CN
   real (dp) :: theta_BC
   real (dp) :: theta_CB
   real (dp) :: theta_NC
   real (dp) :: theta_CN

   real (dp) :: epsilon0_CB
   real (dp) :: epsilon3_BC
   real (dp) :: epsilon3_CB
   real (dp) :: epsilon0_CN
   real (dp) :: epsilon3_NC
   real (dp) :: epsilon3_CN

   real (dp) :: c1_0
   real (dp) :: c1_1
   real (dp) :: c1_2
   real (dp) :: c2_0
   real (dp) :: c2_1
   real (dp) :: c2_2
   real (dp) :: c3_0
   real (dp) :: c3_1
   real (dp) :: c3_2
   real (dp) :: c4_0
   real (dp) :: c4_1
   real (dp) :: c4_2
   real (dp) :: c5_0
   real (dp) :: c5_1
   real (dp) :: c5_2
   real (dp) :: c6_0
   real (dp) :: c6_1
   real (dp) :: c6_2
   real (dp) :: c7_0
   real (dp) :: c7_1
   real (dp) :: c7_2
   real (dp) :: c8_0
   real (dp) :: c8_1
   real (dp) :: c8_2
   real (dp) :: c9_0
   real (dp) :: c9_1
   real (dp) :: c9_2
   real (dp) :: c10_0
   real (dp) :: c10_1
   real (dp) :: c10_2
   real (dp) :: c11_0
   real (dp) :: c11_1
   real (dp) :: c11_2
   real (dp) :: c12_0
   real (dp) :: c12_1
   real (dp) :: c12_2
   real (dp) :: c13_0
   real (dp) :: c13_1
   real (dp) :: c13_2
   real (dp) :: c14_0
   real (dp) :: c14_1
   real (dp) :: c14_2

   real (dp) :: tAA1
   real (dp) :: tAA2
   real (dp) :: tAB1
   real (dp) :: tAB2
   real (dp) :: tBA0
   real (dp) :: tBA2
   real (dp) :: tBA5

   real (dp) :: c0AAV0
   real (dp) :: c1AAV0
   real (dp) :: c2AAV0
   real (dp) :: c0ABV0
   real (dp) :: c1ABV0
   real (dp) :: c2ABV0

   real (dp) :: t2KSLGfromGBNA
   real (dp) :: t2KSLGfromGBNB
   real (dp) :: t3KSLGfromGBN
   real (dp) :: t4KSLGfromGBN
   real (dp) :: t5KSLGfromGBNA
   real (dp) :: t5KSLGfromGBNB
   real (dp) :: t6KSLGfromGBNA
   real (dp) :: t6KSLGfromGBNB
   real (dp) :: t7KSLGfromGBN
   real (dp) :: t8KSLGfromGBN
   real (dp) :: t2KSLBNfromGBNA
   real (dp) :: t2KSLBNfromGBNB
   real (dp) :: t3KSLBNfromGBN
   real (dp) :: t4KSLBNfromGBN
   real (dp) :: t5KSLBNfromGBNA
   real (dp) :: t5KSLBNfromGBNB
   real (dp) :: t6KSLBNfromGBNA
   real (dp) :: t6KSLBNfromGBNB
   real (dp) :: t7KSLBNfromGBN
   real (dp) :: t8KSLBNfromGBN

   real (dp) :: signChange

   real (dp) :: tAB

   real (dp) :: aGSrivani
   real (dp) :: t2Value

   integer :: countG1
   integer :: countG2
   integer :: countG3
   integer :: countG4
   integer :: countG5
   integer :: countG6
   integer :: countG7
   integer :: countG8
   integer :: countBN1
   integer :: countBN2
   integer :: countBN3
   integer :: countBN4
   integer :: countBN5
   integer :: countBN6
   integer :: countBN7
   integer :: countBN8

   logical :: singleLayerXYZ, changeLatticeParameterForSrivaniModel
   real (dp) :: Cabd, BfactorCab, z0 ! it should be ok to take the first neighbor as it is effective model with only 3 neighbors
   real (dp) :: CabG, CabBN
   real (dp) :: CabG_global
   logical :: distanceDependentEffectiveModel, findThetasGeometrically

   logical :: deactivateInterlayer, removeF2G2Flag, switchV3Sign, DCT, oldParameterSet, KoshinoIntralayer, forceBilayerF2G2Intralayer, MayouIntralayer
   logical :: deactivateInterlayer12, deactivateInterlayer23, deactivateInterlayer34
   logical :: deactivateInterlayerBG
   logical :: deactivateInterlayert2GBN1to2, deactivateInterlayert2GBN2to3
   logical :: deactivateInterlayert2BG2to3
   logical :: deactivateInterlayert3BG1to2, deactivateInterlayert3BG2to3, deactivateInterlayert3BG3to4, deactivateInterlayert3BG4to5, deactivateInterlayert3BG5to6
   logical :: useOnlyVAB, addExponentialDecayForDihedral, newFittingFunctions, onlyV0
   logical :: deactivateInterlayerTwisted, renormalizeHoppings, deactivateV6, deactivateV3
   real(dp) :: epsFactor, renormalizeHoppingFactorAAp, renormalizeHoppingFactorABp, renormalizeHoppingFactorBAp, renormalizeHoppingFactorBBp

   logical :: tBGSwitchDxDy
   logical :: renormalizeCoupling, KoshinoSR, writeData, addDisplacements, realStrain, onlyFirstNeighborRealStrain, BernalReadXYZ
   real(dp) :: couplingFactor, couplingFactor2
   logical :: differentCouplings, deactivateIntrasublattice, deactivateIntersublattice, deactivateIntraSublatticeForC

   real(dp) :: lc, rc, qsigma, qpi, Fc, aval, bval, avalC1B, avalC2B, avalC1N , avalC2N, bvalC1B, bvalC2B, bvalC1N , bvalC2N
   real(dp) :: GBNAngle, tBGAngle

   real(dp) :: aBN1
   real(dp) :: aBNR, aBN1R ! B to N distance
   real(dp) :: aGR, accR ! B to N distance
   real(dp) :: rotationAngle, epstBG

   integer :: counter1, nLayers
   logical :: readRigidXYZ, oppositedxdy, bulk, HaldaneOppositePhase
   logical :: onlyvppsigma, setnLayersToZero, HaldaneBothLayers_deprecated
   logical :: corrugatedInterlayerTwoCenter   ! optional: interlayer 2-centre with local p_z axes
   integer :: pzNGhost                         ! size of the MPI halo region of pzNormal
   integer :: pzNfallback                       ! # atoms that fell back to z_hat in ComputePzNormals
   ! Local p_z orbital axis per atom (+MPI halo). SAVEd only so MIO_Allocate can
   ! reuse the buffer across repeated HamHopping calls; scoped to this routine.
   real(dp), pointer :: pzNormal(:,:) => NULL()

   logical :: threeLayerShort
   real(dp) :: tA1B2, tA1A2, tB1B2, tA1A3, tB1A3
   real(dp) :: tA1B2_2, tA1A2_2, tB1B2_2, tA1A3_2, tB1A3_2
   real(dp) :: tB1A2, tA1B3
   real(dp) :: tA1A1_1, tB1B1_1, tA2A2_1, tA1B3_1, tB1A2_1
   real(dp) :: tA1A1_2, tB1B1_2, tA2A2_2, tA1B3_2, tB1A2_2

   logical :: Frank
   integer :: n1, n2, lllll, nnnnn, mB
   real(dp) :: diffx, diffy, mmphi

#ifdef DEBUG
   call MIO_Debug('HamHopping',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('ham')
#endif /* TIMER */

   !call MIO_InputParameter('readRigidXYZ',readRigidXYZ,.false.)
   !if (readRigidXYZ) then
   !   print*, Rat(:,1), RatInit(:,1)
   !   print*, "we replace Rat by RatInit... be careful, only tested for the specific case of the effective distant-dependent continuum model" !   Rat = RatInit
   !   print*, Rat(:,1), RatInit(:,1)
   !end if

   !i = 1000
   !j = 1
   !print*, (NeighD(:,j,i))
   !print*, (Rat(:,i))
   !print*, (Rat(:,NList(j,i)))

   !j = 2
   !print*, (NeighD(:,j,i))
   !print*, (Rat(:,i))
   !print*, (Rat(:,NList(j,i)))

   !j = 3
   !print*, (NeighD(:,j,i))
   !print*, (Rat(:,i))
   !print*, (Rat(:,NList(j,i)))

   counter1 = 0
   hopp = 0.0_dp
   numberOfInterlayerHoppings = 0
   flux = Bmag*pi/fluxq
   !print*, "for ", Bmag, " we have a flux of ", flux
   call MIO_InputParameter('MoireStrain',moireStrain,.false.)
   call MIO_InputParameter('FrankMagneticField',Frank,.false.)
   call MIO_InputParameter('MagField.Integer',mB,1)
   if (moireStrain) then
      if (.not. frac) call AtomsSetFrac()
      !$OMP PARALLEL DO PRIVATE(v1,v2)
      do i=1,nAt
         do j=1,Nneigh(i)
            h1 = 0.50_dp*hStr*(sin(twopi*Rat(1,i)*sCell)+sin(twopi*Rat(2,i)*sCell))
            h2 = 0.50_dp*hStr*(sin(twopi*Rat(1,NList(j,i))*sCell)+sin(twopi*Rat(2,NList(j,i))*sCell))
            d = sqrt(dot_product(NeighD(:,j,i),NeighD(:,j,i))) ! is this correct, shouldn't it be 1:2
            d = sqrt(d**2+(h2-h1)**2)
            hopp(j,i) = exp(-strFactor*(d/dCC-1.0_dp))
            if (Species(i)==Species(NList(j,i))) hopp(j,i) = -hopp(j,i)
         end do
      end do
      !$OMP END PARALLEL DO
   else

         numberOfDel1 = 0
         numberOfDel2 = 0
         numberOfDel3 = 0
         numberOfHAA1 = 0
         numberOfHAA2 = 0
         numberOfHAB1 = 0
         numberOfHAB2 = 0
         numberOfHBA1 = 0
         numberOfHBA2 = 0

         call MIO_InputParameter('changeLatticeParameterForSrivaniModel',changeLatticeParameterForSrivaniModel,.false.)
         !call MIO_InputParameter('changeLatticeParameterForSrivaniModel',changeLatticeParameterForSrivaniModel,.false.)
         if (changeLatticeParameterForSrivaniModel) then
            aGSrivani = 2.4389777651302801_dp
            call MIO_Print('We change aG for the calculation of rbar into 2.43 instead of 2.46 to agree with Srivanis fitting','ham')
         else
            aGSrivani = aG
         end if

         call MIO_InputParameter('removeF2G2Flag',removeF2G2Flag,.false.)
         if (removeF2G2Flag) then
            F2G2Model = .false.
         else
            call MIO_InputParameter('F2G2Model',F2G2Model,.true.)
         end if
         call MIO_InputParameter('threeLayerShort',threeLayerShort,.false.)
         call MIO_InputParameter('GBNtwoLayersF2G2s',GBNtwoLayersF2G2s,.false.)
         call MIO_InputParameter('GBNtwoLayers',GBNtwoLayers,.false.)
         call MIO_InputParameter('BNBNtwoLayers',BNBNtwoLayers,.false.)
         call MIO_InputParameter('tBGOffDiag',tBGOffDiag,.false.)
         call MIO_InputParameter('tBGOffDiagPRB',tBGOffDiagPRB,.false.)
         call MIO_InputParameter('Latticepercent',eps,-0.018181818181818_dp)
         call MIO_InputParameter('MoirePotCab',Cab,0.01134_dp)
         call MIO_InputParameter('MoirePotCabBN',CabBN,0.004418_dp)
         Cab = Cab/g0
         CabBN = CabBN/g0
         call MIO_InputParameter('BilayerOneParameter',BilayerOneParameter,.false.)
         call MIO_InputParameter('BilayerThreeParameters',BilayerThreeParameters,.false.)
         call MIO_InputParameter('BfactorCab',BfactorCab,3.3_dp)
         call MIO_InputParameter('InterlayerDistance',z0,3.35_dp)
         if (BilayerOneParameter) then
             call MIO_InputParameter('BilayertAB1',tAB1,0.361_dp)
             tAB1 = -tAB1/g0 ! Add minus sign to compensate for intrinsic minus sign
         else if (BilayerThreeParameters) then
             call MIO_InputParameter('BilayertAB1',tAB1,0.361_dp)
             tAB1 = -tAB1/g0
             call MIO_InputParameter('BilayertAB3',tAB3,0.283_dp)
             tAB3 = -tAB3/g0
             call MIO_InputParameter('BilayertAB4',tAB4,0.138_dp)
             tAB4 = -tAB4/g0
         else if (F2G2Model) then
             call MIO_InputParameter('BilayertAA1',tAA1,0.09244_dp)
             tAA1 = -tAA1/g0
             call MIO_InputParameter('BilayertAA3',tAA2,-0.02299_dp)
             tAA2 = -tAA2/g0
             call MIO_InputParameter('BilayertAB1',tAB1,0.1391_dp)
             tAB1 = -tAB1/g0
             call MIO_InputParameter('BilayertAB3',tAB2,-0.07211_dp)
             tAB2 = -tAB2/g0
             call MIO_InputParameter('BilayertBA0',tBA0,0.331_dp)
             tBA0 = -tBA0/g0
             call MIO_InputParameter('BilayertBA0',tBA2,-0.01016_dp)
             tBA2 = -tBA2/g0
             call MIO_InputParameter('BilayertBA0',tBA5,0.0001_dp)
             tBA5 = -tBA5/g0
             !print*, "tAA1, tAA2, tAB1, tAB2, tBA0, tBA2, tBA5:"
             !print*, tAA1, tAA2, tAB1, tAB2, tBA0, tBA2, tBA5
         !else
         !    print*, "You didnt specify any Bernal stacked parameters, is that correct?"
         end if
         call MIO_InputParameter('TrilayerFanZhang',FanZhang,.false.)
         if (FanZhang) then
             call MIO_InputParameter('TrilayerGamma3',tr1,0.0_dp)
             call MIO_InputParameter('TrilayerGamma3',tr3,0.0_dp)
             call MIO_InputParameter('TrilayerGamma4',tr4,0.0_dp)
             call MIO_InputParameter('TrilayerGamma2',tr2,0.0_dp)
             call MIO_InputParameter('TrilayerDelta',trDelta,0.0_dp)
             tr1 = tr1/g0
             tr2 = tr2/g0
             tr3 = tr3/g0
             tr4 = tr4/g0
             !tAB1 = -tAB1/g0 ! Add minus sign to compensate for intrinsic minus sign
         !else
         !    print*, "You didnt specify any trilayer parameters, is that correct?"
         end if
         call MIO_InputParameter('MoirePotPhiab',Phiab,0.342084533390889_dp)
         if (tBGOffDiag) then
            !if (tBGOffDiagPRB) then
            !   call MIO_InputParameter('MoirePotPhiabG',PhiabG,0.0_dp)
            !else
               PhiabG = 0.0_dp
            !end if
         else
            call MIO_InputParameter('MoirePotPhiabG',PhiabG,3.5_dp)
         end if
         PhiabG = PhiabG*pi/180.0_dp
         call MIO_InputParameter('MoirePotPhiabBN',PhiabBN,26.1_dp)
         PhiabBN = PhiabBN*pi/180.0_dp

         call MIO_InputParameter('TypeOfSystem',str,'Graphene')
         call MIO_InputParameter('TypeOfBL',BilayerModel,'None')
         call MIO_InputParameter('TypeOfSL',SinglelayerModel,'None') ! then we don't have to change the code and keep the bilayer parts even for the single layer
         call MIO_InputParameter('addExponentialDecayForDihedral',addExponentialDecayForDihedral,.false.)
         call MIO_InputParameter('MoireBilayerTopAngle',MoireBilayerTopAngle,0.0_dp)
         MoireBilayerTopAngleGrad = MoireBilayerTopAngle*pi/180.0_dp
         call MIO_InputParameter('MoireBilayerBottomAngle',MoireBilayerBottomAngle,0.0_dp)
         MoireBilayerBottomAngleGrad = MoireBilayerBottomAngle*pi/180.0_dp
         call MIO_InputParameter('MoireBLDeactivateUpperLayer',deactivateUpperLayer,.false.)
         call MIO_InputParameter('MoiretDBLDeactivateUpperLayers',deactivateUpperLayers,.false.)
         call MIO_InputParameter('MoireTwisted',twisted,.false.)
         call MIO_InputParameter('MoireTwistAngle',twistAngle,0.0_dp)
         twistAngleGrad = twistAngle*pi/180.0_dp
         call MIO_InputParameter('MoireAddSecondMoire',zz,.false.)
         call MIO_InputParameter('basedOnMoireCellParameters',ll,.false.)
         call MIO_InputParameter('MoireFirstMoireMassFactor',sign1,1.0_dp)
         call MIO_InputParameter('WriteDataFiles',writeData,.false.)
         if (writeData) then
            open(587,FILE='HAB')
            open(586,FILE='BfactorExp.dat')
         end if
         if (ll) then
             call MIO_InputParameter('MoireCellParameters',nnn,[0,0,0,0])
             ggg = nnn(1)**2 + nnn(2)**2 + nnn(1)*nnn(2)
             delta = sqrt(real(nnn(3)**2 + nnn(4)**2 + nnn(3)*nnn(4))/ggg)
             phiForEffectiveModel = acos((2.0_dp*nnn(1)*nnn(3)+2.0_dp*nnn(2)*nnn(4) + nnn(1)*nnn(4) + nnn(2)*nnn(3))/(2.0_dp*delta*ggg))
             !print*, "phiForEffectiveModel: ", phiForEffectiveModel/pi*180.0_dp
             twistAngleGrad = phiForEffectiveModel
             !ggg = nnn(1)**2 + nnn(2)**2 + nnn(1)*nnn(2)
             !eps = sqrt(real(nnn(3)**2 + nnn(4)**2 + nnn(3)*nnn(4))/ggg)
             !eps = -0.02463661_dp
             !eps2 = eps
         else
             call MIO_InputParameter('MoireTwistAngle',twistAngle,0.0_dp)
             call MIO_InputParameter('MoireTwistAngle2',twistAngle2,0.0_dp)
             twistAngleGrad = twistAngle*pi/180.0_dp
             !print*, "phiForEffectiveModel: ", twistAngleGrad/pi*180.0_dp
             twistAngleGrad2 = twistAngle2*pi/180.0_dp
         end if
         call MIO_InputParameter('MoireTwisted2',twisted2,.false.)
         call MIO_InputParameter('MoireLayerShift1',shift1,0.0_dp)
         call MIO_InputParameter('MoireLayerShift2',shift2,0.0_dp)
         call MIO_InputParameter('MoireSecondMoireRotateFirst',rotateFirst,.true.)
         call MIO_InputParameter('minDelta',minDelta,1.0_dp)
         phi2M= atan((1.0_dp+eps)*sin(twistAngleGrad2)/((1.0_dp+eps)*cos(twistAngleGrad2)-1.0_dp))
         shift2_x = shift2*cos(phi2M)
         shift2_y = shift2*sin(phi2M)


         !call MIO_InputParameter('twistedBilayerAngle',twistedBilayerAngle,0.0_dp)
         call MIO_InputParameter('MoireCellParameters',n,[0,0,0,0])
         gAngle = n(1)**2 + n(2)**2 + n(1)*n(2)
         deltaAngle = sqrt(real(n(3)**2 + n(4)**2 + n(3)*n(4))/gAngle)
         twistedBilayerAngle = acos((2.0_dp*n(1)*n(3)+2.0_dp*n(2)*n(4) + n(1)*n(4) + n(2)*n(3))/(2.0_dp*deltaAngle*gAngle))
         !print*, "angle coming from Cell Parameters: ", twistedBilayerAngle/(pi/180.0_dp)
         twistedBilayerAngleGrad = twistedBilayerAngle
         !print*, "angle in grad coming from Cell Parameters: ", twistedBilayerAngleGrad
         call MIO_InputParameter('twistedBLtbt',tbt,0.113_dp)
         tbt = tbt/g0
         !print*, "tbt (in gamma0) = ", tbt

         !call MIO_InputParameter('BLdelta',BLdelta,0.184*aG)
         ! Koshino
         aCC = aG/sqrt(3.0_dp)
         ! vpppi0 is one parameter with two roles: it scales the pi part of the
         ! two-centre INTERLAYER hopping in every model, and with
         ! KoshinoIntralayer it is also the INTRALAYER nearest-neighbour hopping.
         ! Its default therefore follows the intralayer model: 2.7 eV, the
         ! original Moon-Koshino value, with KoshinoIntralayer; 3.5 eV, the
         ! calibration that goes with the F2G2-type intralayer models, otherwise.
         call MIO_InputParameter('KoshinoIntralayer',KoshinoIntralayer,.false.)
         if (KoshinoIntralayer) then
            call MIO_InputParameter('vpppi0',vpppi0,2.7_dp)
            call MIO_Print('Two-centre Vpppi0 = '//trim(num2str(vpppi0,4))// &
              ' eV (intralayer and interlayer; default 2.7 with KoshinoIntralayer)','ham')
         else
            call MIO_InputParameter('vpppi0',vpppi0,3.5_dp)
            if (.not. MIO_StringComp(BilayerModel,'None')) then
               call MIO_Print('Two-centre Vpppi0 = '//trim(num2str(vpppi0,4))// &
                 ' eV (interlayer pi term only; default 3.5 without KoshinoIntralayer)','ham')
            end if
         end if
         vpppi0 = vpppi0/g0
         call MIO_InputParameter('vppsigma0',vppsigma0,0.48_dp)
         call MIO_InputParameter('BLdelta',BLdelta,0.184_dp*aG)
         vppsigma0 = -vppsigma0/g0
         ! Optional: interlayer two-center hopping with local p_z orbital axes
         ! (corrugated sheets). Default .false. keeps the legacy global-z
         ! expression with the same operands; radial Vpppi/Vppsigma/BLdelta unchanged.
         call MIO_InputParameter('corrugatedInterlayerTwoCenter',corrugatedInterlayerTwoCenter,.false.)
         ! Mayou
         qpi = aCC*2.218_dp
         lc = 0.265_dp
         rc = 6.14_dp

         call MIO_InputParameter('InterlayerDistance',interlayerdistance,3.22_dp)
         qsigma = interlayerdistance * 2.218_dp
         call MIO_InputParameter('distanceDependentEffectiveModel',distanceDependentEffectiveModel,.false.)
         call MIO_InputParameter('F2G2Model',F2G2Model,.true.) ! default true
         call MIO_InputParameter('useOldGrapheneF2G2',useOldGrapheneF2G2,.false.) ! default true
         call MIO_InputParameter('deactivateInterlayer',deactivateInterlayer,.false.)
         call MIO_InputParameter('deactivateInterlayer12',deactivateInterlayer12,.false.)
         call MIO_InputParameter('deactivateInterlayer23',deactivateInterlayer23,.false.)
         call MIO_InputParameter('deactivateInterlayer34',deactivateInterlayer34,.false.)
         call MIO_InputParameter('deactivateInterlayert2GBN1to2',deactivateInterlayert2GBN1to2,.false.)
         call MIO_InputParameter('deactivateInterlayert2GBN2to3',deactivateInterlayert2GBN2to3,.false.)
         call MIO_InputParameter('deactivateInterlayert2BG2to3',deactivateInterlayert2BG2to3,.false.)
         call MIO_InputParameter('deactivateInterlayert3BG1to2',deactivateInterlayert3BG1to2,.false.)
         call MIO_InputParameter('deactivateInterlayert3BG2to3',deactivateInterlayert3BG2to3,.false.)
         call MIO_InputParameter('deactivateInterlayert3BG3to4',deactivateInterlayert3BG3to4,.false.)
         call MIO_InputParameter('deactivateInterlayert3BG4to5',deactivateInterlayert3BG4to5,.false.)
         call MIO_InputParameter('deactivateInterlayert3BG5to6',deactivateInterlayert3BG5to6,.false.)
         call MIO_InputParameter('deactivateInterlayerBG',deactivateInterlayerBG,.false.)
         call MIO_InputParameter('deactivateInterlayerTwisted',deactivateInterlayerTwisted,.false.)
         !call MIO_Print('Note that for bilayer systems, we use the F2G2 model by default for interlayer (unless it is a twisted system) and intralayer terms','ham')

         call MIO_InputParameter('periodicStrainPeriod',nPeriod,1)
         call MIO_InputParameter('strainedMoireMaxDisplacement',umax,0.5_dp)
         call MIO_InputParameter('strainedMoire',strainedMoire,.false.)
         call MIO_InputParameter('SuperCell',sCell,60)
         LMoire = norm(ucell(:,1))/sCell

         call MIO_InputParameter('twistedBLAddShift',addShift,.false.)
         call MIO_InputParameter('BernalShift',Bernal,.false.)
         call MIO_InputParameter('bridgeShift',bridge,.false.)
         if (addShift) then
            if (Bernal) then
               xShift = aG/2.0_dp
               yShift = aG/sqrt(3.0_dp)/2.0_dp
               !print*, "yShift= ", yShift, aG/3.0_dp
            else if (bridge) then
               xShift = aG/2.0_dp/2.0_dp
               yShift = aG/sqrt(3.0_dp)/2.0_dp/2.0_dp
            end if
         else
            xShift = 0.0_dp
            yShift = 0.0_dp
            !print*, "yShift= ", yShift, aG/3.0_dp
         end if

         if (frac) call AtomsSetCart()
         if (GBNtwoLayersF2G2s) then
            call MIO_Print('defining the GBNtwoLayersF2G2 parameters','ham')
            t1K = g0/g0
            call MIO_InputParameter('SingleLayert2KSL',t2KSLGfromGBNA,-0.24498_dp)
            call MIO_InputParameter('SingleLayert2KSL',t2KSLGfromGBNB,-0.24523_dp)
            call MIO_InputParameter('SingleLayert3K',t3KSLGfromGBN,0.19334_dp)
            call MIO_InputParameter('SingleLayert4K',t4KSLGfromGBN,-0.0_dp)
            call MIO_InputParameter('SingleLayert5KSL',t5KSLGfromGBNA,-0.06618_dp)
            call MIO_InputParameter('SingleLayert5KSL',t5KSLGfromGBNB,-0.06624_dp)
            call MIO_InputParameter('BilayertK6A',t6KSLGfromGBNA,0.0_dp)
            call MIO_InputParameter('BilayertK6A',t6KSLGfromGBNB,0.0_dp)
            call MIO_InputParameter('BilayertK7A',t7KSLGfromGBN,0.0_dp)
            call MIO_InputParameter('BilayertK8A',t8KSLGfromGBN,0.0_dp)
            !call MIO_InputParameter('SingleLayert2KSL',t2KSLGfromGBNA,-0.24463_dp)
            !call MIO_InputParameter('SingleLayert2KSL',t2KSLGfromGBNB,-0.24466_dp)
            !call MIO_InputParameter('SingleLayert3K',t3KSLGfromGBN,0.30358_dp)
            !call MIO_InputParameter('SingleLayert4K',t4KSLGfromGBN,-0.021313_dp)
            !call MIO_InputParameter('SingleLayert5KSL',t5KSLGfromGBNA,-0.059532_dp)
            !call MIO_InputParameter('SingleLayert5KSL',t5KSLGfromGBNB,-0.059609_dp)
            !call MIO_InputParameter('BilayertK6A',t6KSLGfromGBNA,0.021509_dp)
            !call MIO_InputParameter('BilayertK6A',t6KSLGfromGBNB,0.022454_dp)
            !call MIO_InputParameter('BilayertK7A',t7KSLGfromGBN,0.015464_dp)
            !call MIO_InputParameter('BilayertK8A',t8KSLGfromGBN,0.023004_dp)
            t2KSLGfromGBNA = t2KSLGfromGBNA/g0
            t2KSLGfromGBNB = t2KSLGfromGBNB/g0
            t3KSLGfromGBN = t3KSLGfromGBN/g0
            t4KSLGfromGBN = t4KSLGfromGBN/g0
            t5KSLGfromGBNA = t5KSLGfromGBNA/g0
            t5KSLGfromGBNB = t5KSLGfromGBNB/g0
            t6KSLGfromGBNA = t6KSLGfromGBNA/g0
            t6KSLGfromGBNB = t6KSLGfromGBNB/g0
            t7KSLGfromGBN = t7KSLGfromGBN/g0
            t8KSLGfromGBN = t8KSLGfromGBN/g0
            call MIO_InputParameter('SingleLayert2KSL',t2KSLBNfromGBNA,-0.081055_dp)
            call MIO_InputParameter('SingleLayert2KSL',t2KSLBNfromGBNB,-0.24562_dp)
            call MIO_InputParameter('SingleLayert3K',t3KSLBNfromGBN,0.15399_dp)
            call MIO_InputParameter('SingleLayert4K',t4KSLBNfromGBN,-0.0_dp)
            call MIO_InputParameter('SingleLayert5KSL',t5KSLBNfromGBNA,-0.065654_dp)
            call MIO_InputParameter('SingleLayert5KSL',t5KSLBNfromGBNB,-0.04892_dp)
            call MIO_InputParameter('BilayertK6A',t6KSLBNfromGBNA,0.0_dp)
            call MIO_InputParameter('BilayertK6A',t6KSLBNfromGBNB,0.0_dp)
            call MIO_InputParameter('BilayertK7A',t7KSLBNfromGBN,0.0_dp)
            call MIO_InputParameter('BilayertK8A',t8KSLBNfromGBN,0.0_dp)
            !call MIO_InputParameter('SingleLayert2KSL',t2KSLBNfromGBNA,-0.080464_dp)
            !call MIO_InputParameter('SingleLayert2KSL',t2KSLBNfromGBNB,-0.24597_dp)
            !call MIO_InputParameter('SingleLayert3K',t3KSLBNfromGBN,0.26486_dp)
            !call MIO_InputParameter('SingleLayert4K',t4KSLBNfromGBN,-0.049096_dp)
            !call MIO_InputParameter('SingleLayert5KSL',t5KSLBNfromGBNA,-0.051535_dp)
            !call MIO_InputParameter('SingleLayert5KSL',t5KSLBNfromGBNB,-0.040432_dp)
            !call MIO_InputParameter('BilayertK6A',t6KSLBNfromGBNA,0.035441_dp)
            !call MIO_InputParameter('BilayertK6A',t6KSLBNfromGBNB,0.023648_dp)
            !call MIO_InputParameter('BilayertK7A',t7KSLBNfromGBN,0.007505_dp)
            !call MIO_InputParameter('BilayertK8A',t8KSLBNfromGBN,0.018199_dp)
            !call MIO_InputParameter('BilayertK7A',t7KSLBNfromGBN,0.0_dp)
            !call MIO_InputParameter('BilayertK8A',t8KSLBNfromGBN,0.0_dp)
            t2KSLBNfromGBNA = t2KSLBNfromGBNA/g0
            t2KSLBNfromGBNB = t2KSLBNfromGBNB/g0
            t3KSLBNfromGBN = t3KSLBNfromGBN/g0
            t4KSLBNfromGBN = t4KSLBNfromGBN/g0
            t5KSLBNfromGBNA = t5KSLBNfromGBNA/g0
            t5KSLBNfromGBNB = t5KSLBNfromGBNB/g0
            t6KSLBNfromGBNA = t6KSLBNfromGBNA/g0
            t6KSLBNfromGBNB = t6KSLBNfromGBNB/g0
            t7KSLBNfromGBN = t7KSLBNfromGBN/g0
            t8KSLBNfromGBN = t8KSLBNfromGBN/g0
         else if (F2G2Model) then
            t1K = g0/g0
            if (useOldGrapheneF2G2) then
                call MIO_InputParameter('SingleLayert2KSL',t2KSL,-0.21264_dp)
                call MIO_InputParameter('SingleLayert3K',t3KSL,0.23442_dp)
                call MIO_InputParameter('SingleLayert4K',t4KSL,-0.05350_dp)
                !call MIO_InputParameter('SingleLayert4K',t4KSL,0.0_dp)
                call MIO_InputParameter('SingleLayert5KSL',t5KSL,-0.07326_dp)
                call MIO_InputParameter('BilayertK6A',t6K,0.0_dp)
                call MIO_InputParameter('BilayertK7A',t7K,0.0_dp)
                call MIO_InputParameter('BilayertK8A',t8K,0.0_dp)
            !else if (threeLayerShort) then
            !    call MIO_InputParameter('SingleLayert2KSL',t2KSL,-0.2238_dp)
            !    call MIO_InputParameter('SingleLayert3K',t3KSL,0.1859_dp)
            !    call MIO_InputParameter('SingleLayert4K',t4KSL,0.0_dp)
            !    call MIO_InputParameter('SingleLayert5KSL',t5KSL,-0.0446_dp)
            !    call MIO_InputParameter('BilayertK6A',t6K,0.0_dp)
            !    call MIO_InputParameter('BilayertK7A',t7K,0.0_dp)
            !    call MIO_InputParameter('BilayertK8A',t8K,0.0_dp)
            else
                call MIO_InputParameter('SingleLayert2KSL',t2KSL,-0.2354_dp)
                call MIO_InputParameter('SingleLayert3K',t3KSL,0.1877_dp)
                call MIO_InputParameter('SingleLayert4K',t4KSL,0.0_dp)
                call MIO_InputParameter('SingleLayert5KSL',t5KSL,-0.0633_dp)
                call MIO_InputParameter('BilayertK6A',t6K,0.0_dp)
                call MIO_InputParameter('BilayertK7A',t7K,0.0_dp)
                call MIO_InputParameter('BilayertK8A',t8K,0.0_dp)
            end if
            t2KSL = t2KSL/g0
            t3KSL = t3KSL/g0
            t4KSL = t4KSL/g0
            t5KSL = t5KSL/g0
            t6K = t6K/g0
            t7K = t7K/g0
            t8K = t8K/g0
            ! BNBN F2G2
            call MIO_InputParameter('SingleLayert2KSL',t2KSLBN_B,-0.0542_dp)
            call MIO_InputParameter('SingleLayert2KSL',t2KSLBN_N,-0.2228_dp)
            call MIO_InputParameter('SingleLayert3K',t3KSLBN,0.1329_dp)
            call MIO_InputParameter('SingleLayert4K',t4KSLBN,0.0_dp)
            !call MIO_InputParameter('SingleLayert4K',t4KSL,0.0_dp)
            call MIO_InputParameter('SingleLayert5KSL',t5KSLBN_B,-0.0566_dp)
            call MIO_InputParameter('SingleLayert5KSL',t5KSLBN_N,-0.0429_dp)
            call MIO_InputParameter('BilayertK6A',t6KBN,0.0_dp)
            call MIO_InputParameter('BilayertK7A',t7KBN,0.0_dp)
            call MIO_InputParameter('BilayertK8A',t8KBN,0.0_dp)
            t2KSLBN_B = t2KSLBN_B/g0
            t2KSLBN_N = t2KSLBN_N/g0
            t3KSLBN = t3KSLBN/g0
            t4KSLBN = t4KSLBN/g0
            t5KSLBN_B = t5KSLBN_B/g0
            t5KSLBN_N = t5KSLBN_N/g0
            t6KBN = t6KBN/g0
            t7KBN = t7KBN/g0
            t8KBN = t8KBN/g0
         else
            t1K = g0/g0
            t2KSL = -0.2425_dp/g0
            t3KSL = 0.2656_dp/g0
            t4KSL = -0.0235_dp/g0
            t5KSL = -0.0524_dp/g0
            t6K = 0.0209_dp/g0
            t7K = 0.0148_dp/g0
            t8K = 0.0211_dp/g0
            t2KA = t2KSL
            t2KB = t2KSL
            t5KA = t5KSL
            t5KB = t5KSL
         end if
         call MIO_InputParameter('forceBilayerF2G2Intralayer',forceBilayerF2G2Intralayer,.false.)
         if (forceBilayerF2G2Intralayer) then
             t1K = g0/g0
             call MIO_InputParameter('Bilayert2KA',t2KA,-0.2235_dp)
             call MIO_InputParameter('Bilayert2KB',t2KB,-0.2260_dp)
             call MIO_InputParameter('Bilayert3K',t3K,0.1984_dp)
             call MIO_InputParameter('Bilayert4K',t4K,0.0_dp)
             call MIO_InputParameter('Bilayert5KA',t5KA,-0.04016_dp)
             call MIO_InputParameter('Bilayert5KB',t5KB,-0.0404_dp)
             call MIO_InputParameter('BilayertK6A',t6K,0.0_dp)
             call MIO_InputParameter('BilayertK7A',t7K,0.0_dp)
             call MIO_InputParameter('BilayertK8A',t8K,0.0_dp)
             t2KA = t2KA/g0
             t2KB = t2KB/g0
             t3K = t3K/g0
             t4K = t4K/g0
             t5KA = t5KA/g0
             t5KB = t5KB/g0
             t6K = t6K/g0
             t7K = t7K/g0
             t8K = t8K/g0
         end if
         call MIO_InputParameter('t3BG',t3BG,.false.)
         if (t3BG) then
             t1K = g0/g0
             call MIO_InputParameter('Bilayert2KA',t2KA,-0.2235_dp)
             call MIO_InputParameter('Bilayert2KB',t2KB,-0.2260_dp)
             call MIO_InputParameter('Bilayert3K',t3K,0.1984_dp)
             call MIO_InputParameter('Bilayert4K',t4K,0.0_dp)
             call MIO_InputParameter('Bilayert5KA',t5KA,-0.04016_dp)
             call MIO_InputParameter('Bilayert5KB',t5KB,-0.0404_dp)
             call MIO_InputParameter('BilayertK6A',t6K,0.0_dp)
             call MIO_InputParameter('BilayertK7A',t7K,0.0_dp)
             call MIO_InputParameter('BilayertK8A',t8K,0.0_dp)
             t2KA = t2KA/g0
             t2KB = t2KB/g0
             t3K = t3K/g0
             t4K = t4K/g0
             t5KA = t5KA/g0
             t5KB = t5KB/g0
             t6K = t6K/g0
             t7K = t7K/g0
             t8K = t8K/g0
             !print*, "you are inside t3BG", t2KA, t2KB
         end if
         if (MIO_StringComp(BilayerModel,'BLKaxiras')) then
            call MIO_InputParameter('useBNGKaxiras',useBNGKaxiras,.false.)
            tAB = 0.29_dp
            if (useBNGKaxiras) then
                lambda0_CB = 0.3905_dp/g0
                epsilon0_CB = 1.5426_dp
                !x0_CB = 0.0_dp
                kappa0_CB = 1.8229_dp
                lambda3_BC = -0.0588_dp/g0
                epsilon3_BC = 3.0827_dp
                x3_BC = 0.6085_dp
                !kappa3_BC = 0.0_dp
                lambda3_CB = -0.0651_dp/g0
                epsilon3_CB = 3.7998_dp
                x3_CB = 0.6341_dp
                !kappa3_CB = 0.0_dp
                xi0_CB = epsilon0_CB
                xi3_BC = epsilon3_BC
                xi3_CB = epsilon3_CB
                lambda0_CN = 0.2517_dp/g0
                epsilon0_CN = 1.6061_dp
                !x0_CN = 0.0_dp
                kappa0_CN = 2.1909_dp
                lambda3_NC = -0.0606_dp/g0
                epsilon3_NC = 3.3502_dp
                x3_NC = 0.5142_dp
                !kappa3_NC = 0.0_dp
                lambda3_CN = -0.0465_dp/g0
                epsilon3_CN = 3.0464_dp
                x3_CN = 0.5264_dp
                !kappa3_CN = 0.0_dp
                xi0_CN = epsilon0_CN
                xi3_NC = epsilon3_NC
                xi3_CN = epsilon3_CN
                t1K = g0/g0
                t2KB = 0.0594_dp/g0
                t2KN = 0.2276_dp/g0
                t3K = -0.2163_dp/g0
            else ! just use GG
                call MIO_InputParameter('addPressureDependence',addPressureDependence,.true.)
                !call MIO_InputParameter('addPressureDependence',DCT,.true.)
                if (addPressureDependence) then
                    !if (DCT) then
                    !
                    !
                    !else
                    c1_0 = 0.310_dp
                    c1_1 = -1.882_dp
                    c1_2 = 7.741_dp
                    c2_0 = 1.750_dp
                    c2_1 = -1.618_dp
                    c2_2 = 1.848_dp
                    c3_0 = 1.990_dp
                    c3_1 = 1.007_dp
                    c3_2 = 2.427_dp
                    c4_0 = -0.068_dp
                    c4_1 = 0.399_dp
                    c4_2 = -1.739_dp
                    c5_0 = 3.286_dp
                    c5_1 = -0.914_dp
                    c5_2 = 12.011_dp
                    c6_0 = 0.5_dp
                    c6_1 = 0.322_dp
                    c6_2 = 0.908_dp
                    c7_0 = -0.008_dp
                    c7_1 = 0.046_dp
                    c7_2 = -0.183_dp
                    c8_0 = 2.272_dp
                    c8_1 = -0.721_dp
                    c8_2 = -4.414_dp
                    c9_0 = 1.217_dp
                    c9_1 = 0.027_dp
                    c9_2 = -0.658_dp
                    c10_0 = 1.562_dp
                    c10_1 = -0.371_dp
                    c10_2 = -0.134_dp
                    !end if
                    !if (F2G2Model) then
                    !   t1K = g0/g0
                    !   call MIO_InputParameter('Bilayert2KA',t2KA,-0.2235_dp)
                    !   call MIO_InputParameter('Bilayert2KB',t2KB,-0.2260_dp)
                    !   call MIO_InputParameter('Bilayert3K',t3K,0.1984_dp)
                    !   call MIO_InputParameter('Bilayert4K',t4K,0.0_dp)
                    !   call MIO_InputParameter('Bilayert5KA',t5KA,-0.04016_dp)
                    !   call MIO_InputParameter('Bilayert5KB',t5KB,-0.0404_dp)
                    !   call MIO_InputParameter('Bilayert6K',t6K,0.0_dp)
                    !   call MIO_InputParameter('Bilayert7K',t7K,0.0_dp)
                    !   call MIO_InputParameter('Bilayert8K',t8K,0.0_dp)
                    !   t2KA = t2KA/g0
                    !   t2KB = t2KB/g0
                    !   t3K = t3K/g0
                    !   t4K = t4K/g0
                    !   t5KA = t5KA/g0
                    !   t5KB = t5KB/g0
                    !   t6K = t6K/g0
                    !   t7K = t7K/g0
                    !   t8K = t8K/g0
                    !   t1K = g0/g0
                    !else
                    !   t2K = -0.2425_dp/g0
                    !   t3K = 0.2656_dp/g0
                    !   t4K = -0.0235_dp/g0
                    !   t5K = -0.0524_dp/g0
                    !   t6K = 0.0209_dp/g0
                    !   t7K = 0.0148_dp/g0
                    !   t8K = 0.0211_dp/g0
                    !   t2KA = t2K
                    !   t2KB = t2K
                    !   t5KA = t5K
                    !   t5KB = t5K
                    !end if
                else
                    lambda0 = 0.3155_dp/g0
                    epsilon0 = 1.7543_dp
                    x0 = 0.0_dp
                    kappa0 = 2.0010_dp
                    lambda3 = -0.0688_dp/g0
                    epsilon3 = 3.4692_dp
                    x3 = 0.5212_dp
                    kappa3 = 0.0_dp
                    lambda6 = -0.008300_dp/g0
                    epsilon6 = 2.876400_dp
                    x6 = 1.52060_dp
                    kappa6 = 1.57310_dp
                    xi0 = epsilon0
                    xi3 = epsilon3
                    xi6 = epsilon6
                    !print*, "lambda0 =", lambda0
                    !print*, "lambda3 =", lambda3
                    !print*, "lambda6 =", lambda6
                    !print*, "kappa0 =", kappa0
                    !print*, "kappa3 =", kappa3
                    !print*, "kappa6 =", kappa6
                    !print*, "xi0 =", xi0
                    !print*, "xi3 =", xi3
                    !print*, "xi6 =", xi6
                    !print*, "x0 =", x0
                    !print*, "x3 =", x3
                    !print*, "x6 =", 6
               end if
            end if
         else if (MIO_StringComp(BilayerModel,'BLSrivani')) then
            call MIO_InputParameter('useBNGSrivani',useBNGSrivani,.false.)
            call MIO_InputParameter('oppositedxdy',oppositedxdy,.false.)
            call MIO_InputParameter('sublatticeDependent',sublatticeDependent,.false.)
            call MIO_InputParameter('useThetaIJ',useThetaIJ,.false.)
            call MIO_InputParameter('useTheta',useTheta,.false.)
            call MIO_InputParameter('sublatticeIndependent',sublatticeIndependent,.false.)
            if (useBNGSrivani) then
                lambda0_CB = 0.3905_dp/g0
                epsilon0_CB = 1.5426_dp
                !x0_CB = 0.0_dp
                kappa0_CB = 1.8229_dp
                lambda3_BC = -0.0588_dp/g0
                epsilon3_BC = 3.0827_dp
                x3_BC = 0.6085_dp
                !kappa3_BC = 0.0_dp
                lambda3_CB = -0.0651_dp/g0
                epsilon3_CB = 3.7998_dp
                x3_CB = 0.6341_dp
                !kappa3_CB = 0.0_dp
                xi0_CB = epsilon0_CB
                xi3_BC = epsilon3_BC
                xi3_CB = epsilon3_CB
                lambda0_CN = 0.2517_dp/g0
                epsilon0_CN = 1.6061_dp
                !x0_CN = 0.0_dp
                kappa0_CN = 2.1909_dp
                lambda3_NC = -0.0606_dp/g0
                epsilon3_NC = 3.3502_dp
                x3_NC = 0.5142_dp
                !kappa3_NC = 0.0_dp
                lambda3_CN = -0.0465_dp/g0
                epsilon3_CN = 3.0464_dp
                x3_CN = 0.5264_dp
                !kappa3_CN = 0.0_dp
                xi0_CN = epsilon0_CN
                xi3_NC = epsilon3_NC
                xi3_CN = epsilon3_CN
                t1K = g0/g0
                t2KB = 0.0594_dp/g0
                t2KN = 0.2276_dp/g0
                t3K = -0.2163_dp/g0
            else
                call MIO_InputParameter('addPressureDependence',addPressureDependence,.true.)
                call MIO_InputParameter('oldParameterSet',oldParameterSet,.false.)
                call MIO_InputParameter('useOnlyVAB',useOnlyVAB,.false.)
                if (addPressureDependence) then
                    call MIO_Print('Adding pressure (distance) dependent Srivani parameters','ham')
                    if (oldParameterSet) then
                        tAB = 0.355_dp
                        c1_0 = 0.3571385063263838_dp
                        c1_1 = -1.891884708936455_dp
                        c1_2 = 6.652946573261291_dp
                        c2_0 = 1.8867190568791063_dp
                        c2_1 = -0.41000467960493125_dp
                        c2_2 = 4.578192878264396_dp
                        c3_0 = 1.8270139657989923_dp
                        c3_1 = -0.40382997041549323_dp
                        c3_2 = 0.8339980356453551_dp
                        c4_0 = 0.07516409521688965_dp
                        c4_1 = -0.4003623778798789_dp
                        c4_2 = 1.8994745927765773_dp
                        c5_0 = 3.6275668776932046_dp
                        c5_1 = -3.6911506279524273_dp
                        c5_2 = -10.77882608119433_dp
                        c6_0 = 0.5315846479951798_dp
                        c6_1 = 0.005642504054093722_dp
                        c6_2 = -2.1672957606642598_dp
                        c7_0 = -0.009733220417244328_dp
                        c7_1 = 0.0508999779194876_dp
                        c7_2 = -0.1994906574479795_dp
                        c8_0 = 2.694749513617616_dp
                        c8_1 = -0.7312487102958507_dp
                        c8_2 = 8.375059780569929_dp
                        c9_0 = 1.5483753173366424_dp
                        c9_1 = 0.1578717903624763_dp
                        c9_2 = -0.4478833970593441_dp
                        c10_0 =1.5906430728955332_dp
                        c10_1 =-0.26841222887036603_dp
                        c10_2 =0.6153669604354582_dp
                    else
                        !c1_0 = 0.3360066919                          2270793_dp
                        !c1_1 = -1.929949468                          9057044_dp
                        !c1_2 = 5.9426366687                          61513_dp
                        !c2_0 = 1.8020183122                          69422_dp
                        !c2_1 = -1.297300578                          7982854_dp
                        !c2_2 = 2.2929720244                          985514_dp
                        !c3_0 = 1.7585214228                          921264_dp
                        !c3_1 = -0.738930236                          5991159_dp
                        !c3_2 = 0.7382964144                          771087_dp
                        !c4_0 = 0.0684965252                          4022183_dp
                        !c4_1 = -0.441876050                          9738568_dp
                        !c4_2 = 0.6355652920                          68162_dp
                        !c5_0 = 3.1765595945                          8736_dp
                        !c5_1 = -8.668955264                          395734_dp
                        !c5_2 = -4.565706935                          115334_dp
                        !c6_0 = 0.39776800205933815_dp
                        !c6_1 = -0.9045887266045054_dp
                        !c6_2 = -1.3002086877459056_dp
                        !c7_0 = 0.011706116712585913_dp
                        !c7_1 = -0.057290384984340935_dp
                        !c7_2 = 1.1876301867657304_dp
                        !c8_0 = 1.2665913145651069_dp
                        !c8_1 = 1.4879363006310338_dp
                        !c8_2 = -5.7011812463480185_dp
                        !c9_0 = 1.003729197086963_dp
                        !c9_1 = -0.7170670900835974_dp
                        !c9_2 = -9.969709323763947_dp
                        !c10_0 =-2.1171606094526427_dp
                        !c10_1 =-1.1346126431332875_dp
                        !c10_2 =-6.1966938524073_dp
                        !c11_0 =0.030521037568468558_dp
                        !c11_1 =-0.051532907774878865_dp
                        !c11_2 =2.185952444870147_dp
                        !c12_0 =0.4604653597397456_dp
                        !c12_1 =0.08328035937371601_dp
                        !c12_2 =6.600227553330765_dp
                        !c13_0 =-0.5619164303945731_dp
                        !c13_1 =-2.497339711639136_dp
                        !c13_2 =8.65566068091937_dp
                        !c14_0 =2.190603782512041_dp
                        !c14_1 =0.6273043778418543_dp
                        !c14_2 =1.2043001956176094_dp
                        tAB = 0.3357_dp
                        !c1_0 = 0.4332707088245484
                        !c1_1 = -2.325309314458914
                        !c1_2 = 3.1914228801054376
                        !c2_0 = 1.8536617705776488
                        !c2_1 = -0.1117253191783712
                        !c2_2 = -24.766546630000487
                        !c3_0 = 1.7913844513064086
                        !c3_1 = -0.5099517871187137
                        !c3_2 = -5.375985788072558
                        !c4_0 = 0.10425992611308767
                        !c4_1 = -1.045722730568572
                        !c4_2 = 5.735298412084255
                        !c5_0 = 3.292881484137035
                        !c5_1 = 2.661330783775027
                        !c5_2 = -128.28901878128022
                        !c6_0 = 0.3726057881342914
                        !c6_1 = 1.44600764984851
                        !c6_2 = -21.52615676392712
                        !c7_0 = -0.019074365151675893
                        !c7_1 = 0.22125761262807156
                        !c7_2 = -1.9730605318579442
                        !c8_0 = 1.1784542997555516
                        !c8_1 = 2.3669704479327196
                        !c8_2 = -27.825121319945847
                        !c9_0 = 0.9460232858013488
                        !c9_1 = 0.6209696680905099
                        !c9_2 = -10.117818809881792
                        !c10_0 =2.0747368801425905
                        !c10_1 =0.020378997950056467
                        !c10_2 =15.464068143926946
                        !c11_0 =0.05706633732282606
                        !c11_1 =-0.7001457195512684
                        !c11_2 =6.415621197423687
                        !c12_0 =0.47163687834310253
                        !c12_1 =0.013383162167131547
                        !c12_2 =-3.706093128922178e-06
                        !c13_0 =-0.6175549975628934
                        !c13_1 =-0.2124493950429618
                        !c13_2 =7.404171266791098
                        !c14_0 =2.1700121321486385
                        !c14_1 =-0.013403161485675559
                        !c14_2 =-1.0604054159972579e-05
                        if (useOnlyVAB) then
                            c1_0 =   0.3356095024843812
                            c1_1 =   -1.9555147793478367
                            c1_2 =   5.542631379367501
                            c2_0 =   1.718698822030771
                            c2_1 =   -1.989378866684143
                            c2_2 =   6.0025023827468775
                            c3_0 =   1.7543868745706959
                            c3_1 =   -0.7739949690181664
                            c3_2 =   1.903243378532328
                            c4_0 =   0.06851906419484396
                            c4_1 =   -0.6872433666263775
                            c4_2 =   3.769207297230312
                            c5_0 =   3.1733357766403234
                            c5_1 =   2.5647115615881324
                            c5_2 =   -123.63154378822028
                            c6_0 =   0.39765806464770465
                            c6_1 =   1.543230705089051
                            c6_2 =   -22.973479855570684
                            c7_0 =   -0.012256483205855817
                            c7_1 =   0.16990943196521704
                            c7_2 =   -2.1647529202020253
                            c8_0 =   1.1716245804866041
                            c8_1 =   4.487359824823402
                            c8_2 =   -58.75949608817599
                            c9_0 =   0.8917021392108377
                            c9_1 =   1.0187648341101367
                            c9_2 =   -16.229901468075237
                            c10_0 =  2.09879167259061
                            c10_1 =  0.29714623365951714
                            c10_2 =  12.863626996121283
                            c11_0 =  0.04071154326934091
                            c11_1 =  -0.5873777395380401
                            c11_2 =  7.654992092022496
                            c12_0 =  0.4488364197392876
                            c12_1 =  0.013400104807997755
                            c12_2 =  -2.142393269209168e-06
                            c13_0 =  -0.7078365546490298
                            c13_1 =  -0.06723178248496989
                            c13_2 =  2.0037377859863668
                            c14_0 =  2.109173635556271
                            c14_1 =  -0.14140771101700622
                            c14_2 =  7.571534745370173
                        else
                            c1_0 =  0.3575872537033165_dp
                            c1_1 =  -1.9347697178845817_dp
                            c1_2 =  3.8706979283760394_dp
                            c2_0 =  1.9106957393116117_dp
                            c2_1 =  -0.09989684567144147_dp
                            c2_2 =  -8.720932407008418_dp
                            c3_0 =  1.8228523809344535_dp
                            c3_1 =  -0.3667529691024049_dp
                            c3_2 =  -2.3050530340090782_dp
                            c4_0 =  0.071115719665038_dp
                            c4_0 = -c4_0
                            !print*, "changing the sign of V3 to match with Kaxiras model"
                            c4_1 =  -0.42093264500073585_dp
                            c4_2 =  0.9178867554728614_dp
                            c5_0 =  3.9416508464100466_dp
                            c5_1 =  -2.214147419045211_dp
                            c5_2 =  41.76842643035146_dp
                            c6_0 =  0.5481096417108225_dp
                            c6_1 =  0.04235827800799908_dp
                            c6_2 =  3.973507493054687_dp
                            c7_0 =  -0.031087946815144728_dp
                            c7_1 =  0.19481521423521497_dp
                            c7_2 =  -0.7770045478761431_dp
                            c8_0 =  1.1147111796813856_dp
                            c8_1 =  -0.013399935958530826_dp
                            c8_2 =  -1.4302593524023142e-06_dp
                            c9_0 =  0.9353113982829261_dp
                            c9_1 =  0.20106572786860608_dp
                            c9_2 =  -1.13259450482717_dp
                            c10_0 = 1.3970867695033276_dp
                            c10_1 = 0.032853758467485825_dp
                            c10_2 = -0.6913705164987503_dp
                            c11_0 = 0.08985168760341603_dp
                            c11_1 = -0.5898948298140482_dp
                            c11_2 = 2.5050264708064574_dp
                            c12_0 = 0.49260894402823696_dp
                            c12_1 = 0.013400115682404422_dp
                            c12_2 = -2.3561065136309557e-06_dp
                            c13_0 = -0.5744255984608703_dp
                            c13_1 = 0.16416466983562875_dp
                            c13_2 = -0.00014391022969531333_dp
                            c14_0 = 1.3345101551849814_dp
                            c14_1 = -0.01942463416793111_dp
                            c14_2 = 0.22425121204108625_dp
                            !c1_0 = 0.33559774158126227
                            !c1_1 = -1.8011108015126933
                            !c1_2 = 2.4719748458020905
                            !c2_0 = 1.799761302761298
                            !c2_1 = -0.108474483023478
                            !c2_2 = -24.046432757807015
                            !c3_0 = 1.7581554112735707
                            !c3_1 = -0.5004921970060755
                            !c3_2 = -5.2762717963848145
                            !c4_0 = 0.06851906419484396
                            !c4_1 = -0.6872433666263775
                            !c4_2 = 3.769207297230312
                            !c5_0 = 3.1733357766403234
                            !c5_1 = 2.5647115615881324
                            !c5_2 = -123.63154378822028
                            !c6_0 = 0.39765806464770465
                            !c6_1 = 1.543230705089051
                            !c6_2 = -22.973479855570684
                            !c7_0 = -0.013084941395536264
                            !c7_1 = 0.1517818740226915
                            !c7_2 = -1.3535119564196156
                            !c8_0 = 1.2284838218753473
                            !c8_1 = 2.467456656444414
                            !c8_2 = -29.00639547379991
                            !c9_0 = 0.9544482537082066
                            !c9_1 = 0.6265000202261579
                            !c9_2 = -10.207929682758467
                            !c10_0 =2.1056269285153353
                            !c10_1 =0.02069373935912141
                            !c10_2 =15.694040617299631
                            !c11_0 =0.038429929746244776
                            !c11_1 =-0.47149601494384413
                            !c11_2 =4.320443175164844
                            !c12_0 =0.47223611698000756
                            !c12_1 =0.0134001659317156
                            !c12_2 =-3.705808808521222e-06
                            !c13_0 =-0.6115664635903063
                            !c13_1 =-0.2103895430576807
                            !c13_2 =7.332377481695472
                            !c14_0 =2.1694119672988874
                            !c14_1 =-0.013399912532871622
                            !c14_2 =-2.535705314295427e-06
                            !print*, "checking the double precision: ", c14_2, c14_1
                        end if
                    end if
                    !if (F2G2Model) then
                    !   t1K = g0/g0
                    !   call MIO_InputParameter('Bilayert2KA',t2KA,-0.2235_dp)
                    !   call MIO_InputParameter('Bilayert2KB',t2KB,-0.2260_dp)
                    !   call MIO_InputParameter('Bilayert3K',t3K,0.1984_dp)
                    !   call MIO_InputParameter('Bilayert4K',t4K,0.0_dp)
                    !   call MIO_InputParameter('Bilayert5KA',t5KA,-0.04016_dp)
                    !   call MIO_InputParameter('Bilayert5KB',t5KB,-0.0404_dp)
                    !   call MIO_InputParameter('BilayertK6A',t6K,0.0_dp)
                    !   call MIO_InputParameter('BilayertK7A',t7K,0.0_dp)
                    !   call MIO_InputParameter('BilayertK8A',t8K,0.0_dp)
                    !   t2KA = t2KA/g0
                    !   t2KB = t2KB/g0
                    !   t3K = t3K/g0
                    !   t4K = t4K/g0
                    !   t5KA = t5KA/g0
                    !   t5KB = t5KB/g0
                    !   t6K = t6K/g0
                    !   t7K = t7K/g0
                    !   t8K = t8K/g0
                    !else
                    !   t1K = g0/g0
                    !   t2K = -0.2425_dp/g0
                    !   t3K = 0.2656_dp/g0
                    !   t4K = -0.0235_dp/g0
                    !   t5K = -0.0524_dp/g0
                    !   t6K = 0.0209_dp/g0
                    !   t7K = 0.0148_dp/g0
                    !   t8K = 0.0211_dp/g0
                    !   t2KA = t2K
                    !   t2KB = t2K
                    !   t5KA = t5K
                    !   t5KB = t5K
                    !end if
                else !  New parameters from 2020 by Srivani
                    !print*, "Using the new parameters from 2020 for 3.35 interlayer distance, V6 expression probably wrong"
                    lambda0 =  0.3575896463595361_dp/g0
                    epsilon0 = 1.9108224313643323_dp
                    x0 = 0.0_dp
                    kappa0 = 1.8227815545297545_dp
                    lambda3 = 0.07111484887484568_dp/g0
                    epsilon3 = 3.9428998931290584_dp
                    x3 = 0.5482006335511792_dp
                    kappa3 = 0.0_dp
                    lambda6 = -0.013721066052788262_dp/g0
                    epsilon6 = 1.25349967_dp
                    x6 = 1.09324627_dp
                    kappa6 = 1.7524165109575078_dp
                    lambda6b = 0.029028452050848814_dp/g0
                    epsilon6b = 0.4749792471502847_dp
                    x6b = 1.9276014_dp
                    kappa6b = -0.5685377300000001_dp
                    !lambda0 = 0.3155_dp/g0
                    !epsilon0 = 1.7543_dp
                    !x0 = 0.0_dp
                    !kappa0 = 2.0010_dp
                    !lambda3 = -0.0688_dp/g0
                    !epsilon3 = 3.4692_dp
                    !x3 = 0.5212_dp
                    !kappa3 = 0.0_dp
                    !lambda6 = -0.008300_dp/g0
                    !epsilon6 = 2.876400_dp
                    !x6 = 1.52060_dp
                    !kappa6 = 1.57310_dp
                    !xi0 = epsilon0
                    !xi3 = epsilon3
                    !xi6 = epsilon6
                    xi0 = epsilon0
                    xi3 = epsilon3
                    xi6 = epsilon6
                    xi6b = epsilon6b
                    !if (sublatticeDependent) then
                    !   print*, "Actually, lets use the sublattice dependent parameters"
                       lambda0AAp     = 0.3777452105728493_dp/g0
                       xi0AAp         =  1.673723757628059_dp
                       kappa0AAp      = -1.9375254860531055_dp
                       lambda0bAAp     = -0.03297124774024884_dp/g0
                       xi0bAAp         =  1.0731649983324736_dp
                       kappa0bAAp      = 2.1774201575458982_dp

                       lambda0ABp     = 0.3339824357142346_dp/g0
                       xi0ABp         =  1.6086712193494612_dp
                       kappa0ABp      = -1.8083599334947213_dp
                       lambda0bABp     = -0.017946450497759416_dp/g0
                       xi0bABp         =  0.798265793648077_dp
                       kappa0bABp      = 2.093789675583678_dp

                       lambda0BAp     = 0.3339824357142346_dp/g0
                       xi0BAp         =  1.6086712193494612_dp
                       kappa0BAp      = -1.8083599334947213_dp
                       lambda0bBAp     = -0.017946450497759416_dp/g0
                       xi0bBAp         =  0.798265793648077_dp
                       kappa0bBAp      = 2.093789675583678_dp

                       lambda3ABp     =  0.22358533630187713_dp/g0
                       xi3ABp         =  2.515609650785262_dp
                       x3ABp          =  0.6615199942493694_dp
                       lambda3bABp     =  -0.13763943627632735_dp/g0
                       xi3bABp         =  2.5155817180683107_dp
                       x3bABp        =  0.6618620509164502_dp

                       lambda3BAp     =  0.22358533630187713_dp/g0
                       xi3BAp         =  2.515609650785262_dp
                       x3BAp          =  0.6615199942493694_dp
                       lambda3bBAp     =  -0.13763943627632735_dp/g0
                       xi3bBAp         =  2.5155817180683107_dp
                       x3bBAp        =  0.6618620509164502_dp
                    !end if
                    !lambda0 = 0.3155_dp/g0
                    !epsilon0 = 1.7543_dp
                    !x0 = 0.0_dp
                    !kappa0 = 2.0010_dp
                    !lambda3 = -0.0688_dp/g0
                    !epsilon3 = 3.4692_dp
                    !x3 = 0.5212_dp
                    !kappa3 = 0.0_dp
                    !lambda6 = -0.008300_dp/g0
                    !epsilon6 = 2.876400_dp
                    !x6 = 1.52060_dp
                    !kappa6 = 1.57310_dp
                    !xi0 = epsilon0
                    !xi3 = epsilon3
                    !xi6 = epsilon6
                    !if (F2G2Model) then
                    !   t1K = g0/g0
                    !   call MIO_InputParameter('Bilayert2KA',t2KA,-0.2235_dp)
                    !   call MIO_InputParameter('Bilayert2KB',t2KB,-0.2260_dp)
                    !   call MIO_InputParameter('Bilayert3K',t3K,0.1984_dp)
                    !   call MIO_InputParameter('Bilayert4K',t4K,0.0_dp)
                    !   call MIO_InputParameter('Bilayert5KA',t5KA,-0.04016_dp)
                    !   call MIO_InputParameter('Bilayert5KB',t5KB,-0.0404_dp)
                    !   call MIO_InputParameter('BilayertK6A',t6K,0.0_dp)
                    !   call MIO_InputParameter('BilayertK7A',t7K,0.0_dp)
                    !   call MIO_InputParameter('BilayertK8A',t8K,0.0_dp)
                    !   t2KA = t2KA/g0
                    !   t2KB = t2KB/g0
                    !   t3K = t3K/g0
                    !   t4K = t4K/g0
                    !   t5KA = t5KA/g0
                    !   t5KB = t5KB/g0
                    !   t6K = t6K/g0
                    !   t7K = t7K/g0
                    !   t8K = t8K/g0
                    !else
                    !   t1K = g0/g0
                    !   t2K = -0.2425_dp/g0
                    !   t3K = 0.2656_dp/g0
                    !   t4K = -0.0235_dp/g0
                    !   t5K = -0.0524_dp/g0
                    !   t6K = 0.0209_dp/g0
                    !   t7K = 0.0148_dp/g0
                    !   t8K = 0.0211_dp/g0
                    !   t2KA = t2K
                    !   t2KB = t2K
                    !   t5KA = t5K
                    !   t5KB = t5K
                    !end if
                    !print*, "lambda0 =", lambda0
                    !print*, "lambda3 =", lambda3
                    !print*, "lambda6 =", lambda6
                    !print*, "kappa0 =", kappa0
                    !print*, "kappa3 =", kappa3
                    !print*, "kappa6 =", kappa6
                    !print*, "xi0 =", xi0
                    !print*, "xi3 =", xi3
                    !print*, "xi6 =", xi6
                    !print*, "x0 =", x0
                    !print*, "x3 =", x3
                    !print*, "x6 =", x6
               end if
            end if
         end if
         call MIO_InputParameter('Neigh.LayerDistFactor',distFact,1.0_dp)
         call MIO_InputParameter('KaxirasCutoff',KaxirasCutoff,1.0_dp)
         !KaxirasCutoff = 1.5_dp
         KaxirasCutoff2 = (KaxirasCutoff)**2.0_dp
         !print*, "Kaxiras Cutoff = ", KaxirasCutoff
         !SrivaniCutoff = (aG/sqrt(3.0_dp)*(distFact-1.2_dp))**2.0_dp

         call MIO_InputParameter('MoireOffDiag',l,.false.)
         call MIO_InputParameter('GBNOffDiag',GBNOffDiag,.false.)
         call MIO_InputParameter('addDisplacements',addDisplacements,.false.)

         call MIO_InputParameter('singleLayerXYZ',singleLayerXYZ,.false.)
         call MIO_InputParameter('addPressureDependence',addPressureDependence,.true.)
         call MIO_InputParameter('switchV3Sign',switchV3Sign,.false.)
         call MIO_InputParameter('oldParameterSet',oldParameterSet,.false.)
         call MIO_InputParameter('deactivateV6',deactivateV6,.false.)
         call MIO_InputParameter('deactivateV3',deactivateV3,.false.)
         call MIO_InputParameter('KoshinoIntralayer',KoshinoIntralayer,.false.)
         call MIO_InputParameter('MayouIntralayer',MayouIntralayer,.false.)
         call MIO_InputParameter('useBNGKaxiras',useBNGKaxiras,.false.)
         call MIO_InputParameter('twoLayers',twoLayers,.false.)
         call MIO_InputParameter('oneLayer',oneLayer,.false.)
         call MIO_InputParameter('GBNtwoLayers',GBNtwoLayers,.false.)
         call MIO_InputParameter('BNBNtwoLayers',BNBNtwoLayers,.false.)
         call MIO_InputParameter('encapsulatedThreeLayers',encapsulatedThreeLayers,.false.)
         call MIO_InputParameter('removeTopMoireInL2',removeTopMoireInL2,.false.)
         call MIO_InputParameter('encapsulatedFourLayers',encapsulatedFourLayers,.false.)
         call MIO_InputParameter('encapsulatedFiveLayers',encapsulatedFiveLayers,.false.)
         call MIO_InputParameter('encapsulatedSixLayers',encapsulatedSixLayers,.false.)
         call MIO_InputParameter('encapsulatedSevenLayers',encapsulatedSevenLayers,.false.)
         call MIO_InputParameter('t3GwithBN',t3GwithBN,.false.)
         call MIO_InputParameter('BNt2GBN',BNt2GBN,.false.)
         call MIO_InputParameter('t2GBN',t2GBN,.false.)
         call MIO_InputParameter('encapsulatedFourLayersF2G2',encapsulatedFourLayersF2G2,.false.)
         call MIO_InputParameter('GBNuseDisplacementFile',GBNuseDisplacementFile,.false.)
         call MIO_InputParameter('tBGuseDisplacementFile',tBGuseDisplacementFile,.false.)
         call MIO_InputParameter('BNBNtwoLayers',BNBNtwoLayers,.false.)
         call MIO_InputParameter('onlyV0',onlyV0,.false.)
         call MIO_InputParameter('GBNAngle',GBNAngle,0.0_dp)
         GBNAngle = GBNAngle*pi/180.0_dp
         call MIO_InputParameter('tBGAngle',tBGAngle,0.0_dp)
         tBGAngle = tBGAngle*pi/180.0_dp
         call MIO_InputParameter('rotationAngle',rotationAngle,0.0_dp)
         rotationAngle = rotationAngle*pi/180.0_dp
         call MIO_InputParameter('GlobalTwist',Globalang,0.0_dp)
         call MIO_InputParameter('GlobalTwist2',Globalang2,0.0_dp)
         call MIO_InputParameter('GlobalPhiL1',GlobalPhiL1,0.0_dp)
         call MIO_InputParameter('GlobalPhiL2',GlobalPhiL2,0.0_dp)
         call MIO_InputParameter('GlobalPhiL2a',GlobalPhiL2a,0.0_dp)
         call MIO_InputParameter('GlobalPhiL2b',GlobalPhiL2b,0.0_dp)
         call MIO_InputParameter('GlobalPhiL3',GlobalPhiL3,0.0_dp)
         epstBG = 0.0_dp
         !if (tBGOffDiag) then
         !    rotationAngle = rotationAngle*pi/180.0_dp
         !    epstBG = 0.0_dp
         !    GlobalPhi = atan((1.0_dp+epstBG)*sin(rotationAngle)/((1.0_dp+epstBG)*cos(rotationAngle)-1.0_dp))
         !else
         !    GBNAngle = GBNAngle*pi/180.0_dp
             GlobalPhi = Globalang*pi/180.0_dp
             GlobalPhi2 = Globalang2*pi/180.0_dp
             GlobalPhiL1 = GlobalPhiL1*pi/180.0_dp
             GlobalPhiL2 = GlobalPhiL2*pi/180.0_dp
             GlobalPhiL2a = GlobalPhiL2a*pi/180.0_dp
             GlobalPhiL2b = GlobalPhiL2b*pi/180.0_dp
             GlobalPhiL3 = GlobalPhiL3*pi/180.0_dp
         !end if
         call MIO_InputParameter('threeLayers',threeLayers,.false.)
         call MIO_InputParameter('fourLayersSandwiched',fourLayersSandwiched,.false.)
         call MIO_InputParameter('helicalTwistedMBM',helicalTwistedMBM,.false.)
         call MIO_InputParameter('fiveLayersSandwiched',fiveLayersSandwiched,.false.)
         call MIO_InputParameter('sixLayersSandwiched',sixLayersSandwiched,.false.)
         call MIO_InputParameter('sevenLayersSandwiched',sevenLayersSandwiched,.false.)
         call MIO_InputParameter('eightLayersSandwiched',eightLayersSandwiched,.false.)
         call MIO_InputParameter('tenLayersSandwiched',tenLayersSandwiched,.false.)
         call MIO_InputParameter('twentyLayersSandwiched',twentyLayersSandwiched,.false.)
         call MIO_InputParameter('middleTwist',middleTwist,.false.)
         call MIO_InputParameter('fourLayers',fourLayers,.false.)
         call MIO_InputParameter('findThetasGeometrically',findThetasGeometrically,.false.)
         !call MIO_InputParameter('renormalizeHoppings',renormalizeHoppings,.false.)
         call MIO_InputParameter('renormalizeHoppingFactorAAp',renormalizeHoppingFactorAAp,1.0_dp)
         call MIO_InputParameter('renormalizeHoppingFactorBBp',renormalizeHoppingFactorBBp,1.0_dp)
         call MIO_InputParameter('renormalizeHoppingFactorABp',renormalizeHoppingFactorABp,1.0_dp)
         call MIO_InputParameter('renormalizeHoppingFactorBAp',renormalizeHoppingFactorBAp,1.0_dp)
         call MIO_InputParameter('newFittingFunctions',newFittingFunctions,.false.)
         call MIO_InputParameter('renormalizeCoupling',renormalizeCoupling,.false.)
         call MIO_InputParameter('couplingFactor',couplingFactor,1.0_dp)
         call MIO_InputParameter('couplingFactor2',couplingFactor2,1.0_dp)
         call MIO_InputParameter('differentCouplings',differentCouplings,.false.)
         call MIO_InputParameter('readRigidXYZ',readRigidXYZ,.false.)
         !call MIO_Print('By default, we replace the Koshino Intralayer terms with the F2G2 model, make sure this is what you want','ham')

         call MIO_InputParameter('KoshinoSR',KoshinoSR,.false.)
         call MIO_InputParameter('onlyvppsigma',onlyvppsigma,.false.)

         call MIO_InputParameter('Bulk',bulk,.false.)
         call MIO_InputParameter('BernalReadXYZ',BernalReadXYZ,.false.)
         call MIO_InputParameter('MoireSecondMoireMassFactor',sign2,1.0_dp)
         call MIO_InputParameter('deactivateIntrasublattice',deactivateIntrasublattice,.false.)
         call MIO_InputParameter('deactivateIntersublattice',deactivateIntersublattice,.false.)
         call MIO_InputParameter('deactivateIntraSublatticeForC',deactivateIntraSublatticeForC,.false.)
         call MIO_InputParameter('WriteDataFiles',writeData,.false.)
         call MIO_InputParameter('t2Value',t2Value,0.0083_dp)
         ! f1
         call MIO_InputParameter('tA1B2',tA1B2,-0.1418_dp)
         tA1B2 = tA1B2/g0
         call MIO_InputParameter('tA1A2',tA1A2,-0.0883_dp)
         tA1A2 = tA1A2/g0
         call MIO_InputParameter('tB1B2',tB1B2,-0.08965_dp)
         tB1B2 = tB1B2/g0
         call MIO_InputParameter('tA1A3',tA1A3,-0.00541_dp)
         tA1A3 = tA1A3/g0
         call MIO_InputParameter('tB1A3',tB1A3,-0.00853_dp)
         tB1A3 = tB1A3/g0
         ! f2
         call MIO_InputParameter('tA1B2_2',tA1B2_2,0.0755_dp)
         tA1B2_2 = tA1B2_2/g0
         call MIO_InputParameter('tA1A2_2',tA1A2_2,0.0277_dp)
         tA1A2_2 = tA1A2_2/g0
         call MIO_InputParameter('tB1B2_2',tB1B2_2,0.0275_dp)
         tB1B2_2 = tB1B2_2/g0
         call MIO_InputParameter('tA1A3_2',tA1A3_2,0.0009_dp)
         tA1A3_2 = tA1A3_2/g0
         call MIO_InputParameter('tB1A3_2',tB1A3_2,0.0028_dp)
         tB1A3_2 = tB1A3_2/g0
         ! g0 ! other terms are onsite energy terms or the t2value term
         call MIO_InputParameter('tA1B3',tA1B3,-0.007551_dp)
         tA1B3 = tA1B3/g0
         call MIO_InputParameter('tB1A2',tB1A2,-0.3372_dp)
         tB1A2 = tB1A2/g0
         ! g1
         call MIO_InputParameter('tA1A1_1',tA1A1_1,-0.2238_dp)
         tA1A1_1 = tA1A1_1/g0
         call MIO_InputParameter('tB1B1_1',tB1B1_1,-0.2269_dp)
         tB1B1_1 = tB1B1_1/g0
         call MIO_InputParameter('tA2A2_1',tA2A2_1,-0.2264_dp)
         tA2A2_1 = tA2A2_1/g0
         call MIO_InputParameter('tA1B3_1',tA1B3_1,-0.001738_dp)
         tA1B3_1 = tA1B3_1/g0
         call MIO_InputParameter('tB1A2_1',tB1A2_1,0.00786_dp)
         tB1A2_1 = tB1A2_1/g0
         ! g2
         call MIO_InputParameter('tA1A1_2',tA1A1_2,-0.0446_dp)
         tA1A1_2 = tA1A1_2/g0
         call MIO_InputParameter('tB1B1_2',tB1B1_2,-0.04606_dp)
         tB1B1_2 = tB1B1_2/g0
         call MIO_InputParameter('tA2A2_2',tA2A2_2,-0.044_dp)
         tA2A2_2 = tA2A2_2/g0
         call MIO_InputParameter('tA1B3_2',tA1B3_2,0.001772_dp)
         tA1B3_2 = tA1B3_2/g0
         call MIO_InputParameter('tB1A2_2',tB1A2_2,0.00078_dp)
         tB1A2_2 = tB1A2_2/g0

         if (bulk) then
            if (twoLayers) then
                nLayers = 2
            else if (threeLayers) then
                nLayers = 3
            else if (fourLayersSandwiched .or. helicalTwistedMBM) then
                nLayers = 4
            else if (fiveLayersSandwiched) then
                nLayers = 5
            else if (sixLayersSandwiched) then
                nLayers = 6
            else if (sevenLayersSandwiched) then
                nLayers = 7
            else if (eightLayersSandwiched) then
                nLayers = 8
            else if (tenLayersSandwiched) then
                nLayers = 10
            else if (twentyLayersSandwiched) then
                nLayers = 20
            else
               print*, "is this a new type of bulk system?"
            end if
            call MIO_Print('We set the number of layers for the bulk','ham')
         else
            nLayers = 0
         end if
         call MIO_InputParameter('setnLayersToZero',setnLayersToZero,.false.)
         if (setnLayersToZero) then ! useful when using the small bulk neighbor but system is not bulk
            nLayers = 0
         end if

         !aval = 1.34_dp
         !bval = 3.25_dp
         !avalC1B = 1.34_dp
         !bvalC2B = 3.25_dp
         !avalC1N = 1.34_dp
         !bvalC2N = 3.25_dp

         countG1  = 0
         countG2  = 0
         countG3  = 0
         countG4  = 0
         countG5  = 0
         countG6  = 0
         countG7  = 0
         countG8  = 0
         countBN1 = 0
         countBN2 = 0
         countBN3 = 0
         countBN4 = 0
         countBN5 = 0
         countBN6 = 0
         countBN7 = 0
         countBN8 = 0

         call MIO_InputParameter('MoirePotCabG',CabG_global,0.002235_dp)
         call MIO_InputParameter('tBGSwitchDxDy',tBGSwitchDxDy,.false.)

         call MIO_Allocate(HABreal,[inode1],[inode2],'H0','ham')
         call MIO_Allocate(HABimag,[inode1],[inode2],'H0','ham')
         if (frac) call AtomsSetCart()

         ! --- Local p_z orbital axes for the optional corrugated interlayer -----
         ! two-center model. pzNormal is allocated and zeroed unconditionally so
         ! interlayerTwoCenter() can be called branch-free at every site; it is
         ! only populated (owned atoms + MPI halo) when the flag is set.
#ifdef MPI
         pzNGhost = 0
         if (associated(rcvList)) pzNGhost = size(rcvList)
#else
         pzNGhost = 0
#endif /* MPI */
         call MIO_Allocate(pzNormal,[1,inode1],[3,inode2+pzNGhost],'pzNormal','ham')
         pzNormal = 0.0_dp
         pzNfallback = 0
         if (corrugatedInterlayerTwoCenter) then
            call MIO_Print('Interlayer two-center: local p_z normals from same-layer first-shell neighbours','ham')
            call ComputePzNormals(pzNormal, pzNfallback)
            call HamUpdateNormals(pzNormal)
#ifdef DEBUG
            call HamPrintNormalStats(pzNormal, pzNfallback)
#endif /* DEBUG */
         end if
         ! ---------------------------------------------------------------------
         !print*, "WARNING WARNING, remove next line, just temporary for now"
         !call AtomsSetFrac() ! REMOVE THIS, jut temporary
         !!!$OMP& PRIVATE (expFactor, dist, yDistance, DBShiftRatio,boundaryWidth, limit2, limit1, limit0), &
         !!!$OMP& PRIVATE (nnnn), &
         !$OMP PARALLEL DO PRIVATE(d,nlay,delta,del,realH,imagH,Habjj,dx,dy,dxTemp,dyTemp,HBL,HAA,HAB,HBA,dxi,dxj), &
         !$OMP& PRIVATE(realH_b,imagH_b,Habjj_b,dx_b,dy_b,dxTemp_b,dyTemp_b), &
         !$OMP& PRIVATE(realH_t,imagH_t,Habjj_t,dx_t,dy_t,dxTemp_t,dyTemp_t), &
         !$OMP& PRIVATE (dyi,dyj,vpppi,vppsigma,posOrNeg,jj,jj1,jj2,maxDist,k,firstNN), &
         !$OMP& PRIVATE (rbar, V0, V3, V6, theta, d122,d232, d132, d232F), &
         !$OMP& PRIVATE (Cabd, Fc, renormalizeHoppingFactorAAp, renormalizeHoppingFactorBBp, renormalizeHoppingFactorABp, renormalizeHoppingFactorBAp), &
         !$OMP& PRIVATE (uvec, vvec, vvecR, epsKax), &
         !$OMP& PRIVATE (lambda0, lambda3, lambda6, lambda6b, xi0, xi3, xi6, xi6b, x3, x6, x6b, kappa0, kappa6, kappa6b), &
         !$OMP& PRIVATE (signChange), &
         !$OMP& PRIVATE (a0AA, b0AA, c0AA, d0AA, h0AA, j0AA, k0AA), &
         !$OMP& PRIVATE (a3AA, b3AA, c3AA, d3AA), &
         !$OMP& PRIVATE (a6AA, b6AA, c6AA, d6AA), &
         !$OMP& PRIVATE (a0BA, b0BA, c0BA, d0BA, h0BA, j0BA, k0BA), &
         !$OMP& PRIVATE (a3BA, b3BA, c3BA, d3BA), &
         !$OMP& PRIVATE (a6BA, b6BA, c6BA, d6BA), &
         !$OMP& PRIVATE (a0AB, b0AB, c0AB, d0AB, h0AB, j0AB, k0AB), &
         !$OMP& PRIVATE (a3AB, b3AB, c3AB, d3AB), &
         !$OMP& PRIVATE (a6AB, b6AB, c6AB, d6AB), &
         !$OMP& PRIVATE (a0BB, b0BB, c0BB, d0BB, h0BB, j0BB, k0BB), &
         !$OMP& PRIVATE (a3BB, b3BB, c3BB, d3BB), &
         !$OMP& PRIVATE (a6BB, b6BB, c6BB, d6BB), &
         !$OMP& PRIVATE (V0AA, V3AA, V6AA), &
         !$OMP& PRIVATE (V0AB, V3AB, V6AB), &
         !$OMP& PRIVATE (V0AAp, V3AAp, V6AAp), &
         !$OMP& PRIVATE (V0ABp, V3ABp, V6ABp), &
         !$OMP& PRIVATE (V0BAp, V3BAp, V6BAp), &
         !$OMP& PRIVATE (V0BBp, V3BBp, V6BBp), &
         !$OMP& PRIVATE (checking, theta12, theta21, mmm,mm,kk,kkk),&
         !$OMP& PRIVATE (aval, bval),&
         !$OMP& PRIVATE (CabG),&
         !$OMP& PRIVATE (expFactor, dist, distXY, distXYZ), &
         !$OMP& REDUCTION (+:numberOfDel1,numberOfDel2,numberOfDel3,numberOfInterlayerHoppings), &
         !$OMP& REDUCTION (+:numberOfHAA1, numberOfHAA2, numberOfHAB1, numberOfHAB2, numberOfHBA1, numberOfHBA2), &
         !$OMP& REDUCTION (+:countG1, countG2, countG3, countG4, countG5, countG6, countG7, countG8), &
         !$OMP& REDUCTION (+:countBN1, countBN2, countBN3, countBN4, countBN5, countBN6, countBN7, countBN8)
         do i=1,nAt
            delta = 0.0
            nlay = (Species(i)-1)/2 + 1
            del(1) = 0.0_dp
            del(2) = 0.0_dp
            del(3) = 0.0_dp
! ---  Assign the hopping modifications for each lattice i ---
            if (l) then ! MoireOffDiag flag (used for effective moire systems), next one is tBGOffDiag and GBNOffDiag
                if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell') .or. (MIO_StringComp(str,'ReadXYZ') .and. .not.(singleLayerXYZ))) then
                   !if (Rat(3,i).lt.20.0_dp) then
                   if (layerIndex(i).eq.1) then
                   !if (i <= nAtC1) then
                      !dx = (Rat(1,i)+ tauX1)*eps
                      !dy = (Rat(2,i)+ tauY1)*eps
                      dx = ((1.0_dp+eps) * cos(MoireBilayerBottomAngleGrad) - 1) * (Rat(1,i)) &
                                - ((1.0_dp+eps) * sin(MoireBilayerBottomAngleGrad) * (Rat(2,i)))
                      dy = ((1.0_dp+eps) * cos(MoireBilayerBottomAngleGrad) - 1) * (Rat(2,i)) &
                                + ((1.0_dp+eps) * sin(MoireBilayerBottomAngleGrad) * (Rat(1,i)))
                      dx = dx + tauX1
                      dy = dy + tauY1
                   else
                      !dx = ((1.0_dp+eps) * cos(MoireBilayerTopAngleGrad) - 1) * (Rat(1,i)+tauX2) &
                      !          - ((1.0_dp+eps) * sin(MoireBilayerTopAngleGrad) * (Rat(2,i)+tauY2))
                      !dy = ((1.0_dp+eps) * cos(MoireBilayerTopAngleGrad) - 1) * (Rat(2,i)+tauY2) &
                      !          + ((1.0_dp+eps) * sin(MoireBilayerTopAngleGrad) * (Rat(1,i)+tauX2))
                      if (rotateFirst) then
                         dx = ((1.0_dp+eps) * cos(MoireBilayerTopAngleGrad) - 1) * (Rat(1,i)) &
                                   - ((1.0_dp+eps) * sin(MoireBilayerTopAngleGrad) * (Rat(2,i)))
                         dy = ((1.0_dp+eps) * cos(MoireBilayerTopAngleGrad) - 1) * (Rat(2,i)) &
                                   + ((1.0_dp+eps) * sin(MoireBilayerTopAngleGrad) * (Rat(1,i)))
                         dx = dx + tauX2
                         dy = dy + tauY2 !- aG/sqrt(3.0_dp)  ! removed this, as in reality it is already at the right position for two BN layer right on top of each other
                      else
                         dx = ((1.0_dp+eps) * cos(MoireBilayerTopAngleGrad) - 1) * (Rat(1,i)) &
                                   - ((1.0_dp+eps) * sin(MoireBilayerTopAngleGrad) * (Rat(2,i) - tauY2/eps))
                         dy = ((1.0_dp+eps) * cos(MoireBilayerTopAngleGrad) - 1) * (Rat(2,i) - tauY2/eps) &
                                   + ((1.0_dp+eps) * sin(MoireBilayerTopAngleGrad) * (Rat(1,i)))
                      end if
                   end if
                !else if (MIO_StringComp(str,'TwistedBilayer')) then
                !        ! add interlayer hoppings
                else ! Konda, modify hopping terms
                   if (twisted) then
                       if (MIO_StringComp(BilayerModel,'Jeil') ) then
                          if (layerIndex(i).eq.1) then
                             twistAngleGrad = twistedBilayerAngleGrad
                             dx = Rat(1,i) * ((1.0_dp+eps) * cos(twistAngleGrad) - 1.0_dp) &
                                  - Rat(2,i) * (1.0_dp+eps) * sin(twistAngleGrad)
                             dy = Rat(2,i) * ((1.0_dp+eps) * cos(twistAngleGrad) - 1.0_dp) &
                                  + Rat(1,i) * (1.0_dp+eps) * sin(twistAngleGrad)
                          else
                             twistAngleGrad = -twistedBilayerAngleGrad
                             dx = Rat(1,i) * ((1.0_dp+eps) * cos(twistAngleGrad) - 1.0_dp) &
                                  - Rat(2,i) * (1.0_dp+eps) * sin(twistAngleGrad)
                             dy = Rat(2,i) * ((1.0_dp+eps) * cos(twistAngleGrad) - 1.0_dp) &
                                  + Rat(1,i) * (1.0_dp+eps) * sin(twistAngleGrad)
                          end if
                       else
                          dx = Rat(1,i) * ((1.0_dp+eps) * cos(twistAngleGrad) - 1.0_dp) &
                                  - Rat(2,i) * (1.0_dp+eps) * sin(twistAngleGrad)
                          dy = Rat(2,i) * ((1.0_dp+eps) * cos(twistAngleGrad) - 1.0_dp) &
                                  + Rat(1,i) * (1.0_dp+eps) * sin(twistAngleGrad)
                       end if
                   else
                       if (readRigidXYZ) then
                          dx = RatInit(1,i)*eps
                          dy = RatInit(2,i)*eps
                       else
                          dx = Rat(1,i)*eps
                          dy = Rat(2,i)*eps
                       end if
                   end if
                end if
                dy = dy + shift1
                if (layerIndex(i).eq.1) then
                    if (sign1.lt.0.0) then ! Go from BN to NB orientation if necessary
                        dx = -dx
                        dy = -dy
                    end if
                else if (layerIndex(i).eq.2) then
                    if (sign2.lt.0.0) then ! Go from BN to NB orientation for second layer of bilayer graphene
                        dx = -dx
                        dy = -dy
                    end if
                end if
                if (addDisplacements) then
                    dx = dx + displacements(1,i)
                    dy = dy + displacements(2,i)
                end if
                if (distanceDependentEffectiveModel) then
                   call distanceDependentC(Cabd, Cab, BfactorCab, interlayerDistances(i), z0)
                else
                   Cabd = Cab
                end if
                if (layerIndex(i).eq.2 .and. deactivateUpperLayer) then
                   Habjj = 0.0_dp
                else if (layerIndex(i).ne.1 .and. deactivateUpperLayers) then
                    Habjj = 0.0_dp
                else if (FanZhang) then
                    if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                       Habjj = 0.0_dp
                    else
                       call offdiago(Habjj,dx,dy,Cabd,Phiab)
                    end if
                else
                    if (MIO_StringComp(BilayerModel,'Jeil')) then
                       call offdiago2(Habjj,dx,dy,Cabd,Phiab)
                    else
                       call offdiago(Habjj,dx,dy,Cabd,Phiab)
                    end if
                end if
                realH = real(Habjj)
                imagH = imag(Habjj)
                !if (sign1.lt.0.0) then
                !  imagH = -imagH
                !end if
                HABreal(i) = realH
                HABimag(i) = imagH

                if (Species(i).eq.1) then   ! ---  B sublattice

                  del(1) =  2.0_dp*realH/3.0_dp
                  del(2) = -(  realH + sqrt(3.0_dp)*imagH )/3.0_dp
                  del(3) =  ( -realH + sqrt(3.0_dp)*imagH )/3.0_dp

                end if

                if (Species(i).eq.2) then   ! ---  A sublattice

                   del(1) = 2.0_dp*realH/3.0_dp
                   del(2) = ( -realH + sqrt(3.0_dp)*imagH )/3.0_dp
                   del(3) = -(  realH + sqrt(3.0_dp)*imagH )/3.0_dp
                 ! del(2) and del(3) are inverted to take into account difference
                 ! between A and B
                end if
            end if
            if (tBGOffDiag) then
              if (tBGOffDiag) then
                  if (tBGOffDiagPRB) then
                     CabG = CabG_global
                  else
                     call fitFunSrivani(0.0134_dp, -0.0977_dp,0.1790_dp, abs(interlayerDistances(i)), CabG)
                  end if
              else
                  call MIO_InputParameter('MoirePotCabG',CabG,0.001987_dp)
              end if
              CabG = CabG/g0
              !if (i.eq. 1) print*, "adding the offdiagonal intralyer terms for tBG"
              if (tBGuseDisplacementFile) then
                 if (layerIndex(i).eq.1) then
                     dx = displacements(1,i)*cos(GlobalPhi2)+displacements(2,i)*sin(GlobalPhi2)
                     dy = displacements(2,i)*cos(GlobalPhi2)-displacements(1,i)*sin(GlobalPhi2)
                 else if (layerIndex(i).eq.2) then
                     dx = displacements(1,i)*cos(GlobalPhi)+displacements(2,i)*sin(GlobalPhi)
                     dy = displacements(2,i)*cos(GlobalPhi)-displacements(1,i)*sin(GlobalPhi)
                 end if
                 if (layerIndex(i).eq.1) then
                    dxTemp = -dx/aG
                    dyTemp = -dy/aG
                 else if (layerIndex(i).eq.2) then
                    dxTemp = dx/aG
                    dyTemp = dy/aG
                 end if
              else
                 dx = ((1.0_dp+eps) * cos(tBGAngle) - 1.0_dp) * (Rat(1,i)) &
                           - ((1.0_dp+eps) * sin(tBGAngle) * (Rat(2,i)))
                 dy = ((1.0_dp+eps) * cos(tBGAngle) - 1.0_dp) * (Rat(2,i)) &
                           + ((1.0_dp+eps) * sin(tBGAngle) * (Rat(1,i)))
                 if (layerIndex(i).eq.1) then
                    dxTemp = dx/aG
                    dyTemp = dy/aG
                 else if (layerIndex(i).eq.2) then
                    dxTemp = -dx/aG
                    dyTemp = -dy/aG
                 end if
              end if
              if (tBGSwitchDxDy) then
                 if (layerIndex(i).eq.1) then ! graxhene
                    call offdiagotBG(Habjj,dyTemp,dxTemp,CabG,PhiabG)
                 else if (layerIndex(i).eq.2) then
                    call offdiagotBG(Habjj,dyTemp,dxTemp,CabG,PhiabG)
                 end if
              else
                 if (layerIndex(i).eq.1) then ! graphene
                    call offdiagotBG(Habjj,dxTemp,dyTemp,CabG,PhiabG)
                 else if (layerIndex(i).eq.2) then
                    call offdiagotBG(Habjj,dxTemp,dyTemp,CabG,PhiabG)
                 end if
              end if
              realH = real(Habjj)
              imagH = imag(Habjj)
              !print*, "realH", realH ! temp

              if (Species(i).eq.2) then   ! ---  C2 or N
              !if (Species(i).eq.1 .or. Species(i).eq.3) then   ! ---  C2 or N

                 del(1) = del(1) + 2.0_dp*realH/3.0_dp
                 del(2) = del(2) + ( -realH + sqrt(3.0_dp)*imagH )/3.0_dp
                 del(3) = del(3)  -(  realH + sqrt(3.0_dp)*imagH )/3.0_dp
              end if

              if (Species(i).eq.1) then   ! ---  C1 or B (note the conventions from the book chapter putting A sublattice at (0,0) while we put the B-sublattice there when using in-house code as well as the xyz files shared by Jiaqi)
              !if (Species(i).eq.2 .or. Species(i).eq.4) then   ! ---  C1 or B

                del(1) = del(1) + 2.0_dp*realH/3.0_dp
                del(2) = del(2) -(  realH + sqrt(3.0_dp)*imagH )/3.0_dp
                del(3) = del(3) + ( -realH + sqrt(3.0_dp)*imagH )/3.0_dp

              end if
            end if
            if (GBNOffDiag) then
                !if (i.eq. 1) print*, "adding the offdiagonal intralyer terms for GBN"
                if (GBNuseDisplacementFile) then
                   if (encapsulatedThreeLayers .and. layerIndex(i).eq.2) then
                      dx_b = displacements_b(1,i)*cos(GlobalPhiL2a)+displacements_b(2,i)*sin(GlobalPhiL2a)
                      dy_b = displacements_b(2,i)*cos(GlobalPhiL2a)-displacements_b(1,i)*sin(GlobalPhiL2a)
                      dx_t = displacements_t(1,i)*cos(GlobalPhiL2b)+displacements_t(2,i)*sin(GlobalPhiL2b)
                      dy_t = displacements_t(2,i)*cos(GlobalPhiL2b)-displacements_t(1,i)*sin(GlobalPhiL2b)
                   else if (encapsulatedThreeLayers .and. layerIndex(i).eq.1) then
                      dx = displacements(1,i)*cos(GlobalPhiL1)+displacements(2,i)*sin(GlobalPhiL1)
                      dy = displacements(2,i)*cos(GlobalPhiL1)-displacements(1,i)*sin(GlobalPhiL1)
                   else if (encapsulatedThreeLayers .and. layerIndex(i).eq.3) then
                      dx = displacements(1,i)*cos(GlobalPhiL3)+displacements(2,i)*sin(GlobalPhiL3)
                      dy = displacements(2,i)*cos(GlobalPhiL3)-displacements(1,i)*sin(GlobalPhiL3)
                   else
                      dx = displacements(1,i)*cos(GlobalPhi)+displacements(2,i)*sin(GlobalPhi)
                      dy = displacements(2,i)*cos(GlobalPhi)-displacements(1,i)*sin(GlobalPhi)
                   end if
                   !dx = displacements(1,i)
                   !dy = displacements(2,i)
                   if (GBNtwoLayers) then
                      if (layerIndex(i).eq.1) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.2) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   !else if (BNBNtwoLayers) then
                   !   if (layerIndex(i).eq.1) then
                   !      dxTemp = dx/aBN
                   !      dyTemp = dy/aBN
                   !   else if (layerIndex(i).eq.2) then
                   !      dxTemp = -dx/aBN
                   !      dyTemp = -dy/aBN
                   !   end if
                   else if (encapsulatedThreeLayers) then
                      if (layerIndex(i).eq.2) then
                         dxTemp_b = dx_b/aG
                         dyTemp_b = dy_b/aG
                         dxTemp_t = dx_t/aG
                         dyTemp_t = dy_t/aG
                      else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   else if (encapsulatedFourLayers) then
                      if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   else if (encapsulatedFiveLayers) then
                      if (layerIndex(i).eq.2 .or. layerIndex(i).eq.4) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.5) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   else if (t3GwithBN) then
                      if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.4) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.1) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   else if (BNt2GBN) then
                      if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   else if (t2GBN) then
                      if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then ! G layers
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.1) then ! BN layers
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   end if
                else
                   dx = ((1.0_dp+eps) * cos(GBNAngle) - 1.0_dp) * (Rat(1,i)) &
                             - ((1.0_dp+eps) * sin(GBNAngle) * (Rat(2,i)))
                   dy = ((1.0_dp+eps) * cos(GBNAngle) - 1.0_dp) * (Rat(2,i)) &
                             + ((1.0_dp+eps) * sin(GBNAngle) * (Rat(1,i)))
                   if (GBNtwoLayers) then
                      if (layerIndex(i).eq.1) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.2) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   !else if (BNBNtwoLayers) then
                   !   if (layerIndex(i).eq.1) then
                   !      dxTemp = dx/aG
                   !      dyTemp = dy/aG
                   !   else if (layerIndex(i).eq.2) then
                   !      dxTemp = -dx/aBN
                   !      dyTemp = -dy/aBN
                   !   end if
                   else if (encapsulatedThreeLayers) then
                      if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   else if (encapsulatedFourLayers) then
                      if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   else if (encapsulatedFiveLayers) then
                      if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   else if (t3GwithBN) then
                      if (layerIndex(i).eq.4 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.2) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.1) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   else if (BNt2GBN) then
                      if (layerIndex(i).eq.3 .or. layerIndex(i).eq.2) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   else if (t2GBN) then
                      if (layerIndex(i).eq.3 .or. layerIndex(i).eq.2) then
                         dxTemp = dx/aG
                         dyTemp = dy/aG
                      else if (layerIndex(i).eq.1) then
                         dxTemp = -dx/aBN
                         dyTemp = -dy/aBN
                      end if
                   end if
                end if
                if (GBNtwoLayers) then
                   if (layerIndex(i).eq.1) then ! graphene
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabG,PhiabG)
                   else if (layerIndex(i).eq.2) then
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabBN,PhiabBN)
                   end if
                !else if (BNBNtwoLayers) then
                !   if (layerIndex(i).eq.1) then ! graphene
                !      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabBN,PhiabBN)
                !   else if (layerIndex(i).eq.2) then
                !      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabBN,PhiabBN)
                !   end if
                else if (encapsulatedThreeLayers) then
                   if (layerIndex(i).eq.2) then ! graphene
                      call offdiagoGBN(Habjj_b,dxTemp_b,dyTemp_b,CabG,PhiabG)
                      call offdiagoGBN(Habjj_t,dxTemp_t,dyTemp_t,CabG,PhiabG)
                   else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3) then
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabBN,PhiabBN)
                   end if
                else if (encapsulatedFourLayers) then
                   if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then ! graphene
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabG,PhiabG)
                   else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabBN,PhiabBN)
                   end if
                else if (encapsulatedFiveLayers) then
                   if (layerIndex(i).eq.2 .or. layerIndex(i).eq.4) then ! graphene
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabG,PhiabG)
                   else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.5) then
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabBN,PhiabBN)
                   end if
                else if (t3GwithBN) then
                   if (layerIndex(i).eq.4 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.2) then ! graphene
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabG,PhiabG)
                   else if (layerIndex(i).eq.1) then
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabBN,PhiabBN)
                   end if
                else if (BNt2GBN) then
                   if (layerIndex(i).eq.3 .or. layerIndex(i).eq.2) then ! graphene
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabG,PhiabG)
                   else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabBN,PhiabBN)
                   end if
                else if (t2GBN) then
                   if (layerIndex(i).eq.3 .or. layerIndex(i).eq.2) then ! graphene
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabG,PhiabG)
                   else if (layerIndex(i).eq.1) then
                      call offdiagoGBN(Habjj,dxTemp,dyTemp,CabBN,PhiabBN)
                   end if
                end if
                if (encapsulatedThreeLayers .and. layerIndex(i).eq.2) then
                   realH_b = real(Habjj_b)
                   imagH_b = imag(Habjj_b)
                   realH_t = real(Habjj_t)
                   imagH_t = imag(Habjj_t)
                else
                   realH = real(Habjj)
                   imagH = imag(Habjj)
                end if
                !print*, "realH", realH ! temp

                if (Species(i).eq.2 .or. Species(i).eq.4) then   ! ---  C2 or N
                !if (Species(i).eq.1 .or. Species(i).eq.3) then   ! ---  C2 or N
                   if (encapsulatedThreeLayers .and. layerIndex(i).eq.2) then
                      del(1) = del(1) + 2.0_dp*realH_b/3.0_dp
                      del(2) = del(2) + ( -realH_b + sqrt(3.0_dp)*imagH_b )/3.0_dp
                      del(3) = del(3)  -(  realH_b + sqrt(3.0_dp)*imagH_b )/3.0_dp
                      if (removeTopMoireInL2) then
                         del(1) = del(1)
                         del(2) = del(2)
                         del(3) = del(3)
                      else
                         del(1) = del(1) + 2.0_dp*realH_t/3.0_dp
                         del(2) = del(2) + ( -realH_t + sqrt(3.0_dp)*imagH_t )/3.0_dp
                         del(3) = del(3)  -(  realH_t + sqrt(3.0_dp)*imagH_t )/3.0_dp
                      end if
                   else if (GBNtwoLayers) then
                      del(1) = 2.0_dp*realH/3.0_dp
                      del(2) = ( -realH + sqrt(3.0_dp)*imagH )/3.0_dp
                      del(3) = -(  realH + sqrt(3.0_dp)*imagH )/3.0_dp
                   else
                      del(1) = del(1) + 2.0_dp*realH/3.0_dp
                      del(2) = del(2) + ( -realH + sqrt(3.0_dp)*imagH )/3.0_dp
                      del(3) = del(3)  -(  realH + sqrt(3.0_dp)*imagH )/3.0_dp
                   end if
                end if

                if (Species(i).eq.1 .or. Species(i).eq.3) then   ! ---  C1 or B (note the conventions from the book chapter putting A sublattice at (0,0) while we put the B-sublattice there when using in-house code as well as the xyz files shared by Jiaqi)
                !if (Species(i).eq.2 .or. Species(i).eq.4) then   ! ---  C1 or B
                  if (encapsulatedThreeLayers .and. layerIndex(i).eq.2) then
                     del(1) = del(1) + 2.0_dp*realH_b/3.0_dp
                     del(2) = del(2)  -(  realH_b + sqrt(3.0_dp)*imagH_b )/3.0_dp
                     del(3) = del(3) + ( -realH_b + sqrt(3.0_dp)*imagH_b )/3.0_dp
                     if (removeTopMoireInL2) then
                        del(1) = del(1)
                        del(2) = del(2)
                        del(3) = del(3)
                     else
                        del(1) = del(1) + 2.0_dp*realH_t/3.0_dp
                        del(2) = del(2)  -(  realH_t + sqrt(3.0_dp)*imagH_t )/3.0_dp
                        del(3) = del(3) + ( -realH_t + sqrt(3.0_dp)*imagH_t )/3.0_dp
                     end if
                  else if (GBNtwoLayers) then
                     del(1) =  2.0_dp*realH/3.0_dp
                     del(2) = -(  realH + sqrt(3.0_dp)*imagH )/3.0_dp
                     del(3) =  ( -realH + sqrt(3.0_dp)*imagH )/3.0_dp
                  else
                     del(1) = del(1) + 2.0_dp*realH/3.0_dp
                     del(2) = del(2)  -(  realH + sqrt(3.0_dp)*imagH )/3.0_dp
                     del(3) = del(3) + ( -realH + sqrt(3.0_dp)*imagH )/3.0_dp
                  end if
                end if
            !else if (tBGOffDiag) then
                  !del(1) = 0.0_dp
                  !del(2) = 0.0_dp
                  !del(3) = 0.0_dp
            end if
            if (MIO_StringComp(BilayerModel,'Jeil')) then
                jj1 = 0
                jj2 = 0
                maxDist = 1000.0_dp
                do j=1,Nneigh(i)
                   if(layerIndex(i).ne.layerIndex(NList(j,i))) then
                      d = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                      if ((d.lt.maxDist) .and. d.lt.(0.25_dp*aCC)) then
                         jj1 = j
                         maxDist = d
                      end if
                   end if
                end do
                if (jj1.eq.0) then
                  do j=1,Nneigh(i)
                     if(layerIndex(i).ne.layerIndex(NList(j,i))) then
                        d = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                        if (jj1.eq.0) then
                           if (d.lt.(0.75_dp*aCC)) then
                              jj1 = j
                           end if
                        else
                           if (d.lt.(0.75_dp*aCC)) then
                              jj2 = j
                           end if
                        end if
                     end if
                  end do
                  if (jj1.ne.0 .and. jj2.eq.0) then
                     jj1 = 0
                     jj2 = 0
                  end if
                end if
                !if (jj.eq.0) then
                !  do j=1,Nneigh(i)
                !    hopp(j,i) = 0.0_dp
                !  end do
                !  cycle
                !end if
                firstNN = 0
                k = 0
                if (jj1.ne.0 .and. jj2.eq.0) then
                   do j=1,Nneigh(jj1)
                      if(layerIndex(jj1).eq.layerIndex(NList(j,jj1))) then
                         k = k+1
                         firstNN(k) = NList(j,jj1)
                      end if
                   end do
                end if
                !if (k.ne.3) print*, "not enough first neighbors"
            end if
            do j=1,Nneigh(i)
               !if (abs(NeighD(3,j,i))<0.01_dp) then ! INPLANE/INTRALAYER
               jj = NList(j,i)
               if (layerIndex(i) .eq. layerIndex(jj)) then
                  d = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                  do ilvl=1,tbnn
                     if (d < Nradii(ilvl,nlay)) then ! Konda, assign del1, del2, del3 to delta
                     !if (d <  ((4.0_dp*aCC)**2) * 1.2_dp**2) then
                        !hopp(j,i) = cmplx(gn(Species(i),Species(NList(j,i)),ilvl))
                        if (d<1.7_dp) then
                           if (abs(NeighD(1,j,i)) < aG/4.0_dp) then
                              delta = del(1)
                              numberOfDel1 = numberOfDel1+1
                           else if (NeighD(1,j,i) < -aG/4.0_dp ) then
                              !if (sign1.lt.0.0) then
                              !    delta = del(3)
                              !else
                                 delta = del(3)
                              !end if
                              numberOfDel2 = numberOfDel2+1
                           else if (NeighD(1,j,i) > aG/4.0_dp) then
                              !if (sign1.lt.0.0) then
                              !    delta = del(2)
                              !else
                                 delta = del(2)
                              !end if
                              numberOfDel3 = numberOfDel3+1
                           endif
                        end if
                        !print*, delta
                        if (GBNtwoLayers) then
                           dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                           if (GBNtwoLayersF2G2s) then
                               if (layerIndex(i) .eq. 1) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLGfromGBN
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLGfromGBN
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLGfromGBN
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLGfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 2) then ! For hBN layer
                                   !if (i.eq.1) print*, "aBN=", aBN
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBNfromGBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBNfromGBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLBNfromGBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLBNfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                           else if (F2G2Model) then ! use the single layer F2G2
                               if (layerIndex(i) .eq. 1) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSL
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSL
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6K
                                      else
                                         hopp(j,i) = hopp(j,i) + t6K
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7K
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8K
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 2) then ! For hBN layer
                                   !if (i.eq.1) print*, "aBN=", aBN
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                              ! Add GBNtwoLayers F2G2 hopping terms
                           end if
                           if (d > 1.7_dp) then
                               delta = 0.0_dp
                           end if
                           !print*, "we are adding the gn values"
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp) ! Konda, use the chosen value of delta
                           !if (i.eq.1) print*, "GBNF2G2", hopp(j,i)
                        else if (BNBNtwoLayers) then
                           dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                           !if (GBNtwoLayersF2G2s) then
                           !    if (layerIndex(i) .eq. 1) then ! For graphene layer
                           !        if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                           !           countG1 = countG1+1
                           !           hopp(j,i) = hopp(j,i)
                           !        else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                           !           countG2 = countG2+1
                           !           if (Species(i).eq.1) then
                           !              hopp(j,i) = hopp(j,i) + t2KSLGfromGBNA
                           !           else
                           !              hopp(j,i) = hopp(j,i) + t2KSLGfromGBNB
                           !           end if
                           !        else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                           !           countG3 = countG3+1
                           !           hopp(j,i) = hopp(j,i) + t3KSLGfromGBN
                           !        else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                           !           countG4 = countG4+1
                           !           hopp(j,i) = hopp(j,i) + t4KSLGfromGBN
                           !        else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                           !           countG5 = countG5+1
                           !           if (Species(i).eq.1) then
                           !              hopp(j,i) = hopp(j,i) + t5KSLGfromGBNA
                           !           else
                           !              hopp(j,i) = hopp(j,i) + t5KSLGfromGBNB
                           !           end if
                           !        else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                           !           countG6 = countG6+1
                           !           if (Species(i).eq.1) then
                           !              hopp(j,i) = hopp(j,i) + t6KSLGfromGBNA
                           !           else
                           !              hopp(j,i) = hopp(j,i) + t6KSLGfromGBNB
                           !           end if
                           !        else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                           !           countG7 = countG7+1
                           !           hopp(j,i) = hopp(j,i) + t7KSLGfromGBN
                           !        else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                           !           countG8 = countG8+1
                           !           hopp(j,i) = hopp(j,i) + t8KSLGfromGBN
                           !        else
                           !             hopp(j,i) = hopp(j,i) + 0.0_dp
                           !        end if
                           !    else if (layerIndex(i) .eq. 2) then ! For hBN layer
                           !        !if (i.eq.1) print*, "aBN=", aBN
                           !        aBN1 = aBN/sqrt(3.0_dp)
                           !        if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                           !           countBN1 = countBN1+1
                           !           hopp(j,i) = hopp(j,i)
                           !        else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                           !           countBN2 = countBN2+1
                           !           if (Species(i).eq.3) then
                           !              hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNA
                           !           else
                           !              hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNB
                           !           end if
                           !        else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                           !           countBN3 = countBN3+1
                           !           hopp(j,i) = hopp(j,i) + t3KSLBNfromGBN
                           !        else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                           !           countBN4 = countBN4+1
                           !           hopp(j,i) = hopp(j,i) + t4KSLBNfromGBN
                           !        else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                           !           countBN5 = countBN5+1
                           !           if (Species(i).eq.3) then
                           !              hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNA
                           !           else
                           !              hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNB
                           !           end if
                           !        else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                           !           countBN6 = countBN6+1
                           !           if (Species(i).eq.3) then
                           !              hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNA
                           !           else
                           !              hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNB
                           !           end if
                           !        else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                           !           countBN7 = countBN7+1
                           !           hopp(j,i) = hopp(j,i) + t7KSLBNfromGBN
                           !        else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                           !           countBN8 = countBN8+1
                           !           hopp(j,i) = hopp(j,i) + t8KSLBNfromGBN
                           !        else
                           !             hopp(j,i) = hopp(j,i) + 0.0_dp
                           !        end if
                           !    end if
                           !else if (F2G2Model) then ! use the single layer F2G2
                               !if (layerIndex(i) .eq. 1) then ! For graphene layer
                               !Same model for both layers
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSL
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSL
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6K
                                      else
                                         hopp(j,i) = hopp(j,i) + t6K
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7K
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8K
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               !else if (layerIndex(i) .eq. 2) then ! For hBN layer
                               !    !if (i.eq.1) print*, "aBN=", aBN
                               !    aBN1 = aBN/sqrt(3.0_dp)
                               !    if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                               !       countBN1 = countBN1+1
                               !       hopp(j,i) = hopp(j,i)
                               !    else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                               !       countBN2 = countBN2+1
                               !       if (Species(i).eq.3) then
                               !          hopp(j,i) = hopp(j,i) + t2KSLBN_B
                               !       else
                               !          hopp(j,i) = hopp(j,i) + t2KSLBN_N
                               !       end if
                               !    else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                               !       countBN3 = countBN3+1
                               !       hopp(j,i) = hopp(j,i) + t3KSLBN
                               !    else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                               !       countBN4 = countBN4+1
                               !       hopp(j,i) = hopp(j,i) + t4KSLBN
                               !    else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                               !       countBN5 = countBN5+1
                               !       if (Species(i).eq.3) then
                               !          hopp(j,i) = hopp(j,i) + t5KSLBN_B
                               !       else
                               !          hopp(j,i) = hopp(j,i) + t5KSLBN_N
                               !       end if
                               !    else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                               !       countBN6 = countBN6+1
                               !       if (Species(i).eq.3) then
                               !          hopp(j,i) = hopp(j,i) + t6KBN
                               !       else
                               !          hopp(j,i) = hopp(j,i) + t6KBN
                               !       end if
                               !    else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                               !       countBN7 = countBN7+1
                               !       hopp(j,i) = hopp(j,i) + t7KBN
                               !    else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                               !       countBN8 = countBN8+1
                               !       hopp(j,i) = hopp(j,i) + t8KBN
                               !    else
                               !         hopp(j,i) = hopp(j,i) + 0.0_dp
                               !    end if
                               !end if
                              ! Add GBNtwoLayers F2G2 hopping terms
                           !end if
                           if (d > 1.7_dp) then
                               delta = 0.0_dp
                           end if
                           !print*, "we are adding the gn values"
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp) ! Konda, use the chosen value of delta
                           !if (i.eq.1) print*, "GBNF2G2", hopp(j,i)
                        else if (encapsulatedThreeLayers) then
                           dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                           ! Add GBNtwoLayers F2G2 hopping terms
                           if (GBNtwoLayersF2G2s) then
                               if (layerIndex(i) .eq. 2) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNB
                                      end if
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (layerIndex(i).eq.2) then
                                      !        if (Species(i).eq.1) then
                                      !             hopp(j,i) = hopp(j,i) + t2KA
                                      !        else if (Species(i).eq.2) then
                                      !             hopp(j,i) = hopp(j,i) + t2KB
                                      !        end if
                                      !    else if (layerIndex(i).eq.3) then
                                      !        if (Species(i).eq.2) then
                                      !             hopp(j,i) = hopp(j,i) + t2KA
                                      !        else if (Species(i).eq.1) then
                                      !             hopp(j,i) = hopp(j,i) + t2KB
                                      !        end if
                                      !    end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLGfromGBN
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (Species(i).eq.1) then
                                      !         hopp(j,i) = hopp(j,i) + t3K
                                      !    else if (Species(i).eq.2) then
                                      !         hopp(j,i) = hopp(j,i) + t3K
                                      !    end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLGfromGBN
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (Species(i).eq.1) then
                                      !         hopp(j,i) = hopp(j,i) + t4K
                                      !    else if (Species(i).eq.2) then
                                      !         hopp(j,i) = hopp(j,i) + t4K
                                      !    end if
                                      !end if
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNB
                                      end if
                                      !if (encapsulatedFourLayersF2G2) then
                                      !   if (layerIndex(i).eq.2) then
                                      !       if (Species(i).eq.1) then
                                      !            hopp(j,i) = hopp(j,i) + t5KA
                                      !       else if (Species(i).eq.2) then
                                      !            hopp(j,i) = hopp(j,i) + t5KB
                                      !       end if
                                      !   else if (layerIndex(i).eq.3) then
                                      !       if (Species(i).eq.2) then
                                      !            hopp(j,i) = hopp(j,i) + t5KA
                                      !       else if (Species(i).eq.1) then
                                      !            hopp(j,i) = hopp(j,i) + t5KB
                                      !       end if
                                      !   end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLGfromGBN
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLGfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 3) then ! For hBN layer
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBNfromGBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then  ! to account for the fact we are using aBN here
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBNfromGBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLBNfromGBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLBNfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                           else if (F2G2Model) then ! use the single layer F2G2
                               if (layerIndex(i) .eq. 2) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSL
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSL
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6K
                                      else
                                         hopp(j,i) = hopp(j,i) + t6K
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7K
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8K
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 1 .or. layerIndex(i).eq. 3) then ! For hBN layer
                                   !if (i.eq.1) print*, "aBN=", aBN
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                           end if
                           if (d > 1.7_dp) then
                               delta = 0.0_dp
                           end if
                           !print*, "we are adding the gn values"
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp) ! Konda, use the chosen value of delta
                           !if (i.eq.1) print*, "GBNF2G2", hopp(j,i)
                        else if (encapsulatedFourLayers) then
                           dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                           if (GBNtwoLayersF2G2s) then
                            ! Add GBNtwoLayers F2G2 hopping terms
                               if (layerIndex(i) .eq. 2 .or. layerIndex(i).eq.3) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNB
                                      end if
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (layerIndex(i).eq.2) then
                                      !        if (Species(i).eq.1) then
                                      !             hopp(j,i) = hopp(j,i) + t2KA
                                      !        else if (Species(i).eq.2) then
                                      !             hopp(j,i) = hopp(j,i) + t2KB
                                      !        end if
                                      !    else if (layerIndex(i).eq.3) then
                                      !        if (Species(i).eq.2) then
                                      !             hopp(j,i) = hopp(j,i) + t2KA
                                      !        else if (Species(i).eq.1) then
                                      !             hopp(j,i) = hopp(j,i) + t2KB
                                      !        end if
                                      !    end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLGfromGBN
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (Species(i).eq.1) then
                                      !         hopp(j,i) = hopp(j,i) + t3K
                                      !    else if (Species(i).eq.2) then
                                      !         hopp(j,i) = hopp(j,i) + t3K
                                      !    end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLGfromGBN
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (Species(i).eq.1) then
                                      !         hopp(j,i) = hopp(j,i) + t4K
                                      !    else if (Species(i).eq.2) then
                                      !         hopp(j,i) = hopp(j,i) + t4K
                                      !    end if
                                      !end if
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNB
                                      end if
                                      !if (encapsulatedFourLayersF2G2) then
                                      !   if (layerIndex(i).eq.2) then
                                      !       if (Species(i).eq.1) then
                                      !            hopp(j,i) = hopp(j,i) + t5KA
                                      !       else if (Species(i).eq.2) then
                                      !            hopp(j,i) = hopp(j,i) + t5KB
                                      !       end if
                                      !   else if (layerIndex(i).eq.3) then
                                      !       if (Species(i).eq.2) then
                                      !            hopp(j,i) = hopp(j,i) + t5KA
                                      !       else if (Species(i).eq.1) then
                                      !            hopp(j,i) = hopp(j,i) + t5KB
                                      !       end if
                                      !   end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLGfromGBN
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLGfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 4) then ! For hBN layer
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBNfromGBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then  ! to account for the fact we are using aBN here
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBNfromGBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLBNfromGBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLBNfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                           else if (F2G2Model) then ! use the single layer F2G2
                               if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSL
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSL
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6K
                                      else
                                         hopp(j,i) = hopp(j,i) + t6K
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7K
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8K
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 1 .or. layerIndex(i).eq. 4) then ! For hBN layer
                                   !if (i.eq.1) print*, "aBN=", aBN
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                           end if
                           if (d > 1.7_dp) then
                               delta = 0.0_dp
                           end if
                           !print*, "we are adding the gn values"
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp) ! Konda, use the chosen value of delta
                           !if (i.eq.1) print*, "GBNF2G2", hopp(j,i)
                        else if (encapsulatedFiveLayers) then
                           dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                           if (GBNtwoLayersF2G2s) then
                            ! Add GBNtwoLayers F2G2 hopping terms
                               if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3 .or. layerIndex(i).eq.4) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNB
                                      end if
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (layerIndex(i).eq.2) then
                                      !        if (Species(i).eq.1) then
                                      !             hopp(j,i) = hopp(j,i) + t2KA
                                      !        else if (Species(i).eq.2) then
                                      !             hopp(j,i) = hopp(j,i) + t2KB
                                      !        end if
                                      !    else if (layerIndex(i).eq.3) then
                                      !        if (Species(i).eq.2) then
                                      !             hopp(j,i) = hopp(j,i) + t2KA
                                      !        else if (Species(i).eq.1) then
                                      !             hopp(j,i) = hopp(j,i) + t2KB
                                      !        end if
                                      !    end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLGfromGBN
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (Species(i).eq.1) then
                                      !         hopp(j,i) = hopp(j,i) + t3K
                                      !    else if (Species(i).eq.2) then
                                      !         hopp(j,i) = hopp(j,i) + t3K
                                      !    end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLGfromGBN
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (Species(i).eq.1) then
                                      !         hopp(j,i) = hopp(j,i) + t4K
                                      !    else if (Species(i).eq.2) then
                                      !         hopp(j,i) = hopp(j,i) + t4K
                                      !    end if
                                      !end if
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNB
                                      end if
                                      !if (encapsulatedFourLayersF2G2) then
                                      !   if (layerIndex(i).eq.2) then
                                      !       if (Species(i).eq.1) then
                                      !            hopp(j,i) = hopp(j,i) + t5KA
                                      !       else if (Species(i).eq.2) then
                                      !            hopp(j,i) = hopp(j,i) + t5KB
                                      !       end if
                                      !   else if (layerIndex(i).eq.3) then
                                      !       if (Species(i).eq.2) then
                                      !            hopp(j,i) = hopp(j,i) + t5KA
                                      !       else if (Species(i).eq.1) then
                                      !            hopp(j,i) = hopp(j,i) + t5KB
                                      !       end if
                                      !   end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLGfromGBN
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLGfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 5) then ! For hBN layer
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBNfromGBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then  ! to account for the fact we are using aBN here
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBNfromGBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLBNfromGBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLBNfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                           else if (F2G2Model) then ! use the single layer F2G2
                               if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3 .or. layerIndex(i) .eq. 4) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSL
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSL
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6K
                                      else
                                         hopp(j,i) = hopp(j,i) + t6K
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7K
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8K
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 1 .or. layerIndex(i).eq. 5) then ! For hBN layer
                                   !if (i.eq.1) print*, "aBN=", aBN
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                           end if
                           if (d > 1.7_dp) then
                               delta = 0.0_dp
                           end if
                           !print*, "we are adding the gn values"
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp) ! Konda, use the chosen value of delta
                           !if (i.eq.1) print*, "GBNF2G2", hopp(j,i)
                        else if (encapsulatedSixLayers) then
                           dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                           if (GBNtwoLayersF2G2s) then
                            ! Add GBNtwoLayers F2G2 hopping terms
                               if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3 .or. layerIndex(i).eq.4 .or. layerIndex(i).eq. 5) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNB
                                      end if
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (layerIndex(i).eq.2) then
                                      !        if (Species(i).eq.1) then
                                      !             hopp(j,i) = hopp(j,i) + t2KA
                                      !        else if (Species(i).eq.2) then
                                      !             hopp(j,i) = hopp(j,i) + t2KB
                                      !        end if
                                      !    else if (layerIndex(i).eq.3) then
                                      !        if (Species(i).eq.2) then
                                      !             hopp(j,i) = hopp(j,i) + t2KA
                                      !        else if (Species(i).eq.1) then
                                      !             hopp(j,i) = hopp(j,i) + t2KB
                                      !        end if
                                      !    end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLGfromGBN
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (Species(i).eq.1) then
                                      !         hopp(j,i) = hopp(j,i) + t3K
                                      !    else if (Species(i).eq.2) then
                                      !         hopp(j,i) = hopp(j,i) + t3K
                                      !    end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLGfromGBN
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (Species(i).eq.1) then
                                      !         hopp(j,i) = hopp(j,i) + t4K
                                      !    else if (Species(i).eq.2) then
                                      !         hopp(j,i) = hopp(j,i) + t4K
                                      !    end if
                                      !end if
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNB
                                      end if
                                      !if (encapsulatedFourLayersF2G2) then
                                      !   if (layerIndex(i).eq.2) then
                                      !       if (Species(i).eq.1) then
                                      !            hopp(j,i) = hopp(j,i) + t5KA
                                      !       else if (Species(i).eq.2) then
                                      !            hopp(j,i) = hopp(j,i) + t5KB
                                      !       end if
                                      !   else if (layerIndex(i).eq.3) then
                                      !       if (Species(i).eq.2) then
                                      !            hopp(j,i) = hopp(j,i) + t5KA
                                      !       else if (Species(i).eq.1) then
                                      !            hopp(j,i) = hopp(j,i) + t5KB
                                      !       end if
                                      !   end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLGfromGBN
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLGfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 6) then ! For hBN layer
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBNfromGBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then  ! to account for the fact we are using aBN here
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBNfromGBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLBNfromGBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLBNfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                           else if (F2G2Model) then ! use the single layer F2G2
                               if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3 .or. layerIndex(i) .eq. 4 .or. layerIndex(i) .eq. 5) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSL
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSL
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6K
                                      else
                                         hopp(j,i) = hopp(j,i) + t6K
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7K
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8K
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 1 .or. layerIndex(i).eq. 6) then ! For hBN layer
                                   !if (i.eq.1) print*, "aBN=", aBN
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                           end if
                           if (d > 1.7_dp) then
                               delta = 0.0_dp
                           end if
                           !print*, "we are adding the gn values"
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp) ! Konda, use the chosen value of delta
                           !if (i.eq.1) print*, "GBNF2G2", hopp(j,i)
                        else if (encapsulatedSevenLayers) then
                           dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                           if (GBNtwoLayersF2G2s) then
                            ! Add GBNtwoLayers F2G2 hopping terms
                               if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3 .or. layerIndex(i).eq.4 .or. layerIndex(i).eq.5 .or. layerIndex(i).eq.6 ) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNB
                                      end if
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (layerIndex(i).eq.2) then
                                      !        if (Species(i).eq.1) then
                                      !             hopp(j,i) = hopp(j,i) + t2KA
                                      !        else if (Species(i).eq.2) then
                                      !             hopp(j,i) = hopp(j,i) + t2KB
                                      !        end if
                                      !    else if (layerIndex(i).eq.3) then
                                      !        if (Species(i).eq.2) then
                                      !             hopp(j,i) = hopp(j,i) + t2KA
                                      !        else if (Species(i).eq.1) then
                                      !             hopp(j,i) = hopp(j,i) + t2KB
                                      !        end if
                                      !    end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLGfromGBN
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (Species(i).eq.1) then
                                      !         hopp(j,i) = hopp(j,i) + t3K
                                      !    else if (Species(i).eq.2) then
                                      !         hopp(j,i) = hopp(j,i) + t3K
                                      !    end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLGfromGBN
                                      !if (encapsulatedFourLayersF2G2) then
                                      !    if (Species(i).eq.1) then
                                      !         hopp(j,i) = hopp(j,i) + t4K
                                      !    else if (Species(i).eq.2) then
                                      !         hopp(j,i) = hopp(j,i) + t4K
                                      !    end if
                                      !end if
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNB
                                      end if
                                      !if (encapsulatedFourLayersF2G2) then
                                      !   if (layerIndex(i).eq.2) then
                                      !       if (Species(i).eq.1) then
                                      !            hopp(j,i) = hopp(j,i) + t5KA
                                      !       else if (Species(i).eq.2) then
                                      !            hopp(j,i) = hopp(j,i) + t5KB
                                      !       end if
                                      !   else if (layerIndex(i).eq.3) then
                                      !       if (Species(i).eq.2) then
                                      !            hopp(j,i) = hopp(j,i) + t5KA
                                      !       else if (Species(i).eq.1) then
                                      !            hopp(j,i) = hopp(j,i) + t5KB
                                      !       end if
                                      !   end if
                                      !end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLGfromGBN
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLGfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 1 .or. layerIndex(i) .eq. 7) then ! For hBN layer
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBNfromGBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then  ! to account for the fact we are using aBN here
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBNfromGBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLBNfromGBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLBNfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                           else if (F2G2Model) then ! use the single layer F2G2
                               if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3 .or. layerIndex(i) .eq. 4 .or. layerIndex(i).eq.5 .or. layerIndex(i).eq.6) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSL
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSL
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6K
                                      else
                                         hopp(j,i) = hopp(j,i) + t6K
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7K
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8K
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 1 .or. layerIndex(i).eq. 7) then ! For hBN layer
                                   !if (i.eq.1) print*, "aBN=", aBN
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                           end if
                           if (d > 1.7_dp) then
                               delta = 0.0_dp
                           end if
                           !print*, "we are adding the gn values"
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp) ! Konda, use the chosen value of delta
                           !if (i.eq.1) print*, "GBNF2G2", hopp(j,i)
                        else if (t3GwithBN) then
                           dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                            ! Add GBNtwoLayers F2G2 hopping terms
                           if (layerIndex(i).eq.2) then ! For graphene layer (the one contacting with hBN)
                               if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                  countG1 = countG1+1
                                  hopp(j,i) = hopp(j,i)
                               else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                  countG2 = countG2+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t2KSL !GfromGBNA
                                  else
                                     hopp(j,i) = hopp(j,i) + t2KSL !GfromGBNB
                                  end if
                               else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                  countG3 = countG3+1
                                  hopp(j,i) = hopp(j,i) + t3KSL !GfromGBN
                               else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                  countG4 = countG4+1
                                  hopp(j,i) = hopp(j,i) + t4KSL !GfromGBN
                               else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                  countG5 = countG5+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t5KSL !GfromGBNA
                                  else
                                     hopp(j,i) = hopp(j,i) + t5KSL !GfromGBNB
                                  end if
                               else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                  countG6 = countG6+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t6KSL !GfromGBNA
                                  else
                                     hopp(j,i) = hopp(j,i) + t6KSL !GfromGBNB
                                  end if
                               else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                  countG7 = countG7+1
                                  hopp(j,i) = hopp(j,i) + t7KSL !GfromGBN
                               else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                  countG8 = countG8+1
                                  hopp(j,i) = hopp(j,i) + t8KSL !GfromGBN
                               else
                                    hopp(j,i) = hopp(j,i) + 0.0_dp
                               end if
                           else if (layerIndex(i).eq.3) then ! For graphene layer using single layer F2G2
                               if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                  countG1 = countG1+1
                                  hopp(j,i) = hopp(j,i)
                               else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                  countG2 = countG2+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t2KSL
                                  else
                                     hopp(j,i) = hopp(j,i) + t2KSL
                                  end if
                               else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                  countG3 = countG3+1
                                  hopp(j,i) = hopp(j,i) + t3KSL
                               else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                  countG4 = countG4+1
                                  hopp(j,i) = hopp(j,i) + t4KSL
                               else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                  countG5 = countG5+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t5KSL
                                  else
                                     hopp(j,i) = hopp(j,i) + t5KSL
                                  end if
                               else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                  countG6 = countG6+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t6KSL
                                  else
                                     hopp(j,i) = hopp(j,i) + t6KSL
                                  end if
                               else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                  countG7 = countG7+1
                                  hopp(j,i) = hopp(j,i) + t7KSL
                               else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                  countG8 = countG8+1
                                  hopp(j,i) = hopp(j,i) + t8KSL
                               else
                                    hopp(j,i) = hopp(j,i) + 0.0_dp
                               end if
                           else if (layerIndex(i).eq.4) then ! For graphene layer using single layer F2G2
                               if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                  countG1 = countG1+1
                                  hopp(j,i) = hopp(j,i)
                               else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                  countG2 = countG2+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t2KSL
                                  else
                                     hopp(j,i) = hopp(j,i) + t2KSL
                                  end if
                               else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                  countG3 = countG3+1
                                  hopp(j,i) = hopp(j,i) + t3KSL
                               else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                  countG4 = countG4+1
                                  hopp(j,i) = hopp(j,i) + t4KSL
                               else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                  countG5 = countG5+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t5KSL
                                  else
                                     hopp(j,i) = hopp(j,i) + t5KSL
                                  end if
                               else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                  countG6 = countG6+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t6KSL
                                  else
                                     hopp(j,i) = hopp(j,i) + t6KSL
                                  end if
                               else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                  countG7 = countG7+1
                                  hopp(j,i) = hopp(j,i) + t7KSL
                               else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                  countG8 = countG8+1
                                  hopp(j,i) = hopp(j,i) + t8KSL
                               else
                                    hopp(j,i) = hopp(j,i) + 0.0_dp
                               end if
                           else if (layerIndex(i) .eq. 1) then ! For hBN layer
                               aBN1 = aBN/sqrt(3.0_dp)
                               if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                  countBN1 = countBN1+1
                                  hopp(j,i) = hopp(j,i)
                               else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                  countBN2 = countBN2+1
                                  if (Species(i).eq.3) then
                                     hopp(j,i) = hopp(j,i) + t2KSLBN_B !fromGBNA
                                  else
                                     hopp(j,i) = hopp(j,i) + t2KSLBN_N !fromGBNB
                                  end if
                               else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                  countBN3 = countBN3+1
                                  hopp(j,i) = hopp(j,i) + t3KSLBN !fromGBN
                               else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then  ! to account for the fact we are using aBN here
                                  countBN4 = countBN4+1
                                  hopp(j,i) = hopp(j,i) + t4KSLBN !fromGBN
                               else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                  countBN5 = countBN5+1
                                  if (Species(i).eq.3) then
                                     hopp(j,i) = hopp(j,i) + t5KSLBN_B !fromGBNA
                                  else
                                     hopp(j,i) = hopp(j,i) + t5KSLBN_N !fromGBNB
                                  end if
                               else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                  countBN6 = countBN6+1
                                  if (Species(i).eq.3) then
                                     hopp(j,i) = hopp(j,i) + t6KBN !fromGBNA
                                  else
                                     hopp(j,i) = hopp(j,i) + t6KBN !fromGBNB
                                  end if
                               else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                  countBN7 = countBN7+1
                                  hopp(j,i) = hopp(j,i) + t7KBN !fromGBN
                               else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                  countBN8 = countBN8+1
                                  hopp(j,i) = hopp(j,i) + t8KBN !fromGBN
                               else
                                    hopp(j,i) = hopp(j,i) + 0.0_dp
                               end if
                           end if
                           if (d > 1.7_dp) then
                               delta = 0.0_dp
                           end if
                           !print*, "we are adding the gn values"
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp) ! Konda, use the chosen value of delta
                           !if (i.eq.1) print*, "GBNF2G2", hopp(j,i)
                           !print*, "F2G2", hopp(j,i)
                        !else if (BNBNtwoLayers) then
                        !   ! Add BNBNtwoLayers F2G2 hopping terms
                        !   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                        !      hopp(j,i) = hopp(j,i)
                        !   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                        !      hopp(j,i) = hopp(j,i) + t2KSLBN
                        !   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                        !      hopp(j,i) = hopp(j,i) + t3KSLBN
                        !   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                        !      hopp(j,i) = hopp(j,i) + t4KSLBN
                        !   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                        !      hopp(j,i) = hopp(j,i) + t5KSLBN
                        !   else if (dist .gt. 2.0_dp*aG*0.9_dp .and. dist .lt. 2.0_dp*aG*1.1_dp) then
                        !        hopp(j,i) = hopp(j,i) + t6KBN
                        !   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                        !        hopp(j,i) = hopp(j,i) + t7KBN
                        !   else if (dist .gt. 4.0_dp*aCC*0.9_dp .and. dist .lt. 4.0_dp*aCC*1.1_dp) then
                        !        hopp(j,i) = hopp(j,i) + t8KBN
                        !   else
                        !        hopp(j,i) = hopp(j,i) + 0.0_dp
                        !   end if
                        else if (BNt2GBN) then
                           dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                            ! Add GBNtwoLayers F2G2 hopping terms
                           if (layerIndex(i).eq.2) then ! For graphene layer (the one contacting with hBN)
                               if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                  countG1 = countG1+1
                                  hopp(j,i) = hopp(j,i)
                               else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                  countG2 = countG2+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t2KSL !GfromGBNA
                                  else
                                     hopp(j,i) = hopp(j,i) + t2KSL !GfromGBNB
                                  end if
                               else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                  countG3 = countG3+1
                                  hopp(j,i) = hopp(j,i) + t3KSL !GfromGBN
                               else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                  countG4 = countG4+1
                                  hopp(j,i) = hopp(j,i) + t4KSL !GfromGBN
                               else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                  countG5 = countG5+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t5KSL !GfromGBNA
                                  else
                                     hopp(j,i) = hopp(j,i) + t5KSL !GfromGBNB
                                  end if
                               else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                  countG6 = countG6+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t6KSL !GfromGBNA
                                  else
                                     hopp(j,i) = hopp(j,i) + t6KSL !GfromGBNB
                                  end if
                               else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                  countG7 = countG7+1
                                  hopp(j,i) = hopp(j,i) + t7KSL !GfromGBN
                               else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                  countG8 = countG8+1
                                  hopp(j,i) = hopp(j,i) + t8KSL !GfromGBN
                               else
                                    hopp(j,i) = hopp(j,i) + 0.0_dp
                               end if
                           else if (layerIndex(i).eq.3) then ! For graphene layer using single layer F2G2
                               if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                  countG1 = countG1+1
                                  hopp(j,i) = hopp(j,i)
                               else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                  countG2 = countG2+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t2KSL
                                  else
                                     hopp(j,i) = hopp(j,i) + t2KSL
                                  end if
                               else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                  countG3 = countG3+1
                                  hopp(j,i) = hopp(j,i) + t3KSL
                               else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                  countG4 = countG4+1
                                  hopp(j,i) = hopp(j,i) + t4KSL
                               else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                  countG5 = countG5+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t5KSL
                                  else
                                     hopp(j,i) = hopp(j,i) + t5KSL
                                  end if
                               else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                  countG6 = countG6+1
                                  if (Species(i).eq.1) then
                                     hopp(j,i) = hopp(j,i) + t6KSL
                                  else
                                     hopp(j,i) = hopp(j,i) + t6KSL
                                  end if
                               else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                  countG7 = countG7+1
                                  hopp(j,i) = hopp(j,i) + t7KSL
                               else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                  countG8 = countG8+1
                                  hopp(j,i) = hopp(j,i) + t8KSL
                               else
                                    hopp(j,i) = hopp(j,i) + 0.0_dp
                               end if
                           else if (layerIndex(i) .eq. 1 .or. layerIndex(i).eq. 4) then ! For hBN layer
                               aBN1 = aBN/sqrt(3.0_dp)
                               if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                  countBN1 = countBN1+1
                                  hopp(j,i) = hopp(j,i)
                               else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                  countBN2 = countBN2+1
                                  if (Species(i).eq.3) then
                                     hopp(j,i) = hopp(j,i) + t2KSLBN_B !fromGBNA
                                  else
                                     hopp(j,i) = hopp(j,i) + t2KSLBN_N !fromGBNB
                                  end if
                               else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                  countBN3 = countBN3+1
                                  hopp(j,i) = hopp(j,i) + t3KSLBN !fromGBN
                               else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then  ! to account for the fact we are using aBN here
                                  countBN4 = countBN4+1
                                  hopp(j,i) = hopp(j,i) + t4KSLBN !fromGBN
                               else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                  countBN5 = countBN5+1
                                  if (Species(i).eq.3) then
                                     hopp(j,i) = hopp(j,i) + t5KSLBN_B !fromGBNA
                                  else
                                     hopp(j,i) = hopp(j,i) + t5KSLBN_N !fromGBNB
                                  end if
                               else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                  countBN6 = countBN6+1
                                  if (Species(i).eq.3) then
                                     hopp(j,i) = hopp(j,i) + t6KBN !fromGBNA
                                  else
                                     hopp(j,i) = hopp(j,i) + t6KBN !fromGBNB
                                  end if
                               else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                  countBN7 = countBN7+1
                                  hopp(j,i) = hopp(j,i) + t7KBN !fromGBN
                               else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                  countBN8 = countBN8+1
                                  hopp(j,i) = hopp(j,i) + t8KBN !fromGBN
                               else
                                    hopp(j,i) = hopp(j,i) + 0.0_dp
                               end if
                           end if
                           if (d > 1.7_dp) then
                               delta = 0.0_dp
                           end if
                           !print*, "we are adding the gn values"
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp) ! Konda, use the chosen value of delta
                           !if (i.eq.1) print*, "GBNF2G2", hopp(j,i)
                           !print*, "F2G2", hopp(j,i)
                        !else if (BNBNtwoLayers) then
                        !   ! Add BNBNtwoLayers F2G2 hopping terms
                        !   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                        !      hopp(j,i) = hopp(j,i)
                        !   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                        !      hopp(j,i) = hopp(j,i) + t2KSLBN
                        !   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                        !      hopp(j,i) = hopp(j,i) + t3KSLBN
                        !   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                        !      hopp(j,i) = hopp(j,i) + t4KSLBN
                        !   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                        !      hopp(j,i) = hopp(j,i) + t5KSLBN
                        !   else if (dist .gt. 2.0_dp*aG*0.9_dp .and. dist .lt. 2.0_dp*aG*1.1_dp) then
                        !        hopp(j,i) = hopp(j,i) + t6KBN
                        !   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                        !        hopp(j,i) = hopp(j,i) + t7KBN
                        !   else if (dist .gt. 4.0_dp*aCC*0.9_dp .and. dist .lt. 4.0_dp*aCC*1.1_dp) then
                        !        hopp(j,i) = hopp(j,i) + t8KBN
                        !   else
                        !        hopp(j,i) = hopp(j,i) + 0.0_dp
                        !   end if
                        else if (t2GBN) then
                           dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                           if (F2G2Model) then ! use the single layer F2G2
                               if (layerIndex(i) .eq. 2 .or. layerIndex(i) .eq. 3) then ! For graphene layer
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSL
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSL
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6K
                                      else
                                         hopp(j,i) = hopp(j,i) + t6K
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7K
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8K
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 1) then ! For hBN layer
                                   !if (i.eq.1) print*, "aBN=", aBN
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_B
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBN_N
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KBN
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                              ! Add GBNtwoLayers F2G2 hopping terms
                           else if (GBNtwoLayersF2G2s) then
                            ! Add GBNtwoLayers F2G2 hopping terms
                               if (layerIndex(i).eq.2) then ! For graphene layer (the one contacting with hBN)
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLGfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLGfromGBN
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLGfromGBN
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLGfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLGfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLGfromGBN
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLGfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i).eq.3) then ! For graphene layer using single layer F2G2
                                   if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      countG1 = countG1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      countG2 = countG2+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      countG3 = countG3+1
                                      hopp(j,i) = hopp(j,i) + t3KSL
                                   else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      countG4 = countG4+1
                                      hopp(j,i) = hopp(j,i) + t4KSL
                                   else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      countG5 = countG5+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSL
                                      end if
                                   else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then
                                      countG6 = countG6+1
                                      if (Species(i).eq.1) then
                                         hopp(j,i) = hopp(j,i) + t6KSL
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSL
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      countG7 = countG7+1
                                      hopp(j,i) = hopp(j,i) + t7KSL
                                   else if (dist .gt. 4.0_dp*aCC*0.97_dp .and. dist .lt. 4.0_dp*aCC*1.03_dp) then
                                      countG8 = countG8+1
                                      hopp(j,i) = hopp(j,i) + t8KSL
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               else if (layerIndex(i) .eq. 1) then ! For hBN layer
                                   aBN1 = aBN/sqrt(3.0_dp)
                                   if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                                      countBN1 = countBN1+1
                                      hopp(j,i) = hopp(j,i)
                                   else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                                      countBN2 = countBN2+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t2KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                                      countBN3 = countBN3+1
                                      hopp(j,i) = hopp(j,i) + t3KSLBNfromGBN
                                   else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then  ! to account for the fact we are using aBN here
                                      countBN4 = countBN4+1
                                      hopp(j,i) = hopp(j,i) + t4KSLBNfromGBN
                                   else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                                      countBN5 = countBN5+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t5KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. 2.0_dp*aBN*0.97_dp .and. dist .lt. 2.0_dp*aBN*1.03_dp) then
                                      countBN6 = countBN6+1
                                      if (Species(i).eq.3) then
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNA
                                      else
                                         hopp(j,i) = hopp(j,i) + t6KSLBNfromGBNB
                                      end if
                                   else if (dist .gt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aBN)**2.0_dp + (aBN1**2.0_dp)) *1.1_dp) then
                                      countBN7 = countBN7+1
                                      hopp(j,i) = hopp(j,i) + t7KSLBNfromGBN
                                   else if (dist .gt. 4.0_dp*aBN1*0.97_dp .and. dist .lt. 4.0_dp*aBN1*1.03_dp) then
                                      countBN8 = countBN8+1
                                      hopp(j,i) = hopp(j,i) + t8KSLBNfromGBN
                                   else
                                        hopp(j,i) = hopp(j,i) + 0.0_dp
                                   end if
                               end if
                           end if
                           if (d > 1.7_dp) then
                               delta = 0.0_dp
                           end if
                           !print*, "we are adding the gn values"
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp) ! Konda, use the chosen value of delta
                           !if (i.eq.1) print*, "GBNF2G2", hopp(j,i)
                        else if (MIO_StringComp(BilayerModel,'Koshino') .or. MIO_StringComp(BilayerModel,'Mayou') .or. MIO_StringComp(BilayerModel,'HTC') .or. MIO_StringComp(SinglelayerModel,'HTC') ) then
                           if (frac) call AtomsSetCart()
                           if (KoshinoIntralayer .or. MayouIntralayer) then
                              ! The interlayer renormalisation factor (KoshinoSR) is not
                              ! meant for intralayer pairs and may not be set yet here.
                              renormalizeHoppingFactorAAp = 1.0_dp
                              if (MIO_StringComp(BilayerModel,'Koshino')) then
                                 !print*, "adding the inplane Koshino terms"
                                 if (frac) call AtomsSetCart()
                                 if (strainedMoire) then
                                    dist = sqrt((NeighD(1,j,i)+ umax* sin(2.0_dp*pi*nPeriod*Rat(1,i)/LMoire))**2.0_dp  &
                                                +NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                                 else
                                    dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                                 end if
                                 !if (KoshinoSR) then
                                 !    renormalizeHoppingFactorAAp  = exp((abs(NeighD(3,j,i))-bval)/aval)
                                 !end if
                                 vppsigma = vppsigma0 * exp(-(dist-interlayerdistance)/(BLdelta)) * renormalizeHoppingFactorAAp
                                 vpppi = vpppi0 * exp(-(dist-aCC)/(BLdelta)) * renormalizeHoppingFactorAAp
                                 hopp(j,i) = hopp(j,i) + (vpppi * (1.0_dp - ((NeighD(3,j,i))/(dist))**2.0_dp) &
                                                   + vppsigma * ((NeighD(3,j,i))/(dist))**2.0_dp)
                              else if (MIO_StringComp(BilayerModel,'Mayou')) then
                                 !print*, "adding the out of plane Mayou terms"
                                 if (frac) call AtomsSetCart()
                                 if (strainedMoire) then
                                    dist = sqrt((NeighD(1,j,i)+ umax* sin(2.0_dp*pi*nPeriod*Rat(1,i)/LMoire))**2.0_dp  &
                                                +NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                                 else
                                    dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                                 end if
                                 Fc = 1.0_dp/(1.0_dp+exp((dist-rc)/lc))
                                 vppsigma = vppsigma0 * exp(qsigma*(1.0_dp-dist/interlayerdistance)) * Fc  * renormalizeHoppingFactorAAp
                                 vpppi = vpppi0 * exp(qpi*(1.0_dp - dist/aCC)) * Fc  * renormalizeHoppingFactorAAp
                                 hopp(j,i) = hopp(j,i) + (vpppi * (1.0_dp - ((NeighD(3,j,i))/(dist))**2.0_dp) &
                                                   + vppsigma * ((NeighD(3,j,i))/(dist))**2.0_dp)
                                 !hopp(j,i) = 0.0_dp
                                 !hopp(j,i) = - vppsigma * ((NeighD(3,j,i))/(dist))**2.0_dp
                              end if
                              !dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                              !vpppi = vpppi0 * exp(-(dist-aCC)/(BLdelta)) ! inplane
                              !vppsigma = vppsigma0 * exp(-(dist-interlayerdistance)/(BLdelta))
                              !hopp(j,i) =  (vpppi * (1.0_dp - ((NeighD(3,j,i))/(dist))**2.0_dp) &
                              !                + vppsigma * ((NeighD(3,j,i))/(dist))**2.0_dp)
                              !hopp(j,i) = -vpppi * (1.0_dp - ((NeighD(3,j,i))/(dist))**2.0_dp)
                           else
                              dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                              aCC = aG/sqrt(3.0_dp)
                              if (useBNGKaxiras) then
                                 if (dist .gt. 0.1_dp .and. dist .lt. aCC*1.1_dp) then
                                      hopp(j,i) = hopp(j,i)
                                 else if (dist .gt. aCC*1.1_dp .and. dist .lt. aG*1.1_dp) then
                                      if (Species(i).eq.3) then
                                          hopp(j,i) = hopp(j,i) + t2KB
                                      else if (Species(i).eq.4) then
                                          hopp(j,i) = hopp(j,i) + t2KN
                                      end if
                                 else if (dist .gt. aG*1.1_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      hopp(j,i) = hopp(j,i) + t3K
                                 else
                                      hopp(j,i) = hopp(j,i) + 0.0_dp
                                 end if
                              else
                                 if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                      !if (i.eq.1) print*, "assigning intralayer F2G2 terms 1"
                                      hopp(j,i) = hopp(j,i)
                                 else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                      !if (i.eq.1) print*, "assigning intralayer F2G2 terms 2"
                                      if (threeLayers .or. fourLayersSandwiched .or. helicalTwistedMBM .or. fiveLayersSandwiched .or. sixLayersSandwiched .or. sevenLayersSandwiched .or. eightLayersSandwiched .or. tenLayersSandwiched .or. twentyLayersSandwiched) then
                                           if ((layerIndex(i).eq.3) .or. (middleTwist)) then ! middleTwist MUST be included for all the above systems, otherwise, we are incorrectly doing bilayer F2G2
                                              hopp(j,i) = hopp(j,i) + t2KSL
                                           else
                                              if (layerIndex(i).eq.1) then
                                                  if (Species(i).eq.1) then
                                                       hopp(j,i) = hopp(j,i) + t2KA
                                                  else if (Species(i).eq.2) then
                                                       hopp(j,i) = hopp(j,i) + t2KB
                                                  end if
                                              else if (layerIndex(i).eq.2) then
                                                  if (Species(i).eq.2) then
                                                       hopp(j,i) = hopp(j,i) + t2KA
                                                  else if (Species(i).eq.1) then
                                                       hopp(j,i) = hopp(j,i) + t2KB
                                                  end if
                                              end if
                                           end if
                                      else if (t3BG) then
                                           if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.5) then
                                               if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t2KA
                                               else if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t2KB
                                               end if
                                           else if (layerIndex(i).eq.2 .or. layerIndex(i).eq.4 .or. layerIndex(i).eq.6) then
                                               if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t2KA
                                               else if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t2KB
                                               end if
                                           end if
                                      else if (twoLayers .or. oneLayer) then
                                           if (forceBilayerF2G2Intralayer) then
                                              if (layerIndex(i).eq.1) then
                                                  if (Species(i).eq.1) then
                                                       hopp(j,i) = hopp(j,i) + t2KA
                                                  else if (Species(i).eq.2) then
                                                       hopp(j,i) = hopp(j,i) + t2KB
                                                  end if
                                              else if (layerIndex(i).eq.2) then
                                                  if (Species(i).eq.2) then
                                                       hopp(j,i) = hopp(j,i) + t2KA
                                                  else if (Species(i).eq.1) then
                                                       hopp(j,i) = hopp(j,i) + t2KB
                                                  end if
                                              end if
                                           else
                                              hopp(j,i) = hopp(j,i) + t2KSL
                                           end if
                                      else ! fourLayers
                                           if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3) then
                                               if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t2KA
                                               else if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t2KB
                                               end if
                                           else if (layerIndex(i).eq.2 .or. layerIndex(i).eq.4) then
                                               if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t2KA
                                               else if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t2KB
                                               end if
                                           end if
                                      end if
                                 else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                      !if (i.eq.1) print*, "assigning intralayer F2G2 terms 3"
                                      if (threeLayers .or. fourLayersSandwiched .or. helicalTwistedMBM .or. fiveLayersSandwiched .or. sixLayersSandwiched .or. sevenLayersSandwiched .or. eightLayersSandwiched .or. tenLayersSandwiched .or. twentyLayersSandwiched) then
                                           if ((layerIndex(i).eq.3) .or. (middleTwist)) then
                                              hopp(j,i) = hopp(j,i) + t3KSL
                                           else
                                              if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                                   hopp(j,i) = hopp(j,i) + t3K
                                              else if (Species(i).eq.2) then
                                                   hopp(j,i) = hopp(j,i) + t3K
                                              end if
                                           end if
                                      else if (t3BG) then
                                           if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                                hopp(j,i) = hopp(j,i) + t3K
                                           else if (Species(i).eq.2) then
                                                hopp(j,i) = hopp(j,i) + t3K
                                           end if
                                      else if (twoLayers .or. oneLayer) then
                                           if (forceBilayerF2G2Intralayer) then
                                              if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                                   hopp(j,i) = hopp(j,i) + t3K
                                              else if (Species(i).eq.2) then
                                                   hopp(j,i) = hopp(j,i) + t3K
                                              end if
                                           else
                                              hopp(j,i) = hopp(j,i) + t3KSL
                                           end if
                                      else ! fourLayers
                                           if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                                hopp(j,i) = hopp(j,i) + t3K
                                           else if (Species(i).eq.2) then
                                                hopp(j,i) = hopp(j,i) + t3K
                                           end if
                                      end if
                                 else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                      !if (i.eq.1) print*, "assigning intralayer F2G2 terms 4"
                                      if (threeLayers .or. fourLayersSandwiched .or. helicalTwistedMBM .or. fiveLayersSandwiched .or. sixLayersSandwiched .or. sevenLayersSandwiched .or. eightLayersSandwiched .or. tenLayersSandwiched .or. twentyLayersSandwiched) then
                                           if ((layerIndex(i).eq.3) .or. (middleTwist)) then
                                              hopp(j,i) = hopp(j,i) + t4KSL
                                           else
                                              if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                                   hopp(j,i) = hopp(j,i) + t4K
                                              else if (Species(i).eq.2) then
                                                   hopp(j,i) = hopp(j,i) + t4K
                                              end if
                                           end if
                                      else if (t3BG) then
                                           if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                                hopp(j,i) = hopp(j,i) + t4K
                                           else if (Species(i).eq.2) then
                                                hopp(j,i) = hopp(j,i) + t4K
                                           end if
                                      else if (twoLayers .or. oneLayer) then
                                           if (forceBilayerF2G2Intralayer) then
                                              if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                                   hopp(j,i) = hopp(j,i) + t4K
                                              else if (Species(i).eq.2) then
                                                   hopp(j,i) = hopp(j,i) + t4K
                                              end if
                                           else
                                              hopp(j,i) = hopp(j,i) + t4KSL
                                           end if
                                      else ! fourLayers
                                           if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                                hopp(j,i) = hopp(j,i) + t4K
                                           else if (Species(i).eq.2) then
                                                hopp(j,i) = hopp(j,i) + t4K
                                           end if
                                      end if
                                 else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                      !if (i.eq.1) print*, "assigning intralayer F2G2 terms 5", delta
                                      if (threeLayers .or. fourLayersSandwiched .or. helicalTwistedMBM .or. fiveLayersSandwiched .or. sixLayersSandwiched .or. sevenLayersSandwiched .or. eightLayersSandwiched .or. tenLayersSandwiched .or. twentyLayersSandwiched) then
                                           if ((layerIndex(i).eq.3) .or. (middleTwist)) then
                                              hopp(j,i) = hopp(j,i) + t5KSL
                                           else
                                              if (layerIndex(i).eq.1) then
                                                  if (Species(i).eq.1) then
                                                       hopp(j,i) = hopp(j,i) + t5KA
                                                  else if (Species(i).eq.2) then
                                                       hopp(j,i) = hopp(j,i) + t5KB
                                                  end if
                                              else if (layerIndex(i).eq.2) then
                                                  if (Species(i).eq.2) then
                                                       hopp(j,i) = hopp(j,i) + t5KA
                                                  else if (Species(i).eq.1) then
                                                       hopp(j,i) = hopp(j,i) + t5KB
                                                  end if
                                              end if
                                           end if
                                      else if (t3BG) then
                                           if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.5) then
                                               if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t5KA
                                               else if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t5KB
                                               end if
                                           else if (layerIndex(i).eq.2 .or. layerIndex(i).eq.4 .or. layerIndex(i).eq.6) then
                                               if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t5KA
                                               else if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t5KB
                                               end if
                                           end if
                                      else if (twoLayers .or. oneLayer) then
                                           if (forceBilayerF2G2Intralayer) then
                                              if (layerIndex(i).eq.1) then
                                                  if (Species(i).eq.1) then
                                                       hopp(j,i) = hopp(j,i) + t5KA
                                                  else if (Species(i).eq.2) then
                                                       hopp(j,i) = hopp(j,i) + t5KB
                                                  end if
                                              else if (layerIndex(i).eq.2) then
                                                  if (Species(i).eq.2) then
                                                       hopp(j,i) = hopp(j,i) + t5KA
                                                  else if (Species(i).eq.1) then
                                                       hopp(j,i) = hopp(j,i) + t5KB
                                                  end if
                                              end if
                                           else
                                              hopp(j,i) = hopp(j,i) + t5KSL
                                           end if
                                      else ! four layers
                                           if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3) then
                                               if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t5KA
                                               else if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t5KB
                                               end if
                                           else if (layerIndex(i).eq.2 .or. layerIndex(i).eq.4) then
                                               if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t5KA
                                               else if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t5KB
                                               end if
                                           end if
                                           if (Species(i).eq.1) then
                                                hopp(j,i) = hopp(j,i) + t5KA
                                           else if (Species(i).eq.2) then
                                                hopp(j,i) = hopp(j,i) + t5KB
                                           end if
                                      end if
                                 else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then ! narrow: the seventh shell is only 4 % further out
                                      !if (i.eq.1) print*, "assigning intralayer F2G2 terms 6", delta
                                      hopp(j,i) = hopp(j,i) + t6K
                                 else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                      !if (i.eq.1) print*, "assigning intralayer F2G2 terms 7", delta
                                      hopp(j,i) = hopp(j,i) + t7K
                                 else if (dist .gt. 4.0_dp*aCC*0.9_dp .and. dist .lt. 4.0_dp*aCC*1.1_dp) then
                                      !if (i.eq.1) print*, "assigning intralayer F2G2 terms 8", delta
                                      hopp(j,i) = hopp(j,i) + t8K
                                 else
                                      !if (i.eq.1) print*, "assigning intralayer F2G2 terms beyond", dist
                                      hopp(j,i) = hopp(j,i) + 0.0_dp
                                 end if
                                 !if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                 !     hopp(j,i) = hopp(j,i)
                                 !else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                 !     if (Species(i).eq.1) then
                                 !          hopp(j,i) = hopp(j,i) + t2KA
                                 !     else if (Species(i).eq.2) then
                                 !          hopp(j,i) = hopp(j,i) + t2KB
                                 !     end if
                                 !else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                 !     hopp(j,i) = hopp(j,i) + t3K
                                 !else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                 !     hopp(j,i) = hopp(j,i) + t4K
                                 !else if (dist .gt. 4.0_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                 !     if (Species(i).eq.1) then
                                 !          hopp(j,i) = hopp(j,i) + t5KA
                                 !     else if (Species(i).eq.2) then
                                 !          hopp(j,i) = hopp(j,i) + t5KB
                                 !     end if
                                 !else if (dist .gt. 2.0_dp*aG*0.9_dp .and. dist .lt. 2.0_dp*aG*1.1_dp) then
                                 !     hopp(j,i) = hopp(j,i) + t6K
                                 !else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                 !     hopp(j,i) = hopp(j,i) + t7K
                                 !else if (dist .gt. 4.0_dp*aCC*0.9_dp .and. dist .lt. 4.0_dp*aCC*1.1_dp) then
                                 !     hopp(j,i) = hopp(j,i) + t8K
                                 !end if
                                 !     hopp(j,i) = hopp(j,i) + 0.0_dp
                              end if
                              if (d > 1.5_dp) then
                                  delta = 0.0_dp
                              end if
                              !print*, "we are adding the gn values"
                              hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp) ! Konda, use the chosen value of delta
                           end if
                        else if (MIO_StringComp(BilayerModel,'BLKaxiras') .or. MIO_StringComp(BilayerModel,'BLSrivani')) then ! this is where we assign the intralayer F2G2
                           dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                           aCC = aG/sqrt(3.0_dp)
                           if (useBNGKaxiras) then
                              if (dist .gt. 0.1_dp .and. dist .lt. aCC*1.1_dp) then
                                   hopp(j,i) = hopp(j,i)
                              else if (dist .gt. aCC*1.1_dp .and. dist .lt. aG*1.1_dp) then
                                   if (Species(i).eq.3) then
                                       hopp(j,i) = hopp(j,i) + t2KB
                                   else if (Species(i).eq.4) then
                                       hopp(j,i) = hopp(j,i) + t2KN
                                   end if
                              else if (dist .gt. aG*1.1_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                   hopp(j,i) = hopp(j,i) + t3K
                              else
                                   hopp(j,i) = hopp(j,i) + 0.0_dp
                              end if
                           else
                              !if (i.eq.1) counter1 = counter1 + 1
                              !if (i.eq.1) print*, counter1
                              if (dist .gt. aCC*0.9_dp .and. dist .lt. aCC*1.1_dp) then
                                   !if (i.eq.1) print*, "assigning intralayer F2G2 terms 1"
                                   hopp(j,i) = hopp(j,i)
                              else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                                   !if (i.eq.1) print*, "assigning intralayer F2G2 terms 2"
                                   if (threeLayers .or. fourLayersSandwiched .or. helicalTwistedMBM .or. fiveLayersSandwiched .or. sixLayersSandwiched .or. sevenLayersSandwiched .or. eightLayersSandwiched .or. tenLayersSandwiched .or. twentyLayersSandwiched) then
                                        if ((layerIndex(i).eq.3) .or. (middleTwist)) then
                                           hopp(j,i) = hopp(j,i) + t2KSL
                                        else
                                           if (layerIndex(i).eq.1) then
                                               if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t2KA
                                               else if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t2KB
                                               end if
                                           else if (layerIndex(i).eq.2) then
                                               if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t2KA
                                               else if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t2KB
                                               end if
                                           end if
                                        end if
                                   else if (twoLayers .or. oneLayer) then
                                        if (forceBilayerF2G2Intralayer) then
                                           if (layerIndex(i).eq.1) then
                                               if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t2KA
                                               else if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t2KB
                                               end if
                                           else if (layerIndex(i).eq.2) then
                                               if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t2KA
                                               else if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t2KB
                                               end if
                                           end if
                                        else
                                           hopp(j,i) = hopp(j,i) + t2KSL
                                        end if
                                   else ! fourLayers
                                        if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3) then
                                            if (Species(i).eq.1) then
                                                 hopp(j,i) = hopp(j,i) + t2KA
                                            else if (Species(i).eq.2) then
                                                 hopp(j,i) = hopp(j,i) + t2KB
                                            end if
                                        else if (layerIndex(i).eq.2 .or. layerIndex(i).eq.4) then
                                            if (Species(i).eq.2) then
                                                 hopp(j,i) = hopp(j,i) + t2KA
                                            else if (Species(i).eq.1) then
                                                 hopp(j,i) = hopp(j,i) + t2KB
                                            end if
                                        end if
                                   end if
                              else if (dist .gt. 2.0_dp*aCC*0.9_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                   !if (i.eq.1) print*, "assigning intralayer F2G2 terms 3"
                                   if (threeLayers .or. fourLayersSandwiched .or. helicalTwistedMBM .or. fiveLayersSandwiched .or. sixLayersSandwiched .or. sevenLayersSandwiched .or. eightLayersSandwiched .or. tenLayersSandwiched .or. twentyLayersSandwiched) then
                                        if ((layerIndex(i).eq.3) .or. (middleTwist)) then
                                           hopp(j,i) = hopp(j,i) + t3KSL
                                        else
                                           if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                                hopp(j,i) = hopp(j,i) + t3K
                                           else if (Species(i).eq.2) then
                                                hopp(j,i) = hopp(j,i) + t3K
                                           end if
                                        end if
                                   else if (twoLayers .or. oneLayer) then
                                        if (forceBilayerF2G2Intralayer) then
                                           if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                                hopp(j,i) = hopp(j,i) + t3K
                                           else if (Species(i).eq.2) then
                                                hopp(j,i) = hopp(j,i) + t3K
                                           end if
                                        else
                                           hopp(j,i) = hopp(j,i) + t3KSL
                                        end if
                                   else ! fourLayers
                                        if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                             hopp(j,i) = hopp(j,i) + t3K
                                        else if (Species(i).eq.2) then
                                             hopp(j,i) = hopp(j,i) + t3K
                                        end if
                                   end if
                              else if (dist .gt. 2.0_dp*aCC*1.1_dp .and. dist .lt. 4.0_dp) then
                                   !if (i.eq.1) print*, "assigning intralayer F2G2 terms 4"
                                   if (threeLayers .or. fourLayersSandwiched .or. helicalTwistedMBM .or. fiveLayersSandwiched .or. sixLayersSandwiched .or. sevenLayersSandwiched .or. eightLayersSandwiched .or. tenLayersSandwiched .or. twentyLayersSandwiched) then
                                        if ((layerIndex(i).eq.3) .or. (middleTwist)) then
                                           hopp(j,i) = hopp(j,i) + t4KSL
                                        else
                                           if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                                hopp(j,i) = hopp(j,i) + t4K
                                           else if (Species(i).eq.2) then
                                                hopp(j,i) = hopp(j,i) + t4K
                                           end if
                                        end if
                                   else if (twoLayers .or. oneLayer) then
                                        if (forceBilayerF2G2Intralayer) then
                                           if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                                hopp(j,i) = hopp(j,i) + t4K
                                           else if (Species(i).eq.2) then
                                                hopp(j,i) = hopp(j,i) + t4K
                                           end if
                                        else
                                           hopp(j,i) = hopp(j,i) + t4KSL
                                        end if
                                   else ! fourLayers
                                        if (Species(i).eq.1) then ! my A sublattice seems to have index 2
                                             hopp(j,i) = hopp(j,i) + t4K
                                        else if (Species(i).eq.2) then
                                             hopp(j,i) = hopp(j,i) + t4K
                                        end if
                                   end if
                              else if (dist .gt. (aCC*3.0_dp)*0.9_dp .and. dist .lt. (aCC*3.0_dp)*1.1_dp) then
                                   !if (i.eq.1) print*, "assigning intralayer F2G2 terms 5", delta
                                   if (threeLayers .or. fourLayersSandwiched .or. helicalTwistedMBM .or. fiveLayersSandwiched .or. sixLayersSandwiched .or. sevenLayersSandwiched .or. eightLayersSandwiched .or. tenLayersSandwiched .or. twentyLayersSandwiched) then
                                        if ((layerIndex(i).eq.3) .or. (middleTwist)) then
                                           hopp(j,i) = hopp(j,i) + t5KSL
                                        else
                                           if (layerIndex(i).eq.1) then
                                               if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t5KA
                                               else if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t5KB
                                               end if
                                           else if (layerIndex(i).eq.2) then
                                               if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t5KA
                                               else if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t5KB
                                               end if
                                           end if
                                        end if
                                   else if (twoLayers .or. oneLayer) then
                                        if (forceBilayerF2G2Intralayer) then
                                           if (layerIndex(i).eq.1) then
                                               if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t5KA
                                               else if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t5KB
                                               end if
                                           else if (layerIndex(i).eq.2) then
                                               if (Species(i).eq.2) then
                                                    hopp(j,i) = hopp(j,i) + t5KA
                                               else if (Species(i).eq.1) then
                                                    hopp(j,i) = hopp(j,i) + t5KB
                                               end if
                                           end if
                                        else
                                           hopp(j,i) = hopp(j,i) + t5KSL
                                        end if
                                   else ! four layers
                                        if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3) then
                                            if (Species(i).eq.1) then
                                                 hopp(j,i) = hopp(j,i) + t5KA
                                            else if (Species(i).eq.2) then
                                                 hopp(j,i) = hopp(j,i) + t5KB
                                            end if
                                        else if (layerIndex(i).eq.2 .or. layerIndex(i).eq.4) then
                                            if (Species(i).eq.2) then
                                                 hopp(j,i) = hopp(j,i) + t5KA
                                            else if (Species(i).eq.1) then
                                                 hopp(j,i) = hopp(j,i) + t5KB
                                            end if
                                        end if
                                        if (Species(i).eq.1) then
                                             hopp(j,i) = hopp(j,i) + t5KA
                                        else if (Species(i).eq.2) then
                                             hopp(j,i) = hopp(j,i) + t5KB
                                        end if
                                   end if
                              else if (dist .gt. 2.0_dp*aG*0.97_dp .and. dist .lt. 2.0_dp*aG*1.03_dp) then ! narrow: the seventh shell is only 4 % further out
                                   !if (i.eq.1) print*, "assigning intralayer F2G2 terms 6", delta
                                   hopp(j,i) = hopp(j,i) + t6K
                              else if (dist .gt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *0.9_dp .and. dist .lt. dsqrt((2.0_dp*aG)**2.0_dp + (aCC**2.0_dp)) *1.1_dp) then
                                   !if (i.eq.1) print*, "assigning intralayer F2G2 terms 7", delta
                                   hopp(j,i) = hopp(j,i) + t7K
                              else if (dist .gt. 4.0_dp*aCC*0.9_dp .and. dist .lt. 4.0_dp*aCC*1.1_dp) then
                                   !if (i.eq.1) print*, "assigning intralayer F2G2 terms 8", delta
                                   hopp(j,i) = hopp(j,i) + t8K
                              else
                                   !if (i.eq.1) print*, "assigning intralayer F2G2 terms beyond", dist
                                   hopp(j,i) = hopp(j,i) + 0.0_dp
                              end if
                           end if
                           if (d > 1.5_dp) then
                               delta = 0.0_dp
                           end if
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp)
                        else if (MIO_StringComp(BilayerModel,'Srivani')) then
                           !print*, "adding Srivani in plane terms"
                           !if (frac) call AtomsSetCart()
                           !dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                           !rbar = dist/aG
                           !theta = atan(NeighD(2,j,i)/NeighD(1,j,i))
                           !V0 = lambda0 * exp(-xi0*rbar**2.0_dp) * cos(kappa0 * rbar)
                           !V3 = lambda3 * rbar**2.0_dp * exp(-xi3*(rbar-x3)**2.0_dp)
                           !V6 = lambda6 * exp(-xi6*(rbar-x6)**2.0_dp) * sin(kappa6 * rbar)
                           !hopp(j,i) = V0 + V3 * cos(3.0_dp*theta) + V6*cos(6.0_dp*theta)
                           !hopp(j,i) = hopp(j,i)
                           !hopp(j,i) = cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta )
                           dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                           aCC = aG/sqrt(3.0_dp)
                           if (dist .gt. 0.1_dp .and. dist .lt. aCC*1.1_dp) then
                                hopp(j,i) = hopp(j,i) + t1K
                           else if (dist .gt. aCC*1.1_dp .and. dist .lt. aG*1.1_dp) then
                                hopp(j,i) = hopp(j,i) + t2K
                           else if (dist .gt. aG*1.1_dp .and. dist .lt. aCC*2.0_dp*1.1_dp) then
                                hopp(j,i) = hopp(j,i) + t3K
                           else
                                hopp(j,i) = hopp(j,i) + 0.0_dp
                           end if
                        else if (MIO_StringComp(BilayerModel,'Jeil')) then
                           if (d > 1.5_dp) then
                               delta = 0.0_dp
                           end if
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) + 1.0_dp*delta ,kind=dp)
                        else
                           !print*, i, j, cmplx(gn(Species(i),Species(NList(j,i)),ilvl))
                           if (d > 1.5_dp) then
                               delta = 0.0_dp
                           end if
                           hopp(j,i) = hopp(j,i) + cmplx(gn(Species(i),Species(NList(j,i)),ilvl) - 1.0_dp*delta,kind=dp)
                        end if
                        exit
                     end if
                  end do
                  !if ((Species(i)==1 .and. Species(NList(j,i))==2) .or. &
                  !  (Species(i)==2 .and. Species(NList(j,i))==1)) then
                  !   hopp(j,i) = cmplx(1.0_dp)
                  !else if ((Species(i)==3 .and. Species(NList(j,i))==4) .or. &
                  !  (Species(i)==4 .and. Species(NList(j,i))==3)) then
                  !   hopp(j,i) = cmplx(gBN0)
                  !end if
               else if (((layerIndex(i) .eq. layerIndex(jj) + 2 .or. layerIndex(i).eq. layerIndex(jj) - 2)) .and. encapsulatedSevenLayers .and. layerIndex(i) .ne. 1 .and. layerIndex(i) .ne. 7 .and. layerIndex(jj).ne.1 .and. layerIndex(jj).ne. 7) then ! next nearest interlayer
                   dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                   if (dist < 0.1_dp) then
                      hopp(j,i) = hopp(j,i) + tA1B3 ! t2 hopping term
                   else if (dist>aCC*0.9_dp .and. dist<aCC*1.1_dp) then ! First nearest neighbors surrounding the atom in the middle of the hexagon
                      hopp(j,i) = hopp(j,i) + tB1A3
                   else if (dist>aG*0.9_dp .and. dist<aG*1.1_dp) then ! g1 term
                      hopp(j,i) = hopp(j,i) + tA1B3_1
                   else if (dist>aCC*2.0_dp*0.9_dp .and. dist<aCC*2.0_dp*1.1_dp) then ! f2 term
                      hopp(j,i) = hopp(j,i) + tB1A3_2
                   else if (dist>(aCC*3.0_dp)*0.9_dp .and. dist<(aCC*3.0_dp)*1.1_dp) then ! g2 term
                      hopp(j,i) = hopp(j,i) + tA1B3_2
                   else
                      hopp(j,i) = hopp(j,i) + 0.0_dp
                   end if
               else if (((layerIndex(i) .eq. layerIndex(jj) + 2 .or. layerIndex(i).eq. layerIndex(jj) - 2)) .and. encapsulatedSixLayers .and. layerIndex(i) .ne. 1 .and. layerIndex(i) .ne. 6 .and. layerIndex(jj).ne.1 .and. layerIndex(jj).ne. 6) then ! next nearest interlayer
                   dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                   if (dist < 0.1_dp) then
                      hopp(j,i) = hopp(j,i) + tA1B3 ! t2 hopping term
                   else if (dist>aCC*0.9_dp .and. dist<aCC*1.1_dp) then ! First nearest neighbors surrounding the atom in the middle of the hexagon
                      hopp(j,i) = hopp(j,i) + tB1A3
                   else if (dist>aG*0.9_dp .and. dist<aG*1.1_dp) then ! g1 term
                      hopp(j,i) = hopp(j,i) + tA1B3_1
                   else if (dist>aCC*2.0_dp*0.9_dp .and. dist<aCC*2.0_dp*1.1_dp) then ! f2 term
                      hopp(j,i) = hopp(j,i) + tB1A3_2
                   else if (dist>(aCC*3.0_dp)*0.9_dp .and. dist<(aCC*3.0_dp)*1.1_dp) then ! g2 term
                      hopp(j,i) = hopp(j,i) + tA1B3_2
                   else
                      hopp(j,i) = hopp(j,i) + 0.0_dp
                   end if
               else if (((layerIndex(i) .eq. layerIndex(jj) + 2 .or. layerIndex(i).eq. layerIndex(jj) - 2)) .and. encapsulatedFiveLayers .and. layerIndex(i) .ne. 1 .and. layerIndex(i) .ne. 5 .and. layerIndex(jj).ne.1 .and. layerIndex(jj).ne. 5) then ! next nearest interlayer
                   dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                   if (dist < 0.1_dp) then
                      hopp(j,i) = hopp(j,i) + tA1B3 ! t2 hopping term
                   else if (dist>aCC*0.9_dp .and. dist<aCC*1.1_dp) then ! First nearest neighbors surrounding the atom in the middle of the hexagon
                      hopp(j,i) = hopp(j,i) + tB1A3
                   else if (dist>aG*0.9_dp .and. dist<aG*1.1_dp) then ! g1 term
                      hopp(j,i) = hopp(j,i) + tA1B3_1
                   else if (dist>aCC*2.0_dp*0.9_dp .and. dist<aCC*2.0_dp*1.1_dp) then ! f2 term
                      hopp(j,i) = hopp(j,i) + tB1A3_2
                   else if (dist>(aCC*3.0_dp)*0.9_dp .and. dist<(aCC*3.0_dp)*1.1_dp) then ! g2 term
                      hopp(j,i) = hopp(j,i) + tA1B3_2
                   else
                      hopp(j,i) = hopp(j,i) + 0.0_dp
                   end if
               else if (layerIndex(i) .eq. layerIndex(jj) + 1 .or. layerIndex(i).eq. layerIndex(jj)-1  .or. layerIndex(i)+nLayers .eq. layerIndex(jj) + 1 .or. layerIndex(i)-nLayers .eq. layerIndex(jj)-1) then ! INTERLAYER ! nLayers is to account for bulk, make sure this variable is correctly defined
                  !print*, "doing interlayer for ", i, jj, nLayers
                  if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell') .or. BernalReadXYZ) then
                    if (frac) call AtomsSetCart()
                    dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                    if (abs(NeighD(1,j,i)) < 0.01_dp .and. abs(NeighD(2,j,i)) < 0.01_dp) then
                       hopp(j,i) = hopp(j,i) -tAB1
                       numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                       print*, "here1"
                       !hopp(j,i) = interlayerBL(HBL, dx, dy, tAB)
                    else if (dist>aCC*0.95_dp .and. dist<aCC*1.05_dp .and. BilayerThreeParameters) then
                       if (Species(i).eq.Species(NList(j,i))) then
                          hopp(j,i) = hopp(j,i) -tAB4
                          numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                       else
                          hopp(j,i) = hopp(j,i) -tAB3
                          numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                       end if
                    else if (dist>aCC*0.95_dp*2.0_dp .and. dist<aCC*1.05_dp*2.0_dp .and. BilayerThreeParameters) then
                       if (Species(i).eq.Species(NList(j,i))) then
                          hopp(j,i) = hopp(j,i) -tAB4
                          numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                       else
                          hopp(j,i) = hopp(j,i) -tAB3
                          numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                       end if
                    else
                       hopp(j,i) = hopp(j,i) + 0.0_dp
                    end if
                  else if (MIO_StringComp(str,'TrilayerBasedOnMoireCell')) then
                    !!print*, "doing the trilayer for ", str
                    if (frac) call AtomsSetCart()
                    distXY = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                    distXYZ = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp + NeighD(3,j,i)**2.0_dp)
                    if (distXY < 0.1_dp .and. abs(NeighD(3,j,i)) < 3.5_dp) then
                        hopp(j,i) = hopp(j,i) + tr1
                        numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                        print*, "here2"
                    else if (distXY < 0.1_dp .and. abs(NeighD(3,j,i)) > 3.5_dp) then
                        hopp(j,i) = hopp(j,i) + tr2/2.0_dp
                        numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                    else if (NeighD(3,j,i) < 3.5_dp) then
                       if (Species(i).eq.Species(NList(j,i))) then
                         hopp(j,i) = hopp(j,i) + tr4
                         numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                       else
                         hopp(j,i) = hopp(j,i) + tr3
                         numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                       end if
                    !if (abs(NeighD(1,j,i)) < 0.01_dp .and. abs(NeighD(2,j,i)) < 0.01_dp) then
                    !   hopp(j,i) = -tAB1
                    !   numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                    !   !hopp(j,i) = interlayerBL(HBL, dx, dy, tAB)
                    !else if (dist>aCC*0.95_dp .and. dist<aCC*1.05_dp .and. BilayerThreeParameters) then
                    !   if (Species(i).eq.Species(NList(j,i))) then
                    !      hopp(j,i) = -tAB4
                    !      numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                    !   else
                    !      hopp(j,i) = -tAB3
                    !      numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                    !   end if
                    !else if (dist>aCC*0.95_dp*2.0_dp .and. dist<aCC*1.05_dp*2.0_dp .and. BilayerThreeParameters) then
                    !   if (Species(i).eq.Species(NList(j,i))) then
                    !      hopp(j,i) = -tAB4
                    !      numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                    !   else
                    !      hopp(j,i) = -tAB3
                    !      numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                    !   end if
                    else
                       hopp(j,i) = hopp(j,i)+ 0.0_dp
                    end if
                    hopp(j,i) = -hopp(j,i)
                    print*, "if this shows up, check the code... this might be wrong if you have a moire lattice... the moire part already has the sign change"
                  else if (MIO_StringComp(str,'TwistedBilayer') .or. (MIO_StringComp(str,'ReadXYZ') .and. .not.(singleLayerXYZ)) .or. MIO_StringComp(str,'TwistedBilayerBasedOnMoireCell')) then ! add intralayer hopping terms (see Jeil's paper)
                   if(MIO_StringComp(BilayerModel,'Jeil')) then
                      !if (jj.eq.0) then
                      !   hopp(j,i) = 0.0_dp
                      !else if ((j.eq.jj .and. Species(j)==Species(jj))  &
                      !  .or. ((j.eq.firstNN(1) .or. j.eq.firstNN(2) .or. j.eq.firstNN(3)) .and. Species(j).ne.Species(firstNN(1)))) then ! only add closest Neighbor
                      !if ((j.eq.jj1) .or. (j.eq.jj2)  &
                      !  .or. (j.eq.firstNN(1) .or. j.eq.firstNN(2) .or. j.eq.firstNN(3))) then ! only add closest Neighbor
                         !eps = 0.0_dp
                         !dxi = ((1.0_dp+eps) * cos(twistedBilayerAngleGrad) - 1.0_dp) * Rat(1,i) - ((1.0_dp+eps) * sin(twistedBilayerAngleGrad) * Rat(2,i))
                         !dyi = ((1.0_dp+eps) * cos(twistedBilayerAngleGrad) - 1.0_dp) * Rat(2,i) + ((1.0_dp+eps) * sin(twistedBilayerAngleGrad) * Rat(1,i))
                         !dxj = ((1.0_dp+eps) * cos(twistedBilayerAngleGrad) - 1.0_dp) * Rat(1,NList(j,i)) - ((1.0_dp+eps) * sin(twistedBilayerAngleGrad) * Rat(2,NList(j,i)))
                         !dyj = ((1.0_dp+eps) * cos(twistedBilayerAngleGrad) - 1.0_dp) * Rat(2,NList(j,i)) + ((1.0_dp+eps) * sin(twistedBilayerAngleGrad) * Rat(1,NList(j,i)))
                         !dx = Rat(1,i) * ((1.0_dp+eps) * cos(twistAngleGrad) - 1.0_dp) - Rat(2,i) * (1.0_dp+eps) * sin(twistAngleGrad)
                         !dy = Rat(2,i) * ((1.0_dp+eps) * cos(twistAngleGrad) - 1.0_dp) + Rat(1,i) * (1.0_dp+eps) * sin(twistAngleGrad)
                         !print*, twistedBilayerAngleGrad, Rat(1,i), Rat(2,i)
                         dxi = ((1.0_dp) * cos(twistedBilayerAngleGrad) - 1.0_dp) * Rat(1,i) &
                                     - ((1.0_dp) * sin(twistedBilayerAngleGrad) * Rat(2,i)) + xShift
                         dyi = ((1.0_dp) * cos(twistedBilayerAngleGrad) - 1.0_dp) * Rat(2,i) &
                                     + ((1.0_dp) * sin(twistedBilayerAngleGrad) * Rat(1,i)) + yShift
                         dxj = ((1.0_dp) * cos(twistedBilayerAngleGrad) - 1.0_dp) * Rat(1,NList(j,i)) &
                                     - ((1.0_dp) * sin(twistedBilayerAngleGrad) * Rat(2,NList(j,i))) + xShift
                         dyj = ((1.0_dp) * cos(twistedBilayerAngleGrad) - 1.0_dp) * Rat(2,NList(j,i)) &
                                     + ((1.0_dp) * sin(twistedBilayerAngleGrad) * Rat(1,NList(j,i))) + yShift
                         !dxi = - twistedBilayerAngleGrad * Rat(2,i)
                         !dyi =  twistedBilayerAngleGrad * Rat(1,i)
                         !dxj = - twistedBilayerAngleGrad * Rat(2,NList(j,i))
                         !dyj =  twistedBilayerAngleGrad * Rat(1,NList(j,i))
                         !call MIO_Print(' d_i '//trim(num2str(dxi))//' '//trim(num2str(dyi)),'HamHopping')
                         !call MIO_Print(' d_j '//trim(num2str(dxj))//' '//trim(num2str(dyj)),'HamHopping')
                         !if(i==1) then
                         !   print*, layerIndex(i)
                         !end if
                     !    print*, "layerindices", layerIndex(i), layerIndex(NList(j,i))
                         dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                         expFactor = exp(-dist/10.0_dp) * ((NeighD(3,j,i))/(dist))**2.0_dp
                         !expFactor = 1.0_dp * ((NeighD(3,j,i))/(dist))**2.0_dp
                         if (Species(i)==Species(NList(j,i))) then
                              !either AA' or BB' (Habjj,dx,dy,Cab,Phiab
                              if (layerIndex(i)==1) then
                                  posOrNeg = 1.0_dp
                                  call interlayerBLAA(HAA,dxi,dyi,tbt,posOrNeg)
                                  hopp(j,i) =  hopp(j,i) + HAA*expFactor
                                  !if(NList(j,i)==1) then
                                  !    print*, "HAA1 = ", HAA, i
                                  !end if
                                  numberOfHAA1 = numberOfHAA1 + 1
                              else if(layerIndex(i)==2) then
                                  posOrNeg = -1.0_dp
                                  call interlayerBLAA(HAA,dxj,dyj,tbt,posOrNeg)
                                  hopp(j,i) = hopp(j,i) +  HAA*expFactor
                                  !if(i==1) then
                                  !    print*, "HAA2 = ", HAA, i
                                  !end if
                                  numberOfHAA2 = numberOfHAA2 + 1
                              !else
                              !    print*, "not attributed error1"
                              end if
                         else if (Species(i) < Species(NList(j,i))) then  ! i -> A, j -> B
                              if (layerIndex(i)==1) then
                                  ! AB'
                                  posOrNeg = 1.0_dp
                                  call interlayerBLAB(HAB,dxi,dyi,tbt,posOrNeg)
                                  hopp(j,i) = hopp(j,i) + HAB*expFactor
                                  !if(NList(j,i)==1) then
                                  !    print*, "HAB1 = ", HAB, i
                                  !end if
                                  numberOfHAB1 = numberOfHAB1 + 1
                              else if(layerIndex(i)==2) then
                                  ! BA'
                                  posOrNeg = -1.0_dp
                                  call interlayerBLBA(HBA,dxj,dyj,tbt,posOrNeg)
                                  hopp(j,i) = hopp(j,i) +  HBA*expFactor
                                  !if (i==1) then
                                  !    print*, "HBA1 = ", HBA, i
                                  !end if
                                  numberOfHBA2 = numberOfHBA2 + 1
                              !else
                              !    print*, "not attributed error2"
                              end if
                         else if (Species(i) > Species(NList(j,i))) then  ! i -> B, j -> A
                              if (layerIndex(i)==1) then
                                  ! BA'
                                  posOrNeg = 1.0_dp
                                  call interlayerBLBA(HBA,dxi,dyi,tbt,posOrNeg)
                                  hopp(j,i) = hopp(j,i) +  HBA*expFactor
                                  !if(NList(j,i)==1) then
                                  !     print*, "HBA2 = ", HBA, i
                                  !end if
                                  numberOfHBA1 = numberOfHBA1 + 1
                              else if (layerIndex(i)==2) then
                                  ! AB'
                                  posOrNeg = -1.0_dp
                                  call interlayerBLAB(HAB,dxj,dyj,tbt,posOrNeg)
                                  hopp(j,i) = hopp(j,i) + HAB*expFactor
                                  !if (i==1) then
                                  !     print*, "HAB2 = ", HAB, i
                                  !end if
                                  numberOfHAB2 = numberOfHAB2 + 1
                              !else
                              !    print*, "not attributed error3"
                              end if
                         !else
                         !    print*, "not attributed error4"
                         end if
                      !else
                      !   hopp(j,i) = 0.0_dp
                      !end if
                   else if (GBNtwoLayers .or. encapsulatedThreeLayers .or. encapsulatedFourLayers .or. encapsulatedFiveLayers .or. encapsulatedSixLayers .or. encapsulatedSevenLayers .or. t3GwithBN .or. BNt2GBN .or. t2GBN) then ! Adds the interlayer TC model
                      if (deactivateInterlayert2GBN1to2 .and. ((layerIndex(i).eq. 1 .and. layerIndex(jj).eq.2) .or. (layerIndex(i).eq.2 .and. layerIndex(jj).eq.1))) then
                          cycle
                      else if (deactivateInterlayert2GBN2to3 .and. ((layerIndex(i).eq. 2 .and. layerIndex(jj).eq.3) .or. (layerIndex(i).eq.3 .and. layerIndex(jj).eq.2))) then
                          cycle
                      else if (deactivateInterlayer) then
                          cycle
                      end if
                      if (frac) call AtomsSetCart()
                      if (threeLayerShort) then
                          dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                          if (dist<0.2_dp) then ! on top of each other
                             hopp(j,i) = hopp(j,i) + tB1A2
                             numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             print*, "here3"
                          else if (dist>aCC*0.9_dp .and. dist<aCC*1.1_dp) then ! First nearest neighbors surrounding the atom in the middle of the hexagon
                             if (Species(i).eq.Species(NList(j,i))) then ! same sublattice
                                hopp(j,i) = hopp(j,i) + tA1A2
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             else ! Different sublattice
                                hopp(j,i) = hopp(j,i) + tA1B2
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             end if
                          else if (dist>aG*0.9_dp .and. dist<aG*1.1_dp) then ! g1 term
                             if (Species(i).eq.Species(NList(j,i))) then
                                hopp(j,i) = hopp(j,i) + tA1A1_1
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             else
                                hopp(j,i) = hopp(j,i) + tB1A2_1
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             end if
                          else if (dist>aCC*2.0_dp*0.9_dp .and. dist<aCC*2.0_dp*1.1_dp) then ! f2 term
                             if (Species(i).eq.Species(NList(j,i))) then
                                hopp(j,i) = hopp(j,i) + tA1A2_2
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             else
                                hopp(j,i) = hopp(j,i) + tA1B2_2
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             end if
                          else if (dist>(aCC*3.0_dp)*0.9_dp .and. dist<(aCC*3.0_dp)*1.1_dp) then ! g2 term
                             if (Species(i).ne.Species(NList(j,i))) then
                                hopp(j,i) = hopp(j,i) + tA1A1_2
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             else
                                hopp(j,i) = hopp(j,i) + tB1A2_2
                             end if
                          else
                             hopp(j,i) = hopp(j,i) + 0.0_dp
                          end if
                      else if (((layerIndex(i) .eq. 2 .and. layerIndex(jj) .eq. 3) .or. (layerIndex(i) .eq. 3 .and. layerIndex(jj) .eq. 2)) .and. encapsulatedFourLayersF2G2) then
                          dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                          if (dist<0.2_dp) then
                             hopp(j,i) = hopp(j,i) + tBA0
                             numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             print*, "here3"
                          else if (dist>aCC*0.9_dp .and. dist<aCC*1.1_dp) then ! First nearest neighbors surrounding the atom in the middle of the hexagon
                             if (Species(i).eq.Species(NList(j,i))) then ! same sublattice
                                hopp(j,i) = hopp(j,i) + tAA1
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             else ! Different sublattice
                                hopp(j,i) = hopp(j,i) + tAB1
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             end if
                          else if (dist>aG*0.9_dp .and. dist<aG*1.1_dp) then ! BA' (corresponds to g1, not g2)
                             if (Species(i).eq.Species(NList(j,i))) then
                                hopp(j,i) = hopp(j,i) + tBA2
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             else
                                hopp(j,i) = hopp(j,i) + tBA2
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             end if
                          else if (dist>aCC*2.0_dp*0.9_dp .and. dist<aCC*2.0_dp*1.1_dp) then ! Second neirest neighbor from central atom in middle of hexagon
                             if (Species(i).eq.Species(NList(j,i))) then
                                hopp(j,i) = hopp(j,i) + tAA2
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             else
                                hopp(j,i) = hopp(j,i) + tAB2
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             end if
                          else if (dist>(aCC*3.0_dp)*0.9_dp .and. dist<(aCC*3.0_dp)*1.1_dp) then ! BA' (corresponds to g2, not g5)
                             if (Species(i).ne.Species(NList(j,i))) then
                                hopp(j,i) = hopp(j,i) + tBA5
                                numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                             else
                                hopp(j,i) = hopp(j,i)
                             end if
                          else
                             hopp(j,i) = hopp(j,i) + 0.0_dp
                          end if
                      else if ((layerIndex(i) .eq. layerIndex(jj) + 1 .or. layerIndex(i).eq. layerIndex(jj)-1)) then ! INTERLAYER
                          if (deactivateInterlayerBG .and. ((layerIndex(i).eq.2 .and. layerIndex(jj).eq.3) .or. (layerIndex(i).eq.3 .and. layerIndex(jj).eq.2))) then
                              cycle
                          else if (deactivateInterlayer) then
                              cycle
                          else if (deactivateInterlayer12 .and. ((layerIndex(i).eq.1 .and. layerIndex(jj).eq.2) .or. (layerIndex(i).eq.2 .and. layerIndex(jj).eq.1))) then
                              cycle
                          else if (deactivateInterlayer23 .and. ((layerIndex(i).eq.2 .and. layerIndex(jj).eq.3) .or. (layerIndex(i).eq.3 .and. layerIndex(jj).eq.2))) then
                              cycle
                          else if (deactivateInterlayer34 .and. ((layerIndex(i).eq.3 .and. layerIndex(jj).eq.4) .or. (layerIndex(i).eq.4 .and. layerIndex(jj).eq.3))) then
                              cycle
                          end if
                          if (deactivateIntraSublatticeForC .and. (Species(i).lt.3 .and. Species(jj).lt.3) .and. (Species(i).eq.Species(jj))) then
                              cycle
                          end if
                          ! add the hopping terms depending if it is the C-B, C-N or C-C
                          ! interaction
                          if (i.eq.1) call MIO_Print('Adding the interlayer terms for the HTC model for G and BN interactions','ham')
                          if (Species(i).eq.1 .and. Species(NList(j,i)).eq.3) then ! C to B
                            aval = 1.02796_dp
                            bval = 3.047657_dp
                          else if (Species(i).eq.3 .and. Species(NList(j,i)).eq.1) then ! C to B
                            aval = 1.02796_dp
                            bval = 3.047657_dp
                          else if (Species(i).eq.1 .and. Species(NList(j,i)).eq.4) then ! C to N
                            aval = 1.19218_dp
                            bval = 3.47553_dp
                          else if (Species(i).eq.4 .and. Species(NList(j,i)).eq.1) then ! C to N
                            aval = 1.19218_dp
                            bval = 3.47553_dp
                          else if (Species(i).eq.2 .and. Species(NList(j,i)).eq.3) then ! C to B
                            aval = 1.049223_dp
                            bval = 3.0646389_dp
                          else if (Species(i).eq.3 .and. Species(NList(j,i)).eq.2) then ! C to B
                            aval = 1.049223_dp
                            bval = 3.0646389_dp
                          else if (Species(i).eq.2 .and. Species(NList(j,i)).eq.4) then ! C to N
                            aval = 1.16223_dp
                            bval = 3.4799_dp
                          else if (Species(i).eq.4 .and. Species(NList(j,i)).eq.2) then ! C to N
                            aval = 1.16223_dp
                            bval = 3.4799_dp
                          else if (Species(i).eq.1 .and. Species(NList(j,i)).eq.2) then ! C to C
                            aval = 1.34_dp
                            bval = 3.25_dp
                          else if (Species(i).eq.2 .and. Species(NList(j,i)).eq.1) then ! C to C
                            aval = 1.34_dp
                            bval = 3.25_dp
                          else if (Species(i).eq.1 .and. Species(NList(j,i)).eq.1) then ! C to C
                            aval = 1.34_dp
                            bval = 3.25_dp
                          else if (Species(i).eq.2 .and. Species(NList(j,i)).eq.2) then ! C to C
                            aval = 1.34_dp
                            bval = 3.25_dp
                          end if
                          renormalizeHoppingFactorAAp  = exp((abs(NeighD(3,j,i))-bval)/aval)
                          !renormalizeHoppingFactorAAp  = 1.0_dp
                          dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                          if ((Species(i).eq.1 .and. Species(jj) .eq. 2) .or. (Species(i).eq.2 .and. Species(jj) .eq. 1) .or. (Species(i).eq.1 .and. Species(jj) .eq. 1) .or. (Species(i).eq.2 .and. Species(jj) .eq. 2)) then
                             BLdelta = 0.184_dp*aG
                             vppsigma = vppsigma0 * exp(-(dist-interlayerdistance)/(BLdelta)) * renormalizeHoppingFactorAAp
                             vpppi = vpppi0 * exp(-(dist-aCC)/(BLdelta)) * renormalizeHoppingFactorAAp ! inplane
                          else
                             BLdelta = 0.184_dp*2.438977_dp
                             vppsigma = vppsigma0 * exp(-(dist-interlayerdistance)/(BLdelta)) * renormalizeHoppingFactorAAp
                             vpppi = vpppi0 * exp(-(dist-2.438977_dp/sqrt(3.0))/(BLdelta)) * renormalizeHoppingFactorAAp
                          end if
                          hopp(j,i) = hopp(j,i) + interlayerTwoCenter(vpppi, vppsigma, NeighD(:,j,i), dist, &
                                               pzNormal(:,i), pzNormal(:,NList2(j,i)), corrugatedInterlayerTwoCenter)
                      end if
                          !Rz = R(3);
                          !Rd = sqrt(R(1)^2+R(2)^2+R(3)^2);
                          !
                          !Vpppi0 = -2.7;
                          !Vppsigma0 = 0.48;
                          !
                          !a = 2.438977765130281;
                          !a0 = a/sqrt(3);
                          !d0 = 3.35;
                          !r0 = 0.184*a;
                          !
                          !Vpppi = Vpppi0*exp(-(Rd-a0)/r0);
                          !Vppsigma = Vppsigma0*exp(-(Rd-d0)/r0);
                          !
                          !t = Vpppi*(1-(Rz/Rd)^2)+Vppsigma*((Rz/Rd)^2);
                      !print*, "TC", hopp(j,i)
                   !else if (BNBNtwoLayers) then
                   !   if (deactivateInterlayer) then
                   !       cycle
                   !   end if
                   !   if (Species(i).eq.3 .and. Species(NList(j,i).eq.3) then ! C to B
                   !     renormalizeHoppingFactorAAp  = exp((abs(NeighD(3,j,i))-bvalCB)/avalCB)
                   !   else if (Species(i).eq.3 .and. Species(NList(j,i).eq.4) then ! C to N
                   !     renormalizeHoppingFactorAAp  = exp((abs(NeighD(3,j,i))-bvalCN)/avalCN)
                   !   else if (Species(i).eq.4 .and. Species(NList(j,i).eq.3) then ! C to B
                   !     renormalizeHoppingFactorAAp = 1.0
                   !   else if (Species(i).eq.4 .and. Species(NList(j,i).eq.4) then ! C to N
                   !     renormalizeHoppingFactorAAp = 1.0
                   !   end if
                   !   vppsigma = vppsigma0 * exp(-(dist-interlayerdistance)/(BLdelta)) * renormalizeHoppingFactorAAp
                   !   vpppi = vpppi0 * exp(-(dist-aCC)/(BLdelta)) * renormalizeHoppingFactorAAp
                   !   hopp(j,i) = hopp(j,i) + (vpppi * (1.0_dp - ((NeighD(3,j,i))/(dist))**2.0_dp) &
                   !                        + vppsigma * ((NeighD(3,j,i))/(dist))**2.0_dp)
                   else if (BNBNtwoLayers) then
                      if ((layerIndex(i) .eq. layerIndex(jj) + 1 .or. layerIndex(i).eq. layerIndex(jj)-1)) then ! INTERLAYER
                          !if (deactivateInterlayerBG .and. ((layerIndex(i).eq.2 .and. layerIndex(jj).eq.3) .or. (layerIndex(i).eq.3 .and. layerIndex(jj).eq.2))) then
                          !    cycle
                          !else if (deactivateInterlayer) then
                          !    cycle
                          !else if (deactivateInterlayer12 .and. ((layerIndex(i).eq.1 .and. layerIndex(jj).eq.2) .or. (layerIndex(i).eq.2 .and. layerIndex(jj).eq.1))) then
                          !    cycle
                          !else if (deactivateInterlayer23 .and. ((layerIndex(i).eq.2 .and. layerIndex(jj).eq.3) .or. (layerIndex(i).eq.3 .and. layerIndex(jj).eq.2))) then
                          !    cycle
                          !else if (deactivateInterlayer34 .and. ((layerIndex(i).eq.3 .and. layerIndex(jj).eq.4) .or. (layerIndex(i).eq.4 .and. layerIndex(jj).eq.3))) then
                          !    cycle
                          !end if
                          ! add the hopping terms depending if it is the C-B, C-N or C-C
                          ! interaction
                          if (i.eq.1) call MIO_Print('Adding the interlayer terms for the HTC model for BN and BN interactions based on Fengpings parametrization','ham')
                          !if (Species(i).eq.1 .and. Species(NList(j,i)).eq.3) then ! C to B
                          !  aval = 1.02796_dp
                          !  bval = 3.047657_dp
                      end if
                   else if (MIO_StringComp(BilayerModel,'Koshino').or. MIO_StringComp(BilayerModel,'HTC') .or. MIO_StringComp(SinglelayerModel,'HTC')) then
                      if (i.eq.1 .and. j.eq.1) call MIO_Print('We add interlayer interactions for sandwiched systems, tBG, etc','ham')
                      if (deactivateInterlayer) then
                          cycle
                      else if (deactivateInterlayer12 .and. ((layerIndex(i).eq.1 .and. layerIndex(jj).eq.2) .or. (layerIndex(i).eq.2 .and. layerIndex(jj).eq.1))) then
                          cycle
                      else if (deactivateInterlayer23 .and. ((layerIndex(i).eq.2 .and. layerIndex(jj).eq.3) .or. (layerIndex(i).eq.3 .and. layerIndex(jj).eq.2))) then
                          cycle
                      else if (deactivateInterlayer34 .and. ((layerIndex(i).eq.3 .and. layerIndex(jj).eq.4) .or. (layerIndex(i).eq.4 .and. layerIndex(jj).eq.3))) then
                          cycle
                      else if (deactivateInterlayert2BG2to3 .and. ((layerIndex(i).eq. 2 .and. layerIndex(jj).eq.3) .or. (layerIndex(i).eq.3 .and. layerIndex(jj).eq.2))) then
                          cycle
                      else if (deactivateInterlayert3BG1to2 .and. ((layerIndex(i).eq. 1 .and. layerIndex(jj).eq.2) .or. (layerIndex(i).eq.2 .and. layerIndex(jj).eq.1))) then
                          print*, "deactivated layer t3BG1to2"
                          cycle
                      else if (deactivateInterlayert3BG2to3 .and. ((layerIndex(i).eq. 2 .and. layerIndex(jj).eq.3) .or. (layerIndex(i).eq.3 .and. layerIndex(jj).eq.2))) then
                          print*, "deactivated layer t3BG2to3"
                          cycle
                      else if (deactivateInterlayert3BG3to4 .and. ((layerIndex(i).eq. 3 .and. layerIndex(jj).eq.4) .or. (layerIndex(i).eq.4 .and. layerIndex(jj).eq.3))) then
                          print*, "deactivated layer t3BG3to4"
                          cycle
                      else if (deactivateInterlayert3BG4to5 .and. ((layerIndex(i).eq. 4 .and. layerIndex(jj).eq.5) .or. (layerIndex(i).eq.5 .and. layerIndex(jj).eq.4))) then
                          print*, "deactivated layer t3BG4to5"
                          cycle
                      else if (deactivateInterlayert3BG5to6 .and. ((layerIndex(i).eq. 5 .and. layerIndex(jj).eq.6) .or. (layerIndex(i).eq.6 .and. layerIndex(jj).eq.5))) then
                          print*, "deactivated layer t3BG5to6"
                          cycle
                      end if
                      !print*, "adding the out of plane Koshino terms"
                      if (frac) call AtomsSetCart()
                        if ((.not.((layerIndex(i).eq.2 .and. layerIndex(jj).eq.3) .or. (layerIndex(i).eq.3 .and. layerIndex(jj).eq.2) .or. (layerIndex(i).eq.5 .and. layerIndex(jj).eq.4) .or. (layerIndex(i).eq.4 .and. layerIndex(jj).eq.5)) ) .and. (fourLayers .or. t3BG)) then
                          if ((layerIndex(i) .eq. layerIndex(jj) + 1 .or. layerIndex(i).eq. layerIndex(jj)-1) .and. BilayerThreeParameters) then ! INTERLAYER
                            !print*, "doing interlayer for ", i, jj
                            !if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell')) then
                              !if (frac) call AtomsSetCart()
                              dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                              if (abs(NeighD(1,j,i)) < 0.2_dp .and. abs(NeighD(2,j,i)) < 0.2_dp) then
                                 hopp(j,i) = hopp(j,i) -tAB1
                                 numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 print*, "here4"
                                 !hopp(j,i) = interlayerBL(HBL, dx, dy, tAB)
                              else if (dist>aCC*0.95_dp .and. dist<aCC*1.05_dp .and. BilayerThreeParameters) then
                                 if (Species(i).eq.Species(NList(j,i))) then
                                    hopp(j,i) = hopp(j,i) -tAB4
                                    numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 else
                                    hopp(j,i) = hopp(j,i) -tAB3
                                    numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 end if
                              else if (dist>aCC*0.95_dp*2.0_dp .and. dist<aCC*1.05_dp*2.0_dp .and. BilayerThreeParameters) then
                                 if (Species(i).eq.Species(NList(j,i))) then
                                    hopp(j,i) = hopp(j,i) -tAB4
                                    numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 else
                                    hopp(j,i) = hopp(j,i) -tAB3
                                    numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 end if
                              else
                                 hopp(j,i) = hopp(j,i) + 0.0_dp
                              end if
                          else if ((layerIndex(i) .eq. layerIndex(jj) + 1 .or. layerIndex(i).eq. layerIndex(jj)-1) .and. BilayerOneParameter) then ! INTERLAYER
                            !print*, "doing interlayer for ", i, jj
                            !if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell')) then
                              !if (frac) call AtomsSetCart()
                              dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                              if (abs(NeighD(1,j,i)) < 0.2_dp .and. abs(NeighD(2,j,i)) < 0.2_dp) then
                                 hopp(j,i) = hopp(j,i) -tAB1
                                 numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 print*, "here5"
                                 !hopp(j,i) = interlayerBL(HBL, dx, dy, tAB)
                              else
                                 hopp(j,i) = hopp(j,i) + 0.0_dp
                              end if
                          else if ((layerIndex(i) .eq. layerIndex(jj) + 1 .or. layerIndex(i).eq. layerIndex(jj)-1) .and. F2G2Model) then ! INTERLAYER
                              dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                              !if (i.eq.1) then
                              !    print*, "we are applygin the F2G2 model for the interlayer terms in four layer system"
                              !    print*, "aG"
                              !    print*, dist
                              !end if
                              !if (abs(NeighD(1,j,i)) < 0.2_dp .and. abs(NeighD(2,j,i)) < 0.2_dp) then ! On top of each other
                              if (dist<0.2_dp) then
                                 !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 0", delta
                                 !print*, "2: assigning inerlayer F2G2 terms BA 0", i
                                 hopp(j,i) = hopp(j,i) + tBA0
                                 numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 !hopp(j,i) = interlayerBL(HBL, dx, dy, tAB)
                              else if (dist>aCC*0.9_dp .and. dist<aCC*1.1_dp) then ! First nearest neighbors surrounding the atom in the middle of the hexagon
                                 if (Species(i).eq.Species(NList(j,i))) then ! same sublattice
                                    !if (i.eq.1) print*, "assigning inerlayer F2G2 terms 1 AA", delta
                                    !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms 1 AA", delta
                                    hopp(j,i) = hopp(j,i) + tAA1
                                    numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 else ! Different sublattice
                                    !if (i.eq.1) print*, "assigning inerlayer F2G2 terms 1 AB", delta
                                    !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms 1 AB", delta
                                    hopp(j,i) = hopp(j,i) + tAB1
                                    numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 end if
                              else if (dist>aG*0.9_dp .and. dist<aG*1.1_dp) then ! BA' (corresponds to g1, not g2)
                                 if (Species(i).eq.Species(NList(j,i))) then
                                    !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 2", delta
                                    !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms BA 2", delta
                                    hopp(j,i) = hopp(j,i) + tBA2
                                    numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 else
                                    !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 2 bis", delta
                                    !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms BA 2 bis", delta
                                    hopp(j,i) = hopp(j,i) + tBA2
                                    numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 end if
                              else if (dist>aCC*2.0_dp*0.9_dp .and. dist<aCC*2.0_dp*1.1_dp) then ! Second neirest neighbor from central atom in middle of hexagon
                                 if (Species(i).eq.Species(NList(j,i))) then
                                    !if (i.eq.1) print*, "assigning inerlayer F2G2 terms AA 2 ", delta
                                    !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms AA 2", delta
                                    hopp(j,i) = hopp(j,i) + tAA2
                                    numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 else
                                    !if (i.eq.1) print*, "assigning inerlayer F2G2 terms AB 2 ", delta
                                    !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms AB 2", delta
                                    hopp(j,i) = hopp(j,i) + tAB2
                                    numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 end if
                              else if (dist>(aCC*3.0_dp)*0.9_dp .and. dist<(aCC*3.0_dp)*1.1_dp) then ! BA' (corresponds to g2, not g5)
                                 if (Species(i).ne.Species(NList(j,i))) then
                                    !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 5 ", delta
                                    !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms BA 5", delta
                                    hopp(j,i) = hopp(j,i) + tBA5
                                    numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                                 else
                                    !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 5 bis ", delta
                                    !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms BA 5 bis", delta
                                    hopp(j,i) = hopp(j,i)
                                 end if
                              else
                                 !if (i.eq.1) print*, "assigning interlayer F2G2 terms BA 5 bis ", delta
                                 hopp(j,i) = hopp(j,i) + 0.0_dp
                              end if
                          end if
                        else ! among others, this is used for the sandwiched systems
                          !print*, "layerindices: ", layerIndex(i), layerIndex(NList(j,i))
                          if (strainedMoire) then
                             dist = sqrt((NeighD(1,j,i)+ umax* sin(2.0_dp*pi*nPeriod*Rat(1,i)/LMoire))**2.0_dp  &
                                         +NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                          else
                             dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                          end if
                          if (KoshinoSR) then
                              aval = 1.34_dp
                              bval = 3.25_dp
                              renormalizeHoppingFactorAAp  = exp((abs(NeighD(3,j,i))-bval)/aval)
                          else
                              renormalizeHoppingFactorAAp  = 1.0_dp
                          end if
                          if (Species(i).eq.Species(NList(j,i)) .and. deactivateIntrasublattice) then
                              hopp(j,i) = 0.0d0
                          !print*, renormalizeHoppingFactorAAp
                          else if (Species(i).ne.Species(NList(j,i)) .and. deactivateIntersublattice) then
                              hopp(j,i) = 0.0d0
                          !print*, renormalizeHoppingFactorAAp
                          else
                              vppsigma = vppsigma0 * exp(-(dist-interlayerdistance)/(BLdelta)) * renormalizeHoppingFactorAAp
                              vpppi = vpppi0 * exp(-(dist-aCC)/(BLdelta)) * renormalizeHoppingFactorAAp
                              ! onlyvppsigma: drop the pi term (Vpppi -> 0); the helper handles the rest.
                              hopp(j,i) = hopp(j,i) + interlayerTwoCenter(merge(0.0_dp,vpppi,onlyvppsigma), vppsigma, &
                                                 NeighD(:,j,i), dist, pzNormal(:,i), pzNormal(:,NList2(j,i)), &
                                                 corrugatedInterlayerTwoCenter)
                              !hopp(j,i) = hopp(j,i) + (vpppi * (1.0_dp - ((NeighD(3,j,i))/(dist))**2.0_dp) &
                              !                  + vppsigma * ((NeighD(3,j,i))/(dist))**2.0_dp)
                          end if
                        end if
                   else if (MIO_StringComp(BilayerModel,'Mayou')) then
                      if (deactivateInterlayer) then
                          cycle
                      end if
                      !print*, "adding the out of plane Mayou terms"
                      if (frac) call AtomsSetCart()
                      if (strainedMoire) then
                         dist = sqrt((NeighD(1,j,i)+ umax* sin(2.0_dp*pi*nPeriod*Rat(1,i)/LMoire))**2.0_dp  &
                                     +NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                      else
                         dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                      end if
                      Fc = 1.0_dp/(1.0_dp+exp((dist-rc)/lc))
                      vppsigma = vppsigma0 * exp(qsigma*(1.0_dp-dist/interlayerdistance)) * Fc * renormalizeHoppingFactorAAp
                      vpppi = vpppi0 * exp(qpi*(1.0_dp - dist/aCC)) * Fc * renormalizeHoppingFactorAAp
                      hopp(j,i) = hopp(j,i) + interlayerTwoCenter(vpppi, vppsigma, NeighD(:,j,i), dist, &
                                        pzNormal(:,i), pzNormal(:,NList2(j,i)), corrugatedInterlayerTwoCenter)
                      !hopp(j,i) = 0.0_dp
                      !hopp(j,i) = - vppsigma * ((NeighD(3,j,i))/(dist))**2.0_dp
                   else if (MIO_StringComp(BilayerModel,'BLKaxiras') .or. MIO_StringComp(BilayerModel,'BLSrivani')) then
                      if (deactivateInterlayer) then
                          cycle
                      end if
                      !print*, "doing BLKaxiras"
                      !theta12 =  acos((Rat(i,1)*Rat(1,NList(j,i)) + Rat(i,2)*Rat(2,NList(j,i))) /
                      !First, theta12
                      !d122 = NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp
                      !index1 = i
                      !index2 = NList(j,i)
                      !index3 = NList(1,index2)
                      !if (index3 .ne. index1) then
                      !  d232 = NeighD(1,1,index2)**2.0_dp+NeighD(2,1,index2)**2.0_dp
                      !else
                      !  d232 = NeighD(1,2,index2)**2.0_dp+NeighD(2,1,index2)**2.0_dp
                      !end if
                      !if (index3 .eq. index1) index3 = NList(2,index2)
                      !rmax = 10.0_dp
                      !do ix=-1,1; do iy=-1,1
                      !   if (index1==index3 .and. ix==0 .and. iy==0) cycle
                      !   ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
                      !   v = Rat(:,index1) - Rat(:,index3) - matmul(ucell,ncell(:,1))
                      !   d = norm(v)
                      !   if (d<rmax) then
                      !      NeighD13_x = -v(1)
                      !      NeighD13_y = -v(2)
                      !   end if
                      !end do
                      !d132 = NeighD13_x**2.0_dp+NeighD13_y**2.0_dp

                      jj = NList(j,i)
                      !if (frac) call AtomsSetCart()
                      if ((.not.((layerIndex(i).eq.2 .and. layerIndex(jj).eq.3) .or. (layerIndex(i).eq.3 .and. layerIndex(jj).eq.2))) .and. fourLayers) then
                        if ((layerIndex(i) .eq. layerIndex(jj) + 1 .or. layerIndex(i).eq. layerIndex(jj)-1) .and. BilayerThreeParameters) then ! INTERLAYER
                          !print*, "doing interlayer for ", i, jj
                          !if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell')) then
                            !if (frac) call AtomsSetCart()
                            dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                            if (abs(NeighD(1,j,i)) < 0.2_dp .and. abs(NeighD(2,j,i)) < 0.2_dp) then
                               hopp(j,i) = hopp(j,i) -tAB1
                               numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               print*, "here7"
                               !hopp(j,i) = interlayerBL(HBL, dx, dy, tAB)
                            else if (dist>aCC*0.95_dp .and. dist<aCC*1.05_dp .and. BilayerThreeParameters) then
                               if (Species(i).eq.Species(NList(j,i))) then
                                  hopp(j,i) = hopp(j,i) -tAB4
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  hopp(j,i) = hopp(j,i) -tAB3
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else if (dist>aCC*0.95_dp*2.0_dp .and. dist<aCC*1.05_dp*2.0_dp .and. BilayerThreeParameters) then
                               if (Species(i).eq.Species(NList(j,i))) then
                                  hopp(j,i) = hopp(j,i) -tAB4
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  hopp(j,i) = hopp(j,i) -tAB3
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else
                               hopp(j,i) = hopp(j,i) + 0.0_dp
                            end if
                        else if ((layerIndex(i) .eq. layerIndex(jj) + 1 .or. layerIndex(i).eq. layerIndex(jj)-1) .and. BilayerOneParameter) then ! INTERLAYER
                          !print*, "doing interlayer for ", i, jj
                          !if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell')) then
                            !if (frac) call AtomsSetCart()
                            dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                            if (abs(NeighD(1,j,i)) < 0.2_dp .and. abs(NeighD(2,j,i)) < 0.2_dp) then
                               hopp(j,i) = hopp(j,i) -tAB1
                               numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               print*, "here8"
                               !hopp(j,i) = interlayerBL(HBL, dx, dy, tAB)
                            else
                               hopp(j,i) = hopp(j,i) + 0.0_dp
                            end if
                        else if ((layerIndex(i) .eq. layerIndex(jj) + 1 .or. layerIndex(i).eq. layerIndex(jj)-1) .and. F2G2Model) then ! INTERLAYER
                            dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                            !if (i.eq.1) then
                            !    print*, "we are applygin the F2G2 model for the interlayer terms in four layer system"
                            !    print*, "aG"
                            !    print*, dist
                            !end if
                            !if (abs(NeighD(1,j,i)) < 0.2_dp .and. abs(NeighD(2,j,i)) < 0.2_dp) then ! On top of each other
                            if (dist<0.2_dp) then
                               !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 0", delta
                               !print*, "2: assigning inerlayer F2G2 terms BA 0", i
                               hopp(j,i) = hopp(j,i) + tBA0
                               numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               print*, "here9"
                               !hopp(j,i) = interlayerBL(HBL, dx, dy, tAB)
                            else if (dist>aCC*0.9_dp .and. dist<aCC*1.1_dp) then ! First nearest neighbors surrounding the atom in the middle of the hexagon
                               if (Species(i).eq.Species(NList(j,i))) then ! same sublattice
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms 1 AA", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms 1 AA", delta
                                  hopp(j,i) = hopp(j,i) + tAA1
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else ! Different sublattice
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms 1 AB", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms 1 AB", delta
                                  hopp(j,i) = hopp(j,i) + tAB1
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else if (dist>aG*0.9_dp .and. dist<aG*1.1_dp) then ! BA' (corresponds to g1, not g2)
                               if (Species(i).eq.Species(NList(j,i))) then
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 2", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms BA 2", delta
                                  hopp(j,i) = hopp(j,i) + tBA2
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 2 bis", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms BA 2 bis", delta
                                  hopp(j,i) = hopp(j,i) + tBA2
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else if (dist>aCC*2.0_dp*0.9_dp .and. dist<aCC*2.0_dp*1.1_dp) then ! Second neirest neighbor from central atom in middle of hexagon
                               if (Species(i).eq.Species(NList(j,i))) then
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms AA 2 ", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms AA 2", delta
                                  hopp(j,i) = hopp(j,i) + tAA2
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms AB 2 ", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms AB 2", delta
                                  hopp(j,i) = hopp(j,i) + tAB2
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else if (dist>(aCC*3.0_dp)*0.9_dp .and. dist<(aCC*3.0_dp)*1.1_dp) then ! BA' (corresponds to g2, not g5)
                               if (Species(i).ne.Species(NList(j,i))) then
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 5 ", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms BA 5", delta
                                  hopp(j,i) = hopp(j,i) + tBA5
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 5 bis ", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms BA 5 bis", delta
                                  hopp(j,i) = hopp(j,i)
                               end if
                            else
                               !if (i.eq.1) print*, "assigning interlayer F2G2 terms BA 5 bis ", delta
                               hopp(j,i) = hopp(j,i) + 0.0_dp
                            end if
                        else
                          if (strainedMoire) then
                             dist = sqrt((NeighD(1,j,i)+ umax* sin(2.0_dp*pi*nPeriod*Rat(1,i)/LMoire))**2.0_dp  &
                                         +NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                          else
                             dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                          end if
                          vppsigma = vppsigma0 * exp(-(dist-interlayerdistance)/(BLdelta))
                 !epsx    x(j,i) = -u0 * 2.0_dp * pi / LMoire * cos(2.0_dp * pi * nPeriod * Rat(1,i) / LMoire)
                          vpppi = vpppi0 * exp(-(dist-aCC)/(BLdelta)) ! inplane
                          hopp(j,i) = hopp(j,i) + interlayerTwoCenter(vpppi, vppsigma, NeighD(:,j,i), dist, &
                                            pzNormal(:,i), pzNormal(:,NList2(j,i)), corrugatedInterlayerTwoCenter)
                        end if
                      else if (((layerIndex(i).eq.1 .and. layerIndex(jj).eq.2) .or. (layerIndex(i).eq.2 .and. layerIndex(jj).eq.1)) .and. threeLayers .and. (.not. middleTwist)) then
                        if ((layerIndex(i) .eq. layerIndex(jj) + 1 .or. layerIndex(i).eq. layerIndex(jj)-1) .and. BilayerThreeParameters) then ! INTERLAYER
                          !print*, "doing interlayer for ", i, jj
                          !if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell')) then
                            !if (frac) call AtomsSetCart()
                            dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                            if (abs(NeighD(1,j,i)) < 0.2_dp .and. abs(NeighD(2,j,i)) < 0.2_dp) then
                               hopp(j,i) = hopp(j,i) -tAB1
                               numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               print*, "here10"
                               !hopp(j,i) = interlayerBL(HBL, dx, dy, tAB)
                            else if (dist>aCC*0.95_dp .and. dist<aCC*1.05_dp .and. BilayerThreeParameters) then
                               if (Species(i).eq.Species(NList(j,i))) then
                                  hopp(j,i) = hopp(j,i) -tAB4
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  hopp(j,i) = hopp(j,i) -tAB3
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else if (dist>aCC*0.95_dp*2.0_dp .and. dist<aCC*1.05_dp*2.0_dp .and. BilayerThreeParameters) then
                               if (Species(i).eq.Species(NList(j,i))) then
                                  hopp(j,i) = hopp(j,i) -tAB4
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  hopp(j,i) = hopp(j,i) -tAB3
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else
                               hopp(j,i) = hopp(j,i) + 0.0_dp
                            end if
                        else if ((layerIndex(i) .eq. layerIndex(jj) + 1 .or. layerIndex(i).eq. layerIndex(jj)-1) .and. F2G2Model) then ! INTERLAYER USING F2G2 model for bottom 2 layers
                            dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                            !if (i.eq.1) then
                            !    print*, "we are applygin the F2G2 model for the interlayer terms in three layer system"
                            !    print*, "aG"
                            !    print*, dist
                            !end if
                            !if (abs(NeighD(1,j,i)) < 0.2_dp .and. abs(NeighD(2,j,i)) < 0.2_dp) then ! On top of each other
                            if (dist<0.2_dp) then
                               !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 0", delta
                               !print*, "2: assigning inerlayer F2G2 terms BA 0", i
                               hopp(j,i) = hopp(j,i) + tBA0
                               numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               print*, "here11"
                               !hopp(j,i) = interlayerBL(HBL, dx, dy, tAB)
                            else if (dist>aCC*0.9_dp .and. dist<aCC*1.1_dp) then ! First nearest neighbors surrounding the atom in the middle of the hexagon
                               if (Species(i).eq.Species(NList(j,i))) then ! same sublattice
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms 1 AA", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms 1 AA", delta
                                  hopp(j,i) = hopp(j,i) + tAA1
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else ! Different sublattice
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms 1 AB", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms 1 AB", delta
                                  hopp(j,i) = hopp(j,i) + tAB1
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else if (dist>aG*0.9_dp .and. dist<aG*1.1_dp) then ! BA' (corresponds to g1, not g2)
                               if (Species(i).eq.Species(NList(j,i))) then
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 2", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms BA 2", delta
                                  hopp(j,i) = hopp(j,i) + tBA2
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 2 bis", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms BA 2 bis", delta
                                  hopp(j,i) = hopp(j,i) + tBA2
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else if (dist>aCC*2.0_dp*0.9_dp .and. dist<aCC*2.0_dp*1.1_dp) then ! Second neirest neighbor from central atom in middle of hexagon
                               if (Species(i).eq.Species(NList(j,i))) then
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms AA 2 ", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms AA 2", delta
                                  hopp(j,i) = hopp(j,i) + tAA2
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms AB 2 ", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms AB 2", delta
                                  hopp(j,i) = hopp(j,i) + tAB2
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               end if
                            else if (dist>(aCC*3.0_dp)*0.9_dp .and. dist<(aCC*3.0_dp)*1.1_dp) then ! BA' (corresponds to g2, not g5)
                               if (Species(i).ne.Species(NList(j,i))) then
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 5 ", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms BA 5", delta
                                  hopp(j,i) = hopp(j,i) + tBA5
                                  numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               else
                                  !if (i.eq.1) print*, "assigning inerlayer F2G2 terms BA 5 bis ", delta
                                  !if (i.eq.2) print*, "2: assigning inerlayer F2G2 terms BA 5 bis", delta
                                  hopp(j,i) = hopp(j,i)
                               end if
                            else
                               !if (i.eq.1) print*, "assigning interlayer F2G2 terms BA 5 bis ", delta
                               hopp(j,i) = hopp(j,i) + 0.0_dp
                            end if
                        else if ((layerIndex(i) .eq. layerIndex(jj) + 1 .or. layerIndex(i).eq. layerIndex(jj)-1) .and. BilayerOneParameter) then ! INTERLAYER
                          !print*, "doing interlayer for ", i, jj
                          !if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell')) then
                            !if (frac) call AtomsSetCart()
                            dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp)
                            if (abs(NeighD(1,j,i)) < 0.2_dp .and. abs(NeighD(2,j,i)) < 0.2_dp) then
                               hopp(j,i) = hopp(j,i) -tAB1
                               numberOfInterlayerHoppings = numberOfInterlayerHoppings + 1
                               print*, "here12"
                               !hopp(j,i) = interlayerBL(HBL, dx, dy, tAB)
                            else
                               hopp(j,i) = hopp(j,i) + 0.0_dp
                            end if
                        else
                          !print*, "using Koshino for outer layers"
                          if (strainedMoire) then
                             dist = sqrt((NeighD(1,j,i)+ umax* sin(2.0_dp*pi*nPeriod*Rat(1,i)/LMoire))**2.0_dp  &
                                         +NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                          else
                             dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                          end if
                          vppsigma = vppsigma0 * exp(-(dist-interlayerdistance)/(BLdelta))
                 !epsx    x(j,i) = -u0 * 2.0_dp * pi / LMoire * cos(2.0_dp * pi * nPeriod * Rat(1,i) / LMoire)
                          vpppi = vpppi0 * exp(-(dist-aCC)/(BLdelta)) ! inplane
                          hopp(j,i) = hopp(j,i) + interlayerTwoCenter(vpppi, vppsigma, NeighD(:,j,i), dist, &
                                            pzNormal(:,i), pzNormal(:,NList2(j,i)), corrugatedInterlayerTwoCenter)
                        end if




                      else if ((layerIndex(i) .eq. layerIndex(jj) +1) .or. (layerIndex(i) .eq. layerIndex(jj) - 1)) then ! adding the kaxiras or srivani hopping terms between layers

                          if (deactivateInterlayerTwisted) then
                              cycle
                          end if

                          d122 = NeighD(1,j,i)**2.0_dp + NeighD(2,j,i)**2.0_dp
                          !if (d122.eq.0.0_dp) then
                          !    print*, "adding the hopping when exactly on top of each other"
                          !    hopp(j,i) = hopp(j,i) - tAB/g0
                          !    cycle
                          !end if
                          !print*, "d122 =", d122
                          if (d122.lt.KaxirasCutoff2 .and. layerIndex(i) .ne. layerIndex(jj)) then ! probably layer condition not important alreayd in interlayer part of code
                              if (d122.eq.0) then
                                 theta12 = 0.0_dp
                                 theta21 = 0.0_dp
                              else if (findThetasGeometrically) then
                                 !if (i.eq.1 .and. j.eq.1) print*, "we are assigning the thetas geometrically"
                                 !print*, "we are assigning the thetas geometrically, check uvec and vvec initialization"
                                 uvec(1) = NeighD(1,j,i)
                                 uvec(2) = NeighD(2,j,i)
                                 uvec(1) = 0.0_dp
                                 vvec(2) = 0.0_dp
                                 vvec(3) = 0.0_dp
                                 if (layerIndex(i).gt.layerIndex(jj)) then ! Attom 1 in top layer (seems like top layer is the unrotated one)
                                     if (Species(i).eq.Species(jj)) then
                                        if (Species(i).eq.1) then
                                             vvec(2) = -aCC ! AAp
                                        else
                                             vvec(2) = aCC ! BBp
                                        end if
                                        vvecR(1) = vvec(1)*cos(twistAngleGrad) - vvec(2)*sin(twistAngleGrad)
                                        vvecR(2) = vvec(1)*sin(twistAngleGrad) + vvec(2)*cos(twistAngleGrad)
                                        vvecR(3) = vvec(3)
                                        theta12 = atan2(norm2(cross(-uvec, vvecR)), dot_product(-uvec, vvecR))
                                        theta21 = atan2(norm2(cross(uvec, vvec)), dot_product(uvec, vvec)) ! Same for AAp and BBp
                                     else
                                        if (Species(i).eq.1) then
                                             vvec(2) = aCC ! ABp
                                        else
                                             vvec(2) = -aCC ! BAp
                                        end if
                                        vvecR(1) = vvec(1)*cos(twistAngleGrad) - vvec(2)*sin(twistAngleGrad)
                                        vvecR(2) = vvec(1)*sin(twistAngleGrad) + vvec(2)*cos(twistAngleGrad)
                                        vvecR(3) = vvec(3)
                                        theta12 = atan2(norm2(cross(-uvec, vvecR)), dot_product(-uvec, vvecR))
                                        theta21 = atan2(norm2(cross(uvec, -vvec)), dot_product(uvec, -vvec))! inverted for ABp and BAp
                                     end if
                                 else if (layerIndex(i).lt.layerIndex(jj)) then ! Atom 1 in bottom layer
                                     if (Species(i).eq.Species(jj)) then
                                        if (Species(i).eq.1) then
                                             vvec(2) = -aCC ! AAp
                                        else
                                             vvec(2) = aCC ! BBp
                                        end if
                                        vvecR(1) = vvec(1)*cos(twistAngleGrad) - vvec(2)*sin(twistAngleGrad)
                                        vvecR(2) = vvec(1)*sin(twistAngleGrad) + vvec(2)*cos(twistAngleGrad)
                                        vvecR(3) = vvec(3)
                                        theta12 = atan2(norm2(cross(-uvec, vvec)), dot_product(-uvec, vvec))
                                        theta21 = atan2(norm2(cross(uvec, vvecR)), dot_product(uvec, vvecR)) ! Same for AAp and BBp
                                     else
                                        if (Species(i).eq.1) then
                                             vvec(2) = aCC ! ABp
                                        else
                                             vvec(2) = -aCC ! BAp
                                        end if
                                        vvecR(1) = vvec(1)*cos(twistAngleGrad) - vvec(2)*sin(twistAngleGrad)
                                        vvecR(2) = vvec(1)*sin(twistAngleGrad) + vvec(2)*cos(twistAngleGrad)
                                        vvecR(3) = vvec(3)
                                        theta12 = atan2(norm2(cross(-uvec, vvec)), dot_product(-uvec, vvec))
                                        theta21 = atan2(norm2(cross(uvec, -vvecR)), dot_product(uvec, -vvecR)) ! inverted for ABp and BAp
                                     end if
                                     theta21 = theta21 !- twistAngleGrad
                                 end if
                                 !print*, "theta12, theta21 with method1: ", theta12/pi*180.0_dp, theta21/pi*180.0_dp, Species(i), Species(jj), layerIndex(i), layerIndex(jj)

                              else
                                 checking = .false.
                                 d232 = 0.0_dp
                                 d132 = 0.0_dp
                                 d232F = 100.0_dp
                                 ! for theta12
                                 do kk=1, Nneigh(jj)
                                     mm = NList(kk,jj)
                                     if(layerIndex(mm).eq.layerIndex(jj)) then
                                        !checking = .true.
                                        d232 = NeighD(1,kk,jj)**2.0_dp + NeighD(2,kk,jj)**2.0_dp
                                        if (d232 .lt. d232F) then
                                           d232F = d232
                                           !print*, "here", d232F, d232
                                           if (d232F .lt. 1.43_dp**2.0_dp .and. d232F .gt. 0.0_dp) then
                                                checking = .true.
                                                exit
                                           end if
                                        end if
                                     end if
                                 end do
                                 !if (checking .eq. .false.) print*, "check1 not good"
                                 checking = .false.

                                 do kkk=1,Nneigh(i)
                                     mmm = NList(kkk,i)
                                     if (mmm.eq.mm) then
                                        checking = .true.
                                        d132 = NeighD(1,kkk,i)**2.0_dp + NeighD(2,kkk,i)**2.0_dp
                                        !print*, "here2", d232F, d232
                                        !if (layerIndex(mmm).eq.layerIndex(i) .or. d132.eq.0.0_dp) print*, "WARNING2", d132, layerIndex(mmm), layerIndex(i)
                                        exit
                                     end if
                                 end do
                                 !if (checking .eq. .false.) print*, "check2 not good"
                                 checking = .false.

                                 if (d232F.lt.0.1_dp .or. d132.le.0.0_dp) print*, "WARNING3", d232F, d132
                                 !if (abs(d122 + d232F - d132) .lt. 0.00001_dp) then
                                 !    theta21 = pi
                                 !else
                                     theta12 = acos(max(-1.0_dp, min(1.0_dp,(d122 + d232F - d132) / (2.0_dp * sqrt(d122) * sqrt(d232F)))))
                                 !end if
                                 !if (isnan(theta12)) print*, "theta12", d122, d232F, d132

                                 ! for theta21
                                 d232 = 0.0_dp
                                 d132 = 0.0_dp
                                 d232F = 100.0_dp
                                 do kkk=1,Nneigh(i)
                                     mmm = NList(kkk,i)
                                     if(layerIndex(mmm).eq.layerIndex(i)) then
                                        !checking = .true.
                                        !print*, "here3", d232F, d232
                                        d132 = NeighD(1,kkk,i)**2.0_dp + NeighD(2,kkk,i)**2.0_dp
                                        if (d132 .lt. 1.43_dp**2.0_dp .and. d132 .gt. 0.0_dp) then
                                            checking = .true.
                                            exit
                                        end if
                                     end if
                                 end do
                                 !if (checking .eq. .false.) print*, "check3 not good"
                                 checking = .false.

                                 do kk=1, Nneigh(jj)
                                     mm = NList(kk,jj)
                                     if(mm.eq.mmm) then
                                        checking = .true.
                                        d232 = NeighD(1,kk,jj)**2.0_dp + NeighD(2,kk,jj)**2.0_dp
                                        !print*, "here4", d232F, d232
                                        !if (layerIndex(mm).eq.layerIndex(jj) .or. d232.eq.0.0_dp) print*, "WARNING5", d232, layerIndex(mm), layerIndex(jj)
                                        !if (d232 .lt. d232F) then
                                        !d232F = d232
                                        !end if
                                        exit
                                     end if
                                 end do
                                 !if (checking .eq. .false.) print*, "check4 not good"
                                 checking = .false.

                                 !if (d232F.lt.0.1_dp .or. d132.le.0.0_dp) print*, "WARNING4"
                                 !if (abs(d122 + d132 - d232) .lt. 0.00001_dp) then
                                 !    theta21 = pi
                                 !else
                                     theta21 = acos(max(-1.0_dp, min(1.0_dp,(d122 + d132 - d232) / (2.0_dp * sqrt(d122) * sqrt(d132)))))
                                 !end if
                                 !if (isnan(theta21)) print*, "theta21", d122, d132, d232
                                 !print*, "theta12, theta21 with method2: ", theta12/pi*180.0_dp, theta21/pi*180.0_dp
                                 !print*, "--------------------------------------------------------------------"

                              end if

                              !if (layerIndex(i).eq.2) then
                              !    dx = NeighD(1,j,i)
                              !    dy = NeighD(2,j,i)
                              !else
                              !    do kk=1, Nneigh(NList(j,i))
                              !        if(NList(kk,NList(j,i)).eq.i) then
                              !           dx = NeighD(1,kk,NList(j,i))
                              !           dy = NeighD(2,kk,NList(j,i))
                              !        end if
                              !    end do
                              !end if
                              !dist = sqrt(dx**2.0_dp+dy**2.0_dp)! +NeighD(3,j,i)**2.0_dp)
                              !if (isnan(theta12) .or. isnan(theta21))  then
                              !       hopp(j,i) = 0.0_dp
                              !       exit
                              !end if
                              dist = sqrt(d122)
                              !print*, "dist =", dist
                              !rbar = dist/2.4389777651302801_dp ! Need to use 2.46 instead? Srivani uses this value that she got from Jeil's code
                              rbar = dist /aGSrivani
                              !print*, "rbar =", rbar
                              !if (MIO_StringComp(BilayerModel,'Srivani')) then
                              !     V0AA = a1AA * exp(-((rbar-b1AA)/c1AA)**4.0_dp) * cos(d1AA*rbar)
                              !     V3AA = a2AA * rbar**2.0_dp *exp(-b2AA*(rbar-c2AA)**2.0_dp)
                              !     V6AA = a3AA * exp(-b3AA * (rbar + c3AA)**2.0_dp)
                              !     V0AB = a1AB * exp(-((rbar-b1AB)/c1AB)**4.0_dp) * cos(d1AB*rbar)
                              !     V3AB = a2AB * rbar**2.0_dp *exp(-b2AB*(rbar-c2AB)**2.0_dp)
                              !     V6AB = a3AB * exp(-b3AB * (rbar + c3AB)**2.0_dp)
                              !     V0 = (V0AA + V0AB)/2.0_dp
                              !     V3 = V3AB/2.0_dp
                              !     V6 = (V6AA + V6AB)/4.0_dp
                              !     hopp(j,i) = hopp(j,i) -(V0 + V3 * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21)))
                              if (MIO_StringComp(BilayerModel,'BLKaxiras') .or. MIO_StringComp(BilayerModel,'BLSrivani')) then
                                   if (useBNGSrivani) then
                                        !if (addPressureDependence) then
                                            epsKax = -(1.0_dp - abs((NeighD(3,j,i)))/(interlayerdistance))
                                            ! AA'
                                            p1a0AA = 7.639_dp
                                            p1b0AA = -6.525_dp
                                            p1c0AA = 4.270_dp
                                            p1d0AA = 0.388_dp
                                            p1h0AA = -2.046_dp
                                            p2a0AA = -2.240_dp
                                            p2b0AA = -1.808_dp
                                            p2c0AA = 1.471_dp
                                            p2d0AA = -0.122_dp
                                            p2h0AA = -0.455_dp
                                            p3a0AA = 0.585_dp
                                            p3b0AA = -1.092_dp
                                            p3c0AA = 1.725_dp
                                            p3d0AA = -2.094_dp
                                            p3h0AA = -0.242_dp
                                            a0AA = p1a0AA * epsKax**2.0_dp + p2a0AA * epsKax + p3a0AA
                                            b0AA = p1b0AA * epsKax**2.0_dp + p2b0AA * epsKax + p3b0AA
                                            c0AA = p1c0AA * epsKax**2.0_dp + p2c0AA * epsKax + p3c0AA
                                            d0AA = p1d0AA * epsKax**2.0_dp + p2d0AA * epsKax + p3d0AA
                                            h0AA = p1h0AA * epsKax**2.0_dp + p2h0AA * epsKax + p3h0AA
                                            a0AA = a0AA/g0
                                            p1a3AA = 0.143_dp
                                            p1b3AA = 130.9_dp
                                            p1c3AA = 9.964_dp
                                            p1d3AA = 4.066_dp
                                            p2a3AA = 0.008_dp
                                            p2b3AA = -5.614_dp
                                            p2c3AA = -0.266_dp
                                            p2d3AA = -0.241_dp
                                            p3a3AA = -0.009_dp
                                            p3b3AA = 3.860_dp
                                            p3c3AA = 1.288_dp
                                            p3d3AA = -3.225_dp
                                            a3AA = p1a3AA * epsKax**2.0_dp + p2a3AA * epsKax + p3a3AA
                                            b3AA = p1b3AA * epsKax**2.0_dp + p2b3AA * epsKax + p3b3AA
                                            c3AA = p1c3AA * epsKax**2.0_dp + p2c3AA * epsKax + p3c3AA
                                            d3AA = p1d3AA * epsKax**2.0_dp + p2d3AA * epsKax + p3d3AA
                                            a3AA = a3AA/g0
                                            p1a6AA = -0.253_dp
                                            p1b6AA = 33.48_dp
                                            p1c6AA = 2.609_dp
                                            p1d6AA = 1.040_dp
                                            p2a6AA = 0.076_dp
                                            p2b6AA = -3.107_dp
                                            p2c6AA = -0.513_dp
                                            p2d6AA = -0.462_dp
                                            p3a6AA = -0.016_dp
                                            p3b6AA = 2.350_dp
                                            p3c6AA = 1.398_dp
                                            p3d6AA = 1.490_dp
                                            a6AA = p1a6AA * epsKax**2.0_dp + p2a6AA * epsKax + p3a6AA
                                            b6AA = p1b6AA * epsKax**2.0_dp + p2b6AA * epsKax + p3b6AA
                                            c6AA = p1c6AA * epsKax**2.0_dp + p2c6AA * epsKax + p3c6AA
                                            d6AA = p1d6AA * epsKax**2.0_dp + p2d6AA * epsKax + p3d6AA
                                            a6AA = a6AA/g0
                                            ! AB'
                                            p1a0AB = 8.098_dp
                                            p1b0AB = -35.33_dp
                                            p1c0AB = 26.12_dp
                                            p1d0AB = 0.859_dp
                                            p1h0AB = 1.328_dp
                                            p2a0AB = -1.135_dp
                                            p2b0AB = -0.989_dp
                                            p2c0AB = 0.712_dp
                                            p2d0AB = 1.083_dp
                                            p2h0AB = -0.281_dp
                                            p3a0AB = 0.289_dp
                                            p3b0AB = -1.000_dp
                                            p3c0AB = 1.664_dp
                                            p3d0AB = -2.271_dp
                                            p3h0AB = -0.235_dp
                                            a0AB = p1a0AB * epsKax**2.0_dp + p2a0AB * epsKax + p3a0AB
                                            b0AB = p1b0AB * epsKax**2.0_dp + p2b0AB * epsKax + p3b0AB
                                            c0AB = p1c0AB * epsKax**2.0_dp + p2c0AB * epsKax + p3c0AB
                                            d0AB = p1d0AB * epsKax**2.0_dp + p2d0AB * epsKax + p3d0AB
                                            h0AB = p1h0AB * epsKax**2.0_dp + p2h0AB * epsKax + p3h0AB
                                            a0AB = a0AB/g0
                                            p1a3AB = 0.821_dp
                                            p1b3AB = 3.425_dp
                                            p1c3AB = 0.961_dp
                                            p1d3AB = -2.317_dp
                                            p2a3AB = -0.370_dp
                                            p2b3AB = 0.161_dp
                                            p2c3AB = 0.404_dp
                                            p2d3AB = 0.102_dp
                                            p3a3AB = 0.081_dp
                                            p3b3AB = 1.908_dp
                                            p3c3AB = 0.633_dp
                                            p3d3AB = -2.006_dp
                                            a3AB = p1a3AB * epsKax**2.0_dp + p2a3AB * epsKax + p3a3AB
                                            b3AB = p1b3AB * epsKax**2.0_dp + p2b3AB * epsKax + p3b3AB
                                            c3AB = p1c3AB * epsKax**2.0_dp + p2c3AB * epsKax + p3c3AB
                                            d3AB = p1d3AB * epsKax**2.0_dp + p2d3AB * epsKax + p3d3AB
                                            a3AB = a3AB/g0
                                            p1a6AB = -0.005_dp
                                            p1b6AB = -22.81_dp
                                            p1c6AB = -2.983_dp
                                            p1d6AB = -0.121_dp
                                            p2a6AB = 0.041_dp
                                            p2b6AB = -3.370_dp
                                            p2c6AB = -0.221_dp
                                            p2d6AB = -0.129_dp
                                            p3a6AB = -0.008_dp
                                            p3b6AB = 2.271_dp
                                            p3c6AB = 1.441_dp
                                            p3d6AB = 1.687_dp
                                            a6AB = p1a6AB * epsKax**2.0_dp + p2a6AB * epsKax + p3a6AB
                                            b6AB = p1b6AB * epsKax**2.0_dp + p2b6AB * epsKax + p3b6AB
                                            c6AB = p1c6AB * epsKax**2.0_dp + p2c6AB * epsKax + p3c6AB
                                            d6AB = p1d6AB * epsKax**2.0_dp + p2d6AB * epsKax + p3d6AB
                                            a6AB = a6AB/g0
                                            ! BA'
                                            p1a0BA =2.856_dp
                                            p1b0BA =1.963_dp
                                            p1c0BA =-1.491_dp
                                            p1d0BA =4.168_dp
                                            p1h0BA =1.479_dp
                                            p2a0BA =-1.758_dp
                                            p2b0BA =-3.274_dp
                                            p2c0BA =2.527_dp
                                            p2d0BA =-0.2887_dp
                                            p2h0BA =-0.8621_dp
                                            p3a0BA =0.541_dp
                                            p3b0BA =-1.270_dp
                                            p3c0BA =1.895_dp
                                            p3d0BA =-1.981_dp
                                            p3h0BA =-0.2946_dp
                                            a0BA = p1a0BA * epsKax**2.0_dp + p2a0BA * epsKax + p3a0BA
                                            b0BA = p1b0BA * epsKax**2.0_dp + p2b0BA * epsKax + p3b0BA
                                            c0BA = p1c0BA * epsKax**2.0_dp + p2c0BA * epsKax + p3c0BA
                                            d0BA = p1d0BA * epsKax**2.0_dp + p2d0BA * epsKax + p3d0BA
                                            h0BA = p1h0BA * epsKax**2.0_dp + p2h0BA * epsKax + p3h0BA
                                            a0BA = a0BA/g0
                                            p1a3BA =-1.870_dp
                                            p1b3BA =-0.442_dp
                                            p1c3BA =-1.202_dp
                                            p1d3BA =-0.529_dp
                                            p2a3BA =0.540_dp
                                            p2b3BA =-1.665_dp
                                            p2c3BA =-0.024_dp
                                            p2d3BA =0.534_dp
                                            p3a3BA =-0.112_dp
                                            p3b3BA =2.074_dp
                                            p3c3BA =0.627_dp
                                            p3d3BA =-1.846_dp
                                            a3BA = p1a3BA * epsKax**2.0_dp + p2a3BA * epsKax + p3a3BA
                                            b3BA = p1b3BA * epsKax**2.0_dp + p2b3BA * epsKax + p3b3BA
                                            c3BA = p1c3BA * epsKax**2.0_dp + p2c3BA * epsKax + p3c3BA
                                            d3BA = p1d3BA * epsKax**2.0_dp + p2d3BA * epsKax + p3d3BA
                                            a3BA = a3BA/g0
                                            p1a6BA = -0.080_dp
                                            p1b6BA = 14.51_dp
                                            p1c6BA = 1.584_dp
                                            p1d6BA = 0.840_dp
                                            p2a6BA = 0.041_dp
                                            p2b6BA = -6.092_dp
                                            p2c6BA = -0.674_dp
                                            p2d6BA = 0.055_dp
                                            p3a6BA = -0.007_dp
                                            p3b6BA = 2.302_dp
                                            p3c6BA = 1.526_dp
                                            p3d6BA = 1.711_dp
                                            a6BA = p1a6BA * epsKax**2.0_dp + p2a6BA * epsKax + p3a6BA
                                            b6BA = p1b6BA * epsKax**2.0_dp + p2b6BA * epsKax + p3b6BA
                                            c6BA = p1c6BA * epsKax**2.0_dp + p2c6BA * epsKax + p3c6BA
                                            d6BA = p1d6BA * epsKax**2.0_dp + p2d6BA * epsKax + p3d6BA
                                            a6BA = a6BA/g0
                                            ! BA'
                                            p1a0BB = 4.504_dp
                                            p1b0BB = -16.91_dp
                                            p1c0BB =12.81_dp
                                            p1d0BB =-2.317_dp
                                            p1h0BB =0.647_dp
                                            p2a0BB =-1.331_dp
                                            p2b0BB =-0.373_dp
                                            p2c0BB =0.291_dp
                                            p2d0BB =1.158_dp
                                            p2h0BB =0.023_dp
                                            p3a0BB =0.325_dp
                                            p3b0BB =-0.795_dp
                                            p3c0BB =1.497_dp
                                            p3d0BB =-2.325_dp
                                            p3h0BB =-0.148_dp
                                            a0BB = p1a0BB * epsKax**2.0_dp + p2a0BB * epsKax + p3a0BB
                                            b0BB = p1b0BB * epsKax**2.0_dp + p2b0BB * epsKax + p3b0BB
                                            c0BB = p1c0BB * epsKax**2.0_dp + p2c0BB * epsKax + p3c0BB
                                            d0BB = p1d0BB * epsKax**2.0_dp + p2d0BB * epsKax + p3d0BB
                                            h0BB = p1h0BB * epsKax**2.0_dp + p2h0BB * epsKax + p3h0BB
                                            a0BB = a0BB/g0
                                            p1a3BB = -0.236_dp
                                            p1b3BB =484.4_dp
                                            p1c3BB =185.6_dp
                                            p1d3BB =24.74_dp
                                            p2a3BB =0.030_dp
                                            p2b3BB =9.329_dp
                                            p2c3BB =10.78_dp
                                            p2d3BB =-5.209_dp
                                            p3a3BB =-0.001_dp
                                            p3b3BB =0.844_dp
                                            p3c3BB =0.818_dp
                                            p3d3BB =-4.035_dp
                                            a3BB = p1a3BB * epsKax**2.0_dp + p2a3BB * epsKax + p3a3BB
                                            b3BB = p1b3BB * epsKax**2.0_dp + p2b3BB * epsKax + p3b3BB
                                            c3BB = p1c3BB * epsKax**2.0_dp + p2c3BB * epsKax + p3c3BB
                                            d3BB = p1d3BB * epsKax**2.0_dp + p2d3BB * epsKax + p3d3BB
                                            a3BB = a3BB/g0
                                            p1a6BB = -0.116_dp
                                            p1b6BB = -21.16_dp
                                            p1c6BB = -1.550_dp
                                            p1d6BB = -1.011_dp
                                            p2a6BB = 0.077_dp
                                            p2b6BB = 0.365_dp
                                            p2c6BB = 0.218_dp
                                            p2d6BB = -0.383_dp
                                            p3a6BB = -0.016_dp
                                            p3b6BB = 2.577_dp
                                            p3c6BB = 1.369_dp
                                            p3d6BB = 1.496_dp
                                            a6BB = p1a6BB * epsKax**2.0_dp + p2a6BB * epsKax + p3a6BB
                                            b6BB = p1b6BB * epsKax**2.0_dp + p2b6BB * epsKax + p3b6BB
                                            c6BB = p1c6BB * epsKax**2.0_dp + p2c6BB * epsKax + p3c6BB
                                            d6BB = p1d6BB * epsKax**2.0_dp + p2d6BB * epsKax + p3d6BB
                                            a6BB = a6BB/g0
                                        !else
                                        !    V0_CB = lambda0_CB * exp(-xi0_CB*rbar**2.0_dp) * cos(kappa0_CB * rbar)
                                        !    V3_BC = lambda3_BC * rbar**2.0_dp * exp(-xi3_BC*(rbar-x3_BC)**2.0_dp)
                                        !    V3_CB = lambda3_CB * rbar**2.0_dp * exp(-xi3_CB*(rbar-x3_CB)**2.0_dp)
                                        !    V0_CN = lambda0_CN * exp(-xi0_CN*rbar**2.0_dp) * cos(kappa0_CN * rbar)
                                        !    V3_NC = lambda3_NC * rbar**2.0_dp * exp(-xi3_NC*(rbar-x3_NC)**2.0_dp)
                                        !    V3_NC = lambda3_CN * rbar**2.0_dp * exp(-xi3_CN*(rbar-x3_CN)**2.0_dp)
                                        !    theta_BC = theta12
                                        !    theta_CB = theta21
                                        !    theta_NC = theta12
                                        !    theta_CN = theta21
                                        !end if
                                            V0AAp = a0AA*exp(-((rbar-b0AA)/c0AA)**4) * cos(d0AA*rbar - h0AA)
                                            V0ABp = a0AB*exp(-((rbar-b0AB)/c0AB)**4) * cos(d0AB*rbar - h0AB)
                                            V0BAp = a0BA*exp(-((rbar-b0BA)/c0BA)**4) * cos(d0BA*rbar - h0BA)
                                            V0BBp = a0BB*exp(-((rbar-b0BB)/c0BB)**4) * cos(d0BB*rbar - h0BB)

                                            if (onlyV0) then
                                                V3AAp = 0.0_dp
                                                V3ABp = 0.0_dp
                                                V3BAp = 0.0_dp
                                                V3BBp = 0.0_dp

                                                V6AAp = 0.0_dp
                                                V6ABp = 0.0_dp
                                                V6BAp = 0.0_dp
                                                V6BBp = 0.0_dp
                                            else
                                                V3AAp = a3AA*rbar*exp(-b3AA*(rbar-c3AA)**2) * sin(d3AA*rbar)
                                                V3ABp = a3AB*rbar*exp(-b3AB*(rbar-c3AB)**2) * sin(d3AB*rbar)
                                                V3BAp = a3BA*rbar*exp(-b3BA*(rbar-c3BA)**2) * sin(d3BA*rbar)
                                                V3BBp = a3BB*rbar*exp(-b3BB*(rbar-c3BB)**2) * sin(d3BB*rbar)

                                                V6AAp = a6AA*rbar*exp(-b6AA*(rbar-c6AA)**2) * sin(d6AA*rbar)
                                                V6ABp = a6AB*rbar*exp(-b6AB*(rbar-c6AB)**2) * sin(d6AB*rbar)
                                                V6BAp = a6BA*rbar*exp(-b6BA*(rbar-c6BA)**2) * sin(d6BA*rbar)
                                                V6BBp = a6BB*rbar*exp(-b6BB*(rbar-c6BB)**2) * sin(d6BB*rbar)
                                            end if
                                            signChange = 1.0_dp
                                            !print*, "adding interlayer terms for GBNtwoLayers"
                                            if (layerIndex(i).eq.1) then ! graphene layer
                                                if (oppositedxdy) then
                                                   theta = atan2(NeighD(1,j,i),NeighD(2,j,i))
                                                else
                                                   theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i))
                                                end if
                                                if (Species(i).eq.1 .and. Species(jj).eq.3) then ! AAp
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = atan2(NeighD(1,j,i),NeighD(2,j,i))
                                                     if (useThetaIJ) then
                                                         hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta12)) + V6AAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorAAp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta)) + V6AAp*(cos(6.0_dp*theta))) * renormalizeHoppingFactorAAp
                                                     end if
                                                else if (Species(i).eq.2 .and. Species(jj).eq.4) then ! BBp
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = atan2(NeighD(1,j,i),NeighD(2,j,i))
                                                     if (useThetaIJ) then
                                                         hopp(j,i) = hopp(j,i) -(V0BBp + V3BBp * (cos(3.0_dp*theta12)) + V6BBp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBBp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0BBp + V3BBp * (cos(3.0_dp*theta)) + V6BBp*(cos(6.0_dp*theta) )) * renormalizeHoppingFactorBBp
                                                     end if
                                                else if (Species(i).eq.1 .and. Species(jj).eq.4) then ! ABp
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = atan2(NeighD(1,j,i),NeighD(2,j,i))
                                                     if (useThetaIJ) then
                                                         hopp(j,i) = hopp(j,i) -(V0ABp + signChange * V3ABp * (cos(3.0_dp*theta12)) + V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorABp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0ABp + signChange * V3ABp * (cos(3.0_dp*theta)) + V6ABp*(cos(6.0_dp*theta) )) * renormalizeHoppingFactorABp
                                                     end if
                                                else if (Species(i).eq.2 .and. Species(jj).eq.3) then ! BAp
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = atan2(NeighD(1,j,i),NeighD(2,j,i))
                                                     if (useThetaIJ) then
                                                         hopp(j,i) = hopp(j,i) -(V0BAp + signChange * V3BAp * (cos(3.0_dp*theta12)) + V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBAp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0BAp + signChange * V3BAp * (cos(3.0_dp*theta)) + V6BAp*(cos(6.0_dp*theta) )) * renormalizeHoppingFactorBAp
                                                     end if
                                                end if
                                            else if (layerIndex(i).eq.2) then ! hBN layer
                                                if (oppositedxdy) then
                                                   theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i))
                                                else
                                                   theta = atan2(NeighD(1,j,i),NeighD(2,j,i)) ! AAp
                                                end if
                                                if (Species(i).eq.3 .and. Species(jj).eq.1) then
                                                     !theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i)) ! AAp
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i)) ! AAp
                                                     !theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i)) ! AAp
                                                     if (useThetaIJ) then
                                                         hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta21)) + V6AAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorAAp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta)) + V6AAp*(cos(6.0_dp*theta) )) * renormalizeHoppingFactorAAp
                                                     end if
                                                else if (Species(i).eq.4 .and. Species(jj).eq.2) then ! BBp
                                                     !theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i)) ! AAp
                                                     !theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i)) ! AAp
                                                     if (useThetaIJ) then
                                                         hopp(j,i) = hopp(j,i) -(V0BBp + V3BBp * (cos(3.0_dp*theta21)) + V6BBp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBBp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0BBp + V3BBp * (cos(3.0_dp*theta)) + V6BBp*(cos(6.0_dp*theta) )) * renormalizeHoppingFactorBBp
                                                     end if
                                                else if (Species(i).eq.4 .and. Species(jj).eq.1) then ! ABp
                                                     !theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i)) ! AAp
                                                     !theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i)) ! AAp
                                                     if (useThetaIJ) then
                                                         hopp(j,i) = hopp(j,i) -(V0ABp + signChange * V3ABp * (cos(3.0_dp*theta21)) + V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorABp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0ABp + signChange * V3ABp * (cos(3.0_dp*theta)) + V6ABp*(cos(6.0_dp*theta) )) * renormalizeHoppingFactorABp
                                                     end if
                                                else if (Species(i).eq.3 .and. Species(jj).eq.2) then ! BAp
                                                     !theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i)) ! AAp
                                                     !theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i)) ! AAp
                                                     if (useThetaIJ) then
                                                         hopp(j,i) = hopp(j,i) -(V0BAp + signChange * V3BAp * (cos(3.0_dp*theta21)) + V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBAp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0BAp + signChange * V3BAp * (cos(3.0_dp*theta)) + V6BAp*(cos(6.0_dp*theta) )) * renormalizeHoppingFactorBAp
                                                     end if
                                                end if
                                            end if
                                        !if (Species(i) .eq. 3 .or. Species(jj) .eq. 3) then
                                        !    hopp(j,i) = hopp(j,i) -(V0_CB + V3_BC * (cos(3.0_dp*theta_BC)) + V3_CB*(cos(3.0_dp*theta_CB)))
                                        !else if (Species(i) .eq. 4 .or. Species(jj) .eq. 4) then
                                        !    hopp(j,i) = hopp(j,i) -(V0_CN + V3_NC * (cos(3.0_dp*theta_NC)) + V3_CN*(cos(3.0_dp*theta_CN)))
                                        !end if
                                   else
                                        epsKax = -(1.0_dp - abs((NeighD(3,j,i)))/(interlayerdistance))
                                        if (addPressureDependence) then
                                            lambda0 = c1_0 + c1_1 * epsKax + c1_2 * epsKax**2.0_dp
                                            lambda3 = c4_0 + c4_1 * epsKax + c4_2 * epsKax**2.0_dp
                                            lambda6 = c7_0 + c7_1 * epsKax + c7_2 * epsKax**2.0_dp
                                            lambda6b = c11_0 + c11_1 * epsKax + c11_2 * epsKax**2.0_dp
                                            lambda0 = lambda0/g0
                                            if (switchV3Sign) then
                                                lambda3 = -lambda3/g0
                                            else
                                                lambda3 = lambda3/g0
                                            end if
                                            lambda6 = lambda6/g0
                                            lambda6b = lambda6b/g0
                                            xi0 = c2_0 + c2_1 * epsKax + c2_2 * epsKax**2.0_dp
                                            xi3 = c5_0 + c5_1 * epsKax + c5_2 * epsKax**2.0_dp
                                            xi6 = c8_0 + c8_1 * epsKax + c8_2 * epsKax**2.0_dp
                                            xi6b = c12_0 + c12_1 * epsKax + c12_2 * epsKax**2.0_dp
                                            x3 = c6_0 + c6_1 * epsKax + c6_2 * epsKax**2.0_dp
                                            x6 = c9_0 + c9_1 * epsKax + c9_2 * epsKax**2.0_dp
                                            x6b = c13_0 + c13_1 * epsKax + c13_2 * epsKax**2.0_dp
                                            kappa0 = c3_0 + c3_1 * epsKax + c3_2 * epsKax**2.0_dp
                                            kappa6 = c10_0 + c10_1 * epsKax + c10_2 * epsKax**2.0_dp
                                            kappa6b = c14_0 + c14_1 * epsKax + c14_2 * epsKax**2.0_dp
                                            !p1a0AA = 11.28_dp
                                            !p1b0AA = -11.33_dp
                                            !p1c0AA = 7.593_dp
                                            !p1d0AA = 2.259_dp
                                            !p1h0AA = 2.074_dp
                                            !p1j0AA = -0.04202_dp
                                            !p1k0AA = -1.03_dp
                                            !p2a0AA = -2.638_dp
                                            !p2b0AA = -0.4476_dp
                                            !p2c0AA = 0.4556_dp
                                            !p2d0AA = -0.4138_dp
                                            !p2h0AA = 0.08359_dp
                                            !p2j0AA = 0.01403_dp
                                            !p2k0AA = -0.1982_dp
                                            !p3a0AA = 0.474_dp
                                            !p3b0AA =  -1.238_dp
                                            !p3c0AA = 1.828_dp
                                            !p3d0AA = 2.348_dp
                                            !p3h0AA = 0.2894_dp
                                            !p3j0AA = -0.003067_dp
                                            !p3k0AA = 2.576_dp
                                            p1a0AA = 9.636_dp
                                            p1b0AA = -7.911_dp
                                            p1c0AA = 5.295_dp
                                            p1d0AA = 1.530_dp
                                            p1h0AA = 1.788_dp
                                            p2a0AA = -2.504_dp
                                            p2b0AA = -0.247_dp
                                            p2c0AA = 0.289_dp
                                            p2d0AA = -0.484_dp
                                            p2h0AA = 0.065_dp
                                            p3a0AA = 0.447_dp
                                            p3b0AA = -1.004_dp
                                            p3c0AA = 1.629_dp
                                            p3d0AA = 2.222_dp
                                            p3h0AA = 0.238
                                            a0AA = p1a0AA * epsKax**2.0_dp + p2a0AA * epsKax + p3a0AA
                                            b0AA = p1b0AA * epsKax**2.0_dp + p2b0AA * epsKax + p3b0AA
                                            c0AA = p1c0AA * epsKax**2.0_dp + p2c0AA * epsKax + p3c0AA
                                            d0AA = p1d0AA * epsKax**2.0_dp + p2d0AA * epsKax + p3d0AA
                                            h0AA = p1h0AA * epsKax**2.0_dp + p2h0AA * epsKax + p3h0AA
                                            j0AA = p1j0AA * epsKax**2.0_dp + p2j0AA * epsKax + p3j0AA
                                            k0AA = p1k0AA * epsKax**2.0_dp + p2k0AA * epsKax + p3k0AA
                                            !p1a0AB = 5.815_dp
                                            !p1b0AB = -2.896_dp
                                            !p1c0AB = 2.021_dp
                                            !p1d0AB = 0.4907_dp
                                            !p1h0AB = 0.2609_dp
                                            !p1j0AB = -0.05488_dp
                                            !p1k0AB = 1.705_dp
                                            !p2a0AB = -2.106_dp
                                            !p2b0AB = -2.119_dp
                                            !p2c0AB = 1.516_dp
                                            !p2d0AB = 0.181_dp
                                            !p2h0AB = -0.4669_dp
                                            !p2j0AB = 0.01493_dp
                                            !p2k0AB = 0.6104_dp
                                            !p3a0AB = 0.4734_dp
                                            !p3b0AB = -1.514_dp
                                            !p3c0AB = 2.041_dp
                                            !p3d0AB = -2.245_dp
                                            !p3h0AB = -0.354_dp
                                            !p3j0AB = -0.002003_dp
                                            !p3k0AB = 2.605_dp
                                            p1a0AB = 5.498_dp
                                            p1b0AB = -1.997_dp
                                            p1c0AB = 1.457_dp
                                            p1d0AB = 0.670_dp
                                            p1h0AB = 0.404_dp
                                            p2a0AB = -2.026_dp
                                            p2b0AB = -1.961_dp
                                            p2c0AB = 1.468_dp
                                            p2d0AB = 0.044_dp
                                            p2h0AB = -0.521_dp
                                            p3a0AB = 0.440_dp
                                            p3b0AB = -1.277_dp
                                            p3c0AB = 1.849_dp
                                            p3d0AB = -2.148_dp
                                            p3h0AB = -0.315_dp
                                            a0AB = p1a0AB * epsKax**2.0_dp + p2a0AB * epsKax + p3a0AB
                                            b0AB = p1b0AB * epsKax**2.0_dp + p2b0AB * epsKax + p3b0AB
                                            c0AB = p1c0AB * epsKax**2.0_dp + p2c0AB * epsKax + p3c0AB
                                            d0AB = p1d0AB * epsKax**2.0_dp + p2d0AB * epsKax + p3d0AB
                                            h0AB = p1h0AB * epsKax**2.0_dp + p2h0AB * epsKax + p3h0AB
                                            j0AB = p1j0AB * epsKax**2.0_dp + p2j0AB * epsKax + p3j0AB
                                            k0AB = p1k0AB * epsKax**2.0_dp + p2k0AB * epsKax + p3k0AB
                                            a0AA = a0AA/g0
                                            j0AA = j0AA/g0
                                            a0AB = a0AB/g0
                                            j0AB = j0AB/g0
                                            c0AAV0 = 0.791_dp
                                            c1AAV0 = -4.663_dp
                                            c2AAV0 = 13.775_dp
                                            c0ABV0 = 0.948_dp
                                            c1ABV0 = -5.57_dp
                                            c2ABV0 = 16.093_dp

                                            aval = 2.48_dp
                                            bval = 3.19_dp
                                            renormalizeHoppingFactorAAp  = exp((abs(NeighD(3,j,i))-bval)/aval)
                                            renormalizeHoppingFactorBBp  = exp((abs(NeighD(3,j,i))-bval)/aval)
                                            aval = 18.17_dp
                                            bval = 2.73_dp
                                            renormalizeHoppingFactorABp  = exp((abs(NeighD(3,j,i))-bval)/aval)
                                            renormalizeHoppingFactorBAp  = exp((abs(NeighD(3,j,i))-bval)/aval)
                                            !renormalizeHoppingFactorAAp  = c2AAV0 * epsKax**2.0_dp + c1AAV0 * epsKax + c0AAV0
                                            !renormalizeHoppingFactorABp  = c2ABV0 * epsKax**2.0_dp + c1ABV0 * epsKax + c0ABV0
                                            !renormalizeHoppingFactorBAp  = c2ABV0 * epsKax**2.0_dp + c1ABV0 * epsKax + c0ABV0
                                            !print*, "renormalizeHoppingFactors", renormalizeHoppingFactorAAp, renormalizeHoppingFactorABp

                                            !print*, "lambda0 =", lambda0
                                            !print*, "lambda3 =", lambda3
                                            !print*, "lambda6 =", lambda6
                                            !print*, "kappa0 =", kappa0
                                            !print*, "kappa3 =", kappa3
                                            !print*, "kappa6 =", kappa6
                                            !print*, "xi0 =", xi0
                                            !print*, "xi3 =", xi3
                                            !print*, "xi6 =", xi6
                                            !print*, "x0 =", x0
                                            !print*, "x3 =", x3
                                            !print*, "x6 =", x6
                                        end if
                                        !if (useBNGSrivani) then
                                        !    V0AAp = a0AA*exp(-((rbar-b0AA)/c0AA)**4) * cos(d0AA*rbar - h0AA)
                                        !    V0ABp = a0AB*exp(-((rbar-b0AB)/c0AB)**4) * cos(d0AB*rbar - h0AB)
                                        !    V0BAp = a0BA*exp(-((rbar-b0BA)/c0BA)**4) * cos(d0BA*rbar - h0BA)
                                        !    V0BBp = a0BB*exp(-((rbar-b0BB)/c0BB)**4) * cos(d0BB*rbar - h0BB)
                                        !    V3AAp = a3AA*rbar*exp(-b3AA*(rbar-c3AA)**2) * sin(d3AA*rbar)
                                        !    V3ABp = a3AB*rbar*exp(-b3AB*(rbar-c3AB)**2) * sin(d3AB*rbar)
                                        !    V3BAp = a3BA*rbar*exp(-b3BA*(rbar-c3BA)**2) * sin(d3BA*rbar)
                                        !    V3BBp = a3BB*rbar*exp(-b3BB*(rbar-c3BB)**2) * sin(d3BB*rbar)
                                        !    V6AAp = a6AA*rbar*exp(-b6AA*(rbar-c6AA)**2) * sin(d6AA*rbar)
                                        !    V6ABp = a6AB*rbar*exp(-b6AB*(rbar-c6AB)**2) * sin(d6AB*rbar)
                                        !    V6BAp = a6BA*rbar*exp(-b6BA*(rbar-c6BA)**2) * sin(d6BA*rbar)
                                        !    V6BBp = a6BB*rbar*exp(-b6BB*(rbar-c6BB)**2) * sin(d6BB*rbar)
                                        !else
                                            V0 = lambda0 * exp(-xi0*rbar**2.0_dp) * cos(kappa0 * rbar)
                                            V3 = lambda3 * rbar**2.0_dp * exp(-xi3*(rbar-x3)**2.0_dp)
                                            if (newFittingFunctions) then
                                               V0AAp = a0AA*exp(-((rbar-b0AA)/c0AA)**4) * cos(d0AA*rbar - h0AA)
                                               V0ABp = a0AB*exp(-((rbar-b0AB)/c0AB)**4) * cos(d0AB*rbar - h0AB)
                                               V0BAp = a0AB*exp(-((rbar-b0AB)/c0AB)**4) * cos(d0AB*rbar - h0AB)
                                            else
                                               V0AAp = lambda0AAp * exp(-xi0AAp*rbar**2.0_dp)*cos(kappa0AAp*rbar) + lambda0bAAp*rbar**3.0_dp*exp(-xi0bAAp*rbar**2.0_dp)*cos(kappa0bAAp*rbar)
                                               V0ABp = lambda0ABp * exp(-xi0ABp*rbar**2.0_dp)*cos(kappa0ABp*rbar) + lambda0bABp*rbar**2.0_dp*exp(-xi0bABp*rbar**2.0_dp)*cos(kappa0bABp*rbar)
                                               V0BAp = lambda0BAp * exp(-xi0BAp*rbar**2.0_dp)*cos(kappa0BAp*rbar) + lambda0bBAp*rbar**2.0_dp*exp(-xi0bBAp*rbar**2.0_dp)*cos(kappa0bBAp*rbar)
                                            end if
                                            if (deactivateV3) then
                                               !print*, "deactivating V3"
                                               V3 = 0.0_dp
                                               V3AAp = 0.0_dp
                                               V3ABp = 0.0_dp
                                               V3BAp = 0.0_dp
                                            else
                                               V3AAp = 0.0_dp
                                               V3ABp = lambda3ABp*rbar**2.0_dp * exp(-xi3ABp*(rbar-x3ABp)**2.0_dp)  + lambda3bABp *rbar**3.0_dp * exp(-xi3bABp*(rbar-x3bABp)**2.0_dp)
                                               V3BAp = lambda3BAp*rbar**2.0_dp * exp(-xi3BAp*(rbar-x3BAp)**2.0_dp)  + lambda3bBAp *rbar**3.0_dp * exp(-xi3bBAp*(rbar-x3bBAp)**2.0_dp)
                                            end if
                                            !print*, V3ABp, V3BAp
                                            V6AAp = 0.0_dp
                                            V6ABp = 0.0_dp
                                            V6BAp = 0.0_dp
                                            !if (DCT) then
                                            !    call cosine_transform_inverse( )
                                            !if (pol) then
                                            !    call V6  =
                                            !else
                                            if (oldParameterSet) then
                                                V6 = lambda6 * exp(-xi6*(rbar-x6)**2.0_dp) * sin(kappa6 * rbar)
                                            else if (deactivateV6) then
                                                V6 = lambda6 * exp(-xi6*(rbar-x6)**2.0_dp) * sin(kappa6 * rbar)
                                            else
                                                V6 = lambda6 * rbar * exp(-xi6*(rbar-x6)**2.0_dp) * sin(kappa6 * rbar) + lambda6b * rbar**2.0_dp*exp(-xi6b*(rbar-x6b)**2.0_dp) * sin(kappa6b * rbar)
                                            end if
                                            if (deactivateV6) then
                                               V6 = 0.0_dp
                                            end if
                                        !end if
                                        !print*, "V0, V3, V6, theta12, theta21", V0, V3, V6, theta12, theta21
                                        !print*, "here5", d232F, d232
                                        if (sublatticeDependent) then
                                          if (fourLayers) then
                                            if (layerIndex(i).eq.2) then
                                                !if (Species(i).eq.Species(jj)) then
                                                !      hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21))) !+ V6AAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) ! The minus sign is added at the end here
                                                !else if (Species(i).eq.1 .and. Species(jj).eq.2) then
                                                !      hopp(j,i) = hopp(j,i) -(V0ABp - V3ABp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21))) !+ V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) ! The minus sign is added at the end here
                                                !else if (Species(i).eq.2 .and. Species(jj).eq.1) then
                                                      hopp(j,i) = hopp(j,i) -(V0BAp + V3BAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21))) !+ V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) ! The minus sign is added at the end here
                                                !end if
                                            else if (layerIndex(i).eq.3) then
                                                !if (Species(i).eq.Species(jj)) then
                                                !      hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21))) !+ V6AAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) ! The minus sign is added at the end here
                                                !else if (Species(i).eq.2 .and. Species(jj).eq.1) then
                                                !      hopp(j,i) = hopp(j,i) -(V0ABp - V3ABp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21))) !+ V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) ! The minus sign is added at the end here
                                                !else if (Species(i).eq.1 .and. Species(jj).eq.2) then
                                                      hopp(j,i) = hopp(j,i) -(V0BAp + V3BAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21))) !+ V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) ! The minus sign is added at the end here
                                                !end if
                                            end if
                                          else if (twoLayers) then
                                            !print*, "holaaaaaaa"
                                            !if (layerIndex(i).eq.1) then
                                            !    if (Species(i).eq.Species(jj)) then
                                            !          !print*, "holaaaaaaa1"
                                            !          hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6AAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorAAp ! The minus sign is added at the end here
                                            !    else if (Species(i).eq.1 .and. Species(jj).eq.2) then
                                            !          !print*, "holaaaaaaa2"
                                            !          hopp(j,i) = hopp(j,i) -(V0ABp - V3ABp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorABp ! The minus sign is added at the end here
                                            !    else if (Species(i).eq.2 .and. Species(jj).eq.1) then
                                            !          !print*, "holaaaaaaa3"
                                            !         hopp(j,i) = hopp(j,i) -(V0BAp + V3BAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBAp ! The minus sign is added at the end here
                                            !    end if
                                            !else if (layerIndex(i).eq.2) then
                                            !    if (Species(i).eq.Species(jj)) then
                                            !          !print*, "holaaaaaaa4"
                                            !          hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6AAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorAAp ! The minus sign is added at the end here
                                            !    else if (Species(i).eq.2 .and. Species(jj).eq.1) then
                                            !         !print*, "holaaaaaaa5"
                                            !         hopp(j,i) = hopp(j,i) -(V0ABp - V3ABp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorABp  ! The minus sign is added at the end here
                                            !    else if (Species(i).eq.1 .and. Species(jj).eq.2) then
                                            !          !print*, "holaaaaaaa6"
                                            !         hopp(j,i) = hopp(j,i) -(V0BAp + V3BAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBAp ! The minus sign is added at the end here
                                            !    end if
                                            !end if
                                            signChange = 1.0_dp
                                            !if (Rat(1,i) + Rat(2,i).lt.1.0_dp) signChange = -1.0_dp
                                            !theta12 = pi-theta12
                                            !theta21 = pi-theta21
                                            if (layerIndex(i).eq.1) then
                                                if (Species(i).eq.Species(jj)) then
                                                     !theta = atan(-NeighD(1,j,i)/-NeighD(2,j,i))
                                                     !theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i))
                                                     !if (NeighD(1,j,i).gt.0.0) theta = theta + 2*pi
                                                     !theta = theta + 2*pi
                                                     !theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = pi/2.0_dp - atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = modulo(theta,2.0_dp*pi)
                                                     !print*, "theta1 =", theta/pi*180.0_dp
                                                     !print*, "1: ", theta/pi*180.0_dp, theta12/pi*180.0_dp, theta21/pi*180.0_dp
                                                     !print*, "1: ", cos(theta), cos(theta12), cos(theta21)
                                                     !theta = atan2(NeighD(1,j,i),NeighD(2,j,i))
                                                     !print*, "theta1 =", theta/pi*180.0_dp
                                                     if (useThetaIJ) then
                                                         !hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6AAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorAAp
                                                         hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta12)) + V6AAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorAAp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta)) + V6AAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorAAp
                                                     end if
                                                else if (Species(i).eq.1 .and. Species(jj).eq.2) then
                                                     !if (neighD(2,j,i) .eq. 0.0_dp .and. neighD(1,j,i) .eq. 0.0_dp) then
                                                     !    theta = 0.0_dp
                                                     !else
                                                         !theta = atan(NeighD(2,j,i)/NeighD(1,j,i))
                                                     !theta = atan(-NeighD(1,j,i)/-NeighD(2,j,i))
                                                     !if (NeighD(2,j,i).lt.0.0) theta = theta + pi
                                                     !theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i))
                                                     !if (NeighD(1,j,i).gt.0.0) theta = theta + 2*pi
                                                     !theta = theta + 2*pi
                                                     !theta = pi/2.0_dp - atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = pi/2.0_dp - atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = modulo(theta,2.0_dp*pi)
                                                     !print*, "theta2 =", theta/pi*180.0_dp
                                                     !print*, "2: ", theta/pi*180.0_dp, theta12/pi*180.0_dp, theta21/pi*180.0_dp
                                                     !print*, "2: ", cos(theta), cos(theta12), cos(theta21)
                                                     !theta = atan2(NeighD(1,j,i),NeighD(2,j,i))
                                                     !end if
                                                     !if (neighD(2,j,i) .lt. 0.0_dp) then
                                                     !    theta = theta + pi
                                                     !end if
                                                     !print*, "theta2 =", theta/pi*180.0_dp
                                                     if (useThetaIJ) then
                                                         !hopp(j,i) = hopp(j,i) -(V0ABp - V3ABp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorABp
                                                         hopp(j,i) = hopp(j,i) -(V0ABp - signChange * V3ABp * (cos(3.0_dp*theta12)) + V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorABp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0ABp - signChange * V3ABp * (cos(3.0_dp*theta)) + V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorABp
                                                         !print*, "ehck"
                                                     end if
                                                else if (Species(i).eq.2 .and. Species(jj).eq.1) then
                                                     !if (neighD(2,j,i) .eq. 0.0_dp .and. neighD(1,j,i) .eq. 0.0_dp) then
                                                     !    theta = 0.0_dp
                                                     !else
                                                         !theta = atan(NeighD(2,j,i)/NeighD(1,j,i))
                                                     !    theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = atan(-NeighD(1,j,i)/-NeighD(2,j,i))
                                                     !if (NeighD(2,j,i).lt.0.0) theta = theta + pi
                                                     !theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i))
                                                     !if (NeighD(1,j,i).gt.0.0) theta = theta + 2*pi
                                                     !theta = theta + 2*pi
                                                     !theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = pi/2.0_dp - atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = modulo(theta,2.0_dp*pi)
                                                     !print*, "theta3 =", theta/pi*180.0_dp
                                                     !print*, "3: ", theta/pi*180.0_dp, theta12/pi*180.0_dp, theta21/pi*180.0_dp
                                                     !print*, "3: ", cos(theta), cos(theta12), cos(theta21)
                                                     !theta = atan2(NeighD(1,j,i),NeighD(2,j,i))
                                                     !end if
                                                     !print*, "theta3 =", theta/pi*180.0_dp
                                                     if (useThetaIJ) then
                                                         !hopp(j,i) = hopp(j,i) -(V0BAp + V3BAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBAp
                                                         hopp(j,i) = hopp(j,i) -(V0BAp + signChange * V3BAp * (cos(3.0_dp*theta12)) + V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBAp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0BAp + signChange * V3BAp * (cos(3.0_dp*theta)) + V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBAp
                                                     end if
                                                end if
                                            else if (layerIndex(i).eq.2) then
                                                if (Species(i).eq.Species(jj)) then
                                                     !theta = atan(NeighD(1,j,i)/NeighD(2,j,i))
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = atan(NeighD(1,j,i)/NeighD(2,j,i))
                                                     !if (NeighD(2,j,i).lt.0.0) theta = theta + pi
                                                     !theta = atan2(NeighD(1,j,i),NeighD(2,j,i))
                                                     !if (NeighD(1,j,i).lt.0.0) theta = theta + 2*pi
                                                     !theta = theta + 2*pi
                                                     !theta = pi/2.0_dp - atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = pi/2.0_dp - atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = modulo(theta,2.0_dp*pi)
                                                     !print*, "theta4 =", theta/pi*180.0_dp
                                                     !print*, "4: ", theta/pi*180.0_dp, theta12/pi*180.0_dp, theta21/pi*180.0_dp
                                                     !theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i))
                                                     !print*, "theta4 =", theta/pi*180.0_dp
                                                     if (useThetaIJ) then
                                                         !hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6AAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorAAp
                                                         hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta21)) + V6AAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorAAp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0AAp + V3AAp * (cos(3.0_dp*theta)) + V6AAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorAAp
                                                     end if
                                                else if (Species(i).eq.2 .and. Species(jj).eq.1) then
                                                     !if (neighD(2,j,i) .eq. 0.0_dp .and. neighD(1,j,i) .eq. 0.0_dp) then
                                                     !    theta = 0.0_dp
                                                     !else
                                                     !    theta = atan(NeighD(1,j,i)/NeighD(2,j,i))
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = atan(NeighD(1,j,i)/NeighD(2,j,i))
                                                     !if (NeighD(2,j,i).lt.0.0) theta = theta + pi
                                                     !theta = atan2(NeighD(1,j,i),NeighD(2,j,i))
                                                     !if (NeighD(1,j,i).lt.0.0) theta = theta + 2*pi
                                                     !theta = theta + 2*pi
                                                     !theta = pi/2.0_dp - atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = pi/2.0_dp - atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = modulo(theta,2.0_dp*pi)
                                                     !print*, "theta5 =", theta/pi*180.0_dp
                                                     !print*, "5: ", theta/pi*180.0_dp, theta12/pi*180.0_dp, theta21/pi*180.0_dp
                                                     !print*, "5: ", cos(theta), cos(theta12), cos(theta21)
                                                     !theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i))
                                                     !end if
                                                     !if (neighD(2,j,i) .lt. 0.0_dp) then
                                                     !    theta = theta + pi
                                                     !end if
                                                     !print*, "theta5 =", theta/pi*180.0_dp
                                                     if (useThetaIJ) then
                                                         !hopp(j,i) = hopp(j,i) -(V0ABp - V3ABp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorABp
                                                         hopp(j,i) = hopp(j,i) -(V0ABp - signChange * V3ABp * (cos(3.0_dp*theta21)) + V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorABp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0ABp - signChange * V3ABp * (cos(3.0_dp*theta)) + V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorABp
                                                     end if
                                                else if (Species(i).eq.1 .and. Species(jj).eq.2) then
                                                     !if (neighD(2,j,i) .eq. 0.0_dp .and. neighD(1,j,i) .eq. 0.0_dp) then
                                                     !    theta = 0.0_dp
                                                     !else
                                                     !    theta = atan(NeighD(1,j,i)/NeighD(2,j,i))
                                                     !theta = atan2(NeighD(1,j,i),NeighD(2,j,i))
                                                     !if (NeighD(1,j,i).lt.0.0) theta = theta + 2*pi
                                                     !theta = theta + 2*pi
                                                     !theta = atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     theta = atan2(-NeighD(2,j,i),-NeighD(1,j,i))
                                                     !theta = pi/2.0_dp - atan2(NeighD(2,j,i),NeighD(1,j,i))
                                                     !theta = modulo(theta,2.0_dp*pi)
                                                     !print*, "theta6 =", theta/pi*180.0_dp
                                                     !print*, "6: ", theta/pi*180.0_dp, theta12/pi*180.0_dp, theta21/pi*180.0_dp
                                                     !print*, "6: ", cos(theta), cos(theta12), cos(theta21)
                                                     !theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i))
                                                     !end if
                                                     !print*, "theta6 =", theta/pi*180.0_dp
                                                     if (useThetaIJ) then
                                                         !hopp(j,i) = hopp(j,i) -(V0BAp + V3BAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBAp
                                                         hopp(j,i) = hopp(j,i) -(V0BAp + signChange * V3BAp * (cos(3.0_dp*theta21)) + V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBAp
                                                     else
                                                         hopp(j,i) = hopp(j,i) -(V0BAp + signChange * V3BAp * (cos(3.0_dp*theta)) + V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBAp
                                                     end if
                                                end if
                                            end if
                                          end if
                                        else
                                          !if (i.eq.1 .and. j.eq.1) print*, "doing sublattice independent"
                                          !hopp(j,i) = hopp(j,i) -(V0 + V3 * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) ! The minus sign is added at the end here
                                          if (layerIndex(i).eq.1) then
                                              theta = atan2(-NeighD(1,j,i),-NeighD(2,j,i))
                                          else if (layerIndex(i).eq.2) then
                                              theta = atan2(NeighD(1,j,i),NeighD(2,j,i))
                                          end if
                                          if (useTheta) then
                                              !print*, "doing sublattice independent using thetaIJ"
                                              !hopp(j,i) = hopp(j,i) -(V0 + V3 * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6ABp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) ! The minus sign is added at the end here
                                              !hopp(j,i) = hopp(j,i) -(V0BAp + V3BAp * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) !* renormalizeHoppingFactorBAp
                                              !hopp(j,i) = hopp(j,i) -(V0 + V3 * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorAAp
                                              hopp(j,i) = hopp(j,i) -(V0 + V3 * (cos(3.0_dp*theta)) + V6*(cos(6.0_dp*theta)))  ! The minus sign is added at the end here
                                          else ! Kaxiras model
                                              !print*, "doing sublattice independent using theta"
                                              !hopp(j,i) = hopp(j,i) -(V0 + V3 * (cos(3.0_dp*theta)) + V6*(cos(6.0_dp*theta)))  ! The minus sign is added at the end here
                                              !hopp(j,i) = hopp(j,i) -(V0BAp + V3BAp * (cos(3.0_dp*theta12)) + V6BAp*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) * renormalizeHoppingFactorBAp
                                              hopp(j,i) = hopp(j,i) -(V0 + V3 * (cos(3.0_dp*theta12) + cos(3.0_dp*theta21)) + V6*(cos(6.0_dp*theta12) + cos(6.0_dp*theta21))) !* renormalizeHoppingFactorAAp
                                          end if
                                        end if
                                        if (addExponentialDecayForDihedral) then
                                            dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                                            hopp(j,i) = hopp(j,i) * exp(-(dist-interlayerdistance)/(BLdelta))
                                        end if
                                        !print*, "hopp1 =", hopp(j,i)
                                        !hopp(j,i) = float(int(hopp(j,i) * 100000000.0_dp + 0.5_dp)) / 100000000.0_dp ! to remove some precision issues leading to non hermitian H
                                        !print*, "hopp2 =", hopp(j,i)
                                   end if
                              end if
                          !else if (d122.eq.0.0_dp) then
                          !    print*, "adding the hopping when exactly on top of each other"
                          !    hopp(j,i) = hopp(j,i) - tAB/g0
                          else
                              hopp(j,i) = hopp(j,i) + 0.0_dp
                          end if
                      else
                          hopp(j,i) = hopp(j,i)
                      end if
                   end if
                   if (renormalizeCoupling) then
                      !print*, "we are renormalizing the couplings"
                      if (differentCouplings) then ! This works for tBG on hBN, maybe need to update for other layered combinations
                         if (t3GwithBN) then
                            if ((layerIndex(i) .eq.1 .and. layerIndex(jj) .eq. 2) .or. (layerIndex(i) .eq.2 .and. layerIndex(jj) .eq. 1)) then
                               hopp(j,i) = couplingFactor2 * hopp(j,i)
                            else
                               hopp(j,i) = couplingFactor * hopp(j,i)
                            end if
                         else if (BNt2GBN) then
                            if ((layerIndex(i) .eq.1 .and. layerIndex(jj) .eq. 2) .or. (layerIndex(i) .eq.2 .and. layerIndex(jj) .eq. 1)) then
                               hopp(j,i) = couplingFactor2 * hopp(j,i)
                            else if ((layerIndex(i) .eq.3 .and. layerIndex(jj) .eq. 4) .or. (layerIndex(i) .eq.4 .and. layerIndex(jj) .eq. 3)) then
                               hopp(j,i) = couplingFactor2 * hopp(j,i)
                            else
                               hopp(j,i) = couplingFactor * hopp(j,i)
                            end if
                         else
                            if ((layerIndex(i) .eq.2 .and. layerIndex(jj) .eq. 3) .or. (layerIndex(i) .eq.3 .and. layerIndex(jj) .eq. 2)) then
                               hopp(j,i) = couplingFactor2 * hopp(j,i)
                            else
                               hopp(j,i) = couplingFactor * hopp(j,i)
                            end if
                         end if
                      else
                         hopp(j,i) = couplingFactor * hopp(j,i)
                      end if
                   end if
                    !else
                    !   hopp(j,i) = 0.0_dp
                    !end if
                  else
                    d = sqrt(dot_product(NeighD(1:2,j,i),NeighD(1:2,j,i)))
                    hopp(j,i) = hopp(j,i) + cmplx(gIntLay*exp(-d/dIntLay),kind=dp)
                  end if
               else ! end of INTERLAYER
                  hopp(j,i) = hopp(j,i) + 0.0_dp
               end if
            end do
        end do
        !!$OMP END PARALLEL DO
        !print*, "counts for GBN"
        !print*, countG1, countG2, countG3, countG4, countG5, countG6, countG7, countG8
        !print*, countBN1, countBN2, countBN3, countBN4, countBN5, countBN6, countBN7, countBN8
        if (MIO_StringComp(BilayerModel,'Jeil')) then ! to add the missing neighbors because of fact*aCC condition
          !$OMP PARALLEL DO PRIVATE(i,j,iii,jj)
          do i=1,nAt
            do j=1,Nneigh(i)
                 if(hopp(j,i) .eq. 0.0_dp) then
                    iii = NList(j,i)
                    do jj=1,Nneigh(iii)
                       if(NList(jj,iii).eq.i) then
                           hopp(j,i) = hopp(j,i) + cmplx(hopp(jj,iii),kind=dp)
                       end if
                    end do
                 end if
            end do
          end do
          !$OMP END PARALLEL DO
        end if
        call MIO_Print(' HAA '//trim(num2str(numberOfHAA1))//' '//trim(num2str(numberOfHAA2)),'HamHopping')
        call MIO_Print(' HAB '//trim(num2str(numberOfHAB1))//' '//trim(num2str(numberOfHAB2)),'HamHopping')
        call MIO_Print(' HBA '//trim(num2str(numberOfHBA1))//' '//trim(num2str(numberOfHBA2)),'HamHopping')

        if (zz) then
          call MIO_InputParameter('LatticepercentFactor',epsFactor,1.0_dp)
! ---  Repeat one more time - This time, add only a term delta... ---
          !$OMP PARALLEL DO PRIVATE(d,nlay,delta,del,realH,imagH,Habjj,dx,dy,HBL), &
          !$OMP& PRIVATE(Cabd), &
          !$OMP& REDUCTION (+:numberOfDel1,numberOfDel2,numberOfDel3,numberOfInterlayerHoppings)
          do i=1,nAt
              delta = 0.0
              nlay = (Species(i)-1)/2 + 1
              if (l) then
                if (ll) then
                   dx = Rat(1,i) * ((1.0_dp+eps2) * cos(twistAngleGrad2) - 1.0_dp) - (Rat(2,i)) * (1.0_dp+eps2) * sin(twistAngleGrad2)
                   dy = (Rat(2,i)) * ((1.0_dp+eps2) * cos(twistAngleGrad2) - 1.0_dp) + Rat(1,i) * (1.0_dp+eps2) * sin(twistAngleGrad2)
                else
                   !dx = Rat(1,i) * ((1.0_dp+eps) * cos(twistAngleGrad2) - 1.0_dp) - (Rat(2,i)) * (1.0_dp+eps) * sin(twistAngleGrad2)
                   !dy = (Rat(2,i)) * ((1.0_dp+eps) * cos(twistAngleGrad2) - 1.0_dp) + Rat(1,i) * (1.0_dp+eps) * sin(twistAngleGrad2)
                   if (rotateFirst) then
                       dx = Rat(1,i) * ((1.0_dp+eps*epsFactor) * cos(twistAngleGrad2) - 1.0_dp) &
                              - (Rat(2,i)) * (1.0_dp+eps*epsFactor) * sin(twistAngleGrad2)
                       dy = (Rat(2,i)) * ((1.0_dp+eps*epsFactor) * cos(twistAngleGrad2) - 1.0_dp) &
                              + Rat(1,i) * (1.0_dp+eps*epsFactor) * sin(twistAngleGrad2)
                   else
                       dx = Rat(1,i)  * ((1.0_dp+eps*epsFactor) * cos(twistAngleGrad2) - 1.0_dp) &
                              - (Rat(2,i) - shift2/eps/epsFactor) * (1.0_dp+eps) * sin(twistAngleGrad2)
                       dy = (Rat(2,i) - shift2/eps/epsFactor) * ((1.0_dp+eps) * cos(twistAngleGrad2) - 1.0_dp) &
                              + Rat(1,i) * (1.0_dp+eps/epsFactor) * sin(twistAngleGrad2)
                   end if
                end if
                if (addDisplacements) then
                    dx = dx + displacements(1,i)
                    dy = dy + displacements(2,i)
                end if
                if (sign2.lt.0.0) then
                    dx = -dx
                    dy = -dy
                end if
                if (rotateFirst) then ! This is the default behavior
                   dx = dx + shift2_x
                   dy = dy + shift2_y
                end if
                ! --- theta is in radian so just multiplied without conversion
                if (distanceDependentEffectiveModel) then
                   call distanceDependentC(Cabd, Cab, BfactorCab, interlayerDistances(i), z0) ! it should be ok to take the first neighbor as it is effective model with only 3 neighbors
                else
                   Cabd = Cab
                end if
                call offdiago(Habjj,dx,dy,Cabd,Phiab)
                if (writeData) then
                   write(587,*) Habjj
                end if
                realH = real(Habjj)
                imagH = imag(Habjj)
                if (sign2.lt.0.0) then
                  imagH = -imagH
                end if
                !print*, "Habjj", realH, imagH
                !if(sign2.gt.0.0) then

                    if (Species(i).eq.2) then   ! ---  B sublattice

                       del(1) = 2.0_dp*realH/3.0_dp
                       del(2) = ( -realH + sqrt(3.0_dp)*imagH )/3.0_dp
                       del(3) = -(  realH + sqrt(3.0_dp)*imagH )/3.0_dp
                     ! del(2) and del(3) are inverted to take into account
                     ! difference
                     ! between A and B
                    end if

                    if (Species(i).eq.1) then   ! ---  A sublattice

                      del(1) =  2.0_dp*realH/3.0_dp
                      del(2) = -(  realH + sqrt(3.0_dp)*imagH )/3.0_dp
                      del(3) =  ( -realH + sqrt(3.0_dp)*imagH )/3.0_dp

                    end if

                !else ! See discussions with Jeil, both GBN and GNB have same
                !virtual strain effect
                !    if (Species(i).eq.1) then   ! ---  A sublattice

                !      del(1) =  2.0_dp*realH/3.0_dp
                !      del(2) = -(  realH + sqrt(3.0_dp)*imagH )/3.0_dp
                !      del(3) =  ( -realH + sqrt(3.0_dp)*imagH )/3.0_dp

                !    end if


                !    if (Species(i).eq.2) then   ! ---  B sublattice

                !       del(1) = 2.0_dp*realH/3.0_dp
                !       del(3) = -(  realH + sqrt(3.0_dp)*imagH )/3.0_dp
                !       del(2) = ( -realH + sqrt(3.0_dp)*imagH )/3.0_dp
                !     ! del(2) and del(3) are inverted to take into account
                !     ! difference
                !     ! between A and B
                !    end if
                !end if
              else
                    del(1) = 0.0_dp
                    del(2) = 0.0_dp
                    del(3) = 0.0_dp
              end if
              do j=1,Nneigh(i)
                 !if (abs(NeighD(3,j,i))<0.01_dp) then
                 if (layerIndex(i) .eq. layerIndex(jj)) then
                    d = sqrt(NeighD(1,j,i)**2+NeighD(2,j,i)**2)
                    do ilvl=1,tbnn
                       if (d < Nradii(ilvl,nlay)) then
                          !hopp(j,i) =
                          !cmplx(gn(Species(i),Species(NList(j,i)),ilvl))
                          if (abs(NeighD(1,j,i)) < 0.01_dp) then
                             delta = del(1)
                             numberOfdel1 = numberOfDel1+1
                          !else if (NeighD(1,j,i) < -0.02_dp ) then
                          !   delta = del(2)
                          !   numberOfdel2 = numberOfDel2+1
                          !else if (NeighD(1,j,i) > 0.02_dp) then
                          !   delta = del(3)
                          !   numberOfdel3 = numberOfDel3+1
                          !endif
                          else if (NeighD(1,j,i) < -0.02_dp ) then
                             if (sign2.lt.0.0) then
                                 delta = del(2)
                             else
                                 delta = del(3)
                             end if
                             numberOfdel2 = numberOfDel2+1
                          else if (NeighD(1,j,i) > 0.02_dp) then
                             if (sign2.lt.0.0) then
                                 delta = del(3)
                             else
                                 delta = del(2)
                             end if
                             numberOfdel3 = numberOfDel3+1
                          endif
                          !print*, delta
                          if (d > 1.7_dp) then ! only use delta !=0 for first NN hoppings
                                               ! useful to separate when other
                                               ! routine asks for more distant
                                               ! neighbors (realStrain for
                                               ! instance)
                              delta = 0.0_dp
                          end if
                          hopp(j,i) = hopp(j,i) - cmplx(1.0_dp*delta ,kind=dp)
                          exit
                       end if
                    end do
                    !if ((Species(i)==1 .and. Species(NList(j,i))==2) .or. &
                    !  (Species(i)==2 .and. Species(NList(j,i))==1)) then
                    !   hopp(j,i) = cmplx(1.0_dp)
                    !else if ((Species(i)==3 .and. Species(NList(j,i))==4) .or. &
                    !  (Species(i)==4 .and. Species(NList(j,i))==3)) then
                    !   hopp(j,i) = cmplx(gBN0)
                    !end if
                 else
                    d = sqrt(dot_product(NeighD(1:2,j,i),NeighD(1:2,j,i)))
                    hopp(j,i) = hopp(j,i) + cmplx(gIntLay*exp(-d/dIntLay),kind=dp)
                 end if
              end do
          end do
        !$OMP END PARALLEL DO
        if (writeData) then
            close(586)
        end if
        end if

! --- Repeatation ended
        call MIO_Print(' All 3 numbers should be equal '//trim(num2str(numberOfDel1))//' '// &
                      trim(num2str(numberOfDel2))//'  '//trim(num2str(numberOfDel3)),'HamHopping')
        call MIO_Print(' Number of interlayer hoppings '//trim(num2str(numberOfInterlayerHoppings)),'HamHopping')
   end if

   call MIO_InputParameter('RandomStrain',randomStrain,.false.)
   if (randomStrain) then
      if (frac) call AtomsSetCart()
      !call MIO_Allocate(fprimex,nAt,'fprimex','ham')
      !call MIO_Allocate(fprimey,nAt,'fprimey','ham')
      !call MIO_Allocate(epsxx,nAt,'epsxx','ham')
      !call MIO_Allocate(epsyy,nAt,'epsyy','ham')
      !call MIO_Allocate(epsxy,nAt,'epsxy','ham')
      call MIO_Allocate(epsxy,[1,inode1],[maxNeigh,inode2],'epsxy','ham')
      call MIO_Allocate(epsxx,[1,inode1],[maxNeigh,inode2],'epsxx','ham')
      call MIO_Allocate(epsyy,[1,inode1],[maxNeigh,inode2],'epsyy','ham')
      call MIO_Allocate(fprimex,[1,inode1],[maxNeigh,inode2],'fprimex','ham')
      call MIO_Allocate(fprimey,[1,inode1],[maxNeigh,inode2],'fprimey','ham')
      !$OMP PARALLEL DO
      do i=1,nAt
         do j=1,Nneigh(i)
           if (NeighD(3,j,i) .lt. 0.01_dp) then
             !if (abs(NeighD(1,j,i)-(aG/2.0_dp)<0.01_dp) .and. HERE ) then
             !    if (NeighD(1,j,i) > 0.0_dp) then
             !        !fprimex(i) = NeighD(3,j,i)/(aG/2.0_dp)
             !        fprimex(i) = (Rat(3,NList(j,i)) - Rat(3,i))/(aG/2.0_dp)
             !    !else
             !    !    fprimex(i) = -(Rat(3,NList(j,i)) - Rat(3,i))/(aG/2.0_dp)
             !    end if
             !end if
             !if (abs(NeighD(2,j,i)-(aG/sqrt(3.0_dp))<0.01_dp)) then
             !    if (NeighD(2,j,i) > 0.0_dp) then
             !        fprimey(i) = (Rat(3,NList(j,i)) - Rat(3,i))/(aG/sqrt(3.0_dp))
             !    else
             !        fprimey(i) = -(Rat(3,NList(j,i)) - Rat(3,i))/(aG/sqrt(3.0_dp))
             !    end if
             !end if
             if (NeighD(1,j,i) > 0.0_dp) then
                 fprimex(j,i) = (Rat(3,NList(j,i)) - Rat(3,i))/(aG/2.0_dp)
             else
                 fprimex(j,i) = -(Rat(3,NList(j,i)) - Rat(3,i))/(aG/2.0_dp)
             end if
             if (NeighD(2,j,i) > 0.0_dp) then
                  fprimey(j,i) = (Rat(3,NList(j,i)) - Rat(3,i))/(aG/sqrt(3.0_dp))
             else
                  fprimey(j,i) = -(Rat(3,NList(j,i)) - Rat(3,i))/(aG/sqrt(3.0_dp))
             end if
             epsxx(j,i) = 0.5_dp * fprimex(j,i)**2.0_dp
             epsyy(j,i) = 0.5_dp * fprimey(j,i)**2.0_dp
             epsxy(j,i) = 0.5_dp * fprimex(j,i)*fprimey(j,i)
           end if
         end do
      end do
      !$OMP END PARALLEL DO
      !print*, 'eps and fprime'
      !print*, epsxx(1), epsyy(1), epsxy(1)
      !print*, epsxx(36), epsyy(36), epsxy(36)
      !print*, fprimex(1), fprimey(1)
      !print*, fprimex(36), fprimey(36)
      !d_lij   = round(1e6 * 1/1.42 * ( epsxx .*xij.^2 + epsyy.*yij.^2 +2*epsxy.*xij.*yij)) * 1e-6;
      !d_tij   = d_tij + t*exp(-3.37*((d_lij+1.42)/1.42-1))-t;
      !$OMP PARALLEL DO PRIVATE(dlij)
      do i=1,nAt
         do j=1,Nneigh(i)
           if (NeighD(3,j,i) .lt. 0.01_dp) then
             !dlij = NINT(10**6 * 1.0_dp/(aG/sqrt(3.0_dp)) * (epsxx(i) * NeighD(1,j,i)**2.0_dp + epsyy(i) * NeighD(2,j,i)**2.0_dp + 2.0_dp * epsxy(i) * NeighD(1,j,i)*NeighD(2,j,i))) * 10**(-6)
             dlij = 1.0_dp/(aG/sqrt(3.0_dp)) * (epsxx(j,i) * NeighD(1,j,i)**2.0_dp &
                  + epsyy(j,i) * NeighD(2,j,i)**2.0_dp + 2.0_dp * epsxy(j,i) * NeighD(1,j,i)*NeighD(2,j,i))
             !print*, dlij
             !if (abs(real(hopp(j,i)))<0.001) then
             !    hopp(j,i) = exp(-3.37_dp*((dlij + aG/sqrt(3.0_dp))/(aG/sqrt(3.0_dp))-1.0_dp))-1.0_dp
             !else
                 hopp(j,i) = hopp(j,i) * exp(-3.37_dp*((dlij + aG/sqrt(3.0_dp))/(aG/sqrt(3.0_dp))-1.0_dp))
             !end if
           end if
         end do
      end do
      !$OMP END PARALLEL DO
   end if

   ! Based on PHYSICAL REVIEW B 80, 045401 2009
   call MIO_InputParameter('realStrain',realStrain,.false.)
   call MIO_InputParameter('onlyFirstNeighborRealStrain',onlyFirstNeighborRealStrain,.false.)
   if (realStrain) then
      call MIO_InputParameter('realStrainReferenceLatticeConstant',aGR,aG) !
      accR = aGR/sqrt(3.0_dp)
      call MIO_InputParameter('realStrainReferenceLatticeConstantBN',aBNR,aBN) !
      aBN1R = aBNR/sqrt(3.0_dp)
! aCCR --> aCC
! aGR --> aG
! aBN1R --> aBN1 = aBN/sqrt(3)
! aBNR --> aBN
      !$OMP PARALLEL DO PRIVATE(dist, refDist)
      do i=1,nAt
         do j=1,Nneigh(i)
            if (layerIndex(i).eq.layerIndex(NList(j,i))) then
               dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
               if (GBNtwoLayers) then
                  if (layerIndex(i).eq.1) then
                     if (dist .gt. acc*0.9_dp .and. dist .lt. acc*1.1_dp) then
                         refDist = aGR/sqrt(3.0_dp)
                     else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                         refDist = aGR
                     else if (dist .gt. 2.0_dp*acc*0.9_dp .and. dist .lt. acc*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*accR
                     else if (dist .gt. 2.0_dp*acc*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (acc*3.0_dp)*0.9_dp .and. dist .lt. (acc*3.0_dp)*1.1_dp) then
                         refDist = accR*3.0_dp
                     end if
                  else if (layerIndex(i).eq.2) then
                     if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                         refDist = aBNR/sqrt(3.0_dp)
                     else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                         refDist = aBNR
                     else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*aBN1R
                     else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                         refDist = aBN1R*3.0_dp
                     end if
                  end if
               else if (encapsulatedThreeLayers) then
                  if (layerIndex(i).eq.2) then
                     if (dist .gt. acc*0.9_dp .and. dist .lt. acc*1.1_dp) then
                         refDist = aGR/sqrt(3.0_dp)
                     else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                         refDist = aGR
                     else if (dist .gt. 2.0_dp*acc*0.9_dp .and. dist .lt. acc*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*accR
                     else if (dist .gt. 2.0_dp*acc*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (acc*3.0_dp)*0.9_dp .and. dist .lt. (acc*3.0_dp)*1.1_dp) then
                         refDist = accR*3.0_dp
                     end if
                  else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.3) then
                     if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                         refDist = aBNR/sqrt(3.0_dp)
                     else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                         refDist = aBNR
                     else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*aBN1R
                     else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                         refDist = aBN1R*3.0_dp
                     end if
                  end if
               else if (encapsulatedFourLayers) then
                  if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3) then
                     if (dist .gt. acc*0.9_dp .and. dist .lt. acc*1.1_dp) then
                         refDist = aGR/sqrt(3.0_dp)
                     else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                         refDist = aGR
                     else if (dist .gt. 2.0_dp*acc*0.9_dp .and. dist .lt. acc*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*accR
                     else if (dist .gt. 2.0_dp*acc*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (acc*3.0_dp)*0.9_dp .and. dist .lt. (acc*3.0_dp)*1.1_dp) then
                         refDist = accR*3.0_dp
                     end if
                  else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                     if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                         refDist = aBNR/sqrt(3.0_dp)
                     else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                         refDist = aBNR
                     else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*aBN1R
                     else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                         refDist = aBN1R*3.0_dp
                     end if
                  end if
               else if (encapsulatedFiveLayers) then
                  if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.4) then
                     if (dist .gt. acc*0.9_dp .and. dist .lt. acc*1.1_dp) then
                         refDist = aGR/sqrt(3.0_dp)
                     else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                         refDist = aGR
                     else if (dist .gt. 2.0_dp*acc*0.9_dp .and. dist .lt. acc*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*accR
                     else if (dist .gt. 2.0_dp*acc*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (acc*3.0_dp)*0.9_dp .and. dist .lt. (acc*3.0_dp)*1.1_dp) then
                         refDist = accR*3.0_dp
                     end if
                  else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.5) then
                     if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                         refDist = aBNR/sqrt(3.0_dp)
                     else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                         refDist = aBNR
                     else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*aBN1R
                     else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                         refDist = aBN1R*3.0_dp
                     end if
                  end if
               else if (encapsulatedSixLayers) then
                  if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.4 .or. layerIndex(i).eq.5) then
                     if (dist .gt. acc*0.9_dp .and. dist .lt. acc*1.1_dp) then
                         refDist = aGR/sqrt(3.0_dp)
                     else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                         refDist = aGR
                     else if (dist .gt. 2.0_dp*acc*0.9_dp .and. dist .lt. acc*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*accR
                     else if (dist .gt. 2.0_dp*acc*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (acc*3.0_dp)*0.9_dp .and. dist .lt. (acc*3.0_dp)*1.1_dp) then
                         refDist = accR*3.0_dp
                     end if
                  else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.6) then
                     if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                         refDist = aBNR/sqrt(3.0_dp)
                     else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                         refDist = aBNR
                     else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*aBN1R
                     else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                         refDist = aBN1R*3.0_dp
                     end if
                  end if
               else if (encapsulatedSevenLayers) then
                  if (layerIndex(i).eq.2 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.4 .or. layerIndex(i).eq.5 .or. layerIndex(i).eq.6) then
                     if (dist .gt. acc*0.9_dp .and. dist .lt. acc*1.1_dp) then
                         refDist = aGR/sqrt(3.0_dp)
                     else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                         refDist = aGR
                     else if (dist .gt. 2.0_dp*acc*0.9_dp .and. dist .lt. acc*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*accR
                     else if (dist .gt. 2.0_dp*acc*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (acc*3.0_dp)*0.9_dp .and. dist .lt. (acc*3.0_dp)*1.1_dp) then
                         refDist = accR*3.0_dp
                     end if
                  else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.7) then
                     if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                         refDist = aBNR/sqrt(3.0_dp)
                     else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                         refDist = aBNR
                     else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*aBN1R
                     else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                         refDist = aBN1R*3.0_dp
                     end if
                  end if
               else if (t3GwithBN) then
                  if (layerIndex(i).eq.4 .or. layerIndex(i).eq.3 .or. layerIndex(i).eq.2) then
                     if (dist .gt. acc*0.9_dp .and. dist .lt. acc*1.1_dp) then
                         refDist = aGR/sqrt(3.0_dp)
                     else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                         refDist = aGR
                     else if (dist .gt. 2.0_dp*acc*0.9_dp .and. dist .lt. acc*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*accR
                     else if (dist .gt. 2.0_dp*acc*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (acc*3.0_dp)*0.9_dp .and. dist .lt. (acc*3.0_dp)*1.1_dp) then
                         refDist = accR*3.0_dp
                     end if
                  else if (layerIndex(i).eq.1) then
                     if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                         refDist = aBNR/sqrt(3.0_dp)
                     else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                         refDist = aBNR
                     else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*aBN1R
                     else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                         refDist = aBN1R*3.0_dp
                     end if
                  end if
               else if (BNt2GBN) then
                  if (layerIndex(i).eq.3 .or. layerIndex(i).eq.2) then
                     if (dist .gt. acc*0.9_dp .and. dist .lt. acc*1.1_dp) then
                         refDist = aGR/sqrt(3.0_dp)
                     else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                         refDist = aGR
                     else if (dist .gt. 2.0_dp*acc*0.9_dp .and. dist .lt. acc*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*accR
                     else if (dist .gt. 2.0_dp*acc*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (acc*3.0_dp)*0.9_dp .and. dist .lt. (acc*3.0_dp)*1.1_dp) then
                         refDist = accR*3.0_dp
                     end if
                  else if (layerIndex(i).eq.1 .or. layerIndex(i).eq.4) then
                     if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                         refDist = aBNR/sqrt(3.0_dp)
                     else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                         refDist = aBNR
                     else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*aBN1R
                     else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                         refDist = aBN1R*3.0_dp
                     end if
                  end if
               else if (t2GBN) then
                  if (layerIndex(i).eq.3 .or. layerIndex(i).eq.2) then
                     if (dist .gt. acc*0.9_dp .and. dist .lt. acc*1.1_dp) then
                         refDist = aGR/sqrt(3.0_dp)
                     else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                         refDist = aGR
                     else if (dist .gt. 2.0_dp*acc*0.9_dp .and. dist .lt. acc*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*accR
                     else if (dist .gt. 2.0_dp*acc*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (acc*3.0_dp)*0.9_dp .and. dist .lt. (acc*3.0_dp)*1.1_dp) then
                         refDist = accR*3.0_dp
                     end if
                  else if (layerIndex(i).eq.1) then
                     if (dist .gt. aBN1*0.9_dp .and. dist .lt. aBN1*1.1_dp) then
                         refDist = aBNR/sqrt(3.0_dp)
                     else if (dist .gt. aBN*0.9_dp .and. dist .lt. aBN*1.1_dp) then
                         refDist = aBNR
                     else if (dist .gt. 2.0_dp*aBN1*0.9_dp .and. dist .lt. aBN1*2.0_dp*1.1_dp) then
                         refDist = 2.0_dp*aBN1R
                     else if (dist .gt. 2.0_dp*aBN1*1.1_dp .and. dist .lt. 4.0_dp) then
                         refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                     else if (dist .gt. (aBN1*3.0_dp)*0.9_dp .and. dist .lt. (aBN1*3.0_dp)*1.1_dp) then
                         refDist = aBN1R*3.0_dp
                     end if
                  end if
               else
                  if (dist .gt. acc*0.9_dp .and. dist .lt. acc*1.1_dp) then
                      refDist = aGR/sqrt(3.0_dp)
                  else if (dist .gt. aG*0.9_dp .and. dist .lt. aG*1.1_dp) then
                      refDist = aGR
                  else if (dist .gt. 2.0_dp*acc*0.9_dp .and. dist .lt. acc*2.0_dp*1.1_dp) then
                      refDist = 2.0_dp*accR
                  else if (dist .gt. 2.0_dp*acc*1.1_dp .and. dist .lt. 4.0_dp) then
                      refDist = 3.75696_dp ! define this value based on aG and aCC, does it matter of 4th nearest neighbor equal 0 in F2G2 model?
                  else if (dist .gt. (acc*3.0_dp)*0.9_dp .and. dist .lt. (acc*3.0_dp)*1.1_dp) then
                      refDist = accR*3.0_dp
                  end if
               end if
               if (onlyFirstNeighborRealStrain) then
                   if (dist < 1.5_dp) then
                       !print*,  "doing NN1"
                       hopp(j,i) = hopp(j,i) * exp(-3.37_dp*(dist / (refDist)-1.0_dp))
                   end if
               else
                   !print*, "not doing NN1"
                   hopp(j,i) = hopp(j,i) * exp(-3.37_dp*(dist / (refDist)-1.0_dp))
               end if
            end if
         end do
      end do
      !$OMP END PARALLEL DO
   end if
   !print*, hopp

   call MIO_InputParameter('periodicStrain',periodicStrain,.false.)
   if (periodicStrain) then
      call MIO_InputParameter('periodicStrainu0',u0,0.1_dp)
      call MIO_InputParameter('periodicStrainPeriod',nPeriod,1)
      call MIO_InputParameter('SuperCell',sCell,60)
      LMoire = norm(ucell(:,1))/sCell
      !print*, "LMoire= ", LMoire
      if (frac) call AtomsSetCart()
      call MIO_Allocate(epsxy,[1,inode1],[maxNeigh,inode2],'epsxy','ham')
      call MIO_Allocate(epsxx,[1,inode1],[maxNeigh,inode2],'epsxx','ham')
      call MIO_Allocate(epsyy,[1,inode1],[maxNeigh,inode2],'epsyy','ham')
      call MIO_Allocate(fprimex,[1,inode1],[maxNeigh,inode2],'fprimex','ham')
      call MIO_Allocate(fprimey,[1,inode1],[maxNeigh,inode2],'fprimey','ham')
      !$OMP PARALLEL DO
      do i=1,nAt
         do j=1,Nneigh(i)
             if (NeighD(1,j,i) > 0.0_dp) then
                 epsxx(j,i) = -u0 * 2.0_dp * pi / LMoire * cos(2.0_dp * pi * nPeriod * Rat(1,i) / LMoire)
             else
                 epsxx(j,i) = -u0 * 2.0_dp * pi / LMoire * cos(2.0_dp * pi * nPeriod * Rat(1,NList(j,i)) / LMoire)
             end if
             if (NeighD(2,j,i) > 0.0_dp) then
                 epsyy(j,i) = -u0 * 2.0_dp * pi / LMoire * cos(2.0_dp * pi * nPeriod * Rat(2,i) / LMoire)
             else
                 epsyy(j,i) = -u0 * 2.0_dp * pi / LMoire * cos(2.0_dp * pi * nPeriod * Rat(2,NList(j,i)) / LMoire)
             end if
             epsxy(j,i) = 0.0_dp
             !epsxx(j,i) = 0.5_dp * fprimex(j,i)**2.0_dp
             !epsyy(j,i) = 0.5_dp * fprimey(j,i)**2.0_dp
             !epsxy(j,i) = 0.5_dp * fprimex(j,i)*fprimey(j,i)
         end do
      end do
      !$OMP END PARALLEL DO
      !print*, 'eps and fprime'
      !print*, epsxx(1), epsyy(1), epsxy(1)
      !print*, epsxx(36), epsyy(36), epsxy(36)
      !print*, fprimex(1), fprimey(1)
      !print*, fprimex(36), fprimey(36)
      !d_lij   = round(1e6 * 1/1.42 * ( epsxx .*xij.^2 + epsyy.*yij.^2 +2*epsxy.*xij.*yij)) * 1e-6;
      !d_tij   = d_tij + t*exp(-3.37*((d_lij+1.42)/1.42-1))-t;
      !$OMP PARALLEL DO PRIVATE(dlij)
      do i=1,nAt
         do j=1,Nneigh(i)
           if (NeighD(3,j,i) .lt. 0.01_dp) then
             !dlij = NINT(10**6 * 1.0_dp/(aG/sqrt(3.0_dp)) * (epsxx(i) * NeighD(1,j,i)**2.0_dp + epsyy(i) * NeighD(2,j,i)**2.0_dp + 2.0_dp * epsxy(i) * NeighD(1,j,i)*NeighD(2,j,i))) * 10**(-6)
             dlij = 1.0_dp/(aG/sqrt(3.0_dp)) * (epsxx(j,i) * NeighD(1,j,i)**2.0_dp &
                  + epsyy(j,i) * NeighD(2,j,i)**2.0_dp + 2.0_dp * epsxy(j,i) * NeighD(1,j,i)*NeighD(2,j,i))
             !print*, dlij
             !if (abs(real(hopp(j,i)))<0.001) then
             !    hopp(j,i) = exp(-3.37_dp*((dlij + aG/sqrt(3.0_dp))/(aG/sqrt(3.0_dp))-1.0_dp))-1.0_dp
             !else
                 hopp(j,i) = hopp(j,i) * exp(-3.37_dp*((dlij + aG/sqrt(3.0_dp))/(aG/sqrt(3.0_dp))-1.0_dp))
             !end if
           end if
         end do
      end do
      !$OMP END PARALLEL DO
   end if

   call MIO_InputParameter('realisticBubbles',realisticBubbles,.false.)
   if (realisticBubbles) then
      call MIO_InputParameter('BigBubble',bigBubble,.false.)
      call MIO_InputParameter('manyBubbles',manyBubbles,.false.)
      call MIO_InputParameter('bubbleInPlaneStrain',bubbleInPlaneStrain,.false.)
      if (bigBubble) then
         call MIO_Print('We add a single bubble','ham')
         if (.not. frac) call AtomsSetFrac()
         minDiffBubbleX = 0.01_dp
         minDiffBubbleY = 0.01_dp
         totImp = 1
         call MIO_Allocate(indxImp,totImp,'indxImp','gauss')
         call MIO_Allocate(bubbleCenterX,totImp,'bubbleCenterX','gauss')
         call MIO_Allocate(bubbleCenterY,totImp,'bubbleCenterY','gauss')
         !!$OMP PARALLEL DO # not sure parallelization would work here
         !do i=1,nAt
         do i=inode1,inode2
            diffBubbleX = abs(Rat(1,i) - 0.5d0)
            diffBubbleY = abs(Rat(2,i) - 0.5d0)
            if (diffBubbleX.lt.minDiffBubbleX .and. diffBubbleY.lt.minDiffBubbleY) then
                minDiffBubbleX = diffBubbleX
                minDiffBubbleY = diffBubbleY
                indxImp(1) = i
            end if
         end do
         !!$OMP END PARALLEL DO
         !print*, "Center of bubble at (in frac): ", Rat(1,indxImp(1)), Rat(2,indxImp(1))
         if (frac) call AtomsSetCart()
         bubbleCenterX(1) = Rat(1,indxImp(1))
         bubbleCenterY(1) = Rat(2,indxImp(1))
         !print*, "indxImp = ", indxImp
         indexOfBigBubbleCenter = indxImp(1)
      else if (manyBubbles) then
         call MIO_Print('We add a lot of bubbles','ham')
         call MIO_InputParameter('bubblePercentage',per,10.0_dp)
         per = per/100.0_dp
         call MIO_Allocate(def,[inode1],[inode2],'def','gauss')
         def = .false.
         totImp = 0

         call random_seed(size = kk)
         allocate(seed(kk))
         call system_clock(COUNT=clock)
         seed = clock + 37 * (/ (i - 1, i = 1, kk) /)
         call random_seed(PUT = seed)
         deallocate(seed)

         nImp = nint((inode2-inode1+1)*per)
         do i=1,nImp
            do
               !j = nint(RandNum(rng)*(in2-in1)) + in1
               call random_number(rand)
               j = nint(rand*(inode2-inode1)) + inode1
               !if (j>in2 .or. j<in1) write(*,*) 'Gauss err, j:', j
               if (def(j) .eqv. .false.) then
                  def(j) = .true.
                  totImp = totImp + 1
                  exit
               end if
            end do
         end do

         per = totImp*100.0_dp/nAt
         call MIO_Print('')
         call MIO_Print('  final percentage          : '//trim(num2str(per,3))//'%')
         call MIO_Print('  total number of impurities: '//trim(num2str(totImp)))
         call MIO_Allocate(indxImp,totImp,'indxImp','gauss')
         call MIO_Allocate(bubbleCenterX,totImp,'bubbleCenterX','gauss')
         call MIO_Allocate(bubbleCenterY,totImp,'bubbleCenterY','gauss')
         nImp = 0
         if (frac) call AtomsSetCart()
         do i=inode1,inode2
            if (def(i)) then
               nImp = nImp + 1
               indxImp(nImp) = i
               bubbleCenterX(nImp) = Rat(1,indxImp(nImp))
               bubbleCenterY(nImp) = Rat(2,indxImp(nImp))
            end if
         end do
         call MIO_InputParameter('bubbleGaussian',bubbleGauss,.false.)
         if (bubbleGauss) then
            call GaussPotDefinedPositions(H0)
         end if
         !call MIO_InputParameter('WriteDataFiles',l,.false.)
         !if (l) then
         !   open(1,FILE='e')
         !   do i=1,nAt
         !      write(1,*) H0(i)
         !   end do
         !   close(1)
         !end if
      end if

      !call MIO_Allocate(fprimex,nAt,'fprimex','ham')
      !call MIO_Allocate(fprimey,nAt,'fprimey','ham')
      !call MIO_Allocate(epsxx,nAt,'epsxx','ham')
      !call MIO_Allocate(epsyy,nAt,'epsyy','ham')
      !call MIO_Allocate(epsxy,nAt,'epsxy','ham')
      call MIO_Allocate(epsxy,[1,inode1],[maxNeigh,inode2],'epsxy','ham')
      call MIO_Allocate(epsxx,[1,inode1],[maxNeigh,inode2],'epsxx','ham')
      call MIO_Allocate(epsyy,[1,inode1],[maxNeigh,inode2],'epsyy','ham')
      call MIO_InputParameter('bubbleSigmaR',bubbleSigmaR,1.42_dp)
      bubbleSigmaR2 = bubbleSigmaR**2.0d0
      !print*, "bubbleSigmaR2 = ", bubbleSigmaR2
      call MIO_InputParameter('bubbleC',bubbleC,1.42_dp)
      !print*, "bubbleC = ", bubbleC
      !call MIO_Allocate(fprimex,[1,inode1],[maxNeigh,inode2],'fprimex','ham')
      !call MIO_Allocate(fprimey,[1,inode1],[maxNeigh,inode2],'fprimey','ham')

      !print*, "indxImp = ", indxImp
      !indexOfBigBubbleCenter = indxImp(1)
      call MIO_InputParameter('bubbleShift',bubbleShift,0.07_dp)
      call MIO_InputParameter('bubbleRadius',bubbleRadius,100.0_dp)
      call MIO_Allocate(onsiteShift,nAt,'onsiteShift','ham')

      do ii=1,totImp
        !$OMP PARALLEL DO PRIVATE(dx,dy,dist2,bubbleR,bubbleTheta,dlij)
        do i=1,nAt
          !do ic1=-1,1; do ic2=-1,1
           dx = (Rat(1,indxImp(ii)) - Rat(1,i))! + ic1*ucell(1,1) + ic2*ucell(1,2)
           dy = (Rat(2,indxImp(ii)) - Rat(2,i))! + ic1*ucell(2,1) + ic2*ucell(2,2)
           dist2 = dx**2.0d0 + dy**2.0d0 !(the correlation function is in x only
           if (dist2.lt.bubbleSigmaR2) then
             onsiteShift(i) = bubbleShift
           else
             onsiteShift(i) = 0.0d0
           end if
           H0(i) = H0(i) + onsiteShift(i)
           do j=1,Nneigh(i)
           !if (NeighD(3,j,i) .lt. 0.01_dp) then
             !!if (abs(NeighD(1,j,i)-(aG/2.0_dp)<0.01_dp) .and. HERE ) then
             !!    if (NeighD(1,j,i) > 0.0_dp) then
             !!        !fprimex(i) = NeighD(3,j,i)/(aG/2.0_dp)
             !!        fprimex(i) = (Rat(3,NList(j,i)) - Rat(3,i))/(aG/2.0_dp)
             !!    !else
             !!    !    fprimex(i) = -(Rat(3,NList(j,i)) - Rat(3,i))/(aG/2.0_dp)
             !!    end if
             !!end if
             !!if (abs(NeighD(2,j,i)-(aG/sqrt(3.0_dp))<0.01_dp)) then
             !!    if (NeighD(2,j,i) > 0.0_dp) then
             !!        fprimey(i) = (Rat(3,NList(j,i)) - Rat(3,i))/(aG/sqrt(3.0_dp))
             !!    else
             !!        fprimey(i) = -(Rat(3,NList(j,i)) - Rat(3,i))/(aG/sqrt(3.0_dp))
             !!    end if
             !!end if
             !if (NeighD(1,j,i) > 0.0_dp) then
             !    fprimex(j,i) = (Rat(3,NList(j,i)) - Rat(3,i))/(aG/2.0_dp)
             !else
             !    fprimex(j,i) = -(Rat(3,NList(j,i)) - Rat(3,i))/(aG/2.0_dp)
             !end if
             !if (NeighD(2,j,i) > 0.0_dp) then
             !     fprimey(j,i) = (Rat(3,NList(j,i)) - Rat(3,i))/(aG/sqrt(3.0_dp))
             !else
             !     fprimey(j,i) = -(Rat(3,NList(j,i)) - Rat(3,i))/(aG/sqrt(3.0_dp))
             !end if
             !epsxx(j,i) = 0.5_dp * fprimex(j,i)**2.0_dp
             !epsyy(j,i) = 0.5_dp * fprimey(j,i)**2.0_dp
             !epsxy(j,i) = 0.5_dp * fprimex(j,i)*fprimey(j,i)
             bubbleR = sqrt((Rat(1,i)-bubbleCenterX(ii))**2.0d0 + (Rat(2,i)-bubbleCenterY(ii))**2.0d0)
             bubbleTheta = atan2(Rat(2,i)-bubbleCenterY(ii),Rat(1,i)-bubbleCenterX(ii))
             if (bubbleInPlaneStrain) then
                !epsxx   =  2 *C *r.* sin(theta) .*exp(-(r>R).*(r-R).^2/(2*sigmaR^2));
                !epsyy   = -2 *C *r.* sin(theta) .*exp(-(r>R).*(r-R).^2/(2*sigmaR^2));
                !epsxy   =  2 *C *r.* cos(theta) .*exp(-(r>R).*(r-R).^2/(2*sigmaR^2));
                if (bubbleR.gt.bubbleRadius) then
                  epsxx(j,i) = bubbleC*2.0d0 * bubbleR *  sin(bubbleTheta) &
                      * exp(-(bubbleR-bubbleRadius)**2.0d0/(2.0d0*bubbleSigmaR**2.0d0))
                  epsyy(j,i) = -bubbleC*2.0d0 * bubbleR * sin(bubbleTheta) &
                      * exp(-(bubbleR-bubbleRadius)**2.0d0/(2.0d0*bubbleSigmaR**2.0d0))
                  epsxy(j,i) = bubbleC*2.0d0 * bubbleR *  cos(bubbleTheta) &
                      * exp(-(bubbleR-bubbleRadius)**2.0d0/(2.0d0*bubbleSigmaR**2.0d0))
                else
                  epsxx(j,i) = bubbleC*2.0d0 * bubbleR *  sin(bubbleTheta)
                  epsyy(j,i) = -bubbleC*2.0d0 * bubbleR * sin(bubbleTheta)
                  epsxy(j,i) = bubbleC*2.0d0 * bubbleR *  cos(bubbleTheta)
                  !print*, epsxx(j,i), epsyy(j,i), epsxy(j,i)
                end if
             else
                epsxx(j,i) = bubbleC**2.0d0 * bubbleR**2.0d0 &
                  * exp(-2.0d0*bubbleR**2.0d0/(2.0d0*bubbleSigmaR**2.0d0)) * cos(bubbleTheta)**2.0d0/(2.0d0*bubbleSigmaR**4.0d0)
                epsyy(j,i) = bubbleC**2.0d0 * bubbleR**2.0d0 &
                  * exp(-2.0d0*bubbleR**2.0d0/(2.0d0*bubbleSigmaR**2.0d0)) * sin(bubbleTheta)**2.0d0/(2.0d0*bubbleSigmaR**4.0d0)
                epsxy(j,i) = bubbleC**2.0d0 * bubbleR**2.0d0 &
                  * exp(-2.0d0*bubbleR**2.0d0/(2.0d0*bubbleSigmaR**2.0d0)) &
                                   * cos(bubbleTheta)*sin(bubbleTheta)/(2.0d0*bubbleSigmaR**4.0d0)
             end if
             !print*, bubbleR, bubbleTheta, bubbleCenterX, bubbleCenterY, epsxx(j,i), epsyy(j,i), epsxy(j,i)
             !if (NeighD(3,j,i) .lt. 0.01_dp) then
                !dlij = NINT(10**6 * 1.0_dp/(aG/sqrt(3.0_dp)) * (epsxx(i) * NeighD(1,j,i)**2.0_dp + epsyy(i) * NeighD(2,j,i)**2.0_dp + 2.0_dp * epsxy(i) * NeighD(1,j,i)*NeighD(2,j,i))) * 10**(-6)
                dlij = 1.0_dp/(aG/sqrt(3.0_dp)) * (epsxx(j,i) * NeighD(1,j,i)**2.0_dp &
                       + epsyy(j,i) * NeighD(2,j,i)**2.0_dp + 2.0_dp * epsxy(j,i) * NeighD(1,j,i)*NeighD(2,j,i))
                !print*, dlij
                !if (abs(real(hopp(j,i)))<0.001) then
                !    hopp(j,i) = exp(-3.37_dp*((dlij + aG/sqrt(3.0_dp))/(aG/sqrt(3.0_dp))-1.0_dp))-1.0_dp
                !else
                 hopp(j,i) = hopp(j,i) * exp(-3.37_dp*((dlij + aG/sqrt(3.0_dp))/(aG/sqrt(3.0_dp))-1.0_dp))
                !end if
              !end if
           end do
        end do
        !$OMP END PARALLEL DO
      end do


      !print*, 'eps and fprime'
      !print*, epsxx(1), epsyy(1), epsxy(1)
      !print*, epsxx(36), epsyy(36), epsxy(36)
      !print*, fprimex(1), fprimey(1)
      !print*, fprimex(36), fprimey(36)
      !d_lij   = round(1e6 * 1/1.42 * ( epsxx .*xij.^2 + epsyy.*yij.^2 +2*epsxy.*xij.*yij)) * 1e-6;
      !d_tij   = d_tij + t*exp(-3.37*((d_lij+1.42)/1.42-1))-t;

      !do ii=1,totImp
      !   !$OMP PARALLEL DO PRIVATE(dlij)
      !   do i=1,nAt
      !      H0(i) = H0(i) + onsiteShift(i)
      !      do j=1,Nneigh(i)
      !        !if (NeighD(3,j,i) .lt. 0.01_dp) then
      !          !dlij = NINT(10**6 * 1.0_dp/(aG/sqrt(3.0_dp)) * (epsxx(i) * NeighD(1,j,i)**2.0_dp + epsyy(i) * NeighD(2,j,i)**2.0_dp + 2.0_dp * epsxy(i) * NeighD(1,j,i)*NeighD(2,j,i))) * 10**(-6)
      !          dlij = 1.0_dp/(aG/sqrt(3.0_dp)) * (epsxx(j,i) * NeighD(1,j,i)**2.0_dp + epsyy(j,i) * NeighD(2,j,i)**2.0_dp + 2.0_dp * epsxy(j,i) * NeighD(1,j,i)*NeighD(2,j,i))
      !          !print*, dlij
      !          !if (abs(real(hopp(j,i)))<0.001) then
      !          !    hopp(j,i) = exp(-3.37_dp*((dlij + aG/sqrt(3.0_dp))/(aG/sqrt(3.0_dp))-1.0_dp))-1.0_dp
      !          !else
      !           hopp(j,i) = hopp(j,i) * exp(-3.37_dp*((dlij + aG/sqrt(3.0_dp))/(aG/sqrt(3.0_dp))-1.0_dp))
      !          !end if
      !        !end if
      !      end do
      !   end do
      !   !$OMP END PARALLEL DO
      !end do

      call MIO_InputParameter('WriteDataFiles',l,.false.)
      if (l) then
         open(1,FILE='e')
         do i=1,nAt
            write(1,*) H0(i)
         end do
         close(1)
      end if
   end if

   if (realisticBubbles) then
      call MIO_Deallocate(indxImp,'indxImp','gauss')
      call MIO_Deallocate(bubbleCenterX,'bubbleCenterX','gauss')
      call MIO_Deallocate(bubbleCenterY,'bubbleCenterY','gauss')
      call MIO_Deallocate(epsxy,'epsxy','ham')
      call MIO_Deallocate(epsxx,'epsxx','ham')
      call MIO_Deallocate(epsyy,'epsyy','ham')
      call MIO_Deallocate(onsiteShift,'onsiteShift','ham')
      !print*, dlij
   end if

   if (magfield) then
      if (Frank) then
         if (frac) call AtomsSetCart()
         do i=1,nAt
            do j=1,Nneigh(i)
               diffx = neighD(1,j,i)
               diffy = neighD(2,j,i)
               if (abs(diffy).lt.0.1d0 .and. diffx .gt. 0.1_dp) then
                   lllll = 1
                   nnnnn = 3
               else if (abs(diffy).lt.0.1d0 .and. diffx .lt. -0.1_dp) then
                   lllll = 3
                   nnnnn = 1
               else if (diffy.gt.0.1d0 .and. diffx .lt. -0.1_dp) then
                   lllll = 1
                   nnnnn = 2
               else if (diffy.lt.-0.1d0 .and. diffx .gt. 0.1_dp) then
                   lllll = 2
                   nnnnn = 1
               else if (diffy.lt.-0.1d0 .and. diffx .lt. -0.1_dp) then
                   lllll = 1
                   nnnnn = 2
               else if (diffy.lt.0.1d0 .and. diffx .gt. 0.1_dp) then
                   lllll = 2
                   nnnnn = 1
               else if (diffy.gt.0.1d0 .and. diffx .gt. 0.1_dp) then
                   lllll = 3
                   nnnnn = 4
               else if (diffy.lt.-0.1d0 .and. diffx .lt. -0.1_dp) then
                   lllll = 4
                   nnnnn = 3
               else if (diffy.lt.-0.1d0 .and. diffx .gt. 0.1_dp) then
                   lllll = 3
                   nnnnn = 4
               else if (diffy.lt.0.1d0 .and. diffx .lt. -0.1_dp) then
                   lllll = 3
                   nnnnn = 4
               end if
               !n2 = norm(ucell(:,2))/sqrt(3.0_dp)*3.0_dp/aG ! n2 = N_y
               n2 = ucell(2,2)/aG ! n2 = N_y
               n1 = ucell(1,1)/(3.0_dp*aG/sqrt(3.0_dp))/2.0_dp  ! n1 = Nx/2
               if (lllll .eq. 1 .and. nnnnn .eq. 3) mmphi = -1.d0 / n2 * (j-1)
               if (lllll .eq. 1 .and. nnnnn .eq. 2 .and. diffy .ge. 0.d0) mmphi = -1.d0 / n1 * (i-1)
               if (lllll .eq. 1 .and. nnnnn .eq. 2 .and. diffy .lt. 0.d0) mmphi = 0.d0
               if (lllll .eq. 2 .and. nnnnn .eq. 1 .and. diffy .ge. 0.d0) mmphi = 0.d0
               if (lllll .eq. 2 .and. nnnnn .eq. 1 .and. diffy .lt. 0.d0) mmphi = 1.d0 / n1 * (i-1)
               if (lllll .eq. 2 .and. nnnnn .eq. 0) mmphi= 1.d0 / n2 * (j-1)
               if (lllll .eq. 3 .and. nnnnn .eq. 0 .and. diffy .ge. 0.d0) mmphi= 0.d0
               if (lllll .eq. 3 .and. nnnnn .eq. 0 .and. diffy .lt. 0.d0) mmphi= 1.d0 / n1 * (i-1)+1.d0/ n1 / 2.d0
               if (lllll .eq. 3 .and. nnnnn .eq. 1) mmphi= 1.d0 / n2 * (j-1)
               if (lllll .eq. 4 .and. nnnnn .eq. 2) mmphi= -1.d0 / n2 * (j-1)
               if (lllll .eq. 4 .and. nnnnn .eq. 3 .and. diffy .ge. 0.d0) mmphi= -1.d0 / n1 * (i-1)-1.d0/ n1 / 2.d0
               if (lllll .eq. 4 .and. nnnnn .eq. 3 .and. diffy .lt. 0.d0) mmphi= 0.d0
               !flux = Bmag*pi/fluxq
               !print*, "flux, B, pi, fluxq", flux, Bmag, pi, fluxq
               print*, mB, n2, n1
               print*, "Magnetic field using Franks approach: ", dble(mB) / dble(n2) / dble(n1) * 39471.80806616257_dp, "T"
               phase = -mmphi*2.0_dp*pi*mB
               !print*, "phase= ", phase, v1, v2
               hopp(j,i) = hopp(j,i)*exp(cmplx_i*phase)
            end do
         end do
      else
         aCC = aG/sqrt(3.0_dp)
         if (frac) call AtomsSetCart()
         !call MIO_InputParameter('Neigh.fastNN',l,.false.)
         !!$OMP PARALLEL DO PRIVATE(v1,v2,phase)
         do i=1,nAt
            do j=1,Nneigh(i)
               !if (abs(neighCell(1,j,i)).gt.1 .or. abs(neighCell(2,j,i)).gt.1) print*, i, neighCell(1,j,i), neighCell(2,j,i)
               v1 = Rat(:,i) + Rat(:,NList(j,i))
               !dist2 = (Rat(1,i)-Rat(1,j))**2 + (Rat(2,i)-Rat(2,j))**2
               !if (dist2.gt.(aCC*1.2)**2 .and. neighD(1,j,i).lt.0 ) then
               !    n1minm1 = 1.0_dp
               !else if (dist2.gt.(aCC*1.2)**2 .and. neighD(1,j,i).gt.0 ) then
               !    n1minm1 = -1.0_dp
               !else
               !    n1minm1 = 0.0_dp
               !end if
               !if (dist2.gt.(aCC*1.2)**2 .and. neighD(2,j,i).lt.0 ) then
               !    n2minm2 = 1.0_dp
               !else if (dist2.gt.(aCC*1.2)**2 .and. neighD(2,j,i).gt.0 ) then
               !    n2minm2 = -1.0_dp
               !else
               !    n2minm2 = 0.0_dp
               !end if
               !v2 = n1minm1*ucell(:,1) + n2minm2*ucell(:,2)
               v2 = neighCell(1,j,i)*ucell(:,1) + neighCell(2,j,i)*ucell(:,2) ! neighCell(1...) contains n1-m1, neighCell(2...) contains n2-m2 (see Cresti's notes)
               v2 = CrossProd(v1,v2)
               v1 = CrossProd(Rat(:,i),Rat(:,NList(j,i))) + v2 ! Implementation of third expersion in Cresti's notes
               phase = flux*v1(3)/1.0d20
               hopp(j,i) = hopp(j,i)*exp(cmplx_i*phase)
            end do
         end do
         !!$OMP END PARALLEL DO
      end if
   end if

   call MIO_InputParameter('HaldaneNNN',l,.false.)
   if (l) then
      if (frac) call AtomsSetCart()
      call MIO_InputParameter('HaldaneSpecifyFlux',l,.false.)
      call MIO_InputParameter('HaldaneSpecifyPhase',ll,.false.)
      call MIO_InputParameter('HaldaneSpecifyRange',lll,.false.)
      call MIO_InputParameter('paperOrientation',paperOrientation,.false.)
      call MIO_InputParameter('HaldaneT2Complex',t2Complex,0.0_dp)
      call MIO_InputParameter('HaldaneT2',t2, 0.0_dp)
      call MIO_InputParameter('HaldaneOppositePhase',HaldaneOppositePhase,.false.)
      ! Note: HaldaneBothLayers is deprecated. Use HaldaneLayerControl and HaldaneLayers instead.
      ! For backward compatibility, we still check if it's set, but it's ignored if HaldaneLayerControl is enabled.
      call MIO_InputParameter('HaldaneBothLayers',HaldaneBothLayers_deprecated,.false.)
      if (HaldaneBothLayers_deprecated .and. .not. HaldaneLayerControl) then
         call MIO_Print('WARNING: HaldaneBothLayers is deprecated. Use HaldaneLayerControl and HaldaneLayers instead.','ham')
         call MIO_Print('  For backward compatibility, applying Haldane to all layers.','ham')
      end if
      !if (l .eq. ll) then
      !   call MIO_Kill('You must specify either the flux or the phase','Haldane','HamHopping')
      !end if
      if (l) then
          call MIO_InputParameter('HaldaneFlux',flux,0.0_dp)
          call MIO_InputParameter('HaldaneSetFluxQ',l,.false.)
          if (l) then
             HaldanePhase = 2.0_dp*pi*flux/1.0471975512_dp
          else
             HaldanePhase = 2.0_dp*pi*flux/fluxq
          end if
      else if (ll) then
          call MIO_InputParameter('HaldanePhase',HaldanePhase,0.0_dp)
      else if (lll) then
          HaldanePhase = HaldanePhase
      end if
      if (HaldaneOppositePhase) then
        HaldanePhase = -HaldanePhase
      end if
      print*, "Haldane phase = ", HaldanePhase
      print*, "HaldaneT2 = ", t2, " eV"
      print*, "Expected gap = ", 2.0_dp*sqrt(3.0_dp)*3.0_dp*t2*sin(HaldanePhase), " eV (according to PRL 106, 236804)"
      !print*, "flux = ", flux
      !print*, "fluxq = ", fluxq
      !print*, "cmplx_i = ", cmplx_i
      !print*, "pi = ", pi
      print*, Nneigh(1)
      !$OMP PARALLEL DO PRIVATE(d)
      do i=1,nAt
         ! Check if Haldane should be applied to this layer
         if (.not. HaldaneEnabledForLayer(layerIndex(i))) then
            cycle
         else
            do j=1,Nneigh(i)
               if (layerIndex(NList(j,i)).eq.layerIndex(i)) then
                  d = sqrt(dot_product(NeighD(1:2,j,i),NeighD(1:2,j,i)))
                  if ((d-0.1).lt.aG .and. (d+0.1).gt.aG) then
                    if (paperOrientation) then
                         if (Species(i).eq.1) then   ! ---  B sublattice (see convention in my notes)
                            if ((NeighD(1,j,i) .gt. 0.1_dp .and. abs(NeighD(2,j,i)) .lt. 0.01_dp) .or. &
                                (NeighD(2,j,i) .lt. -0.1_dp .and. NeighD(1,j,i) .lt. -0.1_dp) .or. &
                                (NeighD(2,j,i) .gt. 0.1_dp .and. NeighD(1,j,i) .lt. 0.1_dp)) then
                                    !hopp(j,i) = hopp(j,i) + cmplx_i * t2Complex
                                    hopp(j,i) = hopp(j,i)*exp(cmplx_i*HaldanePhase)
                                    !print*, "hi1"
                            else
                                    !hopp(j,i) = hopp(j,i) - cmplx_i * t2Complex
                                    hopp(j,i) = hopp(j,i)*exp(-cmplx_i*HaldanePhase)
                                    !print*, "hi2"
                            end if
                         else if (Species(i).eq.2) then   ! ---  A sublattice
                            if ((NeighD(1,j,i) .lt. -0.1_dp .and. abs(NeighD(2,j,i)) .lt. 0.01_dp) .or. &
                                (NeighD(2,j,i) .gt. 0.1_dp .and. NeighD(1,j,i) .gt. 0.1_dp) .or. &
                                (NeighD(2,j,i) .lt. -0.1_dp .and. NeighD(1,j,i) .gt. -0.1_dp)) then
                                    hopp(j,i) = hopp(j,i)*exp(cmplx_i*HaldanePhase)
                                    !hopp(j,i) = hopp(j,i) + cmplx_i * t2Complex
                                    !print*, "hi3"
                            else
                                    hopp(j,i) = hopp(j,i)*exp(-cmplx_i*HaldanePhase)
                                    !hopp(j,i) = hopp(j,i) - cmplx_i * t2Complex
                                    !print*, "hi4"
                            end if
                         end if
                     else ! when 60 degree rotaton compared to PRL paper
                         if (Species(i).eq.2) then   ! ---  B sublattice (see convention in my notes)
                            if ((NeighD(2,j,i) .gt. 0.1_dp .and. abs(NeighD(1,j,i)) .lt. 0.01_dp) .or. &
                                (NeighD(1,j,i) .lt. -0.1_dp .and. NeighD(2,j,i) .lt. -0.1_dp) .or. &
                                (NeighD(1,j,i) .gt. 0.1_dp .and. NeighD(2,j,i) .lt. -0.1_dp)) then
                                    !hopp(j,i) = hopp(j,i) + cmplx_i * t2Complex
                                    hopp(j,i) = hopp(j,i) + t2 * exp(cmplx_i*HaldanePhase)
                                    !print*, "hi1"
                            else
                                    !hopp(j,i) = hopp(j,i) - cmplx_i * t2Complex
                                    hopp(j,i) = hopp(j,i) + t2 * exp(-cmplx_i*HaldanePhase)
                                    !print*, "hi2"
                            end if
                         else if (Species(i).eq.1) then   ! ---  A sublattice
                            if ((NeighD(2,j,i) .lt. -0.1_dp .and. abs(NeighD(1,j,i)) .lt. 0.01_dp) .or. &
                                (NeighD(1,j,i) .gt. 0.1_dp .and. NeighD(2,j,i) .gt. 0.1_dp) .or. &
                                (NeighD(1,j,i) .lt. -0.1_dp .and. NeighD(2,j,i) .gt. 0.1_dp)) then
                                    hopp(j,i) = hopp(j,i) + t2*exp(cmplx_i*HaldanePhase)
                                    !hopp(j,i) = hopp(j,i) + cmplx_i * t2Complex
                                    !print*, "hi3"
                            else
                                    hopp(j,i) = hopp(j,i)+ t2*exp(-cmplx_i*HaldanePhase)
                                    !hopp(j,i) = hopp(j,i) - cmplx_i * t2Complex
                                    !print*, "hi4"
                            end if
                         end if
                     end if
                  end if
               end if
            end do
         end if
      end do
      !$OMP END PARALLEL DO
   end if

   call MIO_InputParameter('WriteDataFiles',w,.false.)
   if (w) then
      !print*, "prefix is ", prefix
      !print*, "here it is only printing the first 3 neighbor hoppings!!!"
      ! One record holds every hopping of an atom; the default record length
      ! (maxlinel) is too short for that.
      call file%Open(name=trim(prefix)//'.'//'s.mag',serial=.true.,recl=64*maxNeigh)
      !open(1,FILE='s')
      u = file%GetUnit()
      do i=1,nAt
         write(u,*) (hopp(j,i),j=1,Nneigh(i))
         !write(u,*) (hopp(1,i))
         !write(u,*) (hopp(2,i))
         !write(u,*) (hopp(3,i))
      end do
      call file%Close()
      !close(1)
   end if
   call MIO_InputParameter('WriteDataFiles',l,.false.)
   if (l) then
      call file%Open(name=trim(prefix)//'.'//'e',serial=.true.)
      !open(1,FILE='e')
      u = file%GetUnit()
      do i=1,nAt
         write(u,*) H0(i), Species(i), layerIndex(i)
      end do
      call file%Close()
      call file%Open(name=trim(prefix)//'.'//'HABreal',serial=.true.)
      !open(1,FILE='e')
      u = file%GetUnit()
      do i=1,nAt
         write(u,*) HABreal(i)
      end do
      call file%Close()
      call file%Open(name=trim(prefix)//'.'//'HABimag',serial=.true.)
      !open(1,FILE='e')
      u = file%GetUnit()
      do i=1,nAt
         write(u,*) HABimag(i)
      end do
      call file%Close()
      call file%Open(name=trim(prefix)//'.'//'pos',serial=.true.)
      !open(1,FILE='e')
      u = file%GetUnit()
      do i=1,nAt
         write(u,*) (Rat(j,i), j=1,3)
      end do
      call file%Close()
      call file%Open(name=trim(prefix)//'.'//'cell',serial=.true.)
      !open(1,FILE='e')
      u = file%GetUnit()
      do i=1,3
         write(u,*) (ucell(j,i), j=1,3)
      end do
      do i=1,3
         write(u,*) (rcell(j,i), j=1,3)
      end do
      call file%Close()
      call file%Open(name=trim(prefix)//'.'//'bottom.e',serial=.true.)
      !open(1,FILE='e')
      u = file%GetUnit()
      do i=1,nAt
         write(u,*) H0Bottom(i), Species(i)
      end do
      call file%Close()
      call MIO_InputParameter('MoireAddSecondMoire',zz,.false.)
      if (zz) then
          call file%Open(name=trim(prefix)//'.'//'top.e',serial=.true.)
          !open(1,FILE='e')
          u = file%GetUnit()
          do i=1,nAt
             write(u,*) H0Top(i), Species(i)
          end do
          call file%Close()
      end if
      !close(1)
   end if

   ! Zeeman effect is now applied in diag.F90 via ApplySOCtoHamiltonian routine
   ! This allows proper spin-channel-aware application for spin-polarized calculations

   if (magfield) then
      ! Landau level calculations (Peierls phase) go here
      ! (Currently empty - all Landau level code was moved elsewhere)
   end if

#ifdef TIMER
   call MIO_TimerStop('ham')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('HamHopping',1)
#endif /* DEBUG */

end subroutine HamHopping

!> @brief Forward discrete cosine transform of a real vector.
!! @param[in]  n  Vector length
!! @param[in]  d  Input data array
!! @param[out] c  Output cosine coefficients
subroutine cosine_transform_data ( n, d, c )

!*****************************************************************************80
!
!! COSINE_TRANSFORM_DATA does a cosine transform on a vector of data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input, real ( kind = 8 ) D(N), the vector of data.
!
!    Output, real ( kind = 8 ) C(N), the cosine transform coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00

  c(1:n) = 0.0D+00

  do i = 1, n
    do j = 1, n
      c(i) = c(i) + cos ( r8_pi * real ( 2 * j - 1, kind = 8 ) &
        * real ( i - 1, kind = 8 ) / 2.0D+00 / real ( n, kind = 8 ) ) * d(j);
    end do
  end do

  c(1:n) = c(1:n) * sqrt ( 2.0D+00 / real ( n, kind = 8 ) )

  return
end
!> @brief Inverse discrete cosine transform to recover data.
!! @param[in]  n  Vector length
!! @param[in]  c  Cosine coefficients
!! @param[out] d  Reconstructed data array
subroutine cosine_transform_inverse ( n, c, d )

!*****************************************************************************80
!
!! COSINE_TRANSFORM_INVERSE does an inverse cosine transform.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input, real ( kind = 8 ) C(N), the cosine transform coefficients.
!
!    Output, real ( kind = 8 ) D(N), the vector of data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00

  do i = 1, n
    d(i) = c(1) / 2.0D+00
    do j = 2, n
      d(i) = d(i) + cos ( r8_pi * real ( 2 * i - 1, kind = 8 ) &
        * real ( j - 1, kind = 8 ) / 2.0D+00 / real ( n, kind = 8 ) ) * c(j)
    end do
  end do

  d(1:n) = d(1:n) * sqrt ( 2.0D+00 / real ( n, kind = 8 ) )

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

!===============================================================================
! FUNCTION: cross
!
! DESCRIPTION:
!   Computes the cross product of two 3D vectors.
!
! PARAMETERS:
!   a(3): First input vector
!   b(3): Second input vector
!
!> @brief Cross product of two 3D vectors.
!! @param[in]  a First input vector (size 3)
!! @param[in]  b Second input vector (size 3)
!! @return axb Cross product a × b (size 3)
!! @details
!!   axb(1) = a(2)*b(3) - a(3)*b(2)\n
!!   axb(2) = a(3)*b(1) - a(1)*b(3)\n
!!   axb(3) = a(1)*b(2) - a(2)*b(1)
function cross(a,b) result (axb)
   implicit none

   real(dp), dimension(3), intent(in) :: a, b
   real(dp), dimension(3) :: axb

   axb(1) = a(2)*b(3) - a(3)*b(2)
   axb(2) = a(3)*b(1) - a(1)*b(3)
   axb(3) = a(1)*b(2) - a(2)*b(1)

end function cross


!> @brief Check if SOC should be applied to a specific layer
!! @param[in] layerIndex Layer index to check
!! @return True if SOC should be applied to this layer
function SOCEnabledForLayer(layerIndex) result(enabled)
   integer, intent(in) :: layerIndex
   logical :: enabled
   integer :: i

   if (.not. SOCLayerControl) then
      ! If layer control is disabled, apply SOC to all layers
      enabled = .true.
   else
      ! Check if layerIndex is in the SOCLayersArray
      enabled = .false.
      if (allocated(SOCLayersArray) .and. size(SOCLayersArray) > 0) then
         do i = 1, size(SOCLayersArray)
            if (SOCLayersArray(i) == layerIndex) then
               enabled = .true.
               return
            end if
         end do
      end if
   end if
end function SOCEnabledForLayer


!> @brief Check if Haldane should be applied to a specific layer
!! @param[in] layerIndex Layer index to check
!! @return True if Haldane should be applied to this layer
function HaldaneEnabledForLayer(layerIndex) result(enabled)
   integer, intent(in) :: layerIndex
   logical :: enabled
   integer :: i

   if (.not. HaldaneLayerControl) then
      ! If layer control is disabled, apply Haldane to all layers
      enabled = .true.
   else
      ! Check if layerIndex is in the HaldaneLayersArray
      enabled = .false.
      if (allocated(HaldaneLayersArray) .and. size(HaldaneLayersArray) > 0) then
         do i = 1, size(HaldaneLayersArray)
            if (HaldaneLayersArray(i) == layerIndex) then
               enabled = .true.
               return
            end if
         end do
      end if
   end if
end function HaldaneEnabledForLayer


end module ham
