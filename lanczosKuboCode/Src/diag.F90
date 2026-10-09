module diag

   use mio

   implicit none

   PRIVATE

   complex(dp), pointer :: ZWork(:)=>NULL()
   real(dp), pointer :: DWork(:)=>NULL()
   integer :: lwork
   real(dp), save :: Efermi=0.0_dp, Emin=-50.0_dp, Emax=20.0_dp
   real(dp), save :: moireAngle=33.004491598883078_dp
   real(dp), save :: gGridRotationAngle=0.0_dp
   logical, save :: skipGRotation=.false.  ! Skip G-rotation when rcell is already rotated
   integer, save :: tapwNG=180, M_tapw=0
   real(dp), save :: physicalTwistAngle=1.08_dp  ! Physical twist angle between graphene layers (degrees)
     logical, save :: calculateChern=.false.
  integer, save :: nk_chern_x=10, nk_chern_y=10
  real(dp), save :: fermi_energy=0.0_dp
  logical, save :: useHighSymmetryGrid=.false.  ! Use colleague's k-grid definition for high-symmetry points
  logical, save :: addHighSymmetryRefinement=.false.  ! Add extra k-points near high-symmetry points for dispersive bands
  logical, save :: useTriangularTruncation=.false.
  logical, save :: checkTAPWUnitary=.false.  ! Enable TAPW unitary transformation check
  logical, save :: useKprimeValley=.false.  ! Use K' valley instead of K valley for TAPW reference
  logical, save :: tapwDebug=.false.  ! Enable verbose debug output for TAPW calculations
  logical, save :: useRigidPositions=.false.  ! Use rigid reference positions for TAPW X matrix construction
  logical, save :: tapwLowdin=.false.  ! Loewdin-orthonormalise X (X S^-1/2, S = X^H X); relaxed positions make S /= 1
  ! Both valleys in one TAPW basis (opt-in): the G-list is the union of the shells around K and K',
  ! so X^H H X also carries the intervalley block.  Default .false. = the single-valley basis.
  logical, save :: tapwBothValleys=.false.
  logical, save :: tapwValleyDecouple=.false.  ! diagnostic: zero the K-K' block of Hproj (known-answer test)
  integer, save :: tapw_NGvalley1=0            ! G-vectors of the first valley block (G 1..tapw_NGvalley1)
  logical, save :: tbv_announced=.false.
  logical, save :: socDebug=.false.  ! Enable verbose debug output for SOC calculations
  logical, save :: forceBlockTAPW=.false.  ! Force block TAPW path even without Rashba SOC (for testing)
  logical, save :: gammaCentred3D=.true.   ! Gamma-centred i/N grid for 3D TAPW bands (.false. = Monkhorst-Pack)
  logical, save :: calculate3DTAPWBands=.false.  ! Calculate 3D bands (kx, ky, eigenvalue) using Monkhorst-Pack grid in DiagBands
  integer, save :: nk_3D_x=10, nk_3D_y=10    ! Grid size for 3D bands calculation

#ifdef SEMICL
#include "diag_semicl_decl.inc"
#endif
  logical, save :: bf_announced=.false.             ! keep the G-list message to one line
  ! Plane-wave composition of TAPW eigenstates (opt-in diagnostic, Diag.GWeights) -> <prefix>.GWeights:
  ! per k and band in gw_b1..gw_b2, the gw_ntop largest weights sum_label |c_(G,label)|^2.
  logical, save :: gWeightsTAPW=.false.             ! Diag.GWeights
  integer, save :: gw_b1=0, gw_b2=0, gw_ntop=8, gw_unit=-1, gw_ng=0
  logical, save :: gw_open=.false.                  ! newunit= returns NEGATIVE units: test this, not gw_unit
  ! Layer / sublattice composition of TAPW eigenstates (opt-in, Diag.LayerWeights) -> <prefix>.LayerWeights:
  ! per k and band in lw_b1..lw_b2 (0 0 = all), the atom-space weight sum_{i in label} |psi_i|^2, psi = X c.
  logical, save :: layerWeightsTAPW=.false.         ! Diag.LayerWeights
  integer, save :: lw_b1=0, lw_b2=0, lw_unit=-1
  logical, save :: lw_open=.false.
  ! Diag.GWeightsHam (opt-in, needs Diag.GWeights): the projected TAPW Hamiltonian before diagonalisation
  ! -> <prefix>.TAPWHam, for offline first-order (two-level) analysis of avoided crossings.
  logical, save :: gWeightsHam=.false.              ! Diag.GWeightsHam
  logical, save :: gwh_open=.false.
  integer, save :: gwh_unit=-1

  ! TAPW runtime configuration cached outside OpenMP regions
  logical, save :: tapw_cfg_initialized = .false.
  real(dp), save :: tapw_cfg_aG = 2.46019_dp
  real(dp), save :: tapw_cfg_shift = 0.0_dp
  logical, save :: tapw_cfg_useShift = .false.
  logical, save :: tapw_cfg_saveRitz = .false.
  real(dp), save :: tapw_cfg_tol = 0.1_dp

   ! Module variables to store G-vectors from TAPW calculation for Berry curvature
   real(dp), allocatable, save :: saved_Gx(:), saved_Gy(:)
   integer, save :: saved_NG=0, saved_Nlabel=0

   ! Storage for TAPW eigenvectors and Hamiltonians (only when calculateChern=true)
   complex(dp), allocatable, save :: stored_eigenvectors(:,:,:,:)  ! (M, M, nk, nspin)
   real(dp), allocatable, save :: stored_eigenvalues(:,:,:)  ! (M, nk, nspin)
   complex(dp), allocatable, save :: stored_hamiltonians(:,:,:,:)  ! (M, M, nk, nspin)
   integer, save :: stored_M=0, stored_nk=0, stored_nspin=0

   ! PERFORMANCE OPTIMIZATION: Pre-allocated TAPW arrays for reuse across k-points
   complex(dp), allocatable, save, target :: tapw_H_dense(:,:) ! (N, N) - reused across k-points
   complex(dp), allocatable, save, target :: tapw_Hproj(:,:) ! (M, M) - reused across k-points
   real(dp), allocatable, save, target :: tapw_eigvals(:) ! (M) - reused across k-points
   complex(dp), allocatable, save, target :: tapw_ZWorkLoc(:) ! (2*M) - reused across k-points
   real(dp), allocatable, save, target :: tapw_DWorkLoc(:) ! (3*M) - reused across k-points

   ! Storage for TAPW transformation matrix and TB parameters (for Option B Berry curvature)
   complex(dp), allocatable, save :: stored_X_matrix(:,:)  ! (N, M) - TAPW transformation matrix
   integer, save :: stored_N = -1  ! Number of TB orbitals
   real(dp), allocatable, save :: stored_H0(:)  ! On-site energies
   complex(dp), allocatable, save :: stored_hopp(:,:)  ! Hopping parameters
   integer, allocatable, save :: stored_NList(:,:), stored_Nneigh(:), stored_neighCell(:,:,:)
   integer, save :: stored_maxNeigh = -1
   real(dp), allocatable, save :: stored_Kpts(:,:)  ! (3, nk) - k-point coordinates

   ! Cache variables for position difference matrices (module-level)
   real(dp), allocatable, save :: cached_delX(:,:), cached_delY(:,:)
   integer, save :: cached_M = -1, cached_NG = -1, cached_Nlabel = -1
   logical, save :: cache_valid = .false.
   real(dp), allocatable, save :: cached_Gx(:), cached_Gy(:)

   ! Flag to track if any SOC terms are enabled (for optimization)
   logical, save :: anySOCEnabled = .false.

   public :: DiagInit
   public :: DiagDOS, Diag3DBands, DiagPDOS, DiagChern
   public :: DiagBands, DiagBandsG, DiagBandsAroundK, DiagBandsRashba
   public :: DiagSpectralFunction, DiagSpectralFunctionKGrid
   public :: DiagSpectralFunctionKGridInequivalent
   public :: DiagSpectralFunctionKGridInequivalent_v2
   public :: DiagSpectralFunctionKGridInequivalentEnergyCut
   public :: DiagSpectralFunctionKGridInequivalentEnergyCut_v2
   public :: DiagSpectralFunctionKGridInequivalentEnergyCutNickDale
   public :: moireAngle, gGridRotationAngle, skipGRotation, tapwNG, M_tapw, calculateChern, nk_chern_x, nk_chern_y, fermi_energy, useTriangularTruncation, checkTAPWUnitary, physicalTwistAngle, useKprimeValley, tapwDebug, useRigidPositions, tapwLowdin, tapwBothValleys, tapwValleyDecouple, socDebug, forceBlockTAPW
   public :: calculate3DTAPWBands, nk_3D_x, nk_3D_y, gammaCentred3D
#ifdef SEMICL
   public :: berryFluxTAPW, berryBandMin, berryBandMax
#endif
#ifdef SEMICL
   public :: orbMomentTAPW, om_degtol, om_decomp, om_nearcut
#endif
#ifdef SEMICL
   public :: orbSuscTAPW
#endif
   public :: gWeightsTAPW, gw_b1, gw_b2, gw_ntop, gWeightsHam
   public :: layerWeightsTAPW, lw_b1, lw_b2
#ifdef SEMICL
   public :: berryLinksTAPW
#endif

contains

subroutine DiagInit(N)

   use ham, only : Zterm, IntrinsicSOCterm, IsingSOCterm, PIASOCterm, RashbaSOCterm

   integer, intent(in) :: N

   complex(dp) :: A(1,1), OPT(1)
   real(dp) :: W(1), W2(1)
   integer :: INFO

#ifdef DEBUG
   call MIO_Debug('DiagInit',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   ! Set global flag if any SOC term is enabled
   anySOCEnabled = Zterm .or. IntrinsicSOCterm .or. IsingSOCterm .or. PIASOCterm .or. RashbaSOCterm

   call ZHEEV('N','L',N,A,N,W,OPT,-1,W2,INFO)
   if (INFO /= 0) call MIO_Kill('Error in workspace query for diagonalization','diag','DiagInit')
   lwork = int(OPT(1))
   call MIO_Allocate(ZWork,lwork,'ZWork','diag')
   call MIO_Allocate(DWork,3*N-2,'DWork','diag')

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagInit',1)
#endif /* DEBUG */

end subroutine DiagInit

subroutine DiagDOS()

   use cell,                 only : rcell, ucell
   use atoms,                only : nEl, nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : twopi
   use math

   integer, parameter :: intorder=5

   integer :: nk(3), ptot, i1, i2, i3, ik, Epts, u, is, uu
   complex(dp), pointer :: H(:,:,:)=>NULL()
   real(dp), pointer :: Kgrid(:,:)=>NULL(), Eig(:,:)=>NULL(), DOS(:,:)=>NULL(), E(:)=>NULL()
   real(dp), pointer :: EStore(:,:)=>NULL()
   type(cl_file) :: file, file2
   character(len=100) :: flnm, flnm2
   real(dp) :: eps, E1, E2, s, sp, Ep

   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   complex(dp) :: HLoc(nAt, nAt)

#ifdef DEBUG
   call MIO_Debug('DiagDOS',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_Print('Calculating DOS by diagonalization','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1])
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
   call MIO_InputParameter('DOS.Emin',E1,-10.0_dp)
   call MIO_InputParameter('DOS.Emax',E2,10.0_dp)
   call MIO_Allocate(DOS,[Epts,nspin],'DOS','diag')
   DOS = 0.0d0
   call MIO_Allocate(E,Epts,'E','diag')
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(EStore,[ptot,nAt],'DOS','diag')
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      Kgrid(:,ik) = rcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
        rcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + rcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
   end do; end do; end do
   call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
   call MIO_Allocate(Eig,[nAt,nspin],'Eig','diag')
   flnm = trim(prefix)//'.diag.DOS'
   call file%Open(name=flnm,serial=.true.)
   u = file%GetUnit()
   Emax = -huge(1.0_dp)
   Emin = huge(1.0_dp)
   do i2=1,Epts
        E(i2) = E1 + (E2-E1)*(i2-1)/(Epts-1)
   end do
   do is=1,nspin
      if (nspin > 1) then
         call MIO_Print('Diagonalizing spin '//trim(num2str(is))//' of '//trim(num2str(nspin))//': Processing '//trim(num2str(ptot))//' k-points','diag')
      else
         call MIO_Print('Diagonalizing: Processing '//trim(num2str(ptot))//' k-points','diag')
      end if
      !$OMP PARALLEL DO PRIVATE(ELoc, KptsLoc, ik, HLoc), REDUCTION(min:Emin), REDUCTION(max:Emax), &
      !$OMP& SHARED(E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, ptot)
      do ik=1,ptot
         ! Progress reporting every 10% completion
         if (ptot > 1) then
            if (ik == ptot .or. (ik > 0 .and. int(10.0_dp*(ik-1)/ptot) < int(10.0_dp*ik/ptot))) then
               !$OMP CRITICAL
               call MIO_Print('  k-point '//trim(num2str(ik))//' of '//trim(num2str(ptot))//' ('//trim(num2str(int(100.0_dp*ik/ptot)))//'%)','diag')
               !$OMP END CRITICAL
            end if
         end if
         KptsLoc = Kgrid(:,ik)
         ELoc = 0.0_dp
         HLoc = 0.0_dp
         !if (modulo(ik,int(ptot/10)).eq.0) print*, "progress is: ", ik/int(ptot/10)*10, "percent"
         call DiagHam(nAt,nspin,is,HLoc,ELoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
         !Eig(:,is) = Eig(:,is)*g0
         ELoc = ELoc*g0
         !      DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-Eig(i1,is))**2/(2.0_dp*eps**2))
         Emin = min(Emin,ELoc(1))
         Emax = max(Emax,ELoc(nAt))
         !      DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-ELoc(i1))**2/(2.0_dp*eps**2))
         do i1=1,nAt
            EStore(ik,i1) = ELoc(i1)
         end do
      end do
      !$OMP END PARALLEL DO
      call MIO_Print('  Completed diagonalization for spin '//trim(num2str(is))//' of '//trim(num2str(nspin)),'diag')
   end do
   call MIO_Print('Accumulating DOS from eigenvalues...','diag')
   do is=1,nspin
      if (nspin > 1) then
         call MIO_Print('  Processing spin '//trim(num2str(is))//' of '//trim(num2str(nspin)),'diag')
      end if
      do ik = 1,ptot
         ! Progress reporting every 10% completion
         if (ptot > 1) then
            if (ik == ptot .or. (ik > 0 .and. int(10.0_dp*(ik-1)/ptot) < int(10.0_dp*ik/ptot))) then
               call MIO_Print('    Accumulating k-point '//trim(num2str(ik))//' of '//trim(num2str(ptot))//' ('//trim(num2str(int(100.0_dp*ik/ptot)))//'%)','diag')
            end if
         end if
         do i1=1,nAt
            do i2=1,Epts
               DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
            end do
         end do
      end do
   end do
   call MIO_Print('Completed DOS accumulation','diag')
   if (nspin==1) then
      DOS = 2.0_dp*DOS/(eps*sqrt(twopi)*ptot)
   else
      DOS = DOS/(eps*sqrt(twopi)*ptot)
   end if
   do i1=1,Epts
      write(u,*) E(i1), (DOS(i1,is),is=1,nspin)
   end do
   call file%Close()
   call MIO_Deallocate(Eig,'Eig','diag')
   call MIO_Deallocate(H,'H','diag')
   call MIO_Print('Emin: '//trim(num2str(Emin,4))//', Emax: '//trim(num2str(Emax,4)),'diag')
   eps = (E2-E1)/(Epts-1)
   sp = 0.0_dp
   Ep = Emin
   do ik=1,Epts-2*intorder
      s = 0.0_dp
      do is=1,nspin
         s = s + TrapezoidalInt(DOS(:intorder*2+ik,is),intorder*2+ik,eps,intorder)
      end do
      if (s>=nEl) then
         Efermi = (E(intorder*2+ik)+Ep)/2.0_dp
         exit
      else
         if (s/=sp) then
            Ep = E(intorder*2+ik)
            sp = s
         end if
      end if
   end do
   call MIO_Print('Efermi: '//trim(num2str(Efermi,5)),'diag')
   call MIO_Print('')

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagDOS',1)
#endif /* DEBUG */
   !

end subroutine DiagDOS

subroutine DiagPDOS()

   use cell,                 only : rcell, ucell
   use atoms,                only : nEl, nAt, layerIndex
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : twopi
   use math

   integer, parameter :: intorder=5

   integer :: nk(3), ptot, i1, i2, i3, ik, Epts, u, is, uu, i1bis, ivec, ivec2, PDOSLayerIndex, numberOfLayers
   complex(dp), pointer :: H(:,:,:)=>NULL()
   real(dp), pointer :: Kgrid(:,:)=>NULL(), Eig(:,:)=>NULL()
   real(dp), pointer :: DOS(:,:)=>NULL(), E(:)=>NULL()
   real(dp), pointer :: DOS_thread(:,:,:)=>NULL()
   real(dp), pointer :: EStore(:,:)=>NULL()
   type(cl_file) :: file, file2
   character(len=100) :: flnm, flnm2
   character(len=2) :: fmt, x1
   real(dp) :: eps, E1, E2, s, sp, Ep

   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   complex(dp) :: HLoc(nAt, nAt)
   real(dp) :: pipj, vectormultip
   complex(dp), pointer :: HStore(:,:,:)=>NULL()

   integer :: omp_get_thread_num, omp_get_max_threads, index_ii

#ifdef DEBUG
   call MIO_Debug('DiagPDOS',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_Print('Calculating PDOS by diagonalization','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1])
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
   call MIO_InputParameter('DOS.Emin',E1,-10.0_dp)
   call MIO_InputParameter('DOS.Emax',E2,10.0_dp)
   call MIO_Allocate(E,Epts,'E','diag')
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(EStore,[ptot,nAt],'DOS','diag')
   call MIO_Allocate(HStore,[ptot,nAt,nAt],'DOS','diag')
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      Kgrid(:,ik) = rcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
        rcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + rcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
   end do; end do; end do
   call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
   call MIO_Allocate(Eig,[nAt,nspin],'Eig','diag')
   Emax = -huge(1.0_dp)
   Emin = huge(1.0_dp)
   do i2=1,Epts
        E(i2) = E1 + (E2-E1)*(i2-1)/(Epts-1)
   end do
   do is=1,nspin
      !$OMP PARALLEL DO PRIVATE(ELoc, KptsLoc, ik, HLoc), REDUCTION(min:Emin), REDUCTION(max:Emax), &
      !$OMP& SHARED(E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
      do ik=1,ptot
         KptsLoc = Kgrid(:,ik)
         ELoc = 0.0_dp
         HLoc = 0.0_dp
         !if (modulo(ik,int(ptot/10)).eq.0) print*, "progress is: ", ik/int(ptot/10)*10, "percent"
         call DiagHamPDOS(nAt,nspin,is,HLoc,ELoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
         !Eig(:,is) = Eig(:,is)*g0
         ELoc = ELoc*g0
         !      DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-Eig(i1,is))**2/(2.0_dp*eps**2))
         Emin = min(Emin,ELoc(1))
         Emax = max(Emax,ELoc(nAt))
         !      DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-ELoc(i1))**2/(2.0_dp*eps**2))
         do i1=1,nAt
            EStore(ik,i1) = ELoc(i1)
            do i1bis=1,nAt
              HStore(ik,i1,i1bis)=HLoc(i1,i1bis)
            end do
         end do
      end do
      !$OMP END PARALLEL DO
   end do
   call MIO_InputParameter('numberOfLayers',numberOfLayers,2)
   do PDOSLayerIndex=1,numberOfLayers
      flnm = trim(prefix)//'.diag.DOS.Layer'//trim(num2str(PDOSLayerIndex))
      call file%Open(name=flnm,serial=.true.)
      u = file%GetUnit()
      call MIO_Allocate(DOS,[Epts,nspin],'DOS','diag')
      call MIO_Allocate(DOS_thread,[Epts,nspin,omp_get_max_threads()],'DOS','diag')
      DOS = 0.0_dp
      DOS_thread = 0.0_dp
      do is=1,nspin
         do ik = 1,ptot ! Loop over k
            do i1=1,nAt ! Loop over atoms (same as number of bands here)
               if (layerIndex(i1).eq.PDOSLayerIndex) then
                  do i2=1,Epts
                     DOS(i2,is) = 0.0_dp
                     do ivec=1,nAt ! add this for the vector multiplication projection operator
                         !$OMP PARALLEL DO PRIVATE(vectormultip, pipj), &
                         !$OMP& SHARED(E, EStore, eps, DOS_thread)
                         do ivec2=1,nAt ! add this for the vector multiplication projection operator
                     ! outer product between vector multiplication of two eigenvectors
                     ! and the eigenenergy vector. The outer product is thus over the
                     ! eigenenergies
                           vectormultip=conjg(HStore(ik,ivec,i1))*HStore(ik,ivec2,i1) !
                           pipj=EStore(ik,i1)*vectormultip
                           DOS_thread(i2,is,omp_get_thread_num()+1) = DOS_thread(i2,is,omp_get_thread_num()+1) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2)) * pipj
                           !DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2)) * pipj
                         end do
                         !$OMP END PARALLEL DO
                     end do
                  end do
               end if
            end do
         end do
      end do
      ! TO DO, add OMP parallelisaton here if it works without
      DO index_ii = 1,omp_get_max_threads()
         do is=1,nspin
            do i2=1,Epts
               DOS(i2,is) = DOS(i2,is) + DOS_thread(i2,is,index_ii)
            enddo
         enddo
      ENDDO
      call MIO_Deallocate(DOS_thread,'DOS','diag')
      if (nspin==1) then
         DOS = 2.0_dp*DOS/(eps*sqrt(twopi)*ptot)
      else
         DOS = DOS/(eps*sqrt(twopi)*ptot)
      end if
      do i1=1,Epts
         write(u,*) E(i1), (DOS(i1,is),is=1,nspin)
      end do
      call file%Close()
      call MIO_Deallocate(Eig,'Eig','diag')
      call MIO_Deallocate(H,'H','diag')
      call MIO_Print('Emin: '//trim(num2str(Emin,4))//', Emax: '//trim(num2str(Emax,4)),'diag')
      eps = (E2-E1)/(Epts-1)
      sp = 0.0_dp
      Ep = Emin
      do ik=1,Epts-2*intorder
         s = 0.0_dp
         do is=1,nspin
            s = s + TrapezoidalInt(DOS(:intorder*2+ik,is),intorder*2+ik,eps,intorder)
         end do
         if (s>=nEl) then
            Efermi = (E(intorder*2+ik)+Ep)/2.0_dp
            exit
         else
            if (s/=sp) then
               Ep = E(intorder*2+ik)
               sp = s
            end if
         end if
      end do
      call MIO_Print('Efermi: '//trim(num2str(Efermi,5)),'diag')
      call MIO_Print('')
      call MIO_Deallocate(DOS,'DOS','diag')
   end do

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagPDOS',1)
#endif /* DEBUG */
   !

end subroutine DiagPDOS

subroutine Diag3DBands()

   use cell,                 only : rcell, ucell
   use atoms,                only : nEl, nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : twopi
   use math

   integer, parameter :: intorder=5

   integer :: nk(3), ptot, i1, i2, i3, ik, Epts, u, is, uu, i
   complex(dp), pointer :: H(:,:,:)=>NULL()
   real(dp), pointer :: Kgrid(:,:)=>NULL(), Eig(:,:)=>NULL(), DOS(:,:)=>NULL(), E(:)=>NULL()
   type(cl_file) :: file
   character(len=100) :: flnm, flnm2
   real(dp) :: eps, E1, E2, s, sp, Ep

#ifdef DEBUG
   call MIO_Debug('Diag3DBands',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_Print('Calculating 3D Bands by diagonalization','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1])
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
   call MIO_InputParameter('DOS.Emin',E1,-10.0_dp)
   call MIO_InputParameter('DOS.Emax',E2,10.0_dp)
   call MIO_Allocate(DOS,[Epts,nspin],'DOS','diag')
   call MIO_Allocate(E,Epts,'E','diag')
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      Kgrid(:,ik) = rcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
        rcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + rcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
   end do; end do; end do
   call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
   call MIO_Allocate(Eig,[nAt,nspin],'Eig','diag')
   !call file%Open(name=flnm,serial=.true.)

   flnm2 = trim(prefix)//'.3Dbands'
   uu=98
   open(uu,FILE=flnm2,STATUS='replace')
   write(uu,'(3i8)') nAt, nspin, ptot
   do is=1,nspin
      do ik=1,ptot
         if (modulo(ik,int(ptot/10)).eq.0) call MIO_Print('progress is: '//trim(num2str((ik/ptot/10.0_dp*10.0_dp),4))//' percent','diag')
         call DiagHam(nAt,nspin,is,H(:,:,is),Eig(:,is),Kgrid(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
         !hv = max(hv,maxval(E(:,is),mask=E(:,is)<=Efermi/g0))
         !lc = min(lc,minval(E(:,is),mask=E(:,is)>Efermi/g0))
         write(uu,'(f12.6,10f14.6,/,(10x,10f14.6))') Kgrid(:,ik),(Eig(i,is)*g0, i=1,nAt)

         Eig(:,is) = Eig(:,is)*g0
         Emin = min(Emin,Eig(1,is))
         Emax = max(Emax,Eig(nAt,is))
         !      DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-Eig(i1,is))**2/(2.0_dp*eps**2))
      end do
   end do
   !call file%Close()
   call MIO_Deallocate(Eig,'Eig','diag')
   call MIO_Deallocate(H,'H','diag')

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('Diag3DBands',1)
#endif /* DEBUG */
   !

end subroutine Diag3DBands

subroutine DiagBands()

   use cell,                 only : rcell, ucell, aG
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin, RashbaSOCterm, IsingSOCterm
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, uu, uuu, uuuu, neig, i_eig
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: H(:,:,:)=>NULL()
   ! Block Hamiltonian arrays for Rashba
   complex(dp), allocatable :: HBlock(:,:)
   real(dp), allocatable :: ELocBlock(:)
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS, useSameNumberOfPoints
   real(dp) :: theta, volume

   real(dp) :: gcell(3,3), grcell(3,3)

   real(dp) :: KptsLoc(3)
   real(dp), allocatable :: ELoc(:)  ! Made allocatable for dynamic TAPW sizing
   ! Allocated only when a dense-H path runs (DiagHam / DiagHamWF). As an automatic nAt x nAt array it was
   ! created, and zeroed at every k, in TAPW and sparse runs that never read it (11 TB at 845 k atoms).
   complex(dp), allocatable :: HLoc(:,:)
   logical :: useDifferentLatticeVectors, keepWaveFunction, sparseDiagSolver, useTAPW, useDenseMatrixTAPW
   real(dp) :: rcell_rotated(3,3)  ! For TAPW G-grid rotation consistency
   real(dp) :: rcell_original(3,3) ! Backup of original rcell for TAPW
   character(len=5) :: useTAPW_str, forceBlockTAPW_str  ! Debug output variables

#ifdef DEBUG
   call MIO_Debug('DiagBands',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_InputParameter('Bands.NumPoints',nPts0,100)
   call MIO_InputParameter('Bands.SparseNeig',neig,100)
   call MIO_InputParameter('Bands.UseDifferentLatticeVectors',useDifferentLatticeVectors,.false.)
   call MIO_InputParameter('Bands.UseSameNumberOfPoints',useSameNumberOfPoints,.false.)
   call MIO_InputParameter('keepWaveFunction',keepWaveFunction,.false.)
    call MIO_InputParameter('sparseDiagSolver',sparseDiagSolver,.false.)
    call MIO_InputParameter('useTAPW',useTAPW,.false.)
    call MIO_InputParameter('useDenseMatrixTAPW',useDenseMatrixTAPW,.false.)
    ! The Zeeman, Ising and Rashba terms are implemented in the TAPW path only.
    if (.not. useTAPW) then
       block
          logical :: zeeman, pzeeman
          call MIO_InputParameter('ZeemanTerm',zeeman,.false.)
          call MIO_InputParameter('PseudoZeemanTerm',pzeeman,.false.)
          if (zeeman .or. pzeeman .or. IsingSOCterm) then
             call MIO_Print('WARNING: ZeemanTerm, PseudoZeemanTerm and IsingSOCterm are implemented for TAPW '// &
               'calculations only (useTAPW .true.); they have NO effect on this calculation.','diag')
          end if
          if (RashbaSOCterm .and. nspin == 1) then
             call MIO_Kill('RashbaSOCterm is implemented for TAPW calculations only (useTAPW .true.).', &
               'diag','DiagBands')
          end if
       end block
    end if

   if (useDifferentLatticeVectors) then
       ucell(1,2) = 0.0_dp
       rcell(2,1) = 0.0_dp
   end if

   if (MIO_InputFindBlock('Bands.Path',nPath)) then
      call MIO_Print('Band calculation','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_InputBlock('Bands.Path',path)

      ! Store fractional coordinates for debugging
      open(unit=98, file='kpath_debug_fractional', status='replace')
      write(98, '(A)') '# Fractional k-path coordinates'
      do ip=1,nPath
         write(98, '(3F16.8)') path(1,ip), path(2,ip), path(3,ip)
      end do
      close(98)

      do ip=1,nPath
         ! Convert fractional coordinates to absolute k-space coordinates
         path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
      end do

      ! Store absolute coordinates for debugging
      open(unit=98, file='kpath_debug_absolute', status='replace')
      write(98, '(A)') '# Absolute k-path coordinates'
      do ip=1,nPath
         write(98, '(3F16.8)') path(1,ip), path(2,ip), path(3,ip)
      end do
      close(98)

      print *, "K-path debug data written to kpath_debug_* files"
         !   !print*, "theta=", theta
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1 ! coz 4 points, 3 segments
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               if (useSameNumberOfPoints) then
                  nPts(ip) = nPts0
                  ptsTot = ptsTot + nPts0
               else
                  nPts(ip) = nint(real(d*nPts0)/real(d0))
                  ptsTot = ptsTot + nPts(ip)
               end if
            end do
            nPts(ip) = nPts(ip) + 1 ! Add the missing point at the end of last segment
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot+1],'Kpts','diag') ! add missing point
      Kpts(:,1) = path(:,1)
      ip = 1
      d = 0.0_dp
      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            Kpts(:,ip) = path(:,i) + (j)*(path(:,i+1)-path(:,i))/nPts(i)
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            d = d + sqrt(dot_product(v,v))
         end do
      end do
      ! For TAPW: allocate larger E array to accommodate projected space M = NG × Nlabel
      if (useTAPW) then
         call MIO_Allocate(E,[max(nAt, 15000),nspin,ptsTot+1],'E','diag')
      else
         call MIO_Allocate(E,[nAt,nspin,ptsTot+1],'E','diag')
      end if
      flnm = trim(prefix)//'.ham'
      uu=101
      open(uu,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.prob'
      uuuu=103
      open(uuuu,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.eig'
      uuu=102
      open(uuu,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.bands'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp

      ! Override k-path with k-grid for Chern calculation or 3D TAPW bands calculation
      if ((calculateChern .and. useTAPW) .or. calculate3DTAPWBands) then
         if (calculateChern .and. useTAPW) then
            call MIO_Print('TAPW Chern number calculation mode detected','diag')
            call MIO_Print('Generating '//trim(num2str(nk_chern_x))//'x'//trim(num2str(nk_chern_y))// &
                          ' Monkhorst-Pack k-grid for Chern calculation','diag')
         else if (calculate3DTAPWBands) then
            call MIO_Print('3D TAPW bands calculation mode detected','diag')
            call MIO_Print('Generating '//trim(num2str(nk_3D_x))//'x'//trim(num2str(nk_3D_y))// &
                          merge(' Gamma-centred  ',' Monkhorst-Pack ',gammaCentred3D)// &
                          'k-grid for 3D TAPW bands calculation','diag')
         end if

         ! Deallocate existing k-path arrays
         if (associated(path)) call MIO_Deallocate(path,'path','diag')
         if (associated(nPts)) call MIO_Deallocate(nPts,'nPts','diag')
         if (associated(Kpts)) call MIO_Deallocate(Kpts,'Kpts','diag')
         if (associated(E)) call MIO_Deallocate(E,'E','diag')

         ! Generate Monkhorst-Pack k-grid
         if (calculateChern .and. useTAPW) then
            ptsTot = nk_chern_x * nk_chern_y
         else if (calculate3DTAPWBands) then
            ptsTot = nk_3D_x * nk_3D_y
         end if
         call MIO_Allocate(Kpts,[3,ptsTot+1],'Kpts','diag')  ! +1 to match regular bands allocation

         ! Allocate E array based on calculation type
         if (calculateChern .and. useTAPW) then
            ! For TAPW: allocate larger E array to accommodate projected space M = NG × Nlabel
            ! Use conservative estimate: max(nAt, 15000) to handle large TAPW projections
            call MIO_Allocate(E,[max(nAt, 15000),nspin,ptsTot+1],'E','diag')  ! +1 to match Kpts
         else
            ! For regular 3D bands: use nAt
            call MIO_Allocate(E,[nAt,nspin,ptsTot+1],'E','diag')  ! +1 to match Kpts
         end if

         ! Generate 2D k-grid using global rcell (already rotated for TAPW)
         ip = 0
         if (calculateChern .and. useTAPW) then
            if (useHighSymmetryGrid) then
               ! Colleague's definition: k = (i/Nk) * b1 + (j/Nk) * b2 with i,j = 0,1,2,...Nk-1
               ! This ensures proper sampling of Gamma, K, and M points
               call MIO_Print('Using high-symmetry k-grid definition: k = (i/Nk)*b1 + (j/Nk)*b2','diag')
               call MIO_Print('  This ensures proper sampling of Gamma (0,0), K, and M points','diag')
               if (tapwDebug) call MIO_Print('DEBUG: rcell matrix:','diag')
               call MIO_Print('  b1 = ['//trim(num2str(rcell(1,1),6))//','//trim(num2str(rcell(2,1),6))//','//trim(num2str(rcell(3,1),6))//']','diag')
               call MIO_Print('  b2 = ['//trim(num2str(rcell(1,2),6))//','//trim(num2str(rcell(2,2),6))//','//trim(num2str(rcell(3,2),6))//']','diag')
               call MIO_Print('  b3 = ['//trim(num2str(rcell(1,3),6))//','//trim(num2str(rcell(2,3),6))//','//trim(num2str(rcell(3,3),6))//']','diag')

               ! First pass: regular grid
               do j = 0, nk_chern_y-1
                  do i = 0, nk_chern_x-1
                     ip = ip + 1
                     Kpts(:,ip) = rcell(:,1)*real(i,dp)/real(nk_chern_x,dp) + &
                                 rcell(:,2)*real(j,dp)/real(nk_chern_y,dp) + &
                                 rcell(:,3)*0.0_dp  ! 2D system: kz = 0
                  end do
               end do

               ! Second pass: Add extra points near high-symmetry points for dispersive bands
               if (addHighSymmetryRefinement) then
                  call MIO_Print('Adding extra k-points near high-symmetry points for dispersive bands','diag')
                  call add_high_symmetry_refinement(Kpts, ip, rcell, nk_chern_x, nk_chern_y)
               end if
            else
               ! Original Monkhorst-Pack formula
               call MIO_Print('Using standard Monkhorst-Pack k-grid definition','diag')
               do j = 1, nk_chern_y
                  do i = 1, nk_chern_x
                     ip = ip + 1
                     ! Use exact same MP formula as Diag3DBands (lines 443-444) with global rcell
                     Kpts(:,ip) = rcell(:,1)*(2*i-nk_chern_x-1)/(2.0_dp*nk_chern_x) + &
                                 rcell(:,2)*(2*j-nk_chern_y-1)/(2.0_dp*nk_chern_y) + &
                                 rcell(:,3)*0.0_dp  ! 2D system: kz = 0
                  end do
               end do
            end if
         else if (calculate3DTAPWBands) then
            if (gammaCentred3D) then
               ! Gamma-centred k = (i/N)*b1 + (j/N)*b2, i,j = 0..N-1.
               ! Unlike Monkhorst-Pack, this samples Gamma exactly, and for N
               ! divisible by 6 also M (i=N/2) and K (i=2N/3).  The MP grid puts
               ! all three exactly half a step off, which leaves a band extremum
               ! at Gamma or K resolved only by the cells straddling it.
               ! It is also the convention a PERIODIC marching-squares contour
               ! tracer needs, since index N wraps exactly onto index 0.
               call MIO_Print('3D TAPW bands: Gamma-centred i/N k-grid','diag')
               do j = 0, nk_3D_y-1
                  do i = 0, nk_3D_x-1
                     ip = ip + 1
                     Kpts(:,ip) = rcell(:,1)*real(i,dp)/real(nk_3D_x,dp) + &
                                 rcell(:,2)*real(j,dp)/real(nk_3D_y,dp) + &
                                 rcell(:,3)*0.0_dp
                  end do
               end do
            else
               do j = 1, nk_3D_y
                  do i = 1, nk_3D_x
                     ip = ip + 1
                     ! Use exact same MP formula as Diag3DBands (lines 443-444) with global rcell
                     Kpts(:,ip) = rcell(:,1)*(2*i-nk_3D_x-1)/(2.0_dp*nk_3D_x) + &
                                 rcell(:,2)*(2*j-nk_3D_y-1)/(2.0_dp*nk_3D_y) + &
                                 rcell(:,3)*0.0_dp  ! 2D system: kz = 0
                  end do
               end do
            end if
         end if

         ! Set the extra k-point to avoid uninitialized memory (required by OpenMP loop)
         Kpts(:, ptsTot+1) = Kpts(:, ptsTot)  ! Duplicate last k-point

         if (calculateChern .and. useTAPW) then
            call MIO_Print('Generated '//trim(num2str(ptsTot))//' k-points for Chern calculation','diag')
         else if (calculate3DTAPWBands) then
            call MIO_Print('Generated '//trim(num2str(ptsTot))//' k-points for 3D TAPW bands calculation','diag')
         end if
      end if

      ! Allocate ELoc with appropriate size for TAPW or regular calculations
      if (useTAPW) then
         ! For TAPW: allocate large enough for maximum possible M = NG * Nlabel
         ! Use a conservative estimate based on typical N_G values
         allocate(ELoc(max(nAt, 15000)))  ! Allow up to 15k TAPW states
         ELoc = 0.0_dp  ! Initialize to avoid uninitialized memory issues
         call MIO_Print('Allocated ELoc for TAPW with size: '//trim(num2str(size(ELoc))), 'diag')

      ! PERFORMANCE OPTIMIZATION: Pre-allocate TAPW arrays for reuse across k-points
      ! This avoids expensive allocate/deallocate for each of 20,736 k-points
      if (calculateChern) then
         call MIO_Print('Pre-allocating TAPW arrays for high-performance Chern calculation...', 'diag')
         ! Estimate M based on typical values (will be updated in first call)

         ! Pre-allocate all arrays that are currently allocated per k-point
         if (.not. allocated(tapw_H_dense)) then
            allocate(tapw_H_dense(nAt, nAt))
            call MIO_Print('Pre-allocated H_dense: '//trim(num2str(nAt))//'x'//trim(num2str(nAt))//' = '//trim(num2str(nAt*nAt*8/1024/1024))//' MB', 'diag')
         end if

         if (.not. allocated(tapw_Hproj)) then
            allocate(tapw_Hproj(200, 200))
            call MIO_Print('Pre-allocated Hproj: 200x200 = 0.3 MB', 'diag')
         end if

         if (.not. allocated(tapw_eigvals)) then
            allocate(tapw_eigvals(200))
            call MIO_Print('Pre-allocated eigvals: 200 elements', 'diag')
         end if

         if (.not. allocated(tapw_ZWorkLoc)) then
            allocate(tapw_ZWorkLoc(400))
            call MIO_Print('Pre-allocated ZWorkLoc: 400 elements', 'diag')
         end if

         if (.not. allocated(tapw_DWorkLoc)) then
            allocate(tapw_DWorkLoc(600))
            call MIO_Print('Pre-allocated DWorkLoc: 600 elements', 'diag')
         end if

         call MIO_Print('TAPW memory pool ready for '//trim(num2str(ptsTot))//' k-points', 'diag')
      end if
      else
         ! For regular calculations: use number of atoms
         allocate(ELoc(nAt))
         ELoc = 0.0_dp  ! Initialize
      end if

      if (keepWaveFunction .or. .not. (sparseDiagSolver .or. useTAPW)) allocate(HLoc(nAt,nAt))
      d = 0.0_dp
      hv = -huge(0.0_dp)
      lc = huge(0.0_dp)
      E = 0.0_dp
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !      !if (modulo(ip,int(ptsTot/10)).eq.0) print*, "progress is: ", ip/int(ptsTot/10)*10, "percent"
      !             !write(uu,*) (HLoc(j,i),j=1,nAt)
      !         !call DiagHamSparse(nAt,nspin,is,HLoc,ELoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
      !         !call DiagHamSparse2(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)

      ! Print initial TAPW information
      if (useTAPW) then
         call MIO_Print('Starting TAPW calculations for '//trim(num2str(ptsTot))//' k-points...', 'diag')
      end if

      if (ptsTot == 1) then
         ! Special case: only one k-point, do not use OpenMP here
         do ip = 1, ptsTot
            ! Progress tracking for TAPW
            if (useTAPW) then
               call MIO_Print('TAPW: Processing k-point '//trim(num2str(ip))//'/'//trim(num2str(ptsTot)), 'diag')
            end if

            do is = 1, nspin
               ! Don't zero entire ELoc array - causes issues with large TAPW arrays
               ! ELoc will be properly set by the diagonalization routines
               if (allocated(HLoc) .and. .not. sparseDiagSolver) then
                  HLoc = 0.0_dp
               end if
               KptsLoc = Kpts(:, ip)
               if (keepWaveFunction) then
                  call DiagHamWF(nAt, nspin, is, HLoc, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
                  do i = 1, nAt
                     if (ip == 1) then
                        write(uu, '(20000(F10.5,1X))') (real(HLoc(j, i)), j = 1, nAt)
                        write(uuuu, '(20000(F10.5,1X))') (real(conjg(HLoc(j, i)) * HLoc(j, i)), j = 1, nAt)
                        write(uuu, *) ELoc(i)
                     end if
                  end do
               else if (sparseDiagSolver) then
                  call DiagHamSparse(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig)
               else if (useTAPW .and. RashbaSOCterm .and. is == 1) then
                  ! TAPW + Rashba SOC: Build block Hamiltonian once, then use TAPW diagonalization
                  ! Only Rashba requires block Hamiltonian due to spin-flip terms
                  !$OMP CRITICAL
                  if (.not. allocated(HBlock)) allocate(HBlock(2*nAt, 2*nAt))
                  call BuildBlockHamiltonianOnly(nAt, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, HBlock)
                  call DiagH0TAPW_withBlockH(nAt, nspin, is, ELoc, KptsLoc, ucell, HBlock, maxNeigh, hopp, NList, Nneigh, neighCell, neig, ip)
                  !$OMP END CRITICAL
               else if (useTAPW) then
                  call DiagH0TAPW(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig, ip)
               else if (RashbaSOCterm .and. is == 1) then
                  ! Use block Hamiltonian for Rashba (only call once, gives 2N eigenvalues)
                  !$OMP CRITICAL
                  if (.not. allocated(HBlock)) allocate(HBlock(2*nAt, 2*nAt))
                  if (.not. allocated(ELocBlock)) allocate(ELocBlock(2*nAt))
                  call BuildBlockHamiltonianOnly(nAt, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, HBlock)
                  call DiagBlockHamiltonian(HBlock, ELocBlock, nAt)
                  !$OMP END CRITICAL
                  ! Store eigenvalues - first nAt go to spin-up, next nAt to spin-down
                  E(1:nAt, 1, ip) = ELocBlock(1:nAt)
                  E(1:nAt, 2, ip) = ELocBlock(nAt+1:2*nAt)
                  ! Skip the "Copy appropriate number" section for Rashba since we've already set E
                  cycle
               else if (RashbaSOCterm .and. is == 2) then
                  ! Skip spin-down iteration for Rashba since block Hamiltonian already handled both spins
                  cycle
               else if (.not. RashbaSOCterm) then
                  call DiagHam(nAt, nspin, is, HLoc, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
               end if
               ! Copy appropriate number of eigenvalues based on calculation type
               if (useTAPW) then
                  ! For TAPW: copy M_tapw eigenvalues (set by DiagH0TAPW)
                  if (M_tapw > 0) then
                     E(1:M_tapw, is, ip) = ELoc(1:M_tapw)
                  else
                     ! Fallback if M_tapw not set properly
                     E(1:nAt, is, ip) = ELoc(1:nAt)
                  end if
               else
                  ! For regular calculations: copy nAt eigenvalues
                  E(1:nAt, is, ip) = ELoc(1:nAt)
               end if
            end do
         end do
         do ip=1,ptsTot
            if (ip.eq.1) then
               if (useTAPW) then
                  write(u,'(3i8)') M_tapw, nspin, ptsTot+1  ! TAPW: M eigenvalues
               else
                  write(u,'(3i8)') nAt, nspin, ptsTot+1     ! Regular: nAt eigenvalues
               end if
            end if
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            d = d + sqrt(dot_product(v,v))
            if (useTAPW) then
               ! For TAPW: only output the meaningful eigenvalues (M = NG * Nlabel)
               ! Use M = 80 for now (20 * 4), could be made dynamic
               write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is,ip)*g0, i=1,M_tapw),is=1,nspin)
            else
               ! For other methods: output all nAt eigenvalues
               write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
            end if
            if (sum(nPts(:nPath))==ip-1) then
               nPath = nPath+1
               call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
                 '   '//trim(num2str(d,6)),'diag')
            end if
         end do
      else if (sparseDiagSolver) then
         ! Many k-points, but parallelisation happens inside DiagHamSparse so no OMP here (some obsolete code here, too lazy to
         ! remove the parts that will not be considered because we are only doing the sparseDiagSolver part here
         do ip = 1, ptsTot
            ! Progress tracking for TAPW (show every 10%)
            if (useTAPW .and. modulo(ip, max(1, int(ptsTot/10))) == 0 .and. ip <= ptsTot) then
               call MIO_Print('TAPW progress: '//trim(num2str(nint(100.0_dp*ip/ptsTot)))//'% ('// &
                             trim(num2str(ip))//'/'//trim(num2str(ptsTot))//' k-points)', 'diag')
            end if

            do is = 1, nspin
               ! Don't zero entire ELoc array - causes issues with large TAPW arrays
               ! ELoc will be properly set by the diagonalization routines
               if (allocated(HLoc) .and. .not. sparseDiagSolver) then
                  HLoc = 0.0_dp
               end if
               KptsLoc = Kpts(:, ip)
               if (keepWaveFunction) then
                  call DiagHamWF(nAt, nspin, is, HLoc, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
                  do i = 1, nAt
                     if (ip == 1) then
                        write(uu, '(20000(F10.5,1X))') (real(HLoc(j, i)), j = 1, nAt)
                        write(uuuu, '(20000(F10.5,1X))') (real(conjg(HLoc(j, i)) * HLoc(j, i)), j = 1, nAt)
                        write(uuu, *) ELoc(i)
                     end if
                  end do
               else if (sparseDiagSolver) then
                  call DiagHamSparse(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig)
               else if (useTAPW .and. RashbaSOCterm .and. is == 1) then
                  ! TAPW + Rashba SOC: Build block Hamiltonian once, then use TAPW diagonalization
                  !$OMP CRITICAL
                  if (.not. allocated(HBlock)) allocate(HBlock(2*nAt, 2*nAt))
                  call BuildBlockHamiltonianOnly(nAt, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, HBlock)
                  call DiagH0TAPW_withBlockH(nAt, nspin, is, ELoc, KptsLoc, ucell, HBlock, maxNeigh, hopp, NList, Nneigh, neighCell, neig, ip)
                  !$OMP END CRITICAL
               else if (useTAPW) then
                  call MIO_Print('calling TAPW routine for '//trim(num2str(ip)),'diag')
                  call DiagH0TAPW(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig, ip)
               else if (RashbaSOCterm .and. is == 1) then
                  ! Use block Hamiltonian for Rashba (only call once, gives 2N eigenvalues)
                  !$OMP CRITICAL
                  if (.not. allocated(HBlock)) allocate(HBlock(2*nAt, 2*nAt))
                  if (.not. allocated(ELocBlock)) allocate(ELocBlock(2*nAt))
                  call BuildBlockHamiltonianOnly(nAt, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, HBlock)
                  call DiagBlockHamiltonian(HBlock, ELocBlock, nAt)
                  !$OMP END CRITICAL
                  ! Store eigenvalues - first nAt go to spin-up, next nAt to spin-down
                  E(1:nAt, 1, ip) = ELocBlock(1:nAt)
                  E(1:nAt, 2, ip) = ELocBlock(nAt+1:2*nAt)
                  ! Skip the "Copy appropriate number" section for Rashba since we've already set E
                  cycle
               else if (RashbaSOCterm .and. is == 2) then
                  ! Skip spin-down iteration for Rashba since block Hamiltonian already handled both spins
                  cycle
               else if (.not. RashbaSOCterm) then
                  call DiagHam(nAt, nspin, is, HLoc, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
               end if
               ! Copy appropriate number of eigenvalues based on calculation type
               if (useTAPW) then
                  ! For TAPW: copy M_tapw eigenvalues (set by DiagH0TAPW)
                  if (M_tapw > 0) then
                     E(1:M_tapw, is, ip) = ELoc(1:M_tapw)
                  else
                     ! Fallback if M_tapw not set properly
                     E(1:nAt, is, ip) = ELoc(1:nAt)
                  end if
               else
                  ! For regular calculations: copy nAt eigenvalues
                  E(1:nAt, is, ip) = ELoc(1:nAt)
               end if
            end do
         end do

         ! Write header and loop for single k-point case
         if (calculate3DTAPWBands) then
            ! 3D TAPW bands: write header and loop to ptsTot
            if (useTAPW) then
               write(u,'(3i8)') M_tapw, nspin, ptsTot  ! TAPW: M eigenvalues
            else
               write(u,'(3i8)') nAt, nspin, ptsTot     ! Regular: nAt eigenvalues
            end if

            do ip=1,ptsTot
               ! For 3D TAPW bands: write kx, ky coordinates instead of cumulative distance
               if (useTAPW) then
                  write(u,'(3f12.6,10f14.6,/,(10x,10f14.6))') Kpts(1,ip), Kpts(2,ip), Kpts(3,ip), ((E(i,is,ip)*g0, i=1,M_tapw),is=1,nspin)
               else
                  write(u,'(3f12.6,10f14.6,/,(10x,10f14.6))') Kpts(1,ip), Kpts(2,ip), Kpts(3,ip), ((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
               end if
            end do
         else
            ! Regular bands: write header and loop to ptsTot+1
            if (useTAPW) then
               write(u,'(3i8)') M_tapw, nspin, ptsTot+1  ! TAPW: M eigenvalues
            else
               write(u,'(3i8)') nAt, nspin, ptsTot+1     ! Regular: nAt eigenvalues
            end if

            do ip=1,ptsTot+1
               ! For regular band structure: use cumulative distance
               v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
               d = d + sqrt(dot_product(v,v))
               if (useTAPW) then
                  ! For TAPW: only output the meaningful eigenvalues (M = NG * Nlabel)
                  write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is,ip)*g0, i=1,M_tapw),is=1,nspin)
               else
                  ! For other methods: output all nAt eigenvalues
                  write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
               end if
               if (.not. calculate3DTAPWBands .and. sum(nPts(:nPath))==ip-1) then
                  nPath = nPath+1
                  call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
                    '   '//trim(num2str(d,6)),'diag')
               end if
            end do
         end if
      else
         ! General case: more than one k-point, use OpenMP
#ifdef _OPENMP
         ! TAPW must run with one OpenMP thread: the k-loop below is an OpenMP
         ! loop, and the TAPW routines called from it read input parameters and
         ! call threaded LAPACK, neither of which is safe inside it. The threads
         ! of the linear-algebra library do the parallel work instead.
         if (useTAPW) then
            block
               integer, external :: omp_get_max_threads
               if (omp_get_max_threads() > 1) then
                  call MIO_Kill('TAPW calculations must be run with one OpenMP thread. Set OMP_NUM_THREADS=1 '// &
                    'and give the cores to the linear-algebra library instead (MKL_NUM_THREADS or '// &
                    'OPENBLAS_NUM_THREADS = number of cores), then run again.','diag','DiagBands')
               end if
            end block
         end if
#endif

         ! Cache TAPW-related inputs once outside the OpenMP region to avoid nested timer/input
         if (useTAPW .and. .not. tapw_cfg_initialized) then
            call MIO_InputParameter('TAPW.aG',tapw_cfg_aG,2.46019_dp)
            call MIO_InputParameter('Diag.SparseSetShift',tapw_cfg_shift,0.0_dp)
            call MIO_InputParameter('Diag.SparseUseShift',tapw_cfg_useShift,.false.)
            call MIO_InputParameter('Diag.SparseSaveRitz',tapw_cfg_saveRitz,.false.)
            call MIO_InputParameter('Diag.SparseSetTol',tapw_cfg_tol,0.1_dp)
            tapw_cfg_initialized = .true.
         end if

#ifdef SEMICL
         if (berryFluxTAPW .or. orbMomentTAPW) &
            call TAPWGeomInit(ptsTot, nAt, nspin, Kpts, H0, maxNeigh, &
                              hopp, NList, Nneigh, neighCell, neig, ucell)
#endif

         !$OMP PARALLEL DO PRIVATE(ELoc, KptsLoc, is, ip, HLoc, i, j, HBlock, ELocBlock, i_eig, useTAPW_str, forceBlockTAPW_str), &
         !$OMP& SHARED(E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, Kpts, neig, &
         !$OMP& sparseDiagSolver, keepWaveFunction, useTAPW, uu, uuu, uuuu, ptsTot, RashbaSOCterm, forceBlockTAPW, socDebug, M_tapw)
         do ip = 1, ptsTot + 1
            ! Progress tracking for TAPW (show every 10%, but only from thread 1 to avoid spam)
            !$OMP CRITICAL
            if (useTAPW .and. modulo(ip, max(1, int(ptsTot/10))) == 0 .and. ip <= ptsTot) then
               call MIO_Print('TAPW progress: '//trim(num2str(nint(100.0_dp*ip/ptsTot)))//'% ('// &
                             trim(num2str(ip))//'/'//trim(num2str(ptsTot))//' k-points)', 'diag')
            end if
            !$OMP END CRITICAL

            do is = 1, nspin
               ! Don't zero entire ELoc array in parallel - causes race conditions with TAPW
               ! ELoc will be properly set by the diagonalization routines
               if (allocated(HLoc) .and. .not. sparseDiagSolver) then
                  HLoc = 0.0_dp
               end if
               KptsLoc = Kpts(:, ip)
               if (keepWaveFunction) then
                  call DiagHamWF(nAt, nspin, is, HLoc, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
                  do i = 1, nAt
                     if (ip == 1) then
                        write(uu, '(20000(F10.5,1X))') (real(HLoc(j, i)), j = 1, nAt)
                        write(uuuu, '(20000(F10.5,1X))') (real(conjg(HLoc(j, i)) * HLoc(j, i)), j = 1, nAt)
                        write(uuu, *) ELoc(i)
                     end if
                  end do
               else if (sparseDiagSolver) then
                  call DiagHamSparse(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig)
               else if (useTAPW .and. (RashbaSOCterm .or. forceBlockTAPW) .and. is == 1) then
                  ! Debug: Verify conditions are met
                  if (forceBlockTAPW .and. ip <= 2) then
                     if (useTAPW) then
                        useTAPW_str = 'true '
                     else
                        useTAPW_str = 'false'
                     end if
                     if (forceBlockTAPW) then
                        forceBlockTAPW_str = 'true '
                     else
                        forceBlockTAPW_str = 'false'
                     end if
                     call MIO_Print('DEBUG: forceBlockTAPW flag check - useTAPW='//trim(useTAPW_str)//', forceBlockTAPW='//trim(forceBlockTAPW_str)//', is='//trim(num2str(is))//', ip='//trim(num2str(ip)), 'diag')
                  end if
                  ! TAPW + Rashba SOC (or forced block path for testing): Build block Hamiltonian once, then use TAPW diagonalization
                  ! Only Rashba requires block Hamiltonian due to spin-flip terms, but forceBlockTAPW allows testing without SOC
                  !$OMP CRITICAL
                  if (forceBlockTAPW) then
                     call MIO_Print('TAPW: Using forced block path (forceBlockTAPW=true, no SOC)', 'diag')
                     call MIO_Print('About to call BuildBlockHamiltonianOnly for TAPW (forced block path, no SOC)', 'diag')
                  else if (socDebug) then
                     call MIO_Print('About to call BuildBlockHamiltonianOnly for TAPW+Rashba', 'diag')
                  end if
                  if (.not. allocated(HBlock)) allocate(HBlock(2*nAt, 2*nAt))
                  call BuildBlockHamiltonianOnly(nAt, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, HBlock)
                  if (forceBlockTAPW) then
                     call MIO_Print('BuildBlockHamiltonianOnly completed, calling DiagH0TAPW_withBlockH', 'diag')
                  else if (socDebug) then
                     call MIO_Print('BuildBlockHamiltonianOnly completed, calling DiagH0TAPW_withBlockH', 'diag')
                  end if
                  call DiagH0TAPW_withBlockH(nAt, nspin, is, ELoc, KptsLoc, ucell, HBlock, maxNeigh, hopp, NList, Nneigh, neighCell, neig, ip)
                  if (forceBlockTAPW) then
                     call MIO_Print('DiagH0TAPW_withBlockH completed', 'diag')
                  else if (socDebug) then
                     call MIO_Print('DiagH0TAPW_withBlockH completed', 'diag')
                  end if
                  !$OMP END CRITICAL
                  ! For Rashba SOC with TAPW, we now have 2M_tapw eigenvalues (spin-mixed from block projection)
                  ! For forceBlockTAPW (no SOC), we have 2M_tapw eigenvalues in degenerate pairs (λ₁,λ₁), (λ₂,λ₂), ...
                  ! Store eigenvalues for both spins (they're spin-mixed, but we need to populate E array)
                  ! ELoc contains min(2*M_tapw, nAt) eigenvalues, so we need to be careful with bounds
                  if (M_tapw > 0) then
                     if (forceBlockTAPW) then
                        ! For no-SOC case: eigenvalues are in degenerate pairs, extract only M unique ones
                        ! Pattern: (λ₁,λ₁), (λ₂,λ₂), ..., so take eigvals(2i-1) for i=1..M
                        ! Extract one eigenvalue from each degenerate pair
                        do i_eig = 1, min(M_tapw, size(E,1))
                           if (2*i_eig - 1 <= size(ELoc)) then
                              E(i_eig, 1, ip) = ELoc(2*i_eig - 1)  ! Take first of each pair
                              E(i_eig, 2, ip) = ELoc(2*i_eig - 1)  ! Same for both spins (degenerate)
                           else
                              E(i_eig, 1, ip) = 0.0_dp
                              E(i_eig, 2, ip) = 0.0_dp
                           end if
                        end do
                        ! Zero out remaining entries if M_tapw < size(E,1)
                        if (M_tapw < size(E,1)) then
                           E(M_tapw+1:size(E,1), 1, ip) = 0.0_dp
                           E(M_tapw+1:size(E,1), 2, ip) = 0.0_dp
                        end if
                     else
                        ! For Rashba SOC: use all 2M_tapw eigenvalues (spin-mixed)
                        ! Use min(2*M_tapw, nAt) since ELoc now contains 2M_tapw eigenvalues from block projection
                        E(1:min(2*M_tapw, nAt), 1, ip) = ELoc(1:min(2*M_tapw, nAt))
                        ! For Rashba, spin 2 uses the same eigenvalues (they're spin-mixed)
                        E(1:min(2*M_tapw, nAt), 2, ip) = ELoc(1:min(2*M_tapw, nAt))
                        ! If 2*M_tapw > nAt, zero out the remaining entries
                        if (2*M_tapw > nAt) then
                           E(nAt+1:min(2*M_tapw, size(E,1)), 1, ip) = 0.0_dp
                           E(nAt+1:min(2*M_tapw, size(E,1)), 2, ip) = 0.0_dp
                        end if
                     end if
                  else
                     ! Fallback if M_tapw not set
                     E(1:nAt, 1, ip) = ELoc(1:nAt)
                     E(1:nAt, 2, ip) = ELoc(1:nAt)
                  end if
                  ! Skip the regular eigenvalue copying section
                  cycle
               else if (useTAPW .and. (RashbaSOCterm .or. forceBlockTAPW) .and. is == 2) then
                  ! Skip spin-down iteration for TAPW+Rashba (or forced block path) since block Hamiltonian already handled both spins
                  cycle
               else if (useTAPW) then
#ifdef SEMICL
                  if (allocated(bf_evec) .and. is == 1 .and. ip <= ptsTot) then
                     call DiagH0TAPW(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig, ip, &
                                     bf_evec(:,:,ip))
                  else
#endif
                     call DiagH0TAPW(nAt, nspin, is, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell,neig, ip)
#ifdef SEMICL
                  end if
#endif
               else if (RashbaSOCterm .and. is == 1) then
                  ! Use block Hamiltonian for Rashba (only call once, gives 2N eigenvalues)
                  !$OMP CRITICAL
                  if (socDebug) call MIO_Print('About to call BuildBlockHamiltonianOnly for Rashba', 'diag')
                  if (.not. allocated(HBlock)) allocate(HBlock(2*nAt, 2*nAt))
                  if (.not. allocated(ELocBlock)) allocate(ELocBlock(2*nAt))
                  call BuildBlockHamiltonianOnly(nAt, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, HBlock)
                  if (socDebug) call MIO_Print('BuildBlockHamiltonianOnly completed, calling DiagBlockHamiltonian', 'diag')
                  call DiagBlockHamiltonian(HBlock, ELocBlock, nAt)
                  if (socDebug) call MIO_Print('DiagBlockHamiltonian completed, storing eigenvalues', 'diag')
                  !$OMP END CRITICAL
                  ! Store eigenvalues - first nAt go to spin-up, next nAt to spin-down
                  E(1:nAt, 1, ip) = ELocBlock(1:nAt)
                  E(1:nAt, 2, ip) = ELocBlock(nAt+1:2*nAt)
                  ! Skip the "Copy appropriate number" section for Rashba since we've already set E
                  cycle
               else if (RashbaSOCterm .and. is == 2) then
                  ! Skip spin-down iteration for Rashba since block Hamiltonian already handled both spins
                  cycle
               else if (.not. RashbaSOCterm) then
                  call DiagHam(nAt, nspin, is, HLoc, ELoc, KptsLoc, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
               end if
               ! Copy appropriate number of eigenvalues based on calculation type
               if (useTAPW) then
                  ! For TAPW: copy M_tapw eigenvalues (set by DiagH0TAPW)
                  if (M_tapw > 0) then
                     E(1:M_tapw, is, ip) = ELoc(1:M_tapw)
                  else
                     ! Fallback if M_tapw not set properly
                     E(1:nAt, is, ip) = ELoc(1:nAt)
                  end if
               else
                  ! For regular calculations: copy nAt eigenvalues
                  E(1:nAt, is, ip) = ELoc(1:nAt)
               end if
            end do
         end do
         !$OMP END PARALLEL DO

         ! Print completion message for TAPW
         if (useTAPW) then
            call MIO_Print('TAPW calculations completed for all k-points!', 'diag')
         end if

#ifdef SEMICL
         if (allocated(bf_evec)) call BerryFluxWrite(nk_3D_x, nk_3D_y, ptsTot)
         if (allocated(om_M)) call OrbMomentWrite(nk_3D_x, nk_3D_y, ptsTot, E)
         if (allocated(os_w)) call OrbSuscWrite(nk_3D_x, nk_3D_y, ptsTot, E)
#endif

         ! Skip band structure output for Chern calculation (not needed and causes array access issues)
         if (.not. (calculateChern .and. useTAPW)) then
            ! For 3D grid calculations, only loop to ptsTot (no extra k-point)
            ! For regular path calculations, loop to ptsTot+1 (includes end point)
            if (calculate3DTAPWBands) then
               do ip=1,ptsTot
                  if (ip == 1) then
                     if (useTAPW) then
                        write(u,'(3i8)') M_tapw, nspin, ptsTot  ! TAPW: M eigenvalues
                     else
                        write(u,'(3i8)') nAt, nspin, ptsTot     ! Regular: nAt eigenvalues
                     end if
                  end if

                  ! For 3D TAPW bands: write kx, ky coordinates instead of cumulative distance
                  if (useTAPW) then
                     write(u,'(3f12.6,10f14.6,/,(10x,10f14.6))') Kpts(1,ip), Kpts(2,ip), Kpts(3,ip), ((E(i,is,ip)*g0, i=1,M_tapw),is=1,nspin)
                  else
                     write(u,'(3f12.6,10f14.6,/,(10x,10f14.6))') Kpts(1,ip), Kpts(2,ip), Kpts(3,ip), ((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
                  end if
               end do
            else
               do ip=1,ptsTot + 1
                  if (ip == 1) then
                     if (useTAPW) then
                        write(u,'(3i8)') M_tapw, nspin, ptsTot+1  ! TAPW: M eigenvalues
                     else
                        write(u,'(3i8)') nAt, nspin, ptsTot+1     ! Regular: nAt eigenvalues
                     end if
                  end if

                  ! For regular band structure: use cumulative distance
                  v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
                  d = d + sqrt(dot_product(v,v))
                  if (useTAPW) then
                     write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is,ip)*g0, i=1,M_tapw),is=1,nspin)
                  else
                     write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is,ip)*g0, i=1,nAt),is=1,nspin)
                  end if
                  if (.not. calculate3DTAPWBands .and. sum(nPts(:nPath))==ip-1) then
                     nPath = nPath+1
                     call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
                       '   '//trim(num2str(d,6)),'diag')
                  end if
               end do
            end if
         end if
      end if

      ! Chern number calculation for TAPW
      if (calculateChern .and. useTAPW) then
         call MIO_Print('')
         call MIO_Print('Starting Chern number calculation...','diag')
         call MIO_Print('Note: For Chern calculation, consider setting OMP_NUM_THREADS=1 to avoid MKL/OpenMP conflicts','diag')
         call CalculateChernTAPW(Kpts, ptsTot, E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, RashbaSOCterm)
      end if

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Print('')
      close(u)
   end if

   ! Deallocate ELoc
   if (allocated(ELoc)) then
      deallocate(ELoc)
   end if

   ! Always ensure timer is stopped, regardless of exit path
#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagBands',1)
#endif /* DEBUG */

end subroutine DiagBands

subroutine DiagChern()

   use cell,                 only : rcell, ucell, aG
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, uu, uuu, uuuu
   real(dp), pointer :: path(:,:)=>NULL(), Kgrid(:,:)=>NULL(), Chern(:,:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   integer :: nk(3), i1, i2, i3, ik
   real(dp) :: d0, v(3), d, hv, lc, totalChern
   complex(dp), pointer :: H(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS, useSameNumberOfPoints
   real(dp) :: theta, volume, dArea

   real(dp) :: gcell(3,3), grcell(3,3)

   real(dp) :: KptsLoc(3)
   real(dp) :: ChernLoc(nAt)
   complex(dp) :: HLoc(nAt, nAt)
   logical :: useDifferentLatticeVectors, keepWaveFunction

#ifdef DEBUG
   call MIO_Debug('DiagChern',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_InputParameter('Bands.NumPoints',nPts0,100)
   call MIO_InputParameter('Bands.UseDifferentLatticeVectors',useDifferentLatticeVectors,.false.)
   call MIO_InputParameter('Bands.UseSameNumberOfPoints',useSameNumberOfPoints,.false.)
   if (useDifferentLatticeVectors) then
       ucell(1,2) = 0.0_dp
       rcell(2,1) = 0.0_dp
   end if
      !   !print*, "0: ", path(:,ip)
      !   !print*, "1: ", path(:,ip)
      !   !if (MoireBS) then
      !   !   !print*, "theta=", theta
      !   !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
      !   !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
      !   !end if
      !   !print*, "2: ", path(:,ip)
      !   !call MIO_Allocate(nPts,nPath-1,'nPts','diag')
      !         !ptsTot = ptsTot + nPts(ip)
      !!ip = 0
      !      !Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)

      call MIO_Print('Calculating Chern number by diagonalization','diag')
      call MIO_InputParameter('KGrid',nk,[1,1,1])
      ptsTot = nk(1)*nk(2)*nk(3)
      call MIO_Allocate(Kgrid,[3,ptsTot],'Kgrid','diag')
      ik = 0
      do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
         ik = ik+1
         Kgrid(:,ik) = rcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
           rcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + rcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
      end do; end do; end do

      call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
      call MIO_Allocate(Chern,[nAt,nspin,ptsTot+1],'Chern','diag')
      flnm = trim(prefix)//'.chern'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot
      call MIO_InputParameter('keepWaveFunction',keepWaveFunction,.false.)
      d = 0.0_dp
      hv = -huge(0.0_dp)
      lc = huge(0.0_dp)
      Chern = 0.0_dp
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !$OMP PARALLEL DO PRIVATE(ChernLoc, KptsLoc, is, ip, HLoc), &
      !$OMP& SHARED(Chern, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, Kgrid)
      do ip=1,ptsTot
         do is=1,nspin
            !if (modulo(ip,int(ptsTot/10)).eq.0) print*, "progress is: ", ip/int(ptsTot/10)*10, "percent"
            ChernLoc = 0.0_dp
            KptsLoc = Kgrid(:,ip)
            !       !write(uu,*) (HLoc(j,i),j=1,nAt)
               call DiagHamChern(nAt,nspin,is,HLoc,ChernLoc,KptsLoc,ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            Chern(:,is,ip) = ChernLoc
         end do
      end do
      !$OMP END PARALLEL DO
      do ip=1,ptsTot
         !!!v = Kpts(:,ip+1) - Kpts(:,max(ip,1))
         write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') Kgrid(:,ip),((Chern(i,is,ip), i=1,nAt),is=1,nspin)
      end do
      dArea = norm(CrossProd(rcell(:,1),rcell(:,2)))/(nk(1)*nk(2))
      do i=1,nAt
         totalChern = 0.0_dp
         do ip=1,ptsTot
            totalChern = totalChern+Chern(i,1,ip)*dArea
         end do
         call MIO_Print('Chern number for band '//trim(num2str(i))//' equals '//trim(num2str(totalChern,6)),'diag')
      end do
      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(Chern,'Chern','diag')
      call MIO_Deallocate(H,'H','diag')
      call MIO_Print('')
      close(u)

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagChern',1)
#endif /* DEBUG */

end subroutine DiagChern

subroutine DiagBandsRashba()

   use cell,                 only : rcell, ucell, aG
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: H(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm
   real(dp) :: rcell_rotated(3,3), rotation_angle, cos_rot, sin_rot, temp_vec(3)

   logical :: MoireBS
   real(dp) :: theta, volume

   real(dp) :: gcell(3,3), grcell(3,3)

   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt*2)
   complex(dp) :: HLoc(nAt*2, nAt*2)
   logical :: useDifferentLatticeVectors

#ifdef DEBUG
   call MIO_Debug('DiagBandsRashba',0)
#endif /* DEBUG */
#ifdef TIMER
#endif /* TIMER */

   call MIO_InputParameter('Bands.NumPoints',nPts0,100)
   call MIO_InputParameter('Bands.UseDifferentLatticeVectors',useDifferentLatticeVectors,.false.)
   if (useDifferentLatticeVectors) then
       ucell(1,2) = 0.0_dp
       rcell(2,1) = 0.0_dp
   end if
   if (MIO_InputFindBlock('Bands.Path',nPath)) then
      call MIO_Print('Band calculation','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_InputBlock('Bands.Path',path)
      do ip=1,nPath
         path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
         !   !print*, "theta=", theta
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1 ! coz 4 points, 3 segments
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
            nPts(ip) = nPts(ip) + 1 ! Add the missing point at the end of last segment
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot+1],'Kpts','diag') ! add missing point
      Kpts(:,1) = path(:,1)
      ip = 1
      d = 0.0_dp
      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            Kpts(:,ip) = path(:,i) + (j)*(path(:,i+1)-path(:,i))/nPts(i)
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            d = d + sqrt(dot_product(v,v))
         end do
      end do
      call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
      call MIO_Allocate(E,[nAt*2,nspin,ptsTot+1],'E','diag')
      flnm = trim(prefix)//'.bands'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt*2, nspin, ptsTot+1
      d = 0.0_dp
      hv = -huge(0.0_dp)
      lc = huge(0.0_dp)
      E = 0.0_dp
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      !$OMP PARALLEL DO PRIVATE(ELoc, KptsLoc, is, ip, HLoc), &
      !$OMP& SHARED(E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, Kpts)
      do ip=1,ptsTot + 1
         do is=1,nspin
            !if (modulo(ip,int(ptsTot/10)).eq.0) print*, "progress is: ", ip/int(ptsTot/10)*10, "percent"
            ELoc = 0.0_dp
            HLoc = 0.0_dp
            KptsLoc = Kpts(:,ip)
            ! DiagHamRashba removed - use BuildBlockHamiltonianOnly + DiagBlockHamiltonian instead
            call MIO_Kill('DiagBandsRashba: This routine needs to be updated to use BuildBlockHamiltonianOnly + DiagBlockHamiltonian', 'diag', 'DiagBandsRashba')
            E(:,is,ip) = ELoc
         end do
      end do
      !$OMP END PARALLEL DO
      do ip=1,ptsTot + 1
         v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
         d = d + sqrt(dot_product(v,v))
         write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is,ip)*g0, i=1,nAt*2),is=1,nspin)
         if (sum(nPts(:nPath))==ip-1) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
      end do
      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(H,'H','diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagBandsRashba',1)
#endif /* DEBUG */

end subroutine DiagBandsRashba

subroutine DiagBandsAroundK()

   use cell,                 only : rcell, ucell, aG
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: H(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS
   real(dp) :: theta, volume, vn(3), grapheneK(3)

   real(dp) :: gcell(3,3), grcell(3,3)

#ifdef DEBUG
   call MIO_Debug('DiagBandsAroundK',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   gcell(:,1) = (/aG,0.0_dp,0.0_dp/)
   gcell(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp/)
   gcell(:,3) = (/0.0_dp,0.0_dp,40.0_dp/)

   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)

   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   call MIO_InputParameter('Bands.NumPoints',nPts0,100)
   if (MIO_InputFindBlock('Bands.Path',nPath)) then
      call MIO_Print('Band calculation','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_InputBlock('Bands.Path',path)
      do ip=1,nPath
         path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
         grapheneK = 2.0/3.0*grcell(:,1) + 1.0/3.0*grcell(:,2) + 0.0*grcell(:,3)
         path(:,ip) = path(:,ip) + grapheneK
         !   !print*, "theta=", theta
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1 ! coz 4 points, 3 segments
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
            nPts(ip) = nPts(ip) + 1 ! Add the missing point at the end of last segment
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot+1],'Kpts','diag') ! add missing point
      Kpts(:,1) = path(:,1)
      ip = 1
      d = 0.0_dp
      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            Kpts(:,ip) = path(:,i) + (j)*(path(:,i+1)-path(:,i))/nPts(i)
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            d = d + sqrt(dot_product(v,v))
         end do
      end do
      call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.bands'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot+1
      d = 0.0_dp
      hv = -huge(0.0_dp)
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      do ip=1,ptsTot + 1
         do is=1,nspin
            call DiagHam(nAt,nspin,is,H(:,:,is),E(:,is),Kpts(:,ip),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !!end do
            hv = max(hv,maxval(E(:,is),mask=E(:,is)<=Efermi/g0))
            lc = min(lc,minval(E(:,is),mask=E(:,is)>Efermi/g0))
         end do
         v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
         d = d + sqrt(dot_product(v,v))
         write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)
         if (sum(nPts(:nPath))==ip-1) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
      end do
      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(H,'H','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagBandsAroundK',1)
#endif /* DEBUG */

end subroutine DiagBandsAroundK

subroutine DiagBandsG()

   use cell,                 only : rcell, ucell, aG
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: H(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS
   real(dp) :: theta, volume

   real(dp) :: gcell(3,3), grcell(3,3), vn(3)

#ifdef DEBUG
   call MIO_Debug('DiagBandsG',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   gcell(:,1) = [aG,0.0_dp,0.0_dp]
   gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
   gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]

   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)

   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   call MIO_InputParameter('Bands.NumPoints',nPts0,100)
   if (MIO_InputFindBlock('Bands.Path',nPath)) then
      call MIO_Print('Band calculation','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_InputBlock('Bands.Path',path)
      do ip=1,nPath
         path(:,ip) = path(1,ip)*grcell(:,1) + path(2,ip)*grcell(:,2) + path(3,ip)*grcell(:,3)
         !   !print*, "theta=", theta
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1 ! coz 4 points, 3 segments
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
            nPts(ip) = nPts(ip) + 1 ! Add the missing point at the end of last segment
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot+1],'Kpts','diag') ! add missing point
      Kpts(:,1) = path(:,1)
      ip = 1
      d = 0.0_dp
      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            Kpts(:,ip) = path(:,i) + (j)*(path(:,i+1)-path(:,i))/nPts(i)
            v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
            d = d + sqrt(dot_product(v,v))
         end do
      end do
      call MIO_Allocate(H,[nAt,nAt,nspin],'H','diag')
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.bands'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot+1
      d = 0.0_dp
      hv = -huge(0.0_dp)
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      do ip=1,ptsTot + 1
         do is=1,nspin
            call DiagHam(nAt,nspin,is,H(:,:,is),E(:,is),Kpts(:,ip),gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !!end do
            hv = max(hv,maxval(E(:,is),mask=E(:,is)<=Efermi/g0))
            lc = min(lc,minval(E(:,is),mask=E(:,is)>Efermi/g0))
         end do
         v = Kpts(:,ip) - Kpts(:,max(ip-1,1))
         d = d + sqrt(dot_product(v,v))
         write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)
         if (sum(nPts(:nPath))==ip-1) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ip))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
      end do
      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(H,'H','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagBandsG',1)
#endif /* DEBUG */

end subroutine DiagBandsG

subroutine DiagSpectralFunction()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   integer :: ik, iee, ie
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS, GaussConv
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()

   real(dp) :: GVec(3) !, G1(3), G2(3)

   real(dp) :: eps, factor
   integer :: i1, i2

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunction',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)

   if (MIO_InputFindBlock('Spectral.Path',nPath)) then
      call MIO_Print('Spectral function calculation','diag')
      call MIO_Print('Based on PRB 95, 085420 (2017)','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_InputBlock('Spectral.Path',path)
      do ip=1,nPath
         path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
         !   !print*, "theta=", theta
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptsTot],'KptsG','diag')
      Kpts(:,1) = path(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      KptsG(:,1) = Kpts(:,1) + GVec

      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            Kpts(:,ip) = path(:,i) + (j-1)*(path(:,i+1)-path(:,i))/nPts(i)
            KptsG(:,ip) = Kpts(:,ip) + GVec ! Kpts is in SC, Kpts is for Graphene (PC)
         end do
      end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectral'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot
      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      call MIO_Allocate(Pkc,[1,1],[ptsTot,2],'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      Epts2 = CEILING(Epts/factor)
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptsTot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptsTot,Epts],'Ake','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      is = 1
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      Pkc = 0.0_dp ! spectral weight PkscI(k)
      !$OMP PARALLEL DO PRIVATE (is, ik, nPath, iee, ie)
      do ik=1,ptsTot ! k loop
            call DiagSpectralWeightWeiKu(nAt,nspin,is,Pkc(ik,:),E(:,is),Kpts(:,ik),KptsG(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !!end do
            !hv = max(hv,maxval(E(:,is),mask=E(:,is)<=Efermi/g0)) ! mask restrict search for E smaller than Efermi
            !lc = min(lc,minval(E(:,is),mask=E(:,is)>Efermi/g0))
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                   is = 1
                   !   !definitionDOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
                   if(abs(E(ie,is) - Energy(iee)).lt.(0.005/g0)) then
                        Ake(ik,iee) = Ake(ik,iee) + abs(Pkc(ik,ie))**2
                   end if
                end do
            end do
            !   !if(abs(E(ie,is) - Energy(iee)).lt.(0.1/g0)) then
            AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
            !      !DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-Eig(i1,is))**2/(2.0_dp*eps**2))
            !      !Ake(ik,i2) = Ake(ik,i2) + Pkc(ik,i1)
         !       ! write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)
         !       ! float of lenght 12 with 6 after the comma
         !       ! repeat float of lengt 14 with 6 after the commq 10 times
         !       ! go to next line
         !       ! 10 empty spaces
         !       ! repeat float of lengt 14 with 6 after the commq 10 times, as many times as needed because of ()
         !       ! Ill make it simpler, but maybe bigger file
      end do
      !$OMP END PARALLEL DO
      do ik=1,ptsTot ! k loop
         v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
         d = d + sqrt(dot_product(v,v))
         if (sum(nPts(:nPath))==ik) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ik))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake(ik,iee))
            end do
         end if
      end do

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunction',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunction

subroutine DiagSpectralFunctionKGrid()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is
   integer :: ik, iee, ie
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm

   logical :: MoireBS, GaussConv
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()

   real(dp) :: GVec(3) , G1(3), G2(3)

   real(dp) :: eps, factor
   integer :: i1, i2

   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)

   integer :: cellSize

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGrid',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)

   call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
   gcell(:,1) = [aG,0.0_dp,0.0_dp]
   gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
   gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]

   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   if (MIO_InputFindBlock('Spectral.Path',nPath)) then
      call MIO_Print('Spectral function calculation','diag')
      call MIO_Print('Based on PRB 95, 085420 (2017)','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_Allocate(pathG,[3,nPath],'path','diag')
      call MIO_InputBlock('Spectral.Path',path)
      call MIO_InputBlock('Spectral.Path',pathG)
      do ip=1,nPath
         path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
         pathG(:,ip) = pathG(1,ip)*grcell(:,1) + pathG(2,ip)*grcell(:,2) + pathG(3,ip)*grcell(:,3)
         !   !print*, "theta=", theta
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptsTot],'KptsG','diag')
      KptsG(:,1) = path(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      call MIO_InputParameter('CellSize', cellSize, 1)
      Kpts(1,1) = KptsG(1,1)/cellSize
      Kpts(2,1) = KptsG(2,1)/cellSize
      Kpts(3,1) = KptsG(3,1)

      G1 = matmul(rcell,[1,1,0])
      G2 = matmul(rcell,[1,1,0])
      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            KptsG(:,ip) = pathG(:,i) + (j-1)*(pathG(:,i+1)-pathG(:,i))/nPts(i)

            Kpts(1,ip) = KptsG(1,ip)/cellSize
            Kpts(2,ip) = KptsG(2,ip)/cellSize
            Kpts(3,ip) = KptsG(3,ip)
         end do
      end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectral'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot
      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptsTot,nAt,2],'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      Epts2 = CEILING(Epts/factor)
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptsTot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptsTot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(Ake1,[ptsTot,Epts],'Ake1','diag')
      call MIO_Allocate(AkeGaussian1,[ptsTot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(Ake2,[ptsTot,Epts],'Ake2','diag')
      call MIO_Allocate(AkeGaussian2,[ptsTot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      Ake1 = 0.0_dp
      AkeGaussian1 = 0.0_dp
      Ake2 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      is = 1
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      Pkc = 0.0_dp ! spectral weight PkscI(k)
      do ik=1,ptsTot ! k loop
            call DiagSpectralWeightWeiKu(nAt,nspin,is,Pkc(ik,:,:),E(:,is),Kpts(:,ik),KptsG(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            !!end do
            !hv = max(hv,maxval(E(:,is),mask=E(:,is)<=Efermi/g0)) ! mask restrict search for E smaller than Efermi
            !lc = min(lc,minval(E(:,is),mask=E(:,is)>Efermi/g0))
            !       !if (ie.eq.1) then
            !       !end if
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                   is = 1
                       if(abs(E(ie,is) - Energy(iee)).lt.(0.005/g0)) then
                            Ake1(ik,iee) = Ake1(ik,iee) + abs(Pkc(ik,ie,1))**2
                            Ake2(ik,iee) = Ake2(ik,iee) + abs(Pkc(ik,ie,2))**2
                       end if
                end do
            end do
            Ake = Ake1 + Ake2
            !   !if(abs(E(ie,is) - Energy(iee)).lt.(0.1/g0)) then
            AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
            !      !DOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-Eig(i1,is))**2/(2.0_dp*eps**2))
            !      !Ake(ik,i2) = Ake(ik,i2) + Pkc(ik,i1)
         !       ! write(u,'(f12.6,10f14.6,/,(10x,10f14.6))') d,((E(i,is)*g0, i=1,nAt),is=1,nspin)
         !       ! float of lenght 12 with 6 after the comma
         !       ! repeat float of lengt 14 with 6 after the commq 10 times
         !       ! go to next line
         !       ! 10 empty spaces
         !       ! repeat float of lengt 14 with 6 after the commq 10 times, as many times as needed because of ()
         !       ! Ill make it simpler, but maybe bigger file
         v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
         d = d + sqrt(dot_product(v,v))
         if (sum(nPts(:nPath))==ik) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ik))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake(ik,iee))
            end do
         end if
      end do

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGrid',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGrid

subroutine DiagSpectralFunctionKGridInequivalent()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt, frac, layerIndex
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, u1, u2, u3, u4, uu1,uu2,uu3,uu4,uuu2,uuuu2,uuuuu2,uuuuuu2,uuuuuuu2,uuuuuuuu2
   integer :: ik, iee, ie
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   real(dp), pointer :: KptsGFrac(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm
   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)

   logical :: MoireBS, GaussConv
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake1Loc(:)=>NULL()
   complex(dp), pointer :: Ake2Loc(:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: Ake3(:,:)=>NULL()
   complex(dp), pointer :: Ake4(:,:)=>NULL()
   complex(dp), pointer :: Ake1B(:,:)=>NULL()
   complex(dp), pointer :: Ake2B(:,:)=>NULL()
   complex(dp), pointer :: Ake2C(:,:)=>NULL()
   complex(dp), pointer :: Ake2D(:,:)=>NULL()
   complex(dp), pointer :: Ake2E(:,:)=>NULL()
   complex(dp), pointer :: Ake2F(:,:)=>NULL()
   complex(dp), pointer :: Ake2G(:,:)=>NULL()
   complex(dp), pointer :: Ake2H(:,:)=>NULL()
   complex(dp), pointer :: Ake3B(:,:)=>NULL()
   complex(dp), pointer :: Ake4B(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()

   complex(dp) :: PkcLocA(nAt,4)
   complex(dp) :: PkcLocB(nAt,4)
   complex(dp) :: PkcLocC(nAt,4)
   complex(dp) :: PkcLocD(nAt,4)
   complex(dp) :: PkcLocE(nAt,4)
   complex(dp) :: PkcLocF(nAt,4)
   complex(dp) :: PkcLocG(nAt,4)
   complex(dp) :: PkcLocH(nAt,4)

   real(dp) :: GVec(3) , G01(3), G10(3), G11(3), unfoldedK(3)

   real(dp) :: eps, factor, energyGridResolution
   integer :: i1, i2

   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)
   real(dp) :: gcell1(3,3), gcell2(3,3)
   real(dp) :: rot(3,3)

   integer :: cellSize

   real(dp) :: rcellInv(3,3)

   real(dp) :: ll, kk

   integer :: mmm(4)

   character(len=80) :: line
   integer :: id
   real(dp) :: ELoc1(nAt),ELoc2(nAt)

   real(dp) :: delta, phi, gg, aa, alignmentAngle

   logical :: rotateRefSystem, alignRefSystem, foldByOne, useCoordinates
   logical :: foldByZero

   real(dp) :: topBottomRatio
   logical :: WeiKu, Nishi, Lee, WeiKuOld, useGaussianBroadening

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalent',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)

   call MIO_InputParameter('Spectral.WeiKu',WeiKu,.false.)
   call MIO_InputParameter('Spectral.UseGaussianBroadening',useGaussianBroadening,.false.)
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('Spectral.WeiKuOld',WeiKuOld,.false.)
   call MIO_InputParameter('Spectral.Nishi',Nishi,.false.)
   call MIO_InputParameter('Spectral.Lee',Lee,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByOne,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByZero,.false.)
   call MIO_InputParameter('Spectral.UseCoordinates',useCoordinates,.false.)
   if (MIO_InputSearchLabel('MoireCellParameters',line,id)) then
       call MIO_InputParameter('MoireCellParameters',mmm,[0,0,0,0])
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
       call MIO_InputParameter('Spectral.RotateReferenceSystem',rotateRefSystem,.false.)
       gcell1 = gcell
       if (rotateRefSystem .eqv. .true.) then
           gg = mmm(1)**2 + mmm(2)**2 + mmm(1)*mmm(2)
           delta = sqrt(real(mmm(3)**2 + mmm(4)**2 + mmm(3)*mmm(4))/gg)
           phi = acos((2.0_dp*mmm(1)*mmm(3)+2.0_dp*mmm(2)*mmm(4) + mmm(1)*mmm(4) + mmm(2)*mmm(3))/(2.0_dp*delta*gg))
           phi = -phi*180.0_dp/pi
           aa = phi*pi/180.0_dp
           rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
           rot(:,2) = [sin(aa),cos(aa),0.0_dp]
           rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
           gcell = matmul(rot,gcell)
       end if
   else
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
   end if
   call MIO_InputParameter('Spectral.AlignReferenceSystem',alignRefSystem,.false.)
   if (alignRefSystem .eqv. .true.) then
       call MIO_InputParameter('Spectral.AlignmentAngle',alignmentAngle,0.0d0)
       phi = alignmentAngle
       aa = -phi*pi/180.0_dp
       rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
       rot(:,2) = [sin(aa),cos(aa),0.0_dp]
       rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
       gcell = matmul(rot,gcell)
   end if

   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume
   gcell2 = gcell

   print*, grcell(:,1)
   print*, grcell(:,2)
   print*, grcell(:,3)

   if (MIO_InputFindBlock('Spectral.Path',nPath)) then
      call MIO_Print('Spectral function calculation','diag')
      call MIO_Print('Based on PRB 95, 085420 (2017)','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_Allocate(pathG,[3,nPath],'path','diag')
      call MIO_InputBlock('Spectral.Path',path)
      call MIO_InputBlock('Spectral.Path',pathG)
      do ip=1,nPath
         if (useCoordinates) then
            path(:,ip) = path(:,ip)
            pathG(:,ip) = pathG(:,ip)
         else
            path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
            pathG(:,ip) = pathG(1,ip)*grcell(:,1) + pathG(2,ip)*grcell(:,2) + pathG(3,ip)*grcell(:,3)
         end if
      end do
      print*, "The k-points for the reference system equal in cartesian coordinates:", pathG
      print*, "The k-points for the moire system equal in cartesian coordinates:", path
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptsTot],'KptsG','diag')
      call MIO_Allocate(KptsGFrac,[3,ptsTot],'KptsGFrac','diag')
      KptsG(:,1) = pathG(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      G10 = matmul(rcell,[1,0,0]) ! Lattice vectors of moire cell
      G01 = matmul(rcell,[0,1,0])
      G11 = matmul(rcell,[1,1,0])
      ll = (KptsG(1,1)*G10(2)/G10(1) - KptsG(2,1)) / (G01(1)*G10(2)/G10(1) - G01(2))
      kk = (KptsG(1,1) - ll * G01(1)) / G10(1)
      if (kk.gt.0) then
          kk = floor(kk)
      else
          kk = ceiling(kk)
      end if
      if (ll.gt.0) then
          ll = floor(ll)
      else
          ll = ceiling(ll)
      end if

      if (foldByOne) then
         kpts(1,1) = KptsG(1,1) - G10(1) - G01(1)
         kpts(2,1) = KptsG(2,1) - G10(2) - G01(2)
      else if (foldByZero) then
         kpts(1,1) = KptsG(1,1)
         kpts(2,1) = KptsG(2,1)
      else
         kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1)
         kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2)
      end if

      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            KptsG(:,ip) = pathG(:,i) + (j-1)*(pathG(:,i+1)-pathG(:,i))/nPts(i)
            ll = (KptsG(1,ip)*G10(2)/G10(1) - KptsG(2,ip)) / (G01(1)*G10(2)/G10(1) - G01(2))
            kk = (KptsG(1,ip) - ll * G01(1)) / G10(1)
            if (kk.gt.0) then
                kk = floor(kk)
            else
                kk = ceiling(kk)
            end if
            if (ll.gt.0) then
                ll = floor(ll)
            else
                ll = ceiling(ll)
            end if

            if (foldByOne) then
               kpts(1,ip) = KptsG(1,ip) - G10(1) - G01(1)
               kpts(2,ip) = KptsG(2,ip) - G10(2) - G01(2)
            else if (foldByZero) then
               kpts(1,ip) = KptsG(1,ip)
               kpts(2,ip) = KptsG(2,ip)
            else
               kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1)
               kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2)
            end if

         end do
      end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectralA'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot

      flnm = trim(prefix)//'.spectral1A'
      u1=101
      open(u1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral1B'
      uu1=201
      open(uu1,FILE=flnm,STATUS='replace')

      flnm = trim(prefix)//'.spectral2A'
      u2=102
      open(u2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2B'
      uu2=202
      open(uu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2C'
      uuu2=302
      open(uuu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2D'
      uuuu2=402
      open(uuuu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2E'
      uuuuu2=502
      open(uuuuu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2F'
      uuuuuu2=602
      open(uuuuuu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2G'
      uuuuuuu2=702
      open(uuuuuuu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2H'
      uuuuuuuu2=802
      open(uuuuuuuu2,FILE=flnm,STATUS='replace')

      flnm = trim(prefix)//'.spectral3A'
      u3=103
      open(u3,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3B'
      uu3=203
      open(uu3,FILE=flnm,STATUS='replace')

      flnm = trim(prefix)//'.spectral4A'
      u4=104
      open(u4,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral4B'
      uu4=204
      open(uu4,FILE=flnm,STATUS='replace')

      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptsTot,nAt,2],'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      Epts2 = CEILING(Epts/factor)
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptsTot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptsTot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(AkeGaussian1,[ptsTot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(AkeGaussian2,[ptsTot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      AkeGaussian1 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      is = 1
      call MIO_InputParameter('Spectral.topBottomRatio',topBottomRatio,1.0_dp)
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      call MIO_InputParameter('Spectral.energyGridResolution',energyGridResolution,0.005_dp)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      call MIO_Allocate(Ake1,[ptsTot,Epts],'Ake1','diag')
      call MIO_Allocate(Ake1B,[ptsTot,Epts],'Ake1B','diag')
      call MIO_Allocate(Ake2,[ptsTot,Epts],'Ake2','diag')
      call MIO_Allocate(Ake2B,[ptsTot,Epts],'Ake2B','diag')
      call MIO_Allocate(Ake2C,[ptsTot,Epts],'Ake2C','diag')
      call MIO_Allocate(Ake2D,[ptsTot,Epts],'Ake2D','diag')
      call MIO_Allocate(Ake2E,[ptsTot,Epts],'Ake2E','diag')
      call MIO_Allocate(Ake2F,[ptsTot,Epts],'Ake2F','diag')
      call MIO_Allocate(Ake2G,[ptsTot,Epts],'Ake2G','diag')
      call MIO_Allocate(Ake2H,[ptsTot,Epts],'Ake2H','diag')
      call MIO_Allocate(Ake3,[ptsTot,Epts],'Ake3','diag')
      call MIO_Allocate(Ake3B,[ptsTot,Epts],'Ake3B','diag')
      call MIO_Allocate(Ake4,[ptsTot,Epts],'Ake4','diag')
      call MIO_Allocate(Ake4B,[ptsTot,Epts],'Ake4B','diag')
      Ake1  = 0.0_dp
      Ake1B = 0.0_dp
      Ake2  = 0.0_dp
      Ake2B = 0.0_dp
      Ake2C = 0.0_dp
      Ake2D = 0.0_dp
      Ake2E = 0.0_dp
      Ake2F = 0.0_dp
      Ake2G = 0.0_dp
      Ake2H = 0.0_dp
      Ake3  = 0.0_dp
      Ake3B = 0.0_dp
      Ake4  = 0.0_dp
      Ake4B = 0.0_dp
!HERE
      !$OMP PARALLEL DO PRIVATE(iee, ie, unfoldedK, ELoc, PkcLocA, PkcLocB, PkcLocC,PkcLocD,PkcLocE,PkcLocF,PkcLocG,PkcLocH, KptsLoc), &
      !$OMP& SHARED(KptsG, Kpts, AkeGaussian1, AkeGaussian2, AkeGaussian, nAt, nspin, is, ucell, gcell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, gaussian, Epts, Ake)
      do ik=1,ptsTot ! K loop
         ELoc = 0.0_dp
         ELoc1 = 0.0_dp
         ELoc2 = 0.0_dp
         unfoldedK = KptsG(:,ik)
         PkcLocA = 0.0_dp
         PkcLocB = 0.0_dp
         PkcLocB = 0.0_dp
         PkcLocC = 0.0_dp
         PkcLocD = 0.0_dp
         PkcLocE = 0.0_dp
         PkcLocF = 0.0_dp
         PkcLocG = 0.0_dp
         PkcLocH = 0.0_dp
         KptsLoc = Kpts(:,ik)
            if (WeiKu) then
                if (WeiKuOld) then
                    call DiagSpectralWeightWeiKuInequivalentOld(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                else
                    call DiagSpectralWeightWeiKuInequivalentMoreOrbitals(nAt,nspin,is,PkcLocA,PkcLocB,PkcLocC,PkcLocD,PkcLocE,PkcLocF,PkcLocG,PkcLocH,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                end if
            else if (Lee) then
                call DiagSpectralWeightWeiKuInequivalentLee(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
            else if (Nishi) then
                call DiagSpectralWeightWeiKuInequivalentNishi(nAt,nspin,is,PkcLocA,ELoc1,ELoc2,KptsLoc,unfoldedK,ucell,gcell1,gcell2,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            end if
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                       if (WeiKu) then
                          if (useGaussianBroadening) then
                             !definitionDOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
                             Ake1(ik,iee) = Ake1(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,1))**2
                             Ake2(ik,iee) = Ake2(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,2))**2
                             Ake3(ik,iee) = Ake3(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,3))**2
                             Ake4(ik,iee) = Ake4(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,4))**2
                          else
                             if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                                 ! A sublattice
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2
                                 Ake4(ik,iee) = Ake4(ik,iee) + abs(PkcLocA(ie,4))**2

                                 ! B sublattice
                                 Ake1B(ik,iee) = Ake1B(ik,iee) + abs(PkcLocB(ie,1))**2
                                 Ake2B(ik,iee) = Ake2B(ik,iee) + abs(PkcLocB(ie,2))**2
                                 Ake3B(ik,iee) = Ake3B(ik,iee) + abs(PkcLocB(ie,3))**2
                                 Ake4B(ik,iee) = Ake4B(ik,iee) + abs(PkcLocB(ie,4))**2
                                 ! C sublattice
                                 Ake2C(ik,iee) = Ake2C(ik,iee) + abs(PkcLocC(ie,2))**2
                                 ! D sublattice
                                 Ake2D(ik,iee) = Ake2D(ik,iee) + abs(PkcLocD(ie,2))**2
                                 ! E sublattice
                                 Ake2E(ik,iee) = Ake2E(ik,iee) + abs(PkcLocE(ie,2))**2
                                 ! F sublattice
                                 Ake2F(ik,iee) = Ake2F(ik,iee) + abs(PkcLocF(ie,2))**2
                                 ! G sublattice
                                 Ake2G(ik,iee) = Ake2G(ik,iee) + abs(PkcLocG(ie,2))**2
                                 ! H sublattice
                                 Ake2H(ik,iee) = Ake2H(ik,iee) + abs(PkcLocH(ie,2))**2
                             end if
                          end if
                       else if (Lee) then
                          if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + real(PkcLocA(ie,1))
                             Ake2(ik,iee) = Ake2(ik,iee) + real(PkcLocA(ie,2))
                             Ake3(ik,iee) = Ake3(ik,iee) + real(PkcLocA(ie,3))
                             Ake4(ik,iee) = Ake4(ik,iee) + real(PkcLocA(ie,4))
                          end if
                       else if (Nishi) then
                          if(abs(ELoc1(ie) - Energy(iee)).lt.(energyGridResolution/g0).or.abs(ELoc2(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + PkcLocA(ie,1)
                             Ake2(ik,iee) = Ake2(ik,iee) + PkcLocA(ie,2)
                             Ake3(ik,iee) = Ake3(ik,iee) + PkcLocA(ie,3)
                             Ake4(ik,iee) = Ake4(ik,iee) + PkcLocA(ie,4)
                          end if
                       end if
                end do
            end do
         Ake(ik,:) = Ake1(ik,:) + Ake2(ik,:) + Ake3(ik,:) + Ake4(ik,:)
         AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
      end do
      !$OMP END PARALLEL DO
      do ik=1,ptsTot ! K loop
         v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
         d = d + sqrt(dot_product(v,v))
         if (sum(nPts(:nPath))==ik) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ik))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                   ! A sublattice
                   write(u1,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake1(ik,iee))
                   write(u2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2(ik,iee))
                   write(u3,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake3(ik,iee))
                   write(u4,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake4(ik,iee))
                   ! B sublattice
                   write(uu1,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake1B(ik,iee))
                   write(uu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2B(ik,iee))
                   write(uu3,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake3B(ik,iee))
                   write(uu4,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake4B(ik,iee))
                   ! C sublattice
                   write(uuu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2C(ik,iee))
                   ! D sublattice
                   write(uuuu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2D(ik,iee))
                   ! E sublattice
                   write(uuuuu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2E(ik,iee))
                   ! F sublattice
                   write(uuuuuu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2F(ik,iee))
                   ! G sublattice
                   write(uuuuuuu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2G(ik,iee))
                   ! H sublattice
                   write(uuuuuuuu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2H(ik,iee))
            end do
         end if
      end do

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalent',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGridInequivalent

subroutine DiagSpectralFunctionKGridInequivalent_v2()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt, frac, layerIndex
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, u1, u2, u3, u4, uu1,uu2,uu3,uu4
   integer :: ik, iee, ie
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   real(dp), pointer :: KptsGFrac(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm
   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)

   logical :: MoireBS, GaussConv, DirectBand
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: EAke(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake1Loc(:)=>NULL()
   complex(dp), pointer :: Ake2Loc(:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: Ake3(:,:)=>NULL()
   complex(dp), pointer :: Ake4(:,:)=>NULL()
   complex(dp), pointer :: Ake1B(:,:)=>NULL()
   complex(dp), pointer :: Ake2B(:,:)=>NULL()
   complex(dp), pointer :: Ake3B(:,:)=>NULL()
   complex(dp), pointer :: Ake4B(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()

   complex(dp) :: PkcLocA(nAt,4)
   complex(dp) :: PkcLocB(nAt,4)

   real(dp) :: GVec(3) , G01(3), G10(3), G11(3), unfoldedK(3)

   real(dp) :: eps, factor, energyGridResolution
   integer :: i1, i2

   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)
   real(dp) :: gcell1(3,3), gcell2(3,3)
   real(dp) :: rot(3,3)

   integer :: cellSize

   real(dp) :: rcellInv(3,3)

   real(dp) :: ll, kk

   integer :: mmm(4)

   character(len=80) :: line
   integer :: id
   real(dp) :: ELoc1(nAt),ELoc2(nAt)

   real(dp) :: delta, phi, gg, aa, alignmentAngle

   logical :: rotateRefSystem, alignRefSystem, foldByOne, useCoordinates
   logical :: foldByZero

   real(dp) :: topBottomRatio
   logical :: WeiKu, Nishi, Lee, WeiKuOld, useGaussianBroadening

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalent_v2',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)
   print*, rcell(:,1)
   print*, rcell(:,2)
   print*, rcell(:,3)

   call MIO_InputParameter('Spectral.WeiKu',WeiKu,.false.)
   call MIO_InputParameter('Spectral.UseGaussianBroadening',useGaussianBroadening,.false.)
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('Spectral.WeiKuOld',WeiKuOld,.false.)
   call MIO_InputParameter('Spectral.Nishi',Nishi,.false.)
   call MIO_InputParameter('Spectral.Lee',Lee,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByOne,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByZero,.false.)
   call MIO_InputParameter('Spectral.UseCoordinates',useCoordinates,.false.)
   if (MIO_InputSearchLabel('MoireCellParameters',line,id)) then
       call MIO_InputParameter('MoireCellParameters',mmm,[0,0,0,0])
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
       call MIO_InputParameter('Spectral.RotateReferenceSystem',rotateRefSystem,.false.)
       gcell1 = gcell
       if (rotateRefSystem .eqv. .true.) then
           gg = mmm(1)**2 + mmm(2)**2 + mmm(1)*mmm(2)
           delta = sqrt(real(mmm(3)**2 + mmm(4)**2 + mmm(3)*mmm(4))/gg)
           phi = acos((2.0_dp*mmm(1)*mmm(3)+2.0_dp*mmm(2)*mmm(4) + mmm(1)*mmm(4) + mmm(2)*mmm(3))/(2.0_dp*delta*gg))
           phi = -phi*180.0_dp/pi
           aa = phi*pi/180.0_dp
           rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
           rot(:,2) = [sin(aa),cos(aa),0.0_dp]
           rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
           gcell = matmul(rot,gcell)
       end if
   else
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
   end if
   call MIO_InputParameter('Spectral.AlignReferenceSystem',alignRefSystem,.false.)
   if (alignRefSystem .eqv. .true.) then
       call MIO_InputParameter('Spectral.AlignmentAngle',alignmentAngle,0.0d0)
       phi = alignmentAngle
       aa = -phi*pi/180.0_dp
       rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
       rot(:,2) = [sin(aa),cos(aa),0.0_dp]
       rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
       gcell = matmul(rot,gcell)
   end if

   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume
   gcell2 = gcell

   print*, grcell(:,1)
   print*, grcell(:,2)
   print*, grcell(:,3)

   if (MIO_InputFindBlock('Spectral.Path',nPath)) then
      call MIO_Print('Spectral function calculation','diag')
      call MIO_Print('Based on PRB 95, 085420 (2017)','diag')
      call MIO_Allocate(path,[3,nPath],'path','diag')
      call MIO_Allocate(pathG,[3,nPath],'path','diag')
      call MIO_InputBlock('Spectral.Path',path)
      call MIO_InputBlock('Spectral.Path',pathG)
      do ip=1,nPath
         print*, "0: ", path(:,ip)
         if (useCoordinates) then
            path(:,ip) = path(:,ip)
            pathG(:,ip) = pathG(:,ip)
         else
            path(:,ip) = path(1,ip)*rcell(:,1) + path(2,ip)*rcell(:,2) + path(3,ip)*rcell(:,3)
            pathG(:,ip) = pathG(1,ip)*grcell(:,1) + pathG(2,ip)*grcell(:,2) + pathG(3,ip)*grcell(:,3)
         end if
         print*, "1: ", path(:,ip)
         print*, "2: ", path(:,ip)
      end do
      if (nPath==1) then
         call MIO_Allocate(nPts,1,'nPts','diag')
         nPts(1) = 1
         ptsTot = 1
      else
         call MIO_Allocate(nPts,nPath-1,'nPts','diag')
         nPts(1) = nPts0
         ptsTot = nPts0
         if (nPath > 2) then
            v = path(:,2) - path(:,1)
            d0 = sqrt(dot_product(v,v))
            do ip=2,nPath-1
               v = path(:,ip+1) - path(:,ip)
               d = sqrt(dot_product(v,v))
               nPts(ip) = nint(real(d*nPts0)/real(d0))
               ptsTot = ptsTot + nPts(ip)
            end do
         end if
      end if
      call MIO_Allocate(Kpts,[3,ptsTot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptsTot],'KptsG','diag')
      call MIO_Allocate(KptsGFrac,[3,ptsTot],'KptsGFrac','diag')
      KptsG(:,1) = pathG(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      print*, matmul(rcell,[1,1,0]), matmul(rcell,[1,0,0]), matmul(rcell,[0,1,0])
      G10 = matmul(rcell,[1,0,0]) ! Lattice vectors of moire cell
      G01 = matmul(rcell,[0,1,0])
      G11 = matmul(rcell,[1,1,0])
      ll = (KptsG(1,1)*G10(2)/G10(1) - KptsG(2,1)) / (G01(1)*G10(2)/G10(1) - G01(2))
      kk = (KptsG(1,1) - ll * G01(1)) / G10(1)
      if (kk.gt.0) then
          kk = floor(kk)
      else
          kk = ceiling(kk)
      end if
      if (ll.gt.0) then
          ll = floor(ll)
      else
          ll = ceiling(ll)
      end if

      if (foldByOne) then
         kpts(1,1) = KptsG(1,1) - G10(1) - G01(1)
         kpts(2,1) = KptsG(2,1) - G10(2) - G01(2)
      else if (foldByZero) then
         kpts(1,1) = KptsG(1,1)
         kpts(2,1) = KptsG(2,1)
      else
         kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1)
         kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2)
      end if

      do i=1,nPath-1
         do j=1,nPts(i)
            ip = ip + 1
            KptsG(:,ip) = pathG(:,i) + (j-1)*(pathG(:,i+1)-pathG(:,i))/nPts(i)
            ll = (KptsG(1,ip)*G10(2)/G10(1) - KptsG(2,ip)) / (G01(1)*G10(2)/G10(1) - G01(2))
            kk = (KptsG(1,ip) - ll * G01(1)) / G10(1)
            if (kk.gt.0) then
                kk = floor(kk)
            else
                kk = ceiling(kk)
            end if
            if (ll.gt.0) then
                ll = floor(ll)
            else
                ll = ceiling(ll)
            end if

            if (foldByOne) then
               kpts(1,ip) = KptsG(1,ip) - G10(1) - G01(1)
               kpts(2,ip) = KptsG(2,ip) - G10(2) - G01(2)
            else if (foldByZero) then
               kpts(1,ip) = KptsG(1,ip)
               kpts(2,ip) = KptsG(2,ip)
            else
               kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1)
               kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2)
            end if

         end do
      end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectral_v2'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptsTot

      flnm = trim(prefix)//'.spectral1A'
      u1=101
      open(u1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral1B'
      uu1=201
      open(uu1,FILE=flnm,STATUS='replace')

      flnm = trim(prefix)//'.spectral2A'
      u2=102
      open(u2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2B'
      uu2=202
      open(uu2,FILE=flnm,STATUS='replace')

      flnm = trim(prefix)//'.spectral3A'
      u3=103
      open(u3,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3B'
      uu3=203
      open(uu3,FILE=flnm,STATUS='replace')

      flnm = trim(prefix)//'.spectral4A'
      u4=104
      open(u4,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral4B'
      uu4=204
      open(uu4,FILE=flnm,STATUS='replace')

      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptsTot,nAt,2],'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      Epts2 = CEILING(Epts/factor)
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptsTot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptsTot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(AkeGaussian1,[ptsTot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(AkeGaussian2,[ptsTot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      AkeGaussian1 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      print*, "nspin ", nspin
      is = 1
      call MIO_InputParameter('Spectral.topBottomRatio',topBottomRatio,1.0_dp)
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      call MIO_InputParameter('Spectral.energyGridResolution',energyGridResolution,0.005_dp)
      call MIO_InputParameter('Spectral.DirectBandVal',DirectBand,.false.)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      call MIO_Allocate(Ake1,[ptsTot,Epts],'Ake1','diag')
      call MIO_Allocate(Ake1B,[ptsTot,Epts],'Ake1B','diag')
      call MIO_Allocate(Ake2,[ptsTot,Epts],'Ake2','diag')
      call MIO_Allocate(Ake2B,[ptsTot,Epts],'Ake2B','diag')
      call MIO_Allocate(Ake3,[ptsTot,Epts],'Ake3','diag')
      call MIO_Allocate(Ake3B,[ptsTot,Epts],'Ake3B','diag')
      call MIO_Allocate(Ake4,[ptsTot,Epts],'Ake4','diag')
      call MIO_Allocate(Ake4B,[ptsTot,Epts],'Ake4B','diag')
      Ake1  = 0.0_dp
      Ake1B = 0.0_dp
      Ake2  = 0.0_dp
      Ake2B = 0.0_dp
      Ake3  = 0.0_dp
      Ake3B = 0.0_dp
      Ake4  = 0.0_dp
      Ake4B = 0.0_dp
      if (DirectBand) then
        call MIO_Allocate(EAke,[ptsTot,Epts],'EAke','diag')
        EAke = 0.0_dp
      end if
!HERE
      !$OMP PARALLEL DO PRIVATE(iee, ie, unfoldedK, ELoc, PkcLocA, PkcLocB, KptsLoc), &
      !$OMP& SHARED(KptsG, Kpts, AkeGaussian1, AkeGaussian2, AkeGaussian, nAt, nspin, is, ucell, gcell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, gaussian, Epts, Ake)
      do ik=1,ptsTot ! K loop
         ELoc = 0.0_dp
         ELoc1 = 0.0_dp
         ELoc2 = 0.0_dp
         unfoldedK = KptsG(:,ik)
         PkcLocA = 0.0_dp
         PkcLocB = 0.0_dp
         KptsLoc = Kpts(:,ik)
            if (WeiKu) then
                if (WeiKuOld) then
                    call DiagSpectralWeightWeiKuInequivalentOld(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                else
                    call DiagSpectralWeightWeiKuInequivalent(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                end if
            else if (Lee) then
                call DiagSpectralWeightWeiKuInequivalentLee(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
            else if (Nishi) then
                call DiagSpectralWeightWeiKuInequivalentNishi(nAt,nspin,is,PkcLocA,ELoc1,ELoc2,KptsLoc,unfoldedK,ucell,gcell1,gcell2,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            end if
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                       if (WeiKu) then
                          if (useGaussianBroadening) then
                             !definitionDOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
                             Ake1(ik,iee) = Ake1(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,1))**2
                             Ake2(ik,iee) = Ake2(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,2))**2
                             Ake3(ik,iee) = Ake3(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,3))**2
                          else if (DirectBand) then
                             if (iee == int(ie-(nAt/2-Epts/2))) then
                                 EAke(ik,iee) = ELoc(ie)

                                 ! A sublattice
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2 !* topBottomRatio
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2 !* topBottomRatio
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2 !* topBottomRatio
                                 Ake4(ik,iee) = Ake4(ik,iee) + abs(PkcLocA(ie,4))**2 !* topBottomRatio

                                 ! B sublattice
                                 Ake1B(ik,iee) = Ake1B(ik,iee) + abs(PkcLocB(ie,1))**2 !* topBottomRatio
                                 Ake2B(ik,iee) = Ake2B(ik,iee) + abs(PkcLocB(ie,2))**2 !* topBottomRatio
                                 Ake3B(ik,iee) = Ake3B(ik,iee) + abs(PkcLocB(ie,3))**2 !* topBottomRatio
                                 Ake4B(ik,iee) = Ake4B(ik,iee) + abs(PkcLocB(ie,4))**2 !* topBottomRatio
                             end if
                          else
                             if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                                 ! A sublattice
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2
                                 Ake4(ik,iee) = Ake4(ik,iee) + abs(PkcLocA(ie,4))**2

                                 ! B sublattice
                                 Ake1B(ik,iee) = Ake1B(ik,iee) + abs(PkcLocB(ie,1))**2
                                 Ake2B(ik,iee) = Ake2B(ik,iee) + abs(PkcLocB(ie,2))**2
                                 Ake3B(ik,iee) = Ake3B(ik,iee) + abs(PkcLocB(ie,3))**2
                                 Ake4B(ik,iee) = Ake4B(ik,iee) + abs(PkcLocB(ie,4))**2
                             end if
                          end if
                       else if (Lee) then
                          if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + real(PkcLocA(ie,1))
                             Ake2(ik,iee) = Ake2(ik,iee) + real(PkcLocA(ie,2))
                             Ake3(ik,iee) = Ake3(ik,iee) + real(PkcLocA(ie,3))
                          end if
                       else if (Nishi) then
                          if(abs(ELoc1(ie) - Energy(iee)).lt.(energyGridResolution/g0).or.abs(ELoc2(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + PkcLocA(ie,1)
                             Ake2(ik,iee) = Ake2(ik,iee) + PkcLocA(ie,2)
                             Ake3(ik,iee) = Ake3(ik,iee) + PkcLocA(ie,3)
                          end if
                       end if
                end do
            end do
         Ake(ik,:) = Ake1(ik,:) + Ake2(ik,:) + Ake3(ik,:)
         AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
      end do
      !$OMP END PARALLEL DO
      do ik=1,ptsTot ! K loop
         v = KptsG(:,ik) - KptsG(:,max(ik-1,1))
         d = d + sqrt(dot_product(v,v))
         if (sum(nPts(:nPath))==ik) then
            nPath = nPath+1
            call MIO_Print('Point '//trim(num2str(nPath))//': '//trim(num2str(ik))// &
              '   '//trim(num2str(d,6)),'diag')
         end if
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else if (DirectBand) then
            do iee=1,Epts
                write(u,'(f12.6,f12.6,10f12.6)') d,  REAL(EAke(ik,iee)),           REAL(Ake1(ik,iee)), REAL(Ake1B(ik,iee)),&
                                                                                 & REAL(Ake2(ik,iee)), REAL(Ake2B(ik,iee)),&
                                                                                 & REAL(Ake3(ik,iee)), REAL(Ake3B(ik,iee)),&
                                                                                 & REAL(Ake4(ik,iee)), REAL(Ake4B(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                   ! A sublattice
                   write(u1,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake1(ik,iee))
                   write(u2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2(ik,iee))
                   write(u3,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake3(ik,iee))
                   write(u4,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake4(ik,iee))
                   ! B sublattice
                   write(uu1,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake1B(ik,iee))
                   write(uu2,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake2B(ik,iee))
                   write(uu3,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake3B(ik,iee))
                   write(uu4,'(f12.6,f12.6,f22.6)') d, Energy(iee), REAL(Ake4B(ik,iee))
            end do
         end if
      end do

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)
   end if

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalent_v2',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGridInequivalent_v2

subroutine DiagSpectralFunctionKGridInequivalentEnergyCut()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt, frac, layerIndex
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, u1, u2, u3, u4, uu1, uu2, uu3, uu4
   integer :: ik, iee, ie
   integer :: nk(3), ptot,  i1, i2, i3 !, ik, Epts, u, is, uu
   real(dp), pointer :: Kgrid(:,:)=>NULL()
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   real(dp), pointer :: KptsGFrac(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm
   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   real(dp) :: ELoc1(nAt),ELoc2(nAt)

   logical :: MoireBS, GaussConv
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake1Loc(:)=>NULL()
   complex(dp), pointer :: Ake2Loc(:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: Ake3(:,:)=>NULL()
   complex(dp), pointer :: Ake4(:,:)=>NULL()
   complex(dp), pointer :: Ake1B(:,:)=>NULL()
   complex(dp), pointer :: Ake2B(:,:)=>NULL()
   complex(dp), pointer :: Ake3B(:,:)=>NULL()
   complex(dp), pointer :: Ake4B(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()

   complex(dp) :: PkcLocA(nAt,4)
   complex(dp) :: PkcLocB(nAt,4)

   real(dp) :: GVec(3) , G01(3), G10(3), G11(3), unfoldedK(3)

   real(dp) :: eps, factor, energyGridResolution

   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)
   real(dp) :: gcell1(3,3), gcell2(3,3)
   real(dp) :: rot(3,3)

   integer :: cellSize

   real(dp) :: rcellInv(3,3)

   real(dp) :: ll, kk

   integer :: mmm(4)

   character(len=80) :: line
   integer :: id

   real(dp) :: delta, phi, gg, aa, alignmentAngle

   logical :: rotateRefSystem, alignRefSystem, rotateOpposite, foldByOne
   logical :: foldByZero

   real(dp) :: K1(3)

   real(dp) :: gridCut
   real(dp) :: K1x1
   real(dp) :: K1x2
   real(dp) :: deltaKx
   real(dp) :: K1y1
   real(dp) :: K1y2
   real(dp) :: deltaKy
   real(dp) :: K1z1
   real(dp) :: K1z2
   real(dp) :: deltaKz
   real(dp) :: topBottomRatio

   logical :: WeiKu, Nishi, useCoordinates, WeiKuOld, useGaussianBroadening

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)

   call MIO_InputParameter('Spectral.FoldByOne',foldByOne,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByZero,.false.)
   if (MIO_InputSearchLabel('MoireCellParameters',line,id)) then
       call MIO_InputParameter('MoireCellParameters',mmm,[0,0,0,0])
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
       call MIO_InputParameter('Spectral.RotateReferenceSystem',rotateRefSystem,.false.)
       call MIO_InputParameter('Spectral.RotateOpposite',rotateOpposite,.false.)
       gcell1 = gcell
       if (rotateRefSystem .eqv. .true.) then
           gg = mmm(1)**2 + mmm(2)**2 + mmm(1)*mmm(2)
           delta = sqrt(real(mmm(3)**2 + mmm(4)**2 + mmm(3)*mmm(4))/gg)
           phi = acos((2.0_dp*mmm(1)*mmm(3)+2.0_dp*mmm(2)*mmm(4) + mmm(1)*mmm(4) + mmm(2)*mmm(3))/(2.0_dp*delta*gg))
           if (rotateOpposite) then
               phi = phi*180.0_dp/pi
           else
               phi = -phi*180.0_dp/pi
           end if
           aa = phi*pi/180.0_dp
           rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
           rot(:,2) = [sin(aa),cos(aa),0.0_dp]
           rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
           gcell = matmul(rot,gcell)
       end if
   else
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
   end if
   call MIO_InputParameter('Spectral.AlignReferenceSystem',alignRefSystem,.false.)
   if (alignRefSystem .eqv. .true.) then
       call MIO_InputParameter('Spectral.AlignmentAngle',alignmentAngle,0.0d0)
       phi = alignmentAngle
       aa = -phi*pi/180.0_dp
       rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
       rot(:,2) = [sin(aa),cos(aa),0.0_dp]
       rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
       gcell = matmul(rot,gcell)
   end if

   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   gcell2 = gcell

   call MIO_InputParameter('Spectral.WeiKu',WeiKu,.false.)
   call MIO_InputParameter('Spectral.UseGaussianBroadening',useGaussianBroadening,.false.)
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('Spectral.WeiKuOld',WeiKuOld,.false.)
   call MIO_InputParameter('Spectral.Nishi',Nishi,.false.)

   call MIO_Print('Calculating Spectral function around K1 (2/3,1/3)','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1])
   call MIO_InputParameter('KGridCut',gridCut,0.1_dp)
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   call MIO_InputParameter('Spectral.UseCoordinates',useCoordinates,.false.)
   if (useCoordinates) then
     if (MIO_InputFindBlock('Spectral.Path',nPath)) then
        call MIO_Allocate(path,[3,nPath],'path','diag')
        call MIO_InputBlock('Spectral.Path',path)
        K1 = path(:,1)
     end if
   else
      K1 = grcell(:,1)*2.0_dp/3.0_dp + grcell(:,2)*1.0_dp/3.0_dp + grcell(:,3)*0.0  ! Macro K valley of graphene
   end if
   print*, "K1: ", K1
   K1x1 = K1(1) - gridCut
   K1x2 = K1(1) + gridCut
   deltaKx = (K1x2 - K1x1)/nk(1)
   K1y1 = K1(2) - gridCut
   K1y2 = K1(2) + gridCut
   deltaKy = (K1y2 - K1y1)/nk(2)
   K1z1 = K1(3) - gridCut
   K1z2 = K1(3) + gridCut
   deltaKz = (K1z2 - K1z1)/nk(3)
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      Kgrid(1,ik) = (K1x1 + (i1 * deltaKx))
      Kgrid(2,ik) = (K1y1 + (i2 * deltaKy))
      Kgrid(3,ik) = (K1z1 + (i3 * deltaKz))
   end do; end do; end do

      !   !if (MoireBS) then
      !   !   !print*, "theta=", theta
      !   !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
      !   !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
      !   !end if
      call MIO_Allocate(Kpts,[3,ptot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptot],'KptsG','diag')
      call MIO_Allocate(KptsGFrac,[3,ptot],'KptsGFrac','diag')
      KptsG(:,1) = Kgrid(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      G10 = matmul(rcell,[1,0,0])
      G01 = matmul(rcell,[0,1,0])
      G11 = matmul(rcell,[1,1,0])
      ll = (KptsG(1,1)*G10(2)/G10(1) - KptsG(2,1)) / (G01(1)*G10(2)/G10(1) - G01(2))
      kk = (KptsG(1,1) - ll * G01(1)) / G10(1)
      if (kk.gt.0) then
          kk = floor(kk)
      else
          kk = ceiling(kk)
      end if
      if (ll.gt.0) then
          ll = floor(ll)
      else
          ll = ceiling(ll)
      end if

      if (foldByOne) then
         kpts(1,1) = KptsG(1,1) - G10(1) - G01(1)
         kpts(2,1) = KptsG(2,1) - G10(2) - G01(2)
      else if (foldByZero) then
         kpts(1,1) = KptsG(1,1)
         kpts(2,1) = KptsG(2,1)
      else
         kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1)
         kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2)
      end if

         do ip=1,ptot
            KptsG(:,ip) = Kgrid(:,ip)
            ! 2 equations, 2 unknowns. Bring point back to SC reciprocal cell
            ! using G10 and G01.
            ll = (KptsG(1,ip)*G10(2)/G10(1) - KptsG(2,ip)) / (G01(1)*G10(2)/G10(1) - G01(2))
            kk = (KptsG(1,ip) - ll * G01(1)) / G10(1)
            if (kk.gt.0) then
                kk = floor(kk)
            else
                kk = ceiling(kk)
            end if
            if (ll.gt.0) then
                ll = floor(ll)
            else
                ll = ceiling(ll)
            end if

            if (foldByOne) then
               kpts(1,ip) = KptsG(1,ip) - G10(1) - G01(1)
               kpts(2,ip) = KptsG(2,ip) - G10(2) - G01(2)
            else if (foldByZero) then
               kpts(1,ip) = KptsG(1,ip)
               kpts(2,ip) = KptsG(2,ip)
            else
               kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1)
               kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2)
            end if

            !!!Rat(:,i) = Rat(1,i)*ucell(:,1) + Rat(2,i)*ucell(:,2) + Rat(3,i)*ucell(:,3)

         end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectralA'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptot

      flnm = trim(prefix)//'.spectral1A'
      u1=101
      open(u1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2A'
      u2=102
      open(u2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3A'
      u3=103
      open(u3,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral4A'
      u4=104
      open(u4,FILE=flnm,STATUS='replace')

      flnm = trim(prefix)//'.spectral1B'
      uu1=201
      open(uu1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2B'
      uu2=202
      open(uu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3B'
      uu3=203
      open(uu3,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral4B'
      uu4=204
      open(uu4,FILE=flnm,STATUS='replace')

      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptot,nAt,2],'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      Epts2 = CEILING(Epts/factor)
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(AkeGaussian1,[ptot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(AkeGaussian2,[ptot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      AkeGaussian1 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      is = 1
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      call MIO_InputParameter('Spectral.energyGridResolution',energyGridResolution,0.005_dp)
      call MIO_InputParameter('Spectral.topBottomRatio',topBottomRatio,1.0_dp)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      call MIO_Allocate(Ake1,[ptot,Epts],'Ake1','diag')
      call MIO_Allocate(Ake2,[ptot,Epts],'Ake2','diag')
      call MIO_Allocate(Ake3,[ptot,Epts],'Ake3','diag')
      call MIO_Allocate(Ake4,[ptot,Epts],'Ake4','diag')
      call MIO_Allocate(Ake1B,[ptot,Epts],'Ake1B','diag')
      call MIO_Allocate(Ake2B,[ptot,Epts],'Ake2B','diag')
      call MIO_Allocate(Ake3B,[ptot,Epts],'Ake3B','diag')
      call MIO_Allocate(Ake4B,[ptot,Epts],'Ake4B','diag')
      Ake1 = 0.0_dp
      Ake2 = 0.0_dp
      Ake3 = 0.0_dp
      Ake4 = 0.0_dp
      Ake1B = 0.0_dp
      Ake2B = 0.0_dp
      Ake3B = 0.0_dp
      Ake4B = 0.0_dp
!HERE
      !$OMP PARALLEL DO PRIVATE(iee, ie, unfoldedK, ELoc, PkcLocA, PkcLocB, KptsLoc), &
      !$OMP& SHARED(KptsG, Kpts, AkeGaussian1, AkeGaussian2, AkeGaussian, nAt, nspin, is, ucell, gcell, gcell1, gcell2, H0, maxNeigh, hopp, NList, Nneigh, neighCell, gaussian, Epts, Ake)
      do ik=1,ptot ! K loop
         ELoc = 0.0_dp
         ELoc1 = 0.0_dp
         ELoc2 = 0.0_dp
         unfoldedK = KptsG(:,ik)
         PkcLocA = 0.0_dp
         PkcLocB = 0.0_dp
         KptsLoc = Kpts(:,ik)
            if (WeiKu) then
                if (WeiKuOld) then
                   call DiagSpectralWeightWeiKuInequivalentOld(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                else
                   call DiagSpectralWeightWeiKuInequivalent(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                end if
            else if (Nishi) then
                call DiagSpectralWeightWeiKuInequivalentNishi(nAt,nspin,is,PkcLocA,ELoc1,ELoc2,KptsLoc,unfoldedK,ucell,gcell1,gcell2,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            end if
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                       if (WeiKu) then
                          if (useGaussianBroadening) then
                             !definitionDOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
                             Ake1(ik,iee) = Ake1(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,1))**2
                             Ake2(ik,iee) = Ake2(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,2))**2
                             Ake3(ik,iee) = Ake3(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,3))**2
                             Ake4(ik,iee) = Ake4(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,4))**2
                          else
                             if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                                 ! A sublattice
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2 !* topBottomRatio
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2 !* topBottomRatio
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2 !* topBottomRatio
                                 Ake4(ik,iee) = Ake4(ik,iee) + abs(PkcLocA(ie,4))**2 !* topBottomRatio

                                 ! B sublattice
                                 Ake1B(ik,iee) = Ake1B(ik,iee) + abs(PkcLocB(ie,1))**2 !* topBottomRatio
                                 Ake2B(ik,iee) = Ake2B(ik,iee) + abs(PkcLocB(ie,2))**2 !* topBottomRatio
                                 Ake3B(ik,iee) = Ake3B(ik,iee) + abs(PkcLocB(ie,3))**2 !* topBottomRatio
                                 Ake4B(ik,iee) = Ake4B(ik,iee) + abs(PkcLocB(ie,4))**2 !* topBottomRatio
                             end if
                          end if
                       else if (Nishi) then
                          if(abs(ELoc1(ie) - Energy(iee)).lt.(energyGridResolution/g0).or.abs(ELoc2(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + PkcLocA(ie,1)
                             Ake2(ik,iee) = Ake2(ik,iee) + PkcLocA(ie,2)
                             Ake3(ik,iee) = Ake3(ik,iee) + PkcLocA(ie,3)
                             Ake4(ik,iee) = Ake4(ik,iee) + PkcLocA(ie,4)
                          end if
                       end if
                end do
            end do
         Ake(ik,:) = Ake1(ik,:) + Ake2(ik,:) + Ake3(ik,:) + Ake4(ik,:)
         AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
      end do
      !$OMP END PARALLEL DO
      do ik=1,ptot ! K loop
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                write(u1,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake1(ik,iee))
                write(u2,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake2(ik,iee))
                write(u3,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake3(ik,iee))
                write(u4,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake4(ik,iee))
                write(uu1,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake1B(ik,iee))
                write(uu2,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake2B(ik,iee))
                write(uu3,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake3B(ik,iee))
                write(uu4,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake4B(ik,iee))
            end do
         end if
      end do

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Kgrid,'Kgrid','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGridInequivalentEnergyCut

subroutine DiagSpectralFunctionKGridInequivalentEnergyCut_v2()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt, frac, layerIndex
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, u1, u2, u3, u4, uu1, uu2, uu3, uu4
   integer :: ik, iee, ie
   integer :: nk(3), ptot,  i1, i2, i3 !, ik, Epts, u, is, uu
   real(dp), pointer :: Kgrid(:,:)=>NULL()
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   real(dp), pointer :: KptsGFrac(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm
   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   real(dp) :: ELoc1(nAt),ELoc2(nAt)

   logical :: MoireBS, GaussConv, DirectBand
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake1Loc(:)=>NULL()
   complex(dp), pointer :: Ake2Loc(:)=>NULL()
   complex(dp), pointer :: EAke(:,:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: Ake3(:,:)=>NULL()
   complex(dp), pointer :: Ake4(:,:)=>NULL()
   complex(dp), pointer :: Ake1B(:,:)=>NULL()
   complex(dp), pointer :: Ake2B(:,:)=>NULL()
   complex(dp), pointer :: Ake3B(:,:)=>NULL()
   complex(dp), pointer :: Ake4B(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()

   complex(dp) :: PkcLocA(nAt,4)
   complex(dp) :: PkcLocB(nAt,4)

   real(dp) :: GVec(3) , G01(3), G10(3), G11(3), unfoldedK(3)

   real(dp) :: eps, factor, energyGridResolution

   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)
   real(dp) :: gcell1(3,3), gcell2(3,3)
   real(dp) :: rot(3,3)

   integer :: cellSize

   real(dp) :: rcellInv(3,3)

   real(dp) :: ll, kk

   integer :: mmm(4)

   character(len=80) :: line
   integer :: id

   real(dp) :: delta, phi, gg, aa, alignmentAngle

   logical :: rotateRefSystem, alignRefSystem, rotateOpposite, foldByOne
   logical :: foldByZero

   real(dp) :: K1(3)

   real(dp) :: gridCut
   real(dp) :: K1x1
   real(dp) :: K1x2
   real(dp) :: deltaKx
   real(dp) :: K1y1
   real(dp) :: K1y2
   real(dp) :: deltaKy
   real(dp) :: K1z1
   real(dp) :: K1z2
   real(dp) :: deltaKz
   real(dp) :: topBottomRatio

   logical :: WeiKu, Nishi, useCoordinates, WeiKuOld, useGaussianBroadening

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut_v2',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)
   print*, rcell(:,1)
   print*, rcell(:,2)
   print*, rcell(:,3)

   call MIO_InputParameter('Spectral.FoldByOne',foldByOne,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByZero,.false.)
   if (MIO_InputSearchLabel('MoireCellParameters',line,id)) then
       call MIO_InputParameter('MoireCellParameters',mmm,[0,0,0,0])
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
       call MIO_InputParameter('Spectral.RotateReferenceSystem',rotateRefSystem,.false.)
       call MIO_InputParameter('Spectral.RotateOpposite',rotateOpposite,.false.)
       gcell1 = gcell
       print*,"gcell before rotation ", gcell(:,1)
       print*,"gcell before rotation ", gcell(:,2)
       print*,"gcell before rotation ", gcell(:,3)
       if (rotateRefSystem .eqv. .true.) then
           gg = mmm(1)**2 + mmm(2)**2 + mmm(1)*mmm(2)
           delta = sqrt(real(mmm(3)**2 + mmm(4)**2 + mmm(3)*mmm(4))/gg)
           phi = acos((2.0_dp*mmm(1)*mmm(3)+2.0_dp*mmm(2)*mmm(4) + mmm(1)*mmm(4) + mmm(2)*mmm(3))/(2.0_dp*delta*gg))
           if (rotateOpposite) then
               phi = phi*180.0_dp/pi
           else
               phi = -phi*180.0_dp/pi
           end if
           aa = phi*pi/180.0_dp
           print*, "phi =", phi
           rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
           rot(:,2) = [sin(aa),cos(aa),0.0_dp]
           rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
           gcell = matmul(rot,gcell)
       end if
   else
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
   end if
   call MIO_InputParameter('Spectral.AlignReferenceSystem',alignRefSystem,.false.)
   if (alignRefSystem .eqv. .true.) then
       call MIO_InputParameter('Spectral.AlignmentAngle',alignmentAngle,0.0d0)
       phi = alignmentAngle
       aa = -phi*pi/180.0_dp
       rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
       rot(:,2) = [sin(aa),cos(aa),0.0_dp]
       rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
       gcell = matmul(rot,gcell)
   end if

   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   gcell2 = gcell
   print*,"gcell after rotation ", gcell(:,1)
   print*,"gcell after rotation ", gcell(:,2)
   print*,"gcell after rotation ", gcell(:,3)
   print*,"grcell ", grcell(:,1)
   print*,"grcell ", grcell(:,2)
   print*,"grcell ", grcell(:,3)
   print*,"ucell ", ucell(:,1)
   print*,"ucell ", ucell(:,2)
   print*,"ucell ", ucell(:,3)
   print*,"rcell ", rcell(:,1)
   print*,"rcell ", rcell(:,2)
   print*,"rcell ", rcell(:,3)

   call MIO_InputParameter('Spectral.WeiKu',WeiKu,.false.)
   call MIO_InputParameter('Spectral.UseGaussianBroadening',useGaussianBroadening,.false.)
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('Spectral.WeiKuOld',WeiKuOld,.false.)
   call MIO_InputParameter('Spectral.Nishi',Nishi,.false.)

   call MIO_Print('Calculating Spectral function around K1 (2/3,1/3)','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1])
   call MIO_InputParameter('KGridCut',gridCut,0.1_dp)
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   call MIO_InputParameter('Spectral.UseCoordinates',useCoordinates,.false.)
   if (useCoordinates) then
     if (MIO_InputFindBlock('Spectral.Path',nPath)) then
        call MIO_Allocate(path,[3,nPath],'path','diag')
        call MIO_InputBlock('Spectral.Path',path)
        K1 = path(:,1)
     end if
   else
      K1 = grcell(:,1)*2.0_dp/3.0_dp + grcell(:,2)*1.0_dp/3.0_dp + grcell(:,3)*0.0  ! Macro K valley of graphene
   end if
   print*, "K1: ", K1
   K1x1 = K1(1) - gridCut
   K1x2 = K1(1) + gridCut
   deltaKx = (K1x2 - K1x1)/nk(1)
   K1y1 = K1(2) - gridCut
   K1y2 = K1(2) + gridCut
   deltaKy = (K1y2 - K1y1)/nk(2)
   K1z1 = K1(3) - gridCut
   K1z2 = K1(3) + gridCut
   deltaKz = (K1z2 - K1z1)/nk(3)
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      Kgrid(1,ik) = (K1x1 + (i1 * deltaKx))
      Kgrid(2,ik) = (K1y1 + (i2 * deltaKy))
      Kgrid(3,ik) = (K1z1 + (i3 * deltaKz))
   end do; end do; end do

      !   !if (MoireBS) then
      !   !   !print*, "theta=", theta
      !   !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
      !   !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
      !   !end if
      call MIO_Allocate(Kpts,[3,ptot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptot],'KptsG','diag')
      call MIO_Allocate(KptsGFrac,[3,ptot],'KptsGFrac','diag')
      KptsG(:,1) = Kgrid(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      print*, matmul(rcell,[1,1,0]), matmul(rcell,[1,0,0]), matmul(rcell,[0,1,0])
      G10 = matmul(rcell,[1,0,0])
      G01 = matmul(rcell,[0,1,0])
      G11 = matmul(rcell,[1,1,0])
      ll = (KptsG(1,1)*G10(2)/G10(1) - KptsG(2,1)) / (G01(1)*G10(2)/G10(1) - G01(2))
      kk = (KptsG(1,1) - ll * G01(1)) / G10(1)
      if (kk.gt.0) then
          kk = floor(kk)
      else
          kk = ceiling(kk)
      end if
      if (ll.gt.0) then
          ll = floor(ll)
      else
          ll = ceiling(ll)
      end if

      if (foldByOne) then
         kpts(1,1) = KptsG(1,1) - G10(1) - G01(1)
         kpts(2,1) = KptsG(2,1) - G10(2) - G01(2)
      else if (foldByZero) then
         kpts(1,1) = KptsG(1,1)
         kpts(2,1) = KptsG(2,1)
      else
         kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1)
         kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2)
      end if

         do ip=1,ptot
            KptsG(:,ip) = Kgrid(:,ip)
            ! 2 equations, 2 unknowns. Bring point back to SC reciprocal cell
            ! using G10 and G01.
            ll = (KptsG(1,ip)*G10(2)/G10(1) - KptsG(2,ip)) / (G01(1)*G10(2)/G10(1) - G01(2))
            kk = (KptsG(1,ip) - ll * G01(1)) / G10(1)
            if (kk.gt.0) then
                kk = floor(kk)
            else
                kk = ceiling(kk)
            end if
            if (ll.gt.0) then
                ll = floor(ll)
            else
                ll = ceiling(ll)
            end if

            if (foldByOne) then
               kpts(1,ip) = KptsG(1,ip) - G10(1) - G01(1)
               kpts(2,ip) = KptsG(2,ip) - G10(2) - G01(2)
            else if (foldByZero) then
               kpts(1,ip) = KptsG(1,ip)
               kpts(2,ip) = KptsG(2,ip)
            else
               kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1)
               kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2)
            end if

            !!!Rat(:,i) = Rat(1,i)*ucell(:,1) + Rat(2,i)*ucell(:,2) + Rat(3,i)*ucell(:,3)

         end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectral_v2'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptot

      flnm = trim(prefix)//'.spectral1A'
      u1=101
      open(u1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2A'
      u2=102
      open(u2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3A'
      u3=103
      open(u3,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral4A'
      u4=104
      open(u4,FILE=flnm,STATUS='replace')

      flnm = trim(prefix)//'.spectral1B'
      uu1=201
      open(uu1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2B'
      uu2=202
      open(uu2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3B'
      uu3=203
      open(uu3,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral4B'
      uu4=204
      open(uu4,FILE=flnm,STATUS='replace')

      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptot,nAt,2],'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      Epts2 = CEILING(Epts/factor)
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(AkeGaussian1,[ptot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(AkeGaussian2,[ptot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      AkeGaussian1 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      print*, "nspin ", nspin
      is = 1
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      call MIO_InputParameter('Spectral.energyGridResolution',energyGridResolution,0.005_dp)
      call MIO_InputParameter('Spectral.topBottomRatio',topBottomRatio,1.0_dp)
      call MIO_InputParameter('Spectral.DirectBandVal',DirectBand,.false.)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      call MIO_Allocate(Ake1,[ptot,Epts],'Ake1','diag')
      call MIO_Allocate(Ake2,[ptot,Epts],'Ake2','diag')
      call MIO_Allocate(Ake3,[ptot,Epts],'Ake3','diag')
      call MIO_Allocate(Ake4,[ptot,Epts],'Ake4','diag')
      call MIO_Allocate(Ake1B,[ptot,Epts],'Ake1B','diag')
      call MIO_Allocate(Ake2B,[ptot,Epts],'Ake2B','diag')
      call MIO_Allocate(Ake3B,[ptot,Epts],'Ake3B','diag')
      call MIO_Allocate(Ake4B,[ptot,Epts],'Ake4B','diag')
      Ake1 = 0.0_dp
      Ake2 = 0.0_dp
      Ake3 = 0.0_dp
      Ake4 = 0.0_dp
      Ake1B = 0.0_dp
      Ake2B = 0.0_dp
      Ake3B = 0.0_dp
      Ake4B = 0.0_dp

      if (DirectBand) then
        call MIO_Allocate(EAke,[ptot,Epts],'EAke','diag')
        EAke = 0.0_dp
      end if

!HERE
      !$OMP PARALLEL DO PRIVATE(iee, ie, unfoldedK, ELoc, PkcLocA, PkcLocB, KptsLoc), &
      !$OMP& SHARED(KptsG, Kpts, AkeGaussian1, AkeGaussian2, AkeGaussian, nAt, nspin, is, ucell, gcell, gcell1, gcell2, H0, maxNeigh, hopp, NList, Nneigh, neighCell, gaussian, Epts, Ake)
      do ik=1,ptot ! K loop
         ELoc = 0.0_dp
         ELoc1 = 0.0_dp
         ELoc2 = 0.0_dp
         unfoldedK = KptsG(:,ik)
         PkcLocA = 0.0_dp
         PkcLocB = 0.0_dp
         KptsLoc = Kpts(:,ik)
            if (WeiKu) then
                if (WeiKuOld) then
                   call DiagSpectralWeightWeiKuInequivalentOld(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                else
                   call DiagSpectralWeightWeiKuInequivalent(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                end if
            else if (Nishi) then
                call DiagSpectralWeightWeiKuInequivalentNishi(nAt,nspin,is,PkcLocA,ELoc1,ELoc2,KptsLoc,unfoldedK,ucell,gcell1,gcell2,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            end if
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                       if (WeiKu) then
                          if (useGaussianBroadening) then
                             !definitionDOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
                             Ake1(ik,iee) = Ake1(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,1))**2
                             Ake2(ik,iee) = Ake2(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,2))**2
                             Ake3(ik,iee) = Ake3(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,3))**2
                             Ake4(ik,iee) = Ake4(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,4))**2
                          else if (DirectBand) then
                             if (iee == int(ie-(nAt/2-Epts/2))) then
                                 EAke(ik,iee) = ELoc(ie)

                                 ! A sublattice
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2 !* topBottomRatio
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2 !* topBottomRatio
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2 !* topBottomRatio
                                 Ake4(ik,iee) = Ake4(ik,iee) + abs(PkcLocA(ie,4))**2 !* topBottomRatio

                                 ! B sublattice
                                 Ake1B(ik,iee) = Ake1B(ik,iee) + abs(PkcLocB(ie,1))**2 !* topBottomRatio
                                 Ake2B(ik,iee) = Ake2B(ik,iee) + abs(PkcLocB(ie,2))**2 !* topBottomRatio
                                 Ake3B(ik,iee) = Ake3B(ik,iee) + abs(PkcLocB(ie,3))**2 !* topBottomRatio
                                 Ake4B(ik,iee) = Ake4B(ik,iee) + abs(PkcLocB(ie,4))**2 !* topBottomRatio
                             end if
                          else
                             if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                                 ! A sublattice
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2 !* topBottomRatio
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2 !* topBottomRatio
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2 !* topBottomRatio
                                 Ake4(ik,iee) = Ake4(ik,iee) + abs(PkcLocA(ie,4))**2 !* topBottomRatio

                                 ! B sublattice
                                 Ake1B(ik,iee) = Ake1B(ik,iee) + abs(PkcLocB(ie,1))**2 !* topBottomRatio
                                 Ake2B(ik,iee) = Ake2B(ik,iee) + abs(PkcLocB(ie,2))**2 !* topBottomRatio
                                 Ake3B(ik,iee) = Ake3B(ik,iee) + abs(PkcLocB(ie,3))**2 !* topBottomRatio
                                 Ake4B(ik,iee) = Ake4B(ik,iee) + abs(PkcLocB(ie,4))**2 !* topBottomRatio
                             end if
                          end if
                       else if (Nishi) then
                          if(abs(ELoc1(ie) - Energy(iee)).lt.(energyGridResolution/g0).or.abs(ELoc2(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + PkcLocA(ie,1)
                             Ake2(ik,iee) = Ake2(ik,iee) + PkcLocA(ie,2)
                             Ake3(ik,iee) = Ake3(ik,iee) + PkcLocA(ie,3)
                             Ake4(ik,iee) = Ake4(ik,iee) + PkcLocA(ie,4)
                          end if
                       end if
                end do
            end do
         Ake(ik,:) = Ake1(ik,:) + Ake2(ik,:) + Ake3(ik,:)
         AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
      end do
      !$OMP END PARALLEL DO
      print*, "lets write it all out"
      do ik=1,ptot ! K loop
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else if (DirectBand) then
            do iee=1,Epts
                write(u,'(3f12.6,f12.6,10f12.6)') KptsG(:,ik), REAL(EAke(ik,iee)), REAL(Ake1(ik,iee)), REAL(Ake1B(ik,iee)),&
                                                                                 & REAL(Ake2(ik,iee)), REAL(Ake2B(ik,iee)),&
                                                                                 & REAL(Ake3(ik,iee)), REAL(Ake3B(ik,iee)),&
                                                                                 & REAL(Ake4(ik,iee)), REAL(Ake4B(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                write(u1,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake1(ik,iee))
                write(u2,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake2(ik,iee))
                write(u3,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake3(ik,iee))
                write(u4,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake4(ik,iee))
                write(uu1,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake1B(ik,iee))
                write(uu2,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake2B(ik,iee))
                write(uu3,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake3B(ik,iee))
                write(uu4,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake4B(ik,iee))
            end do
         end if
      end do

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Kgrid,'Kgrid','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut_v2',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGridInequivalentEnergyCut_v2

subroutine DiagSpectralFunctionKGridInequivalentEnergyCutNickDale()

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt, frac, layerIndex
   use ham,                  only : H0, hopp, nspin
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh
   use name,                 only : prefix
   use tbpar,                only : g0
   use constants,            only : pi, twopi
   use math

   integer :: nPts0, nPath, ip, ptsTot, i, j, u, is, u1, u2, u3
   integer :: ik, iee, ie
   integer :: nk(3), ptot,  i1, i2, i3 !, ik, Epts, u, is, uu
   real(dp), pointer :: Kgrid(:,:)=>NULL()
   real(dp), pointer :: path(:,:)=>NULL(), Kpts(:,:)=>NULL(), E(:,:)=>NULL()
   real(dp), pointer :: pathG(:,:)=>NULL()
   real(dp), pointer :: KptsG(:,:)=>NULL()
   real(dp), pointer :: KptsGFrac(:,:)=>NULL()
   integer, pointer :: nPts(:)=>NULL()
   real(dp) :: d0, v(3), d, hv, lc
   complex(dp), pointer :: Hts(:,:,:)=>NULL()
   complex(dp), pointer :: Htsp(:,:,:)=>NULL()
   !type(cl_file) :: file
   character(len=100) :: flnm
   real(dp) :: KptsLoc(3)
   real(dp) :: ELoc(nAt)
   real(dp) :: ELoc1(nAt),ELoc2(nAt)

   logical :: MoireBS, GaussConv
   real(dp) :: theta

   integer :: Epts, Epts2
   real(dp) :: E1, E2
   real(dp), pointer :: Energy(:)=>NULL()
   real(dp), pointer :: gaussian(:)=>NULL()

   complex(dp), pointer :: Pkc(:,:,:)=>NULL()
   complex(dp), pointer :: Ake(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian1(:,:)=>NULL()
   complex(dp), pointer :: Ake1Loc(:)=>NULL()
   complex(dp), pointer :: Ake2Loc(:)=>NULL()
   complex(dp), pointer :: Ake1(:,:)=>NULL()
   complex(dp), pointer :: Ake2(:,:)=>NULL()
   complex(dp), pointer :: Ake3(:,:)=>NULL()
   complex(dp), pointer :: AkeGaussian2(:,:)=>NULL()

   complex(dp) :: PkcLocA(nAt,3)
   complex(dp) :: PkcLocB(nAt,3)

   real(dp) :: GVec(3) , G01(3), G10(3), G11(3), unfoldedK(3)

   real(dp) :: eps, factor, energyGridResolution

   real(dp) :: area, volume, grcell(3,3), aG
   real(dp) :: gcell(3,3), vn(3)
   real(dp) :: gcell1(3,3), gcell2(3,3)
   real(dp) :: rot(3,3)

   integer :: cellSize

   real(dp) :: rcellInv(3,3)

   real(dp) :: ll, kk

   integer :: mmm(4)

   character(len=80) :: line
   integer :: id

   real(dp) :: delta, phi, gg, aa, alignmentAngle

   logical :: rotateRefSystem, alignRefSystem, rotateOpposite, foldByOne
   logical :: foldByZero, lowerGridHalf

   real(dp) :: K1(3)

   real(dp) :: gridCut
   real(dp) :: gridCutX
   real(dp) :: gridCutY
   real(dp) :: K1x1
   real(dp) :: K1x2
   real(dp) :: deltaKx
   real(dp) :: K1y1
   real(dp) :: K1y2
   real(dp) :: deltaKy
   real(dp) :: K1z1
   real(dp) :: K1z2
   real(dp) :: deltaKz
   real(dp) :: topBottomRatio

   logical :: WeiKu, Nishi, useCoordinates, WeiKuOld, useGaussianBroadening

#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('diag')
#endif /* TIMER */

   call MIO_InputParameter('Spectral.NumPoints',nPts0,100)

   call MIO_InputParameter('Spectral.FoldByOne',foldByOne,.false.)
   call MIO_InputParameter('Spectral.FoldByOne',foldByZero,.false.)
   if (MIO_InputSearchLabel('MoireCellParameters',line,id)) then
       call MIO_InputParameter('MoireCellParameters',mmm,[0,0,0,0])
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
       call MIO_InputParameter('Spectral.RotateReferenceSystem',rotateRefSystem,.false.)
       call MIO_InputParameter('Spectral.RotateOpposite',rotateOpposite,.false.)
       gcell1 = gcell
       if (rotateRefSystem .eqv. .true.) then
           gg = mmm(1)**2 + mmm(2)**2 + mmm(1)*mmm(2)
           delta = sqrt(real(mmm(3)**2 + mmm(4)**2 + mmm(3)*mmm(4))/gg)
           phi = acos((2.0_dp*mmm(1)*mmm(3)+2.0_dp*mmm(2)*mmm(4) + mmm(1)*mmm(4) + mmm(2)*mmm(3))/(2.0_dp*delta*gg))
           if (rotateOpposite) then
               phi = phi*180.0_dp/pi
           else
               phi = -phi*180.0_dp/pi
           end if
           aa = phi*pi/180.0_dp
           rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
           rot(:,2) = [sin(aa),cos(aa),0.0_dp]
           rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
           gcell = matmul(rot,gcell)
       end if
   else
       call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
       gcell(:,1) = [aG,0.0_dp,0.0_dp]
       gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
       gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]
   end if
   call MIO_InputParameter('Spectral.AlignReferenceSystem',alignRefSystem,.false.)
   if (alignRefSystem .eqv. .true.) then
       call MIO_InputParameter('Spectral.AlignmentAngle',alignmentAngle,0.0d0)
       phi = alignmentAngle
       aa = -phi*pi/180.0_dp
       rot(:,1) = [cos(aa),-sin(aa),0.0_dp]
       rot(:,2) = [sin(aa),cos(aa),0.0_dp]
       rot(:,3) = [0.0_dp,0.0_dp,1.0_dp]
       gcell = matmul(rot,gcell)
   end if

   vn = CrossProd(gcell(:,1),gcell(:,2))
   volume = dot_product(gcell(:,3),vn)
   area = norm(vn)
   grcell(:,1) = twopi*CrossProd(gcell(:,2),gcell(:,3))/volume
   grcell(:,2) = twopi*CrossProd(gcell(:,3),gcell(:,1))/volume
   grcell(:,3) = twopi*CrossProd(gcell(:,1),gcell(:,2))/volume

   gcell2 = gcell

   call MIO_InputParameter('Spectral.WeiKu',WeiKu,.false.)
   call MIO_InputParameter('Spectral.UseGaussianBroadening',useGaussianBroadening,.false.)
   call MIO_InputParameter('Epsilon',eps,0.01_dp)
   call MIO_InputParameter('Spectral.WeiKuOld',WeiKuOld,.false.)
   call MIO_InputParameter('Spectral.Nishi',Nishi,.false.)

   call MIO_Print('Calculating Spectral function around K1 (2/3,1/3)','diag')
   call MIO_InputParameter('KGrid',nk,[1,1,1])
   call MIO_InputParameter('KGridCut',gridCut,0.1_dp)
   call MIO_InputParameter('KGridCutX',gridCutX,0.1_dp)
   call MIO_InputParameter('KGridCutY',gridCutY,0.1_dp)
   call MIO_InputParameter('KGridLowerGridHalf',lowerGridHalf,.false.)
   ptot = nk(1)*nk(2)*nk(3)
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','diag')
   ik = 0
   call MIO_InputParameter('Spectral.UseCoordinates',useCoordinates,.false.)
   if (useCoordinates) then
     if (MIO_InputFindBlock('Spectral.Path',nPath)) then
        call MIO_Allocate(path,[3,nPath],'path','diag')
        call MIO_InputBlock('Spectral.Path',path)
        K1 = path(:,1)
     end if
   else
      K1 = grcell(:,1)*2.0_dp/3.0_dp + grcell(:,2)*1.0_dp/3.0_dp + grcell(:,3)*0.0
   end if
   print*, "K1: ", K1
   if (lowerGridHalf) then
      K1x1 = K1(1) - gridCutX
      K1x2 = K1(1)
   else
      K1x1 = K1(1)
      K1x2 = K1(1) + gridCutX
   end if
   deltaKx = (K1x2 - K1x1)/nk(1)
   K1y1 = K1(2) + gridCutY
   K1z1 = K1(3) - gridCut
   K1z2 = K1(3) + gridCut
   deltaKz = (K1z2 - K1z1)/nk(3)
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      Kgrid(1,ik) = (K1x1 + (i1 * deltaKx))
      Kgrid(2,ik) = K1y1
      Kgrid(3,ik) = (K1z1 + (i3 * deltaKz))
   end do; end do; end do

      !   !if (MoireBS) then
      !   !   !print*, "theta=", theta
      !   !   call MIO_InputParameter('twistedBilayerAngle',theta,0.0_dp) ! Ref. PRB 76, 73103
      !   !   path(:,ip) = path(:,ip)*theta/180.0_dp*pi
      !   !end if
      call MIO_Allocate(Kpts,[3,ptot],'Kpts','diag')
      call MIO_Allocate(KptsG,[3,ptot],'KptsG','diag')
      call MIO_Allocate(KptsGFrac,[3,ptot],'KptsGFrac','diag')
      KptsG(:,1) = Kgrid(:,1)
      ip = 0
      d = 0.0_dp
      GVec = matmul(rcell,[1,0,0]) ! we only want to translate them by one reciprocal lattice vector
      G10 = matmul(rcell,[1,0,0])
      G01 = matmul(rcell,[0,1,0])
      G11 = matmul(rcell,[1,1,0])
      ll = (KptsG(1,1)*G10(2)/G10(1) - KptsG(2,1)) / (G01(1)*G10(2)/G10(1) - G01(2))
      kk = (KptsG(1,1) - ll * G01(1)) / G10(1)
      if (kk.gt.0) then
          kk = floor(kk)
      else
          kk = ceiling(kk)
      end if
      if (ll.gt.0) then
          ll = floor(ll)
      else
          ll = ceiling(ll)
      end if

      if (foldByOne) then
         kpts(1,1) = KptsG(1,1) - G10(1) - G01(1)
         kpts(2,1) = KptsG(2,1) - G10(2) - G01(2)
      else if (foldByZero) then
         kpts(1,1) = KptsG(1,1)
         kpts(2,1) = KptsG(2,1)
      else
         kpts(1,1) = KptsG(1,1) - kk*G10(1) - ll*G01(1)
         kpts(2,1) = KptsG(2,1) - kk*G10(2) - ll*G01(2)
      end if

         do ip=1,ptot
            KptsG(:,ip) = Kgrid(:,ip)
            ! 2 equations, 2 unknowns. Bring point back to SC reciprocal cell
            ! using G10 and G01.
            ll = (KptsG(1,ip)*G10(2)/G10(1) - KptsG(2,ip)) / (G01(1)*G10(2)/G10(1) - G01(2))
            kk = (KptsG(1,ip) - ll * G01(1)) / G10(1)
            if (kk.gt.0) then
                kk = floor(kk)
            else
                kk = ceiling(kk)
            end if
            if (ll.gt.0) then
                ll = floor(ll)
            else
                ll = ceiling(ll)
            end if

            if (foldByOne) then
               kpts(1,ip) = KptsG(1,ip) - G10(1) - G01(1)
               kpts(2,ip) = KptsG(2,ip) - G10(2) - G01(2)
            else if (foldByZero) then
               kpts(1,ip) = KptsG(1,ip)
               kpts(2,ip) = KptsG(2,ip)
            else
               kpts(1,ip) = KptsG(1,ip) - kk*G10(1) - ll*G01(1)
               kpts(2,ip) = KptsG(2,ip) - kk*G10(2) - ll*G01(2)
            end if

            !!!Rat(:,i) = Rat(1,i)*ucell(:,1) + Rat(2,i)*ucell(:,2) + Rat(3,i)*ucell(:,3)

         end do
      call MIO_Allocate(E,[nAt,nspin],'E','diag')
      flnm = trim(prefix)//'.spectral'
      u=99
      open(u,FILE=flnm,STATUS='replace')
      write(u,'(f16.8)') Efermi
      write(u,'(2f16.8)') 0.0_dp, d
      write(u,'(2f16.8)') Emin-2.0_dp, Emax+2.0_dp
      write(u,'(3i8)') nAt, nspin, ptot
      flnm = trim(prefix)//'.spectral1'
      u1=101
      open(u1,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral2'
      u2=102
      open(u2,FILE=flnm,STATUS='replace')
      flnm = trim(prefix)//'.spectral3'
      u3=103
      open(u3,FILE=flnm,STATUS='replace')
      d = 0.0_dp
      hv = -huge(0.0_dp) ! HUGE(X) returns the largest number that is not an infinity in the model of the type of X.
      lc = huge(0.0_dp)
      call MIO_Print('')
      call MIO_Print('Path with '//trim(num2str(nPath))//' points:','diag')
      nPath = 1
      call MIO_Print('Point 1:   1   '//trim(num2str(0.0_dp,6)),'diag')
      call MIO_Allocate(Pkc,[1,1,1],[ptot,nAt,2],'Pkc','diag')
      call MIO_InputParameter('NumberofEnergyPoints',Epts,1000)
      call MIO_InputParameter('Spectral.Emin',E1,-1.0_dp)
      call MIO_InputParameter('Spectral.Emax',E2,1.0_dp)
      call MIO_Allocate(Energy,Epts,'Energy','diag')
      call MIO_InputParameter('Epsilon',eps,0.01_dp)
      factor = (E2-E1)/(6.0*eps)
      Epts2 = CEILING(Epts/factor)
      if (mod(Epts2,2).ne.0) then
         Epts2 = Epts2+1
      end if
      call MIO_Allocate(gaussian,Epts2,'Energy','diag')
      do iee=1,Epts
           Energy(iee) = E1 + (E2-E1)*(iee-1)/(Epts-1)
      end do
      call MIO_Allocate(Ake,[ptot,Epts],'Ake','diag')
      call MIO_Allocate(AkeGaussian,[ptot,Epts],'AkeGaussian','diag')
      call MIO_Allocate(AkeGaussian1,[ptot,Epts],'AkeGaussian1','diag')
      call MIO_Allocate(AkeGaussian2,[ptot,Epts],'AkeGaussian2','diag')
      Ake = 0.0_dp
      AkeGaussian = 0.0_dp
      AkeGaussian1 = 0.0_dp
      AkeGaussian2 = 0.0_dp
      is = 1
      call MIO_InputParameter('Spectral.GaussianConvolution',GaussConv,.false.)
      call MIO_InputParameter('Spectral.energyGridResolution',energyGridResolution,0.005_dp)
      call MIO_InputParameter('Spectral.topBottomRatio',topBottomRatio,1.0_dp)
      do iee=1,Epts2
          gaussian(iee) = exp(-(Energy(iee)-Energy(Epts2/2))**2/(2.0_dp*eps**2))
      end do
      call MIO_Allocate(Ake1,[ptot,Epts],'Ake1','diag')
      call MIO_Allocate(Ake2,[ptot,Epts],'Ake2','diag')
      call MIO_Allocate(Ake3,[ptot,Epts],'Ake3','diag')
      Ake1 = 0.0_dp
      Ake2 = 0.0_dp
      Ake3 = 0.0_dp
!HERE
      !$OMP PARALLEL DO PRIVATE(iee, ie, unfoldedK, ELoc, PkcLocA, PkcLocB, KptsLoc), &
      !$OMP& SHARED(KptsG, Kpts, AkeGaussian1, AkeGaussian2, AkeGaussian, nAt, nspin, is, ucell, gcell, gcell1, gcell2, H0, maxNeigh, hopp, NList, Nneigh, neighCell, gaussian, Epts, Ake)
      do ik=1,ptot ! K loop
         ELoc = 0.0_dp
         ELoc1 = 0.0_dp
         ELoc2 = 0.0_dp
         unfoldedK = KptsG(:,ik)
         PkcLocA = 0.0_dp
         PkcLocB = 0.0_dp
         KptsLoc = Kpts(:,ik)
            if (WeiKu) then
                if (WeiKuOld) then
                   call DiagSpectralWeightWeiKuInequivalentOld(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                else
                   call DiagSpectralWeightWeiKuInequivalent(nAt,nspin,is,PkcLocA,PkcLocB,ELoc,KptsLoc,unfoldedK,ucell,gcell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,topBottomRatio)
                end if
            else if (Nishi) then
                call DiagSpectralWeightWeiKuInequivalentNishi(nAt,nspin,is,PkcLocA,ELoc1,ELoc2,KptsLoc,unfoldedK,ucell,gcell1,gcell2,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
            end if
            do iee=1,Epts  ! epsilon
                do ie=1,nAt   ! epsilonIksc
                       if (WeiKu) then
                          if (useGaussianBroadening) then
                             !definitionDOS(i2,is) = DOS(i2,is) + exp(-(E(i2)-EStore(ik,i1))**2/(2.0_dp*eps**2))
                             Ake1(ik,iee) = Ake1(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,1))**2
                             Ake2(ik,iee) = Ake2(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,2))**2
                             Ake3(ik,iee) = Ake3(ik,iee) + exp(-(ELoc(ie) - Energy(iee))**2.0/(2.0_dp*eps**2)) * abs(PkcLocA(ie,3))**2
                          else
                             if(abs(ELoc(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                                 Ake1(ik,iee) = Ake1(ik,iee) + abs(PkcLocA(ie,1))**2 !* topBottomRatio
                                 Ake2(ik,iee) = Ake2(ik,iee) + abs(PkcLocA(ie,2))**2 !* topBottomRatio
                                 Ake3(ik,iee) = Ake3(ik,iee) + abs(PkcLocA(ie,3))**2 !* topBottomRatio
                             end if
                          end if
                       else if (Nishi) then
                          if(abs(ELoc1(ie) - Energy(iee)).lt.(energyGridResolution/g0).or.abs(ELoc2(ie) - Energy(iee)).lt.(energyGridResolution/g0)) then
                             Ake1(ik,iee) = Ake1(ik,iee) + PkcLocA(ie,1)
                             Ake2(ik,iee) = Ake2(ik,iee) + PkcLocA(ie,2)
                             Ake3(ik,iee) = Ake3(ik,iee) + PkcLocA(ie,3)
                          end if
                       end if
                end do
            end do
         Ake(ik,:) = Ake1(ik,:) + Ake2(ik,:) + Ake3(ik,:)
         AkeGaussian(ik,:) = convolve(real(Ake(ik,:)),gaussian,Epts)
      end do
      !$OMP END PARALLEL DO
      print*, "lets write it all out"
      do ik=1,ptot ! K loop
         if (GaussConv) then
            do iee=1,Epts  ! epsilon
                write(u,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(AkeGaussian(ik,iee))
            end do
         else
            do iee=1,Epts  ! epsilon
                write(u,'(3f12.6,f12.6,f12.6)') KptsG(:,ik), Energy(iee), REAL(Ake(ik,iee))
                write(u1,'(f12.6,f12.6,f22.6)') KptsG(:,ik), Energy(iee), REAL(Ake1(ik,iee))
                write(u2,'(f12.6,f12.6,f22.6)') KptsG(:,ik), Energy(iee), REAL(Ake2(ik,iee))
                write(u3,'(f12.6,f12.6,f22.6)') KptsG(:,ik), Energy(iee), REAL(Ake3(ik,iee))
            end do
         end if
      end do

      call MIO_Print('')
      !call file%Close()
      call MIO_Deallocate(E,'E','diag')
      call MIO_Deallocate(Kpts,'Ktsp','diag')
      call MIO_Deallocate(KptsG,'KtspG','diag')
      call MIO_Deallocate(Kgrid,'Kgrid','diag')
      call MIO_Deallocate(Energy,'Energy','diag')
      call MIO_Print('Band gap: '//trim(num2str(g0*(lc-hv),5)),'diag')
      call MIO_Print('')
      close(u)

#ifdef TIMER
   call MIO_TimerStop('diag')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('DiagSpectralFunctionKGridInequivalentEnergyCut',1)
#endif /* DEBUG */

end subroutine DiagSpectralFunctionKGridInequivalentEnergyCutNickDale

subroutine DiagHam(N,ns,is,HLoc,ELoc,KLoc,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species
   use tbpar,                 only : U

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp), intent(out) :: HLoc(N,N)
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   real(dp) :: R(3), zz

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      ! Add SCF terms if spin-polarized (matches other routines)
      ! Commented out: not doing any SCF calculation for now (matches BuildBlockHamiltonianOnly)
      !   !zz = charge(1,i)*charge(2,i) ! Zch
      ! Apply SOC modifications
      if (anySOCEnabled) call ApplySOCtoHamiltonian(i, is, ns, HLoc)

      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(-cmplx_i*dot_product(KLoc,R))
         ! PIA hopping not yet properly implemented - commented out
      end do
   end do
   if (edgeHopp) then
      do i=1,nQ
         do j=1,nEdgeN(i)
            in = NeI(j,i)
            R = matmul(cell,NedgeCell(:,j,i))
            HLoc(in,edgeIndx(i)) = HLoc(in,edgeIndx(i)) + edgeH(j,i)*exp(-cmplx_i*dot_product(KLoc,R))
         end do
      end do
   end if
   call ZHEEV('N','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   if (info/=0) then
      call MIO_Kill('Error in diagonalization','diag','DiagHam')
   end if

end subroutine DiagHam

subroutine DiagHamSparse(N, ns, is, ELoc, KLoc, cell, H0, maxN, hopp, NList, Nneigh, neighCell,neig)
    use constants, only : cmplx_i
    use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
    use scf, only : charge, Zch
    use atoms, only : Species
    use tbpar, only : U, g0
    use name, only : prefix
    implicit none
    integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is, neig
    real(dp), intent(out) :: ELoc(N)
    real(dp), intent(in) :: KLoc(3), cell(3,3)
    real(dp), intent(inout) :: H0(N)
    complex(dp), intent(in) :: hopp(maxN,N)
    real(dp) :: tol

    integer :: i, j, in, info
    real(dp) :: R(3), zz, resid_norm

    ! ARPACK parameters and variables
    integer :: nev, ncv, lworkl, ido, ierr
    character(1) :: bmat
    character(2) :: which
    complex(dp), allocatable :: resid(:), v(:,:), workd(:), d(:), workev(:)
    complex(dp), allocatable :: workl(:)
    complex(dp), allocatable :: z(:,:)
    double precision, allocatable :: rwork(:)
    integer, allocatable :: iparam(:), ipntr(:)
    logical, allocatable :: select(:)
    integer :: max_iter, nn, iter
    logical :: useShift
    double precision :: shift

    integer, allocatable :: row_ptr(:), col_ind(:)
    complex(dp), allocatable :: values(:)
    real(dp), allocatable :: rand_real(:), rand_imag(:)
    complex(dp) :: sigma

    ! PARDISO variables

    INTEGER :: nnn, nnz

!C.. Internal solver memory pointer for 64-bit architectures
    INTEGER*8 pt(64)
!C.. Internal solver memory pointer for 32-bit architectures
!C.. INTEGER*4 pt(64)

!C.. This is OK in both cases
!    TYPE(MKL_PARDISO_HANDLE) pt(64)
!C.. All other variables
    INTEGER maxfct, mnum, mtype, phase, nrhs, msglvl
    INTEGER iparm(64)
    INTEGER, allocatable :: ia(:) ! row_ptr
    INTEGER, allocatable :: ja(:) ! col_ind
    INTEGER idum(1)
    complex(dp) ::  ddum(1)
    integer error
    complex(dp), allocatable :: a(:) ! values

    DATA nrhs /1/, maxfct /1/, mnum /1/

    complex(dp), allocatable :: EVectors(:,:)
    logical :: saveRitz, symmetric
    integer :: nconv
    complex(dp), allocatable :: ax(:)
    double precision, allocatable :: rd(:,:)
    double precision :: dznrm2, dlapy2

    ! Shift-invert branch (Diag.SparseUseShift)
    logical :: countStates, countOK
    integer :: nbelow(2), nbelow0(1), ichk(2), m5, kk
    real(dp) :: dEchk(2), dE0(1)
    real(dp), allocatable :: esort(:)
    integer, save :: se_unit
    logical, save :: se_open = .false.

    ! Parameters

    ! Variables
    !type(C_PTR) :: mkl_handle
    !!integer(int_kind) :: N, status, nnz, l, i, j, k, row_start, row_end, temp_index
    !!integer(int_kind), allocatable :: row_ptr(:), col_ind(:)
    !type(SPARSE_MATRIX_DESCR) :: descr

    nev=neig
    ncv=nev*10
    ! Shift-invert converges on the levels nearest the shift in a few restarts: a 2*nev+1 basis is enough, and the
    ! regular-mode 10*nev would not fit in memory for the hundreds of levels this mode is meant for.
    call MIO_InputParameter('Diag.SparseUseShift',useShift,.false.)
    if (useShift) ncv = min(N, max(2*nev+1, 20))
    lworkl= 3*NCV**2 + 5*NCV

    max_iter = 10000
    allocate(resid(N), v(N, ncv), workd(3*N), workl(lworkl), rwork(ncv), d(nev+1), iparam(11), ipntr(14), select(ncv), rand_real(N), rand_imag(N))
    allocate(z(N,nev))
    allocate(workev(2*ncv))
    allocate(ax(N))
    allocate(rd(ncv,3))

    ! Ensure correct size of nev and ncv
    if ( (nev < 1) .or. (nev >= ncv) .or. (ncv > N) ) then
        print *, 'Error: invalid parameters - nev=', nev, 'ncv=', ncv, 'N=', N
        error stop 1
    end if

    ! Initialize arrays

    bmat = 'I'

    iparam(1) = 1
    iparam(3) = max_iter
    call MIO_InputParameter('Diag.SparseSetShift',shift,0.0_dp)
    call MIO_InputParameter('Diag.SparseUseShift',useShift,.false.)
    call MIO_InputParameter('Diag.SparseSaveRitz',saveRitz,.false.)
    if (useShift) then
       which = 'LM'
       iparam(7) = 3
       ! Diag.SparseSetShift is in eV, like every other energy in Gendata.in; H is in units of g0 internally
       sigma = cmplx(shift/g0, 0.0_dp, kind=dp)
       call MIO_InputParameter('Diag.SparseCount',countStates,.false.)
       call MIO_Print('DiagHamSparse: shift-invert, '//trim(num2str(nev))//' eigenvalues nearest Diag.SparseSetShift = '// &
          trim(num2str(shift,6))//' eV (PARDISO, complex Hermitian)','diag')
    else
       which = 'SM'
       sigma = (0.0_dp, 0.0_dp)
       iparam(7) = 1
    end if
    ! Initialize ARPACK parameters

    if (useShift) then
       ! the regular-mode default 0.1 is far too loose for a list of interior eigenvalues
       call MIO_InputParameter('Diag.SparseSetTol',tol,1.0e-9_dp)
    else
       call MIO_InputParameter('Diag.SparseSetTol',tol,0.1_dp)
    end if

    ido = 0
    info = 0

    nn = N

    ! Initialize the starting vector resid with random values

    ! Check resid for initial state
    if (size(resid) /= N) then
        print *, 'Error: resid size mismatch: ', size(resid), ' expected: ', N
        error stop 1
    end if

     ! Debug print for resid
    if (any(resid /= resid)) then
        print *, 'Error: resid contains NaN values initially.'
        error stop 1
    end if

    resid_norm = sqrt(sum(abs(resid)**2))

    ! Debug prints

    ! Create CSR sparse matrix storage
    print *, "initialize the sparse matrix and put it in csr format"
    call initialize_sparse_matrix(N, maxN, H0, hopp, NList, Nneigh, neighCell, ns, is, KLoc, cell, row_ptr, col_ind, values,sigma)
    print *, "done"
    symmetric = is_structurally_symmetric(values, row_ptr, col_ind, N)

    print*, "is it symmetric?", symmetric

    ! Debug prints for CSR matrix
    if (any(values /= values)) then
        print *, 'Error: values contains NaN values after initialization.'
        error stop 1
    end if

    if (maxval(abs(values)) > 1e10) then
        print *, 'Warning: values contains extremely large values.'
    end if

    if (useShift) then
       ! values = H - sigma (initialize_sparse_matrix). PARDISO wants the upper triangle of the Hermitian matrix,
       ! columns sorted; the solves below then apply OP = (H - sigma)^(-1).
       call sparse_upper_sorted(N, row_ptr, col_ind, values, ia, ja, a)
       if (countStates) then
          dE0(1) = 0.0_dp
          call sparse_count_below(N, ia, ja, a, 1, dE0, nbelow0, countOK)
       end if
       pt = 0
       iparm = 0
       iparm(1) = 1     ! no solver defaults
       iparm(2) = 3     ! parallel nested dissection
       iparm(10) = 8    ! pivot perturbation 1e-8, used only if a pivot fails (reported in iparm(14))
       iparm(18) = -1   ! report the number of non-zeros of the factor
       iparm(21) = 1    ! Bunch-Kaufman 1x1 and 2x2 pivots
       mtype = -4       ! complex Hermitian indefinite
       msglvl = 0
       nnn = N
       phase = 12       ! analysis + numerical factorisation
       call pardiso(pt, maxfct, mnum, mtype, phase, nnn, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
       if (error /= 0) then
          print *, 'DiagHamSparse: PARDISO factorisation error ', error
          error stop 1
       end if
       call MIO_Print('DiagHamSparse: factor of H - shift: '//trim(num2str(iparm(18)))//' non-zeros, '// &
          trim(num2str(iparm(14)))//' perturbed pivots','diag')
    end if
    call znaupd(ido, bmat, nn, which, nev, tol, resid, ncv, v, nn, iparam, ipntr, workd, workl, lworkl, rwork, info)

    resid_norm = sqrt(sum(abs(resid)**2))

    ! Debug prints

    ! Check for convergence and errors
    if (info == -5) then
        print *, 'Maximum number of iterations reached, INFO = ', info
        print *, 'Number of converged eigenvalues:', iparam(5)
        print *, 'iparam: ', iparam
        print *, 'ipntr: ', ipntr
        print *, 'ido: ', ido
        print *, 'bmat: ', bmat
        print *, 'nn: ', nn
        print *, 'which: ', which
        print *, 'nev: ', nev
        print *, 'tol: ', tol
    else if (info /= 0) then
        print *, 'Error with znaupd, INFO = ', info
        print *, 'iparam: ', iparam
        print *, 'ipntr: ', ipntr
        error stop 1
    else
        print *, 'znaupd converged successfully'
    endif

    !    ! Sort each row
    !                    ! Swap indices
    !                    ! Swap values
    !    ! Set matrix descriptor
    !    descr.type = SPARSE_MATRIX_TYPE_GENERAL
    !    descr.mode = SPARSE_FILL_MODE_LOWER
    !    descr.diag = SPARSE_DIAG_NON_UNIT

    !    ! Create the matrix handle

    !    ! Perform the analysis phase

    !    ! Optimize the matrix structure
    !!   ! PARDISO initialization
    !!       nnn = nn
    !!       !nrhs = 1
    !!       a = values
    !!       ia = row_ptr
    !!       ja = col_ind
    !!       maxfct = 1
    !!       mnum = 1
    !!       msglvl = 0
    !!       error = 0
    !!
    !!       !call pardisoinit(pt, mtype, iparm)
    !!          iparm(i) = 0
    !!          pt(i) = 0
    !!       END DO
    !!       iparm(1) = 1 ! no solver default
    !!       iparm(2) = 0 ! fill-in reordering from METIS
    !!       iparm(3) = 1 ! numbers of processors
    !!       iparm(4) = 0 ! no iterative-direct algorithm
    !!       iparm(5) = 0 ! no user fill-in reducing permutation
    !!       iparm(6) = 0 ! =0 solution on the first n components of x
    !!       iparm(7) = 0 ! not in use
    !!       iparm(8) = 9 ! numbers of iterative refinement steps
    !!       iparm(9) = 0 ! not in use
    !!       iparm(10) = 13 ! perturb the pivot elements with 1E-13
    !!       !iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS ! try changing this
    !!       iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS ! try changing this
    !!       iparm(12) = 0 ! not in use
    !!       !iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric) ! try changin this
    !!       iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric) ! try changin this
    !!       iparm(14) = 0 ! Output: number of perturbed pivots
    !!       iparm(15) = 0 ! not in use
    !!       iparm(16) = 0 ! not in use
    !!       iparm(17) = 0 ! not in use
    !!       iparm(18) = -1 ! Output: number of nonzeros in the factor LU
    !!       iparm(19) = -1 ! Output: Mflops for LU factorization
    !!       iparm(20) = 0 ! Output: Numbers of CG Iterations
    !!       error = 0 ! initialize error flag
    !!       msglvl = 0 ! print statistical information
    !!       !mtype = 11 ! real unsymmetric
    !!       !mtype = 13 ! complex and structurally nonsymmetric
    !!       mtype = 6 ! complex and structurally symmetric

    !!       !DO i = 1, 64
    !!       !  pt(i)%DUMMY = 0
    !!       !END DO
    !!
    !!       phase = 11  ! Reordering and Symbolic Factorization
    !!           stop
    !!       end if
    !!
    !!       phase = 22  ! Numerical factorization
    !!           stop
    !!       end if

    iter = 0

    if (useShift) then
       do while (ido /= 99)
           if (ido == -1 .or. ido == 1) then
               ! y = (H - sigma)^(-1) x with the factor computed above
               phase = 33
               call pardiso(pt, maxfct, mnum, mtype, phase, nnn, a, ia, ja, idum, nrhs, iparm, msglvl, &
                            workd(ipntr(1)), workd(ipntr(2)), error)
               if (error /= 0) then
                  print *, 'DiagHamSparse: PARDISO solve error ', error
                  error stop 1
               end if
               !! Solve the linear system (A - sigma*I) * y = x using PARDISO
               !!phase = 33  ! Back substitution and iterative refinement
               !!    stop
               !!end if
           else
               print *, 'Error: ido has unexpected value ', ido
               error stop 1
           end if

           ! ARPACK iteration
           call znaupd(ido, bmat, nn, which, nev, tol, resid, ncv, v, nn, iparam, ipntr, workd, workl, lworkl, rwork, info)
           if (info /= 0) then
               print *, 'Error with znaupd during iteration, info = ', info
               error stop 1
           end if
!          ! resid_norm = sqrt(sum(abs(workd(ipntr(2):ipntr(2) + nn - 1))**2))
       end do
       print*, 'Solve completed ... '
    else
       do while (ido /= 99)
           if (ido == -1 .or. ido == 1) then
               ! Print the input vector for debugging

               ! Perform sparse matrix-vector multiplication
               call sparse_matvec(nn, row_ptr, col_ind, values, &
                                  workd(ipntr(1):ipntr(1)+nn-1), &
                                  workd(ipntr(2):ipntr(2)+nn-1))

               ! Print the output vector for debugging
           else
               print *, 'Error: ido has unexpected value ', ido
               error stop 1
           end if

           ! ARPACK iteration
           call znaupd(ido, bmat, nn, which, nev, tol, resid, ncv, v, nn, iparam, ipntr, workd, workl, lworkl, rwork, info)
           if (info /= 0) then
               print *, 'Error with znaupd during iteration, info = ', info
               error stop 1
           end if
       end do
    end if

    !       !call cgttrs('N', n, 1, dl, dd, du, du2, ipiv, workd(ipntr(2)), n, ierr)

    ! Debug print for extreme values
    if (maxval(abs(resid)) > 1e10) then
        print *, 'Warning: resid contains extremely large values.'
    end if

    if (saveRitz) then
       allocate(EVectors(nn, nev))
       call zneupd(.true.,'A', select, d, z, nn, sigma, workev, bmat,nn, which, nev, tol, resid, ncv, v, nn, iparam, ipntr,workd,workl, lworkl, rwork, info )
       nconv = iparam(5)
       do j=1, nconv
          call av(n, v(1,j), ax)
          call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
          rd(j,1) = dble (d(j))
          rd(j,2) = dimag (d(j))
          rd(j,3) = dznrm2 (n, ax, 1)
          rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
       end do
       call dmout (6, nconv, 3, rd, ncv, -6, 'Ritz values (Real, Imag) and relative residuals')
       ELoc(:nev) = real(d(:nev))
       EVectors = reshape(z, (/nn, nev/))
       ! Open file for writing
       open(unit=10, file='eigenvectors.txt', status='replace')

       ! Write eigenvectors to file
       ! Write eigenvectors to file
       print*, "ucell:", cell
       do i = 1, nev
           do j = 1, nn
               write(10, '(ES20.12,ES20.12)', advance='no') real(EVectors(j, i)), aimag(EVectors(j, i))
           end do
           write(10, *)
       end do
       close(10)
    else
       call zneupd(.false.,'A', select, d, z, nn, sigma, workev, bmat,nn, which, nev, tol, resid, ncv, v, nn, iparam, ipntr,workd,workl, lworkl, rwork, info )

       ELoc(:nev) = real(d(:nev))
    end if
    ! Debug prints after zneupd
    if (any(d /= d)) then
        print *, 'Error: d contains NaN values after zneupd.'
    endif
    if (info /= 0) then
        print *, 'Error with zneupd, ierr = ', info
        error stop 1
    end if
    print *, "eigenvalues: ", d

    if (useShift) then
       phase = -1  ! release the solve factor before the count matrices are factorised
       call pardiso(pt, maxfct, mnum, mtype, phase, nnn, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
       ! ARPACK returns the levels unordered; d already holds E = sigma + 1/theta (zneupd, mode 3)
       nconv = min(iparam(5), nev)
       allocate(esort(nconv))
       esort = real(d(1:nconv))
       call sparse_sort_real(nconv, esort)
       kk = count(esort < real(sigma))        ! returned levels below the shift
       if (countStates) then
          ! Two-way check that no copy of a degenerate level was dropped (implicitly restarted Arnoldi can do that):
          ! the exact count in the widest level gap of each outer fifth must equal the index bookkeeping.
          m5 = max(nconv/5, 2)
          ichk = 0
          if (nconv >= 4) then
             ichk(1) = maxloc(esort(2:m5) - esort(1:m5-1), 1)
             ichk(2) = nconv - m5 + maxloc(esort(nconv-m5+2:nconv) - esort(nconv-m5+1:nconv-1), 1)
             dEchk(1) = 0.5_dp*(esort(ichk(1)) + esort(ichk(1)+1)) - real(sigma)
             dEchk(2) = 0.5_dp*(esort(ichk(2)) + esort(ichk(2)+1)) - real(sigma)
             call sparse_count_below(N, ia, ja, a, 2, dEchk, nbelow, countOK)
             countOK = countOK .and. all(nbelow == nbelow0(1) - kk + ichk)
          end if
          call MIO_Print('DiagHamSparse: '//trim(num2str(nbelow0(1)))//' states below the shift; returned levels have '// &
             'absolute indices '//trim(num2str(nbelow0(1)-kk+1))//' .. '//trim(num2str(nbelow0(1)-kk+nconv)),'diag')
          if (nconv >= 4) then
             if (countOK) then
                call MIO_Print('DiagHamSparse: count check at two energies: OK (complete)','diag')
             else
                call MIO_Print('DiagHamSparse: count check at two energies: MISMATCH (levels missing, or perturbed '// &
                   'pivots) - counts '//trim(num2str(nbelow(1)))//', '//trim(num2str(nbelow(2)))//' vs bookkeeping '// &
                   trim(num2str(nbelow0(1)-kk+ichk(1)))//', '//trim(num2str(nbelow0(1)-kk+ichk(2))),'diag')
             end if
          end if
       end if
       ! <prefix>.SparseEig: one block per k, rows "index E[eV]"; index is absolute (1 = bottom of the spectrum) with
       ! Diag.SparseCount, otherwise counted from the first returned level
       !$OMP CRITICAL (sparseeig_write)
       if (.not. se_open) then
          open(newunit=se_unit, file=trim(prefix)//'.SparseEig', status='replace')
          write(se_unit,'(a)') '# shift-invert eigenvalues (DiagHamSparse). Blocks: "# k kx ky kz  shift[eV]  nconv  '// &
             'nbelow(-1 = not counted)  complete(T/F/-)", then rows: index  E[eV]'
          se_open = .true.
       end if
       if (countStates) then
          write(se_unit,'(a,3f14.8,f16.8,2i10,1x,l1)') '# k ', KLoc, shift, nconv, nbelow0(1), countOK
       else
          write(se_unit,'(a,3f14.8,f16.8,2i10,1x,a)') '# k ', KLoc, shift, nconv, -1, '-'
       end if
       do i = 1, nconv
          if (countStates) then
             write(se_unit,'(i10,f18.9)') nbelow0(1) - kk + i, esort(i)*g0
          else
             write(se_unit,'(i10,f18.9)') i, esort(i)*g0
          end if
       end do
       flush(se_unit)
       !$OMP END CRITICAL (sparseeig_write)
       deallocate(esort, ia, ja, a)
    end if

     deallocate(resid, v, workd, workl, rwork, d, iparam, ipntr, select, row_ptr, col_ind, values, rand_real, rand_imag)

end subroutine DiagHamSparse

!> Upper triangle of the CSR matrix of initialize_sparse_matrix in the form PARDISO needs for a Hermitian matrix:
!> every row starts with its diagonal and has ascending columns. Only the (i, j >= i) entries are read, the way the
!> dense path's ZHEEV reads one triangle (the full matrix carries far-neighbour hoppings present in one direction only,
!> ~1e-5 eV). The diagonal must be present in every row (it is: with Diag.SparseUseShift it is always written).
subroutine sparse_upper_sorted(N, row_ptr, col_ind, values, ia, ja, a)
    implicit none
    integer, intent(in) :: N, row_ptr(:), col_ind(:)
    complex(dp), intent(in) :: values(:)
    integer, allocatable, intent(out) :: ia(:), ja(:)
    complex(dp), allocatable, intent(out) :: a(:)
    integer :: i, k, l, m, nnzU

    nnzU = 0
    do i = 1, N
       do k = row_ptr(i), row_ptr(i+1) - 1
          if (col_ind(k) >= i) nnzU = nnzU + 1
       end do
    end do
    allocate(ia(N+1), ja(nnzU), a(nnzU))
    l = 1
    do i = 1, N
       ia(i) = l
       do k = row_ptr(i), row_ptr(i+1) - 1
          if (col_ind(k) < i) cycle
          m = l   ! insertion into the sorted run ia(i) .. l-1 (rows hold ~100 entries)
          do while (m > ia(i))
             if (ja(m-1) <= col_ind(k)) exit
             ja(m) = ja(m-1)
             a(m) = a(m-1)
             m = m - 1
          end do
          ja(m) = col_ind(k)
          a(m) = values(k)
          l = l + 1
       end do
       if (l == ia(i)) then
          print *, 'sparse_upper_sorted: empty row ', i
          error stop 1
       else if (ja(ia(i)) /= i) then
          print *, 'sparse_upper_sorted: missing diagonal in row ', i
          error stop 1
       end if
       a(ia(i)) = real(a(ia(i)))
    end do
    ia(N+1) = l
end subroutine sparse_upper_sorted

!> Exact number of eigenvalues of the Hermitian matrix (ia, ja, a) below each of the nE energies dE (same units and
!> origin as the matrix: dE = 0 counts the negative eigenvalues), from the inertia of a sparse LDL^T (Sylvester).
!> PARDISO reports the inertia for REAL symmetric matrices only (mtype -2; for the complex Hermitian mtype -4 iparm(22)
!> and iparm(23) come back 0), so a complex H = A + iB is counted through its real symmetric embedding
!> [[A, -B], [B, A]] of size 2N, whose spectrum is that of H with every level doubled. A real H is counted directly.
!> ok = no perturbed pivot and no zero pivot at any energy; the counts are exact only then.
subroutine sparse_count_below(N, ia, ja, a, nE, dE, nbelow, ok)
    implicit none
    integer, intent(in) :: N, ia(:), ja(:), nE
    complex(dp), intent(in) :: a(:)
    real(dp), intent(in) :: dE(nE)
    integer, intent(out) :: nbelow(nE)
    logical, intent(out) :: ok
    integer, allocatable :: ib(:), jb(:), idg(:), lptr(:), lcol(:), lfill(:)
    real(dp), allocatable :: b(:), bw(:), lval(:)
    integer :: i, j, k, l, nU, nB, M, ie
    logical :: cplx
    integer*8 :: ptc(64)
    integer :: iparmc(64), maxfct, mnum, mtype, phase, nrhs, msglvl, error, idum(1)
    real(dp) :: ddum(1)

    nU = ia(N+1) - 1
    cplx = any(aimag(a(1:nU)) /= 0.0_dp)
    if (.not. cplx) then
       M = N
       allocate(ib(N+1), jb(nU), b(nU), idg(N))
       ib = ia(1:N+1)
       jb = ja(1:nU)
       b = real(a(1:nU))
       idg = ia(1:N)
    else
       ! strictly lower part of B = Im H by rows (B is antisymmetric): row j holds B(j,i) = -Im H(i,j), i < j
       allocate(lptr(N+1), lfill(N))
       lptr = 0
       nB = 0
       do i = 1, N
          do k = ia(i) + 1, ia(i+1) - 1
             if (aimag(a(k)) /= 0.0_dp) then
                lptr(ja(k)+1) = lptr(ja(k)+1) + 1
                nB = nB + 1
             end if
          end do
       end do
       lptr(1) = 1
       do i = 1, N
          lptr(i+1) = lptr(i) + lptr(i+1)
       end do
       allocate(lcol(max(nB,1)), lval(max(nB,1)))
       lfill = lptr(1:N)
       do i = 1, N    ! ascending i: every lower row comes out sorted
          do k = ia(i) + 1, ia(i+1) - 1
             if (aimag(a(k)) /= 0.0_dp) then
                j = ja(k)
                lcol(lfill(j)) = i
                lval(lfill(j)) = -aimag(a(k))
                lfill(j) = lfill(j) + 1
             end if
          end do
       end do
       ! upper triangle of [[A, -B], [B, A]]: row i = [A(i, j>=i) | -B(i, all j)], row N+i = [ . | A(i, j>=i)]
       M = 2*N
       allocate(ib(M+1), jb(2*nU + 2*nB), b(2*nU + 2*nB), idg(M))
       l = 1
       do i = 1, N
          ib(i) = l
          idg(i) = l
          do k = ia(i), ia(i+1) - 1
             jb(l) = ja(k)
             b(l) = real(a(k))
             l = l + 1
          end do
          do k = lptr(i), lptr(i+1) - 1      ! j < i: -B(i,j)
             jb(l) = N + lcol(k)
             b(l) = -lval(k)
             l = l + 1
          end do
          do k = ia(i) + 1, ia(i+1) - 1      ! j > i: -B(i,j) = -Im H(i,j)
             if (aimag(a(k)) /= 0.0_dp) then
                jb(l) = N + ja(k)
                b(l) = -aimag(a(k))
                l = l + 1
             end if
          end do
       end do
       do i = 1, N
          ib(N+i) = l
          idg(N+i) = l
          do k = ia(i), ia(i+1) - 1
             jb(l) = N + ja(k)
             b(l) = real(a(k))
             l = l + 1
          end do
       end do
       ib(M+1) = l
       deallocate(lptr, lfill, lcol, lval)
    end if

    allocate(bw(size(b)))
    ptc = 0
    iparmc = 0
    iparmc(1) = 1     ! no solver defaults
    iparmc(2) = 3     ! parallel nested dissection
    iparmc(10) = 8    ! pivot perturbation 1e-8, only if a pivot fails (reported in iparm(14))
    iparmc(11) = 0    ! no scaling and
    iparmc(13) = 0    ! no weighted matching: either would spoil the inertia
    iparmc(18) = -1
    iparmc(21) = 1    ! Bunch-Kaufman 1x1 and 2x2 pivots
    maxfct = 1
    mnum = 1
    nrhs = 1
    msglvl = 0
    mtype = -2
    ok = .true.
    phase = 11
    call pardiso(ptc, maxfct, mnum, mtype, phase, M, b, ib, jb, idum, nrhs, iparmc, msglvl, ddum, ddum, error)
    if (error /= 0) then
       print *, 'sparse_count_below: PARDISO analysis error ', error
       error stop 1
    end if
    do ie = 1, nE
       bw = b
       bw(idg) = bw(idg) - dE(ie)
       phase = 22
       call pardiso(ptc, maxfct, mnum, mtype, phase, M, bw, ib, jb, idum, nrhs, iparmc, msglvl, ddum, ddum, error)
       if (error /= 0) then
          print *, 'sparse_count_below: PARDISO factorisation error ', error
          error stop 1
       end if
       if (iparmc(14) /= 0 .or. iparmc(22) + iparmc(23) /= M) ok = .false.
       if (cplx) then
          if (mod(iparmc(23), 2) /= 0) ok = .false.
          nbelow(ie) = iparmc(23)/2
       else
          nbelow(ie) = iparmc(23)
       end if
    end do
    if (.not. ok) call MIO_Print('sparse_count_below: WARNING perturbed or zero pivots - the count is not exact; '// &
       'move the energy slightly','diag')
    phase = -1
    call pardiso(ptc, maxfct, mnum, mtype, phase, M, bw, ib, jb, idum, nrhs, iparmc, msglvl, ddum, ddum, error)
    deallocate(ib, jb, b, bw, idg)
end subroutine sparse_count_below

subroutine sparse_sort_real(n, x)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(inout) :: x(n)
    integer :: i, j
    real(dp) :: t
    do i = 2, n
       t = x(i)
       j = i - 1
       do while (j >= 1)
          if (x(j) <= t) exit
          x(j+1) = x(j)
          j = j - 1
       end do
       x(j+1) = t
    end do
end subroutine sparse_sort_real

!subroutine DiagH0TAPW(N, ns, is, ELoc, eigvec, KLoc, cell_real, H0, maxN, hopp, NList, Nneigh, neighCell,neig)
subroutine DiagH0TAPW(N, ns, is, ELoc, KLoc, cell_real, H0, maxN, hopp, NList, Nneigh, neighCell,neig, kpoint_index, evecOut)
    use constants, only : cmplx_i
    use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
    use scf, only : charge, Zch
    use atoms, only : Species, RAt, frac, AtomsSetCart, AtomsSetFrac, layerIndex
    use tbpar, only : U, g0
    use cell,                 only : rcell, ucell
    use constants,            only : pi
    use omp_lib,              only : omp_in_parallel
    implicit none
    integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is, neig
    integer, intent(in), optional :: kpoint_index
    logical :: tapw_pre_ok
    real(dp), intent(out) :: ELoc(N)
    ! Optional output of TAPW eigenvectors; unused unless a caller asks for it.
    complex(dp), intent(out), optional :: evecOut(:,:)
    real(dp), intent(in) :: KLoc(3), cell_real(3,3)
    real(dp), intent(inout) :: H0(N)
    complex(dp), intent(in) :: hopp(maxN,N)
    real(dp) :: tol
    real(dp) :: aGtemp
    real(dp) :: tapw_aG

    integer :: nuniq, ii, jj, code
                                    integer, allocatable :: label_raw(:), uniq_codes(:)

    integer :: i, j, in, info
    real(dp) :: R(3), zz, resid_norm

    ! ARPACK parameters and variables
    integer :: nev, ncv, lworkl, ido, ierr
    character(1) :: bmat
    character(2) :: which
    complex(dp), allocatable :: resid(:), v(:,:), workd(:), d(:), workev(:)
    complex(dp), allocatable :: workl(:)
    complex(dp), allocatable :: z(:,:)
    double precision, allocatable :: rwork(:)
    integer, allocatable :: iparam(:), ipntr(:)
    logical, allocatable :: select(:)
    integer :: max_iter, nn, iter
    logical :: useShift
    double precision :: shift

    integer, allocatable :: row_ptr(:), col_ind(:)
    complex(dp), allocatable :: values(:)
    real(dp), allocatable :: rand_real(:), rand_imag(:)
    complex(dp) :: sigma

    ! PARDISO variables

    INTEGER :: nnn, nnz

!C.. Internal solver memory pointer for 64-bit architectures
    INTEGER*8 pt(64)
!C.. Internal solver memory pointer for 32-bit architectures
!C.. INTEGER*4 pt(64)

!C.. This is OK in both cases
!    TYPE(MKL_PARDISO_HANDLE) pt(64)
!C.. All other variables
    INTEGER maxfct, mnum, mtype, phase, nrhs, msglvl
    INTEGER iparm(64)
    INTEGER, allocatable :: ia(:) ! row_ptr
    INTEGER, allocatable :: ja(:) ! col_ind
    INTEGER idum(1)
    complex(dp) ::  ddum(1)
    integer error
    complex(dp), allocatable :: a(:) ! values

    DATA nrhs /1/, maxfct /1/, mnum /1/

    complex(dp), allocatable :: EVectors(:,:)
    logical :: saveRitz, symmetric, useDenseMatrixTAPW, checkTAPWUnitary
    integer :: nconv
    complex(dp), allocatable :: ax(:)
    double precision, allocatable :: rd(:,:)
    double precision :: dznrm2, dlapy2

    ! === TAPW variables ===
    integer :: NG, Nlabel, M, M1
    integer, allocatable :: label(:)
    integer, allocatable :: original_label_codes(:)  ! Store original label codes before remapping
    real(dp), allocatable :: Gx(:), Gy(:)
    complex(dp), allocatable :: XArray(:,:)
    complex(dp), allocatable :: Hproj(:,:)
    real(dp), allocatable :: eigvals(:)
    complex(dp), allocatable :: ZWorkLoc(:)
    real(dp), allocatable :: DWorkLoc(:)
    real(dp), allocatable :: rigid_positions(:,:)  ! Rigid reference positions for X matrix
    complex(dp), allocatable :: Hproj_copy(:,:)  ! For Chern calculation storage
    integer :: lwork
    real(dp) :: temp_positions(3, N)  ! Temporary storage for position swapping
    logical :: original_frac_state     ! Original coordinate system state

    real(dp) :: k_ref(2)
    integer :: NGrange

    ! Variables for average mass term extraction
    integer :: label_A_idx, label_B_idx, label_B_hBN_idx, label_N_hBN_idx, i_atom
    integer :: count_A, count_B
    real(dp) :: V_A_0, V_B_0, V_B_hBN_0, V_N_hBN_0, m0, delta_avg, m0_meV, delta_avg_meV
    real(dp) :: m0_hBN, delta_avg_hBN, m0_hBN_meV, delta_avg_hBN_meV
    logical :: found_A, found_B, found_B_hBN, found_N_hBN

    ! Variables for gap extraction from eigenvalues
    real(dp) :: gap_from_eigenvalues, gap_from_eigenvalues_meV
    integer :: idx_min1, idx_min2
    real(dp), allocatable :: eigval_sorted(:)
    integer, allocatable :: idx_sorted(:)
    real(dp) :: temp_eig
    integer :: temp_idx

    ! Additional variables for K-point calculation
    real(dp) :: aG, refF(2), sGlattice(2,2), rG(2,2), refpoints1(2)
    real(dp) :: rotation_angle, cos_rot, sin_rot, det, sGlattice_inv(2,2)
    real(dp) :: k_ref_from_rcell(2)  ! For comparing with actual rcell
    real(dp) :: angle_rcell_1, angle_rcell_2, angle_rG_1, angle_rG_2
    real(dp) :: angle_diff_1, angle_diff_2
    logical :: useTriangularDistanceOrder  ! Flag to control distance-based ordering

    ! Read TAPW-specific input parameters outside OpenMP; inside OpenMP use cached config
    if (.not. omp_in_parallel()) then
       call MIO_InputParameter('TAPW.aG',tapw_aG,2.46019_dp)
    else
       tapw_aG = tapw_cfg_aG
    end if

    ! Parameters

    ! Variables
    !type(C_PTR) :: mkl_handle
    !!integer(int_kind) :: N, status, nnz, l, i, j, k, row_start, row_end, temp_index
    !!integer(int_kind), allocatable :: row_ptr(:), col_ind(:)
    !type(SPARSE_MATRIX_DESCR) :: descr

    nev=neig
    ncv=nev*10
    lworkl= 3*NCV**2 + 5*NCV

    max_iter = 10000
    allocate(resid(N), v(N, ncv), workd(3*N), workl(lworkl), rwork(ncv), d(nev+1), iparam(11), ipntr(14), select(ncv), rand_real(N), rand_imag(N))
    allocate(z(N,nev))
    allocate(workev(2*ncv))
    allocate(ax(N))
    allocate(rd(ncv,3))

    ! (The work arrays above are a leftover of the sparse solver; this routine does not call it. The
    ! checks of nev, ncv and of the uninitialised start vector that stood here stopped TAPW runs on
    ! cells of fewer than ten times Bands.SparseNeig atoms, and at random when the memory of the
    ! unset start vector happened to hold NaN.)

    ! Initialize arrays

    bmat = 'I'

    iparam(1) = 1
    iparam(3) = max_iter
    if (.not. omp_in_parallel()) then
       call MIO_InputParameter('Diag.SparseSetShift',shift,0.0_dp)
       call MIO_InputParameter('Diag.SparseUseShift',useShift,.false.)
       call MIO_InputParameter('Diag.SparseSaveRitz',saveRitz,.false.)
    else
       shift = tapw_cfg_shift
       useShift = tapw_cfg_useShift
       saveRitz = tapw_cfg_saveRitz
    end if
    if (useShift) then
       which = 'LM'
       iparam(7) = 3
       sigma = cmplx(shift, 0.0_dp)
       print*, "We will find eigenvalues close to ", shift
    else
       which = 'SM'
       sigma = (0.0_dp, 0.0_dp)
       iparam(7) = 1
    end if
    ! Initialize ARPACK parameters

    if (.not. omp_in_parallel()) then
       call MIO_InputParameter('Diag.SparseSetTol',tol,0.1_dp)
    else
       tol = tapw_cfg_tol
    end if

    ido = 0
    info = 0

    nn = N

    ! Initialize the starting vector resid with random values

    ! Check resid for initial state
    if (size(resid) /= N) then
        print *, 'Error: resid size mismatch: ', size(resid), ' expected: ', N
        error stop 1
    end if

    ! Debug prints

    ! === TAPW settings (configurable via input parameters) ===
    call MIO_InputParameter('useDenseMatrixTAPW',useDenseMatrixTAPW,.false.)
    call MIO_InputParameter('Diag.CheckTAPWUnitary',checkTAPWUnitary,.false.)

    ! Only initialize sparse matrix if NOT using dense matrix approach
    if (.not. useDenseMatrixTAPW) then
    ! Create CSR sparse matrix storage
    print *, "initialize the sparse matrix and put it in csr format"
    if (tapwDebug) then
       print *, "DEBUG TAPW: KLoc for sparse matrix =", KLoc
       print *, "DEBUG TAPW: |KLoc| =", sqrt(dot_product(KLoc, KLoc))
       print *, "DEBUG TAPW: cell_real(1,:) =", cell_real(1,:)
       print *, "DEBUG TAPW: cell_real(2,:) =", cell_real(2,:)

       ! Cell is already in correct column format for initialize_sparse_matrix
       ! matmul(cell, neighCell) expects cell(:,k) = k-th lattice vector
       print *, "DEBUG TAPW: Using cell_real in column format (correct for matmul):"
       print *, "TAPW:   a1 = cell_real(:,1) =", cell_real(:,1)
       print *, "TAPW:   a2 = cell_real(:,2) =", cell_real(:,2)
    end if

     call initialize_sparse_matrix(N, maxN, H0, hopp, NList, Nneigh, neighCell, ns, is, KLoc, cell_real, row_ptr, col_ind, values,sigma)
    print *, "done"
    symmetric = is_structurally_symmetric(values, row_ptr, col_ind, N)

    print*, "is it symmetric?", symmetric
    else
        ! Skip sparse matrix initialization when using dense matrix approach
        if (tapwDebug) print *, "Skipping sparse matrix initialization (using dense matrix approach)"
        ! Initialize dummy values to avoid uninitialized variables
        allocate(row_ptr(1), col_ind(1), values(1))
        row_ptr(1) = 1
        col_ind(1) = 1
        values(1) = cmplx(0.0_dp, 0.0_dp)
        symmetric = .true.
    end if

    ! Debug prints for CSR matrix
    if (any(values /= values)) then
        print *, 'Error: values contains NaN values after initialization.'
        error stop 1
    end if

    if (maxval(abs(values)) > 1e10) then
        print *, 'Warning: values contains extremely large values.'
    end if

    ! === TAPW settings (configurable via input parameters) ===
    NG = 20         ! for example (this gets overwritten by actual G-vector count)
    ! Nlabel will be computed dynamically from layer×sublattice combinations

    ! Calculate reference K-point using proper graphene lattice vectors
    ! Reference fractional coordinates depend on valley selection

    ! Graphene lattice constant (Angstroms) - use TAPW-specific value
    aG = tapw_aG  ! TAPW-specific graphene lattice constant from input

    ! Reference fractional coordinates in BZ (valley-dependent)
    if (useKprimeValley) then
       refF = (/ 1.0_dp/3.0_dp, 2.0_dp/3.0_dp /)  ! K' valley: [1/3, 2/3]
       call MIO_Print('TAPW targeting K'' valley with fractional coordinates [1/3, 2/3]','diag')
    else
       refF = (/ 2.0_dp/3.0_dp, 1.0_dp/3.0_dp /)  ! K valley: [2/3, 1/3] (default)
       if (tapwDebug) call MIO_Print('TAPW targeting K valley with fractional coordinates [2/3, 1/3]','diag')
    end if

    ! Build graphene lattice vectors (60-degree structure) - COLUMN storage like rcell
    ! sGlattice(:,1) = first vector, sGlattice(:,2) = second vector
    sGlattice(1,1) = 1.0_dp * aG                    ! First vector: (aG, 0)
    sGlattice(2,1) = 0.0_dp
    sGlattice(1,2) = cos(60.0_dp*pi/180.0_dp) * aG  ! Second vector: (aG*cos60°, aG*sin60°)
    sGlattice(2,2) = sin(60.0_dp*pi/180.0_dp) * aG

    ! Apply rotation if needed (configurable via Diag.MoireAngle)
    rotation_angle = moireAngle * pi / 180.0_dp
    cos_rot = cos(rotation_angle)
    sin_rot = sin(rotation_angle)

    ! Rotate lattice vectors (column storage: rotate each column)
    sGlattice = matmul(reshape((/ cos_rot, sin_rot, -sin_rot, cos_rot /), (/2,2/)), sGlattice)

    ! Convert to reciprocal space: rG = 2π * inverse(lattice.T)
    ! Manual 2x2 matrix inverse for sGlattice
    det = sGlattice(1,1)*sGlattice(2,2) - sGlattice(1,2)*sGlattice(2,1)
    if (abs(det) < 1e-10_dp) then
        print *, "Error: Lattice matrix is singular, det =", det
        error stop 1
    end if

    sGlattice_inv(1,1) =  sGlattice(2,2) / det
    sGlattice_inv(1,2) = -sGlattice(1,2) / det
    sGlattice_inv(2,1) = -sGlattice(2,1) / det
    sGlattice_inv(2,2) =  sGlattice(1,1) / det

    ! With column storage for sGlattice, rG should also use column storage like rcell
    rG = 2.0_dp * pi * transpose(sGlattice_inv)  ! Column storage: rG(:,i) = i-th reciprocal vector
    print *, "A*rG =", matmul(sGlattice, rG)  ! Should give 2π*identity with column storage

    ! Calculate reference K-point in reciprocal space using graphene vectors (by design)
    ! rG stores vectors as columns: rG(:,1) = b1, rG(:,2) = b2 (like rcell)
    k_ref = refF(1)*rG(:,1) + refF(2)*rG(:,2)

    ! Set NGrange based on calculation mode
    if (checkTAPWUnitary) then
        ! For unitary check: need complete first BZ coverage
        call MIO_Print('UNITARY CHECK MODE: Calculating NGrange for complete first BZ coverage','diag')
        call calculate_optimal_NGrange_for_complete_BZ(abs(physicalTwistAngle), tapwNG, NGrange)
    else
        ! For valley separation: use user-specified N_G for reduced basis around K-point
        NGrange = tapwNG
        if (tapwDebug) call MIO_Print('VALLEY SEPARATION MODE: Using user-specified N_G = '//trim(num2str(NGrange)),'diag')
        call MIO_Print('This creates a reduced TAPW basis centered around K-point','diag')
    end if

    ! Debug output for TAPW parameters (always show basic parameters)
    if (tapwDebug) call MIO_Print('TAPW Parameters:','diag')
    call MIO_Print('  N_G (requested G-vectors): '//trim(num2str(tapwNG)),'diag')
    call MIO_Print('  Moire angle: '//trim(num2str(moireAngle,6))//' degrees','diag')
    call MIO_Print('  TAPW graphene lattice constant: '//trim(num2str(tapw_aG,6))//' Angstroms','diag')
    if (useKprimeValley) then
       call MIO_Print('  Valley selection: K'' valley [1/3, 2/3]','diag')
    else
       call MIO_Print('  Valley selection: K valley [2/3, 1/3]','diag')
    end if
    if (useRigidPositions) then
       call MIO_Print('  Position mode: Rigid reference positions','diag')
    else
       call MIO_Print('  Position mode: Current relaxed positions','diag')
    end if

    if (tapwDebug) then
       print *, "  G-grid rotation angle:", gGridRotationAngle, "degrees"
       print *, "  Triangular truncation:", useTriangularTruncation
    end if

    if (tapwDebug) then
       ! Debug: Print lattice vectors and reference point calculation
       print *, "DEBUG: Lattice vectors and reference point:"
       print *, "  Graphene lattice constant aG =", aG
       print *, "  sGlattice(1,:) =", sGlattice(1,:)
       print *, "  sGlattice(2,:) =", sGlattice(2,:)
       print *, "  sGlattice_inv(1,:) =", sGlattice_inv(1,:)
       print *, "  sGlattice_inv(2,:) =", sGlattice_inv(2,:)
       print *, "  rG(1,:) =", rG(1,:)
       print *, "  rG(2,:) =", rG(2,:)
       print *, "  refF (fractional) =", refF
       print *, "  k_ref (Cartesian) =", k_ref
       print *, "  |k_ref| =", sqrt(k_ref(1)**2 + k_ref(2)**2)
    end if

    if (tapwDebug) then
       ! Debug: Compare orientations of graphene vs moiré lattice vectors
       print *, "DEBUG: Orientation comparison between graphene and moiré lattices:"
       print *, "  Moiré rcell(:,1) =", rcell(:,1)
       print *, "  Moiré rcell(:,2) =", rcell(:,2)
    end if
    if (tapwDebug) then
       print *, "  Graphene rG(:,1) =", rG(:,1)
       print *, "  Graphene rG(:,2) =", rG(:,2)

       ! Compare orientations (angles) rather than magnitudes
       ! Both rcell and rG now use column storage: (:,i) = i-th vector
       angle_rcell_1 = atan2(rcell(2,1), rcell(1,1)) * 180.0_dp / pi  ! rcell(:,1) = first vector
       angle_rcell_2 = atan2(rcell(2,2), rcell(1,2)) * 180.0_dp / pi  ! rcell(:,2) = second vector
       angle_rG_1 = atan2(rG(2,1), rG(1,1)) * 180.0_dp / pi          ! rG(:,1) = first vector
       angle_rG_2 = atan2(rG(2,2), rG(1,2)) * 180.0_dp / pi          ! rG(:,2) = second vector

       ! Handle angle wrapping
       angle_diff_1 = angle_rG_1 - angle_rcell_1
       angle_diff_2 = angle_rG_2 - angle_rcell_2
       if (angle_diff_1 > 180.0_dp) angle_diff_1 = angle_diff_1 - 360.0_dp
       if (angle_diff_1 < -180.0_dp) angle_diff_1 = angle_diff_1 + 360.0_dp
       if (angle_diff_2 > 180.0_dp) angle_diff_2 = angle_diff_2 - 360.0_dp
       if (angle_diff_2 < -180.0_dp) angle_diff_2 = angle_diff_2 + 360.0_dp

       print *, "  Orientation angles:"
    end if
    if (tapwDebug) then
       print *, "    Moiré b1 angle =", angle_rcell_1, "degrees"
       print *, "    Moiré b2 angle =", angle_rcell_2, "degrees"
       print *, "    Graphene b1 angle =", angle_rG_1, "degrees"
       print *, "    Graphene b2 angle =", angle_rG_2, "degrees"
       print *, "    Angle difference b1 =", angle_diff_1, "degrees"
       print *, "    Angle difference b2 =", angle_diff_2, "degrees"

       ! Keep k_ref from graphene calculation (by design)
       print *, "  k_ref from graphene vectors =", k_ref
       print *, "  rcell(2,:) =", rcell(2,:)
       print *, "  |rcell(1)| =", sqrt(rcell(1,1)**2 + rcell(1,2)**2)
       print *, "  |rcell(2)| =", sqrt(rcell(2,1)**2 + rcell(2,2)**2)
    end if

    ! Choose G-vector generation method based on input flag
    if (useTriangularTruncation) then
       call MIO_Print('Using triangular G-vector truncation (following Python get_Gvecs_tri)','diag')

       ! Control distance-based ordering with input parameter
       call MIO_InputParameter('TAPW.UseTriangularDistanceOrder', useTriangularDistanceOrder, .true.)
       if (useTriangularDistanceOrder) then
          call MIO_Print('TRIANGULAR: Distance-based ordering ENABLED (default)','diag')
       else
          call MIO_Print('TRIANGULAR: Distance-based ordering DISABLED - using grid order','diag')
       end if

       call generate_triangular_G_list(rcell, k_ref, NGrange, Gx, Gy, NG, rG, useTriangularDistanceOrder)
    else
       call MIO_Print('Using hexagonal shell G-vector truncation (original method)','diag')
       if (checkTAPWUnitary) then
           ! For unitary check, center G-grid around (0,0) to cover entire first BZ
           call MIO_Print('UNITARY CHECK MODE: Centering G-grid around Γ(0,0) for complete BZ coverage','diag')
           call generate_shifted_G_list_with_graphene_BZ(rcell, k_ref, NGrange, Gx, Gy, NG, [0.0_dp, 0.0_dp], rG)
       else
           ! Normal mode: center around K-point for local expansion (Python style)
           call generate_shifted_G_list_reduced(rcell, k_ref, NGrange, Gx, Gy, NG)
       end if
    end if

    ! Both valleys (opt-in, Diag.TAPWBothValleys): append the shells around the
    ! other graphene K point.  G 1..tapw_NGvalley1 stay the valley of refF.
    if (tapwBothValleys) &
       call tapw_append_second_valley(rcell, rG, sGlattice, refF, NGrange, useDenseMatrixTAPW, Gx, Gy, NG)

    ! Save G-vectors for Berry curvature calculation
#ifdef SEMICL
    if (calculateChern .or. berryFluxTAPW) then
#else
    if (calculateChern) then
#endif
       if (allocated(saved_Gx)) deallocate(saved_Gx)
       if (allocated(saved_Gy)) deallocate(saved_Gy)
       allocate(saved_Gx(NG), saved_Gy(NG))
       saved_Gx = Gx
       saved_Gy = Gy
       ! print once: this routine runs at every k-point of the grid
       if (.not. bf_announced) then
          call MIO_Print('Saved '//trim(num2str(NG))//' G-vectors for Berry curvature calculation','diag')
          bf_announced = .true.
       end if
       saved_NG = NG
    end if

    ! Construct proper TAPW labels: combine layer and sublattice information
    allocate(label(N))
    allocate(original_label_codes(N))
    do i = 1, N
       label(i) = 10 * layerIndex(i) + Species(i)  ! e.g., layer 1 sublattice A -> 11, layer 2 sublattice B -> 22
       original_label_codes(i) = label(i)  ! Store original codes
    end do

    ! Compute number of unique labels dynamically
    call compute_unique_labels(label, N, Nlabel)

    ! Remap labels to contiguous indices 1..Nlabel
    call remap_labels_to_contiguous(label, N, Nlabel)

    !! New (one call does all, deterministically):

    !! ensure label(:) exists
    !
    !! raw two-digit codes (e.g., 11,12,21,22)
    !
    !! collect uniques (unsorted first)
    !   ! is 'code' already in uniq_codes(1:nuniq)?
    !
    !! insertion sort uniq_codes(1:nuniq) ascending (no external routine)
    !
    !! map raw -> rank in sorted uniques => labels in 1..Nlabel
    !
    !
    !! guards before build_X
    !
    ! Now label(i) ∈ {1..Nlabel}, stable order, ready for build_X

    ! Save Nlabel for Berry curvature calculation
#ifdef SEMICL
    if (calculateChern .or. berryFluxTAPW) then
#else
    if (calculateChern) then
#endif
       saved_Nlabel = Nlabel
    end if

    M = NG * Nlabel
    M_tapw = M  ! Update M_tapw to actual projected matrix size

    if (tapwDebug) then
       ! Debug: Show TAPW structure
       print *, "TAPW Structure:"
       print *, "  Number of atoms (N):", N
       print *, "  Number of G-vectors (NG):", NG
       print *, "  Number of unique labels (Nlabel):", Nlabel
       print *, "  Projected space dimension (M = NG × Nlabel):", M
       print *, "  Updated M_tapw =", M_tapw, "(now matches actual matrix size)"
    else
       ! Always show basic TAPW dimensions
       if (tapwDebug) call MIO_Print('TAPW matrix: '//trim(num2str(N))//' atoms → '//trim(num2str(M))//' projected states (NG='//trim(num2str(NG))//', Nlabel='//trim(num2str(Nlabel))//')','diag')
    end if

    allocate(XArray(N, M))
    if (frac) call AtomsSetCart()

    if (tapwDebug) then
       ! Debug: Print atomic coordinates and G-vector information
       print *, "DEBUG: Atomic coordinates and X matrix construction:"
       print *, "  Total atoms N =", N
       print *, "  frac flag =", frac
       print *, "  Sample atomic coordinates (first 5 atoms):"
       do i = 1, min(5, N)
          print '(A,I3,A,2F12.6,A,I0)', "    Atom ", i, ": (", Rat(1,i), Rat(2,i), "), label=", label(i)
       end do
       print *, "  Coordinate ranges: x=[", minval(Rat(1,1:N)), ",", maxval(Rat(1,1:N)), "]"
       print *, "  Coordinate ranges: y=[", minval(Rat(2,1:N)), ",", maxval(Rat(2,1:N)), "]"
       print *, "  G-vector ranges: Gx=[", minval(Gx), ",", maxval(Gx), "]"
       print *, "  G-vector ranges: Gy=[", minval(Gy), ",", maxval(Gy), "]"
    end if

    if (tapwDebug) then
       ! Debug: Show phase calculations for multiple atoms and G-vectors
       print *, "DEBUG: Phase calculations G·r for X-matrix construction:"
       block
          real(dp) :: real_phase
          do i = 1, min(3, N)  ! First 3 atoms
             do j = 1, min(3, NG)  ! First 3 G-vectors
                real_phase = Gx(j)*Rat(1,i) + Gy(j)*Rat(2,i)
                print '(A,I0,A,I0,A,F12.6,A,F12.6,A,F12.6,A,2F12.6)', &
                      "    Atom ", i, ", G-vector ", j, ": G=[", Gx(j), ",", Gy(j), &
                      "], G·r=", real_phase, ", exp(iG·r)=", cos(real_phase), sin(real_phase)
             end do
          end do
       end block
    end if

    if (tapwDebug) then
       ! Output atomic coordinates for visualization
       open(unit=98, file='atomic_coords_debug.dat', status='replace')
       write(98, '(A)') '# x, y, label (for matplotlib plotting)'
       do i = 1, N
          write(98, '(2F16.8,I6)') Rat(1,i), Rat(2,i), label(i)
       end do
       close(98)
    end if
    if (tapwDebug) then
       print *, "DEBUG: Atomic coordinates written to atomic_coords_debug.dat for visualization"

       ! Output comprehensive debugging data for TAPW visualization
       print *, "DEBUG: About to call output_brillouin_zones_debug..."
       call output_brillouin_zones_debug(aG, moireAngle, k_ref, rG)

       ! Output moiré BZ data
       print *, "DEBUG: Outputting moiré BZ data..."
       call output_moire_bz_debug(rcell)

       ! Output k-path data if available
       print *, "DEBUG: Outputting k-path data..."
       call output_kpath_debug()
    end if

    if (useRigidPositions) then
       ! Use rigid reference positions for X matrix construction
       ! This ensures perfect TAPW unitarity even with lattice reconstruction
       call MIO_Print('Using rigid reference positions for TAPW X matrix construction','diag')

       ! Read rigid positions from generateInit.xyz
       allocate(rigid_positions(3, N))
       call read_rigid_positions_for_tapw('generateInit.xyz', N, rigid_positions)

       ! Wrap rigid positions to single moiré unit cell (same as before, but for rigid positions)
       ! Note: We need to temporarily modify the coordinate system for wrapping
       ! Store current positions and coordinate system state
       temp_positions = Rat  ! Save current (relaxed) positions
       original_frac_state = frac

       ! Replace current positions with rigid positions for wrapping
       Rat = rigid_positions
    else
       ! Use current positions in Rat directly (default behavior)
       call MIO_Print('Using current atomic positions for TAPW X matrix construction','diag')

       ! Store current state for consistency with rigid position path
       temp_positions = Rat  ! Save current positions (no change needed)
       original_frac_state = frac
    end if

    if (tapwDebug) then
       if (useRigidPositions) then
          if (tapwDebug) call MIO_Print('DEBUG: First 3 atomic positions BEFORE build_X (rigid):','diag')
       else
          if (tapwDebug) call MIO_Print('DEBUG: First 3 atomic positions BEFORE build_X (current):','diag')
       end if
       do i = 1, min(3, N)
           call MIO_Print('  Atom '//trim(num2str(i))//': ['//trim(num2str(Rat(1,i),6))//','//&
                         trim(num2str(Rat(2,i),6))//','//trim(num2str(Rat(3,i),6))//']','diag')
       end do
    end if

    if (useRigidPositions) then
       ! Convert rigid positions to fractional coordinates relative to moiré cell
       if (.not. frac) call AtomsSetFrac()

       ! Wrap fractional coordinates to [0,1) to ensure single unit cell
       do i = 1, N
          Rat(1,i) = Rat(1,i) - floor(Rat(1,i))
          Rat(2,i) = Rat(2,i) - floor(Rat(2,i))
          ! Don't wrap z-coordinate (only x,y are periodic in 2D moiré)
       end do

       ! Convert wrapped rigid positions back to Cartesian for build_X
       call AtomsSetCart()
       call MIO_Print('Rigid positions wrapped and converted to Cartesian for X matrix','diag')
    else
       ! For current positions, ensure we're in the right coordinate system for build_X
       if (frac) call AtomsSetCart()
       call MIO_Print('Using current positions directly for X matrix construction','diag')
    end if

    ! Build X matrix using current positions (rigid or relaxed)
    call build_X(XArray, Rat(1,:), Rat(2,:), label, Gx, Gy, N, NG, Nlabel)
    if (tapwLowdin) call lowdin_orthonormalize_X(XArray, N, M)

    ! RESTORE original atomic positions for TB calculations (if they were modified)
    if (useRigidPositions) then
       Rat = temp_positions
       if (original_frac_state .and. .not. frac) call AtomsSetFrac()
       if (.not. original_frac_state .and. frac) call AtomsSetCart()
       call MIO_Print('Restored original relaxed positions for TB Hamiltonian construction','diag')

       if (tapwDebug) then
          if (tapwDebug) call MIO_Print('DEBUG: First 3 atomic positions AFTER restoration (relaxed):','diag')
          do i = 1, min(3, N)
              call MIO_Print('  Atom '//trim(num2str(i))//': ['//trim(num2str(Rat(1,i),6))//','//&
                            trim(num2str(Rat(2,i),6))//','//trim(num2str(Rat(3,i),6))//']','diag')
          end do
       end if

       ! Clean up rigid positions array
       deallocate(rigid_positions)
    else
       call MIO_Print('Using same positions for both X matrix and TB Hamiltonian construction','diag')
    end if

    if (tapwDebug) then
       ! Comprehensive numerical verification of TAPW setup (without ucell_T - will calculate internally)
       call verify_tapw_numerical_consistency(rcell, rG, k_ref, gGridRotationAngle, &
                                             Gx, Gy, NG, aG)
    end if

    if (tapwDebug) then
       ! Verify X-matrix column norms
       call verify_x_matrix_norms(XArray, N, NG)
    end if

    ! PERFORMANCE OPTIMIZATION: Use pre-allocated arrays if available
    ! size() must not be evaluated for an unallocated array, and .and. does not short-circuit
    tapw_pre_ok = .false.
    if (allocated(tapw_Hproj)) tapw_pre_ok = (size(tapw_Hproj,1) >= M .and. size(tapw_Hproj,2) >= M)
    if (tapw_pre_ok) then
        ! Copy from pre-allocated array (no pointer overhead)
    allocate(Hproj(M, M))
        Hproj(1:M, 1:M) = tapw_Hproj(1:M, 1:M)
        call MIO_Print('Using pre-allocated Hproj array for performance','diag')
    else
        allocate(Hproj(M, M))
        call MIO_Print('Allocated new Hproj array (not pre-allocated)','diag')
    end if

    ! Branch: use dense or sparse matrix approach for TAPW transformation
    if (useDenseMatrixTAPW) then
        call MIO_Print('Using DENSE matrix approach for TAPW transformation','diag')
        call transform_dense_hamiltonian_tapw(N, M, XArray, Hproj, KLoc, cell_real, H0, maxN, hopp, NList, Nneigh, neighCell, ns, is)
    else
        call MIO_Print('Using SPARSE matrix approach for TAPW transformation','diag')
        ! For large systems (>1M atoms), force sparse approach to avoid memory issues
        if (N > 1000000) then
            call MIO_Print('Large system detected (N='//trim(num2str(N))//'), using sparse-only TAPW transformation','diag')
            call transform_sparse_hamiltonian(N, M, row_ptr, col_ind, values, XArray, Hproj)
        else
            call transform_sparse_hamiltonian(N, M, row_ptr, col_ind, values, XArray, Hproj)
        end if
        if (.not. allocated(values)) stop "values not allocated after transform"
    end if

    ! Intervalley block of the two-valley basis: its size is the valley mixing
    ! carried by the lattice; Diag.TAPWValleyDecouple removes it (known-answer test).
    if (tapwBothValleys) then
       M1 = tapw_NGvalley1*Nlabel
       call MIO_Print('TAPW valleys: k = '//trim(num2str(KLoc(1),6))//' '//trim(num2str(KLoc(2),6))// &
            '   max|H(K,K'')| = '//trim(num2str(maxval(abs(Hproj(1:M1,M1+1:M)))*g0*1000.0_dp,6))//' meV','diag')
       if (tapwValleyDecouple) then
          Hproj(1:M1,M1+1:M) = cmplx(0.0_dp, 0.0_dp, kind=dp)
          Hproj(M1+1:M,1:M1) = cmplx(0.0_dp, 0.0_dp, kind=dp)
       end if
    end if

    if (tapwDebug) then
       ! Verify projected Hamiltonian Hermiticity
       call verify_projected_hamiltonian_hermiticity(Hproj, M_tapw)
    end if

    ! NUMERICAL CHECK: Verify TAPW unitary transformation property (if enabled)
    if (checkTAPWUnitary .and. tapwDebug) then
        call MIO_Print('Running TAPW unitary transformation check (enabled via Diag.CheckTAPWUnitary)','diag')
        call MIO_Print('Current k-point for unitary check: ['//trim(num2str(KLoc(1),6))//','//trim(num2str(KLoc(2),6))//']','diag')
        call diagnose_tapw_phase_consistency(KLoc, Gx, Gy, NG, XArray, N, Nlabel)
        call verify_tapw_unitary_transformation(N, M, row_ptr, col_ind, values, XArray, Hproj, NG, Nlabel, Gx, Gy)
    else if (checkTAPWUnitary .and. .not. tapwDebug) then
        call MIO_Print('TAPW unitary check enabled but debug output suppressed (set Diag.TAPWDebug=.true. for details)','diag')
    else if (tapwDebug) then
        call MIO_Print('TAPW unitary check disabled (set Diag.CheckTAPWUnitary=.true. to enable)','diag')
    end if

    ! PERFORMANCE OPTIMIZATION: Use pre-allocated arrays if available
    tapw_pre_ok = .false.
    if (allocated(tapw_eigvals)) tapw_pre_ok = (size(tapw_eigvals) >= M)
    if (tapw_pre_ok) then
    allocate(eigvals(M))
        eigvals(1:M) = tapw_eigvals(1:M)
        call MIO_Print('Using pre-allocated eigvals array for performance','diag')
    else
        allocate(eigvals(M))
        call MIO_Print('Allocated new eigvals array (not pre-allocated)','diag')
    end if

    lwork = 2*M
    tapw_pre_ok = .false.
    if (allocated(tapw_ZWorkLoc)) tapw_pre_ok = (size(tapw_ZWorkLoc) >= 2*M)
    if (tapw_pre_ok) then
        allocate(ZWorkLoc(2*M))
        ZWorkLoc(1:2*M) = tapw_ZWorkLoc(1:2*M)
        call MIO_Print('Using pre-allocated ZWorkLoc array for performance','diag')
    else
        allocate(ZWorkLoc(2*M))
        call MIO_Print('Allocated new ZWorkLoc array (not pre-allocated)','diag')
    end if

    tapw_pre_ok = .false.
    if (allocated(tapw_DWorkLoc)) tapw_pre_ok = (size(tapw_DWorkLoc) >= 3*M)
    if (tapw_pre_ok) then
        allocate(DWorkLoc(3*M))
        DWorkLoc(1:3*M) = tapw_DWorkLoc(1:3*M)
        call MIO_Print('Using pre-allocated DWorkLoc array for performance','diag')
    else
        allocate(DWorkLoc(3*M))
        call MIO_Print('Allocated new DWorkLoc array (not pre-allocated)','diag')
    end if

    ! Extract average mass term from G=0 block BEFORE ZHEEV overwrites Hproj
    !
    ! TAPW basis structure: |G, label⟩ where label combines layerIndex and Species
    ! - Original label code: 10*layerIndex + Species
    ! - For GBNtwoLayers: layer 1 has Species 1 (A sublattice) and 2 (B sublattice)
    ! - Original codes: 11 (layer1-A), 12 (layer1-B), 23 (layer2-B), 24 (layer2-N)
    ! - After remapping: these become contiguous indices 1..Nlabel
    !
    ! TAPW basis ordering (from build_X): col = (i_G - 1) * N_label + i_label
    ! For NG=1: columns are 1, 2, 3, 4 corresponding to |G=0, label=1⟩, |G=0, label=2⟩, etc.
    ! Hproj(i,j) is matrix element between |G=0, label=i⟩ and |G=0, label=j⟩
    !
    ! To extract G=0 mass term: find which TAPW basis index corresponds to layer 1 sublattice A and B
    ! Also extract hBN layer gap for verification
    if (NG == 1 .and. allocated(original_label_codes)) then
       ! Find which remapped label indices correspond to:
       ! - Layer 1 (graphene): Species 1 = A sublattice (code 11), Species 2 = B sublattice (code 12)
       ! - Layer 2 (hBN): Species 3 = B atom (code 23), Species 4 = N atom (code 24)
       found_A = .false.
       found_B = .false.
       found_B_hBN = .false.
       found_N_hBN = .false.
       do i_atom = 1, N
          ! Original code 11 = layer 1 (graphene), Species 1 (A sublattice)
          if (original_label_codes(i_atom) == 11 .and. .not. found_A) then
             label_A_idx = label(i_atom)  ! Get remapped label index (1..Nlabel)
             found_A = .true.
          end if
          ! Original code 12 = layer 1 (graphene), Species 2 (B sublattice)
          if (original_label_codes(i_atom) == 12 .and. .not. found_B) then
             label_B_idx = label(i_atom)  ! Get remapped label index (1..Nlabel)
             found_B = .true.
          end if
          ! Original code 23 = layer 2 (hBN), Species 3 (B atom)
          if (original_label_codes(i_atom) == 23 .and. .not. found_B_hBN) then
             label_B_hBN_idx = label(i_atom)  ! Get remapped label index (1..Nlabel)
             found_B_hBN = .true.
          end if
          ! Original code 24 = layer 2 (hBN), Species 4 (N atom)
          if (original_label_codes(i_atom) == 24 .and. .not. found_N_hBN) then
             label_N_hBN_idx = label(i_atom)  ! Get remapped label index (1..Nlabel)
             found_N_hBN = .true.
          end if
          if (found_A .and. found_B .and. found_B_hBN .and. found_N_hBN) exit
       end do

       ! Extract diagonal elements from Hproj (G=0 block)
       ! Hproj is M×M where M = NG * Nlabel = 1 * Nlabel = Nlabel
       ! For NG=1, the basis is |G=0, label⟩, so Hproj(label, label) gives the G=0 block

       ! Extract graphene layer 1 gap (A and B sublattices)
       ! NOTE: Hproj diagonal elements are normalized averages due to TAPW basis normalization
       ! In build_X, each column is normalized by 1/sqrt(N_atoms_with_label)
       ! So Hproj(label, label) = (1/N_label) * sum_{i,j with label} H(i,j)
       ! This gives an average over all atoms with that label, which may differ from
       ! the direct on-site energy if there are spatial variations in the moiré pattern
       !
       ! IMPORTANT: The normalization factor is the same for rigid and relaxed systems,
       ! but the spatial averaging can give different results:
       ! - Rigid system: atoms with same label have similar on-site energies → average ≈ individual value
       ! - Relaxed system: atoms with same label have varying on-site energies due to moiré pattern
       !   → average can be much smaller than typical individual value if there are regions
       !   with opposite contributions (some positive, some negative mass term contributions)
       ! - The gap in bands (from eigenvalues) may still look reasonable because it includes
       !   off-diagonal contributions and the full structure of Hproj, not just G=0 diagonal
       if (found_A .and. found_B .and. label_A_idx <= M .and. label_B_idx <= M) then
          V_A_0 = real(Hproj(label_A_idx, label_A_idx))
          V_B_0 = real(Hproj(label_B_idx, label_B_idx))
          m0 = (V_A_0 - V_B_0) / 2.0_dp
          delta_avg = 2.0_dp * abs(m0)

          ! Convert to meV: energies are in units of g0, so multiply by g0 (eV) and 1000 (meV/eV)
          m0_meV = m0 * g0 * 1000.0_dp
          delta_avg_meV = delta_avg * g0 * 1000.0_dp

          ! Count atoms with each label for diagnostic purposes
          count_A = 0
          count_B = 0
          do i_atom = 1, N
             if (original_label_codes(i_atom) == 11) count_A = count_A + 1
             if (original_label_codes(i_atom) == 12) count_B = count_B + 1
          end do

          ! Print graphene average mass term for this k-point
          if (present(kpoint_index)) then
             print *, "TAPW Average Mass Term - GRAPHENE Layer 1 (k-point ", kpoint_index, "):"
          else
             print *, "TAPW Average Mass Term - GRAPHENE Layer 1:"
          end if
          print *, "  TAPW basis: label_A_idx = ", label_A_idx, " (layer 1, Species 1 = A, GRAPHENE)"
          print *, "  TAPW basis: label_B_idx = ", label_B_idx, " (layer 1, Species 2 = B, GRAPHENE)"
          print *, "  Number of atoms with label A (code 11): ", count_A
          print *, "  Number of atoms with label B (code 12): ", count_B
          print *, "  TAPW normalization factor for A: 1/sqrt(", count_A, ") = ", 1.0_dp/sqrt(real(count_A, dp))
          print *, "  TAPW normalization factor for B: 1/sqrt(", count_B, ") = ", 1.0_dp/sqrt(real(count_B, dp))
          print *, "  V_A(G=0) = ", V_A_0, " (in units of g0)"
          print *, "  V_B(G=0) = ", V_B_0, " (in units of g0)"
          print *, "  m0 = (V_A - V_B)/2 = ", m0, " (in units of g0) = ", m0_meV, " meV"
          print *, "  Delta_avg = 2|m0| = ", delta_avg, " (in units of g0) = ", delta_avg_meV, " meV"
          print *, "  NOTE: V_A and V_B are TAPW-averaged values (normalized by 1/N_label)"
          print *, "        For relaxed systems, spatial variations in moiré pattern can cause"
          print *, "        on-site energies to vary across atoms with same label, leading to"
          print *, "        smaller averaged mass term than individual atom values"
          print *, "        The gap in bands (from eigenvalues) may still be reasonable because"
          print *, "        it includes off-diagonal contributions and full Hproj structure"
       else
          if (.not. found_A) print *, "Warning: Could not find layer 1 sublattice A (code 11, GRAPHENE) in TAPW basis"
          if (.not. found_B) print *, "Warning: Could not find layer 1 sublattice B (code 12, GRAPHENE) in TAPW basis"
          if (found_A .and. label_A_idx > M) print *, "Warning: label_A_idx = ", label_A_idx, " > M = ", M
          if (found_B .and. label_B_idx > M) print *, "Warning: label_B_idx = ", label_B_idx, " > M = ", M
       end if

       ! Extract hBN layer 2 gap (B and N atoms) for verification
       if (found_B_hBN .and. found_N_hBN .and. label_B_hBN_idx <= M .and. label_N_hBN_idx <= M) then
          V_B_hBN_0 = real(Hproj(label_B_hBN_idx, label_B_hBN_idx))
          V_N_hBN_0 = real(Hproj(label_N_hBN_idx, label_N_hBN_idx))
          m0_hBN = (V_B_hBN_0 - V_N_hBN_0) / 2.0_dp
          delta_avg_hBN = 2.0_dp * abs(m0_hBN)

          ! Convert to meV
          m0_hBN_meV = m0_hBN * g0 * 1000.0_dp
          delta_avg_hBN_meV = delta_avg_hBN * g0 * 1000.0_dp

          ! Print hBN average mass term for this k-point
          if (present(kpoint_index)) then
             print *, "TAPW Average Mass Term - hBN Layer 2 (k-point ", kpoint_index, "):"
          else
             print *, "TAPW Average Mass Term - hBN Layer 2:"
          end if
          print *, "  TAPW basis: label_B_hBN_idx = ", label_B_hBN_idx, " (layer 2, Species 3 = B, hBN)"
          print *, "  TAPW basis: label_N_hBN_idx = ", label_N_hBN_idx, " (layer 2, Species 4 = N, hBN)"
          print *, "  V_B_hBN(G=0) = ", V_B_hBN_0, " (in units of g0)"
          print *, "  V_N_hBN(G=0) = ", V_N_hBN_0, " (in units of g0)"
          print *, "  m0_hBN = (V_B_hBN - V_N_hBN)/2 = ", m0_hBN, " (in units of g0) = ", m0_hBN_meV, " meV"
          print *, "  Delta_avg_hBN = 2|m0_hBN| = ", delta_avg_hBN, " (in units of g0) = ", delta_avg_hBN_meV, " meV"
       else
          if (.not. found_B_hBN) print *, "Warning: Could not find layer 2 B atom (code 23, hBN) in TAPW basis"
          if (.not. found_N_hBN) print *, "Warning: Could not find layer 2 N atom (code 24, hBN) in TAPW basis"
          if (found_B_hBN .and. label_B_hBN_idx > M) print *, "Warning: label_B_hBN_idx = ", label_B_hBN_idx, " > M = ", M
          if (found_N_hBN .and. label_N_hBN_idx > M) print *, "Warning: label_N_hBN_idx = ", label_N_hBN_idx, " > M = ", M
       end if
    end if

    ! Projected Hamiltonian before ZHEEV overwrites it (opt-in diagnostic; default .false. changes nothing).
    if (gWeightsTAPW .and. gWeightsHam .and. is == 1 .and. present(kpoint_index)) &
         call GWeightsHamWrite(M, NG, Nlabel, N, Gx, Gy, Hproj, label, original_label_codes, kpoint_index, KLoc)

    ! Save Hamiltonian copy for Chern calculation before ZHEEV overwrites it
    if (calculateChern .and. present(kpoint_index)) then
       allocate(Hproj_copy(M, M))
       Hproj_copy = Hproj
    end if

    call ZHEEV('V','U', M, Hproj, M, eigvals, ZWorkLoc, lwork, DWorkLoc, info)

    if (info /= 0) then
       print *, "Diagonalization failed: ZHEEV info =", info
       error stop 1
    end if

    ! Extract gap from eigenvalues for comparison with diagonal extraction
    ! This is the actual gap that appears in the bands
    if (NG == 1 .and. M >= 2) then
       ! Find the two lowest eigenvalues (should correspond to graphene Dirac cone)
       ! Allocate arrays for sorting
       allocate(eigval_sorted(M))
       allocate(idx_sorted(M))

       ! Sort eigenvalues to find the two lowest
       eigval_sorted = eigvals
       do i = 1, M
          idx_sorted(i) = i
       end do
       ! Simple bubble sort for small M (NG=1 means M is small, typically 4)
       do i = 1, M-1
          do j = 1, M-i
             if (eigval_sorted(j) > eigval_sorted(j+1)) then
                ! Swap
                temp_eig = eigval_sorted(j)
                temp_idx = idx_sorted(j)
                eigval_sorted(j) = eigval_sorted(j+1)
                idx_sorted(j) = idx_sorted(j+1)
                eigval_sorted(j+1) = temp_eig
                idx_sorted(j+1) = temp_idx
             end if
          end do
       end do

       idx_min1 = idx_sorted(1)
       idx_min2 = idx_sorted(2)
       gap_from_eigenvalues = eigval_sorted(2) - eigval_sorted(1)
       gap_from_eigenvalues_meV = gap_from_eigenvalues * g0 * 1000.0_dp

       if (present(kpoint_index)) then
          print *, "TAPW Gap from Eigenvalues - GRAPHENE Layer 1 (k-point ", kpoint_index, "):"
       else
          print *, "TAPW Gap from Eigenvalues - GRAPHENE Layer 1:"
       end if
       print *, "  Lowest eigenvalue (idx ", idx_min1, "): ", eigval_sorted(1), " (in units of g0)"
       print *, "  2nd lowest eigenvalue (idx ", idx_min2, "): ", eigval_sorted(2), " (in units of g0)"
       print *, "  Gap from eigenvalues = ", gap_from_eigenvalues, " (in units of g0) = ", gap_from_eigenvalues_meV, " meV"
       print *, "  NOTE: This is the actual gap that appears in the bands"
       print *, "        It includes off-diagonal contributions from Hproj"
       print *, "        Compare with Delta_avg from diagonal extraction above"

       ! Deallocate sorting arrays
       deallocate(eigval_sorted, idx_sorted)
    end if

    ! Store eigenvectors and Hamiltonians for Chern calculation (if needed)
    if (calculateChern .and. present(kpoint_index)) then
       ! Initialize storage arrays if this is the first k-point
       if (.not. allocated(stored_eigenvectors)) then
          stored_M = M
          stored_N = N
          stored_nspin = ns  ! Store number of spin channels
          ! Calculate actual number of k-points needed based on Chern grid
          stored_nk = nk_chern_x * nk_chern_y
          call MIO_Print('Allocating TAPW storage for '//trim(num2str(stored_nk))//' k-points ('//trim(num2str(nk_chern_x))//'x'//trim(num2str(nk_chern_y))//' grid) for '//trim(num2str(stored_nspin))//' spin channel(s)','diag')
          allocate(stored_eigenvectors(M, M, stored_nk, stored_nspin))
          allocate(stored_hamiltonians(M, M, stored_nk, stored_nspin))

          ! Store X matrix for Berry curvature calculation (TB parameters are module-level)
          allocate(stored_X_matrix(N, M))
          stored_X_matrix = XArray

          call MIO_Print('Allocated TAPW storage for Chern calculation: '//trim(num2str(M))//'x'//trim(num2str(M))//'x'//trim(num2str(stored_nk))//'x'//trim(num2str(stored_nspin)),'diag')
          call MIO_Print('Stored X matrix ('//trim(num2str(N))//'x'//trim(num2str(M))//') for Berry curvature calculation','diag')

          ! VERIFICATION: Check stored data integrity
          call MIO_Print('=== TAPW STORAGE VERIFICATION ===','diag')
          call MIO_Print('Stored dimensions: N='//trim(num2str(stored_N))//', M='//trim(num2str(stored_M))//', nk='//trim(num2str(stored_nk))//', nspin='//trim(num2str(stored_nspin)),'diag')
          call MIO_Print('XArray dimensions: '//trim(num2str(size(XArray,1)))//'x'//trim(num2str(size(XArray,2))),'diag')
          call MIO_Print('XArray norm: '//trim(num2str(sqrt(sum(abs(XArray)**2)),8)),'diag')
          call MIO_Print('XArray max element: '//trim(num2str(maxval(abs(XArray)),8)),'diag')
          call MIO_Print('XArray min element: '//trim(num2str(minval(abs(XArray)),8)),'diag')
       else
          ! Check for consistency - M should not change between k-points
          if (M /= stored_M) then
             call MIO_Print('ERROR: TAPW matrix size changed between k-points','diag')
             call MIO_Print('  Initial stored_M: '//trim(num2str(stored_M)),'diag')
             call MIO_Print('  Current M: '//trim(num2str(M)),'diag')
             call MIO_Print('  This indicates inconsistent TAPW calculation','diag')
             ! Use the smaller dimension to avoid bounds errors
             M = min(M, stored_M)
             call MIO_Print('  Using M = '//trim(num2str(M))//' (minimum of both)','diag')
          end if
       end if

       ! Store this k-point's results with bounds checking
       if (kpoint_index <= stored_nk .and. is <= stored_nspin) then
          if (M <= size(stored_eigenvectors,1) .and. M <= size(stored_eigenvectors,2)) then
             stored_eigenvectors(1:M,1:M,kpoint_index,is) = Hproj(1:M,1:M)  ! ZHEEV returns eigenvectors in Hproj
             stored_hamiltonians(1:M,1:M,kpoint_index,is) = Hproj_copy(1:M,1:M)  ! Store original Hamiltonian

             ! CRITICAL FIX: Store TAPW eigenvalues to match eigenvectors
             if (.not. allocated(stored_eigenvalues)) then
                allocate(stored_eigenvalues(M, stored_nk, stored_nspin))
                call MIO_Print('Allocated stored_eigenvalues array: '//trim(num2str(M))//'x'//trim(num2str(stored_nk))//'x'//trim(num2str(stored_nspin)),'diag')
             end if
             stored_eigenvalues(1:M,kpoint_index,is) = eigvals(1:M)  ! Store TAPW eigenvalues

             ! DEBUG: Verify eigenvalue storage for first k-point
             if (kpoint_index == 1 .and. is == 1) then
                call MIO_Print('=== TAPW EIGENVALUE STORAGE VERIFICATION (k-point 1, spin 1) ===','diag')
                call MIO_Print('Stored eigenvalues dimensions: '//trim(num2str(size(stored_eigenvalues,1)))//'x'//trim(num2str(size(stored_eigenvalues,2)))//'x'//trim(num2str(size(stored_eigenvalues,3))),'diag')
                call MIO_Print('First 5 TAPW eigenvalues: ['//trim(num2str(eigvals(1),6))//','//trim(num2str(eigvals(2),6))//','//trim(num2str(eigvals(3),6))//','//trim(num2str(eigvals(4),6))//','//trim(num2str(eigvals(5),6))//']','diag')
                call MIO_Print('TAPW eigenvalue range: ['//trim(num2str(minval(eigvals(1:M)),6))//','//trim(num2str(maxval(eigvals(1:M)),6))//']','diag')
             end if

             ! VERIFICATION: Check stored eigenvectors for first k-point
             if (kpoint_index == 1) then
                call MIO_Print('=== EIGENVECTOR STORAGE VERIFICATION (k-point 1) ===','diag')
                call MIO_Print('Stored eigenvectors dimensions: '//trim(num2str(size(stored_eigenvectors,1)))//'x'//trim(num2str(size(stored_eigenvectors,2)))//'x'//trim(num2str(size(stored_eigenvectors,3))),'diag')
                call MIO_Print('Hproj (eigenvectors) norm: '//trim(num2str(sqrt(sum(abs(Hproj)**2)),8)),'diag')
                call MIO_Print('Hproj max element: '//trim(num2str(maxval(abs(Hproj)),8)),'diag')
                call MIO_Print('Hproj min element: '//trim(num2str(minval(abs(Hproj)),8)),'diag')
                ! Check orthogonality of first few eigenvectors
                call MIO_Print('First eigenvector norm: '//trim(num2str(sqrt(sum(abs(Hproj(:,1))**2)),8)),'diag')
                if (M > 1) call MIO_Print('Second eigenvector norm: '//trim(num2str(sqrt(sum(abs(Hproj(:,2))**2)),8)),'diag')
                if (M > 2) call MIO_Print('Third eigenvector norm: '//trim(num2str(sqrt(sum(abs(Hproj(:,3))**2)),8)),'diag')
             end if
          else
             call MIO_Print('ERROR: Cannot store k-point data - array too small','diag')
             call MIO_Print('  M='//trim(num2str(M))//', stored array size='//trim(num2str(size(stored_eigenvectors,1)))//'x'//trim(num2str(size(stored_eigenvectors,2))),'diag')
          end if
       end if

       ! Clean up temporary copy
       if (allocated(Hproj_copy)) deallocate(Hproj_copy)
    end if

    if (tapwDebug) then
       ! Debug: Print TAPW eigenvalues
       if (tapwDebug) print *, "TAPW projected matrix size M =", M
       print *, "TAPW eigenvalues (first 10):", eigvals(1:min(10,M))
       print *, "Max eigenvalue:", maxval(eigvals)
       print *, "Min eigenvalue:", minval(eigvals)
    end if

    if (.not. allocated(eigvals)) stop "eigvals not allocated"

    ! Berry-flux path: hand back the eigenvectors of the requested band window.
    ! ZHEEV('V') leaves them in the columns of Hproj, in ascending-eigenvalue
    ! order, so column b IS sheet b of generate.bands.
#ifdef SEMICL
    if (present(evecOut)) then
       evecOut = cmplx(0.0_dp, 0.0_dp, kind=dp)
       do i = bf_b1, min(bf_b2, M)
          evecOut(1:M, i-bf_b1+1) = Hproj(1:M, i)
       end do
    end if
#endif

    ! Plane-wave (G) composition of the requested sheets (opt-in diagnostic; default .false. changes nothing).
    if (gWeightsTAPW .and. is == 1 .and. present(kpoint_index)) &
         call GWeightsWrite(M, NG, Nlabel, Gx, Gy, Hproj, eigvals, kpoint_index, KLoc)

    ! Layer / sublattice composition of the requested sheets (opt-in; default .false. changes nothing).
    if (layerWeightsTAPW .and. is == 1 .and. present(kpoint_index)) &
         call LayerWeightsWrite(N, M, Nlabel, XArray, Hproj, eigvals, label, original_label_codes, kpoint_index, KLoc)

#ifdef SEMICL
    if (orbMomentTAPW .and. allocated(om_M) .and. is == 1 .and. present(kpoint_index)) then
       if (useDenseMatrixTAPW) then
          if (kpoint_index >= 1 .and. kpoint_index <= size(om_M,2)) then
             call OrbMomentAtK(N, M, XArray, Hproj, eigvals, kpoint_index, KLoc, &
                               cell_real, maxN, hopp, NList, Nneigh, neighCell)
          end if
       else if (.not. om_warned) then
          om_warned = .true.
          call MIO_Print('Diag.OrbMoment needs useDenseMatrixTAPW .true. - '// &
               'the moment was NOT computed','diag')
       end if
    end if
#endif
    if (tapwDebug) print *, "TAPW projected matrix size M = ", M

    ! Handle size mismatch gracefully - TAPW can produce more/fewer states than atoms
    if (M > size(ELoc)) then
        if (tapwDebug) then
           print *, "Warning: More TAPW states (M=", M, ") than output array size (", size(ELoc), ")"
           print *, "Taking first", size(ELoc), "eigenvalues out of", M, "total TAPW states"
        end if
        ELoc = eigvals(1:size(ELoc))  ! Take first size(ELoc) eigenvalues
    else
        ELoc(1:M) = eigvals           ! Assign meaningful TAPW eigenvalues only
        if (M < size(ELoc)) then
            ELoc(M+1:) = 0.0_dp       ! Zero out unused entries
        end if
    end if

    if (tapwDebug) then
       print *, "TAPW: Using", M, "eigenvalues for output"
    end if

    ! Set M_tapw for band output (module variable) and validate consistency
    if (M_tapw > 0 .and. M_tapw /= M) then
       call MIO_Print('WARNING: M_tapw inconsistency detected','diag')
       call MIO_Print('  Previous M_tapw: '//trim(num2str(M_tapw)),'diag')
       call MIO_Print('  Current M: '//trim(num2str(M)),'diag')
       call MIO_Print('  This may cause array bounds issues in subsequent calculations','diag')
    end if
    M_tapw = M

! eigvals now contains eigenvalues of projected H

    if (tapwDebug) print *, "Deallocating arrays..."

    ! PERFORMANCE OPTIMIZATION: Only deallocate if not using pre-allocated arrays
    ! Use a simple approach: always deallocate, pre-allocated arrays will be reused
    if (allocated(eigvals)) then
        deallocate(eigvals)
        if (tapwDebug) print *, "Deallocated eigvals"
    endif

    if (allocated(ZWorkLoc)) then
        deallocate(ZWorkLoc)
        if (tapwDebug) print *, "Deallocated ZWorkLoc"
    endif

    if (allocated(DWorkLoc)) then
        deallocate(DWorkLoc)
        if (tapwDebug) print *, "Deallocated DWorkLoc"
    endif

    if (allocated(XArray)) then
        deallocate(XArray)
        if (tapwDebug) print *, "Deallocated XArray"
    endif

    if (allocated(Hproj)) then
        deallocate(Hproj)
        if (tapwDebug) print *, "Deallocated Hproj"
    endif

    if (allocated(Gx)) then
        deallocate(Gx)
        if (tapwDebug) print *, "Deallocated Gx"
    endif

    if (allocated(Gy)) then
        deallocate(Gy)
        if (tapwDebug) print *, "Deallocated Gy"
    endif

    if (allocated(label)) then
        deallocate(label)
        if (tapwDebug) print *, "Deallocated label"
    endif

    if (allocated(row_ptr)) then
        deallocate(row_ptr)
        if (tapwDebug) print *, "Deallocated row_ptr"
    endif

    if (allocated(col_ind)) then
        deallocate(col_ind)
        if (tapwDebug) print *, "Deallocated col_ind"
    endif

    if (allocated(values)) then
        deallocate(values)
        if (tapwDebug) print *, "Deallocated values"
    endif

    if (allocated(rand_real)) then
        deallocate(rand_real)
        if (tapwDebug) print *, "Deallocated rand_real"
    endif

    if (allocated(rand_imag)) then
        deallocate(rand_imag)
        if (tapwDebug) print *, "Deallocated rand_imag"
    endif

    if (allocated(a)) then
        deallocate(a)
        if (tapwDebug) print *, "Deallocated a"
    endif

    if (allocated(EVectors)) then
        deallocate(EVectors)
        if (tapwDebug) print *, "Deallocated EVectors"
    endif

    if (allocated(ax)) then
        deallocate(ax)
        if (tapwDebug) print *, "Deallocated ax"
    endif

    if (allocated(rd)) then
        deallocate(rd)
        if (tapwDebug) print *, "Deallocated rd"
    endif

    if (allocated(XArray)) then
        deallocate(XArray)
        if (tapwDebug) print *, "Deallocated XArray"
    endif

    if (allocated(Hproj)) then
        deallocate(Hproj)
        if (tapwDebug) print *, "Deallocated Hproj"
    endif

    if (allocated(ZWorkLoc)) then
        deallocate(ZWorkLoc)
        if (tapwDebug) print *, "Deallocated ZWorkLoc"
    endif

    if (allocated(DWorkLoc)) then
        deallocate(DWorkLoc)
        if (tapwDebug) print *, "Deallocated DWorkLoc"
    endif

    if (tapwDebug) print *, "Done deallocating."

end subroutine DiagH0TAPW

!> @brief TAPW diagonalization using pre-built block Hamiltonian with SOC
!! @details Uses TAPW's specialized algorithms on a pre-built 2N×2N block Hamiltonian
!! @param[in]     N         Number of atoms
!! @param[in]     ns        Number of spin channels
!! @param[in]     is        Spin channel index
!! @param[out]    ELoc      Eigenvalues (N)
!! @param[in]     KLoc      k-point vector
!! @param[in]     cell_real Real space cell matrix
!! @param[in]     HBlock    Pre-built block Hamiltonian (2N×2N)
!! @param[in]     maxN      Maximum number of neighbors
!! @param[in]     hopp      Hopping parameters
!! @param[in]     NList     Neighbor list
!! @param[in]     Nneigh    Number of neighbors per atom
!! @param[in]     neighCell Neighbor cell indices
!! @param[in]     neig      Number of eigenvalues
!! @param[in]     kpoint_index k-point index (optional)
subroutine DiagH0TAPW_withBlockH(N, ns, is, ELoc, KLoc, cell_real, HBlock, maxN, hopp, NList, Nneigh, neighCell, neig, kpoint_index)
    use constants, only : cmplx_i
    use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
    use scf, only : charge, Zch
    use atoms, only : Species, RAt, frac, AtomsSetCart, AtomsSetFrac, layerIndex
    use tbpar, only : U
    use neigh, only : neighD
    use cell,                 only : rcell, ucell
    use ham, only : RashbaSOCterm
    use constants,            only : pi
    use omp_lib,              only : omp_in_parallel
    implicit none
    integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is, neig
    integer, intent(in), optional :: kpoint_index
    real(dp), intent(out) :: ELoc(N)
    real(dp), intent(in) :: KLoc(3), cell_real(3,3)
    complex(dp), intent(inout) :: HBlock(2*N, 2*N)  ! Pre-built block Hamiltonian (inout to allow Hermitian completion)
    complex(dp), intent(in) :: hopp(maxN,N)
    real(dp) :: tol
    real(dp) :: aGtemp
    real(dp) :: tapw_aG

    integer :: nuniq, ii, jj, code
    integer, allocatable :: label_raw(:), uniq_codes(:)

    integer :: i, j, in, info
    real(dp) :: R(3), zz, resid_norm

    ! ARPACK parameters and variables
    integer :: nev, ncv, lworkl, ido, ierr
    character(1) :: bmat
    character(2) :: which
    real(dp) :: sigma
    complex(dp), allocatable :: workd(:), workl(:), resid(:), v(:,:)
    logical :: select(neig)
    real(dp), allocatable :: d(:), rwork(:)
    logical :: rvec
    integer :: ldz, lwork, lrwork
    complex(dp), allocatable :: ZWorkLoc(:)
    real(dp), allocatable :: DWorkLoc(:)

    ! TAPW-specific variables (same as original DiagH0TAPW)
    ! Note: M_tapw is module-level variable, not local
    integer :: M, NG, Nlabel, NGrange
    complex(dp), allocatable :: XArray(:,:), Hproj(:,:), Hproj_copy(:,:)
    complex(dp), allocatable :: tmpX(:,:)  ! Temporary matrix for building block-diagonal XArray
    real(dp), allocatable :: eigvals(:)
    complex(dp), allocatable :: eigvec(:,:)
    real(dp), allocatable :: Gx(:), Gy(:)
    integer, allocatable :: label(:)
    ! Note: tapwDebug and calculateChern are module-level variables, not local
    logical :: useDenseMatrixTAPW = .true.

    ! Debug variables for HBlock and Hproj checking
    ! Note: ii and jj are already declared above (line 7370)
    real(dp) :: diag_min, diag_max, proj_diag_min, proj_diag_max
    real(dp) :: eig_full_min, eig_full_max, eig_full_range
    integer :: n_positive, n_negative, n_zero
    complex(dp), allocatable :: diag_temp(:), proj_diag_temp(:)
    ! Debug variables for block projection structure check
    real(dp) :: max_offblock, offblock_val, max_block_diff, block_diff
    real(dp) :: max_hblock_diff, hblock_diff, max_hblock_offdiag
    real(dp) :: eig_pair_diff, max_degen_diff
    integer :: eig_pair_count, max_degen_diff_idx
    integer :: check_i, check_j
    ! Debug variables for unique eigenvalue counting
    real(dp) :: unique_tol
    integer :: n_unique, unique_check_i
    logical :: is_unique
    ! Debug variables for Rashba bond checking
    integer :: test_i, test_j, test_in, test_count
    real(dp) :: max_rashba, test_dist, test_acc
    complex(dp) :: test_rashba
    ! Debug variables for eigenvalue checking
    integer :: n_eig_check
    real(dp) :: eig_min, eig_max, eig_range
    ! Variables for position handling (matching DiagH0TAPW)
    real(dp) :: temp_positions(3, N)  ! Temporary storage for position swapping
    logical :: original_frac_state     ! Original coordinate system state
    real(dp), allocatable :: rigid_positions(:,:)  ! Rigid reference positions for X matrix

    ! TAPW k_ref calculation variables (all declarations must come before executable statements)
    ! Note: tapw_aG is already declared above (line 7369)
    real(dp) :: k_ref(2), refF(2), aG, rotation_angle, cos_rot, sin_rot
    real(dp) :: sGlattice(2,2), rG(2,2), det, sGlattice_inv(2,2)
    logical :: useTriangularDistanceOrder
    ! Variables for Hermitian completion check
    integer :: herm_i, herm_j
    real(dp) :: max_herm_error
    ! Variables for Hproj Hermiticity check
    real(dp) :: herm_error, max_herm_error_proj
    integer :: proj_i, proj_j
    ! Variables for Hproj structure check
    ! Note: test_i, test_j already declared above (line 7408), ii already declared earlier
    real(dp) :: max_offdiag, test_offdiag, trace_Hproj, norm_Hproj

    if (tapwBothValleys) call MIO_Kill('Diag.TAPWBothValleys is not implemented on the block (SOC) TAPW path', &
         'diag','DiagH0TAPW_withBlockH')

    ! TAPW parameters: all read from module-level variables (set via input parameters)
    ! Calculate NGrange (same logic as original DiagH0TAPW)
    if (checkTAPWUnitary) then
        call calculate_optimal_NGrange_for_complete_BZ(abs(physicalTwistAngle), tapwNG, NGrange)
    else
        NGrange = tapwNG
    end if

    ! Read TAPW-specific input parameters outside OpenMP; inside OpenMP use cached config
    ! (same as original DiagH0TAPW)
    if (.not. omp_in_parallel()) then
       call MIO_InputParameter('TAPW.aG', tapw_aG, 2.46019_dp)
    else
       tapw_aG = tapw_cfg_aG
    end if
    aG = tapw_aG

    ! Reference fractional coordinates in BZ (valley-dependent)
    if (useKprimeValley) then
       refF = (/ 1.0_dp/3.0_dp, 2.0_dp/3.0_dp /)  ! K' valley: [1/3, 2/3]
    else
       refF = (/ 2.0_dp/3.0_dp, 1.0_dp/3.0_dp /)  ! K valley: [2/3, 1/3] (default)
    end if

    ! Build graphene lattice vectors (60-degree structure) - COLUMN storage like rcell
    sGlattice(1,1) = 1.0_dp * aG
    sGlattice(2,1) = 0.0_dp
    sGlattice(1,2) = cos(60.0_dp*pi/180.0_dp) * aG
    sGlattice(2,2) = sin(60.0_dp*pi/180.0_dp) * aG

    ! Apply rotation if needed (using module-level moireAngle)
    rotation_angle = moireAngle * pi / 180.0_dp
    cos_rot = cos(rotation_angle)
    sin_rot = sin(rotation_angle)
    sGlattice = matmul(reshape((/ cos_rot, sin_rot, -sin_rot, cos_rot /), (/2,2/)), sGlattice)

    ! Convert to reciprocal space: rG = 2π * inverse(lattice.T)
    det = sGlattice(1,1)*sGlattice(2,2) - sGlattice(1,2)*sGlattice(2,1)
    if (abs(det) < 1e-10_dp) then
        call MIO_Kill('Lattice matrix is singular in DiagH0TAPW_withBlockH', 'diag')
    end if
    sGlattice_inv(1,1) =  sGlattice(2,2) / det
    sGlattice_inv(1,2) = -sGlattice(1,2) / det
    sGlattice_inv(2,1) = -sGlattice(2,1) / det
    sGlattice_inv(2,2) =  sGlattice(1,1) / det
    rG = 2.0_dp * pi * transpose(sGlattice_inv)

    ! Calculate reference K-point in reciprocal space
    k_ref = refF(1)*rG(:,1) + refF(2)*rG(:,2)

    if (tapwDebug) then
       call MIO_Print('DiagH0TAPW_withBlockH: Starting TAPW with block Hamiltonian', 'diag')
       call MIO_Print('  Using N_G = '//trim(num2str(tapwNG))//', NGrange = '//trim(num2str(NGrange)), 'diag')
       call MIO_Print('  Moire angle: '//trim(num2str(moireAngle,6))//' degrees','diag')
       if (useKprimeValley) then
          call MIO_Print('  Valley selection: K'' valley [1/3, 2/3]','diag')
       else
          call MIO_Print('  Valley selection: K valley [2/3, 1/3]','diag')
       end if
    end if

    ! Generate G-vectors (same as original DiagH0TAPW)
    allocate(Gx(1000), Gy(1000))  ! Temporary allocation
    if (useTriangularTruncation) then
       call MIO_Print('Using triangular G-vector truncation','diag')
       call MIO_InputParameter('TAPW.UseTriangularDistanceOrder', useTriangularDistanceOrder, .true.)
       call generate_triangular_G_list(rcell, k_ref, NGrange, Gx, Gy, NG, rG, useTriangularDistanceOrder)
    else if (checkTAPWUnitary) then
        call generate_shifted_G_list_with_graphene_BZ(rcell, k_ref, NGrange, Gx, Gy, NG, [0.0_dp, 0.0_dp], rG)
    else
        call generate_shifted_G_list_reduced(rcell, k_ref, NGrange, Gx, Gy, NG)
    end if

    ! Resize G-vectors to actual size
    deallocate(Gx, Gy)
    allocate(Gx(NG), Gy(NG))
    if (useTriangularTruncation) then
       call generate_triangular_G_list(rcell, k_ref, NGrange, Gx, Gy, NG, rG, useTriangularDistanceOrder)
    else if (checkTAPWUnitary) then
        call generate_shifted_G_list_with_graphene_BZ(rcell, k_ref, NGrange, Gx, Gy, NG, [0.0_dp, 0.0_dp], rG)
    else
        call generate_shifted_G_list_reduced(rcell, k_ref, NGrange, Gx, Gy, NG)
    end if

    ! Generate labels (same as original DiagH0TAPW)
    allocate(label(N))
    do i = 1, N
       label(i) = 10 * layerIndex(i) + Species(i)  ! e.g., layer 1 sublattice A -> 11, layer 2 sublattice B -> 22
    end do

    ! Compute number of unique labels dynamically
    call compute_unique_labels(label, N, Nlabel)

    ! Remap labels to contiguous indices 1..Nlabel
    call remap_labels_to_contiguous(label, N, Nlabel)

    ! Calculate M (spinless TAPW basis size)
    M = NG * Nlabel
    ! For Rashba SOC with block-diagonal projection, we have 2M eigenvalues (spin-split)
    ! For forceBlockTAPW without SOC, we have 2M eigenvalues but they're degenerate pairs
    ! Set M_tapw based on whether we need to collapse duplicates
    if (forceBlockTAPW .and. .not. RashbaSOCterm) then
       ! No SOC: eigenvalues are degenerate pairs, will collapse to M unique ones
       M_tapw = M
    else
       ! Rashba SOC: real spin-split bands, need all 2M eigenvalues
       M_tapw = 2*M
    end if

    if (tapwDebug) then
       call MIO_Print('  Generated NG = '//trim(num2str(NG))//' G-vectors, Nlabel = '//trim(num2str(Nlabel))//', M = '//trim(num2str(M)), 'diag')
    end if

    ! Allocate TAPW arrays
    ! CRITICAL: Use block-diagonal projector (2N × 2M) to preserve spin channels
    ! This prevents collapsing both spin channels into the same M-dimensional space
    ! Structure: XArray = [X_up  0  ]  (rows 1..N: spin-up, rows N+1..2N: spin-down)
    !                    [0    X_up]  (cols 1..M: spin-up subspace, cols M+1..2M: spin-down subspace)
    allocate(XArray(2*N, 2*M))  ! Block-diagonal: 2N×2M for spinful TAPW
    allocate(Hproj(2*M, 2*M))   ! Projected Hamiltonian: 2M×2M (preserves spin structure)
    allocate(eigvals(2*M))
    allocate(eigvec(2*M, 2*M))

    ! Handle atomic positions for X matrix construction (same as DiagH0TAPW)
    if (useRigidPositions) then
       ! Use rigid reference positions for X matrix construction
       ! This ensures perfect TAPW unitarity even with lattice reconstruction
       call MIO_Print('Using rigid reference positions for TAPW X matrix construction','diag')

       ! Allocate and read rigid positions
       allocate(rigid_positions(3, N))
       call read_rigid_positions_for_tapw('generateInit.xyz', N, rigid_positions)

       ! Note: We need to temporarily modify the coordinate system for wrapping
       ! Store current positions and coordinate system state
       temp_positions = Rat  ! Save current (relaxed) positions
       original_frac_state = frac

       ! Replace current positions with rigid positions for wrapping
       Rat = rigid_positions
       deallocate(rigid_positions)

       call MIO_Print('Rigid positions loaded, wrapping to single unit cell','diag')
    else
       ! Use current positions in Rat directly (default behavior)
       call MIO_Print('Using current atomic positions for TAPW X matrix construction','diag')

       ! Store current state for consistency with rigid position path
       temp_positions = Rat  ! Save current positions (no change needed)
       original_frac_state = frac
    end if

    if (tapwDebug) then
       if (useRigidPositions) then
          call MIO_Print('DEBUG: First 3 atomic positions BEFORE build_X (rigid):','diag')
       else
          call MIO_Print('DEBUG: First 3 atomic positions BEFORE build_X (current):','diag')
       end if
       do i = 1, min(3, N)
           call MIO_Print('  Atom '//trim(num2str(i))//': ['//trim(num2str(Rat(1,i),6))//','//&
                         trim(num2str(Rat(2,i),6))//','//trim(num2str(Rat(3,i),6))//']','diag')
       end do
    end if

    if (useRigidPositions) then
       ! Convert rigid positions to fractional coordinates relative to moiré cell
       if (.not. frac) call AtomsSetFrac()

       ! Wrap fractional coordinates to [0,1) to ensure single unit cell
       do i = 1, N
          Rat(1,i) = Rat(1,i) - floor(Rat(1,i))
          Rat(2,i) = Rat(2,i) - floor(Rat(2,i))
          ! Don't wrap z-coordinate (only x,y are periodic in 2D moiré)
       end do

       ! Convert wrapped rigid positions back to Cartesian for build_X
       call AtomsSetCart()
       call MIO_Print('Rigid positions wrapped and converted to Cartesian for X matrix','diag')
    else
       ! For current positions, ensure we're in the right coordinate system for build_X
       if (frac) call AtomsSetCart()
       call MIO_Print('Using current positions directly for X matrix construction','diag')
    end if

    ! Build spinless TAPW basis once
    allocate(tmpX(N, M))
    call build_X(tmpX, Rat(1,:), Rat(2,:), label, Gx, Gy, N, NG, Nlabel)
    if (tapwLowdin) call lowdin_orthonormalize_X(tmpX, N, M)

    ! RESTORE original atomic positions for TB calculations (if they were modified)
    if (useRigidPositions) then
       Rat = temp_positions
       if (original_frac_state .and. .not. frac) call AtomsSetFrac()
       if (.not. original_frac_state .and. frac) call AtomsSetCart()
       call MIO_Print('Restored original relaxed positions for TB Hamiltonian construction','diag')
    end if

    ! Place into block-diagonal structure:
    ! Columns 1..M: spin-up sector (rows 1..N)
    ! Columns M+1..2*M: spin-down sector (rows N+1..2*N)
    XArray = (0.0_dp, 0.0_dp)
    XArray(1:N, 1:M) = tmpX          ! Spin-up block
    XArray(N+1:2*N, M+1:2*M) = tmpX  ! Spin-down block (same basis)
    deallocate(tmpX)

    ! Debug: Show M and NG before projection
    if (tapwDebug .and. present(kpoint_index)) then
       call MIO_Print('DiagH0TAPW_withBlockH: About to project HBlock (2N='//trim(num2str(2*N))//'×'//trim(num2str(2*N))//') to Hproj (2M='//trim(num2str(2*M))//'×'//trim(num2str(2*M))//')','diag')
       call MIO_Print('  Using block-diagonal projector XArray (2N='//trim(num2str(2*N))//'×2M='//trim(num2str(2*M))//') to preserve spin structure','diag')
       call MIO_Print('  This will perform: Hproj = XArray† * HBlock * XArray','diag')
       call MIO_Print('  Matrix multiplication complexity: O(2M × 2N × 2N) ≈ '//trim(num2str(2*M*2*N*2*N))//' operations','diag')
    end if

    ! Debug: Check HBlock structure before projection (for forceBlockTAPW testing)
    if (forceBlockTAPW .and. present(kpoint_index) .and. kpoint_index <= 2) then
       call MIO_Print('=== HBlock STRUCTURE CHECK (before projection) ===','diag')
       call MIO_Print('  HBlock dimensions: '//trim(num2str(2*N))//'×'//trim(num2str(2*N)),'diag')

       ! Check diagonal blocks (should be identical without SOC)
       max_hblock_diff = 0.0_dp
       do check_j = 1, min(N, 10)
          do check_i = 1, min(N, 10)
             hblock_diff = abs(HBlock(check_i, check_j) - HBlock(check_i+N, check_j+N))
             if (hblock_diff > max_hblock_diff) then
                max_hblock_diff = hblock_diff
             end if
          end do
       end do
       call MIO_Print('  Max difference between HBlock diagonal blocks (first 10×10): '//trim(num2str(max_hblock_diff,8)),'diag')

       ! Check off-diagonal blocks (should be zero without SOC/Rashba)
       max_hblock_offdiag = 0.0_dp
       do check_j = 1, min(N, 10)
          do check_i = 1, min(N, 10)
             offblock_val = abs(HBlock(check_i, check_j+N))
             if (offblock_val > max_hblock_offdiag) then
                max_hblock_offdiag = offblock_val
             end if
             offblock_val = abs(HBlock(check_i+N, check_j))
             if (offblock_val > max_hblock_offdiag) then
                max_hblock_offdiag = offblock_val
             end if
          end do
       end do
       call MIO_Print('  Max off-diagonal block element in HBlock (first 10×10): '//trim(num2str(max_hblock_offdiag,8)),'diag')

       ! Show first few diagonal elements from each block
       if (N >= 5) then
          call MIO_Print('  First 5 HBlock diagonal elements (upper block): ['//&
                        trim(num2str(real(HBlock(1,1)),6))//', '//trim(num2str(real(HBlock(2,2)),6))//', '//&
                        trim(num2str(real(HBlock(3,3)),6))//', '//trim(num2str(real(HBlock(4,4)),6))//', '//&
                        trim(num2str(real(HBlock(5,5)),6))//']','diag')
          call MIO_Print('  First 5 HBlock diagonal elements (lower block): ['//&
                        trim(num2str(real(HBlock(N+1,N+1)),6))//', '//trim(num2str(real(HBlock(N+2,N+2)),6))//', '//&
                        trim(num2str(real(HBlock(N+3,N+3)),6))//', '//trim(num2str(real(HBlock(N+4,N+4)),6))//', '//&
                        trim(num2str(real(HBlock(N+5,N+5)),6))//']','diag')
       end if

       if (max_hblock_diff > 1e-6_dp) then
          call MIO_Print('  WARNING: HBlock diagonal blocks differ! This may cause Hproj blocks to differ.','diag')
          call MIO_Print('  Check if SOC terms are enabled (Zterm, IntrinsicSOCterm, PIASOCterm)','diag')
       end if
    end if

    ! Debug: Check HBlock before projection
    if (tapwDebug .and. present(kpoint_index) .and. kpoint_index <= 2) then
       allocate(diag_temp(2*N))

       call MIO_Print('DiagH0TAPW_withBlockH: Checking HBlock before projection (k-point '//trim(num2str(kpoint_index))//')','diag')
       call MIO_Print('  HBlock dimensions: '//trim(num2str(2*N))//'×'//trim(num2str(2*N)),'diag')

       ! Extract diagonal elements
       do ii = 1, 2*N
          diag_temp(ii) = HBlock(ii, ii)
       end do
       diag_min = minval(real(diag_temp))
       diag_max = maxval(real(diag_temp))
       call MIO_Print('  HBlock diagonal range: ['//trim(num2str(diag_min,6))//', '//trim(num2str(diag_max,6))//']','diag')

       ! Check Rashba spin-flip terms by scanning neighbors
       ! Instead of checking arbitrary pairs, check actual neighbor bonds
       test_acc = 2.46_dp / sqrt(3.0_dp)  ! Graphene nearest neighbor cutoff
       max_rashba = 0.0_dp
       test_count = 0

       ! Scan first few atoms and their neighbors to find Rashba terms
       do test_i = 1, min(10, N)
          do test_j = 1, Nneigh(test_i)
             test_in = NList(test_j, test_i)
             ! Check if this should be a Rashba bond (same layer, different species, nearest neighbor)
             test_dist = sqrt(neighD(1, test_j, test_i)**2 + neighD(2, test_j, test_i)**2)
             if (layerIndex(test_i) == layerIndex(test_in) .and. &
                 Species(test_i) /= Species(test_in) .and. &
                 test_dist < test_acc * 1.1_dp .and. &
                 test_i /= test_in) then
                test_rashba = HBlock(test_i, test_in + N)
                max_rashba = max(max_rashba, abs(test_rashba))
                test_count = test_count + 1
                if (test_count <= 5) then
                   call MIO_Print('  Rashba bond: HBlock('//trim(num2str(test_i))//','//trim(num2str(test_in+N))//') = '//trim(num2str(real(test_rashba),6))//' + i*'//trim(num2str(aimag(test_rashba),6)),'diag')
                end if
             end if
          end do
       end do

       call MIO_Print('  Total Rashba bonds checked: '//trim(num2str(test_count))//', max |Rashba| = '//trim(num2str(max_rashba,6)),'diag')

       deallocate(diag_temp)
    end if

    ! Project the block Hamiltonian using TAPW transformation
    ! CRITICAL: HBlock from BuildBlockHamiltonianOnly may only have lower triangle set
    ! ZGEMM requires full matrix, so we need to complete Hermitian structure if needed
    ! However, BuildBlockHamiltonianOnly actually sets HBlock(in, i) for all neighbors,
    ! which may include upper triangle elements. For safety, we check/complete Hermiticity.
    if (tapwDebug .and. present(kpoint_index)) then
       call MIO_Print('DiagH0TAPW_withBlockH: Starting matrix multiplication (this may take a while for large M)...','diag')
       ! Check: verify HBlock is Hermitian (should be naturally Hermitian from BuildBlockHamiltonianOnly)
       ! If not, this indicates a bug that should be fixed at the source
       max_herm_error = 0.0_dp
       do check_i = 1, min(10, 2*N)
          do check_j = check_i+1, min(10, 2*N)
             if (abs(HBlock(check_i, check_j) - conjg(HBlock(check_j, check_i))) > max_herm_error) then
                max_herm_error = abs(HBlock(check_i, check_j) - conjg(HBlock(check_j, check_i)))
             end if
          end do
       end do
       call MIO_Print('  HBlock Hermiticity check (first 10x10): max error = '//trim(num2str(max_herm_error,8)),'diag')
       if (max_herm_error > 1e-10_dp) then
          call MIO_Print('  WARNING: HBlock is not Hermitian! This indicates a bug in BuildBlockHamiltonianOnly.','diag')
       end if
    end if

    ! Note: We do NOT force Hermiticity here. If HBlock is not Hermitian, that's a bug
    ! that should be fixed in BuildBlockHamiltonianOnly, not hidden by patching here.

    ! Project block Hamiltonian using TAPW transformation
    ! where XArray is (2N×M) and HBlock is (2N×2N)
    !
    ! HBlock structure (from BuildBlockHamiltonianOnly):
    !   HBlock(1:N, 1:N) = H_up (spin-up block)
    !   HBlock(N+1:2*N, N+1:2*N) = H_dn (spin-down block)
    !   HBlock(1:N, N+1:2*N) = H_Rashba (spin-flip up→down)
    !   HBlock(N+1:2*N, 1:N) = H_Rashba* (spin-flip down→up, Hermitian conjugate)
    !
    ! XArray structure (from build_X):
    !   XArray(1:N, 1:M) = X_up (TAPW basis for spin-up, exp(i*G·r))
    !   XArray(N+1:2*N, 1:M) = X_up (same TAPW basis for spin-down)
    !
    ! Mathematical expansion:
    !        = X_up†*H_up*X_up + X_up†*H_Rashba*X_up + X_up†*H_Rashba**X_up + X_up†*H_dn*X_up
    !        = X_up†*(H_up + H_dn)*X_up + X_up†*(H_Rashba + H_Rashba*)*X_up
    !
    ! This matches the non-Rashba case when H_Rashba=0, but adds Rashba spin-flip terms.
    ! Project to (2M × 2M) space to preserve spin structure
    call transform_block_hamiltonian_tapw(2*N, 2*M, XArray, HBlock, Hproj)
    if (tapwDebug .and. present(kpoint_index)) then
       call MIO_Print('DiagH0TAPW_withBlockH: Matrix multiplication completed','diag')

       ! Debug: Verify Hproj is Hermitian (should be if HBlock is Hermitian)
       if (kpoint_index <= 2) then
          max_herm_error_proj = 0.0_dp
          do proj_j = 1, min(2*M, 20)
             do proj_i = proj_j+1, min(2*M, 20)
                herm_error = abs(Hproj(proj_i, proj_j) - conjg(Hproj(proj_j, proj_i)))
                if (herm_error > max_herm_error_proj) then
                   max_herm_error_proj = herm_error
                end if
             end do
          end do
          call MIO_Print('  Hproj Hermiticity check (first 20×20): max error = '//trim(num2str(max_herm_error_proj,8)),'diag')

          ! Debug: Compare a few elements of Hproj with expected values
          if (2*M >= 5) then
             call MIO_Print('  Hproj(1:5,1:5) diagonal: ['//trim(num2str(real(Hproj(1,1)),6))//', '//&
                           trim(num2str(real(Hproj(2,2)),6))//', '//trim(num2str(real(Hproj(3,3)),6))//', '//&
                           trim(num2str(real(Hproj(4,4)),6))//', '//trim(num2str(real(Hproj(5,5)),6))//']','diag')
          end if
       end if
    end if

    ! Debug: Check Hproj after projection and compare with expected structure
    ! For no-SOC case (forceBlockTAPW), Hproj should be block-diagonal with identical blocks
    if ((tapwDebug .or. forceBlockTAPW) .and. present(kpoint_index) .and. kpoint_index <= 2) then
       if (forceBlockTAPW) then
          call MIO_Print('=== BLOCK TAPW PROJECTION DEBUG (no SOC) ===','diag')
          call MIO_Print('  Expected: Hproj should be block-diagonal [H_TAPW 0; 0 H_TAPW]','diag')
          call MIO_Print('  Hproj dimensions: '//trim(num2str(2*M))//'×'//trim(num2str(2*M)),'diag')

          ! Check if Hproj is block-diagonal (off-diagonal blocks should be ~0 without SOC)
          max_offblock = 0.0_dp
          ! Check upper-right block (rows 1..M, cols M+1..2*M)
          do check_j = M+1, min(2*M, M+10)
             do check_i = 1, min(M, 10)
                offblock_val = abs(Hproj(check_i, check_j))
                if (offblock_val > max_offblock) then
                   max_offblock = offblock_val
                end if
             end do
          end do
          ! Check lower-left block (rows M+1..2*M, cols 1..M)
          do check_j = 1, min(M, 10)
             do check_i = M+1, min(2*M, M+10)
                offblock_val = abs(Hproj(check_i, check_j))
                if (offblock_val > max_offblock) then
                   max_offblock = offblock_val
                end if
             end do
          end do
          call MIO_Print('  Max off-diagonal block element (first 10×10): '//trim(num2str(max_offblock,8)),'diag')

          ! Compare diagonal blocks (should be identical without SOC)
          max_block_diff = 0.0_dp
          do check_j = 1, min(M, 10)
             do check_i = 1, min(M, 10)
                block_diff = abs(Hproj(check_i, check_j) - Hproj(check_i+M, check_j+M))
                if (block_diff > max_block_diff) then
                   max_block_diff = block_diff
                end if
             end do
          end do
          call MIO_Print('  Max difference between diagonal blocks (first 10×10): '//trim(num2str(max_block_diff,8)),'diag')

          if (max_block_diff > 1e-6_dp) then
             call MIO_Print('  WARNING: Diagonal blocks differ! They should be identical without SOC','diag')
          else
             call MIO_Print('  SUCCESS: Diagonal blocks are identical (projection structure is correct)','diag')
          end if

          ! Show first few eigenvalues from each block
          if (2*M >= 5) then
             call MIO_Print('  First 5 Hproj diagonal elements (upper block): ['//&
                           trim(num2str(real(Hproj(1,1)),6))//', '//trim(num2str(real(Hproj(2,2)),6))//', '//&
                           trim(num2str(real(Hproj(3,3)),6))//', '//trim(num2str(real(Hproj(4,4)),6))//', '//&
                           trim(num2str(real(Hproj(5,5)),6))//']','diag')
             call MIO_Print('  First 5 Hproj diagonal elements (lower block): ['//&
                           trim(num2str(real(Hproj(M+1,M+1)),6))//', '//trim(num2str(real(Hproj(M+2,M+2)),6))//', '//&
                           trim(num2str(real(Hproj(M+3,M+3)),6))//', '//trim(num2str(real(Hproj(M+4,M+4)),6))//', '//&
                           trim(num2str(real(Hproj(M+5,M+5)),6))//']','diag')
          end if

          ! Check if eigenvalues are doubly degenerate (should be for identical blocks)
          ! Note: This check happens AFTER diagonalization, so we need to check eigvals array
          ! For now, we'll check this after ZHEEV is called
       end if
    end if

    if (tapwDebug .and. present(kpoint_index) .and. kpoint_index <= 2) then
       allocate(proj_diag_temp(2*M))

       call MIO_Print('DiagH0TAPW_withBlockH: Checking Hproj after projection','diag')
       call MIO_Print('  Hproj dimensions: '//trim(num2str(2*M))//'×'//trim(num2str(2*M)),'diag')

       ! Extract diagonal elements (now 2M×2M)
       do jj = 1, 2*M
          proj_diag_temp(jj) = Hproj(jj, jj)
       end do
       proj_diag_min = minval(real(proj_diag_temp))
       proj_diag_max = maxval(real(proj_diag_temp))
       call MIO_Print('  Hproj diagonal range: ['//trim(num2str(proj_diag_min,6))//', '//trim(num2str(proj_diag_max,6))//']','diag')

       ! Debug: Check off-diagonal elements to verify Rashba terms are present in Hproj
       if (2*M >= 10) then
          max_offdiag = 0.0_dp
          do test_j = 1, min(10, 2*M)
             do test_i = test_j+1, min(10, 2*M)
                test_offdiag = abs(Hproj(test_i, test_j))
                if (test_offdiag > max_offdiag) then
                   max_offdiag = test_offdiag
                end if
             end do
          end do
          call MIO_Print('  Hproj max off-diagonal (first 10×10): '//trim(num2str(max_offdiag,6)),'diag')

          ! Debug: Show specific off-diagonal elements to check k-dependence
          ! Also check cross-spin terms (should be non-zero with Rashba)
          if (kpoint_index <= 2 .and. 2*M >= 5) then
             call MIO_Print('  Hproj(1,2) = '//trim(num2str(real(Hproj(1,2)),6))//' + i*'//trim(num2str(aimag(Hproj(1,2)),6)),'diag')
             call MIO_Print('  Hproj(2,1) = '//trim(num2str(real(Hproj(2,1)),6))//' + i*'//trim(num2str(aimag(Hproj(2,1)),6)),'diag')
             if (2*M >= M+1) then
                call MIO_Print('  Hproj(1,M+1) [cross-spin] = '//trim(num2str(real(Hproj(1,M+1)),6))//' + i*'//trim(num2str(aimag(Hproj(1,M+1)),6)),'diag')
             end if
          end if
       end if

       ! Debug: Check trace and norm of Hproj to verify it's not all zeros
       trace_Hproj = 0.0_dp
       norm_Hproj = 0.0_dp
       do jj = 1, min(2*M, 100)  ! Check first 100 diagonal elements
          trace_Hproj = trace_Hproj + real(Hproj(jj, jj))
          norm_Hproj = norm_Hproj + abs(Hproj(jj, jj))**2
          do ii = 1, min(2*M, 100)
             if (ii /= jj) then
                norm_Hproj = norm_Hproj + abs(Hproj(ii, jj))**2
             end if
          end do
       end do
       call MIO_Print('  Hproj trace (first 100×100): '//trim(num2str(trace_Hproj,6))//', norm²: '//trim(num2str(norm_Hproj,6)),'diag')

       deallocate(proj_diag_temp)
    end if

    ! Make a copy of Hproj before diagonalization (it gets overwritten by ZHEEV)
    ! This copy is needed for Berry curvature calculation
    if (calculateChern .and. present(kpoint_index)) then
       allocate(Hproj_copy(2*M, 2*M))
       Hproj_copy = Hproj
    end if

    ! Diagonalize projected Hamiltonian (now 2M×2M to preserve spin structure)
    lwork = max(1, 2*(2*M)-1)
    lrwork = max(1, 3*(2*M)-2)
    allocate(ZWorkLoc(lwork))
    allocate(DWorkLoc(lrwork))

    call ZHEEV('V', 'U', 2*M, Hproj, 2*M, eigvals, ZWorkLoc, lwork, DWorkLoc, info)

    if (info /= 0) then
       call MIO_Print('ZHEEV failed in DiagH0TAPW_withBlockH: info = '//trim(num2str(info)), 'diag')
       call MIO_Kill('Error in TAPW block Hamiltonian diagonalization', 'diag', 'DiagH0TAPW_withBlockH')
    end if

    ! Debug: Check eigenvalue degeneracy after diagonalization (for forceBlockTAPW testing)
    if (forceBlockTAPW .and. present(kpoint_index) .and. kpoint_index <= 2) then
       call MIO_Print('=== EIGENVALUE DEGENERACY CHECK (after diagonalization) ===','diag')
       call MIO_Print('  Total eigenvalues: '//trim(num2str(2*M)),'diag')
       call MIO_Print('  Expected: Each eigenvalue should appear twice (degenerate pairs)','diag')
       call MIO_Print('  Note: For identical diagonal blocks, eigvals should be sorted pairs','diag')

       max_degen_diff = 0.0_dp
       max_degen_diff_idx = 0
       eig_pair_count = 0

       ! Check consecutive pairs: eigvals(2i-1) ≈ eigvals(2i) for i=1..M
       ! For block-diagonal matrix with identical blocks, ZHEEV sorts eigenvalues as:
       ! λ₁, λ₁, λ₂, λ₂, ..., λₘ, λₘ (each eigenvalue appears twice consecutively)
       ! NOT as: λ₁...λₘ, λ₁...λₘ (all first block, then all second block)
       call MIO_Print('  Note: For block-diagonal [A 0; 0 A], ZHEEV sorts eigenvalues','diag')
       call MIO_Print('        as consecutive pairs: (λ₁,λ₁), (λ₂,λ₂), ..., (λₘ,λₘ)','diag')

       do check_i = 1, min(M, 10)
          ! Check if eigvals(2i-1) ≈ eigvals(2i) (consecutive degenerate pairs)
          eig_pair_diff = abs(eigvals(2*check_i-1) - eigvals(2*check_i))
          if (eig_pair_diff < 1e-6_dp) then
             eig_pair_count = eig_pair_count + 1
          end if
          if (eig_pair_diff > max_degen_diff) then
             max_degen_diff = eig_pair_diff
             max_degen_diff_idx = check_i
          end if
       end do

       call MIO_Print('  Consecutive pair check (first 10 pairs: 2i-1 and 2i): '//trim(num2str(eig_pair_count))//'/'//trim(num2str(min(M,10)))//' pairs match','diag')
       if (max_degen_diff > 1e-6_dp) then
          call MIO_Print('  WARNING: Max eigenvalue pair difference: '//trim(num2str(max_degen_diff,8))//' at pair '//trim(num2str(max_degen_diff_idx)),'diag')
          call MIO_Print('  First 5 consecutive eigenvalue pairs:','diag')
          do check_i = 1, min(5, M)
             call MIO_Print('    Pair '//trim(num2str(check_i))//': eigvals('//trim(num2str(2*check_i-1))//') = '//trim(num2str(eigvals(2*check_i-1),6))//', eigvals('//trim(num2str(2*check_i))//') = '//trim(num2str(eigvals(2*check_i),6))//', diff = '//trim(num2str(abs(eigvals(2*check_i-1)-eigvals(2*check_i)),8)),'diag')
          end do
       else
          call MIO_Print('  SUCCESS: Eigenvalues form proper consecutive pairs (eigvals(2i-1) ≈ eigvals(2i))','diag')
          call MIO_Print('  First 5 consecutive eigenvalue pairs (should be identical):','diag')
          do check_i = 1, min(5, M)
             call MIO_Print('    Pair '//trim(num2str(check_i))//': '//trim(num2str(eigvals(2*check_i-1),6))//' (appears twice at indices '//trim(num2str(2*check_i-1))//','//trim(num2str(2*check_i))//')','diag')
          end do
    end if

       ! Also check if we have exactly M unique eigenvalues (each appearing twice)
       ! This is a stronger check: count unique eigenvalues (within tolerance)
       unique_tol = 1e-6_dp
       n_unique = 0
       do unique_check_i = 1, 2*M
          is_unique = .true.
          ! Check if this eigenvalue is different from all previous ones
          do check_j = 1, unique_check_i-1
             if (abs(eigvals(unique_check_i) - eigvals(check_j)) < unique_tol) then
                is_unique = .false.
                exit
             end if
          end do
          if (is_unique) then
             n_unique = n_unique + 1
          end if
       end do
       call MIO_Print('  Unique eigenvalue count: '//trim(num2str(n_unique))//' (expected: '//trim(num2str(M))//' unique eigenvalues, each appearing twice)','diag')
       if (n_unique == M) then
          call MIO_Print('  SUCCESS: Exactly M unique eigenvalues (perfect degeneracy)','diag')
       else
          call MIO_Print('  WARNING: Expected M='//trim(num2str(M))//' unique eigenvalues, found '//trim(num2str(n_unique)),'diag')
       end if
    end if

    ! Debug: Check FULL eigenvalue range BEFORE extraction
    ! This is critical to understand if positive eigenvalues exist
    if (tapwDebug .and. present(kpoint_index) .and. kpoint_index <= 2) then
       eig_full_min = minval(eigvals(1:2*M))
       eig_full_max = maxval(eigvals(1:2*M))
       eig_full_range = eig_full_max - eig_full_min
       n_positive = count(eigvals(1:2*M) > 0.0_dp)
       n_negative = count(eigvals(1:2*M) < 0.0_dp)
       n_zero = count(abs(eigvals(1:2*M)) < 1e-10_dp)

       call MIO_Print('=== FULL EIGENVALUE ANALYSIS (2M='//trim(num2str(2*M))//' total eigenvalues) ===','diag')
       call MIO_Print('  Full eigenvalue range: ['//trim(num2str(eig_full_min,6))//', '//trim(num2str(eig_full_max,6))//']','diag')
       call MIO_Print('  Eigenvalue span: '//trim(num2str(eig_full_range,6)),'diag')
       call MIO_Print('  Positive eigenvalues: '//trim(num2str(n_positive))//' / '//trim(num2str(2*M)),'diag')
       call MIO_Print('  Negative eigenvalues: '//trim(num2str(n_negative))//' / '//trim(num2str(2*M)),'diag')
       call MIO_Print('  Zero eigenvalues: '//trim(num2str(n_zero))//' / '//trim(num2str(2*M)),'diag')
       call MIO_Print('  ELoc array size: '//trim(num2str(size(ELoc)))//' (will store '//trim(num2str(min(2*M, size(ELoc))))//' eigenvalues)','diag')

       ! Show first and last few eigenvalues
       if (2*M >= 5) then
          call MIO_Print('  First 5 eigenvalues (lowest energy): ['//trim(num2str(eigvals(1),6))//', '//trim(num2str(eigvals(2),6))//', '//trim(num2str(eigvals(3),6))//', '//trim(num2str(eigvals(4),6))//', '//trim(num2str(eigvals(5),6))//']','diag')
       end if
       if (2*M >= 5) then
          call MIO_Print('  Last 5 eigenvalues (highest energy): ['//trim(num2str(eigvals(2*M-4),6))//', '//trim(num2str(eigvals(2*M-3),6))//', '//trim(num2str(eigvals(2*M-2),6))//', '//trim(num2str(eigvals(2*M-1),6))//', '//trim(num2str(eigvals(2*M),6))//']','diag')
       end if

       ! Check for truncation issue
       if (2*M > size(ELoc)) then
          call MIO_Print('  WARNING: Truncation will occur! 2M='//trim(num2str(2*M))//' > ELoc size='//trim(num2str(size(ELoc))),'diag')
          call MIO_Print('  Will lose '//trim(num2str(2*M - size(ELoc)))//' eigenvalues (highest energy states)','diag')
          if (size(ELoc)+1 <= 2*M .and. eigvals(size(ELoc)+1) > 0.0_dp) then
             call MIO_Print('  CRITICAL: First lost eigenvalue is POSITIVE: '//trim(num2str(eigvals(size(ELoc)+1),6)),'diag')
          end if
       else
          call MIO_Print('  No truncation: All '//trim(num2str(2*M))//' eigenvalues will be stored in ELoc','diag')
       end if
    else
       ! Initialize to avoid uninitialized variable warnings
       n_positive = 0
       n_negative = 0
       n_zero = 0
    end if

    ! Extract eigenvalues - handle size mismatch gracefully
    ! For Rashba SOC: we have 2M spin-split eigenvalues
    ! For forceBlockTAPW without SOC: we have 2M eigenvalues in degenerate pairs, collapse to M unique ones
    if (forceBlockTAPW .and. .not. RashbaSOCterm) then
       ! No SOC case: collapse degenerate pairs to M unique eigenvalues
       ! Pattern: (λ₁,λ₁), (λ₂,λ₂), ..., so take eigvals(2i-1) for i=1..M
       do i = 1, min(M, size(ELoc))
          ELoc(i) = eigvals(2*i - 1)  ! Take one eigenvalue from each degenerate pair
       end do
       ! Zero out remaining entries if M < size(ELoc)
       if (M < size(ELoc)) then
          ELoc(M+1:size(ELoc)) = 0.0_dp
       end if
       if (tapwDebug) then
          call MIO_Print('TAPW with block H (no SOC): Collapsed '//trim(num2str(2*M))//' eigenvalues to '//trim(num2str(min(M, size(ELoc))))//' unique ones','diag')
       end if
    else
       ! Rashba SOC case: use all 2M eigenvalues (spin-split)
       if (2*M > size(ELoc)) then
          if (tapwDebug) then
             call MIO_Print('Warning: More TAPW states (2M='//trim(num2str(2*M))//') than output array size ('//trim(num2str(size(ELoc)))//')','diag')
             call MIO_Print('Taking first '//trim(num2str(size(ELoc)))//' eigenvalues out of '//trim(num2str(2*M))//' total TAPW states','diag')
             call MIO_Print('  This will only include eigenvalues up to: '//trim(num2str(eigvals(size(ELoc)),6)),'diag')
             if (eigvals(size(ELoc)) < 0.0_dp .and. n_positive > 0) then
                call MIO_Print('  CRITICAL: All positive eigenvalues will be lost due to truncation!','diag')
             end if
          end if
          ELoc(1:size(ELoc)) = eigvals(1:size(ELoc))  ! Take first size(ELoc) eigenvalues
       else
          ELoc(1:2*M) = eigvals(1:2*M)  ! Assign all 2M TAPW eigenvalues
          if (2*M < size(ELoc)) then
             ELoc(2*M+1:size(ELoc)) = 0.0_dp  ! Zero out unused entries
          end if
       end if
    end if

    if (tapwDebug) then
       if (forceBlockTAPW .and. .not. RashbaSOCterm) then
          call MIO_Print('TAPW with block H (no SOC): Using '//trim(num2str(min(M_tapw, size(ELoc))))//' eigenvalues for output (M_tapw='//trim(num2str(M_tapw))//', ELoc size='//trim(num2str(size(ELoc)))//')','diag')
       else
          call MIO_Print('TAPW with block H: Using '//trim(num2str(min(2*M, size(ELoc))))//' eigenvalues for output (2M='//trim(num2str(2*M))//', ELoc size='//trim(num2str(size(ELoc)))//')','diag')
       end if
    end if

    ! Store TAPW data for Berry curvature calculation (SOC-aware with 2N×2M X matrix)
    ! For Rashba SOC, we need to store both spin channels from the spin-mixed eigenvalues
    if (calculateChern .and. present(kpoint_index)) then
       ! Initialize storage arrays if this is the first k-point (same as DiagH0TAPW)
       if (.not. allocated(stored_eigenvectors)) then
          stored_M = 2*M  ! For Rashba: 2M eigenvalues/eigenvectors (spin-mixed)
          stored_N = N
          stored_nspin = ns  ! Store number of spin channels
          ! Calculate actual number of k-points needed based on Chern grid
          stored_nk = nk_chern_x * nk_chern_y
          call MIO_Print('Allocating TAPW storage for '//trim(num2str(stored_nk))//' k-points ('//trim(num2str(nk_chern_x))//'x'//trim(num2str(nk_chern_y))//' grid) for '//trim(num2str(stored_nspin))//' spin channel(s)','diag')
          allocate(stored_eigenvectors(2*M, 2*M, stored_nk, stored_nspin))
          allocate(stored_hamiltonians(2*M, 2*M, stored_nk, stored_nspin))

          ! Store X matrix for Berry curvature calculation (TB parameters are module-level)
          allocate(stored_X_matrix(2*N, 2*M))
          stored_X_matrix = XArray

          call MIO_Print('Allocated TAPW storage for Chern calculation: '//trim(num2str(2*M))//'x'//trim(num2str(2*M))//'x'//trim(num2str(stored_nk))//'x'//trim(num2str(stored_nspin)),'diag')
          call MIO_Print('Stored SOC-aware X matrix (2N='//trim(num2str(2*N))//'x 2M='//trim(num2str(2*M))//') for Berry curvature calculation','diag')

          ! VERIFICATION: Check stored data integrity
          call MIO_Print('=== TAPW STORAGE VERIFICATION (Rashba SOC) ===','diag')
          call MIO_Print('Stored dimensions: N='//trim(num2str(stored_N))//', M='//trim(num2str(stored_M))//', nk='//trim(num2str(stored_nk))//', nspin='//trim(num2str(stored_nspin)),'diag')
          call MIO_Print('XArray dimensions: '//trim(num2str(size(XArray,1)))//'x'//trim(num2str(size(XArray,2))),'diag')
          call MIO_Print('XArray norm: '//trim(num2str(sqrt(sum(abs(XArray)**2)),8)),'diag')
          call MIO_Print('XArray max element: '//trim(num2str(maxval(abs(XArray)),8)),'diag')
          call MIO_Print('XArray min element: '//trim(num2str(minval(abs(XArray)),8)),'diag')
       else
          ! Check for consistency - M should not change between k-points
          if (2*M /= stored_M) then
             call MIO_Print('ERROR: M dimension changed between k-points in DiagH0TAPW_withBlockH','diag')
             call MIO_Print('  Stored M: '//trim(num2str(stored_M))//', Current M: '//trim(num2str(2*M)),'diag')
             error stop 1
          end if

          ! Verify X matrix dimensions are consistent
          if (.not. allocated(stored_X_matrix)) then
             call MIO_Print('ERROR: stored_X_matrix not allocated but stored_eigenvectors is allocated','diag')
             error stop 1
          end if
          if (size(stored_X_matrix, 1) /= 2*N .or. size(stored_X_matrix, 2) /= 2*M) then
             call MIO_Print('ERROR: Stored X matrix dimensions inconsistent with current calculation','diag')
             call MIO_Print('  Stored: '//trim(num2str(size(stored_X_matrix,1)))//'x'//trim(num2str(size(stored_X_matrix,2))),'diag')
             call MIO_Print('  Current: '//trim(num2str(2*N))//'x'//trim(num2str(2*M)),'diag')
          end if
       end if

       ! Store eigenvectors and hamiltonians for both spin channels
       ! For Rashba with block-diagonal projection, we now have 2M eigenvalues/eigenvectors (spin-mixed)
       ! Store all 2M states for both spin channels (they're spin-mixed from the block diagonalization)
       if (kpoint_index <= stored_nk) then
          ! Store all 2M eigenvalues/eigenvectors (but only up to array size)
          ! For spin channel 1
          if (2*M <= size(stored_eigenvectors,1) .and. 2*M <= size(stored_eigenvectors,2) .and. stored_nspin >= 1) then
             stored_eigenvectors(1:2*M,1:2*M,kpoint_index,1) = Hproj(1:2*M,1:2*M)  ! ZHEEV returns eigenvectors in Hproj
             stored_hamiltonians(1:2*M,1:2*M,kpoint_index,1) = Hproj_copy(1:2*M,1:2*M)  ! Store original Hamiltonian

             ! CRITICAL: Store TAPW eigenvalues (now 2M)
             if (.not. allocated(stored_eigenvalues)) then
                allocate(stored_eigenvalues(2*M, stored_nk, stored_nspin))
                call MIO_Print('Allocated stored_eigenvalues array: '//trim(num2str(2*M))//'x'//trim(num2str(stored_nk))//'x'//trim(num2str(stored_nspin)),'diag')
             end if
             stored_eigenvalues(1:2*M,kpoint_index,1) = eigvals(1:2*M)

             ! For Rashba SOC, spin 2 uses the same spin-mixed states
             ! In practice, both spin channels see the same eigenvalues/eigenvectors
             if (stored_nspin >= 2) then
                stored_eigenvectors(1:2*M,1:2*M,kpoint_index,2) = Hproj(1:2*M,1:2*M)
                stored_hamiltonians(1:2*M,1:2*M,kpoint_index,2) = Hproj_copy(1:2*M,1:2*M)
                stored_eigenvalues(1:2*M,kpoint_index,2) = eigvals(1:2*M)
             end if

             ! DEBUG: Verify eigenvalue storage for first k-point
             if (kpoint_index == 1 .and. is == 1) then
                call MIO_Print('=== TAPW EIGENVALUE STORAGE VERIFICATION (Rashba, k-point 1, spin-mixed) ===','diag')
                call MIO_Print('Stored eigenvalues dimensions: '//trim(num2str(size(stored_eigenvalues,1)))//'x'//trim(num2str(size(stored_eigenvalues,2)))//'x'//trim(num2str(size(stored_eigenvalues,3))),'diag')
                call MIO_Print('First 5 TAPW eigenvalues (spin-mixed): ['//trim(num2str(eigvals(1),6))//','//trim(num2str(eigvals(2),6))//','//trim(num2str(eigvals(3),6))//','//trim(num2str(eigvals(4),6))//','//trim(num2str(eigvals(5),6))//']','diag')
                call MIO_Print('TAPW eigenvalue range: ['//trim(num2str(minval(eigvals(1:2*M)),6))//','//trim(num2str(maxval(eigvals(1:2*M)),6))//']','diag')
             end if
          end if
       end if

       deallocate(Hproj_copy)
    end if

    ! Clean up
    deallocate(XArray, Hproj, eigvals, eigvec, ZWorkLoc, DWorkLoc, Gx, Gy, label)

    if (tapwDebug) call MIO_Print('DiagH0TAPW_withBlockH: Completed successfully', 'diag')

end subroutine DiagH0TAPW_withBlockH

!> @brief Transform block Hamiltonian using TAPW projection
!! @details Projects 2N×2N block Hamiltonian using TAPW X matrix
!!
!! The block Hamiltonian has structure:
!!   HBlock = [ H_up      H_Rashba   ]
!!            [ H_Rashba*  H_dn      ]
!!
!! The TAPW X matrix is 2N×M where:
!!   - Rows 1 to N: X_up (spin-up TAPW basis)
!!   - Rows N+1 to 2N: X_up (spin-down TAPW basis, same as spin-up)
!!
!! The projection Hproj = X† * HBlock * X yields:
!!   Hproj = X_up†*(H_up + H_dn)*X_up + 2*X_up†*H_Rashba*X_up
!!
!! This preserves the Rashba spin-flip coupling. The factor of 2 comes from
!! the two cross-terms (↑→↓ and ↓→↑) in the block structure.
!!
!! @param[in]     N         Block size (2N for SOC, where N is number of atoms)
!! @param[in]     M         TAPW matrix size
!! @param[in]     X         TAPW projection matrix (2N×M)
!! @param[in]     HBlock    Block Hamiltonian (2N×2N)
!! @param[out]    Hproj     Projected Hamiltonian (M×M)
subroutine transform_block_hamiltonian_tapw(N, M, X, HBlock, Hproj)
    use constants, only : cmplx_i, cmplx_0, cmplx_1
    implicit none

    integer, intent(in) :: N, M
    complex(dp), intent(in) :: X(N, M)
    complex(dp), intent(in) :: HBlock(N, N)
    complex(dp), intent(out) :: Hproj(M, M)

    complex(dp), allocatable :: temp(:,:)
    integer :: i, j

    ! Allocate temporary matrix
    ! For Hproj = X† * HBlock * X where X is (N×M) and HBlock is (N×N)
    ! Step 1: temp = HBlock * X where HBlock is (N×N), X is (N×M), temp is (N×M)
    allocate(temp(N, M))

    ! Project block Hamiltonian: Hproj = X† * HBlock * X
    ! X is (N×M), HBlock is (N×N), so:
    ! First: temp = HBlock * X  where HBlock is (N×N), X is (N×M), temp is (N×M)
    call ZGEMM('N', 'N', N, M, N, cmplx_1, HBlock, N, X, N, cmplx_0, temp, N)

    ! Second: Hproj = X† * temp  where X† is (M×N), temp is (N×M), Hproj is (M×M)
    call ZGEMM('C', 'N', M, M, N, cmplx_1, X, N, temp, N, cmplx_0, Hproj, M)

    deallocate(temp)

end subroutine transform_block_hamiltonian_tapw

subroutine generate_shifted_G_list_reduced(rcell, k_ref, NGrange, Gx, Gy, NG, center_point, apply_BZ_filter)
  implicit none
  double precision, intent(in) :: rcell(3,3), k_ref(2)
  integer, intent(in) :: NGrange
  double precision, allocatable, intent(out) :: Gx(:), Gy(:)
  integer, intent(out) :: NG
  double precision, intent(in), optional :: center_point(2)  ! Optional center point for G-grid
  logical, intent(in), optional :: apply_BZ_filter  ! Optional hexagonal BZ filtering

  ! All variable declarations must be at the top for Intel Fortran
  double precision :: b1(2), b2(2), det, inv_rcell(2,2), G_nn(2)
  double precision :: ucell_T(2,2)  ! Direct lattice transpose for coordinate conversion
  double precision :: G_nn_method1(2), G_nn_method2(2), frac_roundtrip1(2), frac_roundtrip2(2)  ! For verification
  integer :: G_nn_int(2), G_nn_int_transpose(2), i, j, count, nmax, u, v
  integer :: idx, temp_idx, ii
  integer, allocatable :: sorted_indices(:)
  double precision, allocatable :: sorted_distances(:)
  double precision, allocatable :: sel_dist(:), sel_ang(:)
  integer, allocatable :: perm(:)
  double precision, allocatable :: Gpoints(:,:), distances(:)
  integer, allocatable :: inverse(:)
  double precision, parameter :: tol = 1.0e-5_dp
  integer :: num_unique, k
  double precision, allocatable :: unique_distances(:)
  logical :: found_distance
  double precision :: current_dist, temp_dist
  double precision :: rotation_angle, cos_rot, sin_rot, Gx_rot, Gy_rot
  integer :: initial_NG
  double precision :: grid_center(2)  ! Center point for G-grid generation
  logical :: use_BZ_filter  ! Whether to apply hexagonal BZ filtering
  integer :: NG_before_filter, NG_after_filter  ! For counting filtered G-vectors

  ! Handle optional parameters
  if (present(center_point)) then
     grid_center = center_point
  else
     grid_center = k_ref  ! Default: center around k_ref
  end if

  if (present(apply_BZ_filter)) then
     use_BZ_filter = apply_BZ_filter
  else
     use_BZ_filter = .false.  ! Default: no BZ filtering
  end if

  ! Get reciprocal lattice vectors (moiré supercell)
  b1 = rcell(1:2,1); b2 = rcell(1:2,2)

  ! Get direct lattice vectors from reciprocal lattice
  ! rcell = 2π * inv(ucell.T), so ucell.T = 2π * inv(rcell)
  det = b1(1)*b2(2) - b1(2)*b2(1)
  ucell_T = reshape([b2(2), -b2(1), -b1(2), b1(1)], [2,2]) * (2.0_dp * 3.14159265358979323846_dp) / det

  G_nn_method1 = matmul(grid_center, ucell_T) / (2.0_dp * 3.14159265358979323846_dp)
  G_nn_int(1) = nint(G_nn_method1(1))
  G_nn_int(2) = nint(G_nn_method1(2))

  G_nn = G_nn_int(1)*b1 + G_nn_int(2)*b2

  ! Generate (2*NGrange+1)^2 grid around G_nn - exactly like Python
  nmax = (2*NGrange + 1)**2
  allocate(Gpoints(2, nmax), distances(nmax), inverse(nmax))

  count = 0
  ! CRITICAL FIX: Match Python meshgrid ordering exactly
  ! Python: U,V = np.meshgrid(arange(G_nn[0]-N_G, G_nn[0]+N_G+1), arange(G_nn[1]-N_G, G_nn[1]+N_G+1))
  ! Python: U.flatten(), V.flatten() - this orders by varying V first, then U
  do v = -NGrange, NGrange  ! SWAPPED: v outer loop (like Python meshgrid)
     do u = -NGrange, NGrange  ! SWAPPED: u inner loop
        count = count + 1
        ! Calculate G-point: Gpoint = (G_nn_int + [u,v]) * [b1, b2]
        Gpoints(1, count) = (G_nn_int(1) + u)*b1(1) + (G_nn_int(2) + v)*b2(1)
        Gpoints(2, count) = (G_nn_int(1) + u)*b1(2) + (G_nn_int(2) + v)*b2(2)
        distances(count) = sqrt((Gpoints(1, count) - G_nn(1))**2 + (Gpoints(2, count) - G_nn(2))**2)
     end do
  end do

  ! Step 1: Round distances to 5 decimal places (equivalent to Python np.round(..., 5))
  do i = 1, count
     distances(i) = nint(distances(i) * 1.0e5_dp) / 1.0e5_dp
  end do

  ! Step 2: Find unique distances and sort them (like Python np.unique does by default)
  allocate(unique_distances(count))
  num_unique = 0

  ! First pass: collect unique distances in order they appear
  do i = 1, count
     current_dist = distances(i)
     found_distance = .false.

     ! Check if this distance already exists in unique_distances
     do j = 1, num_unique
        if (abs(unique_distances(j) - current_dist) < tol) then
           found_distance = .true.
           exit
        end if
     end do

     ! If not found, add to unique distances
     if (.not. found_distance) then
        num_unique = num_unique + 1
        unique_distances(num_unique) = current_dist
     end if
  end do

  ! Sort unique distances (ascending order, like Python np.unique)
  do i = 1, num_unique-1
     do j = i+1, num_unique
        if (unique_distances(i) > unique_distances(j)) then
           temp_dist = unique_distances(i)
           unique_distances(i) = unique_distances(j)
           unique_distances(j) = temp_dist
        end if
     end do
  end do

  ! Second pass: create inverse mapping based on sorted unique distances
  do i = 1, count
     current_dist = distances(i)
     ! Find which sorted unique distance this point belongs to
     do j = 1, num_unique
        if (abs(unique_distances(j) - current_dist) < tol) then
           inverse(i) = j - 1  ! 0-based indexing like Python
           exit
        end if
     end do
  end do

  ! Step 3: Count how many G-points have inverse < NGrange (equivalent to Python inverse < N_G)
  NG = 0
  do i = 1, count
     if (inverse(i) < NGrange) then
        NG = NG + 1
     end if
  end do

  ! Step 4: Allocate output arrays
  allocate(Gx(NG), Gy(NG))

  ! Step 5: Store selected G-points IN ORIGINAL MESHGRID ORDER
  idx = 0
  do i = 1, count
     if (inverse(i) < NGrange) then
        idx = idx + 1
        Gx(idx) = Gpoints(1, i)
        Gy(idx) = Gpoints(2, i)
     end if
  end do

  ! Output G-vectors for matplotlib visualization (thread-safe)
  ! Use critical section to avoid conflicts in parallel execution
  !$OMP CRITICAL(g_vectors_file)
  close(98, status='keep')
  open(unit=98, file='g_vectors_debug.dat', status='replace')
  write(98, '(A)') '# Gx, Gy (for matplotlib plotting)'
  do i = 1, NG
     write(98, '(2F16.8)') Gx(i), Gy(i)
  end do
  close(98)
  !$OMP END CRITICAL(g_vectors_file)

  ! Cleanup
  deallocate(Gpoints, distances, unique_distances, inverse)

end subroutine generate_shifted_G_list_reduced

!-----------------------------------------------------------------------------
! Diag.TAPWBothValleys: append to the G-list of valley refF the same shells
! around the other graphene K point, so one TAPW basis holds both valleys and
! X^H H X carries the intervalley block exactly.  On return G 1..NG1 are the
! refF valley (tapw_NGvalley1 = NG1), G NG1+1..NG the other one.
! The two blocks are independent columns only if no G of one block equals a G
! of the other modulo the graphene reciprocal lattice (then exp(iG.r) coincide
! on every graphene site); that holds when both shells stay inside |K|/2,
! the valley-purity cap, and is checked here pair by pair.
!-----------------------------------------------------------------------------
subroutine tapw_append_second_valley(rcell, rG, sGlattice, refF, NGrange, useDense, Gx, Gy, NG)
  use ham, only : IsingSOCterm
  use constants, only : pi
  implicit none
  real(dp), intent(in) :: rcell(3,3), rG(2,2), sGlattice(2,2), refF(2)
  integer, intent(in) :: NGrange
  logical, intent(in) :: useDense
  real(dp), allocatable, intent(inout) :: Gx(:), Gy(:)
  integer, intent(inout) :: NG
  real(dp), allocatable :: Gx2(:), Gy2(:), tx(:), ty(:)
  real(dp) :: k1(2), k2(2), f(2), d1, d2, kcap
  integer :: NG1, NG2, i, j, nalias

  if (useTriangularTruncation) call MIO_Kill('Diag.TAPWBothValleys: not implemented with '// &
       'Diag.TriangularTruncation (its triangle is valley-specific)','diag','tapw_append_second_valley')
  if (checkTAPWUnitary) call MIO_Kill('Diag.TAPWBothValleys: Diag.CheckTAPWUnitary already spans '// &
       'the whole BZ','diag','tapw_append_second_valley')
  if (.not. useDense) call MIO_Kill('Diag.TAPWBothValleys needs useDenseMatrixTAPW .true.', &
       'diag','tapw_append_second_valley')
  if (IsingSOCterm) call MIO_Kill('Diag.TAPWBothValleys: IsingSOCterm is applied as a valley sign '// &
       'fixed by Diag.UseKprimeValley','diag','tapw_append_second_valley')

  k1 = refF(1)*rG(:,1) + refF(2)*rG(:,2)
  k2 = refF(2)*rG(:,1) + refF(1)*rG(:,2)          ! (2/3,1/3) <-> (1/3,2/3)
  call generate_shifted_G_list_reduced(rcell, k2, NGrange, Gx2, Gy2, NG2)

  NG1 = NG
  allocate(tx(NG1+NG2), ty(NG1+NG2))
  tx(1:NG1) = Gx(1:NG1);  tx(NG1+1:) = Gx2
  ty(1:NG1) = Gy(1:NG1);  ty(NG1+1:) = Gy2
  call move_alloc(tx, Gx)
  call move_alloc(ty, Gy)
  NG = NG1 + NG2
  tapw_NGvalley1 = NG1

  nalias = 0
  do i = 1, NG1
     do j = NG1+1, NG
        f(1) = ((Gx(i)-Gx(j))*sGlattice(1,1) + (Gy(i)-Gy(j))*sGlattice(2,1))/(2.0_dp*pi)
        f(2) = ((Gx(i)-Gx(j))*sGlattice(1,2) + (Gy(i)-Gy(j))*sGlattice(2,2))/(2.0_dp*pi)
        if (all(abs(f - nint(f)) < 1.0e-6_dp)) nalias = nalias + 1
     end do
  end do
  if (nalias > 0) then
     call MIO_Print('Diag.TAPWBothValleys: '//trim(num2str(nalias))//' G pairs of the two valleys '// &
          'coincide modulo the graphene reciprocal lattice - lower Diag.N_G','diag')
     call MIO_Kill('Diag.TAPWBothValleys: the two valley shells overlap','diag','tapw_append_second_valley')
  end if

  !$OMP CRITICAL(tapw_both_valleys)
  if (.not. tbv_announced) then
     tbv_announced = .true.
     kcap = 0.5_dp*sqrt(sum(k1**2))
     d1 = sqrt(maxval((Gx(1:NG1)-k1(1))**2 + (Gy(1:NG1)-k1(2))**2))
     d2 = sqrt(maxval((Gx(NG1+1:NG)-k2(1))**2 + (Gy(NG1+1:NG)-k2(2))**2))
     call MIO_Print('TAPW BOTH VALLEYS: G 1..'//trim(num2str(NG1))//' around refF ('// &
          trim(num2str(refF(1),4))//','//trim(num2str(refF(2),4))//'), G '//trim(num2str(NG1+1))//'..'// &
          trim(num2str(NG))//' around the other valley','diag')
     call MIO_Print('  max|G-K| = '//trim(num2str(d1,5))//', '//trim(num2str(d2,5))// &
          ' 1/Ang; valley-purity cap |K|/2 = '//trim(num2str(kcap,5))//'; aliased pairs: 0','diag')
     if (max(d1,d2) >= kcap) call MIO_Print('  WARNING: a shell reaches past |K|/2','diag')
  end if
  ! generate_shifted_G_list_reduced rewrites this file at every k: overwrite it with the full list
  open(unit=98, file='g_vectors_debug.dat', status='replace')
  write(98, '(A)') '# Gx, Gy, valley block (1 = refF, 2 = other)'
  do i = 1, NG
     write(98, '(2F16.8,I3)') Gx(i), Gy(i), merge(1, 2, i <= NG1)
  end do
  close(98)
  !$OMP END CRITICAL(tapw_both_valleys)
  deallocate(Gx2, Gy2)
end subroutine tapw_append_second_valley

subroutine generate_shifted_G_list(rcell, k_ref, NGrange, Gx, Gy, NG, center_point, apply_BZ_filter)
  implicit none
  double precision, intent(in) :: rcell(3,3), k_ref(2)
  integer, intent(in) :: NGrange
  double precision, allocatable, intent(out) :: Gx(:), Gy(:)
  integer, intent(out) :: NG
  double precision, intent(in), optional :: center_point(2)  ! Optional center point for G-grid
  logical, intent(in), optional :: apply_BZ_filter  ! Optional hexagonal BZ filtering

  ! All variable declarations must be at the top for Intel Fortran
  double precision :: b1(2), b2(2), det, inv_rcell(2,2), G_nn(2)
  double precision :: ucell_T(2,2)  ! Direct lattice transpose for coordinate conversion
  double precision :: G_nn_method1(2), G_nn_method2(2), frac_roundtrip1(2), frac_roundtrip2(2)  ! For verification
  integer :: G_nn_int(2), G_nn_int_transpose(2), i, j, count, nmax, u, v
  integer :: idx, temp_idx, ii
  integer, allocatable :: sorted_indices(:)
  double precision, allocatable :: sorted_distances(:)
  double precision, allocatable :: sel_dist(:), sel_ang(:)
  integer, allocatable :: perm(:)
  double precision, allocatable :: Gpoints(:,:), distances(:)
  integer, allocatable :: inverse(:)
  double precision, parameter :: tol = 1.0e-5_dp
  integer :: num_unique, k
  double precision, allocatable :: unique_distances(:)
  logical :: found_distance
  double precision :: current_dist, temp_dist
  double precision :: rotation_angle, cos_rot, sin_rot, Gx_rot, Gy_rot
  integer :: initial_NG
  double precision :: grid_center(2)  ! Center point for G-grid generation
  logical :: use_BZ_filter  ! Whether to apply hexagonal BZ filtering
  integer :: NG_before_filter, NG_after_filter  ! For counting filtered G-vectors

  ! Get reciprocal lattice vectors (moiré supercell)
  b1 = rcell(1:2,1); b2 = rcell(1:2,2)

  ! FIXED: Use direct lattice for coordinate conversion (like Python)
  ! Python: G_nn = point @ position.Tvector.T / (2*PI)
  ! We need ucell (direct lattice) not rcell (reciprocal lattice)

  ! Determine the center point for G-grid generation
  if (present(center_point)) then
     grid_center = center_point
     call MIO_Print('Using custom center point for G-grid: ['//&
                    trim(num2str(grid_center(1),6))//','//trim(num2str(grid_center(2),6))//']','diag')
  else
     grid_center = k_ref
     call MIO_Print('Using k_ref as center point for G-grid: ['//&
                    trim(num2str(grid_center(1),6))//','//trim(num2str(grid_center(2),6))//']','diag')
  end if

  ! Determine if hexagonal BZ filtering should be applied
  use_BZ_filter = .false.
  if (present(apply_BZ_filter)) then
     use_BZ_filter = apply_BZ_filter
     if (use_BZ_filter) then
        call MIO_Print('HEXAGONAL BZ FILTERING ENABLED: Will exclude G-vectors outside first BZ','diag')
     end if
  end if

  ! Get direct lattice vectors from reciprocal lattice
  ! rcell = 2π * inv(ucell.T), so ucell.T = 2π * inv(rcell)
  det = b1(1)*b2(2) - b1(2)*b2(1)
  ucell_T = reshape([b2(2), -b2(1), -b1(2), b1(1)], [2,2]) * (2.0_dp * 3.14159265358979323846_dp) / det

  ! Calculate inv_rcell for old debug prints (keeping for reference)
  inv_rcell = reshape([b2(2), -b2(1), -b1(2), b1(1)], [2,2]) / det

  ! Debug: Check matrix calculations
  print *, "DEBUG: OLD vs NEW coordinate conversion methods:"
  print *, "  det =", det
  print *, "  OLD inv_rcell method (reciprocal lattice inverse):"
  print *, "    inv_rcell(1,1) = b2(2)/det =", b2(2), "/", det, "=", b2(2)/det
  print *, "    inv_rcell(1,2) = -b2(1)/det =", -b2(1), "/", det, "=", -b2(1)/det
  print *, "    inv_rcell(2,1) = -b1(2)/det =", -b1(2), "/", det, "=", -b1(2)/det
  print *, "    inv_rcell(2,2) = b1(1)/det =", b1(1), "/", det, "=", b1(1)/det

  ! Debug: Direct lattice calculation details
  print *, "DEBUG: Direct lattice coordinate conversion (FIXED):"
  print *, "  ucell_T (direct lattice transpose):"
  print *, "    ucell_T(1,:) =", ucell_T(1,:)
  print *, "    ucell_T(2,:) =", ucell_T(2,:)

  ! Debug: Matrix multiplication k_ref @ ucell_T / (2π) (Python equivalent)
  print *, "DEBUG: k_ref @ ucell_T / (2π) calculation (Python equivalent):"
  print *, "  k_ref(1) * ucell_T(1,1) + k_ref(2) * ucell_T(2,1) ="
  print *, "  ", k_ref(1), "*", ucell_T(1,1), "+", k_ref(2), "*", ucell_T(2,1), "="
  print *, "  ", (k_ref(1) * ucell_T(1,1) + k_ref(2) * ucell_T(2,1)) / (2.0_dp * 3.14159265358979323846_dp)
  print *, "  k_ref(1) * ucell_T(1,2) + k_ref(2) * ucell_T(2,2) ="
  print *, "  ", k_ref(1), "*", ucell_T(1,2), "+", k_ref(2), "*", ucell_T(2,2), "="
  print *, "  ", (k_ref(1) * ucell_T(1,2) + k_ref(2) * ucell_T(2,2)) / (2.0_dp * 3.14159265358979323846_dp)

  ! CRITICAL CHECK: Test both transpose methods for coordinate conversion
  print *, "DEBUG: Testing coordinate conversion methods:"
  print *, "  Method 1 (current): k_ref * ucell_T =", matmul(k_ref, ucell_T) / (2.0_dp * 3.14159265358979323846_dp)
  print *, "  Method 2 (transpose): k_ref * ucell_T^T =", matmul(k_ref, transpose(ucell_T)) / (2.0_dp * 3.14159265358979323846_dp)

  ! Use standard Fortran nint() for rounding to nearest G-grid center
  ! For triangular lattice (1/3, 2/3 values), nint() gives same results as Python np.round()
  G_nn_method1 = matmul(grid_center, ucell_T) / (2.0_dp * 3.14159265358979323846_dp)
  G_nn_int(1) = nint(G_nn_method1(1))
  G_nn_int(2) = nint(G_nn_method1(2))
  call MIO_Print('Using standard Fortran nint() for G-grid center selection','diag')
  print *, "DEBUG: Fractional coordinates:", G_nn_method1
  print *, "DEBUG: nint()-rounded G_nn_int =", G_nn_int

  ! Also calculate what standard Fortran nint() would give for comparison
  G_nn_int_transpose(1) = nint(G_nn_method1(1))
  G_nn_int_transpose(2) = nint(G_nn_method1(2))
  print *, "DEBUG: Standard nint() would give:", G_nn_int_transpose
  if (any(G_nn_int /= G_nn_int_transpose)) then
     call MIO_Print('WARNING: Python and Fortran rounding differ - using Python for consistency','diag')
  end if

  ! COMPREHENSIVE VERIFICATION: Test both methods by reconstruction
  print *, "DEBUG: === COMPREHENSIVE COORDINATE CONVERSION VERIFICATION ==="

  ! Method 1 reconstruction: G_nn_int -> G_nn -> k_space
  G_nn_method1 = G_nn_int(1)*b1 + G_nn_int(2)*b2
  print *, "Method 1 reconstruction:"
  print *, "  G_nn_int =", G_nn_int
  print *, "  G_nn_method1 = G_nn_int(1)*b1 + G_nn_int(2)*b2 =", G_nn_method1
  print *, "  Distance from grid_center: |G_nn_method1 - grid_center| =", sqrt(sum((G_nn_method1 - grid_center)**2))

  ! Method 2 reconstruction: G_nn_int_transpose -> G_nn -> k_space
  G_nn_method2 = G_nn_int_transpose(1)*b1 + G_nn_int_transpose(2)*b2
  print *, "Method 2 reconstruction:"
  print *, "  G_nn_int_transpose =", G_nn_int_transpose
  print *, "  G_nn_method2 = G_nn_int_transpose(1)*b1 + G_nn_int_transpose(2)*b2 =", G_nn_method2
  print *, "  Distance from grid_center: |G_nn_method2 - grid_center| =", sqrt(sum((G_nn_method2 - grid_center)**2))

  ! Round-trip test: fractional -> integer -> fractional
  frac_roundtrip1 = matmul(G_nn_method1, ucell_T) / (2.0_dp * 3.14159265358979323846_dp)
  frac_roundtrip2 = matmul(G_nn_method2, ucell_T) / (2.0_dp * 3.14159265358979323846_dp)
  print *, "Round-trip test (should recover integer values):"
  print *, "  Method 1: G_nn_method1 -> fractional =", frac_roundtrip1
  print *, "  Method 2: G_nn_method2 -> fractional =", frac_roundtrip2

  ! The correct method should give the nearest G-vector to grid_center
  print *, "VERDICT: Method with smaller distance from grid_center is likely correct"
  if (sqrt(sum((G_nn_method1 - grid_center)**2)) < sqrt(sum((G_nn_method2 - grid_center)**2))) then
print *, "  -> Method 1 (current) appears CORRECT: closer to grid_center"
  else
     print *, "  -> Method 2 (transpose) appears CORRECT: closer to grid_center"
  end if

  G_nn = G_nn_int(1)*b1 + G_nn_int(2)*b2
  print *, "DEBUG: G_nn reconstruction:"
  print *, "  G_nn_int(1)*b1 =", G_nn_int(1), "*", b1, "=", G_nn_int(1)*b1
  print *, "  G_nn_int(2)*b2 =", G_nn_int(2), "*", b2, "=", G_nn_int(2)*b2
  print *, "  G_nn = G_nn_int(1)*b1 + G_nn_int(2)*b2 =", G_nn

  ! Generate (2*NGrange+1)^2 grid around G_nn - exactly like Python
  nmax = (2*NGrange + 1)**2
  allocate(Gpoints(2, nmax), distances(nmax), inverse(nmax))

  count = 0
  ! CRITICAL FIX: Match Python meshgrid ordering exactly
  ! Python: U,V = np.meshgrid(arange(G_nn[0]-N_G, G_nn[0]+N_G+1), arange(G_nn[1]-N_G, G_nn[1]+N_G+1))
  ! Python: U.flatten(), V.flatten() - this orders by varying V first, then U
  do v = -NGrange, NGrange  ! SWAPPED: v outer loop (like Python meshgrid)
     do u = -NGrange, NGrange  ! SWAPPED: u inner loop
        count = count + 1
        ! Calculate G-point: Gpoint = (G_nn_int + [u,v]) * [b1, b2]
        Gpoints(1, count) = (G_nn_int(1) + u)*b1(1) + (G_nn_int(2) + v)*b2(1)
        Gpoints(2, count) = (G_nn_int(1) + u)*b1(2) + (G_nn_int(2) + v)*b2(2)
        distances(count) = sqrt((Gpoints(1, count) - G_nn(1))**2 + (Gpoints(2, count) - G_nn(2))**2)

        ! DEBUG: Check for the center point
        if (u == 0 .and. v == 0) then
           if (tapwDebug) call MIO_Print('DEBUG: Found center point u=0,v=0:','diag')
           call MIO_Print('  Position in array: '//trim(adjustl(num2str(real(count,dp),0))),'diag')
           call MIO_Print('  G_point = ['//trim(adjustl(num2str(Gpoints(1,count),6)))//','//&
                          trim(adjustl(num2str(Gpoints(2,count),6)))//']','diag')
           call MIO_Print('  G_nn = ['//trim(adjustl(num2str(G_nn(1),6)))//','//&
                          trim(adjustl(num2str(G_nn(2),6)))//']','diag')
           call MIO_Print('  Distance = '//trim(adjustl(num2str(distances(count),6))),'diag')
        end if
     end do
  end do

  ! Find unique distances and create inverse mapping - Python np.unique equivalent
  allocate(unique_distances(nmax))
  num_unique = 0

  if (tapwDebug) call MIO_Print('DEBUG: First 15 raw distances from G_nn:','diag')
  do i = 1, min(15, count)
     call MIO_Print('  Point '//trim(adjustl(num2str(real(i,dp),0)))//': distance = '//trim(adjustl(num2str(distances(i),6))),'diag')
  end do

  ! Also check around the center position
  if (tapwDebug) call MIO_Print('DEBUG: Looking for distance=0 in all raw distances:','diag')
  do i = 1, count
     if (abs(distances(i)) < 1.0e-10_dp) then
        call MIO_Print('  Found distance=0 at position '//trim(adjustl(num2str(real(i,dp),0))),'diag')
     end if
  end do

  do i = 1, count
     current_dist = distances(i)
     found_distance = .false.

     ! DEBUG: Track the center point specifically
     if (i == 25) then
        if (tapwDebug) call MIO_Print('DEBUG: Processing center point (position 25):','diag')
        call MIO_Print('  Original distance = '//trim(adjustl(num2str(current_dist,10))),'diag')
     end if

     ! Check if this distance already exists (rounded to 5 decimal places like Python)
     current_dist = nint(current_dist * 1.0e5_dp) / 1.0e5_dp  ! Round to 5 decimal places

     if (i == 25) then
        call MIO_Print('  Rounded distance = '//trim(adjustl(num2str(current_dist,10))),'diag')
     end if

     do j = 1, num_unique
        if (abs(current_dist - unique_distances(j)) < 1.0e-10_dp) then  ! Exact match after rounding
           found_distance = .true.
           inverse(i) = j - 1  ! 0-based indexing like Python
           if (i == 25) then
              call MIO_Print('  Found existing shell: '//trim(adjustl(num2str(real(j-1,dp),0))),'diag')
           end if
           exit
        end if
     end do

     ! If new distance, add it to unique list
     if (.not. found_distance) then
        num_unique = num_unique + 1
        unique_distances(num_unique) = current_dist
        inverse(i) = num_unique - 1  ! 0-based indexing like Python
        if (i == 25) then
           call MIO_Print('  Created new shell: '//trim(adjustl(num2str(real(num_unique-1,dp),0)))//&
                         ' with distance = '//trim(adjustl(num2str(current_dist,10))),'diag')
        end if
     end if
  end do

  ! Sort unique distances to get proper shell ordering
  if (tapwDebug) call MIO_Print('DEBUG: Before sorting - first 5 unique distances:','diag')
  do i = 1, min(5, num_unique)
     call MIO_Print('  Shell '//trim(adjustl(num2str(real(i-1,dp),0)))//' distance = '//trim(adjustl(num2str(unique_distances(i),6))),'diag')
  end do

  call sort_distances_with_inverse_update(unique_distances, num_unique, inverse, count)

  if (tapwDebug) call MIO_Print('DEBUG: After sorting - first 5 unique distances:','diag')
  do i = 1, min(5, num_unique)
     call MIO_Print('  Shell '//trim(adjustl(num2str(real(i-1,dp),0)))//' distance = '//trim(adjustl(num2str(unique_distances(i),6))),'diag')
  end do

  ! Count how many G-vectors we need: inverse < NGrange (Python condition)
  NG = 0
  if (tapwDebug) call MIO_Print('DEBUG: Checking inverse array for selection:','diag')
  do i = 1, count
     if (i <= 10 .or. i == 25) then  ! Show first 10 + center point
        call MIO_Print('  Point '//trim(adjustl(num2str(real(i,dp),0)))//': inverse = '//&
                       trim(adjustl(num2str(real(inverse(i),dp),0)))//', NGrange = '//&
                       trim(adjustl(num2str(real(NGrange,dp),0)))//', selected = '//&
                       merge('YES','NO ',inverse(i) < NGrange),'diag')
     end if
     if (inverse(i) < NGrange) then
        NG = NG + 1
     end if
  end do

  if (tapwDebug) call MIO_Print('DEBUG: G-vector selection summary','diag')
  call MIO_Print('  Total grid points generated: '//trim(adjustl(num2str(real(count,dp),0))),'diag')
  call MIO_Print('  Unique distance shells: '//trim(adjustl(num2str(real(num_unique,dp),0))),'diag')
  call MIO_Print('  NGrange (shells to keep): '//trim(adjustl(num2str(real(NGrange,dp),0))),'diag')
  call MIO_Print('  Selected G-vectors: '//trim(adjustl(num2str(real(NG,dp),0))),'diag')

  if (present(center_point)) then
     ! Unitary TAPW: use current shell-based logic
     if (tapwDebug) call MIO_Print('DEBUG: Using shell-based G-vector selection for unitary TAPW','diag')

     ! Allocate output arrays
     allocate(Gx(NG), Gy(NG))

     ! Copy all G-vectors where inverse < NGrange (current logic)
     k = 0
     NG_before_filter = 0
     do i = 1, count
        if (inverse(i) < NGrange) then
           NG_before_filter = NG_before_filter + 1

           ! Apply hexagonal BZ filter if enabled
           if (use_BZ_filter) then
              if (.not. is_inside_hexagonal_BZ(Gpoints(1, i), Gpoints(2, i), b1, b2)) then
                 ! Skip this G-vector - it's outside the hexagonal BZ
                 cycle
              end if
           end if

           k = k + 1
           Gx(k) = Gpoints(1, i)
           Gy(k) = Gpoints(2, i)
        end if
     end do

     NG_after_filter = k
     NG = NG_after_filter

  else
     ! Normal TAPW: use Python-like direct distance sorting
     if (tapwDebug) call MIO_Print('DEBUG: Using Python-like distance-based G-vector selection for normal TAPW','diag')

     ! Create arrays of indices and distances for sorting
     allocate(sorted_indices(count), sorted_distances(count))

     ! Initialize indices and copy distances
     do i = 1, count
        sorted_indices(i) = i
        sorted_distances(i) = distances(i)
     end do

     ! Sort by distance (simple bubble sort for now)
     do i = 1, count-1
        do j = i+1, count
           if (sorted_distances(i) > sorted_distances(j)) then
              ! Swap distances
              temp_dist = sorted_distances(i)
              sorted_distances(i) = sorted_distances(j)
              sorted_distances(j) = temp_dist
              ! Swap indices
              temp_idx = sorted_indices(i)
              sorted_indices(i) = sorted_indices(j)
              sorted_indices(j) = temp_idx
           end if
        end do
     end do

     ! Select first NG G-vectors (closest to reference point)
     ! NG was calculated earlier based on shell logic, but now we select directly
     allocate(Gx(NG), Gy(NG))

     k = 0
     do i = 1, count
        if (k >= NG) exit  ! Stop when we have enough G-vectors

        idx = sorted_indices(i)
        k = k + 1
        Gx(k) = Gpoints(1, idx)
        Gy(k) = Gpoints(2, idx)

        if (k <= 5) then  ! Debug first 5
           call MIO_Print('  Selected G('//trim(adjustl(num2str(real(k,dp),0)))//') distance='//&
                         trim(adjustl(num2str(sorted_distances(i),6)))//' = ['//&
                         trim(adjustl(num2str(Gx(k),6)))//','//trim(adjustl(num2str(Gy(k),6)))//']','diag')
        end if
     end do

     deallocate(sorted_indices, sorted_distances)
  end if

  ! Report filtering results
  if (use_BZ_filter) then
     call MIO_Print('HEXAGONAL BZ FILTERING RESULTS:','diag')
     call MIO_Print('  G-vectors before BZ filter: '//trim(adjustl(num2str(real(NG_before_filter,dp),0))),'diag')
     call MIO_Print('  G-vectors after BZ filter: '//trim(adjustl(num2str(real(NG_after_filter,dp),0))),'diag')
     call MIO_Print('  G-vectors removed: '//trim(adjustl(num2str(real(NG_before_filter-NG_after_filter,dp),0))),'diag')
     call MIO_Print('  Efficiency gain: '//trim(adjustl(num2str(100.0_dp*(NG_before_filter-NG_after_filter)/NG_before_filter,1)))//'%','diag')
  end if

  ! No need to reallocate - arrays are already the right size and contain correct data
  ! The filtering was done in-place, so Gx(1:NG) and Gy(1:NG) contain the correct G-vectors

  ! CRITICAL: Apply final uniqueness like Python's np.unique(np.vstack(self.G_list),axis=0)
  ! This removes duplicate G-vectors that might come from multiple reference points
  if (tapwDebug) call MIO_Print('DEBUG: About to call remove_duplicate_G_vectors with NG='//trim(adjustl(num2str(real(NG,dp),0))),'diag')
  initial_NG = NG
  call remove_duplicate_G_vectors(Gx, Gy, NG)
  if (tapwDebug) call MIO_Print('DEBUG: After remove_duplicate_G_vectors, NG='//trim(adjustl(num2str(real(NG,dp),0))),'diag')
  if (NG /= initial_NG) then
     call MIO_Print('G-vector uniqueness: '//trim(adjustl(num2str(real(initial_NG,dp),0)))//&
                    ' → '//trim(adjustl(num2str(real(NG,dp),0)))//' (removed '//&
                    trim(adjustl(num2str(real(initial_NG-NG,dp),0)))//' duplicates)','diag')
  end if

  ! === FINAL: enforce Python-like order ===
  ! Sort by distance to G_nn first, then by angle for tie-breaking (matches Python behavior)
  if (tapwDebug) call MIO_Print('DEBUG: Applying final canonical ordering (distance then angle)','diag')

  allocate(sel_dist(NG), sel_ang(NG), perm(NG))
  do ii = 1, NG
     sel_dist(ii) = sqrt((Gx(ii) - G_nn(1))**2 + (Gy(ii) - G_nn(2))**2)
     sel_ang(ii) = atan2(Gy(ii) - G_nn(2), Gx(ii) - G_nn(1))
     perm(ii) = ii
  end do

  call argsort_by_dist_then_angle(sel_dist, sel_ang, perm, NG)
  call apply_permutation_in_place(Gx, Gy, perm, NG)

  if (tapwDebug) call MIO_Print('DEBUG: First 5 distances after final sort:','diag')
  do ii = 1, min(5, NG)
     call MIO_Print('  '//trim(adjustl(num2str(sel_dist(perm(ii)),6)))//'  ang='//trim(adjustl(num2str(sel_ang(perm(ii)),6))),'diag')
  end do

  deallocate(sel_dist, sel_ang, perm)
  ! === END FINAL ORDER ===

  ! Apply G-grid rotation if specified (configurable via Diag.GGridRotationAngle)
  rotation_angle = gGridRotationAngle * 3.14159265358979323846_dp / 180.0_dp
  if (abs(gGridRotationAngle) > 1.0e-10_dp) then
     cos_rot = cos(rotation_angle)
     sin_rot = sin(rotation_angle)
     call MIO_Print("G-grid rotation: "//trim(adjustl(num2str(gGridRotationAngle,3)))//" degrees",'diag')

     do i = 1, NG
        Gx_rot = Gx(i)
        Gy_rot = Gy(i)
        ! Apply rotation matrix: [cos -sin; sin cos] * [Gx; Gy]
        Gx(i) = cos_rot * Gx_rot - sin_rot * Gy_rot
        Gy(i) = sin_rot * Gx_rot + cos_rot * Gy_rot
  end do
  end if

  deallocate(Gpoints, distances, inverse, unique_distances)

  ! Essential info: G-vector generation summary
  call MIO_Print('G-grid: Generated NG='//trim(adjustl(num2str(real(NG,dp),0)))//' vectors, NGrange='//trim(adjustl(num2str(real(NGrange,dp),0))),'diag')

  ! Output G-vectors for matplotlib visualization (thread-safe)
  ! Use critical section to avoid conflicts in parallel execution
  !$OMP CRITICAL(g_vectors_file)
  close(98, status='keep')
  open(unit=98, file='g_vectors_debug.dat', status='replace')
  write(98, '(A)') '# Gx, Gy (for matplotlib plotting)'
  do i = 1, NG
     write(98, '(2F16.8)') Gx(i), Gy(i)
  end do
  close(98)
  !$OMP END CRITICAL(g_vectors_file)

end subroutine generate_shifted_G_list

! Helper subroutine to sort distances and update inverse mapping to match Python np.unique
subroutine sort_distances_with_inverse_update(distances, n_unique, inverse, n_total)
  implicit none
  integer, intent(in) :: n_unique, n_total
  double precision, intent(inout) :: distances(n_unique)
  integer, intent(inout) :: inverse(n_total)

  integer :: i, j, min_idx, old_shell_idx
  double precision :: temp_dist, old_dist
  double precision, allocatable :: old_distances(:)
  integer, allocatable :: old_inverse(:)

  ! COMPLETELY REWRITTEN: Simple and bulletproof approach
  ! Save original arrays
  allocate(old_distances(n_unique), old_inverse(n_total))
  old_distances = distances
  old_inverse = inverse

  ! Simple selection sort of distances
  do i = 1, n_unique-1
     min_idx = i
     do j = i+1, n_unique
        if (distances(j) < distances(min_idx)) then
           min_idx = j
        end if
     end do

     if (min_idx /= i) then
        ! Swap distances
        temp_dist = distances(i)
        distances(i) = distances(min_idx)
        distances(min_idx) = temp_dist
     end if
  end do

  ! Rebuild inverse mapping: for each point, find which new shell its distance belongs to
  do i = 1, n_total
     old_shell_idx = old_inverse(i) + 1  ! Convert to 1-based
     old_dist = old_distances(old_shell_idx)  ! Get the original distance

     ! Find this distance in the new sorted array
     do j = 1, n_unique
        if (abs(distances(j) - old_dist) < 1.0e-12_dp) then
           inverse(i) = j - 1  ! Convert to 0-based
           exit
        end if
     end do
  end do

  deallocate(old_distances, old_inverse)
end subroutine sort_distances_with_inverse_update

! Remove duplicate G-vectors to match Python's final np.unique step
subroutine remove_duplicate_G_vectors(Gx, Gy, NG)
  implicit none
  integer, intent(inout) :: NG
  double precision, intent(inout) :: Gx(NG), Gy(NG)

  integer :: i, j, unique_count
  double precision, allocatable :: unique_Gx(:), unique_Gy(:)
  logical :: is_duplicate
  double precision, parameter :: tol = 1.0e-10_dp

  allocate(unique_Gx(NG), unique_Gy(NG))
  unique_count = 0

  do i = 1, NG
     is_duplicate = .false.
     do j = 1, unique_count
        if (abs(Gx(i) - unique_Gx(j)) < tol .and. abs(Gy(i) - unique_Gy(j)) < tol) then
           is_duplicate = .true.
           exit
        end if
     end do

     if (.not. is_duplicate) then
        unique_count = unique_count + 1
        unique_Gx(unique_count) = Gx(i)
        unique_Gy(unique_count) = Gy(i)
     end if
  end do

  ! Copy back unique vectors
  NG = unique_count
  Gx(1:NG) = unique_Gx(1:NG)
  Gy(1:NG) = unique_Gy(1:NG)

  deallocate(unique_Gx, unique_Gy)
end subroutine remove_duplicate_G_vectors

! Sort permutation array by distance first, then angle for tie-breaking
subroutine argsort_by_dist_then_angle(dist, ang, perm, n)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: dist(n), ang(n)
  integer, intent(inout) :: perm(n)
  integer :: i, j, best, tmpi
  ! Simple selection sort: stable enough for modest NG
  do i = 1, n-1
     best = i
     do j = i+1, n
        if (dist(perm(j)) < dist(perm(best))) then
           best = j
        else if (abs(dist(perm(j)) - dist(perm(best))) <= 1.0d-12) then
           if (ang(perm(j)) < ang(perm(best))) best = j
        end if
     end do
     if (best /= i) then
        tmpi = perm(i)
        perm(i) = perm(best)
        perm(best) = tmpi
     end if
  end do
end subroutine argsort_by_dist_then_angle

! Apply permutation to reorder Gx, Gy arrays
subroutine apply_permutation_in_place(Gx, Gy, perm, n)
  implicit none
  integer, intent(in) :: n
  double precision, intent(inout) :: Gx(n), Gy(n)
  integer, intent(in) :: perm(n)
  double precision, allocatable :: tmpx(:), tmpy(:)
  integer :: i
  allocate(tmpx(n), tmpy(n))
  do i = 1, n
     tmpx(i) = Gx(perm(i))
     tmpy(i) = Gy(perm(i))
  end do
  Gx = tmpx
  Gy = tmpy
  deallocate(tmpx, tmpy)
end subroutine apply_permutation_in_place

logical function is_inside_hexagonal_BZ(Gx, Gy, b1, b2)
    ! Check if G-vector (Gx, Gy) is inside the hexagonal first Brillouin Zone
    ! defined by reciprocal lattice vectors b1, b2
    implicit none
    double precision, intent(in) :: Gx, Gy, b1(2), b2(2)

    ! Local variables
    double precision :: G(2)
    double precision :: dot_products(6)
    double precision :: BZ_normals(6,2)
    double precision :: BZ_distances(6)
    integer :: i

    G = [Gx, Gy]

    ! Define the 6 face normals of the hexagonal BZ
    ! For a hexagon with vertices at ±b1, ±b2, ±(b1+b2), the face normals are:
    BZ_normals(1,:) = b1 / sqrt(dot_product(b1, b1))           ! +b1 direction
    BZ_normals(2,:) = -b1 / sqrt(dot_product(b1, b1))          ! -b1 direction
    BZ_normals(3,:) = b2 / sqrt(dot_product(b2, b2))           ! +b2 direction
    BZ_normals(4,:) = -b2 / sqrt(dot_product(b2, b2))          ! -b2 direction
    BZ_normals(5,:) = (b1 + b2) / sqrt(dot_product(b1 + b2, b1 + b2))   ! +(b1+b2) direction
    BZ_normals(6,:) = -(b1 + b2) / sqrt(dot_product(b1 + b2, b1 + b2))  ! -(b1+b2) direction

    ! Distance from origin to each face (half the distance between opposite faces)
    BZ_distances(1) = 0.5 * sqrt(dot_product(b1, b1))
    BZ_distances(2) = 0.5 * sqrt(dot_product(b1, b1))
    BZ_distances(3) = 0.5 * sqrt(dot_product(b2, b2))
    BZ_distances(4) = 0.5 * sqrt(dot_product(b2, b2))
    BZ_distances(5) = 0.5 * sqrt(dot_product(b1 + b2, b1 + b2))
    BZ_distances(6) = 0.5 * sqrt(dot_product(b1 + b2, b1 + b2))

    ! Check if G is inside all 6 half-spaces (inside hexagon)
    is_inside_hexagonal_BZ = .true.
    do i = 1, 6
        dot_products(i) = dot_product(G, BZ_normals(i,:))
        if (dot_products(i) > BZ_distances(i)) then
            is_inside_hexagonal_BZ = .false.
            return
        end if
    end do

end function is_inside_hexagonal_BZ

subroutine generate_shifted_G_list_with_graphene_BZ(rcell, k_ref, NGrange, Gx, Gy, NG, center_point, rG_graphene)
  ! Generate G-vectors with hexagonal BZ filtering using graphene reciprocal lattice
  implicit none
  double precision, intent(in) :: rcell(3,3), k_ref(2)
  integer, intent(in) :: NGrange
  double precision, allocatable, intent(out) :: Gx(:), Gy(:)
  integer, intent(out) :: NG
  double precision, intent(in) :: center_point(2)  ! Center point for G-grid
  double precision, intent(in) :: rG_graphene(2,2)  ! Graphene reciprocal lattice for BZ filtering

  ! Local variables - same as generate_shifted_G_list but with graphene BZ filtering
  double precision :: b1_moire(2), b2_moire(2), b1_graphene(2), b2_graphene(2)
  integer :: NG_before_filter, NG_after_filter, k, i

  call MIO_Print('GRAPHENE BZ FILTERING: Using graphene reciprocal lattice for BZ boundaries','diag')

  ! Get moiré reciprocal lattice vectors for G-vector generation
  b1_moire = rcell(1:2,1); b2_moire = rcell(1:2,2)

  ! Get graphene reciprocal lattice vectors for BZ filtering
  b1_graphene = rG_graphene(:,1); b2_graphene = rG_graphene(:,2)

  call MIO_Print('Moiré lattice vectors (for G-grid): b1=['//&
                 trim(num2str(b1_moire(1),4))//','//trim(num2str(b1_moire(2),4))//&
                 '], b2=['//trim(num2str(b2_moire(1),4))//','//trim(num2str(b2_moire(2),4))//']','diag')
  call MIO_Print('Graphene lattice vectors (for BZ filter): b1=['//&
                 trim(num2str(b1_graphene(1),4))//','//trim(num2str(b1_graphene(2),4))//&
                 '], b2=['//trim(num2str(b2_graphene(1),4))//','//trim(num2str(b2_graphene(2),4))//']','diag')

  ! First, generate G-vectors using moiré lattice (no filtering)
  call generate_shifted_G_list(rcell, k_ref, NGrange, Gx, Gy, NG, center_point)
  NG_before_filter = NG

  ! Now apply graphene BZ filtering
  k = 0
  do i = 1, NG_before_filter
     if (is_inside_hexagonal_BZ(Gx(i), Gy(i), b1_graphene, b2_graphene)) then
        k = k + 1
        Gx(k) = Gx(i)
        Gy(k) = Gy(i)
     end if
  end do

  NG_after_filter = k
  NG = NG_after_filter

  ! Report filtering results
  call MIO_Print('GRAPHENE BZ FILTERING RESULTS:','diag')
  call MIO_Print('  G-vectors before graphene BZ filter: '//trim(adjustl(num2str(real(NG_before_filter,dp),0))),'diag')
  call MIO_Print('  G-vectors after graphene BZ filter: '//trim(adjustl(num2str(real(NG_after_filter,dp),0))),'diag')
  call MIO_Print('  G-vectors removed: '//trim(adjustl(num2str(real(NG_before_filter-NG_after_filter,dp),0))),'diag')
  if (NG_before_filter > 0) then
     call MIO_Print('  Efficiency gain: '//trim(adjustl(num2str(100.0_dp*(NG_before_filter-NG_after_filter)/NG_before_filter,1)))//'%','diag')
  end if

  ! Output FILTERED G-vectors for matplotlib visualization (overwrites the unfiltered version)
  call MIO_Print('Writing FILTERED G-vectors to g_vectors_debug.dat','diag')
  ! Use critical section to avoid conflicts in parallel execution
  !$OMP CRITICAL(g_vectors_file)
  close(98, status='keep')
  open(unit=98, file='g_vectors_debug.dat', status='replace')
  write(98, '(A)') '# Filtered Gx, Gy (for matplotlib plotting)'
  do i = 1, NG
     write(98, '(2F16.8)') Gx(i), Gy(i)
  end do
  close(98)
  !$OMP END CRITICAL(g_vectors_file)

end subroutine generate_shifted_G_list_with_graphene_BZ

subroutine generate_triangular_G_list(rcell, k_ref, NGrange, Gx, Gy, NG, rG, use_distance_ordering)
  ! Generate G-vectors using triangular truncation following Python get_Gvecs_tri exactly
  ! This implements the same logic as the Python script to avoid K/K' valley overlap
  implicit none

  double precision, intent(in) :: rcell(3,3), k_ref(2)
  integer, intent(in) :: NGrange
  double precision, allocatable, intent(out) :: Gx(:), Gy(:)
  integer, intent(out) :: NG
  double precision, intent(in) :: rG(2,2)  ! Graphene reciprocal lattice vectors
  logical, intent(in), optional :: use_distance_ordering  ! Flag to enable/disable distance-based ordering

  ! Local variables following Python script structure
  integer :: N1, N2, i, j, idx, total_G, selected_count
  double precision, allocatable :: G_frac(:,:)
  double precision :: ratio, kkDvec_ind(2)
  double precision :: BZ_triangle(3,2)  ! Triangle corners in fractional coordinates
  logical, allocatable :: G_selected(:)
  double precision :: bMvec(2,2)  ! Moiré reciprocal lattice vectors

  ! Variables for distance-based sorting
  double precision, allocatable :: temp_Gx(:), temp_Gy(:), distances(:)
  integer, allocatable :: sorted_indices(:)
  double precision :: temp_dist, dx, dy
  integer :: temp_idx
  logical :: do_distance_sort  ! Local flag to control sorting

  ! Parameters following Python script
  ratio = 1.0_dp  ! Full triangle (corresponds to Python ratio parameter)

  ! Set distance sorting flag (default to true if not provided)
  if (present(use_distance_ordering)) then
     do_distance_sort = use_distance_ordering
  else
     do_distance_sort = .true.  ! Default: use distance-based ordering
  end if

  ! CRITICAL FIX: Convert k_ref from Cartesian to fractional coordinates in moiré lattice
  ! k_ref is the graphene K-point in Cartesian coordinates
  ! We need to express it in fractional coordinates of the moiré lattice for grid selection
  call MIO_Print('k_ref (Cartesian): ['//trim(num2str(k_ref(1),6))//','//trim(num2str(k_ref(2),6))//']','diag')
  call cart_to_frac_single_point(k_ref, bMvec, kkDvec_ind)
  call MIO_Print('k_ref converted to fractional: ['//trim(num2str(kkDvec_ind(1),6))//','//trim(num2str(kkDvec_ind(2),6))//']','diag')

  ! Set grid size following Python nsize logic
  ! NGrange now controls the density of the G-vector grid (like nsize in Python)
  N1 = NGrange  ! Direct correspondence to Python nsize[0]
  N2 = NGrange  ! Direct correspondence to Python nsize[1]

  call MIO_Print('=== Python-style triangular G-vector generation ===','diag')
  call MIO_Print('Grid size: ['//trim(num2str(N1))//','//trim(num2str(N2))//']','diag')
  call MIO_Print('Valley position (frac): ['//trim(num2str(kkDvec_ind(1),6))//','//trim(num2str(kkDvec_ind(2),6))//']','diag')
  call MIO_Print('Triangle ratio: '//trim(num2str(ratio,3)),'diag')

  ! Python: G_frac = np.meshgrid(np.fft.ifftshift(np.arange(-N2,N2)), np.fft.ifftshift(np.arange(-N1,N1)))
  ! np.fft.ifftshift for arange(-N,N) gives: [0, 1, 2, ..., N-1, -N, -N+1, ..., -1]
  total_G = (2*N2) * (2*N1)
  allocate(G_frac(2, total_G))
  allocate(G_selected(total_G))

  ! Generate G-vector grid exactly like Python: np.meshgrid + np.fft.ifftshift
  idx = 0
  ! Python meshgrid: first array varies along columns (j), second along rows (i)
  do i = -N2, N2-1  ! Corresponds to np.fft.ifftshift(np.arange(-N2,N2))
     do j = -N1, N1-1  ! Corresponds to np.fft.ifftshift(np.arange(-N1,N1))
        idx = idx + 1
        ! Python: G_frac = np.vstack(((G_frac[1]).flatten(),(G_frac[0]).flatten()))
        ! This means: G_frac[0] = meshgrid[1] (j values), G_frac[1] = meshgrid[0] (i values)
        G_frac(1, idx) = real(j, dp)  ! j corresponds to first component
        G_frac(2, idx) = real(i, dp)  ! i corresponds to second component
     end do
  end do

  call MIO_Print('Generated '//trim(num2str(total_G))//' G-vectors in fractional coordinates','diag')
  call MIO_Print('G-grid range: ['//trim(num2str(minval(G_frac(1,:)),2))//','//trim(num2str(minval(G_frac(2,:)),2))//'] to ['// &
                 trim(num2str(maxval(G_frac(1,:)),2))//','//trim(num2str(maxval(G_frac(2,:)),2))//']','diag')

  ! Set up moiré reciprocal lattice vectors (bMvec in Python)
  ! rcell uses column storage: rcell(:,i) = i-th vector
  bMvec(:,1) = rcell(1:2, 1)  ! First moiré reciprocal vector
  bMvec(:,2) = rcell(1:2, 2)  ! Second moiré reciprocal vector

  call MIO_Print('Moiré reciprocal lattice vectors:','diag')
  call MIO_Print('  b1 = ['//trim(num2str(bMvec(1,1),6))//','//trim(num2str(bMvec(2,1),6))//']','diag')
  call MIO_Print('  b2 = ['//trim(num2str(bMvec(1,2),6))//','//trim(num2str(bMvec(2,2),6))//']','diag')

  ! CORRECT APPROACH: Use alternating corners of graphene BZ, shifted to K-point
  ! This follows your exact instructions and uses the proper graphene BZ geometry
  call create_triangle_from_graphene_BZ(k_ref, ratio, BZ_triangle, rG, bMvec)

  ! Select G-vectors inside triangular boundary (Python: G_bounds.contains_point)
  call select_G_vectors_in_triangle(G_frac, total_G, BZ_triangle, G_selected)

  ! Count selected G-vectors
  selected_count = count(G_selected)

  call MIO_Print('Triangular selection results:','diag')
  call MIO_Print('  Total G-vectors: '//trim(num2str(total_G)),'diag')
  call MIO_Print('  Selected G-vectors: '//trim(num2str(selected_count)),'diag')
  call MIO_Print('  Selection efficiency: '//trim(num2str(100.0_dp*selected_count/total_G,1))//'%','diag')

  if (selected_count == 0) then
     call MIO_Print('ERROR: No G-vectors selected by triangular truncation!','diag')
     call MIO_Print('This will cause M=0 and matrix dimension errors.','diag')
     call MIO_Print('Check BZ triangle calculation or increase grid size.','diag')
     error stop 1
  end if

  ! OPTIONAL DISTANCE-BASED ORDERING: Match generate_shifted_G_list_reduced approach
  allocate(Gx(selected_count), Gy(selected_count))

  if (do_distance_sort) then
     call MIO_Print('TRIANGULAR: Using distance-based G-vector ordering (matching hexagonal method)','diag')

     ! Allocate temporary arrays for distance-based sorting
     allocate(temp_Gx(selected_count), temp_Gy(selected_count))
     allocate(distances(selected_count), sorted_indices(selected_count))

     ! First pass: collect selected G-vectors and compute distances
     idx = 0
     do i = 1, total_G
        if (G_selected(i)) then
           idx = idx + 1
           ! Convert fractional to Cartesian using moiré reciprocal lattice
           temp_Gx(idx) = G_frac(1,i) * bMvec(1,1) + G_frac(2,i) * bMvec(1,2)
           temp_Gy(idx) = G_frac(1,i) * bMvec(2,1) + G_frac(2,i) * bMvec(2,2)

           ! Calculate distance from k_ref (same as hexagonal method)
           dx = temp_Gx(idx) - k_ref(1)
           dy = temp_Gy(idx) - k_ref(2)
           distances(idx) = sqrt(dx*dx + dy*dy)
           sorted_indices(idx) = idx
        end if
     end do

     ! Sort by distance using bubble sort (same algorithm as hexagonal method)
     do i = 1, selected_count-1
        do j = i+1, selected_count
           if (distances(i) > distances(j)) then
              ! Swap distances
              temp_dist = distances(i)
              distances(i) = distances(j)
              distances(j) = temp_dist
              ! Swap indices
              temp_idx = sorted_indices(i)
              sorted_indices(i) = sorted_indices(j)
              sorted_indices(j) = temp_idx
           end if
        end do
     end do

     ! Copy G-vectors in distance-sorted order
     do i = 1, selected_count
        idx = sorted_indices(i)
        Gx(i) = temp_Gx(idx)
        Gy(i) = temp_Gy(idx)
     end do

     ! Cleanup temporary arrays
     deallocate(temp_Gx, temp_Gy, distances, sorted_indices)

  else
     call MIO_Print('TRIANGULAR: Using original grid-order G-vector ordering (Python-like)','diag')

     ! Original approach: copy G-vectors in grid generation order
     idx = 0
     do i = 1, total_G
        if (G_selected(i)) then
           idx = idx + 1
           ! Convert fractional to Cartesian using moiré reciprocal lattice
           Gx(idx) = G_frac(1,i) * bMvec(1,1) + G_frac(2,i) * bMvec(1,2)
           Gy(idx) = G_frac(1,i) * bMvec(2,1) + G_frac(2,i) * bMvec(2,2)
        end if
     end do
  end if

  NG = selected_count

  if (do_distance_sort) then
     if (tapwDebug) then
        call MIO_Print('Generated '//trim(num2str(NG))//' G-vectors using triangular truncation with distance-based ordering','diag')
        call MIO_Print('First 5 G-vectors (distance-sorted):','diag')
        do i = 1, min(5, NG)
           ! Recalculate distance for debug output
           dx = Gx(i) - k_ref(1)
           dy = Gy(i) - k_ref(2)
           temp_dist = sqrt(dx*dx + dy*dy)
           call MIO_Print('  G('//trim(num2str(i))//') distance='//&
                         trim(num2str(temp_dist,6))//' = ['//trim(num2str(Gx(i),6))//','//trim(num2str(Gy(i),6))//']','diag')
        end do
     end if
  else
     if (tapwDebug) then
        call MIO_Print('Generated '//trim(num2str(NG))//' G-vectors using triangular truncation with grid-order','diag')
        call MIO_Print('First 5 G-vectors (grid-order):','diag')
        do i = 1, min(5, NG)
           call MIO_Print('  G('//trim(num2str(i))//') = ['//trim(num2str(Gx(i),6))//','//trim(num2str(Gy(i),6))//']','diag')
        end do
     end if
  end if

  ! Output G-vectors for matplotlib visualization
  call MIO_Print('Writing triangular G-vectors to g_vectors_triangular.dat','diag')
  open(unit=199, file='g_vectors_triangular.dat', status='replace')
  write(199, '(A)') '# Triangular G-vectors (Gx, Gy) for matplotlib plotting'
  write(199, '(A,I0)') '# Number of G-vectors: ', NG
  write(199, '(A,2I0)') '# Grid size [N1,N2]: ', N1, N2
  write(199, '(A,2F12.6)') '# Valley position (frac): ', kkDvec_ind
  write(199, '(A,F8.3)') '# Triangle ratio: ', ratio
  write(199, '(A)') '# Triangle corners (fractional):'
  do i = 1, 3
     write(199, '(A,I1,A,2F12.6)') '#   Corner ', i, ': ', BZ_triangle(i, :)
  end do
  write(199, '(A)') '# Format: Gx Gy'
  do i = 1, NG
     write(199, '(2F16.8)') Gx(i), Gy(i)
  end do
  close(199)

  ! Cleanup
  deallocate(G_frac, G_selected)

end subroutine generate_triangular_G_list

subroutine select_G_vectors_in_triangle(G_frac, total_G, BZ_corners, G_selected)
  ! Select G-vectors that lie within the triangular boundary
  ! Simplified point-in-triangle test
  implicit none

  integer, intent(in) :: total_G
  real(dp), intent(in) :: G_frac(2, total_G)
  real(dp), intent(in) :: BZ_corners(3, 2)
  logical, intent(out) :: G_selected(total_G)

  integer :: i
  real(dp) :: point(2)

  ! For each G-vector, test if it's inside the triangle
  do i = 1, total_G
     point(1) = G_frac(1, i)
     point(2) = G_frac(2, i)
     G_selected(i) = point_in_triangle(point, BZ_corners)
  end do

end subroutine select_G_vectors_in_triangle

logical function point_in_triangle(point, triangle)
  ! Test if a point is inside a triangle using barycentric coordinates
  implicit none

  real(dp), intent(in) :: point(2)
  real(dp), intent(in) :: triangle(3, 2)

  real(dp) :: v0(2), v1(2), v2(2)
  real(dp) :: dot00, dot01, dot02, dot11, dot12
  real(dp) :: inv_denom, u, v

  ! Triangle vertices
  v0 = triangle(3, :) - triangle(1, :)  ! C - A
  v1 = triangle(2, :) - triangle(1, :)  ! B - A
  v2 = point - triangle(1, :)           ! P - A

  ! Compute dot products
  dot00 = dot_product(v0, v0)
  dot01 = dot_product(v0, v1)
  dot02 = dot_product(v0, v2)
  dot11 = dot_product(v1, v1)
  dot12 = dot_product(v1, v2)

  ! Compute barycentric coordinates
  inv_denom = 1.0_dp / (dot00 * dot11 - dot01 * dot01)
  u = (dot11 * dot02 - dot01 * dot12) * inv_denom
  v = (dot00 * dot12 - dot01 * dot02) * inv_denom

  ! Check if point is in triangle
  point_in_triangle = (u >= 0.0_dp) .and. (v >= 0.0_dp) .and. (u + v <= 1.0_dp)

end function point_in_triangle

subroutine create_triangle_from_graphene_BZ(k_ref_cart, ratio, triangle_frac, rG, bMvec)
  ! Create triangle from alternating corners of graphene BZ, shifted to K-point
  ! This is the CORRECT approach following your instructions exactly
  implicit none

  real(dp), intent(in) :: k_ref_cart(2)     ! K-point in Cartesian coordinates
  real(dp), intent(in) :: ratio             ! Scaling factor
  real(dp), intent(out) :: triangle_frac(3,2) ! Triangle corners in fractional coordinates
  real(dp), intent(in) :: rG(2,2)           ! Graphene reciprocal lattice vectors
  real(dp), intent(in) :: bMvec(2,2)        ! Moiré reciprocal lattice vectors

  ! Local variables
  real(dp) :: b1(2), b2(2)
  real(dp) :: bz_vertices(6,2)              ! All 6 graphene BZ vertices
  real(dp) :: triangle_cart(3,2)            ! Triangle corners in Cartesian
  integer :: alternating_indices(3)
  integer :: i

  call MIO_Print('Creating triangle from alternating graphene BZ corners (CORRECT approach)','diag')

  ! Use graphene reciprocal lattice vectors (same as in the debug file generation)
  b1 = rG(:,1)
  b2 = rG(:,2)

  ! Calculate 6 graphene BZ vertices (exactly like in brillouin_zones_debug_K1.dat)
  bz_vertices(1,:) = (2.0_dp*b1 + b2) / 3.0_dp     ! K point
  bz_vertices(2,:) = (b1 - b2) / 3.0_dp             ! Next vertex clockwise
  bz_vertices(3,:) = -(b1 + 2.0_dp*b2) / 3.0_dp    ! -K' point
  bz_vertices(4,:) = -(2.0_dp*b1 + b2) / 3.0_dp    ! -K point
  bz_vertices(5,:) = -(b1 - b2) / 3.0_dp            ! Next vertex
  bz_vertices(6,:) = (b1 + 2.0_dp*b2) / 3.0_dp     ! K' point

  call MIO_Print('Graphene BZ vertices (Cartesian):','diag')
  do i = 1, 6
     call MIO_Print('  BZ('//trim(num2str(i))//') = ['//trim(num2str(bz_vertices(i,1),6))//','// &
                    trim(num2str(bz_vertices(i,2),6))//']','diag')
  end do

  ! Select alternating corners to form triangle (1/6 of hexagon)
  if (useKprimeValley) then
     alternating_indices = [1, 3, 5]  ! Select alternating vertices for K' valley
     call MIO_Print('Using K'' valley triangular indices: [1, 3, 5]','diag')
  else
     alternating_indices = [2, 4, 6]  ! Skip vertex 1 (K-point), take alternating ones for K valley
     call MIO_Print('Using K valley triangular indices: [2, 4, 6]','diag')
  end if

  do i = 1, 3
     triangle_cart(i, :) = bz_vertices(alternating_indices(i), :) * ratio
  end do

  ! Shift triangle to be centered at K-point
  do i = 1, 3
     triangle_cart(i, 1) = triangle_cart(i, 1) + k_ref_cart(1)
     triangle_cart(i, 2) = triangle_cart(i, 2) + k_ref_cart(2)
  end do

  call MIO_Print('Triangle corners (Cartesian, shifted to K-point):','diag')
  do i = 1, 3
     call MIO_Print('  Triangle('//trim(num2str(i))//') = ['//trim(num2str(triangle_cart(i,1),6))//','// &
                    trim(num2str(triangle_cart(i,2),6))//']','diag')
  end do

  ! Convert triangle corners to fractional coordinates in moiré lattice
  call cart_to_frac_multiple_points(triangle_cart, 3, bMvec, triangle_frac)

  call MIO_Print('Triangle corners (fractional in moiré lattice):','diag')
  do i = 1, 3
     call MIO_Print('  Triangle_frac('//trim(num2str(i))//') = ['//trim(num2str(triangle_frac(i,1),6))//','// &
                    trim(num2str(triangle_frac(i,2),6))//']','diag')
  end do

end subroutine create_triangle_from_graphene_BZ

subroutine cart_to_frac_coords(cart_coords, basis_vectors, frac_coords)
  ! Convert Cartesian coordinates to fractional coordinates (cart2fracArr equivalent)
  ! Implements the coordinate transformation using matrix inversion
  implicit none

  real(dp), intent(in) :: cart_coords(3,2)    ! Cartesian coordinates
  real(dp), intent(in) :: basis_vectors(2,2)  ! Basis vectors (columns)
  real(dp), intent(out) :: frac_coords(3,2)   ! Fractional coordinates

  real(dp) :: inv_basis(2,2), det
  integer :: i

  ! Calculate inverse of basis vectors matrix
  det = basis_vectors(1,1) * basis_vectors(2,2) - basis_vectors(1,2) * basis_vectors(2,1)

  if (abs(det) < 1.0e-12_dp) then
     call MIO_Print('Warning: Singular basis vectors matrix in cart_to_frac_coords','diag')
     frac_coords = cart_coords  ! Fallback
     return
  end if

  inv_basis(1,1) =  basis_vectors(2,2) / det
  inv_basis(1,2) = -basis_vectors(1,2) / det
  inv_basis(2,1) = -basis_vectors(2,1) / det
  inv_basis(2,2) =  basis_vectors(1,1) / det

  ! Convert each point: frac = inv_basis * cart
  do i = 1, 3
     frac_coords(i, 1) = inv_basis(1,1) * cart_coords(i,1) + inv_basis(1,2) * cart_coords(i,2)
     frac_coords(i, 2) = inv_basis(2,1) * cart_coords(i,1) + inv_basis(2,2) * cart_coords(i,2)
  end do

end subroutine cart_to_frac_coords

subroutine cart_to_frac_multiple_points(cart_points, n_points, basis_vectors, frac_points)
  ! Convert multiple Cartesian points to fractional coordinates
  implicit none

  integer, intent(in) :: n_points
  real(dp), intent(in) :: cart_points(n_points,2)    ! Cartesian coordinates
  real(dp), intent(in) :: basis_vectors(2,2)         ! Basis vectors (columns)
  real(dp), intent(out) :: frac_points(n_points,2)   ! Fractional coordinates

  real(dp) :: inv_basis(2,2), det
  integer :: i

  ! Calculate inverse of basis vectors matrix
  det = basis_vectors(1,1) * basis_vectors(2,2) - basis_vectors(1,2) * basis_vectors(2,1)

  if (abs(det) < 1.0e-12_dp) then
     call MIO_Print('Warning: Singular basis vectors matrix in cart_to_frac_multiple_points','diag')
     frac_points = cart_points  ! Fallback
     return
  end if

  inv_basis(1,1) =  basis_vectors(2,2) / det
  inv_basis(1,2) = -basis_vectors(1,2) / det
  inv_basis(2,1) = -basis_vectors(2,1) / det
  inv_basis(2,2) =  basis_vectors(1,1) / det

  ! Convert each point: frac = inv_basis * cart
  do i = 1, n_points
     frac_points(i, 1) = inv_basis(1,1) * cart_points(i,1) + inv_basis(1,2) * cart_points(i,2)
     frac_points(i, 2) = inv_basis(2,1) * cart_points(i,1) + inv_basis(2,2) * cart_points(i,2)
  end do

end subroutine cart_to_frac_multiple_points

subroutine cart_to_frac_single_point(cart_point, basis_vectors, frac_point)
  ! Convert single Cartesian point to fractional coordinates
  implicit none

  real(dp), intent(in) :: cart_point(2)      ! Cartesian coordinates
  real(dp), intent(in) :: basis_vectors(2,2) ! Basis vectors (columns)
  real(dp), intent(out) :: frac_point(2)     ! Fractional coordinates

  real(dp) :: inv_basis(2,2), det

  ! Calculate inverse of basis vectors matrix
  det = basis_vectors(1,1) * basis_vectors(2,2) - basis_vectors(1,2) * basis_vectors(2,1)

  if (abs(det) < 1.0e-12_dp) then
     call MIO_Print('Warning: Singular basis vectors matrix in cart_to_frac_single_point','diag')
     frac_point = cart_point  ! Fallback
     return
  end if

  inv_basis(1,1) =  basis_vectors(2,2) / det
  inv_basis(1,2) = -basis_vectors(1,2) / det
  inv_basis(2,1) = -basis_vectors(2,1) / det
  inv_basis(2,2) =  basis_vectors(1,1) / det

  ! Convert point: frac = inv_basis * cart
  frac_point(1) = inv_basis(1,1) * cart_point(1) + inv_basis(1,2) * cart_point(2)
  frac_point(2) = inv_basis(2,1) * cart_point(1) + inv_basis(2,2) * cart_point(2)

end subroutine cart_to_frac_single_point

! Subroutine to compute number of unique labels
subroutine compute_unique_labels(labels, N, num_unique)
  implicit none
  integer, intent(in) :: N
  integer, intent(in) :: labels(N)
  integer, intent(out) :: num_unique

  integer :: i, j
  integer, allocatable :: unique_labels(:)
  logical :: found

  allocate(unique_labels(N))  ! Maximum possible unique labels
  num_unique = 0

  do i = 1, N
     found = .false.
     do j = 1, num_unique
        if (labels(i) == unique_labels(j)) then
           found = .true.
           exit
        end if
     end do

     if (.not. found) then
        num_unique = num_unique + 1
        unique_labels(num_unique) = labels(i)
     end if
  end do

  deallocate(unique_labels)

  print *, "Found", num_unique, "unique layer×sublattice combinations"
end subroutine compute_unique_labels

! Subroutine to remap labels to contiguous indices 1..Nlabel
subroutine remap_labels_to_contiguous(labels, N, num_unique)
  implicit none
  integer, intent(in) :: N, num_unique
  integer, intent(inout) :: labels(N)

  integer :: i, j, label_idx
  integer, allocatable :: unique_labels(:), label_map(:)
  logical :: found

  allocate(unique_labels(num_unique), label_map(num_unique))

  ! First pass: collect unique labels
  label_idx = 0
  do i = 1, N
     found = .false.
     do j = 1, label_idx
        if (labels(i) == unique_labels(j)) then
           found = .true.
           exit
        end if
     end do

     if (.not. found) then
        label_idx = label_idx + 1
        unique_labels(label_idx) = labels(i)
        label_map(label_idx) = label_idx  ! Map to contiguous index
     end if
  end do

  ! Second pass: remap all labels to contiguous indices
  do i = 1, N
     do j = 1, num_unique
        if (labels(i) == unique_labels(j)) then
           labels(i) = j  ! Remap to contiguous index 1..num_unique
           exit
        end if
     end do
  end do

  ! Debug output
  print *, "Label mapping (original -> contiguous):"
  do i = 1, num_unique
     print *, "  ", unique_labels(i), "->", i
  end do

  deallocate(unique_labels, label_map)
end subroutine remap_labels_to_contiguous

subroutine transform_sparse_hamiltonian(N, M, row_ptr, col_ind, values, X, Hproj)
  use constants, only : cmplx_1, cmplx_0
  use omp_lib,   only : omp_in_parallel
  implicit none

  integer, intent(in) :: N, M
  integer, intent(in) :: row_ptr(N+1), col_ind(:)
  complex(dp), intent(in) :: values(:)
  complex(dp), intent(in) :: X(N,M)
  complex(dp), intent(out) :: Hproj(M,M)

  complex(dp), allocatable :: Y(:,:)
  integer :: i, k, j, nnz, mcol, nthr
  integer, external :: mkl_get_max_threads
  complex(dp) :: tmp, acc

  ! Validate CSR size
  if (size(row_ptr) /= N+1) then
     print *, "ERROR: row_ptr size mismatch!"
     print *, "ERROR: row_ptr has", size(row_ptr), "entries but should have", N+1
     print *, "ERROR: This will cause array bounds violations in sparse multiplication"
     error stop "row_ptr size mismatch in transform_sparse_hamiltonian"
  endif

  nnz = row_ptr(N+1) - 1
  if (size(values) < nnz) then
     print *, "ERROR: size(values) <", nnz
     error stop "CSR format inconsistency: not enough values"
  endif
  if (maxval(col_ind(1:nnz)) > N) then
     print *, "ERROR: max(col_ind) >", N
     error stop "col_ind contains out-of-bounds indices for X"
  endif

  ! Validate matrix dimensions
  if (size(X,1) /= N) then
     print *, "ERROR: Matrix dimension mismatch!"
     print *, "ERROR: X matrix has", size(X,1), "rows but sparse matrix H has", N, "rows"
     print *, "ERROR: This will cause incorrect matrix multiplication"
     error stop "Matrix dimension mismatch in transform_sparse_hamiltonian"
  endif

  ! Allocate Y = H * X
  allocate(Y(N,M))

  ! Index checks, once, outside the multiplication
  do i = 1, N
     do k = row_ptr(i), row_ptr(i+1) - 1
        j = col_ind(k)
        if (j < 1 .or. j > N) then
           print *, "ERROR: j = col_ind(k) = ", j, " out of bounds at i=", i, " k=", k
           error stop 1
        endif
        if (k < 1 .or. k > size(values)) then
           print *, "ERROR: k=", k, " out of bounds (values size=", size(values), ")"
           error stop 1
        endif
     end do
  end do

  ! Sparse multiplication Y = H * X, one column of X at a time, so that X(:,mcol) and Y(:,mcol) are
  ! contiguous. (The row form Y(i,:) = Y(i,:) + values(k)*X(j,:) walks both N x M arrays with stride N:
  ! 6 h for N = 8.4e5, M = 7e3.) Every Y(i,mcol) sums the same terms in the same k order as before.
  ! Threaded over columns with MKL's thread count, because TAPW runs set OMP_NUM_THREADS=1 (k-loop)
  ! and MKL_NUM_THREADS=ncpus; serial when already inside an active parallel region.
  nthr = max(1, mkl_get_max_threads())
  if (omp_in_parallel()) nthr = 1
  !$OMP PARALLEL DO NUM_THREADS(nthr) DEFAULT(SHARED) PRIVATE(mcol, i, k, acc) SCHEDULE(STATIC)
  do mcol = 1, M
     do i = 1, N
        acc = (0.0_dp, 0.0_dp)
        do k = row_ptr(i), row_ptr(i+1) - 1
           acc = acc + values(k) * X(col_ind(k), mcol)
        end do
        Y(i, mcol) = acc
     end do
  end do
  !$OMP END PARALLEL DO

  ! Debug: Check sparse matrix structure
  print *, "Debug: Sparse matrix info:"
  print *, "Debug: N =", N, "M =", M
  print *, "Debug: row_ptr range:", row_ptr(1), "to", row_ptr(N+1)
  print *, "Debug: col_ind range:", minval(col_ind(1:row_ptr(N+1)-1)), "to", maxval(col_ind(1:row_ptr(N+1)-1))
  print *, "Debug: First few row_ptr entries:", (row_ptr(i), i=1,min(10,N+1))
  print *, "Debug: First few col_ind entries:", (col_ind(i), i=1,min(10,row_ptr(N+1)-1))
  print *, "Debug: First few values:", (values(i), i=1,min(10,row_ptr(N+1)-1))

  ! Debug: Check if Y calculation is reasonable
  print *, "Debug: After Y = H * X calculation:"
  print *, "Debug: Max |values| in sparse matrix:", maxval(abs(values))
  print *, "Debug: Max |X| element used in multiplication:", maxval(abs(X))
  print *, "Debug: Expected max |Y| should be around:", maxval(abs(values)) * maxval(abs(X)) * maxval(row_ptr(2:N+1) - row_ptr(1:N))

  ! Debug: Check intermediate result Y = H * X
  print *, "Intermediate Y matrix size:", size(Y,1), "x", size(Y,2)
  print *, "Max |Y|:", maxval(abs(Y))
  print *, "Min |Y|:", minval(abs(Y))
  print *, "Y diagonal elements (first 10):", (abs(Y(i,i)), i=1,min(10,N))

         ! Compute Hproj = X^† * Y
       ! Debug: Check matrix values before multiplication
       print *, "Debug: Original max values - X:", maxval(abs(X)), "Y:", maxval(abs(Y))

       ! Debug: Check matrices before zgemm
       print *, "Debug: Before zgemm call:"
       print *, "Debug: X matrix shape:", shape(X), "leading dimension:", N
       print *, "Debug: Y matrix shape:", shape(Y), "leading dimension:", N
       print *, "Debug: Hproj matrix shape:", shape(Hproj), "leading dimension:", M
       print *, "Debug: zgemm parameters: 'C', 'N',", M, M, N

       ! Matrix multiplication: Hproj = X^† * Y
       ! First try with transpose 'C' (conjugate transpose)
       call zgemm('C', 'N', M, M, N, cmplx_1, X, N, Y, N, cmplx_0, Hproj, M)

       ! Debug: Check result with transpose
       print *, "Debug: After zgemm('C') call:"
       print *, "Debug: Hproj max/min:", maxval(abs(Hproj)), minval(abs(Hproj))

       ! If that fails, try without transpose 'N' (no transpose)
       if (maxval(abs(Hproj)) == 0.0_dp) then
          print *, "Debug: Transpose failed, trying without transpose..."
          call zgemm('N', 'N', M, M, N, cmplx_1, X, N, Y, N, cmplx_0, Hproj, M)
          print *, "Debug: After zgemm('N') call:"
          print *, "Debug: Hproj max/min:", maxval(abs(Hproj)), minval(abs(Hproj))
       endif

       ! If both zgemm calls fail, try manual matrix multiplication
       if (maxval(abs(Hproj)) == 0.0_dp) then
          print *, "Debug: Both zgemm calls failed, trying manual multiplication..."
          Hproj = (0.0_dp, 0.0_dp)
          do i = 1, M
             do j = 1, M
                do k = 1, N
                   Hproj(i,j) = Hproj(i,j) + conjg(X(k,i)) * Y(k,j)
                end do
             end do
          end do
          print *, "Debug: After manual multiplication:"
          print *, "Debug: Hproj max/min:", maxval(abs(Hproj)), minval(abs(Hproj))
       endif

       ! Debug: Check result immediately after zgemm
       print *, "Debug: After zgemm call:"
       print *, "Debug: Hproj max/min:", maxval(abs(Hproj)), minval(abs(Hproj))
       print *, "Debug: Hproj first few diagonal elements:", (abs(Hproj(i,i)), i=1,min(5,M))

  ! Debug: Check final projected Hamiltonian
  print *, "Final Hproj matrix size:", size(Hproj,1), "x", size(Hproj,2)
  print *, "Max |Hproj|:", maxval(abs(Hproj))
  print *, "Min |Hproj|:", minval(abs(Hproj))
  print *, "Hproj diagonal elements (first 10):", (real(Hproj(i,i)), i=1,min(10,M))

  deallocate(Y)
end subroutine transform_sparse_hamiltonian

!> Write the plane-wave composition of TAPW eigenstates gw_b1..gw_b2 at one k to <prefix>.GWeights.
!! Basis column (iG-1)*Nlabel + label (build_X), so weight(iG) = sum_label |U(col, b)|^2 (sums to 1 per band
!! up to the Loewdin mixing of columns, which is small: eig(X^H X) 0.988..1.007).  The G list (Ang^-1) is written
!! once in the header and must not change between k-points.  Energies in eV (eigvals are in g0).
subroutine GWeightsWrite(M, NG, Nlabel, Gx, Gy, U, ev, ik, kvec)
  use name,  only : prefix
  use tbpar, only : g0
  implicit none
  integer, intent(in) :: M, NG, Nlabel, ik
  real(dp), intent(in) :: Gx(NG), Gy(NG), ev(M), kvec(3)
  complex(dp), intent(in) :: U(M, M)
  real(dp) :: w(NG), wtop(gw_ntop)
  integer :: itop(gw_ntop), b, ig, l, t, j
  if (.not. gw_open) then
     open(newunit=gw_unit, file=trim(prefix)//'.GWeights', status='replace')
     gw_open = .true.
     gw_ng = NG
     write(gw_unit,'(a)') '# TAPW plane-wave composition: weight(G) = sum_label |c_(G,label)|^2 of eigenstate b'
     write(gw_unit,'(a,i0,1x,i0,a,i0,a,i0)') '# bands ', gw_b1, gw_b2, '   top ', gw_ntop, '   Nlabel ', Nlabel
     write(gw_unit,'(a,i0)') '# NG ', NG
     do ig = 1, NG
        write(gw_unit,'(a,i6,2(1x,es20.12))') '# G ', ig, Gx(ig), Gy(ig)
     end do
     if (tapwBothValleys) then
        write(gw_unit,'(a,i0,a)') '# valleys: G 1..', tapw_NGvalley1, ' = refF valley, the rest = the other'
        write(gw_unit,'(a)') '# ik kx ky band E[eV] sum_w w_valley1  then (iG weight) x top, largest first'
     else
        write(gw_unit,'(a)') '# ik kx ky band E[eV] sum_w  then (iG weight) x top, largest first'
     end if
  end if
  if (NG /= gw_ng) call MIO_Print('WARNING: Diag.GWeights: NG changed between k-points','diag')
  do b = max(1, gw_b1), min(M, gw_b2)
     w = 0.0_dp
     do ig = 1, NG
        do l = 1, Nlabel
           w(ig) = w(ig) + abs(U((ig-1)*Nlabel + l, b))**2
        end do
     end do
     wtop = -1.0_dp; itop = 0
     do ig = 1, NG                                   ! insertion into the running top list
        if (w(ig) <= wtop(gw_ntop)) cycle
        t = gw_ntop
        do while (t > 1)
           if (w(ig) <= wtop(t-1)) exit
           t = t - 1
        end do
        do j = gw_ntop, t+1, -1
           wtop(j) = wtop(j-1); itop(j) = itop(j-1)
        end do
        wtop(t) = w(ig); itop(t) = ig
     end do
     if (tapwBothValleys) then
        write(gw_unit,'(i6,2(1x,es16.8),1x,i5,1x,f14.8,2(1x,f10.6),*(1x,i6,1x,f10.7))') ik, kvec(1), kvec(2), b, &
             ev(b)*g0, sum(w), sum(w(1:tapw_NGvalley1)), (itop(j), wtop(j), j = 1, gw_ntop)
     else
        write(gw_unit,'(i6,2(1x,es16.8),1x,i5,1x,f14.8,1x,f10.6,*(1x,i6,1x,f10.7))') ik, kvec(1), kvec(2), b, &
             ev(b)*g0, sum(w), (itop(j), wtop(j), j = 1, gw_ntop)
     end if
  end do
  flush(gw_unit)
end subroutine GWeightsWrite

!> Write the layer / sublattice composition of TAPW eigenstates lw_b1..lw_b2 (0 0 = all M) at one k to
!! <prefix>.LayerWeights.  The weight is taken in ATOM space, w(l,b) = sum_{i: label(i)=l} |psi_i|^2 with
!! psi = X U(:,b), so it stays exact when the Loewdin step mixes the (G,label) columns; sum_l w = 1 when
!! X^H X = 1.  Labels are the TAPW ones, code = 10*layerIndex + Species.  Energies in eV (eigvals are in g0).
subroutine LayerWeightsWrite(N, M, Nlabel, X, U, ev, label, labcode, ik, kvec)
  use name,  only : prefix
  use tbpar, only : g0
  use constants, only : cmplx_0, cmplx_1
  implicit none
  integer, intent(in) :: N, M, Nlabel, ik, label(N), labcode(N)
  real(dp), intent(in) :: ev(M), kvec(3)
  complex(dp), intent(in) :: X(N, M), U(M, M)
  integer, parameter :: nchunk = 32
  complex(dp), allocatable :: psi(:,:)
  real(dp) :: w(Nlabel, nchunk)
  integer :: code(Nlabel), b, b1, b2, bb, nb, i, l
  b1 = lw_b1; b2 = lw_b2
  if (b1 < 1) b1 = 1
  if (b2 < 1 .or. b2 > M) b2 = M
  if (.not. lw_open) then
     open(newunit=lw_unit, file=trim(prefix)//'.LayerWeights', status='replace')
     lw_open = .true.
     code = 0
     do i = 1, N
        if (label(i) >= 1 .and. label(i) <= Nlabel) then
           if (code(label(i)) == 0) code(label(i)) = labcode(i)
        end if
     end do
     write(lw_unit,'(a)') '# TAPW layer/sublattice composition: w(label) = sum_{atoms i of that label} |psi_i|^2, psi = X c'
     write(lw_unit,'(a,i0,1x,i0,a,i0,a,i0)') '# bands ', b1, b2, '   M ', M, '   Nlabel ', Nlabel
     write(lw_unit,'(a,*(1x,i0))') '# label codes (10*layer + species):', code
     write(lw_unit,'(a)') '# ik kx ky band E[eV] sum_w  then w(label) in the order of the codes above'
     call MIO_Print('Diag.LayerWeights: layer/sublattice weights of bands '//trim(num2str(b1))//'..'// &
          trim(num2str(b2))//' -> '//trim(prefix)//'.LayerWeights','diag')
  end if
  allocate(psi(N, nchunk))
  do b = b1, b2, nchunk
     nb = min(nchunk, b2 - b + 1)
     call ZGEMM('N', 'N', N, nb, M, cmplx_1, X, N, U(1, b), M, cmplx_0, psi, N)
     w = 0.0_dp
     do bb = 1, nb
        do i = 1, N
           l = label(i)
           if (l >= 1 .and. l <= Nlabel) w(l, bb) = w(l, bb) + real(psi(i, bb))**2 + aimag(psi(i, bb))**2
        end do
     end do
     do bb = 1, nb
        write(lw_unit,'(i6,2(1x,es16.8),1x,i5,1x,f14.8,1x,f10.7,*(1x,f10.7))') ik, kvec(1), kvec(2), b + bb - 1, &
             ev(b + bb - 1)*g0, sum(w(:, bb)), w(:, bb)
     end do
  end do
  deallocate(psi)
  flush(lw_unit)
end subroutine LayerWeightsWrite

!> Write the projected TAPW Hamiltonian H(k) (M x M, before ZHEEV, in eV) to <prefix>.TAPWHam.
!! Basis column (iG-1)*Nlabel + label (build_X).  The header gives the G list (Ang^-1) and, per remapped label,
!! its original code 10*layerIndex + Species (11/12 graphene A/B, 23/24 hBN B/N for GBNtwoLayers), so the
!! reader can downfold labels without guessing.  Per k: a line 'k ik kx ky', then M*M lines 'Re Im' in
!! column-major order (the full matrix: ZHEEV reads only the upper triangle, the reader should too).
subroutine GWeightsHamWrite(M, NG, Nlabel, N, Gx, Gy, H, label, labcode, ik, kvec)
  use name,  only : prefix
  use tbpar, only : g0
  implicit none
  integer, intent(in) :: M, NG, Nlabel, N, ik, label(N), labcode(N)
  real(dp), intent(in) :: Gx(NG), Gy(NG), kvec(3)
  complex(dp), intent(in) :: H(M, M)
  integer :: ig, l, i, j, code(Nlabel)
  if (.not. gwh_open) then
     open(newunit=gwh_unit, file=trim(prefix)//'.TAPWHam', status='replace')
     gwh_open = .true.
     code = 0
     do i = 1, N
        if (label(i) >= 1 .and. label(i) <= Nlabel) then
           if (code(label(i)) == 0) code(label(i)) = labcode(i)
        end if
     end do
     write(gwh_unit,'(a)') '# TAPW projected Hamiltonian H(k) [eV] before diagonalisation; column (iG-1)*Nlabel+label'
     write(gwh_unit,'(a,i0,a,i0,a,i0)') '# M ', M, '   NG ', NG, '   Nlabel ', Nlabel
     write(gwh_unit,'(a,*(1x,i0))') '# labelcodes', (code(l), l = 1, Nlabel)
     do ig = 1, NG
        write(gwh_unit,'(a,i6,2(1x,es20.12))') '# G ', ig, Gx(ig), Gy(ig)
     end do
  end if
  write(gwh_unit,'(a,i6,2(1x,es20.12))') 'k ', ik, kvec(1), kvec(2)
  do j = 1, M
     do i = 1, M
        write(gwh_unit,'(2(1x,es23.15))') real(H(i,j), dp)*g0, aimag(H(i,j))*g0
     end do
  end do
  flush(gwh_unit)
end subroutine GWeightsHamWrite

subroutine build_X(X, xcoord, ycoord, label, Gx, Gy, N_orbit, N_G, N_label)
  implicit none
  integer, intent(in) :: N_orbit, N_G, N_label
  double precision, intent(in) :: xcoord(N_orbit), ycoord(N_orbit)          ! positions
  integer, intent(in) :: label(N_orbit)                 ! orbital label per atom (1..N_label)
  double precision, intent(in) :: Gx(N_G), Gy(N_G)                ! G vectors
  complex(dp), intent(out) :: X(N_orbit, N_G * N_label)

  integer :: i_orb, i_G, i_label, col
  integer :: count(N_label)
  double precision :: norm(N_label)
  double precision :: G_dot_r
  complex(dp) :: phase

  count = 0

  ! First count how many orbitals of each label
  do i_orb = 1, N_orbit
     count(label(i_orb)) = count(label(i_orb)) + 1
  end do

  ! Compute 1/sqrt(N_label) normalization
  do i_label = 1, N_label
     if (count(i_label) > 0) then
        norm(i_label) = 1.0d0 / sqrt(dble(count(i_label)))
     else
        norm(i_label) = 0.0d0
     end if
  end do

  if (tapwDebug) then
     ! Debug: Print label counts and normalizations
     print *, "DEBUG: Label counts and normalizations in build_X:"
     do i_label = 1, N_label
        print *, "  Label ", i_label, ": count =", count(i_label), ", norm = 1/√count =", norm(i_label)
     end do

     ! Debug: Check X-matrix phases for first few elements
     print *, "DEBUG: X-matrix phase factors (first 3 G-vectors, first 3 atoms):"
     do i_G = 1, min(3, N_G)
        do i_orb = 1, min(3, N_orbit)
           if (label(i_orb) == 1) then  ! Only check first label
              G_dot_r = Gx(i_G) * xcoord(i_orb) + Gy(i_G) * ycoord(i_orb)
              phase = cmplx(0.0_dp, G_dot_r, kind=dp)
              print *, "  G(", i_G, ") · r(", i_orb, ") =", G_dot_r
              print *, "  exp(i*G·r) =", exp(phase)
           end if
        end do
     end do
  end if

  ! Build X matrix with improved numerical precision
  do i_G = 1, N_G
     do i_label = 1, N_label
        col = (i_G - 1) * N_label + i_label
        do i_orb = 1, N_orbit
           if (label(i_orb) == i_label) then
              ! Calculate G·r in double precision to minimize rounding errors
              G_dot_r = Gx(i_G) * xcoord(i_orb) + Gy(i_G) * ycoord(i_orb)
              phase = cmplx(0.0_dp, G_dot_r, kind=dp)
              X(i_orb, col) = norm(i_label) * exp(phase)
           else
              X(i_orb, col) = cmplx(0.0_dp, 0.0_dp, kind=dp)
           end if
        end do
     end do
  end do
end subroutine build_X

!> Loewdin orthonormalisation of the TAPW basis, X <- X S^(-1/2) with S = X^H X.
!! build_X is orthonormal only for atoms on lattice sites; with relaxed positions S /= 1 (measured
!! eig(S) 0.988..1.007 on the relaxed 14.5 nm G/hBN cell) and the plain ZHEEV on X^H H X is then not
!! the Rayleigh-Ritz problem of the span of X.  After this call X^H X = 1 and the span is unchanged.
!! The spectrum of S is printed on the first call (the absence of that line = the flag did not act).
subroutine lowdin_orthonormalize_X(X, N, M)
  use constants, only : cmplx_0, cmplx_1
  implicit none
  integer, intent(in) :: N, M
  complex(dp), intent(inout) :: X(N, M)
  complex(dp), allocatable :: S(:,:), T(:,:), P(:,:), Y(:,:), work(:)
  real(dp), allocatable :: ev(:), rwork(:)
  integer :: i, info, lwork
  integer, save :: ncall = 0
  real(dp) :: dmax

  allocate(S(M, M), T(M, M), ev(M), rwork(max(1, 3*M - 2)))
  call ZGEMM('C', 'N', M, M, N, cmplx_1, X, N, X, N, cmplx_0, S, M)
  dmax = 0.0_dp
  do i = 1, M
     S(i, i) = S(i, i) - cmplx_1
     dmax = max(dmax, maxval(abs(S(:, i))))
     S(i, i) = S(i, i) + cmplx_1
  end do
  lwork = -1
  allocate(work(1))
  call ZHEEV('V', 'U', M, S, M, ev, work, lwork, rwork, info)
  lwork = max(1, int(real(work(1))))
  deallocate(work); allocate(work(lwork))
  call ZHEEV('V', 'U', M, S, M, ev, work, lwork, rwork, info)
  if (info /= 0 .or. minval(ev) <= 0.0_dp) then
     call MIO_Kill('Loewdin: X^H X is not positive definite (info '//trim(num2str(info))// &
                   ', min eig '//trim(num2str(minval(ev),6))//')', 'diag', 'lowdin_orthonormalize_X')
  end if
  ncall = ncall + 1
  if (ncall == 1) call MIO_Print('TAPW Loewdin: eig(X^H X) in ['//trim(num2str(minval(ev),8))//', '// &
       trim(num2str(maxval(ev),8))//'], max|X^H X - 1| = '//trim(num2str(dmax,8)), 'diag')
  ! P = S^(-1/2) = U diag(ev^-1/2) U^H  (U in S after ZHEEV); then X <- X P
  do i = 1, M
     T(:, i) = S(:, i) / sqrt(ev(i))
  end do
  allocate(P(M, M), Y(N, M))
  call ZGEMM('N', 'C', M, M, M, cmplx_1, T, M, S, M, cmplx_0, P, M)
  call ZGEMM('N', 'N', N, M, M, cmplx_1, X, N, P, M, cmplx_0, Y, N)
  X = Y
  deallocate(work, S, T, P, Y, ev, rwork)
end subroutine lowdin_orthonormalize_X

logical function is_structurally_symmetric(a, ia, ja, n)
    implicit none
    integer, intent(in) :: ia(:), ja(:)
    complex(dp), intent(in) :: a(:)
    integer, intent(in) :: n
    integer :: i, j, k, m
    logical :: found

    ! Declare arrays to store transpose of the CSR matrix
    integer, allocatable :: trans_ia(:), trans_ja(:)
    integer :: nnz

    ! Initialize the return value
    is_structurally_symmetric = .true.

    ! Number of non-zeros in the matrix
    nnz = size(ja)

    ! Allocate transpose arrays
    allocate(trans_ia(n + 1))
    allocate(trans_ja(nnz))

    ! Initialize trans_ia array
    trans_ia = 0

    ! Count the number of entries in each column of the original matrix
    do i = 1, n
        do k = ia(i), ia(i+1) - 1
            j = ja(k)
            trans_ia(j+1) = trans_ia(j+1) + 1
        end do
    end do

    ! Convert counts to starting indices
    trans_ia(1) = 1
    do i = 2, n + 1
        trans_ia(i) = trans_ia(i) + trans_ia(i-1)
    end do

    ! Fill the trans_ja array
    do i = 1, n
        do k = ia(i), ia(i+1) - 1
            j = ja(k)
            m = trans_ia(j)
            trans_ja(m) = i
            trans_ia(j) = m + 1
        end do
    end do

    ! Restore the trans_ia array to correct starting indices
    do i = n, 1, -1
        trans_ia(i+1) = trans_ia(i)
    end do
    trans_ia(1) = 1

    ! Compare the structure of the original matrix with its transpose
    do i = 1, n
        do k = ia(i), ia(i+1) - 1
            j = ja(k)
            found = .false.
            do m = trans_ia(j), trans_ia(j+1) - 1
                if (trans_ja(m) == i) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                is_structurally_symmetric = .false.
                return
            end if
        end do
    end do

    ! Clean up
    deallocate(trans_ia)
    deallocate(trans_ja)
end function is_structurally_symmetric

real(dp) function norm2(x)
    implicit none
    complex(dp), intent(in) :: x(:)
    norm2 = sqrt(sum(abs(x)**2))
end function norm2

subroutine initialize_sparse_matrix(N, maxN, H0, hopp, NList, Nneigh, neighCell, ns, is, KLoc, cell, row_ptr, col_ind, values,sigma)
    use constants, only : cmplx_i
    use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
    use scf, only : charge, Zch
    use atoms, only : Species, layerIndex
    use tbpar, only : U
    use neigh, only : NeighD
    use ham, only : Zterm, gZeeman, IntrinsicSOCterm, lambdaI, IsingSOCterm, lambdaIsing, PIASOCterm, lambdaPIA, SOCEnabledForLayer
    use magf, only : BmagZeeman
    implicit none
    integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
    real(dp), intent(in) :: KLoc(3), cell(3,3)
    real(dp), intent(inout) :: H0(N)
    complex(dp), intent(in) :: hopp(maxN,N)
    integer, allocatable, intent(out) :: row_ptr(:), col_ind(:)
    complex(dp), allocatable, intent(out) :: values(:)
    integer :: i, j, k, l, nnz, in
    real(dp) :: R(3)
    integer :: nnz_temp
    double precision :: alpha
    logical :: useShift, found
    complex(dp) :: sigma
    real(dp) :: soc_diagonal_contrib

    call MIO_InputParameter('Diag.SparseUseShift',useShift,.false.)

    ! First pass to count non-zero elements
    nnz_temp = 0
    do i = 1, N
        if (useShift) then
           nnz_temp = nnz_temp + 1
        else
           if (H0(i) /= 0.0_dp) nnz_temp = nnz_temp + 1
        end if
        do j = 1, Nneigh(i)
            if (hopp(j,i) /= 0.0_dp) nnz_temp = nnz_temp + 1
        end do
    end do
    if (edgeHopp) then
        do i = 1, nQ
            do j = 1, nEdgeN(i)
                if (edgeH(j,i) /= 0.0_dp) nnz_temp = nnz_temp + 1
            end do
        end do
    end if

    ! Allocate space for CSR arrays
    allocate(row_ptr(N+1), col_ind(nnz_temp + 1000), values(nnz_temp + 1000))

    row_ptr(1) = 1
    l = 1  ! Use a different variable to increment the non-zero element index

    ! Populate CSR arrays directly
    do i = 1, N
        row_ptr(i) = l
        ! Diagonal element - calculate SOC contributions first
        ! Note: SOC terms may be non-zero even if H0(i) == 0, so we need to check
        soc_diagonal_contrib = 0.0_dp

        ! Calculate SOC contributions to diagonal element
        ! Note: This is backward compatible - only applies when SOC is enabled
        if (Zterm .or. IntrinsicSOCterm .or. IsingSOCterm .or. PIASOCterm) then
           ! Apply Zeeman effect (λVZ) - spin-dependent onsite
           if (Zterm .and. (ns==2)) then
              if (is==1) then
                 soc_diagonal_contrib = soc_diagonal_contrib + gZeeman*BmagZeeman  ! Spin-up: +λVZ
              else
                 soc_diagonal_contrib = soc_diagonal_contrib - gZeeman*BmagZeeman  ! Spin-down: -λVZ
              end if
           end if

           ! Apply Intrinsic SOC (λI) - gap-opening onsite term
           ! For graphene, this should be sublattice-dependent to open gap
           if (IntrinsicSOCterm .and. SOCEnabledForLayer(layerIndex(i))) then
              ! Apply +λI to sublattice A, -λI to sublattice B (or vice versa)
              ! This opens a gap of size 2λI at Dirac point
              ! Species(i)=1 is A sublattice, Species(i)=2 is B sublattice
              if (Species(i) == 1) then
                 soc_diagonal_contrib = soc_diagonal_contrib + lambdaI
              else
                 soc_diagonal_contrib = soc_diagonal_contrib - lambdaI
              end if
           end if

           ! Apply Ising SOC (λIsing) - Valley-Zeeman term (τ_z s_z)
           ! This is a valley-dependent spin splitting: τ_z = +1 for K valley, -1 for K' valley
           ! Does not open a global gap, only splits spins differently in K vs K' valleys
           if (IsingSOCterm .and. SOCEnabledForLayer(layerIndex(i)) .and. (ns==2)) then
              ! Valley sign: +1 for K valley, -1 for K' valley
              if (.not. useKprimeValley) then
                 ! K valley: τ_z = +1
                 if (is == 1) then
                    ! Spin-up: +λ for all sites
                    soc_diagonal_contrib = soc_diagonal_contrib + lambdaIsing
                 else
                    ! Spin-down: -λ for all sites
                    soc_diagonal_contrib = soc_diagonal_contrib - lambdaIsing
                 end if
              else
                 ! K' valley: τ_z = -1 (flip sign)
                 if (is == 1) then
                    ! Spin-up: -λ for all sites
                    soc_diagonal_contrib = soc_diagonal_contrib - lambdaIsing
                 else
                    ! Spin-down: +λ for all sites
                    soc_diagonal_contrib = soc_diagonal_contrib + lambdaIsing
                 end if
              end if
           end if

           ! Apply Pseudo-inversion asymmetry (λPIA) - onsite component
           if (PIASOCterm .and. SOCEnabledForLayer(layerIndex(i))) then
              ! PIA is spin-independent onsite term
              soc_diagonal_contrib = soc_diagonal_contrib + lambdaPIA
           end if
        end if

        ! Add SCF terms if spin-polarized (matches DiagHam implementation)
        ! Commented out: not doing any SCF calculation for now (matches BuildBlockHamiltonianOnly)

        ! Add diagonal element if H0(i) is non-zero OR if SOC/SCF contributions are non-zero
        if (useShift .or. H0(i) /= 0.0_dp .or. soc_diagonal_contrib /= 0.0_dp) then
            values(l) = H0(i) - sigma + soc_diagonal_contrib
            col_ind(l) = i
            l = l + 1
        end if

        ! Off-diagonal elements
        !             !found = .false.
        !             !! Check if we already have an entry for this off-diagonal element
        !             !do k = row_ptr(i), l - 1
        !             !    if (col_ind(k) == in) then
        !             !        print*, "entering this loop means you are working on a small system, no?"
        !             !        values(k) = values(k) - hopp(j, i) * exp(cmplx_i * dot_product(KLoc, R))
        !             !        found = .true.
        !             !        exit
        !             !    end if
        !             !end do
        !             !! If no entry exists, add a new one
        !             !if (.not. found) then
        !             !end if
           do j = 1, Nneigh(i)
               in = NList(j,i)
               ! Use actual atomic position difference instead of lattice vector
               ! NeighD(:,j,i) contains the vector from atom i to atom j: (τ_j - τ_i) + T_ij
               ! For consistency with paper formulation, use -NeighD to get distance from i to j
               R(1:3) = 0.0_dp
               R(1:2) = -NeighD(1:2, j, i)  ! Use negative to get distance from i to j
               if (hopp(j,i) /= 0.0_dp) then
                  found = .false.
                  ! Check if we already have an entry for this off-diagonal element
                  do k = row_ptr(i), l - 1
                      if (col_ind(k) == in) then
                          values(k) = values(k) - conjg(hopp(j, i)) * exp(cmplx_i * dot_product(KLoc, R))
                          found = .true.
                          exit
                      end if
                  end do
                  ! If no entry exists, add a new one
                  if (.not. found) then
                      if (l > nnz_temp) then
                         print *, "CSR OVERFLOW: l =", l, " > nnz_temp =", nnz_temp
                         error stop "Sparse matrix allocation overflow in initialize_sparse_matrix"
                      endif
                      ! Debug first few k-dependent phases for K-point analysis (now using NeighD)
                      if (tapwDebug .and. i <= 3 .and. j <= 2) then
                         print *, "DEBUG: atom", i, "neighbor", j, ": R (from NeighD) =", R
                         print *, "DEBUG: k·R =", dot_product(KLoc,R)
                         print *, "DEBUG: exp(-i*k·R) =", exp(-cmplx_i*dot_product(KLoc,R))
                      end if
                      ! CSR row i, column in holds element (i,in). The dense TAPW path writes -hopp(j,i)*exp(-ik.R) at
                      ! (in,i), so (i,in) is its Hermitian conjugate. Writing the unconjugated value here stored H^T = H*,
                      ! which swaps the K and K' valleys in X^H H X (measured 2026-09-30: sparse K = dense K' to 0.001 meV).
                      ! Full-TB sparse eigenvalues are unaffected (H and H* have the same spectrum).
                      values(l) = -conjg(hopp(j, i)) * exp(cmplx_i * dot_product(KLoc, R))
                      col_ind(l) = in
                      l = l + 1
                  end if
                  ! Apply PIA hopping terms if enabled (matches DiagHam implementation)
                  if (PIASOCterm) then
                     ! PIA hopping adds imaginary off-diagonal terms
                     ! Note: This modifies the value we just added, so we need to get the correct index
                     ! Find the index in values array for this off-diagonal element
                     do k = row_ptr(i), l - 1
                        if (col_ind(k) == in) then
                           ! Apply PIA hopping: lambdaPIA * (neighD_y + i*neighD_x)
                           ! Note: NeighD is from i to j, matching ApplyPIAHopping convention
                           values(k) = values(k) + lambdaPIA * cmplx(NeighD(2, j, i), NeighD(1, j, i))
                           exit
                        end if
                     end do
                  end if
               end if
           end do

        ! Edge hopping elements
    end do

    nnz = l - 1
    row_ptr(N+1) = nnz + 1

    ! Trim the allocated arrays to actual size

    ! Debug print for sparse matrix
    print *, 'Sparse matrix row_ptr: ', row_ptr(1:min(N+1,10))
    print *, 'Sparse matrix col_ind: ', col_ind(1:min(nnz_temp,10))
    print *, 'Sparse matrix values: ', values(1:min(nnz,10))
    print *, "CSR max l =", l-1, "allocated nnz_temp =", nnz_temp
end subroutine initialize_sparse_matrix

!subroutine initialize_sparse_matrix(N, maxN, H0, hopp, NList, Nneigh, neighCell, ns, is, KLoc, cell, row_ptr, col_ind, values)
!
!          !zz = charge(1,i)*charge(2,i) ! Zch
!          !print*, "nownow", KLoc
!          !print*, "niwniw", R
!    ! Count the number of non-zero elements
!
!    ! Allocate space for CSR arrays
!
!
!    ! Populate CSR arrays
!    ! Parallelize the outer loop with OpenMP
!    !!$OMP PARALLEL DO PRIVATE(i, j, k_local) SHARED(row_ptr, col_ind, values, HLoc) REDUCTION(+:k)
!    !do i = 1, N
!    !    k_local = k  ! Local copy of k for each thread
!    !    do j = 1, N
!    !        if (HLoc(i, j) /= (0.0, 0.0)) then
!    !            values(k_local) = HLoc(j, i)
!    !            col_ind(k_local) = j
!    !            k_local = k_local + 1
!    !        end if
!    !    end do
!    !    row_ptr(i+1) = k_local
!    !    k = k_local  ! Update global k with the local k value
!    !end do
!    !!$OMP END PARALLEL DO
!
!    ! Debug print for sparse matrix
!    !print *, 'Sparse matrix row_ptr: ', row_ptr(1:min(N+1,10))
!    !print *, 'Sparse matrix col_ind: ', col_ind(1:min(nnz,10))
!    !print *, 'Sparse matrix values: ', values(1:min(nnz,10))

!subroutine sparse_matvec(nn, row_ptr, col_ind, values, x, y)
!

subroutine sparse_matvec(nn, row_ptr, col_ind, values, x, y)
    implicit none
    integer, intent(in) :: nn, row_ptr(:), col_ind(:)
    complex(dp), intent(in) :: values(:), x(:)
    complex(dp), intent(out) :: y(nn)
    integer :: i, j

    y = 0.0_dp

    ! Parallelize the outer loop with OpenMP
    !$OMP PARALLEL DO PRIVATE(i, j) SHARED(row_ptr, col_ind, values, x, y)
    do i = 1, nn
        do j = row_ptr(i), row_ptr(i+1) - 1
            y(i) = y(i) + values(j) * x(col_ind(j))
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine sparse_matvec

subroutine DiagHamWF(N,ns,is,HLoc,ELoc,KLoc,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species
   use tbpar,                 only : U

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp), intent(out) :: HLoc(N,N)
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   real(dp) :: R(3), zz

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      ! Apply SOC modifications
      if (anySOCEnabled) call ApplySOCtoHamiltonian(i, is, ns, HLoc)

      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(-cmplx_i*dot_product(KLoc,R))
         ! PIA hopping not yet properly implemented - commented out
      end do
   end do
   if (edgeHopp) then
      do i=1,nQ
         do j=1,nEdgeN(i)
            in = NeI(j,i)
            R = matmul(cell,NedgeCell(:,j,i))
            HLoc(in,edgeIndx(i)) = HLoc(in,edgeIndx(i)) + edgeH(j,i)*exp(-cmplx_i*dot_product(KLoc,R))
         end do
      end do
   end if
   call ZHEEV('V','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   if (info/=0) then
      call MIO_Kill('Error in diagonalization','diag','DiagHam')
   end if

end subroutine DiagHamWF

subroutine DiagHamChern(N,ns,is,HLoc,ChernLoc,KLoc,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species
   use tbpar,                 only : U

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp), intent(out) :: HLoc(N,N)
   complex(dp) :: HLocM1dx(N,N), HLocM1dy(N,N), H_derivativeDX(N,N), H_derivativeDY(N,N)
   complex(dp) :: dHdx(N,N)
   complex(dp) :: dHdy(N,N)
   complex(dp) :: eigvec(N,N), Vxmn(N,N), Vymn(N,N), temp_matrix(N,N), ones_matrix(N,N), difference_matrix(N,N), eigval_matrix(N,N), eigval_repeated_matrix(N, N)
   real(dp), intent(out) :: ChernLoc(N)
   complex(dp) :: ChernLocSum(N)
   real(dp) :: ELoc(N)
   real(dp) :: eigval(N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   real(dp) :: KLocM1(3)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info, k
   real(dp) :: R(3), zz, eps

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   real(8) :: delta_kx, delta_ky

   delta_kx = 0.01d0
   delta_ky = 0.01d0

   call MIO_InputParameter('Epsilon',eps,0.001_dp)

   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      ! Apply SOC modifications
      call ApplySOCtoHamiltonian(i, is, ns, HLoc)

      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         KLocM1(1) = KLoc(1)-delta_kx
         KLocM1(2) = KLoc(2)
         KLocM1(3) = KLoc(3)
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(-cmplx_i*dot_product(KLocM1,R))
      end do
   end do
   HLocM1dy = HLoc

   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      ! Apply SOC modifications
      call ApplySOCtoHamiltonian(i, is, ns, HLoc)

      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         KLocM1(1) = KLoc(1)
         KLocM1(2) = KLoc(2)-delta_ky
         KLocM1(3) = KLoc(3)
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(-cmplx_i*dot_product(KLocM1,R))
      end do
   end do
   HLocM1dy = HLoc

   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      ! Apply SOC modifications
      if (anySOCEnabled) call ApplySOCtoHamiltonian(i, is, ns, HLoc)

      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(-cmplx_i*dot_product(KLoc,R))
         ! PIA hopping not yet properly implemented - commented out
      end do
   end do

   ! Assuming HLoc and HLocM1 are matrices, you would first initialize a matrix for the derivative
   H_derivativeDX = 0.0_dp
   H_derivativeDY = 0.0_dp

   ! Calculation of Hamiltonians remains the same, so I won't repeat that here.

   ! Now, compute the finite difference derivative
   do i=1,N
      do j=1,N
         H_derivativeDX(i,j) = (HLoc(i,j) - HLocM1dx(i,j)) / (2.0_dp * delta_kx)
         dHdx(i,j) = cmplx(0.0, 1.0) * HLoc(i,j) * H_derivativeDX(i,j)
         H_derivativeDY(i,j) = (HLoc(i,j) - HLocM1dy(i,j)) / (2.0_dp * delta_ky)
         dHdy(i,j) = cmplx(0.0, 1.0) * HLoc(i,j) * H_derivativeDY(i,j)
      end do
   end do

   call ZHEEV('V','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   if (info/=0) then
      call MIO_Kill('Error in diagonalization','diag','DiagHam')
   end if

   ! Calculating Vxmn and Vymn
   eigvec = HLoc
   eigval = ELoc
   print*, "eigval", eigval

   Vxmn = matmul(matmul(transpose(eigvec), dHdx), eigvec)
   Vymn = matmul(matmul(transpose(eigvec), dHdy), eigvec)

   ! Setting diagonals to zero
   do j = 1, N
       Vxmn(j,j) = 0.0_dp
       Vymn(j,j) = 0.0_dp
   end do
   print*, "Vxmn", Vxmn

   !! Calculating Chern
   forall(i=1:N, j=1:N) ones_matrix(i, j) = 1.0_dp

   ! Calculating Chern
   forall(i=1:N, j=1:N) eigval_matrix(i, j) = eigval(j)
   forall(i=1:N, j=1:N) eigval_repeated_matrix(i, j) = eigval(i)
   difference_matrix = eigval_matrix - eigval_repeated_matrix
   temp_matrix = Vxmn * transpose(Vymn) / (difference_matrix**2 + eps)

   do j = 1, N
       do k = 1, N
           ChernLoc = ChernLoc + 2.0 * aimag(temp_matrix(j, k))
       end do
   end do
   ! Chern[i] = 2*np.imag(np.sum(Vxmn*(Vymn.T/((eigval*np.ones((N_orbit,N_orbit))).T-eigval)**2+eps), axis=1))

end subroutine DiagHamChern

subroutine DiagHamPDOS(N,ns,is,HLoc,ELoc,KLoc,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species
   use tbpar,                 only : U

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp), intent(out) :: HLoc(N,N)
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   real(dp) :: R(3), zz

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      if (ns==2) then
         if (is==1) then
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
         else
            HLoc(i,i) = HLoc(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
         end if
      end if
      ! Apply SOC modifications
      if (anySOCEnabled) call ApplySOCtoHamiltonian(i, is, ns, HLoc)

      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(-cmplx_i*dot_product(KLoc,R))
         ! PIA hopping not yet properly implemented - commented out
      end do
   end do
   if (edgeHopp) then
      do i=1,nQ
         do j=1,nEdgeN(i)
            in = NeI(j,i)
            R = matmul(cell,NedgeCell(:,j,i))
            HLoc(in,edgeIndx(i)) = HLoc(in,edgeIndx(i)) + edgeH(j,i)*exp(-cmplx_i*dot_product(KLoc,R))
         end do
      end do
   end if
   call ZHEEV('V','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   if (info/=0) then
      call MIO_Kill('Error in diagonalization','diag','DiagHamPDOS')
   end if

end subroutine DiagHamPDOS

subroutine DiagHamArpack(N,ns,is,HLoc,ELoc,KLoc,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species
   use tbpar,                 only : U

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp), intent(out) :: HLoc(N,N)
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   real(dp) :: R(3), zz

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   HLoc = 0.0_dp
   do i=1,N
      HLoc(i,i) = H0(i)
      ! Add SCF terms if spin-polarized (matches other routines)
      ! Commented out: not doing any SCF calculation for now (matches BuildBlockHamiltonianOnly)
      !   !zz = charge(1,i)*charge(2,i) ! Zch
      ! Apply SOC modifications
      call ApplySOCtoHamiltonian(i, is, ns, HLoc)

      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         HLoc(in,i) = HLoc(in,i) - hopp(j,i)*exp(-cmplx_i*dot_product(KLoc,R))
      end do
   end do
   if (edgeHopp) then
      do i=1,nQ
         do j=1,nEdgeN(i)
            in = NeI(j,i)
            R = matmul(cell,NedgeCell(:,j,i))
            HLoc(in,edgeIndx(i)) = HLoc(in,edgeIndx(i)) + edgeH(j,i)*exp(-cmplx_i*dot_product(KLoc,R))
         end do
      end do
   end if
   call ZHEEV('N','L',N,HLoc,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   if (info/=0) then
      call MIO_Kill('Error in diagonalization','diag','DiagHamArpack')
   end if

end subroutine DiagHamArpack

subroutine DiagSpectralWeightNishi(N,ns,is,Pkc,E,K,KG,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   complex(dp), intent(out) :: Pkc(N)
   complex(dp) :: Pkcaux
   real(dp), intent(out) :: E(N)
   real(dp), intent(in) :: K(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy

   real(dp) :: Nc, nc2

   integer :: nTS
   integer :: kk

   call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
   gcell(:,1) = [aG,0.0_dp,0.0_dp]
   gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
   gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]

   cellts = gcell
   celltsp = gcell

   Hts = 0.0_dp
   Htsp = 0.0_dp
   Pkc = 0.0_dp

   RtsVec(1,:) = [matmul(cellts,[1,0,0])]
   RtsVec(2,:) = [matmul(cellts,[0,1,0])]
   RtsVec(3,:) = [matmul(cellts,[1,1,0])]
   RtsVec(4,:) = [matmul(cellts,[-1,0,0])]
   RtsVec(5,:) = [matmul(cellts,[0,-1,0])]
   RtsVec(6,:) = [matmul(cellts,[-1,-1,0])]

   RtspVec(1,:) = [matmul(cellts,[1,0,0])]
   RtspVec(2,:) = [matmul(cellts,[0,1,0])]
   RtspVec(3,:) = [matmul(cellts,[1,1,0])]
   RtspVec(4,:) = [matmul(cellts,[-1,0,0])]
   RtspVec(5,:) = [matmul(cellts,[0,-1,0])]
   RtspVec(6,:) = [matmul(cellts,[-1,-1,0])]

   call MIO_InputParameter('Spectral.numberOfTS',nTS,1)

   nTS = 6

   do i1=1,nTS!SIZE(RtsVec)
      do i2=1,nTS!SIZE(RtspVec)
         Rts = RtsVec(i1,:)
         Rtsp = RtspVec(i2,:)
         do i=1,N
            Hts(i,i) = H0(i)
            Htsp(i,i) = H0(i)
            do j=1,Nneigh(i)
               in = NList(j,i)
               Tm = matmul(cell,neighCell(:,j,i))
               Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(K,Tm-Rts))
               Htsp(in,i) = Htsp(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(K,Tm-Rtsp))
            end do
         end do
         call ZHEEV('V','L',N,Hts,N,E,ZWork,lwork,DWork,info)
         call ZHEEV('V','L',N,Htsp,N,E,ZWork,lwork,DWork,info)

         do i=1,N
            do in=1,N
               Pkcaux = CONJG(Hts(in,i)) * Htsp(in,i)
               Pkc(i) = Pkc(i) + exp(cmplx_i*dot_product(KG,Rts-Rtsp)) * Pkcaux
            end do
         end do
      end do
   end do
   Nc = 1.0_dp
   nc2 = 5.0_dp*5.0_dp
   Pkc = Nc/nc2 * Pkc

   if (info/=0) then
      call MIO_Kill('Error in diagonalization','diag','DiagSpectralWeight')
   end if

   ! steps

end subroutine DiagSpectralWeightNishi

subroutine DiagSpectralWeightWeiKu(N,ns,is,PkcLoc,E,K,KG,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   complex(dp), intent(out) :: PkcLoc(N,2)
   complex(dp) :: Pkcaux, Pkcaux1, Pkcaux3
   real(dp), intent(out) :: E(N)
   real(dp), intent(in) :: K(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy

   real(dp) :: Nc, nc2

   integer :: nTS
   integer :: kk

   integer :: cellSize

   integer :: at1, at2

   call MIO_InputParameter('LatticeParameter',aG,2.46_dp)
   gcell(:,1) = [aG,0.0_dp,0.0_dp]
   gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
   gcell(:,3) = [0.0_dp,0.0_dp,40.0_dp]

   cellts = gcell
   celltsp = gcell

   Hts = 0.0_dp
   PkcLoc = 0.0_dp

   call MIO_InputParameter('CellSize', cellSize, 1)

   if (.not. frac) call AtomsSetFrac()
   do i=1,N
       RtsVec(i,:) = matmul(cellts,[floor((Rat(1,i)-0.001)*cellSize),floor(Rat(2,i)*cellSize),0]) ! Moire
   end do

         do i=1,N
            Hts(i,i) = H0(i)
            do j=1,Nneigh(i)
               in = NList(j,i)
               Tm = matmul(cell,neighCell(:,j,i))
               Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(K,Tm))
            end do
         end do
         call ZHEEV('V','L',N,Hts,N,E,ZWork,lwork,DWork,info)

         !   !do i2=1,N
         !   !do j=1,Nneigh(i) ! maybe go back to full system
         !      !in = i
         !      !in = NList(j,i)
         !      !Rts = matmul(cellts,[1,0,0]) ! Moire
         !      !Rts = matmul(cellts,neighCell(:,j,i)) ! Moire
         !      !Rtsp = matmul(celltsp,neighCell(:,j,i)) ! graphene
         !      !Rts = matmul(cellts,[floor(Rat(1,in)/RtsVec(1,1)),floor(Rat(2,in)/RtsVec(1,2)),0]) ! Moire
         !      !Rtsp = matmul(celltsp,[floor(Rat(1,i)/RtsVec(1,1)),floor(Rat(2,i)/RtsVec(1,2)),0]) ! Moire
         !      !Pkcaux = CONJG(Hts(i,in)) * Htsp(i,in)
         !      !Pkcaux = Pkcaux + exp(-cmplx_i*dot_product(KG,RtsVec(in,:))) * Hts(in,i)
         !   !end do
         !   !do i2=1,N
         !   !do j=1,Nneigh(i) ! maybe go back to full system
         !      !in = i
         !      !in = NList(j,i)
         !      !Rts = matmul(cellts,[1,0,0]) ! Moire
         !      !Rts = matmul(cellts,neighCell(:,j,i)) ! Moire
         !      !Rtsp = matmul(celltsp,neighCell(:,j,i)) ! graphene
         !      !Rts = matmul(cellts,[floor(Rat(1,in)/RtsVec(1,1)),floor(Rat(2,in)/RtsVec(1,2)),0]) ! Moire
         !      !Rtsp = matmul(celltsp,[floor(Rat(1,i)/RtsVec(1,1)),floor(Rat(2,i)/RtsVec(1,2)),0]) ! Moire
         !      !Pkcaux = CONJG(Hts(i,in)) * Htsp(i,in)
         !      !if (i.eq.in) then
         !      !end if
         !      !Pkcaux1 = Pkcaux1 + exp(-cmplx_i*dot_product(KG,RtsVec(in,:))) * Hts(in,1)
         !      !Pkcaux3 = Pkcaux3 + exp(-cmplx_i*dot_product(KG,RtsVec(in,:))) * Hts(in,3)
         !   !PkcLoc(1) = PkcLoc(1) + Pkcaux1
         !   !PkcLoc(2) = PkcLoc(2) + Pkcaux3
         do j=1,N ! These are the eigenvectors with band index J
             do in=1,N ! NOT the sum over eigenvectors. Pick one eigenvector and then sum over its coefficients. Each coefficient corresponds to one orbital in Wannier (or TB) basis.
                 if (Species(in).eq.1) then
                     PkcLoc(j,1) = PkcLoc(j,1) + exp(-cmplx_i*dot_product(KG, RtsVec(in,:))) * Hts(in,j) ! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
                 else if (Species(in).eq.2) then
                     PkcLoc(j,2) = PkcLoc(j,2) + exp(-cmplx_i*dot_product(KG, RtsVec(in,:))) * Hts(in,j)
                 end if
             end do
         end do
   Nc = 1.0
   nc2 = 1.0
   PkcLoc = Nc/nc2 * PkcLoc

   if (info/=0) then
      call MIO_Kill('Error in diagonalization','diag','DiagSpectralWeight')
   end if

   ! steps

end subroutine DiagSpectralWeightWeiKu

subroutine DiagSpectralWeightWeiKuInequivalentOld(N,ns,is,PkcLocA,PkcLocB,ELoc,KptsLoc,KG,cell,gcell,H0,maxN,hopp,NList,Nneigh,neighCell,topBottomRatio)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac, layerIndex
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   complex(dp), intent(out) :: PkcLocA(N,3), PkcLocB(N,3)
   complex(dp) :: Pkcaux, Pkcaux1, Pkcaux3
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KptsLoc(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy

   real(dp) :: Nc, nc2

   integer :: nTS
   real(dp) :: ll, kk

   integer :: cellSize

   integer :: at1, at2

   character(len=80) :: line
   integer :: id

   real(dp) :: G01(3), G10(3), G11(3)
   real(dp) :: G01p(3), G10p(3), G11p(3)

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   real(dp) :: topBottomRatio

   logical :: changeExpSign

   cellts = gcell ! Lattice vectors of PC (either top or bottom layer)

   Hts = 0.0_dp
   PkcLocA = 0.0_dp
   PkcLocB = 0.0_dp
   if (frac) call AtomsSetCart()
       G10 = matmul(cellts,[1,0,0]) ! here, its in real space!!!
       G01 = matmul(cellts,[0,1,0])
       G11 = matmul(cellts,[1,1,0])
       do i=1,N
           ll = (Rat(1,i)*G10(2)/G10(1) - Rat(2,i)) / (G01(1)*G10(2)/G10(1) - G01(2))
           kk = (Rat(1,i) - ll * G01(1)) / G10(1)
           if (kk.gt.0) then
               kk = floor(kk)
           else
               kk = ceiling(kk)
           end if
           if (ll.gt.0) then
               ll = floor(ll)
           else
               ll = ceiling(ll)
           end if
           RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire
       end do
         do i=1,N
            Hts(i,i) = H0(i)
            do j=1,Nneigh(i)
               in = NList(j,i)
               Tm = matmul(cell,neighCell(:,j,i))
               Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KptsLoc,Tm))
            end do
         end do
         call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
         if (info/=0) then
            call MIO_Kill('Error in diagonalization for spectral function','diag','DiagSpectralWeight')
         end if
         ! Note that I'm using the Wei Ku expression for the spectral function
         ! (this is NOT the NISHI one)
             do j=1,N
                 do in=1,N
                   if (Species(in).eq.1) then
                     if (layerIndex(in).eq.1) then
                         PkcLocA(j,1) = PkcLocA(j,1) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.2) then
                         PkcLocA(j,2) = PkcLocA(j,2) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.3) then
                         PkcLocA(j,3) = PkcLocA(j,3) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
                     end if
                   else if (Species(in).eq.2) then
                     if (layerIndex(in).eq.1) then
                         PkcLocB(j,1) = PkcLocB(j,1) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.2) then
                         PkcLocB(j,2) = PkcLocB(j,2) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.3) then
                         PkcLocB(j,3) = PkcLocB(j,3) - exp(cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j)
                     end if
                   end if
                 end do
             end do
   Nc = 1.0
   nc2 = 1.0
   PkcLocA = Nc/nc2 * PkcLocA
   PkcLocB = Nc/nc2 * PkcLocB

end subroutine DiagSpectralWeightWeiKuInequivalentOld

subroutine DiagSpectralWeightWeiKuInequivalent(N,ns,is,PkcLocA,PkcLocB,ELoc,KptsLoc,KG,cell,gcell,H0,maxN,hopp,NList,Nneigh,neighCell,topBottomRatio)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac, layerIndex
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   complex(dp), intent(out) :: PkcLocA(N,4), PkcLocB(N,4)
   complex(dp) :: Pkcaux, Pkcaux1, Pkcaux3
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KptsLoc(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy

   real(dp) :: Nc, nc2

   integer :: nTS
   real(dp) :: ll, kk

   integer :: cellSize

   integer :: at1, at2

   character(len=80) :: line
   integer :: id

   real(dp) :: G01(3), G10(3), G11(3)
   real(dp) :: G01p(3), G10p(3), G11p(3)

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   real(dp) :: topBottomRatio

   logical :: changeExpSign

   cellts = gcell ! Lattice vectors of PC (either top or bottom layer)

   Hts = 0.0_dp
   PkcLocA = 0.0_dp
   PkcLocB = 0.0_dp
   if (frac) call AtomsSetCart()
       G10 = matmul(cellts,[1,0,0]) ! here, its in real space!!!
       G01 = matmul(cellts,[0,1,0])
       G11 = matmul(cellts,[1,1,0])
       do i=1,N
           ll = (Rat(1,i)*G10(2)/G10(1) - Rat(2,i)) / (G01(1)*G10(2)/G10(1) - G01(2))
           kk = (Rat(1,i) - ll * G01(1)) / G10(1)
           if (kk.gt.0) then
               kk = floor(kk)
           else
               kk = ceiling(kk)
           end if
           if (ll.gt.0) then
               ll = floor(ll)
           else
               ll = ceiling(ll)
           end if
           RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire
       end do
         do i=1,N
            Hts(i,i) = H0(i)
            do j=1,Nneigh(i)
               in = NList(j,i)
               Tm = matmul(cell,neighCell(:,j,i))
               Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KptsLoc,Tm))
            end do
         end do
         call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
         if (info/=0) then
            call MIO_Kill('Error in diagonalization for spectral function','diag','DiagSpectralWeight')
         end if
         ! Note that I'm using the Wei Ku expression for the spectral function
         ! (this is NOT the NISHI one)
             do j=1,N
                 do in=1,N
                   if (Species(in).eq.1) then
                     if (layerIndex(in).eq.1) then
                         PkcLocA(j,1) = PkcLocA(j,1) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.2) then
                         PkcLocA(j,2) = PkcLocA(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.3) then
                         PkcLocA(j,3) = PkcLocA(j,3) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.4) then
                         PkcLocA(j,4) = PkcLocA(j,4) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     end if
                   else if (Species(in).eq.2) then
                     if (layerIndex(in).eq.1) then
                         PkcLocB(j,1) = PkcLocB(j,1) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.2) then
                         PkcLocB(j,2) = PkcLocB(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.3) then
                         PkcLocB(j,3) = PkcLocB(j,3) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.4) then
                         PkcLocB(j,4) = PkcLocB(j,4) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     end if
                   end if
                 end do
             end do
   Nc = 1.0
   nc2 = 1.0
   PkcLocA = Nc/nc2 * PkcLocA
   PkcLocB = Nc/nc2 * PkcLocB

end subroutine DiagSpectralWeightWeiKuInequivalent

subroutine DiagSpectralWeightWeiKuInequivalentMoreOrbitals(N,ns,is,PkcLocA,PkcLocB,PkcLocC,PkcLocD,PkcLocE,PkcLocF,PkcLocG,PkcLocH,ELoc,KptsLoc,KG,cell,gcell,H0,maxN,hopp,NList,Nneigh,neighCell,topBottomRatio)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac, layerIndex
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   complex(dp), intent(out) :: PkcLocA(N,4), PkcLocB(N,4)
   complex(dp), intent(out) :: PkcLocC(N,4), PkcLocD(N,4)
   complex(dp), intent(out) :: PkcLocE(N,4), PkcLocF(N,4)
   complex(dp), intent(out) :: PkcLocG(N,4), PkcLocH(N,4)
   complex(dp) :: Pkcaux, Pkcaux1, Pkcaux3
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KptsLoc(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy

   real(dp) :: Nc, nc2

   integer :: nTS
   real(dp) :: ll, kk

   integer :: cellSize

   integer :: at1, at2

   character(len=80) :: line
   integer :: id

   real(dp) :: G01(3), G10(3), G11(3)
   real(dp) :: G01p(3), G10p(3), G11p(3)

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   real(dp) :: topBottomRatio

   logical :: changeExpSign

   cellts = gcell ! Lattice vectors of PC (either top or bottom layer)

   Hts = 0.0_dp
   PkcLocA = 0.0_dp
   PkcLocB = 0.0_dp
   PkcLocC = 0.0_dp
   PkcLocD = 0.0_dp
   PkcLocE = 0.0_dp
   PkcLocF = 0.0_dp
   PkcLocG = 0.0_dp
   PkcLocH = 0.0_dp
   if (frac) call AtomsSetCart()
       G10 = matmul(cellts,[1,0,0]) ! here, its in real space!!!
       G01 = matmul(cellts,[0,1,0])
       G11 = matmul(cellts,[1,1,0])
       do i=1,N
           ll = (Rat(1,i)*G10(2)/G10(1) - Rat(2,i)) / (G01(1)*G10(2)/G10(1) - G01(2))
           kk = (Rat(1,i) - ll * G01(1)) / G10(1)
           if (kk.gt.0) then
               kk = floor(kk)
           else
               kk = ceiling(kk)
           end if
           if (ll.gt.0) then
               ll = floor(ll)
           else
               ll = ceiling(ll)
           end if
           RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire
       end do
         do i=1,N
            Hts(i,i) = H0(i)
            do j=1,Nneigh(i)
               in = NList(j,i)
               Tm = matmul(cell,neighCell(:,j,i))
               Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KptsLoc,Tm))
            end do
         end do
         call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
         if (info/=0) then
            call MIO_Kill('Error in diagonalization for spectral function','diag','DiagSpectralWeight')
         end if
         ! Note that I'm using the Wei Ku expression for the spectral function
         ! (this is NOT the NISHI one)
             do j=1,N
                 do in=1,N
                   if (Species(in).eq.1) then
                     if (layerIndex(in).eq.1) then
                         PkcLocA(j,1) = PkcLocA(j,1) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.2) then
                         PkcLocA(j,2) = PkcLocA(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.3) then
                         PkcLocA(j,3) = PkcLocA(j,3) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.4) then
                         PkcLocA(j,4) = PkcLocA(j,4) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     end if
                   else if (Species(in).eq.2) then
                     if (layerIndex(in).eq.1) then
                         PkcLocB(j,1) = PkcLocB(j,1) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.2) then
                         PkcLocB(j,2) = PkcLocB(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.3) then
                         PkcLocB(j,3) = PkcLocB(j,3) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     else if (layerIndex(in).eq.4) then
                         PkcLocB(j,4) = PkcLocB(j,4) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                     end if
                   else if (Species(in).eq.3) then
                         PkcLocC(j,2) = PkcLocC(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                   else if (Species(in).eq.4) then
                         PkcLocD(j,2) = PkcLocD(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                   else if (Species(in).eq.5) then
                         PkcLocE(j,2) = PkcLocE(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                   else if (Species(in).eq.6) then
                         PkcLocF(j,2) = PkcLocF(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                   else if (Species(in).eq.7) then
                         PkcLocG(j,2) = PkcLocG(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                   else if (Species(in).eq.8) then
                         PkcLocH(j,2) = PkcLocH(j,2) + exp(-cmplx_i*dot_product(KG, -Rat(:,in))) * Hts(in,j)
                   end if
                 end do
             end do
   Nc = 1.0
   nc2 = 1.0
   PkcLocA = Nc/nc2 * PkcLocA
   PkcLocB = Nc/nc2 * PkcLocB

end subroutine DiagSpectralWeightWeiKuInequivalentMoreOrbitals

subroutine DiagSpectralWeightWeiKuInequivalentLee(N,ns,is,PkcLoc1,PkcLoc2,ELoc,KptsLoc,KG,cell,gcell,H0,maxN,hopp,NList,Nneigh,neighCell,topBottomRatio)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac, layerIndex
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh
   use math

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   complex(dp), intent(out) :: PkcLoc1(N,3)
   complex(dp), intent(out) :: PkcLoc2(N,3)
   complex(dp) :: Pkcaux, Pkcaux1, Pkcaux3
   real(dp), intent(out) :: ELoc(N)
   real(dp), intent(in) :: KptsLoc(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: TmVec(N,3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy
   real(dp) :: RatG(3,2)

   real(dp) :: Nc, nc2

   integer :: nTS
   real(dp) :: ll, kk

   integer :: cellSize

   integer :: at1, at2

   character(len=80) :: line
   integer :: id

   real(dp) :: G01(3), G10(3), G11(3)
   real(dp) :: G01p(3), G10p(3), G11p(3)

   complex(dp) :: ZWorkLoc(lwork)
   real(dp) :: DWorkLoc(3*N-2)

   real(dp) :: topBottomRatio

   logical :: changeExpSign

   integer :: ix, iy, ncell(3,9)
   real(dp) :: rMinRp(3), v(3), dmin

   cellts = cell ! Lattice vectors of PC (either top or bottom layer)
   celltsp = gcell

   Hts = 0.0_dp
   PkcLoc1 = 0.0_dp
   PkcLoc2 = 0.0_dp
   if (frac) call AtomsSetCart()
   G10 = matmul(cellts,[1,0,0]) ! here, its in real space!!!
   G01 = matmul(cellts,[0,1,0])
   G11 = matmul(cellts,[1,1,0])
   G10p = matmul(celltsp,[1,0,0]) ! here, its in real space!!!
   G01p = matmul(celltsp,[0,1,0])
   G11p = matmul(celltsp,[1,1,0])
   do i=1,N
       ll = (Rat(1,i)*G10(2)/G10(1) - Rat(2,i)) / (G01(1)*G10(2)/G10(1) - G01(2))
       kk = (Rat(1,i) - ll * G01(1)) / G10(1)
       if (kk.gt.0) then
           kk = floor(kk)
       else
           kk = ceiling(kk)
       end if
       if (ll.gt.0) then
           ll = floor(ll)
       else
           ll = ceiling(ll)
       end if
       RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire

       ll = (Rat(1,i)*G10p(2)/G10p(1) - Rat(2,i)) / (G01p(1)*G10p(2)/G10p(1) - G01p(2))
       kk = (Rat(1,i) - ll * G01p(1)) / G10p(1)
       if (kk.gt.0) then
           kk = floor(kk)
       else
           kk = ceiling(kk)
       end if
       if (ll.gt.0) then
           ll = floor(ll)
       else
           ll = ceiling(ll)
       end if
       RtspVec(i,:) = matmul(celltsp,[int(kk),int(ll),0]) ! Graphene
   end do
   do i=1,N
      Hts(i,i) = H0(i)
      do j=1,Nneigh(i)
         in = NList(j,i)
         Tm = matmul(cell,neighCell(:,j,i))
         Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KptsLoc,Tm))
      end do
   end do
   call ZHEEV('V','L',N,Hts,N,ELoc,ZWorkLoc,lwork,DWorkLoc,info)
   if (info/=0) then
      call MIO_Kill('Error in diagonalization for spectral function','diag','DiagSpectralWeight')
   end if
   ! Note that I'm using the Wei Ku expression for the spectral function
   ! (this is NOT the NISHI one)
   RatG(:,1) = Rat(:,1)
   RatG(:,2) = Rat(:,2)
   do j=1,N
       do in=1,N
           !do ix=-1,1; do iy=-1,1
           !   !if (i==j .and. ix==0 .and. iy==0) cycle

           if (Species(in).eq.1) then
              rMinRp = Rat(:,1)-Rat(:,in)
           else
              rMinRp = Rat(:,2)-Rat(:,in)
           end if
           if (layerIndex(in).eq.1) then
               if (Species(in).eq.1) then
                  PkcLoc1(j,1) = PkcLoc1(j,1) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               else
                  PkcLoc2(j,1) = PkcLoc2(j,1) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               end if
           else if (layerIndex(in).eq.2) then
               if (Species(in).eq.1) then
                  PkcLoc1(j,2) = PkcLoc1(j,2) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               else
                  PkcLoc2(j,2) = PkcLoc2(j,2) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               end if
           else if (layerIndex(in).eq.3) then
               if (Species(in).eq.1) then
                  PkcLoc1(j,3) = PkcLoc1(j,3) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               else
                  PkcLoc2(j,3) = PkcLoc2(j,3) + exp(cmplx_i*dot_product(KG, rMinRp)) * Hts(in,j) * conjg(Hts(in,j))
               end if
           end if
           !       !PkcLoc1(j,1) = PkcLoc1(j,1) + exp(cmplx_i*dot_product(KG, (RatG(:,1) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
           !       !PkcLoc2(j,1) = PkcLoc2(j,1) + exp(cmplx_i*dot_product(KG, (RatG(:,2) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
           !    !PkcLoc(j,1) = PkcLoc(j,1) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
           !       !PkcLoc1(j,2) = PkcLoc1(j,2) + exp(cmplx_i*dot_product(KG, (RatG(:,1) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
           !       !PkcLoc2(j,2) = PkcLoc2(j,2) + exp(cmplx_i*dot_product(KG, (RatG(:,2) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
           !    !PkcLoc(j,2) = PkcLoc(j,2) + exp(-cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) * conjg(Hts(in,j))
           !    !PkcLoc(j,2) = PkcLoc(j,2) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
           !       !PkcLoc1(j,3) = PkcLoc1(j,3) + exp(cmplx_i*dot_product(KG, (RatG(:,1) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
           !       !PkcLoc2(j,3) = PkcLoc2(j,3) + exp(cmplx_i*dot_product(KG, (RatG(:,2) - Rat(:,in)))) * Hts(in,j) * conjg(Hts(in,j))
           !    !PkcLoc(j,3) = PkcLoc(j,3) + exp(-cmplx_i*dot_product(KG, Rat(:,in))) * Hts(in,j) * conjg(Hts(in,j))
           !    !PkcLoc(j,3) = PkcLoc(j,3) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtsVec(in,:)))) * Hts(in,j) * conjg(Hts(in,j))
       end do
   end do
   Nc = 1.0
   nc2 = 1.0
   PkcLoc1 = Nc/nc2 * PkcLoc1
   PkcLoc2 = Nc/nc2 * PkcLoc2

end subroutine DiagSpectralWeightWeiKuInequivalentLee

subroutine DiagSpectralWeightWeiKuInequivalentNishi(N,ns,is,PkcLoc,ELoc1, ELoc2,KptsLoc,KG,cell,gcell1,gcell2,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i
   use interface,             only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf,                   only : charge, Zch
   use atoms,                 only : Species, Rat, AtomsSetCart, AtomsSetFrac, frac, layerIndex
   use tbpar,                 only : U
   use neigh,                 only : maxNeigh

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
   complex(dp) :: Hts(N,N)
   complex(dp) :: Htsp(N,N)
   complex(dp), intent(out) :: PkcLoc(N,2)
   complex(dp) :: Pkcaux, Pkcaux1, Pkcaux3
   real(dp), intent(out) :: ELoc1(N), ELoc2(N)
   real(dp), intent(in) :: KptsLoc(3), cell(3,3), H0(N)
   real(dp), intent(in) :: KG(3)
   real(dp) :: gcell1(3,3)
   real(dp) :: gcell2(3,3)
   real(dp) :: bncell(3,3)
   real(dp) :: cellts(3,3)
   real(dp) :: celltsp(3,3)
   real(dp) :: aG, aBN
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info, inn
   integer :: i1, i2
   real(dp) :: R(3), zz
   real(dp) :: Rts(3), Rtsp(3), Tm(3)
   real(dp) :: RtsVec(N,3), RtspVec(N,3)
   real(dp) :: Rtsx, Rtsy

   real(dp) :: Nc, nc2

   integer :: nTS
   real(dp) :: ll, kk

   integer :: cellSize

   integer :: at1, at2

   character(len=80) :: line
   integer :: id

   real(dp) :: G01(3), G10(3), G11(3)
   real(dp) :: G01p(3), G10p(3), G11p(3)

   complex(dp) :: ZWorkLoc1(lwork)
   complex(dp) :: ZWorkLoc2(lwork)
   real(dp) :: DWorkLoc1(3*N-2)
   real(dp) :: DWorkLoc2(3*N-2)

   cellts = gcell1
   celltsp = gcell2

   Hts = 0.0_dp
   Htsp = 0.0_dp
   PkcLoc = 0.0_dp
   if (frac) call AtomsSetCart()
   G10 = matmul(cellts,[1,0,0]) ! here, its in real space!!!
   G01 = matmul(cellts,[0,1,0])
   G11 = matmul(cellts,[1,1,0])
   G10p = matmul(celltsp,[1,0,0]) ! here, its in real space!!!
   G01p = matmul(celltsp,[0,1,0])
   G11p = matmul(celltsp,[1,1,0])
   do i=1,N
       ll = (Rat(1,i)*G10(2)/G10(1) - Rat(2,i)) / (G01(1)*G10(2)/G10(1) - G01(2))
       kk = (Rat(1,i) - ll * G01(1)) / G10(1)
       if (kk.gt.0) then
           kk = floor(kk)
       else
           kk = ceiling(kk)
       end if
       if (ll.gt.0) then
           ll = floor(ll)
       else
           ll = ceiling(ll)
       end if
       RtsVec(i,:) = matmul(cellts,[int(kk),int(ll),0]) ! Moire

       ll = (Rat(1,i)*G10p(2)/G10p(1) - Rat(2,i)) / (G01p(1)*G10p(2)/G10p(1) - G01p(2))
       kk = (Rat(1,i) - ll * G01p(1)) / G10p(1)
       if (kk.gt.0) then
           kk = floor(kk)
       else
           kk = ceiling(kk)
       end if
       if (ll.gt.0) then
           ll = floor(ll)
       else
           ll = ceiling(ll)
       end if
       RtspVec(i,:) = matmul(celltsp,[int(kk),int(ll),0]) ! Moire
   end do
   do i=1,N
      Hts(i,i) = H0(i)
      Htsp(i,i) = H0(i)
      do j=1,Nneigh(i)
         in = NList(j,i)
         Tm = matmul(cell,neighCell(:,j,i)) ! Corresponds to L in Nishi equation
         Hts(in,i) = Hts(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KptsLoc,Tm-RtsVec(i,:))) ! RtsVec is the multiple of ts
         Htsp(in,i) = Htsp(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(KptsLoc,Tm-RtspVec(i,:)))
      end do
   end do
   call ZHEEV('V','L',N,Hts,N,ELoc1,ZWorkLoc1,lwork,DWorkLoc1,info)
   call ZHEEV('V','L',N,Htsp,N,ELoc2,ZWorkLoc2,lwork,DWorkLoc2,info)
   if (info/=0) then
      call MIO_Kill('Error in diagonalization for spectral function','diag','DiagSpectralWeight')
   end if
   ! This one is NISHI
   do j=1,N ! These are the eigenvectors with band index J
       do in=1,N ! NOT the sum over eigenvectors. Pick one eigenvector and then sum over its coefficients. Each coefficient corresponds to one orbital in Wannier (or TB) basis.
           do inn=1,N
                  if (layerIndex(in).eq.1 .and. layerIndex(inn).eq.1) then
                      PkcLoc(j,1) = PkcLoc(j,1) + exp(cmplx_i*dot_product(KG, (RtsVec(in,:)-RtsVec(inn,:)))) * conjg(Hts(in,j)) * Hts(inn,j)! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
                  else if (layerIndex(in).eq.2 .and. layerIndex(inn).eq.2) then
                      PkcLoc(j,2) = PkcLoc(j,2) + exp(cmplx_i*dot_product(KG, (RtspVec(in,:)-RtspVec(inn,:)))) * conjg(Htsp(in,j)) * Htsp(inn,j)! j is the eigenvector J, is the coeff of orbital N. Order: (in,j)
                  end if
               !       !print*, "hi1"
               !       !PkcLoc(j,2) = PkcLoc(j,2) + exp(cmplx_i*dot_product(KG, RtsVec(in,:))) * Hts(in,j)
               !       !print*, "hi3"
           end do
       end do
   end do
   Nc = 1.0
   nc2 = 1.0
   PkcLoc = Nc/nc2 * PkcLoc

end subroutine DiagSpectralWeightWeiKuInequivalentNishi

function convolve(x, h, Epts)
    implicit none

    integer :: Epts
    !x is the signal array
    !h is the noise/impulse array
    real(dp), allocatable :: convolve(:), y(:)
    real(dp) :: x(Epts), h(Epts)
    integer :: kernelsize, datasize
    integer :: i,j,k

    datasize = size(x)  ! Ake
    kernelsize = size(h) ! gaussian

    allocate(y(datasize))
    allocate(convolve(datasize))

    !last part
    do i=kernelsize,datasize
        y(i) = 0.0
        j=i
        do k=1,kernelsize
            y(i) = y(i) + x(j)*h(k)
            j = j-1
        end do
    end do

    !first part
    do i=1,kernelsize
        y(i) = 0.0
        j=i
        k=1
        do while (j > 0)
            y(i) = y(i) + x(j)*h(k)
            j = j-1
            k = k+1
        end do
    end do

    convolve = y

end function convolve

pure function matinv3(A) result(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(dp), intent(in) :: A(3,3)   !! Matrix
    real(dp)             :: B(3,3)   !! Inverse matrix
    real(dp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
end function

subroutine av (n, v, w)
      integer  ::         n, j
      complex (dp) :: v(n), w(n), rho,  dd, dl, du, s, h, h2
      complex(dp), parameter ::   one = (1.0D+0, 0.0D+0) , two = (2.0D+0, 0.0D+0)
      common            /convct/ rho

      h = one / dcmplx (n+1)
      h2 = h*h
      s = rho / two
      dd = two / h2
      dl = -one/h2 - s/h
      du = -one/h2 + s/h

      w(1) =  dd*v(1) + du*v(2)
      do 10 j = 2,n-1
         w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1)
 10   continue
      w(n) =  dl*v(n-1) + dd*v(n)
      return
      end

subroutine CalculateChernTAPW(Kpts, nk, E, nAt, nspin, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, useSOC)
   use constants, only : pi, twopi
   use math, only : CrossProd, norm
   use atoms, only : Species, layerIndex
   use tbpar, only : g0
   use cell, only : rcell
   implicit none

   ! Input parameters
   integer, intent(in) :: nk, nAt, nspin, maxNeigh
   real(dp), intent(in) :: Kpts(3,nk), E(nAt,nspin,nk)
   real(dp), intent(in) :: ucell(3,3), H0(nAt)
   complex(dp), intent(in) :: hopp(maxNeigh,nAt)
   integer, intent(in) :: NList(maxNeigh,nAt), Nneigh(nAt), neighCell(3,maxNeigh,nAt)
   logical, intent(in), optional :: useSOC  ! Whether SOC is enabled

   ! Local variables for Berry curvature calculation
   integer :: ik, iband, jband, M, Nlabel, NG
   integer :: i, j, is
   real(dp) :: kx, ky, dk, chern_total, chern_target
   real(dp), allocatable :: berry_curv_bands(:)  ! Per-band Berry curvature
   complex(dp), allocatable :: H_k(:,:)  ! Hamiltonians (no numerical derivatives needed)
   complex(dp), allocatable :: eigvec(:,:)  ! Eigenvectors
   real(dp), allocatable :: eigval(:)  ! Eigenvalues
   complex(dp), allocatable :: Vx(:,:), Vy(:,:)  ! Velocity matrix elements
   real(dp), allocatable :: chern_bands(:)  ! Per-band Chern numbers
   integer :: n_bands_near_fermi, max_bands_near_fermi, target_start, target_end
   integer, allocatable :: band_indices(:)
   integer :: fermi_band_index, closest_band
   real(dp) :: min_energy_diff, energy_diff
   integer :: start_band, end_band
   integer :: max_valid_band  ! Maximum valid band index for Fermi search
   real(dp) :: energy_window  ! Energy window around Fermi level (configurable, default ±0.2 eV)
   real(dp) :: energy_window_eV = 0.2_dp  ! Energy window in eV (configurable)
   real(dp) :: fermi_energy_internal  ! Fermi energy in internal units
   logical :: use_all_occupied_bands = .false.  ! Include all bands below Fermi level
   integer :: n_chern_bands = 30  ! Number of specific bands to include (configurable)
   integer :: n_bands_below_fermi = 30  ! Bands below Fermi for convergence (configurable)
   integer :: n_bands_above_fermi = 30  ! Bands above Fermi for convergence (configurable)
   integer :: n_target_bands = 2  ! Number of bands near Fermi to report Chern numbers for
   real(dp) :: min_energy, max_energy  ! For calculating band energy ranges
   real(dp), allocatable :: fermi_distances(:)  ! Distance from Fermi energy for each band
   integer, allocatable :: sorted_indices(:)  ! Indices sorted by distance from Fermi
   integer :: temp_idx  ! Temporary variable for sorting
   real(dp) :: temp_dist, band_energy_at_gamma  ! Temporary variables
   integer, allocatable :: target_band_indices(:)  ! Actual band indices of target bands
   real(dp), allocatable :: target_chern_values(:)  ! Chern numbers of target bands

   ! Variables for BZ area calculation (following Python approach)
   real(dp) :: k_spacing_x, k_spacing_y, actual_BZ_area
   real(dp) :: areaG  ! Reciprocal lattice area

   ! Option to choose between Option A (TB basis) and Option B (TAPW basis) for Berry curvature
   logical :: use_option_a = .false.

   ! SOC flag
   logical :: soc_enabled

   ! Check if SOC is enabled
   soc_enabled = .false.
   if (present(useSOC)) soc_enabled = useSOC

   call MIO_Print('Implementing TAPW Berry curvature calculation using ANALYTICAL derivatives...','diag')
   if (soc_enabled) call MIO_Print('SOC-enabled Chern calculation','diag')
   call MIO_Print('E array dimensions: '//trim(num2str(size(E,1)))//' x '//trim(num2str(size(E,2)))//' x '//trim(num2str(size(E,3))),'diag')
   call MIO_Print('nAt = '//trim(num2str(nAt))//', nk = '//trim(num2str(nk)),'diag')

   ! Store k-points array for Option B Berry curvature calculation
   if (.not. allocated(stored_Kpts)) then
      allocate(stored_Kpts(3, nk))
      stored_Kpts = Kpts(:, 1:nk)  ! Store only the actual k-points (not the extra +1)
      call MIO_Print('Stored k-points array for Option B Berry curvature calculation','diag')

      ! VERIFICATION: Check stored k-points
      call MIO_Print('=== K-POINTS STORAGE VERIFICATION ===','diag')
      call MIO_Print('Stored k-points dimensions: '//trim(num2str(size(stored_Kpts,1)))//'x'//trim(num2str(size(stored_Kpts,2))),'diag')
      call MIO_Print('First k-point: ['//trim(num2str(stored_Kpts(1,1),6))//','//trim(num2str(stored_Kpts(2,1),6))//','//trim(num2str(stored_Kpts(3,1),6))//']','diag')
      if (nk > 1) call MIO_Print('Second k-point: ['//trim(num2str(stored_Kpts(1,2),6))//','//trim(num2str(stored_Kpts(2,2),6))//','//trim(num2str(stored_Kpts(3,2),6))//']','diag')
      call MIO_Print('Last k-point: ['//trim(num2str(stored_Kpts(1,nk),6))//','//trim(num2str(stored_Kpts(2,nk),6))//','//trim(num2str(stored_Kpts(3,nk),6))//']','diag')
   end if

   ! Make energy window and band selection method configurable
   call MIO_InputParameter('Diag.ChernEnergyWindow', energy_window_eV, 0.2_dp)
   call MIO_InputParameter('Diag.ChernAllOccupiedBands', use_all_occupied_bands, .false.)
   call MIO_InputParameter('Diag.ChernNBands', n_chern_bands, 30)

   call MIO_InputParameter('Diag.ChernUseOptionA', use_option_a, .false.)

   if (use_option_a) then
      call MIO_Print('*** USING OPTION A: Berry curvature in TB basis (colleague''s approach) ***','diag')
   else
      call MIO_Print('*** USING OPTION B: Berry curvature in TAPW basis (default approach) ***','diag')
   end if
   call MIO_InputParameter('Diag.ChernBandsBelowFermi', n_bands_below_fermi, 30)
   call MIO_InputParameter('Diag.ChernBandsAboveFermi', n_bands_above_fermi, 30)
   call MIO_InputParameter('Diag.ChernTargetBands', n_target_bands, 2)

   ! Convert Fermi energy from eV to internal units (g0)
   fermi_energy_internal = fermi_energy / g0
   energy_window = energy_window_eV / g0  ! Convert energy window to internal units
   call MIO_Print('Fermi energy: '//trim(num2str(fermi_energy,6))//' eV = '// &
                 trim(num2str(fermi_energy_internal,6))//' (internal units)','diag')
   call MIO_Print('Energy window: ±'//trim(num2str(energy_window*g0,3))//' eV = ±'// &
                 trim(num2str(energy_window,6))//' (internal units)','diag')

   ! For now, implement a placeholder that identifies bands near Fermi energy
   ! and sets up the structure for Berry curvature calculation

   ! Find bands for Chern number calculation
   n_bands_near_fermi = 0
   max_bands_near_fermi = min(size(E,1), max(n_chern_bands, n_bands_below_fermi + n_bands_above_fermi + 10))  ! Use configured number of bands

   ! Debug: Show eigenvalue range at first k-point
   call MIO_Print('Eigenvalue range at first k-point:','diag')
   call MIO_Print('  First 10 eigenvalues: '//trim(num2str(E(1,1,1)*g0,3))//' to '//trim(num2str(E(10,1,1)*g0,3))//' eV','diag')
   call MIO_Print('  Around band 100: '//trim(num2str(E(max(1,100),1,1)*g0,3))//' eV','diag')
   call MIO_Print('  Around band 500: '//trim(num2str(E(min(size(E,1),500),1,1)*g0,3))//' eV','diag')
   call MIO_Print('  Last few: '//trim(num2str(E(size(E,1)-2,1,1)*g0,3))//' to '//trim(num2str(E(size(E,1),1,1)*g0,3))//' eV','diag')

   ! Allocate band indices array
   allocate(band_indices(max_bands_near_fermi))

   if (use_all_occupied_bands) then
      ! Method 1: Include ALL bands below Fermi level (proper total Chern number)
      call MIO_Print('Using ALL occupied bands (E < E_F) for total Chern number','diag')
      do iband = 1, min(size(E,1), 1000)
         ! Check if this band is below Fermi level at ANY k-point
         do ik = 1, nk
            if (E(iband,1,ik) < fermi_energy_internal) then
               n_bands_near_fermi = n_bands_near_fermi + 1
               if (n_bands_near_fermi <= max_bands_near_fermi) then
                  band_indices(n_bands_near_fermi) = iband
                  call MIO_Print('  Including occupied band '//trim(num2str(iband))// &
                                ': E = '//trim(num2str(E(iband,1,ik)*g0,6))//' eV','diag')
               end if
               exit  ! Found this band, move to next
            end if
         end do
         if (n_bands_near_fermi >= max_bands_near_fermi) then
            call MIO_Print('  Maximum number of bands reached: '//trim(num2str(max_bands_near_fermi)),'diag')
            exit
         end if
      end do
   else if (n_bands_below_fermi > 0 .and. n_bands_above_fermi > 0) then
      ! Method 2: Include N bands below and above Fermi for convergent calculation
      ! This is the CORRECT approach for getting accurate Chern numbers of specific bands
      call MIO_Print('Using '//trim(num2str(n_bands_below_fermi))//' bands below + '// &
                    trim(num2str(n_bands_above_fermi))//' bands above Fermi for convergence','diag')
      call MIO_Print('Will report Chern numbers for '//trim(num2str(n_target_bands))//' bands closest to Fermi','diag')

      ! Find the band closest to Fermi energy at first k-point as reference
      ! CRITICAL FIX: Only search within valid TAPW bands [1, stored_M]
      ! The E array may have zero-filled entries beyond M that would incorrectly be selected
      ! when Fermi energy is 0.0

      min_energy_diff = huge(1.0_dp)
      fermi_band_index = 1
      ! Use stored_M if available, otherwise use size(E,1) as fallback
      if (stored_M > 0) then
         max_valid_band = min(stored_M, size(E,1))
      else
         max_valid_band = size(E,1)
      end if

      call MIO_Print('Searching for Fermi band in valid range [1,'//trim(num2str(max_valid_band))//']','diag')
      do iband = 1, max_valid_band
         energy_diff = abs(E(iband,1,1) - fermi_energy_internal)
         if (energy_diff < min_energy_diff) then
            min_energy_diff = energy_diff
            fermi_band_index = iband
         end if
      end do

      call MIO_Print('Band closest to Fermi energy: '//trim(num2str(fermi_band_index))// &
                    ' (E = '//trim(num2str(E(fermi_band_index,1,1)*g0,6))//' eV)','diag')

      ! Include bands in the range [fermi_band - n_below, fermi_band + n_above]
      ! CRITICAL FIX: Ensure end_band doesn't exceed stored_M (valid TAPW bands)
      start_band = max(1, fermi_band_index - n_bands_below_fermi)
      if (stored_M > 0) then
         end_band = min(stored_M, size(E,1), fermi_band_index + n_bands_above_fermi)
      else
         end_band = min(size(E,1), fermi_band_index + n_bands_above_fermi)
      end if

      ! Debug: Show band range calculation
      if (tapwDebug) call MIO_Print('DEBUG: Band range calculation:','diag')
      call MIO_Print('  Requested: '//trim(num2str(n_bands_below_fermi))//' below + '//trim(num2str(n_bands_above_fermi))//' above = '//trim(num2str(n_bands_below_fermi + n_bands_above_fermi + 1))//' total','diag')
      call MIO_Print('  Fermi band index: '//trim(num2str(fermi_band_index)),'diag')
      call MIO_Print('  Total bands available: '//trim(num2str(size(E,1))),'diag')
      call MIO_Print('  Requested start_band: '//trim(num2str(fermi_band_index - n_bands_below_fermi)),'diag')
      call MIO_Print('  Actual start_band: '//trim(num2str(start_band))//' (constrained by max(1, requested))','diag')
      call MIO_Print('  Requested end_band: '//trim(num2str(fermi_band_index + n_bands_above_fermi)),'diag')
      call MIO_Print('  Actual end_band: '//trim(num2str(end_band))//' (constrained by min(size(E,1), requested))','diag')

      n_bands_near_fermi = end_band - start_band + 1
      call MIO_Print('  Final band count: '//trim(num2str(n_bands_near_fermi)),'diag')

      if (n_bands_near_fermi > max_bands_near_fermi) then
         call MIO_Print('WARNING: Requested '//trim(num2str(n_bands_near_fermi))//' bands exceeds limit '// &
                       trim(num2str(max_bands_near_fermi)),'diag')
         n_bands_near_fermi = max_bands_near_fermi
         end_band = start_band + max_bands_near_fermi - 1
      end if

      do i = 1, n_bands_near_fermi
         iband = start_band + i - 1
         ! Bounds check: ensure iband doesn't exceed available bands
         if (iband > size(E,1)) then
            call MIO_Print('  ERROR: Calculated band index '//trim(num2str(iband))//' exceeds E array size '//trim(num2str(size(E,1))),'diag')
            call MIO_Print('    start_band='//trim(num2str(start_band))//', i='//trim(num2str(i))//', n_bands_near_fermi='//trim(num2str(n_bands_near_fermi)),'diag')
            ! Truncate to valid range
            n_bands_near_fermi = i - 1
            exit
         end if
         if (iband > stored_M .and. stored_M > 0) then
            call MIO_Print('  WARNING: Band index '//trim(num2str(iband))//' exceeds stored_M='//trim(num2str(stored_M))//', will be excluded later','diag')
         end if
         band_indices(i) = iband
         if (i <= 5 .or. i > n_bands_near_fermi - 5) then  ! Show first and last 5
            call MIO_Print('  Including band '//trim(num2str(iband))// &
                          ': E = '//trim(num2str(E(iband,1,1)*g0,6))//' eV','diag')
         else if (i == 6) then
            call MIO_Print('  ... ('//trim(num2str(n_bands_near_fermi-10))//' more bands) ...','diag')
         end if
      end do

      call MIO_Print('Total bands for calculation: '//trim(num2str(n_bands_near_fermi)),'diag')
      call MIO_Print('Target bands for Chern reporting: '//trim(num2str(fermi_band_index-n_target_bands/2))// &
                    ' to '//trim(num2str(fermi_band_index+n_target_bands/2)),'diag')

   else if (n_chern_bands > 0 .and. n_chern_bands <= size(E,1)) then
      ! Method 3: Include specific number of lowest bands (for isolated manifolds)
      call MIO_Print('Using first '//trim(num2str(n_chern_bands))//' bands (isolated manifold)','diag')
      n_bands_near_fermi = n_chern_bands
      do i = 1, n_chern_bands
         band_indices(i) = i
         call MIO_Print('  Including band '//trim(num2str(i))// &
                       ': E = '//trim(num2str(E(i,1,1)*g0,6))//' eV (at Γ)','diag')
      end do
   else
      ! Method 3: Include bands in energy window around Fermi level
      call MIO_Print('Using bands within ±'//trim(num2str(energy_window*g0,3))//' eV of Fermi level','diag')
      do iband = 1, min(size(E,1), 1000)
         ! Check this band at ALL k-points to see if it's ever near Fermi energy
         do ik = 1, nk
            if (abs(E(iband,1,ik) - fermi_energy_internal) < energy_window) then
               n_bands_near_fermi = n_bands_near_fermi + 1
               if (n_bands_near_fermi <= max_bands_near_fermi) then
                  band_indices(n_bands_near_fermi) = iband
                  call MIO_Print('  Found band '//trim(num2str(iband))//' at k-point '//trim(num2str(ik))// &
                                ': E = '//trim(num2str(E(iband,1,ik)*g0,6))//' eV','diag')
               end if
               exit  ! Found this band, move to next band
            end if
         end do
         if (n_bands_near_fermi >= max_bands_near_fermi) then
            call MIO_Print('  Maximum number of bands reached: '//trim(num2str(max_bands_near_fermi)),'diag')
            exit
         end if
      end do
   end if

   call MIO_Print('Found '//trim(num2str(n_bands_near_fermi))//' bands near Fermi energy ('// &
                 trim(num2str(fermi_energy,6))//' eV)','diag')

   if (n_bands_near_fermi >= 1) then
      call MIO_Print('Will calculate Chern numbers for '//trim(num2str(n_bands_near_fermi))//' bands','diag')
      if (n_bands_near_fermi >= 2) then
         call MIO_Print('First few bands: '//trim(num2str(band_indices(1)))//' to '//trim(num2str(band_indices(min(n_bands_near_fermi,5)))),'diag')
      end if

      ! Calculate area element for BZ integration (2D) - following Python approach
      ! Python: G = 2*π * inv(R.T), areaG = |cross(G[0,:2], G[1,:2])|, dkx_dky = areaG / (k_x * k_y)
      ! In Fortran: rcell already includes 2π factor, so areaG = |rcell[1] × rcell[2]|
      areaG = norm(CrossProd(rcell(:,1), rcell(:,2)))
      actual_BZ_area = areaG  ! Set the actual BZ area for sanity check
      dk = areaG / (nk_chern_x * nk_chern_y)  ! Area per k-point in reciprocal space

      call MIO_Print('Using proper reciprocal lattice area for integration:','diag')
      call MIO_Print('  Reciprocal lattice area (areaG) = '//trim(num2str(areaG,8)),'diag')
      call MIO_Print('  dk (area per k-point) = '//trim(num2str(dk,8)),'diag')
      call MIO_Print('  Ratio to (2π)²: '//trim(num2str(areaG/((2.0_dp*pi)**2),6)),'diag')

      ! DEBUG: Check integration area
      if (tapwDebug) call MIO_Print('DEBUG: Integration area calculation:','diag')
      call MIO_Print('  nk_chern_x × nk_chern_y = '//trim(num2str(nk_chern_x))//' × '//trim(num2str(nk_chern_y)),'diag')
      call MIO_Print('  dk = '//trim(num2str(dk,8))//' (k-space area per point)','diag')

      ! Initialize Chern numbers for all selected bands (per spin)
      allocate(chern_bands(n_bands_near_fermi))
      chern_bands = 0.0_dp

      ! Allocate Berry curvature workspace for k-point calculations
      allocate(berry_curv_bands(n_bands_near_fermi))

      ! Loop over spin channels to calculate Chern numbers separately for each spin
      ! For non-SOC calculations (nspin=1), this loop executes once
      ! For SOC calculations (nspin=2), this loop executes twice to calculate separately
      do is = 1, nspin
         call MIO_Print('')
         call MIO_Print('=== Calculating Chern numbers for spin channel '//trim(num2str(is))//' ===','diag')

         ! Reset Chern numbers for this spin channel
         chern_bands = 0.0_dp

         ! Identify target bands (closest to Fermi energy) for this spin channel
      allocate(target_band_indices(n_target_bands))
      if (n_target_bands > 0 .and. n_target_bands <= n_bands_near_fermi) then
         ! Find the actual bands closest to Fermi energy from our selected bands
         allocate(fermi_distances(n_bands_near_fermi))
         allocate(sorted_indices(n_bands_near_fermi))

         ! Calculate distance from Fermi energy for each band (at first k-point as reference)
         do i = 1, n_bands_near_fermi
            iband = band_indices(i)
               band_energy_at_gamma = E(iband,is,1)  ! Energy at first k-point for this spin
            fermi_distances(i) = abs(band_energy_at_gamma - fermi_energy_internal)
            sorted_indices(i) = i
         end do

         ! Simple bubble sort to find closest bands
         do i = 1, n_bands_near_fermi - 1
            do j = i + 1, n_bands_near_fermi
               if (fermi_distances(sorted_indices(j)) < fermi_distances(sorted_indices(i))) then
                  temp_idx = sorted_indices(i)
                  sorted_indices(i) = sorted_indices(j)
                  sorted_indices(j) = temp_idx
               end if
            end do
         end do

         ! Store the indices of the closest bands for file naming
         do i = 1, n_target_bands
            target_band_indices(i) = band_indices(sorted_indices(i))
         end do

         ! Debug: Show target band indices
            call MIO_Print('DEBUG: Target band indices for spin '//trim(num2str(is))//' (file naming):','diag')
         do i = 1, n_target_bands
            call MIO_Print('  target_band_indices('//trim(num2str(i))//') = '//trim(num2str(target_band_indices(i))),'diag')
         end do

         deallocate(fermi_distances, sorted_indices)
      else
         ! Fallback: use first n_target_bands from convergence bands
         do i = 1, n_target_bands
            target_band_indices(i) = band_indices(i)
         end do
      end if

         ! Main Berry curvature calculation loop for this spin channel
         call MIO_Print('Calculating Berry curvature at '//trim(num2str(nk))//' k-points for spin '//trim(num2str(is))//'...','diag')
      call MIO_Print('Processing '//trim(num2str(n_bands_near_fermi))//' bands','diag')

      ! Area element already calculated above

      ! Memory-efficient: calculate Berry curvature one k-point at a time
      do ik = 1, nk
         if (modulo(ik, max(1, nk/4)) == 0) then
            call MIO_Print('Berry curvature progress: '//trim(num2str(nint(100.0_dp*ik/nk)))//'%','diag')
         end if

         ! Calculate Berry curvature for all selected bands at this k-point
         if (ik == 1) then  ! Debug info for first k-point only
            call MIO_Print('Debug: Calling CalculateBerryAtKpointFromStored with:','diag')
            call MIO_Print('  n_bands_near_fermi: '//trim(num2str(n_bands_near_fermi)),'diag')
            call MIO_Print('  band_indices size: '//trim(num2str(size(band_indices(1:n_bands_near_fermi)))),'diag')
            call MIO_Print('  berry_curv_bands size: '//trim(num2str(size(berry_curv_bands))),'diag')
         end if
         ! Choose between Option A and Option B (pass spin index)
         if (use_option_a) then
            ! DEBUG: Verify inputs to Option A
            if (ik == 1 .and. is == 1) then
               call MIO_Print('=== OPTION A INPUT VERIFICATION ===','diag')
               call MIO_Print('kpoint_index: '//trim(num2str(ik)),'diag')
               call MIO_Print('spin_index: '//trim(num2str(is)),'diag')
               call MIO_Print('n_bands_near_fermi: '//trim(num2str(n_bands_near_fermi)),'diag')
               call MIO_Print('band_indices(1:5): ['//trim(num2str(band_indices(1)))//','//trim(num2str(band_indices(2)))//','//trim(num2str(band_indices(3)))//','//trim(num2str(band_indices(4)))//','//trim(num2str(band_indices(5)))//']','diag')
               call MIO_Print('E(:,is,ik) dimensions: '//trim(num2str(size(E,1))),'diag')
               call MIO_Print('E(:,is,ik) first 5 values: ['//trim(num2str(E(1,is,ik),6))//','//trim(num2str(E(2,is,ik),6))//','//trim(num2str(E(3,is,ik),6))//','//trim(num2str(E(4,is,ik),6))//','//trim(num2str(E(5,is,ik),6))//']','diag')
            end if
            if (soc_enabled) then
               call CalculateBerryAtKpointFromStored_OptionA_withSOC(ik, is, band_indices(1:n_bands_near_fermi), E(:,is,ik), berry_curv_bands)
            else
               call CalculateBerryAtKpointFromStored_OptionA(ik, is, band_indices(1:n_bands_near_fermi), E(:,is,ik), berry_curv_bands)
            end if
         else
            call CalculateBerryAtKpointFromStored(ik, is, band_indices(1:n_bands_near_fermi), E(:,is,ik), berry_curv_bands)
         end if

         ! DEBUG: Verify k-point ordering and coordinates (only if tapwDebug enabled)
         if (tapwDebug .and. (ik <= 5 .or. modulo(ik, max(1, nk/10)) == 0)) then
            call MIO_Print('DEBUG k-point '//trim(num2str(ik))//': kx='//trim(num2str(Kpts(1,ik),6))//', ky='//trim(num2str(Kpts(2,ik),6))//', kz='//trim(num2str(Kpts(3,ik),6)),'diag')
         end if

         ! CRITICAL: Verify k-point correspondence between storage and retrieval (only if tapwDebug enabled)
         if (tapwDebug .and. ik <= 3) then
            call MIO_Print('=== K-POINT CORRESPONDENCE CHECK ===','diag')
            call MIO_Print('Chern calculation k-point '//trim(num2str(ik))//':','diag')
            call MIO_Print('  Kpts(:,ik) = ['//trim(num2str(Kpts(1,ik),6))//','//trim(num2str(Kpts(2,ik),6))//','//trim(num2str(Kpts(3,ik),6))//']','diag')
            if (allocated(stored_Kpts) .and. ik <= size(stored_Kpts,2)) then
               call MIO_Print('  stored_Kpts(:,ik) = ['//trim(num2str(stored_Kpts(1,ik),6))//','//trim(num2str(stored_Kpts(2,ik),6))//','//trim(num2str(stored_Kpts(3,ik),6))//']','diag')
               ! Check if they match
               if (abs(Kpts(1,ik) - stored_Kpts(1,ik)) < 1e-10_dp .and. &
                   abs(Kpts(2,ik) - stored_Kpts(2,ik)) < 1e-10_dp .and. &
                   abs(Kpts(3,ik) - stored_Kpts(3,ik)) < 1e-10_dp) then
                  call MIO_Print('  ✓ K-points match perfectly','diag')
               else
                  call MIO_Print('  ✗ K-POINT MISMATCH DETECTED!','diag')
                  call MIO_Print('  ✗ This will cause wrong Berry curvature calculation!','diag')
               end if
            else
               call MIO_Print('  WARNING: stored_Kpts not available or index out of bounds','diag')
            end if
         end if

         ! OUTPUT: Write Berry curvature data for hotspot analysis
         call output_berry_curvature_data(ik, Kpts(:,ik), berry_curv_bands, target_band_indices, n_target_bands)

         ! Accumulate Chern numbers (BZ integration) for all bands
         ! Python formula: C = np.sum(berry, axis=0) * dkx_dky / (2*np.pi)
         ! Now using proper reciprocal lattice area, so we need the /2π factor
         do i = 1, n_bands_near_fermi
            chern_bands(i) = chern_bands(i) + berry_curv_bands(i) * dk / (2.0_dp * pi)
         end do
         end do  ! End of k-point loop

         ! Calculate total and target Chern numbers for this spin channel
      chern_total = sum(chern_bands)

         ! Report Chern numbers for target bands (closest to Fermi) for this spin channel
         call MIO_Print('=== Chern Number Results for Spin Channel '//trim(num2str(is))//' ===','diag')
      if (use_option_a) then
         call MIO_Print('*** RESULTS FROM OPTION A (TB basis approach) ***','diag')
      else
         call MIO_Print('*** RESULTS FROM OPTION B (TAPW basis approach) ***','diag')
      end if
      call MIO_Print('Convergence calculation used '//trim(num2str(n_bands_near_fermi))//' bands','diag')

      ! Find target bands (n_target_bands closest to Fermi energy)
      if (n_target_bands > 0 .and. n_target_bands <= n_bands_near_fermi) then
         ! Find the actual bands closest to Fermi energy from our selected bands
         ! Create array to store energy differences from Fermi level

         allocate(fermi_distances(n_bands_near_fermi))
         allocate(sorted_indices(n_bands_near_fermi))

         ! Calculate distance from Fermi energy for each band (at first k-point as reference) for this spin
         do i = 1, n_bands_near_fermi
            iband = band_indices(i)
            band_energy_at_gamma = E(iband,is,1)  ! Energy at first k-point for this spin
            fermi_distances(i) = abs(band_energy_at_gamma - fermi_energy_internal)
            sorted_indices(i) = i
         end do

         ! Simple bubble sort to find closest bands
         do i = 1, n_bands_near_fermi - 1
            do j = i + 1, n_bands_near_fermi
               if (fermi_distances(sorted_indices(j)) < fermi_distances(sorted_indices(i))) then
                  temp_idx = sorted_indices(i)
                  sorted_indices(i) = sorted_indices(j)
                  sorted_indices(j) = temp_idx
               end if
            end do
         end do

         ! Store the indices of the closest bands for later use
         ! We'll process them individually rather than as a range
         target_start = 1  ! Start from first closest band
         target_end = n_target_bands  ! Process n_target_bands closest bands

         ! Debug: Show which bands are actually closest
         if (tapwDebug) call MIO_Print('DEBUG: Bands closest to Fermi energy (E_F = '//trim(num2str(fermi_energy,6))//' eV) for spin '//trim(num2str(is))//':','diag')
         do i = 1, min(n_target_bands, 5)
            iband = band_indices(sorted_indices(i))
            band_energy_at_gamma = E(iband,is,1) * g0
            call MIO_Print('  Rank '//trim(num2str(i))//': Band '//trim(num2str(iband))// &
                          ', E_gamma = '//trim(num2str(band_energy_at_gamma,6))// &
                          ' eV, |E - E_F| = '//trim(num2str(fermi_distances(sorted_indices(i))*g0,6))//' eV','diag')
         end do

         call MIO_Print('Target bands ('//trim(num2str(n_target_bands))//' closest to Fermi):','diag')
         if (tapwDebug) call MIO_Print('  DEBUG: E array dimensions: '//trim(num2str(size(E,1)))//' x '//trim(num2str(size(E,2)))//' x '//trim(num2str(size(E,3))),'diag')

         ! Calculate min/max energies of target bands across all k-points
         call MIO_Print('  Energy ranges of target bands:','diag')
         do i = 1, n_target_bands
            ! Bounds check: ensure sorted_indices(i) is valid
            if (sorted_indices(i) < 1 .or. sorted_indices(i) > n_bands_near_fermi) then
               call MIO_Print('  ERROR: sorted_indices('//trim(num2str(i))//') = '//trim(num2str(sorted_indices(i)))// &
                            ' out of range [1,'//trim(num2str(n_bands_near_fermi))//']','diag')
               cycle
            end if
            iband = band_indices(sorted_indices(i))
            ! Additional bounds check: ensure iband is valid
            if (iband < 1 .or. iband > stored_M) then
               call MIO_Print('  ERROR: Band index '//trim(num2str(iband))//' out of range [1,'//trim(num2str(stored_M))//'] for target band '//trim(num2str(i)),'diag')
               call MIO_Print('    This indicates band_indices('//trim(num2str(sorted_indices(i)))//') = '//trim(num2str(iband))//' is invalid','diag')
               cycle
            end if
            ! Find min and max energies across all k-points for this band
            ! Extract eigenvalues from stored Hamiltonians instead of using E array
            call GetBandEnergyRange(iband, min_energy, max_energy, is)

            call MIO_Print('  Band '//trim(num2str(iband))//': E_min = '// &
                          trim(num2str(min_energy*g0,6))//' eV, E_max = '// &
                          trim(num2str(max_energy*g0,6))//' eV','diag')
         end do

         chern_target = 0.0_dp
         do i = 1, n_target_bands
            call MIO_Print('  Band '//trim(num2str(band_indices(sorted_indices(i))))//': C = '// &
                          trim(num2str(chern_bands(sorted_indices(i)),6)),'diag')
            chern_target = chern_target + chern_bands(sorted_indices(i))
         end do

         ! Store target band information for file output (target_band_indices already allocated earlier)
         allocate(target_chern_values(n_target_bands))
         do i = 1, n_target_bands
            target_chern_values(i) = chern_bands(sorted_indices(i))
         end do

         deallocate(fermi_distances, sorted_indices)
         call MIO_Print('  Target Chern number: C_target = '//trim(num2str(chern_target,6)),'diag')
      else
         ! Fallback: report first few bands
         call MIO_Print('Selected bands (first '//trim(num2str(min(n_bands_near_fermi,5)))//' shown):','diag')

         ! Calculate min/max energies of selected bands across all k-points
         call MIO_Print('  Energy ranges of selected bands:','diag')
         do i = 1, min(n_bands_near_fermi, 5)
            iband = band_indices(i)
            ! Find min and max energies across all k-points for this band
            ! Extract eigenvalues from stored Hamiltonians instead of using E array
            call GetBandEnergyRange(iband, min_energy, max_energy, is)
            call MIO_Print('  Band '//trim(num2str(iband))//': E_min = '// &
                          trim(num2str(min_energy*g0,6))//' eV, E_max = '// &
                          trim(num2str(max_energy*g0,6))//' eV','diag')
         end do

         do i = 1, min(n_bands_near_fermi, 5)
            call MIO_Print('  Band '//trim(num2str(band_indices(i)))//': C = '// &
                          trim(num2str(chern_bands(i),6)),'diag')
         end do
         if (n_bands_near_fermi > 5) then
            call MIO_Print('  ... ('//trim(num2str(n_bands_near_fermi-5))//' more bands)','diag')
         end if
      end if

         call MIO_Print('  Total Chern number (all bands) for spin '//trim(num2str(is))//': C_total = '//trim(num2str(chern_total,6)),'diag')

      ! Validation: Check for reasonable Chern number values
         call MIO_Print('=== Validation Checks for Spin '//trim(num2str(is))//' ===','diag')
      if (abs(chern_total) > 10.0_dp) then
            call MIO_Print('WARNING: Total Chern number unusually large for spin '//trim(num2str(is))//': '//trim(num2str(chern_total,6)),'diag')
      end if

      ! Check for NaN or infinite values
      do i = 1, n_bands_near_fermi
         if (chern_bands(i) /= chern_bands(i)) then  ! NaN check
               call MIO_Print('ERROR: NaN detected in band '//trim(num2str(band_indices(i)))//' for spin '//trim(num2str(is)),'diag')
         end if
      end do

         call MIO_Print('Multi-band Chern calculation completed successfully for spin '//trim(num2str(is)),'diag')

         ! Write results to output file - use target bands instead of all bands (include spin in filename)
      if (n_target_bands > 0 .and. n_target_bands <= n_bands_near_fermi .and. allocated(target_band_indices)) then
         ! Pass only the target bands (closest to Fermi) to file output
         call write_chern_tapw_results(target_band_indices, target_chern_values, &
                                         chern_total, nk_chern_x, nk_chern_y, fermi_energy, is)
         deallocate(target_band_indices, target_chern_values)
      else
         ! Fallback: use first 2 bands if target bands not properly defined
         call write_chern_tapw_results(band_indices(1:2), chern_bands(1:2), &
                                         chern_total, nk_chern_x, nk_chern_y, fermi_energy, is)
      end if
      end do  ! End of spin loop

      deallocate(chern_bands, berry_curv_bands)
   else
      call MIO_Print('Warning: Not enough bands found near Fermi energy for Chern calculation','diag')
   end if

   call MIO_Print('Chern number calculation completed.','diag')

   ! Clean up stored TAPW data and position difference cache
   call cleanup_TAPW_storage()
   call cleanup_position_differences_cache()

end subroutine CalculateChernTAPW

subroutine CalculateBerryAtKpoint(kpt, band_indices, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, berry_curv)
   ! Calculate Berry curvature at a single k-point using analytical derivatives (following Python dHdk method)
   use constants, only : pi, twopi, cmplx_i
   use math
   use atoms, only : nAt
   implicit none

   ! Input parameters
   real(dp), intent(in) :: kpt(3)
   integer, intent(in) :: band_indices(2)
   real(dp), intent(in) :: ucell(3,3), H0(nAt)
   integer, intent(in) :: maxNeigh
   complex(dp), intent(in) :: hopp(maxNeigh,nAt)
   integer, intent(in) :: NList(maxNeigh,nAt), Nneigh(nAt), neighCell(3,maxNeigh,nAt)
   real(dp), intent(out) :: berry_curv

   ! Local variables
   complex(dp), allocatable :: H_k(:,:)  ! TAPW projected Hamiltonian
   complex(dp), allocatable :: eigvec(:,:)  ! Eigenvectors
   real(dp), allocatable :: eigval(:)  ! Eigenvalues
   complex(dp), allocatable :: Vx(:,:), Vy(:,:)  ! Velocity matrix elements
   complex(dp), allocatable :: dHdkx(:,:), dHdky(:,:)  ! Hamiltonian derivatives (analytical)
   complex(dp), allocatable :: temp_x(:,:), temp_y(:,:)  ! Temporary matrices for ZGEMM
   real(dp), allocatable :: delX(:,:), delY(:,:)  ! Position difference matrices
   integer :: i, j, n, m
   real(dp) :: eps = 1.0e-10_dp  ! Small energy denominator regularization
   complex(dp) :: berry_sum

   ! Get TAPW projected Hamiltonian, eigenvalues, eigenvectors, and position differences
   call GetTAPWHamiltonian(kpt, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, &
                          H_k, eigval, eigvec, M, delX, delY)

   ! Calculate analytical Hamiltonian derivatives following Python dHdk method:
   ! dH/dkx = -i * delX .* H(k), dH/dky = -i * delY .* H(k) (element-wise)
   ! CORRECTED: Use element-wise multiplication to match Python's .multiply() method
   allocate(dHdkx(M,M), dHdky(M,M))
   do i = 1, M
      do j = 1, M
         dHdkx(i,j) = -cmplx_i * delX(i,j) * H_k(i,j)
         dHdky(i,j) = -cmplx_i * delY(i,j) * H_k(i,j)
      end do
   end do

   ! Calculate velocity matrix elements: V^x_nm = <n|dH/dkx|m>, V^y_nm = <n|dH/dky|m>
   ! Need to be careful with matrix dimensions to avoid bounds errors
   allocate(Vx(M,M), Vy(M,M))

   ! Use two-step ZGEMM to avoid temporary array issues:
   ! Step 1: temp = dHdkx * eigvec
   ! Step 2: Vx = eigvec^† * temp
   allocate(temp_x(M,M), temp_y(M,M))

   call zgemm('N', 'N', M, M, M, (1.0_dp,0.0_dp), dHdkx, M, eigvec, M, (0.0_dp,0.0_dp), temp_x, M)
   call zgemm('C', 'N', M, M, M, (1.0_dp,0.0_dp), eigvec, M, temp_x, M, (0.0_dp,0.0_dp), Vx, M)

   call zgemm('N', 'N', M, M, M, (1.0_dp,0.0_dp), dHdky, M, eigvec, M, (0.0_dp,0.0_dp), temp_y, M)
   call zgemm('C', 'N', M, M, M, (1.0_dp,0.0_dp), eigvec, M, temp_y, M, (0.0_dp,0.0_dp), Vy, M)

   deallocate(temp_x, temp_y)

   ! Set diagonal elements to zero (Python: np.fill_diagonal(Vxnm, 0, wrap=True))
   do i = 1, M
      Vx(i,i) = (0.0_dp, 0.0_dp)
      Vy(i,i) = (0.0_dp, 0.0_dp)
   end do

   ! Arrays are correctly dimensioned - proceed with Berry curvature calculation

   ! Calculate Berry curvature using exact Python formula:
   ! Berry[ik] = -2*imag(sum((Vx*(Vy.T))/((eigval[n] - eigval[m])**2 + eps), axis=1))
   berry_sum = (0.0_dp, 0.0_dp)
   ! Use actual array dimensions instead of variable M to avoid any modification issues
   do n = 1, size(Vx,1)
      do m = 1, size(Vx,2)
         if (n /= m) then
            ! CORRECTED Berry curvature formula: Vx(n,m) * Vy(m,n) (no conjg!)
            ! Python uses (Vynm.T) which means transpose, not conjugate transpose
            berry_sum = berry_sum + Vx(n,m) * Vy(m,n) / ((eigval(n) - eigval(m))**2 + eps)
         end if
      end do
   end do

   ! Extract imaginary part and apply -2 factor (exact Python formula)
   berry_curv = -2.0_dp * aimag(berry_sum)

   ! Cleanup
   deallocate(H_k, eigvec, eigval, dHdkx, dHdky, Vx, Vy, delX, delY)

end subroutine CalculateBerryAtKpoint

subroutine CalculateBerryAtKpointFromStored(kpoint_index, spin_index, band_indices, eigval_from_E, berry_curv_bands)
   ! Calculate Berry curvature using stored eigenvectors and Hamiltonians
   ! Now returns per-band Berry curvature for the specified bands (variable number)
   ! Updated to include spin_index for separate spin channel calculations
   use constants, only : pi, twopi, cmplx_i
   use atoms, only : nAt
   use ham, only : H0, hopp
   use neigh, only : NList, Nneigh, neighCell, maxNeigh
   implicit none

   ! Input parameters
   integer, intent(in) :: kpoint_index, spin_index
   integer, intent(in) :: band_indices(:)  ! Variable number of band indices
   real(dp), intent(in) :: eigval_from_E(:)  ! Eigenvalues from E array
   real(dp), intent(out) :: berry_curv_bands(:)  ! Per-band Berry curvature (variable size)

   ! Local variables
   complex(dp), allocatable :: H_k(:,:), eigvec(:,:), Vx(:,:), Vy(:,:)
   complex(dp), allocatable :: dHdkx(:,:), dHdky(:,:), temp_x(:,:), temp_y(:,:)
   real(dp), allocatable :: delX(:,:), delY(:,:), eigval(:)
   integer :: i, j, n, m, M_safe
   real(dp) :: energy_threshold, eta_eV, energy_diff
   real(dp) :: eps = 1.0e-10_dp  ! Match Python regularization parameter
   complex(dp) :: berry_sum

   ! Variables for Option B implementation
   real(dp) :: current_kpt(3)
   complex(dp), allocatable :: dH_TB_dkx(:,:), dH_TB_dky(:,:)

   ! Check if stored data is available
   if (.not. allocated(stored_eigenvectors) .or. .not. allocated(stored_hamiltonians)) then
      call MIO_Print('Error: No stored TAPW data available for Chern calculation','diag')
      berry_curv_bands = 0.0_dp
      return
   end if

   if (kpoint_index > stored_nk .or. kpoint_index < 1) then
      call MIO_Print('Error: k-point index '//trim(num2str(kpoint_index))//' out of range [1,'//trim(num2str(stored_nk))//']','diag')
      berry_curv_bands = 0.0_dp
      return
   end if

   ! Check band indices bounds for all bands
   do i = 1, size(band_indices)
      if (band_indices(i) < 1 .or. band_indices(i) > stored_M) then
         call MIO_Print('Error: Band index '//trim(num2str(band_indices(i)))//' out of range [1,'//trim(num2str(stored_M))//']','diag')
         berry_curv_bands = 0.0_dp
         return
      end if
   end do

   ! Check output array size matches input
   if (size(berry_curv_bands) /= size(band_indices)) then
      call MIO_Print('Error: Output array size mismatch','diag')
      call MIO_Print('  berry_curv_bands size: '//trim(num2str(size(berry_curv_bands))),'diag')
      call MIO_Print('  band_indices size: '//trim(num2str(size(band_indices))),'diag')
      berry_curv_bands = 0.0_dp
      return
   end if

   ! FIXED: Use actual stored array dimensions instead of stored_M to avoid bounds errors
   ! Check if stored arrays exist and get their actual dimensions
   if (.not. allocated(stored_eigenvectors) .or. .not. allocated(stored_hamiltonians)) then
      call MIO_Print('ERROR: Stored TAPW arrays not allocated','diag')
      berry_curv_bands = 0.0_dp
      return
   end if

   ! CRITICAL FIX: Always use actual stored array dimensions, not stored_M
   ! stored_M might be corrupted or inconsistent with actual array allocation
   M = min(size(stored_eigenvectors,1), size(stored_eigenvectors,2))

   ! Debug: Check for stored_M inconsistency
   if (kpoint_index == 1) then
      if (tapwDebug) call MIO_Print('CRITICAL DEBUG: Dimension analysis','diag')
      call MIO_Print('  stored_M value: '//trim(num2str(stored_M)),'diag')
      call MIO_Print('  stored_eigenvectors actual size: '//trim(num2str(size(stored_eigenvectors,1)))//'x'//trim(num2str(size(stored_eigenvectors,2))),'diag')
      call MIO_Print('  stored_hamiltonians actual size: '//trim(num2str(size(stored_hamiltonians,1)))//'x'//trim(num2str(size(stored_hamiltonians,2))),'diag')
      call MIO_Print('  Using safe M = '//trim(num2str(M))//' (actual array dimensions)','diag')

      if (stored_M /= M) then
         call MIO_Print('ERROR: stored_M ('//trim(num2str(stored_M))//') inconsistent with actual array size ('//trim(num2str(M))//')','diag')
         call MIO_Print('This indicates a serious bug in array allocation - using actual array size for safety','diag')
      end if
   end if

   ! Debug: Check stored array dimensions
   if (kpoint_index == 1) then
      call MIO_Print('Debug: Stored array info:','diag')
      call MIO_Print('  stored_M: '//trim(num2str(stored_M)),'diag')
      call MIO_Print('  stored_nk: '//trim(num2str(stored_nk)),'diag')
      call MIO_Print('  stored_eigenvectors shape: '//trim(num2str(size(stored_eigenvectors,1)))//'x'//trim(num2str(size(stored_eigenvectors,2)))//'x'//trim(num2str(size(stored_eigenvectors,3))),'diag')
      call MIO_Print('  eigval_from_E size: '//trim(num2str(size(eigval_from_E))),'diag')
      call MIO_Print('  Max band index requested: '//trim(num2str(maxval(band_indices))),'diag')
      call MIO_Print('  Using M = '//trim(num2str(M))//' (min of stored_M and actual array dimensions)','diag')
   end if

   ! Get stored data for this k-point with safe dimensions
   allocate(H_k(M,M), eigvec(M,M), eigval(M))
   allocate(delX(M,M), delY(M,M))
   allocate(dHdkx(M,M), dHdky(M,M))
   allocate(Vx(M,M), Vy(M,M))
   allocate(temp_x(M,M), temp_y(M,M))

   ! Safe array assignment with bounds checking
   ! Check spin_index bounds
   if (spin_index < 1 .or. spin_index > stored_nspin) then
      call MIO_Print('Error: spin_index '//trim(num2str(spin_index))//' out of range [1,'//trim(num2str(stored_nspin))//']','diag')
      berry_curv_bands = 0.0_dp
      return
   end if

   H_k(1:M,1:M) = stored_hamiltonians(1:M,1:M,kpoint_index,spin_index)
   eigvec(1:M,1:M) = stored_eigenvectors(1:M,1:M,kpoint_index,spin_index)

   ! Debug: Check eigenvalue array bounds
   if (kpoint_index == 1) then
      call MIO_Print('Debug: Eigenvalue assignment:','diag')
      call MIO_Print('  M (stored_M): '//trim(num2str(M)),'diag')
      call MIO_Print('  eigval_from_E size: '//trim(num2str(size(eigval_from_E))),'diag')
      call MIO_Print('  min(M, size(eigval_from_E)): '//trim(num2str(min(M,size(eigval_from_E)))),'diag')
   end if

   ! CRITICAL FIX: Use TAPW eigenvalues instead of TB eigenvalues
   ! The eigenvectors are from TAPW diagonalization, so eigenvalues must be too
   if (allocated(stored_eigenvalues)) then
      if (size(stored_eigenvalues,1) >= M .and. size(stored_eigenvalues,2) >= kpoint_index .and. size(stored_eigenvalues,3) >= spin_index) then
         eigval(1:M) = stored_eigenvalues(1:M,kpoint_index,spin_index)
         call MIO_Print('Using stored TAPW eigenvalues for k-point '//trim(num2str(kpoint_index)),'diag')
         call MIO_Print('  TAPW eigenvalues range: ['//trim(num2str(minval(eigval(1:M)),6))//','//trim(num2str(maxval(eigval(1:M)),6))//']','diag')
      else
         call MIO_Print('ERROR: stored_eigenvalues array too small','diag')
         call MIO_Print('  Array size: '//trim(num2str(size(stored_eigenvalues,1)))//'x'//trim(num2str(size(stored_eigenvalues,2))),'diag')
         call MIO_Print('  Required: '//trim(num2str(M))//'x'//trim(num2str(kpoint_index)),'diag')
         error stop 1
      end if
   else
      call MIO_Print('ERROR: No stored TAPW eigenvalues available - using TB eigenvalues (WRONG!)','diag')
      call MIO_Print('This will cause eigenvector-eigenvalue misalignment!','diag')
      call MIO_Print('  TB eigenvalues range: ['//trim(num2str(minval(eigval_from_E),6))//','//trim(num2str(maxval(eigval_from_E),6))//']','diag')
      ! Fallback to TB eigenvalues (this is wrong but prevents crash)
   if (size(eigval_from_E) >= M) then
      eigval(1:M) = eigval_from_E(1:M)
   else
      eigval(1:size(eigval_from_E)) = eigval_from_E
         eigval(size(eigval_from_E)+1:M) = 0.0_dp
      end if
   end if

   ! VERIFICATION: Check eigenvalues for first k-point
   if (kpoint_index == 1) then
      call MIO_Print('=== EIGENVALUES VERIFICATION (k-point 1) ===','diag')
      call MIO_Print('Eigenvalues dimensions: '//trim(num2str(size(eigval))),'diag')
      call MIO_Print('First 5 eigenvalues: ['//trim(num2str(eigval(1),6))//','//trim(num2str(eigval(2),6))//','//trim(num2str(eigval(3),6))//','//trim(num2str(eigval(4),6))//','//trim(num2str(eigval(5),6))//']','diag')
      call MIO_Print('Eigenvalue range: ['//trim(num2str(minval(eigval),6))//','//trim(num2str(maxval(eigval),6))//']','diag')
   end if

   ! SANITY CHECK: Verify eigenvector-eigenvalue alignment (only if tapwDebug enabled)
   if (tapwDebug) then
      call MIO_Print('=== EIGENVECTOR-EIGENVALUE ALIGNMENT CHECK ===','diag')
      call verify_eigenvector_eigenvalue_alignment(H_k, eigvec, eigval, M, kpoint_index)
   end if

   ! OPTION B: Compute TB Hamiltonian derivatives and project to TAPW space
   ! This is the correct approach: dH^TAPW/dk = X† * (dH^TB/dk) * X

   ! Check if we have stored X matrix (TB parameters are now module-level)
   if (.not. allocated(stored_X_matrix)) then
      call MIO_Print('ERROR: Stored X matrix not available for Option B Berry curvature','diag')
      call MIO_Print('Falling back to old method (incorrect but functional)','diag')

      ! Fallback to old method
      call build_TAPW_position_differences(delX, delY, M, saved_Gx, saved_Gy, saved_NG, saved_Nlabel)
      do i = 1, M
         do j = 1, M
            dHdkx(i,j) = -cmplx_i * delX(i,j) * H_k(i,j)
            dHdky(i,j) = -cmplx_i * delY(i,j) * H_k(i,j)
         end do
      end do
   else
      call MIO_Print('Using Option B: Computing TB derivatives and projecting to TAPW space','diag')

      ! Get actual k-point coordinates from stored k-points array
      if (allocated(stored_Kpts) .and. kpoint_index <= size(stored_Kpts, 2)) then
         current_kpt = stored_Kpts(:, kpoint_index)
         if (kpoint_index == 1) then
            call MIO_Print('Using actual k-point coordinates: ['//trim(num2str(current_kpt(1),6))//','//&
                          trim(num2str(current_kpt(2),6))//','//trim(num2str(current_kpt(3),6))//']','diag')
         end if
      else
         call MIO_Print('WARNING: Cannot access stored k-point, using Gamma point','diag')
         current_kpt = [0.0_dp, 0.0_dp, 0.0_dp]
      end if
      allocate(dH_TB_dkx(nAt, nAt), dH_TB_dky(nAt, nAt))

      ! Step 1: Compute TB Hamiltonian derivatives (use module-level variables directly)
      call compute_TB_hamiltonian_derivatives(dH_TB_dkx, dH_TB_dky, nAt, current_kpt, &
                                            H0, maxNeigh, hopp, NList, Nneigh, neighCell)

      ! Step 2: Project TB derivatives to TAPW space
      call project_TB_derivatives_to_TAPW(dH_TB_dkx, dH_TB_dky, stored_X_matrix, &
                                         nAt, M, dHdkx, dHdky)

      ! DEBUG: Compare magnitudes before and after projection
      if (kpoint_index == 1) then
         if (tapwDebug) call MIO_Print('DEBUG: Derivative magnitudes comparison:','diag')
         call MIO_Print('  |dH_TB/dkx| max = '//trim(num2str(maxval(abs(dH_TB_dkx)),8)),'diag')
         call MIO_Print('  |dH_TB/dky| max = '//trim(num2str(maxval(abs(dH_TB_dky)),8)),'diag')
         call MIO_Print('  |dH_TAPW/dkx| max = '//trim(num2str(maxval(abs(dHdkx)),8)),'diag')
         call MIO_Print('  |dH_TAPW/dky| max = '//trim(num2str(maxval(abs(dHdky)),8)),'diag')
         call MIO_Print('  Amplification factor ≈ '//trim(num2str(maxval(abs(dHdkx))/max(maxval(abs(dH_TB_dkx)),1e-12_dp),2)),'diag')
      end if

      ! Clean up TB derivative matrices
      deallocate(dH_TB_dkx, dH_TB_dky)

      call MIO_Print('Successfully computed Option B derivatives: dH^TAPW/dk = X† * (dH^TB/dk) * X','diag')
   end if

   ! Calculate velocity matrix elements: V = eigvec^† * dH/dk * eigvec
   ! Step 1: temp = dH/dk * eigvec
   call ZGEMM('N', 'N', M, M, M, (1.0_dp,0.0_dp), dHdkx, M, eigvec, M, (0.0_dp,0.0_dp), temp_x, M)
   call ZGEMM('N', 'N', M, M, M, (1.0_dp,0.0_dp), dHdky, M, eigvec, M, (0.0_dp,0.0_dp), temp_y, M)

   ! Step 2: V = eigvec^† * temp
   call ZGEMM('C', 'N', M, M, M, (1.0_dp,0.0_dp), eigvec, M, temp_x, M, (0.0_dp,0.0_dp), Vx, M)
   call ZGEMM('C', 'N', M, M, M, (1.0_dp,0.0_dp), eigvec, M, temp_y, M, (0.0_dp,0.0_dp), Vy, M)

   ! Set diagonal elements to zero
   do i = 1, M
      Vx(i,i) = (0.0_dp, 0.0_dp)
      Vy(i,i) = (0.0_dp, 0.0_dp)
   end do

   ! Calculate Berry curvature for each specified band
   ! Ω_z(n) = -2 Im Σ_{m≠n} [V^x_nm × (V^y_nm)*] / (E_n - E_m)²
   berry_curv_bands = 0.0_dp

   ! Make energy threshold configurable for testing
   call MIO_InputParameter('Diag.ChernEnergyThreshold', energy_threshold, 1.0e-4_dp)

   ! Add proper energy broadening parameter in eV
   call MIO_InputParameter('Diag.ChernEtaEV', eta_eV, 1.0e-6_dp)

   ! Make regularization parameter configurable for testing
   call MIO_InputParameter('Diag.ChernEps', eps, 1.0e-10_dp)
   if (tapwDebug) call MIO_Print('DEBUG: Using eps = '//trim(num2str(eps,12))//' for Berry curvature calculation','diag')

   ! CRITICAL FIX: Use actual array dimensions for loop bounds, not M variable
   ! M might be corrupted or inconsistent - always use actual array sizes
   M_safe = min(size(Vx,1), size(Vx,2), size(Vy,1), size(Vy,2), size(eigval))

   ! Additional safety check before Berry curvature calculation
   if (kpoint_index == 1) then
      call MIO_Print('Debug: About to calculate Berry curvature','diag')
      call MIO_Print('  M variable value: '//trim(num2str(M)),'diag')
      call MIO_Print('  M_safe (actual array bounds): '//trim(num2str(M_safe)),'diag')
      call MIO_Print('  Vx matrix shape: '//trim(num2str(size(Vx,1)))//'x'//trim(num2str(size(Vx,2))),'diag')
      call MIO_Print('  Vy matrix shape: '//trim(num2str(size(Vy,1)))//'x'//trim(num2str(size(Vy,2))),'diag')
      call MIO_Print('  eigval array size: '//trim(num2str(size(eigval))),'diag')
      call MIO_Print('  First band index to access: '//trim(num2str(band_indices(1))),'diag')
      call MIO_Print('  Last band index to access: '//trim(num2str(band_indices(size(band_indices)))),'diag')

      if (M /= M_safe) then
         call MIO_Print('CRITICAL WARNING: M variable ('//trim(num2str(M))//') differs from actual array bounds ('//trim(num2str(M_safe))//')','diag')
         call MIO_Print('This indicates M was corrupted after array allocation - using M_safe for loops','diag')
      end if
   end if

   do i = 1, size(band_indices)  ! For each selected band
      n = band_indices(i)  ! Current band index
      if (n <= M_safe .and. n >= 1) then  ! Use M_safe instead of M
         if (kpoint_index == 1 .and. i <= 3) then  ! Debug first few bands
            call MIO_Print('Debug: Processing band '//trim(num2str(i))//' (index '//trim(num2str(n))//')','diag')
         end if

         berry_sum = (0.0_dp, 0.0_dp)
         do m = 1, M_safe  ! Use M_safe instead of M
            if (n /= m) then
               ! Skip nearly degenerate bands to avoid numerical instability
               if (abs(eigval(n) - eigval(m)) < energy_threshold) then
                  cycle  ! Skip this pair - too close in energy
               end if

               ! Additional bounds checking to prevent array access errors
               if (n <= size(Vx,1) .and. m <= size(Vx,2) .and. m <= size(Vy,1) .and. n <= size(Vy,2) .and. &
                   n <= size(eigval) .and. m <= size(eigval)) then
                  ! CORRECTED: Use Vy(m,n) to match Python's (Vynm.T) - transpose operation (no conjg!)
                  ! Python: (Vxnm*(Vynm.T)) means Vx(n,m) * Vy(m,n) - transpose, not conjugate transpose
                  berry_sum = berry_sum + (Vx(n,m) * Vy(m,n)) / ((eigval(n) - eigval(m))**2 + eps)
               else
                  call MIO_Print('WARNING: Array bounds exceeded in Berry calculation','diag')
                  call MIO_Print('  n='//trim(num2str(n))//', m='//trim(num2str(m))//', M='//trim(num2str(M)),'diag')
                  call MIO_Print('  Vx size: '//trim(num2str(size(Vx,1)))//'x'//trim(num2str(size(Vx,2))),'diag')
                  call MIO_Print('  Vy size: '//trim(num2str(size(Vy,1)))//'x'//trim(num2str(size(Vy,2))),'diag')
                  call MIO_Print('  eigval size: '//trim(num2str(size(eigval))),'diag')
               end if
            end if
         end do
         berry_curv_bands(i) = -2.0_dp * aimag(berry_sum)

         ! DEBUG: Compare velocity matrix elements for first few bands
         if (kpoint_index == 1 .and. i <= 3) then
            if (tapwDebug) call MIO_Print('DEBUG Option B: Velocity matrix elements for band '//trim(num2str(n))//':','diag')
            call MIO_Print('  |Vx('//trim(num2str(n))//','//trim(num2str(n+1))//')| = '//trim(num2str(abs(Vx(n,n+1)),8)),'diag')
            call MIO_Print('  |Vy('//trim(num2str(n))//','//trim(num2str(n+1))//')| = '//trim(num2str(abs(Vy(n,n+1)),8)),'diag')
            call MIO_Print('  Berry curvature = '//trim(num2str(berry_curv_bands(i),8)),'diag')
         end if

         ! Note: Large Berry curvature values are physically expected near band crossings
         ! No need to debug or cap these values - they're part of the physics

         if (kpoint_index == 1 .and. i <= 3) then  ! Debug first few bands
            call MIO_Print('  Berry curvature: '//trim(num2str(berry_curv_bands(i),6)),'diag')
            ! DEBUG: Check velocity matrix elements for first k-point
            if (i == 1) then
               if (tapwDebug) call MIO_Print('  DEBUG: Velocity matrix elements for band '//trim(num2str(n))//':','diag')
               call MIO_Print('    |Vx(n,n+1)| = '//trim(num2str(abs(Vx(n,min(n+1,M_safe))),8)),'diag')
               call MIO_Print('    |Vy(n,n+1)| = '//trim(num2str(abs(Vy(n,min(n+1,M_safe))),8)),'diag')
               call MIO_Print('    Energy gap = '//trim(num2str(abs(eigval(n) - eigval(min(n+1,M_safe))),8)),'diag')
            end if
         end if
      else
         call MIO_Print('WARNING: Band index '//trim(num2str(n))//' out of range [1,'//trim(num2str(M_safe))//']','diag')
         berry_curv_bands(i) = 0.0_dp  ! Band index out of range
      end if
   end do

   ! Clean up
   deallocate(H_k, eigvec, eigval, dHdkx, dHdky, Vx, Vy, temp_x, temp_y)
   if (allocated(delX)) deallocate(delX, delY)  ! Only deallocate if allocated (fallback method)
   if (allocated(dH_TB_dkx)) deallocate(dH_TB_dkx, dH_TB_dky)  ! Option B cleanup

end subroutine CalculateBerryAtKpointFromStored

subroutine GetTAPWHamiltonian(kpt, ucell, H0, maxNeigh, hopp, NList, Nneigh, neighCell, H_proj, eigval, eigvec, M, delX, delY)
   ! Get TAPW projected Hamiltonian, eigenvalues, eigenvectors, and position differences at a k-point
   ! FIXED: Now uses actual TAPW matrix size from M_tapw instead of hardcoded value
   use atoms, only : nAt, Species, layerIndex
   implicit none

   ! Input/Output parameters
   real(dp), intent(in) :: kpt(3)
   real(dp), intent(in) :: ucell(3,3), H0(nAt)
   integer, intent(in) :: maxNeigh
   complex(dp), intent(in) :: hopp(maxNeigh,nAt)
   integer, intent(in) :: NList(maxNeigh,nAt), Nneigh(nAt), neighCell(3,maxNeigh,nAt)
   complex(dp), allocatable, intent(out) :: H_proj(:,:), eigvec(:,:)
   real(dp), allocatable, intent(out) :: eigval(:)
   integer, intent(out) :: M
   real(dp), allocatable, intent(out) :: delX(:,:), delY(:,:)  ! Position differences for analytical derivatives

   ! Local variables
   complex(dp), allocatable :: ZWorkLoc(:)
   real(dp), allocatable :: DWorkLoc(:), ELoc_temp(:), H0_temp(:)
   integer :: info, lwork, i, j, N

   N = nAt

   ! Use actual TAPW matrix size from previous calculation
   if (M_tapw <= 0) then
      call MIO_Print('ERROR: M_tapw not initialized. TAPW calculation must be performed first.','diag')
      call MIO_Print('Current M_tapw = '//trim(num2str(M_tapw)),'diag')
      error stop 1
   end if

   ! CRITICAL FIX: For Berry curvature calculation, use stored_M instead of M_tapw
   ! M_tapw might be different from what was used to allocate stored arrays
   if (.not. allocated(stored_hamiltonians) .or. .not. allocated(stored_eigenvectors)) then
      call MIO_Print('ERROR: No stored TAPW data available','diag')
      call MIO_Print('TAPW bands calculation must be performed before Chern calculation','diag')
      error stop 1
   end if

   ! Use the stored M value to ensure consistency with stored arrays
   M = stored_M
   call MIO_Print('Using M = '//trim(num2str(M))//' from stored TAPW data (stored_M)','diag')

   ! Warn if M_tapw differs from stored_M
   if (M_tapw /= stored_M) then
      call MIO_Print('WARNING: M_tapw ('//trim(num2str(M_tapw))//') differs from stored_M ('//trim(num2str(stored_M))//')','diag')
      call MIO_Print('Using stored_M to maintain consistency with stored arrays','diag')
   end if

   ! Allocate output arrays with stored dimensions
   allocate(H_proj(M,M))
   allocate(eigval(M))
   allocate(eigvec(M,M))
   allocate(delX(M,M))
   allocate(delY(M,M))
   allocate(ELoc_temp(N))
   allocate(H0_temp(N))

   call MIO_Print('Using stored TAPW Hamiltonian and eigenvectors from bands calculation','diag')

   ! Use the first k-point's data as reference (this should be improved for better k-point handling)
   if (stored_nk < 1) then
      call MIO_Print('ERROR: No k-points stored in TAPW data','diag')
      error stop 1
   end if

   ! Extract Hamiltonian and eigenvectors from stored data (use first k-point as reference, spin 1)
   ! For non-SOC calculations, spin_index=1 is the only spin channel
   H_proj = stored_hamiltonians(:,:,1,1)
   eigvec = stored_eigenvectors(:,:,1,1)

   ! Extract eigenvalues from the first k-point's diagonalization (already computed in bands calculation)
   ! The eigenvectors are already stored, so we just need to extract the eigenvalues
   lwork = 2*M
   allocate(ZWorkLoc(lwork))
   allocate(DWorkLoc(3*M-2))

   ! Make a copy of the Hamiltonian for diagonalization to get eigenvalues
   ! (eigenvectors are already available from stored data)
   call ZHEEV('N', 'U', M, H_proj, M, eigval, ZWorkLoc, lwork, DWorkLoc, info)

   if (info /= 0) then
      call MIO_Print('ERROR in GetTAPWHamiltonian: ZHEEV failed with info = '//trim(num2str(info)),'diag')
      error stop 1
   end if

   call MIO_Print('Successfully extracted TAPW eigenvalues and eigenvectors','diag')

   ! Calculate position difference matrices for analytical derivatives using saved G-vectors
   call MIO_Print('Computing analytical derivatives for TAPW projected space','diag')

   if (allocated(saved_Gx) .and. allocated(saved_Gy) .and. saved_NG > 0 .and. saved_Nlabel > 0) then
      call MIO_Print('Building delX/delY matrices from '//trim(num2str(saved_NG))//' G-vectors','diag')
      call build_TAPW_position_differences(delX, delY, M, saved_Gx, saved_Gy, saved_NG, saved_Nlabel)
   else
      call MIO_Print('WARNING: No saved G-vectors found - using zero matrices','diag')
      delX = 0.0_dp
      delY = 0.0_dp
   end if

   ! Cleanup
   deallocate(ZWorkLoc, DWorkLoc, ELoc_temp, H0_temp)

end subroutine GetTAPWHamiltonian

subroutine build_TAPW_position_differences(delX, delY, M, Gx, Gy, NG, Nlabel)
   ! Build position difference matrices for TAPW projected space
   ! OPTIMIZED VERSION: Combines caching, block-structured computation, and pre-computed indices
   ! Expected speedup: 100x+ after first call due to caching, 5-10x on first call due to block structure
   implicit none

   ! Input/Output parameters
   real(dp), intent(out) :: delX(M,M), delY(M,M)
   integer, intent(in) :: M, NG, Nlabel
   real(dp), intent(in) :: Gx(NG), Gy(NG)

   ! Local variables
   integer :: i, j
   logical :: gvectors_changed

   ! Check if cache is valid and G-vectors haven't changed
   gvectors_changed = .false.
   if (cache_valid .and. cached_M == M .and. cached_NG == NG .and. cached_Nlabel == Nlabel) then
      if (allocated(cached_Gx) .and. allocated(cached_Gy)) then
         if (any(abs(cached_Gx - Gx) > 1e-12_dp) .or. any(abs(cached_Gy - Gy) > 1e-12_dp)) then
            gvectors_changed = .true.
         end if
      else
         gvectors_changed = .true.
      end if
   else
      gvectors_changed = .true.
   end if

   ! TEMPORARY FIX: Disable cache to avoid segfault
   ! TODO: Investigate why cache reuse causes segmentation fault
   if (.false.) then  ! Force recomputation every time
      ! Use cached result if available and valid
      if (cache_valid .and. cached_M == M .and. cached_NG == NG .and. cached_Nlabel == Nlabel .and. .not. gvectors_changed) then
         ! Additional safety check: verify cached matrix dimensions
         if (allocated(cached_delX) .and. allocated(cached_delY)) then
            if (size(cached_delX,1) == M .and. size(cached_delX,2) == M .and. &
                size(cached_delY,1) == M .and. size(cached_delY,2) == M) then
               ! Safe to use cached matrices
               delX = cached_delX
               delY = cached_delY
               call MIO_Print('Using cached TAPW position difference matrices','diag')
               return
            else
               call MIO_Print('WARNING: Cached matrix dimensions mismatch - recomputing','diag')
               call MIO_Print('  Expected: '//trim(num2str(M))//'x'//trim(num2str(M)),'diag')
               if (allocated(cached_delX)) then
                  call MIO_Print('  Cached delX: '//trim(num2str(size(cached_delX,1)))//'x'//trim(num2str(size(cached_delX,2))),'diag')
               end if
               if (allocated(cached_delY)) then
                  call MIO_Print('  Cached delY: '//trim(num2str(size(cached_delY,1)))//'x'//trim(num2str(size(cached_delY,2))),'diag')
               end if
            end if
         else
            call MIO_Print('WARNING: Cached matrices not allocated - recomputing','diag')
         end if
      end if
   end if

   call MIO_Print('TEMPORARY: Cache disabled - recomputing position difference matrices','diag')

   call MIO_Print('Computing TAPW position difference matrices (optimized block structure)','diag')
   call MIO_Print('  Matrix size: '//trim(num2str(M))//' x '//trim(num2str(M)),'diag')
   call MIO_Print('  Block structure: '//trim(num2str(Nlabel))//' labels × '//trim(num2str(NG))//' G-vectors','diag')

   ! Compute using optimized block-structured approach
   call compute_position_differences_blocked(delX, delY, M, Gx, Gy, NG, Nlabel)

   ! Update cache
   if (allocated(cached_delX)) deallocate(cached_delX, cached_delY)
   if (allocated(cached_Gx)) deallocate(cached_Gx, cached_Gy)

   allocate(cached_delX(M,M), cached_delY(M,M))
   allocate(cached_Gx(NG), cached_Gy(NG))

   cached_delX = delX
   cached_delY = delY
   cached_Gx = Gx
   cached_Gy = Gy
   cached_M = M
   cached_NG = NG
   cached_Nlabel = Nlabel
   cache_valid = .true.

   call MIO_Print('TAPW position difference matrices cached for future use','diag')

end subroutine build_TAPW_position_differences

subroutine compute_position_differences_blocked(delX, delY, M, Gx, Gy, NG, Nlabel)
   ! Optimized block-structured computation of position difference matrices
   ! This eliminates all modulo/division operations and improves cache locality
   implicit none

   ! Input/Output parameters
   real(dp), intent(out) :: delX(M,M), delY(M,M)
   integer, intent(in) :: M, NG, Nlabel
   real(dp), intent(in) :: Gx(NG), Gy(NG)

   ! Local variables
   integer :: label, i_start, i_end, i, j, i_G, j_G

   ! Initialize matrices to zero
   delX = 0.0_dp
   delY = 0.0_dp

   ! TAPW projected space structure: M = NG × Nlabel
   ! Process each label block separately for optimal performance
   do label = 1, Nlabel
      i_start = (label-1) * NG + 1
      i_end = label * NG

      ! Only compute within this label block (inter-block elements remain zero)
      do i = i_start, i_end
         i_G = i - i_start + 1  ! G-vector index within block: 1..NG

         do j = i_start, i_end
            j_G = j - i_start + 1  ! G-vector index within block: 1..NG
            delX(i,j) = Gx(i_G) - Gx(j_G)
            delY(i,j) = Gy(i_G) - Gy(j_G)
         end do
      end do
   end do

end subroutine compute_position_differences_blocked

subroutine compute_TB_hamiltonian_derivatives(dH_TB_dkx, dH_TB_dky, N, KLoc, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
   ! Compute derivatives of the original tight-binding Hamiltonian
   ! dH^TB/dkx and dH^TB/dky using the same gauge as in transform_dense_hamiltonian_tapw
   use constants, only : cmplx_i, cmplx_0
   use neigh, only : NeighD
   implicit none

   ! Input parameters
   integer, intent(in) :: N, maxNeigh
   real(dp), intent(in) :: KLoc(3), H0(N)
   complex(dp), intent(in) :: hopp(maxNeigh,N)
   integer, intent(in) :: NList(maxNeigh,N), Nneigh(N), neighCell(3,maxNeigh,N)
   complex(dp), intent(out) :: dH_TB_dkx(N,N), dH_TB_dky(N,N)

   ! Local variables
   integer :: i, j, in
   real(dp) :: R(3)
   complex(dp) :: phase_factor, derivative_factor_x, derivative_factor_y

   call MIO_Print('Computing TB Hamiltonian derivatives dH^TB/dk','diag')

   ! DEBUG: Check units and magnitudes
   if (N >= 1 .and. Nneigh(1) >= 1) then
      if (tapwDebug) call MIO_Print('DEBUG: Units check for TB derivatives:','diag')
      call MIO_Print('  KLoc = ['//trim(num2str(KLoc(1),6))//','//trim(num2str(KLoc(2),6))//','//trim(num2str(KLoc(3),6))//'] (1/Angstrom)','diag')
      if (NList(1,1) > 0 .and. NList(1,1) <= N) then
         call MIO_Print('  First NeighD = ['//trim(num2str(NeighD(1,1,1),6))//','//trim(num2str(NeighD(2,1,1),6))//'] (Angstrom)','diag')
         call MIO_Print('  R used in derivatives = ['//trim(num2str(-NeighD(1,1,1),6))//','//trim(num2str(-NeighD(2,1,1),6))//'] (Angstrom)','diag')
         call MIO_Print('  k·R = '//trim(num2str(dot_product(KLoc(1:2), -NeighD(1:2,1,1)),6))//' (dimensionless)','diag')
      end if
   end if

   ! Initialize derivative matrices
   dH_TB_dkx = cmplx_0
   dH_TB_dky = cmplx_0

   ! Diagonal terms have zero derivatives (on-site energies are k-independent)
   ! No contribution from H0(i) terms

   ! Off-diagonal hopping terms derivatives
   ! Original: H_ij = -t_ij * exp(-i k·R_ij)
   ! Derivative: dH_ij/dkx = -t_ij * (-i R_ij^x) * exp(-i k·R_ij) = i * R_ij^x * H_ij
   do i = 1, N
      do j = 1, Nneigh(i)
         in = NList(j, i)
         if (in > 0 .and. in <= N) then
            ! Use correct position difference for derivative calculation
            ! NeighD(:,j,i) contains vector from atom i to atom j: (τ_j - τ_i) + T_ij
            ! For derivatives, use negative sign to match dense Hamiltonian convention
            R(1:3) = 0.0_dp
            R(1:2) = -NeighD(1:2, j, i)

            ! Compute the phase factor (same as in original Hamiltonian)
            phase_factor = exp(-cmplx_i * dot_product(KLoc, R))

            ! Derivatives: d/dkx [exp(-i k·R)] = -i R_x exp(-i k·R)
            derivative_factor_x = -cmplx_i * R(1) * phase_factor
            derivative_factor_y = -cmplx_i * R(2) * phase_factor

            ! Apply to hopping term: dH/dk = -t_ij * d/dk[exp(-i k·R)]
            dH_TB_dkx(i, in) = dH_TB_dkx(i, in) - hopp(j, i) * derivative_factor_x
            dH_TB_dky(i, in) = dH_TB_dky(i, in) - hopp(j, i) * derivative_factor_y
         end if
      end do
   end do

   call MIO_Print('TB Hamiltonian derivatives computed','diag')

end subroutine compute_TB_hamiltonian_derivatives

!> @brief Compute derivatives of block Hamiltonian for SOC
!! @details Computes dH^Block/dkx and dH^Block/dky for 2N×2N block Hamiltonian
!! @param[out]    dH_Block_dkx    Block Hamiltonian derivative w.r.t. kx (2N×2N)
!! @param[out]    dH_Block_dky    Block Hamiltonian derivative w.r.t. ky (2N×2N)
!! @param[in]     N               Number of atoms
!! @param[in]     KLoc            k-point vector
!! @param[in]     H0              Onsite energies
!! @param[in]     maxNeigh        Maximum number of neighbors
!! @param[in]     hopp            Hopping parameters
!! @param[in]     NList           Neighbor list
!! @param[in]     Nneigh          Number of neighbors per atom
!! @param[in]     neighCell       Neighbor cell indices
subroutine compute_block_hamiltonian_derivatives(dH_Block_dkx, dH_Block_dky, N, KLoc, cell, H0, maxNeigh, hopp, NList, Nneigh, neighCell)
   ! Compute derivatives of the block Hamiltonian for SOC
   ! dH^Block/dkx and dH^Block/dky using the same gauge as in BuildBlockHamiltonianOnly
   use constants, only : cmplx_i, cmplx_0
   use neigh, only : NeighD
   use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use atoms, only : Species, layerIndex
   use ham, only : RashbaSOCterm, SOCEnabledForLayer, lambdaR
   implicit none

   ! Input parameters
   integer, intent(in) :: N, maxNeigh
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxNeigh,N)
   integer, intent(in) :: NList(maxNeigh,N), Nneigh(N), neighCell(3,maxNeigh,N)
   complex(dp), intent(out) :: dH_Block_dkx(2*N,2*N), dH_Block_dky(2*N,2*N)

   ! Local variables
   integer :: i, j, in
   real(dp) :: R(3)
   complex(dp) :: phase_factor, derivative_factor_x, derivative_factor_y
   complex(dp) :: RashbaHopp
   real(dp) :: dx, dy, dnorm
   ! Note: lambdaR is already imported from ham module
   complex(dp) :: RashbaDerivative_x, RashbaDerivative_y

   call MIO_Print('Computing block Hamiltonian derivatives dH^Block/dk for SOC','diag')

   ! Initialize derivative matrices
   dH_Block_dkx = cmplx_0
   dH_Block_dky = cmplx_0

   ! Diagonal terms have zero derivatives (on-site energies are k-independent)
   ! No contribution from H0(i) terms

   ! Off-diagonal hopping terms derivatives for both spin blocks
   do i = 1, N
      do j = 1, Nneigh(i)
         in = NList(j, i)
         if (in > 0 .and. in <= N) then
            ! Use correct position difference for derivative calculation
            R(1:3) = 0.0_dp
            R(1:2) = -NeighD(1:2, j, i)

            ! Compute the phase factor (same as in original Hamiltonian)
            phase_factor = exp(-cmplx_i * dot_product(KLoc, R))

            ! Derivatives: d/dkx [exp(-i k·R)] = -i R_x exp(-i k·R)
            derivative_factor_x = -cmplx_i * R(1) * phase_factor
            derivative_factor_y = -cmplx_i * R(2) * phase_factor

            ! Apply to hopping term: dH/dk = -t_ij * d/dk[exp(-i k·R)]
            ! Regular hopping (preserved in both spin blocks)
            ! BuildBlockHamiltonianOnly sets HBlock(in, i), so derivatives should match: dH(in, i)/dk
            dH_Block_dkx(in, i) = dH_Block_dkx(in, i) - hopp(j, i) * derivative_factor_x
            dH_Block_dky(in, i) = dH_Block_dky(in, i) - hopp(j, i) * derivative_factor_y

            dH_Block_dkx(in + N, i + N) = dH_Block_dkx(in + N, i + N) - hopp(j, i) * derivative_factor_x
            dH_Block_dky(in + N, i + N) = dH_Block_dky(in + N, i + N) - hopp(j, i) * derivative_factor_y

            ! Apply Rashba spin-flip terms derivatives if enabled
            if (RashbaSOCterm) then
               ! Check if SOC should be applied to this layer (matches ApplySpinFlipSOC)
               if (.not. SOCEnabledForLayer(layerIndex(i))) then
                  cycle  ! Skip Rashba derivatives for this layer
               end if

               ! Rashba term from ApplySpinFlipSOC:
               ! where dx, dy are normalized direction cosines

               ! Normalize direction cosines (same as in ApplySpinFlipSOC)
               ! ApplySpinFlipSOC uses: dx = neighD(1, j, i), dy = neighD(2, j, i)
               dx = NeighD(1, j, i)  ! x-component (matches ApplySpinFlipSOC line 16871)
               dy = NeighD(2, j, i)  ! y-component (matches ApplySpinFlipSOC line 16872)
               dnorm = sqrt(dx*dx + dy*dy)

               if (dnorm > 1.0e-12_dp) then
                  dx = dx / dnorm
                  dy = dy / dnorm
               end if

               ! Rashba coupling strength (use actual parameter from ham module)
               ! Note: lambdaR should be imported from ham module at subroutine level

               ! Rashba term: RashbaHopp = lambdaR * (dx + i*dy) * exp(-i k·R)
               ! Derivative: dRashbaHopp/dkx = lambdaR * (dx + i*dy) * (-i R_x) * exp(-i k·R)
               RashbaDerivative_x = lambdaR * (dx + cmplx_i * dy) * derivative_factor_x
               RashbaDerivative_y = lambdaR * (dx + cmplx_i * dy) * derivative_factor_y

               ! Apply Rashba derivatives to spin-flip terms
               ! H(i, in+N) = RashbaHopp, so dH(i, in+N)/dk = RashbaDerivative
               ! H(in+N, i) = conjg(RashbaHopp), so dH(in+N, i)/dk = conjg(RashbaDerivative)
               dH_Block_dkx(i, in + N) = dH_Block_dkx(i, in + N) + RashbaDerivative_x
               dH_Block_dky(i, in + N) = dH_Block_dky(i, in + N) + RashbaDerivative_y

               dH_Block_dkx(in + N, i) = dH_Block_dkx(in + N, i) + conjg(RashbaDerivative_x)
               dH_Block_dky(in + N, i) = dH_Block_dky(in + N, i) + conjg(RashbaDerivative_y)
            end if
         end if
      end do
   end do

   ! Edge hopping derivatives if applicable
   ! Edge hopping uses lattice vectors: R = matmul(cell, NedgeCell) (matches BuildBlockHamiltonianOnly line 17145)
   if (edgeHopp) then
      do i = 1, nQ
         do j = 1, nEdgeN(i)
            in = NeI(j, i)
            ! Use same convention as BuildBlockHamiltonianOnly: R = matmul(cell, NedgeCell)
            R = matmul(cell, NedgeCell(:, j, i))

            phase_factor = exp(-cmplx_i * dot_product(KLoc, R))
            derivative_factor_x = -cmplx_i * R(1) * phase_factor
            derivative_factor_y = -cmplx_i * R(2) * phase_factor

            ! Apply to both spin blocks (matches BuildBlockHamiltonianOnly line 17147-17148)
            dH_Block_dkx(in, edgeIndx(i)) = dH_Block_dkx(in, edgeIndx(i)) + edgeH(j, i) * derivative_factor_x
            dH_Block_dky(in, edgeIndx(i)) = dH_Block_dky(in, edgeIndx(i)) + edgeH(j, i) * derivative_factor_y

            dH_Block_dkx(in + N, edgeIndx(i) + N) = dH_Block_dkx(in + N, edgeIndx(i) + N) + edgeH(j, i) * derivative_factor_x
            dH_Block_dky(in + N, edgeIndx(i) + N) = dH_Block_dky(in + N, edgeIndx(i) + N) + edgeH(j, i) * derivative_factor_y
         end do
      end do
   end if

   call MIO_Print('Block Hamiltonian derivatives computed for SOC','diag')

end subroutine compute_block_hamiltonian_derivatives

subroutine project_TB_derivatives_to_TAPW(dH_TB_dkx, dH_TB_dky, X, N, M, dH_TAPW_dkx, dH_TAPW_dky)
   ! Project TB Hamiltonian derivatives to TAPW space: dH^TAPW/dk = X† * (dH^TB/dk) * X
   use constants, only : cmplx_0, cmplx_1
   implicit none

   ! Input parameters
   integer, intent(in) :: N, M
   complex(dp), intent(in) :: dH_TB_dkx(N,N), dH_TB_dky(N,N)
   complex(dp), intent(in) :: X(N,M)
   complex(dp), intent(out) :: dH_TAPW_dkx(M,M), dH_TAPW_dky(M,M)

   ! Local variables
   complex(dp), allocatable :: temp_x(:,:), temp_y(:,:)

   call MIO_Print('Projecting TB derivatives to TAPW space','diag')

   ! Allocate temporary matrices
   allocate(temp_x(N,M), temp_y(N,M))

   ! Project dH^TB/dkx to TAPW space: dH^TAPW/dkx = X† * dH^TB/dkx * X
   ! Step 1: temp_x = dH^TB/dkx * X
   call zgemm('N', 'N', N, M, N, cmplx_1, dH_TB_dkx, N, X, N, cmplx_0, temp_x, N)

   ! Step 2: dH^TAPW/dkx = X† * temp_x
   call zgemm('C', 'N', M, M, N, cmplx_1, X, N, temp_x, N, cmplx_0, dH_TAPW_dkx, M)

   ! Project dH^TB/dky to TAPW space: dH^TAPW/dky = X† * dH^TB/dky * X
   ! Step 1: temp_y = dH^TB/dky * X
   call zgemm('N', 'N', N, M, N, cmplx_1, dH_TB_dky, N, X, N, cmplx_0, temp_y, N)

   ! Step 2: dH^TAPW/dky = X† * temp_y
   call zgemm('C', 'N', M, M, N, cmplx_1, X, N, temp_y, N, cmplx_0, dH_TAPW_dky, M)

   ! Clean up
   deallocate(temp_x, temp_y)

   call MIO_Print('TB derivatives successfully projected to TAPW space','diag')

end subroutine project_TB_derivatives_to_TAPW

!> @brief Project block Hamiltonian derivatives to TAPW space for SOC
!! @details Projects dH^Block/dk derivatives to TAPW space: dH^TAPW/dk = X† * (dH^Block/dk) * X
!! @param[in]     dH_Block_dkx    Block Hamiltonian derivative w.r.t. kx (2N×2N)
!! @param[in]     dH_Block_dky    Block Hamiltonian derivative w.r.t. ky (2N×2N)
!! @param[in]     X               TAPW projection matrix (2N×M)
!! @param[in]     N               Number of atoms (2N for SOC)
!! @param[in]     M               TAPW matrix size
!! @param[out]    dH_TAPW_dkx     Projected derivative w.r.t. kx (M×M)
!! @param[out]    dH_TAPW_dky     Projected derivative w.r.t. ky (M×M)
subroutine project_block_derivatives_to_TAPW(dH_Block_dkx, dH_Block_dky, X, N, M, dH_TAPW_dkx, dH_TAPW_dky)
   ! Project block Hamiltonian derivatives to TAPW space: dH^TAPW/dk = X† * (dH^Block/dk) * X
   use constants, only : cmplx_0, cmplx_1
   implicit none

   ! Input parameters
   integer, intent(in) :: N, M
   complex(dp), intent(in) :: dH_Block_dkx(N,N), dH_Block_dky(N,N)
   complex(dp), intent(in) :: X(N,M)
   complex(dp), intent(out) :: dH_TAPW_dkx(M,M), dH_TAPW_dky(M,M)

   ! Local variables
   complex(dp), allocatable :: temp(:,:)

   ! Allocate temporary matrix
   allocate(temp(M, N))

   ! Project derivatives: dH^TAPW/dk = X† * (dH^Block/dk) * X
   ! First: temp = X† * dH^Block/dk
   call ZGEMM('C', 'N', M, N, N, cmplx_1, X, N, dH_Block_dkx, N, cmplx_0, temp, M)
   ! Second: dH^TAPW/dkx = temp * X
   call ZGEMM('N', 'N', M, M, N, cmplx_1, temp, M, X, N, cmplx_0, dH_TAPW_dkx, M)

   ! Same for ky derivative
   call ZGEMM('C', 'N', M, N, N, cmplx_1, X, N, dH_Block_dky, N, cmplx_0, temp, M)
   call ZGEMM('N', 'N', M, M, N, cmplx_1, temp, M, X, N, cmplx_0, dH_TAPW_dky, M)

   deallocate(temp)

end subroutine project_block_derivatives_to_TAPW

subroutine transform_TAPW_eigenvectors_to_TB(eigvec_TAPW, X, N, M, eigvec_TB)
   ! Transform TAPW eigenvectors back to TB basis: ψ_TB = X * ψ_TAPW
   ! Input: eigvec_TAPW(M,M) - TAPW eigenvectors (columns are eigenvectors)
   ! Input: X(N,M) - transformation matrix from TB to TAPW
   ! Output: eigvec_TB(N,M) - TB eigenvectors (only M vectors in N-dimensional space)
   use constants, only : cmplx_0, cmplx_1
   implicit none

   ! Input parameters
   integer, intent(in) :: N, M
   complex(dp), intent(in) :: eigvec_TAPW(M,M)
   complex(dp), intent(in) :: X(N,M)
   complex(dp), intent(out) :: eigvec_TB(N,M)

   call MIO_Print('Transforming TAPW eigenvectors back to TB basis','diag')
   call MIO_Print('  Input: '//trim(num2str(M))//' TAPW eigenvectors in '//trim(num2str(M))//'-dimensional space','diag')
   call MIO_Print('  Output: '//trim(num2str(M))//' TB eigenvectors in '//trim(num2str(N))//'-dimensional space','diag')

   ! Transform eigenvectors: eigvec_TB = X * eigvec_TAPW
   call zgemm('N', 'N', N, M, M, cmplx_1, X, N, eigvec_TAPW, M, cmplx_0, eigvec_TB, N)

   call MIO_Print('TAPW eigenvectors successfully transformed to TB basis','diag')

end subroutine transform_TAPW_eigenvectors_to_TB

subroutine CalculateBerryAtKpointFromStored_OptionA(kpoint_index, spin_index, band_indices, eigval_from_E, berry_curv_bands)
   ! OPTION A: Transform TAPW eigenvectors to TB basis, then calculate Berry curvature in TB space
   ! This follows your colleague's approach for comparison with Option B
   ! Updated to include spin_index for separate spin channel calculations
   use constants, only : cmplx_i, cmplx_0, cmplx_1, pi
   use atoms, only : nAt
   use ham, only : H0, hopp
   use neigh, only : NList, Nneigh, neighCell, maxNeigh
   implicit none

   ! Input parameters
   integer, intent(in) :: kpoint_index, spin_index
   integer, intent(in) :: band_indices(:)
   real(dp), intent(in) :: eigval_from_E(:)
   real(dp), intent(out) :: berry_curv_bands(:)

   ! Local variables
   integer :: M_local, N_local, i, j, n, m, num_bands
   complex(dp), allocatable :: H_k(:,:), eigvec_TAPW(:,:), eigvec_TB(:,:)
   real(dp), allocatable :: eigval(:)
   complex(dp), allocatable :: dH_TB_dkx(:,:), dH_TB_dky(:,:)
   complex(dp), allocatable :: Vx_TB(:,:), Vy_TB(:,:)
   real(dp) :: current_kpt(3), eps, energy_diff, energy_threshold, eta_eV
   integer :: skipped_pairs, total_pairs, M_safe
   complex(dp) :: berry_sum

   call MIO_Print('=== OPTION A: Berry curvature in TB basis ===','diag')

   ! Check stored data availability (only need X matrix, TB parameters are module-level)
   if (.not. allocated(stored_X_matrix)) then
      call MIO_Print('ERROR: Stored X matrix not available for Option A','diag')
      error stop 1
   end if

   if (.not. allocated(stored_hamiltonians) .or. .not. allocated(stored_eigenvectors)) then
      call MIO_Print('ERROR: Stored TAPW data not available for Option A','diag')
      error stop 1
   end if

   ! Get dimensions (use module-level variables directly)
   N_local = nAt  ! Use module-level nAt instead of stored_N
   M_local = stored_M
   num_bands = size(band_indices)

   call MIO_Print('Option A dimensions: N='//trim(num2str(N_local))//', M='//trim(num2str(M_local))//', bands='//trim(num2str(num_bands)),'diag')

   ! Allocate arrays
   allocate(H_k(M_local,M_local), eigvec_TAPW(M_local,M_local), eigvec_TB(N_local,M_local))
   allocate(eigval(M_local))
   allocate(dH_TB_dkx(N_local,N_local), dH_TB_dky(N_local,N_local))
   allocate(Vx_TB(M_local,M_local), Vy_TB(M_local,M_local))

   ! Check spin_index bounds
   if (spin_index < 1 .or. spin_index > stored_nspin) then
      call MIO_Print('Error: spin_index '//trim(num2str(spin_index))//' out of range [1,'//trim(num2str(stored_nspin))//']','diag')
      berry_curv_bands = 0.0_dp
      return
   end if

   ! Get stored TAPW data for this k-point and spin channel
   H_k(1:M_local,1:M_local) = stored_hamiltonians(1:M_local,1:M_local,kpoint_index,spin_index)
   eigvec_TAPW(1:M_local,1:M_local) = stored_eigenvectors(1:M_local,1:M_local,kpoint_index,spin_index)

   ! CRITICAL FIX: Use TAPW eigenvalues instead of TB eigenvalues
   ! The eigenvectors are from TAPW diagonalization, so eigenvalues must be too
   if (allocated(stored_eigenvalues)) then
      if (size(stored_eigenvalues,1) >= M_local .and. size(stored_eigenvalues,2) >= kpoint_index .and. size(stored_eigenvalues,3) >= spin_index) then
         eigval(1:M_local) = stored_eigenvalues(1:M_local,kpoint_index,spin_index)
         call MIO_Print('Using stored TAPW eigenvalues for k-point '//trim(num2str(kpoint_index)),'diag')
         call MIO_Print('  TAPW eigenvalues range: ['//trim(num2str(minval(eigval(1:M_local)),6))//','//trim(num2str(maxval(eigval(1:M_local)),6))//']','diag')
      else
         call MIO_Print('ERROR: stored_eigenvalues array too small','diag')
         call MIO_Print('  Array size: '//trim(num2str(size(stored_eigenvalues,1)))//'x'//trim(num2str(size(stored_eigenvalues,2))),'diag')
         call MIO_Print('  Required: '//trim(num2str(M_local))//'x'//trim(num2str(kpoint_index)),'diag')
         error stop 1
      end if
   else
      call MIO_Print('ERROR: No stored TAPW eigenvalues available - using TB eigenvalues (WRONG!)','diag')
      call MIO_Print('This will cause eigenvector-eigenvalue misalignment!','diag')
      call MIO_Print('  TB eigenvalues range: ['//trim(num2str(minval(eigval_from_E),6))//','//trim(num2str(maxval(eigval_from_E),6))//']','diag')
      ! Fallback to TB eigenvalues (this is wrong but prevents crash)
   if (size(eigval_from_E) >= M_local) then
      eigval(1:M_local) = eigval_from_E(1:M_local)
   else
      eigval(1:size(eigval_from_E)) = eigval_from_E
      eigval(size(eigval_from_E)+1:M_local) = 0.0_dp
      end if
   end if

   ! VERIFICATION: Check eigenvalues for first k-point
   if (kpoint_index == 1) then
      call MIO_Print('=== EIGENVALUES VERIFICATION (k-point 1) ===','diag')
      call MIO_Print('Eigenvalues dimensions: '//trim(num2str(size(eigval))),'diag')
      call MIO_Print('First 5 eigenvalues: ['//trim(num2str(eigval(1),6))//','//trim(num2str(eigval(2),6))//','//trim(num2str(eigval(3),6))//','//trim(num2str(eigval(4),6))//','//trim(num2str(eigval(5),6))//']','diag')
      call MIO_Print('Eigenvalue range: ['//trim(num2str(minval(eigval),6))//','//trim(num2str(maxval(eigval),6))//']','diag')
   end if

   ! SANITY CHECK: Verify eigenvector-eigenvalue alignment (only if tapwDebug enabled)
   if (tapwDebug) then
      call MIO_Print('=== EIGENVECTOR-EIGENVALUE ALIGNMENT CHECK ===','diag')
      call verify_eigenvector_eigenvalue_alignment(H_k, eigvec_TAPW, eigval, M_local, kpoint_index)
   end if

   ! Get actual k-point coordinates
   current_kpt = stored_Kpts(:, kpoint_index)
   if (tapwDebug) then
   call MIO_Print('Using actual k-point coordinates: ['//trim(num2str(current_kpt(1),6))//','//&
                  trim(num2str(current_kpt(2),6))//','//trim(num2str(current_kpt(3),6))//']','diag')
   end if

   ! Step 1: Compute TB Hamiltonian derivatives (use module-level variables directly)
   call compute_TB_hamiltonian_derivatives(dH_TB_dkx, dH_TB_dky, N_local, current_kpt, &
                                         H0, maxNeigh, hopp, NList, Nneigh, neighCell)

   ! Step 1.5: Check hermiticity of TB Hamiltonian derivatives (only if tapwDebug enabled)
   if (tapwDebug) then
      call MIO_Print('Checking hermiticity of TB Hamiltonian derivatives','diag')
      call check_matrix_hermiticity(dH_TB_dkx, N_local, 'dH_TB_dkx')
      call check_matrix_hermiticity(dH_TB_dky, N_local, 'dH_TB_dky')
   end if

   ! Step 2: Transform TAPW eigenvectors back to TB basis
   call transform_TAPW_eigenvectors_to_TB(eigvec_TAPW, stored_X_matrix, N_local, M_local, eigvec_TB)

   ! Note: Eigenvectors are already properly normalized after transformation

   ! DEBUG: Check eigenvector norms after transformation and renormalization
   if (num_bands >= 1) then
      block
         real(dp) :: norm_tapw, norm_tb
         integer :: debug_band
         debug_band = band_indices(1)
         if (debug_band >= 1 .and. debug_band <= M_local) then
            norm_tapw = real(dot_product(eigvec_TAPW(:,debug_band), eigvec_TAPW(:,debug_band)))
            norm_tb = real(dot_product(eigvec_TB(:,debug_band), eigvec_TB(:,debug_band)))
            if (tapwDebug) call MIO_Print('DEBUG Option A: Band '//trim(num2str(debug_band))//' norms AFTER renormalization:','diag')
            call MIO_Print('  TAPW eigenvector norm: '//trim(num2str(sqrt(norm_tapw),8)),'diag')
            call MIO_Print('  TB eigenvector norm: '//trim(num2str(sqrt(norm_tb),8))//' (should be 1.0)','diag')
         end if
      end block
   end if

   ! Step 3: Calculate velocity matrix elements in TB space
   call MIO_Print('Computing velocity matrix elements in TB space','diag')

   ! Note: eigvec_TB is N×M, so we compute M×M velocity matrices
   ! This is computationally expensive: (N×M)† * (N×N) * (N×M) = M×M

   ! Use temporary array for intermediate result
   block
      complex(dp), allocatable :: temp_TB(:,:)
      allocate(temp_TB(N_local,M_local))

      ! Step 3a: temp = dH_TB_dkx * eigvec_TB
      call zgemm('N', 'N', N_local, M_local, N_local, cmplx_1, dH_TB_dkx, N_local, eigvec_TB, N_local, cmplx_0, temp_TB, N_local)

      ! Step 3b: Vx_TB = eigvec_TB† * temp
      call zgemm('C', 'N', M_local, M_local, N_local, cmplx_1, eigvec_TB, N_local, temp_TB, N_local, cmplx_0, Vx_TB, M_local)

      ! Repeat for y-direction
      call zgemm('N', 'N', N_local, M_local, N_local, cmplx_1, dH_TB_dky, N_local, eigvec_TB, N_local, cmplx_0, temp_TB, N_local)
      call zgemm('C', 'N', M_local, M_local, N_local, cmplx_1, eigvec_TB, N_local, temp_TB, N_local, cmplx_0, Vy_TB, M_local)

      deallocate(temp_TB)
   end block

   ! Step 4: Set diagonal elements to zero (same as Option B)
   do i = 1, M_local
      Vx_TB(i,i) = cmplx_0
      Vy_TB(i,i) = cmplx_0
   end do

   ! Step 5: Check hermiticity of velocity operators (only if tapwDebug enabled)
   if (tapwDebug) then
      call MIO_Print('Checking hermiticity of velocity operators Vx_TB and Vy_TB','diag')
      call check_matrix_hermiticity(Vx_TB, M_local, 'Vx_TB')
      call check_matrix_hermiticity(Vy_TB, M_local, 'Vy_TB')
   end if

   ! DEBUG: Check velocity matrix element magnitudes
   if (num_bands >= 2) then
      block
         integer :: n1, n2
         n1 = band_indices(1)
         n2 = band_indices(2)
         if (n1 >= 1 .and. n1 <= M_local .and. n2 >= 1 .and. n2 <= M_local) then
            if (tapwDebug) call MIO_Print('DEBUG Option A: Velocity matrix elements:','diag')
            call MIO_Print('  |Vx_TB('//trim(num2str(n1))//','//trim(num2str(n2))//')| = '//trim(num2str(abs(Vx_TB(n1,n2)),8)),'diag')
            call MIO_Print('  |Vy_TB('//trim(num2str(n1))//','//trim(num2str(n2))//')| = '//trim(num2str(abs(Vy_TB(n1,n2)),8)),'diag')
         end if
      end block
   end if

   ! Step 5: Calculate Berry curvature for requested bands (same formula as Option B)
   ! Make regularization parameter configurable for testing
   call MIO_InputParameter('Diag.ChernEps', eps, 1.0e-10_dp)
   if (tapwDebug) call MIO_Print('DEBUG: Using eps = '//trim(num2str(eps,12))//' for Berry curvature calculation','diag')
   ! Make energy threshold configurable for testing
   call MIO_InputParameter('Diag.ChernEnergyThreshold', energy_threshold, 1.0e-4_dp)

   ! Add proper energy broadening parameter in eV
   call MIO_InputParameter('Diag.ChernEtaEV', eta_eV, 1.0e-6_dp)
   berry_curv_bands = 0.0_dp

   ! DEBUG: Count how many band pairs are skipped
   skipped_pairs = 0
   total_pairs = 0

   call MIO_Print('Computing Berry curvature for '//trim(num2str(num_bands))//' bands in TB space','diag')

   ! DEBUG: Compare velocity matrix elements between Option A and Option B
   if (kpoint_index == 1) then
      if (tapwDebug) call MIO_Print('DEBUG: Comparing velocity matrices between Option A and Option B','diag')
      call MIO_Print('  Option A uses Vx_TB, Vy_TB (TB basis)','diag')
      call MIO_Print('  Option B uses Vx, Vy (TAPW basis)','diag')
      call MIO_Print('  This comparison will help identify the root cause','diag')
   end if

   ! CRITICAL FIX: Use actual array dimensions for loop bounds, same as Option B
   M_safe = min(size(Vx_TB,1), size(Vx_TB,2), size(Vy_TB,1), size(Vy_TB,2), size(eigval))

   ! Additional safety check before Berry curvature calculation
   if (kpoint_index == 1) then
      call MIO_Print('Debug Option A: About to calculate Berry curvature','diag')
      call MIO_Print('  M_local variable value: '//trim(num2str(M_local)),'diag')
      call MIO_Print('  M_safe (actual array bounds): '//trim(num2str(M_safe)),'diag')
      call MIO_Print('  Vx_TB matrix shape: '//trim(num2str(size(Vx_TB,1)))//'x'//trim(num2str(size(Vx_TB,2))),'diag')
      call MIO_Print('  Vy_TB matrix shape: '//trim(num2str(size(Vy_TB,1)))//'x'//trim(num2str(size(Vy_TB,2))),'diag')
      call MIO_Print('  eigval array size: '//trim(num2str(size(eigval))),'diag')
      call MIO_Print('  First band index to access: '//trim(num2str(band_indices(1))),'diag')
      call MIO_Print('  Last band index to access: '//trim(num2str(band_indices(size(band_indices)))),'diag')

      if (M_local /= M_safe) then
         call MIO_Print('CRITICAL WARNING: M_local ('//trim(num2str(M_local))//') differs from actual array bounds ('//trim(num2str(M_safe))//')','diag')
         call MIO_Print('This indicates M_local was corrupted after array allocation - using M_safe for loops','diag')
      end if
   end if

   do i = 1, size(band_indices)  ! For each selected band (same as Option B)
      n = band_indices(i)  ! Current band index
      if (n <= M_safe .and. n >= 1) then  ! Use M_safe instead of M_local (same as Option B)
         if (kpoint_index == 1 .and. i <= 3) then  ! Debug first few bands (same as Option B)
            call MIO_Print('Debug Option A: Processing band '//trim(num2str(i))//' (index '//trim(num2str(n))//')','diag')
         end if

         berry_sum = (0.0_dp, 0.0_dp)
         do m = 1, M_safe  ! Use M_safe instead of M_local (same as Option B)
            if (n /= m) then
               ! Skip nearly degenerate bands to avoid numerical instability (same as Option B)
               if (abs(eigval(n) - eigval(m)) < energy_threshold) then
                  total_pairs = total_pairs + 1
                  skipped_pairs = skipped_pairs + 1
                  cycle  ! Skip this pair - too close in energy
               end if

               ! Additional bounds checking to prevent array access errors (same as Option B)
               if (n <= size(Vx_TB,1) .and. m <= size(Vx_TB,2) .and. m <= size(Vy_TB,1) .and. n <= size(Vy_TB,2) .and. &
                   n <= size(eigval) .and. m <= size(eigval)) then
                  total_pairs = total_pairs + 1
                  ! CORRECTED: Use Vy(m,n) to match Python's (Vynm.T) - transpose operation (no conjg!) (same as Option B)
                  ! Python: (Vxnm*(Vynm.T)) means Vx(n,m) * Vy(m,n) - transpose, not conjugate transpose
                  berry_sum = berry_sum + (Vx_TB(n,m) * Vy_TB(m,n)) / ((eigval(n) - eigval(m))**2 + eps)
               else
                  call MIO_Print('WARNING: Array bounds exceeded in Berry calculation','diag')
                  call MIO_Print('  n='//trim(num2str(n))//', m='//trim(num2str(m))//', M_local='//trim(num2str(M_local)),'diag')
                  call MIO_Print('  Vx_TB size: '//trim(num2str(size(Vx_TB,1)))//'x'//trim(num2str(size(Vx_TB,2))),'diag')
                  call MIO_Print('  Vy_TB size: '//trim(num2str(size(Vy_TB,1)))//'x'//trim(num2str(size(Vy_TB,2))),'diag')
                  call MIO_Print('  eigval size: '//trim(num2str(size(eigval))),'diag')
               end if
            end if
         end do

         berry_curv_bands(i) = -2.0_dp * aimag(berry_sum)

         ! DEBUG: Check actual energy differences being processed
         if (kpoint_index == 1 .and. i == 1) then
            if (tapwDebug) call MIO_Print('DEBUG: Checking energy differences for first band:','diag')
            block
               real(dp) :: energy_diff
               do m = 1, min(10, M_safe)
                  if (n /= m .and. n <= size(eigval) .and. m <= size(eigval)) then
                     energy_diff = abs(eigval(n) - eigval(m))
                     call MIO_Print('  |E('//trim(num2str(n))//') - E('//trim(num2str(m))//')| = '//trim(num2str(energy_diff,12)),'diag')
                     if (energy_diff < 1.0e-6_dp) then
                        call MIO_Print('    *** VERY SMALL ENERGY DIFFERENCE ***','diag')
                     end if
                  end if
               end do
            end block
         end if

         ! DEBUG: Compare velocity matrix elements for first few bands
         if (kpoint_index == 1 .and. i <= 3) then
            if (tapwDebug) call MIO_Print('DEBUG Option A: Velocity matrix elements for band '//trim(num2str(n))//':','diag')
            call MIO_Print('  |Vx_TB('//trim(num2str(n))//','//trim(num2str(n+1))//')| = '//trim(num2str(abs(Vx_TB(n,n+1)),8)),'diag')
            call MIO_Print('  |Vy_TB('//trim(num2str(n))//','//trim(num2str(n+1))//')| = '//trim(num2str(abs(Vy_TB(n,n+1)),8)),'diag')
            call MIO_Print('  Berry curvature = '//trim(num2str(berry_curv_bands(i),8)),'diag')
         end if
      end if
   end do

   call MIO_Print('Option A Berry curvature calculation completed','diag')

   ! DEBUG: Report statistics
   if (kpoint_index == 1) then
      if (tapwDebug) call MIO_Print('DEBUG Option A: Energy threshold statistics:','diag')
      call MIO_Print('  Total band pairs considered: '//trim(num2str(total_pairs)),'diag')
      call MIO_Print('  Pairs skipped (|ΔE| < '//trim(num2str(energy_threshold,6))//'): '//trim(num2str(skipped_pairs)),'diag')
      call MIO_Print('  Skipped percentage: '//trim(num2str(100.0_dp * skipped_pairs / max(total_pairs,1),2))//'%','diag')
      call MIO_Print('  Grid density: '//trim(num2str(nk_chern_x))//'x'//trim(num2str(nk_chern_y)),'diag')
      call MIO_Print('  Total k-points: '//trim(num2str(nk_chern_x * nk_chern_y)),'diag')
   end if

   ! Debug output for first few bands
   if (num_bands >= 1) then
      call MIO_Print('Option A Berry curvature samples:','diag')
      do i = 1, min(3, num_bands)
         call MIO_Print('  Band '//trim(num2str(band_indices(i)))//': Ω = '//trim(num2str(berry_curv_bands(i),8)),'diag')
      end do
   end if

   ! Clean up
   deallocate(H_k, eigvec_TAPW, eigvec_TB, eigval)
   deallocate(dH_TB_dkx, dH_TB_dky, Vx_TB, Vy_TB)

end subroutine CalculateBerryAtKpointFromStored_OptionA

!> @brief SOC-aware version of CalculateBerryAtKpointFromStored_OptionA
!! @details Uses block Hamiltonian derivatives for SOC calculations
subroutine CalculateBerryAtKpointFromStored_OptionA_withSOC(kpoint_index, spin_index, band_indices, eigval_from_E, berry_curv_bands)
   ! SOC-aware version: Transform TAPW eigenvectors to TB basis, then calculate Berry curvature in TB space
   ! Uses block Hamiltonian derivatives ONLY when Rashba is enabled (requires 2N×M X matrix)
   ! For non-Rashba SOC (Zeeman, Intrinsic), falls back to non-block version (N×M X matrix)
   ! Updated to include spin_index for separate spin channel calculations
   use constants, only : cmplx_i, cmplx_0, cmplx_1, pi
   use atoms, only : nAt
   use ham, only : H0, hopp, RashbaSOCterm
   use neigh, only : NList, Nneigh, neighCell, maxNeigh
   use cell, only : ucell
   implicit none

   ! Input parameters
   integer, intent(in) :: kpoint_index, spin_index
   integer, intent(in) :: band_indices(:)
   real(dp), intent(in) :: eigval_from_E(:)
   real(dp), intent(out) :: berry_curv_bands(:)

   ! Local variables (all declarations must come before executable statements)
   integer :: M_local, N_local, i, j, n, m, num_bands, iband
   integer :: X_dim, skipped_pairs, total_pairs, M_safe
   complex(dp), allocatable :: H_k(:,:), eigvec_TAPW(:,:)
   real(dp), allocatable :: eigval(:)
   complex(dp), allocatable :: dH_Block_dkx(:,:), dH_Block_dky(:,:)  ! Block derivatives for Rashba
   complex(dp), allocatable :: dH_TB_dkx(:,:), dH_TB_dky(:,:)  ! Regular TB derivatives for non-Rashba SOC
   complex(dp), allocatable :: Vx_TB(:,:), Vy_TB(:,:)
   complex(dp), allocatable :: eigvec_TB(:,:)
   complex(dp), allocatable :: temp_TB(:,:)
   real(dp) :: current_kpt(3), eps, energy_diff, energy_threshold, eta_eV
   complex(dp) :: berry_sum
   logical :: use_block_derivatives

   ! Check if Rashba is enabled (determines if we use block or regular derivatives)
   use_block_derivatives = RashbaSOCterm

   if (use_block_derivatives) then
      call MIO_Print('=== SOC-AWARE OPTION A (Rashba): Berry curvature with block derivatives ===','diag')
   else
      call MIO_Print('=== SOC-AWARE OPTION A (non-Rashba): Berry curvature with regular TB derivatives ===','diag')
   end if

   ! Check stored data availability
   if (.not. allocated(stored_X_matrix)) then
      call MIO_Print('ERROR: Stored X matrix not available for SOC Option A','diag')
      error stop 1
   end if

   if (.not. allocated(stored_hamiltonians) .or. .not. allocated(stored_eigenvectors)) then
      call MIO_Print('ERROR: Stored TAPW data not available for SOC Option A','diag')
      error stop 1
   end if

   ! Get dimensions
   N_local = nAt  ! Use module-level nAt
   M_local = stored_M
   num_bands = size(band_indices)

   call MIO_Print('SOC Option A dimensions: N='//trim(num2str(N_local))//', M='//trim(num2str(M_local))//', bands='//trim(num2str(num_bands)),'diag')

   ! Check X matrix dimensions to determine if block or regular version should be used
   ! For Rashba: stored_X_matrix is 2N×M (from DiagH0TAPW_withBlockH)
   ! For non-Rashba SOC: stored_X_matrix is N×M (from DiagH0TAPW)
   if (use_block_derivatives) then
      if (size(stored_X_matrix, 1) /= 2*N_local) then
         call MIO_Print('WARNING: Rashba enabled but X matrix is N×M (not 2N×M). This suggests non-Rashba TAPW was used.','diag')
         call MIO_Print('Falling back to regular TB derivatives (non-block version).','diag')
         use_block_derivatives = .false.
      end if
   else
      if (size(stored_X_matrix, 1) == 2*N_local) then
         call MIO_Print('WARNING: Non-Rashba SOC but X matrix is 2N×M (block version). Using block derivatives.','diag')
         use_block_derivatives = .true.
      end if
   end if

   ! Allocate arrays
   allocate(H_k(M_local,M_local), eigvec_TAPW(M_local,M_local))
   allocate(eigval(M_local))
   allocate(Vx_TB(M_local,M_local), Vy_TB(M_local,M_local))

   if (use_block_derivatives) then
      allocate(dH_Block_dkx(2*N_local,2*N_local), dH_Block_dky(2*N_local,2*N_local))  ! Block derivatives
   else
      allocate(dH_TB_dkx(N_local,N_local), dH_TB_dky(N_local,N_local))  ! Regular TB derivatives
   end if

   ! Get current k-point coordinates
   if (kpoint_index <= size(stored_Kpts, 2)) then
      current_kpt = stored_Kpts(:, kpoint_index)
   else
      call MIO_Print('Warning: kpoint_index exceeds stored k-points, using Γ point','diag')
      current_kpt = [0.0_dp, 0.0_dp, 0.0_dp]
   end if

   call MIO_Print('Using k-point coordinates: ['//trim(num2str(current_kpt(1),6))//','//&
                 trim(num2str(current_kpt(2),6))//','//trim(num2str(current_kpt(3),6))//']','diag')

   ! Step 1: Compute Hamiltonian derivatives (block or regular depending on Rashba)
   if (use_block_derivatives) then
   call compute_block_hamiltonian_derivatives(dH_Block_dkx, dH_Block_dky, N_local, current_kpt, ucell, &
                                             H0, maxNeigh, hopp, NList, Nneigh, neighCell)
   else
      call compute_TB_hamiltonian_derivatives(dH_TB_dkx, dH_TB_dky, N_local, current_kpt, &
                                             H0, maxNeigh, hopp, NList, Nneigh, neighCell)
   end if

   ! Check spin_index bounds
   if (spin_index < 1 .or. spin_index > stored_nspin) then
      call MIO_Print('Error: spin_index '//trim(num2str(spin_index))//' out of range [1,'//trim(num2str(stored_nspin))//']','diag')
      berry_curv_bands = 0.0_dp
      return
   end if

   ! Get stored TAPW data for this k-point and spin channel
   H_k(1:M_local,1:M_local) = stored_hamiltonians(1:M_local,1:M_local,kpoint_index,spin_index)
   eigvec_TAPW(1:M_local,1:M_local) = stored_eigenvectors(1:M_local,1:M_local,kpoint_index,spin_index)

   ! CRITICAL: Use TAPW eigenvalues instead of TB eigenvalues (matches non-SOC Option A)
   ! The eigenvectors are from TAPW diagonalization, so eigenvalues must be too
   if (allocated(stored_eigenvalues)) then
      if (size(stored_eigenvalues,1) >= M_local .and. size(stored_eigenvalues,2) >= kpoint_index .and. size(stored_eigenvalues,3) >= spin_index) then
         eigval(1:M_local) = stored_eigenvalues(1:M_local,kpoint_index,spin_index)
         call MIO_Print('Using stored TAPW eigenvalues for k-point '//trim(num2str(kpoint_index))//', spin '//trim(num2str(spin_index)),'diag')
         call MIO_Print('  TAPW eigenvalues range: ['//trim(num2str(minval(eigval(1:M_local)),6))//','//trim(num2str(maxval(eigval(1:M_local)),6))//']','diag')
      else
         call MIO_Print('ERROR: stored_eigenvalues array too small','diag')
         call MIO_Print('  Array size: '//trim(num2str(size(stored_eigenvalues,1)))//'x'//trim(num2str(size(stored_eigenvalues,2)))//'x'//trim(num2str(size(stored_eigenvalues,3))),'diag')
         call MIO_Print('  Required: '//trim(num2str(M_local))//'x'//trim(num2str(kpoint_index))//'x'//trim(num2str(spin_index)),'diag')
         error stop 1
      end if
   else
      call MIO_Print('ERROR: No stored TAPW eigenvalues available - using TB eigenvalues (WRONG!)','diag')
      call MIO_Print('This will cause eigenvector-eigenvalue misalignment!','diag')
      ! Fallback to eigval_from_E (this is wrong but prevents crash)
      if (size(eigval_from_E) >= M_local) then
         eigval(1:M_local) = eigval_from_E(1:M_local)
      else
         eigval(1:size(eigval_from_E)) = eigval_from_E
         eigval(size(eigval_from_E)+1:M_local) = 0.0_dp
      end if
   end if

   ! Step 2: Transform TAPW eigenvectors back to TB basis
   ! For Rashba: transform to 2N-dimensional block TB basis (2N×M)
   ! For non-Rashba SOC: transform to N-dimensional TB basis (N×M)
   X_dim = size(stored_X_matrix, 1)  ! N for non-Rashba, 2N for Rashba

   allocate(eigvec_TB(X_dim, M_local))
   call transform_TAPW_eigenvectors_to_TB(eigvec_TAPW, stored_X_matrix, X_dim, M_local, eigvec_TB)

   ! Step 3: Calculate velocity matrix elements in TB space (matches non-SOC Option A pattern)
   ! For Rashba: Vx_TB(n,m) = ⟨ψ_TB_n| dH^Block/dkx |ψ_TB_m⟩ using 2N×2N block derivatives
   ! For non-Rashba: Vx_TB(n,m) = ⟨ψ_TB_n| dH^TB/dkx |ψ_TB_m⟩ using N×N regular derivatives

   if (use_block_derivatives) then
      call MIO_Print('Computing velocity matrix elements in TB block space (2N×2N)','diag')
      ! Note: eigvec_TB is 2N×M, dH_Block_dk is 2N×2N, so we compute M×M velocity matrices
   else
      call MIO_Print('Computing velocity matrix elements in TB space (N×N)','diag')
      ! Note: eigvec_TB is N×M, dH_TB_dk is N×N, so we compute M×M velocity matrices
   end if

   ! Use temporary array for intermediate result
   allocate(temp_TB(X_dim, M_local))

   if (use_block_derivatives) then
      ! Step 3a: temp = dH_Block_dkx * eigvec_TB
      call zgemm('N', 'N', 2*N_local, M_local, 2*N_local, cmplx_1, dH_Block_dkx, 2*N_local, eigvec_TB, 2*N_local, cmplx_0, temp_TB, 2*N_local)

      ! Step 3b: Vx_TB = eigvec_TB† * temp
      call zgemm('C', 'N', M_local, M_local, 2*N_local, cmplx_1, eigvec_TB, 2*N_local, temp_TB, 2*N_local, cmplx_0, Vx_TB, M_local)

      ! Repeat for y-direction
      call zgemm('N', 'N', 2*N_local, M_local, 2*N_local, cmplx_1, dH_Block_dky, 2*N_local, eigvec_TB, 2*N_local, cmplx_0, temp_TB, 2*N_local)
      call zgemm('C', 'N', M_local, M_local, 2*N_local, cmplx_1, eigvec_TB, 2*N_local, temp_TB, 2*N_local, cmplx_0, Vy_TB, M_local)
   else
      ! Step 3a: temp = dH_TB_dkx * eigvec_TB
      call zgemm('N', 'N', N_local, M_local, N_local, cmplx_1, dH_TB_dkx, N_local, eigvec_TB, N_local, cmplx_0, temp_TB, N_local)

      ! Step 3b: Vx_TB = eigvec_TB† * temp
      call zgemm('C', 'N', M_local, M_local, N_local, cmplx_1, eigvec_TB, N_local, temp_TB, N_local, cmplx_0, Vx_TB, M_local)

      ! Repeat for y-direction
      call zgemm('N', 'N', N_local, M_local, N_local, cmplx_1, dH_TB_dky, N_local, eigvec_TB, N_local, cmplx_0, temp_TB, N_local)
      call zgemm('C', 'N', M_local, M_local, N_local, cmplx_1, eigvec_TB, N_local, temp_TB, N_local, cmplx_0, Vy_TB, M_local)
   end if

   deallocate(temp_TB, eigvec_TB)

   ! Step 4: Set diagonal elements to zero (matches non-SOC Option A)
   do i = 1, M_local
      Vx_TB(i,i) = cmplx_0
      Vy_TB(i,i) = cmplx_0
   end do

   ! Step 5: Calculate Berry curvature for each band (matches non-SOC Option A)
   ! Ω_n = -2 Im[Σ_{m≠n} (Vx_nm * Vy_mn) / (E_n - E_m)²]
   call MIO_InputParameter('Diag.ChernEps', eps, 1.0e-10_dp)
   call MIO_InputParameter('Diag.ChernEnergyThreshold', energy_threshold, 1.0e-4_dp)
   call MIO_InputParameter('Diag.ChernEtaEV', eta_eV, 1.0e-6_dp)

   berry_curv_bands = 0.0_dp
   skipped_pairs = 0
   total_pairs = 0

   do i = 1, num_bands
      iband = band_indices(i)
      berry_sum = cmplx_0

      do j = 1, M_local
         if (j /= iband) then
            total_pairs = total_pairs + 1
            energy_diff = eigval(iband) - eigval(j)

            ! Skip if energy difference is too small (avoid numerical issues)
            if (abs(energy_diff) < energy_threshold) then
               skipped_pairs = skipped_pairs + 1
               cycle
            end if

               ! Berry curvature contribution: Vx_nm * Vy_mn / (E_n - E_m)²
            berry_sum = berry_sum + (Vx_TB(iband, j) * Vy_TB(j, iband)) / (energy_diff * energy_diff + eps)
         end if
      end do

      ! Berry curvature: Ω_n = -2 Im[sum]
      berry_curv_bands(i) = -2.0_dp * aimag(berry_sum)
   end do

   if (tapwDebug) then
      call MIO_Print('SOC Option A: Skipped '//trim(num2str(skipped_pairs))//'/'//trim(num2str(total_pairs))//' band pairs due to small energy differences','diag')
   end if

   call MIO_Print('SOC Option A: Berry curvature calculation completed for '//trim(num2str(num_bands))//' bands','diag')

   ! Clean up
   deallocate(H_k, eigvec_TAPW, eigval)
   if (use_block_derivatives) then
      deallocate(dH_Block_dkx, dH_Block_dky, Vx_TB, Vy_TB)
   else
      deallocate(dH_TB_dkx, dH_TB_dky, Vx_TB, Vy_TB)
   end if

end subroutine CalculateBerryAtKpointFromStored_OptionA_withSOC

subroutine cleanup_position_differences_cache()
   ! Clean up cached position difference matrices
   implicit none

   if (allocated(cached_delX)) deallocate(cached_delX, cached_delY)
   if (allocated(cached_Gx)) deallocate(cached_Gx, cached_Gy)
   cache_valid = .false.
   cached_M = -1
   cached_NG = -1
   cached_Nlabel = -1

   call MIO_Print('TAPW position difference cache cleared','diag')

end subroutine cleanup_position_differences_cache

! Numerical check for TAPW unitary transformation
! This routine verifies that when G covers all moiré reciprocal lattice vectors
! in the graphene first BZ, the transformation X†TX = T (restores full TB model)
subroutine verify_tapw_unitary_transformation(N, M, row_ptr, col_ind, values, X, Hproj, NG, Nlabel, Gx, Gy)
    use constants, only : cmplx_0, cmplx_1
    use mio
    implicit none

    ! Input parameters
    integer, intent(in) :: N, M, NG, Nlabel
    integer, intent(in) :: row_ptr(N+1), col_ind(:)
    complex(dp), intent(in) :: values(:), X(N,M), Hproj(M,M)
    double precision, intent(in) :: Gx(NG), Gy(NG)

    ! Local variables
    complex(dp), allocatable :: T_full(:,:), X_dagger_X(:,:), X_X_dagger(:,:), reconstruction(:,:), temp(:,:)
    real(dp) :: unitary_error, reconstruction_error, projection_error
    real(dp) :: max_XdX_error
    real(dp) :: Gx_i, Gy_i, Gx_j, Gy_j, mag_i, mag_j, diff_x, diff_y, diff_mag
    integer :: i, j, k, nnz, max_XdX_i, max_XdX_j
    integer :: G_i, G_j, label_i, label_j

    call MIO_Print('=== TAPW Unitary Transformation Check ===','diag')

    ! 1. Check if X†X is close to identity (unitarity check)
    allocate(X_dagger_X(M, M))
    call zgemm('C', 'N', M, M, N, cmplx_1, X, N, X, N, cmplx_0, X_dagger_X, M)

    ! Calculate error from identity
    unitary_error = 0.0_dp
    do i = 1, M
        do j = 1, M
            if (i == j) then
                unitary_error = max(unitary_error, abs(X_dagger_X(i,j) - cmplx_1))
            else
                unitary_error = max(unitary_error, abs(X_dagger_X(i,j)))
            end if
        end do
    end do

    call MIO_Print('X†X unitarity check:','diag')
    call MIO_Print('  Max deviation from identity: '//trim(num2str(unitary_error,8)),'diag')

    if (unitary_error < 1.0e-10_dp) then
        call MIO_Print('  ✓ X is numerically unitary (excellent)','diag')
    else if (unitary_error < 1.0e-6_dp) then
        call MIO_Print('  ✓ X is approximately unitary (good)','diag')
    else if (unitary_error < 1.0e-3_dp) then
        call MIO_Print('  ⚠ X has some non-unitarity (acceptable for truncated G-set)','diag')
    else
        call MIO_Print('  ✗ X is significantly non-unitary (check G-vector completeness)','diag')
    end if

    ! 1b. Also check X*X† (should be projection matrix for rectangular X)
    allocate(X_X_dagger(N, N))
    call zgemm('N', 'C', N, N, M, cmplx_1, X, N, X, N, cmplx_0, X_X_dagger, N)

    ! For a unitary matrix, X*X† should equal identity
    ! For rectangular matrix (M<N), X*X† is a projection matrix
    projection_error = 0.0_dp
    do i = 1, N
        do j = 1, N
            if (i == j) then
                projection_error = max(projection_error, abs(X_X_dagger(i,j) - cmplx_1))
            else
                projection_error = max(projection_error, abs(X_X_dagger(i,j)))
            end if
        end do
    end do

    call MIO_Print('X*X† projection check:','diag')
    call MIO_Print('  Max deviation from identity: '//trim(num2str(projection_error,8)),'diag')

    ! Find where the maximum errors actually occur
    call MIO_Print('Locating maximum unitarity errors:','diag')

    ! Find max error location in X†X
    outer_loop_XdX: do i = 1, M
        do j = 1, M
            if (i == j) then
                if (abs(abs(X_dagger_X(i,j) - cmplx_1) - unitary_error) < 1.0e-12_dp) then
                    call MIO_Print('  X†X max diagonal error at ('//trim(num2str(i))//','//trim(num2str(j))//'): '//&
                                   trim(num2str(real(X_dagger_X(i,j)),8))//' (should be 1.0)','diag')
                    exit outer_loop_XdX
                end if
            else
                if (abs(abs(X_dagger_X(i,j)) - unitary_error) < 1.0e-12_dp) then
                    call MIO_Print('  X†X max off-diagonal error at ('//trim(num2str(i))//','//trim(num2str(j))//'): '//&
                                   trim(num2str(real(X_dagger_X(i,j)),8))//' (should be 0.0)','diag')
                    exit outer_loop_XdX
                end if
            end if
        end do
    end do outer_loop_XdX

    ! Find max error location in X*X† (with tolerance for floating point precision)
    outer_loop_XXd: do i = 1, N
        do j = 1, N
            if (i == j) then
                if (abs(abs(X_X_dagger(i,j) - cmplx_1) - projection_error) < 1.0e-12_dp) then
                    call MIO_Print('  X*X† max diagonal error at ('//trim(num2str(i))//','//trim(num2str(j))//'): '//&
                                   trim(num2str(real(X_X_dagger(i,j)),8))//' (should be 1.0)','diag')
                    exit outer_loop_XXd
                end if
            else
                if (abs(abs(X_X_dagger(i,j)) - projection_error) < 1.0e-12_dp) then
                    call MIO_Print('  X*X† max off-diagonal error at ('//trim(num2str(i))//','//trim(num2str(j))//'): '//&
                                   trim(num2str(real(X_X_dagger(i,j)),8))//' (should be 0.0)','diag')
                    exit outer_loop_XXd
                end if
            end if
        end do
    end do outer_loop_XXd

    ! Decode the problematic G-vector indices and coordinates for better diagnostics
    call MIO_Print('G-vector pair analysis for maximum errors:','diag')

    ! Find the maximum off-diagonal error again to get indices
    max_XdX_error = 0.0_dp
    max_XdX_i = 0
    max_XdX_j = 0
    do i = 1, M
        do j = 1, M
            if (i /= j .and. abs(X_dagger_X(i,j)) > max_XdX_error) then
                max_XdX_error = abs(X_dagger_X(i,j))
                max_XdX_i = i
                max_XdX_j = j
            end if
        end do
    end do

    if (max_XdX_i > 0 .and. max_XdX_j > 0) then
        G_i = (max_XdX_i - 1) / Nlabel + 1
        label_i = mod(max_XdX_i - 1, Nlabel) + 1
        G_j = (max_XdX_j - 1) / Nlabel + 1
        label_j = mod(max_XdX_j - 1, Nlabel) + 1

        call MIO_Print('Problematic G-vector pair:','diag')
        call MIO_Print('  Matrix indices: ('//trim(num2str(max_XdX_i))//', '//trim(num2str(max_XdX_j))//')','diag')
        call MIO_Print('  G-vector indices: G('//trim(num2str(G_i))//'), G('//trim(num2str(G_j))//')','diag')
        call MIO_Print('  Sublattice labels: '//trim(num2str(label_i))//', '//trim(num2str(label_j)),'diag')
        if (G_i <= NG .and. G_j <= NG) then
            ! Use temporary variables to avoid complex expressions in num2str
            Gx_i = Gx(G_i)
            Gy_i = Gy(G_i)
            Gx_j = Gx(G_j)
            Gy_j = Gy(G_j)
            mag_i = sqrt(Gx_i**2 + Gy_i**2)
            mag_j = sqrt(Gx_j**2 + Gy_j**2)
            diff_x = Gx_i - Gx_j
            diff_y = Gy_i - Gy_j
            diff_mag = sqrt(diff_x**2 + diff_y**2)

            ! Use simple print statements to avoid num2str issues
            print *, '  G(', G_i, ') = [', Gx_i, ',', Gy_i, ']'
            print *, '  G(', G_j, ') = [', Gx_j, ',', Gy_j, ']'
            print *, '  |G(', G_i, ')| =', mag_i
            print *, '  |G(', G_j, ')| =', mag_j
            print *, '  G-vector difference: [', diff_x, ',', diff_y, ']'
            print *, '  |ΔG| =', diff_mag
        else
            call MIO_Print('  ERROR: G-vector indices out of bounds!','diag')
        end if
    end if

    ! Output multiple regions of X†X matrix for comparison
    ! Region 1: First 100x100
    open(unit=95, file='XdaggerX_matrix_region1.dat', status='replace')
    write(95, '(A)') '# X†X matrix (real parts) - First 100x100 region'
    write(95, '(A)') '# Should be close to identity matrix'
    do i = 1, min(M, 100)
        write(95, '(*(F12.6))') (real(X_dagger_X(i,j)), j = 1, min(M, 100))
    end do
    close(95)

    ! Region 2: Middle region around M/2
    if (M > 200) then
        open(unit=93, file='XdaggerX_matrix_region2.dat', status='replace')
        write(93, '(A)') '# X†X matrix (real parts) - Middle region'
        write(93, '(A)') '# Should be close to identity matrix'
        do i = max(1, M/2-50), min(M, M/2+49)
            write(93, '(*(F12.6))') (real(X_dagger_X(i,j)), j = max(1, M/2-50), min(M, M/2+49))
        end do
        close(93)
    end if

    ! Region 3: Last 100x100
    if (M > 100) then
        open(unit=92, file='XdaggerX_matrix_region3.dat', status='replace')
        write(92, '(A)') '# X†X matrix (real parts) - Last 100x100 region'
        write(92, '(A)') '# Should be close to identity matrix'
        do i = max(1, M-99), M
            write(92, '(*(F12.6))') (real(X_dagger_X(i,j)), j = max(1, M-99), M)
        end do
        close(92)
    end if

    ! Output X*X† matrix regions
    open(unit=94, file='XXdagger_matrix_region1.dat', status='replace')
    write(94, '(A)') '# X*X† matrix (real parts) - First 100x100 region'
    write(94, '(A)') '# Should be close to identity matrix'
    do i = 1, min(N, 100)
        write(94, '(*(F12.6))') (real(X_X_dagger(i,j)), j = 1, min(N, 100))
    end do
    close(94)

    ! Calculate statistics for different regions
    call MIO_Print('Regional unitarity analysis:','diag')

    ! First 100x100 region statistics
    unitary_error = 0.0_dp
    do i = 1, min(M, 100)
        do j = 1, min(M, 100)
            if (i == j) then
                unitary_error = max(unitary_error, abs(X_dagger_X(i,j) - cmplx_1))
            else
                unitary_error = max(unitary_error, abs(X_dagger_X(i,j)))
            end if
        end do
    end do
    call MIO_Print('  First 100×100 X†X region max error: '//trim(num2str(unitary_error,8)),'diag')

    ! Last 100x100 region statistics (if matrix is large enough)
    if (M > 100) then
        unitary_error = 0.0_dp
        do i = max(1, M-99), M
            do j = max(1, M-99), M
                if (i == j) then
                    unitary_error = max(unitary_error, abs(X_dagger_X(i,j) - cmplx_1))
                else
                    unitary_error = max(unitary_error, abs(X_dagger_X(i,j)))
                end if
            end do
        end do
        call MIO_Print('  Last 100×100 X†X region max error: '//trim(num2str(unitary_error,8)),'diag')
    end if

    call MIO_Print('Matrix visualization files written:','diag')
    call MIO_Print('  XdaggerX_matrix_region1.dat (first 100×100)','diag')
    if (M > 200) call MIO_Print('  XdaggerX_matrix_region2.dat (middle region)','diag')
    if (M > 100) call MIO_Print('  XdaggerX_matrix_region3.dat (last 100×100)','diag')
    call MIO_Print('  XXdagger_matrix_region1.dat (first 100×100)','diag')

    if (M == N) then
        if (projection_error < 1.0e-10_dp) then
            call MIO_Print('  ✓ X*X† is numerically identity (square unitary matrix)','diag')
        else if (projection_error < 1.0e-6_dp) then
            call MIO_Print('  ✓ X*X† is approximately identity (good square matrix)','diag')
        else if (projection_error < 1.0e-3_dp) then
            call MIO_Print('  ⚠ X*X† has some deviation from identity','diag')
        else
            call MIO_Print('  ✗ X*X† significantly deviates from identity','diag')
        end if
    else
        call MIO_Print('  Note: X is rectangular ('//trim(num2str(N))//'×'//trim(num2str(M))//'), X*X† is projection matrix','diag')
    end if

    ! 2. Reconstruct full Hamiltonian: T_reconstruct = X * Hproj * X†
    allocate(T_full(N, N), reconstruction(N, N))

    ! Convert sparse matrix to full matrix for comparison
    T_full = cmplx_0
    do i = 1, N
        do k = row_ptr(i), row_ptr(i+1) - 1
            j = col_ind(k)
            T_full(i, j) = values(k)
        end do
    end do

    ! Compute X * Hproj * X†
    allocate(temp(N, M))

    call zgemm('N', 'N', N, M, M, cmplx_1, X, N, Hproj, M, cmplx_0, temp, N)
    call zgemm('N', 'C', N, N, M, cmplx_1, temp, N, X, N, cmplx_0, reconstruction, N)

    ! Calculate reconstruction error
    reconstruction_error = 0.0_dp
    do i = 1, N
        do j = 1, N
            reconstruction_error = max(reconstruction_error, abs(T_full(i,j) - reconstruction(i,j)))
        end do
    end do

    call MIO_Print('','diag')
    call MIO_Print('Full Hamiltonian reconstruction check:','diag')
    call MIO_Print('  Max |T_original - X*H_proj*X†|: '//trim(num2str(reconstruction_error,8)),'diag')

    if (reconstruction_error < 1.0e-10_dp) then
        call MIO_Print('  ✓ Perfect reconstruction - G-set is complete!','diag')
    else if (reconstruction_error < 1.0e-6_dp) then
        call MIO_Print('  ✓ Excellent reconstruction - G-set is nearly complete','diag')
    else if (reconstruction_error < 1.0e-3_dp) then
        call MIO_Print('  ⚠ Good reconstruction - G-set captures main physics','diag')
    else
        call MIO_Print('  ✗ Poor reconstruction - G-set may be incomplete','diag')
        call MIO_Print('    Consider increasing NG or checking G-vector selection','diag')
    end if

    ! 3. Report G-vector coverage information
    call MIO_Print('','diag')
    call MIO_Print('G-vector coverage analysis:','diag')
    call MIO_Print('  Total G-vectors used: '//trim(num2str(real(NG,dp),0)),'diag')
    call MIO_Print('  Labels (sublattices): '//trim(num2str(real(Nlabel,dp),0)),'diag')
    call MIO_Print('  TAPW subspace dimension: '//trim(num2str(real(M,dp),0))//' (vs full space: '//trim(num2str(real(N,dp),0))//')','diag')
    call MIO_Print('  Compression ratio: '//trim(num2str(real(N,dp)/real(M,dp),2)),'diag')

    ! Cleanup
    deallocate(X_dagger_X, X_X_dagger, T_full, reconstruction, temp)

    call MIO_Print('=== End TAPW Unitary Check ===','diag')

end subroutine verify_tapw_unitary_transformation

subroutine diagnose_tapw_phase_consistency(KLoc, Gx, Gy, NG, XArray, N, Nlabel)
    ! Diagnose potential phase factor inconsistencies between X-matrix and Hamiltonian
    use constants, only : cmplx_i
    use atoms, only : Rat
    implicit none

    ! Input parameters
    real(dp), intent(in) :: KLoc(2), Gx(NG), Gy(NG)
    integer, intent(in) :: NG, N, Nlabel
    complex(dp), intent(in) :: XArray(N, NG*Nlabel)

    ! Local variables
    complex(dp) :: X_phase, H_phase, expected_phase
    real(dp) :: G_dot_r, k_dot_r
    integer :: i_atom, i_G

    call MIO_Print('=== TAPW Phase Consistency Diagnostic ===','diag')
    call MIO_Print('Checking consistency between X-matrix and Hamiltonian phases','diag')
    call MIO_Print('Current k-point: ['//trim(num2str(KLoc(1),6))//','//trim(num2str(KLoc(2),6))//']','diag')

    ! Check first few atoms and G-vectors
    do i_atom = 1, min(3, N)
        do i_G = 1, min(3, NG)
            ! Calculate G·r phase (from X-matrix construction)
            G_dot_r = Gx(i_G) * Rat(1, i_atom) + Gy(i_G) * Rat(2, i_atom)
            X_phase = exp(cmplx_i * G_dot_r)

            ! Calculate k·r phase (what Hamiltonian would see)
            k_dot_r = KLoc(1) * Rat(1, i_atom) + KLoc(2) * Rat(2, i_atom)
            H_phase = exp(-cmplx_i * k_dot_r)  ! Note: negative sign like in initialize_sparse_matrix

            ! Expected combined phase for TAPW
            expected_phase = exp(cmplx_i * (G_dot_r - k_dot_r))

            call MIO_Print('Atom '//trim(num2str(i_atom))//', G-vector '//trim(num2str(i_G))//':','diag')
            call MIO_Print('  Position: ['//trim(num2str(Rat(1,i_atom),4))//','//trim(num2str(Rat(2,i_atom),4))//']','diag')
            call MIO_Print('  G·r = '//trim(num2str(G_dot_r,6))//', exp(iG·r) = ('//&
                          trim(num2str(real(X_phase),6))//','//trim(num2str(aimag(X_phase),6))//')','diag')
            call MIO_Print('  k·r = '//trim(num2str(k_dot_r,6))//', exp(-ik·r) = ('//&
                          trim(num2str(real(H_phase),6))//','//trim(num2str(aimag(H_phase),6))//')','diag')
            call MIO_Print('  Combined: exp(i(G-k)·r) = ('//&
                          trim(num2str(real(expected_phase),6))//','//trim(num2str(aimag(expected_phase),6))//')','diag')

            ! Check if k=0 (Gamma point)
            if (sqrt(KLoc(1)**2 + KLoc(2)**2) < 1e-10) then
                call MIO_Print('  → At Γ-point: k·r ≈ 0, phase ≈ 1','diag')
            end if
        end do
    end do

    call MIO_Print('=== End Phase Consistency Diagnostic ===','diag')

end subroutine diagnose_tapw_phase_consistency

subroutine calculate_optimal_NGrange_for_complete_BZ(twist_angle_deg, user_NGrange, optimal_NGrange)
    use constants, only : pi
    implicit none

    ! Input/Output parameters
    real(dp), intent(in) :: twist_angle_deg    ! Twist angle in degrees
    integer, intent(in) :: user_NGrange        ! User-requested NGrange
    integer, intent(out) :: optimal_NGrange    ! Calculated optimal NGrange

    ! Local variables
    real(dp) :: twist_angle_rad, sin_half_theta
    integer :: NGrange_complete, NGrange_safe

    call MIO_Print('=== TAPW G-Vector Grid Completeness Analysis ===','diag')

    ! Convert to radians
    twist_angle_rad = twist_angle_deg * pi / 180.0_dp
    sin_half_theta = sin(twist_angle_rad / 2.0_dp)

    ! Calculate NGrange needed for complete first BZ coverage
    ! Formula: NGrange_complete = ceil(1 / (2*sin(θ/2)))
    ! This ensures G-vectors span the entire graphene first BZ
    NGrange_complete = ceiling(1.0_dp / (2.0_dp * sin_half_theta))

    ! Add safety margin for numerical precision and edge effects
    NGrange_safe = nint(NGrange_complete * 1.2_dp)  ! 20% safety margin

    call MIO_Print('G-vector grid completeness analysis:','diag')
    call MIO_Print('  Twist angle: '//trim(num2str(twist_angle_deg,3))//' degrees','diag')
    call MIO_Print('  sin(θ/2): '//trim(num2str(sin_half_theta,6)),'diag')
    call MIO_Print('  NGrange for COMPLETE BZ coverage: '//trim(num2str(real(NGrange_complete,dp),0)),'diag')
    call MIO_Print('  NGrange with safety margin: '//trim(num2str(real(NGrange_safe,dp),0)),'diag')
    call MIO_Print('  User requested NGrange: '//trim(num2str(real(user_NGrange,dp),0)),'diag')

    ! Determine which NGrange to use
    if (user_NGrange >= NGrange_safe) then
        optimal_NGrange = user_NGrange
        call MIO_Print('  ✓ Using user NGrange (sufficient for complete BZ)','diag')
    else if (user_NGrange >= NGrange_complete) then
        optimal_NGrange = user_NGrange
        call MIO_Print('  ⚠ Using user NGrange (complete but no safety margin)','diag')
    else
        optimal_NGrange = NGrange_safe
        call MIO_Print('  ⚠ User NGrange TOO SMALL - using calculated safe value','diag')
        call MIO_Print('    This ensures complete first BZ coverage for unitary TAPW','diag')
    end if

    call MIO_Print('  FINAL NGrange: '//trim(num2str(real(optimal_NGrange,dp),0)),'diag')

    ! Provide guidance for different regimes
    if (optimal_NGrange > 100) then
        call MIO_Print('','diag')
        call MIO_Print('NOTE: Large NGrange due to small twist angle.','diag')
        call MIO_Print('      Consider using triangular truncation for efficiency.','diag')
        call MIO_Print('      Set useTriangularTruncation=.true. in input.','diag')
    else if (optimal_NGrange > 50) then
        call MIO_Print('','diag')
        call MIO_Print('NOTE: Moderate NGrange - good balance of completeness vs efficiency.','diag')
    else
        call MIO_Print('','diag')
        call MIO_Print('NOTE: Small NGrange - efficient but check unitary test results.','diag')
    end if

    call MIO_Print('=== End G-Vector Grid Analysis ===','diag')

end subroutine calculate_optimal_NGrange_for_complete_BZ

subroutine output_brillouin_zones_debug(aG, moireAngle, k_ref, rG)
   ! Output simple Brillouin zone vertices using the rG vectors from the main calculation
   use constants, only: pi
   implicit none
   real(dp), intent(in) :: aG, moireAngle, k_ref(2), rG(2,2)

   ! Variable declarations
   real(dp) :: bz_vertices(6,2)  ! Hexagonal BZ vertices
   real(dp) :: b1(2), b2(2)      ! Reciprocal lattice vectors
   integer :: i

   if (tapwDebug) call MIO_Print("DEBUG: Writing Brillouin zone data using rG vectors")

   ! Extract reciprocal lattice vectors (rG uses column storage: rG(:,i) = i-th vector)
   b1 = rG(:,1)  ! First reciprocal lattice vector
   b2 = rG(:,2)  ! Second reciprocal lattice vector

   ! Calculate hexagonal BZ vertices in proper order for connecting
   ! For a hexagonal lattice, the BZ vertices in counterclockwise order are:
   ! Starting from K point and going around the hexagon
   bz_vertices(1,:) = (2.0_dp*b1 + b2) / 3.0_dp     ! K point
   bz_vertices(2,:) = (b1 - b2) / 3.0_dp             ! Next vertex clockwise
   bz_vertices(3,:) = -(b1 + 2.0_dp*b2) / 3.0_dp    ! -K' point
   bz_vertices(4,:) = -(2.0_dp*b1 + b2) / 3.0_dp    ! -K point
   bz_vertices(5,:) = -(b1 - b2) / 3.0_dp            ! Next vertex
   bz_vertices(6,:) = (b1 + 2.0_dp*b2) / 3.0_dp     ! K' point

   ! Layer 1: Original graphene BZ
   open(unit=97, file='brillouin_zones_debug_K1.dat', status='replace')
   write(97, '(A)') '# Layer 1 Brillouin zone data for TAPW debugging'
   write(97, '(A)') '# BZ vertices (x, y):'
   do i = 1, 6
      write(97, '(2F16.8)') bz_vertices(i, 1), bz_vertices(i, 2)
   end do
   write(97, '(A)') '# K-point:'
   write(97, '(2F16.8)') k_ref(1), k_ref(2)
   write(97, '(A)') '# Reciprocal lattice vectors b1, b2:'
   write(97, '(2F16.8)') b1(1), b1(2)
   write(97, '(2F16.8)') b2(1), b2(2)
   close(97)

   if (tapwDebug) call MIO_Print("DEBUG: Brillouin zone data written to brillouin_zones_debug.dat")

end subroutine output_brillouin_zones_debug

subroutine output_moire_bz_debug(rcell)
   ! Output moiré Brillouin zone data for visualization
   implicit none
   real(dp), intent(in) :: rcell(3,3)

   real(dp) :: moire_b1(2), moire_b2(2)
   real(dp) :: bz_vertices(6,2)
   integer :: i

   ! Extract 2D moiré reciprocal lattice vectors
   moire_b1 = rcell(1:2, 1)
   moire_b2 = rcell(1:2, 2)

   ! Calculate hexagonal BZ vertices for moiré lattice
   ! Same logic as graphene BZ but using moiré reciprocal vectors
   bz_vertices(1, :) = (2.0_dp/3.0_dp) * moire_b1 + (1.0_dp/3.0_dp) * moire_b2  ! K
   bz_vertices(2, :) = (1.0_dp/3.0_dp) * moire_b1 + (2.0_dp/3.0_dp) * moire_b2  ! K'
   bz_vertices(3, :) = (-1.0_dp/3.0_dp) * moire_b1 + (1.0_dp/3.0_dp) * moire_b2  ! -K+K'
   bz_vertices(4, :) = (-2.0_dp/3.0_dp) * moire_b1 + (-1.0_dp/3.0_dp) * moire_b2  ! -K
   bz_vertices(5, :) = (-1.0_dp/3.0_dp) * moire_b1 + (-2.0_dp/3.0_dp) * moire_b2  ! -K'
   bz_vertices(6, :) = (1.0_dp/3.0_dp) * moire_b1 + (-1.0_dp/3.0_dp) * moire_b2  ! K-K'

   ! Write to file
   open(unit=96, file='moire_bz_debug.dat', status='replace')
   write(96, '(A)') '# Moiré Brillouin zone data for TAPW debugging'
   write(96, '(A)') '# BZ vertices (x, y):'
   do i = 1, 6
      write(96, '(2F16.8)') bz_vertices(i, 1), bz_vertices(i, 2)
   end do

   write(96, '(A)') '# Moire rcell vectors:'
   write(96, '(2F16.8)') moire_b1(1), moire_b1(2)
   write(96, '(2F16.8)') moire_b2(1), moire_b2(2)
   close(96)

   print *, "Moiré BZ data written to moire_bz_debug.dat"

end subroutine output_moire_bz_debug

subroutine output_kpath_debug()
   ! Output k-path data for visualization
   ! This needs to be called from the main routine where path data is available
   implicit none

   ! This is a placeholder - we'll need to modify the main routine to call this
   ! with the actual path data
   print *, "K-path debug output placeholder - needs path data from main routine"

end subroutine output_kpath_debug

subroutine cleanup_TAPW_storage()
   ! Clean up stored TAPW data after Chern calculation
   implicit none

   if (allocated(stored_eigenvectors)) then
      deallocate(stored_eigenvectors)
      call MIO_Print('Deallocated stored eigenvectors','diag')
   end if

   if (allocated(stored_hamiltonians)) then
      deallocate(stored_hamiltonians)
      call MIO_Print('Deallocated stored Hamiltonians','diag')
   end if

   if (allocated(stored_eigenvalues)) then
      deallocate(stored_eigenvalues)
      call MIO_Print('Deallocated stored eigenvalues','diag')
   end if

   stored_M = 0
   stored_nk = 0
   stored_nspin = 0

end subroutine cleanup_TAPW_storage

subroutine write_chern_tapw_results(band_indices, chern_bands, chern_total, nkx, nky, fermi_eV, spin_index)
   ! Write TAPW Chern number results to output file
   ! Updated to include spin_index for separate output files per spin channel
   use name, only : prefix
   implicit none

   ! Input parameters
   integer, intent(in) :: band_indices(:), nkx, nky, spin_index  ! Variable-sized arrays
   real(dp), intent(in) :: chern_bands(:), chern_total, fermi_eV

   ! Local variables
   character(len=256) :: filename
   integer :: unit_chern, i

   ! Create output filename with spin index
   if (spin_index > 1) then
      filename = trim(prefix)//'.chern_tapw_spin'//trim(num2str(spin_index))
   else
   filename = trim(prefix)//'.chern_tapw'
   end if
   unit_chern = 150  ! Use a unique unit number

   call MIO_Print('Writing TAPW Chern results to: '//trim(filename),'diag')

   open(unit=unit_chern, file=filename, status='replace', action='write')

   ! Write header information
   write(unit_chern, '(A)') '# TAPW Chern Number Calculation Results'
   write(unit_chern, '(A)') '# Generated by GRABNES'
   write(unit_chern, '(A)') '#'
   write(unit_chern, '(A,F12.6,A)') '# Fermi Energy: ', fermi_eV, ' eV'
   write(unit_chern, '(A,I0,A,I0)') '# k-grid: ', nkx, ' x ', nky
   write(unit_chern, '(A,I0)') '# Total k-points: ', nkx * nky
   write(unit_chern, '(A)') '#'
   write(unit_chern, '(A)') '# Band_Index  Chern_Number  Rounded_Chern  Deviation'

   ! Write individual band results (variable number of bands)
   do i = 1, size(band_indices)
      write(unit_chern, '(I10,F15.8,I15,F15.8)') band_indices(i), chern_bands(i), &
                                                 nint(chern_bands(i)), abs(chern_bands(i) - nint(chern_bands(i)))
   end do

   ! Write total
   write(unit_chern, '(A)') '#'
   write(unit_chern, '(A,F15.8)') '# Total Chern number: ', chern_total
   write(unit_chern, '(A,I0)') '# Total Chern (rounded): ', nint(chern_total)
   write(unit_chern, '(A,F15.8)') '# Total deviation from integer: ', abs(chern_total - nint(chern_total))

   close(unit_chern)
   call MIO_Print('TAPW Chern results written successfully','diag')

end subroutine write_chern_tapw_results

!> Comprehensive numerical verification of lattice operations and TAPW setup
subroutine verify_tapw_numerical_consistency(rcell, rG, k_ref, gGridRotationAngle, &
                                            Gx, Gy, NG, aG, path_frac, nPath)
   use cell, only: ucell  ! Import actual direct lattice from cell module
   implicit none
   ! Input parameters
   real(dp), intent(in) :: rcell(3,3)           ! Moiré reciprocal lattice
   real(dp), intent(in) :: rG(2,2)              ! Graphene reciprocal lattice (2x2)
   real(dp), intent(in) :: k_ref(2)             ! Reference K-point
   real(dp), intent(in) :: gGridRotationAngle   ! Rotation angle in degrees
   real(dp), intent(in) :: Gx(:), Gy(:)         ! G-vectors
   integer, intent(in) :: NG                    ! Number of G-vectors
   real(dp), intent(in) :: aG                   ! Graphene lattice constant
   real(dp), intent(in), optional :: path_frac(:,:)  ! Fractional path coordinates
   integer, intent(in), optional :: nPath       ! Number of path points

   ! Local variables
   real(dp) :: rcell_2d(2,2), ucell_T(2,2), ucell_rcell_prod(2,2), identity_dev(2,2)
   real(dp) :: trace_val, off_diag_sum, rotation_matrix(2,2), rot_orthog_dev(2,2)
   real(dp) :: k_ref_frac(2), k_ref_roundtrip(2), path_abs(2), path_roundtrip(2)
   real(dp) :: angle_moire_b1, angle_graphene_b1, relative_angle_diff
   real(dp) :: k_ref_magnitude, expected_k_magnitude, magnitude_ratio
   real(dp) :: rotation_angle, cos_rot, sin_rot, det_rot, det
   real(dp) :: pi = 3.14159265358979323846_dp
   real(dp) :: b1(2), b2(2)
   integer :: i

   if (.not. tapwDebug) return  ! Skip if debug not enabled
   call MIO_Print('=== TAPW NUMERICAL VERIFICATION ===','diag')

   ! DEBUG: Examine full rcell structure
   if (tapwDebug) call MIO_Print('DEBUG: Full rcell matrix structure','diag')
   call MIO_Print('   rcell(1,:) = ['//trim(adjustl(num2str(rcell(1,1),8)))//','//&
                  trim(adjustl(num2str(rcell(1,2),8)))//','//trim(adjustl(num2str(rcell(1,3),8)))//']','diag')
   call MIO_Print('   rcell(2,:) = ['//trim(adjustl(num2str(rcell(2,1),8)))//','//&
                  trim(adjustl(num2str(rcell(2,2),8)))//','//trim(adjustl(num2str(rcell(2,3),8)))//']','diag')
   call MIO_Print('   rcell(3,:) = ['//trim(adjustl(num2str(rcell(3,1),8)))//','//&
                  trim(adjustl(num2str(rcell(3,2),8)))//','//trim(adjustl(num2str(rcell(3,3),8)))//']','diag')

   ! Calculate ucell_T from rcell (same as in generate_shifted_G_list)
   rcell_2d = rcell(1:2,1:2)  ! Extract 2D part
   b1 = rcell_2d(:,1); b2 = rcell_2d(:,2)
   det = b1(1)*b2(2) - b1(2)*b2(1)

   ! Calculate ucell_T from rcell using standard 2D lattice inversion
   ucell_T = reshape([b2(2), -b2(1), -b1(2), b1(1)], [2,2]) * (2.0_dp * pi) / det

   if (tapwDebug) call MIO_Print('DEBUG: Extracted 2D components','diag')
   call MIO_Print('   rcell_2d(1,:) = ['//trim(adjustl(num2str(rcell_2d(1,1),8)))//','//&
                  trim(adjustl(num2str(rcell_2d(1,2),8)))//']','diag')
   call MIO_Print('   rcell_2d(2,:) = ['//trim(adjustl(num2str(rcell_2d(2,1),8)))//','//&
                  trim(adjustl(num2str(rcell_2d(2,2),8)))//']','diag')
   call MIO_Print('   det = '//trim(adjustl(num2str(det,8))),'diag')

   ! 1. Direct–reciprocal duality check: ucell_T^T * rcell_2d should equal 2π * Identity
   ! CRITICAL FIX: We need transpose(ucell_T) * rcell_2d, not ucell_T * rcell_2d
   ucell_rcell_prod = matmul(transpose(ucell_T), rcell_2d)
   identity_dev = ucell_rcell_prod / (2.0_dp * pi)  ! Normalize by 2π to check for identity
   identity_dev(1,1) = identity_dev(1,1) - 1.0_dp  ! Subtract identity
   identity_dev(2,2) = identity_dev(2,2) - 1.0_dp

   trace_val = (ucell_rcell_prod(1,1) + ucell_rcell_prod(2,2)) / (2.0_dp * pi)
   off_diag_sum = (abs(ucell_rcell_prod(1,2)) + abs(ucell_rcell_prod(2,1))) / (2.0_dp * pi)

   call MIO_Print('1. Direct-reciprocal duality: ucell_T^T * rcell_2d = 2π * Identity','diag')
   call MIO_Print('   Normalized trace = '//trim(adjustl(num2str(trace_val,6)))//' (should be 2.0)','diag')
   call MIO_Print('   Normalized off-diag sum = '//trim(adjustl(num2str(off_diag_sum,8)))//' (should be ~0.0)','diag')
   call MIO_Print('   Max identity deviation = '//trim(adjustl(num2str(maxval(abs(identity_dev)),8)))//' (should be ~0.0)','diag')
   if (abs(trace_val - 2.0_dp) < 0.01_dp .and. off_diag_sum < 0.01_dp) then
      call MIO_Print('   ✅ Duality check PASSED: Correct transpose relationship verified','diag')
   else
      call MIO_Print('   ❌ Duality check FAILED: Transpose/frame issue detected','diag')
   end if

   ! 1b. Self-consistency check: rcell should equal 2π * inv(ucell_T)
   call MIO_Print('1b. Self-consistency: rcell vs 2π * inv(ucell_T)','diag')
   ! Calculate inv(ucell_T)
   det_rot = ucell_T(1,1)*ucell_T(2,2) - ucell_T(1,2)*ucell_T(2,1)
   rotation_matrix = reshape([ucell_T(2,2), -ucell_T(1,2), -ucell_T(2,1), ucell_T(1,1)], [2,2]) / det_rot
   rotation_matrix = rotation_matrix * (2.0_dp * pi)
   ! Check deviation
   identity_dev = rcell_2d - rotation_matrix
   call MIO_Print('   Max rcell reconstruction error = '//trim(adjustl(num2str(maxval(abs(identity_dev)),8))),'diag')
   call MIO_Print('   This should be ~0 for consistent lattice relationship','diag')

   ! 2. Rotation orthogonality check (if rotation is applied)
   if (abs(gGridRotationAngle) > 1.0e-10_dp) then
      rotation_angle = gGridRotationAngle * pi / 180.0_dp
      cos_rot = cos(rotation_angle)
      sin_rot = sin(rotation_angle)
      rotation_matrix = reshape([cos_rot, sin_rot, -sin_rot, cos_rot], [2,2])

      rot_orthog_dev = matmul(transpose(rotation_matrix), rotation_matrix)
      rot_orthog_dev(1,1) = rot_orthog_dev(1,1) - 1.0_dp
      rot_orthog_dev(2,2) = rot_orthog_dev(2,2) - 1.0_dp
      det_rot = cos_rot*cos_rot + sin_rot*sin_rot

      call MIO_Print('2. Rotation matrix orthogonality (θ='//trim(adjustl(num2str(gGridRotationAngle,3)))//'°)','diag')
      call MIO_Print('   R^T R - I max deviation = '//trim(adjustl(num2str(maxval(abs(rot_orthog_dev)),8))),'diag')
      call MIO_Print('   det(R) = '//trim(adjustl(num2str(det_rot,8)))//' (should be ~1.0)','diag')
   end if

   ! 3. k_ref in rcell basis (fractional coordinates)
   k_ref_frac = matmul(k_ref, ucell_T) / (2.0_dp * pi)
   call MIO_Print('3. k_ref fractional coordinates in rcell basis','diag')
   call MIO_Print('   k_ref_frac = ['//trim(adjustl(num2str(k_ref_frac(1),6)))//','//&
                  trim(adjustl(num2str(k_ref_frac(2),6)))//']','diag')
   call MIO_Print('   Note: Large values normal for moiré supercell (graphene K in moiré basis)','diag')

   ! 3b. Moiré scale factor analysis
   call MIO_Print('3b. Moiré-graphene scale relationship','diag')
   ! Expected graphene lattice constant
   expected_k_magnitude = 4.0_dp * pi / (3.0_dp * aG)
   ! Moiré lattice vector magnitudes
   magnitude_ratio = sqrt(rcell_2d(1,1)**2 + rcell_2d(2,1)**2)  ! |b1|
   det_rot = sqrt(rcell_2d(1,2)**2 + rcell_2d(2,2)**2)         ! |b2|
   call MIO_Print('   Moiré |b1| = '//trim(adjustl(num2str(magnitude_ratio,6))),'diag')
   call MIO_Print('   Moiré |b2| = '//trim(adjustl(num2str(det_rot,6))),'diag')
   call MIO_Print('   Expected graphene |b| = 4π/(3a) = '//trim(adjustl(num2str(expected_k_magnitude,6))),'diag')
   call MIO_Print('   Scale factors: |b1|/|b_graphene| = '//trim(adjustl(num2str(magnitude_ratio/expected_k_magnitude,3)))//&
                  ', |b2|/|b_graphene| = '//trim(adjustl(num2str(det_rot/expected_k_magnitude,3))),'diag')

   ! 4. Nearest-G rounding sensitivity
   call MIO_Print('4. Nearest-G rounding (Voronoi boundary check)','diag')
   call MIO_Print('   k_ref_frac - nint(k_ref_frac) = ['//&
                  trim(adjustl(num2str(k_ref_frac(1) - nint(k_ref_frac(1)),6)))//','//&
                  trim(adjustl(num2str(k_ref_frac(2) - nint(k_ref_frac(2)),6)))//']','diag')
   call MIO_Print('   Values near ±0.5 indicate Voronoi boundary sensitivity','diag')

   ! Check if we're dangerously close to boundary
   if (abs(abs(k_ref_frac(1) - nint(k_ref_frac(1))) - 0.5_dp) < 0.01_dp .or. &
       abs(abs(k_ref_frac(2) - nint(k_ref_frac(2))) - 0.5_dp) < 0.01_dp) then
      call MIO_Print('   ⚠️  WARNING: Very close to Voronoi boundary! G-vector selection may be unstable.','diag')
      call MIO_Print('   Consider adjusting k_ref slightly or using higher precision.','diag')
   end if

   ! 5. Path round-trip test (if path provided)
   if (present(path_frac) .and. present(nPath) .and. nPath > 0) then
      call MIO_Print('5. Path coordinate round-trip test (first point)','diag')
      ! frac → abs → frac
      path_abs = path_frac(1,1)*rcell_2d(:,1) + path_frac(2,1)*rcell_2d(:,2)
      path_roundtrip = matmul(path_abs, ucell_T) / (2.0_dp * pi)
      call MIO_Print('   Original frac = ['//trim(adjustl(num2str(path_frac(1,1),6)))//','//&
                     trim(adjustl(num2str(path_frac(2,1),6)))//']','diag')
      call MIO_Print('   Round-trip frac = ['//trim(adjustl(num2str(path_roundtrip(1),6)))//','//&
                     trim(adjustl(num2str(path_roundtrip(2),6)))//']','diag')
      call MIO_Print('   Difference = ['//trim(adjustl(num2str(path_frac(1,1)-path_roundtrip(1),8)))//','//&
                     trim(adjustl(num2str(path_frac(2,1)-path_roundtrip(2),8)))//']','diag')
   end if

   ! 6. G-grid frame verification
   call MIO_Print('6. G-grid coordinate frame','diag')
   call MIO_Print('   rcell b1 = ['//trim(adjustl(num2str(rcell(1,1),6)))//','//&
                  trim(adjustl(num2str(rcell(2,1),6)))//']','diag')
   call MIO_Print('   rcell b2 = ['//trim(adjustl(num2str(rcell(1,2),6)))//','//&
                  trim(adjustl(num2str(rcell(2,2),6)))//']','diag')
   call MIO_Print('   gGridRotationAngle = '//trim(adjustl(num2str(gGridRotationAngle,3)))//'°','diag')
   call MIO_Print('   skipGRotation = '//merge('T','F',skipGRotation),'diag')

   ! 7. First few G-vectors
   call MIO_Print('7. First 5 G-vectors and norms','diag')
   call MIO_Print('   Compare with Python: should match G-vector order exactly','diag')
   do i = 1, min(5, NG)
      call MIO_Print('   G('//trim(adjustl(num2str(real(i,dp),0)))//') = ['//&
                     trim(adjustl(num2str(Gx(i),6)))//','//trim(adjustl(num2str(Gy(i),6)))//&
                     '], |G| = '//trim(adjustl(num2str(sqrt(Gx(i)**2 + Gy(i)**2),6))),'diag')
   end do

   ! 7b. G-vector distance ordering verification
   call MIO_Print('7b. G-vector distance ordering (should match Python np.unique)','diag')
   do i = 1, min(5, NG)
      call MIO_Print('   Distance from k_ref: |G('//trim(adjustl(num2str(real(i,dp),0)))//')| = '//&
                     trim(adjustl(num2str(sqrt((Gx(i)-k_ref(1))**2 + (Gy(i)-k_ref(2))**2),6))),'diag')
   end do

   ! 8. k_ref magnitude sanity check
   k_ref_magnitude = sqrt(k_ref(1)**2 + k_ref(2)**2)
   expected_k_magnitude = 4.0_dp * pi / (3.0_dp * aG)  ! |K| for graphene
   magnitude_ratio = k_ref_magnitude / expected_k_magnitude
   call MIO_Print('8. k_ref magnitude sanity','diag')
   call MIO_Print('   |k_ref| = '//trim(adjustl(num2str(k_ref_magnitude,6))),'diag')
   call MIO_Print('   Expected |K| = 4π/(3aG) = '//trim(adjustl(num2str(expected_k_magnitude,6))),'diag')
   call MIO_Print('   Ratio = '//trim(adjustl(num2str(magnitude_ratio,6)))//' (should be ~1.0 for K-point)','diag')

   ! 9. Relative orientation between moiré and graphene lattices
   angle_moire_b1 = atan2(rcell(2,1), rcell(1,1)) * 180.0_dp / pi
   angle_graphene_b1 = atan2(rG(2,1), rG(1,1)) * 180.0_dp / pi
   relative_angle_diff = angle_moire_b1 - angle_graphene_b1
   ! Handle angle wrapping
   if (relative_angle_diff > 180.0_dp) relative_angle_diff = relative_angle_diff - 360.0_dp
   if (relative_angle_diff < -180.0_dp) relative_angle_diff = relative_angle_diff + 360.0_dp

   call MIO_Print('9. Relative lattice orientation','diag')
   call MIO_Print('   Moiré b1 angle = '//trim(adjustl(num2str(angle_moire_b1,3)))//'°','diag')
   call MIO_Print('   Graphene b1 angle = '//trim(adjustl(num2str(angle_graphene_b1,3)))//'°','diag')
   call MIO_Print('   Relative rotation = '//trim(adjustl(num2str(relative_angle_diff,3)))//'°','diag')

   call MIO_Print('=== END VERIFICATION ===','diag')

end subroutine verify_tapw_numerical_consistency

!> Check X-matrix column norms for TAPW projection
subroutine verify_x_matrix_norms(XArray, N, NG)
   implicit none
   complex(dp), intent(in) :: XArray(:,:)
   integer, intent(in) :: N, NG
   real(dp) :: col_norms(5), max_norm, min_norm, avg_norm
   integer :: i, ncols_to_check

   ncols_to_check = min(5, NG)
   if (.not. tapwDebug) return  ! Skip if debug not enabled
   call MIO_Print('=== X-MATRIX VERIFICATION ===','diag')

   ! Check first few column 2-norms
   do i = 1, ncols_to_check
      col_norms(i) = sqrt(real(sum(conjg(XArray(:,i)) * XArray(:,i))))
   end do

   max_norm = maxval(col_norms(1:ncols_to_check))
   min_norm = minval(col_norms(1:ncols_to_check))
   avg_norm = sum(col_norms(1:ncols_to_check)) / real(ncols_to_check)

   call MIO_Print('X-matrix column 2-norms (first '//trim(adjustl(num2str(real(ncols_to_check,dp),0)))//' columns):','diag')
   do i = 1, ncols_to_check
      call MIO_Print('  ||X(:,'//trim(adjustl(num2str(real(i,dp),0)))//')|| = '//trim(adjustl(num2str(col_norms(i),6))),'diag')
   end do
   call MIO_Print('  Range: ['//trim(adjustl(num2str(min_norm,6)))//','//trim(adjustl(num2str(max_norm,6)))//&
                  '], Avg = '//trim(adjustl(num2str(avg_norm,6))),'diag')
   call MIO_Print('  Expected: ~1.0 for proper normalization','diag')

end subroutine verify_x_matrix_norms

!> Check Hermiticity of projected Hamiltonian
subroutine verify_projected_hamiltonian_hermiticity(Hproj, M)
   implicit none
   complex(dp), intent(in) :: Hproj(:,:)
   integer, intent(in) :: M
   real(dp) :: max_hermitian_dev, avg_hermitian_dev
   complex(dp) :: hermitian_dev
   integer :: i, j, count

   if (.not. tapwDebug) return  ! Skip if debug not enabled
   call MIO_Print('=== PROJECTED HAMILTONIAN HERMITICITY ===','diag')

   max_hermitian_dev = 0.0_dp
   avg_hermitian_dev = 0.0_dp
   count = 0

   do i = 1, M
      do j = 1, M
         hermitian_dev = Hproj(i,j) - conjg(Hproj(j,i))
         max_hermitian_dev = max(max_hermitian_dev, abs(hermitian_dev))
         avg_hermitian_dev = avg_hermitian_dev + abs(hermitian_dev)
         count = count + 1
      end do
   end do

   if (count > 0) then
      avg_hermitian_dev = avg_hermitian_dev / real(count)
   else
      avg_hermitian_dev = 0.0_dp
   end if

   call MIO_Print('Projected H Hermiticity check:','diag')
   call MIO_Print('  Matrix size M = '//trim(adjustl(num2str(real(M,dp),0)))//', elements checked = '//trim(adjustl(num2str(real(count,dp),0))),'diag')
   call MIO_Print('  Max |H(i,j) - H*(j,i)| = '//trim(adjustl(num2str(max_hermitian_dev,10))),'diag')
   call MIO_Print('  Avg |H(i,j) - H*(j,i)| = '//trim(adjustl(num2str(avg_hermitian_dev,10))),'diag')
   call MIO_Print('  Expected: ~0.0 for Hermitian matrix','diag')

end subroutine verify_projected_hamiltonian_hermiticity

subroutine build_tapw_labels(layerIndex, species, N, label, Nlabel)
  use iso_fortran_env, only: dp => real64
  implicit none
  integer,           intent(in)  :: N
  integer,           intent(in)  :: layerIndex(N), species(N)
  integer,           intent(out) :: label(N)      ! remapped contiguous labels 1..Nlabel
  integer,           intent(out) :: Nlabel        ! number of unique (layer,sublattice) groups

  integer :: i, j
  integer, allocatable :: raw_label(:)            ! two-digit codes: 10*layer + species
  integer, allocatable :: uniq(:)                 ! unique raw codes (small array)
  integer :: nuniq, code, pos

  ! --- Build raw two-digit codes (e.g., 11,12,21,22) ---
  allocate(raw_label(N))
  do i = 1, N
     if (layerIndex(i) < 0 .or. species(i) < 0) then
        error stop "build_tapw_labels: negative layer/species not allowed"
     end if
     raw_label(i) = 10*layerIndex(i) + species(i)
  end do

  ! --- Collect unique codes in first-seen order (Nlabel is small in practice) ---
  allocate(uniq(N))
  nuniq = 0
  do i = 1, N
     code = raw_label(i)
     ! linear search over current uniques
     pos = 0
     do j = 1, nuniq
        if (code == uniq(j)) then
           pos = j; exit
        end if
     end do
     if (pos == 0) then
        nuniq = nuniq + 1
        uniq(nuniq) = code
     end if
  end do

  ! --- Sort uniques ascending to make the mapping deterministic ---
  call sort_int_ascending(uniq, nuniq)

  ! --- Remap raw_label -> contiguous 1..Nlabel following sorted uniq() ---
  Nlabel = nuniq
  do i = 1, N
     code = raw_label(i)
     pos  = 0
     do j = 1, Nlabel
        if (code == uniq(j)) then
           pos = j; exit
        end if
     end do
     if (pos == 0) stop "build_tapw_labels: internal mapping error"
     label(i) = pos
  end do

  ! (Optional) debug print
  print *, "TAPW labels (sorted):"
  do j = 1, Nlabel
     print '(A,I0,A,I0)', "  code ", uniq(j), " -> contiguous ", j
  end do

  deallocate(raw_label, uniq)
end subroutine build_tapw_labels

subroutine sort_int_ascending(a, n)
  implicit none
  integer, intent(inout) :: a(:)
  integer, intent(in)    :: n
  integer :: i, j, key
  do i = 2, n
     key = a(i); j = i - 1
     do while (j >= 1 .and. a(j) > key)
        a(j+1) = a(j); j = j - 1
     end do
     a(j+1) = key
  end do
end subroutine sort_int_ascending

! Transform Hamiltonian using dense matrix approach for TAPW
! This constructs H directly as a dense matrix and computes Hproj = X† H X
! IMPORTANT: Uses actual atomic position differences (NeighD) instead of lattice vectors
! to match the paper formulation: R_ij = (τ_j - τ_i) + T_ij
subroutine transform_dense_hamiltonian_tapw(N, M, X, Hproj, KLoc, cell, H0, maxN, hopp, NList, Nneigh, neighCell, ns, is)
    use constants, only : cmplx_i, cmplx_0, cmplx_1
    use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
    use scf, only : charge, Zch
    use atoms, only : Rat, Species, layerIndex
    use tbpar, only : U
    use neigh, only : NeighD
    use ham, only : Zterm, gZeeman, IntrinsicSOCterm, lambdaI, IsingSOCterm, lambdaIsing, PIASOCterm, lambdaPIA, SOCEnabledForLayer
    use magf, only : BmagZeeman
    implicit none

    ! Input parameters
    integer, intent(in) :: N, M, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), ns, is
    complex(dp), intent(in) :: X(N,M)
    complex(dp), intent(out) :: Hproj(M,M)
    real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
    complex(dp), intent(in) :: hopp(maxN,N)

    ! Local variables
    complex(dp), allocatable :: H_dense(:,:), temp(:,:)
    integer :: i, j, in
    real(dp) :: R(3)
    real(dp) :: soc_diagonal_contrib
    complex(dp) :: phase
    ! Debug variables for hermiticity check
    real(dp) :: max_debug_error, error
    integer :: max_i, max_j

    if (tapwDebug) call MIO_Print('Building dense Hamiltonian matrix for TAPW transformation','diag')

    if (tapwDebug) then
       ! DEBUG: Check atomic positions used in dense TB construction (should be relaxed)
       if (tapwDebug) call MIO_Print('DEBUG: First 3 atomic positions IN dense TB construction:','diag')
       do i = 1, min(3, N)
           call MIO_Print('  Atom '//trim(num2str(i))//': ['//trim(num2str(Rat(1,i),6))//','//&
                         trim(num2str(Rat(2,i),6))//','//trim(num2str(Rat(3,i),6))//']','diag')
       end do
    end if

    ! For very large systems, warn about memory usage but continue with dense approach
    if (N > 1000000) then
        call MIO_Print('Large system detected (N='//trim(num2str(N))//'), using dense TAPW transformation','diag')
        call MIO_Print('WARNING: This will use significant memory for large systems','diag')
        ! Continue with dense matrix approach
    end if

    ! Allocate dense Hamiltonian matrix
    allocate(H_dense(N, N), temp(N, M))
    H_dense = cmplx_0

    ! Build dense Hamiltonian matrix using the same logic as sparse matrix
    ! Diagonal terms with SOC application
    ! Note: SOC terms may be non-zero even if H0(i) == 0, so we need to calculate SOC contributions first
    do i = 1, N
        soc_diagonal_contrib = 0.0_dp

        ! Calculate SOC contributions to diagonal element
        ! Note: This is backward compatible - only applies when SOC is enabled
        if (Zterm .or. IntrinsicSOCterm .or. IsingSOCterm .or. PIASOCterm) then
           ! Apply Zeeman effect (λVZ) - spin-dependent onsite
           if (Zterm .and. (ns==2)) then
              if (is==1) then
                 soc_diagonal_contrib = soc_diagonal_contrib + gZeeman*BmagZeeman  ! Spin-up: +λVZ
              else
                 soc_diagonal_contrib = soc_diagonal_contrib - gZeeman*BmagZeeman  ! Spin-down: -λVZ
              end if
           end if

           ! Apply Intrinsic SOC (λI) - gap-opening onsite term
           ! For graphene, this should be sublattice-dependent to open gap
           if (IntrinsicSOCterm .and. SOCEnabledForLayer(layerIndex(i))) then
              ! Apply +λI to sublattice A, -λI to sublattice B (or vice versa)
              ! This opens a gap of size 2λI at Dirac point
              ! Species(i)=1 is A sublattice, Species(i)=2 is B sublattice
              if (Species(i) == 1) then
                 soc_diagonal_contrib = soc_diagonal_contrib + lambdaI
              else
                 soc_diagonal_contrib = soc_diagonal_contrib - lambdaI
              end if
           end if

           ! Apply Ising SOC (λIsing) - Valley-Zeeman term (τ_z s_z)
           ! This is a valley-dependent spin splitting: τ_z = +1 for K valley, -1 for K' valley
           ! Does not open a global gap, only splits spins differently in K vs K' valleys
           if (IsingSOCterm .and. SOCEnabledForLayer(layerIndex(i)) .and. (ns==2)) then
              ! Valley sign: +1 for K valley, -1 for K' valley
              if (.not. useKprimeValley) then
                 ! K valley: τ_z = +1
                 if (is == 1) then
                    ! Spin-up: +λ for all sites
                    soc_diagonal_contrib = soc_diagonal_contrib + lambdaIsing
                 else
                    ! Spin-down: -λ for all sites
                    soc_diagonal_contrib = soc_diagonal_contrib - lambdaIsing
                 end if
              else
                 ! K' valley: τ_z = -1 (flip sign)
                 if (is == 1) then
                    ! Spin-up: -λ for all sites
                    soc_diagonal_contrib = soc_diagonal_contrib - lambdaIsing
                 else
                    ! Spin-down: +λ for all sites
                    soc_diagonal_contrib = soc_diagonal_contrib + lambdaIsing
                 end if
              end if
           end if

           ! Apply Pseudo-inversion asymmetry (λPIA) - onsite component
           if (PIASOCterm .and. SOCEnabledForLayer(layerIndex(i))) then
              ! PIA is spin-independent onsite term
              soc_diagonal_contrib = soc_diagonal_contrib + lambdaPIA
           end if
        end if

        ! Add SCF terms if spin-polarized (matches DiagHam implementation)
        ! Commented out: not doing any SCF calculation for now (matches BuildBlockHamiltonianOnly)

        ! Add diagonal element if H0(i) is non-zero OR if SOC/SCF contributions are non-zero
        ! Note: Dense path matches DiagHam exactly (no sigma - sigma is only for sparse shift-and-invert)
        if (H0(i) /= 0.0_dp .or. soc_diagonal_contrib /= 0.0_dp) then
            H_dense(i, i) = H0(i) + soc_diagonal_contrib
        end if
    end do

    ! Off-diagonal hopping terms
    ! Use actual atomic position differences as per paper: R_ij = (τ_j - τ_i) + T_ij
    ! We want distance from i to j, so we use -NeighD which gives (τ_i - τ_j) + T_ji
    do i = 1, N
        do j = 1, Nneigh(i)
            in = NList(j, i)
            if (in > 0 .and. in <= N) then
                ! Use actual atomic position difference instead of lattice vector
                ! NeighD(:,j,i) contains the vector from atom i to atom j: (τ_j - τ_i) + T_ij
                ! For consistency with paper formulation, use -NeighD to get distance from i to j
                R(1:3) = 0.0_dp
                R(1:2) = -NeighD(1:2, j, i)  ! Use negative to get distance from i to j
                ! Note: DiagHam sets HLoc(in, i), but we're building H_dense(i, in) to match dense storage
                ! For consistency, we'll set H_dense(in, i) first, then Hermitian conjugate
                H_dense(in, i) = H_dense(in, i) - hopp(j, i) * exp(-cmplx_i * dot_product(KLoc, R))
                ! Apply PIA hopping terms if enabled (matches DiagHam implementation)
                ! PIA hopping not yet properly implemented - commented out
                !   ! PIA hopping: ApplyPIAHopping sets HLoc(in, i) = HLoc(in, i) + phase
                !   ! where phase = lambdaPIA * cmplx(neighD(2, j, i), neighD(1, j, i))
                !   ! Then it sets HLoc(i, in) = HLoc(i, in) - conjg(phase) for Hermiticity
            end if
        end do
    end do

    ! Add edge hopping terms if present
    ! Note: Edge terms still use lattice vector approach as NeighD may not be available for edge atoms
    if (nQ > 0) then
        do i = 1, nQ
            do j = 1, nEdgeN(i)
                in = NeI(j, i)
                if (in > 0 .and. in <= N .and. edgeIndx(i) > 0 .and. edgeIndx(i) <= N) then
                    R = matmul(cell, NedgeCell(:, j, i))
                    H_dense(edgeIndx(i), in) = H_dense(edgeIndx(i), in) + edgeH(j, i) * exp(-cmplx_i * dot_product(KLoc, R))
                end if
            end do
        end do
    end if

    if (tapwDebug) call MIO_Print('Dense Hamiltonian matrix constructed, performing TAPW transformation','diag')
    if (tapwDebug) call MIO_Print('  Matrix dimensions: H('//trim(num2str(N))//'×'//trim(num2str(N))//'), X('//trim(num2str(N))//'×'//trim(num2str(M))//')','diag')

    ! Check hermiticity of dense Hamiltonian (only if tapwDebug enabled)
    if (tapwDebug) then
       call MIO_Print('Checking hermiticity of dense TB Hamiltonian H_dense','diag')
       call check_matrix_hermiticity(H_dense, N, 'H_dense')
    end if

    ! DEBUG: Find the actual hermiticity violations (only if tapwDebug enabled)
    if (tapwDebug .and. N >= 2) then
        call MIO_Print('DEBUG: Finding hermiticity violations:','diag')
        max_debug_error = 0.0_dp
        max_i = 0
        max_j = 0

        do i = 1, N
            do j = i+1, N
                error = abs(H_dense(i,j) - conjg(H_dense(j,i)))
                if (error > max_debug_error) then
                    max_debug_error = error
                    max_i = i
                    max_j = j
                end if
            end do
        end do

        if (max_debug_error > 1e-12) then
            call MIO_Print('  Maximum violation at H('//trim(num2str(max_i))//','//trim(num2str(max_j))//'):','diag')
            call MIO_Print('  H('//trim(num2str(max_i))//','//trim(num2str(max_j))//') = '//&
                         trim(num2str(real(H_dense(max_i,max_j)),8))//' + i*'//trim(num2str(aimag(H_dense(max_i,max_j)),8)),'diag')
            call MIO_Print('  H('//trim(num2str(max_j))//','//trim(num2str(max_i))//') = '//&
                         trim(num2str(real(H_dense(max_j,max_i)),8))//' + i*'//trim(num2str(aimag(H_dense(max_j,max_i)),8)),'diag')
            call MIO_Print('  conjg(H('//trim(num2str(max_j))//','//trim(num2str(max_i))//')) = '//&
                         trim(num2str(real(conjg(H_dense(max_j,max_i))),8))//' + i*'//trim(num2str(aimag(conjg(H_dense(max_j,max_i))),8)),'diag')
            call MIO_Print('  Error = '//trim(num2str(max_debug_error,8)),'diag')

            ! DEBUG: Check if there are duplicate entries in neighbor list for problematic atoms
            call MIO_Print('DEBUG: Checking neighbor list for atoms '//trim(num2str(max_i))//' and '//trim(num2str(max_j))//':','diag')

            ! Check if max_i has max_j as neighbor
            do j = 1, Nneigh(max_i)
                if (NList(j, max_i) == max_j) then
                    call MIO_Print('  Atom '//trim(num2str(max_i))//' has atom '//trim(num2str(max_j))//' as neighbor '//trim(num2str(j)),'diag')
                    call MIO_Print('  NeighD = ['//trim(num2str(NeighD(1,j,max_i),6))//','//trim(num2str(NeighD(2,j,max_i),6))//']','diag')
                    call MIO_Print('  hopp = '//trim(num2str(real(hopp(j,max_i)),8))//' + i*'//trim(num2str(aimag(hopp(j,max_i)),8)),'diag')
                end if
            end do

            ! Check if max_j has max_i as neighbor
            do j = 1, Nneigh(max_j)
                if (NList(j, max_j) == max_i) then
                    call MIO_Print('  Atom '//trim(num2str(max_j))//' has atom '//trim(num2str(max_i))//' as neighbor '//trim(num2str(j)),'diag')
                    call MIO_Print('  NeighD = ['//trim(num2str(NeighD(1,j,max_j),6))//','//trim(num2str(NeighD(2,j,max_j),6))//']','diag')
                    call MIO_Print('  hopp = '//trim(num2str(real(hopp(j,max_j)),8))//' + i*'//trim(num2str(aimag(hopp(j,max_j)),8)),'diag')
                end if
            end do
        else
            call MIO_Print('  No significant hermiticity violations found','diag')
        end if
    end if

    ! Compute Hproj = X† H X using BLAS
    ! Step 1: temp = H * X
    call zgemm('N', 'N', N, M, N, cmplx_1, H_dense, N, X, N, cmplx_0, temp, N)

    ! Step 2: Hproj = X† * temp
    call zgemm('C', 'N', M, M, N, cmplx_1, X, N, temp, N, cmplx_0, Hproj, M)

    ! Check hermiticity of projected Hamiltonian after TAPW transformation (only if tapwDebug enabled)
    if (tapwDebug) then
       call MIO_Print('Checking hermiticity of projected Hamiltonian Hproj','diag')
       call check_matrix_hermiticity(Hproj, M, 'Hproj')
    end if

    if (tapwDebug) call MIO_Print('Dense TAPW transformation completed','diag')
    if (tapwDebug) call MIO_Print('  Projected matrix dimensions: Hproj('//trim(num2str(M))//'×'//trim(num2str(M))//')','diag')

    deallocate(H_dense, temp)
end subroutine transform_dense_hamiltonian_tapw

! Read rigid reference positions from generateInit.xyz for TAPW X matrix construction
! This ensures perfect TAPW unitarity even with lattice reconstruction
subroutine read_rigid_positions_for_tapw(filename, N, rigid_positions)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: N
    real(dp), intent(out) :: rigid_positions(3, N)

    ! Local variables
    integer :: unit, i, N_file
    real(dp) :: cell_vectors(3,3)
    character(len=256) :: line, element
    real(dp) :: x, y, z

    call MIO_Print('Reading rigid reference positions from: '//trim(filename),'diag')

    ! Open the rigid position file
    unit = 97
    open(unit=unit, file=trim(filename), status='old', action='read', iostat=i)
    if (i /= 0) then
        call MIO_Print('ERROR: Cannot open rigid position file: '//trim(filename),'diag')
        call MIO_Print('Make sure generateInit.xyz exists in the working directory','diag')
        error stop 'Failed to open rigid position file'
    end if

    ! Read cell vectors (3 lines) - we don't need them but must skip them
    do i = 1, 3
        read(unit, '(A)') line
    end do

    ! Read number of atoms
    read(unit, *) N_file
    if (N_file /= N) then
        call MIO_Print('ERROR: Atom count mismatch!','diag')
        call MIO_Print('  Current system has '//trim(num2str(N))//' atoms','diag')
        call MIO_Print('  Rigid file has '//trim(num2str(N_file))//' atoms','diag')
        error stop 'Atom count mismatch between current system and rigid reference'
    end if

    ! Read atomic positions (columns 2, 3, 4 are x, y, z)
    do i = 1, N
        read(unit, *) element, x, y, z
        rigid_positions(1, i) = x
        rigid_positions(2, i) = y
        rigid_positions(3, i) = z
    end do

    close(unit)

    call MIO_Print('Successfully read '//trim(num2str(N))//' rigid reference positions','diag')
    call MIO_Print('These will be used for X matrix construction to ensure TAPW unitarity','diag')

end subroutine read_rigid_positions_for_tapw

subroutine GetBandEnergyRange(band_index, min_energy, max_energy, spin_index)
   ! Extract eigenvalues for a specific band across all k-points from stored Hamiltonians
   ! This provides accurate energy ranges instead of relying on potentially incomplete E array
   ! Updated to include spin_index for separate spin channel calculations
   implicit none

   ! Input/output parameters
   integer, intent(in) :: band_index
   real(dp), intent(out) :: min_energy, max_energy
   integer, intent(in), optional :: spin_index  ! Optional, defaults to 1 for backwards compatibility

   ! Local spin index
   integer :: is_local

   ! Local variables
   complex(dp), allocatable :: H_temp(:,:)
   real(dp), allocatable :: eigval_temp(:)
   complex(dp), allocatable :: work_temp(:)
   real(dp), allocatable :: rwork_temp(:)
   real(dp), allocatable :: band_energies(:)
   integer :: ik, info, lwork, M_safe

   ! Check if stored data is available
   if (.not. allocated(stored_hamiltonians) .or. stored_nk == 0) then
      call MIO_Print('WARNING: No stored Hamiltonian data for energy range calculation','diag')
      min_energy = 0.0_dp
      max_energy = 0.0_dp
      return
   end if

   ! Check band index bounds
   if (band_index < 1 .or. band_index > stored_M) then
      call MIO_Print('WARNING: Band index '//trim(num2str(band_index))//' out of range [1,'//trim(num2str(stored_M))//']','diag')
      min_energy = 0.0_dp
      max_energy = 0.0_dp
      return
   end if

   ! Use safe dimensions
   M_safe = min(stored_M, size(stored_hamiltonians,1), size(stored_hamiltonians,2))

   ! Allocate arrays for eigenvalue extraction
   allocate(H_temp(M_safe,M_safe))
   allocate(eigval_temp(M_safe))
   allocate(band_energies(stored_nk))

   ! Set up workspace for ZHEEV
   lwork = 2*M_safe
   allocate(work_temp(lwork))
   allocate(rwork_temp(3*M_safe-2))

   ! Determine spin index (default to 1 for backwards compatibility)
   is_local = 1
   if (present(spin_index)) then
      is_local = spin_index
   end if

   ! Check spin index bounds
   if (is_local < 1 .or. is_local > stored_nspin) then
      call MIO_Print('WARNING: spin_index '//trim(num2str(is_local))//' out of range [1,'//trim(num2str(stored_nspin))//'], using 1','diag')
      is_local = 1
   end if

   ! Extract eigenvalues for this band from each k-point
   do ik = 1, stored_nk
      ! Copy Hamiltonian for this k-point and spin channel
      H_temp(1:M_safe,1:M_safe) = stored_hamiltonians(1:M_safe,1:M_safe,ik,is_local)

      ! Diagonalize to get eigenvalues
      call ZHEEV('N', 'U', M_safe, H_temp, M_safe, eigval_temp, work_temp, lwork, rwork_temp, info)

      if (info /= 0) then
         call MIO_Print('WARNING: ZHEEV failed for k-point '//trim(num2str(ik))//' in GetBandEnergyRange','diag')
         band_energies(ik) = 0.0_dp
      else
         ! Store eigenvalue for this band at this k-point
         if (band_index <= M_safe) then
            band_energies(ik) = eigval_temp(band_index)
         else
            band_energies(ik) = 0.0_dp
         end if
      end if
   end do

   ! Calculate min and max energies across all k-points
   min_energy = minval(band_energies)
   max_energy = maxval(band_energies)

   ! Clean up
   deallocate(H_temp, eigval_temp, work_temp, rwork_temp, band_energies)

end subroutine GetBandEnergyRange

subroutine check_matrix_hermiticity(matrix, n, matrix_name)
   ! Check if a complex matrix is hermitian
   use constants, only : cmplx_0
   implicit none

   ! Input parameters
   integer, intent(in) :: n
   complex(dp), intent(in) :: matrix(n,n)
   character(len=*), intent(in) :: matrix_name

   ! Local variables
   integer :: i, j
   real(dp) :: max_error, error
   complex(dp) :: diff

   max_error = 0.0_dp

   do i = 1, n
      do j = 1, n
         ! Check if matrix(i,j) = conjg(matrix(j,i))
         diff = matrix(i,j) - conjg(matrix(j,i))
         error = abs(diff)
         if (error > max_error) then
            max_error = error
         end if
      end do
   end do

   call MIO_Print('Hermiticity check for '//trim(matrix_name)//': max error = '//trim(num2str(max_error,8)),'diag')

   if (max_error > 1.0e-10_dp) then
      call MIO_Print('WARNING: '//trim(matrix_name)//' is NOT hermitian! Max error = '//trim(num2str(max_error,8)),'diag')
   else
      call MIO_Print(trim(matrix_name)//' is hermitian (error < 1e-10)','diag')
   end if

end subroutine check_matrix_hermiticity

subroutine verify_eigenvector_eigenvalue_alignment(H, eigvec, eigval, M, kpoint_index)
   ! Verify that eigenvectors and eigenvalues are properly aligned from the same diagonalization
   use constants, only : cmplx_0, cmplx_1
   implicit none

   ! Input parameters
   integer, intent(in) :: M, kpoint_index
   complex(dp), intent(in) :: H(M,M), eigvec(M,M)
   real(dp), intent(in) :: eigval(M)

   ! Local variables
   integer :: i, j, k
   complex(dp) :: residual, max_residual, Hv_j
   real(dp) :: tolerance
   logical :: alignment_ok

   tolerance = 1.0e-6_dp  ! More realistic tolerance for numerical precision
   max_residual = cmplx_0
   alignment_ok = .true.

   ! Check that H * eigvec(:,i) = eigval(i) * eigvec(:,i) for each eigenvector
   do i = 1, M
      ! Compute residual: |H * v_i - λ_i * v_i|
      residual = cmplx_0
      do j = 1, M
         ! CORRECTED: H * v_i means sum(H(j,:) * eigvec(:,i)), not H(j,i)
         ! Compute (H * v_i)_j = sum_k H(j,k) * v_i(k)
         Hv_j = cmplx_0
         do k = 1, M
            Hv_j = Hv_j + H(j,k) * eigvec(k,i)
         end do
         ! Now check: (H * v_i)_j - λ_i * v_i(j)
         residual = residual + (Hv_j - eigval(i) * eigvec(j,i)) * conjg(Hv_j - eigval(i) * eigvec(j,i))
      end do
      residual = sqrt(residual)

      if (abs(residual) > abs(max_residual)) then
         max_residual = residual
      end if

      if (abs(residual) > tolerance) then
         alignment_ok = .false.
         if (i <= 5) then  ! Only print first few violations to avoid spam
            call MIO_Print('  MISALIGNMENT: Band '//trim(num2str(i))//' residual = '//trim(num2str(abs(residual),8)),'diag')
         end if
      end if
   end do

   ! Report results
   if (alignment_ok) then
      call MIO_Print('  ✓ Eigenvector-eigenvalue alignment verified (max residual = '//trim(num2str(abs(max_residual),8))//')','diag')
   else
      call MIO_Print('  ✗ EIGENVECTOR-EIGENVALUE MISALIGNMENT DETECTED!','diag')
      call MIO_Print('  ✗ Maximum residual = '//trim(num2str(abs(max_residual),8))//' (tolerance = '//trim(num2str(tolerance,8))//')','diag')
      call MIO_Print('  ✗ This indicates eigenvectors and eigenvalues are from different diagonalizations!','diag')

      ! Additional diagnostic for first few bands
      if (kpoint_index == 1) then
         call MIO_Print('  DEBUG: First 3 eigenvalues and their residuals:','diag')
         do i = 1, min(3, M)
            residual = cmplx_0
            do j = 1, M
               residual = residual + (H(j,i) - eigval(i) * eigvec(j,i)) * conjg(H(j,i) - eigval(i) * eigvec(j,i))
            end do
            residual = sqrt(residual)
            call MIO_Print('    Band '//trim(num2str(i))//': λ = '//trim(num2str(eigval(i),6))//', residual = '//trim(num2str(abs(residual),8)),'diag')
         end do
      end if
   end if

end subroutine verify_eigenvector_eigenvalue_alignment

subroutine add_high_symmetry_refinement(Kpts, ip, rcell, nk_x, nk_y)
   ! Add extra k-points near high-symmetry points for better Berry curvature sampling
   use constants, only : cmplx_0
   implicit none

   ! Input/Output parameters
   real(dp), intent(inout) :: Kpts(3, *)
   integer, intent(inout) :: ip
   real(dp), intent(in) :: rcell(3,3)
   integer, intent(in) :: nk_x, nk_y

   ! Local variables
   integer :: i, j, n_refine
   real(dp) :: b1(3), b2(3), k_refine(3)
   real(dp) :: refinement_factor, kx_step, ky_step

   ! High-symmetry points in fractional coordinates
   ! Γ = (0, 0), K = (1/3, 1/3), M = (1/2, 0), K' = (2/3, 1/3)
   real(dp), parameter :: gamma(2) = [0.0_dp, 0.0_dp]
   real(dp), parameter :: K_point(2) = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
   real(dp), parameter :: M_point(2) = [0.5_dp, 0.0_dp]
   real(dp), parameter :: K_prime(2) = [2.0_dp/3.0_dp, 1.0_dp/3.0_dp]

   b1 = rcell(:,1)
   b2 = rcell(:,2)

   ! Refinement parameters
   refinement_factor = 4.0_dp  ! 4x finer grid near high-symmetry points
   n_refine = 3  ! 3x3 grid around each high-symmetry point

   kx_step = 1.0_dp / (real(nk_x, dp) * refinement_factor)
   ky_step = 1.0_dp / (real(nk_y, dp) * refinement_factor)

   call MIO_Print('  Adding refinement around Γ point','diag')
   call add_refinement_around_point(Kpts, ip, gamma, b1, b2, kx_step, ky_step, n_refine)

   call MIO_Print('  Adding refinement around K point','diag')
   call add_refinement_around_point(Kpts, ip, K_point, b1, b2, kx_step, ky_step, n_refine)

   call MIO_Print('  Adding refinement around M point','diag')
   call add_refinement_around_point(Kpts, ip, M_point, b1, b2, kx_step, ky_step, n_refine)

   call MIO_Print('  Adding refinement around K'' point','diag')
   call add_refinement_around_point(Kpts, ip, K_prime, b1, b2, kx_step, ky_step, n_refine)

   call MIO_Print('  Total k-points after refinement: '//trim(num2str(ip)),'diag')

end subroutine add_high_symmetry_refinement

subroutine add_refinement_around_point(Kpts, ip, point_frac, b1, b2, kx_step, ky_step, n_refine)
   ! Add refinement grid around a specific high-symmetry point
   use constants, only : cmplx_0
   implicit none

   ! Input/Output parameters
   real(dp), intent(inout) :: Kpts(3, *)
   integer, intent(inout) :: ip
   real(dp), intent(in) :: point_frac(2), b1(3), b2(3)
   real(dp), intent(in) :: kx_step, ky_step
   integer, intent(in) :: n_refine

   ! Local variables
   integer :: i, j
   real(dp) :: kx, ky, k_refine(3)
   real(dp) :: start_x, start_y, end_x, end_y

   ! Define refinement region around the point
   start_x = point_frac(1) - real(n_refine, dp) * kx_step / 2.0_dp
   end_x = point_frac(1) + real(n_refine, dp) * kx_step / 2.0_dp
   start_y = point_frac(2) - real(n_refine, dp) * ky_step / 2.0_dp
   end_y = point_frac(2) + real(n_refine, dp) * ky_step / 2.0_dp

   ! Add refined grid points
   do j = 0, n_refine-1
      do i = 0, n_refine-1
         kx = start_x + real(i, dp) * kx_step
         ky = start_y + real(j, dp) * ky_step

         ! Convert to Cartesian coordinates
         k_refine = kx * b1 + ky * b2
         k_refine(3) = 0.0_dp  ! 2D system

         ! Check if point is within first Brillouin zone (optional)
         if (is_in_first_BZ(kx, ky)) then
            ip = ip + 1
            Kpts(:, ip) = k_refine
         end if
      end do
   end do

end subroutine add_refinement_around_point

function is_in_first_BZ(kx_frac, ky_frac) result(in_BZ)
   ! Check if fractional k-point is within first Brillouin zone
   ! For hexagonal lattice: |kx| ≤ 1/2, |ky| ≤ 1/2, |kx + ky| ≤ 1/2
   implicit none

   real(dp), intent(in) :: kx_frac, ky_frac
   logical :: in_BZ

   ! Simple rectangular BZ check (can be made more sophisticated)
   in_BZ = (abs(kx_frac) <= 0.5_dp) .and. (abs(ky_frac) <= 0.5_dp)

end function is_in_first_BZ

subroutine output_berry_curvature_data(kpoint_index, kpt, berry_curv_bands, band_indices, n_target_bands)
   ! Output Berry curvature data for hotspot analysis
   ! Format: kx, ky, BerryCurvature_band1, BerryCurvature_band2, ...
   implicit none

   ! Input parameters
   integer, intent(in) :: kpoint_index, n_target_bands
   real(dp), intent(in) :: kpt(3)
   real(dp), intent(in) :: berry_curv_bands(:)
   integer, intent(in) :: band_indices(:)

   ! Local variables
   integer :: i, iunit
   character(len=100) :: filename
   logical :: file_exists

   ! Only output for target bands (first n_target_bands)
   if (n_target_bands <= 0 .or. n_target_bands > size(berry_curv_bands)) return

   ! Create filename based on target bands
   ! Since bands may not be consecutive, show all band numbers
   if (n_target_bands == 1) then
      write(filename, '(A,I0,A)') 'berry_curvature_bands_', band_indices(1), '.dat'
   else if (n_target_bands == 2) then
      write(filename, '(A,I0,A,I0,A)') 'berry_curvature_bands_', band_indices(1), '_and_', band_indices(2), '.dat'
   else
      ! For 3+ bands, show first and last with count
      write(filename, '(A,I0,A,I0,A,I0,A)') 'berry_curvature_bands_', band_indices(1), '_to_', band_indices(n_target_bands), '_', n_target_bands, 'bands.dat'
   end if

   ! Open file (append mode for multiple k-points)
   if (kpoint_index == 1) then
      ! First k-point: create new file with header
      open(newunit=iunit, file=trim(filename), status='replace', action='write')
      write(iunit, '(A)') '# Berry curvature data for hotspot analysis'
      write(iunit, '(A)') '# Format: kx, ky, BerryCurvature_band1, BerryCurvature_band2, ...'
      write(iunit, '(A)') '# Target bands: '//trim(num2str(band_indices(1)))//' to '//trim(num2str(band_indices(n_target_bands)))
      write(iunit, '(A)') '#'
   else
      ! Subsequent k-points: append to existing file
      open(newunit=iunit, file=trim(filename), status='old', position='append', action='write')
   end if

   ! Write data: kx, ky, Berry curvature for each target band
   write(iunit, '(F12.6,2X,F12.6)', advance='no') kpt(1), kpt(2)
   do i = 1, n_target_bands
      ! Use scientific notation for large numbers to prevent overflow
      if (abs(berry_curv_bands(i)) > 1e6_dp .or. abs(berry_curv_bands(i)) < 1e-6_dp) then
         write(iunit, '(2X,ES15.6)', advance='no') berry_curv_bands(i)
      else
         write(iunit, '(2X,F15.8)', advance='no') berry_curv_bands(i)
      end if
   end do
   write(iunit, *)  ! New line

   close(iunit)

   ! Print summary for first and last k-points
   if (kpoint_index == 1) then
      call MIO_Print('=== BERRY CURVATURE OUTPUT ===','diag')
      call MIO_Print('Writing Berry curvature data to: '//trim(filename),'diag')
      call MIO_Print('Target bands: '//trim(num2str(band_indices(1)))//' to '//trim(num2str(band_indices(n_target_bands))),'diag')
      call MIO_Print('Format: kx, ky, BerryCurvature_band1, BerryCurvature_band2, ...','diag')
   end if

   ! Print hotspot detection for first k-point
   if (kpoint_index == 1) then
      call MIO_Print('=== HOTSPOT DETECTION (k-point 1) ===','diag')
      do i = 1, n_target_bands
         if (abs(berry_curv_bands(i)) > 1.0_dp) then
            call MIO_Print('  HOTSPOT: Band '//trim(num2str(band_indices(i)))//' |Ω| = '//trim(num2str(abs(berry_curv_bands(i)),6)),'diag')
         else if (abs(berry_curv_bands(i)) > 0.1_dp) then
            call MIO_Print('  Moderate: Band '//trim(num2str(band_indices(i)))//' |Ω| = '//trim(num2str(abs(berry_curv_bands(i)),6)),'diag')
         end if
      end do
   end if

end subroutine output_berry_curvature_data

!> @brief Apply SOC modifications to Hamiltonian matrix
!! @param[in]     i      Atom index
!! @param[in]     is     Spin channel (1=spin-up, 2=spin-down)
!! @param[in]     ns     Number of spin channels
!! @param[inout]  HLoc   Hamiltonian matrix being built
!! @details Applies SOC terms including Zeeman, intrinsic SOC, Rashba, etc.
!!          to the Hamiltonian matrix during construction.
subroutine ApplySOCtoHamiltonian(i, is, ns, HLoc)

   use ham, only : Zterm, gZeeman, IntrinsicSOCterm, lambdaI, IsingSOCterm, lambdaIsing, PIASOCterm, lambdaPIA, SOCEnabledForLayer
   use magf, only : BmagZeeman
   use atoms, only : Species, layerIndex

   integer, intent(in) :: i, is, ns
   complex(dp), intent(inout) :: HLoc(:,:)

   ! Apply Zeeman effect (λVZ) - spin-dependent onsite
   if (Zterm .and. (ns==2)) then
      if (is==1) then
         HLoc(i,i) = HLoc(i,i) + gZeeman*BmagZeeman  ! Spin-up: +λVZ
      else
         HLoc(i,i) = HLoc(i,i) - gZeeman*BmagZeeman  ! Spin-down: -λVZ
      end if
   end if

   ! Apply Intrinsic SOC (λI) - gap-opening onsite term
   ! For graphene, this should be sublattice-dependent to open gap
   if (IntrinsicSOCterm .and. SOCEnabledForLayer(layerIndex(i))) then
      ! Apply +λI to sublattice A, -λI to sublattice B (or vice versa)
      ! This opens a gap of size 2λI at Dirac point
      ! Species(i)=1 is A sublattice, Species(i)=2 is B sublattice
      if (Species(i) == 1) then
         HLoc(i,i) = HLoc(i,i) + lambdaI
      else
         HLoc(i,i) = HLoc(i,i) - lambdaI
      end if
   end if

   ! Apply Ising SOC (λIsing) - Valley-Zeeman term (τ_z s_z)
   ! This is a valley-dependent spin splitting: τ_z = +1 for K valley, -1 for K' valley
   ! Does not open a global gap, only splits spins differently in K vs K' valleys
   if (IsingSOCterm .and. SOCEnabledForLayer(layerIndex(i)) .and. (ns==2)) then
      ! Valley sign: +1 for K valley, -1 for K' valley
      if (.not. useKprimeValley) then
         ! K valley: τ_z = +1
         if (is == 1) then
            ! Spin-up: +λ for all sites
            HLoc(i,i) = HLoc(i,i) + lambdaIsing
         else
            ! Spin-down: -λ for all sites
            HLoc(i,i) = HLoc(i,i) - lambdaIsing
         end if
      else
         ! K' valley: τ_z = -1 (flip sign)
         if (is == 1) then
            ! Spin-up: -λ for all sites
            HLoc(i,i) = HLoc(i,i) - lambdaIsing
         else
            ! Spin-down: +λ for all sites
            HLoc(i,i) = HLoc(i,i) + lambdaIsing
         end if
      end if
   end if

   ! Apply Pseudo-inversion asymmetry (λPIA) - onsite component
   if (PIASOCterm .and. SOCEnabledForLayer(layerIndex(i))) then
      ! PIA is spin-independent onsite term
      HLoc(i,i) = HLoc(i,i) + lambdaPIA
   end if

   ! Note: Rashba SOC (λR) requires block Hamiltonian and is handled separately

end subroutine ApplySOCtoHamiltonian

!> @brief Apply PIA hopping terms to the Hamiltonian
!! @param[in]     i       Atom index
!! @param[in]     j       Neighbor index
!! @param[in]     in      Neighbor atom index
!! @param[in]     N       Number of atoms
!! @param[inout]  HLoc    Hamiltonian (N×N or 2N×2N depending on isBlock)
!! @param[in]     isBlock If true, HLoc is 2N×2N block Hamiltonian
subroutine ApplyPIAHopping(i, j, in, N, HLoc, isBlock)

   use ham, only : PIASOCterm, lambdaPIA
   use neigh, only : neighD
   use atoms, only : Species

   integer, intent(in) :: i, j, in, N
   logical, intent(in) :: isBlock
   complex(dp), intent(inout) :: HLoc(:,:)

   complex(dp) :: phase

   if (.not. PIASOCterm) return

   ! PIA hopping: iλPIA * (τ_y + i*τ_x) where τ is the neighbor vector
   ! This creates imaginary hopping terms
   if (isBlock) then
      ! Block Hamiltonian: apply to both spin blocks
      phase = lambdaPIA * cmplx(neighD(2, j, i), neighD(1, j, i))
      HLoc(in, i) = HLoc(in, i) + phase
      HLoc(i, in) = HLoc(i, in) - conjg(phase) ! Hermitian conjugate
      HLoc(in + N, i + N) = HLoc(in + N, i + N) + phase
      HLoc(i + N, in + N) = HLoc(i + N, in + N) - conjg(phase)
   else
      ! Sequential Hamiltonian: apply normally
      phase = lambdaPIA * cmplx(neighD(2, j, i), neighD(1, j, i))
      HLoc(in, i) = HLoc(in, i) + phase
      HLoc(i, in) = HLoc(i, in) - conjg(phase)
   end if

end subroutine ApplyPIAHopping

!> @brief Apply Rashba SOC (spin-flip) terms to the block Hamiltonian
!! @param[in]     i       Atom index
!! @param[in]     j       Neighbor index
!! @param[in]     in      Neighbor atom index
!! @param[in]     KLoc    k-point vector
!! @param[in]     R       Displacement vector in cell coordinates
!! @param[in]     N       Number of atoms
!! @param[inout]  HLoc    Block Hamiltonian (2N×2N)
subroutine ApplySpinFlipSOC(i, j, in, KLoc, R, N, HLoc)

   use ham, only : lambdaR, RashbaSOCterm, SOCEnabledForLayer
   use neigh, only : neighD
   use atoms, only : Species, layerIndex, frac, AtomsSetCart
   use constants, only : cmplx_i

   integer, intent(in) :: i, j, in, N
   real(dp), intent(in) :: KLoc(3), R(3)
   complex(dp), intent(inout) :: HLoc(2*N, 2*N)

   real(dp) :: dist, acc
   real(dp) :: dx, dy, dnorm
   complex(dp) :: RashbaHopp

   if (.not. RashbaSOCterm) return

   ! Check if SOC should be applied to this layer
   if (.not. SOCEnabledForLayer(layerIndex(i))) then
      ! Debug output when SOC is skipped due to layer control
      if (socDebug .and. i <= 5 .and. j == 1) then
         call MIO_Print('SOC skipped: atom i='//trim(num2str(i))//' is in layer '//trim(num2str(layerIndex(i)))//' (not in SOCLayers list)', 'diag')
      end if
      return
   end if

   ! Ensure we're using Cartesian coordinates
   if (frac) call AtomsSetCart()

   ! Debug: Check why Rashba might not be applied
   if (socDebug .and. i <= 2 .and. j <= 10) then  ! Show first 10 neighbors instead of just 2
      acc = 2.46_dp / sqrt(3.0_dp)
      dist = sqrt(neighD(1, j, i)**2.0_dp + neighD(2, j, i)**2.0_dp)
      call MIO_Print('Rashba Check: i='//trim(num2str(i))//', j='//trim(num2str(j))//', in='//trim(num2str(in))//', layer_i='//trim(num2str(layerIndex(i)))//', layer_in='//trim(num2str(layerIndex(in))), 'diag')
      call MIO_Print('Species_i='//trim(num2str(Species(i)))//', Species_in='//trim(num2str(Species(in)))//', dist='//trim(num2str(dist))//', cutoff='//trim(num2str(acc*1.1_dp)), 'diag')

      ! Check each condition individually
      if (i == in) then
         call MIO_Print('  -> REJECTED: Self-loop (i == in)', 'diag')
      else if (layerIndex(i) /= layerIndex(in)) then
         call MIO_Print('  -> REJECTED: Different layers', 'diag')
      else if (Species(i) == Species(in)) then
         call MIO_Print('  -> REJECTED: Same species', 'diag')
      else if (dist >= acc * 1.1_dp) then
         call MIO_Print('  -> REJECTED: Distance too large (dist='//trim(num2str(dist))//' >= cutoff='//trim(num2str(acc*1.1_dp))//')', 'diag')
      else
         call MIO_Print('  -> ALL CONDITIONS MET! (A-B nearest neighbor)', 'diag')
      end if
   end if

   ! Skip self-loops
   if (i == in) return

   ! Rashba only applies within same layer, between different species, and first nearest neighbors
   acc = 2.46_dp / sqrt(3.0_dp)  ! Graphene lattice constant / sqrt(3) ≈ 1.42 Å
   dist = sqrt(neighD(1, j, i)**2.0_dp + neighD(2, j, i)**2.0_dp)

   if ( (layerIndex(i) == layerIndex(in)) .and. &
     (Species(i) /= Species(in)) .and. &
     (dist < acc * 1.1_dp) ) then

     ! Debug: Rashba is being applied!
     if (socDebug .and. i <= 2 .and. j <= 10) then  ! Show first 10 neighbors instead of just 2
        call MIO_Print('*** RASHBA APPLIED *** i='//trim(num2str(i))//', j='//trim(num2str(j))//', in='//trim(num2str(in))//', dist='//trim(num2str(dist))//', cutoff='//trim(num2str(acc*1.1_dp))//' (A-B nearest neighbor)', 'diag')
     end if

     dx = neighD(1, j, i)
     dy = neighD(2, j, i)
     dnorm = sqrt(dx*dx + dy*dy)
     if (dnorm > 1.0e-12_dp) then
        dx = dx / dnorm
        dy = dy / dnorm
     end if

     RashbaHopp = lambdaR * (dx + cmplx_i * dy) * exp(-cmplx_i * dot_product(KLoc, R))

     ! Debug output for first few bonds
     if (socDebug .and. i <= 2 .and. j <= 10) then  ! Show first 10 neighbors instead of just 2
        call MIO_Print('Rashba Debug: i='//trim(num2str(i))//', j='//trim(num2str(j))//', in='//trim(num2str(in))//', lambdaR='//trim(num2str(lambdaR))//', dx='//trim(num2str(dx))//', dy='//trim(num2str(dy)), 'diag')
        call MIO_Print('RashbaHopp = '//trim(num2str(real(RashbaHopp)))//' + i*'//trim(num2str(aimag(RashbaHopp))), 'diag')
        call MIO_Print('Adding to HLoc('//trim(num2str(i))//','//trim(num2str(in+N))//') = spin-flip term', 'diag')
     end if

     ! Spin-flip couplings (↑→↓ and ↓→↑)
     HLoc(i     , in + N) = HLoc(i     , in + N) + RashbaHopp
     HLoc(in + N, i     ) = HLoc(in + N, i     ) + conjg(RashbaHopp)

     ! Debug: Check if we're accidentally adding diagonal terms (which would cause global shift)
     if (socDebug .and. i <= 2 .and. j <= 10) then
        if (i == in) then
           call MIO_Print('ERROR: Adding Rashba to diagonal term HLoc('//trim(num2str(i))//','//trim(num2str(i))//') - this causes global shift!', 'diag')
        else if (i == in + N) then
           call MIO_Print('ERROR: Adding Rashba to diagonal term HLoc('//trim(num2str(i))//','//trim(num2str(in+N))//') - this causes global shift!', 'diag')
        else
           call MIO_Print('OK: Adding Rashba to off-diagonal spin-flip term HLoc('//trim(num2str(i))//','//trim(num2str(in+N))//')', 'diag')
        end if
     end if

   end if
end subroutine ApplySpinFlipSOC

!> @brief Build block Hamiltonian for spin-flip SOC terms
!! @details Creates a 2N×2N block Hamiltonian:
!!          [ H_up + SOC_diag    H_cross(Rashba)   ]
!!          [ H_cross* (Rashba)   H_dn + SOC_diag   ]
subroutine BuildBlockHamiltonian(N, KLoc, cell, H0, maxN, hopp, NList, Nneigh, neighCell, &
                                  HBlock, EBlock)

   use constants, only : cmplx_i
   use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf, only : charge, Zch
   use atoms, only : Species, layerIndex
   use tbpar, only : U
   use ham, only : RashbaSOCterm, nspin
   use neigh, only : neighD

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)
   complex(dp), intent(out) :: HBlock(2*N, 2*N)
   real(dp), intent(out) :: EBlock(2*N)

   integer :: i, j, in, info
   real(dp) :: R(3)
   integer, parameter :: MAX_BLOCK_SIZE = 20000  ! Max atoms for block Hamiltonian
   integer :: lwork_block
   complex(dp) :: ZWorkLoc(2*(2*MAX_BLOCK_SIZE)-1)  ! Static workspace for safety
   real(dp) :: DWorkLoc(3*(2*MAX_BLOCK_SIZE)-2)

   ! Check that we don't exceed MAX_BLOCK_SIZE
   if (2*N > 2*MAX_BLOCK_SIZE) then
      call MIO_Kill('Block Hamiltonian size exceeds MAX_BLOCK_SIZE', 'diag', 'BuildBlockHamiltonian')
   end if

   ! Workspace sizes for 2*N matrix
   lwork_block = max(1, 2*(2*N)-1)

   ! Initialize block Hamiltonian
   HBlock = 0.0_dp

   ! Build Hamiltonian for both spin blocks
   ! Only set lower triangular elements (row >= col) for 'L' format in ZHEEV
   do i = 1, N
      ! Onsite energies for both spin channels (diagonal elements, always lower triangular)
      HBlock(i, i) = H0(i)         ! Spin-up diagonal
      HBlock(i + N, i + N) = H0(i) ! Spin-down diagonal

      ! Apply SCF potential if spin-polarized and SCF is enabled
      !   ! Use charge data even for non-SCF (charge is still available)

      ! Apply diagonal SOC terms to both blocks
      call ApplySOCtoBlock(i, HBlock)

      ! Build hopping for both spin blocks
      ! Debug: Check neighbor count for first few atoms
      if (socDebug .and. i <= 2) then
         call MIO_Print('BuildBlockHamiltonian: Atom '//trim(num2str(i))//' has '//trim(num2str(Nneigh(i)))//' neighbors', 'diag')
      end if

      do j = 1, Nneigh(i)
         in = NList(j, i)

         ! Debug: Check if we're actually looping through all neighbors
         if (socDebug .and. i <= 2 .and. j <= 5) then
            call MIO_Print('  Neighbor loop: i='//trim(num2str(i))//', j='//trim(num2str(j))//', in='//trim(num2str(in))//', total_neighbors='//trim(num2str(Nneigh(i))), 'diag')
         end if

         ! Only set lower triangular elements (row >= col)
            ! Use actual atomic position difference instead of lattice vector
            ! NeighD(:,j,i) contains the vector from atom i to atom j: (τ_j - τ_i) + T_ij
            ! For consistency with paper formulation, use -NeighD to get distance from i to j
            ! This matches transform_dense_hamiltonian_tapw
            R(1:3) = 0.0_dp
            R(1:2) = -NeighD(1:2, j, i)  ! Use negative to get distance from i to j

            ! Regular hopping (preserved in both blocks)
            HBlock(in, i) = HBlock(in, i) - hopp(j, i) * exp(-cmplx_i * dot_product(KLoc, R))
            HBlock(in + N, i + N) = HBlock(in + N, i + N) - hopp(j, i) * exp(-cmplx_i * dot_product(KLoc, R))

            ! Apply PIA hopping terms if enabled (matches non-TAPW SOC path)
            ! PIA hopping not yet properly implemented - commented out

            ! Apply Rashba spin-flip terms if enabled
            if (RashbaSOCterm) then
               call ApplySpinFlipSOC(i, j, in, KLoc, R, N, HBlock)
            end if
      end do
   end do

   ! Edge hopping if applicable
   if (edgeHopp) then
      do i = 1, nQ
         do j = 1, nEdgeN(i)
            in = NeI(j, i)
            R = matmul(cell, NedgeCell(:, j, i))
            ! Apply to both spin blocks
            HBlock(in, edgeIndx(i)) = HBlock(in, edgeIndx(i)) + edgeH(j, i) * exp(-cmplx_i * dot_product(KLoc, R))
            HBlock(in + N, edgeIndx(i) + N) = HBlock(in + N, edgeIndx(i) + N) + edgeH(j, i) * exp(-cmplx_i * dot_product(KLoc, R))
         end do
      end do
   end if

   ! Diagonalize block Hamiltonian
   ! ZHEEV signature: (jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
   call ZHEEV('N', 'L', 2*N, HBlock, 2*N, EBlock, ZWorkLoc, lwork_block, DWorkLoc, info)
   if (info /= 0) then
      call MIO_Print('ZHEEV returned info = '//trim(num2str(info)), 'diag')
      call MIO_Print('N = '//trim(num2str(N))//', 2*N = '//trim(num2str(2*N)), 'diag')
      call MIO_Print('lwork_block = '//trim(num2str(lwork_block)), 'diag')
      call MIO_Kill('Error in block Hamiltonian diagonalization', 'diag', 'BuildBlockHamiltonian')
   end if

end subroutine BuildBlockHamiltonian

!> @brief Build block Hamiltonian for spin-flip SOC terms (build only, no diagonalization)
!! @details Creates a 2N×2N block Hamiltonian:
!!          [ H_up + SOC_diag    H_cross(Rashba)   ]
!!          [ H_cross* (Rashba)   H_dn + SOC_diag   ]
!! @param[in]     N         Number of atoms
!! @param[in]     KLoc      k-point vector
!! @param[in]     cell      Unit cell matrix
!! @param[in]     H0        Onsite energies
!! @param[in]     maxN      Maximum number of neighbors
!! @param[in]     hopp      Hopping parameters
!! @param[in]     NList     Neighbor list
!! @param[in]     Nneigh    Number of neighbors per atom
!! @param[in]     neighCell Neighbor cell indices
!! @param[out]    HBlock    Block Hamiltonian (2N×2N)
subroutine BuildBlockHamiltonianOnly(N, KLoc, cell, H0, maxN, hopp, NList, Nneigh, neighCell, HBlock)

   use constants, only : cmplx_i
   use interface, only : edgeHopp, nEdgeN, edgeH, nQ, edgeIndx, NeI, NedgeCell
   use scf, only : charge, Zch
   use atoms, only : Species, layerIndex
   use tbpar, only : U
   use ham, only : RashbaSOCterm, nspin
   use neigh, only : neighD

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N)
   real(dp), intent(in) :: KLoc(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)
   complex(dp), intent(out) :: HBlock(2*N, 2*N)

   integer :: i, j, in
   real(dp) :: R(3)
   integer, parameter :: MAX_BLOCK_SIZE = 20000  ! Max atoms for block Hamiltonian

   ! Check that we don't exceed MAX_BLOCK_SIZE
   if (2*N > 2*MAX_BLOCK_SIZE) then
      call MIO_Kill('Block Hamiltonian size exceeds MAX_BLOCK_SIZE', 'diag', 'BuildBlockHamiltonianOnly')
   end if

   ! Initialize block Hamiltonian
   HBlock = 0.0_dp

   ! Build Hamiltonian for both spin blocks
   ! Only set lower triangular elements (row >= col) for 'L' format in ZHEEV
   do i = 1, N
      ! Onsite energies for both spin channels (diagonal elements, always lower triangular)
      HBlock(i, i) = H0(i)         ! Spin-up diagonal
      HBlock(i + N, i + N) = H0(i) ! Spin-down diagonal

      ! Apply SCF potential if spin-polarized and SCF is enabled
      !   ! Use charge data even for non-SCF (charge is still available)

      ! Apply diagonal SOC terms to both blocks
      call ApplySOCtoBlock(i, HBlock)

      ! Build hopping for both spin blocks
      ! Debug: Check neighbor count for first few atoms
      if (socDebug .and. i <= 2) then
         call MIO_Print('BuildBlockHamiltonianOnly: Atom '//trim(num2str(i))//' has '//trim(num2str(Nneigh(i)))//' neighbors', 'diag')
      end if

      do j = 1, Nneigh(i)
         in = NList(j, i)

         ! Debug: Check if we're actually looping through all neighbors
         if (socDebug .and. i <= 2 .and. j <= 5) then
            call MIO_Print('  Neighbor loop: i='//trim(num2str(i))//', j='//trim(num2str(j))//', in='//trim(num2str(in))//', total_neighbors='//trim(num2str(Nneigh(i))), 'diag')
         end if

         ! Debug: Progress indicator for atom 1
         if (socDebug .and. i == 1 .and. modulo(j, 10) == 0) then
            call MIO_Print('BuildBlockHamiltonianOnly: Atom 1 progress: j='//trim(num2str(j))//'/'//trim(num2str(Nneigh(i))), 'diag')
         end if

         ! Only set lower triangular elements (row >= col)
            ! Use actual atomic position difference instead of lattice vector
            ! NeighD(:,j,i) contains the vector from atom i to atom j: (τ_j - τ_i) + T_ij
            ! For consistency with paper formulation, use -NeighD to get distance from i to j
            ! This matches transform_dense_hamiltonian_tapw
            R(1:3) = 0.0_dp
            R(1:2) = -NeighD(1:2, j, i)  ! Use negative to get distance from i to j

            ! Regular hopping (preserved in both blocks)
            HBlock(in, i) = HBlock(in, i) - hopp(j, i) * exp(-cmplx_i * dot_product(KLoc, R))
            HBlock(in + N, i + N) = HBlock(in + N, i + N) - hopp(j, i) * exp(-cmplx_i * dot_product(KLoc, R))

            ! Apply PIA hopping terms if enabled (matches non-TAPW SOC path)
            ! PIA hopping not yet properly implemented - commented out

            ! Apply Rashba spin-flip terms if enabled
            ! Note: ApplySpinFlipSOC uses NeighD internally for bond direction (unit vector)
            ! This is correct - R (for phase) and NeighD (for direction) serve different purposes
            if (RashbaSOCterm) then
               call ApplySpinFlipSOC(i, j, in, KLoc, R, N, HBlock)
            end if
      end do

      ! Debug: Show completion of atom i in BuildBlockHamiltonianOnly
      if (socDebug .and. i <= 5) then
         call MIO_Print('BuildBlockHamiltonianOnly: Finished atom '//trim(num2str(i)), 'diag')
      end if
   end do

   ! Debug: Show completion of BuildBlockHamiltonianOnly
   if (socDebug) then
      call MIO_Print('BuildBlockHamiltonianOnly: Completed for all atoms, moving to edge hopping', 'diag')
   end if

   ! Edge hopping if applicable
   if (edgeHopp) then
      do i = 1, nQ
         do j = 1, nEdgeN(i)
            in = NeI(j, i)
            R = matmul(cell, NedgeCell(:, j, i))
            ! Apply to both spin blocks
            HBlock(in, edgeIndx(i)) = HBlock(in, edgeIndx(i)) + edgeH(j, i) * exp(-cmplx_i * dot_product(KLoc, R))
            HBlock(in + N, edgeIndx(i) + N) = HBlock(in + N, edgeIndx(i) + N) + edgeH(j, i) * exp(-cmplx_i * dot_product(KLoc, R))
         end do
      end do
   end if

end subroutine BuildBlockHamiltonianOnly

!> @brief Diagonalize block Hamiltonian using ZHEEV
!! @param[in]     HBlock    Block Hamiltonian (2N×2N)
!! @param[out]    EBlock    Eigenvalues (2N)
!! @param[in]     N         Number of atoms
subroutine DiagBlockHamiltonian(HBlock, EBlock, N)

   integer, intent(in) :: N
   complex(dp), intent(in) :: HBlock(2*N, 2*N)
   real(dp), intent(out) :: EBlock(2*N)

   integer :: info
   integer, parameter :: MAX_BLOCK_SIZE = 20000  ! Max atoms for block Hamiltonian
   integer :: lwork_block
   complex(dp) :: ZWorkLoc(2*(2*MAX_BLOCK_SIZE)-1)  ! Static workspace for safety
   real(dp) :: DWorkLoc(3*(2*MAX_BLOCK_SIZE)-2)

   ! Check that we don't exceed MAX_BLOCK_SIZE
   if (2*N > 2*MAX_BLOCK_SIZE) then
      call MIO_Kill('Block Hamiltonian size exceeds MAX_BLOCK_SIZE', 'diag', 'DiagBlockHamiltonian')
   end if

   ! Workspace sizes for 2*N matrix
   lwork_block = max(1, 2*(2*N)-1)

   ! Debug: Warn about large matrix diagonalization
   if (socDebug) then
      call MIO_Print('DiagBlockHamiltonian: Starting ZHEEV diagonalization for '//trim(num2str(2*N))//'×'//trim(num2str(2*N))//' matrix (this may take several minutes to hours for large systems)', 'diag')
      call MIO_Print('  Matrix size: N='//trim(num2str(N))//', 2N='//trim(num2str(2*N))//', estimated memory: ~'//trim(num2str(int(2*N*2*N*16.0_dp/1024.0_dp/1024.0_dp/1024.0_dp)))//' GB', 'diag')
   end if

   ! Diagonalize block Hamiltonian
   ! ZHEEV signature: (jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
   call ZHEEV('N', 'L', 2*N, HBlock, 2*N, EBlock, ZWorkLoc, lwork_block, DWorkLoc, info)

   ! Debug: Confirm completion
   if (socDebug) then
      call MIO_Print('DiagBlockHamiltonian: ZHEEV completed successfully', 'diag')
   end if
   if (info /= 0) then
      call MIO_Print('ZHEEV returned info = '//trim(num2str(info)), 'diag')
      call MIO_Print('N = '//trim(num2str(N))//', 2*N = '//trim(num2str(2*N)), 'diag')
      call MIO_Print('lwork_block = '//trim(num2str(lwork_block)), 'diag')
      call MIO_Kill('Error in block Hamiltonian diagonalization', 'diag', 'DiagBlockHamiltonian')
   end if

end subroutine DiagBlockHamiltonian

!> @brief Apply diagonal SOC terms to block Hamiltonian
!! @param[in]     i       Atom index
!! @param[inout]  HBlock  Block Hamiltonian (2N×2N)
subroutine ApplySOCtoBlock(i, HBlock)

   use ham, only : Zterm, gZeeman, IntrinsicSOCterm, lambdaI, IsingSOCterm, lambdaIsing, PIASOCterm, lambdaPIA, nspin, SOCEnabledForLayer
   use magf, only : BmagZeeman
   use atoms, only : nAt, Species, layerIndex

   integer, intent(in) :: i
   complex(dp), intent(inout) :: HBlock(2*nAt, 2*nAt)

   ! Apply Zeeman effect to block Hamiltonian
   if (Zterm .and. nspin == 2) then
      HBlock(i, i) = HBlock(i, i) + gZeeman * BmagZeeman       ! Spin-up: +λVZ
      HBlock(i + nAt, i + nAt) = HBlock(i + nAt, i + nAt) - gZeeman * BmagZeeman  ! Spin-down: -λVZ
   end if

   ! Apply Intrinsic SOC (sublattice-dependent for gap opening)
   if (IntrinsicSOCterm .and. SOCEnabledForLayer(layerIndex(i))) then
      ! Apply +λI to sublattice A, -λI to sublattice B to open gap
      if (Species(i) == 1) then
         HBlock(i, i) = HBlock(i, i) + lambdaI
         HBlock(i + nAt, i + nAt) = HBlock(i + nAt, i + nAt) + lambdaI
      else
         HBlock(i, i) = HBlock(i, i) - lambdaI
         HBlock(i + nAt, i + nAt) = HBlock(i + nAt, i + nAt) - lambdaI
      end if
   end if

   ! Apply Ising SOC - Valley-Zeeman term (τ_z s_z)
   ! This is a valley-dependent spin splitting: τ_z = +1 for K valley, -1 for K' valley
   ! Does not open a global gap, only splits spins differently in K vs K' valleys
   if (IsingSOCterm .and. SOCEnabledForLayer(layerIndex(i)) .and. nspin == 2) then
      ! Valley sign: +1 for K valley, -1 for K' valley
      if (.not. useKprimeValley) then
         ! K valley: τ_z = +1
         ! Spin-up (first block): +λ for all sites
         HBlock(i, i) = HBlock(i, i) + lambdaIsing
         ! Spin-down (second block): -λ for all sites
         HBlock(i + nAt, i + nAt) = HBlock(i + nAt, i + nAt) - lambdaIsing
      else
         ! K' valley: τ_z = -1 (flip sign)
         ! Spin-up (first block): -λ for all sites
         HBlock(i, i) = HBlock(i, i) - lambdaIsing
         ! Spin-down (second block): +λ for all sites
         HBlock(i + nAt, i + nAt) = HBlock(i + nAt, i + nAt) + lambdaIsing
      end if
   end if

   ! Apply Pseudo-inversion asymmetry
   if (PIASOCterm .and. SOCEnabledForLayer(layerIndex(i))) then
      HBlock(i, i) = HBlock(i, i) + lambdaPIA
      HBlock(i + nAt, i + nAt) = HBlock(i + nAt, i + nAt) + lambdaPIA
   end if

end subroutine ApplySOCtoBlock

#ifdef SEMICL
#include "diag_semicl.inc"
#endif

end module diag

#ifdef NO_MKL
! Builds without Intel MKL (compile with -DNO_MKL): the two MKL entry points
! used in this file are replaced here. Shift-invert diagonalisation
! (Diag.SparseUseShift, Diag.SparseCount) needs MKL PARDISO and stops with a
! message; everything else is unaffected.
subroutine pardiso(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)

   use mio

   integer(8) :: pt(*)
   integer :: maxfct, mnum, mtype, phase, n, ia(*), ja(*), perm(*), nrhs, iparm(*), msglvl, error
   complex(8) :: a(*), b(*), x(*)

   error = -1
   call MIO_Kill('Diag.SparseUseShift and Diag.SparseCount need MKL PARDISO; this executable was built '// &
     'without MKL (-DNO_MKL)','diag','pardiso')

end subroutine pardiso

integer function mkl_get_max_threads()

#ifdef _OPENMP
   use omp_lib, only : omp_get_max_threads
   mkl_get_max_threads = omp_get_max_threads()
#else
   mkl_get_max_threads = 1
#endif

end function mkl_get_max_threads
#endif
