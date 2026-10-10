module kubo

   use mio
   use atoms,               only : in1, in2, nAt, inode1, inode2

   implicit none

   PRIVATE

!#ifdef MPI
!#define USE_MPI 1
!#include <sprng_f.h>
!#endif /* MPI */

   public :: KuboInitWF
   public :: KuboInitWFPDOS
   public :: KuboInitWFLayerDOS
   public :: KuboInitWFLayerAndSpeciesDOS
   public :: KuboDOS
   public :: KuboTEvol

contains

subroutine KuboInitWF(Psi)

   use constants,            only : twopi, cmplx_i, cmplx_0, cmplx_1
   use parallel,             only : nDiv, procID

   complex(dp), intent(out) :: Psi(inode1:)

   integer :: i, clock, n
   real(dp) :: r
   complex(dp) :: cnum, s
   !type(rand_t), save :: rng
   integer, pointer :: seed(:)
   logical PDOS,setSeed
   integer PDOSAtomNumber, seedValue

#ifdef DEBUG
   call MIO_Debug('KuboInitWF',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   call random_seed(size = n)
   allocate(seed(n))
   call system_clock(COUNT=clock)
   call MIO_InputParameter('Kubo.SetSeed',setSeed,.false.)
   call MIO_InputParameter('Kubo.SeedValue',seedValue,123456)
   if (setSeed) then
       seed = seedValue
       call MIO_Print('seedValue considered from Gendata','kubo')
   else
       seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   end if
   call random_seed(PUT = seed)

   do i=inode1, inode2
      call random_number(r)
      Psi(i) = exp(twopi*cmplx_i*r)/sqrt(real(nAt))
   end do

   cnum = 0.0_dp
   !$OMP PARALLEL DO REDUCTION(+:cnum)
   do i=1,nAt
      cnum = cnum + Psi(i)*conjg(Psi(i))
   end do
   !$OMP END PARALLEL DO
#ifdef MPI
   call MPIRedSum(cnum,s,1,MPI_DOUBLE_COMPLEX)
   cnum = s
#endif /* MPI */
   call MIO_Print('Initial norm of wave function: '//num2str(sqrt(abs(cnum)),12),'kubo')

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboInitWF',1)
#endif /* DEBUG */

end subroutine KuboInitWF

subroutine KuboInitWFPDOS(Psi,PDOSAtomNumber)

   use constants,            only : twopi, cmplx_i, cmplx_0, cmplx_1
   use parallel,             only : nDiv, procID

   complex(dp), intent(out) :: Psi(inode1:)

   integer :: i, clock, n
   real(dp) :: r
   complex(dp) :: cnum, s
   !type(rand_t), save :: rng
   integer, pointer :: seed(:)
   integer PDOSAtomNumber

#ifdef DEBUG
   call MIO_Debug('KuboInitWFPDOS',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   call random_seed(size = n)
   allocate(seed(n))
   call system_clock(COUNT=clock)
   seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   call random_seed(PUT = seed)

   do i=inode1, inode2
      Psi(i) = cmplx_0
   end do
   Psi(PDOSAtomNumber) = cmplx_1

   cnum = 0.0_dp
   !$OMP PARALLEL DO REDUCTION(+:cnum)
   do i=1,nAt
      cnum = cnum + Psi(i)*conjg(Psi(i))
   end do
   !$OMP END PARALLEL DO
#ifdef MPI
   call MPIRedSum(cnum,s,1,MPI_DOUBLE_COMPLEX)
   cnum = s
#endif /* MPI */
   call MIO_Print('Initial norm of wave function: '//num2str(sqrt(abs(cnum)),12),'kubo')

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboInitWFPDOS',1)
#endif /* DEBUG */

end subroutine KuboInitWFPDOS

subroutine KuboInitWFLayerDOS(Psi,layerNumber,numberOfLayers,numberOfAtomsInLayer)

   use constants,            only : twopi, cmplx_i, cmplx_0, cmplx_1
   use parallel,             only : nDiv, procID
   use atoms,                only : layerIndex

   complex(dp), intent(out) :: Psi(inode1:)

   integer :: i, clock, n
   real(dp) :: r
   complex(dp) :: cnum, s
   !type(rand_t), save :: rng
   integer, pointer :: seed(:)
   integer :: SetSeed
   integer layerNumber, numberOfLayers, numberOfAtomsInLayer

#ifdef DEBUG
   call MIO_Debug('KuboInitWFLayerDOS',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   call random_seed(size = n)
   allocate(seed(n))
   call system_clock(COUNT=clock)
   call MIO_InputParameter('SeedSet',SetSeed,1235)
   call MIO_Print('Seed of the random-phase state: '//trim(num2str(SetSeed)),'kubo')
   seed = SetSeed
   call random_seed(PUT = seed)

   do i=1,nAt
      if (layerIndex(i).eq.layerNumber) then
          call random_number(r)
          Psi(i) = exp(twopi*cmplx_i*r)/sqrt(real(numberOfAtomsInLayer))
      else
          Psi(i) = cmplx_0
      end if
   end do

   cnum = 0.0_dp
   !$OMP PARALLEL DO REDUCTION(+:cnum)
   do i=1,nAt
      cnum = cnum + Psi(i)*conjg(Psi(i))
   end do
   !$OMP END PARALLEL DO
#ifdef MPI
   call MPIRedSum(cnum,s,1,MPI_DOUBLE_COMPLEX)
   cnum = s
#endif /* MPI */
   call MIO_Print('Initial norm of wave function: '//num2str(sqrt(abs(cnum)),12),'kubo')

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboInitWFLayerDOS',1)
#endif /* DEBUG */

end subroutine KuboInitWFLayerDOS

subroutine KuboInitWFLayerAndSpeciesDOS(Psi,layerNumber,speciesNumber,numberOfLayers,numberOfAtomsInLayer)

   use constants,            only : twopi, cmplx_i, cmplx_0, cmplx_1
   use parallel,             only : nDiv, procID
   use atoms,                only : layerIndex, Species

   complex(dp), intent(out) :: Psi(inode1:)

   integer :: i, clock, n
   real(dp) :: r
   complex(dp) :: cnum, s
   !type(rand_t), save :: rng
   integer, pointer :: seed(:)
   integer :: SetSeed
   integer layerNumber, numberOfLayers, numberOfAtomsInLayer, speciesNumber

#ifdef DEBUG
   call MIO_Debug('KuboInitWFLayerAndSpeciesDOS',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   call random_seed(size = n)
   allocate(seed(n))
   call system_clock(COUNT=clock)
   call MIO_InputParameter('SeedSet',SetSeed,1235)
   call MIO_Print('Seed of the random-phase state: '//trim(num2str(SetSeed)),'kubo')
   seed = SetSeed
   call random_seed(PUT = seed)

   do i=1,nAt
      if (layerIndex(i).eq.layerNumber .and. Species(i) .eq. speciesNumber) then
          call random_number(r)
          Psi(i) = exp(twopi*cmplx_i*r)/sqrt(real(numberOfAtomsInLayer))
      else
          Psi(i) = cmplx_0
      end if
   end do

   cnum = 0.0_dp
   !$OMP PARALLEL DO REDUCTION(+:cnum)
   do i=1,nAt
      cnum = cnum + Psi(i)*conjg(Psi(i))
   end do
   !$OMP END PARALLEL DO
#ifdef MPI
   call MPIRedSum(cnum,s,1,MPI_DOUBLE_COMPLEX)
   cnum = s
#endif /* MPI */
   call MIO_Print('Initial norm of wave function: '//num2str(sqrt(abs(cnum)),12),'kubo')

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboInitWFLayerAndSpeciesDOS',1)
#endif /* DEBUG */

end subroutine KuboInitWFLayerAndSpeciesDOS

subroutine KuboDOS(Psi,Psin,Psinm1,a,b,H,H0,hopp,NList,nRecurs,eps,nEn,Emin,Emax)

   use name,                 only : prefix
   use kubosubs,             only : KuboRecursion, KuboFrac, KuboFermi

   complex(dp), intent(out) :: Psi(inode1:),Psin(inode1:),Psinm1(inode1:)
   real(dp), intent(out) :: a(nRecurs), b(nRecurs)
   complex(dp), intent(out) :: H(inode1:)
   real(dp), intent(in) :: H0(inode1:)
   complex(dp), intent(in) :: hopp(:,inode1:)
   integer, intent(in) :: NList(1:,inode1:)
   integer, intent(in) :: nRecurs, nEn
   real(dp), intent(in) :: eps
   real(dp) :: Emin, Emax
   real(dp) :: ac, bc

   character(len=100) :: flnm, flnm2

   logical :: calculateInterval

#ifdef DEBUG
   call MIO_Debug('KuboDOS',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   call MIO_print('Calculating DOS by recursion','kubo')
   call KuboRecursion(Psi,Psin,Psinm1,nRecurs,a,b,H,H0,hopp,nList,.true.,trim(prefix))
   ! DOS for all energy range to obtain Fermi energy
   ! The DOS with more resoulution in the range specified
   flnm = trim(prefix)//'.DOS'
   flnm2 = trim(prefix)//'.AB'
   call KuboFrac(a,b,nRecurs,1.0_dp,eps,nEn,Emin,Emax,flnm,flnm2)
   call MIO_Print('DOS written','kubo')
   call MIO_Print('')

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboDOS',1)
#endif /* DEBUG */

end subroutine KuboDOS

subroutine KuboTEvol(Psi,ZUPsi,Psin,Psinm1,XpnPsi,XpnPsim1,c,dT,nT,az,bz,H,H0,hopp,NList,NeighD, &
                    ac,bc,nEn,Emin,Emax,eps,nRecurs,nPol,nWr)

   use constants,            only : cmplx_0
   use name,                 only : prefix
   use kubosubs,             only : KuboRecursion, KuboFrac, KuboEvol

   complex(dp), intent(inout) :: Psi(inode1:), ZUPsi(inode1:),Psin(inode1:)
   complex(dp), intent(inout) :: Psinm1(inode1:),XpnPsi(inode1:),XpnPsim1(inode1:)
   complex(dp), intent(in) :: c(nPol)
   real(dp), intent(in) :: dT, ac, bc, eps, Emin, Emax
   complex(dp), intent(out) :: H(inode1:)
   real(dp), intent(in) :: H0(inode1:), NeighD(1:,1:,inode1:)
   complex(dp), intent(in) :: hopp(:,inode1:)
   integer, intent(in) :: nT, nRecurs, nPol, nWr, NList(1:,inode1:), nEn
   real(dp), intent(inout) :: az(nRecurs), bz(nRecurs)

   logical :: prnt
   integer :: iT, i
   character(len=80) :: str, str2
   complex(dp) :: cnum, s

#ifdef DEBUG
   call MIO_Debug('KuboTEvol',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

   prnt = .false.
   call MIO_Print('')
   call MIO_Print('Time evolution of wave function','kubo')
   call MIO_Print('')
   do iT=1,nT
      if (mod(iT,nWr)==0) then
         prnt = .true.
         write(str,'(a,i5,a)') '   *** Time Step ', iT, ' ***'
         call MIO_Print(trim(str),'kubo')
      end if
      call KuboEvol(Psi,ZUPsi,Psin,Psinm1,XpnPsi,XpnPsim1,c,H0,hopp,NList,NeighD,ac,bc,nPol,prnt)
      if (prnt) then
         cnum = cmplx_0
         !$OMP PARALLEL DO REDUCTION(+:cnum)
         do i=1,nAt
            cnum = cnum + ZUPsi(i)*conjg(ZUPsi(i))
         enddo
         !$OMP END PARALLEL DO
#ifdef MPI
         call MPIAllRedSum(cnum,s,1,MPI_DOUBLE_COMPLEX)
         cnum = s
#endif /* MPI */
         !$OMP PARALLEL DO
         do i=1,nAt
            ZUPsi(i) = ZUPsi(i)/sqrt(abs(cnum))
         end do
         !$OMP END PARALLEL DO
         call KuboRecursion(ZUPsi,Psin,Psinm1,nRecurs,az,bz,H,H0,hopp,NList,.false.)
         !$OMP PARALLEL DO
         do i=1,nAt
            ZUPsi(i) = ZUPsi(i)*sqrt(abs(cnum))
         end do
         !$OMP END PARALLEL DO
         write(str,'(i0.4)') iT
         str2 = str
         str = trim(prefix)//'.Z2_'//trim(str)//'.dat'
         str2 = trim(prefix)//'.AB_'//trim(str2)//'.dat'
         call KuboFrac(az,bz,nRecurs,abs(cnum),eps,nEn,Emin,Emax,str,str2)
         prnt = .false.
      end if
   end do

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboTEvol',1)
#endif /* DEBUG */

end subroutine KuboTEvol

!!-----------------------------------------------------------------------
!
!      SUBROUTINE intervalle (a, b, ac, bc, NRECURS)
!      !INTEGER NRECURS
!!  Fait appel a tqli.f, DOnc aussi a pythag.f (Numerical Recipes)
!
!
!
!
!
!
!!  Intervalle [ac-2bc,ac+2bc] contenant tout le spectre (avec 10% de
!!  marge) :
!      !print*, "inside subroutine place 5"
!
!

end module kubo
