module gauss

   use mio

   implicit none

   PRIVATE

!#ifdef MPI
!#define USE_MPI 1
!#include <sprng_f.h>
!#endif /* MPI */

   public :: GaussPot, GaussHeight, GaussPotDefinedPositions

contains

subroutine GaussPot(H)
   use random,               only : RandSeedFromInput

   use atoms,                only : inode1, inode2, in1, in2, nAt, Rat
   use atoms,                only : AtomsSetCart, frac
   use cell,                 only : ucell
   use parallel,             only : nDiv, procID

   real(dp), intent(inout) :: H(inode1:)

   real(dp) :: per, sigma, w, rand
   real(dp) :: dist2ref, dx, dy, dist2
   integer :: nImp, i, j, totImp, ic1, ic2, n
   integer, pointer :: indxImp(:)=>NULL()
   logical, pointer :: def(:)=>NULL()
   real(dp), pointer :: Amp(:)=>NULL()
#ifdef MPI
   integer, pointer :: impNd(:)=>NULL(), displ(:)=>NULL()
#endif /* MPI */
   !type(rand_t), save :: rng
   integer, pointer :: seed(:)
   integer :: clock

#ifdef DEBUG
   call MIO_Debug('GaussPot',0)
#endif /* DEBUG */

   if (frac) call AtomsSetCart()
   call MIO_InputParameter('GaussPotPercentage',per,1.0_dp)
   call MIO_InputParameter('GaussPotRange',sigma,1.0_dp)
   call MIO_InputParameter('GaussPotStrength',w,1.0_dp)
   call MIO_Print('')
   call MIO_Print('Gaussian potential','ham')
   call MIO_Print('  percentage: '//trim(num2str(per,3))//'%')
   call MIO_Print('  range     : '//trim(num2str(sigma,3))//' Ang')
   call MIO_Print('  strength  : '//trim(num2str(w,3))//' gamma0')
   per = per/100.0_dp
   call MIO_Allocate(def,[inode1],[inode2],'def','gauss')
   def = .false.
   totImp = 0

   call RandSeedFromInput()

   nImp = nint((inode2-inode1+1)*per)
   do i=1,nImp
      do
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
#ifdef MPI
   call MIO_Allocate(impNd,[0],[nProc-1],'impNd','gauss')
   impNd(Node) = totImp
   call MPIAllGather(MPI_IN_PLACE,0,impNd,1,MPI_INTEGER)
   totImp = sum(impNd)
#endif /* MPI */
   per = totImp*100.0_dp/nAt
   call MIO_Print('')
   call MIO_Print('  final percentage          : '//trim(num2str(per,3))//'%')
   call MIO_Print('  total number of impurities: '//trim(num2str(totImp)))
   call MIO_Allocate(indxImp,totImp,'indxImp','gauss')
   nImp = 0
   do i=inode1,inode2
      if (def(i)) then
         nImp = nImp + 1
         indxImp(nImp) = i
      end if
   end do
#ifdef MPI
   call MIO_Allocate(displ,[0],[nProc-1],'displ','gauss')
   displ(0) = 0
   do i=1,nProc-1
      displ(i) = displ(i-1) + impNd(i-1)
   end do
   call MPIAllGatherV(MPI_IN_PLACE,0,indxImp,impNd,displ,MPI_INTEGER)
   ic1 = sum(impNd(0:Node-1)) + 1
   ic2 = ic1 + impNd(Node) - 1
#else
   ic1 = 1
   ic2 = totImp
#endif /* MPI */
   call MIO_Allocate(Amp,totImp,'Amp','gauss')
   do i=ic1,ic2
      call random_number(rand)
      Amp(i) = rand
   end do
#ifdef MPI
   call MPIAllGatherV(MPI_IN_PLACE,0,Amp,impNd,displ,MPI_DOUBLE_PRECISION)
   call MIO_Deallocate(displ,'displ','gauss')
#endif /* MPI */
   dist2ref = (4.0_dp*sigma)**2
   !$OMP PARALLEL PRIVATE(dx,dy,dist2)
   do i=in1,in2
      do j=1,totImp
         do ic1=-1,1; do ic2=-1,1
            dx = Rat(1,indxImp(j)) - Rat(1,i) + ic1*ucell(1,1) + ic2*ucell(1,2)
            dy = Rat(2,indxImp(j)) - Rat(2,i) + ic1*ucell(2,1) + ic2*ucell(2,2)
            dist2 = dx**2 + dy**2
            if (dist2 <  dist2ref) then
               H(i) = H(i) + w*(Amp(j)-0.5_dp)*exp(-dist2/(2.0_dp*sigma**2))
            end if
         end do; end do
      end do
   end do
   !$OMP END PARALLEL

   call MIO_Deallocate(Amp,'Amp','gauss')
   call MIO_Deallocate(indxImp,'indxImp','gauss')
#ifdef MPI
   call MIO_Deallocate(impNd,'impNd','gauss')
#endif /* MPI */
   call MIO_Deallocate(def,'def','gauss')

#ifdef DEBUG
   call MIO_Debug('GaussPot',1)
#endif /* DEBUG */

end subroutine GaussPot

subroutine GaussPotDefinedPositions(H)
   use random,               only : RandSeedFromInput

   use atoms,                only : inode1, inode2, in1, in2, nAt, Rat
   use atoms,                only : AtomsSetCart, frac
   use cell,                 only : ucell
   use parallel,             only : nDiv, procID

   real(dp), intent(inout) :: H(inode1:)

   real(dp) :: per, sigma, w, rand
   real(dp) :: dist2ref, dx, dy, dist2
   integer :: nImp, i, j, totImp, ic1, ic2, n
   integer, pointer :: indxImp(:)=>NULL()
   logical, pointer :: def(:)=>NULL()
   real(dp), pointer :: Amp(:)=>NULL()
#ifdef MPI
   integer, pointer :: impNd(:)=>NULL(), displ(:)=>NULL()
#endif /* MPI */
   !type(rand_t), save :: rng
   integer, pointer :: seed(:)
   integer :: clock

#ifdef DEBUG
   call MIO_Debug('GaussPot',0)
#endif /* DEBUG */

   if (frac) call AtomsSetCart()
   call MIO_InputParameter('bubbleGaussPotRange',sigma,1.0_dp)
   call MIO_InputParameter('bubbleGaussPotStrength',w,1.0_dp)
   call MIO_Print('')
   call MIO_Print('Gaussian potential','ham')
   call MIO_Print('  range     : '//trim(num2str(sigma,3))//' Ang')
   call MIO_Print('  strength  : '//trim(num2str(w,3))//' gamma0')

   call RandSeedFromInput()

   totImp = size(indxImp)
   ic1 = 1
   ic2 = totImp
   call MIO_Allocate(Amp,totImp,'Amp','gauss')
   do i=ic1,ic2
      call random_number(rand)
      Amp(i) = rand
   end do
   dist2ref = (4.0_dp*sigma)**2
   !$OMP PARALLEL PRIVATE(dx,dy,dist2)
   do i=in1,in2
      do j=1,totImp
         do ic1=-1,1; do ic2=-1,1
            dx = Rat(1,indxImp(j)) - Rat(1,i) + ic1*ucell(1,1) + ic2*ucell(1,2)
            dy = Rat(2,indxImp(j)) - Rat(2,i) + ic1*ucell(2,1) + ic2*ucell(2,2)
            dist2 = dx**2 + dy**2
            if (dist2 <  dist2ref) then
               H(i) = H(i) + w*exp(-dist2/(2.0_dp*sigma**2))
            end if
         end do; end do
      end do
   end do
   !$OMP END PARALLEL

   call MIO_Deallocate(Amp,'Amp','gauss')
!#ifdef MPI
!#endif /* MPI */

#ifdef DEBUG
   call MIO_Debug('GaussPot',1)
#endif /* DEBUG */

end subroutine GaussPotDefinedPositions

subroutine GaussHeight()
   use random,               only : RandSeedFromInput

   use neigh,                only : neighD, NList, Nneigh
   use atoms,                only : inode1, inode2, in1, in2, nAt, Rat
   use atoms,                only : AtomsSetCart, frac
   use cell,                 only : ucell
   use parallel,             only : nDiv, procID

   real(dp) :: per, sigma, w, rand
   real(dp) :: dist2ref, dx, dy, dist2
   integer :: nImp, i, j, totImp, ic1, ic2, n, jj, jjj
   integer, pointer :: indxImp(:)=>NULL()
   logical, pointer :: def(:)=>NULL()
   real(dp), pointer :: Amp(:)=>NULL()
#ifdef MPI
   integer, pointer :: impNd(:)=>NULL(), displ(:)=>NULL()
#endif /* MPI */
   !type(rand_t), save :: rng
   integer, pointer :: seed(:)
   integer :: clock
   character(len=80) :: str

#ifdef DEBUG
   call MIO_Debug('GaussHeight',0)
#endif /* DEBUG */

   if (frac) call AtomsSetCart()
   call MIO_InputParameter('GaussHeightPercentage',per,1.0_dp)
   call MIO_InputParameter('GaussHeightRange',sigma,1.0_dp)
   call MIO_InputParameter('GaussHeightStrength',w,1.0_dp)
   call MIO_Print('')
   call MIO_Print('Gaussian height','ham')
   call MIO_Print('  percentage: '//trim(num2str(per,3))//'%')
   call MIO_Print('  range     : '//trim(num2str(sigma,3))//' Ang')
   call MIO_Print('  strength  : '//trim(num2str(w,3))//' gamma0')
   per = per/100.0_dp
   def = .false.
   totImp = 0

   call RandSeedFromInput()

   nImp = nint((inode2-inode1+1)*per)
   do i=1,nImp
      do
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
#ifdef MPI
   call MIO_Allocate(impNd,[0],[nProc-1],'impNd','gauss')
   impNd(Node) = totImp
   call MPIAllGather(MPI_IN_PLACE,0,impNd,1,MPI_INTEGER)
   totImp = sum(impNd)
#endif /* MPI */
   per = totImp*100.0_dp/nAt
   call MIO_Print('')
   call MIO_Print('  final percentage          : '//trim(num2str(per,3))//'%')
   call MIO_Print('  total number of impurities: '//trim(num2str(totImp)))
   nImp = 0
   do i=inode1,inode2
      if (def(i)) then
         nImp = nImp + 1
         indxImp(nImp) = i
      end if
   end do
#ifdef MPI
   call MIO_Allocate(displ,[0],[nProc-1],'displ','gauss')
   displ(0) = 0
   do i=1,nProc-1
      displ(i) = displ(i-1) + impNd(i-1)
   end do
   call MPIAllGatherV(MPI_IN_PLACE,0,indxImp,impNd,displ,MPI_INTEGER)
   ic1 = sum(impNd(0:Node-1)) + 1
   ic2 = ic1 + impNd(Node) - 1
#else
   ic1 = 1
   ic2 = totImp
#endif /* MPI */
   call MIO_Allocate(Amp,totImp,'Amp','gauss')
   do i=ic1,ic2
      call random_number(rand)
      Amp(i) = rand
   end do
#ifdef MPI
   call MPIAllGatherV(MPI_IN_PLACE,0,Amp,impNd,displ,MPI_DOUBLE_PRECISION)
   call MIO_Deallocate(displ,'displ','gauss')
#endif /* MPI */
   dist2ref = (4.0_dp*sigma)**2   ! gives the distance for which you have to make the changes
   if (MIO_StringComp(str,'MoireEncapsulatedBilayer')) then
     !$OMP PARALLEL PRIVATE(dx,dy,dist2)
     do i=in1,in2
        if (Rat(3,i).lt.20.0_dp) then ! we first start by considering only the atoms in the bottom layer
          do j=1,totImp
             do ic1=-1,1; do ic2=-1,1
                dx = Rat(1,indxImp(j)) - Rat(1,i) + ic1*ucell(1,1) + ic2*ucell(1,2)
                dy = Rat(2,indxImp(j)) - Rat(2,i) + ic1*ucell(2,1) + ic2*ucell(2,2)
                dist2 = dx**2 + dy**2 !(the correlation function is in x only
                if (dist2 <  dist2ref) then
                   Rat(3,i) = Rat(3,i) + w*(Amp(j)-0.5_dp)*exp(-dx**2/(2.0_dp*sigma**2))
                   do jj=1,Nneigh(i)
                      if (abs(neighD(1,jj,i)) .lt. 0.01_dp .and. abs(neighD(2,jj,i)) .lt. 0.01_dp) then    ! move the position of the atoms right above the ones in the bottom layer (atom jj, top layer)
                         Rat(3,NList(jj,i)) = Rat(3,NList(jj,i)) + w*(Amp(j)-0.5_dp)*exp(-dx**2/(2.0_dp*sigma**2))
                         do jjj=1,Nneigh(NList(jj,i))
                            if (abs(neighD(1,NList(jjj,NList(jj,i)),NList(jj,i))) .lt. 0.01_dp .and. &
                                      abs(neighD(3,NList(jjj,NList(jj,i)),NList(jj,i))) .lt. 0.01_dp) then   ! move the atoms that are up from atom jj. This way, every atom jj takes care of one other atom in the top layer.
                               Rat(3,NList(jjj,NList(jj,i))) = Rat(3,NList(jjj,NList(jj,i))) &
                                        + w*(Amp(j)-0.5_dp)*exp(-dx**2/(2.0_dp*sigma**2))
                            end if
                         end do
                      end if
                   end do
                end if
             end do; end do
          end do
        end if
     end do
     !$OMP END PARALLEL
   else
     !$OMP PARALLEL PRIVATE(dx,dy,dist2)
     do i=in1,in2
        do j=1,totImp
           do ic1=-1,1; do ic2=-1,1
              dx = Rat(1,indxImp(j)) - Rat(1,i) + ic1*ucell(1,1) + ic2*ucell(1,2)
              dy = Rat(2,indxImp(j)) - Rat(2,i) + ic1*ucell(2,1) + ic2*ucell(2,2)
              dist2 = dx**2 + dy**2 !(the correlation function is in x only
              if (dist2 <  dist2ref) then
                 Rat(3,i) = Rat(3,i) + w*(Amp(j)-0.5_dp)*exp(-dx**2/(2.0_dp*sigma**2))
              end if
           end do; end do
        end do
     end do
     !$OMP END PARALLEL
   end if

   call MIO_Deallocate(Amp,'Amp','gauss')
   call MIO_Deallocate(indxImp,'indxImp','gauss')
#ifdef MPI
   call MIO_Deallocate(impNd,'impNd','gauss')
#endif /* MPI */
   call MIO_Deallocate(def,'def','gauss')

#ifdef DEBUG
   call MIO_Debug('GaussHeight',1)
#endif /* DEBUG */

end subroutine GaussHeight

end module gauss
