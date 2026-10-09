module scf

   use mio

   implicit none

   PRIVATE

   public :: SCFGetCharge, SCFInit

   complex(dp), pointer :: ZWork(:)=>NULL()
   real(dp), pointer :: DWork(:)=>NULL()
   integer :: lwork

   real(dp), public, pointer :: charge(:,:)
   real(dp), parameter, public :: Zch = 1.0_dp

contains

subroutine SCFGetCharge

   use cell,                 only : rcell, ucell
   use atoms,                only : nAt, Species, layerIndex
   use ham,                  only : H0, hopp, Ho
   use neigh,                only : NList, Nneigh, neighCell,maxNeigh

   integer :: nk(3), ptot, ik, i1, i2, i3, nIt, it, is
   real(dp), pointer :: Kgrid(:,:), Eig(:,:), chargePrev(:,:), chTemp(:,:)
   real(dp) :: qaux, deltaq, tol, alpha, dq
   complex(dp), pointer :: H(:,:,:)
   character(4) :: str
   logical :: af, aflayer, flayer

#ifdef DEBUG
   call MIO_Debug('SCFGetCharge',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('scf')
#endif /* TIMER */

   call MIO_Print('==============================================','scf')
   call MIO_Print('            SCF charge calculation','scf')
   call MIO_Print('==============================================','scf')
   call MIO_Print('')
   call MIO_InputParameter('KGrid',nk,[1,1,1])
   call MIO_InputParameter('SCFMaxIterations',nIt,100)
   call MIO_InputParameter('SCFTolerance',tol,0.0001_dp)
   call MIO_InputParameter('SCFMix',alpha,0.5_dp)
   call MIO_InputParameter('MOCharge',dq,0.0_dp)
   call MIO_InputParameter('AFOrdering',af,.true.)
   call MIO_InputParameter('AFOrderingLayer',aflayer,.true.)
   call MIO_InputParameter('FOrderingLayer',flayer,.true.)
   ptot = nk(1)*nk(2)*nk(3)
   call SCFInit(nAt)
   call MIO_Allocate(charge,[2,nAt],'charge','scf')
   call MIO_Allocate(Kgrid,[3,ptot],'Kgrid','scf')
   call MIO_Allocate(H,[nAt,nAt,2],'H','scf')
   call MIO_Allocate(Eig,[nAt,2],'Eig','scf')
   call MIO_Allocate(chargePrev,[2,nAt],'charge','scf')
   call MIO_Allocate(chTemp,[2,nAt],'chTemp','scf')
   ik = 0
   do i3=1,nk(3); do i2=1,nk(2); do i1=1,nk(1)
      ik = ik+1
      Kgrid(:,ik) = rcell(:,1)*(2*i1-nk(1)-1)/(2.0_dp*nk(1)) + &
        rcell(:,2)*(2*i2-nk(2)-1)/(2.0_dp*nk(2)) + rcell(:,3)*(2*i2-nk(3)-1)/(2.0_dp*nk(3))
   end do; end do; end do
   charge = 0.0_dp
   do ik=1,ptot
      call SCFHam0(nAt,H(:,:,1),Eig(:,1),Kgrid(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell)
      do i1=1,nAt
         qaux = 0.0_dp
         do i2=1,nAt/2
            qaux = qaux + H(i1,i2,1)**2
         enddo
         charge(1,i1) = charge(1,i1) + 2.0_dp*qaux/real(ptot,8)
      enddo
   end do
   charge(2,:) = charge(1,:)
   do i1=1,nAt
      if (af) then
         if (Species(i1)==1) then
            charge(1,i1) = charge(1,i1) + dq
            charge(2,i1) = charge(2,i1) - dq
         else if (Species(i1)==2) then
            charge(1,i1) = charge(1,i1) - dq
            charge(2,i1) = charge(2,i1) + dq
         else if (Species(i1)==3) then
            charge(1,i1) = charge(1,i1) + dq
            charge(2,i1) = charge(2,i1) - dq
         else if (Species(i1)==4) then
            charge(1,i1) = charge(1,i1) - dq
            charge(2,i1) = charge(2,i1) + dq
         end if
      else if (aflayer) then
         if (mod(layerIndex(i1),2).eq.0) then
            charge(1,i1) = charge(1,i1) + dq
            charge(2,i1) = charge(2,i1) - dq
         else
            charge(1,i1) = charge(1,i1) - dq
            charge(2,i1) = charge(2,i1) + dq
         end if
      else if (flayer) then
         if (mod(layerIndex(i1),2).eq.0) then
            charge(1,i1) = charge(1,i1) + dq
            charge(2,i1) = charge(2,i1) + dq
         else
            charge(1,i1) = charge(1,i1) + dq
            charge(2,i1) = charge(2,i1) + dq
         end if
      else
         charge(1,i1) = charge(1,i1) + dq
         charge(2,i1) = charge(2,i1) - dq
      end if
   end do
   chargePrev = 0.0_dp
   open(67,FILE="generate.charge.dump")
   do it=0,nIt
      deltaq = 0.0
      do is=1,2
         do i1=1,nAt
            deltaq = deltaq + abs(charge(is,i1) - chargePrev(is,i1))
         enddo
      enddo
      write(str,'(i4)') it
      call MIO_Print('Iteration no. '//str//'   dq='//trim(num2str(deltaq,5)),'scf')
      if (deltaq<tol) exit
      chargePrev = charge
      chTemp = 0.0_dp
      do ik=1,ptot
         do is=1,2
            call SCFHam(nAt,H(:,:,is),Eig(:,is),charge,is,Kgrid(:,ik),ucell,H0,maxNeigh,hopp,NList,Nneigh,neighCell,Species)
            do i1=1,nAt
               qaux = 0.0_dp
               do i2=1,nAt/2
                  qaux = qaux + abs(H(i1,i2,is))**2
               enddo
               chTemp(is,i1) = chTemp(is,i1) + 2.0_dp*qaux/real(ptot,8)
            enddo
         end do
      end do
      if (it<=1) then
         charge = chTemp
      else
         charge = (1.0_dp-alpha)*chTemp + alpha*chargePrev
      end if
      do i1=1,nAt
         write(67,*), it, charge(1,i1), charge(2,i1)
      end do
   end do
   !H0 = Ho
   open(66,FILE="generate.charge")
   do i1=1,nAt
      write(66,*), charge(1,i1), charge(2,i1)
   end do
   call MIO_Deallocate(chTemp,'chTemp','scf')
   call MIO_Deallocate(chargePrev,'chargePrev','scf')
   call MIO_Deallocate(Eig,'Eig','scf')
   call MIO_Deallocate(H,'H','scf')
   call MIO_Deallocate(Kgrid,'Kgrid','scf')
   call MIO_Deallocate(DWork,'DWork','scf')
   call MIO_Deallocate(ZWork,'ZWork','scf')
   call MIO_Print('')

#ifdef TIMER
   call MIO_TimerStop('scf')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('SCFGetCharge',1)
#endif /* DEBUG */

end subroutine SCFGetCharge

subroutine SCFHam0(N,H,E,K,cell,H0,maxN,hopp,NList,Nneigh,neighCell)

   use constants,             only : cmplx_i

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N)
   complex(dp), intent(out) :: H(N,N)
   real(dp), intent(out) :: E(N)
   real(dp), intent(in) :: K(3), cell(3,3), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   real(dp) :: R(3)

   H = 0.0_dp
   do i=1,N
      H(i,i) = H0(i)
      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         H(in,i) = H(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(K,R))
      end do
   end do
   call ZHEEV('V','L',N,H,N,E,ZWork,lwork,DWork,info)
   if (info/=0) then
      call MIO_Kill('Error in diagonalization','scf','SCFHam0')
   end if

end subroutine SCFHam0

subroutine SCFHam(N,H,E,charge,ns,K,cell,H0,maxN,hopp,NList,Nneigh,neighCell,Species)

   use constants,             only : cmplx_i
   use tbpar,                 only : U

   integer, intent(in) :: N, maxN, NList(maxN,N), Nneigh(N), neighCell(3,maxN,N), Species(N), ns
   complex(dp), intent(out) :: H(N,N)
   real(dp), intent(out) :: E(N)
   real(dp), intent(in) :: K(3), cell(3,3), charge(2,N), H0(N)
   complex(dp), intent(in) :: hopp(maxN,N)

   integer :: i, j, in, info
   real(dp) :: R(3), zz

   H = 0.0_dp
   do i=1,N
      H(i,i) = H0(i)
      !zz = charge(1,i)*charge(2,i) ! Zch
      if (ns==1) then
         H(i,i) = H(i,i) + U(Species(i))*(charge(2,i)-Zch)/2.0_dp
      else
         H(i,i) = H(i,i) + U(Species(i))*(charge(1,i)-Zch)/2.0_dp
      end if
      do j=1,Nneigh(i)
         in = NList(j,i)
         R = matmul(cell,neighCell(:,j,i))
         H(in,i) = H(in,i) - hopp(j,i)*exp(cmplx_i*dot_product(K,R))
      end do
   end do
   call ZHEEV('V','L',N,H,N,E,ZWork,lwork,DWork,info)
   if (info/=0) then
      call MIO_Kill('Error in diagonalization','scf','SCFHam')
   end if

end subroutine SCFHam

subroutine SCFInit(N)

   integer, intent(in) :: N

   complex(dp) :: A(1,1), OPT(1)
   real(dp) :: W(1), W2(1)
   integer :: INFO

   call ZHEEV('N','L',N,A,N,W,OPT,-1,W2,INFO)
   if (INFO /= 0) call MIO_Kill('Error in workspace query for SCF diagonalization','scf','SCFInit')
   lwork = int(OPT(1))
   call MIO_Allocate(ZWork,lwork,'ZWork','scf')
   call MIO_Allocate(DWork,3*N-2,'DWork','scf')

end subroutine SCFInit

end module scf
