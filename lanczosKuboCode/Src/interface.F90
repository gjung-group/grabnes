module interface

   use mio

   implicit none

   PRIVATE

   integer, parameter :: maxEdgeNeigh=30

   integer, save :: NcellInt(2)
   real(dp), save :: lambda(2,4), AmpB(4), AmpN(4), rmaxInt
   real(dp), save :: Ampeh, dampeh

   logical, public, save :: edgeHopp
   integer, public, pointer :: nEdgeN(:), edgeIndx(:), NeI(:,:), NedgeCell(:,:,:)
   real(dp), public, pointer :: edgeH(:,:)
   integer, public, save :: nQ

   public :: InterfacePot1
   public :: InterfacePot2

contains

subroutine InterfacePot1()

   use cell,                 only : ucell
   use tbpar,                only : g0

   integer, parameter :: nmax=200
   real(dp), parameter :: Pmin=0.01_dp, tol=0.0001_dp

   real(dp) :: Amp, f, fp, l
   integer :: i

#ifdef DEBUG
   call MIO_Debug('InterfacePot1',0)
#endif /* DEBUG */

   call MIO_InputParameter('InterfaceClambda',lambda(1,1),6.78_dp)
   call MIO_InputParameter('InterfaceBNlambda',lambda(1,3),12.56_dp)
   call MIO_InputParameter('InterfaceClambda_B',lambda(1,1),lambda(1,1))
   call MIO_InputParameter('InterfaceClambda_N',lambda(2,1),lambda(1,1))
   call MIO_InputParameter('InterfaceBNlambda_B',lambda(1,3),lambda(1,3))
   call MIO_InputParameter('InterfaceBNlambda_N',lambda(2,3),lambda(1,3))
   call MIO_InputParameter('InterfaceAmpN',AmpN(3),0.56_dp)
   call MIO_InputParameter('InterfaceAmpB',AmpB(3),0.56_dp)
   call MIO_InputParameter('InterfaceAmpCN',AmpN(1),AmpN(3))
   call MIO_InputParameter('InterfaceAmpCB',AmpB(1),AmpB(3))
   call MIO_Print('Graphene-BN interface potential','ham')
   if (AmpB(1)/=AmpB(3) .or. AmpN(1)/=AmpN(3)) then
      call MIO_Print('  Amp. B for C:  '//trim(num2str(AmpB(1),3))//' eV*Ang')
      call MIO_Print('  Amp. N for C:  '//trim(num2str(AmpN(1),3))//' eV*Ang')
      call MIO_Print('  Amp. B for BN: '//trim(num2str(AmpB(3),3))//' eV*Ang')
      call MIO_Print('  Amp. N for BN: '//trim(num2str(AmpN(3),3))//' eV*Ang')
   else
      call MIO_Print('  Amp. B: '//trim(num2str(AmpB(1),3))//' eV*Ang')
      call MIO_Print('  Amp. N: '//trim(num2str(AmpN(1),3))//' eV*Ang')
   end if
   AmpB(2) = AmpB(1)
   AmpN(2) = AmpN(1)
   AmpB(4) = AmpB(3)
   AmpN(4) = AmpN(3)
   AmpB = AmpB/g0
   AmpN = AmpN/g0
   lambda(:,2) = lambda(:,1)
   lambda(:,4) = lambda(:,3)
   if (lambda(1,1)==lambda(2,1)) then
      call MIO_Print('  lambda for C:  '//trim(num2str(lambda(1,1),3))//' Ang')
   else
      call MIO_Print('  lambda for C at B edge:  '//trim(num2str(lambda(1,1),3))//' Ang')
      call MIO_Print('  lambda for C at N edge:  '//trim(num2str(lambda(2,1),3))//' Ang')
   end if
   if (lambda(1,3)==lambda(2,3)) then
      call MIO_Print('  lambda for BN: '//trim(num2str(lambda(1,3),3))//' Ang')
   else
      call MIO_Print('  lambda for BN at B edge:  '//trim(num2str(lambda(1,3),3))//' Ang')
      call MIO_Print('  lambda for BN at N edge:  '//trim(num2str(lambda(2,3),3))//' Ang')
   end if
   l = maxval(lambda)
   rmaxInt = l
   Amp = maxval(AmpB)
   f = maxval(AmpN)
   Amp = max(Amp,f)
   ! Newton's method to find the cutoff radius
   do i=1,nmax
      f = Amp*exp(-rmaxInt/l)/rmaxInt-Pmin
      fp = -Amp*(exp(rmaxInt/l)/rmaxInt**2 + exp(-rmaxInt/l)/(l*rmaxInt))
      rmaxInt = rmaxInt - f/fp
      if (f<tol) then
         exit
      end if
   end do
   NcellInt(1) = ceiling(rmaxInt/ucell(1,1))
   NcellInt(2) = ceiling(rmaxInt/ucell(2,2))
   call MIO_InputParameter('EdgeHopping',edgeHopp,.false.)
   call MIO_InputParameter('EdgeHoppingAmp',Ampeh,0.2_dp)
   call MIO_InputParameter('EdgeHoppingDamp',dampeh,0.65_dp)

#ifdef DEBUG
   call MIO_Debug('InterfacePot1',1)
#endif /* DEBUG */

end subroutine InterfacePot1

subroutine InterfacePot2(H,Nneigh,NList,neighCell,neighD,Species,Nradii)

   use atoms,                only : inode1, inode2, in1, in2, nAt, Rat, frac, AtomsSetCart
   use neigh,                only : maxNeigh
   use cell,                 only : ucell, aG
   use constants,            only : pi
   use math

   real(dp), intent(inout) :: H(inode1:)
   integer, intent(in) :: Nneigh(inode1:), NList(1:maxNeigh,inode1:inode2), neighCell(1:3,1:maxNeigh,inode1:inode2)
   integer, intent(in) :: Species(inode1:)
   real(dp), intent(in) :: neighD(:,:,:), Nradii(:,:)

   integer :: i, j, icx, icy, ie, sz
   real(dp), pointer :: Rq(:,:)
   integer, pointer :: s(:)
   real(dp) :: v(3), d, r0

#ifdef DEBUG
   call MIO_Debug('InterfacePot2',0)
#endif /* DEBUG */

   nQ = 0
   r0 = Nradii(1,1)
   if (frac) call AtomsSetCart()
   do i=inode1,inode2
      if (Species(i)<3) then
         do j=1,Nneigh(i)
            if (Species(NList(j,i))>2) then
              d = sqrt(NeighD(1,j,i)**2+NeighD(2,j,i)**2)
              if (d<r0) then
                 nQ = nQ + 1
              end if
            end if
         end do
      end if
   end do
   call MIO_Allocate(Rq,[3,nQ],'Rq','interface')
   call MIO_Allocate(s,nQ,'s','interface')
   if (edgeHopp) then
      call MIO_Allocate(nEdgeN,nQ,'nEdgeN','interface')
      call MIO_Allocate(edgeH,[maxEdgeNeigh,nQ],'edgeH','interface')
      call MIO_Allocate(edgeIndx,nQ,'edgeIndx','interface')
      call MIO_Allocate(NeI,[maxEdgeNeigh,nQ],'NeI','interface')
      call MIO_Allocate(NedgeCell,[3,maxEdgeNeigh,nQ],'NeI','interface')
   end if
   nQ = 0
   do i=inode1,inode2
      if (Species(i)<3) then
         do j=1,Nneigh(i)
            if (Species(NList(j,i))==3) then
               d = sqrt(NeighD(1,j,i)**2+NeighD(2,j,i)**2)
               if (d<r0) then
                  nQ = nQ + 1
                  Rq(:,nQ) = Rat(:,i) + NeighD(:,j,i)/2.0_dp !*0.43_dp
                  s(nQ) = 1
                  if (edgeHopp) then
                     edgeIndx(nQ) = i
                  end if
               end if
            else if (Species(NList(j,i))==4) then
               d = sqrt(NeighD(1,j,i)**2+NeighD(2,j,i)**2)
               if (d<r0) then
                  nQ = nQ + 1
                  Rq(:,nQ) = Rat(:,i) + NeighD(:,j,i)/2.0_dp !*0.56_dp
                  s(nQ) = -1
                  if (edgeHopp) then
                     edgeIndx(nQ) = i
                  end if
               end if
            end if
         end do
      end if
   end do
   call MIO_Print('Number of interface bonds: '//trim(num2str(nQ)),'interface')
   call MIO_Print('')
   do i=inode1,inode2
      do j=1,nQ
         do icy=-NcellInt(2),NcellInt(2); do icx=-NcellInt(1),NcellInt(1)
            v = Rat(:,i) - Rq(:,j) - matmul(ucell,[icx,icy,0])
            d = norm(v)
            if (d < rmaxInt) then
               if (s(j)==1) then
                  H(i) = H(i) + AmpB(Species(i))*exp(-d/lambda(1,Species(i)))/d
               else
                  H(i) = H(i) - AmpN(Species(i))*exp(-d/lambda(2,Species(i)))/d
               end if
            end if
         end do; end do
      end do
   end do
   if (edgeHopp) then
      nEdgeN = 0
      sz = size(Nradii,1)
      do i=1,nQ
         do icy=-5*NcellInt(2),5*NcellInt(2); do icx=-5*NcellInt(1),5*NcellInt(1)
            do j=1,nQ
               v = Rat(:,edgeIndx(j)) - Rat(:,edgeIndx(i)) + matmul(ucell,[icx,icy,0])
               d = norm(v)
               if (d < 20.0_dp  .and. d > Nradii(sz,1)) then
                  nEdgeN(i) = nEdgeN(i) + 1
                  if (nEdgeN(i)>maxEdgeNeigh) call MIO_Kill('Number of edge neighbors > maxEdgeNeigh', &
                    'interface','InterfacePot2')
                  NeI(nEdgeN(i),i) = edgeIndx(j)
                  NedgeCell(:,nEdgeN(i),i) = [icx,icy,0]
                  if (s(i)==1) then
                     edgeH(nEdgeN(i),i) = -Ampeh*s(i)*cos(pi*d/aG)*exp(-(d-aG)/dampeh)
                  else
                     edgeH(nEdgeN(i),i) = 0.0_dp
                  end if
               end if
            end do
         end do; end do
      end do
   end if

   call MIO_Deallocate(Rq,'Rq','interface')
   call MIO_Deallocate(s,'s','interface')

#ifdef DEBUG
   call MIO_Debug('InterfacePot2',1)
#endif /* DEBUG */

end subroutine InterfacePot2

end module interface
