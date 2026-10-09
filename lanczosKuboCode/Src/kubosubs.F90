module kubosubs

   use mio
   use atoms,               only : in1, in2, nAt, inode1, inode2

   implicit none

   PRIVATE

   public :: KuboRecursion
   public :: KuboFrac
   public :: KuboInterval
   public :: KuboCn
   public :: KuboEvol
   public :: KuboFermi

   interface KuboUpdate
      module procedure KuboUpdate_d, KuboUpdate_z
   end interface

   !real(dp), save :: Efermi

contains

subroutine KuboRecursion(Psi,Psin,Psinm1,nRecurs,a,b,H,H0,hopp,NList,wr,name)

   use kuboarrays,           only : tempZ
   use constants,            only : cmplx_0
   use neigh,                only : maxNeigh, Nneigh

   complex(dp), intent(inout) :: Psi(inode1:),Psin(inode1:),Psinm1(inode1:)
   real(dp), intent(out) :: a(nRecurs), b(nRecurs)
   complex(dp), intent(out) :: H(inode1:)
   real(dp), intent(in) :: H0(inode1:)
   complex(dp), intent(in) :: hopp(:,inode1:)
   integer, intent(in) :: NList(1:,inode1:)
   integer, intent(in) :: nRecurs
   logical, intent(in) :: wr
   character(len=*), intent(in), optional :: name

   complex(dp) :: sumc, s, cpsi
   integer :: i, j, iRec, u
   character(len=100) :: filenm
   type(cl_file) :: file

   !integer :: numnum

#ifdef DEBUG
   call MIO_Debug('KuboRecursion',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo::Rec')
#endif /* TIMER */

   !print*, "OMP get the max threads"
   !call OMP_GET_MAX_THREADS()

   !$OMP PARALLEL
   !      write(*,*)'Thread rank: ', OMP_GET_THREAD_NUM()
   !$OMP END PARALLEL

   call KuboUpdate(Psi,tempZ)
   sumc = cmplx_0
   !$OMP PARALLEL DO REDUCTION(+:sumc) PRIVATE(cpsi)
   do i=1,nAt
      Psinm1(i) = Psi(i)
      cpsi = H0(i)*Psi(i)
      !do j=1,maxNeigh
      do j=1,Nneigh(i)
         cpsi = cpsi - hopp(j,i)*Psi(NList(j,i))
      end do
      sumc = sumc + cpsi*conjg(Psi(i))
      H(i) = cpsi
   end do
   !$OMP END PARALLEL DO
#ifdef MPI
   call MPIAllRedSum(sumc,s,1,MPI_DOUBLE_COMPLEX)
   sumc = s
#endif /* MPI */
   a(1) = real(sumc)
   sumc = cmplx_0
   !$OMP PARALLEL DO REDUCTION(+:sumc) PRIVATE(cpsi)
   do i=1,nAt
      cpsi = H(i) - a(1)*Psinm1(i)
      Psin(i) = cpsi
      sumc = sumc + cpsi*conjg(cpsi)
   end do
   !$OMP END PARALLEL DO
#ifdef MPI
   call MPIAllRedSum(sumc,s,1,MPI_DOUBLE_COMPLEX)
   sumc = s
#endif /* MPI */
   b(1) = sqrt(abs(sumc))
   !$OMP PARALLEL DO
   do i=1,nAt
      Psin(i) = Psin(i)/b(1)
   end do
   !$OMP END PARALLEL DO
   do iRec=2,nRecurs
      call KuboUpdate(Psin,tempZ)
      sumc = cmplx_0
      !$OMP PARALLEL DO REDUCTION(+:sumc) PRIVATE(cpsi)
      do i=1,nAt
         cpsi = H0(i)*Psin(i)
         do j=1,Nneigh(i)
         !do j=1,maxNeigh
            cpsi = cpsi - hopp(j,i)*Psin(NList(j,i))
         end do
         sumc = sumc + cpsi*conjg(Psin(i))
         H(i) = cpsi
      end do
      !$OMP END PARALLEL DO
#ifdef MPI
      call MPIAllRedSum(sumc,s,1,MPI_DOUBLE_COMPLEX)
      sumc = s
#endif /* MPI */
      a(iRec) = real(sumc)
      sumc = cmplx_0
      !$OMP PARALLEL DO REDUCTION(+:sumc) PRIVATE(cpsi)
      do i=1,nAt
         cpsi = H(i) - a(iRec)*Psin(i) - b(iRec-1)*Psinm1(i)
         sumc = sumc + cpsi*conjg(cpsi)
         H(i) = cpsi
      end do
      !$OMP END PARALLEL DO
#ifdef MPI
      call MPIAllRedSum(sumc,s,1,MPI_DOUBLE_COMPLEX)
      sumc = s
#endif /* MPI */
      b(iRec) = sqrt(abs(sumc))
      !$OMP PARALLEL DO
      do i=1,nAt
         Psinm1(i) = Psin(i)
         Psin(i) = H(i)/b(iRec)
      end do
      !$OMP END PARALLEL DO
   end do

   ! Write a and b
   if (wr) then
      if (present(name)) then
         filenm = trim(name)//'.b0.dat'
      else
         filenm = 'b0.dat'
      end if
      call file%Open(name=filenm,serial=.true.)
      u = file%GetUnit()
#ifdef MPI
      if (IsPNode) then
#endif /* MPI */
         do iRec=1,nRecurs
            write(u,*) b(iRec)
         end do
#ifdef MPI
      end if
#endif /* MPI */
      call file%Close()
      if (present(name)) then
         filenm = trim(name)//'.a0.dat'
      else
         filenm = 'a0.dat'
      end if
      call file%Open(name=filenm,serial=.true.)
      u =  file%GetUnit()
#ifdef MPI
      if (IsPNode) then
#endif /* MPI */
         do iRec=1,nRecurs
            write(u,*) a(iRec)
         end do
#ifdef MPI
      end if
#endif /* MPI */
      call file%Close()
   end if
#ifdef MPI
   call MPIBarrier()
#endif /* MPI */

#ifdef TIMER
   call MIO_TimerStop('kubo::Rec')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboRecursion',1)
#endif /* DEBUG */

end subroutine KuboRecursion

subroutine KuboFrac(a,b,nRecurs,norm,eps,nEn,Emin,Emax,name,name2)

   use constants,            only : cmplx_i, cmplx_1, pi

   real(dp), intent(in) :: a(nRecurs), b(nRecurs)
   integer, intent(in) :: nRecurs, nEn
   real(dp), intent(in) :: norm, eps, Emin, Emax
   character(len=*), intent(in) :: name, name2

   real(dp) :: dE, E, ainf, binf
   integer :: iE, iR, u, v
   complex(dp) :: z, g
   real(dp) :: D
   type(cl_file) :: file

#ifdef DEBUG
   call MIO_Debug('KuboFrac',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo::Frac')
#endif /* TIMER */

   dE = (Emax-Emin)/nEn
   ainf=a(nRecurs)
   binf=b(nRecurs)
   call file%Open(name=name,serial=.true.)
   u = file%GetUnit()
   call file%Open(name=name2,serial=.true.)
   v = file%GetUnit()
#ifdef MPI
   if (IsPNode) then
#endif /* MPI */
      do iE=1,nEn+1
         E = Emin + (iE-1)*dE
         z = E*cmplx_1 + eps*cmplx_i
         g = (z-ainf-sqrt((z-ainf)**2-4.0_dp*binf**2))/(2.0_dp*binf**2)
         if(aimag(g) .ge. 0) then
            g = (z-ainf+sqrt((z-ainf)**2-4.0_dp*binf**2))/(2.0_dp*binf**2)
         endif
         do iR=nRecurs-1,1,-1
            g = 1.0_dp/(z-a(iR)-b(iR)**2*g)
         enddo
         D = -aimag(g)/pi
         write(u,*) E, D*norm
      end do
      write(v,*) norm
      do iR=1,nRecurs
         write(v,*) a(iR), b(iR)
      enddo
#ifdef MPI
   end if
#endif /* MPI */
   call file%Close()
#ifdef MPI
   call MPIBarrier()
#endif /* MPI */

#ifdef TIMER
   call MIO_TimerStop('kubo::Frac')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboFrac',0)
#endif /* DEBUG */

end subroutine KuboFrac

subroutine KuboInterval(nRecurs,a,b,ac,bc)

   integer, intent(in) :: nRecurs
   real(dp), intent(in) :: a(nRecurs), b(nRecurs)
   real(dp), intent(out) :: ac, bc

   real(dp), pointer :: e(:)=>NULL(), d(:)=>NULL()
   integer :: info, i, n
   real(dp) :: emin, emax

#ifdef DEBUG
   call MIO_Debug('KuboInterval',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

#ifdef MPI
   if (IsPnode) then
#endif /* MPI */
      call MIO_Allocate(e,nRecurs-1,'e','kubosubs')
      call MIO_Allocate(d,nRecurs,'d','kubosubs')
      info = 0
      d = a
      e = b(1:nRecurs-1)
      call DSTERF(nRecurs,d,e,info)
#ifdef MPI
   end if
   call MPIBCast(info,1,MPI_INTEGER)
#endif /* MPI */
   if (info /= 0) then
      call MIO_Kill('Interval: Error in diagonalization','kubo','KuboInterval')
   end if
#ifdef MPI
   if (IsPnode) then
#endif /* MPI */
      emin = d(1)
      emax = d(1)
      do i=2,nRecurs
         emin = min(emin,d(i))
         emax = max(emax,d(i))
      enddo
      ac = (emax+emin)/2.0_dp
      bc = 1.2_dp*(emax-emin)/4.0_dp
      call MIO_Deallocate(e,'e','kubosubs')
      call MIO_Deallocate(d,'d','kubosubs')
#ifdef MPI
   end if
   call MPIBCast(ac,1,MPI_DOUBLE_PRECISION)
   call MPIBCast(bc,1,MPI_DOUBLE_PRECISION)
#endif /* MPI */
   call MIO_Print('ac: '//trim(num2str(ac,12))//', bc: '//trim(num2str(bc,12)),'kubo')

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboInterval',1)
#endif /* DEBUG */

end subroutine KuboInterval

subroutine KuboCn(dT,ac,bc,nPol,c)

   use constants,            only : cmplx_i, cmplx_0
   use name,                 only : prefix

   real(dp), intent(in) :: dT, ac, bc
   integer, intent(in) :: nPol
   complex(dp), intent(out) :: c(nPol)

   real(dp), pointer :: e(:)=>NULL(), d(:)=>NULL(), z(:,:)=>NULL()
   !real(dp), pointer :: w(:)
   integer :: i, j, info, u
   complex(dp) :: sc
   character(len=100) :: filenm
   type(cl_file) :: file

#ifdef DEBUG
   call MIO_Debug('KuboCn',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo')
#endif /* TIMER */

#ifdef MPI
   if (IsPNode) then
#endif /* MPI */
      call MIO_Allocate(e,nPol,'e','kubosubs')
      call MIO_Allocate(d,nPol,'d','kubosubs')
      call MIO_Allocate(z,(/nPol,nPol/),'z','kubosubs')
      !call MIO_Allocate(w,2*nPol-2,'w','kubosubs')
      d = ac
      e = bc
      e(2) = sqrt(2.0_dp)*bc
      do i = 1,nPol
         z(i,i) = 1.0_dp
      end do
      call TQL2(nPol,nPol,d,e,z,info)
      !call DSTEV('V',nPol,d,e,z,nPol,w,info)
      if (info == 0) then
         do i=1,nPol
            sc = cmplx_0
            do j=1,nPol
               sc = sc + z(1,j)*z(i,j)*exp(-cmplx_i*d(j)*dt)
            enddo
            c(i)=sc
         enddo
      end if
      !call MIO_Deallocate(w,'w','kubosubs')
      call MIO_Deallocate(z,'z','kubosubs')
      call MIO_Deallocate(d,'d','kubosubs')
      call MIO_Deallocate(e,'e','kubosubs')
#ifdef MPI
   end if
   call MPIBCast(info,1,MPI_INTEGER)
#endif /* MPI */
   if (info /= 0) then
      call MIO_Kill('Cn: Error in diagonalization','kubo','KuboCn')
   end if
#ifdef MPI
   call MPIBCast(c,nPol,MPI_DOUBLE_COMPLEX)
#endif /* MPI */
   filenm = trim(prefix)//'.cn.dat'
   call file%Open(name=filenm,serial=.true.)
   u = file%GetUnit()
#ifdef MPI
   if (IsPNode) then
#endif /* MPI */
      do i=1,nPol
         write(u,*) real(c(i)), abs(c(i))
      end do
#ifdef MPI
   end if
#endif /* MPI */
   call file%Close()
#ifdef MPI
   call MPIBarrier()
#endif /* MPI */

#ifdef TIMER
   call MIO_TimerStop('kubo')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboCn',1)
#endif /* DEBUG */

end subroutine KuboCn

subroutine KuboEvol(Psi,ZUPsi,Psin,Psinm1,XpnPsi,XpnPsim1,c,H0,hopp,NList,NeighD,ac,bc,nPol,prnt)

   use constants,            only : cmplx_0
   use kuboarrays,           only : tempZ
   use neigh,                only : maxNeigh, Nneigh

   complex(dp), intent(inout) :: Psi(inode1:),ZUPsi(inode1:),Psin(inode1:)
   complex(dp), intent(inout) :: Psinm1(inode1:),XpnPsi(inode1:),XpnPsim1(inode1:)
   complex(dp), intent(in) :: c(nPol)
   real(dp), intent(in) :: ac, bc
   real(dp), intent(in) :: H0(inode1:), NeighD(1:,1:,inode1:)
   complex(dp), intent(in) :: hopp(:,inode1:)
   integer, intent(in) :: nPol, NList(1:,inode1:)
   logical, intent(in) :: prnt

   integer :: i, j, iPol
   complex(dp) :: cnum, cnum2

   logical :: l

#ifdef DEBUG
   call MIO_Debug('KuboEvol',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('kubo::Evol')
#endif /* TIMER */

   !$OMP PARALLEL DO
   do i=1,nAt
      Psinm1(i) = ZUPsi(i)*sqrt(2.0_dp)
   end do
   !$OMP END PARALLEL DO
   call KuboUpdate(Psinm1,tempZ)
   !$OMP PARALLEL DO PRIVATE(cnum)
   do i=1,nAt
      cnum = (H0(i)-ac)*Psinm1(i)
      do j=1,Nneigh(i)
         cnum = cnum - hopp(j,i)*Psinm1(NList(j,i))
      end do
      Psin(i) = cnum/(2.0_dp*bc)
      ZUPsi(i) = ZUPsi(i)*c(1)
   end do
   !$OMP END PARALLEL DO

   do iPol=2,nPol
      call KuboUpdate(Psin,tempZ)
      !$OMP PARALLEL DO PRIVATE(cnum)
      do i=1,nAt
         cnum = (H0(i)-ac)*Psin(i)
         do j=1,Nneigh(i)
            cnum = cnum - hopp(j,i)*Psin(NList(j,i))
         end do
         Psinm1(i) = (cnum-bc*Psinm1(i))/bc
      end do
      !$OMP END PARALLEL DO
      !$OMP BARRIER
      !$OMP PARALLEL DO PRIVATE(cnum)
      do i=1,nAt
         ZUPsi(i) = ZUPsi(i) + c(iPol)*Psin(i)
         cnum = Psin(i)
         Psin(i) = Psinm1(i)
         Psinm1(i) = cnum
      end do
      !$OMP END PARALLEL DO
   end do

   !$OMP PARALLEL DO
   do i=1,nAt
      Psinm1(i) = Psi(i)*sqrt(2.0_dp)
      XpnPsim1(i) = cmplx_0
   end do
   !$OMP END PARALLEL DO
   call KuboUpdate(Psinm1,tempZ)
   call MIO_InputParameter('timeEvolutionInYDirection',l,.false.)
   !$OMP PARALLEL DO PRIVATE(cnum,cnum2)
   do i=1,nAt
      cnum = (H0(i)-ac)*Psinm1(i)
      cnum2 = cmplx_0
      do j=1,Nneigh(i)
         cnum = cnum - hopp(j,i)*Psinm1(NList(j,i))
         if (l) then
              cnum2 = cnum2 - NeighD(2,j,i)*hopp(j,i)*Psinm1(NList(j,i))
         else
              cnum2 = cnum2 - NeighD(1,j,i)*hopp(j,i)*Psinm1(NList(j,i))
         end if
      end do
      Psin(i) = cnum/(2.0_dp*bc)
      XpnPsi(i) = cnum2/(2.0_dp*bc)
      Psi(i) = c(1)*Psi(i)
   end do
   !$OMP END PARALLEL DO

   do iPol=2,nPol
      call KuboUpdate(Psin,tempZ)
      call KuboUpdate(XpnPsi,tempZ)
      !$OMP PARALLEL DO PRIVATE(cnum,cnum2)
      do i=1,nAt
         cnum = (H0(i)-ac)*Psin(i)
         cnum2 = (H0(i)-ac)*XpnPsi(i)
         do j=1,Nneigh(i)
            cnum = cnum - hopp(j,i)*Psin(NList(j,i))
            if (l) then
                cnum2 = cnum2 - NeighD(2,j,i)*hopp(j,i)*Psin(NList(j,i)) - hopp(j,i)*XpnPsi(NList(j,i))
            else
                cnum2 = cnum2 - NeighD(1,j,i)*hopp(j,i)*Psin(NList(j,i)) - hopp(j,i)*XpnPsi(NList(j,i))
            end if
         end do
         Psinm1(i) = (cnum - bc*Psinm1(i))/bc
         XpnPsim1(i) = (cnum2 - bc*XpnPsim1(i))/bc
      end do
      !$OMP END PARALLEL DO
      !$OMP BARRIER
      !$OMP PARALLEL DO PRIVATE(cnum,cnum2)
      do i=1,nAt
         Psi(i) = Psi(i) + c(iPol)*Psin(i)
         ZUPsi(i) = ZUPsi(i) + c(iPol)*XpnPsi(i)
         cnum = Psin(i)
         cnum2 = XpnPsi(i)
         Psin(i) = Psinm1(i)
         XpnPsi(i) = XpnPsim1(i)
         Psinm1(i) = cnum
         XpnPsim1(i) = cnum2
      end do
      !$OMP END PARALLEL DO
   end do
   if (prnt) then
      cnum = cmplx_0
      !$OMP PARALLEL DO REDUCTION(+:cnum)
      do i=1,nAt
         cnum = cnum + Psi(i)*conjg(Psi(i))
      end do
      !$OMP END PARALLEL DO
#ifdef MPI
      call MPIAllRedSum(cnum,cnum2,1,MPI_DOUBLE_COMPLEX)
      cnum = cnum2
#endif /* MPI */
      call MIO_Print('Psi norm: '//num2str(sqrt(abs(cnum)),12))
   end if

#ifdef TIMER
   call MIO_TimerStop('kubo::Evol')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('KuboEvol',1)
#endif /* DEBUG */

end subroutine KuboEvol

subroutine KuboUpdate_d(Array,temp)

#ifdef MPI
   use neigh,                only : rcvIndx, sndList, sndIndx, sndSz, nodeSndRcv
#endif /* MPI */

   real(dp), intent(inout) :: Array(inode1:)
   real(dp), intent(out) :: temp(:)

#ifdef MPI
   integer :: i, ip, in

   if (nProc <2) return

   do i=1,sndSz
      temp(i) = Array(sndList(i))
   end do
   do i=1,nProc-1
      in = mod(Node + i,nProc)
      ip = mod(nProc + Node - i,nProc)
      call MPISendRecv(temp(sndIndx(1,in+1):),nodeSndRcv(Node+1,in+1),in,60, &
        Array(inode2+rcvIndx(1,ip+1):),nodeSndRcv(ip+1,Node+1),ip,60,MPI_DOUBLE_PRECISION)
      call MPIBarrier()
   end do
#endif /* MPI */

end subroutine KuboUpdate_d
subroutine KuboUpdate_z(Array,temp)

#ifdef MPI
   use neigh,                only : rcvIndx, sndList, sndIndx, sndSz, nodeSndRcv
#endif /* MPI */

   complex(dp), intent(inout) :: Array(inode1:)
   complex(dp), intent(out) :: temp(:)

#ifdef MPI
   integer :: i, in, ip

   if (nProc <2) return

   do i=1,sndSz
      temp(i) = Array(sndList(i))
   end do
   call MPIBarrier()
   do i=1,nProc-1
      in = mod(Node + i,nProc)
      ip = mod(nProc + Node - i,nProc)
      call MPISendRecv(temp(sndIndx(1,in+1):),nodeSndRcv(Node+1,in+1),in,60, &
        Array(inode2+rcvIndx(1,ip+1):),nodeSndRcv(ip+1,Node+1),ip,60,MPI_DOUBLE_COMPLEX)
      call MPIBarrier()
   end do
#endif /* MPI */

end subroutine KuboUpdate_z

subroutine KuboFermi(npts)

   use atoms,                only : nEl, nAt
   use math

   integer, parameter :: intorder=5

   integer, intent(in) :: npts

   type(cl_file) :: file
   integer :: u, i
   real(dp) :: Energy(npts), DOS(npts), s, h, Ep, sp, Efermi

#ifdef DEBUG
   call MIO_Debug('KuboFermi',0)
#endif /* DEBUG */

   call file%Open(name='DAT.FERMI',serial=.true.)
   u = file%GetUnit()
   do i=1,npts
      read(u,*) Energy(i), DOS(i)
   end do
   h = Energy(2) - Energy(1)
   s = TrapezoidalInt(DOS,npts,h,intorder)
   DOS = 2.0_dp*nAt*DOS/s
   Ep = 0.0_dp
   sp = 0.0_dp
   do i=1,npts-2*intorder
      s = TrapezoidalInt(DOS(:intorder*2+i),intorder*2+i,h,intorder)
      if (s>=nEl) then
         Efermi = (Energy(intorder*2+i)+Ep)/2.0_dp
         exit
      else
         if (s/=sp) then
            Ep = Energy(intorder*2+i)
            sp = s
         end if
      end if
   end do
   call MIO_Print('Efermi: '//trim(num2str(Efermi,5)),'kubo')
   call MIO_Print('')
   call file%Close()

#ifdef DEBUG
   call MIO_Debug('KuboFermi',1)
#endif /* DEBUG */

end subroutine KuboFermi

subroutine TQL2(NM,N,D,E,Z,IER)
!-------------------------------------------------------------------------
!     QL METHOD TO DETERMINE THE EIGENVALUES AND EIGENVECTORS OF:
!
!       1)  A SYMMETRIC TRIDIAGONAL MATRIX.
!       2)  A FULL SYMMETRIC MATRIX AFTER A PREVIOUS CALL TO TRED2.
!
!     CALLING MODE:
!               CALL TQL2(NM,N,D,E,Z,IER)
!     INPUTSS:
!     NM  (I4)  1ST DIMENSION OF MATRICES A AND Z IN CALLING PROGRAM
!     N   (I4)  SIZE OF Z
!     D  (R*8)  MAIN DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
!     E  (R*8)  SUB-DIAGONAL (N) OF THE TRIDIAGONAL MATRIX
!     Z  (R*8)  TABLE (NM,N) STORING THE UNITY MATRIX IF THE TRIDIAGONAL
!               MATRIX IS DEFINED BY D AND E, CASE #1.
!               FOR CASE #2, IT CONTAINS THE ELEMENTS OF THE TRANSFORMATION
!               MATRIX AFTER A CALL TO TRED2.
!     OUTPUTS:
!     D  (R*8)  EIGENVALUES
!     Z  (R*8)  EIGENVECTORS
!     IER (I4)  ERROR CODE = 0,  CONVERGENCE OK.
!                          = L,  NO CONVERGENCE FOR THE Lth EIGENVALUE
!
!     REFERENCE:
!     J.H.WILKINSON,-C.REINSCH,R.S.MARTIN
!     HANDBOOK FOR AUTOMATIC COMPUTATION, VOL.2, LINEAR ALGEBRA
!     SPRINGER-VERLAG 1971.
!-------------------------------------------------------------------------
      real(dp) :: D(N),E(N),Z(NM,N),B,C,F,G,H,P,R,S,EPS,EPS1
      INTEGER :: I,J,K,L,M,N,NM,JM, IER
      DATA EPS /0.D0/,JM /30/
      IER = 0
      IF (N.EQ.1) GO TO 38
!
!     MACHINE EPSILON
!
      IF (EPS.NE.0.D0) GO TO 12
      EPS = 1.D0
   10 EPS = EPS/2.D0
      EPS1 = 1.D0+EPS
      IF (EPS1.GT.1.D0) GO TO 10
!
   12 DO 14 I = 2,N
   14 E(I-1) = E(I)
      E(N) = 0.D0
      F = 0.D0
      B = 0.D0
!
      DO 28 L = 1,N
      J = 0
      H = EPS*(ABS(D(L))+ABS(E(L)))
      IF (B.LT.H) B = H
!
!     SEEK SMALLEST ELEMENT OF SUBDIAGONAL
!
      DO 16 M = L,N
      IF (ABS(E(M)).LE.B) GO TO 18
   16 CONTINUE
   18 IF (M.EQ.L) GO TO 26

!     START ITERATION

   20 IF (J.EQ.JM) GO TO 36
      J = J+1

!     SHIFT

      G = D(L)
      P = (D(L+1)-G)/(2.D0*E(L))
      R = SQRT(P*P+1.D0)
      D(L) = E(L)/(P+SIGN(R,P))
      H = G-D(L)
      DO 22 I = L+1,N
   22 D(I) = D(I)-H
      F = F+H

!     QL TRANSFORMATION

      P = D(M)
      C = 1.D0
      S = 0.D0
      DO 24 I = M-1,L,-1
      G = C*E(I)
      H = C*P
      IF (ABS(P).GE.ABS(E(I))) THEN
      C = E(I)/P
      R = SQRT(C*C+1.D0)
      E(I+1) = S*P*R
      S = C/R
      C = 1.D0/R
      ELSE
      C = P/E(I)
      R = SQRT(C*C+1.D0)
      E(I+1) = S*E(I)*R
      S = 1.D0/R
      C = C*S
      ENDIF
      P = C*D(I)-S*G
      D(I+1) = H+S*(C*G+S*D(I))

!     ELEMENTS OF EIGENVECTORS

      DO 24 K = 1,N
      H = Z(K,I+1)
      Z(K,I+1) = S*Z(K,I)+C*H
      Z(K,I) = Z(K,I)*C-S*H
   24 CONTINUE
      E(L) = S*P
      D(L) = C*P
      IF (ABS(E(L)).GT.B) GO TO 20

!     CONVERGENCE

   26 D(L) = D(L)+F
   28 CONTINUE

!     SORT EIGENVALUES AND EIGENVECTORS
!     IN ASVENDING ORDER

      DO 34 L = 2,N
      I = L-1
      K = I
      P = D(I)
      DO 30 J = L,N
      IF (D(J).GE.P) GO TO 30
      K = J
      P = D(J)
   30 CONTINUE
      IF (K.EQ.I) GO TO 34
      D(K) = D(I)
      D(I) = P
      DO 32 J = 1,N
      P = Z(J,I)
      Z(J,I) = Z(J,K)
   32 Z(J,K) = P
   34 CONTINUE
      GO TO 38

!     NO CONVERGENCE

   36 IER = L
   38 RETURN

end subroutine TQL2

end module kubosubs
