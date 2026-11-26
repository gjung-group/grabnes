program GenMoire

   implicit none

   !double precision :: phideg = 21.051724d0 ! Rotation angle in degrees
   !double precision :: phideg = 6.01d0 ! Rotation angle in degrees
   double precision :: phideg = 0.0d0 ! Rotation angle in degrees
   double precision :: degtol=1.0d0, lattol = 1.0d0 ! Tolerance
   double precision, parameter :: h = 3.3d0, dz = 35.0d0
   !double precision, parameter :: aG = 2.46d0, aBN0 = 6.5796819680864944d0!2.49d0 ! Lattice parameters
   double precision, parameter :: aG = 2.46d0
   double precision :: aBN0 = 2.46d0 ! Lattice parameters
   character(*), parameter :: select = 'delta'
   logical, parameter :: fdfbuild = .true.
   ! Search (MoireFind): primitive cells with at most nCellMax unit cells of the first lattice, and
   ! nMoire1 <= N(alpha - beta) <= nMoire2 (N(alpha - beta) = det(M - M'): 1 for an aligned G/hBN cell (m,0,m-1,0),
   ! (a-b)^2 for a twisted (a,b,b,a) cell). anySign: the sign of phideg is not imposed (only |angle| is compared, as
   ! before); otherwise phideg is the signed (counterclockwise) rotation of the second lattice.
   integer :: nCellMax = 10000, nMoire1 = 1, nMoire2 = 45
   logical :: anySign = .true.
   ! Exact geometry instead of a search: nGiven = m1 m2 m3 m4 (nonzero) is used as given, or replaced by its
   ! primitive cell with reduceGiven (non-primitive cells are useful for symmetry-broken states: keep .false. then).
   integer :: nGiven(4) = [0, 0, 0, 0]
   logical :: reduceGiven = .false.

   double precision, parameter :: basis(2,2) = reshape((/1.0d0/3.0d0,1.0d0/3.0d0, 2.0d0/3.0d0,2.0d0/3.0d0/),(/2,2/))
   double precision, parameter :: pi = acos(-1.0d0)
   character(*), parameter :: grSp(2) = (/'C','C'/), bnSp(2) =(/'B','N'/)
   integer, parameter :: grN(2) = (/6,6/), bnN(2) = (/5,7/)

   double precision :: delta, deltol, phi, g, lambda
   integer :: n(4)
   double precision, pointer :: Xg(:,:), Xbn(:,:)
   integer :: ncell, i
   character(len=32) :: arg

   ! The values above can also be given on the command line (no recompilation):
   !   GenMoire search phideg aBN0 degtol lattol [nCellMax nMoire1 nMoire2 [signed]]
   !   GenMoire cell m1 m2 m3 m4 [reduce]
   if (command_argument_count() > 0) then
      call get_command_argument(1, arg)
      if (arg == 'search') then
         call RealArg(2, phideg); call RealArg(3, aBN0); call RealArg(4, degtol); call RealArg(5, lattol)
         if (command_argument_count() >= 8) then
            call IntArg(6, nCellMax); call IntArg(7, nMoire1); call IntArg(8, nMoire2)
         end if
         if (command_argument_count() >= 9) anySign = .false.
      else if (arg == 'cell') then
         do i=1,4
            call IntArg(i+1, nGiven(i))
         end do
         reduceGiven = command_argument_count() >= 6
      else
         stop 'usage: GenMoire [search phideg aBN0 degtol lattol [nCellMax nMoire1 nMoire2 [signed]] | cell m1 m2 m3 m4 [reduce]]'
      end if
   end if

   delta = aBN0/aG - 1.0d0
   deltol = (aBN0 + lattol)/aG - 1.0d0 - delta
   if (any(nGiven /= 0)) then
      n = nGiven
      call MoireReport(n, reduceGiven)
   else
      call MoireFind(nCellMax,nMoire1,nMoire2,phideg,delta,degtol,deltol,select,anySign,n)
   end if
   write(*,'(a,4(1x,a))') 'Gendata.in: Structure.MoireCellParameters', (trim(IntIn(n(i))), i=1,4)
   if (fdfbuild) then
      g = n(1)**2 + n(2)**2 + n(1)*n(2)
      delta = sqrt(real(n(3)**2 + n(4)**2 + n(3)*n(4))/g)
      phi = acos((2.0d0*n(1)*n(3)+2.0d0*n(2)*n(4) + n(1)*n(4) + n(2)*n(3))/(2.0d0*delta*g))
      if (delta<1.0d0) then
         delta = 1.0d0/delta
      else
         i = n(1)
         n(1) = n(3)
         n(3) = i
         i = n(2)
         n(2) = n(4)
         n(4) = i
      end if
      write(*,*)
      write(*,'(a,f13.6)') 'Final angle:', phi*180.0d0/pi
      write(*,'(a,f13.6,a)') 'Rotation of the second layer:', SignedAngle(n), ' deg (signed, counterclockwise)' 
      write(*,'(a,f13.6)') 'Final BN lattice parameter:', aG*delta
      ncell = n(1)**2 + n(1)*n(2) + n(2)**2
      allocate(Xg(2,ncell*2))
      write(*,'(a,i8)') 'Atoms for graphene:', ncell*2
      ncell = n(3)**2 + n(3)*n(4) + n(4)**2
      allocate(Xbn(2,ncell*2))
      write(*,'(a,i8)') 'Atoms for BN:', ncell*2
      write(*,*)
      call MoireConstruct(Xg,n(1:2),basis,0.0d0,1.0d0)
      call MoireConstruct(Xbn,n(3:4),basis,-phi*180.0d0/pi,delta)
      call MoireWrite(Xg,Xbn,n,aG,-h,dz,'moire.fdf',grSp,bnSp,grN,bnN)
      delta = delta - 1.0d0
      lambda = (1.0d0+delta)*aG/sqrt(2.0d0*(1.0d0+delta)*(1.0d0-cos(phi))+delta**2)
      write(*,*) 'lambda1:', lambda
      lambda = aG*sqrt(1.0d0*n(1)**2+n(2)**2+n(1)*n(2))
      !lambda = aG*sqrt(1.0d0*n(1)**2+n(2)**2+n(1)*n(2))
      write(*,*) 'lambda2:', lambda
      !lambda = 3.0d0*aG**2/2.0d0
      !write(*,*) 'lambda3:', lambda
   end if

contains

subroutine RealArg(k, x)
   integer, intent(in) :: k
   double precision, intent(inout) :: x
   character(len=32) :: a
   call get_command_argument(k, a)
   read(a,*) x
end subroutine RealArg

subroutine IntArg(k, n)
   integer, intent(in) :: k
   integer, intent(inout) :: n
   character(len=32) :: a
   call get_command_argument(k, a)
   read(a,*) n
end subroutine IntArg

subroutine MoireFind(nCellMax,nM1,nM2,phiin,delta,degtol,lattol,select,anySign,mout)

   ! Search of the commensurate cells of two hexagonal lattices with Hermann's indices, in exact integer arithmetic.
   ! A cell is a pair of Eisenstein integers alpha = m1 + m2 tau (first lattice), beta = m3 + m4 tau (second
   ! lattice), tau = exp(i pi/3): the same supercell vector on both lattices. The second lattice is the first one
   ! scaled by |alpha/beta| = 1 + delta and turned by theta = arg(alpha/beta). Every cell with a rotation in
   ! (-180, 180] and any lattice ratio is (u alpha, u beta) for exactly one unit u = tau^k with beta in the sector
   ! m3 > 0, m4 >= 0, so the loop below over that sector, and over all alpha (signed indices) near (1+delta) exp(i phi)
   ! beta, meets each cell once. Only PRIMITIVE cells are kept: gcd(alpha, beta) a unit in Z[tau] (a common integer
   ! factor of the four indices is not the only way to be non-primitive: (4,1,1,4) holds 3 primitive cells), so every
   ! candidate is the smallest cell of its geometry. Among the candidates within the tolerances the selected one has the
   ! smallest error of the chosen kind; equal errors go to the smaller cell. The smallest candidate is printed too.

   implicit none

   double precision, parameter :: pi = acos(-1.0d0)

   integer, intent(in) :: nCellMax, nM1, nM2
   double precision, intent(in) :: phiin, delta, degtol, lattol
   character(*), intent(in) :: select
   logical, intent(in) :: anySign
   integer, intent(out) :: mout(4)

   integer :: m3, m4, a, b, amin, amax, bmin, bmax, nfound, isg, nsg
   integer :: c(4), best(4,3), small(4)
   double precision :: zb(2), z(2), r, rmax, rmin, th, ratio, ep, ea, e, phi, tol
   double precision :: err0(3), nsmall
   double precision :: sgn(2)

   phi = phiin*pi/180.0d0
   tol = degtol*pi/180.0d0
   err0 = huge(1.0d0)
   best = 0
   small = 0
   nsmall = huge(1.0d0)
   nfound = 0
   sgn = [1.0d0, -1.0d0]
   nsg = 1
   if (anySign .and. abs(phi) > 0.0d0) nsg = 2
   write(*,'(a,f16.9)') 'Angle:', phiin
   write(*,'(a,f16.9)') 'Delta:', delta
   rmax = 1.0d0 + delta + lattol
   rmin = max(1.0d0 + delta - lattol, 1.0d-12)
   do m3=1,nint(sqrt(4.0d0*nCellMax/3.0d0)/rmin)+2; do m4=0,nint(sqrt(4.0d0*nCellMax/3.0d0)/rmin)+2
      if (rmin**2*Norm2E([m3,m4]) > nCellMax) cycle
      zb = Cart([m3,m4])
      do isg=1,nsg
         ! the alpha to try: a box around the annular sector r in [rmin, rmax], angle in [phi - tol, phi + tol]
         z = rmax*sqrt(dot_product(zb,zb))*[cos(sgn(isg)*phi+atan2(zb(2),zb(1))), sin(sgn(isg)*phi+atan2(zb(2),zb(1)))]
         r = sqrt(dot_product(zb,zb))*(rmax - rmin + 2.0d0*rmax*tol) + 2.0d0
         bmin = floor((z(2) - r)/(sqrt(3.0d0)/2.0d0)); bmax = ceiling((z(2) + r)/(sqrt(3.0d0)/2.0d0))
         do b=bmin,bmax
            amin = floor(z(1) - r - 0.5d0*b); amax = ceiling(z(1) + r - 0.5d0*b)
            do a=amin,amax
               c = [a, b, m3, m4]
               if (Norm2E(c(1:2)) == 0 .or. Norm2E(c(1:2)) > nCellMax) cycle
               if (Norm2E(c(1:2)-c(3:4)) < nM1 .or. Norm2E(c(1:2)-c(3:4)) > nM2) cycle
               ratio = sqrt(dble(Norm2E(c(1:2)))/dble(Norm2E(c(3:4))))
               th = SignedAngle(c)*pi/180.0d0
               ep = abs(ratio - 1.0d0 - delta)
               if (anySign) then
                  ea = abs(abs(th) - abs(phi))
               else
                  ea = abs(th - phi)
               end if
               if (ep >= lattol .or. ea >= tol) cycle
               if (GcdNorm(c) /= 1) cycle                ! not primitive: its primitive cell is met elsewhere
               nfound = nfound + 1
               e = sqrt(ep**2/max(delta**2,epsilon(1.0d0)) + ea**2/max(phi**2,epsilon(1.0d0)))
               call Keep(err0(1), best(:,1), ea, c)
               call Keep(err0(2), best(:,2), ep, c)
               call Keep(err0(3), best(:,3), e, c)
               if (Norm2E(c(1:2)) < nsmall) then
                  nsmall = Norm2E(c(1:2))
                  small = c
               end if
            end do
         end do
      end do
   end do; end do
   write(*,*)
   write(*,'(a,i0,a,i0,a)') 'Primitive cells within the tolerances: ', nfound, ' (at most ', nCellMax, &
      ' cells of the first lattice)'
   if (nfound == 0) stop 'No commensurate cell within the tolerances: widen degtol / lattol or raise nCellMax'
   write(*,'(a)') "                a     b     a'    b'       Delta            Angle   N(a-b)   #1     #2"
   call Show('For angle:  ', best(:,1))
   call Show('For delta:  ', best(:,2))
   call Show('For both:   ', best(:,3))
   call Show('Smallest:   ', small)
   if (select=='angle' .or. select=='ANGLE') then
      mout = best(:,1)
   else if (select=='delta' .or. select=='DELTA') then
      mout = best(:,2)
   else if (select=='both' .or. select=='BOTH') then
      mout = best(:,3)
   else
      stop 'select must be angle, delta or both'
   end if

end subroutine MoireFind

subroutine Keep(err0, best, err, cc)
   ! keep cc when its error is smaller, or equal with a smaller cell
   double precision, intent(inout) :: err0
   integer, intent(inout) :: best(4)
   double precision, intent(in) :: err
   integer, intent(in) :: cc(4)
   if (err < err0 - 1.0d-12 .or. (abs(err - err0) <= 1.0d-12 .and. Norm2E(cc(1:2)) < Norm2E(best(1:2)))) then
      err0 = err
      best = cc
   end if
end subroutine Keep

subroutine Show(label, cc)
   character(*), intent(in) :: label
   integer, intent(in) :: cc(4)
   double precision :: p
   p = sqrt(dble(Norm2E(cc(1:2)))/dble(Norm2E(cc(3:4))))
   write(*,'(a,4i6,2f16.9,i7,2i7)') label, cc, max(p,1.0d0/p) - 1.0d0, SignedAngle(cc), Norm2E(cc(1:2)-cc(3:4)), &
      Norm2E(cc(1:2)), Norm2E(cc(3:4))
end subroutine Show

subroutine MoireReport(n, reduce)

   ! Indices given by hand: report whether the cell is primitive and, with reduce, replace it by its primitive cell
   ! (same twist, same lattice ratio).

   implicit none

   integer, intent(inout) :: n(4)
   logical, intent(in) :: reduce

   integer(8) :: al(2), be(2), gg(2)

   if (Norm2E(n(1:2)) == 0 .or. Norm2E(n(3:4)) == 0) stop 'nGiven: (m1, m2) and (m3, m4) must both be nonzero'
   al = n(1:2)
   be = n(3:4)
   gg = EisGcd(al, be)
   write(*,'(a,4i6,a,f14.8,a)') 'Given cell:', n, '   rotation of the second layer', SignedAngle(n), ' deg'
   if (EisNorm(gg) == 1) then
      write(*,'(a)') 'The cell is primitive'
   else
      write(*,'(a,i0,a)') 'The cell is NOT primitive: it holds ', EisNorm(gg), ' primitive cells'
      al = EisDivExact(al, gg)
      be = EisDivExact(be, gg)
      write(*,'(a,4i6)') 'Primitive cell (same twist and lattice ratio):', int(al), int(be)
      if (reduce) n = [int(al), int(be)]
   end if

end subroutine MoireReport

! ------------------------------------------------------------ Eisenstein integers x(1) + x(2) tau, tau = exp(i pi/3)
integer function Norm2E(x)
   integer, intent(in) :: x(2)
   Norm2E = x(1)**2 + x(1)*x(2) + x(2)**2
end function Norm2E

function Cart(x) result(z)
   integer, intent(in) :: x(2)
   double precision :: z(2)
   z = [x(1) + 0.5d0*x(2), sqrt(3.0d0)/2.0d0*x(2)]
end function Cart

double precision function SignedAngle(c)
   ! arg(alpha/beta) in degrees: rotation of the (m3, m4) lattice with respect to the (m1, m2) lattice
   integer, intent(in) :: c(4)
   double precision :: za(2), zb(2)
   za = Cart(c(1:2))
   zb = Cart(c(3:4))
   SignedAngle = atan2(za(2)*zb(1) - za(1)*zb(2), za(1)*zb(1) + za(2)*zb(2))*180.0d0/acos(-1.0d0)
end function SignedAngle

integer function GcdNorm(c)
   integer, intent(in) :: c(4)
   integer(8) :: al(2), be(2)
   al = c(1:2)
   be = c(3:4)
   GcdNorm = int(EisNorm(EisGcd(al, be)))
end function GcdNorm

integer(8) function EisNorm(x)
   integer(8), intent(in) :: x(2)
   EisNorm = x(1)*x(1) + x(1)*x(2) + x(2)*x(2)
end function EisNorm

function EisMul(x, y) result(z)
   integer(8), intent(in) :: x(2), y(2)
   integer(8) :: z(2)
   z = [x(1)*y(1) - x(2)*y(2), x(1)*y(2) + x(2)*y(1) + x(2)*y(2)]
end function EisMul

subroutine EisDivMod(x, y, q, r)
   ! x = q y + r with N(r) < N(y)
   integer(8), intent(in) :: x(2), y(2)
   integer(8), intent(out) :: q(2), r(2)
   integer(8) :: q0(2), qq(2), rr(2)
   integer :: da, db
   q0 = nint(dble(EisMul(x, [y(1) + y(2), -y(2)]))/dble(EisNorm(y)), 8)
   q = q0
   r = x - EisMul(q0, y)
   do da=-1,1; do db=-1,1
      qq = q0 + [int(da,8), int(db,8)]
      rr = x - EisMul(qq, y)
      if (EisNorm(rr) < EisNorm(r)) then
         q = qq
         r = rr
      end if
   end do; end do
end subroutine EisDivMod

function EisDivExact(x, y) result(q)
   integer(8), intent(in) :: x(2), y(2)
   integer(8) :: q(2), r(2)
   call EisDivMod(x, y, q, r)
   if (any(r /= 0)) stop 'EisDivExact: not a divisor'
end function EisDivExact

function EisGcd(x, y) result(g)
   integer(8), intent(in) :: x(2), y(2)
   integer(8) :: g(2), a(2), b(2), q(2), r(2)
   a = x
   b = y
   do while (EisNorm(b) /= 0)
      call EisDivMod(a, b, q, r)
      a = b
      b = r
   end do
   g = a
end function EisGcd

function IntIn(i) result(s)
   ! an integer as GRABNES reads it: a negative value is written 0-n
   integer, intent(in) :: i
   character(len=16) :: s
   if (i < 0) then
      write(s,'(a,i0)') '0-', -i
   else
      write(s,'(i0)') i
   end if
end function IntIn

subroutine MoireConstruct(X,mm,bss,angle,p)

   implicit none

   double precision, parameter :: pi = acos(-1.0d0)

   double precision, intent(out) :: X(:,:)
   integer, intent(in) :: mm(2)
   double precision, intent(in) :: bss(:,:), angle, p

   integer :: ncell, nbasis, i1, i2, j, k, n1, n2, mint(2,2)
   double precision :: ucell(2,2), Rcell(2,2), a, m(2,2)
   double precision :: f(2), sq
   double precision, pointer :: basis(:,:)
   integer, save :: nlayer=0

   ncell = mm(1)**2 + mm(1)*mm(2) + mm(2)**2
   nbasis = size(bss,2)
   allocate(basis(2,nbasis))
   if (size(X) /= 2*ncell*nbasis) then
      stop 'Wrong size of X array'
   end if
   a = angle*pi/180.0d0
   ucell(:,1) = (/cos(a),sin(a)/)
   ucell(:,2) = (/cos(a+pi/3.0d0),sin(a+pi/3)/)
   ucell = ucell*p
   m(:,1) = (/mm(1),mm(2)/)
   m(:,2) = (/-mm(2),mm(1)+mm(2)/)
   Rcell = matmul(ucell,m)
   basis = matmul(ucell,bss)
   sq = sqrt(real(ncell))
   mint = m
   n1 = gcd(mint(1,1),mint(2,1))
   n2 = ncell/n1
   nlayer = nlayer  + 1
   write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'Size of layer ',nlayer,':   ',n1,' x ',n2, '   (',mm(1),' x ',mm(2),')'
   !write(*,*) ucell
   k = 0
   !write(*,*) sqrt(Rcell(1,1)**2+Rcell(2,1)**2)/sqrt(ucell(1,1)**2+ucell(2,1)**2)
   do i1 =1,n1; do i2=1,n2
      do j=1,nbasis
         k = k + 1
         X(:,k) = (i1-1)*ucell(:,1) + (i2-1)*ucell(:,2) + basis(:,j)
      end do
   end do; end do
   do i1=1,ncell*nbasis
      !print*, Rcell(1,1), X(2,i1), X(1,i1)
      ! Cramer's rule (the first cell vector may point along y for signed indices)
      f(1) = (X(1,i1)*Rcell(2,2)-X(2,i1)*Rcell(1,2))/(Rcell(1,1)*Rcell(2,2)-Rcell(1,2)*Rcell(2,1))
      f(2) = (X(2,i1)*Rcell(1,1)-X(1,i1)*Rcell(2,1))/(Rcell(1,1)*Rcell(2,2)-Rcell(1,2)*Rcell(2,1))
      X(:,i1) = mod(f,1.0d0)
      if (X(1,i1)<0.0d0) X(1,i1) = X(1,i1) + 1.0d0
      if (X(2,i1)<0.0d0) X(2,i1) = X(2,i1) + 1.0d0
   end do
   deallocate(basis)

end subroutine MoireConstruct

subroutine MoireWrite(Xs,Xo,mm,a,h,z,name,sp1,sp2,N1,N2)

   implicit none

   integer, parameter :: u=10

   double precision, intent(in) :: Xs(:,:), Xo(:,:), a, h, z
   integer, intent(in) :: mm(4), N1(:), N2(:)
   character(*), intent(in) :: name, sp1(:), sp2(:)

   double precision :: ucell(2,2), Rcell(3,3), m1(2,2), m2(2,2)
   double precision :: rz
   integer :: sz1, sz2, i, j, k, nb1, nb2
   character(60), pointer :: Sp(:)
   integer, pointer :: Nel(:), num1(:), num2(:)

   open(u,FILE=name,STATUS='replace')
   ucell(:,1) = (/a,0.0d0/)
   ucell(:,2) = (/a/2.0d0,sqrt(3.0d0)*a/2.0d0/)
   m1(:,1) = (/mm(1),mm(2)/)
   m1(:,2) = (/-mm(2),mm(1)+mm(2)/)
   m2(:,1) = (/mm(3),mm(4)/)
   m2(:,2) = (/-mm(4),mm(3)+mm(4)/)
   Rcell = 0.0d0
   Rcell(1:2,1:2) = matmul(ucell,m1)
   Rcell(3,3) = z
   sz1 = size(Xs,2)
   sz2 = size(Xo,2)
   nb1 = size(sp1)
   nb2 = size(sp2)
   allocate(Sp(nb1+nb2))
   allocate(Nel(nb1+nb2))
   allocate(num1(nb1))
   allocate(num2(nb2))
   Sp(1) = sp1(1)
   Nel(1) = N1(1)
   k = 1
   num1(1) = 1
out1:do i=2,nb1
      do j=1,k
         if (sp1(i)==Sp(j)) then
            num1(i) = j
            cycle out1
         end if
      end do
      k = k+1
      num1(i) = k
      Sp(k) = sp1(i)
      Nel(k) = N1(i)
   end do out1
out2:do i=1,nb2
      do j=1,k
         if (sp2(i)==Sp(j)) then
            num2(i) = j
            cycle out2
         end if
      end do
      k = k+1
      num2(i) = k
      Sp(k) = sp2(i)
      Nel(k) = N2(i)
   end do out2
   write(u,'(a,i0)') 'NumberOfSpecies', k
   write(u,'(a)') '%block ChemicalSpeciesLabel'
   do i=1,k
      write(u,'(2i8,4x,a)') i, Nel(i), Sp(i)
   end do
   write(u,'(a)') '%endblock ChemicalSpeciesLabel'
   write(u,'(a,i0)') 'NumberOfAtoms  ', sz1+sz2
   write(u,'(a)') 'LatticeConstant  1.00 Ang'
   write(u,'(a)') '%block LatticeVectors'
   write(u,'(3f16.9)') ((Rcell(j,i),j=1,3),i=1,3)
   write(u,'(a)') '%endblock LatticeVectors'
   write(u,*)
   write(u,'(a)') 'AtomicCoordinatesFormat     Fractional'
   write(u,'(a)') '%block AtomicCoordinatesAndAtomicSpecies'
   rz = (z/2.0d0 - h/2.0d0)/z
   do i=1,sz1
      write(u,'(3f16.9,i8)') Xs(:,i), rz, num1(mod(i+1,nb1)+1)
   end do
   rz = (z/2.0d0 + h/2.0d0)/z
   do i=1,sz2
      write(u,'(3f16.9,i8)') Xo(:,i), rz, num2(mod(i+1,nb1)+1)
   end do
   write(u,'(a)') '%endblock AtomicCoordinatesAndAtomicSpecies'

   close(u)

end subroutine MoireWrite

function gcd(i1,i2)

    integer :: gcd
    integer, intent(in) :: i1,i2
    integer :: a,b

    a = abs(i1)     ! non-negative for signed indices
    b = abs(i2)
    if (a > b) then
        gcd = a
        a = b
        b = gcd
    end if
    if (a == 0) then
        gcd = b
        return
    end if

    do
      gcd = mod(a, b)
      if (gcd == 0) exit
      a = b
      b = gcd
    end do

    gcd = b

end function gcd

end program GenMoire
