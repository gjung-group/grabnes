module hybrid

   use mio

   implicit none

   PRIVATE

   public :: HybridGen

contains

subroutine HybridGen(nAt,Rat,Species,in1,in2)

   use cell,                 only : ucell, aG, area
   use constants,            only : pi

   ! Only hexagons of BN are introduced in the lattice.
   ! BN zones are centered at one hexagon, from where an area with radius r
   ! extends with BN. The radius is defined r = 3*n*a0/2, n an integer, a0 C-C distance.

   integer, intent(in) :: nAt, in1, in2
   real(dp), intent(in) :: Rat(3,nAt)
   integer, intent(inout) :: Species(nAt)

   real(dp) :: r, per, Gcell(2,2), v(2), d, dx, dy, rand
   integer :: n, i, sz, sCell, ix, iy, j, ic1, ic2, nBN
   real(dp), pointer :: c(:,:)
   integer, pointer :: site(:)
   logical :: l
   !type(rand_t), save :: rng
   integer :: clock
   integer, pointer :: seed(:)

#ifdef DEBUG
   call MIO_Debug('HybridGen',0)
#endif /* DEBUG */

   call random_seed(size = n)
   allocate(seed(n))
   call system_clock(COUNT=clock)
   seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   call random_seed(PUT = seed)
   deallocate(seed)
   call MIO_InputParameter('HybridRadius',r,10.0_dp)
   n = nint(2.0_dp*r/(sqrt(3.0_dp)*aG))
   r = sqrt(3.0_dp)*aG*n/2.0_dp
   call MIO_Print('Radius: '//trim(num2str(r,4))//' Ang','hybrid')
   call MIO_InputParameter('HybridPercentage',per,1.0_dp)
   per = per/100.0_dp
   n = nint(area*per/(pi*r**2))
   call MIO_Print('Number of hexagonal islands: '//trim(num2str(n)),'hybrid')
   call MIO_Allocate(c,[2,n],'c','hybrid')
   call MIO_InputParameter('SuperCell',sCell,1)
   call MIO_InputParameter('CellSize',sz,50)
   sz = sCell*sz
   Gcell(:,1) = [0.0_dp, aG]
   Gcell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp]
   do i=1,n
loop: do
         call random_number(rand)
         ix=int(rand*sz)
         call random_number(rand)
         iy=int(rand*sz)
         c(:,i) = ix*Gcell(:,1) + iy*Gcell(:,2)
         do j=1,i-1
            v = c(:,i) - c(:,j)
            d = sqrt(dot_product(v,v))
            if (d < r*1.1_dp) cycle loop
         end do
         exit loop
      end do loop
   end do
   call MIO_InputParameter('HybridRandomSites',l,.false.)
   if (l) then
      call MIO_Allocate(site,n,'site','hybrid')
      do i=1,n
         call random_number(rand)
         site(i) = int(2.0_dp*rand)
      end do
   end if
   nBN = 0
   !$OMP PARALLEL PRIVATE(dx,dy,d) REDUCTION(+:nBN)
   do i=in1,in2
      do j=1,n
         do ic1=-1,1; do ic2=-1,1
            dx = c(1,j) - Rat(1,i) + ic1*ucell(1,1) + ic2*ucell(1,2)
            dy = c(2,j) - Rat(2,i) + ic1*ucell(2,1) + ic2*ucell(2,2)
            d = sqrt(dx**2+dy**2)
            if (d <  r) then
               if (l) then
                  Species(i) = mod(Species(i) + site(i) + 1,2) + 3
               else
                  if (Species(i)==1) then
                     Species(i)= 3
                  else if (Species(i)==2) then
                     Species(i)= 4
                  else
                     call MIO_Kill('Wrong species to change.','hybrid','HybridGen')
                  end if
               end if
               nBN = nBN + 1
            end if
         end do; end do
      end do
   end do
   !$OMP END PARALLEL
   per = 100.0_dp*nBN/real(nAt,dp)
   call MIO_Print('Final percentage: '//trim(num2str(per,2))//'%','hybrid')
   if (l) call MIO_Deallocate(site)
   call MIO_Deallocate(c)

#ifdef DEBUG
   call MIO_Debug('HybridGen',1)
#endif /* DEBUG */

end subroutine HybridGen

end module hybrid
