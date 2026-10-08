module neigh

   use mio

   implicit none

   PRIVATE

   integer, public, pointer :: NList(:,:), Nneigh(:), neighCell(:,:,:), NList2(:,:)
   real(dp), public, pointer :: NeighD(:,:,:), Nradii(:,:)
   integer, public :: maxNeigh=3
#ifdef MPI
   ! The halo lists are only allocated for nProc>1; keep them disassociated
   ! (rather than undefined) otherwise so that associated() can be tested.
   integer, public, pointer :: rcvList(:)=>NULL(), rcvIndx(:,:)=>NULL()
   integer, public, pointer :: sndList(:)=>NULL(), sndIndx(:,:)=>NULL()
   integer, public, pointer :: nodeSndRcv(:,:)=>NULL()
   integer, public :: rcvSz, sndSz
#endif /* MPI */

   public :: NeighList, fastNNnotsquare, fastNNnotsquareSmall, fastNN, fastNNnotsquareNotRectangle, fastNNnotsquareBulk, fastNNnotsquareBulkSmall

contains

!> @brief Stop cleanly when an atom has more neighbours inside the cutoffs than
!!        the neighbour arrays (sized with maxnn) can hold, instead of writing
!!        past their end.
subroutine NeighCheckCount(NNcount,maxnn)

   integer, intent(in) :: NNcount, maxnn

   if (NNcount > maxnn) then
      call MIO_Kill('An atom has more neighbours within the cutoffs than the '// &
        'neighbour lists can hold ('//trim(num2str(maxnn))//'). Increase '// &
        'TB.NeighLevels and/or Neigh.LayerNeighbors, or reduce the cutoffs.', &
        'neigh','NeighCheckCount')
   end if

end subroutine NeighCheckCount

subroutine NeighList()

   use atoms,                only : frac, Rat, in1, in2, inode1, inode2, nAt, Species
   use atoms,                only : AtomsSetCart, indxDiv, indxNode, secAt, procAt
   use parallel,             only : xDiv, yDiv, xSubdiv, ySubdiv
   use cell,                 only : ucell, aG, aBN
   use constants,            only : pi
   use tbpar,                only : tbnn
   use math

   integer, parameter :: cellneigh(2,9) = reshape([-1,0, 1,0, 0,-1, 0,1, -1,-1, 1,-1, &
                                                  -1,1, 1,1, 0,0],[2,9])
   integer, parameter :: numN(8) = [3,6,3,6,6,6,6,3]
   real(dp), parameter :: rad(8) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp),3.0_dp/sqrt(3.0_dp),3.0_dp/sqrt(3.0_dp),2.0_dp,4.0_dp/sqrt(3.0_dp),4.0_dp/sqrt(3.0_dp)]

   real(dp) :: distFact

   real(dp) :: d, v(3), rmax, vt(3), r0, vTemp(3), rmaxTemp
   integer :: i, j, isec, insec(2,9), ncell(3,9), c1(3), c2(3)
   integer :: ic, ix, iy, k, np, inplaneNeigh, outplaneNeigh, iopl
   logical :: prnt, small, reptdCell(9), l, tester, BLInPlaneInteractionRadius
   logical :: ultrasmall
   character(len=50) :: mess, str
   real(dp) :: dx, dy
   integer :: expectedNumber, ii, n

   integer :: num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,num10,num11,num12,num13,num14
   integer :: NNcount
#ifdef MPI
   integer, pointer :: countN(:,:), cnt(:,:)
#endif /* MPI */

#ifdef DEBUG
   call MIO_Debug('NeighList',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('neigh')
#endif /* TIMER */

#ifdef TIMER
   call MIO_TimerCount('neigh::p1')
#endif /* TIMER */

   call MIO_InputParameter('TypeOfSystem',str,'Graphene')
   inplaneNeigh=sum(numN(1:tbnn))
   !print*, "inplane", inplaneNeigh
   if (MIO_StringComp(str,'Graphene_Over_BN')) then
      call MIO_InputParameter('InterlayerDistance',d,3.22_dp)
      call MIO_InputParameter('Neigh.LayerNeighbors',outplaneNeigh,2)
      maxNeigh = inplaneNeigh + outplaneNeigh
      call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')
   else if (MIO_StringComp(str,'BilayerGraphene_Over_BN')) then ! optimize this
      call MIO_InputParameter('InterlayerDistance',d,3.22_dp)
      call MIO_InputParameter('Neigh.LayerNeighbors',outplaneNeigh,2)
      maxNeigh = inplaneNeigh + outplaneNeigh*2.0_dp
      call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')
   else if (MIO_StringComp(str,'TrilayerBasedOnMoireCell')) then ! optimize this
      call MIO_InputParameter('InterlayerDistance',d,3.22_dp)
      call MIO_InputParameter('Neigh.LayerNeighbors',outplaneNeigh,2)
      maxNeigh = inplaneNeigh + outplaneNeigh*2.0_dp
      call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')
   else if (MIO_StringComp(str,'TwistedBilayerBasedOnMoireCell') .or. MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell')) then
      call MIO_InputParameter('InterlayerDistance',d,3.22_dp)
      call MIO_InputParameter('Neigh.LayerNeighbors',outplaneNeigh,2)
      maxNeigh = inplaneNeigh + outplaneNeigh
      call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')
   else if (MIO_StringComp(str,'ReadXYZ')) then
      call MIO_InputParameter('InterlayerDistance',d,3.22_dp)
      call MIO_InputParameter('Neigh.LayerNeighbors',outplaneNeigh,2)
      maxNeigh = inplaneNeigh + outplaneNeigh
      call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')
   else if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer')) then
      call MIO_InputParameter('InterlayerDistance',d,3.22_dp)
      call MIO_InputParameter('Neigh.LayerNeighbors',outplaneNeigh,1)
      maxNeigh = inplaneNeigh + outplaneNeigh
      call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')
   else
      maxNeigh = inplaneNeigh
      outplaneNeigh = 0
      call MIO_Allocate(Nradii,[tbnn,2],'Nradii','neigh')
   end if
   !print*, maxNeigh
   do i=1,tbnn
      Nradii(i,1) = rad(i)*aG*1.1_dp
   end do
   if (MIO_StringComp(str,'BoronNitride')) then
      Nradii(:,2) = Nradii(:,1)
   else
      do i=1,tbnn
         Nradii(i,2) = rad(i)*aBN*1.1_dp
      end do
   end if
   call MIO_InputParameter('Neigh.LayerDistFactor',distFact,1.0_dp) ! to increase the cutoff for outerlayer neighbors
   if (outplaneNeigh/=0) then
      Nradii(tbnn+1,:) = sqrt((aBN*1.1_dp)**2/3.0_dp + d**2)! * distFact
   end if
   rmax = maxval(Nradii)
   call MIO_InputParameter('Neigh.CompareAll',small,.false.)
   call MIO_InputParameter('Neigh.CompareAll.UltraSmall',ultrasmall,.false.) ! Use this if less than 36 atoms in the system 2*3*3*2 (BL)
   call MIO_InputParameter('Neigh.BLInPlaneInteractionRadius',BLInPlaneInteractionRadius,.false.)
   call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
   call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
   call MIO_Allocate(neighCell,[1,1,inode1],[3,maxNeigh,inode2],'neighCell','neigh')
   call MIO_Allocate(NeighD,[1,1,inode1],[3,maxNeigh,inode2],'NeighD','neigh')
   if (frac) call AtomsSetCart()
   call MIO_Print('')
   call MIO_Print('Generating neighbor list','neigh')
   if (ultrasmall) then
      !print*, "gege", r0, distFact
      !do i=1,nAt
      !!$OMP PARALLEL DO PRIVATE (v,d,ncell,i,j,ix,iy,vTemp,rmaxTemp,r0)
      do i=1,nAt
      r0 = Nradii(tbnn,(Species(i)-1)/2 + 1) ! The cut radius. Taken as the maximum distance that the
      !do i=in1,in2
         do j=1,nAt
            do ix=-3,3; do iy=-3,3
               if (i==j .and. ix==0 .and. iy==0) cycle
            !do ix=0; do iy=0
               !ix=0
               !iy=0
               !if (i==j) cycle
               ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
                                                                  ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
               if (BLInPlaneInteractionRadius) then
                  v = Rat(:,i) - Rat(:,j) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition.
                  vTemp = v
                  vTemp(3) = 0.0_dp
                  rmaxTemp = aG/sqrt(3.0_dp)*1.1_dp*distFact
               else
                  !v = Rat(:,k) + matmul(ucell,ncell(:,ic)) - Rat(:,i)
                  v = Rat(:,i) - Rat(:,j) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition.
                  vTemp = v
                  rmaxTemp = rmax
               end if
               d = norm(vTemp)
               if (d<r0 .and. d.gt.0.01_dp .and. abs(Rat(3,i)-Rat(3,j))<0.01_dp) then
                  !print*, r0, d, i, j
                  Nneigh(i) = Nneigh(i) + 1
                  !if (Nneigh(i)>maxNeigh)  then
                  !   !print*, "Number of neighbors: ", i, Nneigh(i)
                  !   call MIO_Kill('More neighbors than expected found',&
                  !     'neigh','NeighList')
                  !end if
                  NList(Nneigh(i),i) = j
                  NeighD(:,Nneigh(i),i) = -v
                  neighCell(:,Nneigh(i),i) = ncell(:,1)
               !else if (d<rmaxTemp .and. d.gt.0.0001_dp .and. abs(Rat(3,i)-Rat(3,j))>0.01_dp) then
               else if (d<rmaxTemp .and. abs(Rat(3,i)-Rat(3,j))>0.01_dp) then
                  !print*, rmax, d, i, j
                  Nneigh(i) = Nneigh(i) + 1
                  !if (Nneigh(i)>maxNeigh)  then
                  !   !print*, "Number of neighbors: ", i, Nneigh(i)
                  !   call MIO_Kill('More neighbors than expected found',&
                  !     'neigh','NeighList')
                  !end if
                  NList(Nneigh(i),i) = j
                  NeighD(:,Nneigh(i),i) = -v
                  neighCell(:,Nneigh(i),i) = ncell(:,1)
               end if
            end do; end do
         end do

      end do
      !!$OMP END PARALLEL DO
      call MIO_Print('Finished ultrasmall routine','neigh')

      num0 = 0
      num1 = 0
      num2 = 0
      num3 = 0
      num4 = 0
      num5 = 0
      num6 = 0
      num7 = 0
      num8 = 0
      num9 = 0
      num10= 0
      num11= 0
      num12= 0
      num13= 0
      num14= 0

      do i=1,nAt
         NNcount = Nneigh(i)
         IF(NNcount .eq. 0) num0 = num0 + 1
         !IF(NNcount .eq. 0) print*, "NN counter = 0", n
         IF(NNcount .eq. 1) num1 = num1 + 1
         IF(NNcount .eq. 2) num2 = num2 + 1
         !IF(NNcount .eq. 2) print*, "NN counter = 2", n
         IF(NNcount .eq. 3) num3 = num3 + 1
         IF(NNcount .eq. 4) num4 = num4 + 1
         IF(NNcount .eq. 5) num5 = num5 + 1
         IF(NNcount .eq. 6) num6 = num6 + 1
         IF(NNcount .eq. 7) num7 = num7 + 1
         IF(NNcount .eq. 8) num8 = num8 + 1
         IF(NNcount .eq. 9) num9 = num9 + 1
         IF(NNcount .eq. 10) num10 = num10 + 1
         IF(NNcount .eq. 11) num11 = num11 + 1
         IF(NNcount .eq. 12) num12 = num12 + 1
         IF(NNcount .eq. 13) num13 = num13 + 1
         IF(NNcount .eq. 14) num14 = num14 + 1

      end do

      PRINT*,' 0 neighbors: ',num0
      PRINT*,' 1 neighbors: ',num1
      PRINT*,' 2 neighbors: ',num2
      PRINT*,' 3 neighbors: ',num3
      PRINT*,' 4 neighbors: ',num4
      PRINT*,' 5 neighbors: ',num5
      PRINT*,' 6 neighbors: ',num6
      PRINT*,' 7 neighbors: ',num7
      PRINT*,' 8 neighbors: ',num8
      PRINT*,' 9 neighbors: ',num9
      PRINT*,' 10 neighbors: ',num10
      PRINT*,' 11 neighbors: ',num11
      PRINT*,' 12 neighbors: ',num12
      PRINT*,' 13 neighbors: ',num13
      PRINT*,' 14 neighbors: ',num14
      PRINT*,'>14 neighbors: ',nAt-(num0+num1+num2+num3+num4+num5+num6+num7+num8+num9+num10+num11+num12+num13+num14)

   else if (small) then
      !$OMP PARALLEL PRIVATE (v, d,ncell,i,j,ix,iy)
      do i=in1,in2
         do j=1,nAt
            do ix=-1,1; do iy=-1,1
               if (i==j .and. ix==0 .and. iy==0) cycle
               ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
               v = Rat(:,i) - Rat(:,j) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition.
                                                                  ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
               d = norm(v)
               if (d<rmax) then
                  Nneigh(i) = Nneigh(i) + 1
                  !if (Nneigh(i)>maxNeigh)  then
                  !   !print*, "Number of neighbors: ", i, Nneigh(i)
                  !   call MIO_Kill('More neighbors than expected found',&
                  !     'neigh','NeighList')
                  !end if
                  NList(Nneigh(i),i) = j
                  NeighD(:,Nneigh(i),i) = -v
                  neighCell(:,Nneigh(i),i) = ncell(:,1)
               end if
            end do; end do
         end do
      end do
      !$OMP END PARALLEL
   else
      !$OMP PARALLEL PRIVATE (isec,insec,ncell,d,v,iopl,c1,c2,ix,iy,vt,r0,rmaxTemp,vTemp)
      isec = 0
      ! Cycle over all the atoms present in this thread (and node)
      do i=in1,in2
         !call MIO_Print(trim(num2str(i))//'    '//trim(num2str(rmax)),'neigh')
         !iopl = 0  ! Number of neighbots out of plane (used for bilayers)
         iopl = 0  ! Number of neighbots out of plane (used for bilayers)
         r0 = Nradii(tbnn,(Species(i)-1)/2 + 1) ! The cut radius. Taken as the maximum distance that the
                                                !  last shell of neighbors can have
         ! Find the index of each of the 9 neighbor sections
         if (isec/=secAt(i)) then  ! If the section is the same as in the previous cycle, there's no need to look
                                   !  for the sections again
            isec = secAt(i) ! Update the section of the control index
            reptdCell = .false.
            do ic=1,9
               call NeighIndxCell(isec,cellneigh(:,ic),insec(:,ic),ncell(:,ic), &
                 xDiv,yDiv,xSubdiv,ySubdiv,procAt(i),secAt)
               do j=1,ic-1
                  if (insec(1,ic)==insec(1,j)) then
                     reptdCell(ic) = .true.
                     exit
                  end if
               end do
            end do
         end if
         do ic=1,9
            if (reptdCell(ic)) cycle
            do k=insec(1,ic),insec(2,ic)
               if (ic==9 .and. i==k) cycle
               if (BLInPlaneInteractionRadius) then
                  v = Rat(:,k) + matmul(ucell,ncell(:,ic)) - Rat(:,i)
                  vTemp = v
                  vTemp(3) = 0.0_dp
                  !rmaxTemp = r0*distFact
                  rmaxTemp = aG/sqrt(3.0_dp)*1.1_dp*distFact
               else
                  v = Rat(:,k) + matmul(ucell,ncell(:,ic)) - Rat(:,i)
                  vTemp = v
                  rmaxTemp = rmax
               end if
               d = norm(vTemp)
               !print*, d, rmax, r0
               if (d<(rmax)) then
                  if (abs(v(3))<0.01_dp) then
                     if (d<(r0+0.1_dp)) then
                        Nneigh(i) = Nneigh(i) + 1
                        if (Nneigh(i)>maxNeigh)  then
                           !print*, "Number of neighbors: ", i, Nneigh(i)
                           call MIO_Kill('More neighbors than expected found',&
                             'neigh','NeighList')
                        end if
                        NList(Nneigh(i),i) = k
                        neighCell(:,Nneigh(i),i) = ncell(:,ic)
                        NeighD(:,Nneigh(i),i) = v
                     end if
                  else
                  !**************************************************************
                  !**************************************************************
                   if (d<(rmaxTemp)) then
                     if (iopl>=outplaneNeigh) then
                        c1 = ncell(:,ic)
                        ix = k
                        do j=1,maxNeigh
                           if (abs(NeighD(3,j,i))>0.01_dp) then
                              if(norm(NeighD(:,j,i))>d) then
                                 iy = NList(j,i)
                                 c2 = neighCell(:,j,i)
                                 vt = -NeighD(:,j,i)
                                 NList(j,i) = ix
                                 neighCell(:,j,i) = c1
                                 NeighD(:,j,i) = v
                                 d = norm(vt)
                                 v = vt
                                 c1 = c2
                                 ix = iy
                              end if
                           end if
                        end do
                     else
                        Nneigh(i) = Nneigh(i) + 1
                        NList(Nneigh(i),i) = k
                        neighCell(:,Nneigh(i),i) = ncell(:,ic)
                        NeighD(:,Nneigh(i),i) = v
                        iopl = iopl + 1
                     end if
                   end if
                  !**************************************************************
                  !**************************************************************
                  end if
               end if
            end do
         end do
      end do
      !$OMP END PARALLEL
   end if
   !print*, "Hi I'm here again"
   call MIO_InputParameter('neighborSafetyCheck',l,.false.)
   call MIO_InputParameter('neighborExpectedNumber',expectedNumber,3)
   !print*, "Hi I'm here again3"
   if (l) then
      call MIO_Print('Performing some safetycheck','neigh')
      if (frac) call AtomsSetCart()
      !$OMP PARALLEL PRIVATE (i, j, ii, v, d, tester)
      do i=in1,in2
        if(Nneigh(i).lt.expectedNumber) then
           call MIO_Print(trim(num2str(i))//'    '//trim(num2str(Nneigh(i))),'neigh')
           do j=1,nAt
             if (abs(Rat(3,i)-Rat(3,j))<0.1_dp) then ! only consider inplane atoms for this check
                do ix=-5,5; do iy=-5,5
                   if (i==j .and. ix==0 .and. iy==0) cycle ! this means, go to next part in iteration
                   ncell(:,1) = [ix,iy,0]
                   v = Rat(:,i) - Rat(:,j) - matmul(ucell,ncell(:,1))
                   d = norm(v)
                   if (d<(aG/sqrt(3.0_dp))+0.1) then
                      !print*, "dist = ", d
                      tester = .true.
                      do ii=1,Nneigh(i)
                        !print*, NList(ii,i), j
                        if (NList(ii,i) .eq. j) then
                           !print*, "already there: ", j
                           tester = .false.
                        end if
                      end do
                      if (tester) then
                        !print*, "We add atom ", j, " as neighbor to ", i
                        Nneigh(i) = Nneigh(i) + 1
                        if (Nneigh(i)>maxNeigh)  then
                           !print*, "Number of neighbors: ", i, Nneigh(i)
                           call MIO_Kill('More neighbors than expected found',&
                             'neigh','NeighList')
                        end if
                        !print*, "Number of neighbors: ", i, Nneigh(i)
                        NList(Nneigh(i),i) = j
                        NeighD(:,Nneigh(i),i) = -v
                        neighCell(:,Nneigh(i),i) = ncell(:,1)
                        tester = .true.
                      end if
                   end if
                end do; end do
             end if
           end do
        end if
      end do
      !$OMP END PARALLEL
   end if


   !print*, "Hi I'm here again2"
#ifdef TIMER
   call MIO_TimerStop('neigh::p1')
#endif /* TIMER */
   !print*, "hi I'm here"

   call MIO_InputParameter('Neigh.Print',prnt,.false.)
   if (prnt) then
      call MIO_Print('Neighbor list','neigh')
      call MIO_Print('  ____________________________________________','neigh')
      call MIO_Print(' |         |         |         |              |','neigh')
      call MIO_Print(' |  Atom 1 |  Atom 2 | # neigh |     Cell     |','neigh')
      call MIO_Print(' |_________|_________|_________|______________|','neigh')
      call MIO_Print(' |         |         |         |              |','neigh')
      do i=inode1,inode2
         do j=1,Nneigh(i)
            write(mess,'(a,i4,a,i4,a,i4,a,3i4,a)') ' |   ',i,'  |   ', &
              NList(j,i),'  |   ',j, '  |', &
              neighCell(:,j,i),'  |'
            if (Nneigh(i) /= 3) mess = trim(mess)//' *'
            call MIO_Print(mess,'neigh')
         end do
      end do
      call MIO_Print(' |_________|_________|_________|______________|','neigh')
   end if
   call MIO_Print('')

#ifdef TIMER
   call MIO_TimerCount('neigh::p2')
#endif /* TIMER */
   call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
   NList2 = NList
#ifdef MPI
   if (nProc>1) then
      call MIO_Allocate(countN,[nProc,numThreads],'countN','neigh')
      countN = 0
      !$OMP PARALLEL
      do i=in1,in2
         do j=1,maxNeigh
            if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
               do k=nProc,1,-1
                  if(NList2(j,i)>=indxNode(k)) then
                     countN(k,nThread+1) = countN(k,nThread+1)+1
                     exit
                  end if
               end do
            end if
         end do
      end do
      !$OMP END PARALLEL
      rcvSz = sum(countN)
      call MIO_Allocate(rcvIndx,[2,nProc],'rcvIndx','neigh')
      call MIO_Allocate(rcvList,rcvSz,'rcvList','neigh')
      np = 1
      do i=1,nProc
         rcvIndx(1,i) = np
         np = np + sum(countN(i,:))
         rcvIndx(2,i) = np - 1
      end do
      call MIO_Allocate(cnt,[nProc,numThreads],'cnt','neigh')
      !$OMP PARALLEL PRIVATE(np)
      do i=in1,in2
         do j=1,maxNeigh
            if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
               do k=nProc,1,-1
                  if(NList2(j,i)>=indxNode(k)) then
                     np = rcvIndx(1,k) + sum(countN(k,1:nThread)) + cnt(k,nThread+1)
                     cnt(k,nThread+1) = cnt(k,nThread+1) + 1
                     exit
                  end if
               end do
               rcvList(np) = NList2(j,i)
               NList2(j,i) = np + inode2
            end if
         end do
      end do
      !$OMP END PARALLEL
      call MIO_Allocate(nodeSndRcv,[nProc,nProc],'nodeSndRcv','neigh')
      do i=1,nProc
         nodeSndRcv(i,Node+1) = rcvIndx(2,i)-rcvIndx(1,i) + 1
      end do
      call MPIAllGather(MPI_IN_PLACE,0,nodeSndRcv,nProc,MPI_INTEGER)
      sndSz = sum(nodeSndRcv(Node+1,:))
      call MIO_Allocate(sndList,sndSz,'sndList','neigh')
      call MIO_Allocate(sndIndx,(/2,nProc/),'sndIndx','neigh')
      np = 1
      do i=1,nProc
         sndIndx(1,i) = np
         np = np + nodeSndRcv(Node+1,i)
         sndIndx(2,i) = np - 1
      end do
      do i=1,nProc-1
         j = mod(Node + i,nProc)
         k = mod(nProc + Node - i,nProc)
         call MPISendRecv(rcvList(rcvIndx(1,j+1):rcvIndx(2,j+1)),nodeSndRcv(j+1,Node+1),j,50, &
           sndList(sndIndx(1,k+1):sndIndx(2,k+1)),nodeSndRcv(Node+1,k+1),k,50,MPI_INTEGER)
         call MPISendRecv(rcvList(rcvIndx(1,k+1):rcvIndx(2,k+1)),nodeSndRcv(k+1,Node+1),k,51, &
           sndList(sndIndx(1,j+1):sndIndx(2,j+1)),nodeSndRcv(Node+1,j+1),j,51,MPI_INTEGER)
         call MPIBarrier()
      end do
      call MIO_Deallocate(cnt,'cnt','neigh')
      call MIO_Deallocate(countN,'countN','neigh')
   end if
#endif /* MPI */
#ifdef TIMER
   call MIO_TimerStop('neigh::p2')
#endif /* TIMER */
   call MIO_InputParameter('WriteDataFiles',prnt,.false.)
   if (prnt) then
      open(1,FILE='v')
      open(2,FILE='dx')
      open(3,FILE='dy')
      open(33,FILE='dz')
      open(4,FILE='pos')
      do i=1,nAt
         write(1,*) Nneigh(i)
         write(1,'(500(I8,1X))') (NList(j,i),j=1,Nneigh(i))
         write(2,'(500(F10.5,1X))') (NeighD(1,j,i),j=1,Nneigh(i))
         write(3,'(500(F10.5,1X))') (NeighD(2,j,i),j=1,Nneigh(i))
         write(33,'(500(F10.5,1X))') (NeighD(3,j,i),j=1,Nneigh(i))
         !write(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
         !write(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
         !write(33,*) (NeighD(3,j,i),j=1,Nneigh(i))
         if (Species(i)==3) then
            write(4,*) "B    ", Rat(1,i),Rat(2,i),Rat(3,i)
         else if (Species(i)==4) then
            write(4,*) "N    ", Rat(1,i),Rat(2,i),Rat(3,i)
         else
            write(4,*) "C    ", Rat(1,i),Rat(2,i),Rat(3,i)
         end if
      end do
      close(1)
      close(2)
      close(3)
      close(33)
      close(4)
   end if

#ifdef TIMER
   call MIO_TimerStop('neigh')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('NeighList',1)
#endif /* DEBUG */

end subroutine NeighList

subroutine NeighListOld()

   use atoms,                only : Rat, in1, in2, inode1, inode2, frac, nAt
   use atoms,                only : AtomsSetCart, indxDiv, indxNode
   use parallel,             only : xDiv, yDiv, procID
   use cell,                 only : ucell
   use math

   integer, pointer :: pos(:)
   real(dp) :: alatt, r0, normal(3,4), x0(3,4), d, v(3)
   integer :: i, j, ix, iy, ixn, iyn, np, k, ncell(3)
   logical :: prnt, small
   character(len=50) :: mess
#ifdef MPI
   integer, pointer :: countN(:,:), cnt(:,:)
#endif /* MPI */

#ifdef DEBUG
   call MIO_Debug('NeighListOld',0)
#endif /* DEBUG */

   call MIO_InputParameter('LatticeParameter',alatt,2.46_dp)
   r0 = alatt/sqrt(3.0_dp)
   call MIO_InputParameter('Neigh.Distance',r0,r0)
   r0 = r0*1.1_dp
   call MIO_InputParameter('Neigh.CompareAll',small,.false.)

   ! Check if atoms are near the edge of the subcells
   call MIO_Allocate(NList,(/1,inode1/),(/3,inode2/),'NList','neigh')
   call MIO_Allocate(Nneigh,(/inode1/),(/inode2/),'Nneigh','neigh')
   call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
   call MIO_Allocate(NeighD,(/1,1,inode1/),(/3,3,inode2/),'NeighD','neigh')
   if (frac) call AtomsSetCart()
   call MIO_Allocate(pos,(/inode1/),(/inode2/),'pos','neigh')
   call MIO_Print('')
   call MIO_Print('Generating neighbor list','neigh')
   if (small) then
      do i=inode1,inode2
         do j=1,nAt
            do ix=-1,1; do iy=-1,1
               if (i==j .and. ix==0 .and. iy==0) cycle
               ncell = [ix,iy,0]
               v = Rat(:,i) - Rat(:,j) - matmul(ucell,ncell)
               d = norm(v)
               if (d<r0) then
                  Nneigh(i) = Nneigh(i) + 1
                  NList(Nneigh(i),i) = j
                  NeighD(:,Nneigh(i),i) = -v
                  neighCell(:,Nneigh(i),i) = ncell
               end if
            end do; end do
         end do
      end do
   else
      normal(:,1) = CrossProd(ucell(:,1),ucell(:,3))
      normal(:,1) = normal(:,1)/norm(normal(:,1))
      normal(:,2) = CrossProd(ucell(:,2),ucell(:,3))
      normal(:,2) = normal(:,2)/norm(normal(:,2))
      normal(:,3) = normal(:,1)
      normal(:,4) = normal(:,2)
      !$OMP PARALLEL PRIVATE(ix,iy,x0,d,v,ixn,iyn,np,ncell)
      ix = mod(procID-1,xDiv)
      iy = (procID-1)/yDiv
      x0(:,1) = ucell(:,2)*iy/yDiv
      x0(:,2) = ucell(:,1)*ix/xDiv
      x0(:,3) = ucell(:,2)*(iy+1)/yDiv
      x0(:,4) = ucell(:,1)*(ix+1)/xDiv
      do i=in1,in2
         do j=1,4
            d = abs(dot_product(normal(:,j),x0(:,j)-Rat(:,i)))
            if (d<r0) then
               if (pos(i)/=0) then
                  if (j==4 .and. pos(i)==1) then
                     pos(i) = 8
                  else
                     pos(i) = pos(i) + 4
                  end if
                  exit
               else
                  pos(i) = j
               end if
            end if
         end do
      end do
      !$OMP BARRIER
      do i=in1,in2
         do j=i+1,in2
            v = Rat(:,i)-Rat(:,j)
            d = norm(v)
            if (d<r0) then
               Nneigh(i) = Nneigh(i) + 1
               NList(Nneigh(i),i) = j
               NeighD(:,Nneigh(i),i) = -v
               Nneigh(j) = Nneigh(j) + 1
               NList(Nneigh(j),j) = i
               NeighD(:,Nneigh(j),j) = v
            end if
         end do
         if (pos(i)/=0 .and. Nneigh(i)/=3) then
            do j=0,1
               if (pos(i)<5 .and. j==1) exit
               ncell = 0
               np = mod(pos(i)-1+j,4) + 1
               if (np==1) then
                  ixn = ix
                  iyn = iy - 1
               else if(np==2) then
                  ixn = ix - 1
                  iyn = iy
               else if(np==3) then
                  ixn = ix
                  iyn = iy + 1
               else if(np==4) then
                  ixn = ix + 1
                  iyn = iy
               end if
               if (ixn<0) then
                  ixn = xDiv-1
                  ncell(1) = -1
               else if (ixn>=xDiv) then
                  ixn = 0
                  ncell(1) = 1
               end if
               if (iyn<0) then
                  iyn = yDiv-1
                  ncell(2) = -1
               else if (iyn>=yDiv) then
                  iyn = 0
                  ncell(2) = 1
               end if
               if (ixn<0) then
                  ixn = xDiv-1
                  ncell(1) = -1
               else if (ixn>=xDiv) then
                  ixn = 0
                  ncell(1) = 1
               end if
               if (iyn<0) then
                  iyn = yDiv-1
                  ncell(2) = -1
               else if (iyn>=yDiv) then
                  iyn = 0
                  ncell(2) = 1
               end if
               np = iyn*xDiv + ixn + 1
               ixn = indxDiv(np)
               iyn = indxDiv(np+1) - 1
               do k=ixn,iyn
                  v = Rat(:,i) - Rat(:,k) - matmul(ucell,ncell)
                  d = norm(v)
                  if (d<r0) then
                     Nneigh(i) = Nneigh(i) + 1
                     NList(Nneigh(i),i) = k
                     neighCell(:,Nneigh(i),i) = ncell
                     NeighD(:,Nneigh(i),i) = -v
                  end if
               end do
            end do
            if (pos(i)>4 .and. Nneigh(i)/=3) then
               np = pos(i)
               if (np==5) then
                  ixn = ix - 1
                  iyn = iy - 1
               else if(np==6) then
                  ixn = ix - 1
                  iyn = iy + 1
               else if(np==7) then
                  ixn = ix + 1
                  iyn = iy + 1
               else if(np==8) then
                  ixn = ix + 1
                  iyn = iy - 1
               end if
               if (ixn<0) then
                  ixn = xDiv-1
                  ncell(1) = -1
               else if (ixn>=xDiv) then
                  ixn = 0
                  ncell(1) = 1
               end if
               if (iyn<0) then
                  iyn = yDiv-1
                  ncell(2) = -1
               else if (iyn>=yDiv) then
                  iyn = 0
                  ncell(2) = 1
               end if
               np = iyn*xDiv + ixn + 1
               ixn = indxDiv(np)
               iyn = indxDiv(np+1) - 1
               do k=ixn,iyn
                  v = Rat(:,i)-Rat(:,k) - matmul(ucell,ncell)
                  d = norm(v)
                  if (d<r0) then
                     Nneigh(i) = Nneigh(i) + 1
                     NList(Nneigh(i),i) = k
                     neighCell(:,Nneigh(i),i) = ncell
                     NeighD(:,Nneigh(i),i) = -v
                  end if
               end do
            end if
         end if
      end do
      !$OMP END PARALLEL
   end if

   call MIO_InputParameter('Neigh.Print',prnt,.false.)
   if (prnt) then
      call MIO_Print('Neighbor list','neigh')
      call MIO_Print('  ____________________________________________','neigh')
      call MIO_Print(' |         |         |         |              |','neigh')
      call MIO_Print(' |  Atom 1 |  Atom 2 | # neigh |     Cell     |','neigh')
      call MIO_Print(' |_________|_________|_________|______________|','neigh')
      call MIO_Print(' |         |         |         |              |','neigh')
      do i=inode1,inode2
         do j=1,Nneigh(i)
            write(mess,'(a,i4,a,i4,a,i4,a,3i4,a)') ' |   ',i,'  |   ', &
              NList(j,i),'  |   ',j, '  |', &
              neighCell(:,j,i),'  |'
            if (Nneigh(i) /= 3) mess = trim(mess)//' *'
            call MIO_Print(mess,'neigh')
         end do
      end do
      call MIO_Print(' |_________|_________|_________|______________|','neigh')
   end if
   call MIO_Print('')
   call MIO_Deallocate(pos,'pos','neigh')

#ifdef MPI
   if (nProc>1) then
      call MIO_Allocate(countN,(/nProc,numThreads/),'countN','neigh')
      !$OMP PARALLEL
      do i=in1,in2
         do j=1,3
            if (NList(j,i)<inode1 .or. NList(j,i)>inode2) then
               do k=nProc,1,-1
                  if(NList(j,i)>=indxNode(k)) then
                     countN(k,nThread+1) = countN(k,nThread+1)+1
                     exit
                  end if
               end do
            end if
         end do
      end do
      !$OMP END PARALLEL
      rcvSz = sum(countN)
      call MIO_Allocate(rcvIndx,(/2,nProc/),'rcvIndx','neigh')
      call MIO_Allocate(rcvList,rcvSz,'rcvList','neigh')
      np = 1
      do i=1,nProc
         rcvIndx(1,i) = np
         np = np + sum(countN(i,:))
         rcvIndx(2,i) = np - 1
      end do
      call MIO_Allocate(cnt,(/nProc,numThreads/),'cnt','neigh')
      !$OMP PARALLEL PRIVATE(np)
      do i=in1,in2
         do j=1,3
            if (NList(j,i)<inode1 .or. NList(j,i)>inode2) then
               do k=nProc,1,-1
                  if(NList(j,i)>=indxNode(k)) then
                     np = rcvIndx(1,k) + sum(countN(k,1:nThread)) + cnt(k,nThread+1)
                     cnt(k,nThread+1) = cnt(k,nThread+1) + 1
                     exit
                  end if
               end do
               rcvList(np) = NList(j,i)
               NList(j,i) = np + inode2
            end if
         end do
      end do
      !$OMP END PARALLEL
      call MIO_Allocate(nodeSndRcv,[nProc,nProc],'nodeSndRcv','neigh')
      do i=1,nProc
         nodeSndRcv(i,Node+1) = rcvIndx(2,i)-rcvIndx(1,i) + 1
      end do
      call MPIAllGather(MPI_IN_PLACE,0,nodeSndRcv,nProc,MPI_INTEGER)
      sndSz = sum(nodeSndRcv(Node+1,:))
      call MIO_Allocate(sndList,sndSz,'sndList','neigh')
      call MIO_Allocate(sndIndx,(/2,nProc/),'sndIndx','neigh')
      np = 1
      do i=1,nProc
         sndIndx(1,i) = np
         np = np + nodeSndRcv(Node+1,i)
         sndIndx(2,i) = np - 1
      end do
      do i=1,nProc-1
         j = mod(Node + i,nProc)
         k = mod(nProc + Node - i,nProc)
         call MPISendRecv(rcvList(rcvIndx(1,j+1):rcvIndx(2,j+1)),nodeSndRcv(j+1,Node+1),j,50, &
           sndList(sndIndx(1,k+1):sndIndx(2,k+1)),nodeSndRcv(Node+1,k+1),k,50,MPI_INTEGER)
         call MPISendRecv(rcvList(rcvIndx(1,k+1):rcvIndx(2,k+1)),nodeSndRcv(k+1,Node+1),k,51, &
           sndList(sndIndx(1,j+1):sndIndx(2,j+1)),nodeSndRcv(Node+1,j+1),j,51,MPI_INTEGER)
         call MPIBarrier()
      end do
      call MIO_Deallocate(cnt,'cnt','neigh')
      call MIO_Deallocate(countN,'countN','neigh')
   end if
#endif /* MPI */

#ifdef DEBUG
   call MIO_Debug('NeighListOld',1)
#endif /* DEBUG */

end subroutine NeighListOld

subroutine NeighIndxCell(isec,secd,insec,ncell, &
                 xDiv,yDiv,xSubdiv,ySubdiv,proc,secAt)

   implicit none

   integer, intent(in) :: isec, secd(2), xDiv, yDiv, xSubdiv, ySubdiv, proc, secAt(:)
   integer, intent(out) :: insec(2), ncell(3)

   integer :: div(2), is, nsub, n, indx, nAt, i

   ncell = 0
   div = 0
   nsub = xSubdiv*ySubdiv
   if (secd(1)==-1 .and. mod(isec,xSubdiv)==1) then
      div(1) = -1
      if (mod(proc-1,xDiv)+1==1) ncell(1) = -1
   else if (secd(1)==1 .and. mod(isec,xSubdiv)==0) then
      div(1) = 1
      if (mod(proc,xDiv)==0) ncell(1) = 1
   end if
   if (secd(2)==-1 .and. mod(isec-1,nsub)/xSubdiv==0) then
      div(2) = -1
      if ((proc-1)/xDiv==0) ncell(2) = -1
   else if (secd(2)==1 .and. mod(isec-1,nsub)/xSubdiv==ySubdiv-1) then
      div(2) = 1
      if ((proc-1)/xDiv==yDiv-1) ncell(2) = 1
   end if
   is = isec
   if (div(1)/=0) then
      is = is + secd(1)*nsub - secd(1)*(xSubdiv-1)
   else
      is = is + secd(1)
   end if
   if (div(2)/=0) then
      is = is + nsub*(xDiv - 1)*secd(2) + xSubdiv*secd(2)
   else
      is = is + secd(2)*xSubdiv
   end if
   if (ncell(1)/=0) is = is - nsub*xDiv*ncell(1)
   is = mod(is+nsub*xDiv*yDiv-1,nsub*xDiv*yDiv)+1
   nAt = size(secAt)
   n = nAt/(nsub*xDiv*yDiv)
   indx = is*n
   i = indx
   do
      if (secAt(i)==is) then
         if (i==1) then
            insec(1) = i
            exit
         else if (secAt(i-1)==is-1) then
            insec(1) = i
            exit
         else
            if (secAt(i-1)==is) then
               i = i-1
            else
               call MIO_Kill('Problem with atom indices for neighbors. (1)','neigh','NeighIndxCell')
            end if
         end if
      else
         if (secAt(i)>is) then
            i = i - 1
         else
            i = i + 1
         end if
         if (i < 1 .or. i > nAt) then
            call MIO_Kill('Problem with atom indices for neighbors. (2)', &
              'neigh','NeighIndxCell')
         end if
      end if
   end do
   i = indx
   do
      if (secAt(i)==is) then
         if (i==nAt) then
            insec(2) = i
            exit
         else if (secAt(i+1)==is+1) then
            insec(2) = i
            exit
         else
            if (secAt(i+1)==is) then
               i = i+1
            else
               call MIO_Kill('Problem with atom indices for neighbors. (3)','neigh','NeighIndxCell')
            end if
         end if
      else
         if (secAt(i)>is) then
            i = i - 1
         else
            i = i + 1
         end if
         if (i < 1 .or. i > nAt) call MIO_Kill('Problem with atom indices for neighbors. (4)', &
           'neigh','NeighIndxCell')
      end if
   end do

end subroutine NeighIndxCell

subroutine fastNNnotsquare(natoms,x,y,z,aCC,cutoff2,cutoff2bis,A1,A2,maxnn)

   use atoms,                only : inode1, inode2, in1, in2, indxNode, nAt, Species, Rat
   use tbpar,                only : tbnn
   use cell,                 only : aG, aBN,ucell
   use atoms,                only : AtomsSetCart, frac
   use math
   implicit none
   integer, intent(in) :: natoms,maxnn
   real(dp), intent(in) :: x(natoms),y(natoms),z(natoms)
   integer :: nn(natoms,maxnn)
   real(dp), intent(in) :: aCC,cutoff2,cutoff2bis,A1(3),A2(3)
   integer :: num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,num10,num11,num12,num13,num14
   integer :: i,j,k,l,m,n
   real(dp) :: dx(natoms,maxnn),dy(natoms,maxnn),dz(natoms,maxnn)
   !real(dp), intent(out) :: dr(natoms,maxnn)

   real(dp) :: dx_t,dy_t,dz_t,d2_t, d2_t2
   real(dp) :: xmin,xmax,ymin,ymax,xcell,ycell,rho

   real(dp) :: xnew(9*natoms),ynew(9*natoms),znew(9*natoms)
   integer :: ixs(9*natoms),iys(9*natoms),vecino

   integer :: Nx,Ny,ix,iy
   integer, allocatable :: cells(:,:,:)
   integer :: NNcount,addx,addy

   !real(dp), parameter :: rad(3) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp)]
   real(dp), parameter :: rad(5) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp),3.0_dp/sqrt(3.0_dp),3.0_dp/sqrt(3.0_dp)]
   character(len=50) :: str

   integer :: ncell(3,9)
   real(dp) :: d, v(3), vtemp(3)
   real(dp) :: rmax, r0, rmaxTemp

   integer, pointer :: countN(:,:), cnt(:,:)
   integer :: np

   integer :: safetycounter
   logical :: prnt
   logical :: BLInPlaneInteractionRadius

   real(dp) :: distFact, dist

   integer :: sCell

   logical :: readDataFiles, readNeighborDetails
   integer :: jj


   call MIO_InputParameter('TypeOfSystem',str,'Graphene')
   call MIO_InputParameter('ReadDataFiles',readDataFiles,.false.)
   !if (frac) call AtomsSetCart()
   !call AtomsSetFrac()
   if (readDataFiles) then
      call MIO_Print('Reading in the neighbors from v instead of using the fastNNnotsquare routine','neigh')
      call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')

      do i=1,tbnn
         Nradii(i,1) = rad(i)*aG*1.1_dp
      end do
      if (MIO_StringComp(str,'BoronNitride')) then
         Nradii(:,2) = Nradii(:,1)
      else
         do i=1,tbnn
            Nradii(i,2) = rad(i)*aBN*1.1_dp
         end do
      end if
      rmax = maxval(Nradii)

      maxNeigh = maxnn

      !print*, "allocating NList with maxNeigh= ", maxNeigh
      call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
      !call MIO_Allocate(NList,[1,inode1],[maxnn,inode2],'NList','neigh')
      call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
      call MIO_Allocate(NeighD,[1,1,inode1],[3,maxnn,inode2],'NeighD','neigh')
      !call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
      call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,maxnn,inode2/),'neighCell','neigh')
      open(1,FILE='v')
      open(2,FILE='dx')
      open(3,FILE='dy')
      !open(4,FILE='pos')
      !print*, "nAt =", nAt
      do i=1,nAt
         read(1,*) Nneigh(i)
         read(1,*) (NList(j,i),j=1,Nneigh(i))
         !read(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
         !read(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
      end do
      close(1)
      close(2)
      close(3)
      close(4)

      call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
      NList2 = NList

!`d + idx_1 * direction_1 + idx_2 * direction_2` is smaller than `d`, where
!`idx_1` and `idx_2` are both in `[-1,0,1] and `direction_1` and `direction_2`
!are the real space directions of the translational symmetry.`

      call MIO_InputParameter('SuperCell',sCell,1)
      call MIO_InputParameter('readNeighborDetails',readNeighborDetails,.false.)
      if (readNeighborDetails) then ! Fengping spectral functions
         open(21,FILE='neighCell.dat')
         open(22,FILE='neighD.dat')
         do i=1,nAt
            do j=1,Nneigh(i)
                read(21,*) (neighCell(jj,j,i),jj=1,3)
                read(22,*) (neighD(jj,j,i),jj=1,3)
            end do
         end do
         close(21)
         close(22)
      else
         do i=1,nAt
            !r0 = Nradii(tbnn,(Species(i)-1)/2 + 1) ! The cut radius. Taken as the maximum distance that the
            do j=1,Nneigh(i) ! different from ultraSmall because here we already know the neighbors
               do ix=-1,1; do iy=-1,1
                  !print*, "ixiy", ix, iy
                  if (i==NList(j,i) .and. ix==0 .and. iy==0) cycle
                  ncell(:,1) = [ix,iy,0]
                  !print*, "matmul: " , matmul(ucell,ncell(:,1))
                  v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1))
                  d = norm(v)
                  !print*, "d: ", d
                  if (d < aG/sqrt(3.0_dp)*1.1_dp) then
                     neighCell(:,j,i) = ncell(:,1)
                     NeighD(:,j,i) = -v
                  end if
               end do; end do
            end do
         end do
      end if
   else

       call MIO_InputParameter('Neigh.LayerDistFactor',distFact,1.0_dp) ! to increase the cutoff for outerlayer neighbors


       dx=0.0_dp;dy=0.0_dp;dz=0.0_dp !;dr=0.0_dp

       PRINT*,''
       PRINT*,'Building nearest neighbor data...'
       num0 = 0
       num1 = 0
       num2 = 0
       num3 = 0
       num4 = 0
       num5 = 0
       num6 = 0
       num7 = 0
       num8 = 0
       num9 = 0
       num10= 0
       num11= 0
       num12= 0
       num13= 0
       num14= 0

       ! Divide geometry into rectangular cells of size xcell times ycell (units of Angstroms)
       ! A_celdaunitaria=3.0*dsqrt(3.0)*aCC*aCC/2.0=2.6*aCC*aCC=5.24Ang => rho=2/A_celdaunitaria
       !if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer')) then
       !   rho   = 8.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)
       !else
       rho   = 4.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)
       !end if
       xcell = 3.0_dp
       ycell = 3.0_dp
       !PRINT*,'Size of binning cell', xcell, 'X', ycell
       !============================================================
       !call MIO_Allocate(Species,nAt,'Species','atoms')
       !call MIO_Allocate(NeighD,[1,1,inode1],[3,maxNeigh,inode2],'NeighD','neigh')
       !call MIO_Allocate(xnew,[1],[natoms*9],'xnew','neigh')
       !call MIO_Allocate(ynew,[1],[natoms*9],'ynew','neigh')
       !call MIO_Allocate(znew,[1],[natoms*9],'znew','neigh')
       !allocate(xnew(9*natoms),ynew(9*natoms),znew(9*natoms))
       !print*, "numberOfAtoms: ", natoms, nAt
       !$OMP PARALLEL DO
       do i=1,natoms
         xnew(i)=x(i); ynew(i)=y(i); znew(i)=z(i)
       end do
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(n)
       do i=1,natoms
         n=i+natoms
         xnew(n)=x(i)+A1(1)
         ynew(n)=y(i)+A1(2)
         znew(n)=z(i)
       end do
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(n)
       do i=1,natoms
         n=i+2*natoms
         xnew(n)=x(i)+A1(1)+A2(1)
         ynew(n)=y(i)+A1(2)+A2(2)
         znew(n)=z(i)
       end do
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(n)
       do i=1,natoms
         n=i+3*natoms
         xnew(n)=x(i)+A2(1)
         ynew(n)=y(i)+A2(2)
         znew(n)=z(i)
       end do
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(n)
       do i=1,natoms
         n=i+4*natoms
         xnew(n)=x(i)-A1(1)+A2(1)
         ynew(n)=y(i)-A1(2)+A2(2)
         znew(n)=z(i)
       end do
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(n)
       do i=1,natoms
         n=i+5*natoms
         xnew(n)=x(i)-A1(1)
         ynew(n)=y(i)-A1(2)
         znew(n)=z(i)
       end do
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(n)
       do i=1,natoms
         n=i+6*natoms
         xnew(n)=x(i)-A1(1)-A2(1)
         ynew(n)=y(i)-A1(2)-A2(2)
         znew(n)=z(i)
       end do
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(n)
       do i=1,natoms
         n=i+7*natoms
         xnew(n)=x(i)-A2(1)
         ynew(n)=y(i)-A2(2)
         znew(n)=z(i)
       end do
       !$OMP END PARALLEL DO
       !$OMP PARALLEL DO PRIVATE(n)
       do i=1,natoms
         n=i+8*natoms
         xnew(n)=x(i)+A1(1)-A2(1)
         ynew(n)=y(i)+A1(2)-A2(2)
         znew(n)=z(i)
       end do
       !$OMP END PARALLEL DO
       !============================================================
       xmin = minval(xnew); xmax = maxval(xnew)
       ymin = minval(ynew); ymax = maxval(ynew)

       Nx = CEILING((xmax-xmin)/xcell)
       Ny = CEILING((ymax-ymin)/ycell)

       !write(*,*) xmin,xmax,Nx
       !write(*,*) ymin,ymax,Ny
       ALLOCATE(cells(Nx, Ny, 30*CEILING(xcell*ycell*rho))) ! 2 instead of 4 if not bilayer
       !write(*,*) 'Inicia llenado de cajas'
       !write(*,*) 'First dimension of cells is equal to ', Nx
       !write(*,*) 'Second dimension of cells is equal to ', Ny
       !write(*,*) 'Third dimension of cells is equal to ', 2*CEILING(xcell*ycell*rho)
        cells = 0
       !print*, cells
       !!$OMP PARALLEL DO PRIVATE(ix, iy, i) (make sure about the i incrementation first)
       do n = 1,9*natoms

         ix     = INT(1+(xnew(n)-xmin)/xcell)
         iy     = INT(1+(ynew(n)-ymin)/ycell)
         ixs(n) = ix
         iys(n) = iy

         i = 1
         !print*, ix, iy, i
         do while (cells(ix,iy,i) .ne. 0)
           !print*, i, cells(ix,iy,i)
           i = i+1
         end do
         cells(ix,iy,i) = n
         ! HERE
         !ncell(:,1) = [ix,iy,0]
         !neighCell(:,Nneigh(i),i) = ncell(:,1)
       end do
       !!$OMP END PARALLEL DO


       call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')

       do i=1,tbnn
          Nradii(i,1) = rad(i)*aG*1.1_dp
       end do
       if (MIO_StringComp(str,'BoronNitride')) then
          Nradii(:,2) = Nradii(:,1)
       else
          do i=1,tbnn
             Nradii(i,2) = rad(i)*aBN*1.1_dp
          end do
       end if
       rmax = maxval(Nradii)

       maxNeigh = maxnn
       print*, "inode1, inode2, in1, in2, maxnn", inode1, inode2, 1, natoms, in1, in2, maxnn
       call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
       !call MIO_Allocate(NList,[1,inode1],[maxnn,inode2],'NList','neigh')
       call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
       call MIO_Allocate(NeighD,[1,1,inode1],[3,maxnn,inode2],'NeighD','neigh')
       !call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
       call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,maxnn,inode2/),'neighCell','neigh')
       !ALLOCATE(cells(Nx, Ny, 2*CEILING(xcell*ycell*rho)))
       !call MIO_Allocate(2*CEILING(xcell*ycell*rho),
       ! HERE

       write(*,*) 'Inicia busqueda'
       ! Find NNs of each atom by only searching nearby cells
       safetycounter = 0
       !$OMP PARALLEL DO REDUCTION(+:num0, num1, num2, num3, num4, num5, num6, num7, num8, num9, num10, num11, num12, num13, num14) PRIVATE(ix, iy, NNcount, addx, addy, i, dx_t, dy_t, dz_t, m, d2_t, d2_t2, vecino)
       DO n = 1,natoms
         NNcount = 0

         DO addx = -2,2
         DO addy = -2,2
           ix = ixs(n) + addx
           iy = iys(n) + addy

           ! Periodic boundary conditions here
           IF(ix .gt. Nx) ix = ix-Nx
           IF(ix .lt. 1)  ix = ix+Nx
           IF(iy .gt. Ny) iy = iy-Ny
           IF(iy .lt. 1)  iy = iy+Ny

           i = 1
           DO WHILE(cells(ix, iy, i) .ne. 0)
             m = cells(ix, iy, i)
             IF (m .ne. n) THEN
               !write(*,*) 'hi2'
               dx_t = xnew(n) - xnew(m)
               dy_t = ynew(n) - ynew(m)
               dz_t = znew(n) - znew(m)
               ! Periodic boundary conditions here as well
               ! IF(dx_t .gt.  A1(1)-1.2_dp*aCC) dx_t = dx_t - A1(1)
               ! IF(dx_t .lt. -A1(1)+1.2_dp*aCC) dx_t = dx_t + A1(1)
               ! IF(dy_t .gt.  A2(2)-1.2_dp*aCC) dy_t = dy_t - A2(2)
               ! IF(dy_t .lt. -A2(2)+1.2_dp*aCC) dy_t = dy_t + A2(2)
               d2_t = dx_t**2.0_dp + dy_t**2.0_dp + dz_t**2.0_dp
               d2_t2 = dx_t**2.0_dp + dy_t**2.0_dp
                 !IF (n==1) print*, "distances from n=1 ", d2_t
                 !IF (n==3) print*, "distances from n=1 ", d2_t
               IF (abs(dz_t)<1.50_dp) THEN
                 !if (n .eq. 1) print*, d2_t2, cutoff2
                 IF (d2_t2 .lt. cutoff2*1.1_dp) THEN
                   NNcount = NNcount + 1
                   call NeighCheckCount(NNcount,maxnn)
                   vecino=mod(m,natoms)
                   if (vecino==0) vecino=natoms
                   nn(n,NNcount) = vecino
                   !NList(Nneigh(i),i) = j
                   dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
                   !dr(n,NNcount)=dsqrt(d2_t)
                   NList(NNcount,n) = vecino
                   NeighD(1,NNcount,n) = -dx(n,NNcount)
                   NeighD(2,NNcount,n) = -dy(n,NNcount)
                   NeighD(3,NNcount,n) = -dz(n,NNcount)
                   Nneigh(n) = NNcount
                   !neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
                   !neighCell(2,NNcount,n) = iys(m) - iys(n)
                   !neighCell(3,NNcount,n) = 0
       !             if(n==4800) write(*,*) m,NNcount,vecino
                 ENDIf
               ELSE IF (abs(dz_t)>1.50_dp .and. abs(dz_t)<4.50_dp) THEN
                 IF (d2_t2 .lt. cutoff2bis) THEN
                   NNcount = NNcount + 1
                   call NeighCheckCount(NNcount,maxnn)
                   vecino=mod(m,natoms)
                   if (vecino==0) vecino=natoms
                   nn(n,NNcount) = vecino
                   !NList(Nneigh(i),i) = j
                   dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
                   !dr(n,NNcount)=dsqrt(d2_t)
                   NList(NNcount,n) = vecino
                   NeighD(1,NNcount,n) = -dx(n,NNcount)
                   NeighD(2,NNcount,n) = -dy(n,NNcount)
                   NeighD(3,NNcount,n) = -dz(n,NNcount)
                   Nneigh(n) = NNcount
                   !neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
                   !neighCell(2,NNcount,n) = iys(m) - iys(n)
                   !neighCell(3,NNcount,n) = 0
       !             if(n==4800) write(*,*) m,NNcount,vecino
                 ENDIf
               ENDIF
             ENDIF
             i = i+1
           ENDDO

         ENDDO
         ENDDO

         !IF(NNcount .lt. 3) THEN
         !  print*, "Are you sure you wanted to enter this safety routine?"
         !  NNcount = 0
         !  DO m = 1,natoms
         !    IF (n.ne.m) THEN
         !      dx_t = x(n) - x(m)
         !      dy_t = y(n) - y(m)
         !      dz_t = z(n) - z(m)
         !       !Periodic boundary conditions here as well
         !       IF(dx_t .gt.  A1(1)-1.2_dp*aCC) dx_t = dx_t - A1(1)
         !       IF(dx_t .lt. -A1(1)+1.2_dp*aCC) dx_t = dx_t + A1(1)
         !       IF(dy_t .gt.  A2(2)-1.2_dp*aCC) dy_t = dy_t - A2(2)
         !       IF(dy_t .lt. -A2(2)+1.2_dp*aCC) dy_t = dy_t + A2(2)
         !       d2_t = dx_t**2.0_dp + dy_t**2.0_dp + dz_t**2.0_dp
         !       IF (abs(dz_t)<0.01_dp) THEN
         !          IF (d2_t .lt. cutoff2*1.1_dp) THEN
         !               NNcount = NNcount + 1
         !               !vecino=mod(m,natoms)
         !               !if (vecino==0) vecino=natoms
         !               !nn(n,NNcount) = vecino
         !               !NList(Nneigh(i),i) = j
         !               dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
         !               !dr(n,NNcount)=dsqrt(d2_t)
         !               NList(NNcount,n) = m
         !               NeighD(1,NNcount,n) = dx(n,NNcount)
         !               NeighD(2,NNcount,n) = dy(n,NNcount)
         !               NeighD(3,NNcount,n) = dz(n,NNcount)
         !               Nneigh(n) = NNcount
         !

         !               !ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
         !               !neighCell(:,Nneigh(i),i) = ncell(:,1)
         !

         !               !neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
         !               !neighCell(2,NNcount,n) = iys(m) - iys(n)
         !               !neighCell(3,NNcount,n) = 0
       ! !               if(n==4800) write(*,*) m,NNcount,vecino
         !          ENDIf
         !       END IF
         !     END IF
         !  END DO
         !  safetycounter = safetycounter + 1
         !  if (safetycounter.gt.100) print*, "Carefull, you are probably going too many times through this safety loop. &
         ! I originally implemented this because the fastNNnotsquared routine was missing neighbors for 4 atoms only..."
         !ENDIF


       ! Count the number of atoms with 1, 2, 3, 4, 5 nearest neighbors
         IF(NNcount .eq. 0) num0 = num0 + 1
         !IF(NNcount .eq. 0) print*, "NN counter = 0", n
         IF(NNcount .eq. 1) num1 = num1 + 1
         IF(NNcount .eq. 2) num2 = num2 + 1
         !IF(NNcount .eq. 2) print*, "NN counter = 2", n
         IF(NNcount .eq. 3) num3 = num3 + 1
         IF(NNcount .eq. 4) num4 = num4 + 1
         IF(NNcount .eq. 5) num5 = num5 + 1
         IF(NNcount .eq. 6) num6 = num6 + 1
         IF(NNcount .eq. 7) num7 = num7 + 1
         IF(NNcount .eq. 8) num8 = num8 + 1
         IF(NNcount .eq. 9) num9 = num9 + 1
         IF(NNcount .eq. 10) num10 = num10 + 1
         IF(NNcount .eq. 11) num11 = num11 + 1
         IF(NNcount .eq. 12) num12 = num12 + 1
         IF(NNcount .eq. 13) num13 = num13 + 1
         IF(NNcount .eq. 14) num14 = num14 + 1



       ENDDO
       !$OMP END PARALLEL DO

       call MIO_InputParameter('SuperCell',sCell,1)
       !if (sCell==1) then
          !!$OMP PARALLEL PRIVATE (v, d,ncell,i,j,ix,iy,dist)
          !do i=inode1,inode2
          !$OMP PARALLEL DO PRIVATE(v, d, ncell, i,j, dist, ix, iy)
          do i=1,nAt
             !r0 = Nradii(tbnn,(Species(i)-1)/2 + 1) ! The cut radius. Taken as the maximum distance that the
             do j=1,Nneigh(i) ! different from ultraSmall because here we already know the neighbors
                do ix=-5,5; do iy=-5,5
                   if (i==NList(j,i) .and. ix==0 .and. iy==0) cycle
                   ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
                   v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1))
                   d = norm(v)
                   dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                   if (abs(d-dist).lt.0.1_dp) then
                      !print*, "muak"
                      neighCell(:,j,i) = ncell(:,1)
                   end if
                   !v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition.
                                                                      ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
                   !d = norm(v)
                   !if (d**2.0_dp .lt. cutoff2*1.1_dp) then
                   !   !Nneigh(i) = Nneigh(i) + 1
                   !   !if (Nneigh(i)>maxNeigh)  then
                   !   !   !print*, "Number of neighbors: ", i, Nneigh(i)
                   !   !   call MIO_Kill('More neighbors than expected found',&
                   !   !     'neigh','NeighList')
                   !   !end if
                   !   !NList(Nneigh(i),i) = j
                   !   !NeighD(:,Nneigh(i),i) = -v
                   !   neighCell(:,j,i) = ncell(:,1)
                   !!else
                   !!    print*, "Are you sure this was a neighbor?", i, NList(j,i)
                   !end if
                end do; end do
             end do
          end do
          !!$OMP END PARALLEL
       !end if

       !write(*,*) 'hi4'
       !neighCell = cells
       DEALLOCATE(cells)
       !write(*,*) 'hi6'
       !call MIO_Deallocate(xnew,'xnew','neigh')
       !call MIO_Deallocate(ynew,'ynew','neigh')
       !call MIO_Deallocate(znew,'znew','neigh')
       !deallocate(xnew,ynew,znew)
       !write(*,*) 'hi5'


       PRINT*,' 0 neighbors: ',num0
       PRINT*,' 1 neighbors: ',num1
       PRINT*,' 2 neighbors: ',num2
       PRINT*,' 3 neighbors: ',num3
       PRINT*,' 4 neighbors: ',num4
       PRINT*,' 5 neighbors: ',num5
       PRINT*,' 6 neighbors: ',num6
       PRINT*,' 7 neighbors: ',num7
       PRINT*,' 8 neighbors: ',num8
       PRINT*,' 9 neighbors: ',num9
       PRINT*,' 10 neighbors: ',num10
       PRINT*,' 11 neighbors: ',num11
       PRINT*,' 12 neighbors: ',num12
       PRINT*,' 13 neighbors: ',num13
       PRINT*,' 14 neighbors: ',num14
       PRINT*,'>10 neighbors: ',natoms-(num0+num1+num2+num3+num4+num5+num6+num7+num8+num9+num10+num11+num12+num13+num14)

       PRINT*,''
       PRINT*,'...done'
       PRINT*,''

       call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
       NList2 = NList
#ifdef MPI
       !print*, "nProc =", nProc
       !if (nProc>1) then
       !   call MIO_Allocate(countN,[nProc,numThreads],'countN','neigh')
       !   countN = 0
       !   !do i=1,natoms
       !   !$OMP PARALLEL
       !   do i=in1,in2
       !      do j=1,maxNeigh
       !         if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
       !            do k=nProc,1,-1
       !               if(NList2(j,i)>=indxNode(k)) then
       !                  countN(k,nThread+1) = countN(k,nThread+1)+1
       !                  exit
       !               end if
       !            end do
       !         end if
       !      end do
       !   end do
       !   !$OMP END PARALLEL
       !   rcvSz = sum(countN)
       !   call MIO_Allocate(rcvIndx,[2,nProc],'rcvIndx','neigh')
       !   call MIO_Allocate(rcvList,rcvSz,'rcvList','neigh')
       !   np = 1
       !   do i=1,nProc
       !      rcvIndx(1,i) = np
       !      np = np + sum(countN(i,:))
       !      rcvIndx(2,i) = np - 1
       !   end do
       !   call MIO_Allocate(cnt,[nProc,numThreads],'cnt','neigh')
       !   !do i=1,natoms
       !   !$OMP PARALLEL PRIVATE(np)
       !   do i=in1,in2
       !      do j=1,maxNeigh
       !         if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
       !            do k=nProc,1,-1
       !               if(NList2(j,i)>=indxNode(k)) then
       !                  np = rcvIndx(1,k) + sum(countN(k,1:nThread)) + cnt(k,nThread+1)
       !                  cnt(k,nThread+1) = cnt(k,nThread+1) + 1
       !                  exit
       !               end if
       !            end do
       !            rcvList(np) = NList2(j,i)
       !            NList2(j,i) = np + inode2
       !         end if
       !      end do
       !   end do
       !   !$OMP END PARALLEL
       !   call MIO_Allocate(nodeSndRcv,[nProc,nProc],'nodeSndRcv','neigh')
       !   do i=1,nProc
       !      nodeSndRcv(i,Node+1) = rcvIndx(2,i)-rcvIndx(1,i) + 1
       !   end do
       !   call MPIAllGather(MPI_IN_PLACE,0,nodeSndRcv,nProc,MPI_INTEGER)
       !   sndSz = sum(nodeSndRcv(Node+1,:))
       !   call MIO_Allocate(sndList,sndSz,'sndList','neigh')
       !   call MIO_Allocate(sndIndx,(/2,nProc/),'sndIndx','neigh')
       !   np = 1
       !   do i=1,nProc
       !      sndIndx(1,i) = np
       !      np = np + nodeSndRcv(Node+1,i)
       !      sndIndx(2,i) = np - 1
       !   end do
       !   do i=1,nProc-1
       !      j = mod(Node + i,nProc)
       !      k = mod(nProc + Node - i,nProc)
       !      call MPISendRecv(rcvList(rcvIndx(1,j+1):rcvIndx(2,j+1)),nodeSndRcv(j+1,Node+1),j,50, &
       !        sndList(sndIndx(1,k+1):sndIndx(2,k+1)),nodeSndRcv(Node+1,k+1),k,50,MPI_INTEGER)
       !      call MPISendRecv(rcvList(rcvIndx(1,k+1):rcvIndx(2,k+1)),nodeSndRcv(k+1,Node+1),k,51, &
       !        sndList(sndIndx(1,j+1):sndIndx(2,j+1)),nodeSndRcv(Node+1,j+1),j,51,MPI_INTEGER)
       !      call MPIBarrier()
       !   end do
       !   call MIO_Deallocate(cnt,'cnt','neigh')
       !   call MIO_Deallocate(countN,'countN','neigh')
       !end if
#endif /* MPI */
       call MIO_InputParameter('WriteDataFiles',prnt,.false.)
       !call AtomsSetFrac()
       if (prnt) then
          open(1,FILE='v')
          open(2,FILE='dx')
          open(3,FILE='dy')
          open(33,FILE='dz')
          open(4,FILE='pos')
          do i=1,nAt
             write(1,*) Nneigh(i)
             !write(1,'(300(I8,1X))') (NList(j,i),j=1,Nneigh(i))
             !write(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
             !write(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
             !write(33,*) (NeighD(3,j,i),j=1,Nneigh(i))
             write(1,'(500(I8,1X))') (NList(j,i),j=1,Nneigh(i))
             write(2,'(500(F10.5,1X))') (NeighD(1,j,i),j=1,Nneigh(i))
             write(3,'(500(F10.5,1X))') (NeighD(2,j,i),j=1,Nneigh(i))
             write(33,'(500(F10.5,1X))') (NeighD(3,j,i),j=1,Nneigh(i))
             if (Species(i)==3) then
                write(4,*) "B    ", Rat(1,i),Rat(2,i),Rat(3,i)
             else if (Species(i)==4) then
                write(4,*) "N    ", Rat(1,i),Rat(2,i),Rat(3,i)
             else
                write(4,*) "C    ", Rat(1,i),Rat(2,i),Rat(3,i)
             end if
          end do
          close(1)
          close(2)
          close(3)
          close(33)
          close(4)
       end if
   end if


end subroutine fastNNnotsquare

subroutine fastNNnotsquareSmall(natoms,x,y,z,aCC,cutoff2,cutoff2bis,A1,A2,maxnn)

   use atoms,                only : inode1, inode2, in1, in2, indxNode, nAt, Species, Rat
   use tbpar,                only : tbnn
   use cell,                 only : aG, aBN,ucell
   use atoms,                only : AtomsSetCart, frac
   use math
   implicit none
   integer, intent(in) :: natoms,maxnn
   real(dp), intent(in) :: x(natoms),y(natoms),z(natoms)
   integer :: nn(natoms,maxnn)
   real(dp), intent(in) :: aCC,cutoff2,cutoff2bis,A1(3),A2(3)
   integer :: num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,num10,num11,num12,num13,num14
   integer :: i,j,k,l,m,n
   real(dp) :: dx(natoms,maxnn),dy(natoms,maxnn),dz(natoms,maxnn)
   !real(dp), intent(out) :: dr(natoms,maxnn)

   real(dp) :: dx_t,dy_t,dz_t,d2_t, d2_t2
   real(dp) :: xmin,xmax,ymin,ymax,xcell,ycell,rho

   !real(dp) :: xnew(9*natoms),ynew(9*natoms),znew(9*natoms)
   !integer :: ixs(9*natoms),iys(9*natoms),vecino
   real(dp) :: xnew(11*11*natoms*3),ynew(11*11*natoms*3),znew(11*11*natoms*3)
   integer :: ixs(11*11*natoms*3),iys(11*11*natoms*3),vecino
   !integer :: izs(11*11*natoms*3)
   integer :: nInteger(11*11*natoms*3)
   integer :: mInteger(11*11*natoms*3)

   integer :: Nx,Ny,ix,iy
   integer, allocatable :: cells(:,:,:)
   integer :: NNcount,addx,addy

   !real(dp), parameter :: rad(3) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp)]
   real(dp), parameter :: rad(5) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp),3.0_dp/sqrt(3.0_dp),3.0_dp/sqrt(3.0_dp)]
   character(len=50) :: str

   integer :: ncell(3,9)
   real(dp) :: d, v(3), vtemp(3)
   real(dp) :: rmax, r0, rmaxTemp

   integer, pointer :: countN(:,:), cnt(:,:)
   integer :: np

   integer :: safetycounter
   logical :: prnt
   logical :: BLInPlaneInteractionRadius

   real(dp) :: distFact, dist

   integer :: sCell

   logical :: readDataFiles, only000Cell

   integer :: uX, uY, uZ, nmax, unitsX, unitsY


   call MIO_InputParameter('TypeOfSystem',str,'Graphene')
   call MIO_InputParameter('ReadDataFiles',readDataFiles,.false.)
   !if (frac) call AtomsSetCart()
   !call AtomsSetFrac()
   if (readDataFiles) then
      call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')

      do i=1,tbnn
         Nradii(i,1) = rad(i)*aG*1.1_dp
      end do
      if (MIO_StringComp(str,'BoronNitride')) then
         Nradii(:,2) = Nradii(:,1)
      else
         do i=1,tbnn
            Nradii(i,2) = rad(i)*aBN*1.1_dp
         end do
      end if
      rmax = maxval(Nradii)

      maxNeigh = maxnn

      !print*, "allocating NList with maxNeigh= ", maxNeigh
      call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
      !call MIO_Allocate(NList,[1,inode1],[maxnn,inode2],'NList','neigh')
      call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
      call MIO_Allocate(NeighD,[1,1,inode1],[3,maxnn,inode2],'NeighD','neigh')
      !call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
      call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,maxnn,inode2/),'neighCell','neigh')
      open(1,FILE='v')
      open(2,FILE='dx')
      open(3,FILE='dy')
      !open(4,FILE='pos')
      !print*, "nAt =", nAt
      do i=1,nAt
         read(1,*) Nneigh(i)
         read(1,*) (NList(j,i),j=1,Nneigh(i))
         !read(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
         !read(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
      end do
      close(1)
      close(2)
      close(3)
      close(4)

      call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
      NList2 = NList

!`d + idx_1 * direction_1 + idx_2 * direction_2` is smaller than `d`, where
!`idx_1` and `idx_2` are both in `[-1,0,1] and `direction_1` and `direction_2`
!are the real space directions of the translational symmetry.`

      call MIO_InputParameter('SuperCell',sCell,1)
      call MIO_InputParameter('only000Cell',only000Cell,.false.)
      !print*, "ucell", ucell
      !print*, Rat
      !print*, "jjhhh"
      !print*, Rat(:,1)
      !print*, Rat(:,NList(1,1))
      !print*, "mmmm"
      !if (sCell==1) then
         !!$OMP PARALLEL PRIVATE (v, d,ncell,i,j,ix,iy,dist)
         !do i=inode1,inode2
         do i=1,nAt
            !r0 = Nradii(tbnn,(Species(i)-1)/2 + 1) ! The cut radius. Taken as the maximum distance that the
            do j=1,Nneigh(i) ! different from ultraSmall because here we already know the neighbors
               if (only000Cell) then
                  !do ix=-1,1; do iy=-1,1
                  ix = 0
                  iy = 0
                  !print*, "ixiy", ix, iy
                  if (i==NList(j,i) .and. ix==0 .and. iy==0) cycle
                  ncell(:,1) = [ix,iy,0]
                  !print*, "matmul: " , matmul(ucell,ncell(:,1))
                  v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1))
                  d = norm(v)
                  !print*, "d: ", d
                  if (d < aG/sqrt(3.0_dp)*1.1_dp) then
                     neighCell(:,j,i) = ncell(:,1)
                     NeighD(:,j,i) = -v
                  end if
                  !dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                  !dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                  !if (abs(d-dist).lt.0.1_dp) then
                  !   !print*, "muak"
                  !   neighCell(:,j,i) = ncell(:,1)
                  !end if

                  !v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition.
                                                                    ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
                  !d = norm(v)
                  !if (d**2.0_dp .lt. cutoff2*1.1_dp) then
                     !Nneigh(i) = Nneigh(i) + 1
                     !if (Nneigh(i)>maxNeigh)  then
                     !   !print*, "Number of neighbors: ", i, Nneigh(i)
                     !   call MIO_Kill('More neighbors than expected found',&
                     !     'neigh','NeighList')
                     !end if
                     !NList(Nneigh(i),i) = j
                     !NeighD(:,Nneigh(i),i) = -v
                     !neighCell(:,j,i) = ncell(:,1)
                  !else
                  !    print*, "Are you sure this was a neighbor?", i, NList(j,i)
                  !end if
               else
                  do ix=-1,1; do iy=-1,1
                     !print*, "ixiy", ix, iy
                     if (i==NList(j,i) .and. ix==0 .and. iy==0) cycle
                     ncell(:,1) = [ix,iy,0]
                     !print*, "matmul: " , matmul(ucell,ncell(:,1))
                     v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1))
                     d = norm(v)
                     !print*, "d: ", d
                     if (d < aG/sqrt(3.0_dp)*1.1_dp) then
                        neighCell(:,j,i) = ncell(:,1)
                        NeighD(:,j,i) = -v
                     end if
                     !dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                     !dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                     !if (abs(d-dist).lt.0.1_dp) then
                     !   !print*, "muak"
                     !   neighCell(:,j,i) = ncell(:,1)
                     !end if

                     !v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition.
                                                                       ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
                     !d = norm(v)
                     !if (d**2.0_dp .lt. cutoff2*1.1_dp) then
                        !Nneigh(i) = Nneigh(i) + 1
                        !if (Nneigh(i)>maxNeigh)  then
                        !   !print*, "Number of neighbors: ", i, Nneigh(i)
                        !   call MIO_Kill('More neighbors than expected found',&
                        !     'neigh','NeighList')
                        !end if
                        !NList(Nneigh(i),i) = j
                        !NeighD(:,Nneigh(i),i) = -v
                        !neighCell(:,j,i) = ncell(:,1)
                     !else
                     !    print*, "Are you sure this was a neighbor?", i, NList(j,i)
                     !end if
                  end do; end do
               end if
            end do
         end do
         !!$OMP END PARALLEL
      !end if
   else

       call MIO_InputParameter('Neigh.LayerDistFactor',distFact,1.0_dp) ! to increase the cutoff for outerlayer neighbors


       dx=0.0_dp;dy=0.0_dp;dz=0.0_dp !;dr=0.0_dp

       PRINT*,''
       PRINT*,'Building nearest neighbor data...'
       num0 = 0
       num1 = 0
       num2 = 0
       num3 = 0
       num4 = 0
       num5 = 0
       num6 = 0
       num7 = 0
       num8 = 0
       num9 = 0
       num10= 0
       num11= 0
       num12= 0
       num13= 0
       num14= 0

       ! Divide geometry into rectangular cells of size xcell times ycell (units of Angstroms)
       ! A_celdaunitaria=3.0*dsqrt(3.0)*aCC*aCC/2.0=2.6*aCC*aCC=5.24Ang => rho=2/A_celdaunitaria
       !if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer')) then
       !   rho   = 8.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)
       !else
       !end if
       xcell = 3.0_dp
       ycell = 3.0_dp
       !PRINT*,'Size of binning cell', xcell, 'X', ycell
       !============================================================
       !call MIO_Allocate(Species,nAt,'Species','atoms')
       !call MIO_Allocate(NeighD,[1,1,inode1],[3,maxNeigh,inode2],'NeighD','neigh')
       !call MIO_Allocate(xnew,[1],[natoms*9],'xnew','neigh')
       !call MIO_Allocate(ynew,[1],[natoms*9],'ynew','neigh')
       !call MIO_Allocate(znew,[1],[natoms*9],'znew','neigh')
       !allocate(xnew(9*natoms),ynew(9*natoms),znew(9*natoms))
       !print*, "numberOfAtoms: ", natoms, nAt
       unitsX = 10
       unitsY = 10
       n = 0
       do i=1,natoms
         n=n+1
         uX = 0
         uY = 0
         !uZ = 0
         xnew(n)=x(i)+uX*A1(1)+uY*A2(1)
         ynew(n)=y(i)+uX*A1(2)+uY*A2(2)
         !znew(n)=z(i)+uZ*A3
         nInteger(n) = uX
         mInteger(n) = uY
       end do
       do uX=-unitsX/2,unitsX/2
          do uY=-unitsY/2,unitsY/2
             do uZ=-1,1
                do i=1,natoms
                  if (uX==0 .and. uY==0 .and. uZ==0) cycle
                  n=n+1
                  xnew(n)=x(i)+uX*A1(1)+uY*A2(1)
                  ynew(n)=y(i)+uX*A1(2)+uY*A2(2)
                  !znew(n)=z(i)+uZ*A3
                  nInteger(n) = uX
                  mInteger(n) = uY
                end do
             end do
          end do
       end do
       nmax = n
       !print*, "nmax", nmax
       !============================================================
       xmin = minval(xnew); xmax = maxval(xnew)
       ymin = minval(ynew); ymax = maxval(ynew)

       Nx = CEILING((xmax-xmin)/xcell)
       Ny = CEILING((ymax-ymin)/ycell)

       !write(*,*) xmin,xmax,Nx
       !write(*,*) ymin,ymax,Ny
       rho   = 4.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)
       ALLOCATE(cells(Nx, Ny, 4*CEILING(xcell*ycell*rho*unitsX*unitsY))) ! 2 instead of 4 if not bilayer
       !write(*,*) 'Inicia llenado de cajas'
       !write(*,*) 'First dimension of cells is equal to ', Nx
       !write(*,*) 'Second dimension of cells is equal to ', Ny
       !write(*,*) 'Third dimension of cells is equal to ', 2*CEILING(xcell*ycell*rho)
        cells = 0
       !print*, cells
       !!$OMP PARALLEL DO PRIVATE(ix, iy, i) (make sure about the i incrementation first)
       do n = 1,nmax

         ix     = INT(1+(xnew(n)-xmin)/xcell)
         iy     = INT(1+(ynew(n)-ymin)/ycell)
         ixs(n) = ix
         iys(n) = iy

         i = 1
         !print*, ix, iy, i
         do while (cells(ix,iy,i) .ne. 0)
           !print*, i, cells(ix,iy,i)
           i = i+1
         end do
         cells(ix,iy,i) = n
         ! HERE
         !ncell(:,1) = [ix,iy,0]
         !neighCell(:,Nneigh(i),i) = ncell(:,1)
       end do
       !!$OMP END PARALLEL DO


       call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')

       do i=1,tbnn
          Nradii(i,1) = rad(i)*aG*1.1_dp
       end do
       if (MIO_StringComp(str,'BoronNitride')) then
          Nradii(:,2) = Nradii(:,1)
       else
          do i=1,tbnn
             Nradii(i,2) = rad(i)*aBN*1.1_dp
          end do
       end if
       rmax = maxval(Nradii)

       maxNeigh = maxnn
       !print*, "inode1, inode2, in1, in2, maxnn", inode1, inode2, 1, natoms, in1, in2, maxnn
       call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
       !call MIO_Allocate(NList,[1,inode1],[maxnn,inode2],'NList','neigh')
       call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
       call MIO_Allocate(NeighD,[1,1,inode1],[3,maxnn,inode2],'NeighD','neigh')
       !call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
       call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,maxnn,inode2/),'neighCell','neigh')
       !ALLOCATE(cells(Nx, Ny, 2*CEILING(xcell*ycell*rho)))
       !call MIO_Allocate(2*CEILING(xcell*ycell*rho),
       ! HERE

       write(*,*) 'Inicia busqueda'
       ! Find NNs of each atom by only searching nearby cells
       safetycounter = 0
       !!$OMP PARALLEL DO REDUCTION(+:num0, num1, num2, num3, num4, num5, num6, num7, num8, num9, num10, num11, num12, num13, num14) PRIVATE(ix, iy, NNcount, addx, addy, i, dx_t, dy_t, dz_t, m, d2_t, d2_t2, vecino)
       DO n = 1,natoms
         NNcount = 0

         DO addx = -2,2
         DO addy = -2,2
           ix = ixs(n) + addx
           iy = iys(n) + addy

           ! Periodic boundary conditions here
           IF(ix .gt. Nx) ix = ix-Nx
           IF(ix .lt. 1)  ix = ix+Nx
           IF(iy .gt. Ny) iy = iy-Ny
           IF(iy .lt. 1)  iy = iy+Ny

           i = 1
           DO WHILE(cells(ix, iy, i) .ne. 0)
             m = cells(ix, iy, i)
             IF (m .ne. n) THEN
               !write(*,*) 'hi2'
               dx_t = xnew(n) - xnew(m)
               dy_t = ynew(n) - ynew(m)
               dz_t = znew(n) - znew(m)
               ! Periodic boundary conditions here as well
               ! IF(dx_t .gt.  A1(1)-1.2_dp*aCC) dx_t = dx_t - A1(1)
               ! IF(dx_t .lt. -A1(1)+1.2_dp*aCC) dx_t = dx_t + A1(1)
               ! IF(dy_t .gt.  A2(2)-1.2_dp*aCC) dy_t = dy_t - A2(2)
               ! IF(dy_t .lt. -A2(2)+1.2_dp*aCC) dy_t = dy_t + A2(2)
               d2_t = dx_t**2.0_dp + dy_t**2.0_dp + dz_t**2.0_dp
               d2_t2 = dx_t**2.0_dp + dy_t**2.0_dp
                 !IF (n==1) print*, "distances from n=1 ", d2_t
                 !IF (n==3) print*, "distances from n=1 ", d2_t
               IF (abs(dz_t)<1.50_dp) THEN
                 !if (n .eq. 1) print*, d2_t2, cutoff2
                 IF (d2_t2 .lt. cutoff2*1.1_dp) THEN
                   NNcount = NNcount + 1
                   call NeighCheckCount(NNcount,maxnn)
                   vecino=mod(m,natoms)
                   if (vecino==0) vecino=natoms
                   nn(n,NNcount) = vecino
                   !NList(Nneigh(i),i) = j
                   dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
                   !dr(n,NNcount)=dsqrt(d2_t)
                   NList(NNcount,n) = vecino
                   NeighD(1,NNcount,n) = -dx(n,NNcount)
                   NeighD(2,NNcount,n) = -dy(n,NNcount)
                   NeighD(3,NNcount,n) = -dz(n,NNcount)
                   Nneigh(n) = NNcount
                   neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
                   neighCell(2,NNcount,n) = iys(m) - iys(n)
                   neighCell(3,NNcount,n) = 0
                   !neighCell(1,NNcount,n) = nInteger(m) - nInteger(n)
                   !neighCell(2,NNcount,n) = mInteger(m) - mInteger(n)
                   !neighCell(3,NNcount,n) = izs(m) - izs(n)
       !             if(n==4800) write(*,*) m,NNcount,vecino
                 ENDIf
               ELSE IF (abs(dz_t)>1.50_dp .and. abs(dz_t)<4.50_dp) THEN
                 IF (d2_t2 .lt. cutoff2bis) THEN
                   NNcount = NNcount + 1
                   call NeighCheckCount(NNcount,maxnn)
                   vecino=mod(m,natoms)
                   if (vecino==0) vecino=natoms
                   nn(n,NNcount) = vecino
                   !NList(Nneigh(i),i) = j
                   dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
                   !dr(n,NNcount)=dsqrt(d2_t)
                   NList(NNcount,n) = vecino
                   NeighD(1,NNcount,n) = -dx(n,NNcount)
                   NeighD(2,NNcount,n) = -dy(n,NNcount)
                   NeighD(3,NNcount,n) = -dz(n,NNcount)
                   Nneigh(n) = NNcount
                   neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
                   neighCell(2,NNcount,n) = iys(m) - iys(n)
                   neighCell(3,NNcount,n) = 0
                   !neighCell(1,NNcount,n) = nInteger(m) - nInteger(n)
                   !neighCell(2,NNcount,n) = mInteger(m) - mInteger(n)
                   !neighCell(3,NNcount,n) = izs(m) - izs(n)
       !             if(n==4800) write(*,*) m,NNcount,vecino
                 ENDIf
               ENDIF
             ENDIF
             i = i+1
           ENDDO

         ENDDO
         ENDDO

         !IF(NNcount .lt. 3) THEN
         !  print*, "Are you sure you wanted to enter this safety routine?"
         !  NNcount = 0
         !  DO m = 1,natoms
         !    IF (n.ne.m) THEN
         !      dx_t = x(n) - x(m)
         !      dy_t = y(n) - y(m)
         !      dz_t = z(n) - z(m)
         !       !Periodic boundary conditions here as well
         !       IF(dx_t .gt.  A1(1)-1.2_dp*aCC) dx_t = dx_t - A1(1)
         !       IF(dx_t .lt. -A1(1)+1.2_dp*aCC) dx_t = dx_t + A1(1)
         !       IF(dy_t .gt.  A2(2)-1.2_dp*aCC) dy_t = dy_t - A2(2)
         !       IF(dy_t .lt. -A2(2)+1.2_dp*aCC) dy_t = dy_t + A2(2)
         !       d2_t = dx_t**2.0_dp + dy_t**2.0_dp + dz_t**2.0_dp
         !       IF (abs(dz_t)<0.01_dp) THEN
         !          IF (d2_t .lt. cutoff2*1.1_dp) THEN
         !               NNcount = NNcount + 1
         !               !vecino=mod(m,natoms)
         !               !if (vecino==0) vecino=natoms
         !               !nn(n,NNcount) = vecino
         !               !NList(Nneigh(i),i) = j
         !               dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
         !               !dr(n,NNcount)=dsqrt(d2_t)
         !               NList(NNcount,n) = m
         !               NeighD(1,NNcount,n) = dx(n,NNcount)
         !               NeighD(2,NNcount,n) = dy(n,NNcount)
         !               NeighD(3,NNcount,n) = dz(n,NNcount)
         !               Nneigh(n) = NNcount
         !

         !               !ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
         !               !neighCell(:,Nneigh(i),i) = ncell(:,1)
         !

         !               !neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
         !               !neighCell(2,NNcount,n) = iys(m) - iys(n)
         !               !neighCell(3,NNcount,n) = 0
       ! !               if(n==4800) write(*,*) m,NNcount,vecino
         !          ENDIf
         !       END IF
         !     END IF
         !  END DO
         !  safetycounter = safetycounter + 1
         !  if (safetycounter.gt.100) print*, "Carefull, you are probably going too many times through this safety loop. &
         ! I originally implemented this because the fastNNnotsquared routine was missing neighbors for 4 atoms only..."
         !ENDIF


       ! Count the number of atoms with 1, 2, 3, 4, 5 nearest neighbors
         IF(NNcount .eq. 0) num0 = num0 + 1
         !IF(NNcount .eq. 0) print*, "NN counter = 0", n
         IF(NNcount .eq. 1) num1 = num1 + 1
         IF(NNcount .eq. 2) num2 = num2 + 1
         !IF(NNcount .eq. 2) print*, "NN counter = 2", n
         IF(NNcount .eq. 3) num3 = num3 + 1
         IF(NNcount .eq. 4) num4 = num4 + 1
         IF(NNcount .eq. 5) num5 = num5 + 1
         IF(NNcount .eq. 6) num6 = num6 + 1
         IF(NNcount .eq. 7) num7 = num7 + 1
         IF(NNcount .eq. 8) num8 = num8 + 1
         IF(NNcount .eq. 9) num9 = num9 + 1
         IF(NNcount .eq. 10) num10 = num10 + 1
         IF(NNcount .eq. 11) num11 = num11 + 1
         IF(NNcount .eq. 12) num12 = num12 + 1
         IF(NNcount .eq. 13) num13 = num13 + 1
         IF(NNcount .eq. 14) num14 = num14 + 1



       ENDDO
       !!$OMP END PARALLEL DO

       call MIO_InputParameter('SuperCell',sCell,1)
       !if (sCell==1) then
          !!$OMP PARALLEL PRIVATE (v, d,ncell,i,j,ix,iy,dist)
          !do i=inode1,inode2
          !$OMP PARALLEL DO PRIVATE(v, d, ncell, i,j, dist, ix, iy)
          do i=1,nAt
             !r0 = Nradii(tbnn,(Species(i)-1)/2 + 1) ! The cut radius. Taken as the maximum distance that the
             do j=1,Nneigh(i) ! different from ultraSmall because here we already know the neighbors
                do ix=-5,5; do iy=-5,5
                   if (i==NList(j,i) .and. ix==0 .and. iy==0) cycle
                   ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
                   v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1))
                   d = norm(v)
                   dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                   if (abs(d-dist).lt.0.1_dp) then
                      !print*, "muak"
                      neighCell(:,j,i) = ncell(:,1)
                   end if
                   !v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition.
                                                                      ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
                   !d = norm(v)
                   !if (d**2.0_dp .lt. cutoff2*1.1_dp) then
                   !   !Nneigh(i) = Nneigh(i) + 1
                   !   !if (Nneigh(i)>maxNeigh)  then
                   !   !   !print*, "Number of neighbors: ", i, Nneigh(i)
                   !   !   call MIO_Kill('More neighbors than expected found',&
                   !   !     'neigh','NeighList')
                   !   !end if
                   !   !NList(Nneigh(i),i) = j
                   !   !NeighD(:,Nneigh(i),i) = -v
                   !   neighCell(:,j,i) = ncell(:,1)
                   !!else
                   !!    print*, "Are you sure this was a neighbor?", i, NList(j,i)
                   !end if
                end do; end do
             end do
          end do
          !!$OMP END PARALLEL
       !end if

       !write(*,*) 'hi4'
       !neighCell = cells
       DEALLOCATE(cells)
       !write(*,*) 'hi6'
       !call MIO_Deallocate(xnew,'xnew','neigh')
       !call MIO_Deallocate(ynew,'ynew','neigh')
       !call MIO_Deallocate(znew,'znew','neigh')
       !deallocate(xnew,ynew,znew)
       !write(*,*) 'hi5'


       PRINT*,' 0 neighbors: ',num0
       PRINT*,' 1 neighbors: ',num1
       PRINT*,' 2 neighbors: ',num2
       PRINT*,' 3 neighbors: ',num3
       PRINT*,' 4 neighbors: ',num4
       PRINT*,' 5 neighbors: ',num5
       PRINT*,' 6 neighbors: ',num6
       PRINT*,' 7 neighbors: ',num7
       PRINT*,' 8 neighbors: ',num8
       PRINT*,' 9 neighbors: ',num9
       PRINT*,' 10 neighbors: ',num10
       PRINT*,' 11 neighbors: ',num11
       PRINT*,' 12 neighbors: ',num12
       PRINT*,' 13 neighbors: ',num13
       PRINT*,' 14 neighbors: ',num14
       PRINT*,'>10 neighbors: ',natoms-(num0+num1+num2+num3+num4+num5+num6+num7+num8+num9+num10+num11+num12+num13+num14)

       PRINT*,''
       PRINT*,'...done'
       PRINT*,''

       call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
       NList2 = NList
#ifdef MPI
       !print*, "nProc =", nProc
       !if (nProc>1) then
       !   call MIO_Allocate(countN,[nProc,numThreads],'countN','neigh')
       !   countN = 0
       !   !do i=1,natoms
       !   !$OMP PARALLEL
       !   do i=in1,in2
       !      do j=1,maxNeigh
       !         if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
       !            do k=nProc,1,-1
       !               if(NList2(j,i)>=indxNode(k)) then
       !                  countN(k,nThread+1) = countN(k,nThread+1)+1
       !                  exit
       !               end if
       !            end do
       !         end if
       !      end do
       !   end do
       !   !$OMP END PARALLEL
       !   rcvSz = sum(countN)
       !   call MIO_Allocate(rcvIndx,[2,nProc],'rcvIndx','neigh')
       !   call MIO_Allocate(rcvList,rcvSz,'rcvList','neigh')
       !   np = 1
       !   do i=1,nProc
       !      rcvIndx(1,i) = np
       !      np = np + sum(countN(i,:))
       !      rcvIndx(2,i) = np - 1
       !   end do
       !   call MIO_Allocate(cnt,[nProc,numThreads],'cnt','neigh')
       !   !do i=1,natoms
       !   !$OMP PARALLEL PRIVATE(np)
       !   do i=in1,in2
       !      do j=1,maxNeigh
       !         if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
       !            do k=nProc,1,-1
       !               if(NList2(j,i)>=indxNode(k)) then
       !                  np = rcvIndx(1,k) + sum(countN(k,1:nThread)) + cnt(k,nThread+1)
       !                  cnt(k,nThread+1) = cnt(k,nThread+1) + 1
       !                  exit
       !               end if
       !            end do
       !            rcvList(np) = NList2(j,i)
       !            NList2(j,i) = np + inode2
       !         end if
       !      end do
       !   end do
       !   !$OMP END PARALLEL
       !   call MIO_Allocate(nodeSndRcv,[nProc,nProc],'nodeSndRcv','neigh')
       !   do i=1,nProc
       !      nodeSndRcv(i,Node+1) = rcvIndx(2,i)-rcvIndx(1,i) + 1
       !   end do
       !   call MPIAllGather(MPI_IN_PLACE,0,nodeSndRcv,nProc,MPI_INTEGER)
       !   sndSz = sum(nodeSndRcv(Node+1,:))
       !   call MIO_Allocate(sndList,sndSz,'sndList','neigh')
       !   call MIO_Allocate(sndIndx,(/2,nProc/),'sndIndx','neigh')
       !   np = 1
       !   do i=1,nProc
       !      sndIndx(1,i) = np
       !      np = np + nodeSndRcv(Node+1,i)
       !      sndIndx(2,i) = np - 1
       !   end do
       !   do i=1,nProc-1
       !      j = mod(Node + i,nProc)
       !      k = mod(nProc + Node - i,nProc)
       !      call MPISendRecv(rcvList(rcvIndx(1,j+1):rcvIndx(2,j+1)),nodeSndRcv(j+1,Node+1),j,50, &
       !        sndList(sndIndx(1,k+1):sndIndx(2,k+1)),nodeSndRcv(Node+1,k+1),k,50,MPI_INTEGER)
       !      call MPISendRecv(rcvList(rcvIndx(1,k+1):rcvIndx(2,k+1)),nodeSndRcv(k+1,Node+1),k,51, &
       !        sndList(sndIndx(1,j+1):sndIndx(2,j+1)),nodeSndRcv(Node+1,j+1),j,51,MPI_INTEGER)
       !      call MPIBarrier()
       !   end do
       !   call MIO_Deallocate(cnt,'cnt','neigh')
       !   call MIO_Deallocate(countN,'countN','neigh')
       !end if
#endif /* MPI */
       call MIO_InputParameter('WriteDataFiles',prnt,.false.)
       !call AtomsSetFrac()
       if (prnt) then
          open(1,FILE='v')
          open(2,FILE='dx')
          open(3,FILE='dy')
          open(33,FILE='dz')
          open(4,FILE='pos')
          do i=1,nAt
             write(1,*) Nneigh(i)
             write(1,'(300(I8,1X))') (NList(j,i),j=1,Nneigh(i))
             write(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
             write(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
             write(33,*) (NeighD(3,j,i),j=1,Nneigh(i))
             if (Species(i)==3) then
                write(4,*) "B    ", Rat(1,i),Rat(2,i),Rat(3,i)
             else if (Species(i)==4) then
                write(4,*) "N    ", Rat(1,i),Rat(2,i),Rat(3,i)
             else
                write(4,*) "C    ", Rat(1,i),Rat(2,i),Rat(3,i)
             end if
          end do
          close(1)
          close(2)
          close(3)
          close(33)
          close(4)
       end if
   end if


end subroutine fastNNnotsquareSmall

subroutine fastNNnotsquareBulk(natoms,x,y,z,aCC,cutoff2,cutoff2bis,A1,A2,A3,maxnn)

   use atoms,                only : inode1, inode2, in1, in2, indxNode, nAt, Species, Rat
   use tbpar,                only : tbnn
   use cell,                 only : aG, aBN,ucell
   use atoms,                only : AtomsSetCart, frac
   use math
   implicit none
   integer, intent(in) :: natoms,maxnn
   real(dp), intent(in) :: x(natoms),y(natoms),z(natoms)
   integer :: nn(natoms,maxnn)
   real(dp), intent(in) :: aCC,cutoff2,cutoff2bis,A1(3),A2(3),A3
   integer :: num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,num10,num11,num12,num13,num14
   integer :: i,j,k,l,m,n
   real(dp) :: dx(natoms,maxnn),dy(natoms,maxnn),dz(natoms,maxnn)
   !real(dp), intent(out) :: dr(natoms,maxnn)

   real(dp) :: dx_t,dy_t,dz_t,d2_t, d2_t2
   real(dp) :: xmin,xmax,ymin,ymax,xcell,ycell,rho,zmin,zmax
   real(dp) :: zcell

   real(dp) :: xnew(9*natoms*3),ynew(9*natoms*3),znew(9*natoms*3)
   integer :: ixs(9*natoms*3),iys(9*natoms*3),vecino
   integer :: izs(9*natoms*3)
   integer :: nInteger(9*natoms*3)
   integer :: mInteger(9*natoms*3)

   integer :: Nx,Ny,ix,iy,iz
   integer :: Nz
   integer, allocatable :: cells(:,:,:)
   integer :: NNcount,addx,addy,addz

   !real(dp), parameter :: rad(3) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp)]
   real(dp), parameter :: rad(5) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp),3.0_dp/sqrt(3.0_dp),3.0_dp/sqrt(3.0_dp)]
   character(len=50) :: str

   integer :: ncell(3,9)
   real(dp) :: d, v(3), vtemp(3)
   real(dp) :: rmax, r0, rmaxTemp

   integer, pointer :: countN(:,:), cnt(:,:)
   integer :: np

   integer :: safetycounter
   logical :: prnt
   logical :: BLInPlaneInteractionRadius

   real(dp) :: distFact, dist

   integer :: sCell

   logical :: readDataFiles


   call MIO_InputParameter('TypeOfSystem',str,'Graphene')
   call MIO_InputParameter('ReadDataFiles',readDataFiles,.false.)
   !if (frac) call AtomsSetCart()
   !call AtomsSetFrac()
   if (readDataFiles) then
      call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')

      do i=1,tbnn
         Nradii(i,1) = rad(i)*aG*1.1_dp
      end do
      if (MIO_StringComp(str,'BoronNitride')) then
         Nradii(:,2) = Nradii(:,1)
      else
         do i=1,tbnn
            Nradii(i,2) = rad(i)*aBN*1.1_dp
         end do
      end if
      rmax = maxval(Nradii)

      maxNeigh = maxnn

      !print*, "allocating NList with maxNeigh= ", maxNeigh
      call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
      !call MIO_Allocate(NList,[1,inode1],[maxnn,inode2],'NList','neigh')
      call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
      call MIO_Allocate(NeighD,[1,1,inode1],[3,maxnn,inode2],'NeighD','neigh')
      !call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
      call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,maxnn,inode2/),'neighCell','neigh')
      open(1,FILE='v')
      open(2,FILE='dx')
      open(3,FILE='dy')
      !open(4,FILE='pos')
      !print*, "nAt =", nAt
      do i=1,nAt
         read(1,*) Nneigh(i)
         read(1,*) (NList(j,i),j=1,Nneigh(i))
         !read(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
         !read(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
      end do
      close(1)
      close(2)
      close(3)
      close(4)

      call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
      NList2 = NList

!`d + idx_1 * direction_1 + idx_2 * direction_2` is smaller than `d`, where
!`idx_1` and `idx_2` are both in `[-1,0,1] and `direction_1` and `direction_2`
!are the real space directions of the translational symmetry.`

      call MIO_InputParameter('SuperCell',sCell,1)
      !print*, "ucell", ucell
      !print*, Rat
      !print*, "jjhhh"
      !print*, Rat(:,1)
      !print*, Rat(:,NList(1,1))
      !print*, "mmmm"
      !if (sCell==1) then
         !!$OMP PARALLEL PRIVATE (v, d,ncell,i,j,ix,iy,dist)
         !do i=inode1,inode2
         do i=1,nAt
            !r0 = Nradii(tbnn,(Species(i)-1)/2 + 1) ! The cut radius. Taken as the maximum distance that the
            do j=1,Nneigh(i) ! different from ultraSmall because here we already know the neighbors
               do ix=-1,1; do iy=-1,1
                  !print*, "ixiy", ix, iy
                  if (i==NList(j,i) .and. ix==0 .and. iy==0) cycle
                  ncell(:,1) = [ix,iy,0]
                  !print*, "matmul: " , matmul(ucell,ncell(:,1))
                  v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1))
                  d = norm(v)
                  !print*, "d: ", d
                  if (d < aG/sqrt(3.0_dp)*1.1_dp) then
                     neighCell(:,j,i) = ncell(:,1)
                     NeighD(:,j,i) = -v
                  end if
                  !dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                  !dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                  !if (abs(d-dist).lt.0.1_dp) then
                  !   !print*, "muak"
                  !   neighCell(:,j,i) = ncell(:,1)
                  !end if

                  !v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition.
                                                                    ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
                  !d = norm(v)
                  !if (d**2.0_dp .lt. cutoff2*1.1_dp) then
                     !Nneigh(i) = Nneigh(i) + 1
                     !if (Nneigh(i)>maxNeigh)  then
                     !   !print*, "Number of neighbors: ", i, Nneigh(i)
                     !   call MIO_Kill('More neighbors than expected found',&
                     !     'neigh','NeighList')
                     !end if
                     !NList(Nneigh(i),i) = j
                     !NeighD(:,Nneigh(i),i) = -v
                     !neighCell(:,j,i) = ncell(:,1)
                  !else
                  !    print*, "Are you sure this was a neighbor?", i, NList(j,i)
                  !end if
               end do; end do
            end do
         end do
         !!$OMP END PARALLEL
      !end if
   else

       call MIO_InputParameter('Neigh.LayerDistFactor',distFact,1.0_dp) ! to increase the cutoff for outerlayer neighbors


       dx=0.0_dp;dy=0.0_dp;dz=0.0_dp !;dr=0.0_dp

       PRINT*,''
       PRINT*,'Building nearest neighbor data...'
       num0 = 0
       num1 = 0
       num2 = 0
       num3 = 0
       num4 = 0
       num5 = 0
       num6 = 0
       num7 = 0
       num8 = 0
       num9 = 0
       num10= 0
       num11= 0
       num12= 0
       num13= 0
       num14= 0

       ! Divide geometry into rectangular cells of size xcell times ycell (units of Angstroms)
       ! A_celdaunitaria=3.0*dsqrt(3.0)*aCC*aCC/2.0=2.6*aCC*aCC=5.24Ang => rho=2/A_celdaunitaria
       !if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer')) then
       !   rho   = 8.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)
       !else
       !end if
       xcell = 5.0_dp
       ycell = 5.0_dp
       zcell = A3
       !zcell = 5.0_dp
       !PRINT*,'Size of binning cell', xcell, 'X', ycell
       !============================================================
       !call MIO_Allocate(Species,nAt,'Species','atoms')
       !call MIO_Allocate(NeighD,[1,1,inode1],[3,maxNeigh,inode2],'NeighD','neigh')
       !call MIO_Allocate(xnew,[1],[natoms*9],'xnew','neigh')
       !call MIO_Allocate(ynew,[1],[natoms*9],'ynew','neigh')
       !call MIO_Allocate(znew,[1],[natoms*9],'znew','neigh')
       !allocate(xnew(9*natoms),ynew(9*natoms),znew(9*natoms))
       !print*, "numberOfAtoms: ", natoms, nAt
       do i=1,natoms
         xnew(i)=x(i); ynew(i)=y(i); znew(i)=z(i)
         nInteger(i) = 0
         mInteger(i) = 0
       end do
       do i=1,natoms
         n=i+natoms
         xnew(n)=x(i)+A1(1)
         ynew(n)=y(i)+A1(2)
         znew(n)=z(i)
         nInteger(n) = 1
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+2*natoms
         xnew(n)=x(i)+A1(1)+A2(1)
         ynew(n)=y(i)+A1(2)+A2(2)
         znew(n)=z(i)
         nInteger(n) = 1
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+3*natoms
         xnew(n)=x(i)+A2(1)
         ynew(n)=y(i)+A2(2)
         znew(n)=z(i)
         nInteger(n) = 0
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+4*natoms
         xnew(n)=x(i)-A1(1)+A2(1)
         ynew(n)=y(i)-A1(2)+A2(2)
         znew(n)=z(i)
         nInteger(n) = -1
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+5*natoms
         xnew(n)=x(i)-A1(1)
         ynew(n)=y(i)-A1(2)
         znew(n)=z(i)
         nInteger(n) = -1
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+6*natoms
         xnew(n)=x(i)-A1(1)-A2(1)
         ynew(n)=y(i)-A1(2)-A2(2)
         znew(n)=z(i)
         nInteger(n) = -1
         mInteger(n) = -1
       end do
       do i=1,natoms
         n=i+7*natoms
         xnew(n)=x(i)-A2(1)
         ynew(n)=y(i)-A2(2)
         znew(n)=z(i)
         nInteger(n) = 0
         mInteger(n) = -1
       end do
       do i=1,natoms
         n=i+8*natoms
         xnew(n)=x(i)+A1(1)-A2(1)
         ynew(n)=y(i)+A1(2)-A2(2)
         znew(n)=z(i)
         nInteger(n) = 1
         mInteger(n) = -1
       end do
       ! add on top for bulk
       !print*, "A3 =", A3
       !print*, "A1 =", A1
       do i=1,natoms
         n=i + 9*natoms
         xnew(n)=x(i); ynew(n)=y(i); znew(n)=z(i)+A3
         nInteger(n) = 0
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+natoms + 9*natoms
         xnew(n)=x(i)+A1(1)
         ynew(n)=y(i)+A1(2)
         znew(n)=z(i)+A3
         nInteger(n) = 1
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+2*natoms + 9*natoms
         xnew(n)=x(i)+A1(1)+A2(1)
         ynew(n)=y(i)+A1(2)+A2(2)
         znew(n)=z(i)+A3
         nInteger(n) = 1
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+3*natoms + 9*natoms
         xnew(n)=x(i)+A2(1)
         ynew(n)=y(i)+A2(2)
         znew(n)=z(i)+A3
         nInteger(n) = 0
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+4*natoms + 9*natoms
         xnew(n)=x(i)-A1(1)+A2(1)
         ynew(n)=y(i)-A1(2)+A2(2)
         znew(n)=z(i)+A3
         nInteger(n) = -1
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+5*natoms + 9*natoms
         xnew(n)=x(i)-A1(1)
         ynew(n)=y(i)-A1(2)
         znew(n)=z(i)+A3
         nInteger(n) = -1
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+6*natoms + 9*natoms
         xnew(n)=x(i)-A1(1)-A2(1)
         ynew(n)=y(i)-A1(2)-A2(2)
         znew(n)=z(i)+A3
         nInteger(n) = -1
         mInteger(n) = -1
       end do
       do i=1,natoms
         n=i+7*natoms + 9*natoms
         xnew(n)=x(i)-A2(1)
         ynew(n)=y(i)-A2(2)
         znew(n)=z(i)+A3
         nInteger(n) = 0
         mInteger(n) = -1
       end do
       do i=1,natoms
         n=i+8*natoms + 9*natoms
         xnew(n)=x(i)+A1(1)-A2(1)
         ynew(n)=y(i)+A1(2)-A2(2)
         znew(n)=z(i)+A3
         nInteger(n) = 1
         mInteger(n) = -1
       end do
       ! add below for bulk
       do i=1,natoms
         n=i + 9*natoms + 9*natoms
         xnew(n)=x(i); ynew(n)=y(i); znew(n)=z(i)-A3
         nInteger(n) = 0
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)+A1(1)
         ynew(n)=y(i)+A1(2)
         znew(n)=z(i)-A3
         nInteger(n) = 1
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+2*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)+A1(1)+A2(1)
         ynew(n)=y(i)+A1(2)+A2(2)
         znew(n)=z(i)-A3
         nInteger(n) = 1
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+3*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)+A2(1)
         ynew(n)=y(i)+A2(2)
         znew(n)=z(i)-A3
         nInteger(n) = 0
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+4*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)-A1(1)+A2(1)
         ynew(n)=y(i)-A1(2)+A2(2)
         znew(n)=z(i)-A3
         nInteger(n) = -1
         mInteger(n) = 1
       end do
       do i=1,natoms
         n=i+5*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)-A1(1)
         ynew(n)=y(i)-A1(2)
         znew(n)=z(i)-A3
         nInteger(n) = -1
         mInteger(n) = 0
       end do
       do i=1,natoms
         n=i+6*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)-A1(1)-A2(1)
         ynew(n)=y(i)-A1(2)-A2(2)
         znew(n)=z(i)-A3
         nInteger(n) = -1
         mInteger(n) = -1
       end do
       do i=1,natoms
         n=i+7*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)-A2(1)
         ynew(n)=y(i)-A2(2)
         znew(n)=z(i)-A3
         nInteger(n) = 0
         mInteger(n) = -1
       end do
       do i=1,natoms
         n=i+8*natoms + 9*natoms + 9*natoms
         xnew(n)=x(i)+A1(1)-A2(1)
         ynew(n)=y(i)+A1(2)-A2(2)
         znew(n)=z(i)-A3
         nInteger(n) = 1
         mInteger(n) = -1
       end do
       !============================================================
       xmin = minval(xnew); xmax = maxval(xnew)
       ymin = minval(ynew); ymax = maxval(ynew)
       zmin = minval(znew); zmax = maxval(znew)

       Nx = CEILING((xmax-xmin)/xcell)
       Ny = CEILING((ymax-ymin)/ycell)
       Nz = CEILING((zmax-zmin)/zcell)
       !Nz = CEILING((zmax-zmin)/zcell)

       !write(*,*) xmin,xmax,Nx
       !write(*,*) ymin,ymax,Ny
       !write(*,*) zmin,zmax,zy
       rho   = 4.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)*3.0d0
       ALLOCATE(cells(Nx, Ny, 16*CEILING(xcell*ycell*rho))) ! 2 instead of 4 if not bilayer
       !write(*,*) 'Inicia llenado de cajas'
       !write(*,*) 'First dimension of cells is equal to ', Nx
       !write(*,*) 'Second dimension of cells is equal to ', Ny
       !write(*,*) 'Third dimension of cells is equal to ', 2*CEILING(xcell*ycell*rho)
        cells = 0
       !print*, cells
       do n = 1,9*natoms*3 ! Find the index ix and iy of a binning cell and fill it with all the atoms that are inside there

         ix     = INT(1+(xnew(n)-xmin)/xcell)
         iy     = INT(1+(ynew(n)-ymin)/ycell)
         ixs(n) = ix
         iys(n) = iy
         izs(n) = INT(1+(znew(n)-zmin)/zcell)

         i = 1
         !print*, ix, iy, i
         do while (cells(ix,iy,i) .ne. 0)
           !print*, i, cells(ix,iy,i)
           i = i+1
         end do
         cells(ix,iy,i) = n
         ! HERE
         !ncell(:,1) = [ix,iy,0]
         !neighCell(:,Nneigh(i),i) = ncell(:,1)
       end do
       !print*, izs


       call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')

       do i=1,tbnn
          Nradii(i,1) = rad(i)*aG*1.1_dp
       end do
       if (MIO_StringComp(str,'BoronNitride')) then
          Nradii(:,2) = Nradii(:,1)
       else
          do i=1,tbnn
             Nradii(i,2) = rad(i)*aBN*1.1_dp
          end do
       end if
       rmax = maxval(Nradii)

       maxNeigh = maxnn
       !print*, "inode1, inode2, in1, in2, maxnn", inode1, inode2, 1, natoms, in1, in2, maxnn
       call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
       !call MIO_Allocate(NList,[1,inode1],[maxnn,inode2],'NList','neigh')
       call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
       call MIO_Allocate(NeighD,[1,1,inode1],[3,maxnn,inode2],'NeighD','neigh')
       !call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
       call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,maxnn,inode2/),'neighCell','neigh')
       !ALLOCATE(cells(Nx, Ny, 2*CEILING(xcell*ycell*rho)))
       !call MIO_Allocate(2*CEILING(xcell*ycell*rho),
       ! HERE

       write(*,*) 'Inicia busqueda'
       ! Find NNs of each atom by only searching nearby cells
       safetycounter = 0
       DO n = 1,natoms
         NNcount = 0

         DO addx = -2,2
         DO addy = -2,2
           ix = ixs(n) + addx
           iy = iys(n) + addy
           !iz = 0 + addz

           ! Periodic boundary conditions here
           IF(ix .gt. Nx) ix = ix-Nx
           IF(ix .lt. 1)  ix = ix+Nx
           IF(iy .gt. Ny) iy = iy-Ny
           IF(iy .lt. 1)  iy = iy+Ny
           !IF(iz .gt. Nz) iz = iz-Nz
           !IF(iz .lt. 1)  iz = iz+Nz

           i = 1
           DO WHILE(cells(ix, iy, i) .ne. 0)
             m = cells(ix, iy, i)
             IF (m .ne. n) THEN
               !write(*,*) 'hi2'
               dx_t = xnew(n) - xnew(m)
               dy_t = ynew(n) - ynew(m)
               dz_t = znew(n) - znew(m)
               ! Periodic boundary conditions here as well
               ! IF(dx_t .gt.  A1(1)-1.2_dp*aCC) dx_t = dx_t - A1(1)
               ! IF(dx_t .lt. -A1(1)+1.2_dp*aCC) dx_t = dx_t + A1(1)
               ! IF(dy_t .gt.  A2(2)-1.2_dp*aCC) dy_t = dy_t - A2(2)
               ! IF(dy_t .lt. -A2(2)+1.2_dp*aCC) dy_t = dy_t + A2(2)
               d2_t = dx_t**2.0_dp + dy_t**2.0_dp + dz_t**2.0_dp
               d2_t2 = dx_t**2.0_dp + dy_t**2.0_dp
                 !IF (n==1) print*, "distances from n=1 ", d2_t
                 !IF (n==3) print*, "distances from n=1 ", d2_t
               IF (abs(dz_t)<1.50_dp) THEN
                 !if (n .eq. 1) print*, d2_t2, cutoff2
                 IF (d2_t2 .lt. cutoff2*1.1_dp) THEN
                   NNcount = NNcount + 1
                   call NeighCheckCount(NNcount,maxnn)
                   vecino=mod(m,natoms)
                   if (vecino==0) vecino=natoms
                   nn(n,NNcount) = vecino
                   !NList(Nneigh(i),i) = j
                   dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
                   !dr(n,NNcount)=dsqrt(d2_t)
                   NList(NNcount,n) = vecino
                   NeighD(1,NNcount,n) = -dx(n,NNcount)
                   NeighD(2,NNcount,n) = -dy(n,NNcount)
                   NeighD(3,NNcount,n) = -dz(n,NNcount)
                   Nneigh(n) = NNcount
                   !neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
                   !neighCell(2,NNcount,n) = iys(m) - iys(n)
                   neighCell(1,NNcount,n) = nInteger(m) - nInteger(n)
                   neighCell(2,NNcount,n) = mInteger(m) - mInteger(n)
                   neighCell(3,NNcount,n) = izs(m) - izs(n)
                   !neighCell(3,NNcount,n) = izs(n) - izs(m)
       !             if(n==4800) write(*,*) m,NNcount,vecino
                 ENDIf
               ELSE IF (abs(dz_t)>1.50_dp .and. abs(dz_t)<4.50_dp) THEN
                 IF (d2_t2 .lt. cutoff2bis) THEN
                   NNcount = NNcount + 1
                   call NeighCheckCount(NNcount,maxnn)
                   vecino=mod(m,natoms)
                   if (vecino==0) vecino=natoms
                   nn(n,NNcount) = vecino
                   !NList(Nneigh(i),i) = j
                   dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
                   !dr(n,NNcount)=dsqrt(d2_t)
                   NList(NNcount,n) = vecino
                   NeighD(1,NNcount,n) = -dx(n,NNcount)
                   NeighD(2,NNcount,n) = -dy(n,NNcount)
                   NeighD(3,NNcount,n) = -dz(n,NNcount)
                   Nneigh(n) = NNcount
                   !neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
                   !neighCell(2,NNcount,n) = iys(m) - iys(n)
                   !neighCell(3,NNcount,n) = 0
                   neighCell(1,NNcount,n) = nInteger(m) - nInteger(n)
                   neighCell(2,NNcount,n) = mInteger(m) - mInteger(n)
                   neighCell(3,NNcount,n) = izs(m) - izs(n)
                   !neighCell(3,NNcount,n) = izs(n) - izs(m)
                   !neighCell(3,NNcount,n) =
       !             if(n==4800) write(*,*) m,NNcount,vecino
                 ENDIf
               ENDIF
             ENDIF
             i = i+1
           ENDDO

         ENDDO
         ENDDO

         !IF(NNcount .lt. 3) THEN
         !  print*, "Are you sure you wanted to enter this safety routine?"
         !  NNcount = 0
         !  DO m = 1,natoms
         !    IF (n.ne.m) THEN
         !      dx_t = x(n) - x(m)
         !      dy_t = y(n) - y(m)
         !      dz_t = z(n) - z(m)
         !       !Periodic boundary conditions here as well
         !       IF(dx_t .gt.  A1(1)-1.2_dp*aCC) dx_t = dx_t - A1(1)
         !       IF(dx_t .lt. -A1(1)+1.2_dp*aCC) dx_t = dx_t + A1(1)
         !       IF(dy_t .gt.  A2(2)-1.2_dp*aCC) dy_t = dy_t - A2(2)
         !       IF(dy_t .lt. -A2(2)+1.2_dp*aCC) dy_t = dy_t + A2(2)
         !       d2_t = dx_t**2.0_dp + dy_t**2.0_dp + dz_t**2.0_dp
         !       IF (abs(dz_t)<0.01_dp) THEN
         !          IF (d2_t .lt. cutoff2*1.1_dp) THEN
         !               NNcount = NNcount + 1
         !               !vecino=mod(m,natoms)
         !               !if (vecino==0) vecino=natoms
         !               !nn(n,NNcount) = vecino
         !               !NList(Nneigh(i),i) = j
         !               dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
         !               !dr(n,NNcount)=dsqrt(d2_t)
         !               NList(NNcount,n) = m
         !               NeighD(1,NNcount,n) = dx(n,NNcount)
         !               NeighD(2,NNcount,n) = dy(n,NNcount)
         !               NeighD(3,NNcount,n) = dz(n,NNcount)
         !               Nneigh(n) = NNcount
         !

         !               !ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
         !               !neighCell(:,Nneigh(i),i) = ncell(:,1)
         !

         !               !neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
         !               !neighCell(2,NNcount,n) = iys(m) - iys(n)
         !               !neighCell(3,NNcount,n) = 0
       ! !               if(n==4800) write(*,*) m,NNcount,vecino
         !          ENDIf
         !       END IF
         !     END IF
         !  END DO
         !  safetycounter = safetycounter + 1
         !  if (safetycounter.gt.100) print*, "Carefull, you are probably going too many times through this safety loop. &
         ! I originally implemented this because the fastNNnotsquared routine was missing neighbors for 4 atoms only..."
         !ENDIF


       ! Count the number of atoms with 1, 2, 3, 4, 5 nearest neighbors
         IF(NNcount .eq. 0) num0 = num0 + 1
         !IF(NNcount .eq. 0) print*, "NN counter = 0", n
         IF(NNcount .eq. 1) num1 = num1 + 1
         IF(NNcount .eq. 2) num2 = num2 + 1
         !IF(NNcount .eq. 2) print*, "NN counter = 2", n
         IF(NNcount .eq. 3) num3 = num3 + 1
         IF(NNcount .eq. 4) num4 = num4 + 1
         IF(NNcount .eq. 5) num5 = num5 + 1
         IF(NNcount .eq. 6) num6 = num6 + 1
         IF(NNcount .eq. 7) num7 = num7 + 1
         IF(NNcount .eq. 8) num8 = num8 + 1
         IF(NNcount .eq. 9) num9 = num9 + 1
         IF(NNcount .eq. 10) num10 = num10 + 1
         IF(NNcount .eq. 11) num11 = num11 + 1
         IF(NNcount .eq. 12) num12 = num12 + 1
         IF(NNcount .eq. 13) num13 = num13 + 1
         IF(NNcount .eq. 14) num14 = num14 + 1



       ENDDO

       !if (frac) call AtomsSetCart()
       !call MIO_InputParameter('SuperCell',sCell,1)
       !!if (sCell==1) then
       !   !!$OMP PARALLEL PRIVATE (v, d,ncell,i,j,ix,iy,dist)
       !   !do i=inode1,inode2
       !   do i=1,nAt
       !      !r0 = Nradii(tbnn,(Species(i)-1)/2 + 1) ! The cut radius. Taken as the maximum distance that the
       !      do j=1,Nneigh(i) ! different from ultraSmall because here we already know the neighbors
       !         do ix=-3,3; do iy=-3,3!; do iz=-1,1
       !            if (i==NList(j,i) .and. ix==0 .and. iy==0) cycle! .and. iz==0) cycle
       !            ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
       !            v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1))
       !            d = norm(v)
       !            dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
       !            if (abs(d-dist).lt.0.000001_dp) then
       !               !print*, "muak"
       !               !print*, i, j, NList(j,i), "iz", iz
       !               neighCell(1:2,j,i) = ncell(1:2,1)
       !            end if
       !            !v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition.
       !                                                               ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
       !            !d = norm(v)
       !            !if (d**2.0_dp .lt. cutoff2*1.1_dp) then
       !            !   !Nneigh(i) = Nneigh(i) + 1
       !            !   !if (Nneigh(i)>maxNeigh)  then
       !            !   !   !print*, "Number of neighbors: ", i, Nneigh(i)
       !            !   !   call MIO_Kill('More neighbors than expected found',&
       !            !   !     'neigh','NeighList')
       !            !   !end if
       !            !   !NList(Nneigh(i),i) = j
       !            !   !NeighD(:,Nneigh(i),i) = -v
       !            !   neighCell(:,j,i) = ncell(:,1)
       !            !!else
       !            !!    print*, "Are you sure this was a neighbor?", i, NList(j,i)
       !            !end if
       !         end do; end do!; end do
       !      end do
       !   end do
       !   !!$OMP END PARALLEL
       !!end if

       !write(*,*) 'hi4'
       !neighCell = cells
       DEALLOCATE(cells)
       !write(*,*) 'hi6'
       !call MIO_Deallocate(xnew,'xnew','neigh')
       !call MIO_Deallocate(ynew,'ynew','neigh')
       !call MIO_Deallocate(znew,'znew','neigh')
       !deallocate(xnew,ynew,znew)
       !write(*,*) 'hi5'


       PRINT*,' 0 neighbors: ',num0
       PRINT*,' 1 neighbors: ',num1
       PRINT*,' 2 neighbors: ',num2
       PRINT*,' 3 neighbors: ',num3
       PRINT*,' 4 neighbors: ',num4
       PRINT*,' 5 neighbors: ',num5
       PRINT*,' 6 neighbors: ',num6
       PRINT*,' 7 neighbors: ',num7
       PRINT*,' 8 neighbors: ',num8
       PRINT*,' 9 neighbors: ',num9
       PRINT*,' 10 neighbors: ',num10
       PRINT*,' 11 neighbors: ',num11
       PRINT*,' 12 neighbors: ',num12
       PRINT*,' 13 neighbors: ',num13
       PRINT*,' 14 neighbors: ',num14
       PRINT*,'>10 neighbors: ',natoms-(num0+num1+num2+num3+num4+num5+num6+num7+num8+num9+num10+num11+num12+num13+num14)

       PRINT*,''
       PRINT*,'...done'
       PRINT*,''

       call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
       NList2 = NList
#ifdef MPI
       !print*, "nProc =", nProc
       !if (nProc>1) then
       !   call MIO_Allocate(countN,[nProc,numThreads],'countN','neigh')
       !   countN = 0
       !   !do i=1,natoms
       !   !$OMP PARALLEL
       !   do i=in1,in2
       !      do j=1,maxNeigh
       !         if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
       !            do k=nProc,1,-1
       !               if(NList2(j,i)>=indxNode(k)) then
       !                  countN(k,nThread+1) = countN(k,nThread+1)+1
       !                  exit
       !               end if
       !            end do
       !         end if
       !      end do
       !   end do
       !   !$OMP END PARALLEL
       !   rcvSz = sum(countN)
       !   call MIO_Allocate(rcvIndx,[2,nProc],'rcvIndx','neigh')
       !   call MIO_Allocate(rcvList,rcvSz,'rcvList','neigh')
       !   np = 1
       !   do i=1,nProc
       !      rcvIndx(1,i) = np
       !      np = np + sum(countN(i,:))
       !      rcvIndx(2,i) = np - 1
       !   end do
       !   call MIO_Allocate(cnt,[nProc,numThreads],'cnt','neigh')
       !   !do i=1,natoms
       !   !$OMP PARALLEL PRIVATE(np)
       !   do i=in1,in2
       !      do j=1,maxNeigh
       !         if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
       !            do k=nProc,1,-1
       !               if(NList2(j,i)>=indxNode(k)) then
       !                  np = rcvIndx(1,k) + sum(countN(k,1:nThread)) + cnt(k,nThread+1)
       !                  cnt(k,nThread+1) = cnt(k,nThread+1) + 1
       !                  exit
       !               end if
       !            end do
       !            rcvList(np) = NList2(j,i)
       !            NList2(j,i) = np + inode2
       !         end if
       !      end do
       !   end do
       !   !$OMP END PARALLEL
       !   call MIO_Allocate(nodeSndRcv,[nProc,nProc],'nodeSndRcv','neigh')
       !   do i=1,nProc
       !      nodeSndRcv(i,Node+1) = rcvIndx(2,i)-rcvIndx(1,i) + 1
       !   end do
       !   call MPIAllGather(MPI_IN_PLACE,0,nodeSndRcv,nProc,MPI_INTEGER)
       !   sndSz = sum(nodeSndRcv(Node+1,:))
       !   call MIO_Allocate(sndList,sndSz,'sndList','neigh')
       !   call MIO_Allocate(sndIndx,(/2,nProc/),'sndIndx','neigh')
       !   np = 1
       !   do i=1,nProc
       !      sndIndx(1,i) = np
       !      np = np + nodeSndRcv(Node+1,i)
       !      sndIndx(2,i) = np - 1
       !   end do
       !   do i=1,nProc-1
       !      j = mod(Node + i,nProc)
       !      k = mod(nProc + Node - i,nProc)
       !      call MPISendRecv(rcvList(rcvIndx(1,j+1):rcvIndx(2,j+1)),nodeSndRcv(j+1,Node+1),j,50, &
       !        sndList(sndIndx(1,k+1):sndIndx(2,k+1)),nodeSndRcv(Node+1,k+1),k,50,MPI_INTEGER)
       !      call MPISendRecv(rcvList(rcvIndx(1,k+1):rcvIndx(2,k+1)),nodeSndRcv(k+1,Node+1),k,51, &
       !        sndList(sndIndx(1,j+1):sndIndx(2,j+1)),nodeSndRcv(Node+1,j+1),j,51,MPI_INTEGER)
       !      call MPIBarrier()
       !   end do
       !   call MIO_Deallocate(cnt,'cnt','neigh')
       !   call MIO_Deallocate(countN,'countN','neigh')
       !end if
#endif /* MPI */
       call MIO_InputParameter('WriteDataFiles',prnt,.false.)
       !call AtomsSetFrac()
       if (prnt) then
          open(1,FILE='v')
          open(2,FILE='dx')
          open(3,FILE='dy')
          open(33,FILE='dz')
          open(4,FILE='pos')
          do i=1,nAt
             write(1,*) Nneigh(i)
             !write(1,'(300(I8,1X))') (NList(j,i),j=1,Nneigh(i))
             !write(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
             !write(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
             !write(33,*) (NeighD(3,j,i),j=1,Nneigh(i))
             write(1,'(500(I8,1X))') (NList(j,i),j=1,Nneigh(i))
             write(2,'(500(F10.5,1X))') (NeighD(1,j,i),j=1,Nneigh(i))
             write(3,'(500(F10.5,1X))') (NeighD(2,j,i),j=1,Nneigh(i))
             write(33,'(500(F10.5,1X))') (NeighD(3,j,i),j=1,Nneigh(i))
             if (Species(i)==3) then
                write(4,*) "B    ", Rat(1,i),Rat(2,i),Rat(3,i)
             else if (Species(i)==4) then
                write(4,*) "N    ", Rat(1,i),Rat(2,i),Rat(3,i)
             else
                write(4,*) "C    ", Rat(1,i),Rat(2,i),Rat(3,i)
             end if
          end do
          close(1)
          close(2)
          close(3)
          close(33)
          close(4)
       end if
   end if


end subroutine fastNNnotsquareBulk

subroutine fastNNnotsquareBulkSmall(natoms,x,y,z,aCC,cutoff2,cutoff2bis,A1,A2,A3,maxnn)

   use atoms,                only : inode1, inode2, in1, in2, indxNode, nAt, Species, Rat
   use tbpar,                only : tbnn
   use cell,                 only : aG, aBN,ucell
   use atoms,                only : AtomsSetCart, frac
   use math
   implicit none
   integer, intent(in) :: natoms,maxnn
   real(dp), intent(in) :: x(natoms),y(natoms),z(natoms)
   integer :: nn(natoms,maxnn)
   real(dp), intent(in) :: aCC,cutoff2,cutoff2bis,A1(3),A2(3),A3
   integer :: num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,num10,num11,num12,num13,num14
   integer :: i,j,k,l,m,n
   real(dp) :: dx(natoms,maxnn),dy(natoms,maxnn),dz(natoms,maxnn)
   !real(dp), intent(out) :: dr(natoms,maxnn)

   real(dp) :: dx_t,dy_t,dz_t,d2_t, d2_t2
   real(dp) :: xmin,xmax,ymin,ymax,xcell,ycell,rho,zmin,zmax
   real(dp) :: zcell

   real(dp) :: xnew(11*11*natoms*3),ynew(11*11*natoms*3),znew(11*11*natoms*3)
   integer :: ixs(11*11*natoms*3),iys(11*11*natoms*3),vecino
   integer :: izs(11*11*natoms*3)
   integer :: nInteger(11*11*natoms*3)
   integer :: mInteger(11*11*natoms*3)

   integer :: Nx,Ny,ix,iy,iz
   integer :: Nz
   integer, allocatable :: cells(:,:,:)
   integer :: NNcount,addx,addy,addz

   !real(dp), parameter :: rad(3) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp)]
   real(dp), parameter :: rad(5) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp),3.0_dp/sqrt(3.0_dp),3.0_dp/sqrt(3.0_dp)]
   character(len=50) :: str

   integer :: ncell(3,9)
   real(dp) :: d, v(3), vtemp(3)
   real(dp) :: rmax, r0, rmaxTemp

   integer, pointer :: countN(:,:), cnt(:,:)
   integer :: np

   integer :: safetycounter
   logical :: prnt
   logical :: BLInPlaneInteractionRadius

   real(dp) :: distFact, dist

   integer :: sCell

   logical :: readDataFiles, addSecondLayerInteractions

   integer :: uX, uY, uZ, nmax, unitsX, unitsY


   call MIO_InputParameter('TypeOfSystem',str,'Graphene')
   call MIO_InputParameter('ReadDataFiles',readDataFiles,.false.)
   !if (frac) call AtomsSetCart()
   !call AtomsSetFrac()
   if (readDataFiles) then
      call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')

      do i=1,tbnn
         Nradii(i,1) = rad(i)*aG*1.1_dp
      end do
      if (MIO_StringComp(str,'BoronNitride')) then
         Nradii(:,2) = Nradii(:,1)
      else
         do i=1,tbnn
            Nradii(i,2) = rad(i)*aBN*1.1_dp
         end do
      end if
      rmax = maxval(Nradii)

      maxNeigh = maxnn

      !print*, "allocating NList with maxNeigh= ", maxNeigh
      call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
      !call MIO_Allocate(NList,[1,inode1],[maxnn,inode2],'NList','neigh')
      call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
      call MIO_Allocate(NeighD,[1,1,inode1],[3,maxnn,inode2],'NeighD','neigh')
      !call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
      call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,maxnn,inode2/),'neighCell','neigh')
      open(1,FILE='v')
      open(2,FILE='dx')
      open(3,FILE='dy')
      !open(4,FILE='pos')
      !print*, "nAt =", nAt
      do i=1,nAt
         read(1,*) Nneigh(i)
         read(1,*) (NList(j,i),j=1,Nneigh(i))
         !read(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
         !read(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
         !read(4,*) Species(i), Rat(1,i),Rat(2,i),Rat(3,i)
      end do
      close(1)
      close(2)
      close(3)
      close(4)

      call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
      NList2 = NList

!`d + idx_1 * direction_1 + idx_2 * direction_2` is smaller than `d`, where
!`idx_1` and `idx_2` are both in `[-1,0,1] and `direction_1` and `direction_2`
!are the real space directions of the translational symmetry.`

      call MIO_InputParameter('SuperCell',sCell,1)
      !print*, "ucell", ucell
      !print*, Rat
      !print*, "jjhhh"
      !print*, Rat(:,1)
      !print*, Rat(:,NList(1,1))
      !print*, "mmmm"
      !if (sCell==1) then
         !!$OMP PARALLEL PRIVATE (v, d,ncell,i,j,ix,iy,dist)
         !do i=inode1,inode2
         do i=1,nAt
            !r0 = Nradii(tbnn,(Species(i)-1)/2 + 1) ! The cut radius. Taken as the maximum distance that the
            do j=1,Nneigh(i) ! different from ultraSmall because here we already know the neighbors
               do ix=-1,1; do iy=-1,1
                  !print*, "ixiy", ix, iy
                  if (i==NList(j,i) .and. ix==0 .and. iy==0) cycle
                  ncell(:,1) = [ix,iy,0]
                  !print*, "matmul: " , matmul(ucell,ncell(:,1))
                  v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1))
                  d = norm(v)
                  !print*, "d: ", d
                  if (d < aG/sqrt(3.0_dp)*1.1_dp) then
                     neighCell(:,j,i) = ncell(:,1)
                     NeighD(:,j,i) = -v
                  end if
                  !dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                  !dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
                  !if (abs(d-dist).lt.0.1_dp) then
                  !   !print*, "muak"
                  !   neighCell(:,j,i) = ncell(:,1)
                  !end if

                  !v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition.
                                                                    ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
                  !d = norm(v)
                  !if (d**2.0_dp .lt. cutoff2*1.1_dp) then
                     !Nneigh(i) = Nneigh(i) + 1
                     !if (Nneigh(i)>maxNeigh)  then
                     !   !print*, "Number of neighbors: ", i, Nneigh(i)
                     !   call MIO_Kill('More neighbors than expected found',&
                     !     'neigh','NeighList')
                     !end if
                     !NList(Nneigh(i),i) = j
                     !NeighD(:,Nneigh(i),i) = -v
                     !neighCell(:,j,i) = ncell(:,1)
                  !else
                  !    print*, "Are you sure this was a neighbor?", i, NList(j,i)
                  !end if
               end do; end do
            end do
         end do
         !!$OMP END PARALLEL
      !end if
   else

       call MIO_InputParameter('Neigh.LayerDistFactor',distFact,1.0_dp) ! to increase the cutoff for outerlayer neighbors
       call MIO_InputParameter('addSecondLayerInteractions',addSecondLayerInteractions,.false.) ! to increase the cutoff for outerlayer neighbors


       dx=0.0_dp;dy=0.0_dp;dz=0.0_dp !;dr=0.0_dp

       PRINT*,''
       PRINT*,'Building nearest neighbor data...'
       num0 = 0
       num1 = 0
       num2 = 0
       num3 = 0
       num4 = 0
       num5 = 0
       num6 = 0
       num7 = 0
       num8 = 0
       num9 = 0
       num10= 0
       num11= 0
       num12= 0
       num13= 0
       num14= 0

       ! Divide geometry into rectangular cells of size xcell times ycell (units of Angstroms)
       ! A_celdaunitaria=3.0*dsqrt(3.0)*aCC*aCC/2.0=2.6*aCC*aCC=5.24Ang => rho=2/A_celdaunitaria
       !if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer')) then
       !   rho   = 8.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)
       !else
       !end if
       xcell = 5.0_dp
       ycell = 5.0_dp
       zcell = A3
       !zcell = 5.0_dp
       !PRINT*,'Size of binning cell', xcell, 'X', ycell
       !============================================================
       !call MIO_Allocate(Species,nAt,'Species','atoms')
       !call MIO_Allocate(NeighD,[1,1,inode1],[3,maxNeigh,inode2],'NeighD','neigh')
       !call MIO_Allocate(xnew,[1],[natoms*9],'xnew','neigh')
       !call MIO_Allocate(ynew,[1],[natoms*9],'ynew','neigh')
       !call MIO_Allocate(znew,[1],[natoms*9],'znew','neigh')
       !allocate(xnew(9*natoms),ynew(9*natoms),znew(9*natoms))
       !print*, "numberOfAtoms: ", natoms, nAt
       unitsX = 10
       unitsY = 10
       n = 0
       do i=1,natoms
         n=n+1
         uX = 0
         uY = 0
         uZ = 0
         xnew(n)=x(i)+uX*A1(1)+uY*A2(1)
         ynew(n)=y(i)+uX*A1(2)+uY*A2(2)
         znew(n)=z(i)+uZ*A3
         nInteger(n) = uX
         mInteger(n) = uY
       end do
       do uX=-unitsX/2,unitsX/2
          do uY=-unitsY/2,unitsY/2
             do uZ=-1,1
                do i=1,natoms
                  if (uX==0 .and. uY==0 .and. uZ==0) cycle
                  n=n+1
                  xnew(n)=x(i)+uX*A1(1)+uY*A2(1)
                  ynew(n)=y(i)+uX*A1(2)+uY*A2(2)
                  znew(n)=z(i)+uZ*A3
                  nInteger(n) = uX
                  mInteger(n) = uY
                end do
             end do
          end do
       end do
       nmax = n
       !print*, "nmax", nmax
       !============================================================
       xmin = minval(xnew); xmax = maxval(xnew)
       ymin = minval(ynew); ymax = maxval(ynew)
       zmin = minval(znew); zmax = maxval(znew)

       Nx = CEILING((xmax-xmin)/xcell)
       Ny = CEILING((ymax-ymin)/ycell)
       Nz = CEILING((zmax-zmin)/zcell)
       !Nz = CEILING((zmax-zmin)/zcell)

       !write(*,*) xmin,xmax,Nx
       !write(*,*) ymin,ymax,Ny
       !write(*,*) zmin,zmax,zy
       rho   = 4.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)*3.0d0
       ALLOCATE(cells(Nx, Ny, 4*CEILING(xcell*ycell*rho*unitsX*unitsY))) ! 2 instead of 4 if not bilayer
       !write(*,*) 'Inicia llenado de cajas'
       !write(*,*) 'First dimension of cells is equal to ', Nx
       !write(*,*) 'Second dimension of cells is equal to ', Ny
       !write(*,*) 'Third dimension of cells is equal to ', 2*CEILING(xcell*ycell*rho)
        cells = 0
       !print*, cells
       do n = 1,nmax ! Find the index ix and iy of a binning cell and fill it with all the atoms that are inside there

         ix     = INT(1+(xnew(n)-xmin)/xcell)
         iy     = INT(1+(ynew(n)-ymin)/ycell)
         ixs(n) = ix
         iys(n) = iy
         izs(n) = INT(1+(znew(n)-zmin)/zcell)

         i = 1
         !print*, ix, iy, i
         do while (cells(ix,iy,i) .ne. 0)
           !print*, i, cells(ix,iy,i)
           i = i+1
         end do
         cells(ix,iy,i) = n
         ! HERE
         !ncell(:,1) = [ix,iy,0]
         !neighCell(:,Nneigh(i),i) = ncell(:,1)
       end do
       !print*, izs


       call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')

       do i=1,tbnn
          Nradii(i,1) = rad(i)*aG*1.1_dp
       end do
       if (MIO_StringComp(str,'BoronNitride')) then
          Nradii(:,2) = Nradii(:,1)
       else
          do i=1,tbnn
             Nradii(i,2) = rad(i)*aBN*1.1_dp
          end do
       end if
       rmax = maxval(Nradii)

       maxNeigh = maxnn
       !print*, "inode1, inode2, in1, in2, maxnn", inode1, inode2, 1, natoms, in1, in2, maxnn
       call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
       !call MIO_Allocate(NList,[1,inode1],[maxnn,inode2],'NList','neigh')
       call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
       call MIO_Allocate(NeighD,[1,1,inode1],[3,maxnn,inode2],'NeighD','neigh')
       !call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
       call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,maxnn,inode2/),'neighCell','neigh')
       !ALLOCATE(cells(Nx, Ny, 2*CEILING(xcell*ycell*rho)))
       !call MIO_Allocate(2*CEILING(xcell*ycell*rho),
       ! HERE

       write(*,*) 'Inicia busqueda'
       ! Find NNs of each atom by only searching nearby cells
       safetycounter = 0
       DO n = 1,natoms
         NNcount = 0

         DO addx = -2,2
         DO addy = -2,2
           ix = ixs(n) + addx
           iy = iys(n) + addy
           !iz = 0 + addz

           ! Periodic boundary conditions here
           IF(ix .gt. Nx) ix = ix-Nx
           IF(ix .lt. 1)  ix = ix+Nx
           IF(iy .gt. Ny) iy = iy-Ny
           IF(iy .lt. 1)  iy = iy+Ny
           !IF(iz .gt. Nz) iz = iz-Nz
           !IF(iz .lt. 1)  iz = iz+Nz

           i = 1
           DO WHILE(cells(ix, iy, i) .ne. 0)
             m = cells(ix, iy, i)
             IF (m .ne. n) THEN
               !write(*,*) 'hi2'
               dx_t = xnew(n) - xnew(m)
               dy_t = ynew(n) - ynew(m)
               dz_t = znew(n) - znew(m)
               ! Periodic boundary conditions here as well
               ! IF(dx_t .gt.  A1(1)-1.2_dp*aCC) dx_t = dx_t - A1(1)
               ! IF(dx_t .lt. -A1(1)+1.2_dp*aCC) dx_t = dx_t + A1(1)
               ! IF(dy_t .gt.  A2(2)-1.2_dp*aCC) dy_t = dy_t - A2(2)
               ! IF(dy_t .lt. -A2(2)+1.2_dp*aCC) dy_t = dy_t + A2(2)
               d2_t = dx_t**2.0_dp + dy_t**2.0_dp + dz_t**2.0_dp
               d2_t2 = dx_t**2.0_dp + dy_t**2.0_dp
                 !IF (n==1) print*, "distances from n=1 ", d2_t
                 !IF (n==3) print*, "distances from n=1 ", d2_t
               IF (abs(dz_t)<1.50_dp) THEN
                 !if (n .eq. 1) print*, d2_t2, cutoff2
                 IF (d2_t2 .lt. cutoff2*1.1_dp) THEN
                   NNcount = NNcount + 1
                   call NeighCheckCount(NNcount,maxnn)
                   vecino=mod(m,natoms)
                   if (vecino==0) vecino=natoms
                   nn(n,NNcount) = vecino
                   !NList(Nneigh(i),i) = j
                   dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
                   !dr(n,NNcount)=dsqrt(d2_t)
                   NList(NNcount,n) = vecino
                   NeighD(1,NNcount,n) = -dx(n,NNcount)
                   NeighD(2,NNcount,n) = -dy(n,NNcount)
                   NeighD(3,NNcount,n) = -dz(n,NNcount)
                   Nneigh(n) = NNcount
                   !neighCell(1,NNcount,n) = ixs(n) - ixs(m) ! check if it's not n-m
                   !neighCell(2,NNcount,n) = iys(n) - iys(m)
                   neighCell(1,NNcount,n) = nInteger(m) - nInteger(n)
                   neighCell(2,NNcount,n) = mInteger(m) - mInteger(n)
                   neighCell(3,NNcount,n) = izs(m) - izs(n)
                   !neighCell(3,NNcount,n) = izs(n) - izs(m)
       !             if(n==4800) write(*,*) m,NNcount,vecino
                 ENDIf
               ELSE IF (abs(dz_t)>1.50_dp .and. abs(dz_t)<4.50_dp) THEN
                 IF (d2_t2 .lt. cutoff2bis) THEN
                   NNcount = NNcount + 1
                   call NeighCheckCount(NNcount,maxnn)
                   vecino=mod(m,natoms)
                   if (vecino==0) vecino=natoms
                   nn(n,NNcount) = vecino
                   !NList(Nneigh(i),i) = j
                   dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
                   !dr(n,NNcount)=dsqrt(d2_t)
                   NList(NNcount,n) = vecino
                   NeighD(1,NNcount,n) = -dx(n,NNcount)
                   NeighD(2,NNcount,n) = -dy(n,NNcount)
                   NeighD(3,NNcount,n) = -dz(n,NNcount)
                   Nneigh(n) = NNcount
                   !neighCell(1,NNcount,n) = ixs(n) - ixs(m) ! check if it's not n-m
                   !neighCell(2,NNcount,n) = iys(n) - iys(m)
                   !neighCell(3,NNcount,n) = 0
                   neighCell(1,NNcount,n) = nInteger(m) - nInteger(n)
                   neighCell(2,NNcount,n) = mInteger(m) - mInteger(n)
                   neighCell(3,NNcount,n) = izs(m) - izs(n)
                   !neighCell(3,NNcount,n) = izs(n) - izs(m)
                   !neighCell(3,NNcount,n) =
       !             if(n==4800) write(*,*) m,NNcount,vecino
                 ENDIf
               ELSE IF (abs(dz_t)>4.50_dp .and. abs(dz_t)<7.50_dp .and. addSecondLayerInteractions) THEN
                 IF (d2_t2 .lt. cutoff2bis) THEN
                   NNcount = NNcount + 1
                   call NeighCheckCount(NNcount,maxnn)
                   vecino=mod(m,natoms)
                   if (vecino==0) vecino=natoms
                   nn(n,NNcount) = vecino
                   !NList(Nneigh(i),i) = j
                   dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
                   !dr(n,NNcount)=dsqrt(d2_t)
                   NList(NNcount,n) = vecino
                   NeighD(1,NNcount,n) = -dx(n,NNcount)
                   NeighD(2,NNcount,n) = -dy(n,NNcount)
                   NeighD(3,NNcount,n) = -dz(n,NNcount)
                   Nneigh(n) = NNcount
                   !neighCell(1,NNcount,n) = ixs(n) - ixs(m) ! check if it's not n-m
                   !neighCell(2,NNcount,n) = iys(n) - iys(m)
                   !neighCell(3,NNcount,n) = 0
                   neighCell(1,NNcount,n) = nInteger(m) - nInteger(n)
                   neighCell(2,NNcount,n) = mInteger(m) - mInteger(n)
                   neighCell(3,NNcount,n) = izs(m) - izs(n)
                   !neighCell(3,NNcount,n) = izs(n) - izs(m)
                   !neighCell(3,NNcount,n) =
       !             if(n==4800) write(*,*) m,NNcount,vecino
                 ENDIf
               ENDIF
             ENDIF
             i = i+1
           ENDDO

         ENDDO
         ENDDO

         !IF(NNcount .lt. 3) THEN
         !  print*, "Are you sure you wanted to enter this safety routine?"
         !  NNcount = 0
         !  DO m = 1,natoms
         !    IF (n.ne.m) THEN
         !      dx_t = x(n) - x(m)
         !      dy_t = y(n) - y(m)
         !      dz_t = z(n) - z(m)
         !       !Periodic boundary conditions here as well
         !       IF(dx_t .gt.  A1(1)-1.2_dp*aCC) dx_t = dx_t - A1(1)
         !       IF(dx_t .lt. -A1(1)+1.2_dp*aCC) dx_t = dx_t + A1(1)
         !       IF(dy_t .gt.  A2(2)-1.2_dp*aCC) dy_t = dy_t - A2(2)
         !       IF(dy_t .lt. -A2(2)+1.2_dp*aCC) dy_t = dy_t + A2(2)
         !       d2_t = dx_t**2.0_dp + dy_t**2.0_dp + dz_t**2.0_dp
         !       IF (abs(dz_t)<0.01_dp) THEN
         !          IF (d2_t .lt. cutoff2*1.1_dp) THEN
         !               NNcount = NNcount + 1
         !               !vecino=mod(m,natoms)
         !               !if (vecino==0) vecino=natoms
         !               !nn(n,NNcount) = vecino
         !               !NList(Nneigh(i),i) = j
         !               dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
         !               !dr(n,NNcount)=dsqrt(d2_t)
         !               NList(NNcount,n) = m
         !               NeighD(1,NNcount,n) = dx(n,NNcount)
         !               NeighD(2,NNcount,n) = dy(n,NNcount)
         !               NeighD(3,NNcount,n) = dz(n,NNcount)
         !               Nneigh(n) = NNcount
         !

         !               !ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
         !               !neighCell(:,Nneigh(i),i) = ncell(:,1)
         !

         !               !neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
         !               !neighCell(2,NNcount,n) = iys(m) - iys(n)
         !               !neighCell(3,NNcount,n) = 0
       ! !               if(n==4800) write(*,*) m,NNcount,vecino
         !          ENDIf
         !       END IF
         !     END IF
         !  END DO
         !  safetycounter = safetycounter + 1
         !  if (safetycounter.gt.100) print*, "Carefull, you are probably going too many times through this safety loop. &
         ! I originally implemented this because the fastNNnotsquared routine was missing neighbors for 4 atoms only..."
         !ENDIF


       ! Count the number of atoms with 1, 2, 3, 4, 5 nearest neighbors
         IF(NNcount .eq. 0) num0 = num0 + 1
         !IF(NNcount .eq. 0) print*, "NN counter = 0", n
         IF(NNcount .eq. 1) num1 = num1 + 1
         IF(NNcount .eq. 2) num2 = num2 + 1
         !IF(NNcount .eq. 2) print*, "NN counter = 2", n
         IF(NNcount .eq. 3) num3 = num3 + 1
         IF(NNcount .eq. 4) num4 = num4 + 1
         IF(NNcount .eq. 5) num5 = num5 + 1
         IF(NNcount .eq. 6) num6 = num6 + 1
         IF(NNcount .eq. 7) num7 = num7 + 1
         IF(NNcount .eq. 8) num8 = num8 + 1
         IF(NNcount .eq. 9) num9 = num9 + 1
         IF(NNcount .eq. 10) num10 = num10 + 1
         IF(NNcount .eq. 11) num11 = num11 + 1
         IF(NNcount .eq. 12) num12 = num12 + 1
         IF(NNcount .eq. 13) num13 = num13 + 1
         IF(NNcount .eq. 14) num14 = num14 + 1



       ENDDO

       !if (frac) call AtomsSetCart()
       !call MIO_InputParameter('SuperCell',sCell,1)
       !!if (sCell==1) then
       !   !!$OMP PARALLEL PRIVATE (v, d,ncell,i,j,ix,iy,dist)
       !   !do i=inode1,inode2
       !   do i=1,nAt
       !      !r0 = Nradii(tbnn,(Species(i)-1)/2 + 1) ! The cut radius. Taken as the maximum distance that the
       !      do j=1,Nneigh(i) ! different from ultraSmall because here we already know the neighbors
       !         do ix=-5,5; do iy=-5,5!; do iz=-1,1
       !            if (i==NList(j,i) .and. ix==0 .and. iy==0 .and. iz==0) cycle
       !            ncell(:,1) = [ix,iy,iz] ! We only use one of the nine columns if small, the other columns are used for the other case
       !            v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1))
       !            d = norm(v)
       !            dist = sqrt(NeighD(1,j,i)**2.0_dp+NeighD(2,j,i)**2.0_dp +NeighD(3,j,i)**2.0_dp)
       !            if (abs(d-dist).lt.0.000001_dp) then
       !               !print*, "muak"
       !               !print*, i, j, NList(j,i), "iz", iz
       !               neighCell(1:2,j,i) = ncell(1:2,1)
       !            end if
       !            !v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition.
       !                                                               ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
       !            !d = norm(v)
       !            !if (d**2.0_dp .lt. cutoff2*1.1_dp) then
       !            !   !Nneigh(i) = Nneigh(i) + 1
       !            !   !if (Nneigh(i)>maxNeigh)  then
       !            !   !   !print*, "Number of neighbors: ", i, Nneigh(i)
       !            !   !   call MIO_Kill('More neighbors than expected found',&
       !            !   !     'neigh','NeighList')
       !            !   !end if
       !            !   !NList(Nneigh(i),i) = j
       !            !   !NeighD(:,Nneigh(i),i) = -v
       !            !   neighCell(:,j,i) = ncell(:,1)
       !            !!else
       !            !!    print*, "Are you sure this was a neighbor?", i, NList(j,i)
       !            !end if
       !         end do; end do!; end do
       !      end do
       !   end do
       !   !!$OMP END PARALLEL
       !!end if

       !write(*,*) 'hi4'
       !neighCell = cells
       DEALLOCATE(cells)
       !write(*,*) 'hi6'
       !call MIO_Deallocate(xnew,'xnew','neigh')
       !call MIO_Deallocate(ynew,'ynew','neigh')
       !call MIO_Deallocate(znew,'znew','neigh')
       !deallocate(xnew,ynew,znew)
       !write(*,*) 'hi5'


       PRINT*,' 0 neighbors: ',num0
       PRINT*,' 1 neighbors: ',num1
       PRINT*,' 2 neighbors: ',num2
       PRINT*,' 3 neighbors: ',num3
       PRINT*,' 4 neighbors: ',num4
       PRINT*,' 5 neighbors: ',num5
       PRINT*,' 6 neighbors: ',num6
       PRINT*,' 7 neighbors: ',num7
       PRINT*,' 8 neighbors: ',num8
       PRINT*,' 9 neighbors: ',num9
       PRINT*,' 10 neighbors: ',num10
       PRINT*,' 11 neighbors: ',num11
       PRINT*,' 12 neighbors: ',num12
       PRINT*,' 13 neighbors: ',num13
       PRINT*,' 14 neighbors: ',num14
       PRINT*,'>10 neighbors: ',natoms-(num0+num1+num2+num3+num4+num5+num6+num7+num8+num9+num10+num11+num12+num13+num14)

       PRINT*,''
       PRINT*,'...done'
       PRINT*,''

       call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
       NList2 = NList
#ifdef MPI
       !print*, "nProc =", nProc
       !if (nProc>1) then
       !   call MIO_Allocate(countN,[nProc,numThreads],'countN','neigh')
       !   countN = 0
       !   !do i=1,natoms
       !   !$OMP PARALLEL
       !   do i=in1,in2
       !      do j=1,maxNeigh
       !         if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
       !            do k=nProc,1,-1
       !               if(NList2(j,i)>=indxNode(k)) then
       !                  countN(k,nThread+1) = countN(k,nThread+1)+1
       !                  exit
       !               end if
       !            end do
       !         end if
       !      end do
       !   end do
       !   !$OMP END PARALLEL
       !   rcvSz = sum(countN)
       !   call MIO_Allocate(rcvIndx,[2,nProc],'rcvIndx','neigh')
       !   call MIO_Allocate(rcvList,rcvSz,'rcvList','neigh')
       !   np = 1
       !   do i=1,nProc
       !      rcvIndx(1,i) = np
       !      np = np + sum(countN(i,:))
       !      rcvIndx(2,i) = np - 1
       !   end do
       !   call MIO_Allocate(cnt,[nProc,numThreads],'cnt','neigh')
       !   !do i=1,natoms
       !   !$OMP PARALLEL PRIVATE(np)
       !   do i=in1,in2
       !      do j=1,maxNeigh
       !         if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
       !            do k=nProc,1,-1
       !               if(NList2(j,i)>=indxNode(k)) then
       !                  np = rcvIndx(1,k) + sum(countN(k,1:nThread)) + cnt(k,nThread+1)
       !                  cnt(k,nThread+1) = cnt(k,nThread+1) + 1
       !                  exit
       !               end if
       !            end do
       !            rcvList(np) = NList2(j,i)
       !            NList2(j,i) = np + inode2
       !         end if
       !      end do
       !   end do
       !   !$OMP END PARALLEL
       !   call MIO_Allocate(nodeSndRcv,[nProc,nProc],'nodeSndRcv','neigh')
       !   do i=1,nProc
       !      nodeSndRcv(i,Node+1) = rcvIndx(2,i)-rcvIndx(1,i) + 1
       !   end do
       !   call MPIAllGather(MPI_IN_PLACE,0,nodeSndRcv,nProc,MPI_INTEGER)
       !   sndSz = sum(nodeSndRcv(Node+1,:))
       !   call MIO_Allocate(sndList,sndSz,'sndList','neigh')
       !   call MIO_Allocate(sndIndx,(/2,nProc/),'sndIndx','neigh')
       !   np = 1
       !   do i=1,nProc
       !      sndIndx(1,i) = np
       !      np = np + nodeSndRcv(Node+1,i)
       !      sndIndx(2,i) = np - 1
       !   end do
       !   do i=1,nProc-1
       !      j = mod(Node + i,nProc)
       !      k = mod(nProc + Node - i,nProc)
       !      call MPISendRecv(rcvList(rcvIndx(1,j+1):rcvIndx(2,j+1)),nodeSndRcv(j+1,Node+1),j,50, &
       !        sndList(sndIndx(1,k+1):sndIndx(2,k+1)),nodeSndRcv(Node+1,k+1),k,50,MPI_INTEGER)
       !      call MPISendRecv(rcvList(rcvIndx(1,k+1):rcvIndx(2,k+1)),nodeSndRcv(k+1,Node+1),k,51, &
       !        sndList(sndIndx(1,j+1):sndIndx(2,j+1)),nodeSndRcv(Node+1,j+1),j,51,MPI_INTEGER)
       !      call MPIBarrier()
       !   end do
       !   call MIO_Deallocate(cnt,'cnt','neigh')
       !   call MIO_Deallocate(countN,'countN','neigh')
       !end if
#endif /* MPI */
       call MIO_InputParameter('WriteDataFiles',prnt,.false.)
       !call AtomsSetFrac()
       if (prnt) then
          open(1,FILE='v')
          open(2,FILE='dx')
          open(3,FILE='dy')
          open(33,FILE='dz')
          open(4,FILE='pos')
          do i=1,nAt
             write(1,*) Nneigh(i)
             !write(1,'(300(I8,1X))') (NList(j,i),j=1,Nneigh(i))
             !write(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
             !write(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
             !write(33,*) (NeighD(3,j,i),j=1,Nneigh(i))
             write(1,'(500(I8,1X))') (NList(j,i),j=1,Nneigh(i))
             write(2,'(500(F10.5,1X))') (NeighD(1,j,i),j=1,Nneigh(i))
             write(3,'(500(F10.5,1X))') (NeighD(2,j,i),j=1,Nneigh(i))
             write(33,'(500(F10.5,1X))') (NeighD(3,j,i),j=1,Nneigh(i))
             if (Species(i)==3) then
                write(4,*) "B    ", Rat(1,i),Rat(2,i),Rat(3,i)
             else if (Species(i)==4) then
                write(4,*) "N    ", Rat(1,i),Rat(2,i),Rat(3,i)
             else
                write(4,*) "C    ", Rat(1,i),Rat(2,i),Rat(3,i)
             end if
          end do
          close(1)
          close(2)
          close(3)
          close(33)
          close(4)
       end if
   end if


end subroutine fastNNnotsquareBulkSmall

subroutine fastNNnotsquareNotRectangle(natoms,x,y,z,aCC,cutoff2,cutoff2bis,A1,A2,maxnn)

   use atoms,                only : inode1, inode2, in1, in2, indxNode, nAt, Species, Rat
   use atoms,                only : frac, AtomsSetCart
   use tbpar,                only : tbnn
   use cell,                 only : aG, aBN,ucell
   use math
   implicit none
   integer, intent(in) :: natoms,maxnn
   real(dp), intent(in) :: x(natoms),y(natoms),z(natoms)
   integer :: nn(natoms,maxnn)
   real(dp) :: aCC,cutoff2,cutoff2bis,A1(3),A2(3)
   integer :: num0,num1,num2,num3,num4,num5,num6,num7,num8,num9,num10
   integer :: i,j,k,l,m,n
   real(dp) :: dx(natoms,maxnn),dy(natoms,maxnn),dz(natoms,maxnn)
   !real(dp), intent(out) :: dr(natoms,maxnn)

   real(dp) :: dx_t,dy_t,dz_t,d2_t, d2_t2
   real(dp) :: xmin,xmax,ymin,ymax,xcell,ycell,rho

   real(dp) :: xnew(9*natoms),ynew(9*natoms),znew(9*natoms)
   integer :: ixs(9*natoms),iys(9*natoms),vecino

   integer :: Nx,Ny,ix,iy
   integer, allocatable :: cells(:,:,:)
   integer :: NNcount,addx,addy

   real(dp), parameter :: rad(3) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp)]
   character(len=50) :: str

   integer :: ncell(3,9)
   real(dp) :: d, v(3)
   real(dp) :: rmax

   integer, pointer :: countN(:,:), cnt(:,:)
   integer :: np

   integer :: safetycounter
   logical :: prnt

   real(dp) :: fracFactor, fracFactor2


   dx=0.0_dp;dy=0.0_dp;dz=0.0_dp !;dr=0.0_dp

   PRINT*,''
   PRINT*,'Building nearest neighbor data...'
   num0 = 0
   num1 = 0
   num2 = 0
   num3 = 0
   num4 = 0
   num5 = 0
   num6 = 0
   num7 = 0
   num8 = 0
   num9 = 0
   num10= 0

   ! Divide geometry into rectangular cells of size xcell times ycell (units of Angstroms)
   ! A_celdaunitaria=3.0*dsqrt(3.0)*aCC*aCC/2.0=2.6*aCC*aCC=5.24Ang => rho=2/A_celdaunitaria
   !if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer')) then
   !   rho   = 8.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)
   !else
   !rho   = 4.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)
   !!end if
   !xcell = 5.0_dp
   !ycell = 5.0_dp
   !PRINT*,'Size of binning cell', xcell, 'X', ycell
   fracFactor = 3.0*aG
   fracFactor2 = fracFactor**2.0
   cutoff2 = cutoff2/fracFactor2
   cutoff2bis = cutoff2bis/fracFactor2
   rho   = 4.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC/fracFactor*aCC/fracFactor)
   !end if
   xcell = 3.0_dp/fracFactor
   ycell = 3.0_dp/fracFactor
   !PRINT*,'Size of binning cell', xcell, 'X', ycell
   !============================================================
   !call MIO_Allocate(Species,nAt,'Species','atoms')
   !call MIO_Allocate(NeighD,[1,1,inode1],[3,maxNeigh,inode2],'NeighD','neigh')
   !call MIO_Allocate(xnew,[1],[natoms*9],'xnew','neigh')
   !call MIO_Allocate(ynew,[1],[natoms*9],'ynew','neigh')
   !call MIO_Allocate(znew,[1],[natoms*9],'znew','neigh')
   !allocate(xnew(9*natoms),ynew(9*natoms),znew(9*natoms))
   !print*, "numberOfAtoms: ", natoms, nAt
   do i=1,natoms
     xnew(i)=x(i); ynew(i)=y(i); znew(i)=z(i)
   end do
   do i=1,natoms
     n=i+natoms
     xnew(n)=x(i)+1.0
     ynew(n)=y(i)+0.0
     znew(n)=z(i)
   end do
   do i=1,natoms
     n=i+2*natoms
     xnew(n)=x(i)+1.0
     ynew(n)=y(i)+1.0
     znew(n)=z(i)
   end do
   do i=1,natoms
     n=i+3*natoms
     xnew(n)=x(i)+0.0
     ynew(n)=y(i)+1.0
     znew(n)=z(i)
   end do
   do i=1,natoms
     n=i+4*natoms
     xnew(n)=x(i)-1.0
     ynew(n)=y(i)+1.0
     znew(n)=z(i)
   end do
   do i=1,natoms
     n=i+5*natoms
     xnew(n)=x(i)-1.0
     ynew(n)=y(i)-0.0
     znew(n)=z(i)
   end do
   do i=1,natoms
     n=i+6*natoms
     xnew(n)=x(i)-1.0
     ynew(n)=y(i)-1.0
     znew(n)=z(i)
   end do
   do i=1,natoms
     n=i+7*natoms
     xnew(n)=x(i)-0.0
     ynew(n)=y(i)-1.0
     znew(n)=z(i)
   end do
   do i=1,natoms
     n=i+8*natoms
     xnew(n)=x(i)+1.0
     ynew(n)=y(i)-1.0
     znew(n)=z(i)
   end do
   !============================================================
   xmin = minval(xnew); xmax = maxval(xnew)
   ymin = minval(ynew); ymax = maxval(ynew)

   Nx = CEILING((xmax-xmin)/xcell)
   Ny = CEILING((ymax-ymin)/ycell)

   !write(*,*) xmin,xmax,Nx
   !write(*,*) ymin,ymax,Ny
   ALLOCATE(cells(Nx, Ny, 4*CEILING(xcell*ycell*rho))) ! 2 instead of 4 if not bilayer
   !write(*,*) 'Inicia llenado de cajas'
   !write(*,*) 'First dimension of cells is equal to ', Nx
   !write(*,*) 'Second dimension of cells is equal to ', Ny
   !write(*,*) 'Third dimension of cells is equal to ', 4*CEILING(xcell*ycell*rho)
    cells = 0
   !print*, cells
   do n = 1,9*natoms

     ix     = INT(1+(xnew(n)-xmin)/xcell)
     iy     = INT(1+(ynew(n)-ymin)/ycell)
     ixs(n) = ix
     iys(n) = iy

     i = 1
     !print*, ix, iy, i
     do while (cells(ix,iy,i) .ne. 0)
       !print*, i, cells(ix,iy,i)
       i = i+1
     end do
     cells(ix,iy,i) = n
     ! HERE
     !ncell(:,1) = [ix,iy,0]
     !neighCell(:,Nneigh(i),i) = ncell(:,1)
   end do
   !print*, "got out"



   maxNeigh = maxnn
   call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
   !call MIO_Allocate(NList,[1,inode1],[maxnn,inode2],'NList','neigh')
   call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
   call MIO_Allocate(NeighD,[1,1,inode1],[3,maxnn,inode2],'NeighD','neigh')
   !call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
   call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,maxnn,inode2/),'neighCell','neigh')
   !ALLOCATE(cells(Nx, Ny, 2*CEILING(xcell*ycell*rho)))
   !call MIO_Allocate(2*CEILING(xcell*ycell*rho),
   ! HERE

   write(*,*) 'Inicia busqueda'
   ! Find NNs of each atom by only searching nearby cells
   safetycounter = 0
   DO n = 1,natoms
     NNcount = 0

     DO addx = -2,2
     DO addy = -2,2
       ix = ixs(n) + addx
       iy = iys(n) + addy

       ! Periodic boundary conditions here
       IF(ix .gt. Nx) ix = ix-Nx
       IF(ix .lt. 1)  ix = ix+Nx
       IF(iy .gt. Ny) iy = iy-Ny
       IF(iy .lt. 1)  iy = iy+Ny

       i = 1
       DO WHILE(cells(ix, iy, i) .ne. 0)
         m = cells(ix, iy, i)
         IF (m .ne. n) THEN
           !write(*,*) 'hi2'
           dx_t = xnew(n) - xnew(m)
           dy_t = ynew(n) - ynew(m)
           dz_t = znew(n) - znew(m)
           ! Periodic boundary conditions here as well
           ! IF(dx_t .gt.  A1(1)-1.2_dp*aCC) dx_t = dx_t - A1(1)
           ! IF(dx_t .lt. -A1(1)+1.2_dp*aCC) dx_t = dx_t + A1(1)
           ! IF(dy_t .gt.  A2(2)-1.2_dp*aCC) dy_t = dy_t - A2(2)
           ! IF(dy_t .lt. -A2(2)+1.2_dp*aCC) dy_t = dy_t + A2(2)
           d2_t = dx_t**2.0_dp + dy_t**2.0_dp + dz_t**2.0_dp
           d2_t2 = dx_t**2.0_dp + dy_t**2.0_dp
             !IF (n==1) print*, "distances from n=1 ", d2_t
             !IF (n==3) print*, "distances from n=1 ", d2_t
           IF (abs(dz_t)<0.001_dp) THEN
             !print*, cutoff2
             IF (d2_t .lt. cutoff2*1.1_dp) THEN
               NNcount = NNcount + 1
               call NeighCheckCount(NNcount,maxnn)
               vecino=mod(m,natoms)
               if (vecino==0) vecino=natoms
               nn(n,NNcount) = vecino
               !NList(Nneigh(i),i) = j
               dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
               !dr(n,NNcount)=dsqrt(d2_t)
               NList(NNcount,n) = vecino
               NeighD(1,NNcount,n) = dx(n,NNcount)
               NeighD(2,NNcount,n) = dy(n,NNcount)
               NeighD(3,NNcount,n) = dz(n,NNcount)
               Nneigh(n) = NNcount
               neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
               neighCell(2,NNcount,n) = iys(m) - iys(n)
               neighCell(3,NNcount,n) = 0
   !             if(n==4800) write(*,*) m,NNcount,vecino
             ENDIf
           ELSE
             !print*, "cutoffffff", cutoff2bis
             !print*, "dt_t2", d2_t2
             IF (d2_t2 .lt. cutoff2bis) THEN
               NNcount = NNcount + 1
               call NeighCheckCount(NNcount,maxnn)
               vecino=mod(m,natoms)
               if (vecino==0) vecino=natoms
               nn(n,NNcount) = vecino
               !NList(Nneigh(i),i) = j
               dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
               !dr(n,NNcount)=dsqrt(d2_t)
               NList(NNcount,n) = vecino
               NeighD(1,NNcount,n) = dx(n,NNcount)
               NeighD(2,NNcount,n) = dy(n,NNcount)
               NeighD(3,NNcount,n) = dz(n,NNcount)
               Nneigh(n) = NNcount
               neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
               neighCell(2,NNcount,n) = iys(m) - iys(n)
               neighCell(3,NNcount,n) = 0
   !             if(n==4800) write(*,*) m,NNcount,vecino
             ENDIf
           ENDIF
         ENDIF
         i = i+1
       ENDDO

     ENDDO
     ENDDO

     IF(NNcount .lt. 3) THEN
       !print*, "Are you sure you wanted to enter this safety routine?"
       NNcount = 0
       DO m = 1,natoms
         IF (n.ne.m) THEN
           dx_t = x(n) - x(m)
           dy_t = y(n) - y(m)
           dz_t = z(n) - z(m)
            !Periodic boundary conditions here as well
            IF(dx_t .gt.  A1(1)-1.2_dp*aCC) dx_t = dx_t - A1(1)
            IF(dx_t .lt. -A1(1)+1.2_dp*aCC) dx_t = dx_t + A1(1)
            IF(dy_t .gt.  A2(2)-1.2_dp*aCC) dy_t = dy_t - A2(2)
            IF(dy_t .lt. -A2(2)+1.2_dp*aCC) dy_t = dy_t + A2(2)
            d2_t = dx_t**2.0_dp + dy_t**2.0_dp + dz_t**2.0_dp
            IF (abs(dz_t)<0.01_dp) THEN
               IF (d2_t .lt. cutoff2*1.1_dp) THEN
                    NNcount = NNcount + 1
                    call NeighCheckCount(NNcount,maxnn)
                    !vecino=mod(m,natoms)
                    !if (vecino==0) vecino=natoms
                    !nn(n,NNcount) = vecino
                    !NList(Nneigh(i),i) = j
                    dx(n,NNcount)=dx_t; dy(n,NNcount)=dy_t; dz(n,NNcount)=dz_t
                    !dr(n,NNcount)=dsqrt(d2_t)
                    NList(NNcount,n) = m
                    NeighD(1,NNcount,n) = dx(n,NNcount)
                    NeighD(2,NNcount,n) = dy(n,NNcount)
                    NeighD(3,NNcount,n) = dz(n,NNcount)
                    Nneigh(n) = NNcount


                    !ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
                    !neighCell(:,Nneigh(i),i) = ncell(:,1)


                    !neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
                    !neighCell(2,NNcount,n) = iys(m) - iys(n)
                    !neighCell(3,NNcount,n) = 0
   !                if(n==4800) write(*,*) m,NNcount,vecino
               ENDIf
            END IF
          END IF
       END DO
       safetycounter = safetycounter + 1
       if (safetycounter.gt.100) print*, "Carefull, you are probably going too many times through this safety loop. &
      I originally implemented this because the fastNNnotsquared routine was missing neighbors for 4 atoms only..."
     ENDIF


   ! Count the number of atoms with 1, 2, 3, 4, 5 nearest neighbors
     IF(NNcount .eq. 0) num0 = num0 + 1
     !IF(NNcount .eq. 0) print*, "NN counter = 0", n
     IF(NNcount .eq. 1) num1 = num1 + 1
     IF(NNcount .eq. 2) num2 = num2 + 1
     !IF(NNcount .eq. 2) print*, "NN counter = 2", n
     IF(NNcount .eq. 3) num3 = num3 + 1
     IF(NNcount .eq. 4) num4 = num4 + 1
     IF(NNcount .eq. 5) num5 = num5 + 1
     IF(NNcount .eq. 6) num6 = num6 + 1
     IF(NNcount .eq. 7) num7 = num7 + 1
     IF(NNcount .eq. 8) num8 = num8 + 1
     IF(NNcount .eq. 9) num9 = num9 + 1
     IF(NNcount .eq. 10) num10 = num10 + 1



   ENDDO

   if (frac) call AtomsSetCart()

   call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')

   do i=1,tbnn
      Nradii(i,1) = rad(i)*aG*1.1_dp
   end do
   call MIO_InputParameter('TypeOfSystem',str,'Graphene')
   if (MIO_StringComp(str,'BoronNitride')) then
      Nradii(:,2) = Nradii(:,1)
   else
      do i=1,tbnn
         Nradii(i,2) = rad(i)*aBN*1.1_dp
      end do
   end if
   rmax = maxval(Nradii)

   !$OMP PARALLEL PRIVATE (v, d,ncell,i,j,ix,iy)
   do i=inode1,inode2
      do j=1,Nneigh(i)
         do ix=-1,1; do iy=-1,1
            if (i==NList(j,i) .and. ix==0 .and. iy==0) cycle
            ncell(:,1) = [ix,iy,0] ! We only use one of the nine columns if small, the other columns are used for the other case
            v = Rat(:,i) - Rat(:,NList(j,i)) - matmul(ucell,ncell(:,1)) ! If they are from different unit cells, this extra term will make v very big, and it will not satisfy the distance condition.
                                                               ! Only if it is same unit cell (ix = 0, iy = 0) or neighor cell (e.g. 0 and 1), it might be satisfied (unless the bins are very small)
            d = norm(v)
            if (d<rmax) then
               !Nneigh(i) = Nneigh(i) + 1
               !if (Nneigh(i)>maxNeigh)  then
               !   !print*, "Number of neighbors: ", i, Nneigh(i)
               !   call MIO_Kill('More neighbors than expected found',&
               !     'neigh','NeighList')
               !end if
               !NList(Nneigh(i),i) = j
               !NeighD(:,Nneigh(i),i) = -v
               neighCell(:,j,i) = ncell(:,1)
            !else
            !    print*, "Are you sure this was a neighbor?", i, NList(j,i)
            end if
         end do; end do
      end do
   end do
   !$OMP END PARALLEL

   !write(*,*) 'hi4'
   !neighCell = cells
   DEALLOCATE(cells)
   !write(*,*) 'hi6'
   !call MIO_Deallocate(xnew,'xnew','neigh')
   !call MIO_Deallocate(ynew,'ynew','neigh')
   !call MIO_Deallocate(znew,'znew','neigh')
   !deallocate(xnew,ynew,znew)
   !write(*,*) 'hi5'


   PRINT*,' 0 neighbors: ',num0
   PRINT*,' 1 neighbors: ',num1
   PRINT*,' 2 neighbors: ',num2
   PRINT*,' 3 neighbors: ',num3
   PRINT*,' 4 neighbors: ',num4
   PRINT*,' 5 neighbors: ',num5
   PRINT*,' 6 neighbors: ',num6
   PRINT*,' 7 neighbors: ',num7
   PRINT*,' 8 neighbors: ',num8
   PRINT*,' 9 neighbors: ',num9
   PRINT*,' 10 neighbors: ',num10
   PRINT*,'>10 neighbors: ',natoms-(num0+num1+num2+num3+num4+num5+num6+num7+num8+num9+num10)

   PRINT*,''
   PRINT*,'...done'
   PRINT*,''

   call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
   NList2 = NList
#ifdef MPI
   if (nProc>1) then
      call MIO_Allocate(countN,[nProc,numThreads],'countN','neigh')
      countN = 0
      !$OMP PARALLEL
      do i=in1,in2
         do j=1,maxNeigh
            if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
               do k=nProc,1,-1
                  if(NList2(j,i)>=indxNode(k)) then
                     countN(k,nThread+1) = countN(k,nThread+1)+1
                     exit
                  end if
               end do
            end if
         end do
      end do
      !$OMP END PARALLEL
      rcvSz = sum(countN)
      call MIO_Allocate(rcvIndx,[2,nProc],'rcvIndx','neigh')
      call MIO_Allocate(rcvList,rcvSz,'rcvList','neigh')
      np = 1
      do i=1,nProc
         rcvIndx(1,i) = np
         np = np + sum(countN(i,:))
         rcvIndx(2,i) = np - 1
      end do
      call MIO_Allocate(cnt,[nProc,numThreads],'cnt','neigh')
      !$OMP PARALLEL PRIVATE(np)
      do i=in1,in2
         do j=1,maxNeigh
            if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
               do k=nProc,1,-1
                  if(NList2(j,i)>=indxNode(k)) then
                     np = rcvIndx(1,k) + sum(countN(k,1:nThread)) + cnt(k,nThread+1)
                     cnt(k,nThread+1) = cnt(k,nThread+1) + 1
                     exit
                  end if
               end do
               rcvList(np) = NList2(j,i)
               NList2(j,i) = np + inode2
            end if
         end do
      end do
      !$OMP END PARALLEL
      call MIO_Allocate(nodeSndRcv,[nProc,nProc],'nodeSndRcv','neigh')
      do i=1,nProc
         nodeSndRcv(i,Node+1) = rcvIndx(2,i)-rcvIndx(1,i) + 1
      end do
      call MPIAllGather(MPI_IN_PLACE,0,nodeSndRcv,nProc,MPI_INTEGER)
      sndSz = sum(nodeSndRcv(Node+1,:))
      call MIO_Allocate(sndList,sndSz,'sndList','neigh')
      call MIO_Allocate(sndIndx,(/2,nProc/),'sndIndx','neigh')
      np = 1
      do i=1,nProc
         sndIndx(1,i) = np
         np = np + nodeSndRcv(Node+1,i)
         sndIndx(2,i) = np - 1
      end do
      do i=1,nProc-1
         j = mod(Node + i,nProc)
         k = mod(nProc + Node - i,nProc)
         call MPISendRecv(rcvList(rcvIndx(1,j+1):rcvIndx(2,j+1)),nodeSndRcv(j+1,Node+1),j,50, &
           sndList(sndIndx(1,k+1):sndIndx(2,k+1)),nodeSndRcv(Node+1,k+1),k,50,MPI_INTEGER)
         call MPISendRecv(rcvList(rcvIndx(1,k+1):rcvIndx(2,k+1)),nodeSndRcv(k+1,Node+1),k,51, &
           sndList(sndIndx(1,j+1):sndIndx(2,j+1)),nodeSndRcv(Node+1,j+1),j,51,MPI_INTEGER)
         call MPIBarrier()
      end do
      call MIO_Deallocate(cnt,'cnt','neigh')
      call MIO_Deallocate(countN,'countN','neigh')
   end if
#endif /* MPI */
   call MIO_InputParameter('WriteDataFiles',prnt,.false.)
   !call AtomsSetFrac()
   if (prnt) then
      open(1,FILE='v')
      open(2,FILE='dx')
      open(3,FILE='dy')
      open(33,FILE='dz')
      open(4,FILE='pos')
      do i=1,nAt
         write(1,*) Nneigh(i)
         !write(1,'(300(I8,1X))') (NList(j,i),j=1,Nneigh(i))
         !write(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
         !write(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
         !write(33,*) (NeighD(3,j,i),j=1,Nneigh(i))
         write(1,'(500(I8,1X))') (NList(j,i),j=1,Nneigh(i))
         write(2,'(500(F10.5,1X))') (NeighD(1,j,i),j=1,Nneigh(i))
         write(3,'(500(F10.5,1X))') (NeighD(2,j,i),j=1,Nneigh(i))
         write(33,'(500(F10.5,1X))') (NeighD(3,j,i),j=1,Nneigh(i))
         if (Species(i)==3) then
            write(4,*) "B    ", Rat(1,i),Rat(2,i),Rat(3,i)
         else if (Species(i)==4) then
            write(4,*) "N    ", Rat(1,i),Rat(2,i),Rat(3,i)
         else
            write(4,*) "C    ", Rat(1,i),Rat(2,i),Rat(3,i)
         end if
      end do
      close(1)
      close(2)
      close(3)
      close(33)
      close(4)
   end if



end subroutine fastNNnotsquareNotRectangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine: fastNN                                                                        !
! Description: calculate nearest neighbors of a sample with a fast method that uses binning !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE fastNN(natoms, x,y,z, aCC,cutoff, A1,A2, maxnn)

   use atoms,                only : inode1, inode2, in1, in2, indxNode, nAt, Species, Rat
   use tbpar,                only : tbnn
   use cell,                 only : aG, aBN
   use atoms,                only : AtomsSetFrac

   IMPLICIT NONE

   ! Subroutine inputs
   INTEGER              :: natoms
   DOUBLE PRECISION     :: x(natoms),y(natoms),z(natoms), aCC,cutoff
   real(dp), intent(in) :: A1(3),A2(3)

   ! Subroutine outputs
   INTEGER              :: maxnn, nn(natoms,maxnn), near(natoms)

   ! For subdividing sample into cells
   INTEGER              :: Nx,Ny, ix,iy, ixs(natoms),iys(natoms)
   INTEGER, ALLOCATABLE :: cells(:,:,:)
   DOUBLE PRECISION     :: rho,xcell,ycell, xmin,xmax,ymin,ymax, cutoff2

   ! For finding and counting nearest neighbors
   INTEGER              :: num0,num1,num2,num3,num4,num5
   INTEGER              :: NNcount, addx,addy
   DOUBLE PRECISION     :: dx,dy,dz,dist2

   ! Loop index variables
   INTEGER              :: i,m,n,j,k

   integer, pointer :: countN(:,:), cnt(:,:)
   integer :: np

   real(dp), parameter :: rad(3) = [1.0_dp/sqrt(3.0_dp),1.0_dp,2.0_dp/sqrt(3.0_dp)]
   character(len=50) :: str
   logical :: prnt


   real(dp) :: distFact

   call MIO_Allocate(Nradii,[tbnn+1,2],'Nradii','neigh')

   do i=1,tbnn
      Nradii(i,1) = rad(i)*aG*1.1_dp
   end do
   call MIO_InputParameter('TypeOfSystem',str,'Graphene')
   if (MIO_StringComp(str,'BoronNitride')) then
      Nradii(:,2) = Nradii(:,1)
   else
      do i=1,tbnn
         Nradii(i,2) = rad(i)*aBN*1.1_dp
      end do
   end if



   maxNeigh = maxnn
   call MIO_Allocate(NList,[1,inode1],[maxNeigh,inode2],'NList','neigh')
   !call MIO_Allocate(NList,[1,inode1],[maxnn,inode2],'NList','neigh')
   call MIO_Allocate(Nneigh,[inode1],[inode2],'Nneigh','neigh')
   call MIO_Allocate(NeighD,[1,1,inode1],[3,maxnn,inode2],'NeighD','neigh')
   !call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,3,inode2/),'neighCell','neigh')
   call MIO_Allocate(neighCell,(/1,1,inode1/),(/3,maxnn,inode2/),'neighCell','neigh')

   ! Initialize maxnn, set cutoff to cutoff² to avoid using sqrt()
   num0    = 0
   num1    = 0
   num2    = 0
   num3    = 0
   num4    = 0
   num5    = 0
   maxnn   = 0
   cutoff2 = cutoff*cutoff*aCC*aCC


   ! Divide geometry into rectangular cells of size xcell times ycell (units of Angstroms)
   if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer') &
           .or. MIO_StringComp(str,'TwistedBilayerBasedOnMoireCell') .or. MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell')) then
      rho   = 4.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)
   else
      rho   = 4.0_dp / (3.0_dp*DSQRT(3.0_dp)*aCC*aCC)
   end if
   xcell = 3.0_dp
   ycell = 3.0_dp

   xmin = MINVAL(x)
   xmax = MAXVAL(x)
   ymin = MINVAL(y)
   ymax = MAXVAL(y)


   Nx = CEILING((xmax-xmin)/xcell)
   Ny = CEILING((ymax-ymin)/ycell)
   ALLOCATE(cells(Nx, Ny, 4*CEILING(xcell*ycell*rho)))

   cells = 0
   DO n = 1,natoms

     ix     = INT(1+(x(n)-xmin)/xcell)
     iy     = INT(1+(y(n)-ymin)/ycell)
     ixs(n) = ix
     iys(n) = iy

     i = 1
     DO WHILE (cells(ix,iy,i) .ne. 0)
       i = i+1
     ENDDO
     cells(ix,iy,i) = n

   ENDDO

   call MIO_InputParameter('Neigh.LayerDistFactor',distFact,1.0_dp) ! to increase the cutoff for outerlayer neighbors


   ! Find nearest neighbors by searching atoms in current and adjacent cells
   DO n = 1,natoms
     NNcount = 0

     DO addx = -2,2
     DO addy = -2,2
       ix = ixs(n) + addx
       iy = iys(n) + addy

       ! Periodic boundary conditions here
       IF(ix .gt. Nx) ix = ix-Nx
       IF(ix .lt. 1)  ix = ix+Nx
       IF(iy .gt. Ny) iy = iy-Ny
       IF(iy .lt. 1)  iy = iy+Ny

       i = 1
       DO WHILE(cells(ix, iy, i) .ne. 0)
         m = cells(ix, iy, i)
         IF (m .ne. n) THEN
           dx = x(n) - x(m)
           dy = y(n) - y(m)
           dz = z(n) - z(m)

           ! Periodic boundary conditions here as well
           IF(dx .gt.  A1(1)/2.0_dp) dx = dx - A1(1) - A2(2)
           IF(dx .lt. -A1(1)/2.0_dp) dx = dx + A1(1) + A2(1)
           IF(dy .gt.  A2(2)/2.0_dp) dy = dy - A2(2) - A1(2)
           IF(dy .lt. -A2(2)/2.0_dp) dy = dy + A2(2) + A1(2)

           dist2 = dx**2.0_dp + dy**2.0_dp + dz**2.0_dp
           !cutoff2 = aG/sqrt(3.0_dp)*1.1_dp*distFact
           IF (dist2 .lt. cutoff2*2.5) THEN
             NNcount = NNcount + 1
             call NeighCheckCount(NNcount,maxnn)
             nn(n,NNcount) = m
             NList(NNcount,n) = m
             NeighD(1,NNcount,n) = dx
             NeighD(2,NNcount,n) = dy
             NeighD(3,NNcount,n) = dz
             Nneigh(n) = NNcount
             neighCell(1,NNcount,n) = ixs(m) - ixs(n) ! check if it's not n-m
             neighCell(2,NNcount,n) = iys(m) - iys(n)
             neighCell(3,NNcount,n) = 0
           ENDIf
         ENDIF
         i = i+1
       ENDDO

     ENDDO
     ENDDO

     near(n) = NNcount
     IF(NNcount .gt. maxnn) maxnn = NNcount


   ! Count the number of atoms with 1, 2, 3, 4, 5 nearest neighbors
     IF(NNcount .eq. 0) num0 = num0 + 1
     IF(NNcount .eq. 1) num1 = num1 + 1
     IF(NNcount .eq. 2) num2 = num2 + 1
     IF(NNcount .eq. 3) num3 = num3 + 1
     IF(NNcount .eq. 4) num4 = num4 + 1
     IF(NNcount .eq. 5) num5 = num5 + 1



   ENDDO
   DEALLOCATE(cells)


   PRINT*,'  0 neighbors: ',num0
   PRINT*,'  1 neighbors: ',num1
   PRINT*,'  2 neighbors: ',num2
   PRINT*,'  3 neighbors: ',num3
   PRINT*,'  4 neighbors: ',num4
   PRINT*,'  5 neighbors: ',num5
   PRINT*,' >5 neighbors: ',natoms-(num0+num1+num2+num3+num4+num5)


   call MIO_Allocate(NList2,[1,inode1],[maxNeigh,inode2],'NList','neigh')
   NList2 = NList
#ifdef MPI
   if (nProc>1) then
      call MIO_Allocate(countN,[nProc,numThreads],'countN','neigh')
      countN = 0
      !$OMP PARALLEL
      do i=in1,in2
         do j=1,maxNeigh
            if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
               do k=nProc,1,-1
                  if(NList2(j,i)>=indxNode(k)) then
                     countN(k,nThread+1) = countN(k,nThread+1)+1
                     exit
                  end if
               end do
            end if
         end do
      end do
      !$OMP END PARALLEL
      rcvSz = sum(countN)
      call MIO_Allocate(rcvIndx,[2,nProc],'rcvIndx','neigh')
      call MIO_Allocate(rcvList,rcvSz,'rcvList','neigh')
      np = 1
      do i=1,nProc
         rcvIndx(1,i) = np
         np = np + sum(countN(i,:))
         rcvIndx(2,i) = np - 1
      end do
      call MIO_Allocate(cnt,[nProc,numThreads],'cnt','neigh')
      !$OMP PARALLEL PRIVATE(np)
      do i=in1,in2
         do j=1,maxNeigh
            if (NList2(j,i)<inode1 .or. NList2(j,i)>inode2) then
               do k=nProc,1,-1
                  if(NList2(j,i)>=indxNode(k)) then
                     np = rcvIndx(1,k) + sum(countN(k,1:nThread)) + cnt(k,nThread+1)
                     cnt(k,nThread+1) = cnt(k,nThread+1) + 1
                     exit
                  end if
               end do
               rcvList(np) = NList2(j,i)
               NList2(j,i) = np + inode2
            end if
         end do
      end do
      !$OMP END PARALLEL
      call MIO_Allocate(nodeSndRcv,[nProc,nProc],'nodeSndRcv','neigh')
      do i=1,nProc
         nodeSndRcv(i,Node+1) = rcvIndx(2,i)-rcvIndx(1,i) + 1
      end do
      call MPIAllGather(MPI_IN_PLACE,0,nodeSndRcv,nProc,MPI_INTEGER)
      sndSz = sum(nodeSndRcv(Node+1,:))
      call MIO_Allocate(sndList,sndSz,'sndList','neigh')
      call MIO_Allocate(sndIndx,(/2,nProc/),'sndIndx','neigh')
      np = 1
      do i=1,nProc
         sndIndx(1,i) = np
         np = np + nodeSndRcv(Node+1,i)
         sndIndx(2,i) = np - 1
      end do
      do i=1,nProc-1
         j = mod(Node + i,nProc)
         k = mod(nProc + Node - i,nProc)
         call MPISendRecv(rcvList(rcvIndx(1,j+1):rcvIndx(2,j+1)),nodeSndRcv(j+1,Node+1),j,50, &
           sndList(sndIndx(1,k+1):sndIndx(2,k+1)),nodeSndRcv(Node+1,k+1),k,50,MPI_INTEGER)
         call MPISendRecv(rcvList(rcvIndx(1,k+1):rcvIndx(2,k+1)),nodeSndRcv(k+1,Node+1),k,51, &
           sndList(sndIndx(1,j+1):sndIndx(2,j+1)),nodeSndRcv(Node+1,j+1),j,51,MPI_INTEGER)
         call MPIBarrier()
      end do
      call MIO_Deallocate(cnt,'cnt','neigh')
      call MIO_Deallocate(countN,'countN','neigh')
   end if
#endif /* MPI */
   call MIO_InputParameter('WriteDataFiles',prnt,.false.)
   !call AtomsSetFrac()
   if (prnt) then
      open(1,FILE='v')
      open(2,FILE='dx')
      open(3,FILE='dy')
      open(33,FILE='dz')
      open(4,FILE='pos')
      do i=1,nAt
         write(1,*) Nneigh(i)
         !write(1,'(300(I8,1X))') (NList(j,i),j=1,Nneigh(i))
         !write(2,*) (NeighD(1,j,i),j=1,Nneigh(i))
         !write(3,*) (NeighD(2,j,i),j=1,Nneigh(i))
         !write(33,*) (NeighD(3,j,i),j=1,Nneigh(i))
         write(1,'(500(I8,1X))') (NList(j,i),j=1,Nneigh(i))
         write(2,'(500(F10.5,1X))') (NeighD(1,j,i),j=1,Nneigh(i))
         write(3,'(500(F10.5,1X))') (NeighD(2,j,i),j=1,Nneigh(i))
         write(33,'(500(F10.5,1X))') (NeighD(3,j,i),j=1,Nneigh(i))
         if (Species(i)==3) then
            write(4,*) "B    ", Rat(1,i),Rat(2,i),Rat(3,i)
         else if (Species(i)==4) then
            write(4,*) "N    ", Rat(1,i),Rat(2,i),Rat(3,i)
         else
            write(4,*) "C    ", Rat(1,i),Rat(2,i),Rat(3,i)
         end if
      end do
      close(1)
      close(2)
      close(3)
      close(33)
      close(4)
   end if


   RETURN
END SUBROUTINE fastNN


end module neigh
