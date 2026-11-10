module cell

   use mio

   implicit none

   PRIVATE

   real(dp), public, save :: ucell(3,3), rcell(3,3), volume, area, aG, aBN
   integer, public, save :: sCell, sCell2

   public :: CellGet

contains

subroutine CellGet()

   use math,                 only : CrossProd, norm
   use constants,            only : twopi, pi
   use name,                 only : prefix


   ! Switches and options read from the input
   logical :: l, ll
   character(len=80) :: str

   ! Model parameters read from the input
   integer :: m(4), n, nBN, nC
   real(dp) :: h

   ! Other local variables
   integer :: i, j, m1(2,2)
   real(dp) :: delta, g, phi, uc(2,2), vn(3)
   character(len=120) :: dbgline


#ifdef DEBUG
   call MIO_Debug('CellGet',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('cell')
#endif /* TIMER */

   call MIO_InputParameter('Structure.LatticeParameter',aG,2.46_dp)
     !! - `Structure.LatticeParameter`: Graphene lattice constant aG [Angstrom]; the cells built by the code are multiples of it.
   call MIO_InputParameter('Structure.LatticeParameterBN',aBN,2.505_dp)
     !! - `Structure.LatticeParameterBN`: hBN lattice constant [Angstrom] (cells built from Structure.MoireCellParameters set it
     !!   from the indices).
   call MIO_InputParameter('Run.TypeOfSystem',str,'Graphene')
     !! - `Run.TypeOfSystem`: Kind of structure: built by the code (Graphene, BoronNitride, Ribbons,
     !!   TwistedBilayerBasedOnMoireCell, ...) or read from a file (ReadXYZ). It selects the cell, the atoms and parts of the
     !!   Hamiltonian.
   call MIO_Section('Cell','cell')
   call MIO_Print('Type of system: '//trim(str),'cell')
   call MIO_InputParameter('Structure.CellHeight',h,40.0_dp)
     !! - `Structure.CellHeight`: Length of the third (out-of-plane) cell vector [Angstrom] for the cells built by the code.
   call MIO_InputParameter('Structure.SuperCell',sCell,1)
     !! - `Structure.SuperCell`: Repeat the cell N x N in the plane (the atoms are copied, no reordering: copy 0 holds the original
     !!   atoms). Equivalent to an N x N k-sampling of the original cell for local quantities.
   ! Default of the interlayer distance: 3.35 Ang everywhere since October 2026 (3.22 Ang before in the
   ! generators, the neighbor search and parts of the Hamiltonian; it matters for structures read from a file too)
   if (.not. (MIO_StringComp(str,'Graphene') .or. MIO_StringComp(str,'BoronNitride'))) then
      if (.not. MIO_InputSearchLabel('Structure.InterlayerDistance')) then
         call MIO_Print('Structure.InterlayerDistance is not given: 3.35 Ang is used (the default was '// &
           '3.22 Ang before October 2026)','cell')
      end if
   end if
   call MIO_InputParameter('Structure.SuperCellAsymmetric',l,.false.)
     !! - `Structure.SuperCellAsymmetric`: Repeat the cell Structure.SuperCell times along a1 and Structure.SuperCellY times along
     !!   a2.
   call MIO_InputParameter('Structure.SuperCellY',sCell2,1)
     !! - `Structure.SuperCellY`: Structure.SuperCellAsymmetric: number of copies along the second lattice vector.
   if (MIO_StringComp(str,'Graphene') .or. MIO_StringComp(str,'BoronNitride') &
    .or. MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer')) then
      call MIO_InputParameter('Structure.BasedOnMoireCellParameters',ll,.false.)
        !! - `Structure.BasedOnMoireCellParameters`: Graphene / BoronNitride / TwistedBilayer cells: build the commensurate cell
        !!   from the four indices of Structure.MoireCellParameters instead of CellSize x CellSize graphene cells.
      if (ll) then
          uc(:,1) = (/aG,0.0_dp/)
          uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
          call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
            !! - `Structure.MoireCellParameters`: Four integers m1 m2 m3 m4 of the commensurate moire cell: a1 = m1 a + m2 b on the
            !!   first lattice, (m3, m4) on the second; the twist angle and lattice ratio follow from them and are printed.
          call MIO_Print('Moire cell indices (m1 m2 m3 m4): '//trim(num2str(m(1)))//&
            '  '//trim(num2str(m(2)))//'  '//trim(num2str(m(3)))//'  '//trim(num2str(m(4))),&
            'cell')
          g = m(1)**2 + m(2)**2 + m(1)*m(2)
          delta = sqrt(real(m(3)**2 + m(4)**2 + m(3)*m(4))/g)
          phi = acos((2.0_dp*m(1)*m(3)+2.0_dp*m(2)*m(4) + m(1)*m(4) + m(2)*m(3))/(2.0_dp*delta*g))
          if (delta<1.0_dp) then
             delta = 1.0_dp/delta
          else
             i = m(1)
             m(1) = m(3)
             m(3) = i
             i = m(2)
             m(2) = m(4)
             m(4) = i
          end if
          aBN = aG*delta
          call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' deg','cell')
          call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4))//' Ang','cell')
          call MIO_Print('BN lattice parameter: '//trim(num2str(aBN,4))//' Ang','cell')
          m1(:,1) = [m(1),m(2)]
          m1(:,2) = [-m(2),m(1)+m(2)]
          ucell = 0.0_dp
          ucell(1:2,1:2) = matmul(uc,m1)
          ucell(3,3) = h
          g = norm(ucell(:,1))
          call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
            ' Ang','cell')
          ucell(:,:2) = ucell(:,:2)*sCell
          aBN = aG ! added this, not sure if necessary
      else
          call MIO_InputParameter('Structure.CellSize',n,50)
            !! - `Structure.CellSize`: Number of graphene cells along each lattice vector of a cell built by the code (TypeOfSystem
            !!   Graphene, BoronNitride); the cell is CellSize*SuperCell graphene cells on a side.
          call MIO_Print('Length of unit cell: '//trim(num2str(n*aG,5))// &
            ' Ang','cell')
             n = n*sCell
          ucell(:,1) = [aG,0.0_dp,0.0_dp]
          ucell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
            ucell = ucell*n
          ! SuperCellAsymmetric: the atoms are built on CellSize*SuperCell x CellSize*SuperCellY (AtomsConstructAsymm)
          if (l) ucell(:,2) = ucell(:,2)/sCell*sCell2
          ucell(:,3) = [0.0_dp,0.0_dp,h]
          aBN = aG
      end if
   else if (MIO_StringComp(str,'Graphene_Over_BN')) then
      uc(:,1) = (/aG,0.0_dp/)
      uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
      call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
      call MIO_Print('Moire cell indices (m1 m2 m3 m4): '//trim(num2str(m(1)))//&
        '  '//trim(num2str(m(2)))//'  '//trim(num2str(m(3)))//'  '//trim(num2str(m(4))),&
        'cell')
      g = m(1)**2 + m(2)**2 + m(1)*m(2)
      delta = sqrt(real(m(3)**2 + m(4)**2 + m(3)*m(4))/g)
      phi = acos((2.0_dp*m(1)*m(3)+2.0_dp*m(2)*m(4) + m(1)*m(4) + m(2)*m(3))/(2.0_dp*delta*g))
      if (delta<1.0_dp) then
         delta = 1.0_dp/delta
      else
         i = m(1)
         m(1) = m(3)
         m(3) = i
         i = m(2)
         m(2) = m(4)
         m(4) = i
      end if
      aBN = aG*delta
      call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' deg','cell')
      call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4))//' Ang','cell')
      call MIO_Print('BN lattice parameter: '//trim(num2str(aBN,4))//' Ang','cell')
      m1(:,1) = [m(1),m(2)]
      m1(:,2) = [-m(2),m(1)+m(2)]
      ucell = 0.0_dp
      ucell(1:2,1:2) = matmul(uc,m1)
      ucell(3,3) = h
      g = norm(ucell(:,1))
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
      ucell(:,:2) = ucell(:,:2)*sCell
   else if (MIO_StringComp(str,'BilayerGraphene_Over_BN')) then
      uc(:,1) = (/aG,0.0_dp/)
      uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
      call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
      call MIO_Print('Moire cell indices (m1 m2 m3 m4): '//trim(num2str(m(1)))//&
        '  '//trim(num2str(m(2)))//'  '//trim(num2str(m(3)))//'  '//trim(num2str(m(4))),&
        'cell')
      g = m(1)**2 + m(2)**2 + m(1)*m(2)
      delta = sqrt(real(m(3)**2 + m(4)**2 + m(3)*m(4))/g)
      phi = acos((2.0_dp*m(1)*m(3)+2.0_dp*m(2)*m(4) + m(1)*m(4) + m(2)*m(3))/(2.0_dp*delta*g))
      if (delta<1.0_dp) then
         delta = 1.0_dp/delta
      else
         i = m(1)
         m(1) = m(3)
         m(3) = i
         i = m(2)
         m(2) = m(4)
         m(4) = i
      end if
      aBN = aG*delta
      call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' deg','cell')
      call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4))//' Ang','cell')
      call MIO_Print('BN lattice parameter: '//trim(num2str(aBN,4))//' Ang','cell')
      m1(:,1) = [m(1),m(2)]
      m1(:,2) = [-m(2),m(1)+m(2)]
      ucell = 0.0_dp
      ucell(1:2,1:2) = matmul(uc,m1)
      ucell(3,3) = h
      g = norm(ucell(:,1))
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
      ucell(:,:2) = ucell(:,:2)*sCell
   else if (MIO_StringComp(str,'TwistedBilayerBasedOnMoireCell')) then
      uc(:,1) = (/aG,0.0_dp/)
      uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
      call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
      call MIO_Print('Moire cell indices (m1 m2 m3 m4): '//trim(num2str(m(1)))//&
        '  '//trim(num2str(m(2)))//'  '//trim(num2str(m(3)))//'  '//trim(num2str(m(4))),&
        'cell')
      g = m(1)**2 + m(2)**2 + m(1)*m(2)
      delta = sqrt(real(m(3)**2 + m(4)**2 + m(3)*m(4))/g)
      phi = acos((2.0_dp*m(1)*m(3)+2.0_dp*m(2)*m(4) + m(1)*m(4) + m(2)*m(3))/(2.0_dp*delta*g))
      if (delta<1.0_dp) then
         delta = 1.0_dp/delta
      else
         i = m(1)
         m(1) = m(3)
         m(3) = i
         i = m(2)
         m(2) = m(4)
         m(4) = i
      end if
      aBN = aG*delta
      call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' deg','cell')
      call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4))//' Ang','cell')
      call MIO_Print('BN lattice parameter: '//trim(num2str(aBN,4))//' Ang','cell')
      m1(:,1) = [m(1),m(2)]
      m1(:,2) = [-m(2),m(1)+m(2)]
      ucell = 0.0_dp
      ucell(1:2,1:2) = matmul(uc,m1)
      ucell(3,3) = h
      g = norm(ucell(:,1))
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
      ucell(:,:2) = ucell(:,:2)*sCell
   else if (MIO_StringComp(str,'TrilayerBasedOnMoireCell')) then
      uc(:,1) = (/aG,0.0_dp/)
      uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
      call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
      call MIO_Print('Moire cell indices (m1 m2 m3 m4): '//trim(num2str(m(1)))//&
        '  '//trim(num2str(m(2)))//'  '//trim(num2str(m(3)))//'  '//trim(num2str(m(4))),&
        'cell')
      g = m(1)**2 + m(2)**2 + m(1)*m(2)
      delta = sqrt(real(m(3)**2 + m(4)**2 + m(3)*m(4))/g)
      phi = acos((2.0_dp*m(1)*m(3)+2.0_dp*m(2)*m(4) + m(1)*m(4) + m(2)*m(3))/(2.0_dp*delta*g))
      if (delta<1.0_dp) then
         delta = 1.0_dp/delta
      else
         i = m(1)
         m(1) = m(3)
         m(3) = i
         i = m(2)
         m(2) = m(4)
         m(4) = i
      end if
      aBN = aG*delta
      call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' deg','cell')
      call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4))//' Ang','cell')
      call MIO_Print('BN lattice parameter: '//trim(num2str(aBN,4))//' Ang','cell')
      m1(:,1) = [m(1),m(2)]
      m1(:,2) = [-m(2),m(1)+m(2)]
      ucell = 0.0_dp
      ucell(1:2,1:2) = matmul(uc,m1)
      ucell(3,3) = h
      g = norm(ucell(:,1))
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
      ucell(:,:2) = ucell(:,:2)*sCell
   else if (MIO_StringComp(str,'TwistedBilayerBasedOnMoireCellRectangular')) then
      uc(:,1) = (/aG,0.0_dp/)
      uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
      call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
      call MIO_Print('Moire cell indices (m1 m2 m3 m4): '//trim(num2str(m(1)))//&
        '  '//trim(num2str(m(2)))//'  '//trim(num2str(m(3)))//'  '//trim(num2str(m(4))),&
        'cell')
      g = m(1)**2 + m(2)**2 + m(1)*m(2)
      delta = sqrt(real(m(3)**2 + m(4)**2 + m(3)*m(4))/g)
      phi = acos((2.0_dp*m(1)*m(3)+2.0_dp*m(2)*m(4) + m(1)*m(4) + m(2)*m(3))/(2.0_dp*delta*g))
      if (delta<1.0_dp) then
         delta = 1.0_dp/delta
      else
         i = m(1)
         m(1) = m(3)
         m(3) = i
         i = m(2)
         m(2) = m(4)
         m(4) = i
      end if
      aBN = aG*delta
      call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' deg','cell')
      call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4))//' Ang','cell')
      call MIO_Print('BN lattice parameter: '//trim(num2str(aBN,4))//' Ang','cell')
      m1(:,1) = [m(1),m(2)]
      m1(:,2) = [-m(2),m(1)+m(2)]
      ucell = 0.0_dp
      ucell(1:2,1:2) = matmul(uc,m1)
      ucell(3,3) = h
      g = norm(ucell(:,1))
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
      ucell(:,:2) = ucell(:,:2)*sCell
   else if (MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell')) then
      uc(:,1) = (/aG,0.0_dp/)
      uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
      call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
      call MIO_Print('Moire cell indices (m1 m2 m3 m4): '//trim(num2str(m(1)))//&
        '  '//trim(num2str(m(2)))//'  '//trim(num2str(m(3)))//'  '//trim(num2str(m(4))),&
        'cell')
      g = m(1)**2 + m(2)**2 + m(1)*m(2)
      delta = sqrt(real(m(3)**2 + m(4)**2 + m(3)*m(4))/g)
      phi = acos((2.0_dp*m(1)*m(3)+2.0_dp*m(2)*m(4) + m(1)*m(4) + m(2)*m(3))/(2.0_dp*delta*g))
      if (delta<1.0_dp) then
         delta = 1.0_dp/delta
      else
         i = m(1)
         m(1) = m(3)
         m(3) = i
         i = m(2)
         m(2) = m(4)
         m(4) = i
      end if
      aBN = aG*delta
      call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' deg','cell')
      call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4))//' Ang','cell')
      call MIO_Print('BN lattice parameter: '//trim(num2str(aBN,4))//' Ang','cell')
      m1(:,1) = [m(1),m(2)]
      m1(:,2) = [-m(2),m(1)+m(2)]
      ucell = 0.0_dp
      ucell(1:2,1:2) = matmul(uc,m1)
      ucell(3,3) = h
      g = norm(ucell(:,1))
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
      ucell(:,:2) = ucell(:,:2)*sCell
   else if (MIO_StringComp(str,'Ribbons')) then
      call MIO_InputParameter('Structure.RibbonType',str,'Zigzag')
        !! - `Structure.RibbonType`: TypeOfSystem Ribbons: edge of the ribbons, Zigzag or Armchair.
      call MIO_InputParameter('Structure.GrapheneWidth',nC,9)
        !! - `Structure.GrapheneWidth`: TypeOfSystem Ribbons: width of the graphene ribbon, in rows of the ribbon cell.
      call MIO_InputParameter('Structure.BNWidth',nBN,9)
        !! - `Structure.BNWidth`: TypeOfSystem Ribbons: width of the hBN ribbon, in rows of the ribbon cell.
      aBN = aG
      if (MIO_StringComp(str,'Zigzag')) then
         ucell(:,1) = [aG,0.0_dp,0.0_dp]
         ucell(:,2) = [0.0_dp,aG*(nC+nBN)*sqrt(3.0_dp)/2.0_dp,0.0_dp]
         ucell(:,3) = [0.0_dp,0.0_dp,h]
      else if (MIO_StringComp(str,'Armchair')) then
         ucell(:,1) = [aG*sqrt(3.0_dp),0.0_dp,0.0_dp]
         ucell(:,2) = [0.0_dp,aG*(nC+nBN)/2.0_dp,0.0_dp]
         ucell(:,3) = [0.0_dp,0.0_dp,h]
      else
         call MIO_Kill("Type of system not recognized. Check 'RibbonType'", &
           'cell','CellGet')
      end if
   else if (MIO_StringComp(str,'ReadXYZ')) then
      call MIO_InputParameter('Structure.XYZFile',str,trim(prefix)//'.xyz')
        !! - `Structure.XYZFile`: TypeOfSystem ReadXYZ: structure file. Its first three lines are the lattice vectors [Angstrom],
        !!   the atoms follow.
      open(1,FILE=str,STATUS='old')
      do i=1,3
         read(1,*) (ucell(j,i),j=1,3)
      end do
      do i=1,3
         write(dbgline,'(a,i1,a,3f16.8)') 'lattice vector ',i,' read from the structure file [Ang]:',ucell(:,i)
         call MIO_PrintDebug(dbgline,'cell')
      end do
      close(1)
      call MIO_InputParameter('Structure.LatticeParameterBN',aBN,2.505_dp)
      g = norm(ucell(:,1))
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
   else
      if (MIO_StringComp(str,'Hybrid') .or. MIO_StringComp(str,'BLtoSLYoungju')) then
         call MIO_Kill('TypeOfSystem '//trim(str)//' is not supported at present: parts of the code refer to it, '// &
           'but no cell is defined for it.','cell','CellGet')
      end if
      call MIO_Kill('Type of system not recognized: '//trim(str),'cell','CellGet')
   end if
   if (sCell /= 1) then
      if (l) then
         call MIO_Print('Length of supercell: '//trim(num2str(norm(ucell(:,1)),5))// &
            ' Ang x '//trim(num2str(norm(ucell(:,2)),5))//' Ang','cell')
      else
         call MIO_Print('Length of supercell: '//trim(num2str(norm(ucell(:,1)),5))// &
            ' Ang','cell')
      end if
   end if

   vn = CrossProd(ucell(:,1),ucell(:,2))
   volume = dot_product(ucell(:,3),vn)
   area = norm(vn)
   rcell(:,1) = twopi*CrossProd(ucell(:,2),ucell(:,3))/volume
   rcell(:,2) = twopi*CrossProd(ucell(:,3),ucell(:,1))/volume
   rcell(:,3) = twopi*CrossProd(ucell(:,1),ucell(:,2))/volume

#ifdef TIMER
   call MIO_TimerStop('cell')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('CellGet',1)
#endif /* DEBUG */

end subroutine CellGet

end module cell
