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
            !!   first lattice, (m3, m4) on the second; the twist angle and lattice ratio follow from them and are printed. Any
            !!   signs, both pairs nonzero; a non-primitive cell is built as given
            !!   (docs/user-guide/commensurate-cells.md).
          call MIO_Print('Moire cell indices (m1 m2 m3 m4): '//trim(num2str(m(1)))//&
            '  '//trim(num2str(m(2)))//'  '//trim(num2str(m(3)))//'  '//trim(num2str(m(4))),&
            'cell')
          ! the atoms of MoireEncapsulatedBilayer / TwistedBilayer are built on CellSize x CellSize cells (AtomsPos),
          ! not in this cell: the layers did not fill it
          if (MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer')) &
            call MIO_Kill('Structure.BasedOnMoireCellParameters is not implemented for TypeOfSystem '//trim(str)// &
              ': its atoms are built on Structure.CellSize cells, not in the cell of Structure.MoireCellParameters. '// &
              'Use TypeOfSystem TwistedBilayerBasedOnMoireCell or MoireEncapsulatedBilayerBasedOnMoireCell',&
              'cell','CellGet')
          call CellMoireReport(m,.false.)
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
      call CellMoireReport(m,.true.)
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
      call CellMoireReport(m,.true.)
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
      call CellMoireReport(m,.true.)
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
      call CellMoireReport(m,.true.)
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
      call CellMoireReport(m,.true.)
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
      call CellMoireReport(m,.true.)
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

!> Check and describe the commensurate cell of Structure.MoireCellParameters m1 m2 m3 m4 (Hermann's four indices).
!! With tau = exp(i pi/3), the supercell vector A1 = m1 a1 + m2 a2 of the first lattice is the Eisenstein integer
!! alpha = m1 + m2 tau (in units of that lattice constant), and the same vector on the second lattice is
!! beta = m3 + m4 tau; A2 is A1 turned by +60 degrees. Any nonzero alpha, beta (signs and zeros allowed) define a
!! commensurate pair: the second lattice is the first one scaled by |alpha/beta| and turned by theta = arg(alpha/beta).
!! The coincidence lattice of the two layers is (alpha/gcd(alpha,beta)) Z[tau], gcd in the Eisenstein integers, so the
!! cell is primitive iff gcd(alpha,beta) is a unit; otherwise it holds N(gcd) = |gcd|^2 primitive cells. A common
!! integer factor of the four indices is only one case: (4,1,1,4) has none and holds 3 primitive cells.
!! The cell is reported, never reduced: a non-primitive cell is built as given (supercells for symmetry breaking).
!! twoLayers: both pairs build atoms (a zero pair is refused); otherwise one layer is built and a zero pair is noted.
subroutine CellMoireReport(m, twoLayers)

   use constants,            only : pi

   ! Arguments
   integer, intent(in) :: m(4)
   logical, intent(in) :: twoLayers

   ! Named constants
   integer, parameter :: ik = selected_int_kind(18)
   integer(ik), parameter :: units(2,6) = reshape([1_ik,0_ik, 0_ik,1_ik, -1_ik,1_ik, -1_ik,0_ik, 0_ik,-1_ik, 1_ik,-1_ik], &
                                                  [2,6])

   ! Other local variables
   integer :: k
   integer(ik) :: al(2), be(2), gc(2), ap(2), bp(2), x(2)
   real(dp) :: theta, za(2), zb(2)
   logical :: swapped

   al = int(m(1:2),ik)
   be = int(m(3:4),ik)
   if (EisNorm(al) == 0 .or. EisNorm(be) == 0) then
      if (twoLayers) call MIO_Kill('Structure.MoireCellParameters: (m1, m2) and (m3, m4) must both be nonzero, '// &
         'a zero pair defines no cell','cell','CellMoireReport')
      call MIO_Print('Hermann cell: one index pair is zero; the twist angle and the lattice ratio are undefined',&
         'cell')
      return
   end if
   ! Cartesian components of alpha, beta in units of the lattice constant; theta = arg(alpha conj(beta))
   za = [real(al(1),dp) + 0.5_dp*real(al(2),dp), sqrt(3.0_dp)/2.0_dp*real(al(2),dp)]
   zb = [real(be(1),dp) + 0.5_dp*real(be(2),dp), sqrt(3.0_dp)/2.0_dp*real(be(2),dp)]
   theta = atan2(za(2)*zb(1) - za(1)*zb(2), za(1)*zb(1) + za(2)*zb(2))*180.0_dp/pi
   call MIO_Print('Hermann cell: the (m3, m4) lattice is the (m1, m2) lattice turned by '//trim(num2str(theta,6))// &
      ' deg (signed, counterclockwise)','cell')
   ! The layers are built with the pair of the larger norm first (the pairs are exchanged when N(beta) >= N(alpha)):
   ! the layer built second is turned by theta, or by -theta after the exchange
   swapped = EisNorm(be) >= EisNorm(al)
   if (twoLayers) then
      if (swapped) then
         call MIO_Print('Hermann cell: layer 1 (upper) is built from (m3, m4), layer 2 from (m1, m2), turned by '// &
            trim(num2str(-theta,6))//' deg with respect to layer 1','cell')
      else
         call MIO_Print('Hermann cell: layer 1 (upper) is built from (m1, m2), layer 2 from (m3, m4), turned by '// &
            trim(num2str(theta,6))//' deg with respect to layer 1','cell')
      end if
   end if
   ! Terms that take the angle of the cell use the unsigned angle, i.e. the sense of the (a, b, b, a) cells with a > b
   if ((swapped .and. theta > 1.0e-9_dp) .or. (.not. swapped .and. theta < -1.0e-9_dp)) then
      call MIO_Warning('Structure.MoireCellParameters: layer 2 of this cell is turned clockwise. The structure is '// &
         'built as the indices say, but the terms that take the angle of the cell, if used (moire effective model '// &
         'with Structure.BasedOnMoireCellParameters, Spectral.RotateReferenceSystem, the analytic twisted-bilayer '// &
         'terms), use the unsigned angle, the sense of the (a, b, b, a) cells with a > b. The mirror cell (m1+m2, -m2, '// &
         'm3+m4, -m4) = '//trim(num2str(m(1)+m(2)))//' '//trim(num2str(-m(2)))//' '//trim(num2str(m(3)+m(4)))//' '// &
         trim(num2str(-m(4)))//' is turned counterclockwise','cell')
   end if
   gc = EisGcd(al,be)
   if (EisNorm(gc) == 1) then
      call MIO_Print('Hermann cell: primitive (the layers have no smaller common cell)','cell')
   else
      ap = EisDivExact(al,gc)
      bp = EisDivExact(be,gc)
      ! the unit multiple with (m1, m2) in the 60-degree sector m1 > 0, m2 >= 0 (the same cell)
      do k=1,6
         x = EisMul(units(:,k),ap)
         if (x(1) > 0 .and. x(2) >= 0) exit
      end do
      bp = EisMul(units(:,k),bp)
      ap = x
      call MIO_Print('Hermann cell: NOT primitive, it holds '//trim(num2str(int(EisNorm(gc))))//' primitive cells '// &
         '(built as given). Primitive cell with the same twist: '//trim(num2str(int(ap(1))))//' '// &
         trim(num2str(int(ap(2))))//' '//trim(num2str(int(bp(1))))//' '//trim(num2str(int(bp(2)))),'cell')
   end if

end subroutine CellMoireReport

!> Norm |a + b tau|^2 = a^2 + a b + b^2 of an Eisenstein integer (tau = exp(i pi/3)).
pure function EisNorm(x) result(n)
   integer, parameter :: ik = selected_int_kind(18)
   integer(ik), intent(in) :: x(2)
   integer(ik) :: n
   n = x(1)*x(1) + x(1)*x(2) + x(2)*x(2)
end function EisNorm

!> Product (a + b tau)(c + d tau), with tau^2 = tau - 1.
pure function EisMul(x,y) result(z)
   integer, parameter :: ik = selected_int_kind(18)
   integer(ik), intent(in) :: x(2), y(2)
   integer(ik) :: z(2)
   z = [x(1)*y(1) - x(2)*y(2), x(1)*y(2) + x(2)*y(1) + x(2)*y(2)]
end function EisMul

!> Euclidean division x = q y + r, N(r) < N(y): q is the lattice point nearest to x/y = x conj(y)/N(y).
pure subroutine EisDivMod(x,y,q,r)
   integer, parameter :: ik = selected_int_kind(18)
   integer(ik), intent(in) :: x(2), y(2)
   integer(ik), intent(out) :: q(2), r(2)
   integer(ik) :: num(2), q0(2), qq(2), rr(2), ny
   integer :: da, db
   num = EisMul(x,[y(1)+y(2),-y(2)])
   ny = EisNorm(y)
   q0 = nint(real(num,dp)/real(ny,dp),ik)
   q = q0
   r = x - EisMul(q0,y)
   ! rounding in the skew basis is not always the nearest point: try the neighbours
   do da=-1,1; do db=-1,1
      qq = q0 + [int(da,ik),int(db,ik)]
      rr = x - EisMul(qq,y)
      if (EisNorm(rr) < EisNorm(r)) then
         q = qq
         r = rr
      end if
   end do; end do
end subroutine EisDivMod

!> Quotient of an exact division in the Eisenstein integers.
pure function EisDivExact(x,y) result(q)
   integer, parameter :: ik = selected_int_kind(18)
   integer(ik), intent(in) :: x(2), y(2)
   integer(ik) :: q(2), r(2)
   call EisDivMod(x,y,q,r)
end function EisDivExact

!> Greatest common divisor in the Eisenstein integers (Euclid's algorithm; defined up to a unit).
pure function EisGcd(x,y) result(g)
   integer, parameter :: ik = selected_int_kind(18)
   integer(ik), intent(in) :: x(2), y(2)
   integer(ik) :: g(2), a(2), b(2), q(2), r(2)
   a = x
   b = y
   do while (EisNorm(b) /= 0)
      call EisDivMod(a,b,q,r)
      a = b
      b = r
   end do
   g = a
end function EisGcd

end module cell
