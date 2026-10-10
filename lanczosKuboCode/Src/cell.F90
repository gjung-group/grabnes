module cell

   use mio

   implicit none

   PRIVATE

   real(dp), public, save :: ucell(3,3), rcell(3,3), volume, area, aG, aBN
   integer, public, save :: sCell, sCell2

   public :: CellGet

contains

subroutine CellGet()

   use math
   use constants,            only : twopi, pi
   use name,                 only : prefix

   real(dp) :: h, vn(3), uc(2,2), g, delta, phi
   integer :: superCellX, superCellY
   character(len=80) :: str
   integer :: n, n2,  m(4), m1(2,2), i, j, nC, nBN
   logical l, ll

#ifdef DEBUG
   call MIO_Debug('CellGet',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('cell')
#endif /* TIMER */

   call MIO_InputParameter('Structure.LatticeParameter',aG,2.46_dp)
   call MIO_InputParameter('Structure.LatticeParameterBN',aBN,2.505_dp)
   call MIO_InputParameter('Run.TypeOfSystem',str,'Graphene')
   call MIO_InputParameter('Structure.CellHeight',h,40.0_dp)
   call MIO_InputParameter('Structure.SuperCell',sCell,1)
   ! Default of the interlayer distance: 3.35 Ang everywhere since October 2026 (3.22 Ang before in the
   ! generators, the neighbor search and parts of the Hamiltonian; it matters for structures read from a file too)
   if (.not. (MIO_StringComp(str,'Graphene') .or. MIO_StringComp(str,'BoronNitride'))) then
      if (.not. MIO_InputSearchLabel('Structure.InterlayerDistance')) then
         call MIO_Print('NOTE: InterlayerDistance is not given: 3.35 Ang is used (the default was 3.22 Ang '// &
           'before October 2026)','cell')
      end if
   end if
   call MIO_InputParameter('Structure.SuperCellAsymmetric',l,.false.)
   call MIO_InputParameter('Structure.SuperCellY',sCell2,1)
   if (MIO_StringComp(str,'Graphene') .or. MIO_StringComp(str,'BoronNitride') &
    .or. MIO_StringComp(str,'MoireEncapsulatedBilayer') .or. MIO_StringComp(str,'TwistedBilayer')) then
      call MIO_InputParameter('Structure.BasedOnMoireCellParameters',ll,.false.)
      if (ll) then
          uc(:,1) = (/aG,0.0_dp/)
          uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
          call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
          call MIO_Print('Parameters for moiré generation: '//trim(num2str(m(1)))//&
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
          call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' degrees','cell')
          call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4)),'cell')
          call MIO_Print('BN lattice parameter:       '//trim(num2str(aBN,4)),'cell')
          m1(:,1) = [m(1),m(2)]
          m1(:,2) = [-m(2),m(1)+m(2)]
          ucell = 0.0_dp
          ucell(1:2,1:2) = matmul(uc,m1)
          ucell(3,3) = h
          g = norm(ucell(:,1))
          call MIO_Print('')
          call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
            ' Ang','cell')
          ucell(:,:2) = ucell(:,:2)*sCell
          aBN = aG ! added this, not sure if necessary
      else
          call MIO_InputParameter('Structure.CellSize',n,50)
          call MIO_Print('Length of unit cell: '//trim(num2str(n*aG,5))// &
            ' Ang','cell')
             n = n*sCell
          ucell(:,1) = [aG,0.0_dp,0.0_dp]
          ucell(:,2) = [aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp,0.0_dp]
            ucell = ucell*n
          ucell(:,3) = [0.0_dp,0.0_dp,h]
          aBN = aG
      end if
   else if (MIO_StringComp(str,'Graphene_Over_BN')) then
      uc(:,1) = (/aG,0.0_dp/)
      uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
      call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
      call MIO_Print('Parameters for moiré generation: '//trim(num2str(m(1)))//&
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
      call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' degrees','cell')
      call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4)),'cell')
      call MIO_Print('BN lattice parameter:       '//trim(num2str(aBN,4)),'cell')
      m1(:,1) = [m(1),m(2)]
      m1(:,2) = [-m(2),m(1)+m(2)]
      ucell = 0.0_dp
      ucell(1:2,1:2) = matmul(uc,m1)
      ucell(3,3) = h
      g = norm(ucell(:,1))
      call MIO_Print('')
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
      ucell(:,:2) = ucell(:,:2)*sCell
   else if (MIO_StringComp(str,'BilayerGraphene_Over_BN')) then
      uc(:,1) = (/aG,0.0_dp/)
      uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
      call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
      call MIO_Print('Parameters for moiré generation: '//trim(num2str(m(1)))//&
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
      call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' degrees','cell')
      call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4)),'cell')
      call MIO_Print('BN lattice parameter:       '//trim(num2str(aBN,4)),'cell')
      m1(:,1) = [m(1),m(2)]
      m1(:,2) = [-m(2),m(1)+m(2)]
      ucell = 0.0_dp
      ucell(1:2,1:2) = matmul(uc,m1)
      ucell(3,3) = h
      g = norm(ucell(:,1))
      call MIO_Print('')
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
      ucell(:,:2) = ucell(:,:2)*sCell
   else if (MIO_StringComp(str,'TwistedBilayerBasedOnMoireCell')) then
      uc(:,1) = (/aG,0.0_dp/)
      uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
      call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
      call MIO_Print('Parameters for moiré generation: '//trim(num2str(m(1)))//&
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
      call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' degrees','cell')
      call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4)),'cell')
      call MIO_Print('BN lattice parameter:       '//trim(num2str(aBN,4)),'cell')
      m1(:,1) = [m(1),m(2)]
      m1(:,2) = [-m(2),m(1)+m(2)]
      ucell = 0.0_dp
      ucell(1:2,1:2) = matmul(uc,m1)
      ucell(3,3) = h
      g = norm(ucell(:,1))
      call MIO_Print('')
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
      ucell(:,:2) = ucell(:,:2)*sCell
   else if (MIO_StringComp(str,'TrilayerBasedOnMoireCell')) then
      uc(:,1) = (/aG,0.0_dp/)
      uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
      call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
      call MIO_Print('Parameters for moiré generation: '//trim(num2str(m(1)))//&
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
      call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' degrees','cell')
      call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4)),'cell')
      call MIO_Print('BN lattice parameter:       '//trim(num2str(aBN,4)),'cell')
      m1(:,1) = [m(1),m(2)]
      m1(:,2) = [-m(2),m(1)+m(2)]
      ucell = 0.0_dp
      ucell(1:2,1:2) = matmul(uc,m1)
      ucell(3,3) = h
      g = norm(ucell(:,1))
      call MIO_Print('')
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
      ucell(:,:2) = ucell(:,:2)*sCell
   else if (MIO_StringComp(str,'TwistedBilayerBasedOnMoireCellRectangular')) then
      uc(:,1) = (/aG,0.0_dp/)
      uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
      call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
      call MIO_Print('Parameters for moiré generation: '//trim(num2str(m(1)))//&
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
      call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' degrees','cell')
      call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4)),'cell')
      call MIO_Print('BN lattice parameter:       '//trim(num2str(aBN,4)),'cell')
      m1(:,1) = [m(1),m(2)]
      m1(:,2) = [-m(2),m(1)+m(2)]
      ucell = 0.0_dp
      ucell(1:2,1:2) = matmul(uc,m1)
      ucell(3,3) = h
      g = norm(ucell(:,1))
      call MIO_Print('')
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
      ucell(:,:2) = ucell(:,:2)*sCell
   else if (MIO_StringComp(str,'MoireEncapsulatedBilayerBasedOnMoireCell')) then
      uc(:,1) = (/aG,0.0_dp/)
      uc(:,2) = (/aG/2.0_dp,sqrt(3.0_dp)*aG/2.0_dp/)
      call MIO_InputParameter('Structure.MoireCellParameters',m,[0,0,0,0])
      call MIO_Print('Parameters for moiré generation: '//trim(num2str(m(1)))//&
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
      call MIO_Print('Angle: '//trim(num2str(phi*180.0_dp/pi,4))//' degrees','cell')
      call MIO_Print('Graphene lattice parameter: '//trim(num2str(aG,4)),'cell')
      call MIO_Print('BN lattice parameter:       '//trim(num2str(aBN,4)),'cell')
      m1(:,1) = [m(1),m(2)]
      m1(:,2) = [-m(2),m(1)+m(2)]
      ucell = 0.0_dp
      ucell(1:2,1:2) = matmul(uc,m1)
      ucell(3,3) = h
      g = norm(ucell(:,1))
      call MIO_Print('')
      call MIO_Print('Length of unit cell: '//trim(num2str(g,5))// &
        ' Ang','cell')
      ucell(:,:2) = ucell(:,:2)*sCell
   else if (MIO_StringComp(str,'Ribbons')) then
      call MIO_InputParameter('Structure.RibbonType',str,'Zigzag')
      call MIO_InputParameter('Structure.GrapheneWidth',nC,9)
      call MIO_InputParameter('Structure.BNWidth',nBN,9)
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
      open(1,FILE=str,STATUS='old')
      do i=1,3
         read(1,*) (ucell(j,i),j=1,3)
      end do
#ifdef DEBUG
      print*, "ucell was read from ReadXYZ as: "
      print*, ucell(:,1)
      print*, ucell(:,2)
      print*, ucell(:,3)
#endif /* DEBUG */
      close(1)
      call MIO_InputParameter('Structure.LatticeParameterBN',aBN,2.505_dp)
      g = norm(ucell(:,1))
      call MIO_Print('')
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
            ' Ang'//' X '//trim(num2str(norm(ucell(:,2)),5))//' Ang','cell')
      else
         call MIO_Print('Length of supercell: '//trim(num2str(norm(ucell(:,1)),5))// &
            ' Ang','cell')
      end if
   end if
   call MIO_Print('')

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
