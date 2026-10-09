module magf

   use mio

   implicit none

   PRIVATE

   real(dp), save, public :: Bmag=0.0_dp , BmagZeeman=0.0_dp, HaldanePhase
   integer, save, public :: mBi=0, mBf=0, mStep=1 , mPhii=0, mPhif=0, mPhiStep=0
   logical, save, public :: magfield=.false. , HaldPhase=.false.

   public :: MagfInit
   public :: MagfValue
   public :: HaldPhaseValue
   public :: HaldPhaseInit

contains

subroutine MagfInit()

   use cell,                 only : area
   use constants,            only : fluxq

   character(len=80) :: line
   integer :: id
   real(dp) ::BStep

#ifdef DEBUG
   call MIO_Debug('MagfInit',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('magf')
#endif /* TIMER */

   if (MIO_InputSearchLabel('MagField',line,id)) then
      call MIO_InputParameter('MagField',Bmag,0.0_dp)
      !mBi = nint(area*1.0d-20*Bmag/(2.0_dp*fluxq))
      mBi = nint(area*1.0d-20*Bmag/(fluxq))
      mBf = mBi
      print*, "chosen magnetic field is ", Bmag, " with index ", mBf
   else if (MIO_InputSearchLabel('MagField.Integer',line,id)) then
      call MIO_InputParameter('MagField.Integer',mBi,0)
      mBf = mBi
   else
      if (MIO_InputSearchLabel('MagField.InitialB',line,id)) then
         call MIO_InputParameter('MagField.InitialB',Bmag,0.0_dp)
         !mBi = nint(area*1.0d-20*Bmag/(2.0_dp*fluxq))
         mBi = nint(area*1.0d-20*Bmag/(fluxq))
         call MIO_InputParameter('MagField.FinalB',Bmag,0.0_dp)
         !mBf = nint(area*1.0d-20*Bmag/(2.0_dp*fluxq))
         mBf = nint(area*1.0d-20*Bmag/(fluxq))
         if (mBi>mBf) then
            mBi = 0
            mBf = 0
         end if
      else if (MIO_InputSearchLabel('MagField.InitialInteger',line,id)) then
         call MIO_InputParameter('MagField.InitialInteger',mBi,0)
         call MIO_InputParameter('MagField.FinalInteger',mBf,0)
         if (mBi>mBf) then
            mBi = 0
            mBf = 0
         end if
      else
         mBi = 0
         mBf = 0
      end if
   end if
   if (mBi==0 .and. mBf==0) then
      magfield = .false.
      mStep = 1
   else
      magfield = .true.
      if (MIO_InputSearchLabel('MagField.Step',line,id)) then
         call MIO_InputParameter('MagField.Step',BStep,1.0_dp)
         !mStep = nint(area*1.0d-20*BStep/(2.0_dp*fluxq))
         mStep = nint(area*1.0d-20*BStep/(fluxq))
      else
         call MIO_InputParameter('MagField.IntegerStep',mStep,1)
      end if
      if (mStep<=0) mStep = 1
   end if

#ifdef TIMER
   call MIO_TimerStop('magf')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('MagfInit',1)
#endif /* DEBUG */

end subroutine MagfInit

subroutine MagfValue(mB)

   use cell,                 only : area, sCell
   use constants,            only : fluxq
   use name,                 only : prefix, sysname

   integer, intent(in) :: mB
   character(len=80) :: str

#ifdef DEBUG
   call MIO_Debug('MagfValue',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('magf')
#endif /* TIMER */

   ! Read Zeeman magnetic field parameter (separate from Landau level Bmag)
   call MIO_InputParameter('BmagZeeman',BmagZeeman,0.0_dp)

   call MIO_InputParameter('TypeOfSystem',str,'Graphene')
   if (magfield) then
      !if ((MIO_StringComp(str,'ReadXYZ')) .and. sCell.ne.1) then ! for my recursion calculations on sandwiched systems for instance
      !    Bmag = 2.0_dp*mB*fluxq/(area*1.0d-20*sCell*sCell)
      !else
          !Bmag = 2.0_dp*mB*fluxq/(area*1.0d-20)
          Bmag = mB*fluxq/(area*1.0d-20)
      !end if
      if (mBi /= mBf) then
         write(prefix,'(a,I5.5,a,f0.4)') 'm',mB,'_B', Bmag
         call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'magf')
         prefix = './'//trim(prefix)
         call system('mkdir '//trim(prefix))
         prefix = trim(prefix)//'/'//trim(sysname)
      end if
      call MIO_Print('Magnetic field: '//trim(num2str(Bmag,3)),'magf')
      call MIO_Print('m: '//trim(num2str(mB)),'magf')
   else
      Bmag = 0.0_dp
   end if


#ifdef TIMER
   call MIO_TimerStop('magf')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('MagfValue',1)
#endif /* DEBUG */

end subroutine MagfValue


subroutine HaldPhaseInit()

   use cell,                 only : area
   use constants,            only : fluxq

   character(len=80) :: line
   integer :: id
   real(dp) ::BStep

   logical :: l, ll, lll

#ifdef DEBUG
   call MIO_Debug('HaldPhaseInit',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('magf')
#endif /* TIMER */

   !if (MIO_InputSearchLabel('HaldaneSpecifyFlux',line,id)) then
   call MIO_InputParameter('HaldaneSpecifyFlux',l,.false.)
   call MIO_InputParameter('HaldaneSpecifyPhase',ll,.false.)
   call MIO_InputParameter('HaldaneSpecifyRange',lll,.false.)
   if (l) then
      mPhii = 0
      mPhif = mPhii
   else if (ll) then
      mPhii = 0
      mPhif = mPhii
   else
      if (lll) then
         call MIO_InputParameter('HaldPhase.InitialInteger',mPhii,0)
         call MIO_InputParameter('HaldPhase.FinalInteger',mPhif,0)
         if (mPhii>mPhif) then
            mPhii = 0
            mPhif = 0
         end if
      end if
   end if
   if (mPhii==0 .and. mPhif==0) then
      HaldPhase = .false.
      mPhiStep = 1
   else
      HaldPhase = .true.
      call MIO_InputParameter('HaldPhase.Step',mPhiStep,1)
      if (mPhiStep<=0) mPhiStep = 1
   end if

#ifdef TIMER
   call MIO_TimerStop('magf')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('HaldPhaseInit',1)
#endif /* DEBUG */

end subroutine HaldPhaseInit


subroutine HaldPhaseValue(mPhi)

   use cell,                 only : area
   use constants,            only : fluxq, pi
   use name,                 only : prefix, sysname

   integer, intent(in) :: mPhi

#ifdef DEBUG
   call MIO_Debug('HaldPhaseValue',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('magf')
#endif /* TIMER */

   if (HaldPhase) then
      HaldanePhase = 1.0_dp/(mPhif-mPhii)*mPhi*2.0_dp*pi
      if (mPhii /= mPhif) then
         write(prefix,'(a,I5.5,a,f0.4)') 'm',mPhi,'_B', HaldanePhase
         call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'magf')
         prefix = './'//trim(prefix)
         call system('mkdir '//trim(prefix))
         prefix = trim(prefix)//'/'//trim(sysname)
      end if
      call MIO_Print('Haldane phase: '//trim(num2str(HaldanePhase,3)),'HaldPhase')
      call MIO_Print('m: '//trim(num2str(mPhi)),'magf')
    end if

#ifdef TIMER
   call MIO_TimerStop('magf')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('HaldPhaseValue',1)
#endif /* DEBUG */

end subroutine HaldPhaseValue

end module magf
