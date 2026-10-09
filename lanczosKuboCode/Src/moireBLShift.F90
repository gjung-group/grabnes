module moireBLShift

   use mio

   implicit none

   PRIVATE

   integer, save, public :: mSi=0, mSf=0, mSStep=1
   real(dp), save, public :: tauX1, tauY1, tauX2, tauY2

   public :: moireBLShiftInit
   public :: moireBLShiftValue

contains

subroutine moireBLShiftInit()

   use cell,                 only : area

   logical :: u, v

#ifdef DEBUG
   call MIO_Debug('moireBLShiftInit',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('moireBLShift')
#endif /* TIMER */

   call MIO_InputParameter('moireBLShiftDefinedValues',u,.false.)
   call MIO_InputParameter('moireBLShiftRange',v,.false.)
   if (u) then
      mSi = 0
      mSf = mSi
      mSStep = 1
      print*, "hereclup", mSi, mSf
   else if (v) then
      call MIO_InputParameter('moireBLShiftRangeSteps',mSStep,10)
      mSi = 0
      mSf = mSi + mSStep
      print*, "hereclap", mSi, mSf
   end if

#ifdef TIMER
   call MIO_TimerStop('moireBLShift')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('moireBLShiftInit',1)
#endif /* DEBUG */

end subroutine moireBLShiftInit

subroutine moireBLShiftValue(mS)

   use name,                 only : prefix, sysname

   integer, intent(in) :: mS
   logical :: u, v
   real(dp) :: xmin1, xmin2, xmax1, xmax2, ymin1, ymin2, ymax1, ymax2

#ifdef DEBUG
   call MIO_Debug('moireBLShiftValue',0)
#endif /* DEBUG */
#ifdef TIMER
   call MIO_TimerCount('moireBLShift')
#endif /* TIMER */

   call MIO_InputParameter('moireBLShiftDefinedValues',u,.false.)
   call MIO_InputParameter('moireBLShiftRange',v,.false.)
   if (u) then
      call MIO_InputParameter('MoireBilayerBottomXShift',tauX1,0.0_dp)
      call MIO_InputParameter('MoireBilayerBottomYShift',tauY1,0.0_dp)
      call MIO_InputParameter('MoireBilayerTopXShift',tauX2,0.0_dp)
      call MIO_InputParameter('MoireBilayerTopYShift',tauY2,0.0_dp)
   else if (v) then
      call MIO_InputParameter('MoireBilayerBottomXMin',xmin1,0.0_dp)
      call MIO_InputParameter('MoireBilayerBottomXMax',xmax1,0.0_dp)
      call MIO_InputParameter('MoireBilayerBottomYMin',ymin1,0.0_dp)
      call MIO_InputParameter('MoireBilayerBottomYMax',ymax1,0.0_dp)
      call MIO_InputParameter('MoireBilayerTopXMin',xmin2,0.0_dp)
      call MIO_InputParameter('MoireBilayerTopXMax',xmax2,0.0_dp)
      call MIO_InputParameter('MoireBilayerTopYMin',ymin2,0.0_dp)
      call MIO_InputParameter('MoireBilayerTopYMax',ymax2,0.0_dp)
      tauX1= xmin1 + mS*(xmax1-xmin1)/mSStep
      tauY1= ymin1 + mS*(ymax1-ymin1)/mSStep
      tauX2= xmin2 + mS*(xmax2-xmin2)/mSStep
      tauY2= ymin2 + mS*(ymax2-ymin2)/mSStep
      if (mSi /= mSf) then
         write(prefix,'(a,f0.4,a,f0.4,a,f0.4,a,f0.4)') 'X1_',tauX1,'_Y1_', tauY1,'_X2_',tauX2,'_Y2_',tauY2
         call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'moireBLShift')
         prefix = './'//trim(prefix)
         call system('mkdir '//trim(prefix))
         prefix = trim(prefix)//'/'//trim(sysname)
      end if
   end if

#ifdef TIMER
   call MIO_TimerStop('moireBLShift')
#endif /* TIMER */
#ifdef DEBUG
   call MIO_Debug('moireBLShiftValue',1)
#endif /* DEBUG */

end subroutine moireBLShiftValue

end module moireBLShift
