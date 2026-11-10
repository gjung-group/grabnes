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

   call MIO_InputParameter('Moire.BLShiftDefinedValues',u,.false.)
     !! - `Moire.BLShiftDefinedValues`: moireBLShift.Calc: one calculation with the shifts Moire.Bilayer{Bottom,Top}{X,Y}Shift.
   call MIO_InputParameter('Moire.BLShiftRange',v,.false.)
     !! - `Moire.BLShiftRange`: moireBLShift.Calc: sweep each shift from its Moire.Bilayer*Min to *Max in Moire.BLShiftRangeSteps
     !!   steps.
   if (u) then
      mSi = 0
      mSf = mSi
      mSStep = 1
   else if (v) then
      call MIO_InputParameter('Moire.BLShiftRangeSteps',mSStep,10)
        !! - `Moire.BLShiftRangeSteps`: Moire.BLShiftRange: number of steps of the shift sweep.
      mSi = 0
      mSf = mSi + mSStep
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

   call MIO_InputParameter('Moire.BLShiftDefinedValues',u,.false.)
   call MIO_InputParameter('Moire.BLShiftRange',v,.false.)
   if (u) then
      call MIO_InputParameter('Moire.BilayerBottomXShift',tauX1,0.0_dp)
        !! - `Moire.BilayerBottomXShift`: Moire.BLShiftDefinedValues: x shift of the bottom layer (written into the output prefix).
      call MIO_InputParameter('Moire.BilayerBottomYShift',tauY1,0.0_dp)
        !! - `Moire.BilayerBottomYShift`: Moire.BLShiftDefinedValues: y shift of the bottom layer.
      call MIO_InputParameter('Moire.BilayerTopXShift',tauX2,0.0_dp)
        !! - `Moire.BilayerTopXShift`: Moire.BLShiftDefinedValues: x shift of the top layer.
      call MIO_InputParameter('Moire.BilayerTopYShift',tauY2,0.0_dp)
        !! - `Moire.BilayerTopYShift`: Moire.BLShiftDefinedValues: y shift of the top layer.
   else if (v) then
      call MIO_InputParameter('Moire.BilayerBottomXMin',xmin1,0.0_dp)
        !! - `Moire.BilayerBottomXMin`: Moire.BLShiftRange: first x shift of the bottom layer.
      call MIO_InputParameter('Moire.BilayerBottomXMax',xmax1,0.0_dp)
        !! - `Moire.BilayerBottomXMax`: Moire.BLShiftRange: last x shift of the bottom layer.
      call MIO_InputParameter('Moire.BilayerBottomYMin',ymin1,0.0_dp)
        !! - `Moire.BilayerBottomYMin`: Moire.BLShiftRange: first y shift of the bottom layer.
      call MIO_InputParameter('Moire.BilayerBottomYMax',ymax1,0.0_dp)
        !! - `Moire.BilayerBottomYMax`: Moire.BLShiftRange: last y shift of the bottom layer.
      call MIO_InputParameter('Moire.BilayerTopXMin',xmin2,0.0_dp)
        !! - `Moire.BilayerTopXMin`: Moire.BLShiftRange: first x shift of the top layer.
      call MIO_InputParameter('Moire.BilayerTopXMax',xmax2,0.0_dp)
        !! - `Moire.BilayerTopXMax`: Moire.BLShiftRange: last x shift of the top layer.
      call MIO_InputParameter('Moire.BilayerTopYMin',ymin2,0.0_dp)
        !! - `Moire.BilayerTopYMin`: Moire.BLShiftRange: first y shift of the top layer.
      call MIO_InputParameter('Moire.BilayerTopYMax',ymax2,0.0_dp)
        !! - `Moire.BilayerTopYMax`: Moire.BLShiftRange: last y shift of the top layer.
      tauX1= xmin1 + mS*(xmax1-xmin1)/mSStep
      tauY1= ymin1 + mS*(ymax1-ymin1)/mSStep
      tauX2= xmin2 + mS*(xmax2-xmin2)/mSStep
      tauY2= ymin2 + mS*(ymax2-ymin2)/mSStep
      if (mSi /= mSf) then
         write(prefix,'(a,f0.4,a,f0.4,a,f0.4,a,f0.4)') 'X1_',tauX1,'_Y1_', tauY1,'_X2_',tauX2,'_Y2_',tauY2
         call MIO_Section('Bilayer shift '//trim(num2str(mS))//' of '//trim(num2str(mSi))//'..'//trim(num2str(mSf)),'moire')
         call MIO_Print("Writing files in folder '"//trim(prefix)//"'",'moire')
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
