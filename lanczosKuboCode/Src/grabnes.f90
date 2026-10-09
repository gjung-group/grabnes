program grabnes
! testing sphinx

   use grabnes_start,         only : GrabnesStart
   use grabnes_end,           only : GrabnesEnd
   use atoms,                 only : AtomsPos
   use cell,                  only : CellGet
   use calc,                  only : CalcSelect
   use parallel,              only : ParallelDiv

   implicit none

   logical :: cont

   call GrabnesStart(cont)

   if (cont) then
      call CellGet()
      call ParallelDiv()
      call AtomsPos()
      call CalcSelect()
   end if
   call GrabnesEnd()

end program grabnes
