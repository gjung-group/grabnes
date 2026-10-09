module grabnes_end

   implicit none

   PRIVATE

   public :: GrabnesEnd

contains

subroutine GrabnesEnd

   use mio

   call MIO_Finalize()

end subroutine GrabnesEnd

end module grabnes_end
