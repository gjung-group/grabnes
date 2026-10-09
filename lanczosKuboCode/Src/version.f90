module version

   implicit none

   PRIVATE

   character(len=*), parameter :: ver = "KUBO_VERSION"
   character(len=*), parameter :: update = "KUBO_UPDATE"

   public :: VersionGet
   public :: VersionUpdate

contains

!****** Function: VersionGet **************************************************
!******************************************************************************
!
!  Get the version number in variable 'ver', updated at compilation time.
!
!******************************************************************************
function VersionGet() result(res)

   implicit none

   character(len=len(ver)) :: res

   res = ver

end function VersionGet
!****** End function: VersionGet **********************************************
!******************************************************************************


!****** Function: VersionUpdate ***********************************************
!******************************************************************************
!
!  Get the date of the last update in variable 'update', updated at compilation
! time.
!
!******************************************************************************
function VersionUpdate() result(res)

   implicit none

   character(len=len(update)) :: res

   res = update

end function VersionUpdate
!****** End function: VersionUpdate *******************************************
!******************************************************************************

end module version
