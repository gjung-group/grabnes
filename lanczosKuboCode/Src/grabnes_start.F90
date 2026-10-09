module grabnes_start

   implicit none

   PRIVATE

   public :: GrabnesStart

contains

subroutine GrabnesStart(success)

   use mio
   use name
   !use random,               only : RandTest, RandSeed, rand_t

   logical, intent(out) :: success
   !type(rand_t), save :: rng
   logical :: rtest
   integer :: ntest
   character(len=81) :: desc
   integer :: rr

   call MIO_Initialize('grabnes',success)
   if (success) then
      call StartHeader()
      call MIO_InputParameter('Description',desc,'')
      if (desc/='') then
         call MIO_Print('---- Description of the system ----','start')
         call MIO_Print(desc,'start')
         call MIO_Print('-------------------------------------------------------'&
           //'-----------','start')
      end if
      call MIO_InputParameter('Prefix',sysname,'grabnes')
      prefix = sysname
   end if
   call MIO_Print('')
   !call MIO_InputParameter('RandomTest',rtest,.false.)
   !call MIO_InputParameter('RandomTestNumber',ntest,100)
   !if (rtest) then
   !!$OMP PARALLEL
   !   call RandSeed(rng,nThread)
   !   call RandTest(rng,ntest)
   !!$OMP END PARALLEL
   !end if

end subroutine GrabnesStart


!****** Subroutine: StartHeader ***********************************************
!******************************************************************************
!
!  Print initial header
!
!******************************************************************************
subroutine StartHeader()

   use mio
   !use version,              only : VersionGet, VersionUpdate

   character(len=65) :: mess
   character(len=10) :: proc

   !mess = 'GraBNES version: '//VersionGet()
   !call MIO_Print(mess)
   !mess = 'Last updated: '//VersionUpdate()
   !call MIO_Print(mess)
#ifdef MPI
   write(proc,'(i8)') nProc
   if (nProc > 1) then
      mess = 'Running MPI version on '//trim(adjustl(proc))//&
        ' processors'
   elseif (nProc == 1) then
      mess = 'Running MPI version on 1 processor'
   end if
   call MIO_Print(mess)
#else
   call MIO_Print('Running in serial mode')
#endif
   call MIO_Print('')
   call MIO_Print('    +----------------------------------------------------+')
   call MIO_Print('    | +------------------------------------------------+ |')
   call MIO_Print('    | |                                                | |')
   call MIO_Print('    | |                    GraBNES                     | |')
   call MIO_Print('    | |      Graphene and BN Electronic Structure      | |')
   call MIO_Print('    | |      ------------------------------------      | |')
   call MIO_Print('    | |          Quantum Transport Calculator          | |')
   call MIO_Print('    | |                                                | |')
   call MIO_Print('    | +------------------------------------------------+ |')
   call MIO_Print('    +----------------------------------------------------+')
   call MIO_Print('')

end subroutine StartHeader
!****** End subroutine: StartHeader *******************************************
!******************************************************************************


end module grabnes_start
