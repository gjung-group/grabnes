module grabnes_start

   implicit none

   PRIVATE

   public :: GrabnesStart

contains

subroutine GrabnesStart(success)

   use mio
   use name,                 only : prefix, sysname
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
      call MIO_InputParameter('Run.Description',desc,'')
      if (desc/='') call MIO_Print('Description: '//trim(desc),'start')
      call MIO_InputParameter('Run.Prefix',sysname,'grabnes')
      prefix = sysname
      call MIO_Print('Output prefix: '//trim(prefix),'start')
   end if
   !call MIO_InputParameter('Random.Test',rtest,.false.)
   !call MIO_InputParameter('Random.TestNumber',ntest,100)
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

   character(len=10) :: proc

   !mess = 'GraBNES version: '//VersionGet()
   !call MIO_Print(mess)
   !mess = 'Last updated: '//VersionUpdate()
   !call MIO_Print(mess)
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
   call MIO_Section('Run set-up','start')
#ifdef MPI
   write(proc,'(i0)') nProc
   call MIO_Print('MPI processes: '//trim(proc),'start')
#else
   call MIO_Print('Serial build (no MPI)','start')
#endif
   write(proc,'(i0)') numThreads
   call MIO_Print('OpenMP threads: '//trim(proc),'start')

end subroutine StartHeader
!****** End subroutine: StartHeader *******************************************
!******************************************************************************


end module grabnes_start
