! Alternative names of input keys, compiled into the input library: one line per key,
!
!    call InputAddAlias('Section.Name','formerName')
!
! The file is included in InputAliasInit (input.F90). It is empty until the names proposed in
! docs/development/input-key-renaming-proposal.md are settled.

! Misspelt key that cell.F90 read until October 2026; inputs that contain it keep working.
call InputAddAlias('basedOnMoireCellParameters','basedOnMoireCellParamters')
