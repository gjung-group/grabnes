!---------------------------------------------------------!
!             MEMORY, INPUT & OUTPUT LIBRARY              !
!                                                         !
!     Copyright (C) 2012 Rafael Martinez-Gordillo         !
!                                                         !
!   This file is distributed under the terms of the       !
!   GNU General Public License. See the file `LICENSE'    !
!   in the root directory of the present distribution,    !
!   or http://www.gnu.org/copyleft/gpl.txt                !
!                                                         !
!---------------------------------------------------------!

module mem

   use precision,            only : dp, sp
   use sys,                  only : SysIOErr, SysKill

   implicit none

   PRIVATE

   integer, public, parameter :: Isz = [ISZ]
   integer, public, parameter :: Rsz = [RSZ]
   integer, public, parameter :: Dsz = [DSZ]
   integer, public, parameter :: Csz = [CSZ]
   integer, public, parameter :: Zsz = [ZSZ]
   integer, public, parameter :: Ssz = [SSZ]
   integer, public, parameter :: Lsz = [LSZ]

   integer, save :: totalMem = 0, peakMem=0
   integer :: istat

   public :: MemCount
   public :: MemPrint
   public :: MemGetName

   [WRITE_INTERFACE]

contains

!****** Subroutine: MemCount **************************************************
!******************************************************************************
!
!  Counts the memory in use
!
! INPUT -----------------------------------------------------------------------
!
! integer sz                      : Size in bytes of the array allocated
!                                   (if sz>0) or deallocated (if sz<0)
! character arrname               : Name of the array
! character rname                 : Name of the routine that ask for memory
!
!******************************************************************************
subroutine MemCount(sz,arrname,rname)

    implicit none

    integer, intent(in) :: sz
    character(len=*), intent(in) :: arrname, rname

    totalMem = totalMem + sz
    if (totalMem > peakMem) then
        peakMem = totalMem
    end if

end subroutine MemCount
!****** End subroutine: MemCount **********************************************
!******************************************************************************


!****** Subroutine: MemPrint **************************************************
!******************************************************************************
!
!  Print the memory usage
!
!
!******************************************************************************
subroutine MemPrint()

    use format,              only : num2bytes
    use sys,                 only : SysPrint

    implicit none

    character(len=50) :: mess
    character(len=200) :: line
    character(len=80) :: mess2
    integer :: u, ios
    real(kind(1.0d0)) :: kb

    call SysPrint('')
    mess = 'Peak memory usage: '//trim(num2bytes(peakMem))
    call SysPrint(mess,'mem')
    ! The line above counts only the arrays allocated through this library. The peak resident
    ! memory of the process, as the operating system reports it (Linux), is the number to size jobs with.
    open(newunit=u,file='/proc/self/status',status='old',action='read',iostat=ios)
    if (ios == 0) then
       do
          read(u,'(a)',iostat=ios) line
          if (ios /= 0) exit
          if (line(1:6) == 'VmHWM:') then
             read(line(7:),*,iostat=ios) kb
             if (ios == 0) then
                ! (the lines of this file must stay below 100 columns: the generator of mem.f90 cuts them)
                if (kb >= 1048576.0d0) then
                   write(mess2,'(a,f0.2,a)') 'Peak memory of the process (system): ', &
                      kb/1048576.0d0, ' GB'
                else
                   write(mess2,'(a,f0.1,a)') 'Peak memory of the process (system): ', &
                      kb/1024.0d0, ' MB'
                end if
                write(*,'(a)') 'mem: '//trim(mess2)
             end if
             exit
          end if
       end do
       close(u)
    end if

end subroutine MemPrint
!****** End subroutine: MemPrint **********************************************
!******************************************************************************

   [WRITE_SUBS]
!****** Subroutine: Alloc_is***************************************************
subroutine Alloc_is(array,sz,name,routine,save)

   implicit none

   integer, pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_is
!******************************************************************************
!****** Subroutine: Alloc_rs **************************************************
subroutine Alloc_rs(array,sz,name,routine,save)

   implicit none

   real(sp), pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_rs
!******************************************************************************
!****** Subroutine: Alloc_ds **************************************************
subroutine Alloc_ds(array,sz,name,routine,save)

   implicit none

   real(dp), pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_ds
!******************************************************************************
!****** Subroutine: Alloc_cs **************************************************
subroutine Alloc_cs(array,sz,name,routine,save)

   implicit none

   complex(sp), pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_cs
!******************************************************************************
!****** Subroutine: Alloc_zs **************************************************
subroutine Alloc_zs(array,sz,name,routine,save)

   implicit none

   complex(dp), pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_zs
!******************************************************************************
!****** Subroutine: Alloc_ls **************************************************
subroutine Alloc_ls(array,sz,name,routine,save)

   implicit none

   logical, pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_ls
!******************************************************************************
!****** Subroutine: Alloc_ss **************************************************
subroutine Alloc_ss(array,sz,name,routine,save)

   implicit none

   character(*), pointer, intent(inout) :: array(:)
   integer, intent(in) :: sz
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szu(1), szl(1)

   szu = sz
   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc_ss
!******************************************************************************

!****** Subroutine: MemGetName ************************************************
!******************************************************************************
!
!  Get names for the variable and the routine
!
! INPUT -----------------------------------------------------------------------
!
! character name                 : Name of the variable
! character routine              : Calling routine
!
! OUTPUT ----------------------------------------------------------------------
!
! character arrname              : Final name of the array
! character rname                : Final name of the routine
!
!******************************************************************************
subroutine MemGetName(arrname,rname,name,routine)

   implicit none

   character(len=32), intent(out) :: arrname, rname
   character(len=*), intent(in), optional :: name, routine

   if (present(name)) then
      arrname = name
   else
      arrname = 'unknown'
   end if
   if (present(routine)) then
      rname = routine
   else
      rname = 'unknown'
   end if

end subroutine MemGetName
!****** End subroutine: MemGetName ********************************************
!******************************************************************************


end module mem
