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


!****** Subroutine: Alloc[_TYPE] **********************************************
subroutine Alloc[_TYPE]_c(array,szl,szu,name,routine,save)

   implicit none

   [TYPE], pointer :: array([RANK_DEC])
   integer, intent(in) :: szl([RANK]), szu([RANK])
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   [TYPETMP], pointer :: tmp([RANK_DEC])
   character(len=32) :: arrname, rname
   logical :: copy, dealloc, alloc
   integer :: arrSz, oldSz, i
   integer :: prevBounds([RANK]*2), cpBnd([RANK]*2), indx([RANK]*2)

   call MemGetName(arrname,rname,name,routine)
   copy = .false.
   indx(1:[RANK]) = szl
   indx([RANK]+1:2*[RANK]) = szu
   if (associated(array)) then
      oldSz = size(array)*[TYPESZ]
      if (present(save)) then
         copy = save
      end if
      prevBounds(1:[RANK]) = lbound(array)
      prevBounds([RANK]+1:2*[RANK]) = ubound(array)
      if (all(prevBounds==indx)) then
         dealloc = .false.
         alloc = .false.
      else
         dealloc = .true.
         alloc = .true.
      end if
   else
      dealloc = .false.
      alloc = .true.
   end if
   if (copy) tmp => array
   if (dealloc .and. .not. copy) then
      call MemDealloc(array,arrname,rname)
   end if
   if (alloc) then
      allocate(array([BOUNDS]),STAT=istat)
      call SysIOErr(istat)
      arrSz = size(array)*[TYPESZ]
      call MemCount(arrSz,arrname,rname)
      array = [INIT]
   end if
   if (copy) then
      do i=1,[RANK]
         cpBnd(i) = max(prevBounds(i),szl(i))
      end do
      do i=1,[RANK]
         cpBnd(i+[RANK]) = min(prevBounds(i+[RANK]),szu(i))
      end do
      array([COPY]) = &
        tmp([COPY])
      if (alloc) then
         call MemDealloc(tmp,'tmp'//trim(arrname),rname)
      end if
   end if

end subroutine Alloc[_TYPE]_c
subroutine Alloc[_TYPE]_s(array,szu,name,routine,save)

   implicit none

   [TYPE], pointer :: array([RANK_DEC])
   integer, intent(in) :: szu([RANK])
   character(len=*), intent(in), optional :: routine, name
   logical, intent(in), optional :: save

   integer :: szl([RANK])

   szl = 1
   call MemAlloc(array,szl,szu,name,routine,save)

end subroutine Alloc[_TYPE]_s
!******************************************************************************
!****** Subroutine: Dealloc[_TYPE] ********************************************
subroutine Dealloc[_TYPE](array,name,routine)

   implicit none

   [TYPE], pointer :: array([RANK_DEC])
   character(len=*), intent(in), optional :: routine, name

   character(len=32) :: arrname, rname
   integer :: arrSz

   call MemGetName(arrname,rname,name,routine)
   if (associated(array)) then
      arrSz = size(array)*[TYPESZ]
      deallocate(array,STAT=istat)
      call SysIOErr(istat)
      call MemCount(-arrSz,arrname,rname)
   end if

end subroutine Dealloc[_TYPE]
!******************************************************************************
