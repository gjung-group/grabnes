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

module input

   use file,                 only : cl_file
   use io,                   only : maxlinel, maxfnlen
#ifdef _OPENMP
   use omp_lib,              only : omp_in_parallel, omp_get_thread_num
#endif
#ifdef MPI
   use mpi
#endif /* MPI */

   implicit none

   PRIVATE

   type(cl_file) :: InFl
   type(cl_file), save :: InFile, InLns
   integer, save :: nlines, maxrecl, uln, uin
   integer :: iostat

   logical, save :: warned_parallel_once = .false.

   public :: InputInit
   public :: InputClose
   public :: InputParameter
   public :: InputFindBlock
   public :: InputBlock
   public :: InputSearchLabel

   interface InputParameter
      module procedure InputPar_i, InputPar_d, InputPar_r, InputPar_l, InputPar_s
      module procedure InputPar_is, InputPar_ds, InputPar_rs, InputPar_ls
   end interface

   public :: InputParameterVariable
   interface InputParameterVariable
      module procedure InputPar_i_variable
   end interface

   interface InputFindBlock
      module procedure InputFindBl, InputFindBl_i, InputFindBl_s
   end interface

   interface InputBlock
      module procedure InputBl_i, InputBl_r, InputBl_d, InputBl_l, InputBl_s
      module procedure InputBl_im, InputBl_rm, InputBl_dm, InputBl_lm
      module procedure InputBl_imi, InputBl_rmi, InputBl_dmi, InputBl_lmi
   end interface

contains


!****** Subroutine: InputInit *************************************************
!******************************************************************************
!
!  Initialize the input information
!
!******************************************************************************
subroutine InputInit(success)

   use options,              only : OptionsInputName
   use io,                   only : maxfnlen
   use sys,                  only : SysWarning
   use mem,                  only : Ssz, Isz
#ifdef DEBUG
   use debug,                only : DebugPrint
#endif /* DEBUG */
#ifdef MPI
   use mpi
#endif /* MPI */
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   logical, intent(out) :: success

   integer :: maxl, uif
   character(len=maxfnlen) :: inputfile

#ifdef DEBUG
   call DebugPrint('MIO:InputInit',0)
#endif /* DEBUG */
#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */
#ifdef _OPENMP
   if (omp_in_parallel()) then
      if (.not. warned_parallel_once) then
         warned_parallel_once = .true.
         write(*, '(A,I0,A)') 'input: WARNING: input read inside OpenMP parallel region (thread=', omp_get_thread_num(), ')'
      end if
   end if
#endif
#ifdef _OPENMP
   if (omp_in_parallel()) then
      if (.not. warned_parallel_once) then
         warned_parallel_once = .true.
         write(*, '(A,I0,A)') 'input: WARNING: input read inside OpenMP parallel region (thread=', omp_get_thread_num(), ')'
      end if
   end if
#endif
#ifdef _OPENMP
   if (omp_in_parallel()) then
      if (.not. warned_parallel_once) then
         warned_parallel_once = .true.
         write(*, '(A,I0,A)') 'input: WARNING: input read inside OpenMP parallel region (thread=', omp_get_thread_num(), ')'
      end if
   end if
#endif

   success = .true.
   inputfile = OptionsInputName()
   if (inputfile=='') then
      success = .false.
      call SysWarning('No input file specified','input','InputInit')
      return
   end if
   call InFl%Open(name=inputfile,action='READ',status='OLD',serial=.true.)
   call InFile%Open(action='WRITE',status='REPLACE',serial=.true.)
   call InLns%Open(action='WRITE',status='REPLACE',form='UNFORMATTED', &
     recl=maxfnlen*Ssz+Isz,serial=.true.)
   uif = InFl%GetUnit()
   uin = InFile%GetUnit()
   uln = InLns%GetUnit()
   maxrecl = 0
   call InputProcess(uif,uin,uln,inputfile,maxrecl,nlines)
#ifdef MPI
   call MPIBCast(nlines,1,MPI_INTEGER)
   call MPIBCast(maxrecl,1,MPI_INTEGER)
#endif /* MPI */
   call InFl%Close()
   call InFile%Close(.false.)
   call InLns%Close(.false.)
   call InFile%SetDelete(.true.)
   call InLns%SetDelete(.true.)
   call InFile%Open(action='READ',status='OLD',serial=.false.)
   call InLns%Open(action='READ',status='OLD',form='UNFORMATTED',serial=.true.)
#ifdef MPI
   maxrecl = maxrecl+1
   uin = InFile%GetFH()
#else
   uin = InFile%GetUnit()
#endif /* MPI */
   uln = InLns%GetUnit()
   call InputCheckDuplicates()

#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */
#ifdef DEBUG
   call DebugPrint('MIO:InputInit',1)
#endif /* DEBUG */

end subroutine InputInit
!****** End subroutine: InputInit *********************************************
!******************************************************************************


!****** Subroutine: InputCheckDuplicates **************************************
!******************************************************************************
!
!  Warn about keys that appear more than once in the input file. A key is
!  searched from the top of the file and its first occurrence is used, so a
!  second value given further down is silently ignored. The contents of
!  blocks (&begin ... &end) are data, not keys.
!
!******************************************************************************
subroutine InputCheckDuplicates()

   use string,               only : StringComp
   use sys,                  only : SysIOErr
   use mem,                  only : Ssz
#ifdef MPI
   use sys,                  only : SysIOErrMPI
#endif /* MPI */

   implicit none

   integer :: il, i, j, n, l, nrep
#ifdef MPI
   integer(kind=MPI_OFFSET_KIND) :: nbytes
#endif /* MPI */
   character(len=maxrecl) :: line
   character(len=maxrecl), allocatable :: keys(:)
   logical, allocatable :: done(:)
   logical :: inblock

   if (nlines < 2) return
   allocate(keys(nlines), done(nlines))
   n = 0
   inblock = .false.
#ifdef MPI
   nbytes = maxrecl*Ssz
#else
   rewind(uin)
#endif /* MPI */
   do il=0,nlines-1
#ifdef MPI
      call MPI_File_Read_At_All(uin,il*nbytes,line,maxrecl,MPI_CHARACTER,MPI_STATUS_IGNORE,iostat)
      call SysIOErrMPI(iostat,'input','InputCheckDuplicates')
      if (ichar(line(maxrecl:maxrecl))==10) then
         line(maxrecl:maxrecl) = ''
      end if
#else
      read(uin,'(a)',IOSTAT=iostat) line
      call SysIOErr(iostat,'input','InputCheckDuplicates')
#endif /* MPI */
      line = adjustl(line)
      if (len_trim(line)==0) cycle
      if (line(1:1)=='#' .or. line(1:1)=='!') cycle
      if (line(1:1)=='&') then
         inblock = StringComp(line(1:6),'&begin')
         cycle
      end if
      if (inblock) cycle
      l = scan(line,' =:')
      if (l > 1) line = line(1:l-1)
      n = n + 1
      keys(n) = line
   end do
   done(1:n) = .false.
   do i=1,n-1
      if (done(i)) cycle
      nrep = 1
      do j=i+1,n
         if (done(j)) cycle
         if (len_trim(keys(j)) /= len_trim(keys(i))) cycle
         if (StringComp(keys(j),trim(keys(i)))) then
            nrep = nrep + 1
            done(j) = .true.
         end if
      end do
      if (nrep > 1) then
         write(*,'(A,I0,A)') 'input: WARNING: the key "'//trim(keys(i))//'" appears ', nrep, &
           ' times in the input file; only its FIRST value is used.'
      end if
   end do
   deallocate(keys, done)

end subroutine InputCheckDuplicates
!****** End subroutine: InputCheckDuplicates **********************************
!******************************************************************************


!****** Subroutine: InputClose ************************************************
!******************************************************************************
!
!  Close input module
!
!******************************************************************************
subroutine InputClose()

#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */
   call InFile%Close()
   call InLns%Close()

#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end subroutine InputClose
!****** End subroutine: InputClose ********************************************
!******************************************************************************


!****** Subroutine: InputParameter ********************************************
!******************************************************************************
!
!  Returns the parameter with the specified label from input file
!
! INPUT -----------------------------------------------------------------------
!
! character label                : Label of the parameter
! * def(:)                       : Default value
! character(*) defstr            : Default string if it's required
!
! OUTPUT ----------------------------------------------------------------------
!
! * value                        : Value for the parameter
! character(*) str               : Optional string found after the values
!
!******************************************************************************
subroutine InputPar_i(label,value,def,str,defstr)

   use sys,                  only : SysKill
   use parser,               only : ParserCheck
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   integer, intent(out) :: value(:)
   integer, intent(in) :: def(:)
   character(*), intent(out), optional :: str
   character(*), intent(in), optional :: defstr

   integer :: sz, i, ierr, nl, id
   character(len=maxrecl) :: line

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */
#ifdef _OPENMP
   if (omp_in_parallel()) then
      if (.not. warned_parallel_once) then
         warned_parallel_once = .true.
         write(*, '(A,I0,A)') 'input: WARNING: input read inside OpenMP parallel region (thread=', omp_get_thread_num(), ')'
      end if
   end if
#endif
   sz = size(value)
   if (sz /= size(def)) call SysKill('Incorrect size of default value'// &
     " for parameter '"//trim(label)//"'",'input','InputParameter')
   if (InputSearchLabel(label, line, id)) then
      do i=1,sz
         call ParserCheck(value(i),line,i,ierr)
         call InputParErr(label,ierr,id)
      end do
   else
      value = def
   end if
   if (present(str)) then
      if (present(defstr)) then
         if (line=='') then
            str = defstr
         else
            call ParserCheck(str,line,sz+1,ierr)
         end if
      else
         call SysKill('Default value for string must be present for parameter '&
           //trim(label),'input','InputParameter')
      end if
   end if
#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end subroutine InputPar_i
!******************************************************************************
!> @brief Read variable-length integer array from input (stops when elements run out)
!! @param[in] label Input parameter label
!! @param[inout] value Array to fill (will be filled from start, elements beyond count unchanged)
!! @param[out] count Number of elements actually read
!! @param[in] maxCount Maximum number of elements to read (size of value array)
!! @details Reads elements until end-of-record (error -1) or maxCount reached
subroutine InputPar_i_variable(label,value,count,maxCount)

   use parser,               only : ParserCheck
   use sys,                  only : SysKill
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   integer, intent(out) :: value(:)
   integer, intent(out) :: count
   integer, intent(in) :: maxCount

   integer :: i, ierr, id
   integer :: sz
   character(len=maxrecl) :: line

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */

#ifdef _OPENMP
   if (omp_in_parallel()) then
      if (.not. warned_parallel_once) then
         warned_parallel_once = .true.
         write(*, '(A,I0,A)') 'input: WARNING: input read inside OpenMP parallel region (thread=', omp_get_thread_num(), ')'
      end if
   end if
#endif
   sz = size(value)
   if (maxCount > sz) then
      call TimerStop('input')
      call SysKill('maxCount exceeds array size for parameter '//trim(label),'input','InputParameterVariable')
   end if

   count = 0
   if (InputSearchLabel(label, line, id)) then
      ! Read elements one by one until we run out (error -1 means element not found)
      do i = 1, maxCount
         call ParserCheck(value(i), line, i, ierr)
         if (ierr == 0) then
            ! Successfully read element i
            count = count + 1
         else if (ierr == -1) then
            ! No more elements found - this is normal, stop reading
            exit
         else
            ! Some other error occurred - call InputParErr and stop timer before SysKill
            call TimerStop('input')
            call InputParErr(label, ierr, id)
            return
         end if
      end do
   else
      ! Parameter not found - count stays 0, values unchanged
      count = 0
   end if
#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end subroutine InputPar_i_variable
!******************************************************************************
subroutine InputPar_is(label,value,def,str,defstr)

   implicit none

   character(*), intent(in) :: label
   integer, intent(out) :: value
   integer, intent(in) :: def
   character(*), intent(out), optional :: str
   character(*), intent(in), optional :: defstr

   integer :: v(1)

   call InputParameter(label,v,(/def/),str,defstr)
   value = v(1)

end subroutine InputPar_is
!******************************************************************************
subroutine InputPar_d(label,value,def,str,defstr)

   use sys,                  only : SysKill
   use parser,               only : ParserCheck
   use precision,            only : dp
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   real(dp), intent(out) :: value(:)
   real(dp), intent(in) :: def(:)
   character(*), intent(out), optional :: str
   character(*), intent(in), optional :: defstr

   integer :: sz, i, ierr, nl, id
   character(len=maxrecl) :: line

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */
#ifdef _OPENMP
   if (omp_in_parallel()) then
      if (.not. warned_parallel_once) then
         warned_parallel_once = .true.
         write(0, '(A,A,A,I0)') 'input: WARNING: input read inside OpenMP parallel region for label "', trim(label), '" (thread=', omp_get_thread_num(), ')'
      end if
   end if
#endif
   sz = size(value)
   if (sz /= size(def)) call SysKill('Incorrect size of default value'// &
     " for parameter '"//trim(label)//"'",'input','InputParameter')
   if (InputSearchLabel(label, line, id)) then
      do i=1,sz
         call ParserCheck(value(i),line,i,ierr)
         call InputParErr(label,ierr, id)
      end do
   else
      value = def
   end if
   if (present(str)) then
      if (present(defstr)) then
         if (line=='') then
            str = defstr
         else
            call ParserCheck(str,line,sz+1,ierr)
         end if
      else
         call SysKill('Default value for string must be present for parameter '&
           //trim(label),'input','InputParameter')
      end if
   end if
#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end subroutine InputPar_d
!******************************************************************************
subroutine InputPar_ds(label,value,def,str,defstr)

   use precision,            only : dp

   implicit none

   character(*), intent(in) :: label
   real(dp), intent(out) :: value
   real(dp), intent(in) :: def
   character(*), intent(out), optional :: str
   character(*), intent(in), optional :: defstr

   real(dp) :: v(1)

   call InputParameter(label,v,(/def/),str,defstr)
   value = v(1)

end subroutine InputPar_ds
!******************************************************************************
subroutine InputPar_r(label,value,def,str,defstr)

   use sys,                  only : SysKill
   use parser,               only : ParserCheck
   use precision,            only : sp
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   real(sp), intent(out) :: value(:)
   real(sp), intent(in) :: def(:)
   character(*), intent(out), optional :: str
   character(*), intent(in), optional :: defstr

   integer :: sz, i, ierr, nl, id
   character(len=maxrecl) :: line

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */
   sz = size(value)
   if (sz /= size(def)) call SysKill('Incorrect size of default value'// &
     " for parameter '"//trim(label)//"'",'input','InputParameter')
   if (InputSearchLabel(label, line, id)) then
      do i=1,sz
         call ParserCheck(value(i),line,i,ierr)
         call InputParErr(label,ierr, id)
      end do
   else
      value = def
   end if
   if (present(str)) then
      if (present(defstr)) then
         if (line=='') then
            str = defstr
         else
            call ParserCheck(str,line,sz+1,ierr)
         end if
      else
         call SysKill('Default value for string must be present for parameter '&
           //trim(label),'input','InputParameter')
      end if
   end if
#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end subroutine InputPar_r
!******************************************************************************
subroutine InputPar_rs(label,value,def,str,defstr)

   use precision,            only : sp

   implicit none

   character(*), intent(in) :: label
   real(sp), intent(out) :: value
   real(sp), intent(in) :: def
   character(*), intent(out), optional :: str
   character(*), intent(in), optional :: defstr

   real(sp) :: v(1)

   call InputParameter(label,v,(/def/),str,defstr)
   value = v(1)

end subroutine InputPar_rs
!******************************************************************************
subroutine InputPar_l(label,value,def,str,defstr)

   use sys,                  only : SysKill
   use parser,               only : ParserCheck
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   logical, intent(out) :: value(:)
   logical, intent(in) :: def(:)
   character(*), intent(out), optional :: str
   character(*), intent(in), optional :: defstr

   integer :: sz, i, ierr, nl, id
   character(len=maxrecl) :: line

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */
   sz = size(value)
   if (sz /= size(def)) call SysKill('Incorrect size of default value'// &
     " for parameter '"//trim(label)//"'",'input','InputParameter')
   if (InputSearchLabel(label, line, id)) then
      do i=1,sz
         call ParserCheck(value(i),line,i,ierr)
         call InputParErr(label,ierr,id)
      end do
   else
      value = def
   end if
   if (present(str)) then
      if (present(defstr)) then
         if (line=='') then
            str = defstr
         else
            call ParserCheck(str,line,sz+1,ierr)
         end if
      else
         call SysKill('Default value for string must be present for parameter '&
           //trim(label),'input','InputParameter')
      end if
   end if
#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end subroutine InputPar_l
!******************************************************************************
subroutine InputPar_ls(label,value,def,str,defstr)

   implicit none

   character(*), intent(in) :: label
   logical, intent(out) :: value
   logical, intent(in) :: def
   character(*), intent(out), optional :: str
   character(*), intent(in), optional :: defstr

   logical :: v(1)

   call InputParameter(label,v,(/def/),str,defstr)
   value = v(1)

end subroutine InputPar_ls
!******************************************************************************
subroutine InputPar_s(label,value,def)

   use sys,                  only : SysKill
   use parser,               only : ParserCheck
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   character(*), intent(out) :: value
   character(*), intent(in) :: def

   integer :: sz, i, ierr, nl, id
   character(len=maxrecl) :: line

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */
   if (InputSearchLabel(label, line, id)) then
      value = line
   else
      value = def
   end if
#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end subroutine InputPar_s
!****** End subroutine: InputParameter ****************************************
!******************************************************************************

!****** Function: InputFindBlock **********************************************
!******************************************************************************
!
!  Finds a block
!
! INPUT -----------------------------------------------------------------------
!
! character label                : Label of the block
!
! OUTPUT ----------------------------------------------------------------------
!
! logical found                  : True if the block is found
!
!******************************************************************************
function InputFindBl(label) result(found)

   use sys,                  only : SysKill
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   logical :: found

   character(len=maxrecl) :: line
   integer :: blockid, endid
   character(maxfnlen) :: fname
   integer :: ln
   character(6) :: num

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */

   if (InputSearchLabel('&end '//label, line, endid)) then
      if (InputSearchLabel('&begin '//label, line, blockid)) then
         if (blockid < endid) then
            found = .true.
         else
            call InputGetLineNum(blockid, fname, ln)
            write(num,'(i0)') ln
            call SysKill('Block "'//trim(label)//'" not properly constructed'//&
              ' in input file '//trim(fname)//' at line '//trim(num), &
              'input','InputFindBlock')
         end if
      else
         call InputGetLineNum(endid, fname, ln)
         write(num,'(i0)') ln
         call SysKill('Found end of block "'//trim(label)//'"  at line '// &
           trim(num)//' of file '//trim(fname)//', but not the beginning', &
           'input','InputFindBlock')
      end if
   else
      found = .false.
   end if

#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end function InputFindBl
!******************************************************************************
function InputFindBl_i(label,value) result(found)

   use sys,                  only : SysKill
   use parser,               only : ParserCheck
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   integer, intent(out) :: value
   logical :: found

   character(len=maxrecl) :: line
   integer :: blockid, endid, ierr
   character(maxfnlen) :: fname
   integer :: ln
   character(6) :: num

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */

   if (InputSearchLabel('&end '//label, line, endid)) then
      if (InputSearchLabel('&begin '//label, line, blockid)) then
         if (blockid < endid) then
            if (line=='') then
               call InputGetLineNum(blockid, fname, ln)
               write(num,'(i0)') ln
               call SysKill('An integer must be specified for block "'// &
                 trim(label)//'". Wrong block format in input file '// &
                 trim(fname)//' at line '//trim(num), 'input', &
                 'InputFindBlock')
            end if
            call ParserCheck(value,line,1,ierr)
            call InputBlockErr(label,ierr,blockid)
            found = .true.
         else
            call InputGetLineNum(blockid, fname, ln)
            write(num,'(i0)') ln
            call SysKill('Block "'//trim(label)//'" not properly constructed'//&
              ' in input file '//trim(fname)//' at line '//trim(num), &
              'input','InputFindBlock')
         end if
      else
         call InputGetLineNum(endid, fname, ln)
         write(num,'(i0)') ln
         call SysKill('Found end of block "'//trim(label)//'"  at line '// &
           trim(num)//' of file '//trim(fname)//', but not the beginning', &
           'input','InputFindBlock')
      end if
   else
      found = .false.
   end if

#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end function InputFindBl_i
!******************************************************************************
function InputFindBl_s(label,value) result(found)

   use sys,                  only : SysKill
   use parser,               only : ParserCheck
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   character(*), intent(out) :: value
   logical :: found

   character(len=maxrecl) :: line
   integer :: blockid, endid, ierr
   character(maxfnlen) :: fname
   integer :: ln
   character(6) :: num

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */

   if (InputSearchLabel('&end '//label, line, endid)) then
      if (InputSearchLabel('&begin '//label, line, blockid)) then
         if (blockid < endid) then
            if (line=='') then
               call InputGetLineNum(blockid, fname, ln)
               write(num,'(i0)') ln
               call SysKill('A string must be specified for block "'// &
                 trim(label)//'". Wrong block format in input file '// &
                 trim(fname)//' at line '//trim(num), 'input', &
                 'InputFindBlock')
            end if
            value = line
            found = .true.
         else
            call InputGetLineNum(endid, fname, ln)
            write(num,'(i0)') ln
            call SysKill('Block "'//trim(label)//'" not properly constructed'//&
              ' in input file '//trim(fname)//', found end of block at line ' &
              //trim(num), 'input','InputFindBlock')
         end if
      else
         call InputGetLineNum(endid, fname, ln)
         write(num,'(i0)') ln
         call SysKill('Found end of block "'//trim(label)//'"  at line '// &
           trim(num)//' of file '//trim(fname)//', but not the beginning', &
           'input','InputFindBlock')
      end if
   else
      found = .false.
   end if

#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end function InputFindBl_s
!****** End function: InputFindBlock ******************************************
!******************************************************************************

!****** Subroutine: InputBlock ************************************************
!******************************************************************************
!
!  Get elements from a block
!
! INPUT -----------------------------------------------------------------------
!
! character label                : Label of the block
! integer column                 : Column number
!
! OUTPUT ----------------------------------------------------------------------
!
! * value                        : Elements obtained from the block
!
!******************************************************************************
subroutine InputBl_i(label,column,value)

   use sys,                  only : SysKill, SysWarning
   use parser,               only : ParserCheck
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   integer, intent(in) :: column
   integer, intent(out) :: value(:)

   logical :: found
   character(len=maxrecl) :: line
   integer :: blockid, endid
   character(maxfnlen) :: fname
   integer :: ln, sz, ierr
   character(6) :: num

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */

   if (InputSearchLabel('&end '//label, line, endid)) then
      if (InputSearchLabel('&begin '//label, line, blockid)) then
         if (blockid < endid) then
            found = .true.
         else
            call InputGetLineNum(endid, fname, ln)
            write(num,'(i0)') ln
            call SysKill('Block "'//trim(label)//'" not properly constructed'//&
              ' in input file '//trim(fname)//', found end of block at line ' &
              //trim(num), 'input','InputFindBlock')
         end if
      else
         call InputGetLineNum(endid, fname, ln)
         write(num,'(i0)') ln
         call SysKill('Found end of block "'//trim(label)//'"  at line '// &
           trim(num)//' of file '//trim(fname)//', but not the beginning', &
           'input','InputFindBlock')
      end if
   else
      found = .false.
   end if
   if (.not. found) then
      call SysWarning('Block "'//trim(label)//'" not found. This block'// &
        ' will not be read', 'input', 'InputBlock')
      return
   end if
   sz = size(value)
   if (sz/=endid-blockid-1) then
      write(fname,'(a,i0,a,i0)') 'Expected size of block "'//trim(label)//'" is ', sz, &
        ', found ', endid-blockid-1
      call SysKill(fname,'input','InputBlock')
   end if
   do ln=1,sz
      call InputGetLine(blockid+ln,line)
      call ParserCheck(value(ln),line,column,ierr)
      call InputBlockErr(label,ierr,blockid+ln)
   end do

#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end subroutine InputBl_i
!******************************************************************************
subroutine InputBl_im(label,value)

   implicit none

   character(*), intent(in) :: label
   integer, intent(out) :: value(:,:)

   integer :: i, sz

   sz = size(value,1)
   do i=1,sz
      call InputBlock(label,i,value(i,:))
   end do

end subroutine InputBl_im
!******************************************************************************
subroutine InputBl_imi(label,column,value)

   implicit none

   character(*), intent(in) :: label
   integer, intent(out) :: value(:,:)
   integer, intent(in) :: column

   integer :: i, sz

   sz = size(value,1)
   do i=1,sz
      call InputBlock(label,i+column-1,value(i,:))
   end do

end subroutine InputBl_imi
!******************************************************************************
subroutine InputBl_r(label,column,value)

   use sys,                  only : SysKill, SysWarning
   use parser,               only : ParserCheck
   use precision,            only : sp
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   integer, intent(in) :: column
   real(sp), intent(out) :: value(:)

   logical :: found
   character(len=maxrecl) :: line
   integer :: blockid, endid
   character(maxfnlen) :: fname
   integer :: ln, sz, ierr
   character(6) :: num

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */

   if (InputSearchLabel('&end '//label, line, endid)) then
      if (InputSearchLabel('&begin '//label, line, blockid)) then
         if (blockid < endid) then
            found = .true.
         else
            call InputGetLineNum(endid, fname, ln)
            write(num,'(i0)') ln
            call SysKill('Block "'//trim(label)//'" not properly constructed'//&
              ' in input file '//trim(fname)//', found end of block at line ' &
              //trim(num), 'input','InputFindBlock')
         end if
      else
         call InputGetLineNum(endid, fname, ln)
         write(num,'(i0)') ln
         call SysKill('Found end of block "'//trim(label)//'"  at line '// &
           trim(num)//' of file '//trim(fname)//', but not the beginning', &
           'input','InputFindBlock')
      end if
   else
      found = .false.
   end if
   if (.not. found) then
      call SysWarning('Block "'//trim(label)//'" not found. This block'// &
        ' will not be read', 'input', 'InputBlock')
      return
   end if
   sz = size(value)
   if (sz/=endid-blockid-1) then
      write(fname,'(a,i0,a,i0)') 'Expected size of block "'//trim(label)//'" is ', sz, &
        ', found ', endid-blockid-1
      call SysKill(fname,'input','InputBlock')
   end if
   do ln=1,sz
      call InputGetLine(blockid+ln,line)
      call ParserCheck(value(ln),line,column,ierr)
      call InputBlockErr(label,ierr,blockid+ln)
   end do

#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end subroutine InputBl_r
!******************************************************************************
subroutine InputBl_rm(label,value)

   use precision,            only : sp

   implicit none

   character(*), intent(in) :: label
   real(sp), intent(out) :: value(:,:)

   integer :: i, sz

   sz = size(value,1)
   do i=1,sz
      call InputBlock(label,i,value(i,:))
   end do

end subroutine InputBl_rm
!******************************************************************************
subroutine InputBl_rmi(label,column,value)

   use precision,            only : sp

   implicit none

   character(*), intent(in) :: label
   real(sp), intent(out) :: value(:,:)
   integer, intent(in) :: column

   integer :: i, sz

   sz = size(value,1)
   do i=1,sz
      call InputBlock(label,i+column-1,value(i,:))
   end do

end subroutine InputBl_rmi
!******************************************************************************
subroutine InputBl_d(label,column,value)

   use sys,                  only : SysKill, SysWarning
   use parser,               only : ParserCheck
   use precision,            only : dp
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   integer, intent(in) :: column
   real(dp), intent(out) :: value(:)

   logical :: found
   character(len=maxrecl) :: line
   integer :: blockid, endid
   character(maxfnlen) :: fname
   integer :: ln, sz, ierr
   character(6) :: num

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */

   if (InputSearchLabel('&end '//label, line, endid)) then
      if (InputSearchLabel('&begin '//label, line, blockid)) then
         if (blockid < endid) then
            found = .true.
         else
            call InputGetLineNum(endid, fname, ln)
            write(num,'(i0)') ln
            call SysKill('Block "'//trim(label)//'" not properly constructed'//&
              ' in input file '//trim(fname)//', found end of block at line ' &
              //trim(num), 'input','InputFindBlock')
         end if
      else
         call InputGetLineNum(endid, fname, ln)
         write(num,'(i0)') ln
         call SysKill('Found end of block "'//trim(label)//'"  at line '// &
           trim(num)//' of file '//trim(fname)//', but not the beginning', &
           'input','InputFindBlock')
      end if
   else
      found = .false.
   end if
   if (.not. found) then
      call SysWarning('Block "'//trim(label)//'" not found. This block'// &
        ' will not be read', 'input', 'InputBlock')
      return
   end if
   sz = size(value)
   call InputGetLineNum(blockid, fname, ln)
   call InputGetLineNum(endid, fname, ln)
   if (sz/=endid-blockid-1) then
      write(fname,'(a,i0,a,i0)') 'Expected size of block "'//trim(label)//'" is ', sz, &
        ', found ', endid-blockid-1
      call SysKill(fname,'input','InputBlock')
   end if
   do ln=1,sz
      call InputGetLine(blockid+ln,line)
      call ParserCheck(value(ln),line,column,ierr)
      call InputBlockErr(label,ierr,blockid+ln)
   end do

#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end subroutine InputBl_d
!******************************************************************************
subroutine InputBl_dm(label,value)

   use precision,            only : dp

   implicit none

   character(*), intent(in) :: label
   real(dp), intent(out) :: value(:,:)

   integer :: i, sz

   sz = size(value,1)
   do i=1,sz
      call InputBlock(label,i,value(i,:))
   end do

end subroutine InputBl_dm
!******************************************************************************
subroutine InputBl_dmi(label,column,value)

   use precision,            only : dp

   implicit none

   character(*), intent(in) :: label
   real(dp), intent(out) :: value(:,:)
   integer, intent(in) :: column

   integer :: i, sz

   sz = size(value,1)
   do i=1,sz
      call InputBlock(label,i+column-1,value(i,:))
   end do

end subroutine InputBl_dmi
!******************************************************************************
subroutine InputBl_l(label,column,value)

   use sys,                  only : SysKill, SysWarning
   use parser,               only : ParserCheck
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   integer, intent(in) :: column
   logical, intent(out) :: value(:)

   logical :: found
   character(len=maxrecl) :: line
   integer :: blockid, endid
   character(maxfnlen) :: fname
   integer :: ln, sz, ierr
   character(6) :: num

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */

   if (InputSearchLabel('&end '//label, line, endid)) then
      if (InputSearchLabel('&begin '//label, line, blockid)) then
         if (blockid < endid) then
            found = .true.
         else
            call InputGetLineNum(endid, fname, ln)
            write(num,'(i0)') ln
            call SysKill('Block "'//trim(label)//'" not properly constructed'//&
              ' in input file '//trim(fname)//', found end of block at line ' &
              //trim(num), 'input','InputFindBlock')
         end if
      else
         call InputGetLineNum(endid, fname, ln)
         write(num,'(i0)') ln
         call SysKill('Found end of block "'//trim(label)//'"  at line '// &
           trim(num)//' of file '//trim(fname)//', but not the beginning', &
           'input','InputFindBlock')
      end if
   else
      found = .false.
   end if
   if (.not. found) then
      call SysWarning('Block "'//trim(label)//'" not found. This block'// &
        ' will not be read', 'input', 'InputBlock')
      return
   end if
   sz = size(value)
   if (sz/=endid-blockid-1) then
      write(fname,'(a,i0,a,i0)') 'Expected size of block "'//trim(label)//'" is ', sz, &
        ', found ', endid-blockid-1
      call SysKill(fname,'input','InputBlock')
   end if
   do ln=1,sz
      call InputGetLine(blockid+ln,line)
      call ParserCheck(value(ln),line,column,ierr)
      call InputBlockErr(label,ierr,blockid+ln)
   end do

#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end subroutine InputBl_l
!******************************************************************************
subroutine InputBl_lm(label,value)

   implicit none

   character(*), intent(in) :: label
   logical, intent(out) :: value(:,:)

   integer :: i, sz

   sz = size(value,1)
   do i=1,sz
      call InputBlock(label,i,value(i,:))
   end do

end subroutine InputBl_lm
!******************************************************************************
subroutine InputBl_lmi(label,column,value)

   implicit none

   character(*), intent(in) :: label
   logical, intent(out) :: value(:,:)
   integer, intent(in) :: column

   integer :: i, sz

   sz = size(value,1)
   do i=1,sz
      call InputBlock(label,i+column-1,value(i,:))
   end do

end subroutine InputBl_lmi
!******************************************************************************
subroutine InputBl_s(label,column,value,multiword)

   use sys,                  only : SysKill, SysWarning
   use parser,               only : ParserCheck
#ifdef TIMER
   use timer,                only : TimerCount, TimerStop
#endif /* TIMER */

   implicit none

   character(*), intent(in) :: label
   integer, intent(in) :: column
   character(len=*), intent(out) :: value(:)
   logical, intent(in), optional :: multiword

   logical :: found
   character(len=maxrecl) :: line
   integer :: blockid, endid
   character(maxfnlen) :: fname
   integer :: ln, sz, ierr
   character(6) :: num

#ifdef TIMER
   call TimerCount('input')
#endif /* TIMER */

   if (InputSearchLabel('&end '//label, line, endid)) then
      if (InputSearchLabel('&begin '//label, line, blockid)) then
         if (blockid < endid) then
            found = .true.
         else
            call InputGetLineNum(endid, fname, ln)
            write(num,'(i0)') ln
            call SysKill('Block "'//trim(label)//'" not properly constructed'//&
              ' in input file '//trim(fname)//', found end of block at line ' &
              //trim(num), 'input','InputFindBlock')
         end if
      else
         call InputGetLineNum(endid, fname, ln)
         write(num,'(i0)') ln
         call SysKill('Found end of block "'//trim(label)//'"  at line '// &
           trim(num)//' of file '//trim(fname)//', but not the beginning', &
           'input','InputFindBlock')
      end if
   else
      found = .false.
   end if
   if (.not. found) then
      call SysWarning('Block "'//trim(label)//'" not found. This block'// &
        ' will not be read', 'input', 'InputBlock')
      return
   end if
   sz = size(value)
   if (sz/=endid-blockid-1) then
      write(fname,'(a,i0,a,i0)') 'Expected size of block "'//trim(label)//'" is ', sz, &
        ', found ', endid-blockid-1
      call SysKill(fname,'input','InputBlock')
   end if
   do ln=1,sz
      call InputGetLine(blockid+ln,line)
      call ParserCheck(value(ln),line,column,ierr,multiword)
      call InputBlockErr(label,ierr,blockid+ln)
   end do

#ifdef TIMER
   call TimerStop('input')
#endif /* TIMER */

end subroutine InputBl_s
!****** End subroutine: InputBlock ********************************************
!******************************************************************************


!****** Function: InputSearchLabel ********************************************
!******************************************************************************
!
!  Looks for label in input file, copy the line (without the label) in string.
! The function returns true if the label is found, false otherwise.
!
! INPUT -----------------------------------------------------------------------
!
! character label                : Label to look for
!
! OUTPUT ----------------------------------------------------------------------
!
! character str                  : The text found in the line with the
!                                  specified label
! integer lineid                 : ID of the line (not the actual line number)
! logical found                  : True if the label is found, false if not
!
!******************************************************************************
function InputSearchLabel(label,str,lineid) result(found)

   use string,               only : StringComp
   use sys,                  only : SysIOErr
   use mem,                  only : Ssz
#ifdef MPI
   use sys,                  only : SysIOErrMPI
#endif /* MPI */

   implicit none

   character(*), intent(in) :: label
   character(*), intent(out), optional :: str
   integer, intent(out), optional :: lineid
   logical :: found

   integer :: l, il
#ifdef MPI
   integer(kind=MPI_OFFSET_KIND) :: nbytes
#endif /* MPI */
   character(len=maxrecl) :: line
   logical :: equal, id, wr

   found = .false.
   if (present(lineid)) then
      id = .true.
      lineid = 0
   else
      id = .false.
   end if
   if (present(str)) then
      wr = .true.
      str = ''
   else
      wr = .false.
   end if
#ifdef MPI
   nbytes = maxrecl*Ssz
#else
   rewind(uin)
#endif /* MPI */
   do il=0,nlines-1
#ifdef MPI
      call MPI_File_Read_At_All(uin,il*nbytes,line,maxrecl,MPI_CHARACTER,MPI_STATUS_IGNORE,iostat)
      call SysIOErrMPI(iostat,'input','InputSearchLabel')
      if (ichar(line(maxrecl:maxrecl))==10) then
         line(maxrecl:maxrecl) = ''
      end if
#else
      read(uin,'(a)',IOSTAT=iostat) line
      call SysIOErr(iostat,'input','InputSearchLabel')
#endif /* MPI */
      equal = StringComp(line,label)
      if (.not. equal) cycle
      if (wr) then
         l = len_trim(label)+1
         line = adjustl(line(l:))
         if (line(1:1)=='=' .or. line(1:1)==':') then
            line = adjustl(line(2:))
         end if
         str = line
      end if
      found = .true.
      if (id) lineid = il
      exit
   end do

end function InputSearchLabel
!****** End function: InputSearchLabel ****************************************
!******************************************************************************

!****** Subroutine: InputProcess **********************************************
!******************************************************************************
!
!  Read input file(s) and prepare a single file
!
! INPUT -----------------------------------------------------------------------
!
! integer u                      : Unit of input file
! character fname                : File name
! integer maxl                   : Maximum length of line found
!
! OUTPUT ----------------------------------------------------------------------
!
! integer maxl                   : Maximum length of line found
!
!******************************************************************************
recursive subroutine InputProcess(u,uw,ul,fname,maxl,numlines,wr)

   use file,                 only : cl_file
   use sys,                  only : SysIOErr

   implicit none

   integer, intent(in) :: u, uw, ul
   character(maxfnlen), intent(in) :: fname
   integer, intent(inout) :: maxl, numlines
   logical, intent(in), optional :: wr

   integer ::  nl, l, il
   integer :: pos1, pos2
   integer :: nlinep, un
   logical :: hold, wrline, wrflag
   character(len=maxlinel) :: line, line0
   type(cl_file) :: newFile

   if (present(wr)) then
      wrflag = wr
   else
      wrflag = .true.
   end if
   if (wrflag .and. maxl == 0) then
      call InputProcess(u,uw,ul,fname,maxl,numlines,.false.)
   end if
#ifdef MPI
   if (IsPNode) then
#endif /* MPI */
      nl = 0
      line = ''
      hold = .false.
      do
         wrline = .false.
         read(u,'(a)', IOSTAT=iostat, END=20) line0
         call SysIOErr(iostat,'input','InputProcess')
         nl = nl + 1
         line0 = adjustl(line0)
         if(line0(1:1)/='!'  .and. line0(1:1)/='#' .and. line0(1:1)/='*') then
            if(len_trim(line0) > 0 ) then
               line = trim(line)//' '//line0
               line = adjustl(line)
               pos1 = index(line,'!')
               pos2 = index(line,'#')
               if(pos1>0 .and. pos2>0 ) then
                  line = line(:min(pos1,pos2)-1)
               elseif(pos1/=0 .or. pos2/=0 ) then
                  line = line(:max(pos1,pos2)-1)
               end if
               l = len_trim(line)
               do il=1,l
                  if (ichar(line(il:il))==9) line(il:il) = ' '
                  if (line(il:il) == ' ') line(il+1:) = adjustl(line(il+1:))
               end do
               l = len_trim(line)
               if (.not. hold) nlinep = nl
               if (line(l:l)=='\') then
                  hold = .true.
                  line(l:l)=''
                  cycle
               end if
               if (line(1:1)=='$') then
                  if (line(1:9)=='$INCLUDE ') then
                     pos1 = index(line,' ') + 1
                     pos2 = index(line(pos1:),' ') + pos1 - 2
                     call newFile%Open(name=line(pos1:pos2),status='OLD', &
                       action='READ',serial=.true.)
                     un = newFile%GetUnit()
                     call InputProcess(un,uw,ul,line(pos1:pos2),maxl,numlines,wrflag)
                     call newFile%Close()
                     line = ''
                     hold = .false.
                  end if
                  cycle
               else
                  wrline = .true.
               end if
               if (l > maxl) maxl = l
            else
               if (line/='') wrline = .true.
            end if
         end if
         if (wrline .and. wrflag) then
            write(uw,'(a)',IOSTAT=iostat) line(1:maxl)
            call SysIOErr(iostat,'input','InputInit')
            write(ul,IOSTAT=iostat) fname,nlinep
            call SysIOErr(iostat,'input','InputInit')
            line = ''
            hold = .false.
            numlines = numlines+1
         elseif (wrline) then
            line = ''
            hold = .false.
         end if
      end do
20    continue
      rewind(u,IOSTAT=iostat)
      call SysIOErr(iostat,'input','InputInit')
#ifdef MPI
   end if
#endif /* MPI */

end subroutine InputProcess
!****** End subroutine: InputProcess ******************************************
!******************************************************************************

!****** Subroutine: InputParErr ***********************************************
!******************************************************************************
!
!  Manages the errors passed by the parser
!
! INPUT -----------------------------------------------------------------------
!
! integer err                    : Error code
! integer id                     : ID for the processed line
!
!******************************************************************************
subroutine InputParErr(label,err,id)

   use sys,                  only : SysKill, SysWarning, SysPrint, SysIOErr

   implicit none

   character(*), intent(in) :: label
   integer, intent(in) :: err, id

   character(maxfnlen) :: fname
   integer :: ln
   character(6) :: num

   if (err/=0) then
      call InputGetLineNum(id, fname, ln)
      write(num,'(i0)') ln
      if (err==-2) then
         call SysWarning('Taking only integer part for tag "'//trim(label)// &
           '" in file '//trim(fname)//' at line '//trim(num),'input', &
           'InputParErr')
      else
         call SysPrint('Error in input file '//trim(fname)//' for tag "' &
           //trim(label)//'" at line '//trim(num),'input')
         call SysIOErr(err,'input','InputParErr')
      end if
   end if

end subroutine InputParErr
!******************************************************************************
subroutine InputBlockErr(label,err,id)

   use sys,                  only : SysKill, SysWarning, SysPrint, SysIOErr

   implicit none

   character(*), intent(in) :: label
   integer, intent(in) :: err, id

   character(maxfnlen) :: fname
   integer :: ln
   character(6) :: num

   if (err/=0) then
      call InputGetLineNum(id, fname, ln)
      write(num,'(i0)') ln
      if (err==-2) then
         call SysWarning('Taking only integer part for block "'//trim(label)//&
           '" in file '//trim(fname)//' at line '//trim(num),'input', &
           'InputParErr')
      else
         call SysPrint('Error in input file '//trim(fname)//' for block "' &
           //trim(label)//'" at line '//trim(num),'input')
         call SysIOErr(err,'input','InputParErr')
      end if
   end if

end subroutine InputBlockErr
!****** End subroutine: InputParErr *******************************************
!******************************************************************************


!****** Subroutine: InputGetLineNum *******************************************
!******************************************************************************
!
!  Returns the line number for the corresponding line
!
! INPUT -----------------------------------------------------------------------
!
! integer id                     : ID for the processed line
!
! OUTPUT ----------------------------------------------------------------------
!
! character fname                : File name
! integer ln                     : Line number
!
!******************************************************************************
subroutine InputGetLineNum(id,fname,ln)

   use sys,                  only : SysIOErr

   implicit none

   integer, intent(in) :: id
   character(maxfnlen), intent(out) :: fname
   integer, intent(out) :: ln

   integer :: i

#ifdef MPI
   if (IsPNode) then
#endif /* MPI */
      rewind(uln)
#ifdef MPI
   end if
#endif /* MPI */
   do i=0,id
#ifdef MPI
      if (IsPNode) then
#endif /* MPI */
         read(uln,IOSTAT=iostat) fname, ln
#ifdef MPI
      end if
      call MPIBCast(iostat,1,MPI_INTEGER)
#endif /* MPI */
      call SysIOErr(iostat,'input','InputGetLine')
   end do
#ifdef MPI
   call MPIBCast(fname,maxfnlen,MPI_CHARACTER)
   call MPIBCast(ln,1,MPI_INTEGER)
#endif /* MPI */

end subroutine InputGetLineNum
!****** End subroutine: InputGetLineNum ***************************************
!******************************************************************************

!****** Subroutine: InputGetLine **********************************************
!******************************************************************************
!
!  Returns the contents of line with specified id
!
! INPUT -----------------------------------------------------------------------
!
! integer lineid                 : ID of line
!
! OUTPUT ----------------------------------------------------------------------
!
! character str                  : The text found in the line with the
!                                  specified id
!
!******************************************************************************
subroutine InputGetLine(lineid,str)

   use mem,                  only : Ssz
#ifdef MPI
   use sys,                  only : SysIOErrMPI
#else
   use sys,                  only : SysIOErr
#endif /* MPI */

   implicit none

   integer, intent(in) :: lineid
   character(*), intent(out) :: str

#ifdef MPI
   integer(kind=MPI_OFFSET_KIND) :: nbytes
#endif /* MPI */
   character(len=maxrecl) line
   integer :: i

#ifdef MPI
   nbytes = maxrecl*Ssz
   call MPI_File_Read_At_All(uin,lineid*nbytes,line,maxrecl,MPI_CHARACTER, &
     MPI_STATUS_IGNORE,iostat)
   call SysIOErrMPI(iostat,'input','InputGetLine')
   if (ichar(line(maxrecl:maxrecl))==10) then
      line(maxrecl:maxrecl) = ''
   end if
#else
   rewind(uin)
   do i=0,lineid
      read(uin,'(a)',IOSTAT=iostat) line
      call SysIOErr(iostat,'input','InputGetLine')
   end do
#endif /* MPI */
   str = line

end subroutine InputGetLine
!****** End subroutine: InputGetLine ******************************************
!******************************************************************************

end module input
