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

module file

   use io,                 only : maxfnlen, maxlinel
#ifdef MPI
   use mpi
#endif
   implicit none

   PRIVATE

   type, public :: cl_file
      PRIVATE
      character(len=maxfnlen) :: name = 'unknown'
      integer :: unit
      character(len=11) :: form='FORMATTED'
      character(len=7) :: status='UNKNOWN'
      character(len=11) :: access='SEQUENTIAL'
      character(len=6) :: position='ASIS'
      character(len=9) :: action='READWRITE'
      character(len=6) :: delete='KEEP'
      integer :: recl=maxlinel
      integer :: iostat=0
      logical :: opened=.false.
      logical :: init=.false.
      class(cl_file), pointer :: prev => NULL()
      class(cl_file), pointer :: next => NULL()
#ifdef MPI
      logical :: serial=.false.
      integer :: comm=MPI_COMM_WORLD
      integer :: amode
      integer :: disp
      integer :: etype
      integer :: filetype
      character(len=MPI_MAX_DATAREP_STRING) :: datarep
#endif /* MPI */
   contains
      PRIVATE
      procedure, public :: Initialize
      procedure, public :: Destructor
      procedure, public :: Open
      procedure, public :: SetForm
      procedure, public :: SetStatus
      procedure, public :: SetAccess
      procedure, public :: SetPosition
      procedure, public :: SetAction
      procedure, public :: SetRecl
      procedure, public :: SetDelete
      procedure, public :: Close
      procedure, public :: GetUnit
      procedure, public :: GetName
      procedure :: FindUnit
      procedure :: Error
      procedure, public :: IsOpen
      procedure, public :: IsInit
      procedure, public :: Rewind
#ifdef MPI
      procedure :: SetAmode
      procedure, public :: SetSerial
      procedure, public :: GetFH
      procedure, public :: SetView
#endif /* MPI */
   end type cl_file

   public :: CloseAll

   class(cl_file), pointer, save :: headFile => NULL()

contains


!****** Subroutine: Initialize ************************************************
!******************************************************************************
!
!  Initialization routine for the file
!
!******************************************************************************
subroutine Initialize(this)

   use sys,                  only : SysWarning

   class(cl_file), intent(inout), target :: this

   if (this%init) then
      call SysWarning( &
        'Attempt to initialize an already initialized file avoided','file')
      return
   end if
   if (associated(headFile)) then
      headFile%next => this
   end if
   this%prev => headFile
   this%next => NULL()
   headFile => this
   this%init = .true.

end subroutine Initialize
!****** End subroutine: Initialize ********************************************
!******************************************************************************


!****** Subroutine: Destructor ************************************************
!******************************************************************************
!
!  Destructor routine for the file
!
!******************************************************************************
subroutine Destructor(this)

   use sys,                  only : SysWarning

   class(cl_file), intent(inout) :: this

   if (.not. this%init) then
      call SysWarning('Attempt to destroy an uninitialized file avoided',&
        'file')
      return
   end if
   if (this%opened)  then
      call SysWarning('Attempt to destroy an opened file avoided',&
        'file')
      return
   end if
   this%init = .false.
   if (associated(this%prev)) then
      this%prev%next => this%next
   end if
   if (associated(this%next)) then
      this%next%prev => this%prev
   else
      headFile => this%prev
   end if

end subroutine Destructor
!****** End subroutine: Destructor ********************************************
!******************************************************************************


!****** Subroutine: Open ******************************************************
!******************************************************************************
!
!  Open file
!
!******************************************************************************
subroutine Open(this,name,form,status,access,position,action,recl,unit,serial)

#ifdef MPI
   use mpi,                  only : IsPNode
#endif /* MPI */
   use sys,                  only : SysKill
   use format,               only : num2string

   class(cl_file), intent(inout) :: this
   character(len=*), intent(in), optional :: name
   character(len=*), intent(in), optional :: form
   character(len=*), intent(in), optional :: status
   character(len=*), intent(in), optional :: access
   character(len=*), intent(in), optional :: position
   character(len=*), intent(in), optional :: action
   integer, intent(in), optional :: recl, unit
   logical, intent(in), optional :: serial

#ifdef MPI
   integer :: amode
#endif /* MPI */
   logical :: opened

   if (.not. this%init) then
      call this%Initialize()
   end if
   if (present(name)) then
      this%name = name
   end if
#ifdef MPI
   if (present(serial)) then
      this%serial=serial
   else
      this%serial=.false.
   end if
#endif /* MPI */
   if (present(unit)) then
      inquire(unit, OPENED=opened)
      if (opened) then
         call SysKill("File unit not availabe for file '"//trim(this%name)// &
           "'", 'file')
      endif
      this%unit = unit
#ifdef MPI
      this%serial = .true.
#endif /* MPI */
   else
#ifdef MPI
      if (this%serial) then
#endif /* MPI */
         call this%FindUnit()
#ifdef MPI
      end if
#endif /* MPI */
   end if
   if (.not. present(name)) then
      if(this%name=='unknown') then
         this%name = 'file.'//trim(num2string(this%unit))
      end if
#ifdef MPI
      call MPIBCast(this%name,len(this%name),MPI_CHARACTER)
#endif /* MPI */
   end if
   if (present(form)) then
      call this%SetForm(form)
   end if
   if (present(status)) then
      call this%SetStatus(status)
#ifdef MPI
   else
      if (this%serial) then
         call this%SetStatus('UNKNOWN')
      else
         call this%SetStatus('REPLACE')
      end if
#endif /* MPI */
   end if
   if (present(access)) then
      call this%SetAccess(access)
#ifdef MPI
   else
      if (this%serial) then
         call this%SetAccess('SEQUENTIAL')
      else
         call this%SetAccess('DIRECT')
      end if
#endif /* MPI */
   end if
   if (present(position)) then
      call this%SetPosition(position)
   end if
   if (present(action)) then
      call this%SetAction(action)
   end if
   if (present(recl)) then
      call this%SetRecl(recl)
   end if
#ifdef MPI
   if (this%serial) then
      if (IsPNode) then
#endif /* MPI */
         open(this%unit,IOSTAT=this%iostat, FILE=this%name, STATUS=this%status,&
           ACCESS=this%access,FORM=this%form, RECL=this%recl, &
           POSITION=this%position, ACTION=this%action)
#ifdef MPI
      endif
   else
      call this%SetAmode()
      call MPI_File_Open(MPI_COMM_WORLD,trim(this%name),this%amode,&
        MPI_INFO_NULL,this%unit,this%iostat)
      this%opened = .true.
   endif
#endif /* MPI */
   call this%Error('OPEN')
   this%opened = .true.

end subroutine Open
!****** End subroutine: Open **************************************************
!******************************************************************************


!****** Subroutine: Close *****************************************************
!******************************************************************************
!
!  Close file
!
!******************************************************************************
subroutine Close(this,destroy)

#ifdef MPI
   use mpi
#endif /* MPI */

   class(cl_file), intent(inout) :: this
   logical, optional :: destroy

   character(len=6) :: status

   integer :: io

#ifdef MPI
   if (this%serial) then
      if (IsPNode) then
#endif /* MPI */
         close(this%unit,IOSTAT=this%iostat, &
           STATUS=this%delete)
#ifdef MPI
      endif
   else
      call MPI_File_Close(this%unit,this%iostat)
   endif
#endif /* MPI */
   call this%Error('CLOSE')
   this%opened = .false.
   if (present(destroy)) then
      if (destroy) call this%Destructor()
   else
      call this%Destructor()
   end if

end subroutine Close
!****** End subroutine: Close *************************************************
!******************************************************************************

!****** Subroutine: CloseAll **************************************************
!******************************************************************************
!
!  Close all files
!
!******************************************************************************
subroutine CloseAll()

   do while(associated(headFile))
      call headFile%Close()
   end do

end subroutine CloseAll
!****** End subroutine: CloseAll **********************************************
!******************************************************************************


!****** Subroutine: FindUnit **************************************************
!******************************************************************************
!
!  Find an available unit number
!
!******************************************************************************
subroutine FindUnit(this)

   use io,                   only : unit0

   class(cl_file), intent(inout) :: this

   integer :: u
   logical :: file_opened

   u = unit0
   file_opened = .true.
   do while (file_opened)
      u = u+1
      inquire(u, OPENED=file_opened )
   end do
   this%unit = u

end subroutine FindUnit
!****** End subroutine: FindUnit **********************************************
!******************************************************************************


!****** Subroutine: Error *****************************************************
!******************************************************************************
!
!  Error handler
!
!******************************************************************************
subroutine Error(this,action)

   use sys,                  only : SysPrint, SysKill, SysIOErr
#ifdef MPI
   use sys,                  only : SysIOErrMPI
#endif /* MPI */
   use format,               only : num2string

   class(cl_file), intent(inout) :: this
   character(*), intent(in) :: action

   if (this%iostat==0) return
   if (action=='OPEN') then
      if (this%iostat==-2) then
         call SysKill("OPEN: Found EOR while opening file '"&
           //trim(this%name)//"'",'file')
      elseif (this%iostat==-1) then
         call SysKill("OPEN: Found EOF while opening file '"&
           //trim(this%name)//"'",'file')
      else
         call SysPrint('')
         call SysPrint(trim(action)//": Error in file '"// &
           trim(this%name)//"'",'file')
      end if
   else
      call SysPrint('')
      call SysPrint(trim(action)//": Error in file '"// &
        trim(this%name)//"'",'file')
   end if
#ifdef MPI
   if (this%serial) then
#endif /* MPI*/
      call SysIOErr(this%iostat)
#ifdef MPI
   else
      call SysIOErrMPI(this%iostat)
   end if
#endif /* MPI */

end subroutine Error
!****** End subroutine: Error *************************************************
!******************************************************************************


!****** Subroutines: Set* *****************************************************
!******************************************************************************
!
!  Routines for setting values of 'this'
!
!******************************************************************************
subroutine SetForm(this,form)

   use sys,                  only : SysKill, SysPrint

   class(cl_file), intent(inout) :: this
   character(*), intent(in) :: form

   if (form=='FORMATTED' .or. form=='formatted' .or. form=='f' .or. &
     form=='F') then
      this%form = 'FORMATTED'
   elseif (form=='UNFORMATTED' .or. form=='unformatted' .or. form=='u' .or. &
     form=='U') then
      this%form = 'UNFORMATTED'
   elseif (form=='BINARY' .or. form=='binary' .or. form=='b' .or. &
     form=='B') then
      this%form = 'BINARY'
   else
      call SysPrint("ERROR: Wrong argument for 'form' file specification", &
        'file')
      call SysPrint("       Possible values are: 'FORMATTED', "// &
        "'UNFORMATTED', 'BINARY'",'file')
      call SysKill('Wrong specification for file','file')
   end if

end subroutine SetForm


subroutine SetStatus(this,status)

   use sys,                  only : SysKill, SysPrint

   class(cl_file), intent(inout) :: this
   character(*), intent(in) :: status

   if (status=='NEW' .or. status=='new' .or. status=='n' .or. status=='N') then
      this%status = 'NEW'
   elseif (status=='REPLACE' .or. status=='replace' .or. status=='r' .or. &
     status=='R') then
      this%status = 'REPLACE'
   elseif (status=='SCRATCH' .or. status=='scratch' .or. status=='s' .or. &
     status=='S') then
      this%status = 'REPLACE'
      this%delete = 'DELETE'
   elseif (status=='OLD' .or. status=='old' .or. status=='o' .or. &
     status=='O') then
      this%status = 'OLD'
   elseif (status=='UNKNOWN' .or. status=='unknown' .or. status=='u' .or. &
     status=='U') then
      this%status = 'UNKNOWN'
   else
      call SysPrint("ERROR: Wrong argument for 'status' file specification",&
        'mfile')
      call SysPrint("       Possible values are: 'NEW', 'REPLACE'"// &
        ", 'SCRATCH', 'OLD', 'UNKNOWN'",'mfile')
      call SysKill('Wrong specification for file','file')
   end if

end subroutine SetStatus


subroutine SetAccess(this,access)

   use sys,                  only : SysKill, SysPrint

   class(cl_file), intent(inout) :: this
   character(*), intent(in) :: access

   if (access=='SEQUENTIAL' .or. access=='sequential' .or. access=='s' .or. &
     access=='S') then
      this%access = 'SEQUENTIAL'
   elseif (access=='DIRECT' .or. access=='direct' .or. access=='d' .or. &
     access=='D') then
      this%access = 'DIRECT'
   elseif (access=='TRANSPARENT' .or. access=='transparent' .or. &
     access=='t' .or. access=='T') then
      this%access = 'TRANSPARENT'
   else
      call SysPrint("ERROR: Wrong argument for 'access' file specification",&
        'mfile')
      call SysPrint("       Possible values are: 'SEQUENTIAL', 'DIRECT'"// &
        ", 'TRANSPARENT'",'mfile')
      call SysKill('Wrong specification for file','file')
   end if

end subroutine SetAccess


subroutine SetPosition(this,position)

   use sys,                  only : SysKill, SysPrint

   class(cl_file), intent(inout) :: this
   character(*), intent(in) :: position

   if (position=='REWIND' .or. position=='rewind' .or. position=='r' .or. &
     position=='R') then
      this%position = 'REWIND'
   elseif (position=='APPEND' .or. position=='append' .or. position=='ap' .or. &
     position=='AP') then
      this%position = 'APPEND'
   elseif (position=='ASIS' .or. position=='asis' .or. position=='a' .or. &
     position=='A') then
      this%position = 'ASIS'
   else
      call SysPrint("ERROR: Wrong argument for 'position' file "// &
        "specification",'mfile')
      call SysPrint("       Possible values are: 'REWIND', 'APPEND'"// &
        ", 'AIS'",'mfile')
      call SysKill('Wrong specification for file','file')
   end if

end subroutine SetPosition


subroutine SetAction(this,action)

   use sys,                  only : SysKill, SysPrint

   class(cl_file), intent(inout) :: this
   character(*), intent(in) :: action

   if (action=='READ' .or. action=='read' .or. action=='r' .or. &
     action=='R') then
      this%action = 'READ'
   elseif (action=='WRITE' .or. action=='write' .or. action=='w' .or. &
     action=='W') then
      this%action = 'WRITE'
   elseif (action=='READWRITE' .or. action=='readwrite' .or. action=='rw' &
     .or. action=='RW') then
      this%action = 'READWRITE'
   else
      call SysPrint("ERROR: Wrong argument for 'action' file specification"&
        ,'mfile')
      call SysPrint("       Possible values are: 'READ', 'WRITE'"// &
        ", 'READWRITE'",'mfile')
      call SysKill('Wrong specification for file','file')
   end if

end subroutine SetAction


subroutine SetRecl(this,recl)

   class(cl_file), intent(inout) :: this
   integer, intent(in) :: recl

   this%recl = recl

end subroutine SetRecl


#ifdef MPI
subroutine SetSerial(this,serial)

   class(cl_file), intent(inout) :: this
   logical, intent(in) :: serial

   this%serial = serial

end subroutine SetSerial
#endif /* MPI */


subroutine SetDelete(this,delete)

   class(cl_file), intent(inout) :: this
   logical, intent(in) :: delete

   if (delete) then
      this%delete = 'DELETE'
   else
      this%delete =  'KEEP'
   end if

end subroutine SetDelete
!****** End subroutines: Set* *************************************************
!******************************************************************************


!****** Function: GetUnit *****************************************************
!******************************************************************************
!
!  Get the unit of file
!
!******************************************************************************
function GetUnit(this) result(unit)

   use io,                 only : stderr

   integer :: unit
   class(cl_file), intent(inout) :: this

#ifdef MPI
   if (this%serial) then
#endif /* MPI */
      unit = this%unit
#ifdef MPI
   else
      unit = stderr
   end if
#endif /* MPI */

end function GetUnit
!****** End function: GetUnit *************************************************
!******************************************************************************

#ifdef MPI
!****** Function: GetFH *******************************************************
!******************************************************************************
!
!  Get file handler of file
!
!******************************************************************************
function GetFH(this) result(fh)

   integer :: fh
   class(cl_file), intent(inout) :: this

   fh = this%unit

end function GetFH
!****** End function: GetFH ***************************************************
!******************************************************************************
#endif /* MPI */


!****** Function: GetName *****************************************************
!******************************************************************************
!
!  Get the the name of the file
!
!******************************************************************************
function GetName(this) result(name)

   character(len=maxfnlen) :: name
   class(cl_file), intent(inout) :: this

   name = this%name

end function GetName
!****** End function: GetName *************************************************
!******************************************************************************


!****** function: IsOpen ******************************************************
!******************************************************************************
!
!  True if file is open
!
!******************************************************************************
function IsOpen(this)

   class(cl_file), intent(inout) :: this

   logical :: IsOpen

   IsOpen = this%opened

end function IsOpen
!****** End function: IsOpen **************************************************
!******************************************************************************


!****** function: IsInit ******************************************************
!******************************************************************************
!
!  True if file is initialized
!
!******************************************************************************
function IsInit(this)

   class(cl_file), intent(inout) :: this

   logical :: IsInit

   IsInit = this%init

end function IsInit
!****** End function: IsInit **************************************************
!******************************************************************************


#ifdef MPI
!****** Subroutine: SetAMode **************************************************
!******************************************************************************
!
!  Set the file access mode.
!
!******************************************************************************
subroutine SetAmode(this)

   class(cl_file), intent(inout) :: this

   this%amode = 0

   if (this%action=='READ') then
      this%amode = this%amode + MPI_MODE_RDONLY
   else if (this%action=='WRITE') then
      this%amode = this%amode + MPI_MODE_WRONLY
   else if (this%action=='READWRITE') then
      this%amode = this%amode + MPI_MODE_RDWR
   end if
   if (this%status=='REPLACE') then
      this%amode = this%amode + MPI_MODE_CREATE
   else if (this%status=='NEW') then
      this%amode = this%amode + MPI_MODE_EXCL
   end if
   if (this%delete=='DELETE') then
      this%amode = this%amode + MPI_MODE_DELETE_ON_CLOSE
   end if
   if (this%position=='APPEND') then
      this%amode = this%amode + MPI_MODE_APPEND
   end if
   if (this%access=='SEQUENTIAL') then
      this%amode = this%amode + MPI_MODE_SEQUENTIAL
   end if

end subroutine SetAmode
!****** End subroutine: SetAMode **********************************************
!******************************************************************************
#endif /* MPI */

#ifdef MPI
!****** Subroutine: SetView ***************************************************
!******************************************************************************
!
!  Read input file(s) and prepare a single file
!
! INPUT -----------------------------------------------------------------------
!
! integer disp                   : Position where the view begins
! integer etype                  : The data unit of the file
! integer filetype               : The distribution of the data
!
!******************************************************************************
subroutine SetView(this,disp,etype,filetype,datarep)

   use sys,                  only : SysIOErrMPI

   class(cl_file), intent(inout) :: this
   integer, intent(in) :: disp, etype, filetype
   character(*), intent(in) :: datarep

   this%disp = disp
   this%etype = etype
   this%filetype = filetype
   this%datarep = datarep
   call MPI_File_Set_View(this%unit,disp,etype,filetype,datarep, &
     MPI_INFO_NULL,this%iostat)
   call SysIOErrMPI(this%iostat)

end subroutine SetView
!****** End subroutine: SetView ***********************************************
!******************************************************************************
#endif /* MPI */


!****** Subroutine: Rewind ****************************************************
!******************************************************************************
!
!  Rewinds file to its initial position
!
!******************************************************************************
subroutine Rewind(this)

   use sys,                  only : SysIOErr
#ifdef MPI
   use sys,                  only : SysIOErrMPI
#endif /* MPI */

#ifdef MPI
   integer(kind=MPI_OFFSET_KIND),parameter :: zero=0
#endif /* MPI */

   class(cl_file), intent(inout), target :: this

#ifdef MPI
   if (this%serial) then
#endif /* MPI */
      rewind(this%unit,IOSTAT=this%iostat)
      call SysIOErr(this%iostat,'file','Rewind')
#ifdef MPI
   else
      call MPI_File_Set_View(this%unit,zero,this%etype,this%filetype, &
        this%datarep,MPI_INFO_NULL,this%iostat)
      call SysIOErrMPI(this%iostat,'file','Rewind')
   end if
#endif /* MPI */

end subroutine Rewind
!****** End subroutine: Rewind ************************************************
!******************************************************************************

end module file
