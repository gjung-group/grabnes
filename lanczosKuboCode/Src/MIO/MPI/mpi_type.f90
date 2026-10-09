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

module mpi[_TYPE]

   use mpitime,              only : MPITimer
   use mpi_inc
   use mpikinds

   implicit none

   PRIVATE

   integer :: ierr

   public :: MPIBCast
   interface MPIBCast
      module procedure MPIBCast[_TYPE], MPIBCast[_TYPE]_v
   end interface

   public :: MPIRedSum
   interface MPIRedSum
      module procedure MPIRedSum[_TYPE],  MPIRedSum[_TYPE]_v
   end interface

   public :: MPIAllRedSum
   interface MPIAllRedSum
      module procedure MPIAllRedSum[_TYPE],  MPIAllRedSum[_TYPE]_v
   end interface

   public :: MPIGather
   interface MPIGather
      module procedure MPIGather[_TYPE],  MPIGather[_TYPE]_v
   end interface

   public :: MPIAllGather
   interface MPIAllGather
      module procedure MPIAllGather[_TYPE],  MPIAllGather[_TYPE]_v
      module procedure MPIAllGather2[_TYPE], MPIAllGather2[_TYPE]_v
   end interface

   public :: MPIAllGatherV
   interface MPIAllGatherV
      module procedure MPIAllGatherV[_TYPE], MPIAllGatherV[_TYPE]_v
   end interface

   public :: MPISendRecv
   interface MPISendRecv
      module procedure MPISendRecv[_TYPE]
   end interface

contains

subroutine MPIBCast[_TYPE](buffer,sz,type)
   [TYPE], intent(inout) :: buffer
   integer :: sz
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIBCast[_TYPE]
subroutine MPIBCast[_TYPE]_v(buffer,sz,type)
   [TYPE], intent(inout) :: buffer(*)
   integer :: sz
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_Bcast(buffer,sz,type,rankPNode,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIBCast[_TYPE]_v

subroutine MPIRedSum[_TYPE](sendbuf,recvbuf,sz,type)
   [TYPE], intent(in) :: sendbuf
   [TYPE], intent(out) :: recvbuf
   integer, intent(in) :: sz
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIRedSum[_TYPE]
subroutine MPIRedSum[_TYPE]_v(sendbuf,recvbuf,sz,type)
   [TYPE], intent(in) :: sendbuf(*)
   [TYPE], intent(out) :: recvbuf(*)
   integer, intent(in) :: sz
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_Reduce(sendbuf,recvbuf,sz,type,MPI_SUM,rankPNode,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIRedSum[_TYPE]_v

subroutine MPIAllRedSum[_TYPE](sendbuf,recvbuf,sz,type)
   [TYPE], intent(in) :: sendbuf
   [TYPE], intent(out) :: recvbuf
   integer, intent(in) :: sz
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIAllRedSum[_TYPE]
subroutine MPIAllRedSum[_TYPE]_v(sendbuf,recvbuf,sz,type)
   [TYPE], intent(in) :: sendbuf(*)
   [TYPE], intent(out) :: recvbuf(*)
   integer, intent(in) :: sz
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_AllReduce(sendbuf,recvbuf,sz,type,MPI_SUM,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIAllRedSum[_TYPE]_v

subroutine MPIGather[_TYPE](sendbuf,sendsz,recvbuf,recvsz,type)
   [TYPE], intent(in) :: sendbuf
   [TYPE], intent(out) :: recvbuf(*)
   integer, intent(in) :: sendsz, recvsz
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIGather[_TYPE]
subroutine MPIGather[_TYPE]_v(sendbuf,sendsz,recvbuf,recvsz,type)
   [TYPE], intent(in) :: sendbuf(*)
   [TYPE], intent(out) :: recvbuf(*)
   integer, intent(in) :: sendsz, recvsz
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_Gather(sendbuf,sendsz,type,recvbuf,recvsz,type,rankPNode,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIGather[_TYPE]_v

subroutine MPIAllGather[_TYPE](sendbuf,sendsz,recvbuf,recvsz,type)
   [TYPE], intent(in) :: sendbuf
   [TYPE], intent(out) :: recvbuf(*)
   integer, intent(in) :: sendsz, recvsz
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIAllGather[_TYPE]
subroutine MPIAllGather[_TYPE]_v(sendbuf,sendsz,recvbuf,recvsz,type)
   [TYPE], intent(in) :: sendbuf(*)
   [TYPE], intent(out) :: recvbuf(*)
   integer, intent(in) :: sendsz, recvsz
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIAllGather[_TYPE]_v
subroutine MPIAllGather2[_TYPE](sendbuf,sendsz,recvbuf,recvsz,type)
   [TYPE], intent(in) :: sendbuf
   [TYPE], intent(out) :: recvbuf(:,:)
   integer, intent(in) :: sendsz, recvsz
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIAllGather2[_TYPE]
subroutine MPIAllGather2[_TYPE]_v(sendbuf,sendsz,recvbuf,recvsz,type)
   [TYPE], intent(in) :: sendbuf(*)
   [TYPE], intent(out) :: recvbuf(:,:)
   integer, intent(in) :: sendsz, recvsz
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_AllGather(sendbuf,sendsz,type,recvbuf,recvsz,type,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIAllGather2[_TYPE]_v


subroutine MPIAllGatherV[_TYPE](sendbuf,sendsz,recvbuf,recvsz,displs,type)
   integer, intent(in) :: sendbuf
   [TYPE], intent(out) :: recvbuf(*)
   integer, intent(in) :: sendsz, recvsz(*), displs(*)
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIAllGatherV[_TYPE]
subroutine MPIAllGatherV[_TYPE]_v(sendbuf,sendsz,recvbuf,recvsz,displs,type)
   [TYPE], intent(in) :: sendbuf(*)
   [TYPE], intent(out) :: recvbuf(*)
   integer, intent(in) :: sendsz, recvsz(*), displs(*)
   integer, intent(in) :: type

   call MPITimer(1)
   call MPI_AllGatherV(sendbuf,sendsz,type,recvbuf,recvsz,displs,type,MPI_COMM_WORLD,ierr)
   call MPITimer(0)
end subroutine MPIAllGatherV[_TYPE]_v

subroutine MPISendRecv[_TYPE](sendbuf,sendsz,dest,sendtag,recvbuf,recvsz,source,recvtag,type)
   [TYPE], intent(in) :: sendbuf(*)
   [TYPE], intent(out) :: recvbuf(*)
   integer, intent(in) :: sendsz, recvsz, dest, source
   integer, intent(in) :: sendtag, recvtag, type

   call MPITimer(1)
   call MPI_SendRecv(sendbuf,sendsz,type,dest,sendtag,recvbuf,recvsz,type,source,recvtag,&
     MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
   call MPITimer(0)
end subroutine MPISendRecv[_TYPE]

end module mpi[_TYPE]
