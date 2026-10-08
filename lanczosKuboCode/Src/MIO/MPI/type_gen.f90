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

program type_gen

   implicit none

   character(*), parameter :: types(7) = ['integer     ','real(sp)    ', &
                              'real(dp)    ','character(*)','logical     ', &
                              'complex(sp) ','complex(dp) ']
   character(*), parameter :: label(7) = ['_i','_r','_d','_s','_l','_c','_z']

   character(len=1024) :: fl
   character(len=100) :: line, line2
   logical :: exist
   integer :: lh, il, it, indx

   call getarg(1,fl)
   inquire(FILE=fl,EXIST=exist)
   if (exist) then
      write(*,'(a)') "Generating file 'mpi_types.f90'"
   else
      write(*,'(a)') "Cannot generate file 'mpi_types.f90'"
      stop
   end if
   open(10,FILE=fl)
   open(11,FILE='mpi_types.f90',STATUS='REPLACE')
   read(10,*) line
   lh = 0
   do while(line(1:1)=='!')
      write(11,'(a)') line
      lh = lh+1
      read(10,'(a)') line
   end do
   rewind(10)
   do it=1,size(types)
      do il=1,lh
         read(10,*)
      end do
      do
         read(10,'(a)',END=100) line
         indx = index(line,'[TYPE]')
         do while(indx/=0)
            line2 = line(:indx-1)//trim(types(it))//line(indx+6:)
            line = line2
            indx = index(line,'[TYPE]')
         end do
         indx = index(line,'[_TYPE]')
         do while(indx/=0)
            line2 = line(:indx-1)//trim(label(it))//line(indx+7:)
            line = line2
            indx = index(line,'[_TYPE]')
         end do
         write(11,'(a)') line
      end do
100   continue
      rewind(10)
   end do

   close(10)
   close(11)

end program type_gen
