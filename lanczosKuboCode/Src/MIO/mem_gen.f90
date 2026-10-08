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

program mem_gen

   use precision

   implicit none

   integer, parameter :: Isz = kind(1)
   integer, parameter :: Rsz = kind(1.0_sp)
   integer, parameter :: Dsz = kind(1.0_dp)
   integer, parameter :: Csz = Rsz*2
   integer, parameter :: Zsz = Dsz*2
   integer, parameter :: Ssz = kind('a')
   integer, parameter :: Lsz = kind(.true.)

   integer, parameter :: sz(7) = (/Isz,Rsz,Dsz,Lsz,Csz,Zsz,Ssz/)

   character(*), parameter :: types(7) = ['integer    ','real(sp)   ', &
                              'real(dp)   ','logical    ', &
                              'complex(sp)','complex(dp)','character  ']
   character(*), parameter :: label(7) = ['_i','_r','_d','_l','_c','_z','_s']
   character(*), parameter :: init(7) = ['0              ','0.0_sp         ', &
                              '0.0_dp         ', '.false.        ', &
                              '(0.0_sp,0.0_sp)', '(0.0_dp,0.0_dp)',"''             "]
   integer, parameter :: maxrank=4

   character(len=1024) :: flin1,flin2,flout
   character(len=100) :: line, line2
   character(len=60) :: str
   character :: rnk
   logical :: exist
   integer :: lh, il, it, indx, rank, ir

   call getarg(1,flin1)
   call getarg(2,flin2)
   call getarg(3,flout)
   inquire(FILE=flin1,EXIST=exist)
   if (exist) then
      inquire(FILE=flin2,EXIST=exist)
   end if
   if (exist) then
      write(*,'(a)') "Generating file '"//trim(flout)//"'"
   else
      write(*,'(a)') "Cannot generate file '"//trim(flout)//"'"
      stop
   end if
   open(10,FILE=flin1)
   open(11,FILE=flin2)
   open(12,FILE=flout,STATUS='REPLACE')

   do
      read(10,'(a)',END=100) line
      indx = index(line,'[ISZ]')
      if  (indx/=0) then
         write(line2,'(a,i0,a)') line(:indx-1),Isz,trim(line(indx+5:))
         write(12,'(a)') line2
         cycle
      end if
      indx = index(line,'[RSZ]')
      if  (indx/=0) then
         write(line2,'(a,i0,a)') line(:indx-1),Rsz,trim(line(indx+5:))
         write(12,'(a)') line2
         cycle
      end if
      indx = index(line,'[DSZ]')
      if  (indx/=0) then
         write(line2,'(a,i0,a)') line(:indx-1),Dsz,trim(line(indx+5:))
         write(12,'(a)') line2
         cycle
      end if
      indx = index(line,'[CSZ]')
      if  (indx/=0) then
         write(line2,'(a,i0,a)') line(:indx-1),Csz,trim(line(indx+5:))
         write(12,'(a)') line2
         cycle
      end if
      indx = index(line,'[ZSZ]')
      if  (indx/=0) then
         write(line2,'(a,i0,a)') line(:indx-1),Zsz,trim(line(indx+5:))
         write(12,'(a)') line2
         cycle
      end if
      indx = index(line,'[SSZ]')
      if  (indx/=0) then
         write(line2,'(a,i0,a)') line(:indx-1),Ssz,trim(line(indx+5:))
         write(12,'(a)') line2
         cycle
      end if
      indx = index(line,'[LSZ]')
      if  (indx/=0) then
         write(line2,'(a,i0,a)') line(:indx-1),Lsz,trim(line(indx+5:))
         write(12,'(a)') line2
         cycle
      end if
      indx = index(line,'[WRITE_INTERFACE]')
      if  (indx/=0) then
         write(12,'(a)') '   public :: MemAlloc'
         write(12,'(a)') '   public :: MemDealloc'
         write(12,*)
         write(12,'(a)') '   interface MemAlloc'
         do it=1,size(types)
            do rank=1,maxrank
               write(rnk,'(i0)') rank
               str = 'Alloc'//label(it)//rnk
               line = '      module procedure '//trim(str)//'_s, '//trim(str)//'_c'
               write(12,'(a)') trim(line)
            end do
         end do
         write(12,'(a)') '      module procedure Alloc_is'
         write(12,'(a)') '      module procedure Alloc_rs, Alloc_ds'
         write(12,'(a)') '      module procedure Alloc_cs, Alloc_zs'
         write(12,'(a)') '      module procedure Alloc_ls, Alloc_ss'
         write(12,'(a)') '   end interface'
         write(12,'(a)') '   interface MemDealloc'
         do it=1,size(types)
            do rank=1,maxrank
               write(rnk,'(i0)') rank
               str = 'Dealloc'//label(it)//rnk
               line = '      module procedure '//trim(str)
               write(12,'(a)') trim(line)
            end do
         end do
         write(12,'(a)') '   end interface'
         cycle
      end if
      indx = index(line,'[WRITE_SUBS]')
      if  (indx/=0) then
         lh = 0
         read(11,*) line
         do while(line(1:1)=='!')
            lh = lh+1
            read(11,*) line
         end do
         rewind(11)

         do it=1,size(types)
            do rank=1,maxrank
               write(rnk,'(i0)') rank
               do il=1,lh
                  read(11,*)
               end do
               do
                  read(11,'(a)',END=101) line
                  indx = index(line,'[TYPE]')
                  do while(indx/=0)
                     if (types(it)=='character') then
                        line2 = line(:indx-1)//trim(types(it))//'(*)'//line(indx+6:)
                     else
                        line2 = line(:indx-1)//trim(types(it))//line(indx+6:)
                     endif
                     line = line2
                     indx = index(line,'[TYPE]')
                  end do
                  indx = index(line,'[TYPETMP]')
                  do while(indx/=0)
                     if (types(it)=='character') then
                        line2 = line(:indx-1)//trim(types(it))// &
                                '(len=len(array))'//line(indx+9:)
                     else
                        line2 = line(:indx-1)//trim(types(it))//line(indx+9:)
                     endif
                     line = line2
                     indx = index(line,'[TYPETMP]')
                  end do
                  indx = index(line,'[_TYPE]')
                  do while(indx/=0)
                     line2 = line(:indx-1)//trim(label(it))//rnk//line(indx+7:)
                     line = line2
                     indx = index(line,'[_TYPE]')
                  end do
                  indx = index(line,'[RANK_DEC]')
                  do while(indx/=0)
                     str = ':'
                     do ir=2,rank
                        str = trim(str)//',:'
                     end do
                     line2 = line(:indx-1)//trim(str)//line(indx+10:)
                     line = line2
                     indx = index(line,'[RANK_DEC]')
                  end do
                  indx = index(line,'[RANK]')
                  do while(indx/=0)
                     line2 = line(:indx-1)//rnk//line(indx+6:)
                     line = line2
                     indx = index(line,'[RANK]')
                  end do
                  indx = index(line,'[TYPESZ]')
                  do while(indx/=0)
                     if (types(it)=='character') then
                        write(line2,'(a,i0,a)') line(:indx-1)//'len(array)*',sz(it),trim(line(indx+8:))
                     else
                        write(line2,'(a,i0,a)') line(:indx-1),sz(it),line(indx+8:)
                     endif
                     line = line2
                     indx = index(line,'[TYPESZ]')
                  end do
                  indx = index(line,'[BOUNDS]')
                  do while(indx/=0)
                     line2 = 'szl(1):szu(1)'
                     do ir=2,rank
                        write(str,'(a,i0,a,i0,a)') ',szl(',ir,'):szu(',ir,')'
                        line2 = trim(line2)//str
                     end do
                     line2 = line(:indx-1)//trim(line2)//line(indx+8:)
                     line = line2
                     indx = index(line,'[BOUNDS]')
                  end do
                  indx = index(line,'[COPY]')
                  do while(indx/=0)
                     write(line2,'(a,i0,a,i0,a)') 'cpBnd(',1,'):cpBnd(',rank+1,')'
                     do ir=2,rank
                        write(str,'(a,i0,a,i0,a)') ',cpBnd(',ir,'):cpBnd(',rank+ir,')'
                        line2 = trim(line2)//str
                     end do
                     line2 = line(:indx-1)//trim(line2)//line(indx+6:)
                     line = line2
                     indx = index(line,'[COPY]')
                  end do
                  indx = index(line,'[INIT]')
                  do while(indx/=0)
                     line2 = line(:indx-1)//init(it)//line(indx+6:)
                     line = line2
                     indx = index(line,'[INIT]')
                  end do

                  write(12,'(a)') trim(line)
               end do
101            continue
               rewind(11)
            end do
         end do
         cycle
      end if
      write(12,'(a)') trim(line)
   end do
100 continue
   close(10)
   close(11)
   close(12)

end program mem_gen
