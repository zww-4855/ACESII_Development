C*************************************************************
      subroutine Get1EInt(oneh,ovrlp,buf,ibuf,ldim)
C     
C     Get one electron integrals
C         
C*************************************************************
      implicit double precision (a-h,o-z)
c
      logical FileExist
c
      dimension oneh(ldim),ovrlp(ldim),
     &          buf(600),ibuf(600)
c
      ilnbuf=600
      FileExist=.false.
      !print*,"inside 1e- int"
      inquire(file='IIII',exist=FileExist)
      if (FileExist) then
c         write(*,*) 'IIII file exists'
         open(unit=10,file='IIII',form='UNFORMATTED',
     &        access='SEQUENTIAL')

         rewind 10

         call locate(10,'ONEHAMIL')
         call zero(oneh,ldim)
         nut = 1
         do while (nut.gt.0)
            read(10) buf, ibuf, nut
            do int = 1, nut
               oneh(ibuf(int)) = buf(int)
            end do
         end do
c
        !print*, "one hamil is: ", oneh
        ! print*, "outside first do"
         call locate(10,'OVERLAP ')
         call zero(ovrlp,ldim)
         nut = 1
         do while (nut.gt.0)
            read(10) buf, ibuf, nut
         !   print*, 'nut is: ', nut
            do int = 1, nut
          !     print*,"second loop::: int = ", int
               ovrlp(ibuf(int)) = buf(int)
            end do
         end do
         close(10)
         !print*, "done with 1e-"
      else
         write(*,*) 'IIII file does not exist'
         stop
      end if
      return
      end
