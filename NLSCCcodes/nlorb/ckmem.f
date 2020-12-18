      subroutine ckmem(pointer,max)

C-------------------------------------------------------------------
C  OPERATION: Check memory
C-------------------------------------------------------------------

      implicit none

      integer pointer, max
 
      if (pointer.ge.max) then
         write(*,10) pointer, max
         stop
      end if

 10   format('Asking for too much memory ',I25,' of ',I25)

      return
      end

