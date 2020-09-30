










      SUBROUTINE test(t1,nocc,nvirt)
        integer, intent(in)::nocc,nvirt
        double precision, intent(in)::t1(nocc,nvirt)
        integer :: i
        do i=1,nocc
          do j=1,nvirt
                print*, t1(i,j)
          enddo
        enddo
      END SUBROUTINE test
