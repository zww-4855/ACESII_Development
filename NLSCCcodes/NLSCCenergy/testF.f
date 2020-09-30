










      SUBROUTINE testF(fia,fianew,T1Size,nocc,nvirt) 
        integer, intent(in)::nocc,nvirt
        double precision, intent(in)::fianew(nvirt,nocc)
        double precision, intent(in):: fia(nocc*nvirt)
        integer::iter

        iter=1
        do i=1,nocc
            do a=1,nvirt
                if (abs(fianew(a,i)- fia(iter)).gt.0.00004) then
                print*, fianew(a,i), fia(iter)
                endif
                iter=iter+1
        enddo   
        enddo
      END SUBROUTINE 
