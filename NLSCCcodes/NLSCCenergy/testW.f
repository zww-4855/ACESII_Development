










      SUBROUTINE testW(Wab,Wabnew,T2abSize, nocc,nvirt) 
        integer, intent(in)::nocc,nvirt
        double precision, intent(in)::Wabnew(nvirt,nvirt,nocc,nocc)
        double precision, intent(in):: Wab(T2abSize)
        integer::iter

        iter=1
        do i=1,nocc
          do j=1,nocc
            do a=1,nvirt
              do b=1,nvirt
                if (abs(Wabnew(a,b,i,j)- Wab(iter)).gt.0.00004) then
                print*, Wabnew(a,b,i,j), Wab(iter)
                endif
                iter=iter+1
        enddo   
        enddo
        enddo   
        enddo
      END SUBROUTINE 
