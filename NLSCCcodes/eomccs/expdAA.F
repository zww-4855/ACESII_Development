      SUBROUTINE expdAA(Waa,Waanew,T2aaSize, nocc,nvirt) 
        integer, intent(in)::nocc,nvirt
        double precision, intent(inout)::Waanew(nvirt,nvirt,nocc,nocc)
        double precision, intent(in):: Waa(T2aaSize)
        integer::iter

        iter=1
        Waanew=0.0d0
        do i=1,nocc
          do j=1,nocc
            do a=1,nvirt
              do b=1,nvirt
                if (j<i.and.b<a) then
                        Waanew(a,b,i,j)=Waa(iter)
                        iter=iter+1
                endif
        enddo   
        enddo
        enddo   
        enddo
      END SUBROUTINE 
