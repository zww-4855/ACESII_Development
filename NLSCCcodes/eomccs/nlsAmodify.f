










        subroutine nlsAmodify(CISmat0,nocc,nvirt,acesCISevecs,evecSize,
     &                          nroot)
        integer, intent(in)::nroot,nocc,nvirt,evecSize
        double precision,intent(in)::acesCISevecs(evecSize),
     &                  CISmat0(2*nocc*nvirt,2*nocc*nvirt)

        integer::MATDIM,iter
        double precision::temp,matmulcopy(2*nocc*nvirt)
        double precision::intermed(2*nocc*nvirt)

        MATDIM=2*nocc*nvirt
        iter=1
        do i=1,nroot
           call xgemm('N','N',MATDIM,1,MATDIM,1.0D0,CISmat0,
     &               MATDIM,acesCISevecs(iter),MATDIM,0.0D0,
     &              intermed,MATDIM)

        call dcopy(2*nocc*nvirt,acesCISevecs(iter),1,matmulcopy,1)
        temp=dot_product(matmulcopy, intermed)
        print*, 'Second Way/modifyA:',temp*27.2114
        iter=iter+MATDIM
        enddo

        

        end subroutine
