        subroutine nlsMultiply(acesCISevecs,evecSize,CISmat,nocc,nvirt
     &            ,nbas,roots,NLMOQM1,QM1size,NLMOQM2,QM2size)

        integer, intent(in)::roots,nocc,nvirt,QM1size,QM2size
        integer, intent(in)::NLMOQM1(QM1size),NLMOQM2(QM2size)
        double precision, intent(inout)::acesCISevecs(evecSize),
     &                  CISmat(2*nocc*nvirt,2*nocc*nvirt)

        integer::iter,next,offset,MATDIM,root,i,a
        integer::compareIA(2),QMregIA(2)
        double precision::temp
        double precision::intermed(2*nocc*nvirt)
        double precision::matmulcopy(2*nocc*nvirt)
      !  allocate(intermed(roots*2*nocc*nvirt))
        print*, "INSIDE NLSMULTIPLY"
!#ifdef _DEBUG_LVL0
        print*,'*****************************************************'
        print*,'******           FullCIS matrix                ******'
        print*,'*****************************************************'
        MATDIM=2*nocc*nvirt

        intermed=0.0d0
        iter=1
        next=1
        offset=nocc*nvirt
        do root=1,roots!4
         do i=1,nocc
          do a=nocc+1,nbas
            compareIA=0
            QMregIA=0
            compareIA=(/ i,a /)
            call findQMregion(compareIA,size(compareIA),NLMOQM1,
     &           size(NLMOQM1),NLMOQM2,size(NLMOQM2),QMregIA)
            if (QMregIA(1).eq.2) then
              acesCISevecs(iter)=0.0d0
              acesCISevecs(iter+offset)=0.0d0
            endif
          iter=iter+1
          enddo
        enddo


           call xgemm('N','N',2*nocc*nvirt,1,2*nocc*nvirt,1.0D0,CISmat,
     &          2*nocc*nvirt,acesCISevecs(next),2*nocc*nvirt,0.0D0,
     &              intermed,2*nocc*nvirt)
        call dcopy(2*nocc*nvirt,acesCISevecs(next),1,matmulcopy,1)
        temp=dot_product(matmulcopy, intermed)
        iter=iter+nocc*nvirt
        next=next+2*nocc*nvirt
        print*,'Ajith suggested root', temp*27.2114
        enddo
        end subroutine
