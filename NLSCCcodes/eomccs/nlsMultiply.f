










        subroutine nlsMultiply(acesCISevecs,evecSize,CISmat,nocc,nvirt
     &            ,nbas,roots,NLMOQM1,QM1size,NLMOQM2,QM2size)

        integer, intent(in)::roots,nocc,nvirt,QM1size,QM2size
        integer, intent(in)::NLMOQM1(QM1size),NLMOQM2(QM2size)
        double precision, intent(inout)::acesCISevecs(evecSize),
     &                  CISmat(2*nocc*nvirt,2*nocc*nvirt)

        integer::iter,next,offset,MATDIM,root,i,a
<<<<<<< HEAD
        integer::compareIA(2),QMregIA(2)
        double precision::temp
=======
        integer::compareIA(2),QMregIA(2),test(5)
        double precision::trace,temp
>>>>>>> simplifyCIS
        double precision::intermed(2*nocc*nvirt)
        double precision::matmulcopy(2*nocc*nvirt)
      !  allocate(intermed(roots*2*nocc*nvirt))
        print*, "INSIDE NLSMULTIPLY"
!#ifdef _DEBUG_LVL0
        print*,'*****************************************************'
        print*,'******           FullCIS matrix                ******'
        print*,'*****************************************************'
        MATDIM=2*nocc*nvirt

<<<<<<< HEAD
=======
           call output(acesCISevecs,1,2*nocc*nvirt,1,2*nocc*nvirt,
     &          2*nocc*nvirt,2*nocc*nvirt,1)


        print*,'*****************************************************'
        print*,'*****************************************************'
        print*,'printing NLS-CIS matrix to output file for heatmap...'
        print*,'*****************************************************'
        print*,'*****************************************************'
        print*,'*****************************************************'
        trace=0.0d0
        open(unit=2500, file='NLS_CISmat.txt')
        do i=1,2*nocc*nvirt
          do j=1,2*nocc*nvirt
            if (i.eq.j) then
                trace=trace+CISmat(i,i)
            endif
            write(2500,*) CISmat(j,i)
          enddo
        enddo
        write(2500,*) 'trace: ',trace
        close(2500)

        print*,'*****************************************************'
        print*,'*****************************************************'
        print*,'*****************************************************'
        





        print*,'nbas',nocc,nvirt,nbas
>>>>>>> simplifyCIS
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
