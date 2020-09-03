        subroutine countQM2(NLMOQM1,NLMOQM2,qmsize,QM2list)
        integer, intent(in)::qmsize
        integer,intent(in)::NLMOQM1(qmsize),NLMOQM2(qmsize)
        integer,intent(inout)::QM2list(2*qmsize)
        integer::counter1,i,QMcount
        QMcount=0
        counter1=1
        print*, 'inside ocuntQM2'
        print*,'nlmo1', NLMOQM1
        print*
        print*,'nlmo2',NLMOQM2
        do i=1,qmsize
          if (not(any(NLMOQM2(i).eq.NLMOQM1))) then
            QM2list(counter1)=NLMOQM2(i)
            counter1=counter1+1
            QMcount=QMcount+1
          endif
        enddo
        
        print*, QM2list
        do i=1,QMcount
          QM2list(i+QMcount)=QM2list(i)+nocc
        enddo
        end subroutine        
