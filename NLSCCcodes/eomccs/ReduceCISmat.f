










        subroutine ReduceCISmat(CISmat0,CIS0dim,
     &                  NLMOQM1,NLMOQM2,nbas,nocc,nvirt)

        integer, intent(in)::CIS0dim,nbas,nocc,nvirt
        integer,intent(inout)::NLMOQM1(nbas),NLMOQM2(nbas)
        double precision, intent(in)::CISmat0(CIS0dim,CIS0dim)
        double precision,allocatable::CISreduce(:,:),
     &                                          CISreduceVec(:,:)
        integer :: QMreg(4),compareIA(2),compareJB(2),compareIJAB(4)
        integer :: counter_i, counter_j,i,j,a,b
        integer :: QMregIA(2),QMregJB(2)
        integer::ZEROindx(CIS0dim),ia,jb

        iter=1
        do i=1,CIS0dim
          colCount=sum(CISmat0(i,:))
          rowCount=sum(CISmat0(:,i))
          if ((abs(colCount).lt.0.0000000001) .and.
     &                  (abs(rowCount).lt.0.0000000001)) then
                ZEROindx(iter)=i
                iter=iter+1
          endif
        enddo


        ndim=2*nocc*nvirt-iter
        allocate(CISreduce(ndim,ndim),CISreduceVec(ndim,ndim))
        coi=1
        do ia=1,CIS0dim
          if (any(ia.eq.ZEROindx)) cycle
          coj=1
          do jb=1,CIS0dim
            if (any(jb.eq.ZEROindx)) cycle
            CISreduce(coi,coj)=CISmat0(ia,jb)
            print*,CISmat0(ia,jb)
            coj=coj+1
          enddo 
          coi=coi+1
        enddo
        print*, 'revised reduce CIS mat'        
        call output(CISreduce,1,ndim,1,ndim,ndim,ndim,1)
         call eig(CISreduce,CISreduceVec,100,ndim,1)
        do i=1,ndim
          print*, CISreduce(i,i),CISreduce(i,i)*27.2114
        enddo

        call output(CISreduceVec,1,ndim,1,ndim,ndim,ndim,1)

        deallocate(CISreduce,CISreduceVec)
        end subroutine
