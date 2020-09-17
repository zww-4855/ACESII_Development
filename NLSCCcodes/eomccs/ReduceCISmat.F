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
          if ((abs(colCount).lt.0.0000001) .and.
     &                  (abs(rowCount).lt.0.0000001)) then
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
!            print*,CISmat0(ia,jb)
            coj=coj+1
          enddo 
          coi=coi+1
        enddo

#ifdef _DEBUG_LVL0
        print*, '** Reduced NLS-CIS matrix **'        
        call output(CISreduce,1,ndim,1,ndim,ndim,ndim,1)
#endif
         call eig(CISreduce,CISreduceVec,100,ndim,1)
        print*,'*********************************************'
        print*,'FINAL NLS-CIS EIGENVALUES'
        print*,'*********************************************'
        print*, 'original # of roots for full CIS',CIS0dim 
        print*,'NLS-CIS # of roots: ', ndim
        newindx=1
        do i=1,ndim
          print*, i,CISreduce(i,i)*27.2114
          if (abs(CISreduce(i,i)).lt.0.0000001) then
                newindx=newindx+1
          endif
        enddo
        print*,'*********************************************'
        print*,'*********************************************'
        print*,'*********************************************'
#ifdef _DEBUG_LVL0
        call output(CISreduceVec,1,ndim,1,ndim,ndim,ndim,1)
#endif

        print*,'*********************************************'
        print*,'FINAL NLS-CIS EIGENVECTORS FOR FIRST FEW STATES'
        print*,'*********************************************'
        PRINT*, CISreduceVec(:,newindx)
        print*,'max vec',maxval(CISreduceVec(:,newindx))
        print*,'min vec', minval(CISreduceVec(:,newindx))        
        deallocate(CISreduce,CISreduceVec)
        end subroutine
