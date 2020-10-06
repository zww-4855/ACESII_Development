        subroutine ReduceCISmat(CISmat0,CIS0dim,
     &                  NLMOQM1,NLMOQM2,nbas,nocc,nvirt)

        integer, intent(in)::CIS0dim,nbas,nocc,nvirt
        integer,intent(inout)::NLMOQM1(nbas),NLMOQM2(nbas)
        double precision, intent(in)::CISmat0(CIS0dim,CIS0dim)
        double precision,allocatable::CISreduce(:,:),
     &                                          CISreduceVec(:,:)
        double precision,allocatable:: NLScistemp(:)
        integer :: QMreg(4),compareIA(2),compareJB(2),compareIJAB(4)
        integer :: iterR,iterC,counter_i, counter_j,i,j,a,b
        integer :: QMregIA(2),QMregJB(2)
        integer:: colZEROindx(CIS0dim),rowZEROindx(CIS0dim),ia,jb
        logical:: colCount(CIS0dim),rowCount(CIS0dim)
        iterR=1
        iterC=1
        rowZEROindx=0
        colZEROindx=0
        do i=1,CIS0dim
          colCount=.False.
          rowCount=.False.
          do j=1,CIS0dim
             if (abs(CISmat0(i,j)).lt.0.000000001) colCount(j)=.True.
             if (abs(CISmat0(j,i)).lt.0.000000001) rowCount(j)=.True.
          enddo   
          if (all(colCount)) then !.and. all(rowCount)) then
                colZEROindx(iterC)=i
                iterC=iterC+1
          endif
          if (all(rowCount)) then
                rowZEROindx(iterR)=i
                iterR=iterR+1
          endif
        enddo


        ndim=2*nocc*nvirt-iterR
        print*, 'row dim and col dim for NLS-CIS matrix:', iterR,iterC
        print*, 'row iter:', rowZEROindx
        print*
        print*, 'col iter:', colZEROindx
        print*, 'size of NLS-CIS matrix:', ndim
        do i=1,iterR!CIS0dim
           if ((.not.any(rowZEROindx(i).eq.colZEROindx))) then
                print*,'i found where rowscols are not equivalent',i
                print*, rowZEROindx(i)
           endif
        enddo
        allocate(CISreduce(ndim,ndim),CISreduceVec(ndim,ndim),
     &                  NLScistemp(ndim*ndim))
        CISreduce=0.0d0
c        coi=1
c        do ia=1,CIS0dim
c          if (any(ia.eq.rowZEROindx)) cycle
c          coj=1
c          do jb=i,CIS0dim
c            if (any(jb.eq.colZEROindx)) cycle
c            CISreduce(coi,coj)=CISmat0(ia,jb)
c            CISreduce(coj,coi)=CISmat0(jb,ia)
c            coj=coj+1
c          enddo
c          coi=coi+1
c        enddo

        coi=1
        do ia=1,CIS0dim
          if (any(ia.eq.colZEROindx)) cycle
          coj=1
          do jb=1,CIS0dim
            if (any(jb.eq.rowZEROindx)) cycle
            CISreduce(coi,coj)=CISmat0(ia,jb)
            coj=coj+1
          enddo 
          coi=coi+1
        enddo

        do j=1,ndim
          CISreduce(1,j)=CISreduce(j,1)
        enddo

c        num=1
c        do i=1,CIS0dim
c          if (any(i.eq.colZEROindx)) cycle
c          do j=1,CIS0dim
c            if (any(j.eq.rowZEROindx)) cycle
c            NLScistemp(num)=CISmat0(i,j)
c            num=num+1
c          enddo
c        enddo
!        print*,NLScistemp(2),NLScistemp(2+ndim)
c#ifdef _DEBUG_LVL0
        print*, '** Reduced NLS-CIS matrix **'        
        print*, CISmat0(17,17),CISmat0(18,18)
        print*, CISmat0(18,17),CISmat0(17,18)
        print*, CISmat0(19,17),CISmat0(17,19)
        print*,CISmat0(18,19),CISmat0(19,18)
!        call output(NLScistemp,1,4,1,4,ndim,ndim,1)
        call output(CISreduce,1,4,1,4,ndim,ndim,1)!ndim,1,ndim,ndim,ndim,1)
         call output(CISmat0,17,21,17,21,2*nocc*nvirt,2*nocc*nvirt,1)
c#endif
!*****************************************************
!       *Check is matrix is symmetric*
!*****************************************************
        do i=1,ndim
          do j=1,ndim
            if (abs(CISreduce(i,j)-CISreduce(j,i)).gt.1e-5) then
                print*,'mat not symm at index', i,j
             endif
        enddo
        enddo
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
c        call output(CISreduceVec,1,ndim,1,ndim,ndim,ndim,1)
#endif

        print*,'*********************************************'
        print*,'FINAL NLS-CIS EIGENVECTORS FOR FIRST FEW STATES'
        print*,'*********************************************'
        PRINT*, CISreduceVec(:,newindx)
        print*,'max vec',maxval(CISreduceVec(:,newindx))
        print*,'min vec', minval(CISreduceVec(:,newindx))        
        deallocate(CISreduce,CISreduceVec,NLScistemp)
        end subroutine
