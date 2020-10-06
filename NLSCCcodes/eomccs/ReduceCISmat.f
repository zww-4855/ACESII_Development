










        subroutine ReduceCISmat(CISmat0,CIS0dim,
     &                  NLMOQM1,NLMOQM2,nbas,nocc,nvirt)

        integer, intent(in)::CIS0dim,nbas,nocc,nvirt
        integer,intent(inout)::NLMOQM1(nbas),NLMOQM2(nbas)
        double precision, intent(in)::CISmat0(CIS0dim,CIS0dim)
        double precision,allocatable::CISreduce(:,:),
     &                                          CISreduceVec(:,:)
        integer :: QMreg(4),compareIA(2),compareJB(2),compareIJAB(4)
        integer :: iterR,iterC,counter_i, counter_j,i,j,a,b
        integer :: QMregIA(2),QMregJB(2)
        integer:: colZEROindx(CIS0dim),rowZEROindx(CIS0dim),ia,jb
        logical:: colCount(CIS0dim),rowCount(CIS0dim)
        logical :: keptIndx(CIS0dim)
        integer,allocatable::selectedIndx(:)

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


        ndim=CIS0dim-iterR+1 
        do i=1,iterR
           if ((.not.any(rowZEROindx(i).eq.colZEROindx))) then
                print*,'i found where rowscols are not equivalent',i
                print*, rowZEROindx(i)
           endif
        enddo
        allocate(CISreduce(ndim,ndim),CISreduceVec(ndim,ndim),
     &                  selectedIndx(ndim))
        CISreduce=0.0d0
        selectedIndx=0
        CISreduceVec=0.0d0
c        coi=1
c        do ia=1,CIS0dim
c          if (any(ia.eq.colZEROindx)) cycle
c          coj=1
c          do jb=1,CIS0dim
c            if (any(jb.eq.rowZEROindx)) cycle
c            CISreduce(coi,coj)=CISmat0(ia,jb)
c            coj=coj+1
c          enddo 
c          coi=coi+1
c        enddo
        
c Source:
c  https://stackoverflow.com/questions/48363576/how-to-remove-several-columns-from-a-matrix 

        keptIndx=.true.
        keptIndx(rowZEROindx(1:iterR))=.FALSE.
        selectedIndx=PACK([(i,i=1,CIS0dim)],keptIndx)
        print*, 'selected index', selectedIndx
        CISreduce=CISmat0(selectedIndx,selectedIndx)



!*****************************************************
!       *Diagonalize NLS-CIS matrix*
!*****************************************************
        print*,'*********************************************'
        print*,'FINAL NLS-CIS EIGENVALUES'
        print*,'*********************************************'
        print*
        call eig(CISreduce,CISreduceVec,100,ndim,1)
        print*, '** first five roots **'
        print*, '               Root #  Excitation Energy (eV)'
        do i=1,5!ndim
          print*, i, CISreduce(i,i)*27.2114
        enddo
        print*, 'original # of roots for full CIS',CIS0dim 
        print*,'NLS-CIS # of roots: ', ndim
        newindx=1
        print*,'*********************************************'
        print*,'*********************************************'
        print*,'*********************************************'

        print*,'*********************************************'
        print*,'FINAL NLS-CIS EIGENVECTORS FOR FIRST FEW STATES'
        print*,'*********************************************'
        PRINT*, CISreduceVec(:,newindx)
        print*,'max vec',maxval(CISreduceVec(:,newindx))
        print*,'min vec', minval(CISreduceVec(:,newindx))        
        deallocate(CISreduce,CISreduceVec,selectedIndx)
        end subroutine
