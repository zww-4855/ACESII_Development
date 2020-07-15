subroutine addCaps(coords,bondMat,atomName,hcapList,lineCount,&
                 & revQMcoords,revAtomName,QM2Index,QM1Index,QM1count,scaleH)
                
        integer,intent(in)::lineCount
        real,intent(in)::scaleH
        real,intent(inout):: coords(lineCount,3)
        integer,intent(in)::QM2Index(lineCount)
        integer,intent(in)::bondMat(lineCount,4)
        integer,intent(in)::hcapList(lineCount,4),QM1count
        integer,intent(in)::QM1Index(QM1count)
        real,intent(inout)::revQMcoords(3)
        character (len=*),intent(in)::atomName(lineCount)
        character (len=*),intent(inout)::revAtomName        
        integer::tempArr(4),i,j,ix
        integer:: a,revQM2Index(lineCount),QM2realCount
        integer,allocatable::cutIndex(:)
        revQM2Index=0
! Write all of QM1 to toZMAT first, then QM2 after
        open(1555,file="toZMAT")
        do i=1,QM1count
           write(1555,*) atomName(QM1Index(i)), coords(QM1Index(i),1),&
                         & coords(QM1Index(i),2),coords(QM1Index(i),3)          
        enddo

        QM2realCount=0
        do i=1,lineCount
           if (QM2Index(i) .eq. 0) cycle
           if (any(QM2Index(i) .eq. QM1Index)) then
                cycle
           endif
           write(1555,*) atomName(QM2Index(i)), coords(QM2Index(i),1),&
                         & coords(QM2Index(i),2),coords(QM2Index(i),3)

           QM2realCount=QM2realCount+1
        enddo

        print*,"number of sites just in QM1:", QM1count
        print*,'Number of sites in both QM2/QM1 ',QM1count+QM2realCount
        print*,'Number of sites just in QM2, *not* QM1,w/o caps',QM2realCount


! info on pack:
! https://stackoverflow.com/questions/16588214/is-there-an-easy-way-to-find-index-array-zeros-in-fortran
! cutIndex contains atomIndexs outside QM1/2
        print*,'hcap list',hcapList
        do i=1,lineCount
          if (any(hcapList(i,:).ne.0)) then 
              tempArr=hcapList(i,:)
              print*,'tempArr',tempArr
              a=size(pack([(ix,ix=1,4)],tempArr.ne.0))
              print*,'a is:', a
              print*
              allocate(cutIndex(a))

              cutIndex=(/(0, i=1,a)/)
              print*,'cut index ad=bd temp array before:', tempArr,cutIndex
              cutIndex=pack([(ix,ix=1,4)],tempArr.ne.0)
              print*,'cutIndex',cutIndex
              print*,'pack is', pack([(ix,ix=1,4)],tempArr.ne.0)
              print*
              do j=1,a
                 if (cutIndex(j).eq.0) cycle
                 revQMcoords=0.0d0
                 print*,'atomQM -> atomOut', i, tempArr(cutIndex(j)) 
                 call makeUnitVector(coords,lineCount,i,&
                                & tempArr(cutIndex(j)),revQMcoords)
                 revQMcoords=scaleH*revQMcoords
                 QM2realCount=QM2realCount+1
                 write(1555,*) 'H',revQMcoords(1),revQMcoords(2),revQMcoords(3)

              enddo
              deallocate(cutIndex)
  
          endif
        enddo
        close(1555)
! Write to QMcenters
        open(2555,file='QMcenters')
        write(2555,*) QM1count
        write(2555,*)
        do i=1,QM1count
                write(2555,*) i
        enddo   
        write(2555,*)
        write(2555,*) QM2realCount-QM1count
        write(2555,*)
        do i=QM1count+1,QM2realCount
           ! if (any(QM2Index(i) .eq. QM1Index)) cycle
           ! if (revQM2Index(i).eq.0) cycle 
            write(2555,*) i
        enddo
        write(2555,*)

        close(2555)

        close(2555)
end subroutine

