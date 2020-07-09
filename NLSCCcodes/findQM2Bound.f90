subroutine findQM2Bound(coords,lineCount,QM1coords,QM1count,QMcutoff,&
        &      QM1index,xmin,xmax,ymin,ymax,zmin,zmax,QM2index)

        integer, intent(in)::lineCount,QM1count,xmin,xmax
        integer, intent(in)::ymin,ymax,zmin,zmax,QM1index(QM1count)
        real,intent(in)::coords(lineCount,3),QM1coords(QM1count,3)
        real,intent(in):: QMcutoff
        integer,intent(out):: QM2index(lineCount)
        real::xminCut,xmaxCut,yminCut,ymaxCut,zminCut,zmaxCut
        integer::counter,i

        xminCut=coords(xmin,1)-QMcutoff
        xmaxCut=coords(xmax,1)+QMcutoff

        yminCut=coords(ymin,2)-QMcutoff
        ymaxCut=coords(ymax,2)+QMcutoff

        zminCut=coords(zmin,3)-QMcutoff
        zmaxCut=coords(zmax,3)+QMcutoff

        print*, 'xmin end of QM1: ', coords(xmin,1)
        print*, 'xmax end of QM1: ', coords(xmax,1)

        print*, 'ymin end of QM1: ', coords(ymin,2)
        print*, 'ymax end of QM1: ', coords(ymax,2)

        print*, 'zmin end of QM1: ', coords(zmin,3)
        print*, 'zmax end of QM1: ', coords(zmax,3)

        print*,'x cutoffs: ', xminCut,xmaxCut
        print*,'y cutoffs: ', yminCut,ymaxCut
        print*, 'z cutoffs: ',zminCut,zmaxCut
        counter=1
        do i=1,lineCount
           if (any(QM1index==i)) Continue
           print*, 'i is: ', i   
           if ((coords(i,1).gt.xminCut) .and. (coords(i,1).lt.xmaxCut)&
            & .and. (coords(i,2).gt.yminCut).and.(coords(i,2).lt.ymaxCut)&
            & .and.(coords(i,3).gt.zminCut).and.(coords(i,3).lt.zmaxCut)) then

                QM2index(counter)=i
                counter=counter+1
           endif
        enddo
        print*, 'QM2 region is: ', QM2index
end subroutine findQM2Bound
