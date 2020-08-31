        function findQMregion(indx,indxSize,QMregion,QMregSize)!,output)
           integer,intent(in)::indxSize,QMregSize
           integer,intent(in)::indx(indxSize),QMregion(QMregSize)
           !logical, intent(inout)::findQMregion 
           integer:: i,j
           logical:: finalOut(indxSize) 
           logical:: test(indxSize)
           logical::findQMregion
           test=(/ (.TRUE.,j=1,indxSize)  /)
          do i=1,indxSize
                if (any(indx(i).eq.QMregion)) then
                   finalOut(i)=.TRUE.
                else
                   finalOut(i)=.FALSE.
                endif
          enddo
          print*,'**************************************'
          print*,'** INside findQM  **'
           print*,'test is: ', test
           print*,'compare is: ', finalOut
           print*,'indx is: ', indx
          findQMregion=all(finalOut .eq. test,dim=1) 
          print*, 'value of findQMreigon:', findQMregion
          print*,'**************************************'
          print*,'**************************************'
        end function 
