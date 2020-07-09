function distance(coords,lineCount,i,j)
     integer,intent(in)::lineCount,i,j
     real, intent(in)::coords(lineCount,3)
     real::temp,distance
!       distance() will have info on the unit vector
!       corresponding to coords(). This should
!       hopefully help by adding a direction into 
!       the mix. SECOND ELEMENT IS UNIT VECTOR IN
!       X DIRECTION 

!       ELEMENTS OF distance() CORRESPOND TO ELEMENTS OF
!       QM2index() MEANING THEY ARE EITHER IN QM1/2
        
        temp=(coords(i,1)-coords(j,1))**2 + &
           & (coords(i,2)-coords(j,2))**2 + &
           & (coords(i,3)-coords(j,3))**2 
        print*,'all of coords in dist'
        print*,coords(i,:),coords(j,:)
        distance=sqrt(temp)

end function distance
