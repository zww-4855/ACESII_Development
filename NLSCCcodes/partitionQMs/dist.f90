function dist(coords,lineCount,QM2index,i,j)
     integer,intent(in)::lineCount,i,j
     integer,intent(in)::QM2index(lineCount)
     real, intent(in)::coords(lineCount,3)
     real::temp,dist
!       dist() will have info on the unit vector
!       corresponding to coords(). This should
!       hopefully help by adding a direction into 
!       the mix. SECOND ELEMENT IS UNIT VECTOR IN
!       X DIRECTION 

!       ELEMENTS OF dist() CORRESPOND TO ELEMENTS OF
!       QM2index() MEANING THEY ARE EITHER IN QM1/2
        p=QM2index(i)
        q=QM2index(j)
        
        temp=(coords(p,1)-coords(q,1))**2 + &
           & (coords(p,2)-coords(q,2))**2 + &
           & (coords(p,3)-coords(q,3))**2 

        dist=sqrt(temp)

end function dist
