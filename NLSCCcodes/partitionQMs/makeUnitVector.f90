subroutine makeUnitVector(coords,lineCount,inAtom,outAtom,revQMcoords)
        integer,intent(in)::inAtom,outAtom,lineCount
        real,intent(in):: coords(lineCount,3)
        real,intent(inout)::revQMcoords(3)
        real::norm,xdist,ydist,zdist

        print*,'distance from inAtom->outAtom',inAtom,outAtom
        xdist=coords(outAtom,1)-coords(inAtom,1)
        ydist=coords(outAtom,2)-coords(inAtom,2)
        zdist=coords(outAtom,3)-coords(inAtom,3)

        norm=sqrt(xdist**2 + ydist**2 + zdist**2)
        print*,'distances',xdist,ydist,zdist,norm
        revQMcoords(1)= (1.0d0/norm)*xdist + coords(inAtom,1)
        revQMcoords(2)=(1.0d0/norm)*ydist +coords(inAtom,2)
        revQMcoords(3)=(1.0d0/norm)*zdist +coords(inAtom,3)
        
end subroutine
