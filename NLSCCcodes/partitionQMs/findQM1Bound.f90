subroutine findQM1Bound(coords,QMindex,line,xmin,xmax,ymin,ymax,zmin,zmax)
        real,intent(in)::coords(line,3)
        integer,intent(in)::line,QMindex(line)
        integer,intent(inout)::xmin,xmax,ymin,ymax,zmin,zmax
        real::x1,x2,y1,y2,z1,z2

        x1=coords(1,1)
        x2=coords(1,1)
        y1=coords(1,2)
        y2=coords(1,2)
        z1=coords(1,3)
        z2=coords(1,3)

        xmin=QMindex(1)
        xmax=QMindex(1)
        ymin=QMindex(1)
        ymax=QMindex(1)
        zmin=QMindex(1)
        zmax=QMindex(1)

        print*, "qm index in qm1bound", QMindex
        do i=1,line
          if (coords(i,1).lt.x1) then
             xmin=QMindex(i)
             x1=coords(i,1)
          endif
          if (coords(i,1).gt.x2) then
             xmax=QMindex(i) 
             x2=coords(i,1)
          endif          

          if (coords(i,2).lt.y1) then
             ymin=QMindex(i)
             y1=coords(i,2)
          endif
          if (coords(i,2).gt.y2) then
             ymax=QMindex(i)
             y2=coords(i,2)
            
          endif

          if (coords(i,3).lt.z1) then
             zmin=QMindex(i)
             z1=coords(i,3)
          endif
          if (coords(i,3).gt.z2) then
             zmax=QMindex(i)
             z2=coords(i,3)
          endif
        enddo

        print*, xmin,xmax,ymin,ymax,zmin,zmax
end subroutine findQM1Bound 
