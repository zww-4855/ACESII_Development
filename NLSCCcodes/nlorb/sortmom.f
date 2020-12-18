      subroutine sortmom(centbf,mombf,nshells,shells,nbasal,
     & largel,mxnal,natoms,nbas,angmax,junk,count,count2)

C-------------------------------------------------------------------------
C  Summarize the number of shells on each atom and the number of
C  different type of angular momentum on each shell
C-------------------------------------------------------------------------

      implicit none

      integer angmax 
 
      integer centbf(nbas), mombf(nbas), nshells(natoms),
     & shells(largel+1,natoms), nbas, i, j, k, l,
     & temp, junk(nbas), count(natoms,angmax), degen(angmax),
     & nbasal(largel+1,natoms), mxnal, largel, natoms,
     & count2(natoms,angmax)

C yjin ****************************************
C
C      do i=1,nbas
C       write(*,*) i,centbf(i)
C       write(*,*) i,mombf(i)
C      end do
C
C **********************************************


      mxnal = 0
      temp = 0
      k = 0
      call izerosec(nshells,natoms)
      call izerosec(junk,nbas)
      call sphdeg(degen,angmax) 

      call countmom(natoms,nbas,centbf,mombf,angmax,count) 
      call countsh(natoms,angmax,count,count2,degen) 

      do i = 1, natoms
         do j = 1, nbas 
            if (centbf(j).eq.i) then
               
               junk(j) = mombf(j) + 1

               if (junk(j).ne.temp) then
                  k = k + 1
                  nshells(i) = nshells(i) + 1
                  shells(k,i) = mombf(j)
                  nbasal(k,i) = count2(i,mombf(j) + 1)
                  temp = junk(j)
               end if

            end if  
         end do
 
         do k = 1, nshells(i)
            if (nbasal(k,i).gt.mxnal) then
               mxnal = nbasal(k,i)
            end if
         end do
 
         write(*,10) i, nshells(i) 
         do k = 1, nshells(i)
            write(*,20) k, shells(k,i), nbasal(k,i)
         end do
         write(*,*)

         temp = 0
         k = 0

      end do 
 
      write(*,*)
      write(*,*)

 10   format('Atom ',I5,' has ',I5,' shells')
 20   format('Shell ',I5,' has angular momentum ',I5
     & , ' of dimension ',I5)
 
      return
      end
	  
