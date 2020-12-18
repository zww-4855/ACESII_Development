      subroutine owntrans(trans,angmom,centbf,nbasal,nshells,
     & shells,largel,natoms,angmax,nbas)

C----------------------------------------------------------------------
C  Input:
C    angmom: angular momentum of each basis (s=0,p=1,...)
C    centbf: atomic center corresponds to each basis
C    shells(i,a):Ith l-shell type for atom A (s=0,p=1)
C
C  Output: trans
C----------------------------------------------------------------------

      implicit none
 
      integer nbas, i, angmax, angmom(nbas), degsph(angmax),
     & atom, centbf(nbas), nbasal(largel+1,natoms), largel,
     & natoms, io, j, k, nshells(natoms),
     & shells(largel+1,natoms), l, lsh

      double precision trans(nbas,nbas)

      call sphdeg(degsph,angmax)     
C---------------------------------------------------------------------
C  sphdeg: get the degenercy of each angular momentum (s,p,d...)
C          (s=1,p=3...)
C  degsph=1,3,5,7,...
C--------------------------------------------------------------------


      i = 1
      io = i

 10   continue

      if (angmom(i).lt.1) then

         trans(i,i) = 1.0D0
	 i = i + 1
         io = i

      end if

      if (i.gt.nbas) then
         go to 20
      end if
            
      if (angmom(i).ge.1) then

	 atom = centbf(i)
         do l = 1, largel+1
            if (shells(l,atom).eq.angmom(i)) then
               lsh = l
            end if 
         end do

         do j = 1, degsph(angmom(i)+1)
            k = io+nbasal(lsh,atom)*(j-1)
            trans(i,k) = 1.0D0
            i = i + 1
         end do

         io = io + 1
         if (angmom(i-1).ne.angmom(i)) then
            io = i
         end if

      end if

      if (i.gt.nbas) then
         go to 20
      end if
	
      go to 10

 20   continue

      return
      end
