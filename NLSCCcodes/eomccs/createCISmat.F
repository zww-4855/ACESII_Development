        subroutine createCISmat(diaFockA,W,Wab,nocc,nvirt,CISmat)
          integer, intent(in):: nocc,nvirt
          double precision, intent(in)::diaFockA(nocc+nvirt)
          double precision, intent(in)::W(nocc,nvirt,nocc,nvirt)
          double precision, intent(in)::Wab(nvirt,nocc,nocc,nvirt)
          double precision, intent(inout)::CISmat(2*nocc*nvirt,
     &                                          2*nocc*nvirt)

        integer ::del_ij,del_ab,nbas,offset
        integer :: counter_i, counter_j,i,j,a,b
        real(kind=8) :: term,fab, fij,fe,HFenergy
c        CISmat=0.0d0
!       INSERT ALL ALPHA ALPHA ALPHA ALPHA BLOCK
        print*, "PRINTED 2 E- INTS"
        counter_j=1
        do b=1,nvirt
          do j=1,nocc
            counter_i=1
            do a=1,nvirt
              do i=1,nocc
                CISmat(counter_j,counter_i)=CISmat(counter_j,counter_i)
     &                  -W(i,a,j,b)

!                print*, W(i,a,j,b)
                counter_i=counter_i+1
              enddo
            enddo
            counter_j=counter_j+1
            print*
          enddo
        enddo

!       INSERT ALPHA BETA ALHPA BETA BLOCK
c        offset=nocc*nvirt
c        counter_j=1
c        do b=1,nvirt
c          do i=1,nocc
c            counter_i=counter_i+1
c            do j=1,nocc
c              do a=1,nvirt
c               print*, Wab(a,j,i,b)
c               CISmat(counter_j,offset+counter_i)=Wab(a,j,i,b)
c               CISmat(counter_j+offset,counter_i)=Wab(a,j,i,b) 
c               counter_i=counter_i+1
c              enddo
c            enddo
c            counter_j=counter_j+1
c          enddo 
c        enddo
        print*,'***********************************'
        print*, "inside CIS fxn"
        counter_j=1
        offset=nocc*nvirt
        do i=1,nocc
          do a=1,nvirt
            counter_i=1
            do j=1,nocc
              do b=1,nvirt
                 del_ij = (i == j)
                 del_ab = (a == b)
           CISmat(counter_j,counter_i)=CISmat(counter_j,counter_i)+
     &          del_ij*del_ab*(diaFockA(a+nocc)-diaFockA(i))

           CISmat(offset+counter_j,offset+counter_i)=
     &                  CISmat(counter_j,counter_i)

!           CISmat(counter_j,offset+counter_i)=Wab(a,j,i,b)
!           CISmat(offset+counter_j,counter_i)=
!     &                          CISmat(counter_j,offset+counter_i)
                counter_i=counter_i+1
              enddo
            end do
            counter_j=counter_j+1
          enddo
        enddo 

      print*, "symmetric CIS hamil"
      do i=1,2*nocc*nvirt
        do j=1,2*nocc*nvirt
          if (abs(CISmat(i,j)-CISmat(j,i)).gt.1e-5)then
            print*, 'not symmetric at index: ', i,j
          endif
         enddo
      enddo

        end subroutine 
