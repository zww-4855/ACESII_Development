










        subroutine nlscisZERO(CISmat,CISmat0,nocc,nvirt,QM1atoms,QM1num,
     &                     QM2atoms,QM2num)
          integer, intent(in):: nocc,nvirt,QM1num,QM2num
          integer, intent(in):: QM1atoms(QM1num),QM2atoms(QM2num)
          double precision, intent(in)::CISmat(2*nocc*nvirt,
     &                                          2*nocc*nvirt)
          double precision, intent(inout)::CISmat0(2*nocc*nvirt,
     &                                          2*nocc*nvirt)

        integer :: compareIJAB(4)
        integer ::iter,del_ij,del_ab,nbas,offset
        integer :: counter_i, counter_j,i,j,a,b
!        real(kind=8) :: term,fab, fij,fe,HFenergy

        CISmat0=0.0d0
        counter_j=1
        offset=nocc*nvirt
        do i=1,nocc
          do a=1,nvirt
            counter_i=1
            do j=1,nocc
              do b=1,nvirt
                compareIJAB=0
                compareIJAB=(/ i,j,a,b /)
!                print*, 'compare IJAB', compareIJAB

                counter_i=counter_i+1
              enddo
            end do
            counter_j=counter_j+1
          enddo
        enddo
        end subroutine
