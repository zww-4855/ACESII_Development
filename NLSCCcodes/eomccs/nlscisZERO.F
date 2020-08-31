        subroutine nlscisZ(CISmat,CISmat0,nocc,nvirt,QM1atoms,QM1num,
     &                     QM2atoms,QM2num,NLMOQM1,NLMOQM2,nbas)
          integer, intent(in):: nocc,nvirt,QM1num,QM2num,nbas
          integer, intent(in)::NLMOQM1(nbas), NLMOQM2(nbas)
          integer, intent(in):: QM1atoms(QM1num),QM2atoms(QM2num)
          double precision, intent(in)::CISmat(2*nocc*nvirt,
     &                                          2*nocc*nvirt)
          double precision, intent(inout)::CISmat0(2*nocc*nvirt,
     &                                          2*nocc*nvirt)

        integer :: compareIJAB(4),compareIJ(2),compareAB(2)
        integer ::iter,del_ij,del_ab,offset
        integer :: counter_i, counter_j,i,j,a,b
        logical :: findQMregion
!        real(kind=8) :: term,fab, fij,fe,HFenergy
        print*,'inside nlscisZERO'
        CISmat0=0.0d0
        counter_j=1
        offset=nocc*nvirt
        loca=findloc(NLMOQM1,value=0,dim=1)
        locb=findloc(NLMOQM2,value=0,dim=1)
        loca=loca-1
        locb=locb-1
        print*
        print*
        print*,'nlmoQM1&2', NLMOQM1,NLMOQM2
        do i=1,nocc
          do a=nocc+1,nbas
            counter_i=1
            do j=1,nocc
              do b=nocc+1,nbas
                compareIJAB=0
                compareIJAB=(/ i,j,a,b /)
                compareIJ=(/ i,j /)
                compareAB=(/ a,b /)
! first QM1 -> QM1 ie i,j,a,b in QM1
        if(findQMregion(compareIJAB,size(compareIJAB),NLMOQM1,
     &                   size(NLMOQM1))) then
                print*,'IJAB in QM1 fully' 
                CISmat0(counter_j,counter_i)=
     &                  CISmat(counter_j,counter_i)

                CISmat0(offset+counter_j,counter_i)=
     &                  CISmat(offset+counter_j,counter_i)

                CISmat0(counter_j,offset+counter_i)=
     &                  CISmat(counter_j,offset+counter_i)

                CISmat0(offset+counter_j,offset+counter_i)=
     &                  CISmat(offset+counter_j,offset+counter_i)

        elseif((findQMregion(compareIJ,size(compareIJ),NLMOQM1,
     &                   size(NLMOQM1))) .and.
     &          (findQMregion(compareAB,size(compareAB),NLMOQM2,
     &                  size(NLMOQM2)))) then
                print*,"IJ in QM1 and AB in QM2"
                CISmat0(counter_j,counter_i)=
     &                  CISmat(counter_j,counter_i)

                CISmat0(offset+counter_j,counter_i)=
     &                  CISmat(offset+counter_j,counter_i)

                CISmat0(counter_j,offset+counter_i)=
     &                  CISmat(counter_j,offset+counter_i)

                CISmat0(offset+counter_j,offset+counter_i)=
     &                  CISmat(offset+counter_j,offset+counter_i)
         else
                print*,"no indexes agree"
                CISmat0(counter_j,counter_i)=0.0d0
                CISmat0(offset+counter_j,counter_i)=0.0d0
                CISmat0(counter_j,offset+counter_i)=0.0d0
                CISmat0(offset+counter_j,offset+counter_i)=0.0d0
         endif 

                counter_i=counter_i+1
              enddo
            end do
            counter_j=counter_j+1
          enddo
        enddo
        end subroutine
