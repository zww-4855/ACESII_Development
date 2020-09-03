        subroutine nlscisZ(CISmat,CISmat0,nocc,nvirt,QM1atoms,QM1num,
     &                     QM2atoms,QM2num,NLMOQM1,NLMOQM2,nbas)
          integer, intent(in):: nocc,nvirt,QM1num,QM2num,nbas
          integer, intent(in)::NLMOQM1(nbas), NLMOQM2(nbas)
          integer, intent(in):: QM1atoms(QM1num),QM2atoms(QM2num)
          double precision, intent(in)::CISmat(2*nocc*nvirt,
     &                                          2*nocc*nvirt)
          double precision, intent(inout)::CISmat0(2*nocc*nvirt,
     &                                          2*nocc*nvirt)

        integer :: QMreg(4),compareIJAB(4),compareIJ(2),compareAB(2)
        integer ::iter,del_ij,del_ab,offset
        integer :: counter_i, counter_j,i,j,a,b
!        real(kind=8) :: term,fab, fij,fe,HFenergy
        print*,'inside nlscisZERO'
        CISmat0=CISmat !0.0d0
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
                QMreg=0
                compareIJAB=(/ a,j,i,b /)
        
                call findQMregion(compareIJAB,size(compareIJAB),NLMOQM1,
     &                   size(NLMOQM1),NLMOQM2,size(NLMOQM2),QMreg)
                print*,'value of whether all indices are in QM1'
                print*,i,j
                print*,a,b
                print*,all(QMreg.eq.(/ 1,1,1,1 /))
                print*,QMreg
                print*
                print*,NLMOQM1
                print*
                print*,CISmat(counter_j,counter_i)
                print*,'end'
                print*
                if (all(QMreg.eq.(/ 1,1,1,1 /))) then
                   CISmat0(counter_j,counter_i)=
     &                  CISmat(counter_j,counter_i)

                   CISmat0(offset+counter_j,counter_i)=
     &                  CISmat(offset+counter_j,counter_i)

                   CISmat0(counter_j,offset+counter_i)=
     &                  CISmat(counter_j,offset+counter_i)

                   CISmat0(offset+counter_j,offset+counter_i)=
     &                  CISmat(offset+counter_j,offset+counter_i)

                else if (any(QMreg.eq.2)) then
                   CISmat0(counter_j,counter_i)=0.0d0

                   CISmat0(offset+counter_j,counter_i)=0.0d0

                   CISmat0(counter_j,offset+counter_i)=0.0d0

                   CISmat0(offset+counter_j,offset+counter_i)=0.0d0
                else if (any(QMreg.eq.3)) then
                   if (sum(QMreg).eq.6) then! ie <31||11>
                   CISmat0(counter_j,counter_i)=
     &                  CISmat(counter_j,counter_i)!*0.75d0

                   CISmat0(offset+counter_j,counter_i)=
     &                  CISmat(offset+counter_j,counter_i)!*0.75d0

                   CISmat0(counter_j,offset+counter_i)=
     &                  CISmat(counter_j,offset+counter_i)!*0.75d0

                   CISmat0(offset+counter_j,offset+counter_i)=
     &                  CISmat(offset+counter_j,offset+counter_i)!*0.75d0
               
                    else if (sum(QMreg).eq.8) then! ie <33||11>
                   CISmat0(counter_j,counter_i)=
     &                  CISmat(counter_j,counter_i)!*0.5d0

                   CISmat0(offset+counter_j,counter_i)=
     &                  CISmat(offset+counter_j,counter_i)!*0.5d0

                   CISmat0(counter_j,offset+counter_i)=
     &                  CISmat(counter_j,offset+counter_i)!*0.5d0

                   CISmat0(offset+counter_j,offset+counter_i)=
     &                  CISmat(offset+counter_j,offset+counter_i)!*0.5d0
 

                    else if (sum(QMreg).eq.10) then! ie <33||31>
                   CISmat0(counter_j,counter_i)=
     &                  CISmat(counter_j,counter_i)!*0.25d0

                   CISmat0(offset+counter_j,counter_i)=
     &                  CISmat(offset+counter_j,counter_i)!*0.25d0

                   CISmat0(counter_j,offset+counter_i)=
     &                  CISmat(counter_j,offset+counter_i)!*0.25d0

                   CISmat0(offset+counter_j,offset+counter_i)=
     &                  CISmat(offset+counter_j,offset+counter_i)!*0.25d0
                    endif
!                         all(QMreg.eq.(/ 2,2,2,2 /))
                else
                   CISmat0(counter_j,counter_i)=0.0d0

                   CISmat0(offset+counter_j,counter_i)=0.0d0

                   CISmat0(counter_j,offset+counter_i)=0.0d0

                   CISmat0(offset+counter_j,offset+counter_i)=0.0d0
                endif
        

!        elseif((findQMregion(compareIJ,size(compareIJ),NLMOQM1,
!     &                   size(NLMOQM1))) .and.
!     &          (findQMregion(compareAB,size(compareAB),NLMOQM2,
!     &                  size(NLMOQM2)))) then
!                print*,"IJ in QM1 and AB in QM2"
!                CISmat0(counter_j,counter_i)=
!     &                  CISmat(counter_j,counter_i)
!
!                CISmat0(offset+counter_j,counter_i)=
!     &                  CISmat(offset+counter_j,counter_i)
!
!                CISmat0(counter_j,offset+counter_i)=
!     &                  CISmat(counter_j,offset+counter_i)
!
!                CISmat0(offset+counter_j,offset+counter_i)=
!     &                  CISmat(offset+counter_j,offset+counter_i)
!         else
!                print*,"no indexes agree"
!                CISmat0(counter_j,counter_i)=0.0d0
!                CISmat0(offset+counter_j,counter_i)=0.0d0
!                CISmat0(counter_j,offset+counter_i)=0.0d0
!                CISmat0(offset+counter_j,offset+counter_i)=0.0d0
!         endif 

                counter_i=counter_i+1
              enddo
            end do
            counter_j=counter_j+1
          enddo
        enddo

        do i=1,2*nocc*nvirt
          do j=1,2*nocc*nvirt
            if (abs(CISmat0(i,j)-CISmat0(j,i)).gt.1e-5) then
                print*,'mat not symm at index', i,j
             endif
        enddo
        enddo
        end subroutine
