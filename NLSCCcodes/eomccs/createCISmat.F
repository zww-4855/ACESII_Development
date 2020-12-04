        subroutine createCISmat(diaFockA,W,Wab,nocc,nvirt,CISmat)
          integer, intent(in):: nocc,nvirt
          double precision, intent(in)::diaFockA(nocc+nvirt)
          double precision, intent(in)::W(nocc,nvirt,nocc,nvirt)
          double precision, intent(in)::Wab(nvirt*nocc*nocc*nvirt)
          double precision, intent(inout)::CISmat(2*nocc*nvirt,
     &                                          2*nocc*nvirt)

        integer ::iter,del_ij,del_ab,nbas,offset
        integer :: counter_i, counter_j,i,j,a,b
        real(kind=8) :: term,fab, fij,fe,HFenergy
! ****************************************************************
! ****************************************************************
!     * This subroutine handles Diagonal elements of CIS matrix *
!  * All alpha/alpha/alpha/alpha block and all Beta-Beta-beta-beta block
!   
!             A_ia,jb = delij*delab*(e_a-e_i) + <AJ||IB>+2Fia for all alpha
!             A_ia,jb = delij*delab*(e_a-e_i) + <aj||ib>+2Fia for all beta
!
!   ** NOTE::: Off diagonal alpha/beta/alpha/beta block and
!              beta/alpha/beta/alpha block are setup separately.
!       ie A_ia,jb = <Aj|Ib>
!                  handled, and stored in the CIS matrix 'CISmat'
! ****************************************************************
! ****************************************************************
!
!
! ****************************************************************
! ****************************************************************
!       Step0) Insert mixed alpha/beta/alpha/beta &
!              beta/alpha/beta/alpha blocks first. 
! ****************************************************************
! ****************************************************************
        print*,'Waa integrals',W
        print*, 'Wab integrals', Wab
        CISmat=0.0d0
        offset=nocc*nvirt
        iter=1
        do i=1,nocc*nvirt
          do a=1,nocc*nvirt
!           print*,Wab(iter)
           CISmat(i+offset,a)=Wab(iter)
           CISmat(i,offset+a)=Wab(iter)
           iter=iter+1
          enddo
        enddo

! ****************************************************************
! ****************************************************************
!       Step1) Insert the 2e- integral part from record 23 from its
!              loaded format. 
!               
!               * Only inserts <aj||ib> part
!
!               A_ia,jb = <aj||ib>
! ****************************************************************
! ****************************************************************
        counter_j=1
        do b=1,nvirt
          do j=1,nocc
            counter_i=1
            do a=1,nvirt
              do i=1,nocc
                CISmat(counter_j,counter_i)=CISmat(counter_j,counter_i)
     &                  +W(i,a,j,b)

                counter_i=counter_i+1
              enddo
            enddo
            counter_j=counter_j+1
          enddo
        enddo
! ****************************************************************
! ****************************************************************
!       Step2) Insert Diagonal energy elements of CIS matrix *
!             A_ia,jb = delij*delab*(e_a-e_i) 
! ****************************************************************
! ****************************************************************
        counter_j=1
        offset=nocc*nvirt
        do i=1,nocc
          do a=1,nvirt
            counter_i=1
            do j=1,nocc
              do b=1,nvirt
                 del_ij = (i == j)
                 del_ab = (a == b)
           CISmat(counter_j,counter_i)=CISmat(counter_j,counter_i)!+
c     &          del_ij*del_ab*(diaFockA(a+nocc)-diaFockA(i))

           CISmat(offset+counter_j,offset+counter_i)=
     &                  CISmat(counter_j,counter_i)

                counter_i=counter_i+1
              enddo
            end do
            counter_j=counter_j+1
          enddo
        enddo 

! ****************************************************************
! ****************************************************************
!               * Verify that CIS matrix is symmetric *
! ****************************************************************
! ****************************************************************
      print*, "symmetric CIS matrix"
      do i=1,2*nocc*nvirt
        do j=1,2*nocc*nvirt
          if (abs(CISmat(i,j)-CISmat(j,i)).gt.1e-5)then
            print*, 'not symmetric at index: ', i,j
          endif
         enddo
      enddo
! ****************************************************************
! ****************************************************************
! ****************************************************************
        end subroutine 
