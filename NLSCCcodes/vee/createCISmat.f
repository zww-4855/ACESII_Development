










        subroutine createCISmat(W,Wab,nocc,nvirt,MATDIM)
          integer, intent(in):: nocc,nvirt,MATDIM
          double precision, intent(in)::W(nocc,nvirt,nocc,nvirt)
          double precision, intent(in)::Wab(nvirt*nocc*nocc*nvirt)


          double precision,allocatable:: intermedAC(:),acesCISevecs(:)
          double precision::CISmat(2*nocc*nvirt,2*nocc*nvirt)
        double precision:: output1,temp
        integer ::iter,del_ij,del_ab,nbas,offset
        integer :: counter_i, counter_j,i,j,a,b
        real(kind=8) :: term,fab, fij,fe,HFenergy
        COMMON/CALCINFO/NROOT(8)
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
      print*,'CIS matrix is symmetric'
! ****************************************************************
! ****************************************************************
! ****************************************************************
!#ifdef _DEBUG_LVL0
        print*,'*****************************************************'
        print*,'******           FullCIS matrix                ******'
        print*,'*****************************************************'
        call output(CISmat,1,MATDIM,1,MATDIM,MATDIM,MATDIM,1)
        print*,'reprint cis mat'
        print*,CISmat
        print*,'*****************************************************'
        print*,'NUMBER OF ROOTS: ', NROOT(1)
        allocate(intermedAC(NROOT(1)*2*nocc*nvirt),
     &                  acesCISevecs(NROOT(1)*2*nocc*nvirt))
        acesCISevecs=0.0d0
        call Getlst(acesCISevecs,1,NROOT(1),1,1,94)
        print*,'acesii CIS eigenvecs'
        print*, acesCISevecs

        print*,'*****************************************************'
        print*,'******           FullCIS eigenvectors          ******'
        print*,'*****************************************************'
        call output(acesCISevecs,1,MATDIM,1,6,6,6,1)

!        print*,'MATDIM',MATDIM,nocc,nvirt
!        intermedAC=0.0d0
!        iter=1
!        do i=1,1!NROOT(1)!4!6
!           call xgemm('N','N',2*nocc*nvirt,1,2*nocc*nvirt,1.0D0,CISmat,
!     &          2*nocc*nvirt,acesCISevecs(iter),2*nocc*nvirt,0.0D0,
!     &              intermedAC(iter),2*nocc*nvirt)
!          print*,'intermedAC',intermedAC(iter:2*nocc*nvirt)
!          print*,'aces evecs', acesCISevecs(iter:2*nocc*nvirt)
!!           temp=ddot(2*nocc*nvirt,acesCISevecs(iter),1,intermedAC(iter)
!!     &                                  ,1)
!           temp=sdot(6,acesCISevecs(iter),1,
!     &                  intermedAC(iter),1)
!           print*,'recalculating ACESII vars',temp*27.2114
!           print*, 'test', ddot(6,acesCISevecs(iter:2*nocc*nvirt),1,
!     &                  intermedAC(iter:2*nocc*nvirt),1)
!           print*,'dtprd',27.2114*
!     &          dot_product(acesCISevecs(1:6),intermedAC(1:6))
!          output1=0.0d0
!          do j=1,6
!                output1=output1+acesCISevecs(j)*intermedAC(j)
!          enddo
!          print*,'output',27.2114*output1
!           iter=iter+2*nocc*nvirt
!        enddo



        deallocate(intermedAC,acesCISevecs)
        end subroutine 
