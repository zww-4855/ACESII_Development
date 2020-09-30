      SUBROUTINE totalEnergy(T1aa,T1bb,fia_aa,fia_bb,T2aa,T2ab,T2bb,
     &                   Waa,Wab,Wbb,nocc,nvirt,iuhf,Eccsd) 
        integer, intent(in)::nocc,nvirt,iuhf
        double precision::factor
        double precision, intent(in)::T1aa(nvirt,nocc),T1bb(nvirt,nocc),
     &          fia_aa(nvirt,nocc),fia_bb(nvirt,nocc)

        double precision, intent(in)::T2aa(nvirt,nvirt,nocc,nocc),
     &          T2bb(nvirt,nvirt,nocc,nocc),T2ab(nvirt,nvirt,nocc,nocc)

        double precision, intent(in)::Waa(nvirt,nvirt,nocc,nocc),
     &          Wbb(nvirt,nvirt,nocc,nocc),Wab(nvirt,nvirt,nocc,nocc)

        double precision::TauAA(nvirt,nvirt,nocc,nocc),
     &                  TauAB(nvirt,nvirt,nocc,nocc)
        double precision,intent(inout)::Eccsd(nocc)
        double precision::EccsdAA(nocc), EccsdAB(nocc)
        double precision:: tot,energyAA,energyAb
        integer :: i,j,a,b
        double precision:: Vabtest(nocc*nvirt*nocc*nvirt)
        if (iuhf.eq.1) then
          factor=1.0d0
        else
          factor=2.0d0
        endif
! Build Tau and contract with W

        EccsdAB=0.0d0
        EccsdAA=0.0d0
        Eccsd=0.0d0
        do i=1,nocc
         do j=1,nocc
          do a=1,nvirt
            do b=1,nvirt
              TauAB(a,b,i,j)=T2ab(a,b,i,j)+T1aa(a,i)*T1aa(b,j)
              TauAA(a,b,i,j)=T2aa(a,b,i,j)+T1aa(a,i)*T1aa(b,j)
     &                      - T1aa(b,i)*T1aa(a,j)
              EccsdAB(i)=EccsdAB(i)+TauAB(a,b,i,j)*Wab(a,b,i,j)
              EccsdAA(i)=EccsdAA(i)+TauAA(a,b,i,j)*Waa(a,b,i,j)


              enddo
            enddo
          enddo
        enddo
        energyAA=0.0d0
        energyAb=0.0d0
        do i=1,nocc
c          print*, 'Eccsdab', Eccsdab(i)
c          print*, 'Eccsdaa', Eccsdaa(i)
c          print*
          Eccsd(i)=EccsdAB(i)+2*EccsdAA(i)
          energyAA=energyAA+2*EccsdAA(i)
          do a=1,nvirt
            energyAb=energyAb+2*T1aa(a,i)*fia_aa(a,i)
            Eccsd(i)=Eccsd(i)+2*T1aa(a,i)*fia_aa(a,i)
          enddo
        enddo
C        do i=1,nocc
c          do j=1,nvirt
c                print*, T1aa(i,j)
c          enddo
c        enddo
        tot=0.0d0
        print*,'*****************************************************'
        print*,'*****************************************************'
        print*,'** Correlation Energy decomposed according to NLMO**'
        print*
        do i=1,nocc
          print*, 'orbital #', i,' has e_corr = ', Eccsd(i)
          tot=tot+Eccsd(i)
        enddo
        print*,'*****************************************************'
        print*,'*****************************************************'
        print*
        print*,'** Total CCSD energy from all occupied NLMOs **'
        print*, 'AA correlation contribution is: ', sum(EccsdAA)
        print*, 'AB correlation contribution is: ', sum(EccsdAB)
        print*,'total CCSD correlation energy: ', tot
        print*,'*****************************************************'
      END SUBROUTINE 
