










        subroutine GetNLMOQM2(NLMOQM1,NLMOQM2,nbas,QM2NLMOcount)
        integer, intent(in)::nbas
        integer, intent(inout):: QM2NLMOcount
        integer, intent(in):: NLMOQM1(nbas),NLMOQM2(nbas)
        integer::i,counter
        counter=1
        QM2NLMOcount=0
        print*,'inside getnlmoqm2'
        do i=1,size(NLMOQM2)
          if (NLMOQM2(i).eq.0) cycle
          if (not(any(NLMOQM2(i).eq.NLMOQM1))) then
          QM2NLMOcount=QM2NLMOcount+1
          endif
        enddo




        end subroutine
