










      subroutine Get_Cmatrix(C,X,PrimeC,nbas)
        implicit none
        real(kind=8)::PrimeC(nbas,nbas)
        real(kind=8)::C(nbas,nbas),X(nbas,nbas)
        integer ::nbas
        C=0.0d0
        C=matmul(X,PrimeC)
      end
