      subroutine Get_FockPrime(FockPrime,X,FockMat,nbas)
        implicit none
        integer:: nbas
        real(kind=8)::FockPrime(nbas,nbas),X(nbas,nbas)
        real(kind=8)::FockMat(nbas,nbas)
        FockPrime=0.0d0
        FockPrime=matmul(matmul(transpose(X),FockMat),X)
      end
