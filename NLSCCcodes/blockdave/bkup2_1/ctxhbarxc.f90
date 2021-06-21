      subroutine ctxhbarxc(Hbar,C,thetaNew,eVec,NSize,i,initVec)
      integer,intent(in)::NSize,i,initVec
      double precision, intent(in):: Hbar(Nsize,Nsize),C(Nsize,i)
      double precision, intent(inout)::thetaNew(initVec)
      double precision, intent(inout)::eVec(i,i)
      double precision:: T(i,i),Scr(Nsize,i)
      integer::a,junk

! Purpose: Computes the matrix-vector multiplication of C^t*Hbar*C
!          which is stored in T. T is then diagonalized to form evecs and
!          evalues.
! In general, initVec < i unless we are in first iteration of Davidson loop

      T=0.0d0

! vector-Matrix-vector multiplication

      call dgemm('N','N',NSize,i,Nsize,1.0d0,Hbar,Nsize,C,Nsize,0.0d0,Scr,Nsize)

!call dgemm('T','N',i,i,Nsize,1.0d0,C,i,Scr,i,0.0d0,T,i)

      T=matmul(transpose(C),Scr)

! Diagonalize tridiagonal matrix T, then copy evalues to thetaNew
! This routine works because T is symmetric; otherwise use blas then sort
! the resulting evecs/values

      call eig(T,eVec,junk,i,0)


      do a=1,initVec
         thetaNew(a)=T(a,a)
      enddo
      print*,'thetaNew in ctxhbarxc'
      print*,thetaNew

!  print*,'evals:'
!  print*,T
!  print*,'evecs'
!  print*,eVec
!  print*,'matrix mult evec^t T evec'
!  print*,matmul(transpose(eVec),matmul(T,eVec))

   Return
   end 
