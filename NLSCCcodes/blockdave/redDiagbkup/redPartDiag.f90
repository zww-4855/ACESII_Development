subroutine redPartDiag(CISmat,NSize,TOL)
  integer,intent(in)::Nsize
  double precision,intent(in)::TOL
  double precision,intent(in)::CISmat(NSize,NSize)
  double precision::P(NSize,1),Q(Nsize,NSize-1)


  double precision::projP(NSize,NSize),projQ(NSize,NSize)
  double precision::Pert(NSize,NSize),Denom(NSize,NSize)
integer::nbas,natm,nocc,nvirt,scratch(2)
  integer::order

  double precision::H0(NSize,NSize),tmpHold(NSize-1,NSize-1)
  double precision::V(NSize,NSize),D(NSize-1,NSize-1)
  double precision::V0out(NSize,NSize),pVq(1,NSize-1)
  double precision::qVq(NSize-1,NSize-1)
  double precision::E_0(1,1),E_1(1,1)
  double precision,allocatable::AA(:,:),BB(:,:),CC(:,:)
  double precision::tmp(1,1),initEval(1,1),evaltmp(1,1)
  double precision::guessE(1,1),correction(1,1)
  logical::BWPT,RSPT,Converged
  
!      call getrec(1,'JOBARC','NBASTOT',1,nbas)
!      call getrec(1,'JOBARC','NREALATM',1,natm)
!        call getrec(1,'JOBARC','NOCCORB',2,scratch)
!      nocc=scratch(1)
!      nvirt=nbas-nocc
!print*,'nocc:',nocc
!print*,'nvirt:',nvirt
!print*,'NSize:',NSize,nocc*nvirt
!print*,'Nbas:',nbas

  call checksum("initVec:",cisEvecs,NSize*Nblocks,s)
! ** Pick lowest element of diagonal of CIS mat for defining P and Q spaces
  call definePQ(CISmat,NSize,P,Q,projP,projQ)

! ** Use P and Q defininitions to define::
!       PHP and PHQ R0 QHP = PH (\sumQ |Q><Q|/[Eo-EQ]) HP
!                          Eo-EQ=(ei-ea)+(ej-eb)

! 0th order approximation: PHP=(ei-ea)+<ai||jb>
  call defineH0(CISmat,NSize,P,initEval,H0)
  print*,'b4 Vdef'
  call defineV(CISmat,NSize,V)
  print*,'after Vdef'



  BWPT=.True.
  RSPT=.False.


if (RSPT) then
  call defineD(D,initEval,NSize,H0,Q,guessE,RSPT,BWPT)

  xVq=0.0d0
  qVDx=0.0d0
  qVDq=0.0d0
! RS PT correction:

  tmpHold=0.0d0
  evaltmp=initEval(1,1)
  do order=0,1
    allocate(AA(order,order),BB(order,order),CC(order,1))

    call buildRSPert(V,D,P,Q,NSize,pVq,qVq,tmpHold,order,E_0,evaltmp)
    print*
    print*,'****************************'
    print*,'Order: ', order+1
    print*,'Correction at this order:', E_0
    evaltmp=evaltmp+E_0(1,1)
    print*,'Current excitation energy:',evaltmp
    print*,'****************************'
!    if (order.eq.1) then
!      call buildRSPert(V,D,P,Q,NSize,pVq,qVq,tmpHold,0,E_0)
!      print*,'first order correction (reg):',E_0
!      print*,'first order estimate of lowest root:',E_0+initEval(1,1)
!      evaltmp=initEval(1,1)+E_0(1,1)
!      call buildRSPert(V,D,P,Q,NSize,pVq,qVq,tmpHold,1,E_1)
!      print*,'second order correction(reg):',E_1
!      print*,'second order estimate of lowest root:',evaltmp+E_1
!      tmp=(E_0)**2/(E_0-E_1)
!      print*,'first order correction (Pade):',evaltmp+tmp
!    endif
   



    deallocate(AA,BB,CC)
  enddo

endif

if (BWPT) then
  guessE=initEval(1,1)
  print*,'Initial guess:', guessE
  pVq=0.0d0
  qVq=0.0d0

  qVDx=0.0d0
  qVDq=0.0d0
! RS PT correction:

  tmpHold=0.0d0
  evaltmp=initEval(1,1)
  do order=0,0
    allocate(AA(order,order),BB(order,order),CC(order,1))

    call defineD(D,initEval,NSize,H0,Q,guessE,RSPT,BWPT)
    call buildBWPert(V,D,P,Q,NSize,pVq,qVq,tmpHold,order,E_0)
    print*
    print*,'****************************'
    print*,'Order: ', order+1
    print*,'Correction at this order:', E_0
    evaltmp=evaltmp+E_0(1,1)
    print*,'Current excitation energy:',evaltmp
    print*,'****************************'
    guessE=evaltmp
!    print*,'Writing perturbation to file....'
!    V0out=matmul(V,matmul(D,V))
!    open(1000,file='V1.out')
!    do i=1,NSize
!      do j=1,NSize
!        write(1000,*) V0out(j,i)
!      enddo
!    enddo    
    



    deallocate(AA,BB,CC)
  enddo

!! REDEFINED H0/V AS A CHECK FOR H20
  print*,'Redefining H0/V .....'
  print*
  print*
  call altDefine_H0_V(CISmat,NSize,projP,P,initEval,H0,V)
  guessE=initEval(1,1)

  pVq=0.0d0
  qVq=0.0d0

  qVDx=0.0d0
  qVDq=0.0d0
! RS PT correction:

  tmpHold=0.0d0
  evaltmp=initEval(1,1)
  do order=0,4
    call defineD(D,initEval,NSize,H0,Q,guessE,RSPT,BWPT)
    call buildBWPert(V,D,P,Q,NSize,pVq,qVq,tmpHold,order,E_0)
    print*
    print*,'****************************'
    print*,'Order: ', order+1
    print*,'Correction at this order:', E_0
    evaltmp=evaltmp+E_0(1,1)
    print*,'Current excitation energy:',evaltmp
    print*,'****************************'
    guessE=evaltmp
  enddo

endif
!!
!!
!!
!!
!! ** ITERATIVE EXAMPLE AT EACH ORDER
!  guessE=initEval(1,1)
!
!  pVq=0.0d0
!  qVq=0.0d0
!
!  qVDx=0.0d0
!  qVDq=0.0d0
!! RS PT correction:
!
!  tmpHold=0.0d0
!  evaltmp=initEval(1,1)
!  Converged=.False.
!  do order=0,4
!    correction=0.0d0
!    Converged=.False.
!    call defineD(D,initEval,NSize,H0,Q,guessE,RSPT,BWPT)
!    call buildBWPert(V,D,P,Q,NSize,pVq,qVq,tmpHold,order,E_0)
!    print*
!    print*,'****************************'
!    print*,'Order: ', order+1
!    print*,'Correction at this order:', E_0
!    evaltmp=evaltmp+E_0(1,1)
!    print*,'Current excitation energy:',evaltmp
!    print*,'****************************'
!    guessE=evaltmp
!
!    do while (.not.Converged)
!      call defineD(D,initEval,NSize,H0,Q,guessE,RSPT,BWPT)
!      call buildBWPert(V,D,P,Q,NSize,pVq,qVq,tmpHold,order,E_0)
!      if (abs(abs(correction)-abs(E_0)).gt.0.00000001) then
!        guessE=
!      else
!        Converged=.True.
!      endif
!    enddo
!
!
!
!  enddo
!
!
!endif


end subroutine  
