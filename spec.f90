program main
implicit none
double precision pi
double complex one
integer i, j, n                                                     
parameter (n=50)
double precision x(n), z(n), TTT
double complex Psi_T(n),       Psi_P(n),       Tem(n)               
double complex hat_Psi_T(n+2), hat_Psi_P(n+4), hat_Tem(n+2)
double precision real_Psi_T, imag_Psi_T, real_Psi_P, imag_Psi_P, real_Tem, imag_Tem

pi=acos(-1.d0)
one=(0.d0, 1.d0)

do i=1, n
 read(2,'(7E15.6)') z(i), real_Psi_T, imag_Psi_T, real_Psi_P, imag_Psi_P, real_Tem, imag_Tem
 Psi_T(i)=real_Psi_T+one*imag_Psi_T
 Psi_P(i)=real_Psi_P+one*imag_Psi_P
 Tem(i)=real_Tem+one*imag_Tem
enddo
call phys_to_spec(n,Psi_P,n+4,hat_Psi_P,x)
call phys_to_spec(n,Psi_T,n+2,hat_Psi_T,x)
call phys_to_spec(n,Tem,n+2,hat_Tem,x)

open(1,file='tor',form='formatted')
do j=1, n+2
 write(1,'(I5,E15.6)') j, real(hat_Psi_T(j))**2+imag(hat_Psi_T(j))**2
enddo
close(1)

open(1,file='pol',form='formatted')
do j=1, n+4
 write(1,'(I5,E15.6)') j, real(hat_Psi_P(j))**2+imag(hat_Psi_P(j))**2
enddo
close(1)

open(1,file='tem',form='formatted')
do j=1, n+2
 write(1,'(I5,E15.6)') j, real(hat_Tem(j))**2+imag(hat_Tem(j))**2
enddo
close(1)
end program main

subroutine phys_to_spec(np,phys,ns,spec,x)
implicit none
integer  i, j, np, ns       
double complex phys(np), spec(ns)   
double precision x(np), TTT
do j=2, ns
 spec(j)=(0.d0, 0.d0)
 do i=1, np
  spec(j)=spec(j)+phys(i)*TTT(0,j-1,x(i))
 enddo
 spec(j)=spec(j)*2.d0/np
enddo
spec(1)=(0.d0, 0.d0)
do i=1, np
 spec(1)=spec(1)+phys(i)
enddo
spec(1)=spec(1)/np
end subroutine phys_to_spec

!!!   TTT(K,M,X) = the K-th derivative of Tm(X), the M-th Chebyshev polynomial evaluated at X.
      DOUBLE PRECISION FUNCTION TTT(K,M,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(0:10000), B(0:10000)
      DO 10 I=0,M-1
       A(I)=0.D0
10    CONTINUE
      A(M)=1.D0
      DO 20 J=1,K
       CALL DIFF(A,M)
20    CONTINUE
      B(M+2)=0.D0
      B(M+1)=0.D0
      DO 30 I=M,0,-1
       B(I)=2.D0 * X * B(I+1)  -  B(I+2)  +  A(I)
30    CONTINUE
      TTT=(A(0) + B(0) - B(2))/2.D0
      RETURN
      END
      SUBROUTINE DIFF(A,M)
      DOUBLE PRECISION A(0:10000), C(0:10000)
      C(M+1)=0.D0
      C(M)=0.D0
      DO 10 I=M-1,0,-1
       C(I)=C(I+2)  +  DFLOAT(2*I+2) * A(I+1)
10    CONTINUE
      C(0)=C(0)/2.D0
      DO 20 I=0,M
       A(I)=C(I)
20    CONTINUE
      RETURN
      END