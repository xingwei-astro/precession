! gfortran -o precession precession.f90 nag.f
! z=(1/2)x, D=d/dz=2d/dx, D^2=4d/dx

module globe
implicit none
double precision pi
double complex one
integer n
parameter (n=10)
double precision x(0:n+1)
end module globe

program main
use globe
implicit none
integer i, j
double precision z(0:n+1)
double precision ini_re, ini_im
double complex Psi_T(0:n+1),     Psi_P(0:n+1),     Tem(0:n+1)     ! coefficients about (z,t)
double complex hat_Psi_T(0:n+1), hat_Psi_P(0:n+1), hat_Tem(0:n+1) ! Chebyshev coefficients about t
double complex a1(0:n+1,0:n+1), a2(0:n+1,0:n+1), a3(0:n+1,0:n+1)
double complex a1_inv(0:n+1,0:n+1), a2_inv(0:n+1,0:n+1), a3_inv(0:n+1,0:n+1)
double complex b1(0:n+1), b2(0:n+1), b3(0:n+1)
double precision Ek, Pr, epsilon, R_c, delta_R, k_x, k_y, k_z, k2_perp, k2, time, dt
double precision TTT
integer it, nt
parameter (nt=10000)

pi=acos(-1.d0)
one=(0.d0, 1.d0)
Ek=1.d-4
Pr=1.d0
epsilon=1.d-6
k_x=pi
k_y=pi
k_z=pi
k2_perp=k_x**2+k_y**2
k2=k2_perp+k_z**2
R_c=4.d0*k_z**2/k2_perp*2.d0*Pr/(1.d0+Pr)+2.d0*Ek**2*(1.d0+1.d0/Pr)*k2**3/k2_perp
delta_R=1.d-3

x(0)=-1.d0
x(n+1)=1.d0
do i=1,n
 x(i)=-cos(pi/dfloat(n)*(i-0.5d0))
enddo
do i=0,n+1
 z(i)=x(i)/2.d0
enddo

! give initial condition
call random_seed()
do i=1,n
 call random_number(ini_re)
 call random_number(ini_im)
 Psi_T(i)=ini_re+one*ini_im
 call random_number(ini_re)
 call random_number(ini_im)
 Psi_P(i)=ini_re+one*ini_im
 call random_number(ini_re)
 call random_number(ini_im)
 Tem(i)=ini_re+one*ini_im
enddo
Psi_T(0)=Psi_T(1)
Psi_T(n+1)=Psi_T(n)
Psi_P(0)=(0.d0, 0.d0)
Psi_P(n+1)=(0.d0, 0.d0)
Tem(0)=(0.d0, 0.d0)
Tem(n+1)=(0.d0, 0.d0)
call phys_to_spec(Psi_T,hat_Psi_T)
call phys_to_spec(Psi_P,hat_Psi_T)
call phys_to_spec(Tem,hat_Tem)

stop

! collocate Psi_T equation on inner points
do j=0, n+1
 do i=1, n
  a1(i,j)=TTT(0,j,x(i))
 enddo
enddo
! collocate Psi_T equation on boundary points
do j=0, n+1
 a1(0,j)=TTT(1,j,x(0))
 a1(n+1,j)=TTT(1,j,x(n+1))
enddo
!call mat_inv(n+2,a1,a1_inv)
!call mat_inv(n+2,a2,a2_inv)
!call mat_inv(n+2,a3,a3_inv)

! time stepping
time=0.d0
dt=1.d-1
do it=1, nt
 time=time+it*dt
 call spec_to_phys(hat_Psi_T,Psi_T)
 call spec_to_phys(hat_Psi_P,Psi_P)
 call spec_to_phys(hat_Tem,Tem)
 do i=1, n
  b1(i)=(0.d0, 0.d0)
  do j=0, n+1
   b1(i)=b1(i)+(dt*TTT(2,j,x(i))+TTT(0,j,x(i)))*Psi_T(j)
  enddo   
 enddo
 b1(0)=(0.d0, 0.d0)
 b1(n+1)=(0.d0, 0.d0)
!call mat_mul(n+2,a1_inv,b1,hat_Psi_T)
!call mat_mul(n+2,a2_inv,b2,hat_Psi_P)
!call mat_mul(n+2,a3_inv,b3,hat_Tem)
 write(1,'(4E15.6)') time, hat_Psi_T(1), hat_Psi_P(1), hat_Tem(1)
enddo

! output module of Psi_T, Psi_P, Tem
call spec_to_phys(hat_Psi_T,Psi_T)
call spec_to_phys(hat_Psi_P,Psi_P)
call spec_to_phys(hat_Tem,Tem)
end program main

subroutine spec_to_phys(spec,phys)
use globe
implicit none
integer i, j
double complex phys(0:n+1), spec(0:n+1)
double precision TTT
do i=0,n+1
 phys(i)=0.d0
 do j=0, n+1
  phys(i)=phys(i)+spec(j)*TTT(0,j,x(i))
 enddo
enddo
end subroutine spec_to_phys

subroutine phys_to_spec(phys,spec)
use globe
implicit none
integer i, j
double complex phys(0:n+1), spec(0:n+1)
double precision TTT, a(0:n+1,0:n+1)
do j=0, n+1
 do i=0, n+1
  a(i,j)=TTT(0,j,x(i))
 enddo
enddo
!call r8mat_fs(n+2,a,phys,spec)
end subroutine phys_to_spec

subroutine mat_inv(n,a,b)
implicit none
integer n
double complex a(n,n), b(n,n)
double precision RINT(n), DETR, DETI
integer IDETE, IFAIL
call F03AHF(n,A,n,DETR,DETI,IDETE,RINT,IFAIL)
call F04AKF(n,n,A,n,RINT,b,n)
end subroutine mat_inv

subroutine mat_mul(n,a,b,c)
implicit none
integer i, j, k, n
double complex a(n,n), b(n,n), c(n,n)
do i=1,n
 do j=1,n
  c(i,j)=0.d0
  do k=1,n
   c(i,j)=c(i,j)+a(i,k)*b(k,j)
  enddo
 enddo
enddo
end subroutine mat_mul

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!