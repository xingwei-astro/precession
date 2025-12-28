! gfortran -o precession precession.f90 solve.f90
! z=(1/2)x, D=d/dz=2d/dx, D^2=4d/dx

program main
implicit none
double precision pi
double complex one
integer i, j, n
parameter (n=10)
double precision x(0:n+1), z(0:n+1)
double precision ini_re, ini_im
double complex Psi_T(0:n+1),     Psi_P(0:n+1),     Tem(0:n+1)     ! coefficients about (z,t)
double complex hat_Psi_T(0:n+1), hat_Psi_P(0:n+1), hat_Tem(0:n+1) ! Chebyshev coefficients about t
double complex a(0:n+1,0:n+1), b(0:n+1)
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

! calculate Chebyshev coefficients of initial condition
do j=0, n+1
 do i=0, n+1
  a(i,j)=TTT(0,j,x(i))
 enddo
enddo
do i=1, n
 b(i)=Psi_T(i)   
enddo
b(0)=Psi_T(0)
b(n+1)=Psi_T(n+1)
!call r8mat_fs(n+2,a,b,hat_Psi_T)

stop

! time stepping
dt=1.d-1
do it=1, nt
! collocate equation on inner points
 do j=0, n+1
  do i=1, n
   a(i,j)=TTT(0,j,x(i))
  enddo
 enddo
! collocate equation on boundary points
 do j=0, n+1
  a(0,j)=TTT(0,j,x(0))
  a(n+1,j)=TTT(0,j,x(n+1))
 enddo
! the right-hand-side vector
 do i=1, n
  b(i)=0.d0
  do j=0, n+1
   b(i)=b(i)+(dt*TTT(2,j,x(i))+TTT(0,j,x(i)))*Psi_T(j)
  enddo   
 enddo
 b(0)=Psi_T(0)
 b(n+1)=Psi_T(n+1)
! solve a*hat_Psi_T=b
! call r8mat_fs(n+2,a,b,hat_Psi_T)
 write(1,'(4E15.6)') time, hat_Psi_T(1), hat_Psi_P(1), hat_Tem(1)
enddo

! output module of Psi_T, Psi_P, Tem
do i=0,n+1
 Psi_T(i)=0.d0
 do j=0, n+1
  Psi_T(i)=Psi_T(i)+hat_Psi_T(j)*TTT(0,j,x(i))
 enddo
 write(2,'(4E15.6)') x(i), Psi_T(i), Psi_P(i), Tem(i) 
enddo
end program main

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