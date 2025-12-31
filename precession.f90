! gfortran -o precession precession.f90 nag.f
! z=(1/2)x, D=d/dz=2d/dx, D^2=4d/dx
! Psi_T(x(i))=sum_j hat_Psi_T*TTT(0,j-1,x(i))
! i=1~n: inner points; i=n+1, n+2, n+3 n+4: boundary conditions

program main
implicit none
double precision pi
double complex one
integer i, j, n							    ! i is physical, j is spectral, n is dimension
parameter (n=3)
double precision x(n), z(n), TTT				    ! inner points
double precision Ek, Pr, epsilon, R_c, delta_R, k_x, k_y, k_z, k2_perp, k2, time, dt
double complex Psi_T(n),       Psi_P(n),       Tem(n)    	    ! spectral coefficients about (z,t)
double complex hat_Psi_T(n+2), hat_Psi_P(n+4), hat_Tem(n+2) 	    ! Chebyshev coefficients about t
double precision a1(n+2,n+2), a2(n+4,n+4), a3(n+2,n+2)  	    ! coefficient matrices 
double precision a1_inv(n+2,n+2), a2_inv(n+4,n+4), a3_inv(n+2,n+2)  ! inverse of coefficient matrices
double complex b1(n+2), b2(n+4), b3(n+2)  			    ! right-hand-side
double precision energy1, energy2, energy3			    ! spectral energy sum_j(|hat_Psi_T|^2), ...
double precision ini_re, ini_im
integer it, nt
parameter (nt=100)

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

! inner points
do i=1,n
 x(i)=-cos((2*dfloat(i)-1)/(2*dfloat(n))*pi)
 z(i)=x(i)/2.d0
enddo

! collocate Psi_T equation on inner points
do j=1, n+2
 do i=1, n
  a1(i,j)=TTT(0,j-1,x(i))
 enddo
enddo
! collocate Psi_T equation on boundary points
do j=1, n+2
 a1(n+1,j)=TTT(1,j-1,-1.d0)
 a1(n+2,j)=TTT(1,j-1,1.d0)
enddo
b1(n+1)=(0.d0, 0.d0)
b1(n+2)=(0.d0, 0.d0)
! collocate Psi_P equation on inner points
do j=1, n+4
 do i=1, n
  a2(i,j)=TTT(0,j-1,x(i))
 enddo
enddo
! collocate Psi_P equation on boundary points
do j=1, n+4
 a2(n+1,j)=TTT(0,j-1,-1.d0)
 a2(n+2,j)=TTT(0,j-1,1.d0)
 a2(n+3,j)=TTT(2,j-1,-1.d0)
 a2(n+4,j)=TTT(2,j-1,1.d0)
enddo
b2(n+1)=(0.d0, 0.d0)
b2(n+2)=(0.d0, 0.d0)
b2(n+3)=(0.d0, 0.d0)
b2(n+4)=(0.d0, 0.d0)
! collocate Tem equation on inner points
do j=1, n+2
 do i=1, n
  a3(i,j)=TTT(0,j-1,x(i))
 enddo
enddo
! collocate Tem equation on boundary points
do j=1, n+2
 a3(n+1,j)=TTT(0,j-1,-1.d0)
 a3(n+2,j)=TTT(0,j-1,1.d0)
enddo
b3(n+1)=(0.d0, 0.d0)
b3(n+2)=(0.d0, 0.d0)
!!! test
do i=1, n+2
 write(6,'(5E15.6)') a1(i,:)
enddo
!!!
! calculate inverse of coefficient matrices a1, a2, a3 for time stepping
call mat_inv(n+2,a1,a1_inv)
call mat_inv(n+4,a2,a2_inv)
call mat_inv(n+2,a3,a3_inv)
!!! test
energy1=0.d0
energy2=0.d0
do j=1, n+2
 do i=1, n+2
  energy1=energy1+abs(a1(i,j))**2
  energy2=energy2+abs(a1_inv(i,j))**2
 enddo
enddo
call energy(n+2,b1,energy3)
write(6,*) energy1, energy2, energy3
!!!

! give initial condition of Chebyshev coefficients
call random_seed()
do j=1, n+2
 call random_number(ini_re)
 call random_number(ini_im)
 hat_Psi_T(j)=ini_re+one*ini_im
enddo
do j=1, n+4
 call random_number(ini_re)
 call random_number(ini_im)
 hat_Psi_P(j)=ini_re+one*ini_im
enddo
do j=1, n+2
 call random_number(ini_re)
 call random_number(ini_im)
 hat_Tem(j)=ini_re+one*ini_im
enddo

! time stepping
time=0.d0
dt=1.d-6
do it=1, nt
 time=time+it*dt
 call spec_to_phys(n+2,hat_Psi_T,n,Psi_T,x)
 call spec_to_phys(n+4,hat_Psi_P,n,Psi_P,x)
 call spec_to_phys(n+2,hat_Tem,n,Tem,x)
!!! test
if(it.eq.1) then 
 call energy(n,Psi_T,energy1)
 write(6,*) energy1
endif
!!!
 do i=1, n
  b1(i)=(0.d0, 0.d0)
  b2(i)=(0.d0, 0.d0)
  b3(i)=(0.d0, 0.d0)
  do j=1, n
   b1(i)=b1(i)+TTT(2,j-1,x(i))*hat_Psi_T(j)+TTT(1,j-1,x(i))*hat_Psi_P(j)
   b2(i)=b2(i)+TTT(2,j-1,x(i))**2*hat_Psi_P(j)-TTT(2,j-1,x(i))*hat_Psi_T(j)
   b3(i)=b3(i)+TTT(2,j-1,x(i))*hat_Tem(j)+TTT(1,j-1,x(i))*hat_Psi_P(j)
  enddo
  b1(i)=Psi_T(i)+dt*b1(i)
  b2(i)=Psi_P(i)+dt*b2(i)
  b2(i)=Tem(i)+dt*b3(i)
 enddo
!!! test
if(it.eq.1) then 
 call energy(n+2,b1,energy1)
 write(6,*) energy1
endif
!!!
 call mat_mul(n+2,a1_inv,b1,hat_Psi_T)
 call mat_mul(n+4,a2_inv,b2,hat_Psi_P)
 call mat_mul(n+2,a3_inv,b3,hat_Tem)
 call energy(n+2,hat_Psi_T,energy1)
 call energy(n+4,hat_Psi_P,energy2)
 call energy(n+2,hat_Tem,energy3)
 write(1,'(4E15.6)') time, energy1, energy2, energy3
enddo

! output Psi_T, Psi_P, Tem
do i=1,n
 write(2,'(4E15.6)') z(i), abs(Psi_T(i)), abs(Psi_P(i)), abs(Tem(i))
enddo
end program main

subroutine spec_to_phys(ns,spec,np,phys,x)
implicit none
integer i, j, ns, np
double complex spec(ns), phys(np)
double precision x(np), TTT
do i=1,np
 phys(i)=0.d0
 do j=1, ns
  phys(i)=phys(i)+spec(j)*TTT(0,j,x(i))
 enddo
enddo
end subroutine spec_to_phys

subroutine mat_inv(n,a,c)
implicit none
integer i, j, n
double precision a(n,n), b(n,n), c(n,n)
double precision WKSPCE(n), AA(n,n), BB(n,n)
integer IFAIL
do i=1, n
 do j=1, n
  b(i,j)=0.d0
  if(i.eq.j) b(i,j)=1.d0
 enddo
enddo
call F04AEF(a,n,b,n,n,n,c,n,WKSPCE,AA,n,BB,n,IFAIL)
end subroutine mat_inv

subroutine mat_mul(n,a,b,c)
implicit none
integer i, k, n
double precision a(n,n)
double complex b(n), c(n)
do i=1,n
 c(i)=0.d0
 do k=1,n
  c(i)=c(i)+a(i,k)*b(k)
 enddo
enddo
end subroutine mat_mul

subroutine energy(n,a,x)
implicit none
integer i, n
double complex a(n)
double precision x
x=0.d0
do i=1, n
 x=x+abs(a(i))**2
enddo
end subroutine energy

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