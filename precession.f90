! gfortran -o precession precession.f90 nag.f
! Crank-Nicolson scheme, diffusion terms half-implicit half-explicit, the other terms explicit
! z=(1/2)x, D=d/dz=2d/dx, D^2=4d/dx, D^4=16d/dx
! Psi_T(x(i))=sum_j hat_Psi_T(j)*TTT(0,j-1,x(i)), j from 1 and j-1 from 0
! i=1~n: inner points; i=n+1, n+2, n+3 n+4: boundary conditions

program main
implicit none
double precision pi
double complex one
integer i, j, n							    ! i physical, j spectral, n dimension
parameter (n=50)
double precision x(n), z(n), TTT				    ! inner points
double precision Ek, Pr, epsilon, R_c, delta_R			    ! dimensionless parameters
double precision k_x, k_y, k_z, k2_perp, k2			    ! wavenumbers
double complex Psi_T(n),       Psi_P(n),       Tem(n)    	    ! coefficients about (z,t) -- physical space
double complex hat_Psi_T(n+2), hat_Psi_P(n+4), hat_Tem(n+2) 	    ! Chebyshev coefficients about t -- spectral space
double precision a1(n+2,n+2), a2(n+4,n+4), a3(n+2,n+2)  	    ! coefficient matrices 
double precision a1_inv(n+2,n+2), a2_inv(n+4,n+4), a3_inv(n+2,n+2)  ! inverse of coefficient matrices
integer it, nt							    ! time steps	
double precision dt, time, c1, c2				    ! c1 and c2 are coefficients of precession terms
double complex D1_Psi_T(n), D2_Psi_T(n)				    ! derivatives of Psi_T
double complex D1_Psi_P(n), D2_Psi_P(n), D4_Psi_P(n)		    ! derivatives of Psi_P
double complex D2_Tem(n)					    ! derivatives of Tem
double complex b1(n+2), b2(n+4), b3(n+2)  			    ! right-hand-side terms
double precision energy1_0, energy2_0, energy3_0		    ! spectral energy at the last timestep
double precision energy1_1, energy2_1, energy3_1		    ! spectral energy at the next timestep
double precision sigma_u, sigma_T				    ! growth rate

pi=acos(-1.d0)
one=(0.d0, 1.d0)
Ek=1.d-4
Pr=1.d0
epsilon=0.2d0
k_x=5.748d0*pi		! resonance condition of two inertial waves
k_y=0.d0
k_z=pi
k2_perp=k_x**2+k_y**2
k2=k2_perp+k_z**2
!R_c=4.d0*k_z**2/k2_perp+Ek**2*k2**3/k2_perp*(1.d0+2.d0/Pr)*(1.d0+1.d0/Pr)
R_c=0.d0
delta_R=0.d0
dt=1.d-2
nt=200000
write(6,'(7(A10,E15.6,/))') 'Ek=', Ek, 'R_c=', R_c, 'delta_R=', delta_R, 'epsilon=', epsilon, &
                            'k_z=', k_z, 'k_perp=', sqrt(k2_perp), 'dt=', dt

! inner points
do i=1,n
 x(i)=-cos((2*dfloat(i)-1)/(2*dfloat(n))*pi)
 z(i)=x(i)/2.d0
enddo

! collocate Psi_T equation on inner points
do j=1, n+2
 do i=1, n
  a1(i,j)=1.d0/dt*TTT(0,j-1,x(i))-0.5d0*Ek*(4*TTT(2,j-1,x(i))-k2_perp*TTT(0,j-1,x(i)))
 enddo
enddo
! collocate Psi_T equation on boundary points
do j=1, n+2
 a1(n+1,j)=2.d0*TTT(1,j-1,-1.d0)
 a1(n+2,j)=2.d0*TTT(1,j-1,1.d0)
enddo
b1(n+1)=(0.d0, 0.d0)
b1(n+2)=(0.d0, 0.d0)
! collocate Psi_P equation on inner points
do j=1, n+4
 do i=1, n
  a2(i,j)=1.d0/dt*(4*TTT(2,j-1,x(i))-k2_perp*TTT(0,j-1,x(i))) &
         -0.5d0*Ek*(16*TTT(4,j-1,x(i))-8*k2_perp*TTT(2,j-1,x(i))+k2_perp**2*TTT(0,j-1,x(i)))
 enddo
enddo
! collocate Psi_P equation on boundary points
do j=1, n+4
 a2(n+1,j)=TTT(0,j-1,-1.d0)
 a2(n+2,j)=TTT(0,j-1,1.d0)
 a2(n+3,j)=4.d0*TTT(2,j-1,-1.d0)
 a2(n+4,j)=4.d0*TTT(2,j-1,1.d0)
enddo
b2(n+1)=(0.d0, 0.d0)
b2(n+2)=(0.d0, 0.d0)
b2(n+3)=(0.d0, 0.d0)
b2(n+4)=(0.d0, 0.d0)
! collocate Tem equation on inner points
do j=1, n+2
 do i=1, n
  a3(i,j)=1.d0/dt*TTT(0,j-1,x(i))-0.5d0*Ek/Pr*(4*TTT(2,j-1,x(i))-k2_perp*TTT(0,j-1,x(i)))
 enddo
enddo
! collocate Tem equation on boundary points
do j=1, n+2
 a3(n+1,j)=TTT(0,j-1,-1.d0)
 a3(n+2,j)=TTT(0,j-1,1.d0)
enddo
b3(n+1)=(0.d0, 0.d0)
b3(n+2)=(0.d0, 0.d0)
! calculate inverse of coefficient matrices a1, a2, a3 for time stepping
call mat_inv(n+2,a1,a1_inv)
call mat_inv(n+4,a2,a2_inv)
call mat_inv(n+2,a3,a3_inv)

! initial condition of Chebyshev coefficients
!do j=1, n+2
! hat_Psi_T(j)=1.d-1/dfloat(j)
!enddo
!do j=1, n+4
! hat_Psi_P(j)=1.d-1/dfloat(j)
!enddo
!do j=1, n+2
! hat_Tem(j)=1.d-1/dfloat(j)
!enddo
hat_Psi_T(:) = 0.d0
hat_Tem(:)   = 0.d0
hat_Psi_P(:) = 0.d0
! k_z = pi mode
hat_Psi_P(2) =  1.d-6
hat_Psi_P(4) = -1.d-6
! k_z = 2pi mode
hat_Psi_P(3) =  1.d-6
hat_Psi_P(5) = -1.d-6

! time stepping
do it=1, nt
 time=it*dt
 c1=one*k_x*sin(time)+one*k_y*cos(time)
 c2=one*k_x*sin(time)-one*k_y*cos(time)
 ! calculate spectral energy before update
 call energy(n+2,hat_Psi_T,energy1_0)
 call energy(n+4,hat_Psi_P,energy2_0)
 call energy(n+2,hat_Tem,energy3_0)
 ! rescale by energy to avoid too large value
 !hat_Psi_T=hat_Psi_T/sqrt(energy1_0+energy2_0+energy3_0)
 !hat_Psi_P=hat_Psi_P/sqrt(energy1_0+energy2_0+energy3_0)
 !hat_Tem=hat_Tem/sqrt(energy1_0+energy2_0+energy3_0)
 ! calculate right-hand-side terms
 call spec_to_phys(n+2,hat_Psi_T,n,Psi_T,x,0)
 call spec_to_phys(n+2,hat_Psi_T,n,D1_Psi_T,x,1)
 call spec_to_phys(n+2,hat_Psi_T,n,D2_Psi_T,x,2)
 call spec_to_phys(n+4,hat_Psi_P,n,Psi_P,x,0)
 call spec_to_phys(n+4,hat_Psi_P,n,D1_Psi_P,x,1)
 call spec_to_phys(n+4,hat_Psi_P,n,D2_Psi_P,x,2)
 call spec_to_phys(n+4,hat_Psi_P,n,D4_Psi_P,x,4)
 call spec_to_phys(n+2,hat_Tem,n,Tem,x,0)
 call spec_to_phys(n+2,hat_Tem,n,D2_Tem,x,2)
 do i=1, n
  b1(i)=2*epsilon*(z(i)*c1*Psi_T(i)-2*c2*Psi_P(i))+2*D1_Psi_P(i)+Psi_T(i)/dt &
       +0.5d0*Ek*(D2_Psi_T(i)-k2_perp*Psi_T(i))
  b2(i)=2*epsilon*c1*(Psi_T(i)+z(i)*(D2_Psi_P(i)-k2_perp*Psi_P(i)))-2*D1_Psi_T(i) &
       -(R_c+epsilon*delta_R)*Tem(i)+(D2_Psi_P(i)-k2_perp*Psi_P(i))/dt &
       +0.5d0*Ek*(D4_Psi_P(i)-2*k2_perp*D2_Psi_P(i)+k2_perp**2*Psi_P(i))
  b3(i)=2*epsilon*z(i)*c1*Tem(i)+k2_perp*Psi_P(i)+Tem(i)/dt &
       +0.5d0*Ek/Pr*(D2_Tem(i)-k2_perp*Tem(i))
 enddo
 ! multiply by inverse of coefficient matrices to update spectral coefficients
 call mat_mul(n+2,a1_inv,b1,hat_Psi_T)
 call mat_mul(n+4,a2_inv,b2,hat_Psi_P)
 call mat_mul(n+2,a3_inv,b3,hat_Tem)
 ! calculate spectral energy after update and output growth rate
 call energy(n+2,hat_Psi_T,energy1_1)
 call energy(n+4,hat_Psi_P,energy2_1)
 call energy(n+2,hat_Tem,energy3_1)
 sigma_u=log((energy1_1+energy2_1)/(energy1_0+energy2_0))/(2.d0*dt)
 sigma_T=log(energy3_1/energy3_0)/(2.d0*dt)
 write(1,'(5E15.6)') time, energy2_1, log(energy2_1/energy2_0)/(2.d0*dt), sigma_u, sigma_T
enddo

! output Psi_T, Psi_P, Tem
hat_Psi_T=hat_Psi_T/sqrt(energy1_1)
hat_Psi_P=hat_Psi_P/sqrt(energy2_1)
hat_Tem=hat_Tem/sqrt(energy3_1)
call spec_to_phys(n+2,hat_Psi_T,n,Psi_T,x,0)
call spec_to_phys(n+4,hat_Psi_P,n,Psi_P,x,0)
call spec_to_phys(n+2,hat_Tem,n,Tem,x,0)
do i=1, n
 write(2,'(4E15.6)') z(i), abs(Psi_T(i)), abs(Psi_P(i)), abs(Tem(i))
enddo
end program main

!!! use Chebyshev spectral coefficients to calculate k-th derivative in physical space at points x
subroutine spec_to_phys(ns,spec,np,phys,x,k)
implicit none
integer i, j, ns, np, k
double complex spec(ns), phys(np)
double precision x(np), TTT
do i=1, np
 phys(i)=0.d0
 do j=1, ns
  phys(i)=phys(i)+spec(j)*TTT(k,j-1,x(i))*2.d0**k
 enddo
enddo
end subroutine spec_to_phys

subroutine mat_inv(n,a,c)
implicit none
integer i, j, n
double precision a(n,n), b(n,n), c(n,n)
double precision WKSPCE(n), AA(n,n), BB(n,n)
integer IFAIL
do j=1, n
 do i=1, n
  b(i,j)=0.d0
 enddo
 b(j,j)=1.d0
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
