! gfortran -o precession precession.f90 nag.f
! Crank-Nicolson scheme, diffusion terms semi-implicit, the other terms explicit
! z=(1/2)x, D=d/dz=2d/dx, D^2=4d/dx, D^4=16d/dx
! Psi_T(x(i))=sum_j hat_Psi_T(j)*TTT(0,j-1,x(i)), j from 1 and j-1 from 0
! i=1~n: inner points; i=n+1, n+2, n+3 n+4: boundary conditions
! initial condition: two resonant inertial modes k_z=pi and 2pi, k_perp=5.748pi

program main
implicit none
double precision pi
double complex one
integer i, j, n							    ! i physical, j spectral, n dimension
parameter (n=50)
double precision x(n), z(n), TTT				    ! inner points
double precision Ek, Pr, force, R_c, delta_R			    ! dimensionless parameters
double precision k_x, k_y, k2_perp, k_z, k2			    ! wavenumber
double precision k_z_1, k_z_2, k2_1, k2_2, omega_1, omega_2         ! two resonance waves
double complex Psi_T(n),       Psi_P(n),       Tem(n)    	    ! coefficients about (z,t) -- physical space
double complex hat_Psi_T(n+2), hat_Psi_P(n+4), hat_Tem(n+2) 	    ! Chebyshev coefficients about t -- spectral space
double precision a1(n+2,n+2), a2(n+4,n+4), a3(n+2,n+2)  	    ! coefficient matrices 
double precision a1_inv(n+2,n+2), a2_inv(n+4,n+4), a3_inv(n+2,n+2)  ! inverse of coefficient matrices
integer ini							    ! select initial condition
double precision ini_r, ini_i					    ! random initial condition in physical space							
integer it, nt							    ! time steps	
double precision dt, time, c1, c2, c3				    ! c1 c2 c3 are coefficients of precession terms
double complex D1_Psi_T(n), D2_Psi_T(n)				    ! derivatives of Psi_T
double complex D1_Psi_P(n), D2_Psi_P(n), D4_Psi_P(n)		    ! derivatives of Psi_P
double complex D2_Tem(n)					    ! derivatives of Tem
double complex b1(n+2), b2(n+4), b3(n+2)  			    ! right-hand-side terms
double precision energy1_0, energy2_0, energy3_0		    ! spectral energy at the last timestep
double precision energy1_1, energy2_1, energy3_1		    ! spectral energy at the next timestep
integer ns							    ! mode number in Fourier space
parameter (ns=11)		
double complex ft(ns)						    ! modes in Fourier space

pi=acos(-1.d0)
one=(0.d0, 1.d0)
Ek=1.d-8
Pr=1.d-1
force=1.d-2
k_x=4.064639*pi	! resonance condition of two inertial modes k_z=pi and 2pi
k_y=4.064639*pi	! to satisfy omega_1-omega_2=1, k_perp can be solved =5.748*pi
k2_perp=k_x**2+k_y**2
k_z_1=pi
k_z_2=2.d0*pi
k2_1=k2_perp+k_z_1**2
k2_2=k2_perp+k_z_2**2
omega_1=2.d0*k_z_1/sqrt(k2_1)
omega_2=-2.d0*k_z_2/sqrt(k2_2)
k_z=k_z_1
k2=k2_perp+k_z**2
!R_c=8.d0*k_z**2/k2_perp*Pr/(1.d0+Pr)+2.d0*Ek**2*k2**3/k2_perp*(1.d0+Pr)/Pr
R_c=0.d0
!delta_R=10*R_c
delta_R=0.d0
ini=1
dt=1.d-1
nt=5000

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
if(ini.eq.0) then
 ! two inertial modes k_z=pi and 2pi
 do i=1, n
  Psi_P(i)=1.d-6*sin(k_z_1*(z(i)+0.5))+1.d-6*sin(k_z_2*(z(i)+0.5))
  Psi_T(i)=2.d-6*k_z_1/(one*omega_1+Ek*k2_1)*cos(k_z_1*(z(i)+0.5)) &
          +2.d-6*k_z_2/(one*omega_2+Ek*k2_2)*cos(k_z_2*(z(i)+0.5))
  Tem(i)  =2.d-6*k2_perp/(one*omega_1+Ek/Pr*k2_1)*sin(k_z_1*(z(i)+0.5)) &
          +2.d-6*k2_perp/(one*omega_2+Ek/Pr*k2_2)*sin(k_z_2*(z(i)+0.5))
 enddo  
else
 ! random
 call random_seed()
 do i=1,n
  call random_number(ini_r)
  call random_number(ini_i)
  Psi_P(i)=1.d-6*(ini_r+one*ini_i)
  call random_number(ini_r)
  call random_number(ini_i)
  Psi_T(i)=1.d-6*(ini_r+one*ini_i)
  call random_number(ini_r)
  call random_number(ini_i)
  Tem(i)=1.d-6*(ini_r+one*ini_i)
 enddo
endif
! transform to spectral space and calculate spectral energy
call phys_to_spec(n,Psi_P,n+4,hat_Psi_P,x)
call phys_to_spec(n,Psi_T,n+2,hat_Psi_T,x)
call phys_to_spec(n,Tem,n+2,hat_Tem,x)
call energy(n+2,hat_Psi_T,energy1_0)
call energy(n+4,hat_Psi_P,energy2_0)
call energy(n+2,hat_Tem,energy3_0)
! output initial condition in physical space
open(1,file='ini_phys.dat',form='formatted')
do i=1, n
 write(1,'(7E15.6)') z(i), real(Psi_T(i))/sqrt(energy1_0), imag(Psi_T(i))/sqrt(energy1_0), &
                           real(Psi_P(i))/sqrt(energy2_0), imag(Psi_P(i))/sqrt(energy2_0), &
                           real(Tem(i))/sqrt(energy3_0), imag(Tem(i))/sqrt(energy3_0)
enddo
close(1)
! output initial condition in Chebyshev spectral space
open(1,file='ini_tor.dat',form='formatted')
do j=1, n+2
 write(1,'(I5,3E15.6)') j, real(hat_Psi_T(j))/sqrt(energy1_0), imag(hat_Psi_T(j))/sqrt(energy1_0), &
                        abs(hat_Psi_T(j))**2/energy1_0
enddo
close(1)
open(1,file='ini_pol.dat',form='formatted')
do j=1, n+4
 write(1,'(I5,3E15.6)') j, real(hat_Psi_P(j))/sqrt(energy2_0), imag(hat_Psi_P(j))/sqrt(energy2_0), &
                        abs(hat_Psi_P(j))**2/energy2_0
enddo
close(1)
open(1,file='ini_tem.dat',form='formatted')
do j=1, n+2
 write(1,'(I5,3E15.6)') j, real(hat_Tem(j))/sqrt(energy3_0), imag(hat_Tem(j))/sqrt(energy3_0), &
                        abs(hat_Tem(j))**2/energy3_0
enddo
close(1)
!!! output initial condition in Fourier spectral space
open(1,file='ini_fourier.dat',form='formatted')
call Fourier(n,Psi_P,ns,ft)
do j=1, ns
 write(1,'(I10,3E15.6)') j-1, real(ft(j))/sqrt(energy2_0), imag(ft(j))/sqrt(energy2_0), &
                         abs(ft(j))**2/energy2_0
enddo
close(1)

write(6,'(2(A10,I10,/),13(A10,E15.6,/))') 'n=', n, 'ini=', ini, 'Ek=', Ek, 'Pr=', Pr, &
 'R_c=', R_c, 'delta_R', delta_R, 'force=', force, 'k_perp=', sqrt(k2_perp), &
 'k_z_1=', k_z_1, 'k_z_2=', k_z_2, 'k2_1=', k2_1, 'k2_2=', k2_2, &
 'omega_1=', omega_1, 'omega_2=', omega_2, 'dt=', dt

! time stepping
open(1,file='evolution.dat',form='formatted')
do it=1, nt
 time=it*dt
 c1=one*(k_x*sin(time)+k_y*cos(time))
 c2=one*(k_x*sin(time)-k_y*cos(time))
 c3=one*(k_x*cos(time)-k_y*sin(time))
 ! calculate spectral energy before update
 call energy(n+2,hat_Psi_T,energy1_0)
 call energy(n+4,hat_Psi_P,energy2_0)
 call energy(n+2,hat_Tem,energy3_0)
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
  b1(i)=2*force*(z(i)*c1*Psi_T(i)-2*c2*Psi_P(i))+2*D1_Psi_P(i)+Psi_T(i)/dt &
       +0.5d0*Ek*(D2_Psi_T(i)-k2_perp*Psi_T(i))
  b2(i)=2*force*(-z(i)*c1*D2_Psi_P(i)+2*c2*D1_Psi_P(i)-c3*Psi_T(i))-2*D1_Psi_T(i) &
       -(R_c+force*delta_R)*Tem(i)+(D2_Psi_P(i)-k2_perp*Psi_P(i))/dt &
       +0.5d0*Ek*(D4_Psi_P(i)-2*k2_perp*D2_Psi_P(i)+k2_perp**2*Psi_P(i))
  b3(i)=2*force*z(i)*c1*Tem(i)+k2_perp*Psi_P(i)+Tem(i)/dt &
       +0.5d0*Ek/Pr*(D2_Tem(i)-k2_perp*Tem(i))
 enddo
 ! multiply by inverse of coefficient matrices to update spectral coefficients
 call mat_mul(n+2,a1_inv,b1,hat_Psi_T)
 call mat_mul(n+4,a2_inv,b2,hat_Psi_P)
 call mat_mul(n+2,a3_inv,b3,hat_Tem)
 ! calculate spectral energy after update
 call energy(n+2,hat_Psi_T,energy1_1)
 call energy(n+4,hat_Psi_P,energy2_1)
 call energy(n+2,hat_Tem,energy3_1)
 ! output energy and growth rate
 write(1,'(7E15.6)') time, energy1_1, energy2_1, energy3_1, &
                     log(energy1_1/energy1_0)/(2.d0*dt), log(energy2_1/energy2_0)/(2.d0*dt), &
                     log(energy3_1/energy3_0)/(2.d0*dt)
enddo
close(1)

! output final result in physical space
open(1,file='fin_phys.dat',form='formatted')
call spec_to_phys(n+2,hat_Psi_T,n,Psi_T,x,0)
call spec_to_phys(n+4,hat_Psi_P,n,Psi_P,x,0)
call spec_to_phys(n+2,hat_Tem,n,Tem,x,0)
do i=1, n
 write(1,'(7E15.6)') z(i), real(Psi_T(i))/sqrt(energy1_1), imag(Psi_T(i))/sqrt(energy1_1), &
                           real(Psi_P(i))/sqrt(energy2_1), imag(Psi_P(i))/sqrt(energy2_1), &
                           real(Tem(i))/sqrt(energy3_1), imag(Tem(i))/sqrt(energy3_1)
enddo
close(1)
! output final result in Chebyshev spectral space
open(1,file='fin_tor.dat',form='formatted')
do j=1, n+2
 write(1,'(I5,3E15.6)') j, real(hat_Psi_T(j))/sqrt(energy1_1), imag(hat_Psi_T(j))/sqrt(energy1_1), &
                        abs(hat_Psi_T(j))**2/energy1_1
enddo
close(1)
open(1,file='fin_pol.dat',form='formatted')
do j=1, n+4
 write(1,'(I5,3E15.6)') j, real(hat_Psi_P(j))/sqrt(energy2_1), imag(hat_Psi_P(j))/sqrt(energy2_1), &
                        abs(hat_Psi_P(j))**2/energy2_1
enddo
close(1)
open(1,file='fin_tem.dat',form='formatted')
do j=1, n+2
 write(1,'(I5,3E15.6)') j, real(hat_Tem(j))/sqrt(energy3_1), imag(hat_Tem(j))/sqrt(energy3_1), &
                        abs(hat_Tem(j))**2/energy3_1
enddo
close(1)
!!! output final result in Fourier spectral space
open(1,file='fin_fourier.dat',form='formatted')
call Fourier(n,Psi_P,ns,ft)
do j=1, ns
 write(1,'(I10,3E15.6)') j-1, real(ft(j))/sqrt(energy2_1), imag(ft(j))/sqrt(energy2_1), &
                         abs(ft(j))**2/energy2_1
enddo
close(1)
end program main

!!! use Chebyshev spectral coefficients to calculate k-th derivative in physical space at points x
subroutine spec_to_phys(ns,spec,np,phys,x,k)
implicit none
integer i, j, ns, np, k
double complex spec(ns), phys(np)
double precision x(np), TTT
do i=1, np
 phys(i)=(0.d0, 0.d0)
 do j=1, ns
  phys(i)=phys(i)+spec(j)*TTT(k,j-1,x(i))*2.d0**k
 enddo
enddo
end subroutine spec_to_phys

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

subroutine Fourier(np,phys,ns,spec)
implicit none
double complex one
double precision pi
integer  i, j, np, ns
double complex phys(np), spec(ns)
one=(0.d0, 1.d0)
pi=acos(-1.d0)
if(ns>np) then
 write(6,*) "Warning: Fourier ns > np, will zero-pad in Fourier space"
endif
do j=1, ns
 spec(j)=(0.d0, 0.d0)
 do i=1, np
  spec(j)=spec(j)+phys(i)*exp(-one*2*pi*(j-1)*(i-1)/np)
 enddo
 spec(j)=spec(j)/np
enddo
end subroutine Fourier

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