program main
implicit none
double complex one
integer i, n						
parameter (n=50)
double complex Psi_T(n), Psi_P(n), Tem(n)
double precision z(n), x1, y1, x2, y2, x3, y3
integer j, ns
parameter (ns=n/2+1)		
double complex ft(ns)

one=(0.d0, 1.d0)

open(1,file='ini_phys.dat',form='formatted')
do i=1, n
 read(1,'(7E15.6)') z(i), x1, y1, x2, y2, x3, y3
 Psi_T(i)=x1+one*y1
 Psi_P(i)=x2+one*y2
 Tem(i)=x3+one*y3
enddo
close(1)

open(1,file='ini_fourier.dat',form='formatted')
call Fourier(n,Psi_P,ns,ft)
do j=1, ns
 write(1,'(I10,3E15.6)') j-1, real(ft(j)), imag(ft(j)), abs(ft(j))**2
enddo
close(1)

open(1,file='fin_phys.dat',form='formatted')
do i=1, n
 read(1,'(7E15.6)') z(i), x1, y1, x2, y2, x3, y3
 Psi_T(i)=x1+one*y1
 Psi_P(i)=x2+one*y2
 Tem(i)=x3+one*y3
enddo
close(1)

open(1,file='fin_fourier.dat',form='formatted')
call Fourier(n,Psi_P,ns,ft)
do j=1, ns
 write(1,'(I10,3E15.6)') j-1, real(ft(j)), imag(ft(j)), abs(ft(j))**2
enddo
close(1)
end program main

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